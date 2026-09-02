"""Doubly robust estimation and an identifiability map across decisions (v0.7).

v0.6 estimated one decision with inverse-probability weighting and hit two limits
it could not answer from inside itself.

**Was the estimate leaning on a single model?** IPW is consistent only if the
treatment model is right. Part A adds two more estimators over the same overlap
population - g-computation, which leans entirely on an *outcome* model, and AIPW,
which is consistent if *either* model is right - plus inverse-probability-of-
censoring weights so that dependent censoring is handled rather than assumed away.
Three estimators that disagree would mean model misspecification is driving the
answer; three that agree is the strongest statement this data supports.

**Was chemotherapy unusually hard, or is every decision like this?** Part B runs
the positivity diagnostics of v0.6 across all three decisions METABRIC records -
chemotherapy, endocrine therapy, radiotherapy - and reports which of them
observational data can answer at all. The result is a map of where causal
questions are askable, which is what should drive the K-CURE variable request.

Still a point-treatment emulation, still not the MCTS-vs-NCCN regime comparison.
"""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.causal.decisions import (  # noqa: E402
    DECISIONS,
    DEFAULT_SPEC,
    CovariateSpec,
    balance_table,
    build_cohort,
    design_matrix,
    propensity_and_weights,
    trim_to_overlap,
)
from analysis.causal.ipw import (  # noqa: E402
    aipw_risk,
    discrete_time_rows,
    e_value,
    e_value_for_interval,
    effective_sample_size,
    fit_logistic,
    sigmoid,
    survival_probability,
    weighted_kaplan_meier,
)
from analysis.dynamic.cohort import (  # noqa: E402
    BASE_SEED,
    git_commit,
    input_manifest,
)

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
REPORT_DIR = ROOT / "reports" / "doubly-robust-v0.7"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-27"
PRIMARY_TREATMENT = "chemo"
HORIZON_MONTHS = 60.0
INTERVAL_EDGES = np.array([0.0, 12.0, 24.0, 36.0, 48.0, 60.0])
PRIMARY_TRIM = (0.10, 0.90)     # the trim v0.6 selected for this decision
BOOTSTRAP = 500


def horizon_status(cohort: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """``(outcome, observed)`` at the horizon.

    A patient's five-year status is known if they died within the horizon, or if
    follow-up reached it. Anyone censored earlier contributes no outcome, which
    is what the censoring weights exist to correct for.
    """
    months = cohort["os_months"].to_numpy(dtype=float)
    died = cohort["os_event"].to_numpy(dtype=float) == 1
    outcome = ((months <= HORIZON_MONTHS) & died).astype(float)
    observed = ((months >= HORIZON_MONTHS) | ((months < HORIZON_MONTHS) & died))
    return outcome, observed.astype(float)


def pooled_design(
    cohort: pd.DataFrame,
    treatment: str,
    rows: np.ndarray,
    intervals: np.ndarray,
    force_treatment: float | None = None,
) -> np.ndarray:
    """Person-interval design: interval dummies, treatment, baseline covariates."""
    base = design_matrix(cohort)[rows]
    n_intervals = INTERVAL_EDGES.size - 1
    dummies = np.zeros((rows.size, n_intervals - 1))
    for index in range(1, n_intervals):     # interval 0 is the reference
        dummies[intervals == index, index - 1] = 1.0
    if force_treatment is None:
        exposure = cohort[treatment].to_numpy(dtype=float)[rows]
    else:
        exposure = np.full(rows.size, float(force_treatment))
    return np.column_stack([base, dummies, exposure])


def censoring_weights(cohort: pd.DataFrame, treatment: str) -> np.ndarray:
    """Probability of reaching one's horizon status uncensored, given X and A.

    Estimated as a discrete-time censoring hazard (pooled logistic over yearly
    intervals) rather than assumed constant, so that censoring may depend on the
    same covariates that drive treatment.
    """
    months = cohort["os_months"].to_numpy(dtype=float)
    died = cohort["os_event"].to_numpy(dtype=float) == 1
    censored_early = ((months < HORIZON_MONTHS) & ~died).astype(float)
    rows, intervals, censor_here = discrete_time_rows(
        months, censored_early, INTERVAL_EDGES)
    design = pooled_design(cohort, treatment, rows, intervals)
    hazard = sigmoid(design @ fit_logistic(design, censor_here))
    return np.clip(survival_probability(hazard, rows, len(cohort)), 1e-3, 1.0)


def outcome_predictions(
    cohort: pd.DataFrame,
    treatment: str,
    outcome: np.ndarray,
    observed: np.ndarray,
    uncensored: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """IPCW-weighted outcome model, evaluated with treatment set to 1 and to 0.

    Fitting only on patients whose status is known would over-represent those with
    long follow-up; weighting them by 1/P(uncensored) restores the population they
    stand in for.
    """
    design = np.column_stack([
        design_matrix(cohort), cohort[treatment].to_numpy(dtype=float)])
    weights = observed / uncensored
    beta = fit_logistic(design, outcome, sample_weight=weights)

    def predict(value: float) -> np.ndarray:
        forced = np.column_stack([
            design_matrix(cohort), np.full(len(cohort), value)])
        return sigmoid(forced @ beta)

    return predict(1.0), predict(0.0)


def ipw_risks(cohort: pd.DataFrame, treatment: str, weights: np.ndarray) -> dict:
    """v0.6's estimator: IPT-weighted Kaplan-Meier risk at the horizon."""
    def arm(value: int) -> float:
        mask = (cohort[treatment] == value).to_numpy()
        curve = weighted_kaplan_meier(
            cohort.loc[mask, "os_months"].to_numpy(dtype=float),
            cohort.loc[mask, "os_event"].to_numpy(dtype=float),
            weights[mask])
        return 1.0 - curve.at(HORIZON_MONTHS)

    treated_risk, control_risk = arm(1), arm(0)
    return {
        "risk_treated": treated_risk,
        "risk_control": control_risk,
        "risk_difference": treated_risk - control_risk,
        "risk_ratio": treated_risk / control_risk if control_risk > 0 else float("nan"),
    }


def all_estimators(cohort: pd.DataFrame, treatment: str, weights: np.ndarray,
                   propensity: np.ndarray) -> dict[str, dict]:
    """IPW, g-computation and AIPW over the same population."""
    outcome, observed = horizon_status(cohort)
    uncensored = censoring_weights(cohort, treatment)
    predicted_treated, predicted_control = outcome_predictions(
        cohort, treatment, outcome, observed, uncensored)

    g_treated = float(predicted_treated.mean())
    g_control = float(predicted_control.mean())
    return {
        "ipw_km": ipw_risks(cohort, treatment, weights),
        "g_computation": {
            "risk_treated": g_treated,
            "risk_control": g_control,
            "risk_difference": g_treated - g_control,
            "risk_ratio": g_treated / g_control if g_control > 0 else float("nan"),
        },
        "aipw": aipw_risk(
            cohort[treatment].to_numpy(dtype=float), outcome, observed,
            propensity, uncensored, predicted_treated, predicted_control),
        "censoring": {
            "min_uncensored_probability": float(uncensored.min()),
            "mean_uncensored_probability": float(uncensored.mean()),
            "observed_at_horizon_pct": float(observed.mean() * 100),
        },
    }


def bootstrap_aipw(cohort: pd.DataFrame, treatment: str, seed: int) -> dict:
    """Percentile CI repeating fit, trim, refit and both nuisance models."""
    rng = np.random.default_rng(seed)
    n = len(cohort)
    differences, ratios = [], []
    for _ in range(BOOTSTRAP):
        sample = cohort.iloc[rng.integers(0, n, size=n)].reset_index(drop=True)
        try:
            trimmed = trim_to_overlap(sample, treatment, PRIMARY_TRIM)
            result = all_estimators(
                trimmed["cohort"], treatment, trimmed["weights"],
                trimmed["propensity"])["aipw"]
        except (np.linalg.LinAlgError, ValueError):
            continue
        if np.isfinite(result["risk_difference"]):
            differences.append(result["risk_difference"])
        if np.isfinite(result["risk_ratio"]) and result["risk_ratio"] > 0:
            ratios.append(result["risk_ratio"])
    return {
        "risk_difference": [float(np.quantile(differences, q)) for q in (0.025, 0.975)],
        "risk_ratio": [float(np.quantile(ratios, q)) for q in (0.025, 0.975)],
        "replicates": len(differences),
    }


def overlap_diagnostics(
    raw: pd.DataFrame,
    treatment: str,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> dict:
    """Can observational data answer this decision at all, under this spec?"""
    cohort = build_cohort(raw, treatment, spec)
    treated = cohort[treatment].to_numpy(dtype=float)
    propensity, weights = propensity_and_weights(cohort, treatment, spec)
    untrimmed = balance_table(cohort, treatment, weights, spec)

    row = {
        "decision": treatment,
        "label": DECISIONS[treatment],
        "covariate_spec": spec.label,
        "eligible": int(len(cohort)),
        "treated": int(treated.sum()),
        "treated_pct": float(treated.mean() * 100),
        "median_propensity_control": float(np.median(propensity[treated == 0])),
        "median_propensity_treated": float(np.median(propensity[treated == 1])),
        "propensity_separation": float(
            np.median(propensity[treated == 1]) - np.median(propensity[treated == 0])),
        "untrimmed_worst_abs_smd": float(untrimmed["smd_weighted"].abs().max()),
        "untrimmed_effective_sample_size": effective_sample_size(weights),
    }
    try:
        trimmed = trim_to_overlap(cohort, treatment, PRIMARY_TRIM, spec)
        row.update({
            "trimmed_n": trimmed["n"],
            "retained_pct": trimmed["retained_pct"],
            "trimmed_worst_abs_smd": trimmed["worst_abs_smd"],
            "trimmed_balanced_pct": trimmed["balanced_pct"],
            "trimmed_effective_sample_size": trimmed["effective_sample_size"],
            "identifiable": bool(trimmed["balanced_pct"] == 100.0),
        })
    except ValueError:
        row.update({
            "trimmed_n": 0, "retained_pct": 0.0,
            "trimmed_worst_abs_smd": float("nan"), "trimmed_balanced_pct": 0.0,
            "trimmed_effective_sample_size": 0.0, "identifiable": False,
        })
    return row


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)

    # --- Part A: three estimators over the v0.6 overlap population -----------
    cohort = build_cohort(raw, PRIMARY_TREATMENT)
    trimmed = trim_to_overlap(cohort, PRIMARY_TREATMENT, PRIMARY_TRIM)
    print(f"primary decision = {PRIMARY_TREATMENT}: "
          f"{trimmed['n']} of {len(cohort)} in overlap", flush=True)

    estimators = all_estimators(
        trimmed["cohort"], PRIMARY_TREATMENT, trimmed["weights"],
        trimmed["propensity"])
    for name in ("ipw_km", "g_computation", "aipw"):
        result = estimators[name]
        print(f"  {name:14s} RD {result['risk_difference']:+.4f}  "
              f"RR {result['risk_ratio']:.3f}", flush=True)

    intervals = bootstrap_aipw(cohort, PRIMARY_TREATMENT, BASE_SEED)
    ratio_ci = intervals["risk_ratio"]
    print(f"  AIPW 95% CI [{intervals['risk_difference'][0]:+.4f}, "
          f"{intervals['risk_difference'][1]:+.4f}] "
          f"({intervals['replicates']} replicates)", flush=True)

    spread = max(
        abs(estimators[a]["risk_difference"] - estimators[b]["risk_difference"])
        for a in ("ipw_km", "g_computation", "aipw")
        for b in ("ipw_km", "g_computation", "aipw")
    )

    # --- Part B: which decisions are answerable at all? ----------------------
    # Run under two confounder sets. Radiotherapy is chosen largely by surgery
    # type, which the baseline set omits; if wide overlap survives adding it the
    # decision really is close to unconfounded, and if it collapses the overlap
    # was measuring our ignorance rather than clinical equipoise.
    with_surgery = DEFAULT_SPEC.with_surgery()
    overlap_rows = [
        overlap_diagnostics(raw, decision, spec)
        for spec in (DEFAULT_SPEC, with_surgery)
        for decision in DECISIONS
    ]
    overlap = pd.DataFrame(overlap_rows)
    print("\n=== identifiability map ===")
    print(overlap[[
        "covariate_spec", "decision", "eligible", "treated_pct",
        "propensity_separation", "retained_pct", "trimmed_worst_abs_smd",
        "identifiable",
    ]].round(3).to_string(index=False))

    def find(decision: str, label: str) -> dict:
        return next(row for row in overlap_rows
                    if row["decision"] == decision and row["covariate_spec"] == label)

    radio_shift = {
        "retained_pct_baseline": find("radio", "baseline")["retained_pct"],
        "retained_pct_with_surgery": find("radio", "baseline + surgery")["retained_pct"],
        "worst_abs_smd_baseline": find("radio", "baseline")["trimmed_worst_abs_smd"],
        "worst_abs_smd_with_surgery":
            find("radio", "baseline + surgery")["trimmed_worst_abs_smd"],
        "identifiable_baseline": find("radio", "baseline")["identifiable"],
        "identifiable_with_surgery":
            find("radio", "baseline + surgery")["identifiable"],
    }

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "doubly-robust-v0.7",
        "question": (
            "Do IPW, g-computation and AIPW agree on the chemotherapy effect, and "
            "which treatment decisions can this observational data answer at all?"
        ),
        "estimand": (
            "Five-year all-cause mortality risk difference at the horizon, treated "
            "vs untreated, in the propensity-overlap population "
            f"(propensity in {list(PRIMARY_TRIM)}) of each decision."
        ),
        "scope_warning": (
            "Point-treatment emulation of single decisions. Not the MCTS-vs-NCCN "
            "regime comparison, which needs treatment timing METABRIC lacks."
        ),
        "assumptions": {
            "exchangeability": "conditional on the measured baseline covariates only",
            "positivity": "enforced by trimming to the overlap region, not assumed",
            "censoring": "modelled as a covariate- and treatment-dependent "
                         "discrete-time hazard (IPCW), not assumed independent",
            "consistency": "one version of each treatment; METABRIC records no "
                           "regimen or dose",
        },
        "design": {
            "primary_decision": PRIMARY_TREATMENT,
            "eligible_patients": int(len(cohort)),
            "analysed_patients": trimmed["n"],
            "retained_pct": trimmed["retained_pct"],
            "primary_trim": list(PRIMARY_TRIM),
            "horizon_months": HORIZON_MONTHS,
            "interval_edges_months": INTERVAL_EDGES.tolist(),
            "bootstrap_replicates": intervals["replicates"],
        },
        "estimators": {
            name: estimators[name] for name in ("ipw_km", "g_computation", "aipw")
        },
        "censoring": estimators["censoring"],
        "estimator_spread_risk_difference": spread,
        "aipw_ci95": {
            "risk_difference": intervals["risk_difference"],
            "risk_ratio": ratio_ci,
        },
        "e_value_point": e_value(estimators["aipw"]["risk_ratio"]),
        "e_value_ci_limit": e_value_for_interval(*ratio_ci),
        "identifiability_map": overlap_rows,
        "decisions_identifiable_baseline": [
            row["decision"] for row in overlap_rows
            if row["identifiable"] and row["covariate_spec"] == "baseline"
        ],
        "decisions_identifiable_with_surgery": [
            row["decision"] for row in overlap_rows
            if row["identifiable"] and row["covariate_spec"] == "baseline + surgery"
        ],
        "radiotherapy_overlap_was_missing_confounder": radio_shift,
    }

    pd.DataFrame([
        {"estimator": name, **estimators[name]}
        for name in ("ipw_km", "g_computation", "aipw")
    ]).to_csv(TABLE_DIR / "estimator_comparison.csv", index=False)
    overlap.to_csv(TABLE_DIR / "identifiability_map.csv", index=False)
    outcome, observed = horizon_status(trimmed["cohort"])
    pd.DataFrame({
        "patient_id": trimmed["cohort"]["patient_id"],
        PRIMARY_TREATMENT: trimmed["cohort"][PRIMARY_TREATMENT],
        "outcome_5y_death": outcome,
        "status_observed": observed,
        "uncensored_probability": censoring_weights(
            trimmed["cohort"], PRIMARY_TREATMENT),
    }).to_csv(TABLE_DIR / "censoring_weights.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV}),
        "protocol": "docs/target-trial-protocol.md",
        "entry_point": "analysis/17_run_doubly_robust.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(f"\nestimator spread (max |RD difference|) = {spread:.4f}")
    print(f"E-value: point {metrics['e_value_point']:.2f}, "
          f"CI limit {metrics['e_value_ci_limit']:.2f}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
