"""Target trial emulation, step one: IPW for a single treatment decision (v0.6).

Every number this project has published so far is an *association* that a Cox
model learned from observational data. `docs/target-trial-protocol.md` specifies
what it would take to say "causes" instead. That protocol's full estimand - one
dynamic treatment regime versus another - needs g-methods over time-varying
confounding, which METABRIC cannot support because it records no treatment
timing.

So this script does the piece that is possible now: a **point-treatment**
emulation of one decision node, adjuvant chemotherapy versus none, on five-year
overall survival. It is the machinery test - propensity model, positivity,
balance, weighted survival, sensitivity to unmeasured confounding - that the
sequential version will be built on.

**This is not the MCTS-versus-NCCN causal comparison.** It answers a much smaller
question, and saying so plainly is the point.

The headline difficulty is positivity, not confounding. In the full eligible
cohort chemotherapy is close to deterministic given age, receptor status and
nodes: the median propensity is about 0.04 among controls and 0.69 among the
treated. No amount of weighting fixes that, so the analysis is restricted to the
overlap region where the decision was genuinely uncertain - which happens to be
the population the rest of this project cares about.
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
    balance_table,
    build_cohort,
    covariate_frame,
    propensity_and_weights,
    trim_to_overlap,
)
from analysis.causal.ipw import (  # noqa: E402
    e_value,
    e_value_for_interval,
    effective_sample_size,
    weighted_kaplan_meier,
)
from analysis.dynamic.cohort import (  # noqa: E402
    BASE_SEED,
    git_commit,
    input_manifest,
)

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
REPORT_DIR = ROOT / "reports" / "ipw-target-trial-v0.6"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-27"
TREATMENT = "chemo"
HORIZON_MONTHS = 60.0
BOOTSTRAP = 1000

#: Trims considered, loosest first. The primary analysis is the loosest one that
#: balances every covariate; the others are reported as sensitivity so the
#: dependence of the estimate on this choice is visible rather than buried.
TRIM_GRID = ((0.05, 0.95), (0.10, 0.90), (0.15, 0.85))


def risk_at_horizon(cohort, weights: np.ndarray, treated_value: int) -> float:
    """Weighted five-year risk of death in one arm."""
    mask = (cohort[TREATMENT] == treated_value).to_numpy()
    curve = weighted_kaplan_meier(
        cohort.loc[mask, "os_months"].to_numpy(dtype=float),
        cohort.loc[mask, "os_event"].to_numpy(dtype=float),
        weights[mask],
    )
    return 1.0 - curve.at(HORIZON_MONTHS)


def estimate(cohort: pd.DataFrame, weights: np.ndarray) -> dict[str, float]:
    treated_risk = risk_at_horizon(cohort, weights, 1)
    control_risk = risk_at_horizon(cohort, weights, 0)
    return {
        "risk_treated": treated_risk,
        "risk_control": control_risk,
        "risk_difference": treated_risk - control_risk,
        "risk_ratio": treated_risk / control_risk if control_risk > 0 else float("nan"),
    }


def trimmed_analysis(cohort: pd.DataFrame, bounds: tuple[float, float]) -> dict:
    """Overlap trimming plus the naive and IPW estimates for this decision."""
    result = trim_to_overlap(cohort, TREATMENT, bounds)
    trimmed, weights = result["cohort"], result["weights"]
    result["naive"] = estimate(trimmed, np.ones(len(trimmed)))
    result["ipw"] = estimate(trimmed, weights)
    return result


def bootstrap_interval(
    cohort: pd.DataFrame,
    bounds: tuple[float, float],
    seed: int,
) -> dict:
    """Percentile CI repeating the whole pipeline - fit, trim, refit - per replicate.

    Treating the propensity model or the trim as fixed would ignore where most of
    the uncertainty lives and give intervals that are far too narrow.
    """
    rng = np.random.default_rng(seed)
    n = len(cohort)
    differences, ratios = [], []
    for _ in range(BOOTSTRAP):
        sample = cohort.iloc[rng.integers(0, n, size=n)].reset_index(drop=True)
        try:
            result = trimmed_analysis(sample, bounds)
        except (np.linalg.LinAlgError, ValueError):
            continue
        difference = result["ipw"]["risk_difference"]
        ratio = result["ipw"]["risk_ratio"]
        if np.isfinite(difference):
            differences.append(difference)
        if np.isfinite(ratio) and ratio > 0:
            ratios.append(ratio)
    return {
        "risk_difference": [float(np.quantile(differences, q)) for q in (0.025, 0.975)],
        "risk_ratio": [float(np.quantile(ratios, q)) for q in (0.025, 0.975)],
        "replicates": len(differences),
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)
    cohort = build_cohort(raw, TREATMENT)
    treated = cohort[TREATMENT].to_numpy(dtype=float)
    print(f"eligible = {len(cohort)} "
          f"(treated {int(treated.sum())}, control {int(len(cohort) - treated.sum())})",
          flush=True)

    # --- positivity in the full cohort: the finding that shapes everything ----
    full_propensity, full_weights = propensity_and_weights(cohort, TREATMENT)
    full_balance = balance_table(cohort, TREATMENT, full_weights)
    positivity = {
        "median_propensity_treated": float(np.median(full_propensity[treated == 1])),
        "median_propensity_control": float(np.median(full_propensity[treated == 0])),
        "propensity_min": float(full_propensity.min()),
        "propensity_max": float(full_propensity.max()),
        "untrimmed_worst_abs_smd": float(full_balance["smd_weighted"].abs().max()),
        "untrimmed_balanced_pct": float(full_balance["balanced_after"].mean() * 100),
        "untrimmed_effective_sample_size": effective_sample_size(full_weights),
        "untrimmed_max_weight": float(full_weights.max()),
    }
    print(f"full cohort: median propensity {positivity['median_propensity_control']:.3f} "
          f"(control) vs {positivity['median_propensity_treated']:.3f} (treated); "
          f"worst |SMD| after weighting {positivity['untrimmed_worst_abs_smd']:.3f}",
          flush=True)

    # --- trim grid -----------------------------------------------------------
    analyses = []
    for bounds in TRIM_GRID:
        result = trimmed_analysis(cohort, bounds)
        analyses.append(result)
        print(f"  trim {bounds}: n={result['n']:4d} "
              f"worst |SMD|={result['worst_abs_smd']:.3f} "
              f"balanced={result['balanced_pct']:.0f}% "
              f"RD={result['ipw']['risk_difference']:+.3f}", flush=True)

    balanced = [a for a in analyses if a["balanced_pct"] == 100.0]
    if not balanced:
        raise RuntimeError("no trim in TRIM_GRID balances every covariate")
    primary = balanced[0]          # loosest trim that balances everything
    bounds = tuple(primary["bounds"])
    print(f"primary trim = {bounds}", flush=True)

    intervals = bootstrap_interval(cohort, bounds, BASE_SEED)
    ratio_ci = intervals["risk_ratio"]

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "ipw-target-trial-v0.6",
        "question": (
            "Point-treatment emulation: what is the five-year mortality effect of "
            "adjuvant chemotherapy versus none under IPW adjustment for measured "
            "baseline confounders?"
        ),
        "scope_warning": (
            "ONE decision node, not the MCTS-vs-NCCN regime comparison the target "
            "trial protocol specifies. METABRIC records no treatment timing, so the "
            "sequential estimand is not identifiable here."
        ),
        "estimand": (
            "Five-year all-cause mortality risk difference, chemotherapy vs none, "
            f"in the propensity-overlap population (propensity in {list(bounds)})."
        ),
        "design": {
            "eligible_patients": int(len(cohort)),
            "eligible_treated": int(treated.sum()),
            "analysed_patients": primary["n"],
            "analysed_treated": primary["n_treated"],
            "retained_pct": primary["retained_pct"],
            "horizon_months": HORIZON_MONTHS,
            "confounders": list(covariate_frame(cohort).columns),
            "primary_trim": list(bounds),
            "bootstrap_replicates": intervals["replicates"],
        },
        "positivity": positivity,
        "balance": {
            "worst_abs_smd_crude": float(primary["balance"]["smd_crude"].abs().max()),
            "worst_abs_smd_weighted": primary["worst_abs_smd"],
            "covariates_balanced_after_pct": primary["balanced_pct"],
            "effective_sample_size": primary["effective_sample_size"],
            "max_weight": primary["max_weight"],
        },
        "naive": primary["naive"],
        "ipw": primary["ipw"],
        "ipw_ci95": {
            "risk_difference": intervals["risk_difference"],
            "risk_ratio": ratio_ci,
        },
        "e_value_point": e_value(primary["ipw"]["risk_ratio"]),
        "e_value_ci_limit": e_value_for_interval(*ratio_ci),
        "confounding_by_indication_shift": (
            primary["ipw"]["risk_difference"] - primary["naive"]["risk_difference"]
        ),
        "trim_sensitivity": [
            {
                "bounds": a["bounds"], "n": a["n"], "n_treated": a["n_treated"],
                "worst_abs_smd": a["worst_abs_smd"], "balanced_pct": a["balanced_pct"],
                "effective_sample_size": a["effective_sample_size"],
                "risk_difference": a["ipw"]["risk_difference"],
                "risk_ratio": a["ipw"]["risk_ratio"],
            }
            for a in analyses
        ],
    }

    primary["balance"].to_csv(TABLE_DIR / "covariate_balance.csv", index=False)
    full_balance.to_csv(TABLE_DIR / "covariate_balance_untrimmed.csv", index=False)
    pd.DataFrame([
        {"estimator": "naive (unadjusted, trimmed)", **primary["naive"]},
        {"estimator": "IPW (stabilized, trimmed)", **primary["ipw"]},
    ]).to_csv(TABLE_DIR / "effect_estimates.csv", index=False)
    pd.DataFrame(metrics["trim_sensitivity"]).to_csv(
        TABLE_DIR / "trim_sensitivity.csv", index=False)
    pd.DataFrame({
        "patient_id": cohort["patient_id"],
        TREATMENT: cohort[TREATMENT],
        "propensity_full_cohort": full_propensity,
        "in_primary_overlap": (
            (full_propensity >= bounds[0]) & (full_propensity <= bounds[1])
        ),
    }).to_csv(TABLE_DIR / "propensity_full_cohort.csv", index=False)
    pd.DataFrame({
        "patient_id": primary["cohort"]["patient_id"],
        TREATMENT: primary["cohort"][TREATMENT],
        "propensity": primary["propensity"],
        "stabilized_weight": primary["weights"],
    }).to_csv(TABLE_DIR / "propensity_and_weights_trimmed.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV}),
        "protocol": "docs/target-trial-protocol.md",
        "entry_point": "analysis/15_run_ipw_target_trial.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== primary analysis (overlap population) ===")
    print(f"n = {primary['n']} of {len(cohort)} "
          f"({primary['retained_pct']:.0f}% retained), "
          f"worst |SMD| {primary['balance']['smd_crude'].abs().max():.3f} -> "
          f"{primary['worst_abs_smd']:.3f}")
    for label, result in (("naive", primary["naive"]), ("IPW", primary["ipw"])):
        print(f"{label:>5}: treated {result['risk_treated']:.3f}  "
              f"control {result['risk_control']:.3f}  "
              f"RD {result['risk_difference']:+.3f}  "
              f"RR {result['risk_ratio']:.3f}")
    rd_ci = intervals["risk_difference"]
    print(f"IPW risk difference 95% CI [{rd_ci[0]:+.3f}, {rd_ci[1]:+.3f}] "
          f"({intervals['replicates']} replicates)")
    print(f"E-value: point {metrics['e_value_point']:.2f}, "
          f"CI limit {metrics['e_value_ci_limit']:.2f}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
