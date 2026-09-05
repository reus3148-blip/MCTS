"""Endocrine therapy effect, with a negative control (v1.3).

v0.7 drew a map of which decisions this data can answer and then estimated only
one of them. Endocrine therapy came out the most answerable of the three -
89.0% of patients retained in the overlap region, worst |SMD| 0.056, effective
sample 825 - and was left unestimated. This run finishes it.

Endocrine therapy is also the one decision here with unambiguous randomised
evidence behind it: EBCTCG's meta-analyses put five years of tamoxifen at about
a one-third reduction in breast-cancer mortality in ER-positive disease, and at
essentially nothing in ER-negative disease. That gives this analysis something
none of our other estimates have - **a known answer to check against**.

So the design is a positive and a negative control on the same pipeline:

* **Primary (positive control)** - ER-positive patients, where the drug works.
* **Negative control** - ER-negative patients, where it does not. Any clear
  protective estimate here cannot be the drug; it is confounding by indication,
  and it would bound how much of the primary estimate we can believe.

A third arm re-estimates chemotherapy under the surgery-augmented covariate set.
v0.7 estimated chemotherapy under the baseline set and then, in the same report,
showed that omitting surgery type distorts the overlap diagnostics. Its own
finding therefore applies to its own estimate, and this checks whether it moved.

PRE-SPECIFIED PREDICTIONS, recorded before the run
--------------------------------------------------
1. **Primary** - the endocrine risk difference in ER-positive patients is
   negative (protective), consistent with EBCTCG.
2. **Negative control** - the ER-negative risk difference has a 95% interval
   covering zero. If it is clearly protective, we report the primary estimate as
   contaminated by confounding rather than as an effect.
3. Chemotherapy under the surgery-augmented spec stays inside v0.7's interval
   [-0.1354, +0.0630].

Point-treatment emulation of single decisions, as in v0.6/v0.7. Not the
MCTS-vs-NCCN regime comparison, which still needs treatment timing METABRIC lacks.
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
    DEFAULT_SPEC,
    balance_table,
    build_cohort,
    drop_constant_terms,
    propensity_and_weights,
    trim_to_overlap,
)
from analysis.causal.effects import (  # noqa: E402
    HORIZON_MONTHS,
    INTERVAL_EDGES,
    all_estimators,
    bootstrap_aipw,
    estimator_spread,
)
from analysis.causal.ipw import (  # noqa: E402
    e_value,
    e_value_for_interval,
    effective_sample_size,
)
from analysis.dynamic.cohort import BASE_SEED, git_commit, input_manifest  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
REPORT_DIR = ROOT / "reports" / "endocrine-effect-v1.3"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "doubly-robust-v0.7" / "metrics.json"

RUN_DATE = "2026-09-05"
TRIM = (0.10, 0.90)             # the same overlap window v0.6 selected
BOOTSTRAP = 500

PRESPECIFIED_PREDICTION = {
    "primary": (
        "The endocrine-therapy risk difference in ER-positive patients is "
        "negative (protective), consistent with EBCTCG."
    ),
    "negative_control": (
        "The ER-negative risk difference has a 95% interval covering zero. A "
        "clearly protective estimate there is confounding by indication, not the "
        "drug, and bounds how much of the primary estimate is believable."
    ),
    "chemo_respecification": (
        "Chemotherapy under the surgery-augmented covariate set stays inside "
        "v0.7's interval [-0.1354, +0.0630]."
    ),
    "external_benchmark": (
        "EBCTCG (Lancet 2011): five years of tamoxifen reduces breast-cancer "
        "mortality by about a third in ER-positive disease and is without effect "
        "in ER-negative disease."
    ),
}

#: (key, treatment, restriction, label, prespecified). ``restriction`` is applied
#: to the eligible cohort before any model is fitted.
#:
#: The two radiotherapy arms were added *after* the trimming fix, which is what
#: made the decision identifiable at all - they are marked post-hoc for that
#: reason. They are split by surgery type for the same reason the endocrine arm
#: is split by ER status: post-lumpectomy radiotherapy is standard of care rather
#: than a decision, and pooling it with the post-mastectomy decision mixes the
#: indication into the contrast.
ARMS = (
    ("hormone_er_positive", "hormone", ("er", 1.0), "호르몬치료 (ER 양성)", True),
    ("hormone_er_negative", "hormone", ("er", 0.0), "호르몬치료 (ER 음성, 음성대조)", True),
    ("chemo_with_surgery", "chemo", None, "보조 항암치료 (수술 유형 보정)", True),
    ("radio_mastectomy", "radio", ("surgery", "MASTECTOMY"),
     "방사선치료 (전절제 후)", False),
    ("radio_conserving", "radio", ("surgery", "BREAST CONSERVING"),
     "방사선치료 (유방보존술 후)", False),
)


def restrict(cohort: pd.DataFrame, restriction) -> pd.DataFrame:
    if restriction is None:
        return cohort
    column, value = restriction
    return cohort[cohort[column] == value].reset_index(drop=True)


def analyse(raw: pd.DataFrame, treatment: str, restriction, seed: int) -> dict:
    """Overlap diagnostics plus three estimators and a bootstrap interval.

    The covariate spec is fixed to baseline + surgery for every arm: v0.7 showed
    surgery type matters, and using the same set everywhere is what makes the
    three arms comparable to each other.
    """
    spec = DEFAULT_SPEC.with_surgery()
    eligible = restrict(build_cohort(raw, treatment, spec), restriction)
    spec = drop_constant_terms(eligible, spec)

    propensity, weights = propensity_and_weights(eligible, treatment, spec)
    untrimmed = balance_table(eligible, treatment, weights, spec)
    trimmed = trim_to_overlap(eligible, treatment, TRIM, spec)
    if not trimmed["converged"]:
        raise SystemExit(
            f"{treatment}: trim did not reach a fixed point in "
            f"{trimmed['iterations']} iterations")

    estimators = all_estimators(
        trimmed["cohort"], treatment, trimmed["weights"], trimmed["propensity"],
        spec, HORIZON_MONTHS, INTERVAL_EDGES)
    intervals = bootstrap_aipw(
        eligible, treatment, seed, spec, TRIM, BOOTSTRAP,
        HORIZON_MONTHS, INTERVAL_EDGES)

    aipw = estimators["aipw"]
    low, high = intervals["risk_difference"]
    treated = eligible[treatment].to_numpy(dtype=float)
    return {
        "covariate_spec": spec.label,
        "covariates": spec.names,
        "eligible": int(len(eligible)),
        "treated": int(treated.sum()),
        "treated_pct": float(treated.mean() * 100),
        "median_propensity_control": float(np.median(propensity[treated == 0])),
        "median_propensity_treated": float(np.median(propensity[treated == 1])),
        "untrimmed_worst_abs_smd": float(untrimmed["smd_weighted"].abs().max()),
        "untrimmed_effective_sample_size": effective_sample_size(weights),
        "analysed": trimmed["n"],
        "retained_pct": trimmed["retained_pct"],
        "trim_iterations": trimmed["iterations"],
        "trim_converged": trimmed["converged"],
        "max_weight": trimmed["max_weight"],
        "trimmed_worst_abs_smd": trimmed["worst_abs_smd"],
        "trimmed_balanced_pct": trimmed["balanced_pct"],
        "trimmed_effective_sample_size": trimmed["effective_sample_size"],
        "identifiable": bool(trimmed["balanced_pct"] == 100.0),
        "estimators": {name: estimators[name]
                       for name in ("ipw_km", "g_computation", "aipw")},
        "censoring": estimators["censoring"],
        "estimator_spread_risk_difference": estimator_spread(estimators),
        "aipw_risk_difference": aipw["risk_difference"],
        "aipw_risk_ratio": aipw["risk_ratio"],
        "aipw_ci95": [low, high],
        "aipw_risk_ratio_ci95": intervals["risk_ratio"],
        "bootstrap_replicates": intervals["replicates"],
        "interval_covers_null": bool(low <= 0.0 <= high),
        "e_value_point": e_value(aipw["risk_ratio"]),
        "e_value_ci_limit": e_value_for_interval(*intervals["risk_ratio"]),
        "balance": trimmed["balance"],
    }


def trim_convergence(raw: pd.DataFrame) -> pd.DataFrame:
    """Single-pass vs iterated trimming, for every decision and both specs.

    This is the table that shows the defect rather than asserting it: for each
    cell, what one trim-then-refit pass leaves behind and what the fixed point
    leaves behind. Cheap - no bootstrap, no outcome model.
    """
    rows = []
    cells = [(decision, spec, restriction, label)
             for spec in (DEFAULT_SPEC, DEFAULT_SPEC.with_surgery())
             for decision, restriction, label in (
                 ("chemo", None, "전체"),
                 ("hormone", None, "전체"),
                 ("hormone", ("er", 1.0), "ER 양성"),
                 ("radio", None, "전체"),
                 ("radio", ("surgery", "MASTECTOMY"), "전절제 후"),
             )]
    for decision, spec0, restriction, label in cells:
        eligible = restrict(build_cohort(raw, decision, spec0), restriction)
        spec = drop_constant_terms(eligible, spec0)
        for passes, name in ((1, "single_pass"), (25, "iterated")):
            try:
                trimmed = trim_to_overlap(eligible, decision, TRIM, spec,
                                          max_iterations=passes)
            except ValueError:
                continue
            rows.append({
                "decision": decision,
                "population": label,
                "covariate_spec": spec0.label,
                "trimming": name,
                "iterations": trimmed["iterations"],
                "converged": trimmed["converged"],
                "eligible": int(len(eligible)),
                "n": trimmed["n"],
                "retained_pct": trimmed["retained_pct"],
                "worst_abs_smd": trimmed["worst_abs_smd"],
                "balanced_pct": trimmed["balanced_pct"],
                "effective_sample_size": trimmed["effective_sample_size"],
                "max_weight": trimmed["max_weight"],
                "propensity_inside_bounds_pct":
                    trimmed["propensity_inside_bounds_pct"],
                "identifiable": bool(trimmed["balanced_pct"] == 100.0),
            })
    return pd.DataFrame(rows)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)

    results: dict[str, dict] = {}
    balance_rows: list[pd.DataFrame] = []
    for index, (key, treatment, restriction, label, prespecified) in enumerate(ARMS):
        print(f"[{key}] running...", flush=True)
        result = analyse(raw, treatment, restriction, BASE_SEED + index * 17)
        balance = result.pop("balance")
        balance.insert(0, "arm", key)
        balance_rows.append(balance)
        result["label"] = label
        result["treatment"] = treatment
        result["prespecified"] = prespecified
        results[key] = result
        low, high = result["aipw_ci95"]
        print(f"  n {result['analysed']}/{result['eligible']} "
              f"({result['retained_pct']:.1f}%)  "
              f"AIPW RD {result['aipw_risk_difference']:+.4f} "
              f"[{low:+.4f}, {high:+.4f}]", flush=True)

    print("\n[trim convergence] single pass vs fixed point...", flush=True)
    convergence = trim_convergence(raw)
    convergence.to_csv(TABLE_DIR / "trim_convergence.csv", index=False)

    def cell(decision: str, population: str, spec_label: str, trimming: str) -> dict:
        match = convergence[
            (convergence.decision == decision)
            & (convergence.population == population)
            & (convergence.covariate_spec == spec_label)
            & (convergence.trimming == trimming)]
        return match.iloc[0].to_dict()

    trimming_effect = {
        "decisions_needing_more_than_one_pass": int(
            (convergence[convergence.trimming == "iterated"]["iterations"] > 1).sum()),
        "decision_cells": int((convergence.trimming == "iterated").sum()),
        "single_pass_not_at_fixed_point": int(
            (~convergence[convergence.trimming == "single_pass"]["converged"]).sum()),
        "endocrine_er_positive": {
            "single_pass": cell("hormone", "ER 양성", "baseline + surgery", "single_pass"),
            "iterated": cell("hormone", "ER 양성", "baseline + surgery", "iterated"),
        },
        "radiotherapy_with_surgery": {
            "single_pass": cell("radio", "전체", "baseline + surgery", "single_pass"),
            "iterated": cell("radio", "전체", "baseline + surgery", "iterated"),
        },
    }

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    prior_ci = [-0.1354, 0.0630]

    positive = results["hormone_er_positive"]
    negative = results["hormone_er_negative"]
    chemo = results["chemo_with_surgery"]
    verdict = {
        "primary_prediction_met": bool(positive["aipw_risk_difference"] < 0),
        "primary_interval_excludes_null": bool(not positive["interval_covers_null"]),
        "negative_control_prediction_met": bool(negative["interval_covers_null"]),
        "negative_control_risk_difference": negative["aipw_risk_difference"],
        "negative_control_ci95": negative["aipw_ci95"],
        "chemo_inside_v0_7_interval": bool(
            prior_ci[0] <= chemo["aipw_risk_difference"] <= prior_ci[1]),
        "chemo_v0_7_risk_difference": prior["estimators"]["aipw"]["risk_difference"],
        "chemo_shift_from_v0_7": (chemo["aipw_risk_difference"]
                                  - prior["estimators"]["aipw"]["risk_difference"]),
        "v0_7_interval": prior_ci,
    }

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "endocrine-effect-v1.3",
        "question": (
            "v0.7 found endocrine therapy the most answerable decision in this "
            "data and did not estimate it. What is the effect, and does the "
            "pipeline pass a negative control where the drug is known not to work?"
        ),
        "estimand": (
            "Five-year all-cause mortality risk difference at the horizon, "
            "treated vs untreated, in the propensity-overlap population "
            "(propensity in [0.10, 0.90]) of each arm."
        ),
        "scope_warning": (
            "Point-treatment emulation of single decisions. Not the MCTS-vs-NCCN "
            "regime comparison, which needs treatment timing METABRIC lacks."
        ),
        "assumptions": {
            "exchangeability": "conditional on the measured baseline covariates plus surgery type",
            "positivity": "enforced by trimming to the overlap region, not assumed",
            "censoring": "modelled as a covariate- and treatment-dependent discrete-time hazard (IPCW)",
            "consistency": "one version of each treatment; METABRIC records no agent, dose or duration",
        },
        "prespecified_prediction": PRESPECIFIED_PREDICTION,
        "design": {
            "horizon_months": HORIZON_MONTHS,
            "interval_edges_months": INTERVAL_EDGES.tolist(),
            "trim": list(TRIM),
            "bootstrap_replicates": BOOTSTRAP,
            "covariate_spec": "baseline + surgery, constant terms dropped per arm",
            "arms": [{"key": key, "treatment": treatment,
                      "restriction": None if r is None else {r[0]: r[1]},
                      "label": label, "prespecified": prespecified}
                     for key, treatment, r, label, prespecified in ARMS],
        },
        "verdict": verdict,
        "trimming_effect": trimming_effect,
        "arms": results,
    }

    frame = pd.DataFrame([
        {
            "arm": key,
            "label": result["label"],
            "treatment": result["treatment"],
            "prespecified": result["prespecified"],
            "eligible": result["eligible"],
            "analysed": result["analysed"],
            "retained_pct": result["retained_pct"],
            "treated_pct": result["treated_pct"],
            "trimmed_worst_abs_smd": result["trimmed_worst_abs_smd"],
            "identifiable": result["identifiable"],
            "trim_iterations": result["trim_iterations"],
            "effective_sample_size": result["trimmed_effective_sample_size"],
            "ipw_km": result["estimators"]["ipw_km"]["risk_difference"],
            "g_computation": result["estimators"]["g_computation"]["risk_difference"],
            "aipw": result["aipw_risk_difference"],
            "aipw_ci_low": result["aipw_ci95"][0],
            "aipw_ci_high": result["aipw_ci95"][1],
            "risk_ratio": result["aipw_risk_ratio"],
            "estimator_spread": result["estimator_spread_risk_difference"],
            "e_value_point": result["e_value_point"],
            "e_value_ci_limit": result["e_value_ci_limit"],
        }
        for key, result in results.items()
    ])
    frame.to_csv(TABLE_DIR / "arm_results.csv", index=False)
    pd.concat(balance_rows, ignore_index=True).to_csv(
        TABLE_DIR / "covariate_balance.csv", index=False)
    pd.DataFrame([
        {"arm": key, "estimator": name, **result["estimators"][name]}
        for key, result in results.items()
        for name in ("ipw_km", "g_computation", "aipw")
    ]).to_csv(TABLE_DIR / "estimator_comparison.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV}),
        "entry_point": "analysis/31_run_endocrine_effect.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== pre-specified predictions ===")
    print(json.dumps(verdict, indent=2))
    print("\n=== arms ===")
    print(frame[[
        "arm", "eligible", "analysed", "retained_pct", "trimmed_worst_abs_smd",
        "aipw", "aipw_ci_low", "aipw_ci_high", "estimator_spread",
    ]].round(4).to_string(index=False))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
