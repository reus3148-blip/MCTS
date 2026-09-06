"""Does the MCTS advantage come from confounding in the reward model? (v1.4)

This run connects the two halves of the project, which have never been joined.

The causal half (v0.6-v1.3) established that treatment assignment in METABRIC is
confounded by indication badly enough that a negative control fails. The
simulator half (v0.2-v1.2) rests on a Cox reward model **fitted to that same
observational data**, and MCTS optimises whatever that model says.

Reading the fitted coefficients out shows the problem directly:

===========================  ===================  ==========================
term                         reward model says    our adjusted analysis says
===========================  ===================  ==========================
chemotherapy                 HR 1.042 (harmful)   risk difference -0.023
endocrine therapy            HR 1.046 (harmful)   risk difference -0.048
mastectomy                   HR 1.093 (harmful)   (not estimated)
radiotherapy                 HR 0.939 (helpful)   risk difference -0.053
===========================  ===================  ==========================

Two of the three decisions carry the **wrong sign**. And in the static v0.1
comparison every systematic deviation MCTS makes from NCCN runs in the direction
those coefficients point - less chemotherapy (72.5% vs 83.1%), more radiotherapy
(100% vs 89.8%), more breast conservation (93.3% vs 76.1%).

That is the same defect v0.5 found in the response channel, one level up: an
**undeclared** clinical assumption baked into the environment, benefiting one
policy over the other. v0.5's remedy was to neutralise the undeclared channel and
let the declared assumptions in ``configs/dynamic_v0_5.json`` carry the effect,
where sensitivity analysis can reach it. This run applies the same remedy to the
reward model and measures what the headline does.

PRE-SPECIFIED PREDICTION, recorded before the run
--------------------------------------------------
1. **Primary** - with the reward model's treatment coefficients neutralised, the
   MCTS-minus-NCCN utility gap is **smaller** than with them left in. Tested as a
   seed-paired difference over the same twelve seeds.
2. The neutralised arm's gap stays **positive** - the declared config channels
   (timing, response, toxicity, recurrence) should still favour a searching
   policy somewhat.
3. MCTS's treatment rates move **toward** NCCN's when the coefficients are
   neutralised, because the pull that separated them is gone.

If prediction 1 fails, the advantage does not come from the reward model's
treatment terms and this concern is closed. If it holds, the headline has to be
restated as the neutralised number.
"""

from __future__ import annotations

from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.causal.decisions import (  # noqa: E402
    DEFAULT_SPEC,
    build_cohort,
    drop_constant_terms,
    trim_to_overlap,
)
from analysis.dynamic.cohort import (  # noqa: E402
    BASE_SEED,
    balanced_subtype_sample,
    build_reward_models,
    git_commit,
    input_manifest,
    make_risk_table,
)
from analysis.dynamic.config import DynamicConfig, _validate_probabilities  # noqa: E402
from analysis.dynamic.environment import DynamicBreastCancerEnvironment  # noqa: E402
from analysis.dynamic.evaluation import run_policy_episodes  # noqa: E402
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402
from analysis.mcts.environment import DECISIONS as PLAN_DECISIONS  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "reward-confounding-v1.4"
TABLE_DIR = REPORT_DIR / "tables"
CAUSAL_METRICS = ROOT / "reports" / "endocrine-effect-v1.3" / "metrics.json"
REPLICATION_METRICS = ROOT / "reports" / "cohort-replication-v1.2" / "metrics.json"

RUN_DATE = "2026-09-07"
PER_SUBTYPE = 5                 # the v1.2 cohorts A and B together: 40 patients
N_SEEDS = 12
SIMULATIONS = 1024
EPISODES_PER_POLICY = 40
EXPLORATION_WEIGHT = math.sqrt(2.0)
TRIM = (0.10, 0.90)

#: (key, neutralise, label)
ARMS = (
    ("as_fitted", False, "보상모형 그대로 (v0.2~v1.2)"),
    ("treatment_neutral", True, "치료 계수 중립화"),
)

PRESPECIFIED_PREDICTION = {
    "primary": (
        "Neutralising the reward model's treatment coefficients makes the "
        "MCTS-minus-NCCN utility gap smaller, tested seed-paired."
    ),
    "secondary_still_positive": (
        "The neutralised gap stays positive - the declared config channels "
        "should still favour a searching policy somewhat."
    ),
    "secondary_action_rates": (
        "MCTS's treatment rates move toward NCCN's once the coefficients that "
        "separated them are gone."
    ),
    "why": (
        "The reward model puts chemotherapy at HR 1.042 and endocrine therapy at "
        "HR 1.046 - both nominally harmful - while our own adjusted analyses of "
        "the same data put them on the protective side. MCTS optimises the model, "
        "so it can win by exploiting that sign error rather than by treating well."
    ),
}

ACTION_FIELDS = ("timing", "surgery", "chemo", "endocrine", "radiation")


def load_config() -> DynamicConfig:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        config = DynamicConfig(**json.load(handle))
    _validate_probabilities(config)
    return config


def cohort_of_forty(os_test: pd.DataFrame) -> pd.DataFrame:
    """v1.2's two disjoint cohorts, concatenated - the largest we have."""
    first = balanced_subtype_sample(os_test, PER_SUBTYPE)
    second = balanced_subtype_sample(os_test, PER_SUBTYPE, offset=PER_SUBTYPE)
    if set(first["patient_id"]) & set(second["patient_id"]):
        raise SystemExit("cohorts A and B are not disjoint")
    return pd.concat([first, second], ignore_index=True)


def evaluate(sample, os_model, rfs_model, config) -> dict:
    """Per-seed MCTS and NCCN utilities, plus the action mix each policy chose."""
    gaps, mcts_means, nccn_means = [], [], []
    actions = {"MCTS": [], "NCCN": []}
    per_patient = np.zeros((N_SEEDS, len(sample)), dtype=float)
    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        mcts_util, nccn_util = [], []
        for _, row in sample.iterrows():
            patient = patient_from_row(row)
            environment = DynamicBreastCancerEnvironment(
                patient, make_risk_table(row, os_model, rfs_model), config)
            offset = int(hashlib.sha256(
                patient.patient_id.encode()).hexdigest()[:8], 16) % 100_000
            patient_seed = seed + offset
            policy = CachedMCTSPolicy(
                environment, simulations=SIMULATIONS,
                exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed)
            md = pd.DataFrame(run_policy_episodes(
                environment, policy, EPISODES_PER_POLICY, patient_seed + 10_000))
            nd = pd.DataFrame(run_policy_episodes(
                environment, DynamicNccnPolicy(environment),
                EPISODES_PER_POLICY, patient_seed + 10_000))
            mcts_util.append(md["utility"].mean())
            nccn_util.append(nd["utility"].mean())
            per_patient[seed_index, len(mcts_util) - 1] = (
                md["utility"].mean() - nd["utility"].mean())
            actions["MCTS"].append(md)
            actions["NCCN"].append(nd)
        gaps.append(float(np.mean(mcts_util) - np.mean(nccn_util)))
        mcts_means.append(float(np.mean(mcts_util)))
        nccn_means.append(float(np.mean(nccn_util)))

    mix = {}
    for policy, frames in actions.items():
        joined = pd.concat(frames, ignore_index=True)
        mix[policy] = {
            field: joined[field].astype(str).value_counts(normalize=True).mul(100)
            .round(3).to_dict()
            for field in ACTION_FIELDS if field in joined.columns
        }
    return {
        "per_seed_gap": gaps,
        "per_patient_gap": per_patient.mean(axis=0).tolist(),
        "patient_ids": [str(value) for value in sample["patient_id"]],
        "subtypes": [str(value) for value in sample["subtype"]],
        "utility_gap": float(np.mean(gaps)),
        "utility_gap_sd": float(np.std(gaps, ddof=1)),
        "standard_error": float(np.std(gaps, ddof=1) / math.sqrt(N_SEEDS)),
        "mcts_utility": float(np.mean(mcts_means)),
        "nccn_utility": float(np.mean(nccn_means)),
        "action_mix": mix,
    }


def implied_risk_differences(raw: pd.DataFrame, os_model) -> pd.DataFrame:
    """What five-year risk difference the reward model implies per decision.

    Evaluated on the same propensity-overlap populations v1.3 estimated in, so
    the two numbers answer the same question about the same people. The model was
    fitted on data that includes them - that is the point: this measures what the
    model *encodes*, not how it generalises.
    """
    spec_base = DEFAULT_SPEC.with_surgery()
    arms = [
        ("chemo", None, "항암 (겹침 구간)"),
        ("hormone", ("er", 1.0), "호르몬 (ER 양성)"),
        ("hormone", ("er", 0.0), "호르몬 (ER 음성)"),
        ("radio", ("surgery", "MASTECTOMY"), "방사선 (전절제 후)"),
    ]
    rows = []
    for decision, restriction, label in arms:
        cohort = build_cohort(raw, decision, spec_base)
        if restriction is not None:
            cohort = cohort[cohort[restriction[0]] == restriction[1]].reset_index(
                drop=True)
        spec = drop_constant_terms(cohort, spec_base)
        trimmed = trim_to_overlap(cohort, decision, TRIM, spec)["cohort"]

        differences = []
        for _, row in trimmed.iterrows():
            surgery = ("MAST" if str(row.get("surgery", "")).upper()
                       .startswith("MAST") else "BCS")
            observed = {"surgery": surgery, "chemo": int(row["chemo"]),
                        "hormone": int(row["hormone"]), "radio": int(row["radio"])}
            plans = []
            for value in (1, 0):
                plan = dict(observed)
                plan[decision] = value
                plans.append(tuple(plan[name] for name in PLAN_DECISIONS))
            scores = os_model.score_plans(row, plans)
            differences.append(
                (1.0 - scores[plans[0]]) - (1.0 - scores[plans[1]]))
        rows.append({
            "decision": decision,
            "population": label,
            "n": int(len(trimmed)),
            "reward_model_risk_difference": float(np.mean(differences)),
        })
    return pd.DataFrame(rows)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)
    config = load_config()
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = cohort_of_forty(os_test)
    print(f"reward-model confounding: {len(sample)} patients, {N_SEEDS} seeds, "
          f"budget {SIMULATIONS}", flush=True)

    coefficients = os_model.coefficient_table()
    treatment_terms = os_model.treatment_features()
    coefficients["is_treatment_term"] = coefficients["feature"].isin(treatment_terms)
    coefficients.to_csv(TABLE_DIR / "reward_model_coefficients.csv", index=False)

    print("[implied effects] scoring plans on the v1.3 overlap populations...",
          flush=True)
    implied = implied_risk_differences(raw, os_model)
    causal = json.loads(CAUSAL_METRICS.read_text(encoding="utf-8"))
    causal_by_arm = {
        "항암 (겹침 구간)": causal["arms"]["chemo_with_surgery"],
        "호르몬 (ER 양성)": causal["arms"]["hormone_er_positive"],
        "호르몬 (ER 음성)": causal["arms"]["hormone_er_negative"],
        "방사선 (전절제 후)": causal["arms"]["radio_mastectomy"],
    }
    implied["causal_risk_difference"] = [
        causal_by_arm[label]["aipw_risk_difference"] for label in implied["population"]
    ]
    implied["sign_disagrees"] = (
        np.sign(implied["reward_model_risk_difference"])
        != np.sign(implied["causal_risk_difference"]))
    implied["gap"] = (implied["reward_model_risk_difference"]
                      - implied["causal_risk_difference"])
    implied.to_csv(TABLE_DIR / "implied_vs_causal.csv", index=False)
    print(implied.round(4).to_string(index=False), flush=True)

    results = {}
    started = time.perf_counter()
    for key, neutralise, label in ARMS:
        model = os_model.neutralise_treatment_terms() if neutralise else os_model
        print(f"\n[{key}] running...", flush=True)
        result = evaluate(sample, model, rfs_model, config)
        result["label"] = label
        result["neutralised"] = neutralise
        result["concordance"] = float(model.concordance(os_test))
        results[key] = result
        print(f"  gap {result['utility_gap']:+.4f} "
              f"(SE {result['standard_error']:.4f})  "
              f"MCTS {result['mcts_utility']:.4f} vs NCCN {result['nccn_utility']:.4f}  "
              f"({time.perf_counter() - started:.0f}s)", flush=True)

    fitted = np.array(results["as_fitted"]["per_seed_gap"])
    neutral = np.array(results["treatment_neutral"]["per_seed_gap"])
    paired = neutral - fitted
    paired_se = float(paired.std(ddof=1) / math.sqrt(N_SEEDS))
    verdict = {
        "gap_as_fitted": results["as_fitted"]["utility_gap"],
        "gap_treatment_neutral": results["treatment_neutral"]["utility_gap"],
        "paired_difference": float(paired.mean()),
        "paired_standard_error": paired_se,
        "paired_z": float(paired.mean() / paired_se) if paired_se > 0 else float("nan"),
        "primary_prediction_met": bool(paired.mean() < 0),
        "share_of_gap_from_reward_model": float(
            1.0 - results["treatment_neutral"]["utility_gap"]
            / results["as_fitted"]["utility_gap"])
        if results["as_fitted"]["utility_gap"] != 0 else float("nan"),
        "neutral_gap_still_positive": bool(
            results["treatment_neutral"]["utility_gap"] > 0),
        "concordance_as_fitted": results["as_fitted"]["concordance"],
        "concordance_neutral": results["treatment_neutral"]["concordance"],
    }

    # Subtype standardisation, the other correction on this same headline (v1.2).
    # Reported jointly so the two are not mistaken for alternatives: they are
    # independent and they compound.
    population = json.loads(
        (ROOT / "reports" / "cohort-replication-v1.2"
         / "metrics_posthoc_subtype.json").read_text(encoding="utf-8"))
    shares = population["population_reference"]["subtype_share"]
    standardised = {}
    for key, result in results.items():
        frame = pd.DataFrame({
            "patient_id": result["patient_ids"],
            "subtype": result["subtypes"],
            "gap": result["per_patient_gap"],
        })
        by_subtype = frame.groupby("subtype")["gap"].agg(["mean", "std", "count"])
        by_subtype["share"] = [shares.get(name, 0.0) for name in by_subtype.index]
        by_subtype["standard_error"] = (
            by_subtype["std"] / np.sqrt(by_subtype["count"]))
        standardised[key] = {
            "balanced_mean": float(by_subtype["mean"].mean()),
            "prevalence_standardised_mean": float(
                (by_subtype["share"] * by_subtype["mean"]).sum()),
            "prevalence_standardised_standard_error": float(np.sqrt(
                ((by_subtype["share"] * by_subtype["standard_error"]) ** 2).sum())),
            "by_subtype": by_subtype.reset_index().to_dict(orient="records"),
        }
    pd.DataFrame([
        {"arm": key, **row}
        for key, value in standardised.items()
        for row in value["by_subtype"]
    ]).to_csv(TABLE_DIR / "subtype_gaps.csv", index=False)
    verdict["standardised_gap_as_fitted"] = (
        standardised["as_fitted"]["prevalence_standardised_mean"])
    verdict["standardised_gap_treatment_neutral"] = (
        standardised["treatment_neutral"]["prevalence_standardised_mean"])
    verdict["both_corrections_applied"] = (
        standardised["treatment_neutral"]["prevalence_standardised_mean"])

    replication = json.loads(REPLICATION_METRICS.read_text(encoding="utf-8"))
    pooled = next(row for row in replication["by_cohort"]
                  if row["cohort"] == "pooled")
    verdict["v1_2_pooled_gap"] = pooled["baseline_utility_gap"]
    verdict["reproduces_v1_2"] = bool(abs(
        results["as_fitted"]["utility_gap"] - pooled["baseline_utility_gap"]) < 1e-9)

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "reward-confounding-v1.4",
        "question": (
            "MCTS optimises a Cox reward model fitted to confounded observational "
            "data. How much of the MCTS-minus-NCCN gap comes from that model's "
            "treatment coefficients rather than from the declared assumptions?"
        ),
        "estimand": (
            "Difference between the MCTS-minus-NCCN utility gap with the reward "
            "model's treatment coefficients as fitted and with them set to zero, "
            "over the same 40 patients and the same 12 seeds."
        ),
        "scope_warning": (
            "Sensitivity of our own simulator to an undeclared channel in its "
            "reward model. Not a clinical effect."
        ),
        "prespecified_prediction": PRESPECIFIED_PREDICTION,
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "simulations": SIMULATIONS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "config": "configs/dynamic_v0_5.json",
            "neutralised_terms": treatment_terms,
        },
        "verdict": verdict,
        "subtype_standardisation": standardised,
        "reward_model_vs_causal": implied.to_dict(orient="records"),
        "arms": results,
    }

    pd.DataFrame([
        {"arm": key, "seed_index": index, "utility_gap": value}
        for key, result in results.items()
        for index, value in enumerate(result["per_seed_gap"])
    ]).to_csv(TABLE_DIR / "per_seed_gaps.csv", index=False)
    pd.DataFrame([
        {"arm": key, "policy": policy, "field": field, "action": action,
         "pct": pct}
        for key, result in results.items()
        for policy, fields in result["action_mix"].items()
        for field, counts in fields.items()
        for action, pct in counts.items()
    ]).to_csv(TABLE_DIR / "action_mix.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "assumptions": CONFIG_PATH}),
        "entry_point": "analysis/33_run_reward_model_confounding.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== verdict ===")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
