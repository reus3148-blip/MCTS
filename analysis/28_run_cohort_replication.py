"""Out-of-sample replication of the sensitivity ranking (v1.2).

v1.1 raised the cohort from 8 to 20 patients and the ranking held, but it closed
by naming its own weakness: the three cohorts were **nested**, so the 20-patient
result partly re-counted the same eight people. Agreement between nested cohorts
is not a replication.

This run draws a **disjoint** second cohort of 20 (``offset=5`` in
``balanced_subtype_sample``, the same draw order, the next five patients per
subtype) and asks the only question that settles it: does the ranking found in
cohort A appear in patients cohort A never saw?

PRE-SPECIFIED PREDICTION, recorded before the run
-------------------------------------------------
From v1.1's cohort A (20 patients, seeds 12, budget 1024):

1. **Primary** - in cohort B the top two parameters will be the two value
   judgements (``reward.recurrence_free_year`` and
   ``reward.acute_toxicity_penalty``).
2. The baseline utility gap in cohort B will stay positive.
3. Spearman rank correlation between A and B will be positive.

Stated as a prediction, so a miss is a miss. v1.1 already showed the patient
bootstrap puts prediction 1 at only 69.5%, so this is not a formality: roughly
one draw in three should fail it.

Cohort A is re-run rather than read from v1.1's tables, for two reasons: the
per-(seed, patient) matrices are needed for the pooled analysis and are not
stored, and re-running it is a second reproducibility check. Pooling A and B
gives a 40-patient estimate at no extra search cost.

Still synthetic assumptions - this measures our simulator, not a clinical effect.
"""

from __future__ import annotations

import copy
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
from analysis.dynamic.experiment_utils import (  # noqa: E402
    rank_parameters_by_influence,
    rescale_response_major,
)
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "cohort-replication-v1.2"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "sensitivity-patients-v1.1" / "metrics.json"

RUN_DATE = "2026-09-04"
PER_SUBTYPE = 5                # 20 patients per cohort
N_SEEDS = 12                   # fixed at the v1.0 / v1.1 setting
SIMULATIONS = 1024             # fixed at the v1.0 / v1.1 setting
EPISODES_PER_POLICY = 40
EXPLORATION_WEIGHT = math.sqrt(2.0)

#: (label, offset). Offset 0 is v1.1's cohort; offset 5 shares no patient with it.
COHORTS = (("A", 0), ("B", PER_SUBTYPE))

VALUE_JUDGEMENTS = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}

#: Written before the run; echoed into metrics.json so the record shows it was a
#: prediction rather than a description of what came out.
PRESPECIFIED_PREDICTION = {
    "primary": (
        "In cohort B the top two parameters by |gap delta| are the two value "
        "judgements."
    ),
    "secondary_baseline_positive": "Cohort B's baseline utility gap is positive.",
    "secondary_rank_correlation": (
        "Spearman rank correlation of the six-parameter ranking between cohort A "
        "and cohort B is positive."
    ),
    "expected_failure_rate": (
        "v1.1's patient bootstrap put the primary prediction at 69.5%, so about "
        "one draw in three is expected to miss it."
    ),
}


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def set_response_major(data: dict, intensity: str, major: float) -> None:
    data["response_probabilities"][intensity] = rescale_response_major(
        data["response_probabilities"][intensity], major)


def build_variants(base: dict) -> list[dict]:
    """The v0.3 perturbation set, unchanged so every ranking lines up."""
    variants = [{"name": "baseline", "param": "(none)", "value": "-",
                 "data": copy.deepcopy(base)}]

    def add(name, param, value, mutate):
        data = copy.deepcopy(base)
        mutate(data)
        variants.append({"name": name, "param": param, "value": value, "data": data})

    for v in (0.20, 0.50):
        add(f"intensified_major={v}", "response.intensified.major", v,
            lambda d, v=v: set_response_major(d, "intensified", v))
    for v in (0.85, 1.00):
        add(f"intensified_death_hr={v}", "hazard.chemo.intensified.death", v,
            lambda d, v=v: d["hazard_multipliers"]["chemo"]["intensified"]
            .__setitem__("death", v))
    for v in (0.15, 0.45):
        add(f"intensified_toxicity={v}", "toxicity.chemo.intensified", v,
            lambda d, v=v: d["acute_toxicity_probabilities"]["chemo"]
            .__setitem__("intensified", v))
    for v in (0.05, 0.30):
        add(f"toxicity_penalty={v}", "reward.acute_toxicity_penalty", v,
            lambda d, v=v: d["reward"].__setitem__("acute_toxicity_penalty", v))
    for v in (0.90, 1.00):
        add(f"major_response_death_hr={v}", "hazard.response.major.death", v,
            lambda d, v=v: d["hazard_multipliers"]["response"]["major"]
            .__setitem__("death", v))
    for v in (0.25, 1.00):
        add(f"recurrence_free_reward={v}", "reward.recurrence_free_year", v,
            lambda d, v=v: d["reward"].__setitem__("recurrence_free_year", v))
    return variants


def evaluate(config: DynamicConfig, sample, os_model, rfs_model) -> np.ndarray:
    """``(seed, patient)`` matrix of MCTS-minus-NCCN utility gaps."""
    gaps = np.empty((N_SEEDS, len(sample)), dtype=float)
    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        for patient_index, (_, row) in enumerate(sample.iterrows()):
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
            gaps[seed_index, patient_index] = (
                md["utility"].mean() - nd["utility"].mean())
    return gaps


def influence_table(variants: list[dict],
                    matrices: dict[str, np.ndarray]) -> pd.DataFrame:
    """Largest |Δ vs baseline| per parameter, with the seed-paired SE."""
    baseline = matrices["baseline"].mean(axis=1)
    records = []
    for variant in variants:
        if variant["name"] == "baseline":
            continue
        delta = matrices[variant["name"]].mean(axis=1) - baseline
        stderr = float(delta.std(ddof=1)) / math.sqrt(N_SEEDS)
        records.append({
            "parameter": variant["param"],
            "variant": variant["name"],
            "gap_delta": float(delta.mean()),
            "abs_gap_delta": float(abs(delta.mean())),
            "delta_standard_error": stderr,
            "delta_z": float(delta.mean() / stderr) if stderr > 0 else float("nan"),
        })
    frame = pd.DataFrame(records)
    best = frame.loc[frame.groupby("parameter")["abs_gap_delta"].idxmax()]
    return best.sort_values("abs_gap_delta", ascending=False).reset_index(drop=True)


def pooled_matrices(a: dict[str, np.ndarray],
                    b: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Concatenate the two cohorts along the patient axis.

    Both were run at the same seeds, so column-stacking gives the 40-patient
    result without any extra search.
    """
    return {name: np.concatenate([a[name], b[name]], axis=1) for name in a}


def spearman(order_a: list[str], order_b: list[str]) -> float:
    ranks_a = {name: index for index, name in enumerate(order_a)}
    return float(pd.Series([ranks_a[name] for name in order_b]).corr(
        pd.Series(range(len(order_b))), method="spearman"))


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    variants = build_variants(base)

    samples = {label: balanced_subtype_sample(os_test, PER_SUBTYPE, offset=offset)
               for label, offset in COHORTS}
    ids = {label: set(sample["patient_id"]) for label, sample in samples.items()}
    overlap = ids["A"] & ids["B"]
    if overlap:
        raise SystemExit(f"cohorts are not disjoint: {sorted(overlap)}")
    for label, sample in samples.items():
        counts = sample["subtype"].value_counts()
        if len(sample) != PER_SUBTYPE * 4 or not (counts == PER_SUBTYPE).all():
            raise SystemExit(f"cohort {label} is not subtype-balanced: "
                             f"{counts.to_dict()}")
    print(f"cohort replication: A={len(samples['A'])} B={len(samples['B'])} "
          f"disjoint, {N_SEEDS} seeds, budget {SIMULATIONS}, "
          f"{len(variants)} variants", flush=True)

    detail_rows: list[dict] = []
    per_patient_rows: list[dict] = []
    matrices_by_cohort: dict[str, dict[str, np.ndarray]] = {}
    started = time.perf_counter()
    for label, _ in COHORTS:
        sample = samples[label]
        matrices: dict[str, np.ndarray] = {}
        for variant in variants:
            gaps = evaluate(config_from_dict(variant["data"]), sample,
                            os_model, rfs_model)
            matrices[variant["name"]] = gaps
            per_seed = gaps.mean(axis=1)
            detail_rows.append({
                "cohort": label,
                "variant": variant["name"],
                "parameter": variant["param"],
                "value": variant["value"],
                "utility_gap": float(per_seed.mean()),
                "utility_gap_sd": float(per_seed.std(ddof=1)),
                "standard_error": float(per_seed.std(ddof=1) / math.sqrt(N_SEEDS)),
            })
            print(f"  [{label}] {variant['name']:30s} "
                  f"gap={per_seed.mean():+.4f} "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)
        matrices_by_cohort[label] = matrices
        for patient_id, subtype, gap in zip(
                sample["patient_id"], sample["subtype"],
                matrices["baseline"].mean(axis=0)):
            per_patient_rows.append({
                "cohort": label,
                "patient_id": str(patient_id),
                "subtype": str(subtype),
                "baseline_utility_gap": float(gap),
            })

    matrices_by_cohort["pooled"] = pooled_matrices(
        matrices_by_cohort["A"], matrices_by_cohort["B"])
    for variant in variants:
        per_seed = matrices_by_cohort["pooled"][variant["name"]].mean(axis=1)
        detail_rows.append({
            "cohort": "pooled",
            "variant": variant["name"],
            "parameter": variant["param"],
            "value": variant["value"],
            "utility_gap": float(per_seed.mean()),
            "utility_gap_sd": float(per_seed.std(ddof=1)),
            "standard_error": float(per_seed.std(ddof=1) / math.sqrt(N_SEEDS)),
        })

    detail = pd.DataFrame(detail_rows)
    per_patient = pd.DataFrame(per_patient_rows)
    influence_by_cohort = {}
    for label in ("A", "B", "pooled"):
        frame = influence_table(variants, matrices_by_cohort[label])
        frame.insert(0, "cohort", label)
        influence_by_cohort[label] = frame
    influence = pd.concat(influence_by_cohort.values(), ignore_index=True)

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    prior_largest = prior["by_patient_count"][-1]

    summary = []
    for label in ("A", "B", "pooled"):
        frame = influence_by_cohort[label]
        order = list(frame["parameter"])
        block = detail[detail.cohort == label]
        baseline_row = block[block.variant == "baseline"].iloc[0]
        n_patients = matrices_by_cohort[label]["baseline"].shape[1]
        summary.append({
            "cohort": label,
            "n_patients": int(n_patients),
            "baseline_utility_gap": float(baseline_row["utility_gap"]),
            "baseline_seed_standard_error": float(baseline_row["standard_error"]),
            "top_two_are_value_judgements": bool(set(order[:2]) == VALUE_JUDGEMENTS),
            "ranking": order,
            "parameters_distinguishable_from_zero": int(
                (frame["delta_z"].abs() >= 2).sum()),
            "negative_variant_cells": int((block["utility_gap"] < 0).sum()),
            "variant_cells": int(len(block)),
        })
    by_cohort = {row["cohort"]: row for row in summary}

    order_a = by_cohort["A"]["ranking"]
    order_b = by_cohort["B"]["ranking"]
    replication = {
        "primary_prediction_met": by_cohort["B"]["top_two_are_value_judgements"],
        "baseline_positive_in_b": by_cohort["B"]["baseline_utility_gap"] > 0,
        "spearman_a_vs_b": spearman(order_a, order_b),
        "rank_correlation_positive": spearman(order_a, order_b) > 0,
        "identical_ranking": order_a == order_b,
        "baseline_gap_difference": (by_cohort["B"]["baseline_utility_gap"]
                                    - by_cohort["A"]["baseline_utility_gap"]),
        "cohort_a_matches_v1_1": bool(abs(
            by_cohort["A"]["baseline_utility_gap"]
            - prior_largest["baseline_utility_gap"]) < 1e-9),
        "v1_1_cohort_a_baseline": prior_largest["baseline_utility_gap"],
    }

    per_cohort_sd = {
        label: float(per_patient[per_patient.cohort == label]
                     ["baseline_utility_gap"].std(ddof=1))
        for label in ("A", "B")
    }
    pooled_patient_gaps = per_patient["baseline_utility_gap"]
    heterogeneity = {
        "per_patient_gap_sd": {**per_cohort_sd,
                               "pooled": float(pooled_patient_gaps.std(ddof=1))},
        "pooled_min": float(pooled_patient_gaps.min()),
        "pooled_max": float(pooled_patient_gaps.max()),
        "pooled_median": float(pooled_patient_gaps.median()),
        "patients_with_negative_gap": int((pooled_patient_gaps < 0).sum()),
        "n_patients": int(len(pooled_patient_gaps)),
        "between_patient_standard_error_pooled": float(
            pooled_patient_gaps.std(ddof=1) / math.sqrt(len(pooled_patient_gaps))),
    }

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "cohort-replication-v1.2",
        "question": (
            "v1.1's cohorts were nested, so its agreement could be a recount "
            "rather than a replication. Does the sensitivity ranking appear in "
            "twenty patients the first cohort never saw?"
        ),
        "estimand": (
            "Largest absolute shift in the MCTS-minus-NCCN utility gap that each "
            "assumption produces when moved to either end of its v0.3 range, "
            "estimated separately in two disjoint subtype-balanced cohorts of 20."
        ),
        "scope_warning": (
            "Sensitivity of our own simulator to its declared assumptions. Not a "
            "clinical effect."
        ),
        "prespecified_prediction": PRESPECIFIED_PREDICTION,
        "design": {
            "cohorts": {label: {"offset": offset, "n_patients": len(samples[label])}
                        for label, offset in COHORTS},
            "cohorts_disjoint": True,
            "seeds": N_SEEDS,
            "simulations": SIMULATIONS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "variants": len(variants),
            "config": "configs/dynamic_v0_5.json",
        },
        "replication": replication,
        "by_cohort": summary,
        "patient_heterogeneity": heterogeneity,
    }

    detail.to_csv(TABLE_DIR / "variant_results.csv", index=False)
    influence.to_csv(TABLE_DIR / "parameter_influence.csv", index=False)
    per_patient.to_csv(TABLE_DIR / "per_patient_gaps.csv", index=False)
    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "assumptions": CONFIG_PATH}),
        "entry_point": "analysis/28_run_cohort_replication.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== pre-specified prediction ===")
    print(json.dumps(replication, indent=2))
    print("\n=== ranking by cohort ===")
    for label in ("A", "B", "pooled"):
        print(f"\n[{label}]")
        print(influence_by_cohort[label][[
            "parameter", "abs_gap_delta", "delta_standard_error", "delta_z",
        ]].round(4).to_string(index=False))
    print("\n=== summary ===")
    print(pd.DataFrame(summary)[[
        "cohort", "n_patients", "baseline_utility_gap",
        "top_two_are_value_judgements", "parameters_distinguishable_from_zero",
        "negative_variant_cells",
    ]].to_string(index=False))
    print("\n=== patient heterogeneity ===")
    print(json.dumps(heterogeneity, indent=2))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
