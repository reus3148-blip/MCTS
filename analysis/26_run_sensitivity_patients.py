"""One-at-a-time sensitivity as the cohort grows (v1.1).

v1.0 fixed the two knobs it could reach: twelve seeds and a 1024 budget. It
closed by naming what was left - the eight simulated patients - and that is what
this run varies. Everything else is held at the v1.0 setting so any movement in
the ranking has one candidate explanation.

Three cohort sizes are run rather than two. Two points can only say "it moved";
three can say whether it is still moving. The 8-patient arm also re-runs the
exact v1.0 configuration, so it doubles as a reproducibility check - it must
return v1.0's numbers to the last digit or something in the pipeline has drifted.

``balanced_subtype_sample`` extends rather than replaces as ``per_subtype``
grows, so the 20-patient cohort contains the 8. That makes the comparison
nested: the extra patients are added to the same people, not swapped for
different ones.

Two things beyond the ranking come out of the larger cohort:

* the **between-patient spread** of the utility gap, which is what actually
  determines how many patients this comparison needs. Seed noise and search
  noise were the previous bottlenecks; this is the third.
* a **patient bootstrap** of the ranking. The seed-paired standard error in v1.0
  answered "would other seeds rank these differently?". Resampling patients
  answers the question a reader asks instead: "would other patients?"

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
    patients_for_standard_error,
    rank_parameters_by_influence,
    rescale_response_major,
)
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "sensitivity-patients-v1.1"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "sensitivity-precision-v1.0" / "metrics.json"

RUN_DATE = "2026-09-04"
SUBTYPE_COUNTS = (2, 3, 5)     # 8, 12 and 20 patients; 2 reproduces v1.0
N_SEEDS = 12                   # fixed at the v1.0 setting
SIMULATIONS = 1024             # fixed at the v1.0 setting
EPISODES_PER_POLICY = 40
EXPLORATION_WEIGHT = math.sqrt(2.0)

VALUE_JUDGEMENTS = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}

BOOTSTRAP_REPLICATES = 2_000
#: Standard error the cohort would need for the ranked effects (0.003-0.012) to
#: be separable; declared here so the sample-size answer is not chosen after
#: seeing the spread.
TARGET_STANDARD_ERROR = 0.005


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def set_response_major(data: dict, intensity: str, major: float) -> None:
    data["response_probabilities"][intensity] = rescale_response_major(
        data["response_probabilities"][intensity], major)


def build_variants(base: dict) -> list[dict]:
    """The v0.3 perturbation set, unchanged so the rankings line up."""
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
    """``(seed, patient)`` matrix of MCTS-minus-NCCN utility gaps.

    Keeping both axes rather than a single mean is what lets the same run answer
    "would other seeds rank these differently?" and "would other patients?".
    """
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
    """Largest |Δ vs baseline| per parameter, with the v1.0 seed-paired SE."""
    baseline = matrices["baseline"].mean(axis=1)
    records = []
    for variant in variants:
        if variant["name"] == "baseline":
            continue
        delta = matrices[variant["name"]].mean(axis=1) - baseline
        sd = float(delta.std(ddof=1))
        stderr = sd / math.sqrt(N_SEEDS)
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


def bootstrap_ranking(variants: list[dict], matrices: dict[str, np.ndarray],
                      replicates: int, seed: int) -> dict[str, object]:
    """Resample patients with replacement and re-rank the parameters each time.

    Every variant is evaluated on the same patients under the same seeds, so a
    replicate resamples the *patient index* and applies it to all of them. The
    result answers what the seed-paired standard error cannot: how much of the
    ranking is a property of these particular eight (or twenty) people.
    """
    rng = np.random.default_rng(seed)
    parameters = sorted({v["param"] for v in variants if v["param"] != "(none)"})
    by_parameter: dict[str, list[str]] = {p: [] for p in parameters}
    for parameter in parameters:
        by_parameter[parameter] = [
            v["name"] for v in variants if v["param"] == parameter]

    n_patients = matrices["baseline"].shape[1]
    rank_counts = {p: np.zeros(len(parameters), dtype=int) for p in parameters}
    top_two_hits = 0
    for _ in range(replicates):
        index = rng.integers(0, n_patients, size=n_patients)
        baseline = matrices["baseline"][:, index].mean()
        influence = {
            parameter: max(
                (matrices[name][:, index].mean() - baseline
                 for name in names), key=abs)
            for parameter, names in by_parameter.items()
        }
        order = rank_parameters_by_influence(influence)
        for position, parameter in enumerate(order):
            rank_counts[parameter][position] += 1
        if set(order[:2]) == VALUE_JUDGEMENTS:
            top_two_hits += 1
    return {
        "replicates": replicates,
        "top_two_are_value_judgements_share": top_two_hits / replicates,
        "rank_distribution": {
            parameter: (counts / replicates).round(4).tolist()
            for parameter, counts in rank_counts.items()
        },
        "share_ranked_first": {
            parameter: float(counts[0] / replicates)
            for parameter, counts in rank_counts.items()
        },
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    variants = build_variants(base)

    cohorts = {}
    for per_subtype in SUBTYPE_COUNTS:
        sample = balanced_subtype_sample(os_test, per_subtype)
        cohorts[len(sample)] = sample
    sizes = sorted(cohorts)
    largest = cohorts[sizes[-1]]
    nested = all(
        set(cohorts[n]["patient_id"]).issubset(set(largest["patient_id"]))
        for n in sizes)
    print(f"patient-count sensitivity: cohorts {sizes}, {N_SEEDS} seeds, "
          f"budget {SIMULATIONS}, {len(variants)} variants (nested={nested})",
          flush=True)

    detail_rows: list[dict] = []
    per_patient_rows: list[dict] = []
    influence_by_size: dict[int, pd.DataFrame] = {}
    matrices_by_size: dict[int, dict[str, np.ndarray]] = {}
    started = time.perf_counter()
    for n_patients in sizes:
        sample = cohorts[n_patients]
        matrices: dict[str, np.ndarray] = {}
        for variant in variants:
            gaps = evaluate(config_from_dict(variant["data"]), sample,
                            os_model, rfs_model)
            matrices[variant["name"]] = gaps
            per_seed = gaps.mean(axis=1)
            detail_rows.append({
                "n_patients": n_patients,
                "variant": variant["name"],
                "parameter": variant["param"],
                "value": variant["value"],
                "utility_gap": float(per_seed.mean()),
                "utility_gap_sd": float(per_seed.std(ddof=1)),
                "standard_error": float(per_seed.std(ddof=1) / math.sqrt(N_SEEDS)),
            })
            print(f"  [n={n_patients:2d}] {variant['name']:30s} "
                  f"gap={per_seed.mean():+.4f} "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)
        matrices_by_size[n_patients] = matrices
        influence = influence_table(variants, matrices)
        influence.insert(0, "n_patients", n_patients)
        influence_by_size[n_patients] = influence

        baseline_by_patient = matrices["baseline"].mean(axis=0)
        for patient_id, gap in zip(sample["patient_id"], baseline_by_patient):
            per_patient_rows.append({
                "n_patients": n_patients,
                "patient_id": str(patient_id),
                "baseline_utility_gap": float(gap),
            })

    detail = pd.DataFrame(detail_rows)
    influence = pd.concat(influence_by_size.values(), ignore_index=True)
    per_patient = pd.DataFrame(per_patient_rows)

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    prior_1024 = next(b for b in prior["by_budget"] if b["budget"] == "1024")
    prior_order = list(prior_1024["ranking"])

    by_size = []
    for n_patients in sizes:
        frame = influence_by_size[n_patients]
        order = list(frame["parameter"])
        block = detail[detail.n_patients == n_patients]
        baseline_gap = float(
            block[block.variant == "baseline"]["utility_gap"].iloc[0])
        by_size.append({
            "n_patients": n_patients,
            "baseline_utility_gap": baseline_gap,
            "baseline_seed_standard_error": float(
                block[block.variant == "baseline"]["standard_error"].iloc[0]),
            "top_two_are_value_judgements": bool(set(order[:2]) == VALUE_JUDGEMENTS),
            "ranking": order,
            "matches_v1_0_ranking": order == prior_order,
            "spearman_vs_v1_0": float(pd.Series(
                [prior_order.index(p) for p in order]).corr(
                pd.Series(range(len(order))), method="spearman")),
            "parameters_distinguishable_from_zero": int(
                (frame["delta_z"].abs() >= 2).sum()),
            "negative_variant_cells": int((block["utility_gap"] < 0).sum()),
            "variant_cells": int(len(block)),
        })

    largest_size = sizes[-1]
    largest_baseline = per_patient[per_patient.n_patients == largest_size]
    per_patient_sd = float(largest_baseline["baseline_utility_gap"].std(ddof=1))
    heterogeneity = {
        "n_patients": largest_size,
        "per_patient_gap_sd": per_patient_sd,
        "per_patient_gap_min": float(largest_baseline["baseline_utility_gap"].min()),
        "per_patient_gap_max": float(largest_baseline["baseline_utility_gap"].max()),
        "patients_with_negative_gap": int(
            (largest_baseline["baseline_utility_gap"] < 0).sum()),
        "between_patient_standard_error": {
            str(n): per_patient_sd / math.sqrt(n) for n in sizes},
        "target_standard_error": TARGET_STANDARD_ERROR,
        "patients_needed_for_target": patients_for_standard_error(
            per_patient_sd, TARGET_STANDARD_ERROR),
    }

    bootstrap = {
        str(n): bootstrap_ranking(variants, matrices_by_size[n],
                                  BOOTSTRAP_REPLICATES, BASE_SEED + n)
        for n in sizes
    }

    reproduction = {
        "v1_0_baseline_utility_gap": prior_1024["baseline_utility_gap"],
        "rerun_baseline_utility_gap": by_size[0]["baseline_utility_gap"],
        "absolute_difference": abs(
            prior_1024["baseline_utility_gap"] - by_size[0]["baseline_utility_gap"]),
        "reproduces_v1_0": bool(abs(
            prior_1024["baseline_utility_gap"]
            - by_size[0]["baseline_utility_gap"]) < 1e-9),
        "ranking_matches_v1_0": by_size[0]["ranking"] == prior_order,
    }

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "sensitivity-patients-v1.1",
        "question": (
            "v1.0 settled seeds and budget and named the eight patients as the "
            "next bottleneck. Does the sensitivity ranking hold when the cohort "
            "grows to twenty, and how many patients would it actually need?"
        ),
        "estimand": (
            "Largest absolute shift in the MCTS-minus-NCCN utility gap that each "
            "assumption produces when moved to either end of its v0.3 range, as "
            "the simulated cohort grows from 8 to 20 subtype-balanced patients."
        ),
        "scope_warning": (
            "Sensitivity of our own simulator to its declared assumptions. Not a "
            "clinical effect."
        ),
        "design": {
            "patient_counts": sizes,
            "cohorts_are_nested": bool(nested),
            "seeds": N_SEEDS,
            "simulations": SIMULATIONS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "variants": len(variants),
            "config": "configs/dynamic_v0_5.json",
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        },
        "v1_0_reference": {
            "patients": prior["design"]["patients"],
            "seeds": prior["design"]["seeds"],
            "simulations": 1024,
            "baseline_utility_gap": prior_1024["baseline_utility_gap"],
            "ranking": prior_order,
        },
        "reproduction_check": reproduction,
        "by_patient_count": by_size,
        "patient_heterogeneity": heterogeneity,
        "patient_bootstrap": bootstrap,
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
        "entry_point": "analysis/26_run_sensitivity_patients.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== reproduction of v1.0 (n=8) ===")
    print(json.dumps(reproduction, indent=2))
    print("\n=== influence ranking by cohort size ===")
    for n_patients in sizes:
        print(f"\n[n={n_patients}]")
        print(influence_by_size[n_patients][[
            "parameter", "abs_gap_delta", "delta_standard_error", "delta_z",
        ]].round(4).to_string(index=False))
    print("\n=== by cohort size ===")
    print(pd.DataFrame(by_size)[[
        "n_patients", "baseline_utility_gap", "top_two_are_value_judgements",
        "matches_v1_0_ranking", "parameters_distinguishable_from_zero",
        "negative_variant_cells",
    ]].to_string(index=False))
    print("\n=== patient heterogeneity ===")
    print(json.dumps(heterogeneity, indent=2))
    print("\n=== bootstrap: share where the top two are the value judgements ===")
    for n_patients in sizes:
        share = bootstrap[str(n_patients)]["top_two_are_value_judgements_share"]
        print(f"  n={n_patients:2d}  {share:.1%}")
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
