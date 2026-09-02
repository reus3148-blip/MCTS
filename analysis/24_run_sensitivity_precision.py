"""One-at-a-time sensitivity at higher precision (v1.0).

v0.9 showed that three seeds shift a cell mean by 0.010-0.015 and can manufacture
an interaction that four times the seeds erases. The v0.3 one-at-a-time
sensitivity ran on the same eight patients and the same three seeds, so its
ranking inherits the same problem. Its top two - both value judgements - survived
the v0.5 environment re-run, but nothing has tested whether the *lower* half of
that ranking is real.

Two things differ from v0.3, and they are separated rather than changed together:

* **seeds 3 -> 12**, at the original 256 budget. Compared against
  ``sensitivity-v0.5env`` this isolates what precision alone buys.
* **budget 256 -> 1024** at 12 seeds. v0.4 established that 256 cannot resolve an
  action ordering, so this is the ranking under the current standard.

Reporting both means the reader can see which part of any change came from
precision and which from search resolution. Still synthetic assumptions.
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
from analysis.dynamic.experiment_utils import rescale_response_major  # noqa: E402
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "sensitivity-precision-v1.0"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "sensitivity-v0.5env" / "metrics.json"
PRIOR_INFLUENCE = ROOT / "reports" / "sensitivity-v0.5env" / "tables" / "parameter_influence.csv"

RUN_DATE = "2026-08-28"
PATIENTS_PER_SUBTYPE = 2       # same 8 patients as v0.3, so the ranking is comparable
N_SEEDS = 12                   # v0.3 used 3
EPISODES_PER_POLICY = 40
EXPLORATION_WEIGHT = math.sqrt(2.0)

#: (label, simulations). 256 isolates the seed effect against v0.5env; 1024 is
#: the current standard after v0.4.
BUDGETS = (("256", 256), ("1024", 1024))


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


def evaluate(config: DynamicConfig, sample, os_model, rfs_model,
             simulations: int) -> dict:
    """Per-seed utility gaps for one variant, kept for the paired comparison."""
    gaps = []
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
                environment, simulations=simulations,
                exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed)
            md = pd.DataFrame(run_policy_episodes(
                environment, policy, EPISODES_PER_POLICY, patient_seed + 10_000))
            nd = pd.DataFrame(run_policy_episodes(
                environment, DynamicNccnPolicy(environment),
                EPISODES_PER_POLICY, patient_seed + 10_000))
            mcts_util.append(md["utility"].mean())
            nccn_util.append(nd["utility"].mean())
        gaps.append(float(np.mean(mcts_util) - np.mean(nccn_util)))
    return {"per_seed": gaps, "utility_gap": float(np.mean(gaps)),
            "utility_gap_sd": float(np.std(gaps, ddof=1)),
            "standard_error": float(np.std(gaps, ddof=1) / math.sqrt(N_SEEDS))}


def influence_table(rows: list[dict], per_seed: dict[str, list[float]],
                    baseline_name: str) -> pd.DataFrame:
    """Largest |Δ vs baseline| per parameter, with a paired standard error.

    The variants share their seeds with the baseline, so each Δ can be formed
    seed by seed. That paired standard error is what says whether a rung of the
    ranking is distinguishable from the one below it.
    """
    baseline = np.array(per_seed[baseline_name])
    records = []
    for row in rows:
        if row["variant"] == baseline_name:
            continue
        delta = np.array(per_seed[row["variant"]]) - baseline
        records.append({
            "parameter": row["parameter"],
            "variant": row["variant"],
            "gap_delta": float(delta.mean()),
            "abs_gap_delta": float(abs(delta.mean())),
            "delta_standard_error": float(delta.std(ddof=1) / math.sqrt(N_SEEDS)),
            "delta_z": float(delta.mean() / (delta.std(ddof=1) / math.sqrt(N_SEEDS)))
            if delta.std(ddof=1) > 0 else float("nan"),
        })
    frame = pd.DataFrame(records)
    best = frame.loc[frame.groupby("parameter")["abs_gap_delta"].idxmax()]
    return best.sort_values("abs_gap_delta", ascending=False).reset_index(drop=True)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    variants = build_variants(base)
    print(f"precision sensitivity: {len(sample)} patients, {N_SEEDS} seeds, "
          f"{len(variants)} variants x {len(BUDGETS)} budgets", flush=True)

    detail_rows: list[dict] = []
    influence_by_budget: dict[str, pd.DataFrame] = {}
    started = time.perf_counter()
    for label, simulations in BUDGETS:
        rows: list[dict] = []
        per_seed: dict[str, list[float]] = {}
        for variant in variants:
            result = evaluate(config_from_dict(variant["data"]), sample,
                              os_model, rfs_model, simulations)
            per_seed[variant["name"]] = result["per_seed"]
            row = {
                "budget": label, "variant": variant["name"],
                "parameter": variant["param"], "value": variant["value"],
                "utility_gap": result["utility_gap"],
                "utility_gap_sd": result["utility_gap_sd"],
                "standard_error": result["standard_error"],
            }
            rows.append(row)
            detail_rows.append(row)
            print(f"  [{label}] {variant['name']:30s} "
                  f"gap={result['utility_gap']:+.4f} "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)
        influence = influence_table(rows, per_seed, "baseline")
        influence.insert(0, "budget", label)
        influence_by_budget[label] = influence

    detail = pd.DataFrame(detail_rows)
    influence = pd.concat(influence_by_budget.values(), ignore_index=True)

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    prior_rank = pd.read_csv(PRIOR_INFLUENCE)
    prior_order = list(prior_rank["parameter"])

    comparison = []
    for label, frame in influence_by_budget.items():
        order = list(frame["parameter"])
        comparison.append({
            "budget": label,
            "top_two_are_value_judgements": bool(set(order[:2]) == {
                "reward.recurrence_free_year", "reward.acute_toxicity_penalty"}),
            "ranking": order,
            "matches_v0_5env_ranking": order == prior_order,
            "spearman_vs_v0_5env": float(pd.Series(
                [prior_order.index(p) for p in order]).corr(
                pd.Series(range(len(order))), method="spearman")),
            "parameters_distinguishable_from_zero": int(
                (frame["delta_z"].abs() >= 2).sum()),
            "baseline_utility_gap": float(
                detail[(detail.budget == label)
                       & (detail.variant == "baseline")]["utility_gap"].iloc[0]),
        })

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "sensitivity-precision-v1.0",
        "question": (
            "Does the v0.3 sensitivity ranking survive four times the seeds, and "
            "does it change at the search budget v0.4 established as adequate?"
        ),
        "estimand": (
            "Largest absolute shift in the MCTS-minus-NCCN utility gap that each "
            "assumption produces when moved to either end of its v0.3 range."
        ),
        "scope_warning": (
            "Sensitivity of our own simulator to its declared assumptions. Not a "
            "clinical effect."
        ),
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "variants": len(variants),
            "budgets": [label for label, _ in BUDGETS],
            "config": "configs/dynamic_v0_5.json",
        },
        "v0_5env_reference": {
            "seeds": prior["design"]["n_seeds"],
            "simulations": 256,
            "baseline_utility_gap": prior["baseline_utility_gap"],
            "ranking": prior_order,
        },
        "by_budget": comparison,
    }

    detail.to_csv(TABLE_DIR / "variant_results.csv", index=False)
    influence.to_csv(TABLE_DIR / "parameter_influence.csv", index=False)
    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "assumptions": CONFIG_PATH}),
        "entry_point": "analysis/24_run_sensitivity_precision.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== influence ranking ===")
    for label in influence_by_budget:
        print(f"\n[budget {label}]")
        print(influence_by_budget[label][[
            "parameter", "abs_gap_delta", "delta_standard_error", "delta_z",
        ]].round(4).to_string(index=False))
    print("\n=== vs v0.5env (3 seeds, 256) ===")
    print(pd.DataFrame(comparison)[[
        "budget", "baseline_utility_gap", "top_two_are_value_judgements",
        "matches_v0_5env_ranking", "parameters_distinguishable_from_zero",
    ]].to_string(index=False))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
