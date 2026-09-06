"""What is the remaining MCTS advantage actually made of? (v1.5)

v1.4 neutralised the reward model's treatment coefficients and left a gap of
+0.0199, saying only that whatever remains must come from the *declared* channels
in ``configs/dynamic_v0_5.json``. Nobody has looked at which.

Reading the config makes an uncomfortable pattern visible before anything is run:

=================  ==========================  =========================
treatment level    declared survival benefit   declared cost
=================  ==========================  =========================
chemo standard     **none** (HR 1.00 / 1.00)   toxicity 0.15, burden 0.12
endocrine standard **none** (HR 1.00 / 1.00)   toxicity 0.05, burden 0.03
radiation local    **none** (HR 1.00 / 1.00)   toxicity 0.05, burden 0.05
chemo intensified  HR 0.95 / 0.90              toxicity 0.30, burden 0.20
endocrine extended HR 0.98 / 0.92              toxicity 0.10, burden 0.06
radiation regional HR 1.00 / 0.95              toxicity 0.12, burden 0.09
=================  ==========================  =========================

**Standard-of-care treatment is declared to be pure cost.** Until v1.4 the Cox
reward model's coefficients were the only thing standing in for "does standard
treatment help?", and they said it *hurts*. Now they say nothing, so the config's
zeros are all that is left - and the NCCN policy is the one that prescribes
standard treatment by guideline while MCTS is free to decline it.

This run measures how much of the remaining gap that accounts for, as a 2x2 over
the two declared halves.

===========  =================================  ==========================
arm          declared benefit                   declared cost
===========  =================================  ==========================
A baseline   as configured                      as configured
B cost off   as configured                      toxicity penalty and burden 0
C benefit off treatment hazards all 1.0         as configured
D both off   treatment hazards all 1.0          toxicity penalty and burden 0
===========  =================================  ==========================

The response channel is left alone in every arm. It is mean-neutralised (v0.5) so
it adds no expected benefit, but it is the thing MCTS can *adapt* to - the
closed-loop advantage a sequential policy is supposed to have. Isolating that is
a separate question and stays out of this decomposition.

PRE-SPECIFIED PREDICTIONS, recorded before the run
--------------------------------------------------
1. **Null control** - arm D has |gap| < 0.005. With treatment neither helping nor
   costing anything, a searching policy has nothing to find. If D is far from
   zero, this design is measuring something other than treatment and the rest of
   the run cannot be interpreted.
2. **Primary** - arm B is much smaller than arm A. If the remaining advantage is
   cost-avoidance, removing the cost removes it.
3. Arm C is **larger** than arm A: stripping the small benefits that intensified
   and extended treatment carry leaves even less reason to treat, so MCTS
   declines more and NCCN keeps paying.

Prediction 3 is the one most likely to be wrong, and 1 is the one that has to
hold for anything else to mean anything.
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
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "channel-decomposition-v1.5"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "reward-confounding-v1.4" / "metrics.json"

RUN_DATE = "2026-09-07"
PER_SUBTYPE = 5                 # v1.2's cohorts A and B: 40 patients
N_SEEDS = 12
SIMULATIONS = 1024
EPISODES_PER_POLICY = 40
EXPLORATION_WEIGHT = math.sqrt(2.0)
NULL_CONTROL_TOLERANCE = 0.005

#: Treatment decisions whose hazard multipliers count as "declared benefit".
#: ``response`` is deliberately absent - see the module docstring.
BENEFIT_CHANNELS = ("chemo", "endocrine", "radiation", "timing")

ARMS = (
    ("baseline", False, False, "A. 선언된 채널 그대로"),
    ("cost_off", True, False, "B. 비용 제거 (독성 페널티·치료 부담 0)"),
    ("benefit_off", False, True, "C. 이득 제거 (치료 위험비 전부 1.0)"),
    ("both_off", True, True, "D. 둘 다 제거 (영대조)"),
)

PRESPECIFIED_PREDICTION = {
    "null_control": (
        "Arm D has |gap| < 0.005. With treatment neither helping nor costing, a "
        "searching policy has nothing to find."
    ),
    "primary": (
        "Arm B is much smaller than arm A - if the remaining advantage is "
        "cost-avoidance, removing the cost removes it."
    ),
    "secondary": (
        "Arm C is larger than arm A: stripping the small benefits that "
        "intensified and extended treatment carry leaves even less reason to "
        "treat, so MCTS declines more and NCCN keeps paying."
    ),
    "why": (
        "The config declares zero survival benefit for standard chemotherapy, "
        "endocrine therapy and local radiotherapy while charging toxicity and "
        "burden for all three. NCCN prescribes them; MCTS may decline."
    ),
}

ACTION_FIELDS = ("timing", "surgery", "chemo", "endocrine", "radiation")


def zero_leaves(node: dict) -> None:
    """Set every numeric leaf under ``node`` to zero, in place."""
    for key, value in node.items():
        if isinstance(value, dict):
            zero_leaves(value)
        else:
            node[key] = 0.0


def unit_hazards(node: dict) -> None:
    """Set every numeric leaf under ``node`` to one, in place."""
    for key, value in node.items():
        if isinstance(value, dict):
            unit_hazards(value)
        else:
            node[key] = 1.0


def build_variant(base: dict, cost_off: bool, benefit_off: bool) -> dict:
    data = copy.deepcopy(base)
    if cost_off:
        data["reward"]["acute_toxicity_penalty"] = 0.0
        zero_leaves(data["treatment_burden"])
    if benefit_off:
        for channel in BENEFIT_CHANNELS:
            unit_hazards(data["hazard_multipliers"][channel])
    return data


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def declared_asymmetry(base: dict) -> pd.DataFrame:
    """Declared survival benefit against declared cost, level by level.

    Free to compute and it is the reason this analysis exists, so it ships as a
    table rather than as a claim in prose.
    """
    rows = []
    for channel in ("chemo", "endocrine", "radiation"):
        for level, hazards in base["hazard_multipliers"][channel].items():
            rows.append({
                "channel": channel,
                "level": level,
                "death_hazard": hazards["death"],
                "recurrence_hazard": hazards["recurrence"],
                "declares_any_benefit": bool(
                    hazards["death"] < 1.0 or hazards["recurrence"] < 1.0),
                "acute_toxicity_probability":
                    base["acute_toxicity_probabilities"][channel][level],
                "treatment_burden": base["treatment_burden"][channel][level],
            })
    return pd.DataFrame(rows)


def evaluate(sample, os_model, rfs_model, config) -> dict:
    gaps, mcts_means, nccn_means = [], [], []
    actions = {"MCTS": [], "NCCN": []}
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
            actions["MCTS"].append(md)
            actions["NCCN"].append(nd)
        gaps.append(float(np.mean(mcts_util) - np.mean(nccn_util)))
        mcts_means.append(float(np.mean(mcts_util)))
        nccn_means.append(float(np.mean(nccn_util)))

    mix = {}
    for policy, frames in actions.items():
        joined = pd.concat(frames, ignore_index=True)
        mix[policy] = {
            field: joined[field].astype(str).value_counts(normalize=True)
            .mul(100).round(3).to_dict()
            for field in ACTION_FIELDS if field in joined.columns
        }
    return {
        "per_seed_gap": gaps,
        "utility_gap": float(np.mean(gaps)),
        "standard_error": float(np.std(gaps, ddof=1) / math.sqrt(N_SEEDS)),
        "mcts_utility": float(np.mean(mcts_means)),
        "nccn_utility": float(np.mean(nccn_means)),
        "action_mix": mix,
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    # v1.4: the reward model carries confounded treatment coefficients, so every
    # arm here runs on the neutralised model. Treatment effects come only from
    # the config, which is exactly what this run decomposes.
    neutral_model = os_model.neutralise_treatment_terms()

    first = balanced_subtype_sample(os_test, PER_SUBTYPE)
    second = balanced_subtype_sample(os_test, PER_SUBTYPE, offset=PER_SUBTYPE)
    if set(first["patient_id"]) & set(second["patient_id"]):
        raise SystemExit("cohorts A and B are not disjoint")
    sample = pd.concat([first, second], ignore_index=True)

    asymmetry = declared_asymmetry(base)
    asymmetry.to_csv(TABLE_DIR / "declared_asymmetry.csv", index=False)
    print(asymmetry.to_string(index=False), flush=True)
    pure_cost = asymmetry[
        ~asymmetry["declares_any_benefit"]
        & (asymmetry["acute_toxicity_probability"] > 0)]
    print(f"\n순수 비용인 치료 수준: {len(pure_cost)}개 "
          f"({', '.join(pure_cost['channel'] + '.' + pure_cost['level'])})",
          flush=True)

    results = {}
    started = time.perf_counter()
    for key, cost_off, benefit_off, label in ARMS:
        config = config_from_dict(build_variant(base, cost_off, benefit_off))
        print(f"\n[{key}] running...", flush=True)
        result = evaluate(sample, neutral_model, rfs_model, config)
        result.update({"label": label, "cost_off": cost_off,
                       "benefit_off": benefit_off})
        results[key] = result
        print(f"  gap {result['utility_gap']:+.4f} "
              f"(SE {result['standard_error']:.4f})  "
              f"MCTS {result['mcts_utility']:.4f} vs NCCN {result['nccn_utility']:.4f}  "
              f"({time.perf_counter() - started:.0f}s)", flush=True)

    def paired(a: str, b: str) -> dict:
        difference = (np.array(results[b]["per_seed_gap"])
                      - np.array(results[a]["per_seed_gap"]))
        stderr = float(difference.std(ddof=1) / math.sqrt(N_SEEDS))
        return {
            "difference": float(difference.mean()),
            "standard_error": stderr,
            "z": float(difference.mean() / stderr) if stderr > 0 else float("nan"),
        }

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    baseline_gap = results["baseline"]["utility_gap"]
    verdict = {
        "gap_baseline": baseline_gap,
        "gap_cost_off": results["cost_off"]["utility_gap"],
        "gap_benefit_off": results["benefit_off"]["utility_gap"],
        "gap_both_off": results["both_off"]["utility_gap"],
        "null_control_passed": bool(
            abs(results["both_off"]["utility_gap"]) < NULL_CONTROL_TOLERANCE),
        "null_control_tolerance": NULL_CONTROL_TOLERANCE,
        "cost_off_vs_baseline": paired("baseline", "cost_off"),
        "benefit_off_vs_baseline": paired("baseline", "benefit_off"),
        "share_of_gap_from_declared_cost": float(
            1.0 - results["cost_off"]["utility_gap"] / baseline_gap)
        if baseline_gap else float("nan"),
        "primary_prediction_met": bool(
            results["cost_off"]["utility_gap"] < baseline_gap),
        "secondary_prediction_met": bool(
            results["benefit_off"]["utility_gap"] > baseline_gap),
        "reproduces_v1_4_neutral": bool(abs(
            baseline_gap - prior["verdict"]["gap_treatment_neutral"]) < 1e-9),
        "v1_4_neutral_gap": prior["verdict"]["gap_treatment_neutral"],
        "pure_cost_levels": pure_cost.apply(
            lambda row: f"{row['channel']}.{row['level']}", axis=1).tolist(),
    }

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "channel-decomposition-v1.5",
        "question": (
            "After v1.4 removed the reward model's confounded treatment "
            "coefficients, what is the remaining MCTS-minus-NCCN gap made of?"
        ),
        "estimand": (
            "MCTS-minus-NCCN utility gap over the same 40 patients and 12 seeds "
            "under a 2x2 over the config's declared treatment benefit and "
            "declared treatment cost, with the reward model treatment-neutral."
        ),
        "scope_warning": (
            "Decomposition of our own simulator's declared assumptions. Not a "
            "clinical effect."
        ),
        "prespecified_prediction": PRESPECIFIED_PREDICTION,
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "simulations": SIMULATIONS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "config": "configs/dynamic_v0_5.json",
            "reward_model": "treatment-neutral (v1.4)",
            "benefit_channels_zeroed": list(BENEFIT_CHANNELS),
            "response_channel": "left untouched in every arm",
        },
        "verdict": verdict,
        "declared_asymmetry": asymmetry.to_dict(orient="records"),
        "arms": results,
    }

    pd.DataFrame([
        {"arm": key, "seed_index": index, "utility_gap": value}
        for key, result in results.items()
        for index, value in enumerate(result["per_seed_gap"])
    ]).to_csv(TABLE_DIR / "per_seed_gaps.csv", index=False)
    pd.DataFrame([
        {"arm": key, "policy": policy, "field": field, "action": action, "pct": pct}
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
        "entry_point": "analysis/35_run_channel_decomposition.py",
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
