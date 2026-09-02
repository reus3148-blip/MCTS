"""Confirmatory run for the utility interaction (v0.9).

v0.8 found that the two value judgements reverse each other's direction - raising
the toxicity penalty helps MCTS when the recurrence reward is low and hurts it
when the reward is high - at a corner difference-in-differences of -0.0245 with a
standard error of 0.0115. That is z = -2.13 on three seeds: a signal to follow,
not a result to report.

This run tests the same grid with four times the seeds. Everything else is held
fixed - same patients, same episodes, same 1024 budget, same three-by-three
values - so the only change is precision. Quadrupling seeds should halve the
standard error, which turns a z of about 2 into about 4 if the effect is real and
leaves the estimate near zero if it was not.

Pre-specified before running: the interaction is confirmed if |z| >= 3 with the
sign of v0.8, called noise if |z| < 2, and reported as unresolved in between.
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
REPORT_DIR = ROOT / "reports" / "utility-interaction-v0.9"
TABLE_DIR = REPORT_DIR / "tables"
PRIOR_METRICS = ROOT / "reports" / "interaction-sensitivity-v0.8" / "metrics.json"

RUN_DATE = "2026-08-28"
PATIENTS_PER_SUBTYPE = 2
N_SEEDS = 12                  # v0.8 used 3; four times the seeds halves the SE
EPISODES_PER_POLICY = 30
SIMULATIONS = 1024
EXPLORATION_WEIGHT = math.sqrt(2.0)

RECURRENCE_REWARDS = (0.25, 0.50, 1.00)
TOXICITY_PENALTIES = (0.05, 0.15, 0.30)

#: Decision rule fixed before the run, so the verdict is not chosen after seeing z.
CONFIRM_Z = 3.0
NOISE_Z = 2.0


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def evaluate(config: DynamicConfig, sample, os_model, rfs_model) -> dict:
    """Per-seed utility gaps for one cell, kept individually for the paired test."""
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
                environment, simulations=SIMULATIONS,
                exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed)
            md = pd.DataFrame(run_policy_episodes(
                environment, policy, EPISODES_PER_POLICY, patient_seed + 10_000))
            nd = pd.DataFrame(run_policy_episodes(
                environment, DynamicNccnPolicy(environment),
                EPISODES_PER_POLICY, patient_seed + 10_000))
            mcts_util.append(md["utility"].mean())
            nccn_util.append(nd["utility"].mean())
        gaps.append(float(np.mean(mcts_util) - np.mean(nccn_util)))
    return {"per_seed_gaps": gaps, "utility_gap": float(np.mean(gaps)),
            "utility_gap_sd": float(np.std(gaps, ddof=1))}


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    cells = len(RECURRENCE_REWARDS) * len(TOXICITY_PENALTIES)
    print(f"confirmation grid: {len(sample)} patients, {N_SEEDS} seeds, "
          f"{cells} cells", flush=True)

    rows: list[dict] = []
    per_seed: dict[tuple[float, float], list[float]] = {}
    started = time.perf_counter()
    for reward in RECURRENCE_REWARDS:
        for penalty in TOXICITY_PENALTIES:
            data = copy.deepcopy(base)
            data["reward"]["recurrence_free_year"] = reward
            data["reward"]["acute_toxicity_penalty"] = penalty
            result = evaluate(config_from_dict(data), sample, os_model, rfs_model)
            per_seed[(reward, penalty)] = result["per_seed_gaps"]
            rows.append({
                "recurrence_free_year": reward,
                "acute_toxicity_penalty": penalty,
                "utility_gap": result["utility_gap"],
                "utility_gap_sd": result["utility_gap_sd"],
                "standard_error": result["utility_gap_sd"] / math.sqrt(N_SEEDS),
            })
            print(f"  보상={reward:<5} 페널티={penalty:<5} "
                  f"gap={result['utility_gap']:+.4f} "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)

    detail = pd.DataFrame(rows)

    # The four corner cells share their seeds, so the difference-in-differences
    # can be formed seed by seed and tested as a paired quantity. That removes
    # the seed-level variance the unpaired version of v0.8 had to carry.
    lo_r, hi_r = RECURRENCE_REWARDS[0], RECURRENCE_REWARDS[-1]
    lo_p, hi_p = TOXICITY_PENALTIES[0], TOXICITY_PENALTIES[-1]
    paired = np.array([
        (per_seed[(hi_r, hi_p)][i] - per_seed[(hi_r, lo_p)][i])
        - (per_seed[(lo_r, hi_p)][i] - per_seed[(lo_r, lo_p)][i])
        for i in range(N_SEEDS)
    ])
    did = float(paired.mean())
    did_se = float(paired.std(ddof=1) / math.sqrt(N_SEEDS))
    did_z = did / did_se if did_se > 0 else float("nan")

    # Also the unpaired form, so the number is comparable with v0.8's.
    corner_se = np.array([
        detail[(detail.recurrence_free_year == r)
               & (detail.acute_toxicity_penalty == p)]["standard_error"].iloc[0]
        for r, p in ((lo_r, lo_p), (lo_r, hi_p), (hi_r, lo_p), (hi_r, hi_p))
    ])
    unpaired_se = float(np.sqrt(np.square(corner_se).sum()))
    unpaired_z = did / unpaired_se if unpaired_se > 0 else float("nan")

    # The two simple effects the difference-in-differences is built from: what
    # raising the toxicity penalty does at each end of the reward axis.
    simple_low = (float(np.mean(per_seed[(lo_r, hi_p)]))
                  - float(np.mean(per_seed[(lo_r, lo_p)])))
    simple_high = (float(np.mean(per_seed[(hi_r, hi_p)]))
                   - float(np.mean(per_seed[(hi_r, lo_p)])))

    prior = json.loads(PRIOR_METRICS.read_text(encoding="utf-8"))
    prior_grid = next(g for g in prior["grids"] if g["grid"] == "judgement_x_judgement")

    if abs(did_z) >= CONFIRM_Z and np.sign(did) == np.sign(
            prior_grid["corner_difference_in_differences"]):
        verdict = "confirmed"
    elif abs(did_z) < NOISE_Z:
        verdict = "noise"
    else:
        verdict = "unresolved"

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "utility-interaction-v0.9",
        "question": (
            "Is the v0.8 interaction between the recurrence-free reward and the "
            "toxicity penalty real, or was it seed noise at three seeds?"
        ),
        "estimand": (
            "Corner difference-in-differences of the MCTS-minus-NCCN utility gap "
            "over the two reward weights, in the corrected v0.5 environment."
        ),
        "prespecified_rule": {
            "confirmed_if_abs_z_at_least": CONFIRM_Z,
            "noise_if_abs_z_below": NOISE_Z,
            "sign_must_match_v0_8": True,
        },
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "simulations_per_decision": SIMULATIONS,
            "cells": int(len(rows)),
            "config": "configs/dynamic_v0_5.json",
        },
        "v0_8_reference": {
            "seeds": 3,
            "corner_difference_in_differences":
                prior_grid["corner_difference_in_differences"],
            "standard_error": prior_grid["corner_did_standard_error"],
            "z": prior_grid["corner_did_z"],
        },
        "corner_difference_in_differences": did,
        "paired_standard_error": did_se,
        "paired_z": did_z,
        "unpaired_standard_error": unpaired_se,
        "unpaired_z": unpaired_z,
        "simple_effect_low_reward": simple_low,
        "simple_effect_high_reward": simple_high,
        "direction_reverses": bool(simple_low * simple_high < 0),
        "verdict": verdict,
    }

    detail.to_csv(TABLE_DIR / "confirmation_grid.csv", index=False)
    pd.DataFrame({
        "seed_index": range(N_SEEDS),
        "did_per_seed": paired,
    }).to_csv(TABLE_DIR / "did_per_seed.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "assumptions": CONFIG_PATH}),
        "entry_point": "analysis/22_run_utility_interaction_confirm.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== confirmation ===")
    print(f"v0.8 (3 seeds):  DiD {prior_grid['corner_difference_in_differences']:+.4f} "
          f"± {prior_grid['corner_did_standard_error']:.4f}  "
          f"z={prior_grid['corner_did_z']:+.2f}")
    print(f"v0.9 ({N_SEEDS} seeds): DiD {did:+.4f} ± {did_se:.4f} (paired)  "
          f"z={did_z:+.2f}")
    print(f"                 unpaired SE {unpaired_se:.4f}  z={unpaired_z:+.2f}")
    print(f"simple effects: 낮은 보상 {simple_low:+.4f}, 높은 보상 {simple_high:+.4f} "
          f"(부호 반전 {metrics['direction_reverses']})")
    print(f"verdict = {verdict}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
