"""Two-dimensional sensitivity grids (v0.8): do the assumptions interact?

v0.3 moved one assumption at a time and ranked them by how far each shifted the
MCTS-vs-NCCN utility gap. One-at-a-time analysis answers "which assumption
matters most" but is blind to "does assumption A matter *differently* depending
on B" - and if the two value judgements that topped that ranking reinforce each
other, the honest uncertainty around our headline number is wider than the
one-at-a-time range suggested.

Three grids, each crossing the dominant value judgement with something else:

* **judgement x judgement** - recurrence-free reward against toxicity penalty.
  Both are ours to choose, so any interaction here is pure value-space.
* **judgement x clinical** - recurrence-free reward against intensified-chemo
  toxicity probability, a quantity K-CURE could actually estimate. Interaction
  here would mean the value of measuring that probability depends on how we
  weigh recurrence.
* **judgement x discounting** - recurrence-free reward against the annual
  discount rate, which v0.5 made declarable. Discounting was silently 0 for
  v0.2-v0.4, so this is the first look at what that choice was worth.

Runs on the corrected v0.5 environment at the v0.4 budget (1024). Still synthetic
assumptions: this measures the sensitivity of our own model, not clinical effect.
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
REPORT_DIR = ROOT / "reports" / "interaction-sensitivity-v0.8"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-28"
PATIENTS_PER_SUBTYPE = 2      # 8 patients; 27 grid cells have to stay affordable
N_SEEDS = 3
EPISODES_PER_POLICY = 30
SIMULATIONS = 1024            # v0.4: 256 cannot resolve an action ordering
EXPLORATION_WEIGHT = math.sqrt(2.0)


def set_recurrence_reward(data: dict, value: float) -> None:
    data["reward"]["recurrence_free_year"] = value


def set_toxicity_penalty(data: dict, value: float) -> None:
    data["reward"]["acute_toxicity_penalty"] = value


def set_intensified_toxicity(data: dict, value: float) -> None:
    data["acute_toxicity_probabilities"]["chemo"]["intensified"] = value


def set_discount_rate(data: dict, value: float) -> None:
    data["discount_rate_annual"] = value


#: (label, axis-x spec, axis-y spec). The x axis is the same value judgement in
#: every grid - the one v0.3 ranked first - so the three grids are comparable.
GRIDS = (
    {
        "name": "judgement_x_judgement",
        "title": "무재발 보상 × 독성 페널티",
        "x": ("reward.recurrence_free_year", (0.25, 0.50, 1.00), set_recurrence_reward),
        "y": ("reward.acute_toxicity_penalty", (0.05, 0.15, 0.30), set_toxicity_penalty),
        "kind": "가치판단 × 가치판단",
    },
    {
        "name": "judgement_x_clinical",
        "title": "무재발 보상 × 강화항암 독성확률",
        "x": ("reward.recurrence_free_year", (0.25, 0.50, 1.00), set_recurrence_reward),
        "y": ("toxicity.chemo.intensified", (0.15, 0.30, 0.45), set_intensified_toxicity),
        "kind": "가치판단 × 데이터로 추정 가능",
    },
    {
        "name": "judgement_x_discount",
        "title": "무재발 보상 × 연 할인율",
        "x": ("reward.recurrence_free_year", (0.25, 0.50, 1.00), set_recurrence_reward),
        "y": ("discount_rate_annual", (0.00, 0.03, 0.05), set_discount_rate),
        "kind": "가치판단 × 시간선호",
    },
)


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def evaluate(config: DynamicConfig, sample, os_model, rfs_model) -> dict:
    """Mean utility gap over seeds for one assumption setting."""
    gaps, neo, extended = [], [], []
    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        mcts_util, nccn_util, m_neo, m_ext = [], [], [], []
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
            m_neo.append((md["timing"] == "neoadjuvant").mean())
            m_ext.append((md["endocrine"] == "extended").mean())
        gaps.append(np.mean(mcts_util) - np.mean(nccn_util))
        neo.append(np.mean(m_neo) * 100)
        extended.append(np.mean(m_ext) * 100)
    return {
        "utility_gap": float(np.mean(gaps)),
        "utility_gap_sd": float(np.std(gaps, ddof=1)) if len(gaps) > 1 else 0.0,
        "neoadjuvant_rate_pct": float(np.mean(neo)),
        "extended_endocrine_rate_pct": float(np.mean(extended)),
    }


def interaction_summary(
    surface: np.ndarray,
    seed_sd: np.ndarray,
) -> dict[str, float]:
    """How far the grid departs from "each axis acts on its own".

    ``additive_residual`` is the largest gap between the observed surface and the
    best additive fit (row means plus column means minus the overall mean). A
    surface that is truly additive fits exactly, so any residual is interaction.
    ``corner_difference_in_differences`` reads the same thing off the four
    corners, which is the version that survives being put on a slide.

    Both are differences of noisy cells, so the corner statistic carries its own
    standard error, propagated from the per-cell seed spread. Without it a grid
    this small could show a large "interaction" that is only seed noise.
    """
    grand = surface.mean()
    additive = (surface.mean(axis=1, keepdims=True)
                + surface.mean(axis=0, keepdims=True) - grand)
    residual = surface - additive
    did = ((surface[-1, -1] - surface[-1, 0])
           - (surface[0, -1] - surface[0, 0]))
    corner_se = np.array([
        seed_sd[0, 0], seed_sd[0, -1], seed_sd[-1, 0], seed_sd[-1, -1]
    ]) / math.sqrt(N_SEEDS)
    did_se = float(np.sqrt(np.square(corner_se).sum()))
    return {
        "range": float(surface.max() - surface.min()),
        "max_additive_residual": float(np.abs(residual).max()),
        "corner_difference_in_differences": float(did),
        "corner_did_standard_error": did_se,
        "corner_did_z": float(did / did_se) if did_se > 0 else float("nan"),
        "mean_cell_standard_error": float(seed_sd.mean() / math.sqrt(N_SEEDS)),
        "interaction_share_of_range": float(
            np.abs(residual).max() / (surface.max() - surface.min())
            if surface.max() > surface.min() else float("nan")),
        "sign_flips": bool(surface.min() < 0 < surface.max()),
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    print(f"interaction sample = {len(sample)} patients, "
          f"{sum(len(g['x'][1]) * len(g['y'][1]) for g in GRIDS)} cells", flush=True)

    rows: list[dict] = []
    summaries: list[dict] = []
    started = time.perf_counter()
    for grid in GRIDS:
        x_name, x_values, x_setter = grid["x"]
        y_name, y_values, y_setter = grid["y"]
        surface = np.zeros((len(x_values), len(y_values)))
        seed_sd = np.zeros_like(surface)
        for i, x_value in enumerate(x_values):
            for j, y_value in enumerate(y_values):
                data = copy.deepcopy(base)
                x_setter(data, x_value)
                y_setter(data, y_value)
                result = evaluate(
                    config_from_dict(data), sample, os_model, rfs_model)
                surface[i, j] = result["utility_gap"]
                seed_sd[i, j] = result["utility_gap_sd"]
                rows.append({
                    "grid": grid["name"], "kind": grid["kind"],
                    "x_parameter": x_name, "x_value": x_value,
                    "y_parameter": y_name, "y_value": y_value,
                    **result,
                })
                print(f"  {grid['name']:22s} {x_value:<5} x {y_value:<5} "
                      f"gap={result['utility_gap']:+.4f} "
                      f"({time.perf_counter() - started:.0f}s)", flush=True)
        summaries.append({
            "grid": grid["name"], "title": grid["title"], "kind": grid["kind"],
            "x_parameter": x_name, "y_parameter": y_name,
            **interaction_summary(surface, seed_sd),
        })

    detail = pd.DataFrame(rows)
    summary = pd.DataFrame(summaries)
    print("\n=== interaction summary ===")
    print(summary[[
        "grid", "range", "max_additive_residual",
        "corner_difference_in_differences", "interaction_share_of_range",
        "sign_flips",
    ]].round(4).to_string(index=False))

    strongest = max(summaries, key=lambda row: abs(
        row["corner_difference_in_differences"]))
    one_at_a_time_range = float(
        detail[detail["grid"] == "judgement_x_judgement"]["utility_gap"].max()
        - detail[detail["grid"] == "judgement_x_judgement"]["utility_gap"].min())

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "interaction-sensitivity-v0.8",
        "question": (
            "Do the synthetic assumptions interact, or does one-at-a-time "
            "sensitivity already capture the uncertainty in the utility gap?"
        ),
        "estimand": (
            "Mean MCTS-minus-NCCN utility gap in the corrected v0.5 environment, "
            "as a function of two assumptions varied jointly."
        ),
        "scope_warning": (
            "Sensitivity of our own simulator to its declared assumptions. Not a "
            "clinical effect and not an estimate of any real quantity."
        ),
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "simulations_per_decision": SIMULATIONS,
            "config": "configs/dynamic_v0_5.json",
            "cells": int(len(rows)),
        },
        "grids": summaries,
        "strongest_interaction": strongest,
        "judgement_grid_range": one_at_a_time_range,
    }

    detail.to_csv(TABLE_DIR / "interaction_grid.csv", index=False)
    summary.to_csv(TABLE_DIR / "interaction_summary.csv", index=False)
    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "assumptions": CONFIG_PATH}),
        "entry_point": "analysis/19_run_interaction_sensitivity.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(f"\nstrongest interaction: {strongest['grid']} "
          f"(DiD {strongest['corner_difference_in_differences']:+.4f})")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
