"""How much did the response-channel bias inflate our results? (v0.5)

An environment audit found that the response hazard channel fired for one arm
only. A patient who chose neoadjuvant drew a response and picked up a response
hazard multiplier; a surgery-first patient stayed at ``not_applicable`` = 1.0
forever. Because those multipliers did not average to 1.0, choosing neoadjuvant
bought an expected recurrence-hazard discount of 6.5% (standard chemo) to 8.6%
(intensified) that no assumption ever declared.

This project exists to ask whether the NCCN guideline is optimal, so a channel
that hands one policy an undeclared advantage does not merely weaken the
answer - it voids it. ``configs/dynamic_v0_5.json`` plus the neutralisation in
``DynamicBreastCancerEnvironment`` fix it: the response channel is divided by
its own expectation, which keeps major-beats-none ordering while making the
channel mean-neutral, and any real timing effect now lives in an explicit
``hazard_multipliers.timing`` block that sensitivity analysis can see.

This script measures what the bias was worth. It runs the same patients, seeds
and search budget under the old and the fixed environment, paired within seed,
and reports how much of the MCTS advantage - and of the neoadjuvant preference
in particular - came from the bias rather than from planning.

Still the same synthetic assumptions otherwise: not causal, not clinical.
"""

from __future__ import annotations

import dataclasses
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
from analysis.dynamic.config import load_dynamic_config  # noqa: E402
from analysis.dynamic.environment import DynamicBreastCancerEnvironment  # noqa: E402
from analysis.dynamic.evaluation import run_policy_episodes  # noqa: E402
from analysis.dynamic.experiment_utils import confidence_interval  # noqa: E402
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
FIXED_CONFIG = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "environment-fix-v0.5"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-27"
PATIENTS_PER_SUBTYPE = 3          # same 12 patients as v0.3 / v0.4
N_SEEDS = 10
EPISODES_PER_POLICY = 25
SIMULATIONS = 1024                # the v0.4 finding: 256 cannot resolve actions
EXPLORATION_WEIGHT = math.sqrt(2.0)

#: The biased twin differs from the fixed config in exactly one flag, derived in
#: code rather than kept as a second file, so nothing else can drift between the
#: two arms and the v0.2 config keeps the SHA its published manifests recorded.
ARM_NAMES = ("biased", "fixed")


def patient_offset(patient_id: str) -> int:
    return int(hashlib.sha256(patient_id.encode("utf-8")).hexdigest()[:8], 16)


def build_environments(sample, os_model, rfs_model, config):
    environments = {}
    for _, row in sample.iterrows():
        patient = patient_from_row(row)
        environments[patient.patient_id] = DynamicBreastCancerEnvironment(
            patient, make_risk_table(row, os_model, rfs_model), config
        )
    return environments


def evaluate(environments, seed_index: int) -> dict[str, float]:
    """One seed: mean outcomes over patients for both policies."""
    seed = BASE_SEED + seed_index * 1_000
    mcts, nccn = [], []
    for pid, environment in environments.items():
        patient_seed = seed + patient_offset(pid) % 100_000
        policy = CachedMCTSPolicy(
            environment, simulations=SIMULATIONS,
            exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed,
        )
        mcts.append(pd.DataFrame(run_policy_episodes(
            environment, policy, EPISODES_PER_POLICY, patient_seed + 10_000)))
        nccn.append(pd.DataFrame(run_policy_episodes(
            environment, DynamicNccnPolicy(environment),
            EPISODES_PER_POLICY, patient_seed + 10_000)))
    md = pd.concat(mcts, ignore_index=True)
    nd = pd.concat(nccn, ignore_index=True)
    return {
        "mcts_utility": float(md["utility"].mean()),
        "nccn_utility": float(nd["utility"].mean()),
        "utility_gap": float(md["utility"].mean() - nd["utility"].mean()),
        "mcts_survival_pct": float(md["survived_5y"].mean() * 100),
        "survival_gap_pp": float(
            (md["survived_5y"].mean() - nd["survived_5y"].mean()) * 100),
        "mcts_recurrence_pct": float(md["recurred_by_5y"].mean() * 100),
        "mcts_neoadjuvant_pct": float((md["timing"] == "neoadjuvant").mean() * 100),
        "nccn_neoadjuvant_pct": float((nd["timing"] == "neoadjuvant").mean() * 100),
        "mcts_intensified_pct": float((md["chemo"] == "intensified").mean() * 100),
        "mcts_toxicity": float(md["toxicity_count"].mean()),
    }


def analytic_effect_estimate(config, neoadjuvant_share_pct: float) -> dict:
    """Order-of-magnitude size of the bias, so a null result can be read honestly.

    A paired experiment that finds nothing has only shown that the effect is
    smaller than what it could resolve. This estimates how large the effect
    could have been, so the report can say which of the two happened.

    For a representative patient the bias lowers the annual recurrence
    probability of the neoadjuvant arm from ``p_biased`` to ``p_fixed``. Over the
    horizon that buys extra recurrence-free years, each worth
    ``recurrence_free_year / max_followup_reward`` of normalised utility, and it
    only applies to the share of episodes that actually chose neoadjuvant.
    """
    horizon = int(config.horizon_years)
    weight = float(config.reward["recurrence_free_year"]) / config.max_followup_reward

    def recurrence_free_years(annual_probability: float) -> float:
        return sum((1.0 - annual_probability) ** year
                   for year in range(1, horizon + 1))

    # Representative patient, standard chemo (see reports/environment-fix-v0.5).
    p_neutral = 0.01769
    p_biased = p_neutral * config.response_multiplier_mean("standard", "recurrence")
    share = neoadjuvant_share_pct / 100.0
    delta_years = recurrence_free_years(p_biased) - recurrence_free_years(p_neutral)
    return {
        "representative_annual_recurrence_neutral": p_neutral,
        "representative_annual_recurrence_biased": p_biased,
        "extra_recurrence_free_years_from_bias": delta_years,
        "neoadjuvant_share_pct": neoadjuvant_share_pct,
        "expected_utility_gap_inflation": share * delta_years * weight,
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    print(f"impact sample = {len(sample)} patients", flush=True)

    fixed_config = load_dynamic_config(FIXED_CONFIG)
    biased_config = dataclasses.replace(
        fixed_config, response_channel_neutralised=False)
    arms = {
        "biased": build_environments(sample, os_model, rfs_model, biased_config),
        "fixed": build_environments(sample, os_model, rfs_model, fixed_config),
    }

    rows = []
    started = time.perf_counter()
    for seed_index in range(N_SEEDS):
        for name in ARM_NAMES:
            rows.append({
                "seed_index": seed_index, "arm": name,
                **evaluate(arms[name], seed_index),
            })
        print(f"  seed {seed_index + 1}/{N_SEEDS} "
              f"({time.perf_counter() - started:.0f}s)", flush=True)
    per_seed = pd.DataFrame(rows)

    metrics_of_interest = [
        "utility_gap", "survival_gap_pp", "mcts_utility", "nccn_utility",
        "mcts_survival_pct", "mcts_recurrence_pct", "mcts_neoadjuvant_pct",
        "nccn_neoadjuvant_pct", "mcts_intensified_pct", "mcts_toxicity",
    ]

    summary_rows = []
    for metric in metrics_of_interest:
        wide = per_seed.pivot(index="seed_index", columns="arm", values=metric)
        delta = (wide["fixed"] - wide["biased"]).to_numpy()
        stats = confidence_interval(delta)
        summary_rows.append({
            "metric": metric,
            "biased_mean": float(wide["biased"].mean()),
            "fixed_mean": float(wide["fixed"].mean()),
            "paired_delta_mean": stats["mean"],
            "paired_delta_ci95_low": stats["ci95_low"],
            "paired_delta_ci95_high": stats["ci95_high"],
            "crosses_zero": bool(stats["ci95_low"] <= 0 <= stats["ci95_high"]),
        })
    summary = pd.DataFrame(summary_rows)

    def row_for(metric: str) -> dict:
        return summary[summary["metric"].eq(metric)].iloc[0].to_dict()

    gap = row_for("utility_gap")
    neo = row_for("mcts_neoadjuvant_pct")
    share = (
        abs(gap["paired_delta_mean"]) / abs(gap["biased_mean"]) * 100
        if gap["biased_mean"] else float("nan")
    )

    analytic = analytic_effect_estimate(fixed_config, neo["biased_mean"])
    resolvable = (gap["paired_delta_ci95_high"] - gap["paired_delta_ci95_low"]) / 2
    analytic["smallest_resolvable_shift_ci95_half_width"] = resolvable
    analytic["design_is_powered_for_it"] = bool(
        abs(analytic["expected_utility_gap_inflation"]) > resolvable)

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "environment-fix-v0.5",
        "question": (
            "How much of the measured MCTS advantage came from the asymmetric "
            "response hazard channel rather than from planning?"
        ),
        "interpretation": (
            "Paired within-seed comparison of the same experiment under the "
            "biased and the neutralised environment; not causal, not clinical."
        ),
        "defect": {
            "channel": "hazard_multipliers.response",
            "why": (
                "Only a neoadjuvant patient draws a response, so only that arm "
                "picks up a response hazard multiplier; the multipliers did not "
                "average to 1.0."
            ),
            "expected_recurrence_multiplier_standard": 0.935,
            "expected_recurrence_multiplier_intensified": 0.914,
            "surgery_first_multiplier": 1.0,
        },
        "design": {
            "patients": int(len(sample)),
            "seeds": N_SEEDS,
            "episodes_per_policy": EPISODES_PER_POLICY,
            "simulations_per_decision": SIMULATIONS,
        },
        "utility_gap_biased": gap["biased_mean"],
        "utility_gap_fixed": gap["fixed_mean"],
        "utility_gap_shift": gap["paired_delta_mean"],
        "utility_gap_shift_ci95": [
            gap["paired_delta_ci95_low"], gap["paired_delta_ci95_high"]],
        "utility_gap_share_attributable_to_bias_pct": share,
        "mcts_neoadjuvant_pct_biased": neo["biased_mean"],
        "mcts_neoadjuvant_pct_fixed": neo["fixed_mean"],
        "mcts_neoadjuvant_shift_pp": neo["paired_delta_mean"],
        "mcts_neoadjuvant_shift_ci95": [
            neo["paired_delta_ci95_low"], neo["paired_delta_ci95_high"]],
        "metrics_moved_beyond_noise": [
            row["metric"] for row in summary_rows if not row["crosses_zero"]
        ],
        "analytic_effect_estimate": analytic,
    }

    per_seed.to_csv(TABLE_DIR / "per_seed_by_arm.csv", index=False)
    summary.to_csv(TABLE_DIR / "paired_effect_of_fix.csv", index=False)
    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({
            "data": INPUT_CSV, "assumptions": FIXED_CONFIG,
        }),
        "biased_arm": "same config with response_channel_neutralised=False",
        "entry_point": "analysis/13_run_environment_fix_impact.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== effect of neutralising the response channel ===")
    print(summary.to_string(index=False))
    print(f"\nutility gap {gap['biased_mean']:+.4f} -> {gap['fixed_mean']:+.4f}")
    print(f"neoadjuvant {neo['biased_mean']:.1f}% -> {neo['fixed_mean']:.1f}%")
    print("")
    print("expected inflation from the bias  = "
          f"{analytic['expected_utility_gap_inflation']:.6f}")
    print(f"smallest shift this design resolves = {resolvable:.6f}")
    print(f"powered to detect it = {analytic['design_is_powered_for_it']}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
