"""Search-budget scaling study (v0.4): is the v0.3 instability noise or equipoise?

The v0.3 robustness study found that the MCTS first decision agreed across seeds
only 62.9% of the time. Two very different explanations fit that number:

1. **Search noise** - 256 simulations are too few, so the value estimates are
   noisy and the argmax flips at random. More budget would fix it.
2. **Genuine equipoise** - the competing actions really do have (almost) the
   same expected utility under our assumptions, so no budget can separate them.

They imply opposite next steps, so this script separates them. For every
decision node a patient actually reaches, it re-runs the search ``REPLICATES``
times at each budget and tracks two quantities as the budget grows:

* ``value_noise_sd`` - how much the *estimated* value of an action moves between
  independent replicates. Pure Monte-Carlo noise; it must shrink roughly as
  1/sqrt(budget) if the search is behaving.
* ``value_gap`` - the distance between the best and second-best action's mean
  estimated value. A property of the environment, not of the budget.

If the noise shrinks while agreement stays flat and the gap stays near zero, the
instability is equipoise, and reporting "MCTS recommends X" is unjustifiable no
matter how long we search. Part B then checks whether the headline utility gap
itself moves with budget.

Same synthetic environment as v0.2/v0.3: not causal, not clinical.
"""

from __future__ import annotations

from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import random
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
from analysis.dynamic.experiment_utils import (  # noqa: E402
    confidence_interval,
    summarize_decision_replicates,
)
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import DynamicState, patient_from_row  # noqa: E402
from analysis.dynamic.search import stochastic_mcts_search  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
# The environment audit (v0.5) neutralised the response channel and declared the
# discount rate, so this study now runs on the corrected assumptions. The v0.2
# config and the original report directory are left untouched: those results are
# the published record of what the biased environment produced, reproducible from
# their manifests' git_commit_before_run.
CONFIG_PATH = ROOT / "configs" / "dynamic_v0_5.json"
REPORT_DIR = ROOT / "reports" / "budget-scaling-v0.5env"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-28"
PATIENTS_PER_SUBTYPE = 3           # the same 12 patients as robustness-v0.3
EXPLORATION_WEIGHT = math.sqrt(2.0)

# Part A - decision-level convergence.
BUDGETS = (64, 128, 256, 512, 1024, 2048)
REPLICATES = 10                    # independent searches per (node, budget)
NODE_SAMPLE_SEED = BASE_SEED + 31  # seed of the walk that collects decision nodes
UNSTABLE_AGREEMENT_PCT = 75.0      # a node below this is called unstable

# Part B - policy-level outcome (a subset of budgets; each is a full re-run).
POLICY_BUDGETS = (64, 256, 1024)
REFERENCE_BUDGET = 256             # the v0.2/v0.3 setting every budget is paired against
N_SEEDS = 10
EPISODES_PER_POLICY = 25


def patient_offset(patient_id: str) -> int:
    """Stable per-patient seed offset (``hash`` is salted per process)."""
    return int(hashlib.sha256(patient_id.encode("utf-8")).hexdigest()[:8], 16)


def collect_decision_nodes(
    environment: DynamicBreastCancerEnvironment,
    seed: int,
) -> list[DynamicState]:
    """States with a real choice that this patient reaches under the guideline.

    Walking with the NCCN policy - rather than enumerating every reachable state
    - keeps the node set clinically meaningful and small: these are the forks a
    guideline-treated patient actually arrives at. Single-action states are
    dropped because there is nothing for the search to get wrong there.
    """
    policy = DynamicNccnPolicy(environment)
    rng = random.Random(seed)
    state = environment.initial_state()
    nodes: list[DynamicState] = []
    steps = 0
    while not environment.is_terminal(state):
        if steps >= 32:
            raise RuntimeError("node walk exceeded the environment horizon")
        if len(environment.legal_actions(state)) >= 2:
            nodes.append(state)
        state, _, _ = environment.step(state, policy(state), rng)
        steps += 1
    return nodes


def run_part_a(environments, patient_meta, node_map) -> pd.DataFrame:
    rows: list[dict] = []
    total = sum(len(nodes) for nodes in node_map.values()) * len(BUDGETS)
    done = 0
    started = time.perf_counter()
    for pid, nodes in node_map.items():
        environment = environments[pid]
        offset = patient_offset(pid) % 100_000
        for node_index, state in enumerate(nodes):
            for budget in BUDGETS:
                records = []
                for replicate in range(REPLICATES):
                    seed = (
                        BASE_SEED
                        + offset
                        + node_index * 1_000_003
                        + budget * 7_919
                        + replicate * 101
                    )
                    result = stochastic_mcts_search(
                        environment, state, simulations=budget,
                        exploration_weight=EXPLORATION_WEIGHT, seed=seed,
                    )
                    records.append((result.action, result.action_values))
                rows.append({
                    "patient_id": pid,
                    "subtype": patient_meta[pid]["subtype"],
                    "node_index": node_index,
                    "phase": state.phase,
                    "is_root": node_index == 0 and state.phase == "timing",
                    "budget": budget,
                    **summarize_decision_replicates(records),
                })
                done += 1
        print(f"  part A {done}/{total} node-budget cells "
              f"({time.perf_counter() - started:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def run_part_b(environments) -> pd.DataFrame:
    """Utility gap per budget. The NCCN arm is budget-free, so it runs once."""
    rows: list[dict] = []
    nccn_cache: dict[tuple[str, int], pd.DataFrame] = {}
    started = time.perf_counter()
    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        for budget in POLICY_BUDGETS:
            mcts_util, nccn_util = [], []
            mcts_surv, nccn_surv, mcts_tox = [], [], []
            for pid, environment in environments.items():
                patient_seed = seed + patient_offset(pid) % 100_000
                mcts_policy = CachedMCTSPolicy(
                    environment, simulations=budget,
                    exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed,
                )
                md = pd.DataFrame(run_policy_episodes(
                    environment, mcts_policy, EPISODES_PER_POLICY,
                    patient_seed + 10_000))
                key = (pid, seed_index)
                if key not in nccn_cache:
                    nccn_cache[key] = pd.DataFrame(run_policy_episodes(
                        environment, DynamicNccnPolicy(environment),
                        EPISODES_PER_POLICY, patient_seed + 10_000))
                nd = nccn_cache[key]
                mcts_util.append(md["utility"].mean())
                nccn_util.append(nd["utility"].mean())
                mcts_surv.append(md["survived_5y"].mean())
                nccn_surv.append(nd["survived_5y"].mean())
                mcts_tox.append(md["toxicity_count"].mean())
            rows.append({
                "seed_index": seed_index,
                "budget": budget,
                "mcts_utility": float(np.mean(mcts_util)),
                "nccn_utility": float(np.mean(nccn_util)),
                "utility_gap": float(np.mean(mcts_util) - np.mean(nccn_util)),
                "survival_gap_pp": float(
                    (np.mean(mcts_surv) - np.mean(nccn_surv)) * 100),
                "mcts_toxicity": float(np.mean(mcts_tox)),
            })
        print(f"  part B seed {seed_index + 1}/{N_SEEDS} "
              f"({time.perf_counter() - started:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    if REFERENCE_BUDGET not in POLICY_BUDGETS:
        raise ValueError("REFERENCE_BUDGET must be one of POLICY_BUDGETS")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_dynamic_config(CONFIG_PATH)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    print(f"budget-scaling sample = {len(sample)} patients "
          f"({sample['subtype'].value_counts().to_dict()})", flush=True)

    environments: dict[str, DynamicBreastCancerEnvironment] = {}
    patient_meta: dict[str, dict] = {}
    for _, row in sample.iterrows():
        patient = patient_from_row(row)
        environments[patient.patient_id] = DynamicBreastCancerEnvironment(
            patient, make_risk_table(row, os_model, rfs_model), config
        )
        patient_meta[patient.patient_id] = {
            "subtype": patient.subtype, "stage": patient.stage,
        }

    node_map = {
        pid: collect_decision_nodes(environment, NODE_SAMPLE_SEED)
        for pid, environment in environments.items()
    }
    n_nodes = sum(len(nodes) for nodes in node_map.values())
    print(f"decision nodes with >=2 legal actions: {n_nodes}", flush=True)

    part_a = run_part_a(environments, patient_meta, node_map)
    part_b = run_part_b(environments)

    # --- Part A aggregation -------------------------------------------------
    def mean_abs(series: pd.Series) -> float:
        return float(np.nanmean(np.abs(series.to_numpy(dtype=float))))

    def unstable_pct(series: pd.Series) -> float:
        return float((series < UNSTABLE_AGREEMENT_PCT).mean() * 100)

    by_budget = part_a.groupby("budget").agg(
        mean_agreement_pct=("agreement_pct", "mean"),
        min_agreement_pct=("agreement_pct", "min"),
        mean_value_noise_sd=("value_noise_sd", "mean"),
        mean_abs_value_gap=("value_gap", mean_abs),
        mean_separation=("separation",
                         lambda s: float(np.nanmean(s.to_numpy(dtype=float)))),
        unstable_node_pct=("agreement_pct", unstable_pct),
    ).reset_index()

    root = part_a[part_a["is_root"]]
    root_by_budget = root.groupby("budget").agg(
        mean_agreement_pct=("agreement_pct", "mean"),
        mean_value_noise_sd=("value_noise_sd", "mean"),
        mean_abs_value_gap=("value_gap", mean_abs),
    ).reset_index()

    # Does the noise fall like Monte-Carlo error? Slope of log(sd) on
    # log(budget); -0.5 is the textbook rate.
    noise = by_budget[by_budget["mean_value_noise_sd"] > 0]
    if len(noise) > 1:
        slope = float(np.polyfit(np.log(noise["budget"]),
                                 np.log(noise["mean_value_noise_sd"]), 1)[0])
    else:
        slope = float("nan")
    agreement_slope = float(np.polyfit(
        np.log(by_budget["budget"]), by_budget["mean_agreement_pct"], 1)[0])

    # --- Part B aggregation -------------------------------------------------
    gap_rows: list[dict] = []
    for budget, group in part_b.groupby("budget"):
        stats = confidence_interval(group["utility_gap"].to_numpy())
        surv = confidence_interval(group["survival_gap_pp"].to_numpy())
        gap_rows.append({
            "budget": int(budget),
            "utility_gap_mean": stats["mean"],
            "utility_gap_ci95_low": stats["ci95_low"],
            "utility_gap_ci95_high": stats["ci95_high"],
            "survival_gap_pp_mean": surv["mean"],
            "survival_gap_pp_ci95_low": surv["ci95_low"],
            "survival_gap_pp_ci95_high": surv["ci95_high"],
            "mean_mcts_toxicity": float(group["mcts_toxicity"].mean()),
        })
    gap_by_budget = pd.DataFrame(gap_rows)

    # Every budget runs on the same seeds, so pair within seed: that removes the
    # seed-to-seed variance which otherwise swamps a small budget effect.
    wide = part_b.pivot(index="seed_index", columns="budget", values="utility_gap")
    paired_rows: list[dict] = []
    for budget in POLICY_BUDGETS:
        if budget == REFERENCE_BUDGET:
            continue
        delta = (wide[budget] - wide[REFERENCE_BUDGET]).to_numpy()
        stats = confidence_interval(delta)
        paired_rows.append({
            "budget": int(budget),
            "reference_budget": REFERENCE_BUDGET,
            "paired_gap_delta_mean": stats["mean"],
            "paired_gap_delta_ci95_low": stats["ci95_low"],
            "paired_gap_delta_ci95_high": stats["ci95_high"],
            "crosses_zero": bool(stats["ci95_low"] <= 0 <= stats["ci95_high"]),
        })
    paired = pd.DataFrame(paired_rows)

    lowest, highest = by_budget.iloc[0], by_budget.iloc[-1]
    noise_fell = (
        highest["mean_value_noise_sd"] < lowest["mean_value_noise_sd"] * 0.75
    )
    agreement_flat = (
        highest["mean_agreement_pct"] < lowest["mean_agreement_pct"] + 10.0
    )
    verdict = "equipoise" if (noise_fell and agreement_flat) else "budget_limited"

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "budget-scaling-v0.5env",
        "question": (
            "Is the v0.3 seed instability caused by too small a search budget, "
            "or by genuine equipoise between the competing actions?"
        ),
        "interpretation": (
            "Search-behaviour diagnostic on the same synthetic environment; "
            "not causal, not clinical."
        ),
        "design": {
            "patients": int(len(sample)),
            "decision_nodes": int(n_nodes),
            "budgets": list(BUDGETS),
            "replicates_per_node_budget": REPLICATES,
            "policy_budgets": list(POLICY_BUDGETS),
            "policy_seeds": N_SEEDS,
            "episodes_per_policy": EPISODES_PER_POLICY,
        },
        "verdict": verdict,
        "value_noise_log_log_slope": slope,
        "agreement_pct_per_log_budget": agreement_slope,
        "agreement_at_min_budget_pct": float(lowest["mean_agreement_pct"]),
        "agreement_at_max_budget_pct": float(highest["mean_agreement_pct"]),
        "value_noise_sd_at_min_budget": float(lowest["mean_value_noise_sd"]),
        "value_noise_sd_at_max_budget": float(highest["mean_value_noise_sd"]),
        "mean_abs_value_gap_at_max_budget": float(highest["mean_abs_value_gap"]),
        "unstable_node_pct_at_min_budget": float(lowest["unstable_node_pct"]),
        "unstable_node_pct_at_max_budget": float(highest["unstable_node_pct"]),
        "utility_gap_by_budget": gap_rows,
        "paired_utility_gap_vs_reference": paired_rows,
    }

    part_a.to_csv(TABLE_DIR / "node_budget_detail.csv", index=False)
    by_budget.to_csv(TABLE_DIR / "convergence_by_budget.csv", index=False)
    root_by_budget.to_csv(TABLE_DIR / "root_convergence_by_budget.csv", index=False)
    part_b.to_csv(TABLE_DIR / "policy_per_seed.csv", index=False)
    gap_by_budget.to_csv(TABLE_DIR / "utility_gap_by_budget.csv", index=False)
    paired.to_csv(TABLE_DIR / "utility_gap_paired_vs_256.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({
            "data": INPUT_CSV, "assumptions": CONFIG_PATH,
        }),
        "entry_point": "analysis/12_run_budget_scaling.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== convergence by budget ===")
    print(by_budget.to_string(index=False))
    print("\n=== root-decision convergence by budget ===")
    print(root_by_budget.to_string(index=False))
    print("\n=== utility gap by budget ===")
    print(gap_by_budget.to_string(index=False))
    print(f"\n=== utility gap paired against budget {REFERENCE_BUDGET} ===")
    print(paired.to_string(index=False))
    print(f"\nnoise log-log slope = {slope:+.3f} (Monte-Carlo ideal -0.5)")
    print(f"agreement per e-fold budget = {agreement_slope:+.2f} pp")
    print(f"verdict = {verdict}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
