"""Multi-seed robustness study for the dynamic MCTS policy (v0.3 analysis).

The v0.2 report ran a single fixed seed. A single seed cannot tell us whether
the MCTS-vs-NCCN comparison is stable or an artifact of one lucky random draw.
This script re-runs the dynamic evaluation across many independent seeds and
reports the mean and a 95% confidence interval for each headline metric, plus a
per-patient measure of how often the first MCTS decision agrees across seeds.

Nothing here is a clinical recommendation. The environment still uses the same
synthetic transition assumptions as v0.2; this study only measures *stochastic
stability*, not clinical validity or causal effect.
"""

from __future__ import annotations

from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.dynamic.config import load_dynamic_config  # noqa: E402
from analysis.dynamic.environment import DynamicBreastCancerEnvironment  # noqa: E402
from analysis.dynamic.experiment_utils import confidence_interval  # noqa: E402
from analysis.dynamic.evaluation import run_policy_episodes  # noqa: E402
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import RiskEstimate, patient_from_row  # noqa: E402
from analysis.dynamic.search import stochastic_mcts_search  # noqa: E402
from analysis.mcts.environment import all_plans  # noqa: E402
from analysis.mcts.outcome_model import (  # noqa: E402
    RegularizedCoxRewardModel,
    prepare_model_cohort,
    stratified_train_validation_test_split,
    tune_penalizer,
)

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_poc_v0_2.json"
REPORT_DIR = ROOT / "reports" / "robustness-v0.3"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-24"
BASE_SEED = 20_260_720
PENALIZERS = (0.01, 0.1, 1.0)
PATIENTS_PER_SUBTYPE = 3          # 12 patients keeps the multi-seed cost tractable
N_SEEDS = 20                      # independent replications
EPISODES_PER_POLICY = 50
SIMULATIONS = 256
EXPLORATION_WEIGHT = math.sqrt(2.0)
SUBTYPES = ("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT, capture_output=True, check=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def make_risk_table(patient_row, os_model, rfs_model):
    plans = all_plans()
    os_scores = os_model.score_plans(patient_row, plans, months=60.0)
    rfs_scores = rfs_model.score_plans(patient_row, plans, months=60.0)
    return {
        plan: RiskEstimate(
            five_year_os=os_scores[plan],
            five_year_rfs=rfs_scores[plan],
        )
        for plan in plans
    }


def balanced_sample(test: pd.DataFrame) -> pd.DataFrame:
    required = ["tumor_size_mm", "lymph_pos", "stage", "grade", "er", "pr", "her2"]
    complete = test.dropna(subset=required).copy()
    sampled = []
    for index, subtype in enumerate(SUBTYPES):
        group = complete[complete["subtype"].eq(subtype)]
        count = min(PATIENTS_PER_SUBTYPE, len(group))
        sampled.append(group.sample(n=count, random_state=BASE_SEED + index))
    return (
        pd.concat(sampled, ignore_index=True)
        .sort_values(["subtype", "patient_id"])
        .reset_index(drop=True)
    )


def build_models(raw: pd.DataFrame):
    os_cohort = prepare_model_cohort(raw)
    os_train, os_val, os_test, assignments = (
        stratified_train_validation_test_split(os_cohort, seed=BASE_SEED)
    )
    os_pen, _ = tune_penalizer(os_train, os_val, PENALIZERS)
    os_model = RegularizedCoxRewardModel(os_pen).fit(
        pd.concat([os_train, os_val], ignore_index=True)
    )

    rfs_cohort = prepare_model_cohort(
        raw, time_column="rfs_months", event_column="rfs_event"
    )
    split_map = assignments.set_index("patient_id")["split"]
    rfs_split = rfs_cohort["patient_id"].astype(str).map(split_map)
    rfs_train = rfs_cohort[rfs_split.eq("train")].copy()
    rfs_val = rfs_cohort[rfs_split.eq("validation")].copy()
    rfs_pen, _ = tune_penalizer(
        rfs_train, rfs_val, PENALIZERS,
        time_column="rfs_months", event_column="rfs_event",
    )
    rfs_model = RegularizedCoxRewardModel(
        rfs_pen, time_column="rfs_months", event_column="rfs_event",
    ).fit(pd.concat([rfs_train, rfs_val], ignore_index=True))
    return os_model, rfs_model, os_test


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_dynamic_config(CONFIG_PATH)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_models(raw)

    sample = balanced_sample(os_test)
    print(f"robustness sample = {len(sample)} patients "
          f"({sample['subtype'].value_counts().to_dict()})")

    # Pre-build one environment per patient (risk tables are seed-independent).
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

    per_seed_rows: list[dict] = []
    first_action_rows: list[dict] = []

    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        mcts_util, nccn_util = [], []
        mcts_surv, nccn_surv = [], []
        mcts_recur, nccn_recur = [], []
        mcts_tox, nccn_tox = [], []
        action_flags = {k: [] for k in (
            "neoadjuvant", "bcs", "intensified", "extended", "regional")}

        for pid, environment in environments.items():
            pid_hash = int(hashlib.sha256(pid.encode("utf-8")).hexdigest()[:8], 16)
            patient_seed = seed + pid_hash % 100_000

            # First-decision action under this seed (root search).
            root = stochastic_mcts_search(
                environment, environment.initial_state(),
                simulations=SIMULATIONS, exploration_weight=EXPLORATION_WEIGHT,
                seed=patient_seed + 7_777,
            )
            first_action_rows.append({
                "seed_index": seed_index, "patient_id": pid,
                "subtype": patient_meta[pid]["subtype"],
                "first_action": root.action,
            })

            mcts_policy = CachedMCTSPolicy(
                environment, simulations=SIMULATIONS,
                exploration_weight=EXPLORATION_WEIGHT, seed=patient_seed,
            )
            nccn_policy = DynamicNccnPolicy(environment)

            mcts_eps = run_policy_episodes(
                environment, mcts_policy, EPISODES_PER_POLICY, patient_seed + 10_000)
            nccn_eps = run_policy_episodes(
                environment, nccn_policy, EPISODES_PER_POLICY, patient_seed + 10_000)

            md = pd.DataFrame(mcts_eps)
            nd = pd.DataFrame(nccn_eps)
            mcts_util.append(md["utility"].mean())
            nccn_util.append(nd["utility"].mean())
            mcts_surv.append(md["survived_5y"].mean())
            nccn_surv.append(nd["survived_5y"].mean())
            mcts_recur.append(md["recurred_by_5y"].mean())
            nccn_recur.append(nd["recurred_by_5y"].mean())
            mcts_tox.append(md["toxicity_count"].mean())
            nccn_tox.append(nd["toxicity_count"].mean())
            action_flags["neoadjuvant"].append((md["timing"] == "neoadjuvant").mean())
            action_flags["bcs"].append((md["surgery"] == "BCS").mean())
            action_flags["intensified"].append((md["chemo"] == "intensified").mean())
            action_flags["extended"].append((md["endocrine"] == "extended").mean())
            action_flags["regional"].append((md["radiation"] == "regional").mean())

        per_seed_rows.append({
            "seed_index": seed_index,
            "mcts_utility": float(np.mean(mcts_util)),
            "nccn_utility": float(np.mean(nccn_util)),
            "utility_gap": float(np.mean(mcts_util) - np.mean(nccn_util)),
            "mcts_survival_pct": float(np.mean(mcts_surv) * 100),
            "nccn_survival_pct": float(np.mean(nccn_surv) * 100),
            "survival_gap_pp": float((np.mean(mcts_surv) - np.mean(nccn_surv)) * 100),
            "mcts_recurrence_pct": float(np.mean(mcts_recur) * 100),
            "nccn_recurrence_pct": float(np.mean(nccn_recur) * 100),
            "mcts_toxicity": float(np.mean(mcts_tox)),
            "nccn_toxicity": float(np.mean(nccn_tox)),
            "neoadjuvant_rate_pct": float(np.mean(action_flags["neoadjuvant"]) * 100),
            "bcs_rate_pct": float(np.mean(action_flags["bcs"]) * 100),
            "intensified_rate_pct": float(np.mean(action_flags["intensified"]) * 100),
            "extended_rate_pct": float(np.mean(action_flags["extended"]) * 100),
            "regional_rate_pct": float(np.mean(action_flags["regional"]) * 100),
        })
        print(f"seed {seed_index + 1}/{N_SEEDS} done "
              f"(gap={per_seed_rows[-1]['utility_gap']:+.4f})")

    per_seed = pd.DataFrame(per_seed_rows)

    # 95% CIs for the headline metrics across seeds.
    ci_metrics = [
        "utility_gap", "survival_gap_pp", "mcts_utility", "nccn_utility",
        "mcts_survival_pct", "nccn_survival_pct", "mcts_recurrence_pct",
        "nccn_recurrence_pct", "neoadjuvant_rate_pct", "bcs_rate_pct",
        "intensified_rate_pct", "extended_rate_pct", "regional_rate_pct",
    ]
    ci_rows = []
    for metric in ci_metrics:
        stats = confidence_interval(per_seed[metric].to_numpy())
        ci_rows.append({"metric": metric, **stats})
    ci_table = pd.DataFrame(ci_rows)

    # Per-patient first-action stability across seeds.
    fa = pd.DataFrame(first_action_rows)
    stability_rows = []
    for pid, group in fa.groupby("patient_id"):
        counts = group["first_action"].value_counts()
        modal = counts.index[0]
        stability_rows.append({
            "patient_id": pid,
            "subtype": patient_meta[pid]["subtype"],
            "modal_first_action": modal,
            "modal_agreement_pct": float(counts.iloc[0] / len(group) * 100),
            "distinct_actions": int(group["first_action"].nunique()),
        })
    stability = pd.DataFrame(stability_rows).sort_values("modal_agreement_pct")

    gap_ci = confidence_interval(per_seed["utility_gap"].to_numpy())
    surv_ci = confidence_interval(per_seed["survival_gap_pp"].to_numpy())
    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "robustness-v0.3",
        "interpretation": (
            "Stochastic stability of the v0.2 comparison across seeds; "
            "same synthetic environment, not causal, not clinical."
        ),
        "design": {
            "patients": int(len(sample)),
            "patients_per_subtype": PATIENTS_PER_SUBTYPE,
            "n_seeds": N_SEEDS,
            "episodes_per_policy_per_seed": EPISODES_PER_POLICY,
            "simulations_per_decision": SIMULATIONS,
        },
        "utility_gap_ci95": gap_ci,
        "survival_gap_pp_ci95": surv_ci,
        "utility_gap_positive_in_all_seeds": bool((per_seed["utility_gap"] > 0).all()),
        "mean_modal_first_action_agreement_pct": float(
            stability["modal_agreement_pct"].mean()
        ),
        "min_modal_first_action_agreement_pct": float(
            stability["modal_agreement_pct"].min()
        ),
    }

    per_seed.to_csv(TABLE_DIR / "per_seed_summary.csv", index=False)
    ci_table.to_csv(TABLE_DIR / "metric_confidence_intervals.csv", index=False)
    stability.to_csv(TABLE_DIR / "first_action_stability.csv", index=False)
    fa.to_csv(TABLE_DIR / "first_action_by_seed.csv", index=False)

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": {
            "data": {"path": "data/processed/patients_with_nccn.csv",
                     "sha256": sha256(INPUT_CSV)},
            "assumptions": {"path": "configs/dynamic_poc_v0_2.json",
                            "sha256": sha256(CONFIG_PATH)},
        },
        "entry_point": "analysis/10_run_multiseed_robustness.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== 95% confidence intervals ===")
    print(ci_table.to_string(index=False))
    print("\n=== first-action stability (worst 5) ===")
    print(stability.head().to_string(index=False))
    print(f"\nutility gap 95% CI: "
          f"[{gap_ci['ci95_low']:+.4f}, {gap_ci['ci95_high']:+.4f}]")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
