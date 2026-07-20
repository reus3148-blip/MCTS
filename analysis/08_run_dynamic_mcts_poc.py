"""Run the METABRIC-backed stochastic dynamic MCTS proof of concept.

Observed METABRIC data fit the five-year OS and RFS baseline risk models.
Treatment response, intensity, field, toxicity, and utility weights are
synthetic teaching assumptions loaded from an explicit JSON configuration.
Nothing produced by this script is a clinical treatment recommendation.
"""

from __future__ import annotations

from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.dynamic.config import load_dynamic_config  # noqa: E402
from analysis.dynamic.environment import DynamicBreastCancerEnvironment  # noqa: E402
from analysis.dynamic.evaluation import run_policy_episodes, simulate_episode  # noqa: E402
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
REPORT_DIR = ROOT / "reports" / "dynamic-mcts-poc-v0.2"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-07-11"
SEED = 20_260_711
PENALIZERS = (0.01, 0.1, 1.0)
PATIENTS_PER_SUBTYPE = 10
EPISODES_PER_POLICY = 100
PRIMARY_SIMULATIONS = 256
STABILITY_BUDGETS = (16, 32, 64, 128, 256)
REFERENCE_SIMULATIONS = 1024
EXPLORATION_WEIGHT = math.sqrt(2.0)


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
            cwd=ROOT,
            capture_output=True,
            check=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def split_from_assignments(
    cohort: pd.DataFrame,
    assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    split_map = assignments.set_index("patient_id")["split"]
    split = cohort["patient_id"].astype(str).map(split_map)
    return (
        cohort[split.eq("train")].copy(),
        cohort[split.eq("validation")].copy(),
        cohort[split.eq("test")].copy(),
    )


def balanced_dynamic_sample(test: pd.DataFrame) -> pd.DataFrame:
    required = [
        "tumor_size_mm",
        "lymph_pos",
        "stage",
        "grade",
        "er",
        "pr",
        "her2",
    ]
    complete = test.dropna(subset=required).copy()
    sampled = []
    for index, subtype in enumerate([
        "HR+/HER2-",
        "HR+/HER2+",
        "HR-/HER2+",
        "TNBC",
    ]):
        group = complete[complete["subtype"].eq(subtype)]
        count = min(PATIENTS_PER_SUBTYPE, len(group))
        sampled.append(group.sample(n=count, random_state=SEED + index))
    return (
        pd.concat(sampled, ignore_index=True)
        .sort_values(["subtype", "patient_id"])
        .reset_index(drop=True)
    )


def make_risk_table(
    patient: pd.Series,
    os_model: RegularizedCoxRewardModel,
    rfs_model: RegularizedCoxRewardModel,
) -> dict[tuple[object, ...], RiskEstimate]:
    plans = all_plans()
    os_scores = os_model.score_plans(patient, plans, months=60.0)
    rfs_scores = rfs_model.score_plans(patient, plans, months=60.0)
    return {
        plan: RiskEstimate(
            five_year_os=os_scores[plan],
            five_year_rfs=rfs_scores[plan],
        )
        for plan in plans
    }


def count_max_decision_trajectories(
    environment: DynamicBreastCancerEnvironment,
) -> int:
    patient = environment.patient
    config = environment.config
    endocrine_count = 3 if patient.hr_positive else 1

    def surgery_radiation_paths(tumor_size: float) -> int:
        bcs_allowed = (
            tumor_size <= float(config.eligibility["bcs_max_tumor_mm"])
            and patient.stage <= int(config.eligibility["bcs_max_stage"])
        )
        return 2 + (3 if bcs_allowed else 0)

    surgery_first = 3 * endocrine_count * surgery_radiation_paths(
        patient.tumor_size_mm
    )
    initial_actions = environment.legal_actions(environment.initial_state())
    if "neoadjuvant" not in initial_actions:
        return surgery_first

    neoadjuvant = 0
    for _chemo in ("standard", "intensified"):
        for response, multiplier in config.tumor_size_multipliers.items():
            if response == "not_applicable":
                continue
            tumor_size = patient.tumor_size_mm * float(multiplier)
            neoadjuvant += endocrine_count * surgery_radiation_paths(tumor_size)
    return surgery_first + neoadjuvant


def summarize_episodes(episodes: pd.DataFrame) -> pd.DataFrame:
    return (
        episodes.groupby("policy", as_index=False)
        .agg(
            patients=("patient_id", "nunique"),
            episodes=("episode", "size"),
            mean_utility=("utility", "mean"),
            survived_5y_pct=("survived_5y", lambda values: values.mean() * 100),
            recurred_by_5y_pct=(
                "recurred_by_5y",
                lambda values: values.mean() * 100,
            ),
            mean_toxicity_count=("toxicity_count", "mean"),
            neoadjuvant_rate_pct=(
                "timing",
                lambda values: values.eq("neoadjuvant").mean() * 100,
            ),
            bcs_rate_pct=(
                "surgery",
                lambda values: values.eq("BCS").mean() * 100,
            ),
            intensified_chemo_rate_pct=(
                "chemo",
                lambda values: values.eq("intensified").mean() * 100,
            ),
            extended_endocrine_rate_pct=(
                "endocrine",
                lambda values: values.eq("extended").mean() * 100,
            ),
            regional_radiation_rate_pct=(
                "radiation",
                lambda values: values.eq("regional").mean() * 100,
            ),
        )
    )


def summarize_by_subtype(episodes: pd.DataFrame) -> pd.DataFrame:
    return (
        episodes.groupby(["policy", "subtype"], as_index=False)
        .agg(
            patients=("patient_id", "nunique"),
            episodes=("episode", "size"),
            mean_utility=("utility", "mean"),
            survived_5y_pct=("survived_5y", lambda values: values.mean() * 100),
            recurred_by_5y_pct=(
                "recurred_by_5y",
                lambda values: values.mean() * 100,
            ),
            mean_toxicity_count=("toxicity_count", "mean"),
            neoadjuvant_rate_pct=(
                "timing",
                lambda values: values.eq("neoadjuvant").mean() * 100,
            ),
        )
    )


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_dynamic_config(CONFIG_PATH)
    raw = pd.read_csv(INPUT_CSV)

    os_cohort = prepare_model_cohort(raw)
    os_train, os_validation, os_test, assignments = (
        stratified_train_validation_test_split(os_cohort, seed=SEED)
    )
    os_penalizer, os_tuning = tune_penalizer(
        os_train,
        os_validation,
        PENALIZERS,
    )
    os_train_validation = pd.concat(
        [os_train, os_validation],
        ignore_index=True,
    )
    os_model = RegularizedCoxRewardModel(os_penalizer).fit(os_train_validation)

    rfs_cohort = prepare_model_cohort(
        raw,
        time_column="rfs_months",
        event_column="rfs_event",
    )
    rfs_train, rfs_validation, rfs_test = split_from_assignments(
        rfs_cohort,
        assignments,
    )
    rfs_penalizer, rfs_tuning = tune_penalizer(
        rfs_train,
        rfs_validation,
        PENALIZERS,
        time_column="rfs_months",
        event_column="rfs_event",
    )
    rfs_train_validation = pd.concat(
        [rfs_train, rfs_validation],
        ignore_index=True,
    )
    rfs_model = RegularizedCoxRewardModel(
        rfs_penalizer,
        time_column="rfs_months",
        event_column="rfs_event",
    ).fit(rfs_train_validation)

    sample = balanced_dynamic_sample(os_test)
    print(
        f"dynamic sample={len(sample)} "
        f"({sample['subtype'].value_counts().to_dict()})"
    )

    episode_rows: list[dict[str, object]] = []
    patient_rows: list[dict[str, object]] = []
    stability_rows: list[dict[str, object]] = []
    trace_rows: list[dict[str, object]] = []
    traced_subtypes: set[str] = set()

    for patient_index, patient_row in sample.iterrows():
        patient = patient_from_row(patient_row)
        environment = DynamicBreastCancerEnvironment(
            patient,
            make_risk_table(patient_row, os_model, rfs_model),
            config,
        )
        patient_seed = SEED + patient_index * 100_000
        path_count = count_max_decision_trajectories(environment)
        patient_rows.append({
            "patient_id": patient.patient_id,
            "subtype": patient.subtype,
            "stage": patient.stage,
            "tumor_size_mm": patient.tumor_size_mm,
            "lymph_pos": patient.lymph_pos,
            "max_decision_trajectories": path_count,
        })

        initial_state = environment.initial_state()
        reference = stochastic_mcts_search(
            environment,
            initial_state,
            simulations=REFERENCE_SIMULATIONS,
            exploration_weight=EXPLORATION_WEIGHT,
            seed=patient_seed + 7_777,
        )
        ordered_reference_values = sorted(
            reference.action_values.values(),
            reverse=True,
        )
        reference_margin = (
            ordered_reference_values[0] - ordered_reference_values[1]
            if len(ordered_reference_values) > 1
            else 0.0
        )
        reference_near_tie = reference_margin <= 0.01
        for budget in STABILITY_BUDGETS:
            searched = stochastic_mcts_search(
                environment,
                initial_state,
                simulations=budget,
                exploration_weight=EXPLORATION_WEIGHT,
                seed=patient_seed + 7_777,
            )
            stability_rows.append({
                "patient_id": patient.patient_id,
                "subtype": patient.subtype,
                "simulations": budget,
                "selected_action": searched.action,
                "reference_action": reference.action,
                "reference_match": int(searched.action == reference.action),
                "reference_value_margin": reference_margin,
                "reference_near_tie": int(reference_near_tie),
                "match_or_near_tie": int(
                    searched.action == reference.action or reference_near_tie
                ),
            })

        policies = {
            "NCCN": DynamicNccnPolicy(environment),
            "MCTS": CachedMCTSPolicy(
                environment,
                simulations=PRIMARY_SIMULATIONS,
                exploration_weight=EXPLORATION_WEIGHT,
                seed=patient_seed,
            ),
        }
        for policy_name, policy in policies.items():
            results = run_policy_episodes(
                environment,
                policy,
                episodes=EPISODES_PER_POLICY,
                seed=patient_seed + 10_000,
            )
            for episode, result in enumerate(results, start=1):
                episode_rows.append({
                    "patient_id": patient.patient_id,
                    "subtype": patient.subtype,
                    "policy": policy_name,
                    "episode": episode,
                    **result,
                })

            if patient.subtype not in traced_subtypes:
                _, trace = simulate_episode(
                    environment,
                    policy,
                    seed=patient_seed + 99_999,
                    include_trace=True,
                )
                for row in trace:
                    trace_rows.append({
                        "patient_id": patient.patient_id,
                        "subtype": patient.subtype,
                        "policy": policy_name,
                        **row,
                    })
        traced_subtypes.add(patient.subtype)

    episodes = pd.DataFrame(episode_rows)
    patients = pd.DataFrame(patient_rows)
    stability_detail = pd.DataFrame(stability_rows)
    stability = (
        stability_detail.groupby("simulations", as_index=False)
        .agg(
            patients=("patient_id", "size"),
            reference_action_match_pct=(
                "reference_match",
                lambda values: values.mean() * 100,
            ),
            reference_near_tie_pct=(
                "reference_near_tie",
                lambda values: values.mean() * 100,
            ),
            match_or_near_tie_pct=(
                "match_or_near_tie",
                lambda values: values.mean() * 100,
            ),
        )
    )
    policy_summary = summarize_episodes(episodes)
    subtype_summary = summarize_by_subtype(episodes)
    patient_policy_summary = (
        episodes.groupby(["patient_id", "subtype", "policy"], as_index=False)
        .agg(
            mean_utility=("utility", "mean"),
            survived_5y_pct=("survived_5y", lambda values: values.mean() * 100),
            recurred_by_5y_pct=(
                "recurred_by_5y",
                lambda values: values.mean() * 100,
            ),
            mean_toxicity_count=("toxicity_count", "mean"),
        )
    )

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "dynamic-mcts-poc-v0.2",
        "interpretation": (
            "METABRIC-backed baseline risk plus synthetic transitions; "
            "not causal and not clinical"
        ),
        "cohort": {
            "input_rows": int(len(raw)),
            "os_model_rows": int(len(os_cohort)),
            "rfs_model_rows": int(len(rfs_cohort)),
            "dynamic_test_patients": int(len(sample)),
            "episodes_per_policy_patient": EPISODES_PER_POLICY,
        },
        "models": {
            "os_selected_penalizer": os_penalizer,
            "os_test_c_index": os_model.concordance(os_test),
            "rfs_selected_penalizer": rfs_penalizer,
            "rfs_test_c_index": rfs_model.concordance(rfs_test),
        },
        "environment": {
            "horizon_years": config.horizon_years,
            "minimum_max_decision_trajectories": int(
                patients["max_decision_trajectories"].min()
            ),
            "maximum_max_decision_trajectories": int(
                patients["max_decision_trajectories"].max()
            ),
            "assumption_status": config.assumption_status,
        },
        "search": {
            "primary_simulations_per_decision": PRIMARY_SIMULATIONS,
            "reference_simulations": REFERENCE_SIMULATIONS,
            "exploration_weight": EXPLORATION_WEIGHT,
            "seed": SEED,
        },
        "policy_summary": policy_summary.to_dict(orient="records"),
        "search_stability": stability.to_dict(orient="records"),
    }

    os_tuning.to_csv(TABLE_DIR / "os_penalizer_tuning.csv", index=False)
    rfs_tuning.to_csv(TABLE_DIR / "rfs_penalizer_tuning.csv", index=False)
    os_model.coefficient_table().to_csv(
        TABLE_DIR / "os_cox_coefficients.csv",
        index=False,
    )
    rfs_model.coefficient_table().to_csv(
        TABLE_DIR / "rfs_cox_coefficients.csv",
        index=False,
    )
    patients.to_csv(TABLE_DIR / "dynamic_patients.csv", index=False)
    episodes.to_csv(TABLE_DIR / "policy_episodes.csv", index=False)
    policy_summary.to_csv(TABLE_DIR / "policy_summary.csv", index=False)
    subtype_summary.to_csv(TABLE_DIR / "subtype_summary.csv", index=False)
    patient_policy_summary.to_csv(
        TABLE_DIR / "patient_policy_summary.csv",
        index=False,
    )
    stability_detail.to_csv(
        TABLE_DIR / "search_stability_detail.csv",
        index=False,
    )
    stability.to_csv(TABLE_DIR / "search_stability.csv", index=False)
    pd.DataFrame(trace_rows).to_csv(
        TABLE_DIR / "example_traces.csv",
        index=False,
    )

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    shutil.copy2(CONFIG_PATH, REPORT_DIR / "assumptions_snapshot.json")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": {
            "data": {
                "path": str(INPUT_CSV.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(INPUT_CSV),
            },
            "assumptions": {
                "path": str(CONFIG_PATH.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(CONFIG_PATH),
            },
        },
        "entry_point": "analysis/08_run_dynamic_mcts_poc.py",
        "seed": SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== model metrics ===")
    print(json.dumps(metrics["models"], indent=2))
    print("\n=== environment ===")
    print(json.dumps(metrics["environment"], ensure_ascii=False, indent=2))
    print("\n=== search stability ===")
    print(stability.to_string(index=False))
    print("\n=== policy summary ===")
    print(policy_summary.to_string(index=False))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
