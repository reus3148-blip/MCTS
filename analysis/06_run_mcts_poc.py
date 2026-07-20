"""Train the v0.1 survival reward model and evaluate UCT-MCTS policies.

Outputs are written to ``reports/mcts-poc-v1/tables`` and include the frozen
cohort split, model diagnostics, patient-level policy decisions, convergence
checks against exhaustive search, and subgroup summaries.

Important: the Cox reward is predictive/associational. Counterfactual policy
scores in this PoC are not estimates of causal treatment effects.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from lifelines import KaplanMeierFitter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.mcts.environment import (  # noqa: E402
    DECISIONS,
    TreatmentPlanningEnvironment,
    all_plans,
    feasible_plans_for_subtype,
    plan_to_label,
)
from analysis.mcts.outcome_model import (  # noqa: E402
    EVENT_COLUMN,
    TIME_COLUMN,
    RegularizedCoxRewardModel,
    prepare_model_cohort,
    stratified_train_validation_test_split,
    tune_penalizer,
)
from analysis.mcts.search import exhaustive_best_plan, mcts_plan  # noqa: E402
from analysis.nccn_policy import nccn_plan  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
REPORT_DIR = ROOT / "reports" / "mcts-poc-v1"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-07-11"
SEED = 20_260_711
PENALIZER_CANDIDATES = (0.01, 0.1, 1.0)
SIMULATION_BUDGETS = (16, 32, 64, 128, 256, 512)
PRIMARY_SIMULATION_BUDGET = 512
EXPLORATION_WEIGHT = math.sqrt(2.0)
SURVIVAL_HORIZON_MONTHS = 60.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            capture_output=True,
            check=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def actual_plan(row: pd.Series) -> tuple[object, ...]:
    return (
        str(row["surgery"]),
        int(row["chemo"]),
        int(row["hormone"]),
        int(row["radio"]),
    )


def add_plan_columns(
    output: dict[str, object],
    prefix: str,
    plan: tuple[object, ...] | None,
) -> None:
    for decision, value in zip(DECISIONS, plan or (None,) * 4, strict=True):
        output[f"{prefix}_{decision}"] = value
    output[f"{prefix}_plan"] = plan_to_label(plan) if plan else None


def observed_km_survival(frame: pd.DataFrame, months: float) -> float:
    fitter = KaplanMeierFitter()
    fitter.fit(frame[TIME_COLUMN], event_observed=frame[EVENT_COLUMN])
    return float(fitter.survival_function_at_times(months).iloc[0])


def policy_summary(evaluation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, prefix, prediction_column in [
        ("Actual", "actual", "predicted_5y_os_actual"),
        ("NCCN", "nccn", "predicted_5y_os_nccn"),
        ("MCTS", "mcts", "predicted_5y_os_mcts"),
    ]:
        rows.append({
            "policy": label,
            "n": int(len(evaluation)),
            "mean_predicted_5y_os": float(evaluation[prediction_column].mean()),
            "median_predicted_5y_os": float(evaluation[prediction_column].median()),
            "bcs_rate_pct": float(
                evaluation[f"{prefix}_surgery"].eq("BCS").mean() * 100
            ),
            "chemo_rate_pct": float(
                pd.to_numeric(evaluation[f"{prefix}_chemo"]).mean() * 100
            ),
            "hormone_rate_pct": float(
                pd.to_numeric(evaluation[f"{prefix}_hormone"]).mean() * 100
            ),
            "radio_rate_pct": float(
                pd.to_numeric(evaluation[f"{prefix}_radio"]).mean() * 100
            ),
        })
    return pd.DataFrame(rows)


def policy_agreement(evaluation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, left, right in [
        ("MCTS vs NCCN", "mcts", "nccn"),
        ("MCTS vs Actual", "mcts", "actual"),
        ("NCCN vs Actual", "nccn", "actual"),
    ]:
        decision_matches = []
        for decision in DECISIONS:
            matched = evaluation[f"{left}_{decision}"].eq(
                evaluation[f"{right}_{decision}"]
            )
            decision_matches.append(matched)
            rows.append({
                "comparison": label,
                "decision": decision,
                "n": int(len(evaluation)),
                "agreement_pct": float(matched.mean() * 100),
            })
        all_matched = pd.concat(decision_matches, axis=1).all(axis=1)
        rows.append({
            "comparison": label,
            "decision": "all_four",
            "n": int(len(evaluation)),
            "agreement_pct": float(all_matched.mean() * 100),
        })
    return pd.DataFrame(rows)


def subtype_summary(evaluation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subtype, group in evaluation.groupby("subtype", sort=True):
        row: dict[str, object] = {
            "subtype": subtype,
            "n": int(len(group)),
            "mcts_nccn_all_four_agreement_pct": float(
                group["mcts_plan"].eq(group["nccn_plan"]).mean() * 100
            ),
            "mean_predicted_5y_os_nccn": float(
                group["predicted_5y_os_nccn"].mean()
            ),
            "mean_predicted_5y_os_mcts": float(
                group["predicted_5y_os_mcts"].mean()
            ),
            "mean_predicted_difference_mcts_minus_nccn": float(
                group["predicted_5y_os_mcts_minus_nccn"].mean()
            ),
        }
        for decision in DECISIONS:
            row[f"{decision}_agreement_pct"] = float(
                group[f"mcts_{decision}"].eq(
                    group[f"nccn_{decision}"]
                ).mean() * 100
            )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(INPUT_CSV)
    cohort = prepare_model_cohort(raw)
    train, validation, test, split_assignments = (
        stratified_train_validation_test_split(cohort, seed=SEED)
    )

    print(
        f"cohort={len(cohort)} train={len(train)} "
        f"validation={len(validation)} test={len(test)}"
    )
    selected_penalizer, tuning_scores = tune_penalizer(
        train,
        validation,
        PENALIZER_CANDIDATES,
    )
    train_validation = pd.concat([train, validation], ignore_index=True)
    reward_model = RegularizedCoxRewardModel(selected_penalizer).fit(
        train_validation
    )
    print(f"selected penalizer={selected_penalizer}")

    candidate_plans = all_plans()
    support_counts = Counter(
        actual_plan(row) for _, row in train_validation.iterrows()
    )
    patient_rows: list[dict[str, object]] = []
    convergence_rows: list[dict[str, object]] = []

    for patient_number, (_, patient) in enumerate(test.iterrows()):
        all_rewards = reward_model.score_plans(
            patient,
            candidate_plans,
            months=SURVIVAL_HORIZON_MONTHS,
        )
        feasible_plans = feasible_plans_for_subtype(str(patient["subtype"]))
        rewards = {plan: all_rewards[plan] for plan in feasible_plans}
        environment = TreatmentPlanningEnvironment(rewards)
        exact_plan, exact_reward = exhaustive_best_plan(environment)
        actual = actual_plan(patient)
        guideline = nccn_plan(patient)
        plans_by_budget: dict[int, tuple[object, ...]] = {}

        for budget in SIMULATION_BUDGETS:
            searched = mcts_plan(
                environment,
                simulations_per_step=budget,
                exploration_weight=EXPLORATION_WEIGHT,
                seed=SEED + patient_number * 10_000 + budget,
            )
            plans_by_budget[budget] = searched
            convergence_rows.append({
                "patient_id": patient["patient_id"],
                "subtype": patient["subtype"],
                "simulations_per_step": budget,
                "exact_plan_match": int(searched == exact_plan),
                "regret": float(exact_reward - rewards[searched]),
            })

        searched = plans_by_budget[PRIMARY_SIMULATION_BUDGET]
        output: dict[str, object] = {
            "patient_id": patient["patient_id"],
            "subtype": patient["subtype"],
            "age": patient["age"],
            "stage": patient["stage"],
            "grade": patient["grade"],
            "tumor_size_mm": patient["tumor_size_mm"],
            "lymph_pos": patient["lymph_pos"],
            "nccn_complete": int(guideline is not None),
            "feasible_plan_count": len(feasible_plans),
            "predicted_5y_os_actual": all_rewards[actual],
            "predicted_5y_os_nccn": all_rewards[guideline] if guideline else np.nan,
            "predicted_5y_os_mcts": rewards[searched],
            "predicted_5y_os_exact": exact_reward,
            "predicted_5y_os_mcts_minus_nccn": (
                rewards[searched] - all_rewards[guideline]
                if guideline
                else np.nan
            ),
            "mcts_exact_match": int(searched == exact_plan),
            "mcts_regret": float(exact_reward - rewards[searched]),
            "actual_plan_support_n": int(support_counts[actual]),
            "nccn_plan_support_n": (
                int(support_counts[guideline]) if guideline else np.nan
            ),
            "mcts_plan_support_n": int(support_counts[searched]),
        }
        add_plan_columns(output, "actual", actual)
        add_plan_columns(output, "nccn", guideline)
        add_plan_columns(output, "mcts", searched)
        add_plan_columns(output, "exact", exact_plan)
        patient_rows.append(output)

    patient_decisions = pd.DataFrame(patient_rows)
    convergence_detail = pd.DataFrame(convergence_rows)
    convergence_summary = (
        convergence_detail.groupby("simulations_per_step", as_index=False)
        .agg(
            n=("patient_id", "size"),
            exact_plan_match_pct=("exact_plan_match", lambda values: values.mean() * 100),
            mean_regret=("regret", "mean"),
            max_regret=("regret", "max"),
        )
    )

    evaluation = patient_decisions[patient_decisions["nccn_complete"].eq(1)].copy()
    policy_summary_table = policy_summary(evaluation)
    policy_agreement_table = policy_agreement(evaluation)
    subtype_summary_table = subtype_summary(evaluation)

    predicted_actual = reward_model.predict_survival_at(
        test,
        months=SURVIVAL_HORIZON_MONTHS,
    )
    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "mcts-poc-v1",
        "interpretation": "predictive association only; not causal treatment effect",
        "cohort": {
            "input_rows": int(len(raw)),
            "model_rows": int(len(cohort)),
            "train_rows": int(len(train)),
            "validation_rows": int(len(validation)),
            "test_rows": int(len(test)),
            "test_rows_with_complete_nccn": int(len(evaluation)),
        },
        "reward_model": {
            "type": "L2-regularized Cox proportional hazards",
            "horizon_months": SURVIVAL_HORIZON_MONTHS,
            "penalizer_candidates": list(PENALIZER_CANDIDATES),
            "selected_penalizer": selected_penalizer,
            "train_validation_c_index": reward_model.concordance(train_validation),
            "held_out_test_c_index": reward_model.concordance(test),
            "held_out_observed_km_5y_os": observed_km_survival(
                test,
                SURVIVAL_HORIZON_MONTHS,
            ),
            "held_out_mean_predicted_5y_os_actual_treatment": float(
                predicted_actual.mean()
            ),
            "encoder": reward_model.encoder.metadata(),
        },
        "search": {
            "algorithm": "UCT with random rollouts and receding horizon",
            "exploration_weight": EXPLORATION_WEIGHT,
            "simulation_budgets_per_step": list(SIMULATION_BUDGETS),
            "primary_simulations_per_step": PRIMARY_SIMULATION_BUDGET,
            "max_action_space_size": len(candidate_plans),
            "min_action_space_size": int(patient_decisions["feasible_plan_count"].min()),
            "seed": SEED,
            "primary_exact_plan_match_pct": float(
                patient_decisions["mcts_exact_match"].mean() * 100
            ),
            "primary_mean_regret": float(patient_decisions["mcts_regret"].mean()),
            "primary_max_regret": float(patient_decisions["mcts_regret"].max()),
        },
    }

    split_assignments.sort_values(["split", "patient_id"]).to_csv(
        TABLE_DIR / "cohort_split.csv",
        index=False,
        encoding="utf-8",
    )
    tuning_scores.to_csv(
        TABLE_DIR / "penalizer_tuning.csv",
        index=False,
        encoding="utf-8",
    )
    reward_model.coefficient_table().to_csv(
        TABLE_DIR / "cox_coefficients.csv",
        index=False,
        encoding="utf-8",
    )
    patient_decisions.to_csv(
        TABLE_DIR / "patient_policy_decisions.csv",
        index=False,
        encoding="utf-8",
    )
    convergence_detail.to_csv(
        TABLE_DIR / "search_convergence_detail.csv",
        index=False,
        encoding="utf-8",
    )
    convergence_summary.to_csv(
        TABLE_DIR / "search_convergence.csv",
        index=False,
        encoding="utf-8",
    )
    policy_summary_table.to_csv(
        TABLE_DIR / "policy_summary.csv",
        index=False,
        encoding="utf-8",
    )
    policy_agreement_table.to_csv(
        TABLE_DIR / "policy_agreement.csv",
        index=False,
        encoding="utf-8",
    )
    subtype_summary_table.to_csv(
        TABLE_DIR / "subtype_summary.csv",
        index=False,
        encoding="utf-8",
    )

    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "input": {
            "path": str(INPUT_CSV.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(INPUT_CSV),
        },
        "entry_point": "analysis/06_run_mcts_poc.py",
        "seed": SEED,
        "outputs": sorted(
            str(path.relative_to(ROOT)).replace("\\", "/")
            for path in TABLE_DIR.glob("*.csv")
        ),
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== held-out model metrics ===")
    print(json.dumps(metrics["reward_model"], ensure_ascii=False, indent=2))
    print("\n=== MCTS convergence ===")
    print(convergence_summary.to_string(index=False))
    print("\n=== policy summary (complete NCCN subset) ===")
    print(policy_summary_table.to_string(index=False))
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
