"""Shared cohort, reward-model and run-manifest helpers for the v0.3+ studies.

``analysis/10``, ``analysis/11`` and ``analysis/12`` all need the same three
things before they can run anything: fitted OS/RFS Cox reward models, a
subtype-balanced sample drawn from the *held-out* test split, and per-patient
risk tables. They were copied into each script; keeping one implementation here
means the three studies are guaranteed to evaluate the same patients under the
same models, which is what makes their numbers comparable at all.

The numbered ``analysis/1x_*.py`` scripts cannot be imported (module names may
not start with a digit), so this module is also the only place these helpers can
be unit-tested from.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
from typing import Mapping, Sequence

import pandas as pd

from analysis.mcts.environment import Plan, all_plans
from analysis.mcts.outcome_model import (
    RegularizedCoxRewardModel,
    prepare_model_cohort,
    stratified_train_validation_test_split,
    tune_penalizer,
)
from .schema import RiskEstimate

ROOT = Path(__file__).resolve().parents[2]

#: Seed shared by every v0.3+ study so they draw the same split and sample.
BASE_SEED = 20_260_720
PENALIZERS: tuple[float, ...] = (0.01, 0.1, 1.0)
SUBTYPES: tuple[str, ...] = ("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC")

#: Columns a patient must have before the environment can be built for them.
REQUIRED_PATIENT_COLUMNS = (
    "tumor_size_mm", "lymph_pos", "stage", "grade", "er", "pr", "her2",
)


def sha256(path: str | Path) -> str:
    """SHA-256 of a file, streamed so large inputs do not load into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(cwd: str | Path = ROOT) -> str:
    """Commit the run started from, or ``"unavailable"`` outside a checkout."""
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(cwd), capture_output=True, check=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def build_reward_models(
    raw: pd.DataFrame,
    seed: int = BASE_SEED,
    penalizers: Sequence[float] = PENALIZERS,
) -> tuple[RegularizedCoxRewardModel, RegularizedCoxRewardModel, pd.DataFrame]:
    """Fit the OS and RFS Cox reward models and return the held-out OS test set.

    The RFS model reuses the *OS* train/validation/test assignment so a patient
    never lands in one model's training data and the other's test data.
    """
    os_cohort = prepare_model_cohort(raw)
    os_train, os_val, os_test, assignments = (
        stratified_train_validation_test_split(os_cohort, seed=seed)
    )
    os_penalizer, _ = tune_penalizer(os_train, os_val, penalizers)
    os_model = RegularizedCoxRewardModel(os_penalizer).fit(
        pd.concat([os_train, os_val], ignore_index=True)
    )

    rfs_cohort = prepare_model_cohort(
        raw, time_column="rfs_months", event_column="rfs_event"
    )
    split_map = assignments.set_index("patient_id")["split"]
    rfs_split = rfs_cohort["patient_id"].astype(str).map(split_map)
    rfs_train = rfs_cohort[rfs_split.eq("train")].copy()
    rfs_val = rfs_cohort[rfs_split.eq("validation")].copy()
    rfs_penalizer, _ = tune_penalizer(
        rfs_train, rfs_val, penalizers,
        time_column="rfs_months", event_column="rfs_event",
    )
    rfs_model = RegularizedCoxRewardModel(
        rfs_penalizer, time_column="rfs_months", event_column="rfs_event",
    ).fit(pd.concat([rfs_train, rfs_val], ignore_index=True))
    return os_model, rfs_model, os_test


def balanced_subtype_sample(
    test: pd.DataFrame,
    per_subtype: int,
    seed: int = BASE_SEED,
    subtypes: Sequence[str] = SUBTYPES,
    offset: int = 0,
) -> pd.DataFrame:
    """Draw ``per_subtype`` complete-record patients from each molecular subtype.

    Sampling per subtype (rather than at random) keeps rare subtypes such as
    TNBC represented in the small samples these studies can afford. Each subtype
    uses ``seed + index`` so that raising ``per_subtype`` extends a sample rather
    than replacing it.

    ``offset`` skips that many patients per subtype before taking the sample,
    which is how a *disjoint* second cohort is drawn: ``offset=0`` and
    ``offset=per_subtype`` share the same draw order and therefore share no
    patients. v1.1 could only report nested cohorts, and nested cohorts cannot
    tell a replication from a recount of the same people.
    """
    if per_subtype < 1:
        raise ValueError("per_subtype must be positive")
    if offset < 0:
        raise ValueError("offset must be non-negative")
    complete = test.dropna(subset=list(REQUIRED_PATIENT_COLUMNS)).copy()
    sampled = []
    for index, subtype in enumerate(subtypes):
        group = complete[complete["subtype"].eq(subtype)]
        count = min(per_subtype + offset, len(group))
        if count > offset:
            block = group.sample(n=count, random_state=seed + index)
            sampled.append(block.iloc[offset:])
    if not sampled:
        raise ValueError("no complete patients available for the requested subtypes")
    return (
        pd.concat(sampled, ignore_index=True)
        .sort_values(["subtype", "patient_id"])
        .reset_index(drop=True)
    )


def make_risk_table(
    row: pd.Series,
    os_model: RegularizedCoxRewardModel,
    rfs_model: RegularizedCoxRewardModel,
    months: float = 60.0,
) -> Mapping[Plan, RiskEstimate]:
    """Five-year OS/RFS estimates for every static plan, for one patient."""
    plans = all_plans()
    os_scores = os_model.score_plans(row, plans, months=months)
    rfs_scores = rfs_model.score_plans(row, plans, months=months)
    return {
        plan: RiskEstimate(
            five_year_os=os_scores[plan],
            five_year_rfs=rfs_scores[plan],
        )
        for plan in plans
    }


def input_manifest(paths: Mapping[str, str | Path]) -> dict[str, dict[str, str]]:
    """``{label: {path, sha256}}`` for the run manifests the reports carry."""
    return {
        label: {
            "path": str(Path(path).resolve().relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(path),
        }
        for label, path in paths.items()
    }
