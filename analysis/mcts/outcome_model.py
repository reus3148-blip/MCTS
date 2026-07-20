"""Regularized Cox model used as the v0.1 terminal reward model.

The model estimates association-adjusted five-year overall survival. It does
not identify causal treatment effects from the observational METABRIC cohort.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import numpy as np
import pandas as pd

from analysis.nccn_policy import normalize_actual_surgery
from .environment import Plan, plan_to_dict

TIME_COLUMN = "os_months"
EVENT_COLUMN = "os_event"
TREATMENT_COLUMNS = ("surgery", "chemo", "hormone", "radio")


def prepare_model_cohort(
    df: pd.DataFrame,
    time_column: str = TIME_COLUMN,
    event_column: str = EVENT_COLUMN,
) -> pd.DataFrame:
    """Select rows with survival and all four observed treatment decisions."""
    required = [
        "patient_id",
        "subtype",
        time_column,
        event_column,
        *TREATMENT_COLUMNS,
    ]
    cohort = df.dropna(subset=required).copy()
    cohort["surgery"] = cohort["surgery"].apply(normalize_actual_surgery)
    cohort = cohort.dropna(subset=["surgery"])
    cohort = cohort[cohort[time_column] > 0].copy()
    for column in ("chemo", "hormone", "radio", event_column):
        cohort[column] = cohort[column].astype(int)
    return cohort.reset_index(drop=True)


def stratified_train_validation_test_split(
    cohort: pd.DataFrame,
    seed: int,
    train_fraction: float = 0.60,
    validation_fraction: float = 0.20,
    event_column: str = EVENT_COLUMN,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Split by subtype and event with stable patient-level assignments."""
    if train_fraction <= 0 or validation_fraction <= 0:
        raise ValueError("split fractions must be positive")
    if train_fraction + validation_fraction >= 1:
        raise ValueError("train and validation fractions must sum to < 1")

    rng = np.random.default_rng(seed)
    assignments: list[dict[str, str]] = []
    strata = cohort["subtype"].astype(str) + "|" + cohort[event_column].astype(str)

    for stratum in sorted(strata.unique()):
        indices = cohort.index[strata == stratum].to_numpy(copy=True)
        rng.shuffle(indices)
        count = len(indices)
        train_end = max(1, int(round(count * train_fraction)))
        validation_count = max(1, int(round(count * validation_fraction)))
        validation_end = min(count - 1, train_end + validation_count)
        train_end = min(train_end, validation_end - 1)

        split_indices = {
            "train": indices[:train_end],
            "validation": indices[train_end:validation_end],
            "test": indices[validation_end:],
        }
        for split, selected in split_indices.items():
            for index in selected:
                assignments.append({
                    "patient_id": str(cohort.at[index, "patient_id"]),
                    "split": split,
                    "stratum": stratum,
                })

    assignment_frame = pd.DataFrame(assignments)
    split_map = assignment_frame.set_index("patient_id")["split"]
    patient_ids = cohort["patient_id"].astype(str)
    train = cohort[patient_ids.map(split_map).eq("train")].copy()
    validation = cohort[patient_ids.map(split_map).eq("validation")].copy()
    test = cohort[patient_ids.map(split_map).eq("test")].copy()
    return train, validation, test, assignment_frame


@dataclass
class CoxFeatureEncoder:
    """Fixed, inspectable feature map with train-only imputation/scaling."""

    medians: dict[str, float] | None = None
    means: dict[str, float] | None = None
    standard_deviations: dict[str, float] | None = None

    def fit(self, frame: pd.DataFrame) -> "CoxFeatureEncoder":
        medians = {
            column: float(frame[column].median())
            for column in ("age", "tumor_size_mm", "lymph_pos")
        }
        transformed = self._continuous_values(frame, medians)
        means = {
            column: float(transformed[column].mean())
            for column in transformed
        }
        standard_deviations = {
            column: float(transformed[column].std(ddof=0)) or 1.0
            for column in transformed
        }
        self.medians = medians
        self.means = means
        self.standard_deviations = standard_deviations
        return self

    @staticmethod
    def _continuous_values(
        frame: pd.DataFrame,
        medians: dict[str, float],
    ) -> pd.DataFrame:
        values = pd.DataFrame(index=frame.index)
        values["age"] = frame["age"].fillna(medians["age"]).astype(float)
        values["log_tumor"] = np.log1p(
            frame["tumor_size_mm"].fillna(medians["tumor_size_mm"]).clip(lower=0)
        )
        values["log_lymph"] = np.log1p(
            frame["lymph_pos"].fillna(medians["lymph_pos"]).clip(lower=0)
        )
        return values

    def transform(self, frame: pd.DataFrame) -> pd.DataFrame:
        if self.medians is None or self.means is None:
            raise RuntimeError("encoder must be fitted before transform")
        if self.standard_deviations is None:
            raise RuntimeError("encoder scaling parameters are unavailable")

        continuous = self._continuous_values(frame, self.medians)
        output = pd.DataFrame(index=frame.index)
        for source, destination in [
            ("age", "age_z"),
            ("log_tumor", "tumor_z"),
            ("log_lymph", "lymph_z"),
        ]:
            output[destination] = (
                continuous[source] - self.means[source]
            ) / self.standard_deviations[source]

        output["tumor_missing"] = frame["tumor_size_mm"].isna().astype(float)
        output["lymph_missing"] = frame["lymph_pos"].isna().astype(float)

        subtype = frame["subtype"].astype(str)
        for value, name in [
            ("HR+/HER2+", "subtype_hrpos_her2pos"),
            ("HR-/HER2+", "subtype_hrneg_her2pos"),
            ("TNBC", "subtype_tnbc"),
        ]:
            output[name] = subtype.eq(value).astype(float)

        stage = pd.to_numeric(frame["stage"], errors="coerce")
        for value in (0, 2, 3, 4):
            output[f"stage_{value}"] = stage.eq(value).astype(float)
        output["stage_missing"] = stage.isna().astype(float)

        grade = pd.to_numeric(frame["grade"], errors="coerce")
        for value in (2, 3):
            output[f"grade_{value}"] = grade.eq(value).astype(float)
        output["grade_missing"] = grade.isna().astype(float)
        output["menopause_post"] = (
            frame["menopause"].astype(str).str.lower().eq("post").astype(float)
        )

        surgery_mast = frame["surgery"].astype(str).eq("MAST").astype(float)
        chemo = pd.to_numeric(frame["chemo"], errors="raise").astype(float)
        hormone = pd.to_numeric(frame["hormone"], errors="raise").astype(float)
        radio = pd.to_numeric(frame["radio"], errors="raise").astype(float)
        output["surgery_mast"] = surgery_mast
        output["chemo"] = chemo
        output["hormone"] = hormone
        output["radio"] = radio

        her2_positive = subtype.isin({"HR+/HER2+", "HR-/HER2+"}).astype(float)
        hr_positive = subtype.str.startswith("HR+").astype(float)
        stage_three_plus = stage.ge(3).fillna(False).astype(float)
        output["chemo_x_hrpos_her2neg"] = chemo * subtype.eq("HR+/HER2-")
        output["chemo_x_her2pos"] = chemo * her2_positive
        output["chemo_x_tnbc"] = chemo * subtype.eq("TNBC")
        output["hormone_x_hrpos"] = hormone * hr_positive
        output["radio_x_mast"] = radio * surgery_mast
        output["radio_x_stage3plus"] = radio * stage_three_plus
        output["mast_x_stage3plus"] = surgery_mast * stage_three_plus
        output["mast_x_tumor"] = surgery_mast * output["tumor_z"]

        return output.astype(float)

    def metadata(self) -> dict[str, dict[str, float]]:
        return {
            "medians": self.medians or {},
            "means": self.means or {},
            "standard_deviations": self.standard_deviations or {},
        }


class RegularizedCoxRewardModel:
    """Cox model wrapper that exposes terminal five-year survival rewards."""

    def __init__(
        self,
        penalizer: float,
        time_column: str = TIME_COLUMN,
        event_column: str = EVENT_COLUMN,
    ) -> None:
        self.penalizer = float(penalizer)
        self.time_column = time_column
        self.event_column = event_column
        self.encoder = CoxFeatureEncoder()
        self.model = CoxPHFitter(penalizer=self.penalizer)

    def fit(self, frame: pd.DataFrame) -> "RegularizedCoxRewardModel":
        self.encoder.fit(frame)
        design = self.encoder.transform(frame)
        duration_column = "__duration"
        event_column = "__event"
        design[duration_column] = frame[self.time_column].astype(float).to_numpy()
        design[event_column] = frame[self.event_column].astype(int).to_numpy()
        self.model.fit(
            design,
            duration_col=duration_column,
            event_col=event_column,
            show_progress=False,
        )
        return self

    def concordance(self, frame: pd.DataFrame) -> float:
        features = self.encoder.transform(frame)
        risk = self.model.predict_partial_hazard(features).to_numpy()
        return float(concordance_index(
            frame[self.time_column],
            -risk,
            frame[self.event_column],
        ))

    def predict_survival_at(
        self,
        frame: pd.DataFrame,
        months: float = 60.0,
    ) -> np.ndarray:
        features = self.encoder.transform(frame)
        predictions = self.model.predict_survival_function(
            features,
            times=[float(months)],
        )
        return predictions.iloc[0].to_numpy(dtype=float)

    def score_plans(
        self,
        patient: pd.Series,
        plans: Sequence[Plan],
        months: float = 60.0,
    ) -> dict[Plan, float]:
        rows = []
        base = patient.to_dict()
        for plan in plans:
            row = dict(base)
            row.update(plan_to_dict(plan))
            rows.append(row)
        frame = pd.DataFrame(rows)
        predictions = self.predict_survival_at(frame, months=months)
        return {
            tuple(plan): float(prediction)
            for plan, prediction in zip(plans, predictions, strict=True)
        }

    def coefficient_table(self) -> pd.DataFrame:
        table = self.model.summary.reset_index().rename(
            columns={"covariate": "feature"}
        )
        keep = [
            column
            for column in [
                "feature",
                "coef",
                "exp(coef)",
                "se(coef)",
                "coef lower 95%",
                "coef upper 95%",
                "p",
            ]
            if column in table.columns
        ]
        return table[keep]


def tune_penalizer(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    candidates: Iterable[float],
    time_column: str = TIME_COLUMN,
    event_column: str = EVENT_COLUMN,
) -> tuple[float, pd.DataFrame]:
    """Choose the highest validation C-index, preferring stronger shrinkage."""
    rows = []
    for penalizer in candidates:
        model = RegularizedCoxRewardModel(
            float(penalizer),
            time_column=time_column,
            event_column=event_column,
        ).fit(train)
        rows.append({
            "penalizer": float(penalizer),
            "train_c_index": model.concordance(train),
            "validation_c_index": model.concordance(validation),
        })
    scores = pd.DataFrame(rows).sort_values("penalizer").reset_index(drop=True)
    best = max(
        rows,
        key=lambda row: (row["validation_c_index"], row["penalizer"]),
    )
    return float(best["penalizer"]), scores
