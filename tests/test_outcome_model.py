from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from analysis.mcts.environment import all_plans
from analysis.mcts.outcome_model import (
    CoxFeatureEncoder,
    RegularizedCoxRewardModel,
    prepare_model_cohort,
)


class CoxFeatureEncoderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.frame = pd.DataFrame({
            "age": [50.0, 65.0],
            "tumor_size_mm": [20.0, np.nan],
            "lymph_pos": [0.0, 4.0],
            "stage": [1.0, np.nan],
            "grade": [2.0, 3.0],
            "subtype": ["HR+/HER2-", "TNBC"],
            "menopause": ["Pre", "Post"],
            "surgery": ["BCS", "MAST"],
            "chemo": [0, 1],
            "hormone": [1, 0],
            "radio": [1, 1],
        })

    def test_transform_is_numeric_complete_and_stable(self) -> None:
        encoder = CoxFeatureEncoder().fit(self.frame)
        first = encoder.transform(self.frame)
        second = encoder.transform(self.frame.iloc[[0]])
        self.assertFalse(first.isna().any().any())
        self.assertTrue(all(dtype.kind == "f" for dtype in first.dtypes))
        self.assertEqual(list(first.columns), list(second.columns))

    def test_cohort_preparation_accepts_rfs_columns(self) -> None:
        frame = self.frame.assign(
            patient_id=["A", "B"],
            rfs_months=[60.0, 24.0],
            rfs_event=[0, 1],
        )
        prepared = prepare_model_cohort(
            frame,
            time_column="rfs_months",
            event_column="rfs_event",
        )
        self.assertEqual(len(prepared), 2)
        self.assertEqual(prepared["rfs_event"].dtype.kind, "i")


def synthetic_cohort(n: int = 500, seed: int = 3) -> pd.DataFrame:
    """A cohort large enough to fit the Cox reward model on."""
    rng = np.random.default_rng(seed)
    age = rng.normal(60.0, 10.0, n)
    tumor = rng.normal(25.0, 8.0, n).clip(1.0)
    chemo = rng.integers(0, 2, n)
    # Survival driven by the covariates so the fit has something to find; a
    # cohort of pure noise fails to converge and would not exercise anything.
    risk = 0.03 * (age - 60.0) + 0.02 * (tumor - 25.0) + 0.30 * chemo
    frame = pd.DataFrame({
        "patient_id": [f"P{index}" for index in range(n)],
        "age": age,
        "tumor_size_mm": tumor,
        "lymph_pos": rng.integers(0, 5, n).astype(float),
        "stage": rng.integers(1, 5, n).astype(float),
        "grade": rng.integers(1, 4, n).astype(float),
        "subtype": rng.choice(["HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC"], n),
        "menopause": rng.choice(["Pre", "Post"], n),
        "surgery": rng.choice(["BCS", "MAST"], n),
        "chemo": chemo,
        "hormone": rng.integers(0, 2, n),
        "radio": rng.integers(0, 2, n),
        "os_months": rng.exponential(80.0 * np.exp(-risk)).clip(1.0, 180.0),
        "os_event": rng.integers(0, 2, n),
    })
    # The encoder emits a missing-indicator column per partially observed field
    # and a stage-0 dummy. On fully complete synthetic data those columns are
    # constant, which makes the Cox fit diverge - so the fixture carries the same
    # kinds of gaps the real cohort has.
    frame.loc[frame.index[:20], "tumor_size_mm"] = np.nan
    frame.loc[frame.index[20:40], "lymph_pos"] = np.nan
    frame.loc[frame.index[40:60], "stage"] = np.nan
    frame.loc[frame.index[60:80], "grade"] = np.nan
    frame.loc[frame.index[80:110], "stage"] = 0.0
    return prepare_model_cohort(frame)


class NeutraliseTreatmentTermsTests(unittest.TestCase):
    """The v1.4 remedy: the reward model must stop preferring plans on its own."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.cohort = synthetic_cohort()
        cls.model = RegularizedCoxRewardModel(0.5).fit(cls.cohort)
        cls.neutral = cls.model.neutralise_treatment_terms()
        cls.patient = cls.cohort.iloc[0]

    def test_treatment_features_cover_every_decision_and_interaction(self) -> None:
        features = self.model.treatment_features()
        for expected in ("chemo", "hormone", "radio", "surgery_mast",
                         "chemo_x_her2pos", "hormone_x_hrpos", "radio_x_mast",
                         "mast_x_tumor"):
            self.assertIn(expected, features)

    def test_treatment_features_exclude_patient_characteristics(self) -> None:
        features = self.model.treatment_features()
        for patient_term in ("age_z", "tumor_z", "lymph_z", "grade_2",
                             "menopause_post", "subtype_tnbc"):
            self.assertNotIn(patient_term, features)

    def test_neutralised_model_ranks_every_plan_equally(self) -> None:
        scores = self.neutral.score_plans(self.patient, list(all_plans()))
        self.assertAlmostEqual(max(scores.values()), min(scores.values()), places=12)

    def test_the_original_model_does_separate_plans(self) -> None:
        """Guards against a vacuous pass if the fit were degenerate."""
        scores = self.model.score_plans(self.patient, list(all_plans()))
        self.assertGreater(max(scores.values()) - min(scores.values()), 0.0)

    def test_neutralising_does_not_mutate_the_original(self) -> None:
        self.assertNotEqual(float(self.model.model.params_["chemo"]), 0.0)

    def test_patient_characteristic_coefficients_are_untouched(self) -> None:
        for term in ("age_z", "tumor_z", "lymph_z"):
            self.assertAlmostEqual(
                float(self.model.model.params_[term]),
                float(self.neutral.model.params_[term]), places=12)

    def test_patients_are_still_ranked_by_risk(self) -> None:
        """Neutralising treatment must not turn the model into a constant."""
        features = self.neutral.encoder.transform(self.cohort)
        risk = self.neutral.model.predict_partial_hazard(features).to_numpy()
        self.assertGreater(risk.std(), 0.0)


if __name__ == "__main__":
    unittest.main()
