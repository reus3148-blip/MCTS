from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from analysis.mcts.outcome_model import (
    CoxFeatureEncoder,
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


if __name__ == "__main__":
    unittest.main()
