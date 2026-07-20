from __future__ import annotations

import unittest

import pandas as pd

from analysis.nccn_policy import nccn_plan


class SimplifiedNccnPolicyTests(unittest.TestCase):
    def test_hr_positive_low_risk_plan(self) -> None:
        patient = pd.Series({
            "tumor_size_mm": 15,
            "lymph_pos": 0,
            "stage": 1,
            "grade": 1,
            "subtype": "HR+/HER2-",
            "er": 1,
            "pr": 1,
        })
        self.assertEqual(nccn_plan(patient), ("BCS", 0, 1, 1))

    def test_tnbc_stage_three_plan(self) -> None:
        patient = pd.Series({
            "tumor_size_mm": 55,
            "lymph_pos": 5,
            "stage": 3,
            "grade": 3,
            "subtype": "TNBC",
            "er": 0,
            "pr": 0,
        })
        self.assertEqual(nccn_plan(patient), ("MAST", 1, 0, 1))

    def test_missing_stage_returns_no_complete_plan(self) -> None:
        patient = pd.Series({
            "tumor_size_mm": 20,
            "lymph_pos": 0,
            "stage": pd.NA,
            "grade": 2,
            "subtype": "HR+/HER2-",
            "er": 1,
            "pr": 1,
        })
        self.assertIsNone(nccn_plan(patient))


if __name__ == "__main__":
    unittest.main()

