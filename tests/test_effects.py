from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from analysis.causal.decisions import DEFAULT_SPEC, drop_constant_terms
from analysis.causal.effects import (
    HORIZON_MONTHS,
    estimator_spread,
    horizon_status,
    ipw_risks,
)


def make_cohort(rows: list[tuple[float, int, int]]) -> pd.DataFrame:
    """``(os_months, os_event, treatment)`` rows into the frame the stack wants."""
    return pd.DataFrame({
        "os_months": [row[0] for row in rows],
        "os_event": [row[1] for row in rows],
        "chemo": [row[2] for row in rows],
        "patient_id": [f"P{index}" for index in range(len(rows))],
    })


class HorizonStatusTests(unittest.TestCase):
    def test_death_before_the_horizon_is_an_observed_event(self) -> None:
        outcome, observed = horizon_status(make_cohort([(30.0, 1, 1)]))
        self.assertEqual(outcome[0], 1.0)
        self.assertEqual(observed[0], 1.0)

    def test_followup_past_the_horizon_is_an_observed_survivor(self) -> None:
        outcome, observed = horizon_status(make_cohort([(90.0, 0, 1)]))
        self.assertEqual(outcome[0], 0.0)
        self.assertEqual(observed[0], 1.0)

    def test_censored_before_the_horizon_is_not_observed(self) -> None:
        """The case the censoring weights exist for: status is unknown."""
        outcome, observed = horizon_status(make_cohort([(30.0, 0, 1)]))
        self.assertEqual(outcome[0], 0.0)
        self.assertEqual(observed[0], 0.0)

    def test_death_exactly_at_the_horizon_counts_as_an_event(self) -> None:
        outcome, observed = horizon_status(make_cohort([(HORIZON_MONTHS, 1, 1)]))
        self.assertEqual(outcome[0], 1.0)
        self.assertEqual(observed[0], 1.0)

    def test_a_shorter_horizon_moves_a_death_out_of_the_window(self) -> None:
        cohort = make_cohort([(30.0, 1, 1)])
        outcome, observed = horizon_status(cohort, horizon=24.0)
        self.assertEqual(outcome[0], 0.0)
        self.assertEqual(observed[0], 1.0)


class IpwRisksTests(unittest.TestCase):
    def test_unit_weights_reproduce_the_plain_kaplan_meier_risk(self) -> None:
        cohort = make_cohort([
            (30.0, 1, 1), (90.0, 0, 1), (90.0, 0, 1), (90.0, 0, 1),
            (30.0, 1, 0), (30.0, 1, 0), (90.0, 0, 0), (90.0, 0, 0),
        ])
        result = ipw_risks(cohort, "chemo", np.ones(len(cohort)))
        self.assertAlmostEqual(result["risk_treated"], 0.25, places=6)
        self.assertAlmostEqual(result["risk_control"], 0.50, places=6)
        self.assertAlmostEqual(result["risk_difference"], -0.25, places=6)
        self.assertAlmostEqual(result["risk_ratio"], 0.5, places=6)

    def test_weights_shift_the_risk_toward_the_upweighted_patients(self) -> None:
        cohort = make_cohort([
            (30.0, 1, 1), (90.0, 0, 1),
            (90.0, 0, 0), (90.0, 0, 0),
        ])
        heavy = np.array([3.0, 1.0, 1.0, 1.0])
        result = ipw_risks(cohort, "chemo", heavy)
        self.assertGreater(result["risk_treated"], 0.5)

    def test_a_control_arm_with_no_risk_yields_nan_ratio(self) -> None:
        cohort = make_cohort([(30.0, 1, 1), (90.0, 0, 0)])
        result = ipw_risks(cohort, "chemo", np.ones(len(cohort)))
        self.assertTrue(np.isnan(result["risk_ratio"]))


class EstimatorSpreadTests(unittest.TestCase):
    def test_is_the_range_of_the_risk_differences(self) -> None:
        estimators = {
            "ipw_km": {"risk_difference": -0.036},
            "g_computation": {"risk_difference": -0.031},
            "aipw": {"risk_difference": -0.032},
        }
        self.assertAlmostEqual(estimator_spread(estimators), 0.005, places=9)

    def test_is_zero_when_the_estimators_agree_exactly(self) -> None:
        estimators = {name: {"risk_difference": -0.02}
                      for name in ("ipw_km", "g_computation", "aipw")}
        self.assertEqual(estimator_spread(estimators), 0.0)


class DropConstantTermsTests(unittest.TestCase):
    def frame(self, **overrides) -> pd.DataFrame:
        base = {
            "age": [50.0, 60.0, 70.0, 55.0],
            "tumor_size_mm": [10.0, 20.0, 30.0, 15.0],
            "lymph_pos": [0.0, 1.0, 2.0, 0.0],
            "stage": [1.0, 2.0, 3.0, 2.0],
            "grade": [1.0, 2.0, 3.0, 2.0],
            "er": [1.0, 1.0, 1.0, 1.0],
            "pr": [1.0, 0.0, 1.0, 0.0],
            "her2": [0.0, 1.0, 0.0, 1.0],
            "menopause": ["Pre", "Post", "Pre", "Post"],
            "subtype": ["HR+/HER2-", "HR+/HER2+", "HR+/HER2-", "HR+/HER2+"],
            "surgery": ["MASTECTOMY", "BREAST CONSERVING",
                        "MASTECTOMY", "BREAST CONSERVING"],
        }
        base.update(overrides)
        return pd.DataFrame(base)

    def test_drops_a_binary_covariate_that_is_constant(self) -> None:
        spec = drop_constant_terms(self.frame(), DEFAULT_SPEC)
        self.assertNotIn("er", spec.binary)
        self.assertIn("pr", spec.binary)

    def test_drops_categorical_levels_that_are_absent(self) -> None:
        spec = drop_constant_terms(self.frame(), DEFAULT_SPEC)
        levels = dict(spec.categorical)["subtype"]
        self.assertEqual(levels, ("HR+/HER2+",))

    def test_drops_a_categorical_level_present_in_every_row(self) -> None:
        """A level everyone has is as uninformative as one nobody has."""
        frame = self.frame(subtype=["TNBC"] * 4)
        spec = drop_constant_terms(frame, DEFAULT_SPEC)
        self.assertNotIn("subtype", dict(spec.categorical))

    def test_drops_a_constant_continuous_covariate(self) -> None:
        frame = self.frame(grade=[2.0, 2.0, 2.0, 2.0])
        spec = drop_constant_terms(frame, DEFAULT_SPEC)
        self.assertNotIn("grade", spec.continuous)

    def test_keeps_everything_that_varies(self) -> None:
        frame = self.frame(er=[1.0, 0.0, 1.0, 0.0])
        spec = drop_constant_terms(frame, DEFAULT_SPEC)
        self.assertEqual(spec.continuous, DEFAULT_SPEC.continuous)
        self.assertEqual(spec.binary, DEFAULT_SPEC.binary)

    def test_preserves_the_label(self) -> None:
        spec = drop_constant_terms(self.frame(), DEFAULT_SPEC.with_surgery())
        self.assertEqual(spec.label, "baseline + surgery")

    def test_a_spec_naming_a_missing_column_raises(self) -> None:
        """Silently dropping a confounder the cohort lacks would hide a bug."""
        frame = self.frame().drop(columns=["surgery"])
        with self.assertRaises(KeyError):
            drop_constant_terms(frame, DEFAULT_SPEC.with_surgery())


if __name__ == "__main__":
    unittest.main()
