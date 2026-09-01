"""Unit tests for the IPW building blocks of the target trial emulation.

These are the first causal estimates the project will publish, so each step is
checked against a value that can be derived by hand or against a data-generating
process whose truth is known.
"""

from __future__ import annotations

import math
import unittest

import numpy as np

from analysis.causal.ipw import (
    e_value,
    e_value_for_interval,
    effective_sample_size,
    fit_logistic,
    sigmoid,
    stabilized_weights,
    standardized_mean_difference,
    truncate_weights,
    weighted_kaplan_meier,
)


class SigmoidTests(unittest.TestCase):
    def test_matches_the_definition(self) -> None:
        self.assertAlmostEqual(float(sigmoid(np.array([0.0]))[0]), 0.5)
        self.assertAlmostEqual(float(sigmoid(np.array([2.0]))[0]), 1 / (1 + math.exp(-2)))

    def test_is_stable_at_extremes(self) -> None:
        values = sigmoid(np.array([-800.0, 800.0]))
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertAlmostEqual(float(values[0]), 0.0)
        self.assertAlmostEqual(float(values[1]), 1.0)


class LogisticFitTests(unittest.TestCase):
    def test_recovers_known_coefficients(self) -> None:
        rng = np.random.default_rng(11)
        n = 6000
        x = rng.normal(size=n)
        treated = (rng.random(n) < sigmoid(0.5 + 1.2 * x)).astype(float)
        beta = fit_logistic(np.column_stack([np.ones(n), x]), treated)
        self.assertAlmostEqual(beta[0], 0.5, delta=0.08)
        self.assertAlmostEqual(beta[1], 1.2, delta=0.08)

    def test_stays_finite_under_separation(self) -> None:
        """Perfectly separated data must not blow the fit up."""
        x = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
        y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        beta = fit_logistic(np.column_stack([np.ones(6), x]), y)
        self.assertTrue(np.all(np.isfinite(beta)))

    def test_rejects_a_non_binary_outcome(self) -> None:
        with self.assertRaises(ValueError):
            fit_logistic(np.ones((3, 1)), np.array([0.0, 1.0, 2.0]))


class BalanceTests(unittest.TestCase):
    def test_crude_smd_matches_the_formula(self) -> None:
        values = np.array([1.0, 2.0, 3.0, 7.0, 8.0, 9.0])
        treated = np.array([0, 0, 0, 1, 1, 1])
        pooled = math.sqrt((values[3:].var(ddof=1) + values[:3].var(ddof=1)) / 2)
        self.assertAlmostEqual(
            standardized_mean_difference(values, treated), 6.0 / pooled)

    def test_weighting_can_shrink_an_imbalance(self) -> None:
        # Each arm needs within-arm spread, or the pooled sd is zero and the
        # standardised difference is undefined.
        values = np.array([0.0, 1.0, 1.0, 2.0])
        treated = np.array([0, 0, 1, 1])
        crude = standardized_mean_difference(values, treated)
        self.assertGreater(abs(crude), 1.0)
        # Upweight the arms' overlapping values, pulling the means together.
        weighted = standardized_mean_difference(
            values, treated, np.array([1.0, 3.0, 3.0, 1.0]))
        self.assertLess(abs(weighted), abs(crude))
        self.assertAlmostEqual(abs(weighted), abs(crude) / 2.0)

    def test_constant_covariate_is_perfectly_balanced(self) -> None:
        values = np.full(6, 2.5)
        treated = np.array([0, 0, 0, 1, 1, 1])
        self.assertEqual(standardized_mean_difference(values, treated), 0.0)


class WeightTests(unittest.TestCase):
    def test_stabilized_weights_centre_near_one(self) -> None:
        treated = np.array([1.0, 1.0, 0.0, 0.0])
        propensity = np.array([0.5, 0.5, 0.5, 0.5])
        weights = stabilized_weights(treated, propensity)
        self.assertTrue(np.allclose(weights, 1.0))

    def test_rare_treatment_gets_a_large_weight(self) -> None:
        treated = np.array([1.0, 0.0, 0.0, 0.0])
        propensity = np.array([0.05, 0.5, 0.5, 0.5])
        weights = stabilized_weights(treated, propensity)
        self.assertGreater(weights[0], weights[1])

    def test_rejects_a_degenerate_propensity(self) -> None:
        with self.assertRaises(ValueError):
            stabilized_weights(np.array([1.0, 0.0]), np.array([1.0, 0.5]))

    def test_truncation_clips_both_tails(self) -> None:
        weights = np.array([0.1, 1.0, 1.0, 1.0, 50.0])
        clipped = truncate_weights(weights, 0.2, 0.8)
        self.assertGreater(clipped.min(), 0.1)
        self.assertLess(clipped.max(), 50.0)

    def test_effective_sample_size_penalises_unequal_weights(self) -> None:
        equal = np.ones(100)
        unequal = np.concatenate([np.ones(99), np.array([100.0])])
        self.assertAlmostEqual(effective_sample_size(equal), 100.0)
        self.assertLess(effective_sample_size(unequal), 100.0)


class WeightedKaplanMeierTests(unittest.TestCase):
    def test_unweighted_matches_a_hand_computation(self) -> None:
        durations = np.array([1.0, 2.0, 2.0, 3.0, 4.0])
        events = np.array([1.0, 1.0, 0.0, 1.0, 0.0])
        curve = weighted_kaplan_meier(durations, events)
        # 4/5, then 3/4 of that, then 1/2 of that.
        self.assertTrue(np.allclose(curve.survival, [0.8, 0.6, 0.3]))

    def test_doubling_every_weight_changes_nothing(self) -> None:
        durations = np.array([1.0, 2.0, 3.0])
        events = np.array([1.0, 0.0, 1.0])
        plain = weighted_kaplan_meier(durations, events)
        doubled = weighted_kaplan_meier(durations, events, np.full(3, 2.0))
        self.assertTrue(np.allclose(plain.survival, doubled.survival))

    def test_weights_reproduce_duplicated_rows(self) -> None:
        durations = np.array([1.0, 2.0, 3.0])
        events = np.array([1.0, 1.0, 0.0])
        weighted = weighted_kaplan_meier(durations, events, np.array([2.0, 1.0, 1.0]))
        duplicated = weighted_kaplan_meier(
            np.array([1.0, 1.0, 2.0, 3.0]), np.array([1.0, 1.0, 1.0, 0.0]))
        self.assertTrue(np.allclose(weighted.survival, duplicated.survival))

    def test_survival_before_the_first_event_is_one(self) -> None:
        curve = weighted_kaplan_meier(np.array([5.0]), np.array([1.0]))
        self.assertEqual(curve.at(1.0), 1.0)
        self.assertAlmostEqual(curve.at(5.0), 0.0)

    def test_rejects_empty_input(self) -> None:
        with self.assertRaises(ValueError):
            weighted_kaplan_meier(np.array([]), np.array([]))


class EValueTests(unittest.TestCase):
    def test_null_effect_needs_no_confounding(self) -> None:
        self.assertAlmostEqual(e_value(1.0), 1.0)

    def test_is_symmetric_under_inversion(self) -> None:
        self.assertAlmostEqual(e_value(0.7), e_value(1 / 0.7))

    def test_matches_the_published_formula(self) -> None:
        # VanderWeele & Ding: RR + sqrt(RR * (RR - 1)).
        self.assertAlmostEqual(e_value(2.0), 2 + math.sqrt(2.0), places=6)

    def test_interval_covering_the_null_reports_one(self) -> None:
        self.assertEqual(e_value_for_interval(0.55, 1.25), 1.0)

    def test_interval_uses_the_limit_closest_to_the_null(self) -> None:
        self.assertAlmostEqual(e_value_for_interval(0.60, 0.85), e_value(0.85))
        self.assertAlmostEqual(e_value_for_interval(1.15, 2.00), e_value(1.15))

    def test_rejects_a_reversed_interval(self) -> None:
        with self.assertRaises(ValueError):
            e_value_for_interval(1.5, 0.5)


if __name__ == "__main__":
    unittest.main()
