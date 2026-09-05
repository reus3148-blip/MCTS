from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from analysis.causal.decisions import CovariateSpec, trim_to_overlap

BOUNDS = (0.10, 0.90)
#: Two continuous covariates only; the rest of the default spec is irrelevant here.
SPEC = CovariateSpec(continuous=("age", "tumor_size_mm"), binary=(), categorical=())


def separated_cohort(strength: float, seed: int = 11, n: int = 600) -> pd.DataFrame:
    """Treatment strongly predicted by a heavy-tailed covariate.

    The heavy tail is what makes trimming change the answer: ``design_matrix``
    standardises continuous covariates by the *cohort's* own mean and spread, so
    removing the tails rescales the design and the refitted model is sharper than
    the one that chose the trim.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    x[: n // 10] *= 6.0
    z = rng.normal(size=n)
    propensity = 1.0 / (1.0 + np.exp(-(strength * x + 1.5 * z)))
    return pd.DataFrame({
        "age": x,
        "tumor_size_mm": z,
        "chemo": (rng.random(n) < propensity).astype(int),
    })


class TrimToOverlapTests(unittest.TestCase):
    def test_iterating_reaches_a_population_that_satisfies_its_own_bounds(self) -> None:
        result = trim_to_overlap(separated_cohort(6.0), "chemo", BOUNDS, SPEC)
        self.assertTrue(result["converged"])
        self.assertEqual(result["propensity_inside_bounds_pct"], 100.0)

    def test_a_single_pass_can_leave_scores_outside_the_window(self) -> None:
        """The v0.6/v0.7 defect: one trim-then-refit is not a fixed point."""
        once = trim_to_overlap(
            separated_cohort(6.0), "chemo", BOUNDS, SPEC, max_iterations=1)
        self.assertFalse(once["converged"])
        self.assertLess(once["propensity_inside_bounds_pct"], 100.0)

    def test_iterating_removes_more_patients_than_one_pass(self) -> None:
        cohort = separated_cohort(6.0)
        once = trim_to_overlap(cohort, "chemo", BOUNDS, SPEC, max_iterations=1)
        iterated = trim_to_overlap(cohort, "chemo", BOUNDS, SPEC)
        self.assertLess(iterated["n"], once["n"])
        self.assertGreater(iterated["iterations"], 1)

    def test_an_already_overlapping_cohort_converges_without_trimming(self) -> None:
        cohort = separated_cohort(0.2)
        result = trim_to_overlap(cohort, "chemo", (0.0, 1.0), SPEC)
        self.assertTrue(result["converged"])
        self.assertEqual(result["n"], len(cohort))
        self.assertEqual(result["iterations"], 1)

    def test_reports_the_iteration_count_and_retention(self) -> None:
        result = trim_to_overlap(separated_cohort(6.0), "chemo", BOUNDS, SPEC)
        self.assertGreaterEqual(result["iterations"], 1)
        self.assertLess(result["retained_pct"], 100.0)
        self.assertGreater(result["retained_pct"], 0.0)

    def test_weights_and_balance_describe_the_final_population(self) -> None:
        result = trim_to_overlap(separated_cohort(6.0), "chemo", BOUNDS, SPEC)
        self.assertEqual(len(result["weights"]), result["n"])
        self.assertEqual(len(result["propensity"]), result["n"])
        self.assertEqual(len(result["cohort"]), result["n"])

    def test_a_trim_that_empties_one_arm_raises(self) -> None:
        cohort = separated_cohort(6.0)
        with self.assertRaises(ValueError):
            trim_to_overlap(cohort, "chemo", (0.499, 0.501), SPEC)


if __name__ == "__main__":
    unittest.main()
