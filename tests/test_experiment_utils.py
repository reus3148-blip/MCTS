from __future__ import annotations

import math
import unittest

from analysis.dynamic.experiment_utils import (
    confidence_interval,
    rescale_response_major,
    summarize_decision_replicates,
)


class ConfidenceIntervalTests(unittest.TestCase):
    def test_mean_and_symmetric_interval(self) -> None:
        result = confidence_interval([1.0, 2.0, 3.0, 4.0, 5.0])
        self.assertAlmostEqual(result["mean"], 3.0)
        self.assertAlmostEqual(result["sd"], math.sqrt(2.5))
        # interval is symmetric around the mean
        self.assertAlmostEqual(
            result["mean"] - result["ci95_low"],
            result["ci95_high"] - result["mean"],
        )
        self.assertEqual(result["n_seeds"], 5)

    def test_single_value_has_zero_width(self) -> None:
        result = confidence_interval([0.7])
        self.assertEqual(result["mean"], 0.7)
        self.assertEqual(result["sd"], 0.0)
        self.assertEqual(result["ci95_low"], 0.7)
        self.assertEqual(result["ci95_high"], 0.7)

    def test_empty_raises(self) -> None:
        with self.assertRaises(ValueError):
            confidence_interval([])


class RescaleResponseMajorTests(unittest.TestCase):
    def test_probabilities_still_sum_to_one(self) -> None:
        block = {"major": 0.35, "partial": 0.45, "none": 0.20}
        rescaled = rescale_response_major(block, 0.50)
        self.assertAlmostEqual(rescaled["major"], 0.50)
        self.assertAlmostEqual(sum(rescaled.values()), 1.0)

    def test_partial_none_ratio_is_preserved(self) -> None:
        block = {"major": 0.25, "partial": 0.50, "none": 0.25}  # partial:none = 2:1
        rescaled = rescale_response_major(block, 0.40)
        self.assertAlmostEqual(rescaled["partial"] / rescaled["none"], 2.0)

    def test_rejects_out_of_range_major(self) -> None:
        block = {"major": 0.3, "partial": 0.4, "none": 0.3}
        with self.assertRaises(ValueError):
            rescale_response_major(block, 1.5)

    def test_rejects_degenerate_remainder(self) -> None:
        block = {"major": 1.0, "partial": 0.0, "none": 0.0}
        with self.assertRaises(ValueError):
            rescale_response_major(block, 0.5)


class SummarizeDecisionReplicatesTests(unittest.TestCase):
    """The v0.4 budget study reads equipoise off these three numbers."""

    def test_unanimous_replicates_report_full_agreement(self) -> None:
        records = [("a", {"a": 0.6, "b": 0.4})] * 5
        result = summarize_decision_replicates(records)
        self.assertEqual(result["modal_action"], "a")
        self.assertEqual(result["agreement_pct"], 100.0)
        self.assertEqual(result["distinct_actions_chosen"], 1)
        self.assertAlmostEqual(result["value_gap"], 0.2)
        self.assertEqual(result["value_noise_sd"], 0.0)
        self.assertTrue(math.isnan(result["separation"]))

    def test_split_replicates_report_partial_agreement(self) -> None:
        records = [
            ("a", {"a": 0.50, "b": 0.48}),
            ("a", {"a": 0.52, "b": 0.49}),
            ("b", {"a": 0.47, "b": 0.51}),
            ("b", {"a": 0.48, "b": 0.52}),
        ]
        result = summarize_decision_replicates(records)
        self.assertEqual(result["agreement_pct"], 50.0)
        self.assertEqual(result["distinct_actions_chosen"], 2)
        self.assertEqual(result["n_legal_actions"], 2)

    def test_equipoise_shows_a_gap_below_the_noise(self) -> None:
        """Near-equal actions: separation < 1 is the equipoise signature."""
        records = [
            ("a", {"a": 0.500, "b": 0.499}),
            ("b", {"a": 0.470, "b": 0.520}),
            ("a", {"a": 0.530, "b": 0.480}),
            ("b", {"a": 0.480, "b": 0.510}),
        ]
        result = summarize_decision_replicates(records)
        self.assertLess(abs(result["value_gap"]), result["value_noise_sd"])
        self.assertLess(abs(result["separation"]), 1.0)

    def test_a_clear_winner_separates_well_beyond_the_noise(self) -> None:
        records = [
            ("a", {"a": 0.80, "b": 0.30}),
            ("a", {"a": 0.81, "b": 0.31}),
            ("a", {"a": 0.79, "b": 0.29}),
        ]
        result = summarize_decision_replicates(records)
        self.assertGreater(result["separation"], 10.0)

    def test_single_action_node_has_no_gap(self) -> None:
        records = [("only", {"only": 0.5})] * 3
        result = summarize_decision_replicates(records)
        self.assertTrue(math.isnan(result["value_gap"]))
        self.assertEqual(result["n_legal_actions"], 1)

    def test_ties_break_deterministically_by_action_name(self) -> None:
        forward = summarize_decision_replicates(
            [("b", {"a": 0.5, "b": 0.5}), ("a", {"a": 0.5, "b": 0.5})]
        )
        reverse = summarize_decision_replicates(
            [("a", {"a": 0.5, "b": 0.5}), ("b", {"a": 0.5, "b": 0.5})]
        )
        self.assertEqual(forward["modal_action"], reverse["modal_action"])
        self.assertEqual(forward["modal_action"], "a")

    def test_empty_raises(self) -> None:
        with self.assertRaises(ValueError):
            summarize_decision_replicates([])


if __name__ == "__main__":
    unittest.main()
