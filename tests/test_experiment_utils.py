from __future__ import annotations

import math
import unittest

from analysis.dynamic.experiment_utils import (
    confidence_interval,
    rescale_response_major,
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


if __name__ == "__main__":
    unittest.main()
