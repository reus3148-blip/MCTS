"""Unit tests for the shared cohort helpers used by analysis 10/11/12.

These run without METABRIC: the sampling and manifest helpers only need a small
synthetic frame, which keeps the suite runnable on a machine that has no data/.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.dynamic.cohort import (  # noqa: E402
    BASE_SEED,
    REQUIRED_PATIENT_COLUMNS,
    balanced_subtype_sample,
    input_manifest,
    sha256,
)

SUBTYPES = ("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC")


def make_frame(per_subtype: int = 10) -> pd.DataFrame:
    rows = []
    for subtype_index, subtype in enumerate(SUBTYPES):
        for index in range(per_subtype):
            rows.append({
                "patient_id": f"P{subtype_index}{index:03d}",
                "subtype": subtype,
                "tumor_size_mm": 20.0 + index,
                "lymph_pos": index % 3,
                "stage": 2,
                "grade": 2,
                "er": 1, "pr": 1, "her2": 0,
            })
    return pd.DataFrame(rows)


class BalancedSubtypeSampleTests(unittest.TestCase):
    def test_draws_the_requested_count_from_each_subtype(self) -> None:
        sample = balanced_subtype_sample(make_frame(), per_subtype=3)
        counts = sample["subtype"].value_counts().to_dict()
        self.assertEqual(len(sample), 12)
        self.assertEqual(set(counts), set(SUBTYPES))
        self.assertTrue(all(count == 3 for count in counts.values()))

    def test_is_deterministic_for_a_given_seed(self) -> None:
        frame = make_frame()
        first = balanced_subtype_sample(frame, per_subtype=3)
        second = balanced_subtype_sample(frame, per_subtype=3)
        self.assertEqual(list(first["patient_id"]), list(second["patient_id"]))

    def test_a_larger_request_extends_rather_than_replaces(self) -> None:
        """Raising per_subtype must keep the smaller sample's patients."""
        frame = make_frame()
        small = set(balanced_subtype_sample(frame, per_subtype=2)["patient_id"])
        large = set(balanced_subtype_sample(frame, per_subtype=3)["patient_id"])
        self.assertTrue(small.issubset(large))

    def test_drops_patients_with_missing_required_values(self) -> None:
        frame = make_frame()
        frame.loc[frame.index[:8], "tumor_size_mm"] = None
        sample = balanced_subtype_sample(frame, per_subtype=3)
        self.assertFalse(sample[list(REQUIRED_PATIENT_COLUMNS)].isna().any().any())

    def test_takes_what_exists_when_a_subtype_is_short(self) -> None:
        frame = make_frame(per_subtype=1)
        sample = balanced_subtype_sample(frame, per_subtype=5)
        self.assertEqual(len(sample), len(SUBTYPES))

    def test_rejects_a_non_positive_request(self) -> None:
        with self.assertRaises(ValueError):
            balanced_subtype_sample(make_frame(), per_subtype=0)

    def test_rejects_a_frame_with_no_complete_patients(self) -> None:
        frame = make_frame()
        frame["grade"] = None
        with self.assertRaises(ValueError):
            balanced_subtype_sample(frame, per_subtype=1)

    def test_different_seeds_can_draw_different_patients(self) -> None:
        frame = make_frame()
        default = balanced_subtype_sample(frame, per_subtype=2)
        shifted = balanced_subtype_sample(frame, per_subtype=2, seed=BASE_SEED + 5)
        self.assertNotEqual(list(default["patient_id"]), list(shifted["patient_id"]))


class ManifestHelperTests(unittest.TestCase):
    def test_sha256_matches_hashlib(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            payload = b"patient_id,subtype\nP1,TNBC\n"
            path.write_bytes(payload)
            self.assertEqual(sha256(path), hashlib.sha256(payload).hexdigest())

    def test_input_manifest_reports_repo_relative_posix_paths(self) -> None:
        manifest = input_manifest({
            "assumptions": ROOT / "configs" / "dynamic_poc_v0_2.json",
        })
        entry = manifest["assumptions"]
        self.assertEqual(entry["path"], "configs/dynamic_poc_v0_2.json")
        self.assertEqual(len(entry["sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
