"""Run the first simplified NCCN-policy concordance analysis.

The reusable policy rules live in ``analysis/nccn_policy.py``. This numbered
file remains the command-line entry point so the original analysis sequence is
preserved.
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.nccn_policy import (  # noqa: E402
    apply_nccn,
    concordance_by_subtype,
    concordance_per_decision,
    normalize_actual_surgery,
)

CSV_IN = ROOT / "data" / "processed" / "patients.csv"
CSV_OUT = ROOT / "data" / "processed" / "patients_with_nccn.csv"


def main() -> None:
    df = pd.read_csv(CSV_IN)
    print(f"loaded {len(df)} patients")

    recommendations = df.apply(apply_nccn, axis=1)
    output = pd.concat([df, recommendations], axis=1)
    output.to_csv(CSV_OUT, index=False, encoding="utf-8")
    print(f"saved -> {CSV_OUT}")

    valid = output.dropna(subset=["subtype", "os_event"]).copy()
    print(f"\nanalysis subset: {len(valid)}")

    print("\n=== decision-level concordance ===")
    print(concordance_per_decision(valid).to_string(index=False))

    print("\n=== subtype x decision concordance (%) ===")
    print(concordance_by_subtype(valid).to_string(index=False))

    valid["actual_surgery"] = valid["surgery"].apply(
        normalize_actual_surgery
    )
    all_compared = (
        valid["rec_surgery"].notna()
        & valid["actual_surgery"].notna()
        & valid["rec_chemo"].notna()
        & valid["chemo"].notna()
        & valid["rec_hormone"].notna()
        & valid["hormone"].notna()
        & valid["rec_radio"].notna()
        & valid["radio"].notna()
    )
    all_matched = all_compared & (
        (valid["rec_surgery"] == valid["actual_surgery"])
        & (valid["rec_chemo"] == valid["chemo"])
        & (valid["rec_hormone"] == valid["hormone"])
        & (valid["rec_radio"] == valid["radio"])
    )
    percentage = all_matched.sum() / max(all_compared.sum(), 1) * 100
    print(f"\nall four comparable: {all_compared.sum()}")
    print(f"all four matched: {all_matched.sum()} ({percentage:.1f}%)")


if __name__ == "__main__":
    main()
