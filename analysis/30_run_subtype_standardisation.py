"""Post-hoc: what the subtype-balanced sample does to the headline (v1.2).

NOT PRE-SPECIFIED. This was prompted by the v1.2 per-patient table, where the
utility gap turned out to differ several-fold across molecular subtypes. Every
study since v0.3 has sampled *equal numbers per subtype* - the right choice for
keeping rare subtypes present in a tiny cohort - but that also means the headline
average is taken over a cohort whose subtype mix is nothing like the population
it was drawn from. If the gap varies by subtype, the two averages differ.

This script recomputes the same 40 patients' mean gap under the subtype mix of
the held-out test split, and reports both. It reads the committed per-patient
table and re-runs no search, so it is cheap and cannot drift from v1.2.

Direct standardisation over four strata of ten patients is crude, and the result
is exploratory. It is reported because the direction and size of the difference
change how the headline should be read, not because this design can estimate a
population mean.
"""

from __future__ import annotations

from datetime import datetime
import json
import math
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.dynamic.cohort import (  # noqa: E402
    REQUIRED_PATIENT_COLUMNS,
    build_reward_models,
    git_commit,
    input_manifest,
)

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
REPORT_DIR = ROOT / "reports" / "cohort-replication-v1.2"
TABLE_DIR = REPORT_DIR / "tables"
PER_PATIENT = TABLE_DIR / "per_patient_gaps.csv"

RUN_DATE = "2026-09-04"


def main() -> None:
    per_patient = pd.read_csv(PER_PATIENT)
    raw = pd.read_csv(INPUT_CSV)
    _, _, os_test = build_reward_models(raw)
    eligible = os_test.dropna(subset=list(REQUIRED_PATIENT_COLUMNS))
    prevalence = eligible["subtype"].value_counts()
    prevalence = prevalence / prevalence.sum()

    rows = []
    for subtype, block in per_patient.groupby("subtype"):
        values = block["baseline_utility_gap"]
        rows.append({
            "subtype": subtype,
            "n_sampled": int(len(values)),
            "sample_weight": 1.0 / per_patient["subtype"].nunique(),
            "population_share": float(prevalence.get(subtype, 0.0)),
            "mean_gap": float(values.mean()),
            "sd_gap": float(values.std(ddof=1)),
            "standard_error": float(values.std(ddof=1) / math.sqrt(len(values))),
        })
    table = pd.DataFrame(rows).sort_values("population_share", ascending=False)

    balanced_mean = float((table["sample_weight"] * table["mean_gap"]).sum())
    balanced_se = math.sqrt(float(
        ((table["sample_weight"] * table["standard_error"]) ** 2).sum()))
    standardised_mean = float((table["population_share"] * table["mean_gap"]).sum())
    standardised_se = math.sqrt(float(
        ((table["population_share"] * table["standard_error"]) ** 2).sum()))

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "cohort-replication-v1.2-posthoc-subtype",
        "prespecified": False,
        "question": (
            "Every study since v0.3 averages over an equal number of patients per "
            "subtype. If the gap varies by subtype, what does that averaging do to "
            "the headline?"
        ),
        "estimand": (
            "Mean MCTS-minus-NCCN utility gap over the v1.2 forty patients, under "
            "(a) the equal-per-subtype weights the sample was drawn with and "
            "(b) the subtype mix of the held-out test split it was drawn from."
        ),
        "scope_warning": (
            "Direct standardisation over four strata of ten patients. Exploratory; "
            "not a population estimate, and still a synthetic simulator."
        ),
        "population_reference": {
            "source": "held-out OS test split, complete-record patients",
            "n_patients": int(len(eligible)),
            "subtype_share": {k: float(v) for k, v in prevalence.items()},
        },
        "balanced_sample_mean": balanced_mean,
        "balanced_sample_standard_error": balanced_se,
        "prevalence_standardised_mean": standardised_mean,
        "prevalence_standardised_standard_error": standardised_se,
        "ratio_standardised_to_balanced": standardised_mean / balanced_mean,
        "largest_subtype": {
            "subtype": str(table.iloc[0]["subtype"]),
            "population_share": float(table.iloc[0]["population_share"]),
            "sample_share": float(table.iloc[0]["sample_weight"]),
            "mean_gap": float(table.iloc[0]["mean_gap"]),
        },
        "subtype_gap_spread": {
            "min": float(table["mean_gap"].min()),
            "max": float(table["mean_gap"].max()),
            "ratio": float(table["mean_gap"].max() / table["mean_gap"].min()),
        },
    }

    table.to_csv(TABLE_DIR / "subtype_standardisation.csv", index=False)
    with (REPORT_DIR / "metrics_posthoc_subtype.json").open(
            "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({"data": INPUT_CSV, "per_patient": PER_PATIENT}),
        "entry_point": "analysis/30_run_subtype_standardisation.py",
    }
    with (REPORT_DIR / "run_manifest_posthoc_subtype.json").open(
            "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(table.round(4).to_string(index=False))
    print(f"\nbalanced        {balanced_mean:+.4f} +- {balanced_se:.4f}")
    print(f"standardised    {standardised_mean:+.4f} +- {standardised_se:.4f}")
    print(f"ratio           {metrics['ratio_standardised_to_balanced']:.2f}")
    print(f"\nsaved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
