"""Treatment-decision cohorts and propensity machinery shared by the emulations.

``analysis/15`` estimates the effect of one decision; ``analysis/17`` adds doubly
robust estimation and asks the same positivity question of every decision. Both
need identical eligibility, covariate coding and trimming, or their numbers would
not be comparable - so it lives here rather than in either script.

Nothing here fixes *which* decision is being studied: the treatment column is a
parameter, which is what lets the overlap map cover chemotherapy, endocrine
therapy and radiotherapy on the same footing.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .ipw import (
    effective_sample_size,
    fit_logistic,
    sigmoid,
    stabilized_weights,
    standardized_mean_difference,
)

#: Baseline confounders. These are the variables METABRIC records at diagnosis
#: that a clinician would plausibly weigh when choosing adjuvant treatment.
CONTINUOUS = ["age", "tumor_size_mm", "lymph_pos", "stage", "grade"]
BINARY = ["er", "pr", "her2"]
CATEGORICAL = {
    "menopause": ["Pre"],                            # Post is the reference level
    "subtype": ["HR+/HER2+", "HR-/HER2+", "TNBC"],   # HR+/HER2- is reference
}

BALANCE_THRESHOLD = 0.1        # |SMD| below this is conventionally "balanced"


@dataclass(frozen=True)
class CovariateSpec:
    """Which baseline variables a propensity model adjusts for.

    Made explicit because *which* confounders are included decides whether a
    decision looks answerable. Radiotherapy appears to have excellent overlap
    until surgery type is added, at which point it does not - so the spec has to
    be a reported choice, not a hidden constant.
    """

    continuous: tuple[str, ...] = field(default=tuple(CONTINUOUS))
    binary: tuple[str, ...] = field(default=tuple(BINARY))
    categorical: tuple[tuple[str, tuple[str, ...]], ...] = field(
        default=tuple((column, tuple(levels)) for column, levels in CATEGORICAL.items())
    )
    label: str = "baseline"

    @property
    def columns(self) -> list[str]:
        return list(self.continuous) + list(self.binary) + [
            column for column, _ in self.categorical
        ]

    @property
    def names(self) -> list[str]:
        return list(self.continuous) + list(self.binary) + [
            f"{column}={level}" for column, levels in self.categorical
            for level in levels
        ]

    def with_surgery(self) -> "CovariateSpec":
        """The same spec plus surgery type, the main driver of radiotherapy."""
        return CovariateSpec(
            self.continuous, self.binary,
            self.categorical + (("surgery", ("BREAST CONSERVING",)),),
            label="baseline + surgery",
        )


DEFAULT_SPEC = CovariateSpec()

#: Decisions the METABRIC treatment columns can express, in guideline order.
DECISIONS = {
    "chemo": "보조 항암치료",
    "hormone": "호르몬치료",
    "radio": "방사선치료",
}


def build_cohort(
    raw: pd.DataFrame,
    treatment: str,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> pd.DataFrame:
    """Protocol §1.1 eligibility for one decision, as far as METABRIC allows."""
    required = (
        spec.columns + [treatment, "os_months", "os_event", "patient_id"]
    )
    cohort = raw.dropna(subset=required).copy()
    cohort = cohort[cohort["os_months"] > 0]
    cohort[treatment] = cohort[treatment].astype(int)
    return cohort.reset_index(drop=True)


def covariate_frame(
    cohort: pd.DataFrame,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> pd.DataFrame:
    """Covariates on their natural scale, for balance tables."""
    frame = cohort[list(spec.continuous) + list(spec.binary)].astype(float).copy()
    for column, levels in spec.categorical:
        for level in levels:
            frame[f"{column}={level}"] = (cohort[column] == level).astype(float)
    return frame


def design_matrix(
    cohort: pd.DataFrame,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> np.ndarray:
    """Intercept, standardised continuous terms, and reference-coded factors."""
    columns = [np.ones(len(cohort))]
    for column in spec.continuous:
        values = cohort[column].to_numpy(dtype=float)
        spread = values.std(ddof=1)
        columns.append((values - values.mean()) / (spread if spread else 1.0))
    for column in spec.binary:
        columns.append(cohort[column].to_numpy(dtype=float))
    for column, levels in spec.categorical:
        for level in levels:
            columns.append((cohort[column] == level).to_numpy(dtype=float))
    return np.column_stack(columns)


def propensity_and_weights(
    cohort: pd.DataFrame,
    treatment: str,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> tuple[np.ndarray, np.ndarray]:
    treated = cohort[treatment].to_numpy(dtype=float)
    design = design_matrix(cohort, spec)
    propensity = np.clip(
        sigmoid(design @ fit_logistic(design, treated)), 1e-6, 1 - 1e-6)
    return propensity, stabilized_weights(treated, propensity)


def balance_table(
    cohort: pd.DataFrame,
    treatment: str,
    weights: np.ndarray,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> pd.DataFrame:
    treated = cohort[treatment].to_numpy(dtype=float)
    covariates = covariate_frame(cohort, spec)
    rows = []
    for column in covariates.columns:
        values = covariates[column].to_numpy(dtype=float)
        rows.append({
            "covariate": column,
            "mean_treated": float(values[treated == 1].mean()),
            "mean_control": float(values[treated == 0].mean()),
            "smd_crude": standardized_mean_difference(values, treated),
            "smd_weighted": standardized_mean_difference(values, treated, weights),
        })
    table = pd.DataFrame(rows)
    table["balanced_after"] = table["smd_weighted"].abs() < BALANCE_THRESHOLD
    return table


def trim_to_overlap(
    cohort: pd.DataFrame,
    treatment: str,
    bounds: tuple[float, float],
    spec: CovariateSpec = DEFAULT_SPEC,
) -> dict:
    """Fit, trim to the overlap region, refit inside it, then weight and diagnose.

    Refitting after trimming matters: a model fitted on the full cohort is
    dominated by the near-deterministic tails, and its scores are not the ones
    that balance the overlap population.
    """
    propensity, _ = propensity_and_weights(cohort, treatment, spec)
    keep = (propensity >= bounds[0]) & (propensity <= bounds[1])
    trimmed = cohort[keep].reset_index(drop=True)
    if trimmed[treatment].nunique() < 2:
        raise ValueError(f"trim {bounds} leaves a single arm")

    inner_propensity, weights = propensity_and_weights(trimmed, treatment, spec)
    balance = balance_table(trimmed, treatment, weights, spec)
    return {
        "bounds": list(bounds),
        "cohort": trimmed,
        "propensity": inner_propensity,
        "weights": weights,
        "balance": balance,
        "n": int(len(trimmed)),
        "n_treated": int(trimmed[treatment].sum()),
        "retained_pct": float(len(trimmed) / len(cohort) * 100),
        "worst_abs_smd": float(balance["smd_weighted"].abs().max()),
        "balanced_pct": float(balance["balanced_after"].mean() * 100),
        "effective_sample_size": effective_sample_size(weights),
        "max_weight": float(weights.max()),
    }
