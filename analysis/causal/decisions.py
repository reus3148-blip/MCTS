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
    max_iterations: int = 25,
) -> dict:
    """Trim to the overlap region and refit until the scores stay inside it.

    Refitting after trimming matters: a model fitted on the full cohort is
    dominated by the near-deterministic tails, and its scores are not the ones
    that balance the overlap population.

    But the refit produces *new* scores, and those can fall back outside the
    window - so one trim-then-refit pass does not generally leave a population
    that satisfies its own overlap condition. v0.6 and v0.7 did exactly one pass;
    for the endocrine decision restricted to ER-positive patients that left a
    population whose refitted propensities reached 0.999, worst |SMD| 1.58 and an
    effective sample of 77 out of 817. Iterating to a fixed point - trim, refit,
    trim again until every score is inside the window - gives 564 patients with
    worst |SMD| 0.032 and an effective sample of 500.

    ``max_iterations=1`` reproduces the original single-pass behaviour, which is
    how the pre-v1.3 numbers can still be regenerated.
    """
    current = cohort.reset_index(drop=True)
    iterations = 0
    converged = False
    for iterations in range(1, max_iterations + 1):
        propensity, weights = propensity_and_weights(current, treatment, spec)
        inside = (propensity >= bounds[0]) & (propensity <= bounds[1])
        if inside.all():
            converged = True
            break
        trimmed = current[inside].reset_index(drop=True)
        if trimmed[treatment].nunique() < 2:
            raise ValueError(f"trim {bounds} leaves a single arm")
        if iterations == max_iterations:
            current = trimmed
            propensity, weights = propensity_and_weights(current, treatment, spec)
            break
        current = trimmed

    if current[treatment].nunique() < 2:
        raise ValueError(f"trim {bounds} leaves a single arm")
    balance = balance_table(current, treatment, weights, spec)
    return {
        "bounds": list(bounds),
        "cohort": current,
        "propensity": propensity,
        "weights": weights,
        "balance": balance,
        "n": int(len(current)),
        "n_treated": int(current[treatment].sum()),
        "retained_pct": float(len(current) / len(cohort) * 100),
        "worst_abs_smd": float(balance["smd_weighted"].abs().max()),
        "balanced_pct": float(balance["balanced_after"].mean() * 100),
        "effective_sample_size": effective_sample_size(weights),
        "max_weight": float(weights.max()),
        "iterations": int(iterations),
        "converged": bool(converged),
        "propensity_inside_bounds_pct": float(
            ((propensity >= bounds[0]) & (propensity <= bounds[1])).mean() * 100),
    }


def drop_constant_terms(
    cohort: pd.DataFrame,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> CovariateSpec:
    """The same spec with covariates that do not vary in this cohort removed.

    Subgroup analyses make some covariates constant by construction - ``er`` is
    always 1 inside the ER-positive subgroup, and two of the four subtype levels
    are empty there. A constant cannot confound within the subgroup, but leaving
    it in makes the design matrix collinear with the intercept and turns its
    standardised mean difference into 0/0. Dropping it is the honest fix, and
    doing it here means every model in the analysis drops the same terms.
    """
    continuous = tuple(
        column for column in spec.continuous
        if cohort[column].astype(float).nunique(dropna=True) > 1
    )
    binary = tuple(
        column for column in spec.binary
        if cohort[column].astype(float).nunique(dropna=True) > 1
    )
    categorical = tuple(
        (column, tuple(level for level in levels
                       if 0 < (cohort[column] == level).sum() < len(cohort)))
        for column, levels in spec.categorical
    )
    categorical = tuple((column, levels) for column, levels in categorical if levels)
    return CovariateSpec(continuous, binary, categorical, spec.label)
