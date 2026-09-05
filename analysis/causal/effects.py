"""Point-treatment effect estimation: IPCW, three estimators, bootstrap.

Extracted from ``analysis/17_run_doubly_robust.py`` unchanged in behaviour so a
second decision (v1.3) can reuse the same stack rather than copy it. The
numbered scripts cannot be imported - their module names start with a digit - so
this is also the only place these functions can be unit-tested from.

Two things became parameters during the move, both defaulting to what v0.7 used:

* ``spec`` - v0.7's Part A ran only the baseline covariate set, so the outcome
  and censoring models could hardcode it. A decision whose confounders include
  surgery type needs the *same* spec in all three models, not just in the
  propensity model, or the estimators quietly adjust for different things.
* ``horizon``/``edges``/``trim``/``replicates`` - so a decision with different
  follow-up or overlap can be run without editing the module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .decisions import (
    CovariateSpec,
    DEFAULT_SPEC,
    design_matrix,
    trim_to_overlap,
)
from .ipw import (
    aipw_risk,
    discrete_time_rows,
    fit_logistic,
    sigmoid,
    survival_probability,
    weighted_kaplan_meier,
)

#: v0.7's settings, kept as defaults so its numbers reproduce exactly.
HORIZON_MONTHS = 60.0
INTERVAL_EDGES = np.array([0.0, 12.0, 24.0, 36.0, 48.0, 60.0])
DEFAULT_TRIM = (0.10, 0.90)
DEFAULT_BOOTSTRAP = 500


def horizon_status(
    cohort: pd.DataFrame, horizon: float = HORIZON_MONTHS
) -> tuple[np.ndarray, np.ndarray]:
    """``(outcome, observed)`` at the horizon.

    A patient's status is known if they died within the horizon, or if follow-up
    reached it. Anyone censored earlier contributes no outcome, which is what the
    censoring weights exist to correct for.
    """
    months = cohort["os_months"].to_numpy(dtype=float)
    died = cohort["os_event"].to_numpy(dtype=float) == 1
    outcome = ((months <= horizon) & died).astype(float)
    observed = ((months >= horizon) | ((months < horizon) & died))
    return outcome, observed.astype(float)


def pooled_design(
    cohort: pd.DataFrame,
    treatment: str,
    rows: np.ndarray,
    intervals: np.ndarray,
    force_treatment: float | None = None,
    spec: CovariateSpec = DEFAULT_SPEC,
    edges: np.ndarray = INTERVAL_EDGES,
) -> np.ndarray:
    """Person-interval design: interval dummies, treatment, baseline covariates."""
    base = design_matrix(cohort, spec)[rows]
    n_intervals = edges.size - 1
    dummies = np.zeros((rows.size, n_intervals - 1))
    for index in range(1, n_intervals):     # interval 0 is the reference
        dummies[intervals == index, index - 1] = 1.0
    if force_treatment is None:
        exposure = cohort[treatment].to_numpy(dtype=float)[rows]
    else:
        exposure = np.full(rows.size, float(force_treatment))
    return np.column_stack([base, dummies, exposure])


def censoring_weights(
    cohort: pd.DataFrame,
    treatment: str,
    spec: CovariateSpec = DEFAULT_SPEC,
    horizon: float = HORIZON_MONTHS,
    edges: np.ndarray = INTERVAL_EDGES,
) -> np.ndarray:
    """Probability of reaching one's horizon status uncensored, given X and A.

    Estimated as a discrete-time censoring hazard (pooled logistic over yearly
    intervals) rather than assumed constant, so that censoring may depend on the
    same covariates that drive treatment.
    """
    months = cohort["os_months"].to_numpy(dtype=float)
    died = cohort["os_event"].to_numpy(dtype=float) == 1
    censored_early = ((months < horizon) & ~died).astype(float)
    rows, intervals, censor_here = discrete_time_rows(months, censored_early, edges)
    design = pooled_design(cohort, treatment, rows, intervals, spec=spec, edges=edges)
    hazard = sigmoid(design @ fit_logistic(design, censor_here))
    return np.clip(survival_probability(hazard, rows, len(cohort)), 1e-3, 1.0)


def outcome_predictions(
    cohort: pd.DataFrame,
    treatment: str,
    outcome: np.ndarray,
    observed: np.ndarray,
    uncensored: np.ndarray,
    spec: CovariateSpec = DEFAULT_SPEC,
) -> tuple[np.ndarray, np.ndarray]:
    """IPCW-weighted outcome model, evaluated with treatment set to 1 and to 0.

    Fitting only on patients whose status is known would over-represent those
    with long follow-up; weighting them by 1/P(uncensored) restores the
    population they stand in for.
    """
    design = np.column_stack([
        design_matrix(cohort, spec), cohort[treatment].to_numpy(dtype=float)])
    weights = observed / uncensored
    beta = fit_logistic(design, outcome, sample_weight=weights)

    def predict(value: float) -> np.ndarray:
        forced = np.column_stack([
            design_matrix(cohort, spec), np.full(len(cohort), value)])
        return sigmoid(forced @ beta)

    return predict(1.0), predict(0.0)


def ipw_risks(
    cohort: pd.DataFrame,
    treatment: str,
    weights: np.ndarray,
    horizon: float = HORIZON_MONTHS,
) -> dict:
    """v0.6's estimator: IPT-weighted Kaplan-Meier risk at the horizon."""
    def arm(value: int) -> float:
        mask = (cohort[treatment] == value).to_numpy()
        curve = weighted_kaplan_meier(
            cohort.loc[mask, "os_months"].to_numpy(dtype=float),
            cohort.loc[mask, "os_event"].to_numpy(dtype=float),
            weights[mask])
        return 1.0 - curve.at(horizon)

    treated_risk, control_risk = arm(1), arm(0)
    return {
        "risk_treated": treated_risk,
        "risk_control": control_risk,
        "risk_difference": treated_risk - control_risk,
        "risk_ratio": treated_risk / control_risk if control_risk > 0 else float("nan"),
    }


def all_estimators(
    cohort: pd.DataFrame,
    treatment: str,
    weights: np.ndarray,
    propensity: np.ndarray,
    spec: CovariateSpec = DEFAULT_SPEC,
    horizon: float = HORIZON_MONTHS,
    edges: np.ndarray = INTERVAL_EDGES,
) -> dict[str, dict]:
    """IPW, g-computation and AIPW over the same population.

    All three lean on different things - the treatment model, the outcome model,
    and either-one-of-them respectively - so their agreement (or not) is the
    check on model misspecification.
    """
    outcome, observed = horizon_status(cohort, horizon)
    uncensored = censoring_weights(cohort, treatment, spec, horizon, edges)
    predicted_treated, predicted_control = outcome_predictions(
        cohort, treatment, outcome, observed, uncensored, spec)

    g_treated = float(predicted_treated.mean())
    g_control = float(predicted_control.mean())
    return {
        "ipw_km": ipw_risks(cohort, treatment, weights, horizon),
        "g_computation": {
            "risk_treated": g_treated,
            "risk_control": g_control,
            "risk_difference": g_treated - g_control,
            "risk_ratio": g_treated / g_control if g_control > 0 else float("nan"),
        },
        "aipw": aipw_risk(
            cohort[treatment].to_numpy(dtype=float), outcome, observed,
            propensity, uncensored, predicted_treated, predicted_control),
        "censoring": {
            "min_uncensored_probability": float(uncensored.min()),
            "mean_uncensored_probability": float(uncensored.mean()),
            "observed_at_horizon_pct": float(observed.mean() * 100),
        },
    }


def bootstrap_aipw(
    cohort: pd.DataFrame,
    treatment: str,
    seed: int,
    spec: CovariateSpec = DEFAULT_SPEC,
    trim: tuple[float, float] = DEFAULT_TRIM,
    replicates: int = DEFAULT_BOOTSTRAP,
    horizon: float = HORIZON_MONTHS,
    edges: np.ndarray = INTERVAL_EDGES,
    max_iterations: int = 25,
) -> dict:
    """Percentile CI repeating fit, trim, refit and both nuisance models.

    Resampling *before* trimming is deliberate: the trim is part of the estimator,
    so its own variability has to be inside the interval - which also means
    ``max_iterations`` has to match the point estimate's, or the interval would
    describe a different estimator.
    """
    rng = np.random.default_rng(seed)
    n = len(cohort)
    differences, ratios = [], []
    for _ in range(replicates):
        sample = cohort.iloc[rng.integers(0, n, size=n)].reset_index(drop=True)
        try:
            trimmed = trim_to_overlap(sample, treatment, trim, spec,
                                      max_iterations=max_iterations)
            result = all_estimators(
                trimmed["cohort"], treatment, trimmed["weights"],
                trimmed["propensity"], spec, horizon, edges)["aipw"]
        except (np.linalg.LinAlgError, ValueError):
            continue
        if np.isfinite(result["risk_difference"]):
            differences.append(result["risk_difference"])
        if np.isfinite(result["risk_ratio"]) and result["risk_ratio"] > 0:
            ratios.append(result["risk_ratio"])
    return {
        "risk_difference": [float(np.quantile(differences, q)) for q in (0.025, 0.975)],
        "risk_ratio": [float(np.quantile(ratios, q)) for q in (0.025, 0.975)],
        "replicates": len(differences),
    }


def estimator_spread(estimators: dict[str, dict]) -> float:
    """Largest gap between any two risk differences.

    Reported next to the confidence interval: a spread that is small relative to
    the interval says the answer is limited by sample size, not by which model
    was chosen.
    """
    names = ("ipw_km", "g_computation", "aipw")
    values = [estimators[name]["risk_difference"] for name in names]
    return float(max(values) - min(values))
