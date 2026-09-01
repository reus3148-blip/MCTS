"""Inverse-probability weighting building blocks for the target trial emulation.

Everything the emulation needs - a propensity model, balance diagnostics,
stabilized weights, a weighted survival curve and an unmeasured-confounding
sensitivity value - implemented on numpy alone.

Two reasons for not pulling in scikit-learn or statsmodels. The estimates here
are the first causal numbers this project will publish, so every step should be
readable and unit-testable rather than hidden behind a fit() call. And the repo
so far runs on four packages; a causal prototype is a poor excuse to double that.

Nothing here knows about breast cancer. ``analysis/15_*`` supplies the cohort.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def sigmoid(values: np.ndarray) -> np.ndarray:
    """Numerically stable logistic function."""
    out = np.empty_like(values, dtype=float)
    positive = values >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_negative = np.exp(values[~positive])
    out[~positive] = exp_negative / (1.0 + exp_negative)
    return out


def fit_logistic(
    design: np.ndarray,
    outcome: np.ndarray,
    l2: float = 1e-3,
    max_iter: int = 100,
    tol: float = 1e-9,
) -> np.ndarray:
    """Ridge-penalised logistic regression by iteratively reweighted least squares.

    ``design`` must already carry an intercept in column 0; the penalty skips it,
    so shifting the outcome prevalence does not shrink the intercept. The small
    default penalty exists to keep the fit finite under near-separation, which a
    propensity model with strong indication effects can easily hit.
    """
    if design.ndim != 2:
        raise ValueError("design must be 2-dimensional")
    if design.shape[0] != outcome.shape[0]:
        raise ValueError("design and outcome must have the same number of rows")
    if not np.all(np.isin(outcome, (0, 1))):
        raise ValueError("outcome must be binary 0/1")

    n_features = design.shape[1]
    penalty = np.eye(n_features) * l2
    penalty[0, 0] = 0.0
    beta = np.zeros(n_features)

    for _ in range(max_iter):
        eta = design @ beta
        mu = sigmoid(eta)
        # Clip the IRLS working weights away from zero; without it a
        # near-separated fit divides by ~0 and the step explodes.
        w = np.clip(mu * (1.0 - mu), 1e-8, None)
        z = eta + (outcome - mu) / w
        lhs = design.T @ (design * w[:, None]) + penalty
        rhs = design.T @ (w * z)
        step = np.linalg.solve(lhs, rhs)
        if np.max(np.abs(step - beta)) < tol:
            beta = step
            break
        beta = step
    return beta


def standardized_mean_difference(
    values: np.ndarray,
    treated: np.ndarray,
    weights: np.ndarray | None = None,
) -> float:
    """Standardised difference between arms; |SMD| < 0.1 is the usual target.

    Uses the *unweighted* pooled standard deviation as the denominator in both
    the crude and the weighted version, so the before/after numbers of a Love
    plot are on one scale and the weighting cannot flatter itself by inflating
    the spread it divides by.
    """
    values = np.asarray(values, dtype=float)
    treated = np.asarray(treated).astype(bool)
    if weights is None:
        weights = np.ones_like(values)
    weights = np.asarray(weights, dtype=float)

    def weighted_mean(mask: np.ndarray) -> float:
        total = weights[mask].sum()
        if total <= 0:
            return float("nan")
        return float((values[mask] * weights[mask]).sum() / total)

    pooled = math.sqrt(
        (values[treated].var(ddof=1) + values[~treated].var(ddof=1)) / 2.0
    )
    if pooled == 0:
        return 0.0
    return (weighted_mean(treated) - weighted_mean(~treated)) / pooled


def stabilized_weights(treated: np.ndarray, propensity: np.ndarray) -> np.ndarray:
    """Stabilised IPT weights: marginal treatment probability over ``propensity``.

    Stabilising keeps the weights near 1 and the effective sample size close to
    the real one; unstabilised 1/e weights blow up for the few treated patients
    with a low propensity, which is exactly where this cohort is thin.
    """
    treated = np.asarray(treated).astype(float)
    propensity = np.asarray(propensity, dtype=float)
    if np.any((propensity <= 0) | (propensity >= 1)):
        raise ValueError("propensity must lie strictly inside (0, 1)")
    marginal = treated.mean()
    return np.where(
        treated == 1,
        marginal / propensity,
        (1.0 - marginal) / (1.0 - propensity),
    )


def truncate_weights(
    weights: np.ndarray,
    lower_quantile: float = 0.01,
    upper_quantile: float = 0.99,
) -> np.ndarray:
    """Clip extreme weights, trading a little bias for a lot of variance."""
    if not 0 <= lower_quantile < upper_quantile <= 1:
        raise ValueError("quantiles must satisfy 0 <= lower < upper <= 1")
    weights = np.asarray(weights, dtype=float)
    low, high = np.quantile(weights, [lower_quantile, upper_quantile])
    return np.clip(weights, low, high)


def effective_sample_size(weights: np.ndarray) -> float:
    """Kish's effective sample size - how many patients the weights really buy."""
    weights = np.asarray(weights, dtype=float)
    return float(weights.sum() ** 2 / np.square(weights).sum())


@dataclass(frozen=True)
class SurvivalCurve:
    times: np.ndarray
    survival: np.ndarray

    def at(self, horizon: float) -> float:
        """Survival at ``horizon`` (last value at or before it; 1.0 if before all)."""
        eligible = self.times <= horizon
        if not eligible.any():
            return 1.0
        return float(self.survival[eligible][-1])


def weighted_kaplan_meier(
    durations: np.ndarray,
    events: np.ndarray,
    weights: np.ndarray | None = None,
) -> SurvivalCurve:
    """Kaplan-Meier estimate where each patient counts ``weights`` times.

    Written out rather than delegated so the weighting is visible: at each event
    time the number at risk and the number of events are weighted sums, and the
    product-limit step is unchanged.
    """
    durations = np.asarray(durations, dtype=float)
    events = np.asarray(events, dtype=float)
    if weights is None:
        weights = np.ones_like(durations)
    weights = np.asarray(weights, dtype=float)
    if not (durations.shape == events.shape == weights.shape):
        raise ValueError("durations, events and weights must have equal length")
    if durations.size == 0:
        raise ValueError("weighted_kaplan_meier requires at least one observation")

    order = np.argsort(durations, kind="mergesort")
    durations, events, weights = durations[order], events[order], weights[order]

    times, survival = [], []
    running = 1.0
    at_risk = weights.sum()
    index = 0
    total = durations.size
    while index < total:
        time = durations[index]
        tied = index
        event_weight = 0.0
        exiting_weight = 0.0
        while tied < total and durations[tied] == time:
            exiting_weight += weights[tied]
            if events[tied] == 1:
                event_weight += weights[tied]
            tied += 1
        if event_weight > 0 and at_risk > 0:
            running *= 1.0 - event_weight / at_risk
            times.append(time)
            survival.append(running)
        at_risk -= exiting_weight
        index = tied
    return SurvivalCurve(np.asarray(times), np.asarray(survival))


def e_value(risk_ratio: float) -> float:
    """VanderWeele-Ding E-value: the confounding strength that could explain RR.

    An unmeasured confounder would have to be associated with both treatment and
    outcome by at least this risk ratio, above and beyond the measured covariates,
    to move the observed estimate to the null.
    """
    if risk_ratio <= 0:
        raise ValueError("risk_ratio must be positive")
    ratio = risk_ratio if risk_ratio >= 1 else 1.0 / risk_ratio
    return ratio + math.sqrt(ratio * (ratio - 1.0))


def e_value_for_interval(lower: float, upper: float) -> float:
    """E-value for a confidence interval, following VanderWeele-Ding.

    When the interval already covers the null there is nothing for an unmeasured
    confounder to explain away, and the convention is to report 1. Otherwise the
    E-value is computed for the confidence limit *closest to the null*, which is
    the weaker, more conservative claim.
    """
    if lower > upper:
        raise ValueError("lower must not exceed upper")
    if lower <= 1.0 <= upper:
        return 1.0
    limit = lower if lower > 1.0 else upper
    return e_value(limit)
