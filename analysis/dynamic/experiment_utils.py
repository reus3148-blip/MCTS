"""Small, testable helpers shared by the robustness and sensitivity studies.

Keeping these here (rather than inside the numbered ``analysis/1x_*.py`` scripts,
whose module names start with a digit and cannot be imported) lets the unit
tests exercise them directly.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence


def confidence_interval(values: Sequence[float], z: float = 1.96) -> dict[str, float]:
    """Mean and a normal-approximation confidence interval of the sample mean.

    Uses the standard error of the mean (sd / sqrt(n)); with the default
    ``z = 1.96`` this is a 95% interval. A single value yields a zero-width
    interval rather than an error.
    """
    data = [float(v) for v in values]
    n = len(data)
    if n == 0:
        raise ValueError("confidence_interval requires at least one value")
    mean = sum(data) / n
    if n > 1:
        variance = sum((v - mean) ** 2 for v in data) / (n - 1)
        sd = math.sqrt(variance)
    else:
        sd = 0.0
    stderr = sd / math.sqrt(n)
    half = z * stderr
    return {
        "mean": mean,
        "sd": sd,
        "stderr": stderr,
        "ci95_low": mean - half,
        "ci95_high": mean + half,
        "n_seeds": n,
    }


def rescale_response_major(block: Mapping[str, float], major: float) -> dict[str, float]:
    """Return a response-probability block with ``major`` set to the given value.

    ``partial`` and ``none`` are rescaled proportionally so the three
    probabilities still sum to 1. Raises if ``major`` is outside [0, 1] or the
    remaining mass is degenerate.
    """
    if not 0.0 <= major <= 1.0:
        raise ValueError(f"major must be in [0, 1], got {major}")
    old_rest = float(block["partial"]) + float(block["none"])
    if old_rest <= 0.0:
        raise ValueError("cannot rescale when partial + none is zero")
    remainder = 1.0 - major
    return {
        "major": major,
        "partial": remainder * (float(block["partial"]) / old_rest),
        "none": remainder * (float(block["none"]) / old_rest),
    }
