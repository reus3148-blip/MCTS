"""Load and validate transparent assumptions for the dynamic PoC."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class DynamicConfig:
    label: str
    assumption_status: str
    horizon_years: int
    eligibility: Mapping[str, Any]
    response_probabilities: Mapping[str, Mapping[str, float]]
    tumor_size_multipliers: Mapping[str, float]
    acute_toxicity_probabilities: Mapping[str, Mapping[str, float]]
    treatment_burden: Mapping[str, Mapping[str, float]]
    hazard_multipliers: Mapping[str, Any]
    reward: Mapping[str, float]

    @property
    def max_followup_reward(self) -> float:
        return self.horizon_years * (
            float(self.reward["alive_year"])
            + float(self.reward["recurrence_free_year"])
        )

    def normalized(self, raw_reward: float) -> float:
        return raw_reward / self.max_followup_reward


def _validate_probabilities(config: DynamicConfig) -> None:
    for intensity, probabilities in config.response_probabilities.items():
        total = sum(float(value) for value in probabilities.values())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                f"response probabilities for {intensity!r} sum to {total}"
            )
        if any(not 0 <= float(value) <= 1 for value in probabilities.values()):
            raise ValueError(f"invalid response probability for {intensity!r}")

    for treatment, actions in config.acute_toxicity_probabilities.items():
        if any(not 0 <= float(value) <= 1 for value in actions.values()):
            raise ValueError(f"invalid toxicity probability for {treatment!r}")


def load_dynamic_config(path: str | Path) -> DynamicConfig:
    source = Path(path)
    with source.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    config = DynamicConfig(**raw)
    if config.horizon_years < 1:
        raise ValueError("horizon_years must be positive")
    if config.max_followup_reward <= 0:
        raise ValueError("maximum follow-up reward must be positive")
    _validate_probabilities(config)
    return config

