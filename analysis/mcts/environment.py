"""Finite treatment-path environment used by the MCTS PoC.

The patient context remains fixed in v0.1. A state is therefore the prefix of
decisions already made. The terminal reward is supplied by the fitted survival
model for the completed four-decision plan.
"""

from __future__ import annotations

from itertools import product
from typing import Mapping, Sequence

DECISIONS = ("surgery", "chemo", "hormone", "radio")
ACTION_SPACE = {
    "surgery": ("BCS", "MAST"),
    "chemo": (0, 1),
    "hormone": (0, 1),
    "radio": (0, 1),
}

Action = object
Plan = tuple[Action, ...]


def all_plans() -> tuple[Plan, ...]:
    """Return the 16 complete treatment plans in stable order."""
    spaces = [ACTION_SPACE[decision] for decision in DECISIONS]
    return tuple(tuple(values) for values in product(*spaces))


def feasible_plans_for_subtype(subtype: str) -> tuple[Plan, ...]:
    """Apply the receptor-based hard eligibility rule used in PoC v0.1.

    Endocrine therapy is not a meaningful action for HR-negative disease, so
    those patients are searched over the eight plans with ``hormone=0``.
    Other treatment choices remain available for the search to compare.
    """
    plans = all_plans()
    if subtype in {"HR-/HER2+", "TNBC"}:
        return tuple(plan for plan in plans if plan[2] == 0)
    return plans


def plan_to_dict(plan: Sequence[Action]) -> dict[str, Action]:
    """Map a complete or partial plan to named decisions."""
    return dict(zip(DECISIONS, plan, strict=False))


def plan_to_label(plan: Sequence[Action]) -> str:
    """Create a compact, stable plan label for tables and logs."""
    values = plan_to_dict(plan)
    if len(values) != len(DECISIONS):
        return " / ".join(str(value) for value in plan)
    return (
        f"{values['surgery']} | C{int(values['chemo'])} | "
        f"H{int(values['hormone'])} | R{int(values['radio'])}"
    )


class TreatmentPlanningEnvironment:
    """Deterministic four-step environment with cached terminal rewards."""

    def __init__(self, terminal_rewards: Mapping[Plan, float]) -> None:
        expected = set(all_plans())
        provided = set(terminal_rewards)
        if not provided:
            raise ValueError("at least one feasible terminal plan is required")
        if not provided.issubset(expected):
            raise ValueError(f"unknown terminal plans: {provided - expected}")
        self._terminal_rewards = {
            tuple(plan): float(reward)
            for plan, reward in terminal_rewards.items()
        }
        self._terminal_plans = tuple(
            plan for plan in all_plans() if plan in provided
        )

    @property
    def terminal_plans(self) -> tuple[Plan, ...]:
        return self._terminal_plans

    @staticmethod
    def is_terminal(state: Plan) -> bool:
        return len(state) == len(DECISIONS)

    def legal_actions(self, state: Plan) -> tuple[Action, ...]:
        if len(state) > len(DECISIONS):
            raise ValueError("state is longer than the decision horizon")
        if len(state) == len(DECISIONS):
            return ()
        if not any(plan[:len(state)] == state for plan in self._terminal_plans):
            raise ValueError(f"state is not a feasible plan prefix: {state!r}")
        decision = DECISIONS[len(state)]
        return tuple(
            action
            for action in ACTION_SPACE[decision]
            if any(
                plan[:len(state)] == state and plan[len(state)] == action
                for plan in self._terminal_plans
            )
        )

    def transition(self, state: Plan, action: Action) -> Plan:
        legal = self.legal_actions(state)
        if action not in legal:
            raise ValueError(f"illegal action {action!r}; expected one of {legal!r}")
        return (*state, action)

    def reward(self, state: Plan) -> float:
        if not self.is_terminal(state):
            raise ValueError("reward is defined only for terminal states")
        if state not in self._terminal_rewards:
            raise ValueError(f"terminal state is not feasible: {state!r}")
        return self._terminal_rewards[state]
