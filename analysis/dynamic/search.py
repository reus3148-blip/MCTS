"""Chance-aware UCT search for the stochastic dynamic environment."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import random
from typing import Hashable, Protocol


class StochasticEnvironment(Protocol):
    def is_terminal(self, state: Hashable) -> bool: ...
    def legal_actions(self, state: Hashable) -> tuple[str, ...]: ...
    def step(
        self,
        state: Hashable,
        action: str,
        rng: random.Random,
    ) -> tuple[Hashable, float, dict[str, object]]: ...


@dataclass
class ActionStatistics:
    visits: int = 0
    value_sum: float = 0.0
    outcomes: dict[Hashable, "DecisionNode"] = field(default_factory=dict)

    @property
    def mean_value(self) -> float:
        return self.value_sum / self.visits if self.visits else 0.0


@dataclass
class DecisionNode:
    state: Hashable
    visits: int = 0
    actions: dict[str, ActionStatistics] = field(default_factory=dict)


@dataclass(frozen=True)
class StochasticSearchResult:
    action: str
    action_values: dict[str, float]
    action_visits: dict[str, int]
    simulations: int


def _rollout(
    environment: StochasticEnvironment,
    state: Hashable,
    rng: random.Random,
    max_steps: int = 32,
) -> float:
    total = 0.0
    steps = 0
    while not environment.is_terminal(state):
        if steps >= max_steps:
            raise RuntimeError("rollout exceeded the dynamic environment horizon")
        actions = environment.legal_actions(state)
        action = rng.choice(actions)
        state, reward, _ = environment.step(state, action, rng)
        total += reward
        steps += 1
    return total


def _select_action(
    node: DecisionNode,
    exploration_weight: float,
    rng: random.Random,
) -> str:
    unvisited = [
        action
        for action, statistics in node.actions.items()
        if statistics.visits == 0
    ]
    if unvisited:
        return rng.choice(sorted(unvisited))

    log_parent = math.log(max(node.visits, 1))

    def score(action: str) -> tuple[float, str]:
        statistics = node.actions[action]
        exploration = exploration_weight * math.sqrt(
            log_parent / statistics.visits
        )
        return statistics.mean_value + exploration, action

    return max(node.actions, key=score)


def _simulate(
    environment: StochasticEnvironment,
    node: DecisionNode,
    exploration_weight: float,
    rng: random.Random,
) -> float:
    if environment.is_terminal(node.state):
        return 0.0

    legal = environment.legal_actions(node.state)
    for action in legal:
        node.actions.setdefault(action, ActionStatistics())
    action = _select_action(node, exploration_weight, rng)
    statistics = node.actions[action]
    next_state, immediate_reward, _ = environment.step(node.state, action, rng)

    if environment.is_terminal(next_state):
        future_reward = 0.0
    elif next_state not in statistics.outcomes:
        statistics.outcomes[next_state] = DecisionNode(next_state)
        future_reward = _rollout(environment, next_state, rng)
    else:
        future_reward = _simulate(
            environment,
            statistics.outcomes[next_state],
            exploration_weight,
            rng,
        )

    total_reward = immediate_reward + future_reward
    statistics.visits += 1
    statistics.value_sum += total_reward
    node.visits += 1
    return total_reward


def stochastic_mcts_search(
    environment: StochasticEnvironment,
    state: Hashable,
    simulations: int = 128,
    exploration_weight: float = math.sqrt(2.0),
    seed: int = 0,
) -> StochasticSearchResult:
    """Return the robust-child action and root action diagnostics."""
    if environment.is_terminal(state):
        raise ValueError("cannot search from a terminal state")
    if simulations < 1:
        raise ValueError("simulations must be positive")

    root = DecisionNode(state)
    rng = random.Random(seed)
    for _ in range(simulations):
        _simulate(environment, root, exploration_weight, rng)

    best_action = max(
        root.actions,
        key=lambda action: (
            root.actions[action].visits,
            root.actions[action].mean_value,
            action,
        ),
    )
    return StochasticSearchResult(
        action=best_action,
        action_values={
            action: statistics.mean_value
            for action, statistics in root.actions.items()
        },
        action_visits={
            action: statistics.visits
            for action, statistics in root.actions.items()
        },
        simulations=simulations,
    )

