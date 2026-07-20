"""Seeded UCT search and an exhaustive oracle for the treatment environment."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import random

from .environment import Plan, TreatmentPlanningEnvironment


@dataclass
class Node:
    state: Plan
    parent: "Node | None" = None
    action: object | None = None
    visits: int = 0
    value_sum: float = 0.0
    children: dict[object, "Node"] = field(default_factory=dict)
    untried_actions: list[object] = field(default_factory=list)

    @property
    def mean_value(self) -> float:
        return self.value_sum / self.visits if self.visits else 0.0


def _new_node(
    environment: TreatmentPlanningEnvironment,
    state: Plan,
    parent: Node | None = None,
    action: object | None = None,
) -> Node:
    return Node(
        state=state,
        parent=parent,
        action=action,
        untried_actions=list(environment.legal_actions(state)),
    )


def _uct_child(node: Node, exploration_weight: float) -> Node:
    log_parent = math.log(max(node.visits, 1))

    def score(child: Node) -> tuple[float, str]:
        exploration = exploration_weight * math.sqrt(
            log_parent / max(child.visits, 1)
        )
        return child.mean_value + exploration, repr(child.action)

    return max(node.children.values(), key=score)


def _rollout(
    environment: TreatmentPlanningEnvironment,
    state: Plan,
    rng: random.Random,
) -> Plan:
    while not environment.is_terminal(state):
        action = rng.choice(environment.legal_actions(state))
        state = environment.transition(state, action)
    return state


def search_action(
    environment: TreatmentPlanningEnvironment,
    prefix: Plan,
    simulations: int,
    exploration_weight: float,
    rng: random.Random,
) -> object:
    """Search from one plan prefix and return the robust-child action."""
    if simulations < 1:
        raise ValueError("simulations must be positive")
    if environment.is_terminal(prefix):
        raise ValueError("cannot search from a terminal state")

    root = _new_node(environment, prefix)
    for _ in range(simulations):
        node = root

        while (
            not environment.is_terminal(node.state)
            and not node.untried_actions
        ):
            node = _uct_child(node, exploration_weight)

        if not environment.is_terminal(node.state) and node.untried_actions:
            action_index = rng.randrange(len(node.untried_actions))
            action = node.untried_actions.pop(action_index)
            child_state = environment.transition(node.state, action)
            child = _new_node(environment, child_state, node, action)
            node.children[action] = child
            node = child

        terminal_state = _rollout(environment, node.state, rng)
        reward = environment.reward(terminal_state)
        while node is not None:
            node.visits += 1
            node.value_sum += reward
            node = node.parent

    if not root.children:
        raise RuntimeError("search produced no child nodes")
    best = max(
        root.children.values(),
        key=lambda child: (
            child.visits,
            child.mean_value,
            repr(child.action),
        ),
    )
    return best.action


def mcts_plan(
    environment: TreatmentPlanningEnvironment,
    simulations_per_step: int = 256,
    exploration_weight: float = math.sqrt(2.0),
    seed: int = 0,
) -> Plan:
    """Build a complete plan with receding-horizon UCT searches."""
    rng = random.Random(seed)
    plan: Plan = ()
    while not environment.is_terminal(plan):
        action = search_action(
            environment,
            plan,
            simulations_per_step,
            exploration_weight,
            rng,
        )
        plan = environment.transition(plan, action)
    return plan


def exhaustive_best_plan(
    environment: TreatmentPlanningEnvironment,
) -> tuple[Plan, float]:
    """Return the exact optimum; used only to validate this small PoC."""
    best_plan = max(
        environment.terminal_plans,
        key=lambda plan: (environment.reward(plan), repr(plan)),
    )
    return best_plan, environment.reward(best_plan)
