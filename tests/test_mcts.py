from __future__ import annotations

import unittest

from analysis.mcts.environment import (
    TreatmentPlanningEnvironment,
    all_plans,
    feasible_plans_for_subtype,
)
from analysis.mcts.search import exhaustive_best_plan, mcts_plan


class TreatmentEnvironmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.optimum = ("MAST", 1, 0, 1)
        rewards = {
            plan: 0.25 + index / 1000
            for index, plan in enumerate(all_plans())
        }
        rewards[self.optimum] = 0.99
        self.environment = TreatmentPlanningEnvironment(rewards)

    def test_action_space_has_sixteen_complete_plans(self) -> None:
        plans = all_plans()
        self.assertEqual(len(plans), 16)
        self.assertEqual(len(set(plans)), 16)

    def test_transitions_follow_the_fixed_decision_order(self) -> None:
        state = self.environment.transition((), "BCS")
        self.assertEqual(self.environment.legal_actions(state), (0, 1))
        state = self.environment.transition(state, 1)
        state = self.environment.transition(state, 0)
        state = self.environment.transition(state, 1)
        self.assertTrue(self.environment.is_terminal(state))

    def test_hr_negative_patients_cannot_receive_hormone_action(self) -> None:
        feasible = feasible_plans_for_subtype("TNBC")
        environment = TreatmentPlanningEnvironment({
            plan: 0.5 for plan in feasible
        })
        prefix = ("BCS", 1)
        self.assertEqual(environment.legal_actions(prefix), (0,))
        self.assertEqual(len(environment.terminal_plans), 8)

    def test_mcts_finds_the_exhaustive_optimum(self) -> None:
        exact, _ = exhaustive_best_plan(self.environment)
        searched = mcts_plan(
            self.environment,
            simulations_per_step=512,
            seed=20260711,
        )
        self.assertEqual(exact, self.optimum)
        self.assertEqual(searched, exact)

    def test_mcts_is_reproducible_with_a_fixed_seed(self) -> None:
        first = mcts_plan(self.environment, simulations_per_step=64, seed=17)
        second = mcts_plan(self.environment, simulations_per_step=64, seed=17)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
