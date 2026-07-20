from __future__ import annotations

import random
import unittest
from pathlib import Path

from analysis.dynamic.config import load_dynamic_config
from analysis.dynamic.environment import DynamicBreastCancerEnvironment
from analysis.dynamic.evaluation import simulate_episode
from analysis.dynamic.policies import DynamicNccnPolicy
from analysis.dynamic.schema import DynamicState, PatientProfile, RiskEstimate
from analysis.dynamic.search import stochastic_mcts_search
from analysis.mcts.environment import all_plans

ROOT = Path(__file__).resolve().parent.parent


def make_patient(subtype: str = "HR+/HER2-") -> PatientProfile:
    hr_positive = subtype.startswith("HR+")
    her2_positive = "HER2+" in subtype
    return PatientProfile(
        patient_id="TEST-001",
        age=55.0,
        menopause="Post",
        tumor_size_mm=20.0,
        lymph_pos=0,
        stage=1,
        grade=2,
        subtype=subtype,
        er=int(hr_positive),
        pr=int(hr_positive),
        her2=int(her2_positive),
    )


def make_environment(patient: PatientProfile | None = None):
    config = load_dynamic_config(ROOT / "configs" / "dynamic_poc_v0_2.json")
    risks = {
        plan: RiskEstimate(five_year_os=0.82, five_year_rfs=0.75)
        for plan in all_plans()
    }
    return DynamicBreastCancerEnvironment(patient or make_patient(), risks, config)


class SequenceRandom:
    def __init__(self, values):
        self.values = iter(values)

    def random(self):
        return next(self.values)


class DynamicEnvironmentTests(unittest.TestCase):
    def test_hr_negative_endocrine_action_is_masked(self) -> None:
        environment = make_environment(make_patient("TNBC"))
        state = DynamicState(
            phase="endocrine",
            current_tumor_size_mm=20.0,
            timing="surgery_first",
            surgery="BCS",
            chemo="standard",
        )
        self.assertEqual(environment.legal_actions(state), ("none",))

    def test_neoadjuvant_response_can_open_bcs_action(self) -> None:
        patient = make_patient()
        patient = PatientProfile(**{
            **patient.__dict__,
            "tumor_size_mm": 40.0,
            "stage": 2,
        })
        environment = make_environment(patient)
        state = environment.initial_state()
        state, _, _ = environment.step(
            state,
            "neoadjuvant",
            random.Random(1),
        )
        state, _, info = environment.step(
            state,
            "intensified",
            SequenceRandom([0.99, 0.10]),
        )
        self.assertEqual(info["response"], "major")
        self.assertEqual(state.current_tumor_size_mm, 20.0)
        self.assertIn("BCS", environment.legal_actions(state))

    def test_nccn_episode_reaches_a_terminal_five_year_state(self) -> None:
        environment = make_environment()
        result, _ = simulate_episode(
            environment,
            DynamicNccnPolicy(environment),
            seed=20260711,
        )
        self.assertIn(result["survived_5y"], (0, 1))
        self.assertLessEqual(result["terminal_year"], 5)


class ToyStochasticEnvironment:
    @staticmethod
    def is_terminal(state):
        return state != "root"

    @staticmethod
    def legal_actions(state):
        return ("good", "bad") if state == "root" else ()

    @staticmethod
    def step(state, action, rng):
        reward = 1.0 if action == "good" else 0.0
        return action, reward, {}


class StochasticSearchTests(unittest.TestCase):
    def test_search_prefers_the_higher_expected_reward(self) -> None:
        result = stochastic_mcts_search(
            ToyStochasticEnvironment(),
            "root",
            simulations=64,
            seed=7,
        )
        self.assertEqual(result.action, "good")
        self.assertGreater(result.action_values["good"], result.action_values["bad"])


if __name__ == "__main__":
    unittest.main()

