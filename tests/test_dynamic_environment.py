from __future__ import annotations

import dataclasses
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


def make_environment_with(config, patient: PatientProfile | None = None):
    risks = {
        plan: RiskEstimate(five_year_os=0.82, five_year_rfs=0.75)
        for plan in all_plans()
    }
    return DynamicBreastCancerEnvironment(patient or make_patient(), risks, config)


def replace_discount(config, rate: float):
    return dataclasses.replace(config, discount_rate_annual=rate)


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


class ResponseChannelNeutralityTests(unittest.TestCase):
    """The comparison is void if one policy gets an undeclared head start.

    Only a neoadjuvant patient ever draws a response, so the response hazard
    channel fires for that arm alone. Unless its expectation is 1.0, choosing
    neoadjuvant buys a hazard advantage no assumption ever stated.
    """

    def setUp(self) -> None:
        self.config = load_dynamic_config(
            ROOT / "configs" / "dynamic_poc_v0_2.json")
        self.environment = make_environment()

    def _followup_state(self, timing: str, response: str, chemo: str):
        return DynamicState(
            phase="followup",
            current_tumor_size_mm=20.0,
            timing=timing,
            surgery="MAST",
            chemo=chemo,
            endocrine="standard",
            radiation="none",
            response=response,
            year=0,
        )

    def test_raw_response_multipliers_are_not_neutral(self) -> None:
        """Documents the bias the environment has to correct for."""
        mean = self.config.response_multiplier_mean("standard", "recurrence")
        self.assertLess(mean, 1.0)

    def test_expected_hazards_match_between_timing_arms(self) -> None:
        """Averaged over responses, neoadjuvant must not beat surgery-first."""
        for intensity in ("standard", "intensified"):
            baseline = self.environment.annual_event_probabilities(
                self._followup_state("surgery_first", "not_applicable", intensity))
            probabilities = self.config.response_probabilities[intensity]
            for index, outcome in enumerate(("death", "recurrence")):
                expected = 0.0
                for response, probability in probabilities.items():
                    arm = self.environment.annual_event_probabilities(
                        self._followup_state("neoadjuvant", response, intensity))
                    expected += float(probability) * arm[index]
                # Hazards are mean-neutral by construction; the hazard ->
                # probability transform is convex, so allow a small gap.
                self.assertAlmostEqual(
                    expected, baseline[index], delta=0.002,
                    msg=f"{intensity}/{outcome} arms differ",
                )

    def test_switching_neutralisation_off_reproduces_the_defect(self) -> None:
        """Regression guard: this is the bias the audit found, kept measurable."""
        biased_config = dataclasses.replace(
            self.config, response_channel_neutralised=False)
        biased = make_environment_with(biased_config)
        probabilities = self.config.response_probabilities["standard"]

        def expected_recurrence(environment) -> float:
            return sum(
                float(probability) * environment.annual_event_probabilities(
                    self._followup_state("neoadjuvant", response, "standard"))[1]
                for response, probability in probabilities.items()
            )

        baseline = self.environment.annual_event_probabilities(
            self._followup_state("surgery_first", "not_applicable", "standard"))[1]
        # Off: choosing neoadjuvant buys a recurrence discount nobody declared.
        self.assertLess(expected_recurrence(biased), baseline * 0.97)
        # On: the advantage is gone.
        self.assertAlmostEqual(
            expected_recurrence(self.environment), baseline, delta=0.002)

    def test_relative_ordering_of_responses_survives(self) -> None:
        """Neutralising the mean must not flatten major vs none."""
        major, none = (
            self.environment.annual_event_probabilities(
                self._followup_state("neoadjuvant", response, "standard"))[1]
            for response in ("major", "none")
        )
        self.assertLess(major, none)

    def test_timing_channel_carries_a_declared_effect(self) -> None:
        """A real neoadjuvant effect belongs in an explicit config channel."""
        config = load_dynamic_config(ROOT / "configs" / "dynamic_v0_5.json")
        environment = make_environment_with(config)
        neutral = environment.annual_event_probabilities(
            self._followup_state("surgery_first", "not_applicable", "standard"))
        self.assertGreater(neutral[0], 0.0)
        self.assertIn("timing", config.hazard_multipliers)


class DiscountingTests(unittest.TestCase):
    def test_zero_rate_reproduces_undiscounted_rewards(self) -> None:
        config = load_dynamic_config(ROOT / "configs" / "dynamic_poc_v0_2.json")
        self.assertEqual(config.discount_rate_annual, 0.0)
        self.assertEqual(config.discount_factor(5), 1.0)

    def test_positive_rate_shrinks_later_years(self) -> None:
        config = load_dynamic_config(ROOT / "configs" / "dynamic_poc_v0_2.json")
        discounted = replace_discount(config, 0.03)
        self.assertAlmostEqual(discounted.discount_factor(1), 1 / 1.03)
        self.assertLess(discounted.discount_factor(5), discounted.discount_factor(1))


if __name__ == "__main__":
    unittest.main()
