"""One-at-a-time sensitivity analysis of the synthetic assumptions (v0.3).

The v0.2 report's central caveat is that its outcomes depend on synthetic
transition/toxicity/reward assumptions. This script makes that dependence
explicit: it perturbs each key assumption up and down from the baseline
configuration and measures how the MCTS-vs-NCCN utility gap and the MCTS action
mix respond. Assumptions whose perturbation flips the sign of the gap, or
sharply changes the chosen actions, are the ones that most need real K-CURE
estimates before any clinical reading.

To separate the *assumption* effect from seed noise, every configuration is run
over several seeds and averaged. Still not causal, still not clinical.
"""

from __future__ import annotations

import copy
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.dynamic.cohort import (  # noqa: E402
    BASE_SEED,
    balanced_subtype_sample,
    build_reward_models,
    git_commit,
    input_manifest,
    make_risk_table,
)
from analysis.dynamic.config import DynamicConfig, _validate_probabilities  # noqa: E402
from analysis.dynamic.environment import DynamicBreastCancerEnvironment  # noqa: E402
from analysis.dynamic.experiment_utils import rescale_response_major  # noqa: E402
from analysis.dynamic.evaluation import run_policy_episodes  # noqa: E402
from analysis.dynamic.policies import CachedMCTSPolicy, DynamicNccnPolicy  # noqa: E402
from analysis.dynamic.schema import patient_from_row  # noqa: E402

INPUT_CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
CONFIG_PATH = ROOT / "configs" / "dynamic_poc_v0_2.json"
REPORT_DIR = ROOT / "reports" / "sensitivity-v0.3"
TABLE_DIR = REPORT_DIR / "tables"

RUN_DATE = "2026-08-24"
PATIENTS_PER_SUBTYPE = 2       # 8 patients x variants x seeds -> tractable
N_SEEDS = 3
EPISODES_PER_POLICY = 40
SIMULATIONS = 256
EXPLORATION_WEIGHT = math.sqrt(2.0)


def config_from_dict(data: dict) -> DynamicConfig:
    config = DynamicConfig(**data)
    _validate_probabilities(config)
    return config


def set_response_major(data: dict, intensity: str, major: float) -> None:
    """Set the 'major' response probability and rescale partial/none to sum 1."""
    data["response_probabilities"][intensity] = rescale_response_major(
        data["response_probabilities"][intensity], major
    )


def build_variants(base: dict) -> list[dict]:
    """Return a list of {name, param, value, config_dict} perturbations."""
    variants = [{"name": "baseline", "param": "(none)", "value": "-",
                 "data": copy.deepcopy(base)}]

    def add(name, param, value, mutate):
        data = copy.deepcopy(base)
        mutate(data)
        variants.append({"name": name, "param": param, "value": value, "data": data})

    # 1. Neoadjuvant intensified major-response probability
    for v in (0.20, 0.50):
        add(f"resp_intensified_major={v}", "response.intensified.major", v,
            lambda d, v=v: set_response_major(d, "intensified", v))
    # 2. Intensified chemo survival benefit (death hazard multiplier)
    for v in (0.85, 1.00):
        add(f"intensified_death_hr={v}", "hazard.chemo.intensified.death", v,
            lambda d, v=v: d["hazard_multipliers"]["chemo"]["intensified"]
            .__setitem__("death", v))
    # 3. Intensified chemo toxicity probability
    for v in (0.15, 0.45):
        add(f"intensified_toxicity={v}", "toxicity.chemo.intensified", v,
            lambda d, v=v: d["acute_toxicity_probabilities"]["chemo"]
            .__setitem__("intensified", v))
    # 4. Acute toxicity penalty weight in the reward
    for v in (0.05, 0.30):
        add(f"toxicity_penalty={v}", "reward.acute_toxicity_penalty", v,
            lambda d, v=v: d["reward"].__setitem__("acute_toxicity_penalty", v))
    # 5. Major-response death benefit
    for v in (0.90, 1.00):
        add(f"major_response_death_hr={v}", "hazard.response.major.death", v,
            lambda d, v=v: d["hazard_multipliers"]["response"]["major"]
            .__setitem__("death", v))
    # 6. Recurrence-free reward weight
    for v in (0.25, 1.00):
        add(f"recurrence_free_reward={v}", "reward.recurrence_free_year", v,
            lambda d, v=v: d["reward"].__setitem__("recurrence_free_year", v))
    return variants


def evaluate_config(config: DynamicConfig, sample, os_model, rfs_model) -> dict:
    gaps, neo, inten, tox_gap = [], [], [], []
    for seed_index in range(N_SEEDS):
        seed = BASE_SEED + seed_index * 1_000
        m_util, n_util, m_neo, m_inten = [], [], [], []
        m_tox, n_tox = [], []
        for _, row in sample.iterrows():
            patient = patient_from_row(row)
            env = DynamicBreastCancerEnvironment(
                patient, make_risk_table(row, os_model, rfs_model), config)
            ph = int(hashlib.sha256(patient.patient_id.encode()).hexdigest()[:8], 16)
            ps = seed + ph % 100_000
            mcts = CachedMCTSPolicy(env, simulations=SIMULATIONS,
                                    exploration_weight=EXPLORATION_WEIGHT, seed=ps)
            nccn = DynamicNccnPolicy(env)
            md = pd.DataFrame(run_policy_episodes(env, mcts, EPISODES_PER_POLICY, ps + 10_000))
            nd = pd.DataFrame(run_policy_episodes(env, nccn, EPISODES_PER_POLICY, ps + 10_000))
            m_util.append(md["utility"].mean()); n_util.append(nd["utility"].mean())
            m_neo.append((md["timing"] == "neoadjuvant").mean())
            m_inten.append((md["chemo"] == "intensified").mean())
            m_tox.append(md["toxicity_count"].mean()); n_tox.append(nd["toxicity_count"].mean())
        gaps.append(np.mean(m_util) - np.mean(n_util))
        neo.append(np.mean(m_neo) * 100)
        inten.append(np.mean(m_inten) * 100)
        tox_gap.append(np.mean(m_tox) - np.mean(n_tox))
    return {
        "utility_gap": float(np.mean(gaps)),
        "utility_gap_sd": float(np.std(gaps, ddof=1)) if len(gaps) > 1 else 0.0,
        "neoadjuvant_rate_pct": float(np.mean(neo)),
        "intensified_rate_pct": float(np.mean(inten)),
        "toxicity_gap": float(np.mean(tox_gap)),
    }


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        base = json.load(handle)
    raw = pd.read_csv(INPUT_CSV)
    os_model, rfs_model, os_test = build_reward_models(raw)
    sample = balanced_subtype_sample(os_test, PATIENTS_PER_SUBTYPE)
    print(f"sensitivity sample = {len(sample)} patients")

    variants = build_variants(base)
    rows = []
    baseline_gap = None
    for variant in variants:
        config = config_from_dict(variant["data"])
        result = evaluate_config(config, sample, os_model, rfs_model)
        if variant["name"] == "baseline":
            baseline_gap = result["utility_gap"]
        rows.append({
            "variant": variant["name"],
            "parameter": variant["param"],
            "value": variant["value"],
            **result,
        })
        print(f"  {variant['name']:35s} gap={result['utility_gap']:+.4f} "
              f"neo={result['neoadjuvant_rate_pct']:5.1f}% "
              f"inten={result['intensified_rate_pct']:5.1f}%")

    table = pd.DataFrame(rows)
    table["gap_delta_vs_baseline"] = table["utility_gap"] - baseline_gap
    table["sign_flip"] = (table["utility_gap"] * baseline_gap < 0).astype(int)
    table = table.sort_values("gap_delta_vs_baseline")

    ranked = table[table["variant"] != "baseline"].copy()
    ranked["abs_gap_delta"] = ranked["gap_delta_vs_baseline"].abs()
    influence = (ranked.groupby("parameter", as_index=False)["abs_gap_delta"]
                 .max().sort_values("abs_gap_delta", ascending=False))

    metrics = {
        "run_date": RUN_DATE,
        "analysis_label": "sensitivity-v0.3",
        "interpretation": (
            "How the MCTS-vs-NCCN gap responds to each synthetic assumption; "
            "identifies which assumptions most need real K-CURE estimates."
        ),
        "design": {
            "patients": int(len(sample)), "n_seeds": N_SEEDS,
            "episodes_per_policy_per_seed": EPISODES_PER_POLICY,
            "variants": len(variants),
        },
        "baseline_utility_gap": baseline_gap,
        "any_sign_flip": bool(table["sign_flip"].any()),
        "most_influential_parameters": influence.head(3).to_dict(orient="records"),
    }

    table.to_csv(TABLE_DIR / "sensitivity_grid.csv", index=False)
    influence.to_csv(TABLE_DIR / "parameter_influence.csv", index=False)
    with (REPORT_DIR / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    manifest = {
        "run_date": RUN_DATE,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit_before_run": git_commit(),
        "inputs": input_manifest({
            "data": INPUT_CSV, "assumptions": CONFIG_PATH,
        }),
        "entry_point": "analysis/11_run_sensitivity_analysis.py",
        "base_seed": BASE_SEED,
    }
    with (REPORT_DIR / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\n=== sensitivity grid (sorted by gap change) ===")
    print(table[["variant", "utility_gap", "gap_delta_vs_baseline",
                 "neoadjuvant_rate_pct", "intensified_rate_pct"]].to_string(index=False))
    print("\n=== most influential parameters ===")
    print(influence.to_string(index=False))
    print(f"\nbaseline gap={baseline_gap:+.4f}  any_sign_flip={bool(table['sign_flip'].any())}")
    print(f"saved -> {REPORT_DIR}")


if __name__ == "__main__":
    main()
