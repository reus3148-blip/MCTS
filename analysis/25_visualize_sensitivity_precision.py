"""Figure for the higher-precision sensitivity re-run (fig32).

Three seeds put the baseline gap at -0.002; twelve seeds put it at +0.017 with
nothing else changed. That single shift is larger than every effect the original
ranking was built to compare, so the figure leads with it and only then shows
what the ranking looks like once the noise is out.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REPORT_DIR = ROOT / "reports" / "sensitivity-precision-v1.0"
TABLE_DIR = REPORT_DIR / "tables"
FIGURE_DIR = REPORT_DIR / "figures"
PUBLIC_DIR = ROOT / "public" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PUBLIC_DIR.mkdir(parents=True, exist_ok=True)

matplotlib.rcParams.update({
    "font.family": "Malgun Gothic",
    "axes.unicode_minus": False,
    "figure.facecolor": "#fafaf9",
    "axes.facecolor": "#fafaf9",
    "axes.edgecolor": "#292524",
    "axes.labelcolor": "#292524",
    "axes.titlecolor": "#1c1917",
    "xtick.color": "#57534e",
    "ytick.color": "#57534e",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e7e5e4",
    "grid.linewidth": 0.6,
})

INK = "#1c1917"
BLUE = "#2563eb"
GREEN = "#15803d"
RED = "#dc2626"
GRAY = "#78716c"

METRICS = json.loads((REPORT_DIR / "metrics.json").read_text(encoding="utf-8"))

PRETTY = {
    "reward.recurrence_free_year": "무재발 1년의 가치",
    "reward.acute_toxicity_penalty": "급성 독성 페널티",
    "toxicity.chemo.intensified": "강화항암 독성확률",
    "hazard.chemo.intensified.death": "강화항암 사망 HR",
    "response.intensified.major": "강화항암 major 반응확률",
    "hazard.response.major.death": "major 반응 사망 HR",
}
JUDGEMENTS = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig32_precision() -> None:
    detail = pd.read_csv(TABLE_DIR / "variant_results.csv")
    influence = pd.read_csv(TABLE_DIR / "parameter_influence.csv")
    by_budget = {row["budget"]: row for row in METRICS["by_budget"]}

    # The ranking panels carry long Korean labels, so they need real gutters or
    # the tick text lands on the neighbouring bars.
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8),
                             gridspec_kw={"wspace": 0.55})

    # --- baseline shift ----------------------------------------------------
    ax = axes[0]
    base_rows = [
        ("시드 3\n예산 256", METRICS["v0_5env_reference"]["baseline_utility_gap"], GRAY),
        ("시드 12\n예산 256", by_budget["256"]["baseline_utility_gap"], BLUE),
        ("시드 12\n예산 1024", by_budget["1024"]["baseline_utility_gap"], GREEN),
    ]
    positions = np.arange(len(base_rows))
    bars = ax.bar(positions, [row[1] for row in base_rows],
                  color=[row[2] for row in base_rows], width=0.6)
    for bar, (_, value, _) in zip(bars, base_rows):
        ax.text(bar.get_x() + bar.get_width() / 2,
                value + (0.0012 if value >= 0 else -0.0028),
                f"{value:+.4f}", ha="center", fontsize=10, color=INK)
    ax.axhline(0, color=INK, linewidth=1.1)
    ax.set_xticks(positions)
    ax.set_xticklabels([row[0] for row in base_rows], fontsize=9.5)
    ax.set_ylabel("기준선 효용 격차")
    shift = (by_budget["256"]["baseline_utility_gap"]
             - METRICS["v0_5env_reference"]["baseline_utility_gap"])
    ax.set_title(f"시드만 4배로 늘렸는데 {shift:+.4f} 움직였다",
                 fontsize=12, pad=10, color=RED)

    # --- ranking at each budget --------------------------------------------
    for column, budget in enumerate(("256", "1024"), start=1):
        ax = axes[column]
        block = influence[influence["budget"].astype(str) == budget]
        block = block.sort_values("abs_gap_delta")
        names = [PRETTY.get(p, p) for p in block["parameter"]]
        colors = [RED if p in JUDGEMENTS else BLUE for p in block["parameter"]]
        positions = np.arange(len(block))
        ax.barh(positions, block["abs_gap_delta"], color=colors, height=0.6)
        ax.errorbar(block["abs_gap_delta"], positions,
                    xerr=1.96 * block["delta_standard_error"], fmt="none",
                    ecolor=INK, capsize=3, linewidth=1.1)
        ax.set_yticks(positions)
        ax.set_yticklabels(names, fontsize=9)
        ax.set_xlabel("효용 격차에 준 최대 영향 |Δ|")
        ax.set_xlim(0, float(influence["abs_gap_delta"].max()
                             + 1.96 * influence["delta_standard_error"].max()) * 1.1)
        top_two = by_budget[budget]["top_two_are_value_judgements"]
        verdict = "상위 2개가 가치판단" if top_two else "상위 2개가 가치판단 아님"
        ax.set_title(f"시드 12 · 예산 {budget}\n{verdict}",
                     fontsize=11.5, pad=10,
                     color=GREEN if top_two else RED)

    handles = [plt.Rectangle((0, 0), 1, 1, color=RED),
               plt.Rectangle((0, 0), 1, 1, color=BLUE)]
    axes[2].legend(handles, ["가치판단", "데이터로 추정 가능"],
                   frameon=False, fontsize=8.5, loc="lower right")

    fig.suptitle("v1.0 — v0.3의 순위는 예산 1024에서만 성립한다 (오차막대는 95% CI)",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig32_sensitivity_precision.png")


def main() -> None:
    fig32_precision()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
