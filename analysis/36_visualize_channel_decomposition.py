"""Figure for the declared-channel decomposition (fig40).

The finding is an asymmetry in what we ourselves wrote down, so the figure has to
show the config first: three standard-of-care levels carry a declared cost and a
declared benefit of exactly zero. Then the 2x2 that measures what that does, and
finally the arm-by-arm utilities showing which policy moves.

Every number is read from the committed ``metrics.json`` and ``tables/``.
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

REPORT_DIR = ROOT / "reports" / "channel-decomposition-v1.5"
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
AMBER = "#d97706"
GRAY = "#78716c"

METRICS = json.loads((REPORT_DIR / "metrics.json").read_text(encoding="utf-8"))

PRETTY_LEVEL = {
    ("chemo", "standard"): "표준 항암",
    ("chemo", "intensified"): "강화 항암",
    ("endocrine", "standard"): "표준 호르몬",
    ("endocrine", "extended"): "연장 호르몬",
    ("radiation", "local"): "국소 방사선",
    ("radiation", "regional"): "국소+림프절",
}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig40_decomposition() -> None:
    verdict = METRICS["verdict"]
    arms = METRICS["arms"]
    asymmetry = pd.DataFrame(METRICS["declared_asymmetry"])
    asymmetry = asymmetry[asymmetry["level"] != "none"].reset_index(drop=True)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4),
                             gridspec_kw={"wspace": 0.32,
                                          "width_ratios": [1.25, 1, 1]})

    # --- A: what the config declares ---------------------------------------
    ax = axes[0]
    # Benefit as the recurrence-hazard reduction, cost as toxicity + burden, both
    # on a "how much does this move utility" footing rather than raw units.
    benefit = (1.0 - asymmetry["recurrence_hazard"]) * 100
    cost = (asymmetry["acute_toxicity_probability"]
            + asymmetry["treatment_burden"]) * 100
    positions = np.arange(len(asymmetry))
    width = 0.36
    ax.bar(positions - width / 2, benefit, width, color=GREEN,
           label="선언된 이득 (재발위험 감소 %)")
    ax.bar(positions + width / 2, cost, width, color=RED,
           label="선언된 비용 (독성확률 + 부담) %")
    for position, (b, row) in enumerate(zip(benefit, asymmetry.itertuples())):
        if not row.declares_any_benefit:
            ax.text(position - width / 2, 0.6, "0", ha="center", fontsize=11,
                    color=RED, fontweight="bold")
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [PRETTY_LEVEL.get((row.channel, row.level), f"{row.channel}.{row.level}")
         for row in asymmetry.itertuples()], fontsize=8.5, rotation=20,
        ha="right")
    ax.set_ylabel("퍼센트")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    pure = len(verdict["pure_cost_levels"])
    ax.set_title(f"A. 표준치료 {pure}종은 **순수 비용**이다\n"
                 f"우리가 config에 이득을 0으로 적어 뒀다".replace("**", ""),
                 fontsize=12, pad=10, color=RED)

    # --- B: the 2x2 ---------------------------------------------------------
    ax = axes[1]
    keys = ["baseline", "cost_off", "benefit_off", "both_off"]
    labels = ["A 그대로", "B 비용 제거", "C 이득 제거", "D 둘 다 제거\n(영대조)"]
    colors = [GRAY, GREEN, AMBER, BLUE]
    gaps = [arms[key]["utility_gap"] for key in keys]
    errors = [1.96 * arms[key]["standard_error"] for key in keys]
    positions = np.arange(len(keys))
    ax.bar(positions, gaps, yerr=errors, color=colors, width=0.6,
           error_kw={"ecolor": INK, "capsize": 4, "linewidth": 1.1})
    for x, gap, error in zip(positions, gaps, errors):
        offset = error + 0.0012 if gap >= 0 else -error - 0.0032
        ax.text(x, gap + offset, f"{gap:+.4f}", ha="center", fontsize=10.5,
                color=INK)
    ax.axhline(0, color=INK, linewidth=1.2)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("MCTS - NCCN 효용 격차")
    share = verdict["share_of_gap_from_declared_cost"]
    passed = verdict["null_control_passed"]
    ax.set_title(f"B. 비용을 없애면 격차의 {share:.0%}가 사라진다\n"
                 f"영대조 {'통과' if passed else '실패'} "
                 f"({verdict['gap_both_off']:+.4f})",
                 fontsize=12, pad=10, color=RED if share > 0.5 else INK)

    # --- C: which policy moves ----------------------------------------------
    ax = axes[2]
    positions = np.arange(len(keys))
    width = 0.36
    mcts = [arms[key]["mcts_utility"] for key in keys]
    nccn = [arms[key]["nccn_utility"] for key in keys]
    ax.bar(positions - width / 2, mcts, width, color=BLUE, label="MCTS")
    ax.bar(positions + width / 2, nccn, width, color=AMBER, label="NCCN")
    low = min(mcts + nccn)
    ax.set_ylim(low - 0.02, max(mcts + nccn) + 0.012)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("평균 효용")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title("C. 비용을 없애면 NCCN이 올라온다\n"
                 "MCTS는 애초에 그 비용을 내지 않고 있었다",
                 fontsize=12, pad=10, color=INK)

    fig.suptitle("v1.5 — 남은 격차는 무엇으로 이루어져 있나 (오차막대는 95% CI)",
                 fontsize=13.5, y=1.03, color=INK)
    save(fig, "fig40_channel_decomposition.png")


def main() -> None:
    fig40_decomposition()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
