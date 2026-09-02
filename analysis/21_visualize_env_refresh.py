"""Figure for the v0.5-environment refresh (fig30).

v0.3 and v0.4 ran on the environment that quietly rewarded neoadjuvant therapy.
v0.5 measured that bias as far too small for the design to resolve, but "too
small to see" is a prediction, not a result. Re-running both studies on the
corrected environment tests it, and this figure puts the two sets of headline
numbers side by side so the answer is one glance rather than four files.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REPORT_DIR = ROOT / "reports" / "robustness-v0.5env"
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
GRAY = "#78716c"
GREEN = "#15803d"


def load(label: str) -> dict:
    return json.loads(
        (ROOT / "reports" / label / "metrics.json").read_text(encoding="utf-8"))


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig30_environment_refresh() -> None:
    old_robust, new_robust = load("robustness-v0.3"), load("robustness-v0.5env")
    old_budget, new_budget = load("budget-scaling-v0.4"), load("budget-scaling-v0.5env")
    old_sens, new_sens = load("sensitivity-v0.3"), load("sensitivity-v0.5env")

    def gap_at(metrics: dict, budget: int) -> dict:
        return next(row for row in metrics["utility_gap_by_budget"]
                    if row["budget"] == budget)

    rows = [
        ("효용 격차\n(256 · 20시드)",
         old_robust["utility_gap_ci95"]["mean"], new_robust["utility_gap_ci95"]["mean"],
         (old_robust["utility_gap_ci95"]["ci95_low"],
          old_robust["utility_gap_ci95"]["ci95_high"]),
         (new_robust["utility_gap_ci95"]["ci95_low"],
          new_robust["utility_gap_ci95"]["ci95_high"])),
        ("효용 격차\n(1024 · 10시드)",
         gap_at(old_budget, 1024)["utility_gap_mean"],
         gap_at(new_budget, 1024)["utility_gap_mean"],
         (gap_at(old_budget, 1024)["utility_gap_ci95_low"],
          gap_at(old_budget, 1024)["utility_gap_ci95_high"]),
         (gap_at(new_budget, 1024)["utility_gap_ci95_low"],
          gap_at(new_budget, 1024)["utility_gap_ci95_high"])),
        ("민감도 기준선\n(256 · 8명 · 3시드)",
         old_sens["baseline_utility_gap"], new_sens["baseline_utility_gap"],
         None, None),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.6))

    ax = axes[0]
    positions = np.arange(len(rows))
    for position, (_, old, new, old_ci, new_ci) in zip(positions, rows):
        for value, ci, offset, color, label in (
            (old, old_ci, +0.16, GRAY, "편향 환경 (v0.3/v0.4)"),
            (new, new_ci, -0.16, GREEN, "수정 환경 (v0.5)"),
        ):
            if ci is not None:
                ax.plot(ci, [position + offset] * 2, color=color, linewidth=2.4,
                        alpha=0.55, solid_capstyle="round")
            ax.scatter(value, position + offset, color=color, s=68, zorder=3,
                       label=label if position == 0 else None)
            ax.text(value, position + offset + 0.11, f"{value:+.4f}",
                    ha="center", fontsize=8.8, color=INK)
    ax.axvline(0, color=INK, linewidth=1.1)
    ax.set_yticks(positions)
    ax.set_yticklabels([row[0] for row in rows], fontsize=9.5)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel("효용 격차 (MCTS - NCCN)")
    ax.set_title("헤드라인 수치는 사실상 그대로다", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    ax = axes[1]
    diagnostics = [
        ("첫 결정 일치율 (%)",
         old_robust["mean_modal_first_action_agreement_pct"],
         new_robust["mean_modal_first_action_agreement_pct"]),
        ("예산 2048 일치율 (%)",
         old_budget["agreement_at_max_budget_pct"],
         new_budget["agreement_at_max_budget_pct"]),
        ("2048에서 불안정 지점 (%)",
         old_budget["unstable_node_pct_at_max_budget"],
         new_budget["unstable_node_pct_at_max_budget"]),
    ]
    positions = np.arange(len(diagnostics))
    width = 0.36
    ax.barh(positions + width / 2, [row[1] for row in diagnostics], width,
            color=GRAY, label="편향 환경")
    ax.barh(positions - width / 2, [row[2] for row in diagnostics], width,
            color=GREEN, label="수정 환경")
    for position, (_, old, new) in zip(positions, diagnostics):
        ax.text(old + 1.2, position + width / 2, f"{old:.1f}", va="center",
                fontsize=9.5, color=INK)
        ax.text(new + 1.2, position - width / 2, f"{new:.1f}", va="center",
                fontsize=9.5, color=INK)
    ax.set_yticks(positions)
    ax.set_yticklabels([row[0] for row in diagnostics], fontsize=9.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("%")
    ax.set_title("진단 지표도 결론을 바꾸지 않는다", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    fig.suptitle("v0.5 환경 수정 후 재실행 — 예측대로 결론은 유지됐다",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig30_environment_refresh.png")


def main() -> None:
    fig30_environment_refresh()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
