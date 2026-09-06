"""Figure for the reward-model confounding test (fig38).

Three things have to be visible at once for this result to be readable: that the
reward model disagrees in *sign* with our own adjusted analyses, that removing
those coefficients closes a large share of the headline gap, and that it closes
it by lifting NCCN rather than by lowering MCTS - which is the signature of an
undeclared channel that was docking one policy.

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

REPORT_DIR = ROOT / "reports" / "reward-confounding-v1.4"
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

PRETTY_DECISION = {
    "항암 (겹침 구간)": "항암치료",
    "호르몬 (ER 양성)": "호르몬치료\n(ER 양성)",
    "호르몬 (ER 음성)": "호르몬치료\n(ER 음성)",
    "방사선 (전절제 후)": "방사선치료\n(전절제 후)",
}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig38_reward_confounding() -> None:
    verdict = METRICS["verdict"]
    arms = METRICS["arms"]
    implied = pd.DataFrame(METRICS["reward_model_vs_causal"])
    # The ER-negative arm is v1.3's failed negative control, not an estimate of an
    # effect - putting it on a "what our analysis says" axis would misrepresent it,
    # and its magnitude would swamp the scale. Reported in the text instead.
    implied = implied[implied["population"] != "호르몬 (ER 음성)"].reset_index(drop=True)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.2),
                             gridspec_kw={"wspace": 0.34,
                                          "width_ratios": [1.15, 1, 1]})

    # --- A: the model disagrees in sign with our own analyses ---------------
    ax = axes[0]
    positions = np.arange(len(implied))
    width = 0.36
    ax.bar(positions - width / 2, implied["reward_model_risk_difference"], width,
           color=RED, label="보상모형이 믿는 값")
    ax.bar(positions + width / 2, implied["causal_risk_difference"], width,
           color=BLUE, label="우리 보정 분석 (v1.3)")
    ax.axhline(0, color=INK, linewidth=1.2)
    ax.set_xticks(positions)
    ax.set_xticklabels([PRETTY_DECISION.get(name, name)
                        for name in implied["population"]], fontsize=8.5)
    ax.set_ylabel("5년 사망 위험차 (치료 - 미치료)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    disagreements = int(implied["sign_disagrees"].sum())
    ax.set_title(f"A. 보상모형은 치료를 해롭다고 믿는다\n"
                 f"{len(implied)}개 결정 중 {disagreements}개에서 부호가 반대다",
                 fontsize=12, pad=10, color=RED)
    ax.text(0.02, 0.97, "위쪽 = 해롭다", transform=ax.transAxes,
            fontsize=8.5, color=GRAY, va="top")

    # --- B: what neutralising does to the headline --------------------------
    ax = axes[1]
    gaps = [verdict["gap_as_fitted"], verdict["gap_treatment_neutral"]]
    errors = [1.96 * arms["as_fitted"]["standard_error"],
              1.96 * arms["treatment_neutral"]["standard_error"]]
    positions = np.arange(2)
    ax.bar(positions, gaps, yerr=errors, color=[RED, GREEN], width=0.5,
           error_kw={"ecolor": INK, "capsize": 4, "linewidth": 1.1})
    for x, gap, error in zip(positions, gaps, errors):
        ax.text(x, gap + error + 0.0012, f"{gap:+.4f}", ha="center",
                fontsize=11.5, color=INK)
    ax.axhline(0, color=INK, linewidth=1.1)
    ax.set_xticks(positions)
    ax.set_xticklabels(["보상모형 그대로\n(v0.2~v1.2)", "치료 계수\n중립화"],
                       fontsize=9.5)
    ax.set_ylabel("MCTS - NCCN 효용 격차")
    ax.set_ylim(0, max(g + e for g, e in zip(gaps, errors)) + 0.006)
    share = verdict["share_of_gap_from_reward_model"]
    ax.set_title(f"B. 격차의 {share:.0%}가 사라진다\n"
                 f"짝지은 차이 {verdict['paired_difference']:+.4f} "
                 f"(z = {verdict['paired_z']:.1f})",
                 fontsize=12, pad=10, color=RED)

    # --- C: it closes by lifting NCCN, not by lowering MCTS -----------------
    ax = axes[2]
    labels = ["MCTS", "NCCN"]
    before = [arms["as_fitted"]["mcts_utility"], arms["as_fitted"]["nccn_utility"]]
    after = [arms["treatment_neutral"]["mcts_utility"],
             arms["treatment_neutral"]["nccn_utility"]]
    positions = np.arange(2)
    width = 0.34
    ax.bar(positions - width / 2, before, width, color=RED, label="보상모형 그대로")
    ax.bar(positions + width / 2, after, width, color=GREEN, label="치료 계수 중립화")
    for position, (b, a) in enumerate(zip(before, after)):
        ax.annotate("", xy=(position + width / 2, a),
                    xytext=(position - width / 2, b),
                    arrowprops=dict(arrowstyle="->", color=INK, linewidth=1.4))
        ax.text(position, max(b, a) + 0.004, f"{a - b:+.4f}", ha="center",
                fontsize=10.5, color=GREEN if a > b else RED)
    low = min(before + after)
    ax.set_ylim(low - 0.02, max(before + after) + 0.018)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_ylabel("평균 효용")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title("C. MCTS가 나빠진 게 아니라\nNCCN이 벌점을 그만 받은 것이다",
                 fontsize=12, pad=10, color=INK)

    fig.suptitle("v1.4 — MCTS 우세의 얼마가 보상모형의 교란에서 왔나 "
                 "(오차막대는 95% CI)", fontsize=13.5, y=1.03, color=INK)
    save(fig, "fig38_reward_confounding.png")


def fig39_headline_corrections() -> None:
    """The two independent corrections on the same headline, and their product."""
    verdict = METRICS["verdict"]
    rows = [
        ("지금까지 인용해 온 값\n(아형 균등 · 보상모형 그대로)",
         verdict["gap_as_fitted"], GRAY),
        ("치료 계수 중립화만\n(v1.4)",
         verdict["gap_treatment_neutral"], BLUE),
        ("아형 표준화만\n(v1.2)",
         verdict["standardised_gap_as_fitted"], AMBER),
        ("둘 다 적용\n(현재 최선의 값)",
         verdict["both_corrections_applied"], GREEN),
    ]
    fig, ax = plt.subplots(figsize=(10.0, 5.0))
    positions = np.arange(len(rows))
    bars = ax.bar(positions, [row[1] for row in rows],
                  color=[row[2] for row in rows], width=0.58)
    for bar, (_, value, _) in zip(bars, rows):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.0009,
                f"{value:+.4f}", ha="center", fontsize=12, color=INK)
    ax.axhline(0, color=INK, linewidth=1.2)
    ax.set_xticks(positions)
    ax.set_xticklabels([row[0] for row in rows], fontsize=9)
    ax.set_ylabel("MCTS - NCCN 효용 격차")
    ax.set_ylim(0, max(row[1] for row in rows) * 1.20)
    shrink = 1.0 - rows[3][1] / rows[0][1]
    ax.set_title("두 정정은 서로 다른 것을 고치므로 겹쳐서 적용된다 — "
                 f"합쳐서 {shrink:.0%} 감소\n"
                 "하나는 표집 설계, 하나는 보상모형의 미선언 채널",
                 fontsize=12.5, pad=12, color=INK)
    save(fig, "fig39_headline_corrections.png")


def main() -> None:
    fig38_reward_confounding()
    fig39_headline_corrections()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
