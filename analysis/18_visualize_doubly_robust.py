"""Figures for the doubly robust emulation and identifiability map (fig27 - fig28).

Two claims to make visible: three estimators that lean on different models land in
the same place, and the apparent answerability of a decision depends on which
confounders you bothered to measure.
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

REPORT_DIR = ROOT / "reports" / "doubly-robust-v0.7"
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


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig27_estimator_agreement() -> None:
    """Three estimators leaning on different models, and the CI that dwarfs them."""
    labels = {
        "ipw_km": "IPW 가중 KM\n(치료모형 의존)",
        "g_computation": "g-계산\n(결과모형 의존)",
        "aipw": "AIPW\n(둘 중 하나만 맞으면 됨)",
    }
    colors = {"ipw_km": BLUE, "g_computation": AMBER, "aipw": GREEN}
    names = list(labels)
    values = [METRICS["estimators"][name]["risk_difference"] for name in names]
    low, high = METRICS["aipw_ci95"]["risk_difference"]
    spread = METRICS["estimator_spread_risk_difference"]

    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    ax.axvspan(low, high, color=GREEN, alpha=0.10)
    ax.axvline(0, color=INK, linewidth=1.2)
    for position, name in enumerate(names):
        ax.scatter(values[position], position, color=colors[name], s=110, zorder=3)
        ax.text(values[position], position + 0.24,
                f"{values[position]:+.4f}", ha="center", fontsize=10, color=INK)
    ax.errorbar([METRICS["estimators"]["aipw"]["risk_difference"]], [2],
                xerr=[[METRICS["estimators"]["aipw"]["risk_difference"] - low],
                      [high - METRICS["estimators"]["aipw"]["risk_difference"]]],
                fmt="none", ecolor=GREEN, capsize=5, linewidth=1.8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([labels[name] for name in names], fontsize=10)
    ax.set_ylim(-0.6, 2.7)
    ax.set_xlabel("5년 사망 위험 차이 (항암 - 미시행)")
    ax.set_title("서로 다른 모형에 기대는 세 추정량이 같은 곳에 모인다\n"
                 f"추정량 간 최대 차이 {spread:.4f} vs 95% CI 폭 {high - low:.4f}",
                 fontsize=12.5, pad=12)
    ax.text(low, -0.45, f"AIPW 95% CI [{low:+.3f}, {high:+.3f}]",
            fontsize=9, color=GREEN)
    save(fig, "fig27_estimator_agreement.png")


def fig28_identifiability_map() -> None:
    """How answerable each decision looks, and how that survives one more confounder."""
    table = pd.read_csv(TABLE_DIR / "identifiability_map.csv")
    labels = {"chemo": "보조 항암", "hormone": "호르몬치료", "radio": "방사선치료"}
    order = ["chemo", "hormone", "radio"]
    positions = np.arange(len(order))
    width = 0.36

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

    ax = axes[0]
    for offset, (spec, color, label) in enumerate((
        ("baseline", GRAY, "기본 교란요인"),
        ("baseline + surgery", BLUE, "+ 수술 유형"),
    )):
        values = [
            float(table[(table.decision == d) & (table.covariate_spec == spec)]
                  ["retained_pct"].iloc[0]) for d in order
        ]
        bars = ax.bar(positions + (offset - 0.5) * width, values, width,
                      color=color, label=label)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.5,
                    f"{value:.0f}%", ha="center", fontsize=9.5, color=INK)
    ax.set_xticks(positions)
    ax.set_xticklabels([labels[d] for d in order])
    ax.set_ylim(0, 112)
    ax.set_ylabel("겹침 구간에 남는 환자 비율 (%)")
    ax.set_title("결정마다 답할 수 있는 정도가 다르다", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for offset, (spec, color, label) in enumerate((
        ("baseline", GRAY, "기본 교란요인"),
        ("baseline + surgery", BLUE, "+ 수술 유형"),
    )):
        values = [
            float(table[(table.decision == d) & (table.covariate_spec == spec)]
                  ["trimmed_worst_abs_smd"].iloc[0]) for d in order
        ]
        bars = ax.bar(positions + (offset - 0.5) * width, values, width,
                      color=color, label=label)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.006,
                    f"{value:.3f}", ha="center", fontsize=9.5, color=INK)
    ax.axhline(0.1, color=RED, linestyle="--", linewidth=1.3)
    ax.text(2.35, 0.107, "균형 기준 0.1", color=RED, fontsize=9, ha="right")
    ax.set_xticks(positions)
    ax.set_xticklabels([labels[d] for d in order])
    ax.set_ylabel("트리밍 후 최악 |SMD|")
    ax.set_title("방사선은 교란요인 하나를 더하자 균형에 실패한다",
                 fontsize=12, pad=10, color=RED)
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    fig.suptitle("넓은 겹침은 좋은 소식이 아니라 측정하지 않았다는 뜻일 수 있다",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig28_identifiability_map.png")


def main() -> None:
    fig27_estimator_agreement()
    fig28_identifiability_map()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
