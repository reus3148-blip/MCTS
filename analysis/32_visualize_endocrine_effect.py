"""Figures for the endocrine effect and the trimming fix (fig36, fig37).

fig36 carries the estimates and the negative control - the negative control is
the point of the report, so it gets its own panel on the relative scale where
the comparison with the randomised benchmark is legible.

fig37 carries the methodological finding: one trim-then-refit pass is not a fixed
point, and for five of the ten decision cells the population it left behind was
not balanced at all. Two pairs of cells sit on top of each other in the left
panel - restricting to mastectomy patients makes surgery type constant, so the
two covariate specs collapse to the same model there.

Every number is read from the committed ``metrics.json`` and ``tables/``.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REPORT_DIR = ROOT / "reports" / "endocrine-effect-v1.3"
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
BALANCE_THRESHOLD = 0.1
#: EBCTCG: roughly a third off breast-cancer mortality in ER-positive disease,
#: nothing in ER-negative disease. Drawn as reference lines, not as our estimate.
EBCTCG_ER_POSITIVE_RR = 0.67
NO_EFFECT_RR = 1.0


def plain_log_ticks(axis) -> None:
    """Plain decimals on a log axis.

    Mathtext exponents render the minus as U+2212, which Malgun Gothic has no
    glyph for - the ticks come out as boxes. Same fix as ``analysis/14``.
    """
    axis.set_major_formatter(FuncFormatter(
        lambda value, _: f"{value:g}" if value >= 0.01 else f"{value:.3f}"))
    axis.set_minor_formatter(NullFormatter())


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig36_effects() -> None:
    arms = pd.read_csv(TABLE_DIR / "arm_results.csv")
    order = list(arms["arm"])

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4),
                             gridspec_kw={"wspace": 0.36, "width_ratios": [1.5, 1]})

    # --- left: forest plot of every arm -------------------------------------
    ax = axes[0]
    positions = np.arange(len(arms))[::-1]
    for position, row in zip(positions, arms.itertuples()):
        color = BLUE if row.prespecified else GRAY
        if row.arm == "hormone_er_negative":
            color = RED
        ax.plot([row.aipw_ci_low, row.aipw_ci_high], [position, position],
                color=color, linewidth=2.6, solid_capstyle="round", alpha=0.55)
        ax.scatter(row.aipw, position, color=color, s=85, zorder=3)
        ax.text(row.aipw_ci_high + 0.012, position,
                f"{row.aipw:+.4f}", va="center", fontsize=9, color=INK)
    ax.axvline(0, color=INK, linewidth=1.1)
    ax.set_yticks(positions)
    ax.set_yticklabels([
        f"{row.label}\nn={row.analysed}/{row.eligible} "
        f"({row.retained_pct:.0f}%) · |SMD| {row.trimmed_worst_abs_smd:.3f}"
        for row in arms.itertuples()], fontsize=8.5)
    ax.set_xlabel("5년 전체사망 위험차 (치료 - 미치료), 가로선은 95% CI")
    ax.set_xlim(-0.30, 0.26)
    handles = [plt.Line2D([0], [0], color=BLUE, linewidth=3),
               plt.Line2D([0], [0], color=RED, linewidth=3),
               plt.Line2D([0], [0], color=GRAY, linewidth=3)]
    ax.legend(handles, ["사전 선언", "음성대조 (효과가 없어야 함)", "사후 추가"],
              frameon=False, fontsize=8.5, loc="lower left")
    ax.set_title("결정별 효과 추정 — 겹침이 확보된 모든 결정", fontsize=12.5, pad=10)

    # --- right: the negative control, on the relative scale -----------------
    ax = axes[1]
    positive = arms[arms.arm == "hormone_er_positive"].iloc[0]
    negative = arms[arms.arm == "hormone_er_negative"].iloc[0]
    bars = ax.bar([0, 1], [positive.risk_ratio, negative.risk_ratio],
                  color=[BLUE, RED], width=0.5)
    for bar, value in zip(bars, (positive.risk_ratio, negative.risk_ratio)):
        ax.text(bar.get_x() + bar.get_width() / 2, value - 0.075,
                f"{value:.3f}", ha="center", fontsize=13, color="white")
    ax.axhline(NO_EFFECT_RR, color=INK, linewidth=1.2)
    ax.text(-0.55, NO_EFFECT_RR + 0.02, "효과 없음", fontsize=8.5, color=INK,
            ha="left")
    ax.axhline(EBCTCG_ER_POSITIVE_RR, color=GREEN, linewidth=1.4, linestyle="--")
    ax.text(-0.55, EBCTCG_ER_POSITIVE_RR + 0.03,
            "무작위배정 근거가 말하는 ER 양성 값 (약 1/3 감소)",
            fontsize=8.5, color=GREEN, ha="left")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["ER 양성\n(약이 듣는 곳)", "ER 음성\n(약이 안 듣는 곳)"],
                       fontsize=9.5)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("5년 사망 위험비 (RR)")
    ax.set_title("음성대조가 실패했다 — 두 값이 거의 같다\n"
                 "약이 들을 수 없는 곳에서도 같은 크기의 '효과'가 나온다",
                 fontsize=12, pad=10, color=RED)

    fig.suptitle("v1.3 — 겹침이 확보된 결정의 효과 추정과 음성대조",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig36_endocrine_effect.png")


def fig37_trimming() -> None:
    convergence = pd.read_csv(TABLE_DIR / "trim_convergence.csv")
    single = convergence[convergence.trimming == "single_pass"].set_index(
        ["covariate_spec", "decision", "population"])
    iterated = convergence[convergence.trimming == "iterated"].set_index(
        ["covariate_spec", "decision", "population"])
    trimming = METRICS["trimming_effect"]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0),
                             gridspec_kw={"wspace": 0.30})

    # --- left: balance before and after, all ten cells ----------------------
    ax = axes[0]
    xs = single.loc[iterated.index, "worst_abs_smd"]
    ys = iterated["worst_abs_smd"]
    failed = xs >= BALANCE_THRESHOLD
    ax.scatter(xs[~failed], ys[~failed], s=80, color=BLUE, zorder=3,
               label="한 번만 돌려도 균형은 맞았던 칸")
    ax.scatter(xs[failed], ys[failed], s=110, color=RED, zorder=3,
               label="한 번만 돌리면 균형이 깨진 칸")
    ax.axhline(BALANCE_THRESHOLD, color=GRAY, linewidth=1.1, linestyle=":")
    ax.axvline(BALANCE_THRESHOLD, color=GRAY, linewidth=1.1, linestyle=":")
    ax.text(0.011, BALANCE_THRESHOLD * 1.12, "균형 기준 |SMD| 0.1",
            fontsize=8.5, color=GRAY)
    limits = (0.008, 3.0)
    ax.plot(limits, limits, color=GRAY, linewidth=1.0, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*limits)
    ax.set_ylim(*limits)
    plain_log_ticks(ax.xaxis)
    plain_log_ticks(ax.yaxis)
    ax.set_xlabel("한 번만 트리밍했을 때의 최악 |SMD|")
    ax.set_ylabel("고정점까지 반복했을 때")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title(
        f"{trimming['decision_cells']}칸 전부 한 번으로는 부족했다 — "
        f"그중 {int(failed.sum())}칸은 균형까지 깨져 있었다",
        fontsize=12, pad=10, color=RED)

    # --- right: what it did to the endocrine arm ----------------------------
    ax = axes[1]
    endocrine = trimming["endocrine_er_positive"]
    labels = ["분석 환자 수", "유효표본", "최악 |SMD| × 1000"]
    before = [endocrine["single_pass"]["n"],
              endocrine["single_pass"]["effective_sample_size"],
              endocrine["single_pass"]["worst_abs_smd"] * 1000]
    after = [endocrine["iterated"]["n"],
             endocrine["iterated"]["effective_sample_size"],
             endocrine["iterated"]["worst_abs_smd"] * 1000]
    positions = np.arange(len(labels))
    width = 0.36
    ax.bar(positions - width / 2, before, width, color=RED, label="한 번만 (v0.6~v0.7 방식)")
    ax.bar(positions + width / 2, after, width, color=GREEN, label="고정점까지 반복")
    for position, (b, a) in enumerate(zip(before, after)):
        ax.text(position - width / 2, b * 1.06, f"{b:,.0f}", ha="center",
                fontsize=9.5, color=INK)
        ax.text(position + width / 2, a * 1.06, f"{a:,.0f}", ha="center",
                fontsize=9.5, color=INK)
    ax.set_yscale("log")
    plain_log_ticks(ax.yaxis)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel("환자 수 · 유효표본 · |SMD|×1000 (로그 축)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    ax.set_title("호르몬치료 (ER 양성)에서 무슨 일이 있었나\n"
                 "환자는 줄었지만 유효표본은 6배가 됐다",
                 fontsize=12, pad=10, color=INK)

    fig.suptitle("v1.3 — 트리밍을 한 번만 하면 그 결과가 스스로의 조건을 만족하지 않는다",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, "fig37_trim_convergence.png")


def main() -> None:
    fig36_effects()
    fig37_trimming()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
