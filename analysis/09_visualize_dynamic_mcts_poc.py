"""Create figures for the stochastic dynamic MCTS PoC v0.2."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = ROOT / "reports" / "dynamic-mcts-poc-v0.2"
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
POLICY_COLORS = {"NCCN": BLUE, "MCTS": GREEN}


def save_and_publish(fig: plt.Figure, filename: str) -> None:
    report_path = FIGURE_DIR / filename
    fig.savefig(report_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(report_path, PUBLIC_DIR / filename)
    print(report_path)


def add_box(ax, x, y, width, height, text, color, fontsize=10):
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.015",
        linewidth=1.2,
        edgecolor=color,
        facecolor=color + "18",
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color=INK)


def add_arrow(ax, start, end, text=None, color=GRAY):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={"arrowstyle": "->", "color": color, "lw": 1.4},
    )
    if text:
        ax.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + 0.025,
            text,
            ha="center",
            va="bottom",
            fontsize=8,
            color=GRAY,
        )


def plot_environment_flow() -> None:
    fig, ax = plt.subplots(figsize=(12, 6.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("동적 환경 v0.2: 치료 결과가 다음 선택을 바꾼다", fontsize=15, pad=18)

    add_box(ax, 0.08, 0.72, 0.13, 0.12, "진단 시 상태\n병기·아형·종양크기", BLUE)
    add_box(ax, 0.25, 0.72, 0.12, 0.12, "치료 순서\n선택", BLUE)
    add_box(ax, 0.43, 0.86, 0.14, 0.11, "수술 먼저", AMBER)
    add_box(ax, 0.43, 0.56, 0.14, 0.11, "선행 항암\n표준 / 강화", GREEN)
    add_box(ax, 0.61, 0.56, 0.15, 0.14, "반응 발생\n큰 반응 / 부분 / 없음\n+ 독성", GREEN, fontsize=9)
    add_box(ax, 0.61, 0.86, 0.14, 0.11, "수술\nBCS / MAST", AMBER)
    add_box(ax, 0.78, 0.72, 0.13, 0.15, "보조치료\n항암·호르몬\n방사선", BLUE, fontsize=9)
    add_box(ax, 0.92, 0.72, 0.12, 0.15, "5년 추적\n무사건 / 재발\n사망", RED, fontsize=9)

    add_arrow(ax, (0.145, 0.72), (0.19, 0.72))
    add_arrow(ax, (0.31, 0.75), (0.36, 0.83), "수술 우선")
    add_arrow(ax, (0.31, 0.69), (0.36, 0.59), "약물 우선")
    add_arrow(ax, (0.50, 0.56), (0.535, 0.56))
    add_arrow(ax, (0.61, 0.63), (0.61, 0.80), "종양크기 갱신")
    add_arrow(ax, (0.50, 0.86), (0.54, 0.86))
    add_arrow(ax, (0.68, 0.84), (0.715, 0.76))
    add_arrow(ax, (0.845, 0.72), (0.855, 0.72))

    ax.text(
        0.5,
        0.30,
        "핵심 변화",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        0.5,
        0.21,
        "v0.1은 처음 환자정보가 끝까지 고정됐지만, v0.2는 치료 반응·종양크기·독성·재발 상태가 매 단계 바뀐다.",
        ha="center",
        va="center",
        fontsize=10,
        color=GRAY,
    )
    ax.text(
        0.5,
        0.13,
        "실선 상자의 OS/RFS 바닥 위험은 METABRIC에서 학습하고, 반응·독성 차이는 공개된 합성 가정 파일로 분리한다.",
        ha="center",
        va="center",
        fontsize=9,
        color=RED,
    )
    save_and_publish(fig, "fig14_dynamic_environment_flow.png")


def plot_policy_outcomes(summary: pd.DataFrame) -> None:
    indexed = summary.set_index("policy")
    metrics = [
        ("mean_utility", "평균 utility", 1.0, "{:.3f}"),
        ("survived_5y_pct", "5년 생존", 100.0, "{:.1f}%"),
        ("recurred_by_5y_pct", "5년 내 재발", 12.0, "{:.1f}%"),
        ("mean_toxicity_count", "평균 독성 사건", 0.5, "{:.3f}"),
    ]
    policies = ["NCCN", "MCTS"]
    fig, axes = plt.subplots(1, 4, figsize=(12, 4.2))
    for ax, (column, title, ceiling, formatter) in zip(axes, metrics, strict=True):
        values = [indexed.at[policy, column] for policy in policies]
        bars = ax.bar(
            policies,
            values,
            color=[POLICY_COLORS[policy] for policy in policies],
            width=0.62,
        )
        for bar, value in zip(bars, values, strict=True):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + ceiling * 0.025,
                formatter.format(value),
                ha="center",
                va="bottom",
                fontsize=9,
            )
        ax.set_ylim(0, ceiling)
        ax.set_title(title, fontsize=11)
        ax.grid(axis="y")
        ax.grid(axis="x", visible=False)
    fig.suptitle("동적 환경 4,000회씩 반복한 정책 결과", fontsize=14)
    fig.text(
        0.5,
        -0.02,
        "METABRIC 기반 위험 + 합성 전이 가정의 시뮬레이션 결과이며 임상효과 추정이 아님",
        ha="center",
        color=RED,
        fontsize=9,
    )
    fig.tight_layout()
    save_and_publish(fig, "fig15_dynamic_policy_outcomes.png")


def plot_action_choices(summary: pd.DataFrame) -> None:
    indexed = summary.set_index("policy")
    measures = [
        ("neoadjuvant_rate_pct", "선행치료"),
        ("bcs_rate_pct", "BCS"),
        ("intensified_chemo_rate_pct", "강화 항암"),
        ("extended_endocrine_rate_pct", "연장 호르몬"),
        ("regional_radiation_rate_pct", "광범위 방사선"),
    ]
    policies = ["NCCN", "MCTS"]
    x = np.arange(len(measures))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for index, policy in enumerate(policies):
        values = [indexed.at[policy, column] for column, _ in measures]
        offset = (index - 0.5) * width
        bars = ax.bar(
            x + offset,
            values,
            width,
            color=POLICY_COLORS[policy],
            label=policy,
        )
        for bar, value in zip(bars, values, strict=True):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1.2,
                f"{value:.0f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
    ax.set_xticks(x, [label for _, label in measures])
    ax.set_ylim(0, 105)
    ax.set_ylabel("선택 비율 (%)")
    ax.set_title("선택지가 늘어나자 MCTS가 실제로 사용한 행동")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    save_and_publish(fig, "fig16_dynamic_action_choices.png")


def plot_search_stability(stability: pd.DataFrame) -> None:
    x = stability["simulations"]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.plot(
        x,
        stability["reference_action_match_pct"],
        color=BLUE,
        marker="o",
        linewidth=2,
        label="1,024회 기준과 같은 첫 행동",
    )
    ax.plot(
        x,
        stability["match_or_near_tie_pct"],
        color=GREEN,
        marker="s",
        linewidth=2,
        label="같거나 기준 기대값 차이 0.01 이하",
    )
    ax.set_xscale("log", base=2)
    ax.set_xticks(x, labels=x.astype(int))
    ax.set_ylim(0, 105)
    ax.set_xlabel("결정당 시뮬레이션 수")
    ax.set_ylabel("환자 비율 (%)")
    ax.set_title("확률환경에서의 첫 행동 안정성")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_and_publish(fig, "fig17_dynamic_search_stability.png")


def main() -> None:
    summary = pd.read_csv(TABLE_DIR / "policy_summary.csv")
    stability = pd.read_csv(TABLE_DIR / "search_stability.csv")
    plot_environment_flow()
    plot_policy_outcomes(summary)
    plot_action_choices(summary)
    plot_search_stability(stability)


if __name__ == "__main__":
    main()

