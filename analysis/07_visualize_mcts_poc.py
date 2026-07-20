"""Create publication-style figures for the MCTS v0.1 PoC report."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = ROOT / "reports" / "mcts-poc-v1"
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
POLICY_COLORS = {"Actual": "#78716c", "NCCN": BLUE, "MCTS": GREEN}
DECISION_LABELS = {
    "surgery": "수술",
    "chemo": "항암",
    "hormone": "호르몬",
    "radio": "방사선",
    "all_four": "4개 모두",
}


def save_and_publish(fig: plt.Figure, filename: str) -> Path:
    report_path = FIGURE_DIR / filename
    fig.savefig(report_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(report_path, PUBLIC_DIR / filename)
    print(report_path)
    return report_path


def plot_search_convergence(convergence: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    x = convergence["simulations_per_step"]
    axes[0].plot(
        x,
        convergence["exact_plan_match_pct"],
        color=BLUE,
        marker="o",
        linewidth=2,
    )
    axes[0].set_xscale("log", base=2)
    axes[0].set_xticks(x, labels=x.astype(int))
    axes[0].set_ylim(0, 103)
    axes[0].set_xlabel("결정 단계당 시뮬레이션 수")
    axes[0].set_ylabel("완전탐색 최적 경로 일치율 (%)")
    axes[0].set_title("MCTS 수렴 정확도")

    axes[1].plot(
        x,
        convergence["mean_regret"],
        color=RED,
        marker="o",
        linewidth=2,
    )
    axes[1].set_xscale("log", base=2)
    axes[1].set_xticks(x, labels=x.astype(int))
    axes[1].set_xlabel("결정 단계당 시뮬레이션 수")
    axes[1].set_ylabel("평균 regret (5년 OS 확률 단위)")
    axes[1].set_title("최적 보상과의 차이")
    fig.suptitle("UCT-MCTS 내부 검증: 16개 경로 완전탐색과 비교", fontsize=13)
    fig.tight_layout()
    save_and_publish(fig, "fig10_mcts_search_convergence.png")


def plot_policy_agreement(agreement: pd.DataFrame) -> None:
    view = agreement[agreement["comparison"].isin([
        "MCTS vs NCCN",
        "MCTS vs Actual",
    ])].copy()
    decisions = ["surgery", "chemo", "hormone", "radio", "all_four"]
    comparisons = ["MCTS vs NCCN", "MCTS vs Actual"]
    x = np.arange(len(decisions))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for index, comparison in enumerate(comparisons):
        values = (
            view[view["comparison"].eq(comparison)]
            .set_index("decision")
            .loc[decisions, "agreement_pct"]
            .to_numpy()
        )
        offset = (index - 0.5) * width
        bars = ax.bar(
            x + offset,
            values,
            width,
            color=BLUE if index == 0 else AMBER,
            label="NCCN과 일치" if index == 0 else "실제 치료와 일치",
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
    ax.set_xticks(x, [DECISION_LABELS[d] for d in decisions])
    ax.set_ylim(0, 105)
    ax.set_ylabel("일치율 (%)")
    ax.set_title("MCTS 정책은 무엇과 얼마나 같은가")
    ax.legend(frameon=False, ncol=2, loc="upper center")
    fig.tight_layout()
    save_and_publish(fig, "fig11_mcts_policy_agreement.png")


def plot_treatment_rates(summary: pd.DataFrame) -> None:
    measures = [
        ("bcs_rate_pct", "BCS"),
        ("chemo_rate_pct", "항암"),
        ("hormone_rate_pct", "호르몬"),
        ("radio_rate_pct", "방사선"),
    ]
    policies = ["Actual", "NCCN", "MCTS"]
    x = np.arange(len(measures))
    width = 0.24
    fig, ax = plt.subplots(figsize=(9, 4.8))
    indexed = summary.set_index("policy")
    for index, policy in enumerate(policies):
        values = [indexed.at[policy, column] for column, _ in measures]
        offset = (index - 1) * width
        bars = ax.bar(
            x + offset,
            values,
            width,
            color=POLICY_COLORS[policy],
            label={"Actual": "실제", "NCCN": "NCCN", "MCTS": "MCTS"}[policy],
        )
        for bar, value in zip(bars, values, strict=True):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1.0,
                f"{value:.0f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
    ax.set_xticks(x, [label for _, label in measures])
    ax.set_ylim(0, 105)
    ax.set_ylabel("선택 비율 (%)")
    ax.set_title("정책별 치료 선택 비율")
    ax.legend(frameon=False, ncol=3, loc="upper center")
    fig.tight_layout()
    save_and_publish(fig, "fig12_mcts_treatment_rates.png")


def plot_predicted_survival(decisions: pd.DataFrame) -> None:
    view = decisions[decisions["nccn_complete"].eq(1)]
    values = [
        view["predicted_5y_os_actual"].to_numpy() * 100,
        view["predicted_5y_os_nccn"].to_numpy() * 100,
        view["predicted_5y_os_mcts"].to_numpy() * 100,
    ]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    box = ax.boxplot(
        values,
        tick_labels=["실제 치료", "NCCN", "MCTS"],
        patch_artist=True,
        showfliers=False,
        medianprops={"color": INK, "linewidth": 1.5},
    )
    for patch, color in zip(
        box["boxes"],
        [POLICY_COLORS["Actual"], POLICY_COLORS["NCCN"], POLICY_COLORS["MCTS"]],
        strict=True,
    ):
        patch.set_facecolor(color)
        patch.set_alpha(0.72)
    means = [float(np.mean(group)) for group in values]
    ax.scatter([1, 2, 3], means, color=INK, marker="D", s=28, zorder=3, label="평균")
    ax.set_ylim(0, 105)
    ax.set_ylabel("예측 5년 전체생존 확률 (%)")
    ax.set_title("보상모형이 계산한 정책별 5년 OS")
    ax.text(
        0.5,
        -0.18,
        "관찰자료 기반 모형 예측이며 치료의 인과효과 또는 임상 우월성 추정이 아님",
        transform=ax.transAxes,
        ha="center",
        color=RED,
        fontsize=9,
    )
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_and_publish(fig, "fig13_mcts_predicted_survival.png")


def main() -> None:
    convergence = pd.read_csv(TABLE_DIR / "search_convergence.csv")
    agreement = pd.read_csv(TABLE_DIR / "policy_agreement.csv")
    summary = pd.read_csv(TABLE_DIR / "policy_summary.csv")
    decisions = pd.read_csv(TABLE_DIR / "patient_policy_decisions.csv")
    plot_search_convergence(convergence)
    plot_policy_agreement(agreement)
    plot_treatment_rates(summary)
    plot_predicted_survival(decisions)


if __name__ == "__main__":
    main()

