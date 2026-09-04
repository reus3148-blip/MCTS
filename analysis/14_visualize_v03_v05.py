"""Figures for the v0.3 - v0.5 findings (fig18 - fig23).

v0.3, v0.4 and v0.5 carry the findings that matter most for a talk - that our own
result was overstated, that 256 simulations could not resolve an action ordering,
and that our environment was handing one policy an undeclared advantage - and
none of them had a single figure. Tables are enough for a reader who already
knows the argument; they are not enough for an audience.

Every panel here is drawn from the committed tables of those reports, so a figure
can never drift from the numbers it illustrates. Style matches
``analysis/09_visualize_dynamic_mcts_poc.py`` so the deck reads as one system.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import FuncFormatter, NullFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ROBUSTNESS = ROOT / "reports" / "robustness-v0.3"
SENSITIVITY = ROOT / "reports" / "sensitivity-v0.3"
BUDGET = ROOT / "reports" / "budget-scaling-v0.4"
ENVFIX = ROOT / "reports" / "environment-fix-v0.5"
PUBLIC_DIR = ROOT / "public" / "figures"
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
LIGHT = "#d6d3d1"


def save(fig: plt.Figure, label_dir: Path, filename: str) -> None:
    figure_dir = label_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    path = figure_dir / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


# --------------------------------------------------------------------------
# fig18 - the v0.4 headline: search was fine, we were under-running it
# --------------------------------------------------------------------------
def fig18_budget_scaling() -> None:
    table = pd.read_csv(BUDGET / "tables" / "convergence_by_budget.csv")
    metrics = json.loads((BUDGET / "metrics.json").read_text(encoding="utf-8"))
    slope = metrics["value_noise_log_log_slope"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    ax = axes[0]
    ax.plot(table["budget"], table["mean_agreement_pct"],
            marker="o", color=GREEN, linewidth=2.2, label="평균 일치율")
    ax.plot(table["budget"], table["unstable_node_pct"],
            marker="s", color=AMBER, linewidth=1.8, linestyle="--",
            label="불안정 지점 비율")
    ax.axvline(256, color=RED, linewidth=1.2, linestyle=":")
    ax.text(268, 86, "v0.2~v0.4가\n쓰던 예산", color=RED, fontsize=9, va="top")
    ax.set_xscale("log", base=2)
    ax.set_xticks(table["budget"])
    ax.set_xticklabels(table["budget"])
    ax.set_xlabel("탐색 예산 (시뮬레이션 횟수)")
    ax.set_ylabel("%")
    ax.set_title("예산을 키우면 결정이 안정된다", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.plot(table["budget"], table["mean_value_noise_sd"],
            marker="o", color=BLUE, linewidth=2.2, label="추정 잡음 SD")
    ax.plot(table["budget"], table["mean_abs_value_gap"],
            marker="^", color=GRAY, linewidth=1.8, linestyle="--",
            label="1·2위 행동 가치 차이")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(table["budget"])
    ax.set_xticklabels(table["budget"])
    # Malgun Gothic has no minus glyph, so mathtext exponents (10^-2) break.
    # Plain decimals keep the log axis readable.
    ax.set_yticks([0.01, 0.02, 0.03, 0.05, 0.07])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.02f}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("탐색 예산 (시뮬레이션 횟수)")
    ax.set_ylabel("정규화 효용 (로그 축)")
    ax.set_title(f"잡음은 1/√N로 줄고 (기울기 {slope:+.3f}), 격차는 그대로",
                 fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=9)

    fig.suptitle("v0.4 — 흔들림의 원인은 탐색 부족이었다 (일부는)",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, BUDGET, "fig18_budget_scaling.png")


# --------------------------------------------------------------------------
# fig19 - where the irreducible ties live
# --------------------------------------------------------------------------
def fig19_phase_agreement() -> None:
    detail = pd.read_csv(BUDGET / "tables" / "node_budget_detail.csv")
    top = detail[detail["budget"] == detail["budget"].max()]
    grouped = (top.groupby("phase")["agreement_pct"]
               .agg(["mean", "min", "count"]).sort_values("mean"))
    labels = {
        "chemo": "항암", "surgery": "수술", "timing": "선행/수술우선",
        "radiation": "방사선", "endocrine": "호르몬치료",
    }
    names = [f"{labels.get(p, p)}\n(n={int(c)})"
             for p, c in zip(grouped.index, grouped["count"])]
    colors = [RED if value < 70 else (AMBER if value < 90 else GREEN)
              for value in grouped["mean"]]

    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    bars = ax.barh(names, grouped["mean"], color=colors, height=0.6)
    ax.scatter(grouped["min"], names, color=INK, s=26, zorder=3, label="최소")
    for bar, value in zip(bars, grouped["mean"]):
        ax.text(value + 1.5, bar.get_y() + bar.get_height() / 2,
                f"{value:.1f}%", va="center", fontsize=10, color=INK)
    ax.axvline(75, color=GRAY, linestyle=":", linewidth=1.2)
    ax.text(76, -0.45, "불안정 기준 75%", color=GRAY, fontsize=9)
    ax.set_xlim(0, 108)
    ax.set_xlabel("예산 2048에서의 결정 일치율 (%)")
    ax.set_title("계산을 아무리 해도 안 갈리는 결정이 있다\n"
                 "— 호르몬치료·방사선은 우리 가중치가 승패를 정한다",
                 fontsize=12.5, pad=12)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    save(fig, BUDGET, "fig19_phase_agreement.png")


# --------------------------------------------------------------------------
# fig20 - the v0.3 sensitivity tornado
# --------------------------------------------------------------------------
def fig20_sensitivity_tornado() -> None:
    influence = pd.read_csv(SENSITIVITY / "tables" / "parameter_influence.csv")
    influence = influence.sort_values("abs_gap_delta").tail(8)
    pretty = {
        "reward.recurrence_free_year": "무재발 1년의 가치",
        "reward.acute_toxicity_penalty": "급성 독성 페널티",
        "toxicity.chemo.intensified": "강화항암 독성 확률",
        "response.intensified.major": "강화항암 major 반응확률",
        "hazard.chemo.intensified.death": "강화항암 사망 HR",
        "hazard.response.major.death": "major 반응 사망 HR",
    }
    judgement = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}
    names = [pretty.get(p, p) for p in influence["parameter"]]
    colors = [RED if p in judgement else BLUE for p in influence["parameter"]]

    fig, ax = plt.subplots(figsize=(9.8, 4.6))
    ax.barh(names, influence["abs_gap_delta"], color=colors, height=0.6)
    for index, value in enumerate(influence["abs_gap_delta"]):
        ax.text(value + 0.0004, index, f"{value:.4f}",
                va="center", fontsize=9.5, color=INK)
    ax.set_xlabel("효용 격차에 준 최대 영향 |Δ|")
    ax.set_xlim(0, influence["abs_gap_delta"].max() * 1.25)
    ax.set_title("v0.3 — 결론을 가장 크게 흔든 건 임상 효과가 아니라 가치판단",
                 fontsize=12.5, pad=12)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=RED),
        plt.Rectangle((0, 0), 1, 1, color=BLUE),
    ]
    ax.legend(handles, ["가치판단 (우리가 정한 가중치)", "데이터로 추정 가능"],
              frameon=False, fontsize=9, loc="lower right")
    save(fig, SENSITIVITY, "fig20_sensitivity_tornado.png")


# --------------------------------------------------------------------------
# fig21 - v0.3 robustness: the advantage is small and seed-dependent
# --------------------------------------------------------------------------
def fig21_seed_robustness() -> None:
    per_seed = pd.read_csv(ROBUSTNESS / "tables" / "per_seed_summary.csv")
    stability = pd.read_csv(ROBUSTNESS / "tables" / "first_action_stability.csv")
    ci = pd.read_csv(ROBUSTNESS / "tables" / "metric_confidence_intervals.csv")
    gap = ci[ci["metric"] == "utility_gap"].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    ax = axes[0]
    order = per_seed.sort_values("utility_gap")
    colors = [RED if value < 0 else GREEN for value in order["utility_gap"]]
    ax.bar(range(len(order)), order["utility_gap"], color=colors, width=0.7)
    ax.axhline(0, color=INK, linewidth=1.0)
    ax.axhline(gap["mean"], color=BLUE, linestyle="--", linewidth=1.4)
    ax.fill_between([-0.6, len(order) - 0.4], gap["ci95_low"], gap["ci95_high"],
                    color=BLUE, alpha=0.12)
    ax.text(0.2, gap["ci95_high"] + 0.001,
            f"평균 {gap['mean']:+.3f}  95% CI [{gap['ci95_low']:+.3f}, "
            f"{gap['ci95_high']:+.3f}]", color=BLUE, fontsize=9)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_xlabel("시드 (효용 격차 오름차순)")
    ax.set_ylabel("효용 격차 (MCTS - NCCN)")
    ax.set_title("20개 시드 중 2개는 MCTS가 졌다", fontsize=12, pad=10)
    ax.set_xticks([])

    ax = axes[1]
    order = stability.sort_values("modal_agreement_pct")
    colors = [RED if value < 75 else GREEN
              for value in order["modal_agreement_pct"]]
    ax.barh(range(len(order)), order["modal_agreement_pct"],
            color=colors, height=0.65)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order["subtype"], fontsize=8.5)
    ax.axvline(order["modal_agreement_pct"].mean(), color=INK,
               linestyle="--", linewidth=1.3)
    ax.text(order["modal_agreement_pct"].mean() + 1.5, 0.2,
            f"평균 {order['modal_agreement_pct'].mean():.1f}%",
            fontsize=9, color=INK)
    ax.set_xlim(0, 108)
    ax.set_xlabel("첫 결정의 시드 간 일치율 (%)")
    ax.set_title("환자 절반에서 첫 결정이 난수에 따라 갈렸다",
                 fontsize=12, pad=10)

    fig.suptitle("v0.3 — 우리 결과를 우리가 공격해 보기",
                 fontsize=13.5, y=1.02, color=INK)
    save(fig, ROBUSTNESS, "fig21_seed_robustness.png")


# --------------------------------------------------------------------------
# fig22 - the v0.5 defect, drawn
# --------------------------------------------------------------------------
def fig22_response_channel_bias() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    arms = ["수술우선", "선행치료\n+ 표준항암", "선행치료\n+ 강화항암"]

    ax = axes[0]
    before = [1.000, 0.935, 0.914]
    bars = ax.bar(arms, before, color=[GRAY, RED, RED], width=0.55)
    for bar, value in zip(bars, before):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.004,
                f"{value:.3f}", ha="center", fontsize=10.5, color=INK)
    ax.axhline(1.0, color=INK, linewidth=1.2, linestyle="--")
    ax.set_ylim(0.87, 1.03)
    ax.set_ylabel("E[재발 위험 배수]")
    ax.set_title("수정 전 — 선행치료만 고르면\n재발 위험이 6.5~8.6% 깎였다",
                 fontsize=12, pad=10, color=RED)
    arrow = FancyArrowPatch((1, 0.935), (1, 0.998), arrowstyle="<->",
                            color=RED, linewidth=1.4, mutation_scale=12)
    ax.add_patch(arrow)
    ax.text(1.12, 0.963, "선언된 적 없는\n이득", color=RED, fontsize=9)

    ax = axes[1]
    after = [1.000, 1.000, 1.000]
    bars = ax.bar(arms, after, color=[GRAY, GREEN, GREEN], width=0.55)
    for bar, value in zip(bars, after):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.004,
                f"{value:.3f}", ha="center", fontsize=10.5, color=INK)
    ax.axhline(1.0, color=INK, linewidth=1.2, linestyle="--")
    ax.set_ylim(0.87, 1.03)
    ax.set_title("수정 후 — 채널 평균을 1.0으로 중립화\n"
                 "(major>none 상대 순서는 유지)",
                 fontsize=12, pad=10, color=GREEN)

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.suptitle("v0.5 — 우리 환경이 한쪽 정책에만 점수를 주고 있었다",
                 fontsize=13.5, y=1.0, color=INK)
    save(fig, ENVFIX, "fig22_response_channel_bias.png")


# --------------------------------------------------------------------------
# fig23 - reconciling every headline number we have published
# --------------------------------------------------------------------------
def fig23_results_reconciliation() -> None:
    """Every headline gap we have published, read from its own metrics.json.

    The number moved five times and each move had a different cause. Reading the
    values from the committed metrics rather than typing them in means this
    figure cannot drift away from the reports it summarises.
    """
    def load(path: Path) -> dict:
        return json.loads((path / "metrics.json").read_text(encoding="utf-8"))

    dynamic = load(ROOT / "reports" / "dynamic-mcts-poc-v0.2")
    by_policy = {row["policy"]: row for row in dynamic["policy_summary"]}
    v02_gap = (by_policy["MCTS"]["mean_utility"]
               - by_policy["NCCN"]["mean_utility"])

    robust = load(ROBUSTNESS)
    robust_ci = robust["utility_gap_ci95"]
    sensitivity = load(SENSITIVITY)
    budget = load(BUDGET)
    by_budget = {row["budget"]: row for row in budget["utility_gap_by_budget"]}
    envfix = load(ENVFIX)
    precision = load(ROOT / "reports" / "sensitivity-precision-v1.0")
    precision_1024 = next(row for row in precision["by_budget"]
                          if row["budget"] == "1024")
    patients = load(ROOT / "reports" / "sensitivity-patients-v1.1")
    patients_largest = patients["by_patient_count"][-1]
    replication = load(ROOT / "reports" / "cohort-replication-v1.2")
    pooled = next(row for row in replication["by_cohort"]
                  if row["cohort"] == "pooled")
    standardised = json.loads(
        (ROOT / "reports" / "cohort-replication-v1.2"
         / "metrics_posthoc_subtype.json").read_text(encoding="utf-8"))

    rows = [
        (f"v0.2 단일 시드\n(256 · {dynamic['cohort']['dynamic_test_patients']}명 · 1시드)",
         v02_gap, None, None, GRAY),
        (f"v0.3 강건성\n(256 · {robust['design']['patients']}명 · "
         f"{robust['design']['n_seeds']}시드)",
         robust_ci["mean"], robust_ci["ci95_low"], robust_ci["ci95_high"], BLUE),
        (f"v0.3 민감도 기준선\n(256 · {sensitivity['design']['patients']}명 · "
         f"{sensitivity['design']['n_seeds']}시드)",
         sensitivity["baseline_utility_gap"], None, None, AMBER),
        (f"v0.4 예산 256\n(256 · {budget['design']['patients']}명 · "
         f"{budget['design']['policy_seeds']}시드)",
         by_budget[256]["utility_gap_mean"],
         by_budget[256]["utility_gap_ci95_low"],
         by_budget[256]["utility_gap_ci95_high"], BLUE),
        (f"v0.4 예산 1024\n(1024 · {budget['design']['patients']}명 · "
         f"{budget['design']['policy_seeds']}시드)",
         by_budget[1024]["utility_gap_mean"],
         by_budget[1024]["utility_gap_ci95_low"],
         by_budget[1024]["utility_gap_ci95_high"], GREEN),
        (f"v0.5 편향 수정\n(1024 · {envfix['design']['patients']}명 · "
         f"{envfix['design']['seeds']}시드)",
         envfix["utility_gap_fixed"], None, None, GREEN),
        (f"v1.0 민감도 재실행\n(1024 · {precision['design']['patients']}명 · "
         f"{precision['design']['seeds']}시드)",
         precision_1024["baseline_utility_gap"], None, None, AMBER),
        (f"v1.1 환자 확대\n(1024 · {patients_largest['n_patients']}명 · "
         f"{patients['design']['seeds']}시드)",
         patients_largest["baseline_utility_gap"], None, None, AMBER),
        (f"v1.2 두 코호트 합산\n(1024 · {pooled['n_patients']}명 · "
         f"{replication['design']['seeds']}시드)",
         standardised["balanced_sample_mean"],
         standardised["balanced_sample_mean"]
         - 1.96 * standardised["balanced_sample_standard_error"],
         standardised["balanced_sample_mean"]
         + 1.96 * standardised["balanced_sample_standard_error"], GREEN),
        (f"v1.2 아형 표준화\n(사후 · {pooled['n_patients']}명)",
         standardised["prevalence_standardised_mean"],
         standardised["prevalence_standardised_mean"]
         - 1.96 * standardised["prevalence_standardised_standard_error"],
         standardised["prevalence_standardised_mean"]
         + 1.96 * standardised["prevalence_standardised_standard_error"], RED),
    ]
    labels = [row[0] for row in rows]
    positions = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    for position, (_, value, low, high, color) in zip(positions, rows):
        if low is not None:
            ax.plot([low, high], [position, position], color=color,
                    linewidth=2.4, solid_capstyle="round", alpha=0.55)
        ax.scatter(value, position, color=color, s=70, zorder=3)
        ax.text(value, position + 0.30, f"{value:+.4f}",
                ha="center", fontsize=9.5, color=INK)
    ax.axvline(0, color=INK, linewidth=1.1)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_ylim(len(rows) - 0.35, -0.7)
    ax.set_xlabel("효용 격차 (MCTS - NCCN), 가로선은 95% CI")
    ax.set_title("헤드라인 수치가 달라질 때마다 무엇이 바뀌었나\n"
                 "— 예산·표본/시드·환경, 그리고 마지막엔 표집 설계",
                 fontsize=12.5, pad=12)
    ax.annotate("예산 256 → 1024\n(행동 순서를 분해할 해상도)",
                xy=(by_budget[1024]["utility_gap_mean"], 4),
                xytext=(by_budget[1024]["utility_gap_mean"] + 0.0025, 3.1),
                fontsize=8.8, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, linewidth=1.1))
    ax.annotate("표본 8명·3시드로 축소한 설정",
                xy=(sensitivity["baseline_utility_gap"], 2),
                xytext=(sensitivity["baseline_utility_gap"] + 0.0025, 1.35),
                fontsize=8.8, color=AMBER,
                arrowprops=dict(arrowstyle="->", color=AMBER, linewidth=1.1))
    ax.annotate("같은 8명을 시드 12개로 다시 재니\n민감도 기준선이 여기로 올라왔다",
                xy=(precision_1024["baseline_utility_gap"], 6),
                xytext=(precision_1024["baseline_utility_gap"] - 0.0230, 6.55),
                fontsize=8.8, color=AMBER,
                arrowprops=dict(arrowstyle="->", color=AMBER, linewidth=1.1))
    ax.annotate("아형 구성을 실제 분포로 맞추면\n같은 40명이 여기로 내려온다",
                xy=(standardised["prevalence_standardised_mean"], 9),
                xytext=(standardised["prevalence_standardised_mean"] - 0.0205, 8.35),
                fontsize=8.8, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, linewidth=1.1))
    save(fig, ENVFIX, "fig23_results_reconciliation.png")

    print("  reconciliation:", {row[0].splitlines()[0]: round(row[1], 4)
                                for row in rows})


def main() -> None:
    fig18_budget_scaling()
    fig19_phase_agreement()
    fig20_sensitivity_tornado()
    fig21_seed_robustness()
    fig22_response_channel_bias()
    fig23_results_reconciliation()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
