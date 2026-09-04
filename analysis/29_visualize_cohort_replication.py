"""Figure for the disjoint-cohort replication (fig34).

A replication figure has one job: put the two cohorts side by side so the reader
can see for themselves whether the second one reproduces the first, rather than
being told it did. The parameter-by-parameter scatter does that in one glance -
points on the diagonal replicate, points off it do not.

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

REPORT_DIR = ROOT / "reports" / "cohort-replication-v1.2"
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
SHORT = {
    "reward.recurrence_free_year": "무재발 보상",
    "reward.acute_toxicity_penalty": "독성 페널티",
    "toxicity.chemo.intensified": "독성확률",
    "hazard.chemo.intensified.death": "항암 HR",
    "response.intensified.major": "반응확률",
    "hazard.response.major.death": "반응 HR",
}
JUDGEMENTS = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig34_replication() -> None:
    influence = pd.read_csv(TABLE_DIR / "parameter_influence.csv")
    per_patient = pd.read_csv(TABLE_DIR / "per_patient_gaps.csv")
    by_cohort = {row["cohort"]: row for row in METRICS["by_cohort"]}
    replication = METRICS["replication"]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.6),
                             gridspec_kw={"wspace": 0.30, "hspace": 0.42})

    # --- A: the two rankings, side by side ---------------------------------
    ax = axes[0][0]
    frame_a = influence[influence["cohort"] == "A"].set_index("parameter")
    frame_b = influence[influence["cohort"] == "B"].set_index("parameter")
    order = list(frame_a.sort_values("abs_gap_delta").index)
    positions = np.arange(len(order))
    height = 0.36
    ax.barh(positions + height / 2, [frame_a.loc[p, "abs_gap_delta"] for p in order],
            height=height, color=[RED if p in JUDGEMENTS else BLUE for p in order],
            label="코호트 A (20명)")
    ax.barh(positions - height / 2, [frame_b.loc[p, "abs_gap_delta"] for p in order],
            height=height, color=[RED if p in JUDGEMENTS else BLUE for p in order],
            alpha=0.45, hatch="///", edgecolor="white",
            label="코호트 B (겹치지 않는 20명)")
    ax.set_yticks(positions)
    ax.set_yticklabels([PRETTY[p] for p in order], fontsize=9)
    ax.set_xlabel("효용 격차에 준 최대 영향 |Δ|")
    # Both cohorts' bars start at zero and run most of the axis, so an in-axes
    # legend would sit on top of the data.
    ax.legend(frameon=False, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.17), ncol=2)
    met = replication["primary_prediction_met"]
    ax.set_title(
        "A. 사전 선언한 예측: B의 상위 2개는 가치판단\n"
        f"결과 — {'적중' if met else '빗나감'}",
        fontsize=12, pad=10, color=GREEN if met else RED)

    # --- B: parameter-by-parameter agreement -------------------------------
    ax = axes[0][1]
    xs = [frame_a.loc[p, "abs_gap_delta"] for p in order]
    ys = [frame_b.loc[p, "abs_gap_delta"] for p in order]
    limit = max(max(xs), max(ys)) * 1.25
    ax.plot([0, limit], [0, limit], color=GRAY, linewidth=1.1, linestyle="--",
            label="완전 일치선")
    # Points crowd near the diagonal, so labels go outward from it; neighbours on
    # the same side alternate above and below or their text collides.
    left_side = 0
    for parameter, x, y in zip(order, xs, ys):
        color = RED if parameter in JUDGEMENTS else BLUE
        ax.scatter(x, y, s=95, color=color, zorder=3)
        if x < limit / 2:
            offset, align = (10, 6), "left"
        else:
            offset = (-12, 8) if left_side % 2 == 0 else (-12, -16)
            align = "right"
            left_side += 1
        ax.annotate(SHORT[parameter], (x, y), textcoords="offset points",
                    xytext=offset, ha=align, fontsize=8.5, color=INK)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("코호트 A의 |Δ|")
    ax.set_ylabel("코호트 B의 |Δ|")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title(f"B. 대각선 위면 재현된 것\nSpearman "
                 f"{replication['spearman_a_vs_b']:+.2f}",
                 fontsize=12, pad=10, color=INK)

    # --- C: the headline number in each cohort -----------------------------
    ax = axes[1][0]
    labels = ["코호트 A\n(20명)", "코호트 B\n(20명)", "합산\n(40명)"]
    keys = ["A", "B", "pooled"]
    gaps = [by_cohort[k]["baseline_utility_gap"] for k in keys]
    errors = [1.96 * by_cohort[k]["baseline_seed_standard_error"] for k in keys]
    colors = [BLUE, GREEN, INK]
    positions = np.arange(len(keys))
    ax.bar(positions, gaps, yerr=errors, color=colors, width=0.55,
           error_kw={"ecolor": INK, "capsize": 4, "linewidth": 1.1})
    for x, gap in zip(positions, gaps):
        ax.text(x, gap + max(errors) + 0.0015, f"{gap:+.4f}",
                ha="center", fontsize=10, color=INK)
    ax.axhline(0, color=INK, linewidth=1.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel("기준선 효용 격차")
    ax.set_ylim(0, max(gaps) + max(errors) + 0.008)
    difference = replication["baseline_gap_difference"]
    # ASCII hyphen: Malgun Gothic has no glyph for U+2212.
    ax.set_title(f"C. 헤드라인은 코호트를 갈아도 유지된다\nB - A = {difference:+.4f}",
                 fontsize=12, pad=10, color=INK)

    # --- D: pooled per-patient spread --------------------------------------
    ax = axes[1][1]
    heterogeneity = METRICS["patient_heterogeneity"]
    block = per_patient.sort_values("baseline_utility_gap").reset_index(drop=True)
    bar_colors = [BLUE if row == "A" else GREEN for row in block["cohort"]]
    ax.bar(np.arange(len(block)), block["baseline_utility_gap"],
           color=bar_colors, width=0.8)
    mean = float(block["baseline_utility_gap"].mean())
    ax.axhline(mean, color=INK, linewidth=1.3, linestyle="--",
               label=f"평균 {mean:+.4f}")
    ax.axhline(0, color=INK, linewidth=1.0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE),
               plt.Rectangle((0, 0), 1, 1, color=GREEN)]
    ax.legend(handles + [plt.Line2D([0], [0], color=INK, linestyle="--")],
              ["코호트 A", "코호트 B", f"평균 {mean:+.4f}"],
              frameon=False, fontsize=8.5, loc="upper left")
    ax.set_xticks([])
    ax.set_xlabel(f"환자 {heterogeneity['n_patients']}명 (효용 격차 순)")
    ax.set_ylabel("환자별 효용 격차")
    negative = heterogeneity["patients_with_negative_gap"]
    ax.set_title(
        f"D. 두 코호트가 같은 분포에서 나왔다 (음수 {negative}명)\n"
        f"환자간 표준편차 A {heterogeneity['per_patient_gap_sd']['A']:.4f} · "
        f"B {heterogeneity['per_patient_gap_sd']['B']:.4f}",
        fontsize=12, pad=10, color=INK)

    fig.suptitle("v1.2 — 겹치지 않는 두 코호트로 순위를 재검했다 (오차막대는 95% CI)",
                 fontsize=13.5, y=0.98, color=INK)
    save(fig, "fig34_cohort_replication.png")


def fig35_subtype_standardisation() -> None:
    """Post-hoc: the sampling scheme's effect on the headline.

    Equal-per-subtype sampling keeps rare subtypes present, which is why every
    study since v0.3 used it. It also means the headline averages over a mix that
    is nothing like the cohort it was drawn from - and the gap turns out to vary
    several-fold by subtype, so the two averages are not the same number.
    """
    posthoc = json.loads(
        (REPORT_DIR / "metrics_posthoc_subtype.json").read_text(encoding="utf-8"))
    table = pd.read_csv(TABLE_DIR / "subtype_standardisation.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.9),
                             gridspec_kw={"wspace": 0.30, "width_ratios": [1.35, 1]})

    # --- left: gap by subtype, with how much of each we sampled -------------
    ax = axes[0]
    positions = np.arange(len(table))
    ax.bar(positions, table["mean_gap"],
           yerr=1.96 * table["standard_error"], color=BLUE, width=0.55,
           error_kw={"ecolor": INK, "capsize": 4, "linewidth": 1.1})
    ax.axhline(0, color=INK, linewidth=1.0)
    for x, row in zip(positions, table.itertuples()):
        ax.text(x, row.mean_gap + 1.96 * row.standard_error + 0.0025,
                f"{row.mean_gap:+.4f}", ha="center", fontsize=9.5, color=INK)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{row.subtype}\n표본 {row.sample_weight:.0%} · 실제 {row.population_share:.0%}"
         for row in table.itertuples()], fontsize=9)
    ax.set_ylabel("환자별 효용 격차 평균")
    spread = posthoc["subtype_gap_spread"]
    ax.set_title(
        f"아형에 따라 {spread['ratio']:.1f}배 갈린다 — 그런데 표본은 균등하게 뽑았다\n"
        f"가장 흔한 아형(전체의 "
        f"{posthoc['largest_subtype']['population_share']:.0%})이 표본에서는 25%",
        fontsize=12, pad=10, color=INK)

    # --- right: what that does to the headline ------------------------------
    ax = axes[1]
    values = [posthoc["balanced_sample_mean"],
              posthoc["prevalence_standardised_mean"]]
    errors = [1.96 * posthoc["balanced_sample_standard_error"],
              1.96 * posthoc["prevalence_standardised_standard_error"]]
    positions = np.arange(2)
    ax.bar(positions, values, yerr=errors, color=[BLUE, RED], width=0.5,
           error_kw={"ecolor": INK, "capsize": 4, "linewidth": 1.1})
    for x, value, error in zip(positions, values, errors):
        ax.text(x, value + error + 0.0018, f"{value:+.4f}",
                ha="center", fontsize=10.5, color=INK)
    ax.axhline(0, color=INK, linewidth=1.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(["아형 균등 표본\n(우리가 인용해 온 값)",
                        "실제 아형 구성으로\n표준화"], fontsize=9.5)
    ax.set_ylabel("효용 격차")
    ax.set_ylim(0, max(v + e for v, e in zip(values, errors)) + 0.006)
    ax.set_title(
        f"헤드라인이 {1 - posthoc['ratio_standardised_to_balanced']:.0%} 줄어든다",
        fontsize=12, pad=10, color=RED)

    fig.suptitle("v1.2 사후 분석 — 아형 균등 표집이 헤드라인을 얼마나 키웠나 "
                 "(사전 계획 아님, 오차막대는 95% CI)",
                 fontsize=13, y=1.03, color=INK)
    save(fig, "fig35_subtype_standardisation.png")


def main() -> None:
    fig34_replication()
    fig35_subtype_standardisation()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
