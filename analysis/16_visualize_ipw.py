"""Figures for the IPW target trial emulation (fig24 - fig26).

The argument of v0.6 is visual before it is numerical: two arms that barely
overlap, a weighting that only balances once the non-overlapping patients are
removed, and an effect estimate that changes sign once it is adjusted. Each panel
reads from the committed tables of ``reports/ipw-target-trial-v0.6``.
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

from analysis.causal.ipw import weighted_kaplan_meier  # noqa: E402

REPORT_DIR = ROOT / "reports" / "ipw-target-trial-v0.6"
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
TRIM = METRICS["design"]["primary_trim"]


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig24_overlap() -> None:
    """The positivity problem: two arms that barely share a propensity range."""
    scores = pd.read_csv(TABLE_DIR / "propensity_full_cohort.csv")
    treated = scores[scores["chemo"] == 1]["propensity_full_cohort"]
    control = scores[scores["chemo"] == 0]["propensity_full_cohort"]
    bins = np.linspace(0, 1, 41)

    fig, ax = plt.subplots(figsize=(10.2, 4.6))
    ax.hist(control, bins=bins, color=BLUE, alpha=0.6, label=f"항암 안 함 (n={len(control)})")
    ax.hist(treated, bins=bins, color=RED, alpha=0.6, label=f"항암 함 (n={len(treated)})")
    ax.axvspan(TRIM[0], TRIM[1], color=GREEN, alpha=0.10)
    for bound in TRIM:
        ax.axvline(bound, color=GREEN, linestyle="--", linewidth=1.4)
    retained = METRICS["design"]["retained_pct"]
    ax.text((TRIM[0] + TRIM[1]) / 2, ax.get_ylim()[1] * 0.92,
            f"분석 대상 구간\n{TRIM[0]:.2f} ~ {TRIM[1]:.2f}\n"
            f"(전체의 {retained:.0f}%)",
            ha="center", va="top", fontsize=9.5, color=GREEN)
    ax.axvline(METRICS["positivity"]["median_propensity_control"],
               color=BLUE, linewidth=1.2, linestyle=":")
    ax.axvline(METRICS["positivity"]["median_propensity_treated"],
               color=RED, linewidth=1.2, linestyle=":")
    ax.set_xlabel("추정 항암 시행 확률 (propensity score)")
    ax.set_ylabel("환자 수")
    ax.set_title("두 군이 거의 겹치지 않는다 — positivity 위반\n"
                 f"중앙값 {METRICS['positivity']['median_propensity_control']:.3f} "
                 f"vs {METRICS['positivity']['median_propensity_treated']:.3f}",
                 fontsize=12.5, pad=12)
    ax.legend(frameon=False, fontsize=9.5)
    save(fig, "fig24_propensity_overlap.png")


def fig25_love_plot() -> None:
    """Balance before and after weighting, inside the overlap population."""
    balance = pd.read_csv(TABLE_DIR / "covariate_balance.csv")
    untrimmed = pd.read_csv(TABLE_DIR / "covariate_balance_untrimmed.csv")
    pretty = {
        "age": "나이", "tumor_size_mm": "종양 크기", "lymph_pos": "림프절 전이 수",
        "stage": "병기", "grade": "분화도", "er": "ER", "pr": "PR", "her2": "HER2",
        "menopause=Pre": "폐경 전", "subtype=HR+/HER2+": "HR+/HER2+",
        "subtype=HR-/HER2+": "HR-/HER2+", "subtype=TNBC": "TNBC",
    }
    # All three series describe the FULL eligible cohort's covariates, so the
    # crude and the untrimmed-weighted points come from the untrimmed table;
    # mixing populations across series would make the plot unreadable.
    crude = untrimmed.set_index("covariate")["smd_crude"].abs()
    untrimmed_weighted = untrimmed.set_index("covariate")["smd_weighted"].abs()
    order = balance.reindex(
        balance["covariate"].map(crude).sort_values().index).reset_index(drop=True)
    names = [pretty.get(c, c) for c in order["covariate"]]
    positions = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    ax.scatter([crude.get(c, np.nan) for c in order["covariate"]], positions,
               color=GRAY, s=52, label="가중 전 (전체 코호트)")
    ax.scatter([untrimmed_weighted.get(c, np.nan) for c in order["covariate"]],
               positions, color=AMBER, marker="^", s=50,
               label="가중 후 · 트리밍 없음 (실패)")
    ax.scatter(order["smd_weighted"].abs(), positions, color=GREEN, s=52,
               label="가중 후 · 트리밍 적용")
    ax.axvline(0.1, color=RED, linestyle="--", linewidth=1.3)
    ax.text(0.115, -0.6, "균형 기준 0.1", color=RED, fontsize=9)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=9.5)
    ax.set_xlabel("|표준화 평균 차이|")
    ax.set_xscale("symlog", linthresh=0.05)
    ax.set_xticks([0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0])
    ax.set_xticklabels(["0", "0.05", "0.1", "0.25", "0.5", "1.0", "2.0"])
    ax.set_title("트리밍 없이는 가중치가 균형을 만들지 못한다\n"
                 f"최악 |SMD| {crude.max():.2f} → "
                 f"{METRICS['positivity']['untrimmed_worst_abs_smd']:.2f}(트리밍 없음) → "
                 f"{METRICS['balance']['worst_abs_smd_weighted']:.3f}(트리밍)",
                 fontsize=12, pad=12)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    save(fig, "fig25_covariate_balance.png")


def fig26_weighted_survival() -> None:
    """Naive versus IPW survival, and the estimate flipping sign."""
    weights_table = pd.read_csv(TABLE_DIR / "propensity_and_weights_trimmed.csv")
    raw = pd.read_csv(ROOT / "data" / "processed" / "patients_with_nccn.csv")
    merged = weights_table.merge(
        raw[["patient_id", "os_months", "os_event"]], on="patient_id", how="left")

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

    ax = axes[0]
    for arm, color, label in ((1, RED, "항암 함"), (0, BLUE, "항암 안 함")):
        subset = merged[merged["chemo"] == arm]
        for weights, style, alpha in (
            (np.ones(len(subset)), ":", 0.55),
            (subset["stabilized_weight"].to_numpy(), "-", 1.0),
        ):
            curve = weighted_kaplan_meier(
                subset["os_months"].to_numpy(dtype=float),
                subset["os_event"].to_numpy(dtype=float), weights)
            times = np.concatenate([[0.0], curve.times])
            survival = np.concatenate([[1.0], curve.survival])
            ax.step(times, survival, where="post", color=color,
                    linestyle=style, alpha=alpha, linewidth=2.0 if style == "-" else 1.4,
                    label=f"{label} ({'IPW' if style == '-' else '조정 전'})")
    ax.axvline(60, color=INK, linewidth=1.0, linestyle="--")
    ax.text(62, 0.42, "5년", fontsize=9, color=INK)
    ax.set_xlim(0, 120)
    ax.set_ylim(0.4, 1.0)
    ax.set_xlabel("추적 기간 (개월)")
    ax.set_ylabel("생존 확률")
    ax.set_title("IPW 가중 생존곡선 (겹침 구간)", fontsize=12, pad=10)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")

    ax = axes[1]
    naive = METRICS["naive"]["risk_difference"]
    ipw = METRICS["ipw"]["risk_difference"]
    low, high = METRICS["ipw_ci95"]["risk_difference"]
    ax.errorbar([ipw], [0], xerr=[[ipw - low], [high - ipw]], fmt="o",
                color=GREEN, markersize=9, capsize=5, linewidth=2.0)
    ax.scatter([naive], [1], color=GRAY, s=80)
    ax.axvline(0, color=INK, linewidth=1.1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([f"IPW 조정\n(95% CI)", "조정 전"], fontsize=10)
    ax.set_ylim(-0.6, 1.6)
    ax.text(naive, 1.22, f"{naive:+.3f}", ha="center", fontsize=10, color=GRAY)
    ax.text(ipw, 0.25, f"{ipw:+.3f}  [{low:+.3f}, {high:+.3f}]",
            ha="center", fontsize=10, color=GREEN)
    ax.set_xlabel("5년 사망 위험 차이 (항암 - 미시행)")
    ax.set_title("조정하면 부호가 바뀌지만 CI는 0을 포함한다",
                 fontsize=12, pad=10)
    save(fig, "fig26_ipw_effect.png")


def main() -> None:
    fig24_overlap()
    fig25_love_plot()
    fig26_weighted_survival()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
