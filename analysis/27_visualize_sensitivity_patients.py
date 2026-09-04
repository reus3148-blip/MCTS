"""Figure for the cohort-size sensitivity re-run (fig33).

v1.0's figure led with the baseline moving further under a change of seeds than
under any effect being ranked. The question this run asks is the same one aimed
at the remaining knob, so the figure keeps that shape: what the cohort size does
to the headline number first, then how much of the ranking is a property of
these particular patients.

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

REPORT_DIR = ROOT / "reports" / "sensitivity-patients-v1.1"
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
JUDGEMENTS = {"reward.recurrence_free_year", "reward.acute_toxicity_penalty"}


def save(fig: plt.Figure, filename: str) -> None:
    path = FIGURE_DIR / filename
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(path, PUBLIC_DIR / filename)
    print(path)


def fig33_patients() -> None:
    influence = pd.read_csv(TABLE_DIR / "parameter_influence.csv")
    per_patient = pd.read_csv(TABLE_DIR / "per_patient_gaps.csv")
    by_size = METRICS["by_patient_count"]
    sizes = [row["n_patients"] for row in by_size]
    largest = sizes[-1]
    heterogeneity = METRICS["patient_heterogeneity"]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.4),
                             gridspec_kw={"wspace": 0.30, "hspace": 0.42})

    # --- A: does the headline number move with the cohort? ------------------
    ax = axes[0][0]
    gaps = [row["baseline_utility_gap"] for row in by_size]
    seed_se = [row["baseline_seed_standard_error"] for row in by_size]
    positions = np.arange(len(sizes))
    ax.errorbar(positions, gaps, yerr=[1.96 * s for s in seed_se],
                fmt="o-", color=BLUE, ecolor=INK, capsize=4, linewidth=1.6,
                markersize=7)
    patient_se = heterogeneity["between_patient_standard_error"]
    ax.fill_between(
        positions,
        [g - 1.96 * patient_se[str(n)] for g, n in zip(gaps, sizes)],
        [g + 1.96 * patient_se[str(n)] for g, n in zip(gaps, sizes)],
        color=RED, alpha=0.12, label="환자 재추출 95% 구간")
    for x, gap in zip(positions, gaps):
        ax.annotate(f"{gap:+.4f}", (x, gap), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=10, color=INK)
    ax.axhline(0, color=INK, linewidth=1.0)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{n}명" for n in sizes])
    ax.set_xlim(-0.4, len(sizes) - 0.6)
    ax.set_ylabel("기준선 효용 격차")
    ax.set_xlabel("코호트 크기 (시드 12 · 예산 1024 고정)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    shift = gaps[-1] - gaps[0]
    ax.set_title(f"A. 환자를 8명 → {largest}명으로 늘리면 {shift:+.4f}",
                 fontsize=12, pad=10,
                 color=RED if abs(shift) > 0.005 else GREEN)

    # --- B: where that uncertainty comes from -------------------------------
    ax = axes[0][1]
    block = (per_patient[per_patient.n_patients == largest]
             .sort_values("baseline_utility_gap")
             .reset_index(drop=True))
    colors = [GREEN if v >= 0 else RED for v in block["baseline_utility_gap"]]
    ax.bar(np.arange(len(block)), block["baseline_utility_gap"],
           color=colors, width=0.72)
    mean = float(block["baseline_utility_gap"].mean())
    ax.axhline(mean, color=BLUE, linewidth=1.4, linestyle="--",
               label=f"평균 {mean:+.4f}")
    ax.axhline(0, color=INK, linewidth=1.0)
    ax.set_xticks([])
    ax.set_xlabel(f"환자 {largest}명 (효용 격차 순)")
    ax.set_ylabel("환자별 효용 격차")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    negative = heterogeneity["patients_with_negative_gap"]
    ratio = (heterogeneity["per_patient_gap_max"]
             / heterogeneity["per_patient_gap_min"])
    median = float(block["baseline_utility_gap"].median())
    ax.set_title(
        f"B. 부호는 {largest}명 전원 같다 (음수 {negative}명) — 크기는 {ratio:.0f}배 갈린다\n"
        f"환자간 표준편차 {heterogeneity['per_patient_gap_sd']:.4f} · "
        f"중앙값 {median:+.4f}",
        fontsize=12, pad=10, color=INK)

    # --- C: the ranking at the largest cohort -------------------------------
    ax = axes[1][0]
    frame = (influence[influence["n_patients"] == largest]
             .sort_values("abs_gap_delta"))
    names = [PRETTY.get(p, p) for p in frame["parameter"]]
    colors = [RED if p in JUDGEMENTS else BLUE for p in frame["parameter"]]
    positions = np.arange(len(frame))
    ax.barh(positions, frame["abs_gap_delta"], color=colors, height=0.6)
    ax.errorbar(frame["abs_gap_delta"], positions,
                xerr=1.96 * frame["delta_standard_error"], fmt="none",
                ecolor=INK, capsize=3, linewidth=1.1)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("효용 격차에 준 최대 영향 |Δ|")
    ax.set_xlim(0, float(frame["abs_gap_delta"].max()
                         + 1.96 * frame["delta_standard_error"].max()) * 1.12)
    row = next(r for r in by_size if r["n_patients"] == largest)
    verdict = ("상위 2개가 가치판단" if row["top_two_are_value_judgements"]
               else "상위 2개가 가치판단 아님")
    ax.set_title(f"C. 환자 {largest}명에서의 순위 — {verdict}\n"
                 f"|z| ≥ 2 는 6개 중 "
                 f"{row['parameters_distinguishable_from_zero']}개",
                 fontsize=12, pad=10,
                 color=GREEN if row["top_two_are_value_judgements"] else RED)
    handles = [plt.Rectangle((0, 0), 1, 1, color=RED),
               plt.Rectangle((0, 0), 1, 1, color=BLUE)]
    ax.legend(handles, ["가치판단", "데이터로 추정 가능"],
              frameon=False, fontsize=8.5, loc="lower right")

    # --- D: would other patients rank them the same? ------------------------
    ax = axes[1][1]
    shares = [METRICS["patient_bootstrap"][str(n)]
              ["top_two_are_value_judgements_share"] * 100 for n in sizes]
    positions = np.arange(len(sizes))
    bars = ax.bar(positions, shares, color=BLUE, width=0.55)
    for bar, share in zip(bars, shares):
        ax.text(bar.get_x() + bar.get_width() / 2, share + 1.5,
                f"{share:.0f}%", ha="center", fontsize=10, color=INK)
    ax.axhline(50, color=GRAY, linewidth=1.1, linestyle=":")
    ax.text(len(sizes) - 0.55, 52, "동전 던지기", fontsize=8.5, color=GRAY)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{n}명" for n in sizes])
    ax.set_ylim(0, 105)
    ax.set_ylabel("상위 2개가 가치판단인 비율 (%)")
    ax.set_xlabel(f"환자 부트스트랩 {METRICS['design']['bootstrap_replicates']:,}회")
    ax.set_title("D. 다른 환자였다면 같은 순위가 나왔을까",
                 fontsize=12, pad=10, color=INK)

    fig.suptitle(
        "v1.1 — 시드·예산을 고정하고 환자 수만 늘렸을 때 (오차막대는 95% CI)",
        fontsize=13.5, y=0.98, color=INK)
    save(fig, "fig33_sensitivity_patients.png")


def main() -> None:
    fig33_patients()
    print(f"\ncopied to {PUBLIC_DIR}")


if __name__ == "__main__":
    main()
