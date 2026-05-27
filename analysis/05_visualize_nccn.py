"""NCCN 일치율(concordance) 시각화.

04_nccn_policy.py 의 결과(patients_with_nccn.csv)를 읽어 두 그림을 저장한다:
  fig08_nccn_concordance.png             — 결정 4개 전체 일치율 막대
  fig09_nccn_concordance_by_subtype.png  — subtype × 결정 그룹 막대
"""

from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "data" / "processed" / "patients_with_nccn.csv"
OUT = ROOT / "data" / "processed" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# 사이트 라이트 톤 (03_visualize.py 와 동일)
matplotlib.rcParams.update({
    "font.family": "Malgun Gothic",
    "axes.unicode_minus": False,
    "figure.facecolor": "#fafaf9",
    "axes.facecolor": "#fafaf9",
    "axes.edgecolor": "#1c1917",
    "axes.labelcolor": "#1c1917",
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
SUBTYPE_COLORS = {
    "HR+/HER2-": "#0891b2",
    "HR+/HER2+": "#7c3aed",
    "HR-/HER2+": "#2563eb",
    "TNBC":      "#dc2626",
}

DECISIONS = ["surgery", "chemo", "hormone", "radio"]
DECISION_KOR = {"surgery": "수술", "chemo": "보조항암", "hormone": "호르몬", "radio": "방사선"}


def norm_surgery(v):
    if pd.isna(v):
        return None
    return "BCS" if v == "BREAST CONSERVING" else "MAST"


def main():
    df = pd.read_csv(CSV)
    df = df.dropna(subset=["subtype", "os_event"]).copy()
    df["actual_surgery"] = df["surgery"].apply(norm_surgery)

    # --- 결정별 전체 일치율 -----------------------------------------
    pairs = {
        "surgery": ("rec_surgery", "actual_surgery"),
        "chemo":   ("rec_chemo",   "chemo"),
        "hormone": ("rec_hormone", "hormone"),
        "radio":   ("rec_radio",   "radio"),
    }
    overall = []
    for name in DECISIONS:
        rec_col, act_col = pairs[name]
        rec, act = df[rec_col], df[act_col]
        ok = (~rec.isna()) & (~act.isna())
        match = ok & (rec == act)
        overall.append((name, ok.sum(), match.sum(), match.sum() / max(ok.sum(), 1) * 100))

    fig, ax = plt.subplots(figsize=(8, 4.8))
    labels = [DECISION_KOR[n] for n, _, _, _ in overall]
    pcts = [p for _, _, _, p in overall]
    ns = [(m, n) for _, n, m, _ in overall]
    bars = ax.bar(labels, pcts, color=INK, alpha=0.85, edgecolor="white")
    for b, pct, (m, n) in zip(bars, pcts, ns):
        ax.text(b.get_x() + b.get_width() / 2, pct + 1.5,
                f"{pct:.1f}%\n({m}/{n})", ha="center", va="bottom",
                fontsize=9, color=INK)
    ax.set_ylabel("NCCN 일치율 (%)")
    ax.set_ylim(0, 100)
    ax.set_title("결정별 NCCN 권고-실제 일치율 (METABRIC, n≤1980)", pad=14)
    fig.savefig(OUT / "fig08_nccn_concordance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(OUT / "fig08_nccn_concordance.png")

    # --- subtype × 결정 그룹 막대 ----------------------------------
    subs = ["HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "TNBC"]
    matrix = {}
    counts = {}
    for sub in subs:
        g = df[df["subtype"] == sub]
        counts[sub] = len(g)
        row = []
        for name in DECISIONS:
            rec_col, act_col = pairs[name]
            rec, act = g[rec_col], g[act_col]
            ok = (~rec.isna()) & (~act.isna())
            match = ok & (rec == act)
            row.append(match.sum() / max(ok.sum(), 1) * 100)
        matrix[sub] = row

    x = np.arange(len(DECISIONS))
    width = 0.2
    fig, ax = plt.subplots(figsize=(9.5, 5))
    for i, sub in enumerate(subs):
        offset = (i - 1.5) * width
        bars = ax.bar(
            x + offset, matrix[sub], width,
            color=SUBTYPE_COLORS[sub], edgecolor="white",
            label=f"{sub} (n={counts[sub]})",
        )
        for b, v in zip(bars, matrix[sub]):
            ax.text(b.get_x() + b.get_width() / 2, v + 1.5,
                    f"{v:.0f}", ha="center", va="bottom",
                    fontsize=8, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([DECISION_KOR[n] for n in DECISIONS])
    ax.set_ylabel("NCCN 일치율 (%)")
    ax.set_ylim(0, 100)
    ax.set_title("분자아형 × 결정 노드별 NCCN 일치율", pad=14)
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.savefig(OUT / "fig09_nccn_concordance_by_subtype.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(OUT / "fig09_nccn_concordance_by_subtype.png")


if __name__ == "__main__":
    main()
