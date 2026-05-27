"""METABRIC patients.csv 1차 시각화.

분포 4종 + Kaplan-Meier 생존곡선 3종을 data/processed/figures/ 에 PNG로 저장.

분포:
  fig01_age.png       — 나이 분포 (히스토그램)
  fig02_stage.png     — 병기 분포 (막대)
  fig03_subtype.png   — 분자아형 분포 (막대)
  fig04_treatment.png — 치료조합 상위 10 (가로 막대)

생존곡선 (lifelines KaplanMeierFitter):
  fig05_km_os_subtype.png  — subtype별 OS
  fig06_km_os_stage.png    — 병기별 OS
  fig07_km_rfs_subtype.png — subtype별 RFS

* 임상마커가 일괄 결측인 cohort (~529명)는 시각화에서 제외.
"""

from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "data" / "processed" / "patients.csv"
OUT = ROOT / "data" / "processed" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# --- 스타일 (사이트 라이트 톤과 어울리게) -----------------------------
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
SUB = "#57534e"

# subtype 4종 라이트 친화적 컬러 (사이트와 일관)
SUBTYPE_COLORS = {
    "HR+/HER2-": "#0891b2",  # cyan-700
    "HR+/HER2+": "#7c3aed",  # violet-600
    "HR-/HER2+": "#2563eb",  # blue-600
    "TNBC":      "#dc2626",  # red-600
}

STAGE_COLORS = {
    0: "#94a3b8",
    1: "#22c55e",
    2: "#f59e0b",
    3: "#ef4444",
    4: "#7f1d1d",
}


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --- 데이터 ----------------------------------------------------------
df = pd.read_csv(CSV)
print(f"loaded: {len(df)} rows")

# 임상마커 일괄 결측 cohort 제외 (subtype + os_event 둘 다 결측인 행)
analysis = df.dropna(subset=["subtype", "os_event", "os_months"]).copy()
print(f"after dropping cohort-wide NA: {len(analysis)} rows")

# --- 1) 나이 분포 ----------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.hist(analysis["age"].dropna(), bins=30, color=INK, alpha=0.85, edgecolor="white")
ax.set_xlabel("진단 시 나이 (세)")
ax.set_ylabel("환자 수")
ax.set_title(f"진단 시 나이 분포  (n={analysis['age'].notna().sum()})", pad=14)
median = analysis["age"].median()
ax.axvline(median, color="#dc2626", linestyle="--", linewidth=1.2,
           label=f"중앙값 = {median:.1f}세")
ax.legend(frameon=False)
print(save(fig, "fig01_age.png"))

# --- 2) 병기 분포 ----------------------------------------------------
stage_counts = analysis["stage"].value_counts(dropna=False).sort_index()
fig, ax = plt.subplots(figsize=(8, 4.5))
labels = [("결측" if pd.isna(i) else f"Stage {int(i)}") for i in stage_counts.index]
colors = [STAGE_COLORS.get(int(i), "#cbd5e1") if pd.notna(i) else "#d6d3d1"
          for i in stage_counts.index]
bars = ax.bar(labels, stage_counts.values, color=colors, edgecolor="white")
for b, v in zip(bars, stage_counts.values):
    ax.text(b.get_x() + b.get_width()/2, v + 15, str(v),
            ha="center", va="bottom", fontsize=9, color=INK)
ax.set_ylabel("환자 수")
ax.set_title("병기 (Tumor Stage) 분포", pad=14)
print(save(fig, "fig02_stage.png"))

# --- 3) 분자아형 분포 ------------------------------------------------
subtype_counts = analysis["subtype"].value_counts()
fig, ax = plt.subplots(figsize=(8, 4.5))
colors = [SUBTYPE_COLORS.get(s, "#94a3b8") for s in subtype_counts.index]
bars = ax.bar(subtype_counts.index, subtype_counts.values, color=colors, edgecolor="white")
for b, v in zip(bars, subtype_counts.values):
    pct = v / subtype_counts.sum() * 100
    ax.text(b.get_x() + b.get_width()/2, v + 20,
            f"{v}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=9, color=INK)
ax.set_ylabel("환자 수")
ax.set_title(f"분자아형 분포  (n={subtype_counts.sum()})", pad=14)
print(save(fig, "fig03_subtype.png"))

# --- 4) 치료조합 상위 10 ---------------------------------------------
treat_cols = ["surgery", "chemo", "hormone", "radio"]
combo_df = analysis[treat_cols].copy()
for c in combo_df.columns:
    combo_df[c] = combo_df[c].map(lambda v: "?" if pd.isna(v) else (
        ("MAST" if v == "MASTECTOMY" else "BCS") if c == "surgery"
        else ("Y" if int(v) == 1 else "N")
    ))
combos = (combo_df.assign(
    label=combo_df["surgery"] + " / Chemo " + combo_df["chemo"] +
          " / Hormone " + combo_df["hormone"] + " / Radio " + combo_df["radio"]
)["label"].value_counts().head(10))

fig, ax = plt.subplots(figsize=(9, 5.5))
bars = ax.barh(range(len(combos))[::-1], combos.values, color=INK, alpha=0.85,
               edgecolor="white")
ax.set_yticks(range(len(combos))[::-1])
ax.set_yticklabels(combos.index, fontsize=9)
ax.set_xlabel("환자 수")
ax.set_title("치료조합 상위 10  (수술 / 항암 / 호르몬 / 방사선)", pad=14)
for b, v in zip(bars, combos.values):
    ax.text(v + 4, b.get_y() + b.get_height()/2, str(v),
            va="center", fontsize=9, color=INK)
print(save(fig, "fig04_treatment.png"))


# --- KM 곡선 공통 함수 ----------------------------------------------
def km_by_group(data, time_col, event_col, group_col, color_map, title, fname):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    kmf = KaplanMeierFitter()
    groups = data.dropna(subset=[time_col, event_col, group_col])
    for name, sub in groups.groupby(group_col):
        kmf.fit(sub[time_col], event_observed=sub[event_col], label=f"{name} (n={len(sub)})")
        color = color_map.get(
            name if not isinstance(name, float) else int(name),
            "#94a3b8"
        )
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color, linewidth=1.8)
    ax.set_xlabel("개월")
    ax.set_ylabel("생존 확률")
    ax.set_ylim(0, 1.02)
    ax.set_title(title, pad=14)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    return save(fig, fname)


# --- 5) KM OS by subtype --------------------------------------------
print(km_by_group(
    analysis, "os_months", "os_event", "subtype", SUBTYPE_COLORS,
    "전체 생존 (OS) — 분자아형별",
    "fig05_km_os_subtype.png",
))

# --- 6) KM OS by stage ----------------------------------------------
stage_view = analysis.dropna(subset=["stage"]).copy()
stage_view["stage"] = stage_view["stage"].astype(int)
print(km_by_group(
    stage_view, "os_months", "os_event", "stage",
    {0: "#94a3b8", 1: "#22c55e", 2: "#f59e0b", 3: "#ef4444", 4: "#7f1d1d"},
    "전체 생존 (OS) — 병기별",
    "fig06_km_os_stage.png",
))

# --- 7) KM RFS by subtype -------------------------------------------
print(km_by_group(
    analysis, "rfs_months", "rfs_event", "subtype", SUBTYPE_COLORS,
    "무재발 생존 (RFS) — 분자아형별",
    "fig07_km_rfs_subtype.png",
))


# --- 요약 출력 --------------------------------------------------------
print("\n=== 요약 ===")
print(f"분석 대상: {len(analysis)}명")
print(f"OS event (사망): {int(analysis['os_event'].sum())}명 "
      f"({analysis['os_event'].mean() * 100:.1f}%)")
print(f"RFS event (재발): {int(analysis['rfs_event'].sum())}명 "
      f"({analysis['rfs_event'].mean() * 100:.1f}%)")
print(f"중앙 추적기간 (OS): {analysis['os_months'].median():.1f}개월 "
      f"({analysis['os_months'].median()/12:.1f}년)")
print(f"\n분자아형별 5년 OS:")
for sub, group in analysis.groupby("subtype"):
    kmf = KaplanMeierFitter()
    kmf.fit(group["os_months"], event_observed=group["os_event"])
    s60 = kmf.survival_function_at_times(60).iloc[0]
    print(f"  {sub:<12} n={len(group):<5} 5y OS = {s60*100:.1f}%")
