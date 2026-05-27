"""NCCN 침습성 유방암 가이드라인 정책 A (1차 단순화 버전).

이 모듈은 METABRIC patients.csv 의 환자 상태(state)에 NCCN 룰을 적용해
4개 결정 노드의 권고 치료를 계산하고, 실제 받은 치료와의 일치율을 산출한다.

본 연구 scope:
  - METABRIC 시기(1977-2005) 무관 골격 결정 4개에 한정
  - HER2 표적치료·CDK4/6 inhibitor·면역치료는 행동 공간에서 제외

룰 1차 버전 (사용자 합의 2026-05-27):
  ① 수술  : 종양 ≤ 30mm AND stage ≤ 2  → BCS,  else Mastectomy
  ② 항암  : HER2+ stage≥1 OR TNBC stage≥1 OR
            HR+/HER2- AND (림프절 양성 OR 종양>20mm OR grade=3) → Chemo
  ③ 호르몬: HR+ (ER+ 또는 PR+) → Hormone
  ④ 방사선: BCS → Radio,  Mastectomy AND (림프절≥4 OR stage≥3) → Radio
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CSV_IN = ROOT / "data" / "processed" / "patients.csv"
CSV_OUT = ROOT / "data" / "processed" / "patients_with_nccn.csv"


# --- 정책 함수 -------------------------------------------------------

def nccn_surgery(state: pd.Series) -> Optional[str]:
    """수술 종류 권고. 결정에 필요한 변수가 결측이면 None."""
    size = state.get("tumor_size_mm")
    stage = state.get("stage")
    if pd.isna(size) or pd.isna(stage):
        return None
    if size <= 30 and stage <= 2:
        return "BCS"
    return "MAST"


def nccn_chemo(state: pd.Series) -> Optional[int]:
    """보조 항암 권고."""
    subtype = state.get("subtype")
    stage = state.get("stage")
    if pd.isna(subtype) or pd.isna(stage):
        return None
    if stage < 1:  # in situ (stage 0)
        return 0

    if subtype == "HR-/HER2+":
        return 1
    if subtype == "TNBC":
        return 1
    if subtype == "HR+/HER2+":
        return 1  # 시기 무관 골격에서 HER2 양성은 항암 권고

    # HR+/HER2- 고위험 판단
    if subtype == "HR+/HER2-":
        lymph = state.get("lymph_pos")
        size = state.get("tumor_size_mm")
        grade = state.get("grade")
        high_risk = False
        if pd.notna(lymph) and lymph >= 1:
            high_risk = True
        if pd.notna(size) and size > 20:
            high_risk = True
        if pd.notna(grade) and grade == 3:
            high_risk = True
        return 1 if high_risk else 0

    return None


def nccn_hormone(state: pd.Series) -> Optional[int]:
    """호르몬 치료 권고. ER 또는 PR 양성이면 권고."""
    er = state.get("er")
    pr = state.get("pr")
    if pd.isna(er) and pd.isna(pr):
        return None
    er_pos = (er == 1) if pd.notna(er) else False
    pr_pos = (pr == 1) if pd.notna(pr) else False
    return 1 if (er_pos or pr_pos) else 0


def nccn_radio(state: pd.Series, recommended_surgery: Optional[str]) -> Optional[int]:
    """방사선 권고. 수술 종류 결정 결과에 의존한다."""
    if recommended_surgery is None:
        return None
    if recommended_surgery == "BCS":
        return 1  # BCS 후 방사선은 사실상 필수
    # Mastectomy
    lymph = state.get("lymph_pos")
    stage = state.get("stage")
    pmrt = False
    if pd.notna(lymph) and lymph >= 4:
        pmrt = True
    if pd.notna(stage) and stage >= 3:
        pmrt = True
    return 1 if pmrt else 0


def apply_nccn(row: pd.Series) -> pd.Series:
    surgery = nccn_surgery(row)
    chemo = nccn_chemo(row)
    hormone = nccn_hormone(row)
    radio = nccn_radio(row, surgery)
    return pd.Series({
        "rec_surgery": surgery,
        "rec_chemo": chemo,
        "rec_hormone": hormone,
        "rec_radio": radio,
    })


# --- 일치율 분석 ----------------------------------------------------

def normalize_actual_surgery(v):
    if pd.isna(v):
        return None
    return "BCS" if v == "BREAST CONSERVING" else "MAST"


def concordance_per_decision(df: pd.DataFrame) -> pd.DataFrame:
    """각 결정 노드의 일치/불일치/결측 집계."""
    rows = []
    pairs = [
        ("surgery",  "rec_surgery", df["surgery"].apply(normalize_actual_surgery)),
        ("chemo",    "rec_chemo",   df["chemo"]),
        ("hormone",  "rec_hormone", df["hormone"]),
        ("radio",    "rec_radio",   df["radio"]),
    ]
    for name, rec_col, actual in pairs:
        rec = df[rec_col]
        compared = (~rec.isna()) & (~actual.isna())
        match = compared & (rec == actual)
        rows.append({
            "decision":   name,
            "compared_n": int(compared.sum()),
            "match_n":    int(match.sum()),
            "concordance_pct": round(match.sum() / max(compared.sum(), 1) * 100, 1),
            "missing_n":  int((~compared).sum()),
        })
    return pd.DataFrame(rows)


def concordance_by_subtype(df: pd.DataFrame) -> pd.DataFrame:
    """subtype × decision 매트릭스."""
    df = df.copy()
    df["actual_surgery"] = df["surgery"].apply(normalize_actual_surgery)
    out = []
    for sub, g in df.groupby("subtype"):
        row = {"subtype": sub, "n": len(g)}
        for name, rec_col, actual_col in [
            ("surgery", "rec_surgery", "actual_surgery"),
            ("chemo",   "rec_chemo",   "chemo"),
            ("hormone", "rec_hormone", "hormone"),
            ("radio",   "rec_radio",   "radio"),
        ]:
            rec, act = g[rec_col], g[actual_col]
            ok = (~rec.isna()) & (~act.isna())
            match = ok & (rec == act)
            row[name] = round(match.sum() / max(ok.sum(), 1) * 100, 1)
        out.append(row)
    return pd.DataFrame(out)


def main():
    df = pd.read_csv(CSV_IN)
    print(f"loaded {len(df)} patients")

    rec = df.apply(apply_nccn, axis=1)
    out = pd.concat([df, rec], axis=1)
    out.to_csv(CSV_OUT, index=False, encoding="utf-8")
    print(f"saved → {CSV_OUT}")

    # validation cohort 제외
    valid = out.dropna(subset=["subtype", "os_event"]).copy()
    print(f"\nanalysis subset (cohort-wide NA 제외): {len(valid)}")

    print("\n=== 결정별 일치율 ===")
    per = concordance_per_decision(valid)
    print(per.to_string(index=False))

    print("\n=== Subtype × 결정 일치율 (%) ===")
    by_sub = concordance_by_subtype(valid)
    print(by_sub.to_string(index=False))

    # 전 결정 일치 (4/4 모두 일치한 환자)
    valid["actual_surgery"] = valid["surgery"].apply(normalize_actual_surgery)
    all_compared = (
        valid["rec_surgery"].notna() & valid["actual_surgery"].notna() &
        valid["rec_chemo"].notna()   & valid["chemo"].notna() &
        valid["rec_hormone"].notna() & valid["hormone"].notna() &
        valid["rec_radio"].notna()   & valid["radio"].notna()
    )
    all_match = all_compared & (
        (valid["rec_surgery"] == valid["actual_surgery"]) &
        (valid["rec_chemo"]   == valid["chemo"]) &
        (valid["rec_hormone"] == valid["hormone"]) &
        (valid["rec_radio"]   == valid["radio"])
    )
    print(f"\n4개 결정 모두 비교 가능: {all_compared.sum()}명")
    print(f"4개 결정 모두 일치     : {all_match.sum()}명 "
          f"({all_match.sum() / max(all_compared.sum(), 1) * 100:.1f}%)")


if __name__ == "__main__":
    main()
