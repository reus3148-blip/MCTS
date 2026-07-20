"""Simplified NCCN-style policy used as policy A in this project.

This is an auditable research abstraction, not the complete or current NCCN
Breast Cancer Clinical Practice Guideline and not a clinical recommendation.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd


def nccn_surgery(state: pd.Series) -> Optional[str]:
    """Return the simplified surgery recommendation."""
    size = state.get("tumor_size_mm")
    stage = state.get("stage")
    if pd.isna(size) or pd.isna(stage):
        return None
    return "BCS" if size <= 30 and stage <= 2 else "MAST"


def nccn_chemo(state: pd.Series) -> Optional[int]:
    """Return the simplified adjuvant chemotherapy recommendation."""
    subtype = state.get("subtype")
    stage = state.get("stage")
    if pd.isna(subtype) or pd.isna(stage):
        return None
    if stage < 1:
        return 0
    if subtype in {"HR-/HER2+", "TNBC", "HR+/HER2+"}:
        return 1
    if subtype == "HR+/HER2-":
        lymph = state.get("lymph_pos")
        size = state.get("tumor_size_mm")
        grade = state.get("grade")
        high_risk = (
            (pd.notna(lymph) and lymph >= 1)
            or (pd.notna(size) and size > 20)
            or (pd.notna(grade) and grade == 3)
        )
        return int(high_risk)
    return None


def nccn_hormone(state: pd.Series) -> Optional[int]:
    """Recommend endocrine therapy when ER or PR is positive."""
    er = state.get("er")
    pr = state.get("pr")
    if pd.isna(er) and pd.isna(pr):
        return None
    er_positive = pd.notna(er) and er == 1
    pr_positive = pd.notna(pr) and pr == 1
    return int(er_positive or pr_positive)


def nccn_radio(
    state: pd.Series,
    recommended_surgery: Optional[str],
) -> Optional[int]:
    """Return the simplified radiotherapy recommendation."""
    if recommended_surgery is None:
        return None
    if recommended_surgery == "BCS":
        return 1
    lymph = state.get("lymph_pos")
    stage = state.get("stage")
    pmrt = (
        (pd.notna(lymph) and lymph >= 4)
        or (pd.notna(stage) and stage >= 3)
    )
    return int(pmrt)


def apply_nccn(row: pd.Series) -> pd.Series:
    """Apply all four policy-A decisions to one patient."""
    surgery = nccn_surgery(row)
    return pd.Series({
        "rec_surgery": surgery,
        "rec_chemo": nccn_chemo(row),
        "rec_hormone": nccn_hormone(row),
        "rec_radio": nccn_radio(row, surgery),
    })


def normalize_actual_surgery(value) -> Optional[str]:
    """Map METABRIC surgery labels to the policy action space."""
    if pd.isna(value):
        return None
    if value in {"BCS", "BREAST CONSERVING"}:
        return "BCS"
    if value in {"MAST", "MASTECTOMY"}:
        return "MAST"
    return None


def nccn_plan(row: pd.Series) -> Optional[tuple[object, ...]]:
    """Return a complete plan tuple, or None when any decision is unavailable."""
    recommendations = apply_nccn(row)
    if recommendations.isna().any():
        return None
    return (
        recommendations["rec_surgery"],
        int(recommendations["rec_chemo"]),
        int(recommendations["rec_hormone"]),
        int(recommendations["rec_radio"]),
    )


def concordance_per_decision(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize recommendation/actual agreement at each decision node."""
    rows = []
    pairs = [
        ("surgery", "rec_surgery", df["surgery"].apply(normalize_actual_surgery)),
        ("chemo", "rec_chemo", df["chemo"]),
        ("hormone", "rec_hormone", df["hormone"]),
        ("radio", "rec_radio", df["radio"]),
    ]
    for name, recommendation_column, actual in pairs:
        recommendation = df[recommendation_column]
        compared = recommendation.notna() & actual.notna()
        matched = compared & (recommendation == actual)
        rows.append({
            "decision": name,
            "compared_n": int(compared.sum()),
            "match_n": int(matched.sum()),
            "concordance_pct": round(
                matched.sum() / max(compared.sum(), 1) * 100,
                1,
            ),
            "missing_n": int((~compared).sum()),
        })
    return pd.DataFrame(rows)


def concordance_by_subtype(df: pd.DataFrame) -> pd.DataFrame:
    """Return subtype-by-decision agreement percentages."""
    working = df.copy()
    working["actual_surgery"] = working["surgery"].apply(
        normalize_actual_surgery
    )
    output = []
    for subtype, group in working.groupby("subtype"):
        row = {"subtype": subtype, "n": len(group)}
        for name, recommendation_column, actual_column in [
            ("surgery", "rec_surgery", "actual_surgery"),
            ("chemo", "rec_chemo", "chemo"),
            ("hormone", "rec_hormone", "hormone"),
            ("radio", "rec_radio", "radio"),
        ]:
            recommendation = group[recommendation_column]
            actual = group[actual_column]
            compared = recommendation.notna() & actual.notna()
            matched = compared & (recommendation == actual)
            row[name] = round(
                matched.sum() / max(compared.sum(), 1) * 100,
                1,
            )
        output.append(row)
    return pd.DataFrame(output)

