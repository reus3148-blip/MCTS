"""METABRIC 임상 데이터를 연구용 patients.csv로 정리.

로드맵 1단계 표 컬럼 매핑:
    환자ID       → patient_id
    나이         → age
    폐경         → menopause           (Pre / Post)
    종양크기     → tumor_size_mm
    림프절 양성수→ lymph_pos
    병기         → stage               (0–4)
    ER/PR/HER2   → er / pr / her2      (1=Positive, 0=Negative)
    분자아형     → subtype             (HR+/HER2-, HR+/HER2+, HR-/HER2+, TNBC)
    PAM50        → pam50
    Grade        → grade               (1–3)
    Histology    → histology
    치료         → surgery / chemo / hormone / radio
    생존         → os_months / os_event
    재발         → rfs_months / rfs_event
    NPI          → npi                 (Nottingham Prognostic Index)

Ki67은 METABRIC에 없으므로 PAM50 + 3-gene 분류로 대체한다.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "brca_metabric"
OUT_DIR = ROOT / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def yn_to_int(s: pd.Series) -> pd.Series:
    return s.str.upper().map({"YES": 1, "NO": 0})


def pos_neg_to_int(s: pd.Series) -> pd.Series:
    cleaned = s.str.strip().str.lower().replace({"positve": "positive"})
    return cleaned.map({"positive": 1, "negative": 0})


def parse_event(s: pd.Series) -> pd.Series:
    return s.str.split(":").str[0].astype("float").astype("Int64")


def derive_subtype(er: pd.Series, pr: pd.Series, her2: pd.Series) -> pd.Series:
    hr = (er == 1) | (pr == 1)
    subtype = pd.Series(pd.NA, index=er.index, dtype="object")
    subtype[hr & (her2 == 0)] = "HR+/HER2-"
    subtype[hr & (her2 == 1)] = "HR+/HER2+"
    subtype[(~hr) & (her2 == 1)] = "HR-/HER2+"
    subtype[(~hr) & (her2 == 0)] = "TNBC"
    subtype[er.isna() | pr.isna() | her2.isna()] = pd.NA
    return subtype


def main() -> None:
    patients = pd.read_csv(RAW_DIR / "data_clinical_patient.txt", sep="\t", comment="#", low_memory=False)
    samples = pd.read_csv(RAW_DIR / "data_clinical_sample.txt", sep="\t", comment="#", low_memory=False)
    df = patients.merge(samples, on="PATIENT_ID", how="inner")

    out = pd.DataFrame({
        "patient_id":   df["PATIENT_ID"],
        "age":          df["AGE_AT_DIAGNOSIS"].round(1),
        "menopause":    df["INFERRED_MENOPAUSAL_STATE"],
        "tumor_size_mm":df["TUMOR_SIZE"],
        "lymph_pos":    df["LYMPH_NODES_EXAMINED_POSITIVE"].astype("Int64"),
        "stage":        df["TUMOR_STAGE"].astype("Int64"),
        "grade":        df["GRADE"].astype("Int64"),
        "er":           pos_neg_to_int(df["ER_STATUS"]).astype("Int64"),
        "pr":           pos_neg_to_int(df["PR_STATUS"]).astype("Int64"),
        "her2":         pos_neg_to_int(df["HER2_STATUS"]).astype("Int64"),
        "pam50":        df["CLAUDIN_SUBTYPE"],
        "histology":    df["HISTOLOGICAL_SUBTYPE"],
        "laterality":   df["LATERALITY"],
        "surgery":      df["BREAST_SURGERY"],
        "chemo":        yn_to_int(df["CHEMOTHERAPY"]).astype("Int64"),
        "hormone":      yn_to_int(df["HORMONE_THERAPY"]).astype("Int64"),
        "radio":        yn_to_int(df["RADIO_THERAPY"]).astype("Int64"),
        "os_months":    df["OS_MONTHS"].round(2),
        "os_event":     parse_event(df["OS_STATUS"]),
        "rfs_months":   df["RFS_MONTHS"].round(2),
        "rfs_event":    parse_event(df["RFS_STATUS"]),
        "npi":          df["NPI"],
    })
    out.insert(8, "subtype", derive_subtype(out["er"], out["pr"], out["her2"]))

    out_path = OUT_DIR / "patients.csv"
    out.to_csv(out_path, index=False, encoding="utf-8")

    print(f"saved → {out_path}")
    print(f"rows: {len(out)}  cols: {out.shape[1]}")
    print("\n[missing per column]")
    miss = out.isna().sum()
    pct = (miss / len(out) * 100).round(1)
    print(pd.concat([miss.rename("missing"), pct.rename("%")], axis=1).to_string())

    print("\n[subtype distribution]")
    print(out["subtype"].value_counts(dropna=False).to_string())

    print("\n[stage distribution]")
    print(out["stage"].value_counts(dropna=False).sort_index().to_string())

    print("\n[treatment combos - top 10]")
    combo_df = out[["surgery","chemo","hormone","radio"]].copy()
    for c in combo_df.columns:
        combo_df[c] = combo_df[c].map(lambda v: "NA" if pd.isna(v) else str(v))
    combo = combo_df.agg(" | ".join, axis=1).value_counts().head(10)
    print(combo.to_string())

    print("\n[head]")
    print(out.head(5).to_string())


if __name__ == "__main__":
    main()
