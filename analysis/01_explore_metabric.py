"""METABRIC 임상 데이터 1차 탐색.

cBioPortal datahub의 brca_metabric에서 받은 두 파일을 읽어
컬럼·결측·기초통계를 stdout으로 출력한다.

cBioPortal clinical 포맷 규칙:
- 1~4행: 메타데이터(#로 시작, 표시명/설명/타입/우선순위)
- 5행:   실제 컬럼명(PATIENT_ID, ...)
- 6행~:  데이터
"""

from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "brca_metabric"
PATIENT_FILE = DATA_DIR / "data_clinical_patient.txt"
SAMPLE_FILE = DATA_DIR / "data_clinical_sample.txt"


def load_clinical(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def summarize(name: str, df: pd.DataFrame) -> None:
    print(f"\n=== {name} ===")
    print(f"shape: {df.shape[0]} rows × {df.shape[1]} cols")
    print("\n[columns / dtype / non-null / unique / sample values]")
    for col in df.columns:
        s = df[col]
        non_null = s.notna().sum()
        n_unique = s.nunique(dropna=True)
        sample = s.dropna().head(3).tolist()
        print(f"  {col:<35} {str(s.dtype):<10} non-null={non_null:<5} unique={n_unique:<5} e.g.={sample}")

    num_cols = df.select_dtypes(include="number").columns.tolist()
    if num_cols:
        print("\n[numeric describe]")
        print(df[num_cols].describe().round(2).to_string())


def main() -> None:
    patients = load_clinical(PATIENT_FILE)
    samples = load_clinical(SAMPLE_FILE)
    summarize("PATIENTS", patients)
    summarize("SAMPLES", samples)

    print("\n=== JOIN CHECK ===")
    merged = patients.merge(samples, on="PATIENT_ID", how="outer", indicator=True)
    print(merged["_merge"].value_counts().to_string())


if __name__ == "__main__":
    main()
