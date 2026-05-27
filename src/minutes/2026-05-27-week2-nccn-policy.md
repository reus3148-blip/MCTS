---
date: 2026-05-27
title: 2주차 — NCCN 정책 A 1차 구현 + 일치율 분석
summary: NCCN 침습성 유방암 가이드라인을 4개 결정 노드의 if-then 함수로 코드화하고 METABRIC 1,980명에 적용. 4개 결정 모두 NCCN과 일치한 환자는 11.3%에 그쳤다.
---

# 2주차 — NCCN 정책 A 1차 구현 + 일치율 분석

> **목적**: [[2026-05-26-kickoff]] 의 결정 3에서 정의한 "정책 A (NCCN)"를 코드로 구현하고, METABRIC 환자에게 적용해 실제 받은 치료와의 일치율(concordance)을 산출한다.
> **결과물**: [analysis/04_nccn_policy.py](https://github.com/reus3148-blip/genie/blob/main/analysis/04_nccn_policy.py), [analysis/05_visualize_nccn.py](https://github.com/reus3148-blip/MCTS/blob/main/analysis/05_visualize_nccn.py), `patients_with_nccn.csv`

---

## 이번 주에 실제로 한 일

### ① 룰 4개를 사람이 먼저 합의
의학적 정확성이 중요한 부분이라 코딩 전에 룰 초안을 사람 검토로 통과시켰다. NCCN 침습성 유방암 가이드라인의 핵심을 [[2026-05-26-kickoff]] 의 결정 3·4 에서 정한 scope(시기 무관 골격 4개 결정)에 맞춰 단순화했다.

| 결정 | 룰 |
|---|---|
| 수술 | 종양 ≤ 30mm AND stage ≤ 2 → BCS,  else Mastectomy |
| 항암 | HER2+ stage≥1 OR TNBC stage≥1 OR HR+/HER2- 고위험* → Chemo |
| 호르몬 | ER+ 또는 PR+ → Hormone |
| 방사선 | BCS → Radio,  Mastectomy AND (림프절≥4 OR stage≥3) → Radio |

* HR+/HER2- 고위험 = 림프절 양성 OR 종양 > 20mm OR grade = 3

### ② 함수로 코드화 — [04_nccn_policy.py](https://github.com/reus3148-blip/MCTS/blob/main/analysis/04_nccn_policy.py)

각 결정을 독립 함수로 분리. 결정에 필요한 변수가 결측이면 `None` 반환해 비교에서 제외되도록 했다.

```python
def nccn_surgery(state)  -> Optional[str]   # 'BCS' or 'MAST'
def nccn_chemo(state)    -> Optional[int]   # 0 or 1
def nccn_hormone(state)  -> Optional[int]   # 0 or 1
def nccn_radio(state, recommended_surgery) -> Optional[int]
```

방사선은 수술 결정 결과에 의존하므로 명시적으로 인자로 받게 했다 — **MDP 의 결정 순서 일부가 함수 시그니처에 자연스럽게 드러난다.**

### ③ METABRIC 전체에 적용
1,980명 분석 대상에 룰을 일괄 적용해 권고 치료 컬럼 4개(`rec_surgery`, `rec_chemo`, `rec_hormone`, `rec_radio`)를 추가하고 `patients_with_nccn.csv` 로 저장.

### ④ 일치율 계산
권고 vs 실제 비교를 결정별·subtype별로 산출. 결측이 있는 결정은 비교에서 제외.

---

## 결과

### 결정별 NCCN 일치율
![결정별 NCCN 일치율](/figures/fig08_nccn_concordance.png)

| 결정 | 비교 가능 | 일치 | 일치율 |
|---|---|---|---|
| 호르몬 | 1,980 | 1,406 | **71.0%** |
| 방사선 | 1,459 | 933 | 63.9% |
| 수술 | 1,445 | 824 | 57.0% |
| **항암** | **1,466** | **583** | **39.8%** |

**4개 결정 전부 비교 가능: 1,445명 / 4개 모두 일치: 163명 (11.3%)**

→ METABRIC 환자 10명 중 9명은 현재 NCCN 기준으로 보면 최소 1개 결정에서 다른 치료를 받았다. [[2026-05-27-week1-visualization]] 의 "치료조합 top 4가 모두 항암 없음" 관찰과 정확히 같은 신호 — **항암 일치율이 40% 미만**으로 결정 4개 중 가장 낮다.

### 분자아형 × 결정 일치율
![분자아형별 NCCN 일치율](/figures/fig09_nccn_concordance_by_subtype.png)

가장 흥미로운 칸은 **HR+/HER2+ 환자의 보조항암 일치율 25%**. 이 그룹에 NCCN은 항암을 강하게 권고하지만 METABRIC 시대(1977-2005) 영국 임상은 HER2 양성을 일관되게 인식하지 못해 항암을 적게 적용한 것으로 보인다.

| Subtype | n | 수술 | 항암 | 호르몬 | 방사선 |
|---|---|---|---|---|---|
| HR+/HER2- | 1,413 | 57% | **36%** | 71% | 62% |
| HR+/HER2+ | 113 | 58% | **25%** ⚠️ | 74% | 70% |
| HR-/HER2+ | 134 | 55% | **61%** | 73% | 60% |
| TNBC | 320 | 58% | **52%** | 70% | 71% |

호르몬·방사선·수술은 subtype 간 변동이 작다(±3%p 내외). 항암만 분자아형에 따라 25%~61%로 크게 흔들린다. **이게 본 연구에서 MCTS가 가장 큰 차이를 만들 수 있는 결정 노드**라는 신호.

---

## 정리

| 검증 항목 | 결과 |
|---|---|
| 데이터 시기 mismatch ([[2026-05-27-week1-visualization]] 결정 4) | ✅ HR+/HER2+ 항암 일치 25%로 추가 확정 |
| 항암이 가장 mismatch 큰 결정 | ✅ 40% (다른 결정은 57~71%) |
| 호르몬은 시대 무관 안정적 | ✅ 71% — 모든 subtype에서 70%+ |
| 가이드라인-MCTS 비교가 의미 있는 결정 노드 | **항암** (가장 큰 차이) > 방사선 > 수술 > 호르몬 |

이번 주 결과는 [[2026-05-26-kickoff]] 의 결정 4("프레임을 가이드라인 검증 PoC로")가 정량적으로 정당화된다는 것 외에, **연구의 다음 우선순위가 자연스럽게 항암 결정 노드라는 것**을 보여준다.

---

## 한계 명시

룰 1차 버전이 단순화되어 있다는 점을 분명히 한다.

- **수술**: 환자 선호·BRCA 보유 여부·다발성 종양 등은 미반영. METABRIC에 변수 없음.
- **항암**: HR+/HER2- 고위험 정의가 Ki-67 부재로 단순함. PAM50 LumB 활용 시 더 정밀해질 수 있음.
- **호르몬**: 폐경 상태별 약제 선택(Tamoxifen vs AI) 미구분. 본 연구 scope에서 "호르몬 Y/N" 만 다룸.
- **방사선**: PMRT 권고 기준은 시대에 따라 변동이 컸다. 본 룰은 현재 NCCN 단순화.

이 한계들은 [[2026-05-26-kickoff]] 결정 4의 *"PoC scope 정의"* 와 정합. 후속 연구에서 다듬는다.

---

## 다음 작업 (3주차)

1. **불일치 패턴 심층 분석** — 4개 결정 중 어떤 결정의 불일치가 OS/RFS 결과와 관련 있는가? (Cox regression)
2. **MCTS 환경 골격 정의** — state space, action space, transition probability 추정 1차 버전.
3. **단순 환경(48 상태 × 2 행동)에서 MCTS 시뮬레이션** — NCCN 정책과 동일한 결정을 내리는지부터 확인.

---

*기록: Claude — 회의록 자동 생성 후 사람이 검토함.*
