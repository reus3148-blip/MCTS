# Target Trial 프로토콜 — MCTS 치료전략의 인과효과 명세 (초안 v0.1)

- 작성일: 2026-08-24
- 상태: **방법론 초안**. 데이터 분석 전에 인과 질문을 먼저 고정하기 위한 문서.
- 적용 경계: 아래 명세는 "어떤 인과 질문에 답하려는가"를 규정할 뿐,
  현재 PoC(v0.1/v0.2)가 그 질문에 이미 답했다는 뜻이 아니다. 현재 결과는
  여전히 관찰 연관성이며 인과효과가 아니다.

---

## 왜 이 문서가 필요한가

지금까지의 결과("MCTS 정책이 NCCN보다 예측 5년 생존을 +3.7%p")는 **모델이 학습한
연관성을 최적화한 값**이다. "MCTS가 실제로 생존을 높인다"는 인과 주장으로 넘어가려면,
Hernán & Robins의 **target trial emulation** 틀에 따라 (1) 만약 무작위배정 임상시험을
한다면 그 시험이 어떤 모습일지 먼저 명세하고, (2) 관찰자료(METABRIC → K-CURE)로 그
시험을 "모방(emulate)"해야 한다. 프로토콜 없이 회귀만 돌리면 immortal-time bias,
confounding by indication, 시간의존 교란에 그대로 노출된다.

---

## 1. 이상적 표적시험 (Target Trial)

### 1.1 적격성 기준 (Eligibility)

- 조직학적으로 확진된 **침습성 유방암**, 원격전이 없음 (M0)
- 진단 시 **치료 미시작** 상태
- 4개 골격 결정(수술·항암·내분비·방사선)에 대한 임상 정보가 time zero에 존재
- 제외: 양측성 동시 진단, 이전 유방암력, 진단 전 사망/추적불가

> METABRIC 에뮬레이션에서의 대응: `prepare_model_cohort`의 완전 케이스에
> 병기·수용체·4개 실제 치료가 모두 있는 환자. 현재 코호트 1,955명(OS).

### 1.2 치료 전략 (Treatment Strategies)

핵심은 **단일 시점 치료가 아니라 순차적 치료전략(dynamic treatment regime)의 비교**라는
점이다. 두 전략을 비교한다.

- **전략 A — 가이드라인 정책**: 각 결정 노드에서 단순화 NCCN 규칙이 지정한 치료를 시행.
- **전략 B — 모델 유도 정책(MCTS)**: 각 결정 노드에서, 환자의 관찰된 중간 상태(반응·종양크기·독성)에
  조건부로 MCTS가 선택한 치료를 시행.

두 전략 모두 **동적 규칙**이다. 즉 "처음에 무엇을 배정했나"가 아니라 "매 결정 시점의
상태에 따라 규칙 f(state)→action을 따랐나"로 정의된다. 이 때문에 분석은 시점별 규칙
준수를 다루는 **g-방법(g-formula, IPW of marginal structural model, doubly robust)** 이
필요하며, 단순 baseline 보정 회귀로는 부족하다.

### 1.3 배정 (Assignment)

- 이상시험: 적격 환자를 전략 A/B에 **1:1 무작위배정**, 배정 은폐.
- 관찰 에뮬레이션: 무작위배정이 없으므로 time zero의 측정된 교란요인으로 **교환가능성(exchangeability)** 을
  근사. 미측정 교란은 민감도 분석으로 다룬다(§4.4).

### 1.4 Time Zero

- **진단 후 첫 치료 결정 시점**을 time zero로 고정한다.
- 규칙: 적격성·전략배정·추적 시작이 **동일 시점**에 정렬되어야 한다.
  이 정렬을 어기면 immortal-time bias가 발생한다(예: "항암을 끝까지 받은 사람"으로
  전략을 정의하면, 그 정의 자체가 초기 생존자를 선택함).

### 1.5 추적 기간 (Follow-up)

- time zero부터 **5년(60개월)**, 또는 사망·마지막 추적 중 먼저 도래하는 시점까지.
- 검열: 행정적 검열(자료 마감), 추적소실.

### 1.6 결과 (Outcomes)

| 구분 | 결과 | 정의 |
|---|---|---|
| 1차 | 5년 전체생존(OS) | time zero부터 전원인 사망까지 |
| 2차 | 5년 무재발생존(RFS) | 재발 또는 사망까지 |
| 2차 | 급성 치료독성 | 사전정의 이상반응(등급 기준 확정 필요) |
| 탐색 | 복합 utility | 생존·무재발·독성을 사전 가중한 값 (가중치 사전등록) |

### 1.7 인과 대비 (Causal Contrast) / 추정치 (Estimand)

- **Intention-to-treat 유사**: 전략 배정에 따른 효과.
- **Per-protocol**: 전략을 끝까지 준수했을 때의 효과 (시간의존 교란 보정 필요).
- 주 추정치: **5년 위험 차이(risk difference)** 와 **위험비**, 전략 B vs A.
  예) "전략 B를 따를 때의 5년 사망위험 − 전략 A를 따를 때의 5년 사망위험."

---

## 2. 관찰자료로의 에뮬레이션 (Emulation)

| 표적시험 요소 | METABRIC (현재) | K-CURE (향후) |
|---|---|---|
| 적격성 | 완전 케이스 코호트 | 진단 코호트에서 동일 기준 적용 |
| Time zero | 진단/첫 결정(근사) | 실제 첫 치료 결정일 |
| 전략 A/B | NCCN 규칙 / MCTS 규칙 | 동일 |
| 교환가능성 | baseline 공변량 보정 | 시간의존 교란까지 g-방법 |
| 결과 | OS/RFS Cox | 실제 재발·사망일 기반 |
| 준수 | 관찰 안 됨(단일시점) | 치료 순서·용량으로 준수 정의 |

**핵심 한계(METABRIC):** METABRIC은 시점별 중간상태(반응·독성)와 치료 타이밍을
충분히 담지 않아, 동적 전략의 per-protocol 효과를 제대로 에뮬레이트하기 어렵다.
따라서 v0.2의 동적 환경은 **합성 전이로 채운 시뮬레이션**이고, 진짜 에뮬레이션은
시간해상도가 있는 K-CURE에서 가능해진다.

---

## 3. 분석 계획

1. **기술통계**: 전략별 baseline 공변량, positivity 점검(각 공변량 층에서 두 전략이
   모두 관찰되는가).
2. **교란 보정**: time zero 공변량으로 propensity/IPW; 시간의존 교란은 marginal
   structural model 또는 g-formula.
3. **Doubly robust 추정**: 결과모형과 처치모형 중 하나만 옳아도 일치추정.
4. **생존 분석**: IPW 가중 Kaplan-Meier, 5년 위험차/위험비 + bootstrap 신뢰구간.
5. **예측 vs 인과 분리**: Cox 보상모형은 **예측·시뮬레이션 도구**로만 쓰고, 인과
   주장에는 위 g-방법 추정치를 사용한다.

---

## 4. 편향과 민감도

- **4.1 Confounding by indication**: 더 위험한 환자가 더 강한 치료를 받음 →
  IPW/g-방법 + 미측정 교란 민감도(E-value).
- **4.2 Immortal-time bias**: 전략을 "완주자"로 정의하지 않기. 시점정렬 엄수.
- **4.3 Positivity 위반**: 특정 상태에서 한 전략만 관찰되면 추정 불안정 → 층 병합 또는
  추정 제한.
- **4.4 미측정 교란**: E-value, 음성대조(negative control) 결과.
- **4.5 모델 오지정**: doubly robust + 다중 명세 비교.

---

## 5. 사전등록해야 할 항목 (분석 전 고정)

- 적격성·time zero·전략 정의의 정확한 코드
- 결과·독성 등급 기준
- utility 가중치와 그 출처
- 주 추정치와 신뢰구간 방법
- 하위군(분자아형별)과 민감도 목록

---

## 6. 다음 실행 단계

1. 위 적격성·time zero를 METABRIC 변수로 코드화(`analysis/`에 target-trial 코호트 정의 추가).
2. baseline propensity/IPW 모듈 시제품 구현 및 positivity 진단.
3. K-CURE 변수사전 확보 후 시간의존 교란·per-protocol 준수 정의(→ `docs/k-cure-variable-dictionary.md`).
4. 임상 자문: eligibility·utility·독성 정의 검토.

---

## 참고문헌

- Hernán MA, Robins JM. Using Big Data to Emulate a Target Trial When a Randomized
  Trial Is Not Available. *Am J Epidemiol*. 2016. https://doi.org/10.1093/aje/kwv254
- Hernán MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC, 2020.
- Robins JM. A new approach to causal inference in mortality studies. *Math Modelling*. 1986.
- VanderWeele TJ, Ding P. Sensitivity Analysis: the E-value. *Ann Intern Med*. 2017.
- Collins GS, et al. TRIPOD+AI statement. *BMJ*. 2024. https://doi.org/10.1136/bmj-2023-078378
