# METABRIC에서 K-CURE로 옮기는 계획

## 목표

METABRIC에서 코드와 연구 질문을 먼저 검증하고, K-CURE 접근이 가능해지면
같은 동적 환경의 합성 전이를 실제 한국인 진료 데이터 추정치로 교체합니다.

## 공통 환자 입력 계약

현재 환경이 요구하는 최소 필드입니다.

| 표준 필드 | 의미 | METABRIC | K-CURE 준비 |
|---|---|---|---|
| patient_id | 비식별 환자키 | patient_id | 비식별 연계키 매핑 필요 |
| age | 진단 나이 | age | 진단일·출생연도 계산 |
| menopause | 폐경 상태 | menopause | 직접 변수 또는 나이 기반 정의 필요 |
| tumor_size_mm | 종양 크기 | tumor_size_mm | 병리/영상 기준 우선순위 정의 |
| lymph_pos | 양성 림프절 수 | lymph_pos | 수술 병리 자료 매핑 |
| stage | 병기 | stage | AJCC 판본과 진단시점 고정 |
| grade | 조직학적 등급 | grade | 병리 코드 매핑 |
| subtype | HR/HER2 아형 | 파생 | ER·PR·HER2로 동일 파생 |
| er / pr / her2 | 생체표지자 | 0/1 | 검사일과 양성 기준 통일 |

## 동적 전이를 위해 추가로 필요한 K-CURE 변수

| 범주 | 필요한 변수 | 이유 |
|---|---|---|
| 시간 기준 | 진단일, 치료 시작·종료일 | 치료 순서와 time zero 정의 |
| 수술 | 수술일, 수술 종류, 절제 범위 | 수술 전후 상태 구분 |
| 약물 | 약제, 용량, cycle, 투여일 | 표준/강화 같은 합성 action 교체 |
| 반응 | 영상·병리 반응, 평가일 | transition probability 직접 추정 |
| 방사선 | 시작일, 범위, 선량 | local/regional action 정의 |
| 독성 | 이상반응 코드, 입원, 중단일 | 합성 독성 확률 교체 |
| 결과 | 재발일, 사망일, 마지막 추적일 | OS/RFS와 검열 정의 |
| 교란 | 동반질환, 수행상태, 병원·연도 | 치료 선택 편향 보정 |

## 옮기는 순서

1. K-CURE 변수사전에서 위 필드의 원천 테이블과 코드를 찾습니다.
2. `PatientProfile`로 바꾸는 `k_cure_adapter.py`를 작성합니다.
3. 진단일을 time zero로 고정하고 치료 episode를 시간순으로 만듭니다.
4. 현재 JSON의 반응·독성 확률을 K-CURE 추정치로 하나씩 교체합니다.
5. 병원·연도·환자군별 positivity와 결측을 확인합니다.
6. target trial을 명세하고 인과추론 분석과 예측 분석을 구분합니다.
7. METABRIC에서 정한 코드는 동결하고 K-CURE 외부 검증 결과를 비교합니다.

## 바꾸지 않아도 되는 부분

- `DynamicState`의 상태 전이 인터페이스
- 환자별 action mask 구조
- chance-aware UCT-MCTS 탐색기
- episode 평가와 seed 재현 체계
- Minutes와 report 산출 방식

## 반드시 다시 검토할 부분

- 병기·수용체 정의와 코드북
- eligibility 규칙
- utility 가중치
- 시간의존 교란과 immortal-time bias
- 병원별 치료 패턴 차이
- 외부검증 및 임상 자문

