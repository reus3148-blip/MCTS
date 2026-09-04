# MCTS-ONC

공개 유방암 코호트로 치료 의사결정 환경을 만들고, 단순화한 NCCN 정책과
Monte Carlo Tree Search(MCTS) 정책을 비교하는 학부 연구 프로젝트입니다.

> 현재 단계: **환자 수 민감도 v1.1 완료 (2026-09-04)**  
> 연구용 예측 실험이며 실제 환자의 치료 권고 도구가 아닙니다.

## 지금까지 한 일

1. METABRIC 2,509명을 공통 임상 변수 23개로 전처리했습니다.
2. 분석 가능한 1,980명에서 분포, Kaplan-Meier 곡선, NCCN 정책 일치율을
   계산했습니다.
3. 수술·항암·호르몬·방사선의 4단계, 환자별 8~16개 치료 경로 환경을 구현했습니다.
4. 규제 Cox 생존모형을 5년 전체생존 보상으로 사용해 UCT-MCTS를 실행했습니다.
5. 고정 테스트셋 390명에서 MCTS를 완전탐색과 검증하고 NCCN과 비교했습니다.
6. OS/RFS를 사용해 반응·독성·재발 상태가 변하는 5년 확률환경을 만들었습니다.
7. 환자별 의사결정 궤적을 18~135개로 확장하고 동적 MCTS를 8,000회 평가했습니다.
8. **v0.3**: 20개 시드로 강건성을 검증해 MCTS 우세가 작고(생존 +1.9%p) 개별 결정은
   시드에 민감함을(첫 결정 평균 63% 일치) 확인했습니다.
9. **v0.3**: 합성 가정 민감도 분석으로 결론을 가장 크게 좌우하는 요인이 임상 효과가
   아니라 **보상 가중치(가치판단)** 임을 밝히고, target trial 인과 명세와 K-CURE 변수
   요청 명세를 초안했습니다.
10. **v0.4**: 결정 불안정성의 원인을 탐색 부족과 진짜 동률로 분리했습니다. 잡음은
    이론값대로 1/√N 로 줄었고(기울기 −0.508) 예산을 2048로 올리자 일치율이
    57.8% → 80.4%로 올랐지만, 49개 결정 지점 중 10개는 **어떤 예산에서도** 갈리지
    않았고 그 지점들은 호르몬치료·방사선에 몰려 있었습니다.
11. **v0.5**: MCTS·MDP 구현을 형식적 정의에 대조해 감사했습니다. 알고리즘은 맞았으나
    환경에서 **선행치료를 고르기만 하면 재발 위험이 6.5~8.6% 깎이는 미선언 이득**을
    찾아 중립화했습니다. 할인율도 설정에 명시해 민감도 대상으로 만들었습니다.
12. **v0.6**: 관찰 연관성을 인과추정으로 바꾸는 첫 시도로 IPW 표적시험 에뮬레이션을
    돌렸습니다. 항암 시행 여부가 사실상 결정돼 있어 **positivity가 무너졌고**(가중 후
    최악 |SMD| 1.13 → 2.19, 유효표본 24명), 겹침 구간 38%로 좁혀야 균형이 잡혔습니다.
    그 안에서 조정 전후로 **효과의 부호가 뒤집혔습니다**(+0.057 → −0.036, CI는 0 포함).
13. **v0.7**: 서로 다른 모형에 기대는 세 추정량(IPW·g-계산·AIPW)이 CI 폭의 3% 안에서
    일치해 모형 오설정이 아님을 확인했습니다. 결정별 식별 가능성 지도를 그리다가
    **방사선치료의 넓은 겹침이 교란요인(수술 유형) 누락의 산물**임을 발견했습니다 —
    넣자마자 유지율 96%→69%, 균형 실패.
14. **v0.8**: v0.5의 "편향 영향은 검출 한계보다 43배 작다"는 예측을 v0.3·v0.4 재실행으로
    확인했습니다(효용 격차 소수점 넷째 자리까지 동일). 이어서 2차원 격자에서 **두 가치판단이
    서로의 방향을 뒤집는 것처럼** 보였고, 할인율(0~5%)은 결론을 바꾸지 않음을 확인했습니다.
15. **v0.9**: 그 방향 반전을 시드 12개로 확인했더니 **살아남지 못했습니다** — 상호작용이
    절반(−0.0115)이 되고 z가 잡음 기준 아래로 내려갔으며 부호 반전 자체가 사라졌습니다.
    승자의 저주 사례이며, v0.8의 해당 서술에 정정 표기를 달았습니다.
16. **v1.0**: v0.3 민감도를 시드 12개로 재실행했습니다. **시드만 4배 늘렸는데 기준선이
    0.019 움직여** 순위를 매기려던 효과들보다 컸습니다. 핵심 결론(상위 2개가 가치판단)은
    유지되나 **예산 1024에서만** 성립하고, v0.3이 보고한 **부호 반전은 전부 철회**됐습니다.
17. **v1.1**: 마지막 병목이던 환자 수를 8 → 20명으로 늘렸습니다. **헤드라인은 +0.0004만
    움직여 처음으로 버텼지만**, 환자 부트스트랩에서 "상위 2개가 가치판단"의 재현율이
    **8명 24.6% → 20명 69.5%** 였습니다. 결론은 맞았어도 **8명은 그 근거가 못 됐습니다.**
    임상 파라미터들의 영향은 70~86% 사라졌고, |z|≥2인 파라미터가 처음으로 2개가 되면서
    그 둘이 정확히 두 가치판단이었습니다.

전체 줄기를 한 번에 읽으려면 [연구 이야기](docs/research-story.md),
숫자가 왜 달라졌는지는 [숫자 화해](docs/results-reconciliation.md),
제안서 대비 범위 변화는 [제안서 대비 산출물](docs/proposal-vs-delivered.md)을 보세요.
셋 다 사이트의 [/story](https://mcts.blundermate.app/story)에도 있습니다.

핵심 결과와 해석은 [동적 MCTS PoC v0.2 보고서](reports/dynamic-mcts-poc-v0.2/README.md),
[강건성 v0.3 보고서](reports/robustness-v0.3/README.md),
[민감도 v0.3 보고서](reports/sensitivity-v0.3/README.md),
[예산 스케일링 v0.4 보고서](reports/budget-scaling-v0.4/README.md),
[환경 편향 수정 v0.5 보고서](reports/environment-fix-v0.5/README.md),
[IPW 표적시험 v0.6 보고서](reports/ipw-target-trial-v0.6/README.md),
[이중강건 v0.7 보고서](reports/doubly-robust-v0.7/README.md),
[상호작용 민감도 v0.8 보고서](reports/interaction-sensitivity-v0.8/README.md),
[확인 실험 v0.9 보고서](reports/utility-interaction-v0.9/README.md),
[민감도 재실행 v1.0 보고서](reports/sensitivity-precision-v1.0/README.md),
[환자 수 민감도 v1.1 보고서](reports/sensitivity-patients-v1.1/README.md),
[v0.5 환경 재실행](reports/robustness-v0.5env/README.md),
[Target trial 프로토콜](docs/target-trial-protocol.md),
날짜별 작업 근거는 [프로젝트 타임라인](PROJECT_TIMELINE.md)에 정리되어 있습니다.

## 저장소 구조

```text
analysis/
  01~05_*.py              데이터·시각화·NCCN 분석 실행 파일
  06_run_mcts_poc.py      보상모형 학습과 MCTS 정책 비교
  07_visualize_mcts_poc.py 결과 그림 생성 및 사이트 반영
  08_run_dynamic_mcts_poc.py OS/RFS 기반 동적 정책 실험
  09_visualize_dynamic_mcts_poc.py 동적환경 Figure 14~17
  10_run_multiseed_robustness.py  다중 시드 강건성 (v0.3)
  11_run_sensitivity_analysis.py  합성 가정 민감도 (v0.3)
  12_run_budget_scaling.py        탐색 예산 스케일링 진단 (v0.4)
  13_run_environment_fix_impact.py 환경 편향 수정의 영향 측정 (v0.5)
  14_visualize_v03_v05.py         v0.3~v0.5 Figure 18~23
  15_run_ipw_target_trial.py      IPW 표적시험 에뮬레이션 (v0.6)
  16_visualize_ipw.py             IPW Figure 24~26
  17_run_doubly_robust.py         AIPW·IPCW·결정별 식별 가능성 지도 (v0.7)
  18_visualize_doubly_robust.py   Figure 27~28
  19_run_interaction_sensitivity.py 2차원 상호작용 민감도 (v0.8)
  20_visualize_interaction.py     Figure 29
  21_visualize_env_refresh.py     Figure 30 (환경 재실행 대조)
  22_run_utility_interaction_confirm.py 가치판단 상호작용 확인 (v0.9)
  23_visualize_interaction_confirm.py   Figure 31
  24_run_sensitivity_precision.py 민감도 정밀도 재실행 (v1.0)
  25_visualize_sensitivity_precision.py Figure 32
  26_run_sensitivity_patients.py  환자 수를 늘린 민감도 (v1.1)
  27_visualize_sensitivity_patients.py Figure 33
  causal/ipw.py                   프로펜서티·균형·가중 KM·IPCW·AIPW·E-value (numpy 구현)
  causal/decisions.py             결정별 코호트·교란요인 명세·트리밍
  mcts/                   환경·생존모형·UCT 탐색 모듈
  dynamic/                공통 스키마·확률 전이·stochastic MCTS
  dynamic/cohort.py       10~12가 공유하는 코호트·보상모형·매니페스트
data/
  brca_metabric/          원본 데이터 (Git 제외)
  processed/              전처리 데이터 (Git 제외)
reports/mcts-poc-v1/
  tables/                 고정 분할, 계수, 환자별 결정, 요약표
  figures/                보고서 그림
  metrics.json            모형·탐색 핵심 지표
  run_manifest.json       입력 해시·시드·실행 정보
reports/dynamic-mcts-poc-v0.2/
  tables/                 8,000 episode·trace·정책·안정성 결과
  assumptions_snapshot.json 실행 당시 합성 가정
reports/budget-scaling-v0.4/
  tables/                 결정지점 × 예산 수렴, 예산별 효용 격차
reports/environment-fix-v0.5/
  tables/                 편향/수정 환경의 시드별 결과와 짝지은 차이
reports/ipw-target-trial-v0.6/
  tables/                 프로펜서티·공변량 균형·효과추정·트리밍 민감도
reports/doubly-robust-v0.7/
  tables/                 추정량 비교·결정별 식별 가능성 지도·검열 가중
reports/interaction-sensitivity-v0.8/
  tables/                 2차원 격자 27셀·상호작용 통계
reports/utility-interaction-v0.9/
  tables/                 시드 12개 확인 격자·시드별 차이-의-차이
reports/sensitivity-precision-v1.0/
  tables/                 13변형 × 2예산 결과·짝지은 표준오차
reports/sensitivity-patients-v1.1/
  tables/                 13변형 × 3코호트 결과·환자별 효용 격차
reports/*-v0.5env/        v0.3·v0.4를 수정 환경에서 재실행한 결과
configs/dynamic_v0_5.json 현재 환경 가정 (v0.2는 과거 리포트 재현용으로 보존)
docs/k-cure-adaptation.md K-CURE 이행 데이터 계약
src/minutes/              날짜별 연구 기록과 의사결정 근거
tests/                    정책·환경·MCTS·특징변환 자동 검증
```

## 재현 방법

Python 분석:

```powershell
py -m pip install -r requirements.txt
py analysis\04_nccn_policy.py
py analysis\06_run_mcts_poc.py
py analysis\07_visualize_mcts_poc.py
py analysis\08_run_dynamic_mcts_poc.py
py analysis\09_visualize_dynamic_mcts_poc.py
py analysis\10_run_multiseed_robustness.py
py analysis\11_run_sensitivity_analysis.py
py -m unittest discover -s tests -v
```

연구 블로그:

```powershell
npm install
npm run dev
```

## 해석 경계

- Cox 보상모형은 관찰자료의 **연관성**을 학습합니다. 치료의 인과효과를
  추정하지 않습니다. v0.6에서 같은 자료의 조정 전 관찰 비교가 **부호까지 틀린다**는
  것을 직접 확인했습니다.
- v0.6~v0.7의 인과 추정치는 **단일 시점 결정 하나**에 대한 것이며, 표적시험 프로토콜이
  명세한 MCTS vs NCCN 전략 비교가 아닙니다. 그 비교는 치료 시점 정보가 있어야 식별됩니다.
- **겹침(positivity) 진단은 교란요인 목록에 조건부입니다.** v0.7에서 방사선치료의 넓은
  겹침이 수술 유형 누락의 산물이었음을 확인했습니다. 겹침이 넓다고 안심하면 안 됩니다.
- METABRIC은 오래된 치료 시기의 데이터이므로 최신 임상진료에 바로 적용할 수
  없습니다.
- v0.2는 동적으로 움직이지만 반응·독성·치료강도 차이는 합성 학습용 가정입니다.
  K-CURE 또는 임상시험 자료로 교체하기 전에는 실제 효과로 해석할 수 없습니다.
- NCCN 정책은 연구용 단순화 규칙이며 완전한 최신 가이드라인이 아닙니다.
- v0.4 이후 기본 탐색 예산은 **1024 이상**입니다. 256으로 낸 v0.2·v0.3의 결정 단위
  수치는 탐색 해상도에 미달했으므로 그 단서와 함께 읽어야 합니다.
- 호르몬치료·방사선 결정은 어떤 예산에서도 안정되지 않습니다. 이 결정에 대해서는
  "MCTS가 X를 권한다"가 아니라 **"우리 효용 아래 두 선택이 동률"** 로만 말할 수 있습니다.
- **v0.2~v0.4의 수치는 응답 채널 편향이 있던 환경에서 나왔습니다.** v0.5에서 고쳤고,
  2026-08-28 재실행으로 결론이 바뀌지 않음을 확인했습니다(`reports/*-v0.5env`).
  인용은 재실행 쪽 값을 씁니다.
- **보상 가중치는 각각 결과를 크게 움직입니다.** 두 값을 함께 사전등록해야 합니다.
  (v0.8이 보고한 "서로의 방향을 뒤집는다"는 v0.9 확인 실험에서 기각됐습니다.)
- **민감도 분석의 최소 설계는 시드 12개·예산 1024·환자 20명입니다**(v1.0, v1.1).
  3시드 결과는 인용하지 않습니다 — v0.3의 기준선(−0.0043)과 "부호가 바뀐다"는 서술은
  철회됐습니다. **8명으로 낸 순위도 인용하지 않습니다** — 환자를 다시 뽑으면 우리 결론이
  24.6%만 재현됩니다(v1.1).
- **민감도 순위는 순서가 아니라 그룹으로만 말할 수 있습니다.** 환자 20명에서 상위 2개가
  가치판단이고 |z|≥2도 그 둘뿐이지만, 1·2위의 CI가 겹쳐 순서는 구분되지 않습니다.
- **효용 격차의 평균만 인용하지 않습니다.** 환자 20명에서 부호는 전원 같지만 크기가
  87배(+0.0009 ~ +0.0809) 갈리고 상위 3명이 전체의 33%를 가져갑니다(v1.1).
  "평균 +0.031"은 전형적인 환자의 값이 아닙니다.
