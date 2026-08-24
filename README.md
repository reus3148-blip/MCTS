# MCTS-ONC

공개 유방암 코호트로 치료 의사결정 환경을 만들고, 단순화한 NCCN 정책과
Monte Carlo Tree Search(MCTS) 정책을 비교하는 학부 연구 프로젝트입니다.

> 현재 단계: **강건성·민감도 분석 v0.3 완료 (2026-07-20)**  
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

핵심 결과와 해석은 [동적 MCTS PoC v0.2 보고서](reports/dynamic-mcts-poc-v0.2/README.md),
[강건성 v0.3 보고서](reports/robustness-v0.3/README.md),
[민감도 v0.3 보고서](reports/sensitivity-v0.3/README.md),
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
  mcts/                   환경·생존모형·UCT 탐색 모듈
  dynamic/                공통 스키마·확률 전이·stochastic MCTS
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
  추정하지 않습니다.
- METABRIC은 오래된 치료 시기의 데이터이므로 최신 임상진료에 바로 적용할 수
  없습니다.
- v0.2는 동적으로 움직이지만 반응·독성·치료강도 차이는 합성 학습용 가정입니다.
  K-CURE 또는 임상시험 자료로 교체하기 전에는 실제 효과로 해석할 수 없습니다.
- NCCN 정책은 연구용 단순화 규칙이며 완전한 최신 가이드라인이 아닙니다.
