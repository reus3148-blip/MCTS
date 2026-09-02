# MCTS-ONC 프로젝트 타임라인

나중에 제3자가 연구의 순서와 근거 파일을 빠르게 확인할 수 있도록 날짜별
산출물을 연결한 인덱스입니다. 세부 의사결정은 각 Minutes 원문에 남깁니다.

| 날짜 | 이정표 | 실제 산출물 | 상태 |
|---|---|---|---|
| 2026-05-26 | 연구 킥오프와 범위 확정 | `src/minutes/2026-05-26-kickoff.md` | 완료 |
| 2026-05-27 | METABRIC 전처리·기초 시각화 | `patients.csv`, Figure 01~07, `analysis/02~03` | 완료 |
| 2026-05-27 | 단순화 NCCN 정책 A와 일치율 | `patients_with_nccn.csv`, Figure 08~09, `analysis/04~05` | 완료 |
| 2026-05-27 | 연구 블로그와 Minutes 체계 | `/minutes`, 날짜별 Markdown 기록 | 완료 |
| 2026-07-11 | Cox 보상모형 + UCT-MCTS PoC v0.1 | `analysis/mcts`, `analysis/06~07`, Figure 10~13 | 완료 |
| 2026-07-11 | 확률적 동적환경 + MCTS PoC v0.2 | `analysis/dynamic`, `analysis/08~09`, Figure 14~17 | 완료 |
| 2026-07-11 | K-CURE 이행 데이터 계약 | `docs/k-cure-adaptation.md` | 초안 완료 |
| 2026-08-24 | 다중 시드 강건성 분석 v0.3 | `analysis/10`, `reports/robustness-v0.3` | 완료 |
| 2026-08-24 | 합성 가정 민감도 분석 v0.3 | `analysis/11`, `reports/sensitivity-v0.3` | 완료 |
| 2026-08-24 | Target trial 인과 명세 초안 | `docs/target-trial-protocol.md` | 초안 완료 |
| 2026-08-24 | K-CURE 변수 요청·매핑 명세 | `docs/k-cure-variable-dictionary.md` | 초안 완료 |
| 2026-08-26 | 탐색 예산 스케일링 진단 v0.4 | `analysis/12`, `reports/budget-scaling-v0.4` | 완료 |
| 2026-08-27 | MCTS·MDP 개념 감사와 환경 편향 수정 v0.5 | `analysis/13`, `configs/dynamic_v0_5.json`, `reports/environment-fix-v0.5` | 완료 |
| 2026-08-27 | 발표용 자료 정비 — Figure 18~23·이야기·숫자화해·제안서대조 | `analysis/14`, `docs/research-story.md`, `docs/results-reconciliation.md`, `docs/proposal-vs-delivered.md`, `/story` | 완료 |
| 2026-08-27 | IPW 표적시험 에뮬레이션 v0.6 | `analysis/causal`, `analysis/15~16`, `reports/ipw-target-trial-v0.6`, Figure 24~26 | 완료 |
| 2026-08-27 | 이중강건 추정·결정별 식별 가능성 지도 v0.7 | `analysis/17~18`, `reports/doubly-robust-v0.7`, Figure 27~28 | 완료 |
| 2026-08-28 | v0.5 환경 재실행 검증 | `reports/*-v0.5env`, Figure 30 | 완료 |
| 2026-08-28 | 2차원 상호작용 민감도 v0.8 | `analysis/19~20`, `reports/interaction-sensitivity-v0.8`, Figure 29 | 완료 (일부 정정됨) |
| 2026-08-28 | 가치판단 상호작용 확인 v0.9 | `analysis/22~23`, `reports/utility-interaction-v0.9`, Figure 31 | 완료 |
| 다음 단계 | 인과추정 시제품·상호작용 민감도·K-CURE 확보 | IPW/g-방법, 2D 민감도, utility 사전등록, 코드북 매핑, 임상 검토 | 예정 |

## 2026-07-11 기준 현재 위치

```text
METABRIC 원본
  -> 전처리 환자표
  -> NCCN 정책 A
  -> 생존 보상모형
  -> 4단계 치료환경 (HR 적격성 반영, 환자별 8~16개 경로)
  -> UCT-MCTS 정책 B
  -> 보류 테스트셋 비교와 완전탐색 검증 (v0.1)
  -> OS/RFS 기반 5년 확률환경 (v0.2)
  -> 환자별 18~135개 의사결정 궤적
  -> 반응 후 재계획하는 stochastic MCTS

다음: 합성 전이 -> K-CURE 실제 추정, 관찰 연관성 -> 인과효과 추정
```

MCTS v0.1의 수치·방법·한계는
[기술 보고서](reports/mcts-poc-v1/README.md)와
`src/minutes/2026-07-11-mcts-poc-v1.md`에서 확인할 수 있습니다.

동적 v0.2는 [기술 보고서](reports/dynamic-mcts-poc-v0.2/README.md)와
`src/minutes/2026-07-11-dynamic-environment-v02.md`에서 확인할 수 있습니다.

v0.4 탐색 예산 진단은 [기술 보고서](reports/budget-scaling-v0.4/README.md)와
`src/minutes/2026-08-26-budget-scaling-v04.md`에서 확인할 수 있습니다.

v0.5 환경 감사는 [기술 보고서](reports/environment-fix-v0.5/README.md)와
`src/minutes/2026-08-27-environment-audit-v05.md`에서 확인할 수 있습니다.
**새 실험은 `configs/dynamic_v0_5.json`을 씁니다.** v0.2~v0.4 리포트 수치는
편향이 있던 환경에서 나온 것이며, 아직 v0.5로 재실행하지 않았습니다.
이 진단 이후 **기본 탐색 예산은 1024 이상**을 쓰며, 256으로 낸 v0.2·v0.3의
결정 단위 수치는 탐색 해상도 미달이라는 단서와 함께 읽어야 합니다.
