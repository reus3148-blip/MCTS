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
| 다음 단계 | 민감도·인과추론·K-CURE 매핑 | 여러 seed, target trial, 변수사전, 임상 검토 | 예정 |

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
