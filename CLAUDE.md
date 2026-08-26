# CLAUDE.md

이 파일은 Claude Code가 이 저장소에서 작업할 때 참고하는 지침서입니다.

## 프로젝트 개요

MCTS-ONC (Monte Carlo Tree Search for Oncology) — 경북대학교 의과대학 융합형 의사과학자 양성사업 학부과정 의과학 연구동아리 웹사이트.

연구 주제: 강화학습 기반 몬테카를로 트리 탐색을 활용한 유방암 치료 의사결정 모형 개발 및 표준 가이드라인 기반 정책과의 비교 연구.

## 기술 스택

- **React 18** + **Vite 5** + **react-router-dom 7**
- 순수 CSS (CSS-in-JS / Tailwind 미사용)
- 배포: **Vercel** (GitHub main 브랜치 push 시 자동 배포)

## ⚠️ 호환성 요구사항 (필수)

**모바일과 PC 양쪽 모두에서 완벽하게 작동해야 합니다.**

- 모든 페이지·컴포넌트는 `≤ 640px` (모바일), `≤ 768px` (태블릿/햄버거 메뉴 전환점), `≥ 769px` (데스크탑) 세 구간에서 깨짐 없이 표시되어야 함
- 새 컴포넌트를 만들 때 반드시 모바일 폭(360~414px) 기준에서 가독성·터치 영역·줄바꿈을 검증할 것
- 텍스트 가로 오버플로우, 가로 스크롤, 글자 잘림은 절대 허용 금지
- 그리드/플렉스 레이아웃은 `grid-template-columns: repeat(auto-fit, minmax(...))` 또는 미디어 쿼리로 모바일에서 단일 열로 전환
- 폰트 크기는 `clamp()` 또는 미디어 쿼리로 모바일에서 축소
- Navbar는 768px 이하에서 햄버거 메뉴로 전환

## 폴더 구조

```
src/
├── components/        # 재사용 컴포넌트 (Navbar, Hero, Research, Roadmap, Team, Footer)
├── pages/             # 라우터 페이지 (Home, BreastCancerPage, MCTSPage, ResearchPage, RoadmapPage, TeamPage)
├── App.jsx            # BrowserRouter + Routes
├── main.jsx           # 엔트리포인트
└── index.css          # 전역 변수·section 클래스
```

각 컴포넌트는 동일 이름의 `.css` 파일을 함께 둠 (예: `Hero.jsx` + `Hero.css`).

## 연구 코드 (Python) — 웹사이트와 별개

이 저장소는 **역할이 둘**이다. 위 `src/`(React 블로그, `mcts.blundermate.app`으로 배포)와,
아래 `analysis/`의 **MCTS 연구 파이프라인**(로컬에서 `py`로 실행, 배포되지 않음).

```
analysis/
├── 01~05_*.py       # METABRIC 전처리·시각화·NCCN 정책 일치율
├── 06~07_*.py       # MCTS 정적 PoC v1 (Cox 보상모형 + UCT 탐색)
├── 08~09_*.py       # 동적 확률환경 v0.2 (chance-aware UCT)
├── 10~11_*.py       # 다중 시드 강건성·합성 가정 민감도 v0.3
├── 12_*.py          # 탐색 예산 스케일링 진단 v0.4
├── mcts/            # 치료환경·보상모형·UCT 탐색 모듈
└── dynamic/         # 공통 스키마·확률 전이·stochastic MCTS
                     #   cohort.py: 10~12 공유 코호트·보상모형·매니페스트
                     #   experiment_utils.py: 테스트 가능한 통계 헬퍼
tests/               # 단위 테스트 (py -m unittest discover -s tests -v)
reports/             # 재현 가능한 기술 리포트 (metrics·manifest·표·figure)
```

- **현재 상태·실행법·다음 단계·해석 경계**는 `README.md`, `PROJECT_TIMELINE.md`,
  `reports/*/README.md`, `src/minutes/`(날짜별 회의록)에 정리되어 있음 — 연구 작업 전 반드시 참고.
- **데이터**(`data/`)는 대용량이라 git 제외. 이 PC엔 `data/processed/*.csv`가 있어 바로 실행됨.
  새 환경에서는 `analysis/01~02`로 재생성 필요.
- 연구 코드 초안은 주로 Codex/Claude로 작성했고, 임상 규칙·인과추론·결과 해석은 사람 검토 대상(회의록에 명시).
- **실험 규약**: 번호 붙은 `analysis/1x_*.py`는 모듈로 import할 수 없으므로, 재사용·테스트
  대상 로직은 `analysis/dynamic/`에 두고 스크립트는 얇게 유지한다. 새 실험은
  `reports/<label>/`에 `metrics.json`·`run_manifest.json`·`tables/`·`README.md`를 함께 낸다.
- **탐색 예산**: v0.4 진단 이후 기본값은 **1024 이상**(256은 행동 순서를 분해하지 못함).

## 디자인 시스템

CSS 변수는 `src/index.css`의 `:root`에 정의:
- `--bg-dark`, `--bg-card`, `--bg-card2` — 배경 단계
- `--accent` (#00c9b1 청록), `--accent2` (#3b82f6 파랑) — 포인트
- `--text-primary`, `--text-secondary` — 텍스트
- `--border` — 카드/구분선

공통 섹션 헤더는 `.section-label` (대문자 라벨) + `.section-title` + `.section-desc` 3종 구조.

## 팀 정보 (사이트와 일치 필수)

- **전병우** (대표, 의학과 2020110089) — MDP 설계, 외부 컨택, 논문 의학 파트
- **강홍준** (팀원, 의학과 2021113048) — 코딩, 모델 학습, 실험 설계
- **서재환** (팀원, 의학과 2021110282) — 데이터 전처리, 가상 환자 환경 구축, 결과 시각화

## 비공개 정보 (웹사이트 노출 금지)

`팀_내부_노트.md`는 대외비 문서로 git에 커밋되어 있긴 하지만 다음 정보는 사이트에 반영하지 않음:
- 자문 시도 실패 사실
- 팀 코딩 수준 자평
- 예산 정정 내역
- 특정 임상의/병원과의 사적 관계

## 배포 워크플로우

```bash
git add .
git commit -m "..."
git push   # → Vercel 자동 배포
```

GitHub Pages는 사용하지 않음. `vercel.json`이 SPA 라우팅을 처리.
