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
