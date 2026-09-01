import './Checklist.css'

const sections = [
  {
    label: '0',
    title: '행정 / 제출',
    items: [
      { text: '팀원 확정 (전병우 / 강홍준 / 서재환)', done: true },
      { text: '팀원 개인정보 동의서 작성', done: true },
      { text: '제안서 제출', done: true },
      { text: '중간 보고', done: false },
      { text: '최종 보고 / 발표', done: false },
    ],
  },
  {
    label: '1',
    title: '학습 단계',
    period: '1~3개월차 · 팀원 공통',
    items: [
      { text: 'Google Colab 계정 개설', done: false },
      { text: 'Python 기초 (변수·조건·반복·함수)', done: false },
      { text: 'pandas — 환자 표 다루기', done: false },
      { text: 'NumPy 기본', done: false },
      { text: 'matplotlib / seaborn — 그래프', done: false },
      { text: 'scikit-learn 기초 — 전처리·분류', done: false },
      { text: 'gymnasium — RL 환경 인터페이스', done: false },
      { text: 'PyTorch 튜토리얼 (DQN까지)', done: false },
      { text: 'MCTS 알고리즘 학습 (오픈소스 코드 분석)', done: false },
    ],
  },
  {
    label: '2',
    title: '데이터 확보',
    items: [
      { text: 'TCGA-BRCA 다운로드, 변수 목록 정리', done: false },
      { text: 'METABRIC 다운로드, 변수 목록 정리', done: true },
      { text: 'K-CURE 공통 환자 스키마·이행 변수 설계', done: true },
      { text: 'Target trial 인과 명세 초안', done: true },
      { text: '두 데이터셋 변수 매핑 (공통 컬럼 일치화)', done: false },
      { text: 'K-CURE 데이터 신청서 초안 작성', done: false },
      { text: 'K-CURE 승인 추적', done: false },
    ],
  },
  {
    label: '3',
    title: '임상 자료 / NCCN',
    items: [
      { text: 'NCCN 유방암 가이드라인 정독', done: false },
      { text: '한국유방암학회 권고안 정독', done: false },
      { text: 'NCCN + AI 관련 최신 논문 정독 및 차별점 정리', done: false },
      { text: 'NCCN 규칙을 정책 매트릭스(표/코드)로 정리 — 1차 단순화', done: true },
      { text: '임상 자문 확보 시도', done: false },
    ],
  },
  {
    label: '4',
    title: '연구 파이프라인',
    subsections: [
      {
        title: '① 데이터 전처리',
        period: '1~2개월',
        items: [
          { text: '결측치 처리 / 이상치 제거 / 변수 표준화', done: true },
          { text: '치료명 텍스트 → 카테고리 변수 변환', done: true },
          { text: '분자아형(HR+/HER2-, TNBC, HR+/HER2+, HR-/HER2+) 라벨링', done: true },
          { text: '정제된 환자 데이터셋 산출', done: true },
        ],
      },
      {
        title: '② 가상 환자 환경 구축',
        period: '2~3개월',
        items: [
          { text: '상태(State) 정의 — 환자 특성 + 치료 경로 prefix', done: true },
          { text: '행동(Action) 정의 — 수술·항암·호르몬·방사선', done: true },
          { text: '보상(Reward) 함수 1차 — Cox 기반 5년 OS 예측', done: true },
          { text: '전이 확률 v0.2 — METABRIC OS/RFS 기준위험 + 공개 합성 가정', done: true },
          { text: '정적 4단계 치료환경 구현 (PoC v0.1)', done: true },
          { text: '동적 5년 확률환경 구현 (PoC v0.2)', done: true },
          { text: '완전탐색 대조 sanity check', done: true },
        ],
      },
      {
        title: '③ 두 정책 구현',
        period: '3~5개월',
        items: [
          { text: 'NCCN 정책 — 단순화 가이드라인 규칙 코딩', done: true },
          { text: 'UCT-MCTS 정책 — 정적 v0.1 + 확률적 v0.2', done: true },
          { text: '(옵션) DQN 베이스라인', done: false },
        ],
      },
      {
        title: '④ 시뮬레이션',
        period: '5~7개월',
        items: [
          { text: '가상 환자 1,000명 표본 추출', done: false },
          { text: '보류 테스트셋 NCCN 정책 평가', done: true },
          { text: '보류 테스트셋 MCTS 정책 평가', done: true },
          { text: '동적 환경 40명 × 정책별 100회 episode', done: true },
          { text: '시드 고정·재현성 확인', done: true },
        ],
      },
      {
        title: '⑤ 결과 분석',
        period: '7~9개월',
        items: [
          { text: '합성 동적환경 5년 생존·재발·독성 비교', done: true },
          { text: '정책 결정 일치율 계산', done: true },
          { text: '아형별 서브 분석', done: true },
          { text: '그래프·표 시각화', done: true },
          { text: '다중 시드 강건성 검증 (v0.3)', done: true },
          { text: '합성 가정 민감도 분석 (v0.3)', done: true },
          { text: '탐색 예산 스케일링 진단 — 결정 불안정성의 원인 규명 (v0.4)', done: true },
          { text: 'MCTS·MDP 개념 감사 — 환경 편향 발견·수정 (v0.5)', done: true },
          { text: 'IPW 표적시험 에뮬레이션 — 단일 결정 인과추정 시제품 (v0.6)', done: true },
          { text: 'doubly robust(AIPW)·검열 가중(IPCW)', done: false },
          { text: '2차원(상호작용) 민감도 격자', done: false },
          { text: 'v0.5 환경으로 v0.3·v0.4 재실행', done: false },
        ],
      },
    ],
  },
  {
    label: '5',
    title: '산출물',
    items: [
      { text: 'MCTS PoC v0.1 기술 보고서', done: true },
      { text: '동적 MCTS PoC v0.2 기술 보고서', done: true },
      { text: '강건성·민감도 v0.3 기술 보고서', done: true },
      { text: '탐색 예산 스케일링 v0.4 기술 보고서', done: true },
      { text: '환경 편향 수정 v0.5 기술 보고서', done: true },
      { text: '발표용 Figure 18~23 (v0.3~v0.5)', done: true },
      { text: '연구 이야기·숫자 화해·제안서 대조 문서', done: true },
      { text: 'IPW 표적시험 v0.6 기술 보고서', done: true },
      { text: '중간 보고서 초안', done: false },
      { text: '학회 발표 슬라이드 / 포스터', done: false },
      { text: '논문 / 최종 보고서 초안', done: false },
      { text: 'Limitation 섹션 정리', done: true },
      { text: '임상 검토 의뢰', done: false },
    ],
  },
  {
    label: '6',
    title: '운영 / 인프라',
    items: [
      { text: 'GitHub 저장소 셋업 (코드 / 데이터 / 보고서 / 문서 분리)', done: true },
      { text: '클라우드 환경 결정 (Colab → AWS/GCP / 학교 서버)', done: false },
      { text: '정기 팀 회의 주기 확정', done: false },
      { text: '진행 상황 기록 (Minutes + 프로젝트 타임라인)', done: true },
    ],
  },
]

function flatItems(section) {
  if (section.subsections) {
    return section.subsections.flatMap((s) => s.items)
  }
  return section.items
}

function CheckItem({ text, done }) {
  return (
    <li className={`check-item ${done ? 'done' : ''}`}>
      <span className="check-box" aria-hidden="true">
        {done && (
          <svg viewBox="0 0 16 16" width="12" height="12" fill="none">
            <path d="M3 8.5L6.5 12L13 4.5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        )}
      </span>
      <span className="check-text">{text}</span>
    </li>
  )
}

export default function Checklist() {
  const totalAll = sections.reduce((sum, s) => sum + flatItems(s).length, 0)
  const doneAll = sections.reduce(
    (sum, s) => sum + flatItems(s).filter((i) => i.done).length,
    0
  )
  const percent = Math.round((doneAll / totalAll) * 100)

  return (
    <section id="checklist">
      <div className="section">
        <p className="section-label">Checklist</p>
        <h2 className="section-title">연구 체크리스트</h2>
        <p className="section-desc">
          9개월 연구 프로젝트의 전체 작업 목록입니다. 행정 단계부터 산출물까지 단계별로 정리했습니다.
        </p>

        <div className="checklist-progress">
          <div className="checklist-progress-meta">
            <span className="checklist-progress-label">전체 진행도</span>
            <span className="checklist-progress-count">
              {doneAll} / {totalAll} · {percent}%
            </span>
          </div>
          <div className="checklist-progress-bar">
            <div
              className="checklist-progress-fill"
              style={{ width: `${percent}%` }}
            />
          </div>
        </div>

        <div className="checklist">
          {sections.map((s) => {
            const items = flatItems(s)
            const sectionDone = items.filter((i) => i.done).length
            const sectionTotal = items.length
            return (
              <div key={s.label} className="checklist-section">
                <div className="checklist-section-head">
                  <div className="checklist-section-title-wrap">
                    <span className="checklist-section-label">{s.label}</span>
                    <h3 className="checklist-section-title">{s.title}</h3>
                  </div>
                  <div className="checklist-section-meta">
                    {s.period && <span className="checklist-section-period">{s.period}</span>}
                    <span className="checklist-section-count">
                      {sectionDone} / {sectionTotal}
                    </span>
                  </div>
                </div>

                {s.subsections ? (
                  <div className="checklist-subsections">
                    {s.subsections.map((sub) => (
                      <div key={sub.title} className="checklist-subsection">
                        <div className="checklist-subsection-head">
                          <h4 className="checklist-subsection-title">{sub.title}</h4>
                          {sub.period && (
                            <span className="checklist-subsection-period">{sub.period}</span>
                          )}
                        </div>
                        <ul className="check-list">
                          {sub.items.map((item) => (
                            <CheckItem key={item.text} {...item} />
                          ))}
                        </ul>
                      </div>
                    ))}
                  </div>
                ) : (
                  <ul className="check-list">
                    {s.items.map((item) => (
                      <CheckItem key={item.text} {...item} />
                    ))}
                  </ul>
                )}
              </div>
            )
          })}
        </div>
      </div>
    </section>
  )
}
