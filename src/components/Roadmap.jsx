import './Roadmap.css'

const steps = [
  {
    phase: '1단계',
    period: '1~2개월',
    title: '이론 학습 및 환경 정의',
    items: [
      '강화학습, MDP, MCTS 이론 학습',
      '유방암 진료 가이드라인 및 주요 논문 정독',
      '상태·행동·보상 함수 정의',
      'K-CURE 유방암 데이터 신청서 초안 작성',
    ],
  },
  {
    phase: '2단계',
    period: '3~4개월',
    title: '데이터 수집 및 가상 환자 환경 구축',
    items: [
      'TCGA-BRCA, METABRIC 데이터 수집 및 전처리',
      '임상 변수 및 분자생물학적 변수 정리',
      '단순화된 가상 환자 환경 구축',
      'K-CURE 보완 요청 대응',
    ],
  },
  {
    phase: '3단계',
    period: '5~7개월',
    title: '알고리즘 구현 및 비교 실험',
    items: [
      'RL/MCTS 알고리즘 구현',
      'NCCN 가이드라인 기반 정책 vs RL 기반 정책 비교',
      '생존 지표, 재발 위험, 치료 독성, 의사결정 일치율 분석',
      'K-CURE 접근 시 한국인 데이터 기초 통계 분석',
    ],
  },
  {
    phase: '4단계',
    period: '8~9개월',
    title: '결과 정리 및 발표',
    items: [
      '실험 결과 정리 및 임상 타당성 검토',
      '학회 발표 자료 작성',
      '최종 보고서 작성',
      '논문 초고 작성',
    ],
  },
]

export default function Roadmap() {
  return (
    <section id="roadmap">
      <div className="section">
        <p className="section-label">Roadmap</p>
        <h2 className="section-title">연구 로드맵</h2>
        <p className="section-desc">9개월간의 단계적 연구 활동 계획</p>
        <div className="roadmap">
          {steps.map((s, i) => (
            <div key={s.phase} className="roadmap-item">
              <div className="roadmap-left">
                <div className="roadmap-dot" />
                {i < steps.length - 1 && <div className="roadmap-line" />}
              </div>
              <div className="roadmap-card">
                <div className="roadmap-meta">
                  <span className="roadmap-phase">{s.phase}</span>
                  <span className="roadmap-period">{s.period}</span>
                </div>
                <h3 className="roadmap-title">{s.title}</h3>
                <ul className="roadmap-list">
                  {s.items.map((item) => (
                    <li key={item}>{item}</li>
                  ))}
                </ul>
              </div>
            </div>
          ))}
        </div>
      </div>
    </section>
  )
}
