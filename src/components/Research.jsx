import './Research.css'

const cards = [
  {
    icon: '🎯',
    title: '연구 목적',
    desc: '유방암 치료 의사결정을 마르코프 결정 과정(MDP)으로 정식화하고, 강화학습 기반 탐색 정책과 NCCN 표준 가이드라인 정책을 정량적으로 비교합니다.',
  },
  {
    icon: '🧠',
    title: '핵심 방법론',
    desc: '몬테카를로 트리 탐색(MCTS)과 강화학습(RL)을 결합하여, 병기·호르몬 수용체·HER2 상태 등 임상 변수 기반의 순차적 치료 의사결정을 모델링합니다.',
  },
  {
    icon: '🗄️',
    title: '데이터',
    desc: 'TCGA-BRCA와 METABRIC 공개 유방암 데이터를 활용한 가상 환자 환경을 구축하고, K-CURE 한국인 유방암 데이터 신청을 병행합니다.',
  },
  {
    icon: '📊',
    title: '기대 결과',
    desc: '생존 지표, 재발 위험, 치료 독성, 의사결정 일치율을 기준으로 두 정책을 비교하고 개념 증명(Proof-of-Concept) 연구 결과를 도출합니다.',
  },
]

export default function Research() {
  return (
    <section id="research">
      <div className="section">
        <p className="section-label">Research</p>
        <h2 className="section-title">연구 주제</h2>
        <p className="section-desc">
          강화학습 기반 몬테카를로 트리 탐색을 활용한 유방암 치료 의사결정 모형 개발 및
          표준 가이드라인 기반 정책과의 비교 연구
        </p>
        <div className="research-cards">
          {cards.map((c) => (
            <div key={c.title} className="research-card">
              <div className="research-card-icon">{c.icon}</div>
              <h3 className="research-card-title">{c.title}</h3>
              <p className="research-card-desc">{c.desc}</p>
            </div>
          ))}
        </div>
        <div className="research-thesis">
          <div className="thesis-label">연구 제목 (영문)</div>
          <p className="thesis-text">
            Development of a Reinforcement Learning-Based Monte Carlo Tree Search Model
            for Breast Cancer Treatment Decision-Making and Comparison with Guideline-Based Policies
          </p>
        </div>
      </div>
    </section>
  )
}
