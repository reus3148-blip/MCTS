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

const stateVars = [
  { group: '환자 기본 특성', items: ['나이, 폐경 여부', '종양 크기·위치', '림프절 전이', '병기 (Stage I~IV)'] },
  { group: '분자생물학적 특성', items: ['ER / PR (호르몬 수용체)', 'HER2 양성/음성', 'Ki-67 (증식 속도)', 'Triple Negative 여부', 'BRCA1/2 돌연변이'] },
  { group: '치료 진행 중 변화', items: ['투여 약물 종류·용량', '치료 사이클 수', '부작용 지표 (혈액·간 수치)', '종양 크기 변화 추이'] },
]

const actionVars = [
  { name: '수술', detail: '유방보존술 vs 전절제술 / 선행항암 전후 타이밍' },
  { name: '항암화학요법', detail: '약물 종류 · 용량 조절 · 투여 간격 · 휴약기' },
  { name: '방사선치료', detail: '시행 여부 및 범위' },
  { name: '호르몬치료', detail: 'Tamoxifen, Aromatase Inhibitor 추가 여부' },
  { name: '표적치료', detail: 'Trastuzumab, CDK4/6 inhibitor 추가 여부' },
]

const rewards = [
  { sign: '+', label: '종양 크기 감소', color: '#16a34a' },
  { sign: '+', label: '무재발 생존 기간 연장', color: '#16a34a' },
  { sign: '−', label: '치료 부작용 발생', color: '#dc2626' },
  { sign: '−', label: '약물 내성 발현', color: '#dc2626' },
]

const dataSources = [
  { name: 'TCGA-BRCA', cost: '무료', scale: '약 1,098명', origin: 'NIH (미국)', url: 'portal.gdc.cancer.gov', desc: '유전체·임상 정보 통합. 분자아형 분류의 표준 레퍼런스.' },
  { name: 'METABRIC', cost: '무료', scale: '약 2,509명', origin: '캐나다·영국', url: 'cbioportal.org', desc: '장기 추적 임상 결과 풍부. 예후 분석에 유리.' },
  { name: 'K-CURE', cost: '승인 필요', scale: '약 226만 명 (전체 암)', origin: '한국 보건복지부', url: 'k-cure.mohw.go.kr', desc: '한국인 임상 데이터. 외부 검증 및 인종 차이 분석 용도.' },
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

        {/* === 차별점 섹션 === */}
        <div className="research-section">
          <h3 className="research-h3">연구의 차별점</h3>
          <p className="research-p">
            의료 AI는 이미 다양한 형태로 임상에 도입되고 있다. 하지만 대부분은 <strong>"가이드라인을 더 잘 따르도록 돕는"</strong> 역할에 집중되어 있다.
            본 연구는 한 단계 더 나아가, <strong>가이드라인 자체가 최적인지 정량적으로 검증</strong>하는 것이 목표다.
          </p>
          <div className="compare-grid">
            <div className="compare-col compare-col-existing">
              <div className="compare-tag">기존 의료 AI</div>
              <h4 className="compare-title">가이드라인을 <span className="muted">따른다</span></h4>
              <ul className="compare-list">
                <li>맘모그램 이미지로 위험도 예측 (스크리닝)</li>
                <li>LLM으로 가이드라인 추천 자동화</li>
                <li>EMR 연동으로 누락 검사 알림</li>
                <li><strong>단일 시점</strong>의 의사결정 보조</li>
              </ul>
            </div>
            <div className="compare-arrow">→</div>
            <div className="compare-col compare-col-ours">
              <div className="compare-tag accent">본 연구</div>
              <h4 className="compare-title">가이드라인을 <span className="accent">검증한다</span></h4>
              <ul className="compare-list">
                <li>NCCN 정책을 강화학습 정책과 비교</li>
                <li>탐색 기반 모델이 다른 답을 내놓는지 확인</li>
                <li>특정 환자군에서 추가 연구 가설 도출</li>
                <li><strong>순차적 치료 전략 전체</strong>를 평가</li>
              </ul>
            </div>
          </div>
        </div>

        {/* === MDP 설계 === */}
        <div className="research-section">
          <h3 className="research-h3">우리가 설계한 MDP</h3>
          <p className="research-p">
            마르코프 결정 과정(MDP)은 순차적 의사결정 문제를 <strong>상태(State)</strong>, <strong>행동(Action)</strong>,
            <strong> 보상(Reward)</strong>의 세 축으로 정의한다. 유방암 치료 의사결정을 다음과 같이 정식화한다.
          </p>

          <div className="mdp-block">
            <div className="mdp-header">
              <span className="mdp-tag" style={{ background: 'rgba(8,145,178,0.1)', color: '#0891b2', border: '1px solid rgba(8,145,178,0.2)' }}>STATE</span>
              <h4 className="mdp-title">상태 변수</h4>
            </div>
            <div className="state-groups">
              {stateVars.map((g) => (
                <div key={g.group} className="state-group">
                  <div className="state-group-name">{g.group}</div>
                  <ul className="state-list">
                    {g.items.map((it) => <li key={it}>{it}</li>)}
                  </ul>
                </div>
              ))}
            </div>
          </div>

          <div className="mdp-block">
            <div className="mdp-header">
              <span className="mdp-tag" style={{ background: 'rgba(37,99,235,0.1)', color: '#2563eb', border: '1px solid rgba(37,99,235,0.2)' }}>ACTION</span>
              <h4 className="mdp-title">행동 (치료 선택)</h4>
            </div>
            <div className="action-grid">
              {actionVars.map((a) => (
                <div key={a.name} className="action-card">
                  <div className="action-name">{a.name}</div>
                  <div className="action-detail">{a.detail}</div>
                </div>
              ))}
            </div>
          </div>

          <div className="mdp-block">
            <div className="mdp-header">
              <span className="mdp-tag" style={{ background: 'rgba(124,58,237,0.1)', color: '#7c3aed', border: '1px solid rgba(124,58,237,0.2)' }}>REWARD</span>
              <h4 className="mdp-title">보상 함수</h4>
            </div>
            <div className="reward-row">
              {rewards.map((r) => (
                <div key={r.label} className="reward-item" style={{ borderColor: `${r.color}55` }}>
                  <span className="reward-sign" style={{ color: r.color }}>{r.sign}</span>
                  <span className="reward-label">{r.label}</span>
                </div>
              ))}
            </div>
            <p className="research-p-sub">
              ※ 실제 변수의 임상적 타당성과 데이터셋 가용성은 문헌 검토를 통해 지속적으로 보정한다.
            </p>
          </div>
        </div>

        {/* === 데이터 출처 === */}
        <div className="research-section">
          <h3 className="research-h3">데이터 출처</h3>
          <p className="research-p">
            공개 데이터로 가상 환자 환경을 먼저 구축하고, K-CURE 신청은 연구 초기부터 병행하여
            9개월 연구 기간 내 한국인 데이터 활용 가능성을 최대화한다.
          </p>
          <div className="data-table">
            <div className="data-row data-head">
              <div>데이터셋</div>
              <div>출처</div>
              <div>규모</div>
              <div>접근</div>
            </div>
            {dataSources.map((d) => (
              <div key={d.name} className="data-row">
                <div>
                  <div className="data-name">{d.name}</div>
                  <div className="data-desc">{d.desc}</div>
                </div>
                <div className="data-cell">
                  <span className="data-mob-label">출처</span>{d.origin}
                </div>
                <div className="data-cell">
                  <span className="data-mob-label">규모</span>{d.scale}
                </div>
                <div className="data-cell">
                  <span className="data-mob-label">접근</span>
                  <div className="data-cost">{d.cost}</div>
                  <div className="data-url">{d.url}</div>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>
    </section>
  )
}
