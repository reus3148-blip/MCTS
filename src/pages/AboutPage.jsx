import './AboutPage.css'

const infoItems = [
  { label: '소속 대학', value: '경북대학교 의과대학' },
  { label: '지원 사업', value: '융합형 의사과학자 양성사업 (학부과정)' },
  { label: '연구 분야', value: '의료 인공지능 · 강화학습 · 유방암 임상 의사결정' },
  { label: '활동 기간', value: '2026년 (9개월 과정)' },
  { label: '동아리명', value: 'MCTS-ONC — Monte Carlo Tree Search for Oncology' },
]

const links = [
  { label: 'GitHub 저장소', url: 'https://github.com/reus3148-blip/MCTS', desc: '웹사이트 및 연구 코드' },
  { label: '경북대학교 의과대학', url: 'https://medicine.knu.ac.kr', desc: '소속 기관' },
  { label: '융합형 의사과학자 양성사업', url: 'https://www.khidi.or.kr', desc: '한국보건산업진흥원 (KHIDI)' },
]

export default function AboutPage() {
  return (
    <div className="page about-page">
      <div className="section">
        <p className="section-label">About</p>
        <h2 className="section-title">동아리 소개</h2>
        <p className="section-desc">
          MCTS-ONC는 경북대학교 의과대학 학부생들이 자율적으로 운영하는 의과학 연구동아리입니다.
          유방암 치료 의사결정 문제를 인공지능 강화학습 관점에서 분석하고, 표준 가이드라인 기반
          의사결정과의 비교를 목표로 합니다.
        </p>

        <div className="about-block">
          <h3 className="about-h3">기본 정보</h3>
          <div className="info-table">
            {infoItems.map((i) => (
              <div key={i.label} className="info-row">
                <div className="info-label">{i.label}</div>
                <div className="info-value">{i.value}</div>
              </div>
            ))}
          </div>
        </div>

        <div className="about-block">
          <h3 className="about-h3">동아리 비전</h3>
          <p className="about-p">
            본 동아리는 임상 진료에 직접 적용할 도구를 만드는 것이 아니라, <strong>의과학적 문제를
            데이터 과학·인공지능 알고리즘의 관점에서 구조화하고 검증하는 훈련</strong>에 중점을 둔다.
          </p>
          <p className="about-p">
            참여 학생들은 유방암 임상 가이드라인과 주요 논문을 학습하고, 공개 데이터 전처리, 변수 정의,
            알고리즘 구현, 결과 해석을 직접 수행함으로써 <strong>융합형 의사과학자로서 필요한 기초
            연구 역량</strong>을 기른다.
          </p>
        </div>

        <div className="about-block">
          <h3 className="about-h3">관련 링크</h3>
          <div className="link-grid">
            {links.map((l) => (
              <a
                key={l.label}
                href={l.url}
                target="_blank"
                rel="noopener noreferrer"
                className="link-card"
              >
                <div className="link-label">{l.label}</div>
                <div className="link-desc">{l.desc}</div>
                <div className="link-arrow">↗</div>
              </a>
            ))}
          </div>
        </div>

        <div className="about-footer">
          <p>© 2026 MCTS-ONC. 경북대학교 의과대학.</p>
        </div>
      </div>
    </div>
  )
}
