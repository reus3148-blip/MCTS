import './Team.css'

const members = [
  {
    name: '전병우',
    role: '대표',
    dept: '의학과',
    studentId: '2020110089',
    desc: '동아리 총괄 및 연구 설계, 강화학습 알고리즘 구현 담당',
    initial: '전',
  },
  {
    name: '강홍준',
    role: '팀원',
    dept: '의학과',
    studentId: '2021113048',
    desc: '유방암 임상 가이드라인 분석 및 문헌 리뷰 담당',
    initial: '강',
  },
  {
    name: '서재환',
    role: '팀원',
    dept: '의학과',
    studentId: '2021110282',
    desc: '데이터 전처리 및 가상 환자 환경 구축 담당',
    initial: '서',
  },
]

export default function Team() {
  return (
    <section id="team">
      <div className="section">
        <p className="section-label">Team</p>
        <h2 className="section-title">팀 소개</h2>
        <p className="section-desc">경북대학교 의학과 학부생 3인으로 구성된 연구팀</p>
        <div className="team-cards">
          {members.map((m) => (
            <div key={m.name} className="team-card">
              <div className="team-avatar">{m.initial}</div>
              <div className="team-role-badge">{m.role}</div>
              <h3 className="team-name">{m.name}</h3>
              <p className="team-dept">{m.dept} · {m.studentId}</p>
              <p className="team-desc">{m.desc}</p>
            </div>
          ))}
        </div>
      </div>
    </section>
  )
}
