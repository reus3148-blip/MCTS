import { Link } from 'react-router-dom'
import './Hero.css'

export default function Hero() {
  return (
    <section className="hero">
      <div className="hero-content">
        <div className="hero-badge">경북대학교 의과대학 · 융합형 의사과학자 양성사업</div>
        <h1 className="hero-title">
          <span className="hero-title-accent">MCTS-ONC</span>
        </h1>
        <p className="hero-subtitle">Monte Carlo Tree Search for Oncology</p>
        <p className="hero-desc">
          유방암 치료 의사결정을 강화학습과 몬테카를로 트리 탐색으로 모델링하는<br />
          경북대학교 의과대학 의과학 연구동아리
        </p>
        <div className="hero-tags">
          <span className="tag">강화학습 (RL)</span>
          <span className="tag">몬테카를로 트리 탐색 (MCTS)</span>
          <span className="tag">유방암</span>
          <span className="tag">의사결정 모델링</span>
        </div>
        <div className="hero-cta">
          <Link to="/research" className="btn-primary">연구 보기</Link>
          <Link to="/team" className="btn-secondary">팀 소개</Link>
        </div>
      </div>
    </section>
  )
}
