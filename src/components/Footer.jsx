import './Footer.css'

export default function Footer() {
  return (
    <footer className="footer">
      <div className="footer-inner">
        <div className="footer-brand">
          <span className="footer-logo">MCTS-ONC</span>
          <p className="footer-tagline">Monte Carlo Tree Search for Oncology</p>
        </div>
        <div className="footer-info">
          <p>경북대학교 의과대학</p>
          <p>융합형 의사과학자 양성사업 (학부과정)</p>
          <p>2026년도 의과학 연구동아리</p>
        </div>
        <div className="footer-links">
          <a href="https://github.com/reus3148-blip/MCTS" target="_blank" rel="noopener noreferrer">
            GitHub
          </a>
        </div>
      </div>
      <div className="footer-bottom">
        <p>© 2026 MCTS-ONC. 경북대학교 의과대학.</p>
      </div>
    </footer>
  )
}
