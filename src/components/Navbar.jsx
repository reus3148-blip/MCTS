import { useState } from 'react'
import { NavLink } from 'react-router-dom'
import './Navbar.css'

const links = [
  { to: '/',              label: 'Home' },
  { to: '/breast-cancer', label: '유방암' },
  { to: '/mcts',          label: 'MCTS' },
  { to: '/research',      label: 'Research' },
  { to: '/roadmap',       label: 'Roadmap' },
  { to: '/checklist',     label: 'Checklist' },
  { to: '/team',          label: 'Team' },
  { to: '/about',         label: 'About' },
  { to: '/minutes',       label: 'Minutes' },
]

export default function Navbar() {
  const [open, setOpen] = useState(false)
  const close = () => setOpen(false)

  return (
    <nav className={`navbar ${open ? 'open' : ''}`}>
      <div className="navbar-bar">
        <NavLink to="/" className="navbar-logo" onClick={close}>MCTS-ONC</NavLink>
        <button
          className="navbar-toggle"
          aria-label="메뉴 열기"
          aria-expanded={open}
          onClick={() => setOpen((v) => !v)}
        >
          <span></span>
          <span></span>
          <span></span>
        </button>
        <ul className="navbar-links">
          {links.map((l) => (
            <li key={l.to}>
              <NavLink
                to={l.to}
                end={l.to === '/'}
                onClick={close}
                className={({ isActive }) => isActive ? 'nav-link active' : 'nav-link'}
              >
                {l.label}
              </NavLink>
            </li>
          ))}
        </ul>
      </div>
    </nav>
  )
}
