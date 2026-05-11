import { NavLink } from 'react-router-dom'
import './Navbar.css'

const links = [
  { to: '/',              label: 'Home' },
  { to: '/breast-cancer', label: '유방암' },
  { to: '/mcts',          label: 'MCTS' },
  { to: '/research',      label: 'Research' },
  { to: '/roadmap',       label: 'Roadmap' },
  { to: '/team',          label: 'Team' },
]

export default function Navbar() {
  return (
    <nav className="navbar">
      <NavLink to="/" className="navbar-logo">MCTS-ONC</NavLink>
      <ul className="navbar-links">
        {links.map((l) => (
          <li key={l.to}>
            <NavLink
              to={l.to}
              end={l.to === '/'}
              className={({ isActive }) => isActive ? 'nav-link active' : 'nav-link'}
            >
              {l.label}
            </NavLink>
          </li>
        ))}
      </ul>
    </nav>
  )
}
