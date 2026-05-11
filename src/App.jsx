import { BrowserRouter, Routes, Route } from 'react-router-dom'
import Navbar from './components/Navbar'
import Home from './pages/Home'
import BreastCancerPage from './pages/BreastCancerPage'
import MCTSPage from './pages/MCTSPage'
import ResearchPage from './pages/ResearchPage'
import RoadmapPage from './pages/RoadmapPage'
import TeamPage from './pages/TeamPage'
import AboutPage from './pages/AboutPage'

function App() {
  return (
    <BrowserRouter>
      <Navbar />
      <main className="main-content">
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/breast-cancer" element={<BreastCancerPage />} />
          <Route path="/mcts" element={<MCTSPage />} />
          <Route path="/research" element={<ResearchPage />} />
          <Route path="/roadmap" element={<RoadmapPage />} />
          <Route path="/team" element={<TeamPage />} />
          <Route path="/about" element={<AboutPage />} />
        </Routes>
      </main>
    </BrowserRouter>
  )
}

export default App
