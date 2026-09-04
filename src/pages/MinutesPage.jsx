import { useMemo } from 'react'
import { Link, useParams } from 'react-router-dom'
import { marked } from 'marked'
import './MinutesPage.css'

const rawFiles = import.meta.glob('../minutes/*.md', {
  eager: true,
  query: '?raw',
  import: 'default',
})

function parseFrontmatter(raw) {
  const m = raw.match(/^---\r?\n([\s\S]*?)\r?\n---\r?\n([\s\S]*)$/)
  if (!m) return { meta: {}, body: raw }
  const meta = {}
  m[1].split(/\r?\n/).forEach((line) => {
    const idx = line.indexOf(':')
    if (idx < 0) return
    const key = line.slice(0, idx).trim()
    const val = line.slice(idx + 1).trim()
    if (key) meta[key] = val
  })
  return { meta, body: m[2] }
}

const entries = Object.entries(rawFiles)
  .map(([path, raw]) => {
    const slug = path.split('/').pop().replace(/\.md$/, '')
    const { meta, body } = parseFrontmatter(raw)
    return {
      slug,
      date: meta.date || slug.slice(0, 10),
      // Several entries can share a date. `order` puts them back in the sequence
      // they were written; without it same-day entries fall back to filename
      // order, which showed v0.9 above v1.0.
      order: Number(meta.order) || 0,
      title: meta.title || slug,
      summary: meta.summary || '',
      body,
    }
  })
  .sort((a, b) => {
    if (a.date !== b.date) return a.date < b.date ? 1 : -1
    if (a.order !== b.order) return b.order - a.order
    return a.slug < b.slug ? 1 : -1
  })

marked.setOptions({ gfm: true, breaks: false })

function MinutesIndex() {
  return (
    <div className="page minutes-page">
      <div className="section">
        <p className="section-label">Minutes</p>
        <h2 className="section-title">회의록</h2>
        <p className="section-desc">
          MCTS-ONC 연구의 주요 의사결정을 시간 순으로 기록합니다.
          각 항목은 결정 사항·검토한 옵션·선택 사유·결과를 포함합니다.
        </p>

        <div className="minutes-list">
          {entries.length === 0 && (
            <p className="minutes-empty">아직 작성된 회의록이 없습니다.</p>
          )}
          {entries.map((e) => (
            <Link key={e.slug} to={`/minutes/${e.slug}`} className="minutes-card">
              <div className="minutes-card-date">{e.date}</div>
              <div className="minutes-card-title">{e.title}</div>
              {e.summary && <div className="minutes-card-summary">{e.summary}</div>}
              <div className="minutes-card-arrow">→</div>
            </Link>
          ))}
        </div>
      </div>
    </div>
  )
}

function MinutesDetail({ slug }) {
  const entry = entries.find((e) => e.slug === slug)
  const html = useMemo(() => (entry ? marked.parse(entry.body) : ''), [entry])

  if (!entry) {
    return (
      <div className="page minutes-page">
        <div className="section">
          <p className="section-label">Minutes</p>
          <h2 className="section-title">회의록을 찾을 수 없습니다</h2>
          <p className="section-desc">
            요청하신 항목이 존재하지 않습니다.
          </p>
          <Link to="/minutes" className="minutes-back">← 회의록 목록</Link>
        </div>
      </div>
    )
  }

  return (
    <div className="page minutes-page">
      <div className="section minutes-detail">
        <Link to="/minutes" className="minutes-back">← 회의록 목록</Link>
        <div className="minutes-detail-date">{entry.date}</div>
        <article
          className="markdown-body"
          dangerouslySetInnerHTML={{ __html: html }}
        />
      </div>
    </div>
  )
}

export default function MinutesPage() {
  const { slug } = useParams()
  return slug ? <MinutesDetail slug={slug} /> : <MinutesIndex />
}
