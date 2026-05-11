import './LearnPage.css'

const subtypes = [
  {
    name: 'HR+ / HER2−',
    altName: '루미날 A·B',
    ratio: '72.7%',
    profile: 'ER/PR 양성 · HER2 음성',
    desc: '가장 큰 비중을 차지하며 호르몬치료에 잘 반응한다. Ki-67 수치에 따라 Luminal A(예후 양호)와 Luminal B(증식 빠름)로 세분화된다.',
    color: '#00c9b1',
  },
  {
    name: 'Triple Negative',
    altName: 'TNBC',
    ratio: '12.2%',
    profile: 'ER/PR 음성 · HER2 음성',
    desc: '세 수용체 모두 음성이라 표적이 없다. 항암화학요법이 주된 치료이며 재발 위험과 예후 부담이 가장 크다.',
    color: '#f87171',
  },
  {
    name: 'HR+ / HER2+',
    altName: '루미날 HER2',
    ratio: '10.3%',
    profile: 'ER/PR 양성 · HER2 양성',
    desc: '호르몬치료와 HER2 표적치료(트라스투주맙 등)를 모두 받을 수 있다. 치료 옵션이 가장 다양한 유형.',
    color: '#a78bfa',
  },
  {
    name: 'HR− / HER2+',
    altName: 'HER2-enriched',
    ratio: '4.6%',
    profile: 'ER/PR 음성 · HER2 양성',
    desc: '공격적이지만 HER2 표적치료제 등장 이후 예후가 크게 개선되었다. 표적치료가 치료의 중심.',
    color: '#3b82f6',
  },
]

const treatments = [
  { title: '수술', desc: '유방보존술 또는 전절제술. 가능한 경우 가장 먼저 고려되는 국소 치료.' },
  { title: '방사선치료', desc: '유방보존술 후 잔존 미세 병변 제거. 국소 재발 위험을 낮춘다.' },
  { title: '항암화학요법', desc: '전신을 순환하며 암세포 분열을 차단. 수술 전후 보조요법으로 사용.' },
  { title: '호르몬치료', desc: 'HR 양성 종양에서 에스트로겐 신호를 차단 (Tamoxifen, Aromatase Inhibitor 등).' },
  { title: '표적치료', desc: 'HER2, CDK4/6 등 특정 분자 표적을 정밀 공격 (Trastuzumab, Palbociclib 등).' },
  { title: '면역치료', desc: 'PD-L1 양성 TNBC 등에서 면역관문 억제제로 면역세포가 암을 공격하도록 유도.' },
]

export default function BreastCancerPage() {
  return (
    <div className="page learn-page">
      <div className="section">
        <p className="section-label">Learn · Breast Cancer</p>
        <h2 className="section-title">유방암 이해하기</h2>
        <p className="section-desc">
          유방암은 한국 여성암 발생률 1위 질환이지만, 사실 단일한 병이 아니다.
          종양의 분자생물학적 특성에 따라 전혀 다른 질병처럼 행동하며, 그래서 치료도 달라진다.
        </p>

        <div className="learn-block">
          <h3 className="learn-h3">1. 유방암은 한 가지 병이 아니다</h3>
          <p className="learn-p">
            과거에는 유방암을 모두 같은 병으로 보고 동일한 항암제로 치료했다.
            그러나 1990년대 이후 분자생물학 연구가 발전하면서, 유방암 세포 안에 어떤 <strong>수용체(receptor)</strong>가
            발현되어 있느냐에 따라 종양이 자라는 방식과 약물 반응이 완전히 다르다는 것이 밝혀졌다.
          </p>
          <p className="learn-p">
            대표적으로 보는 세 가지 표지자가 <strong>ER(에스트로겐 수용체)</strong>, <strong>PR(프로게스테론 수용체)</strong>,
            <strong> HER2(인간 표피성장인자 수용체 2)</strong> 이며, 이 조합으로 유방암을 네 가지 분자아형(molecular subtype)으로 나눈다.
          </p>
        </div>

        <div className="subtypes-grid">
          {subtypes.map((s) => (
            <div key={s.name} className="subtype-card" style={{ borderTopColor: s.color }}>
              <div className="subtype-header">
                <div>
                  <h4 className="subtype-name">{s.name}</h4>
                  <p className="subtype-alt">{s.altName}</p>
                </div>
                <span className="subtype-ratio" style={{ color: s.color }}>{s.ratio}</span>
              </div>
              <p className="subtype-profile" style={{ color: s.color }}>{s.profile}</p>
              <p className="subtype-desc">{s.desc}</p>
            </div>
          ))}
        </div>

        <div className="strategy-box">
          <div className="strategy-label">본 동아리의 연구 전략</div>
          <p className="strategy-text">
            데이터 가용성이 가장 높은 <strong style={{ color: '#00c9b1' }}>HR+/HER2− (72.7%)</strong> 아형부터 모델링을 시작하고,
            가이드라인 개선 여지가 큰 <strong style={{ color: '#f87171' }}>Triple Negative (12.2%)</strong>로 확장한다.
            희귀 아형은 충분한 표본을 확보하기 어려우므로 후속 연구로 분리한다.
          </p>
        </div>

        <div className="learn-block">
          <h3 className="learn-h3">2. 병기 분류 — TNM 시스템</h3>
          <p className="learn-p">
            분자아형이 "어떤 종류의 암인가"를 말해준다면, 병기(stage)는 "얼마나 진행되었나"를 말해준다.
            국제 표준인 TNM 시스템은 세 가지 정보를 조합한다.
          </p>
          <div className="tnm-grid">
            <div className="tnm-item">
              <span className="tnm-letter">T</span>
              <div>
                <div className="tnm-title">Tumor</div>
                <div className="tnm-desc">원발 종양의 크기와 침윤 깊이</div>
              </div>
            </div>
            <div className="tnm-item">
              <span className="tnm-letter">N</span>
              <div>
                <div className="tnm-title">Node</div>
                <div className="tnm-desc">주변 림프절로의 전이 정도</div>
              </div>
            </div>
            <div className="tnm-item">
              <span className="tnm-letter">M</span>
              <div>
                <div className="tnm-title">Metastasis</div>
                <div className="tnm-desc">다른 장기로의 원격 전이 여부</div>
              </div>
            </div>
          </div>
          <p className="learn-p">
            이 세 가지를 조합해 1기 ~ 4기로 나누며, 1~3기는 다시 A/B/C 등 하위 단계로 세분화된다.
            병기가 높을수록 일반적으로 예후가 나쁘고 더 공격적인 치료가 필요하다.
          </p>
        </div>

        <div className="learn-block">
          <h3 className="learn-h3">3. 치료 옵션</h3>
          <p className="learn-p">
            유방암 치료는 크게 <strong>국소치료(수술·방사선)</strong> 와 <strong>전신치료(항암·호르몬·표적·면역)</strong> 로 나뉜다.
            대부분의 환자는 여러 치료를 순차적으로 받게 된다.
          </p>
          <div className="treatments-grid">
            {treatments.map((t) => (
              <div key={t.title} className="treatment-card">
                <h4 className="treatment-title">{t.title}</h4>
                <p className="treatment-desc">{t.desc}</p>
              </div>
            ))}
          </div>
        </div>

        <div className="learn-block highlight-block">
          <h3 className="learn-h3">4. 왜 의사결정이 복잡한가?</h3>
          <p className="learn-p">
            유방암 치료는 <strong>한 번의 결정이 아니라 일련의 결정</strong>이다.
            예를 들어 한 환자의 치료 경로를 따라가 보면:
          </p>
          <div className="decision-flow">
            <div className="flow-step">진단 (조직검사 + IHC)</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">분자아형 분류 (ER/PR/HER2)</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">수술 가능 여부 판단</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">수술 전 보조요법 (Neoadjuvant) 여부</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">수술 → 병리 결과로 잔존암 확인</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">수술 후 보조요법 (Adjuvant): 항암 + 호르몬 + 방사선</div>
            <div className="flow-arrow">↓</div>
            <div className="flow-step">5~10년 추적 관찰, 재발 시 다음 치료 결정</div>
          </div>
          <p className="learn-p">
            각 단계의 선택이 다음 단계의 가능한 선택지를 바꾼다. 이런 문제 구조를
            수학적으로 표현한 것이 바로 <strong>마르코프 결정 과정(MDP)</strong>이며,
            본 동아리는 이 문제에 <strong>강화학습과 몬테카를로 트리 탐색</strong>을 적용한다.
          </p>
        </div>

        <div className="learn-sources">
          <p className="sources-title">참고 자료</p>
          <ul>
            <li>국가암정보센터 유방암 페이지 (cancer.go.kr)</li>
            <li>Susan G. Komen — Molecular Subtypes of Breast Cancer</li>
            <li>NCCN Clinical Practice Guidelines in Oncology — Breast Cancer</li>
            <li>제10차 한국유방암 진료권고안 (2023)</li>
          </ul>
        </div>
      </div>
    </div>
  )
}
