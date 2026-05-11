import './LearnPage.css'

const phases = [
  {
    num: '1',
    name: 'Selection (선택)',
    desc: '루트 노드에서 시작해 트리를 따라 내려간다. 각 단계에서 UCB 공식을 사용해 "유망하지만 충분히 탐색되지 않은" 자식 노드를 선택한다.',
    detail: '활용(exploitation)과 탐색(exploration)의 균형을 맞추는 핵심 단계',
  },
  {
    num: '2',
    name: 'Expansion (확장)',
    desc: '선택 단계가 끝난 잎(leaf) 노드에서, 아직 시도되지 않은 행동 중 하나를 골라 새 자식 노드를 트리에 추가한다.',
    detail: '한 번의 반복마다 새로운 노드 1개씩 트리가 자란다',
  },
  {
    num: '3',
    name: 'Simulation (시뮬레이션)',
    desc: '새로 추가된 노드에서 게임이 끝날 때까지 무작위(또는 간단한 규칙) 플레이를 진행한다. 이 빠른 모의실험을 "롤아웃(rollout)"이라 부른다.',
    detail: '값비싼 평가 함수 없이 결과를 추정하는 영리한 방법',
  },
  {
    num: '4',
    name: 'Backpropagation (역전파)',
    desc: '시뮬레이션 결과(승/패, 보상)를 루트까지 거슬러 올라가며 거쳐온 모든 노드의 방문 횟수와 누적 보상을 업데이트한다.',
    detail: '다음 Selection 단계의 의사결정에 이 통계가 사용됨',
  },
]

export default function MCTSPage() {
  return (
    <div className="page learn-page">
      <div className="section">
        <p className="section-label">Learn · MCTS</p>
        <h2 className="section-title">몬테카를로 트리 탐색</h2>
        <p className="section-desc">
          MCTS(Monte Carlo Tree Search)는 가능한 선택지가 너무 많아 전부 따져볼 수 없을 때,
          "무작위 시뮬레이션을 영리하게 반복"하여 최선의 결정을 찾아내는 알고리즘이다.
          DeepMind의 AlphaZero가 단 4시간 학습만으로 세계 챔피언 체스 엔진 Stockfish를 꺾었을 때 사용한 핵심 엔진이 바로 이것이다.
        </p>

        <div className="learn-block">
          <h3 className="learn-h3">1. 왜 필요한가? — 결정 트리의 문제</h3>
          <p className="learn-p">
            체스에서 한 수마다 평균 약 35개의 합법 수가 존재하고, 한 경기는 평균 80수(40턴) 정도로 끝난다.
            가능한 모든 경기 경로의 수는 약 <strong>10<sup>120</sup></strong> — 이를
            <strong> 섀넌 수(Shannon Number)</strong>라 부르며, 관측 가능한 우주의 원자 수(약 10<sup>80</sup>)보다도 압도적으로 많다.
            전부 계산하는 것은 불가능하다.
          </p>
          <p className="learn-p">
            전통적인 체스 엔진(Stockfish 등)은 <strong>"가능한 수를 가지치기하며 깊이 탐색"</strong>하고
            <strong> "사람이 손수 만든 평가 함수"</strong>로 점수를 매겨왔다. 강력하지만, 좋은 평가 함수를 설계하는 데 수십 년의 인간 노하우가 필요했다.
          </p>
          <p className="learn-p">
            <strong>MCTS의 발상</strong>은 다르다 — "유망해 보이는 가지만 골라서, 끝까지 시뮬레이션해보자."
            그리고 그 결과를 통계적으로 누적하면, 자연스럽게 좋은 선택지가 부각된다.
            손으로 만든 규칙 없이 게임의 승패 정보만으로도 최선의 수를 찾을 수 있다.
          </p>
        </div>

        <div className="learn-block">
          <h3 className="learn-h3">2. 알고리즘 4단계</h3>
          <p className="learn-p">
            MCTS는 아래 네 단계를 수천~수백만 번 반복하면서 결정 트리를 점진적으로 키워나간다.
          </p>
          <div className="phases-list">
            {phases.map((p) => (
              <div key={p.num} className="phase-row">
                <div className="phase-num">{p.num}</div>
                <div className="phase-body">
                  <h4 className="phase-name">{p.name}</h4>
                  <p className="phase-desc">{p.desc}</p>
                  <p className="phase-detail">{p.detail}</p>
                </div>
              </div>
            ))}
          </div>
        </div>

        <div className="learn-block">
          <h3 className="learn-h3">3. 핵심 수식 — UCB1</h3>
          <p className="learn-p">
            Selection 단계에서 어떤 자식 노드를 고를까? <strong>UCB1(Upper Confidence Bound)</strong> 공식이 답을 준다.
          </p>
          <div className="formula-box">
            <div className="formula">
              UCB1 = <span className="formula-acc">x̄<sub>i</sub></span>
              <span className="formula-op"> + </span>
              <span className="formula-acc">c</span>
              <span className="formula-frac">
                √( ln N / n<sub>i</sub> )
              </span>
            </div>
            <div className="formula-legend">
              <div><span className="formula-key">x̄ᵢ</span> = 노드 i의 평균 보상 <span className="legend-tag">(활용)</span></div>
              <div><span className="formula-key">N</span> = 부모 노드의 총 방문 횟수</div>
              <div><span className="formula-key">nᵢ</span> = 노드 i의 방문 횟수</div>
              <div><span className="formula-key">c</span> = 탐색 강도 상수 (보통 √2)</div>
            </div>
          </div>
          <p className="learn-p">
            왼쪽 항은 "지금까지 이 노드가 얼마나 좋은 결과를 냈는가"(<strong>활용</strong>), 오른쪽 항은
            "이 노드를 충분히 탐색했는가"(<strong>탐색</strong>)를 나타낸다. 적게 방문한 노드일수록
            오른쪽 항이 커져 한 번쯤 가볼 기회를 얻는다.
          </p>
        </div>

        <div className="learn-block">
          <h3 className="learn-h3">4. AlphaZero는 어떻게 이것을 발전시켰나</h3>
          <p className="learn-p">
            전통적인 MCTS의 약점은 시뮬레이션 단계의 "무작위 플레이"가 너무 단순하다는 점이다.
            DeepMind의 AlphaZero(2017)는 여기에 <strong>두 개의 신경망</strong>을 결합해 체스·쇼기·체스을 모두 정복했다:
          </p>
          <ul className="learn-ul">
            <li><strong>정책망(Policy Network)</strong> — Selection 단계에서 "어떤 수가 그럴듯한가"를 추천해, 무작위가 아닌 똑똑한 탐색이 가능하게 함</li>
            <li><strong>가치망(Value Network)</strong> — Simulation을 끝까지 돌리지 않고도 "이 국면의 승률"을 추정해 시간을 절약</li>
          </ul>
          <p className="learn-p">
            놀라운 점은 AlphaZero가 <strong>체스의 규칙만 알려준 상태에서 자기 자신과의 대국만으로 학습</strong>했다는 것이다.
            단 4시간 만에 세계 챔피언 엔진 Stockfish를 압도했고, 인간이 수백 년간 쌓아온 정석을 뛰어넘는 새로운 전략을 스스로 발견했다.
          </p>
          <p className="learn-p">
            이후 MuZero는 게임의 규칙조차 모른 채로 학습하는 단계까지 발전했고,
            이 기술은 게임을 넘어 단백질 구조 예측(AlphaFold)·로봇 제어·핵융합 플라즈마 제어 등으로 확장되었다.
          </p>
        </div>

        <div className="learn-block highlight-block">
          <h3 className="learn-h3">5. 왜 의료 의사결정에 쓸 수 있는가?</h3>
          <p className="learn-p">
            게임과 의료 치료는 놀랍게도 같은 구조를 가진다:
          </p>
          <div className="parallel-table">
            <div className="parallel-row">
              <div className="parallel-cell parallel-head">체스</div>
              <div className="parallel-cell parallel-head">유방암 치료</div>
            </div>
            <div className="parallel-row">
              <div className="parallel-cell">현재 판 상태</div>
              <div className="parallel-cell">환자 상태 (병기·수용체·림프절…)</div>
            </div>
            <div className="parallel-row">
              <div className="parallel-cell">다음에 둘 수</div>
              <div className="parallel-cell">다음 치료 선택 (수술·항암·표적…)</div>
            </div>
            <div className="parallel-row">
              <div className="parallel-cell">한 수 두면 상대가 응수</div>
              <div className="parallel-cell">치료하면 종양이 반응 (관해·진행·재발)</div>
            </div>
            <div className="parallel-row">
              <div className="parallel-cell">게임 종료 → 승/패</div>
              <div className="parallel-cell">추적 종료 → 생존·재발·독성 결과</div>
            </div>
          </div>
          <p className="learn-p">
            즉, 유방암 치료는 본질적으로 <strong>순차적 의사결정 게임</strong>이며,
            MCTS는 이 "게임"에서 표준 가이드라인이 놓칠 수 있는 더 나은 경로를 시뮬레이션으로 탐색할 수 있다.
            본 동아리의 연구는 바로 이 가능성을 학부생 수준에서 정량적으로 검증하는 것이 목표다.
          </p>
        </div>

        <div className="learn-sources">
          <p className="sources-title">참고 자료</p>
          <ul>
            <li>Browne et al. (2012) — A Survey of Monte Carlo Tree Search Methods</li>
            <li>Silver et al. (2017) — Mastering Chess and Shogi by Self-Play with a General Reinforcement Learning Algorithm (AlphaZero)</li>
            <li>Silver et al. (2018) — A general reinforcement learning algorithm that masters chess, shogi, and Go through self-play (Science)</li>
            <li>Shannon, C. E. (1950) — Programming a Computer for Playing Chess</li>
            <li>Wikipedia — Monte Carlo tree search</li>
          </ul>
        </div>
      </div>
    </div>
  )
}
