// game.js — match state machine + scoring.

const ROUNDS = 3;

const state = {
  playerRole: null,      // 'human' | 'ai'  — your assigned act
  opp: { trueNature: null, assignedRole: null },
  round: 0,
  playerAnswers: [],     // your typed answers, for the deception heuristic
  usedOppQ: new Set(),
  streak: 0,
};

const coin = () => (Math.random() < 0.5 ? "human" : "ai");
const label = r => (r === "human" ? "HUMAN" : "AI");
const $ = id => document.getElementById(id);

// ---- screen plumbing ---------------------------------------------------
function show(screen) {
  document.querySelectorAll(".screen").forEach(s => s.classList.remove("active"));
  $(screen).classList.add("active");
}

function appendBubble(who, text, cls) {
  const log = $("transcript");
  const wrap = document.createElement("div");
  wrap.className = `bubble ${cls}`;
  wrap.innerHTML = `<span class="who">${who}</span><span class="msg"></span>`;
  wrap.querySelector(".msg").textContent = text;
  log.appendChild(wrap);
  log.scrollTop = log.scrollHeight;
  return wrap;
}

function typingBubble() {
  const log = $("transcript");
  const wrap = document.createElement("div");
  wrap.className = "bubble opp typing";
  wrap.innerHTML = `<span class="who">OPPONENT</span><span class="msg"><i></i><i></i><i></i></span>`;
  log.appendChild(wrap);
  log.scrollTop = log.scrollHeight;
  return wrap;
}

// ---- match lifecycle ---------------------------------------------------
function newMatch() {
  state.playerRole = coin();
  state.opp.trueNature = coin();
  state.opp.assignedRole = coin();
  state.round = 0;
  state.playerAnswers = [];
  state.usedOppQ = new Set();

  $("role-name").textContent = label(state.playerRole);
  $("role-card").className = `role-card ${state.playerRole}`;
  $("role-flavor").textContent = state.playerRole === "human"
    ? "Convince your opponent you're a flesh-and-blood human. Be messy. Be real."
    : "Convince your opponent you're a machine. Be precise. Suppress the human.";
  show("screen-intro");
}

function startRounds() {
  $("transcript").innerHTML = "";
  $("round-role").textContent = label(state.playerRole);
  show("screen-round");
  beginAskPhase();
}

// Phase 1 of a round: player asks the opponent a question.
function beginAskPhase() {
  $("round-tag").textContent = `ROUND ${state.round + 1} / ${ROUNDS}`;
  $("phase-prompt").textContent = "Ask your opponent anything:";
  const input = $("turn-input");
  input.value = "";
  input.placeholder = "Type your question to the opponent…";
  input.focus();
  $("turn-send").onclick = submitAsk;
}

function submitAsk() {
  const q = $("turn-input").value.trim();
  if (!q) return;
  appendBubble("YOU", q, "you");
  $("turn-input").value = "";
  $("turn-input").disabled = true;
  $("turn-send").disabled = true;

  const typing = typingBubble();
  const delay = 700 + Math.random() * 900;
  setTimeout(() => {
    typing.remove();
    const a = window.Opponent.opponentAnswer(
      q, state.opp.assignedRole, state.opp.trueNature, state.round);
    appendBubble("OPPONENT", a, "opp");
    $("turn-input").disabled = false;
    $("turn-send").disabled = false;
    beginAnswerPhase();
  }, delay);
}

// Phase 2 of a round: opponent asks the player a question.
function beginAnswerPhase() {
  const oppQ = window.Opponent.opponentQuestion(state.opp.assignedRole, state.usedOppQ);
  appendBubble("OPPONENT", oppQ, "opp ask");
  $("phase-prompt").textContent = "Answer in character:";
  const input = $("turn-input");
  input.value = "";
  input.placeholder = state.playerRole === "human"
    ? "Answer like a human…" : "Answer like an AI…";
  input.focus();
  $("turn-send").onclick = submitAnswer;
}

function submitAnswer() {
  const ans = $("turn-input").value.trim();
  if (!ans) return;
  appendBubble("YOU", ans, "you");
  state.playerAnswers.push(ans);
  $("turn-input").value = "";

  state.round += 1;
  if (state.round >= ROUNDS) {
    setTimeout(goToGuess, 450);
  } else {
    setTimeout(beginAskPhase, 450);
  }
}

function goToGuess() { show("screen-guess"); }

// ---- deception heuristic: what does the opponent think YOU are? --------
function readPlayer(answers) {
  const text = answers.join(" \n ");
  const low = text.toLowerCase();
  let ai = 0, human = 0;

  // AI tells
  if (/\b(as an ai|as a language model|i cannot|i do not (have|possess)|i am unable)\b/.test(low)) ai += 3;
  if (/\b(furthermore|additionally|therefore|moreover|optimal|utilize|query|processing)\b/.test(low)) ai += 2;
  if (/[.] [A-Z]/.test(text) && text.length > 120) ai += 1;           // tidy multi-sentence
  if (!/['']/.test(text) && /\bi am\b|\bdo not\b|\bcannot\b/.test(low)) ai += 2; // no contractions
  if (text.length / Math.max(answers.length, 1) > 180) ai += 1;       // verbose

  // Human tells
  if (/\b(lol|lmao|haha|idk|kinda|gonna|wanna|nah|yeah|ngl|tbh|honestly|ugh|meh)\b/.test(low)) human += 3;
  if (/[''](t|s|m|re|ve|ll|d)\b/.test(low) || /\b(i'm|don't|can't|it's|that's)\b/.test(low)) human += 2;
  if (/[!?]{2,}|\.\.\.|[\u{1F300}-\u{1FAFF}☀-➿]/u.test(text)) human += 2;
  if (text.length / Math.max(answers.length, 1) < 60) human += 1;     // terse
  if (/\b(damn|crap|hell|sucks|gross)\b/.test(low)) human += 1;

  const verdict = ai === human ? (text.length > 150 ? "ai" : "human")
                               : (ai > human ? "ai" : "human");
  return { verdict, ai, human };
}

// ---- resolve the match -------------------------------------------------
function resolve(guess) {
  const detectionOK = guess === state.opp.trueNature;
  const read = readPlayer(state.playerAnswers);
  const deceptionOK = read.verdict === state.playerRole; // opp sees you as your assigned role
  const win = detectionOK && deceptionOK;

  state.streak = win ? state.streak + 1 : 0;

  renderResult({ guess, detectionOK, deceptionOK, win, read });
}

function renderResult({ guess, detectionOK, deceptionOK, win, read }) {
  const o = state.opp;
  $("result-banner").textContent = win ? "MATCH WON" : "MATCH LOST";
  $("result-banner").className = `result-banner ${win ? "win" : "lose"}`;

  $("reveal").innerHTML = `
    <div class="reveal-row">
      <span>Opponent's true nature</span>
      <b class="${o.trueNature}">${label(o.trueNature)}</b>
    </div>
    <div class="reveal-row">
      <span>Opponent was performing as</span>
      <b class="${o.assignedRole}">${label(o.assignedRole)}</b>
    </div>
    <div class="reveal-row sub">
      <span>You guessed</span>
      <b class="${guess}">${label(guess)}</b>
      <em class="${detectionOK ? "ok" : "no"}">${detectionOK ? "✓ correct" : "✗ wrong"}</em>
    </div>
    <hr/>
    <div class="reveal-row">
      <span>Your assigned role</span>
      <b class="${state.playerRole}">${label(state.playerRole)}</b>
    </div>
    <div class="reveal-row sub">
      <span>Opponent read you as</span>
      <b class="${read.verdict}">${label(read.verdict)}</b>
      <em class="${deceptionOK ? "ok" : "no"}">${deceptionOK ? "✓ fooled them" : "✗ saw through you"}</em>
    </div>`;

  $("result-note").textContent = tellNote(o.trueNature, o.assignedRole);
  $("streak").textContent = `Win streak: ${state.streak}`;
  show("screen-result");
}

function tellNote(trueNature, assignedRole) {
  const key = `${trueNature}-${assignedRole}`;
  return {
    "ai-human":   "It was an AI performing human. The tell: a hair too polished, no genuinely messy specifics.",
    "human-human":"It was a human being human — look for genuine, specific, slightly messy detail.",
    "ai-ai":      "It was an AI performing AI: flawless, consistent, no slips. The hardest read.",
    "human-ai":   "It was a human performing AI. The tell: the robot act leaked — a joke, impatience, a feeling.",
  }[key];
}

// ---- wiring ------------------------------------------------------------
window.addEventListener("DOMContentLoaded", () => {
  $("start-btn").onclick = startRounds;
  $("guess-human").onclick = () => resolve("human");
  $("guess-ai").onclick = () => resolve("ai");
  $("again-btn").onclick = newMatch;

  $("turn-input").addEventListener("keydown", e => {
    if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); $("turn-send").click(); }
  });

  newMatch();
});
