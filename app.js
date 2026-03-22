// ==========================================
// BoxScores - App Logic
// All data is hardcoded fake data from data.js
// No user-generated or external content is used
// ==========================================

// === ROUTING ===
function navigate(page) {
  document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-link').forEach(l => l.classList.remove('active'));

  const target = document.getElementById('page-' + page);
  if (target) {
    target.classList.add('active');
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }

  document.querySelectorAll(`.nav-link[data-page="${page}"]`).forEach(l => l.classList.add('active'));

  // Close mobile menu
  document.querySelector('.nav-links')?.classList.remove('mobile-open');

  // Trigger page-specific init
  if (page === 'scores') initScores();
  if (page === 'tools') initTools();
  if (page === 'betting') initBetting();
  if (page === 'stories') initStories();
}

// === TICKER ===
function buildTicker() {
  const ticker = document.getElementById('ticker');
  if (!ticker) return;
  const fragment = document.createDocumentFragment();

  // Build items twice for seamless loop
  for (let i = 0; i < 2; i++) {
    Object.entries(SCORES).forEach(([sport, games]) => {
      games.forEach(game => {
        const item = document.createElement('div');
        item.className = 'ticker-item';
        // #1 Ticker click navigates to scores page
        item.addEventListener('click', () => navigate('scores'));

        const status = document.createElement('span');
        status.className = `ticker-status ${game.status}`;
        status.textContent = game.status === 'live' ? 'LIVE' : game.status === 'final' ? 'FINAL' : game.time;
        item.appendChild(status);

        const teams = document.createElement('span');
        teams.className = 'ticker-teams';
        teams.textContent = `${game.away} vs ${game.home}`;
        item.appendChild(teams);

        if (game.status !== 'upcoming') {
          const score = document.createElement('span');
          score.className = 'ticker-score';
          score.textContent = `${game.awayScore} - ${game.homeScore}`;
          item.appendChild(score);
        }

        if (game.status === 'live') {
          const time = document.createElement('span');
          time.className = 'ticker-time';
          time.textContent = `${game.quarter} ${game.time}`;
          item.appendChild(time);
        }

        fragment.appendChild(item);
      });
    });
  }

  ticker.appendChild(fragment);
}

// === SCORES PAGE ===
let currentSport = 'nba';

function initScores() {
  renderScores(currentSport);
  renderStandings(currentSport);
}

// #2 Yesterday/Tomorrow date switcher (mock)
function switchDate(direction) {
  const header = document.querySelector('#page-scores .section-header');
  if (!header) return;
  const buttons = header.querySelectorAll('.btn');
  buttons.forEach(btn => {
    btn.classList.remove('btn-primary');
    btn.classList.add('btn-secondary');
  });

  if (direction === 'yesterday') {
    buttons[0].classList.remove('btn-secondary');
    buttons[0].classList.add('btn-primary');
  } else if (direction === 'tomorrow') {
    buttons[2].classList.remove('btn-secondary');
    buttons[2].classList.add('btn-primary');
  } else {
    buttons[1].classList.remove('btn-secondary');
    buttons[1].classList.add('btn-primary');
  }

  // Show a brief loading skeleton then re-render the same data (mock)
  const grid = document.getElementById('scores-grid');
  if (!grid) return;
  grid.textContent = '';
  for (let i = 0; i < 4; i++) {
    const skel = document.createElement('div');
    skel.className = 'score-card skeleton';
    skel.style.height = '160px';
    grid.appendChild(skel);
  }

  setTimeout(() => {
    renderScores(currentSport);
  }, 600);
}

function switchSport(sport, context) {
  currentSport = sport;
  const tabContainer = context || document;
  tabContainer.querySelectorAll('.sport-tab').forEach(t => t.classList.remove('active'));
  const clickedTab = tabContainer.querySelector(`.sport-tab[data-sport="${sport}"]`);
  if (clickedTab) clickedTab.classList.add('active');

  if (tabContainer.closest('#page-scores') || tabContainer.closest('.scores-section')) {
    renderScores(sport);
    renderStandings(sport);
  }
  if (tabContainer.closest('#page-betting') || tabContainer.closest('.betting-section')) {
    renderBettingOdds(sport);
  }
}

function renderScores(sport) {
  const grid = document.getElementById('scores-grid');
  if (!grid) return;

  const games = SCORES[sport] || [];
  const teamData = TEAMS[sport] || {};

  // Clear existing content
  grid.textContent = '';

  games.forEach(game => {
    const awayTeam = teamData[game.away] || { name: game.away, color: '#333' };
    const homeTeam = teamData[game.home] || { name: game.home, color: '#333' };
    const isAwayWinning = game.awayScore > game.homeScore;
    const isHomeWinning = game.homeScore > game.awayScore;

    const card = document.createElement('div');
    card.className = `score-card ${game.status === 'live' ? 'live' : ''}`;
    card.addEventListener('click', () => openBoxScore(sport, game.away, game.home));

    // Header
    const header = document.createElement('div');
    header.className = 'score-card-header';

    const sportLabel = document.createElement('span');
    sportLabel.className = 'score-card-sport';
    sportLabel.textContent = sport.toUpperCase();

    const statusBadge = document.createElement('span');
    statusBadge.className = `score-card-status ${game.status}`;
    statusBadge.textContent = game.status === 'live' ? 'LIVE' : game.status === 'final' ? 'FINAL' : game.time;

    header.appendChild(sportLabel);
    header.appendChild(statusBadge);
    card.appendChild(header);

    // Teams container
    const teamsDiv = document.createElement('div');
    teamsDiv.className = 'score-card-teams';

    [{ code: game.away, data: awayTeam, score: game.awayScore, winning: isAwayWinning },
     { code: game.home, data: homeTeam, score: game.homeScore, winning: isHomeWinning }].forEach(t => {
      const row = document.createElement('div');
      row.className = 'score-team';

      const info = document.createElement('div');
      info.className = 'score-team-info';

      const logo = document.createElement('div');
      logo.className = 'team-logo';
      logo.style.background = t.data.color || '#333';
      logo.textContent = t.code;

      const name = document.createElement('span');
      name.className = 'score-team-name';
      name.textContent = `${t.data.city || ''} ${t.data.name || t.code}`;

      info.appendChild(logo);
      info.appendChild(name);

      const scoreEl = document.createElement('span');
      scoreEl.className = `score-team-score ${t.winning && game.status !== 'upcoming' ? 'winning' : ''}`;
      scoreEl.textContent = game.status === 'upcoming' ? '-' : t.score;

      row.appendChild(info);
      row.appendChild(scoreEl);
      teamsDiv.appendChild(row);
    });

    card.appendChild(teamsDiv);

    if (game.status === 'live') {
      const footer = document.createElement('div');
      footer.className = 'score-card-footer';
      footer.textContent = `${game.quarter} ${game.time}`;
      card.appendChild(footer);
    }

    grid.appendChild(card);
  });
}

// #3 Box scores with real players and sport-specific stats
function openBoxScore(sport, away, home) {
  const modal = document.getElementById('modal-boxscore');
  const body = document.getElementById('boxscore-body');
  const teamData = TEAMS[sport] || {};
  const awayTeam = teamData[away] || { name: away };
  const homeTeam = teamData[home] || { name: home };

  body.textContent = '';

  const title = document.createElement('h3');
  title.style.marginBottom = '16px';
  title.textContent = `${awayTeam.name || away} @ ${homeTeam.name || home}`;
  body.appendChild(title);

  // Determine sport-specific stat columns
  let statColumns;
  if (sport === 'nba') {
    statColumns = ['MIN', 'PTS', 'REB', 'AST'];
  } else if (sport === 'nfl') {
    statColumns = ['COMP', 'YDS', 'TD', 'INT'];
  } else if (sport === 'mlb') {
    statColumns = ['AB', 'H', 'RBI', 'AVG'];
  } else if (sport === 'nhl') {
    statColumns = ['TOI', 'G', 'A', 'SOG'];
  } else {
    statColumns = ['MIN', 'G', 'A', 'SH'];
  }

  // Get real roster data from ROSTERS
  const sportRosters = ROSTERS[sport] || {};
  const awayRoster = sportRosters[away] || ['Player 1', 'Player 2', 'Player 3', 'Player 4', 'Player 5'];
  const homeRoster = sportRosters[home] || ['Player 1', 'Player 2', 'Player 3', 'Player 4', 'Player 5'];

  [{ team: awayTeam.name || away, color: 'var(--cyan)', roster: awayRoster },
   { team: homeTeam.name || home, color: 'var(--red)', roster: homeRoster }].forEach(side => {
    const section = document.createElement('div');
    section.style.marginBottom = '20px';

    const heading = document.createElement('h4');
    heading.style.color = side.color;
    heading.style.fontSize = '13px';
    heading.style.marginBottom = '8px';
    heading.textContent = side.team;
    section.appendChild(heading);

    const table = document.createElement('table');
    table.className = 'data-table';
    table.style.fontSize = '12px';

    const headerRow = document.createElement('tr');
    ['Player'].concat(statColumns).forEach(col => {
      const th = document.createElement('th');
      th.textContent = col;
      headerRow.appendChild(th);
    });
    table.appendChild(headerRow);

    side.roster.forEach(playerName => {
      const row = document.createElement('tr');

      // Generate sport-specific fake stats
      let statsArr;
      if (sport === 'nba') {
        const mins = Math.floor(Math.random() * 20) + 15;
        const pts = Math.floor(Math.random() * 25) + 2;
        const reb = Math.floor(Math.random() * 10);
        const ast = Math.floor(Math.random() * 8);
        statsArr = [mins, pts, reb, ast];
      } else if (sport === 'nfl') {
        const comp = Math.floor(Math.random() * 20) + 5;
        const yds = Math.floor(Math.random() * 200) + 20;
        const td = Math.floor(Math.random() * 3);
        const intStat = Math.floor(Math.random() * 2);
        statsArr = [comp, yds, td, intStat];
      } else if (sport === 'mlb') {
        const ab = Math.floor(Math.random() * 4) + 2;
        const h = Math.floor(Math.random() * ab);
        const rbi = Math.floor(Math.random() * 3);
        const avg = ab > 0 ? (h / ab).toFixed(3) : '.000';
        statsArr = [ab, h, rbi, avg];
      } else if (sport === 'nhl') {
        const mins = Math.floor(Math.random() * 12) + 8;
        const toi = mins + ':' + String(Math.floor(Math.random() * 60)).padStart(2, '0');
        const g = Math.floor(Math.random() * 2);
        const a = Math.floor(Math.random() * 3);
        const sog = Math.floor(Math.random() * 6) + 1;
        statsArr = [toi, g, a, sog];
      } else {
        const mins = Math.floor(Math.random() * 70) + 20;
        const g = Math.floor(Math.random() * 2);
        const a = Math.floor(Math.random() * 2);
        const sh = Math.floor(Math.random() * 4);
        statsArr = [mins, g, a, sh];
      }

      [playerName].concat(statsArr).forEach(val => {
        const td = document.createElement('td');
        td.textContent = val;
        row.appendChild(td);
      });
      table.appendChild(row);
    });

    section.appendChild(table);
    body.appendChild(section);
  });

  modal.classList.add('active');
}

// === TOOLS PAGE ===
function initTools() {
  populateCompareDropdowns();
  populateRadarDropdown();
  drawRadarChart(false);
}

function populateCompareDropdowns() {
  const nbaPlayers = PLAYERS.filter(p => p.sport === 'nba');
  const sel1 = document.getElementById('compare-player-1');
  const sel2 = document.getElementById('compare-player-2');
  if (!sel1 || !sel2) return;

  [sel1, sel2].forEach(sel => {
    sel.textContent = '';
    const defaultOpt = document.createElement('option');
    defaultOpt.value = '';
    defaultOpt.textContent = 'Select Player...';
    sel.appendChild(defaultOpt);

    nbaPlayers.forEach(p => {
      const opt = document.createElement('option');
      opt.value = p.id;
      opt.textContent = `${p.name} (${p.team})`;
      sel.appendChild(opt);
    });
  });

  sel1.value = '1'; // LeBron
  sel2.value = '2'; // Tatum
  runComparison();
}

function runComparison() {
  const id1 = parseInt(document.getElementById('compare-player-1')?.value);
  const id2 = parseInt(document.getElementById('compare-player-2')?.value);
  const container = document.getElementById('compare-results');
  if (!container || !id1 || !id2) return;

  const p1 = PLAYERS.find(p => p.id === id1);
  const p2 = PLAYERS.find(p => p.id === id2);
  if (!p1 || !p2) return;

  const stats = [
    { label: 'PPG', key: 'ppg' },
    { label: 'RPG', key: 'rpg' },
    { label: 'APG', key: 'apg' },
    { label: 'GP', key: 'gp' },
  ];

  container.textContent = '';

  stats.forEach(stat => {
    const v1 = parseFloat(p1[stat.key]) || 0;
    const v2 = parseFloat(p2[stat.key]) || 0;
    const max = Math.max(v1, v2) || 1;
    const pct1 = (v1 / max * 100).toFixed(0);
    const pct2 = (v2 / max * 100).toFixed(0);

    const row = document.createElement('div');
    row.className = 'compare-stat-row';

    // Left side
    const leftDiv = document.createElement('div');
    const leftVal = document.createElement('span');
    leftVal.className = 'compare-stat-value left';
    leftVal.textContent = v1;
    const leftBarContainer = document.createElement('div');
    leftBarContainer.className = 'compare-bar-container';
    const leftBar = document.createElement('div');
    leftBar.className = 'compare-bar left';
    leftBar.style.width = '0';
    leftBar.dataset.width = pct1 + '%';
    leftBarContainer.appendChild(leftBar);
    leftDiv.appendChild(leftVal);
    leftDiv.appendChild(leftBarContainer);

    // Label
    const label = document.createElement('span');
    label.className = 'compare-stat-label';
    label.textContent = stat.label;

    // Right side
    const rightDiv = document.createElement('div');
    const rightVal = document.createElement('span');
    rightVal.className = 'compare-stat-value right';
    rightVal.textContent = v2;
    const rightBarContainer = document.createElement('div');
    rightBarContainer.className = 'compare-bar-container';
    const rightBar = document.createElement('div');
    rightBar.className = 'compare-bar right';
    rightBar.style.width = '0';
    rightBar.dataset.width = pct2 + '%';
    rightBarContainer.appendChild(rightBar);
    rightDiv.appendChild(rightVal);
    rightDiv.appendChild(rightBarContainer);

    row.appendChild(leftDiv);
    row.appendChild(label);
    row.appendChild(rightDiv);
    container.appendChild(row);
  });

  // Animate bars
  requestAnimationFrame(() => {
    setTimeout(() => {
      container.querySelectorAll('.compare-bar').forEach(bar => {
        bar.style.width = bar.dataset.width;
      });
    }, 50);
  });
}

// === STAT TABLE ===
let sortColumn = 'ppg';
let sortDir = 'desc';
let tableFilter = '';
let tableSportFilter = 'nba';

function renderStatTable() {
  const tbody = document.getElementById('stat-tbody');
  if (!tbody) return;

  let filtered = PLAYERS.filter(p => p.sport === tableSportFilter);

  if (tableFilter) {
    const f = tableFilter.toLowerCase();
    filtered = filtered.filter(p =>
      p.name.toLowerCase().includes(f) ||
      p.team.toLowerCase().includes(f) ||
      p.pos.toLowerCase().includes(f)
    );
  }

  filtered.sort((a, b) => {
    let va = a[sortColumn], vb = b[sortColumn];
    if (typeof va === 'string') va = parseFloat(va) || va;
    if (typeof vb === 'string') vb = parseFloat(vb) || vb;
    if (sortDir === 'asc') return va > vb ? 1 : -1;
    return va < vb ? 1 : -1;
  });

  tbody.textContent = '';

  // #10 Empty state when no results
  if (filtered.length === 0) {
    const emptyRow = document.createElement('tr');
    const emptyCell = document.createElement('td');
    emptyCell.colSpan = 8;
    const emptyDiv = document.createElement('div');
    emptyDiv.className = 'empty-state';
    emptyDiv.style.cssText = 'text-align:center; padding:40px 20px; color:var(--text-muted); font-size:14px';
    emptyDiv.textContent = 'No players found matching your search.';
    emptyCell.appendChild(emptyDiv);
    emptyRow.appendChild(emptyCell);
    tbody.appendChild(emptyRow);

    // Still update sort indicators
    document.querySelectorAll('#stat-table th').forEach(th => {
      th.classList.remove('sorted');
      const arrow = th.querySelector('.sort-arrow');
      if (arrow) arrow.textContent = '';
      if (th.dataset.col === sortColumn) {
        th.classList.add('sorted');
        if (arrow) arrow.textContent = sortDir === 'asc' ? ' \u25B2' : ' \u25BC';
      }
    });
    return;
  }

  filtered.forEach(p => {
    const row = document.createElement('tr');

    if (tableSportFilter === 'nba') {
      [p.name, p.team, p.pos, p.ppg, p.rpg, p.apg, p.fg, p.gp].forEach((val, i) => {
        const td = document.createElement('td');
        if (i === 0) { td.className = 'player-name'; td.textContent = val; }
        else if (i === 1) {
          const badge = document.createElement('span');
          badge.className = 'team-badge';
          badge.textContent = val;
          td.appendChild(badge);
        }
        else { td.textContent = val; }
        row.appendChild(td);
      });
    } else if (tableSportFilter === 'nfl') {
      [p.name, p.team, p.pos, p.passYds || p.recYds || '-', p.td || '-', p.int || p.rec || '-', p.qbr || p.targets || '-', p.gp].forEach((val, i) => {
        const td = document.createElement('td');
        if (i === 0) { td.className = 'player-name'; td.textContent = val; }
        else if (i === 1) {
          const badge = document.createElement('span');
          badge.className = 'team-badge';
          badge.textContent = val;
          td.appendChild(badge);
        }
        else { td.textContent = val; }
        row.appendChild(td);
      });
    } else if (tableSportFilter === 'mlb') {
      [p.name, p.team, p.pos, p.avg, p.hr, p.rbi, p.ops, p.gp].forEach((val, i) => {
        const td = document.createElement('td');
        if (i === 0) { td.className = 'player-name'; td.textContent = val; }
        else if (i === 1) {
          const badge = document.createElement('span');
          badge.className = 'team-badge';
          badge.textContent = val;
          td.appendChild(badge);
        }
        else { td.textContent = val; }
        row.appendChild(td);
      });
    } else if (tableSportFilter === 'nhl') {
      [p.name, p.team, p.pos, p.goals, p.assists, p.points, p.plusMinus, p.gp].forEach((val, i) => {
        const td = document.createElement('td');
        if (i === 0) { td.className = 'player-name'; td.textContent = val; }
        else if (i === 1) {
          const badge = document.createElement('span');
          badge.className = 'team-badge';
          badge.textContent = val;
          td.appendChild(badge);
        }
        else { td.textContent = val; }
        row.appendChild(td);
      });
    }

    tbody.appendChild(row);
  });

  // Update sort indicators
  document.querySelectorAll('#stat-table th').forEach(th => {
    th.classList.remove('sorted');
    const arrow = th.querySelector('.sort-arrow');
    if (arrow) arrow.textContent = '';
    if (th.dataset.col === sortColumn) {
      th.classList.add('sorted');
      if (arrow) arrow.textContent = sortDir === 'asc' ? ' \u25B2' : ' \u25BC';
    }
  });
}

function sortTable(column) {
  if (sortColumn === column) {
    sortDir = sortDir === 'asc' ? 'desc' : 'asc';
  } else {
    sortColumn = column;
    sortDir = 'desc';
  }
  renderStatTable();
}

function filterTable(value) {
  tableFilter = value;
  renderStatTable();
}

function switchTableSport(sport) {
  tableSportFilter = sport;
  document.querySelectorAll('.table-sport-tab').forEach(t => t.classList.remove('active'));
  document.querySelector(`.table-sport-tab[data-sport="${sport}"]`)?.classList.add('active');

  // Update table headers based on sport
  const thead = document.getElementById('stat-thead');
  if (thead) {
    thead.textContent = '';
    const headerRow = document.createElement('tr');
    let columns = [];

    if (sport === 'nba') columns = [['name','Player'],['team','Team'],['pos','Pos'],['ppg','PPG'],['rpg','RPG'],['apg','APG'],['fg','FG%'],['gp','GP']];
    else if (sport === 'nfl') columns = [['name','Player'],['team','Team'],['pos','Pos'],['passYds','YDS'],['td','TD'],['int','INT/REC'],['qbr','QBR/TGT'],['gp','GP']];
    else if (sport === 'mlb') columns = [['name','Player'],['team','Team'],['pos','Pos'],['avg','AVG'],['hr','HR'],['rbi','RBI'],['ops','OPS'],['gp','GP']];
    else if (sport === 'nhl') columns = [['name','Player'],['team','Team'],['pos','Pos'],['goals','G'],['assists','A'],['points','PTS'],['plusMinus','+/-'],['gp','GP']];

    columns.forEach(([col, label]) => {
      const th = document.createElement('th');
      th.dataset.col = col;
      th.addEventListener('click', () => sortTable(col));
      th.textContent = label;
      const arrow = document.createElement('span');
      arrow.className = 'sort-arrow';
      th.appendChild(arrow);
      headerRow.appendChild(th);
    });

    thead.appendChild(headerRow);
  }

  sortColumn = sport === 'nba' ? 'ppg' : sport === 'nfl' ? 'passYds' : sport === 'mlb' ? 'hr' : 'points';
  sortDir = 'desc';
  renderStatTable();
}

// === QUERY BUILDER ===

// #4 Update query stats per sport
function updateQueryStats() {
  const sportSel = document.getElementById('query-sport');
  const statSel = document.getElementById('query-stat');
  if (!sportSel || !statSel) return;

  const sport = sportSel.value;
  statSel.textContent = '';

  let options = [];
  if (sport === 'nba') {
    options = [['ppg', 'PPG'], ['rpg', 'RPG'], ['apg', 'APG'], ['gp', 'GP']];
  } else if (sport === 'nfl') {
    options = [['passYds', 'Pass YDS'], ['td', 'TD'], ['int', 'INT'], ['qbr', 'QBR']];
  } else if (sport === 'mlb') {
    options = [['hr', 'HR'], ['rbi', 'RBI'], ['avg', 'AVG'], ['ops', 'OPS']];
  } else if (sport === 'nhl') {
    options = [['goals', 'Goals'], ['assists', 'Assists'], ['points', 'Points'], ['plusMinus', '+/-']];
  }

  options.forEach(([val, label]) => {
    const opt = document.createElement('option');
    opt.value = val;
    opt.textContent = label;
    statSel.appendChild(opt);
  });
}

function runQuery() {
  const sport = document.getElementById('query-sport')?.value;
  const stat = document.getElementById('query-stat')?.value;
  const condition = document.getElementById('query-condition')?.value;
  const value = document.getElementById('query-value')?.value;
  const results = document.getElementById('query-results');
  if (!results) return;

  results.textContent = 'Running query...';
  results.style.color = 'var(--text-muted)';

  setTimeout(() => {
    let filtered = PLAYERS.filter(p => p.sport === sport);
    const numVal = parseFloat(value);

    filtered = filtered.filter(p => {
      const pVal = parseFloat(p[stat]);
      if (isNaN(pVal)) return false;
      if (condition === 'gt') return pVal > numVal;
      if (condition === 'lt') return pVal < numVal;
      if (condition === 'eq') return pVal === numVal;
      return true;
    });

    results.style.color = '';

    if (filtered.length === 0) {
      results.textContent = 'No results found. Try adjusting your query.';
      results.style.color = 'var(--yellow)';
      return;
    }

    const condSymbol = condition === 'gt' ? '>' : condition === 'lt' ? '<' : '=';
    const header = `-- Results: ${filtered.length} player(s) found\n-- Query: SELECT * FROM ${sport}_players WHERE ${stat} ${condSymbol} ${value}\n\n`;
    const rows = filtered.map(p => {
      return `${p.name.padEnd(25)} ${p.team.padEnd(5)} ${String(p[stat]).padEnd(8)}`;
    }).join('\n');

    results.textContent = header + rows;
  }, 600);
}

// === BETTING PAGE ===
let bettingSport = 'nba';
// #11 Bet slip tracking
let betSlip = [];

function initBetting() {
  // #6 Use bettingSport variable instead of hardcoded 'nba'
  renderBettingOdds(bettingSport);
  renderInjuries();
  renderPredictions();
  renderTrends();
  renderBetSlip();
  populateWinProbDropdown();
  drawWinProbChart();
}

function switchBettingSport(sport) {
  bettingSport = sport;
  document.querySelectorAll('.betting-sport-tab').forEach(t => t.classList.remove('active'));
  document.querySelector(`.betting-sport-tab[data-sport="${sport}"]`)?.classList.add('active');
  renderBettingOdds(sport);
}

function renderBettingOdds(sport) {
  const tbody = document.getElementById('odds-tbody');
  if (!tbody) return;

  tbody.textContent = '';
  const lines = BETTING[sport] || [];

  lines.forEach(line => {
    const row = document.createElement('tr');

    const gameCell = document.createElement('td');
    gameCell.textContent = `${line.away} @ ${line.home}`;
    row.appendChild(gameCell);

    const oddsEntries = [
      { val: line.spread, label: 'Spread', game: `${line.away} @ ${line.home}` },
      { val: line.ml_away, label: 'ML Away', game: `${line.away} @ ${line.home}` },
      { val: line.ml_home, label: 'ML Home', game: `${line.away} @ ${line.home}` },
      { val: `O/U ${line.ou}`, label: 'O/U', game: `${line.away} @ ${line.home}` },
    ];

    oddsEntries.forEach(entry => {
      const td = document.createElement('td');
      const span = document.createElement('span');
      span.className = 'odds-value';
      span.textContent = entry.val;

      // Check if already in bet slip
      const existingIdx = betSlip.findIndex(b => b.game === entry.game && b.label === entry.label);
      if (existingIdx !== -1) {
        span.classList.add('selected');
      }

      span.addEventListener('click', (e) => {
        e.stopPropagation();
        toggleBetSlipItem(entry.game, entry.label, entry.val, span);
      });
      td.appendChild(span);
      row.appendChild(td);
    });

    const timeCell = document.createElement('td');
    timeCell.style.color = line.time === 'LIVE' ? 'var(--red)' : 'var(--text-muted)';
    timeCell.style.fontWeight = '700';
    timeCell.style.fontSize = '11px';
    timeCell.textContent = line.time;
    row.appendChild(timeCell);

    tbody.appendChild(row);
  });
}

// #11 Bet slip functions
function toggleBetSlipItem(game, label, val, spanEl) {
  const existingIdx = betSlip.findIndex(b => b.game === game && b.label === label);
  if (existingIdx !== -1) {
    betSlip.splice(existingIdx, 1);
    if (spanEl) spanEl.classList.remove('selected');
  } else {
    betSlip.push({ game: game, label: label, odds: val });
    if (spanEl) spanEl.classList.add('selected');
  }
  renderBetSlip();
}

function removeBetSlipItem(index) {
  const removed = betSlip[index];
  betSlip.splice(index, 1);

  // Remove 'selected' class from matching odds-value in the table
  document.querySelectorAll('.odds-value.selected').forEach(el => {
    if (el.textContent === removed.odds) {
      const row = el.closest('tr');
      if (row) {
        const gameCell = row.querySelector('td');
        if (gameCell && gameCell.textContent === removed.game) {
          el.classList.remove('selected');
        }
      }
    }
  });
  renderBetSlip();
}

function renderBetSlip() {
  let panel = document.getElementById('bet-slip');
  if (!panel) {
    // Create the panel if it doesn't exist yet
    panel = document.createElement('div');
    panel.id = 'bet-slip';
    panel.className = 'bet-slip';
    document.body.appendChild(panel);
  }

  panel.textContent = '';

  if (betSlip.length === 0) {
    panel.classList.remove('active');
    return;
  }

  panel.classList.add('active');

  // Header
  const header = document.createElement('div');
  header.className = 'bet-slip-header';
  const h4 = document.createElement('h4');
  h4.textContent = 'Bet Slip';
  const badge = document.createElement('span');
  badge.className = 'bet-slip-count';
  badge.textContent = betSlip.length;
  header.appendChild(h4);
  header.appendChild(badge);
  panel.appendChild(header);

  // Body
  const body = document.createElement('div');
  body.className = 'bet-slip-body';
  betSlip.forEach((bet, idx) => {
    const item = document.createElement('div');
    item.className = 'bet-slip-item';

    const info = document.createElement('div');
    const gameLine = document.createElement('div');
    gameLine.style.cssText = 'font-size:12px; color:var(--text-secondary)';
    gameLine.textContent = bet.game;
    const oddsLine = document.createElement('div');
    oddsLine.style.cssText = 'font-size:13px; font-weight:700; color:var(--text-primary)';
    oddsLine.textContent = `${bet.label}: ${bet.odds}`;
    info.appendChild(gameLine);
    info.appendChild(oddsLine);
    item.appendChild(info);

    const removeBtn = document.createElement('button');
    removeBtn.className = 'bet-slip-remove';
    removeBtn.textContent = '\u00D7';
    removeBtn.addEventListener('click', () => removeBetSlipItem(idx));
    item.appendChild(removeBtn);

    body.appendChild(item);
  });
  panel.appendChild(body);

  // Footer
  const footer = document.createElement('div');
  footer.className = 'bet-slip-footer';

  const total = document.createElement('div');
  total.className = 'bet-slip-total';
  total.textContent = `${betSlip.length} selection(s)`;
  footer.appendChild(total);

  const placeBtn = document.createElement('button');
  placeBtn.className = 'btn btn-primary';
  placeBtn.style.cssText = 'width:100%; margin-top:8px; padding:10px; font-size:13px';
  placeBtn.textContent = 'Place Bet';
  placeBtn.addEventListener('click', () => {
    betSlip.length = 0;
    document.querySelectorAll('.odds-value.selected').forEach(el => el.classList.remove('selected'));
    renderBetSlip();
  });
  footer.appendChild(placeBtn);
  panel.appendChild(footer);
}

// #13 Injury filtering
function filterInjuries(sport) {
  // Update tab styling
  document.querySelectorAll('.injury-sport-tab').forEach(t => t.classList.remove('active'));
  const active = document.querySelector(`.injury-sport-tab[data-sport="${sport}"]`);
  if (active) active.classList.add('active');

  renderInjuries(sport);
}

function renderInjuries(sport) {
  const list = document.getElementById('injury-list');
  if (!list) return;

  list.textContent = '';

  const activeSport = sport || 'all';
  const filtered = activeSport === 'all' ? INJURIES : INJURIES.filter(inj => inj.sport === activeSport);

  filtered.forEach(inj => {
    const statusClass = inj.status.toLowerCase().replace(/\s+/g, '-');

    const item = document.createElement('div');
    item.className = `injury-item ${statusClass}`;

    const player = document.createElement('span');
    player.className = 'injury-player';
    player.textContent = inj.player;

    const team = document.createElement('span');
    team.className = 'injury-team';
    team.textContent = `${inj.team} (${inj.sport.toUpperCase()})`;

    const detail = document.createElement('span');
    detail.className = 'injury-detail';
    detail.textContent = inj.injury;

    const status = document.createElement('span');
    status.className = `injury-status ${statusClass}`;
    status.textContent = inj.status;

    const date = document.createElement('span');
    date.className = 'injury-date';
    date.textContent = inj.updated;

    item.appendChild(player);
    item.appendChild(team);
    item.appendChild(detail);
    item.appendChild(status);
    item.appendChild(date);
    list.appendChild(item);
  });
}

function renderPredictions() {
  const grid = document.getElementById('predictions-grid');
  if (!grid) return;

  grid.textContent = '';

  PREDICTIONS.forEach(pred => {
    const card = document.createElement('div');
    card.className = 'prediction-card';

    const matchup = document.createElement('div');
    matchup.className = 'prediction-matchup';

    const awayTeam = document.createElement('span');
    awayTeam.className = 'prediction-team';
    awayTeam.style.color = 'var(--cyan)';
    awayTeam.textContent = pred.away;

    const vs = document.createElement('span');
    vs.className = 'prediction-vs';
    vs.textContent = 'VS';

    const homeTeam = document.createElement('span');
    homeTeam.className = 'prediction-team';
    homeTeam.style.color = 'var(--red)';
    homeTeam.textContent = pred.home;

    matchup.appendChild(awayTeam);
    matchup.appendChild(vs);
    matchup.appendChild(homeTeam);
    card.appendChild(matchup);

    const barBg = document.createElement('div');
    barBg.className = 'prediction-bar-bg';
    const awayBar = document.createElement('div');
    awayBar.className = 'prediction-bar-fill away';
    awayBar.style.width = '0';
    awayBar.dataset.width = pred.awayWin + '%';
    const homeBar = document.createElement('div');
    homeBar.className = 'prediction-bar-fill home';
    homeBar.style.width = '0';
    homeBar.dataset.width = pred.homeWin + '%';
    barBg.appendChild(awayBar);
    barBg.appendChild(homeBar);
    card.appendChild(barBg);

    const pcts = document.createElement('div');
    pcts.className = 'prediction-percentages';
    const awayPct = document.createElement('span');
    awayPct.style.color = 'var(--cyan)';
    awayPct.textContent = pred.awayWin + '%';
    const homePct = document.createElement('span');
    homePct.style.color = 'var(--red)';
    homePct.textContent = pred.homeWin + '%';
    pcts.appendChild(awayPct);
    pcts.appendChild(homePct);
    card.appendChild(pcts);

    const pick = document.createElement('div');
    pick.className = 'prediction-pick';
    pick.textContent = 'AI Pick: ';
    const strong = document.createElement('strong');
    strong.textContent = `${pred.prediction} (${pred.confidence}% confidence)`;
    pick.appendChild(strong);
    card.appendChild(pick);

    grid.appendChild(card);
  });

  // Animate bars
  setTimeout(() => {
    grid.querySelectorAll('.prediction-bar-fill').forEach(bar => {
      bar.style.width = bar.dataset.width;
    });
  }, 200);
}

function renderTrends() {
  const container = document.getElementById('trends-container');
  if (!container) return;

  container.textContent = '';

  TRENDING_BETS.forEach(trend => {
    const item = document.createElement('div');
    item.className = 'trend-item';

    const game = document.createElement('span');
    game.className = 'trend-game';
    game.textContent = trend.game;
    item.appendChild(game);

    const type = document.createElement('span');
    type.className = 'trend-type';
    type.textContent = trend.type;
    item.appendChild(type);

    const bars = document.createElement('div');
    bars.className = 'trend-bars';

    [{ label: `Public ${trend.public}%`, pct: trend.public, cls: 'public' },
     { label: `Sharp ${trend.sharp}%`, pct: trend.sharp, cls: 'sharp' }].forEach(b => {
      const wrap = document.createElement('div');
      wrap.className = 'trend-bar-wrap';

      const lbl = document.createElement('span');
      lbl.className = 'trend-bar-label';
      lbl.textContent = b.label;

      const bg = document.createElement('div');
      bg.className = 'trend-bar-bg';
      const fill = document.createElement('div');
      fill.className = `trend-bar-fill ${b.cls}`;
      fill.style.width = '0';
      fill.dataset.width = b.pct + '%';
      bg.appendChild(fill);

      wrap.appendChild(lbl);
      wrap.appendChild(bg);
      bars.appendChild(wrap);
    });

    item.appendChild(bars);

    const side = document.createElement('span');
    side.className = 'trend-side';
    side.textContent = trend.side;
    item.appendChild(side);

    container.appendChild(item);
  });

  setTimeout(() => {
    container.querySelectorAll('.trend-bar-fill').forEach(bar => {
      bar.style.width = bar.dataset.width;
    });
  }, 300);
}

// === STORIES PAGE ===
let storyFilter = 'all';

function initStories() {
  renderArticles('all');
}

function filterStories(category) {
  storyFilter = category;
  document.querySelectorAll('.category-filter').forEach(f => f.classList.remove('active'));
  document.querySelector(`.category-filter[data-cat="${category}"]`)?.classList.add('active');
  renderArticles(category);
}

// #16 Article images: use article.image for thumbnails, gradient as fallback
function renderArticles(category) {
  const grid = document.getElementById('articles-grid');
  if (!grid) return;

  const filtered = category === 'all' ? ARTICLES : ARTICLES.filter(a => a.category === category);
  const sportIconClasses = { nba: 'fa-basketball', nfl: 'fa-football', mlb: 'fa-baseball-bat-ball', nhl: 'fa-hockey-puck', mls: 'fa-futbol' };

  grid.textContent = '';

  filtered.forEach(article => {
    const card = document.createElement('div');
    card.className = 'article-card';
    card.addEventListener('click', () => openArticle(article.id));

    const thumb = document.createElement('div');
    thumb.className = 'article-thumb';
    if (article.image) {
      thumb.style.backgroundImage = 'url(' + article.image + ')';
      thumb.style.backgroundSize = 'cover';
      thumb.style.backgroundPosition = 'center';
    } else if (article.gradient) {
      thumb.style.background = article.gradient;
      const icon = document.createElement('i');
      icon.className = 'fa-solid ' + (sportIconClasses[article.sport] || 'fa-trophy');
      thumb.appendChild(icon);
    }

    const catBadge = document.createElement('span');
    catBadge.className = `article-category ${article.category}`;
    catBadge.textContent = article.category.replace('-', ' ');
    thumb.appendChild(catBadge);
    card.appendChild(thumb);

    const body = document.createElement('div');
    body.className = 'article-body';
    const title = document.createElement('h3');
    title.textContent = article.title;
    body.appendChild(title);

    const meta = document.createElement('div');
    meta.className = 'article-meta';
    const author = document.createElement('span');
    author.textContent = article.author;
    const readTime = document.createElement('span');
    readTime.textContent = `${article.readTime} read`;
    meta.appendChild(author);
    meta.appendChild(readTime);
    body.appendChild(meta);
    card.appendChild(body);

    grid.appendChild(card);
  });
}

// #8 openArticle works from any page by ensuring modal exists
function openArticle(id) {
  const article = ARTICLES.find(a => a.id === id);
  if (!article) return;

  // Ensure article modal exists regardless of current page
  let modal = document.getElementById('modal-article');
  if (!modal) {
    modal = document.createElement('div');
    modal.id = 'modal-article';
    modal.className = 'modal-overlay';
    const modalInner = document.createElement('div');
    modalInner.className = 'modal';
    modalInner.style.maxWidth = '640px';
    const closeBtn = document.createElement('button');
    closeBtn.className = 'modal-close';
    closeBtn.textContent = '\u00D7';
    closeBtn.addEventListener('click', () => closeModal('modal-article'));
    const bodyDiv = document.createElement('div');
    bodyDiv.id = 'article-body';
    modalInner.appendChild(closeBtn);
    modalInner.appendChild(bodyDiv);
    modal.appendChild(modalInner);
    document.body.appendChild(modal);
  }

  const body = document.getElementById('article-body');
  const sportIconClasses2 = { nba: 'fa-basketball', nfl: 'fa-football', mlb: 'fa-baseball-bat-ball', nhl: 'fa-hockey-puck', mls: 'fa-futbol' };

  body.textContent = '';

  // Hero image
  const heroDiv = document.createElement('div');
  if (article.image) {
    heroDiv.style.cssText = 'height:200px; border-radius:var(--radius); display:flex; align-items:center; justify-content:center; font-size:64px; margin-bottom:24px; background-size:cover; background-position:center';
    heroDiv.style.backgroundImage = 'url(' + article.image + ')';
  } else {
    heroDiv.style.cssText = `background:${article.gradient || 'var(--bg-card)'}; height:200px; border-radius:var(--radius); display:flex; align-items:center; justify-content:center; font-size:64px; margin-bottom:24px`;
    const heroIcon = document.createElement('i');
    heroIcon.className = 'fa-solid ' + (sportIconClasses2[article.sport] || 'fa-trophy');
    heroDiv.appendChild(heroIcon);
  }
  body.appendChild(heroDiv);

  // Category + title + meta
  const catSpan = document.createElement('span');
  catSpan.className = `article-category ${article.category}`;
  catSpan.style.cssText = 'position:static; margin-bottom:12px; display:inline-block';
  catSpan.textContent = article.category.replace('-', ' ');
  body.appendChild(catSpan);

  const h2 = document.createElement('h2');
  h2.style.cssText = 'font-size:24px; margin-bottom:8px';
  h2.textContent = article.title;
  body.appendChild(h2);

  const metaDiv = document.createElement('div');
  metaDiv.style.cssText = 'color:var(--text-muted); font-size:13px; margin-bottom:20px';
  metaDiv.textContent = `By ${article.author} \u00B7 ${article.readTime} read \u00B7 March 21, 2026`;
  body.appendChild(metaDiv);

  // Blurred preview paragraph with fade
  const previewWrap = document.createElement('div');
  previewWrap.style.cssText = 'position:relative; margin-bottom:0; overflow:hidden; max-height:80px';

  const previewText = document.createElement('p');
  previewText.style.cssText = 'color:var(--text-secondary); line-height:1.8; margin-bottom:16px; filter:blur(3px); user-select:none';
  previewText.textContent = 'Breaking down the numbers, the tape, and everything in between. Our analysts dove deep into the data to bring you this exclusive story. The findings are nothing short of remarkable — from the advanced metrics to the eye-test moments that define greatness in the modern era of sports.';
  previewWrap.appendChild(previewText);

  const fadeOverlay = document.createElement('div');
  fadeOverlay.style.cssText = 'position:absolute; bottom:0; left:0; right:0; height:60px; background:linear-gradient(to bottom, transparent, var(--bg-card))';
  previewWrap.appendChild(fadeOverlay);
  body.appendChild(previewWrap);

  // Members-only gate
  const gate = document.createElement('div');
  gate.style.cssText = 'text-align:center; padding:32px 24px; border:1px solid var(--border); border-radius:var(--radius); background:rgba(255,255,255,0.02); margin-top:16px';

  const lockIcon = document.createElement('div');
  lockIcon.style.cssText = 'font-size:40px; margin-bottom:16px';
  const lockI = document.createElement('i');
  lockI.className = 'fa-solid fa-lock';
  lockI.style.color = 'var(--red)';
  lockIcon.appendChild(lockI);
  gate.appendChild(lockIcon);

  const gateTitle = document.createElement('h3');
  gateTitle.style.cssText = 'font-size:20px; font-weight:800; margin-bottom:8px';
  gateTitle.textContent = 'Members Only';
  gate.appendChild(gateTitle);

  const gateDesc = document.createElement('p');
  gateDesc.style.cssText = 'color:var(--text-muted); font-size:14px; margin-bottom:24px; line-height:1.6';
  gateDesc.textContent = 'This article is exclusive to BoxScores members. Create a free account to unlock all stories, analysis, and hot takes.';
  gate.appendChild(gateDesc);

  const gateButtons = document.createElement('div');
  gateButtons.style.cssText = 'display:flex; gap:12px; justify-content:center; flex-wrap:wrap';

  const signupBtn = document.createElement('button');
  signupBtn.className = 'btn btn-primary';
  signupBtn.textContent = 'Sign Up Free';
  signupBtn.onclick = () => { closeModal('modal-article'); openModal('modal-signup'); };
  gateButtons.appendChild(signupBtn);

  const loginBtn = document.createElement('button');
  loginBtn.className = 'btn btn-secondary';
  loginBtn.textContent = 'Log In';
  loginBtn.onclick = () => { closeModal('modal-article'); openModal('modal-login'); };
  gateButtons.appendChild(loginBtn);

  gate.appendChild(gateButtons);

  const proNote = document.createElement('p');
  proNote.style.cssText = 'color:var(--text-muted); font-size:12px; margin-top:16px';
  const proNoteText = document.createTextNode('Already a member? ');
  const proNoteLink = document.createElement('a');
  proNoteLink.href = '#';
  proNoteLink.style.color = 'var(--red)';
  proNoteLink.textContent = 'Log in here';
  proNoteLink.onclick = (e) => { e.preventDefault(); closeModal('modal-article'); openModal('modal-login'); };
  proNote.appendChild(proNoteText);
  proNote.appendChild(proNoteLink);
  gate.appendChild(proNote);

  body.appendChild(gate);

  modal.classList.add('active');
}

// === MODALS ===
function openModal(id) {
  document.getElementById(id)?.classList.add('active');
}

function closeModal(id) {
  document.getElementById(id)?.classList.remove('active');
}

// Close modal on overlay click
document.addEventListener('click', (e) => {
  if (e.target.classList.contains('modal-overlay')) {
    e.target.classList.remove('active');
  }
});

// Close modal on Escape
document.addEventListener('keydown', (e) => {
  if (e.key === 'Escape') {
    document.querySelectorAll('.modal-overlay.active').forEach(m => m.classList.remove('active'));
  }
});

// #5 Close search dropdown when clicking outside
document.addEventListener('click', (e) => {
  const dropdown = document.getElementById('search-dropdown');
  const searchInput = document.querySelector('.nav-search');
  if (!dropdown) return;
  if (searchInput && searchInput.contains(e.target)) return;
  if (dropdown.contains(e.target)) return;
  dropdown.style.display = 'none';
});

// === ANIMATED COUNTERS ===
function animateCounters() {
  const counters = document.querySelectorAll('.stat-number');
  counters.forEach(counter => {
    const target = counter.dataset.target;
    if (!target) return;

    const num = parseInt(target.replace(/[^0-9]/g, ''));
    if (isNaN(num)) { counter.textContent = target; return; }

    const suffix = target.includes('+') ? '+' : target.includes('%') ? '%' : '';
    const duration = 2000;
    const start = Date.now();

    const tick = () => {
      const elapsed = Date.now() - start;
      const progress = Math.min(elapsed / duration, 1);
      const eased = 1 - Math.pow(1 - progress, 3);
      const current = Math.floor(eased * num);

      counter.textContent = current.toLocaleString() + suffix;

      if (progress < 1) requestAnimationFrame(tick);
    };

    tick();
  });
}

// #9 IntersectionObserver for counters
function observeCounters() {
  const statsSection = document.querySelector('.stats-counter');
  if (!statsSection) {
    // Fallback if no section found
    animateCounters();
    return;
  }

  if (!('IntersectionObserver' in window)) {
    animateCounters();
    return;
  }

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        animateCounters();
        observer.unobserve(entry.target);
      }
    });
  }, { threshold: 0.2 });

  observer.observe(statsSection);
}

// === SEARCH ===
function handleSearch(value) {
  const dropdown = document.getElementById('search-dropdown');
  if (!dropdown) return;

  if (!value || value.length < 2) {
    dropdown.style.display = 'none';
    return;
  }

  const query = value.toLowerCase();
  const results = PLAYERS.filter(p =>
    p.name.toLowerCase().includes(query) ||
    p.team.toLowerCase().includes(query)
  ).slice(0, 5);

  if (results.length === 0) {
    dropdown.style.display = 'none';
    return;
  }

  dropdown.textContent = '';
  results.forEach(p => {
    const item = document.createElement('div');
    item.className = 'search-result';
    item.addEventListener('click', () => {
      navigate('tools');
      dropdown.style.display = 'none';
    });

    const nameEl = document.createElement('strong');
    nameEl.textContent = p.name;
    item.appendChild(nameEl);

    const meta = document.createElement('span');
    meta.style.cssText = 'color:var(--text-muted); font-size:12px';
    meta.textContent = ` ${p.team} \u00B7 ${p.sport.toUpperCase()}`;
    item.appendChild(meta);

    dropdown.appendChild(item);
  });
  dropdown.style.display = 'block';
}

// === HAMBURGER MENU ===
function toggleMenu() {
  document.querySelector('.nav-links')?.classList.toggle('mobile-open');
}

// #12 Sign-up modal success
function handleSignup() {
  const modal = document.getElementById('modal-signup');
  if (!modal) return;
  const inner = modal.querySelector('.modal');
  if (!inner) return;

  // Show loading state
  const btn = inner.querySelector('.btn-primary');
  if (btn) {
    btn.disabled = true;
    btn.textContent = 'Creating Account...';
  }

  setTimeout(() => {
    // Swap modal content to success message
    inner.textContent = '';

    const closeBtn = document.createElement('button');
    closeBtn.className = 'modal-close';
    closeBtn.textContent = '\u00D7';
    closeBtn.addEventListener('click', () => {
      closeModal('modal-signup');
      resetSignupModal();
    });
    inner.appendChild(closeBtn);

    const checkmark = document.createElement('div');
    checkmark.style.cssText = 'text-align:center; font-size:64px; margin-bottom:16px; margin-top:24px';
    checkmark.textContent = '\u2713';
    inner.appendChild(checkmark);

    const heading = document.createElement('h2');
    heading.style.cssText = 'text-align:center; margin-bottom:12px';
    heading.textContent = 'Welcome to BoxScores!';
    inner.appendChild(heading);

    const msg = document.createElement('p');
    msg.style.cssText = 'text-align:center; color:var(--text-secondary); margin-bottom:24px';
    msg.textContent = 'Your free account has been created. Start exploring live scores, data tools, and stories right away.';
    inner.appendChild(msg);

    const doneBtn = document.createElement('button');
    doneBtn.className = 'btn btn-primary';
    doneBtn.style.width = '100%';
    doneBtn.textContent = 'Get Started';
    doneBtn.addEventListener('click', () => {
      closeModal('modal-signup');
      resetSignupModal();
    });
    inner.appendChild(doneBtn);
  }, 1000);
}

function resetSignupModal() {
  // Restore original signup modal content after closing
  const modal = document.getElementById('modal-signup');
  if (!modal) return;
  const inner = modal.querySelector('.modal');
  if (!inner) return;

  inner.textContent = '';

  const closeBtn = document.createElement('button');
  closeBtn.className = 'modal-close';
  closeBtn.textContent = '\u00D7';
  closeBtn.addEventListener('click', () => closeModal('modal-signup'));
  inner.appendChild(closeBtn);

  const h2 = document.createElement('h2');
  h2.textContent = 'Join BoxScores';
  inner.appendChild(h2);

  const p = document.createElement('p');
  p.textContent = 'Create your free account and get access to live scores, stories, and basic tools.';
  inner.appendChild(p);

  const nameInput = document.createElement('input');
  nameInput.type = 'text';
  nameInput.className = 'modal-input';
  nameInput.placeholder = 'Full Name';
  inner.appendChild(nameInput);

  const emailInput = document.createElement('input');
  emailInput.type = 'email';
  emailInput.className = 'modal-input';
  emailInput.placeholder = 'Email Address';
  inner.appendChild(emailInput);

  const passInput = document.createElement('input');
  passInput.type = 'password';
  passInput.className = 'modal-input';
  passInput.placeholder = 'Create Password';
  inner.appendChild(passInput);

  const createBtn = document.createElement('button');
  createBtn.className = 'btn btn-primary';
  createBtn.style.width = '100%';
  createBtn.textContent = 'Create Free Account';
  createBtn.addEventListener('click', () => handleSignup());
  inner.appendChild(createBtn);

  const loginP = document.createElement('p');
  loginP.style.cssText = 'text-align:center; margin-top:16px; font-size:13px; color:var(--text-muted)';
  loginP.textContent = 'Already have an account? ';
  const loginLink = document.createElement('a');
  loginLink.href = '#';
  loginLink.textContent = 'Log In';
  loginLink.addEventListener('click', (e) => {
    e.preventDefault();
    closeModal('modal-signup');
    openModal('modal-login');
  });
  loginP.appendChild(loginLink);
  inner.appendChild(loginP);
}

// #14 Standings
function renderStandings(sport) {
  const container = document.getElementById('standings-container');
  if (!container) return;

  const data = STANDINGS[sport] || [];
  const teamData = TEAMS[sport] || {};

  container.textContent = '';

  if (data.length === 0) {
    const empty = document.createElement('div');
    empty.style.cssText = 'text-align:center; padding:20px; color:var(--text-muted)';
    empty.textContent = 'No standings data available.';
    container.appendChild(empty);
    return;
  }

  const table = document.createElement('table');
  table.className = 'data-table';

  const thead = document.createElement('thead');
  const headerRow = document.createElement('tr');
  ['Rank', 'Team', 'W', 'L', 'PCT', 'GB', 'Streak'].forEach(col => {
    const th = document.createElement('th');
    th.textContent = col;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);
  table.appendChild(thead);

  const tbody = document.createElement('tbody');
  data.forEach((row, idx) => {
    const tr = document.createElement('tr');
    const team = teamData[row.team] || { name: row.team, color: '#333' };

    // Rank
    const rankTd = document.createElement('td');
    rankTd.style.fontWeight = '700';
    rankTd.textContent = idx + 1;
    tr.appendChild(rankTd);

    // Team with color circle
    const teamTd = document.createElement('td');
    const teamWrap = document.createElement('div');
    teamWrap.style.cssText = 'display:flex; align-items:center; gap:8px';
    const colorDot = document.createElement('div');
    colorDot.style.cssText = `width:12px; height:12px; border-radius:50%; background:${team.color || '#333'}; flex-shrink:0`;
    const teamName = document.createElement('span');
    teamName.className = 'player-name';
    teamName.textContent = `${team.city || ''} ${team.name || row.team}`;
    teamWrap.appendChild(colorDot);
    teamWrap.appendChild(teamName);
    teamTd.appendChild(teamWrap);
    tr.appendChild(teamTd);

    // W
    const wTd = document.createElement('td');
    wTd.textContent = row.w;
    tr.appendChild(wTd);

    // L
    const lTd = document.createElement('td');
    lTd.textContent = row.l;
    tr.appendChild(lTd);

    // PCT
    const pctTd = document.createElement('td');
    pctTd.textContent = row.pct;
    tr.appendChild(pctTd);

    // GB
    const gbTd = document.createElement('td');
    gbTd.textContent = row.gb;
    tr.appendChild(gbTd);

    // Streak (colored W/L badge)
    const streakTd = document.createElement('td');
    const streakBadge = document.createElement('span');
    streakBadge.style.cssText = 'padding:2px 8px; border-radius:4px; font-size:11px; font-weight:700';
    if (row.streak.startsWith('W')) {
      streakBadge.style.background = 'rgba(0,200,100,0.15)';
      streakBadge.style.color = 'var(--green)';
    } else {
      streakBadge.style.background = 'rgba(255,62,62,0.15)';
      streakBadge.style.color = 'var(--red)';
    }
    streakBadge.textContent = row.streak;
    streakTd.appendChild(streakBadge);
    tr.appendChild(streakTd);

    tbody.appendChild(tr);
  });

  table.appendChild(tbody);
  container.appendChild(table);
}

// #15 Video section - 3D Carousel
let videoFilter = 'all';
let carouselIndex = 0;
let carouselVideos = [];
let carouselAutoTimer = null;
let carouselHovered = false;

function filterVideos(sport) {
  videoFilter = sport;
  document.querySelectorAll('#video-filters .category-filter').forEach(f => f.classList.remove('active'));
  const active = document.querySelector(`#video-filters .category-filter[data-sport="${sport}"]`);
  if (active) active.classList.add('active');
  carouselIndex = 0;
  renderVideos(sport);
}

function renderVideos(sport) {
  const scene = document.getElementById('carousel-scene');
  const dotsContainer = document.getElementById('carousel-dots');
  if (!scene) return;

  scene.textContent = '';
  if (dotsContainer) dotsContainer.textContent = '';

  const activeSport = sport || videoFilter || 'all';
  carouselVideos = activeSport === 'all' ? VIDEOS.slice() : VIDEOS.filter(v => v.sport === activeSport);

  if (carouselVideos.length === 0) return;
  if (carouselIndex >= carouselVideos.length) carouselIndex = 0;

  carouselVideos.forEach(function(video, i) {
    var card = document.createElement('div');
    card.className = 'carousel-card';
    card.setAttribute('data-carousel-index', i);
    card.addEventListener('click', function() {
      if (i === carouselIndex) {
        openVideoModal(video);
      } else {
        carouselIndex = i;
        updateCarouselPositions();
        resetCarouselAutoPlay();
      }
    });

    // Poster
    var posterDiv = document.createElement('div');
    posterDiv.className = 'carousel-poster';
    posterDiv.style.backgroundImage = 'url(' + video.poster + ')';

    // Play button overlay
    var playBtn = document.createElement('div');
    playBtn.className = 'carousel-play-btn';
    posterDiv.appendChild(playBtn);

    // Duration badge
    var durBadge = document.createElement('span');
    durBadge.className = 'carousel-duration';
    durBadge.textContent = video.duration;
    posterDiv.appendChild(durBadge);

    card.appendChild(posterDiv);

    // Info section
    var info = document.createElement('div');
    info.className = 'carousel-info';
    var title = document.createElement('h3');
    title.textContent = video.title;
    info.appendChild(title);
    var sportTag = document.createElement('span');
    sportTag.className = 'sport-tag';
    sportTag.textContent = video.sport.toUpperCase();
    info.appendChild(sportTag);
    card.appendChild(info);

    scene.appendChild(card);
  });

  // Build dot indicators
  if (dotsContainer) {
    carouselVideos.forEach(function(_v, i) {
      var dot = document.createElement('button');
      dot.className = 'carousel-dot';
      if (i === carouselIndex) dot.classList.add('active');
      dot.addEventListener('click', function() {
        carouselIndex = i;
        updateCarouselPositions();
        resetCarouselAutoPlay();
      });
      dotsContainer.appendChild(dot);
    });
  }

  updateCarouselPositions();
  startCarouselAutoPlay();
}

function updateCarouselPositions() {
  var cards = document.querySelectorAll('.carousel-card');
  var dots = document.querySelectorAll('.carousel-dot');
  var total = carouselVideos.length;

  cards.forEach(function(card) {
    var i = parseInt(card.getAttribute('data-carousel-index'), 10);
    // Remove all position classes
    card.classList.remove('active', 'pos-left-1', 'pos-right-1', 'pos-left-2', 'pos-right-2');

    var diff = i - carouselIndex;
    // Wrap around for circular navigation
    if (diff > total / 2) diff -= total;
    if (diff < -total / 2) diff += total;

    if (diff === 0) {
      card.classList.add('active');
    } else if (diff === -1) {
      card.classList.add('pos-left-1');
    } else if (diff === 1) {
      card.classList.add('pos-right-1');
    } else if (diff === -2) {
      card.classList.add('pos-left-2');
    } else if (diff === 2) {
      card.classList.add('pos-right-2');
    }
    // Cards further away remain hidden (opacity: 0 by default)
  });

  dots.forEach(function(dot, i) {
    if (i === carouselIndex) {
      dot.classList.add('active');
    } else {
      dot.classList.remove('active');
    }
  });
}

function carouselNext() {
  carouselIndex = (carouselIndex + 1) % carouselVideos.length;
  updateCarouselPositions();
}

function carouselPrev() {
  carouselIndex = (carouselIndex - 1 + carouselVideos.length) % carouselVideos.length;
  updateCarouselPositions();
}

function startCarouselAutoPlay() {
  stopCarouselAutoPlay();
  carouselAutoTimer = setInterval(function() {
    if (!carouselHovered && carouselVideos.length > 1) {
      carouselNext();
    }
  }, 5000);
}

function stopCarouselAutoPlay() {
  if (carouselAutoTimer) {
    clearInterval(carouselAutoTimer);
    carouselAutoTimer = null;
  }
}

function resetCarouselAutoPlay() {
  stopCarouselAutoPlay();
  startCarouselAutoPlay();
}

function initCarouselControls() {
  var prevBtn = document.getElementById('carousel-prev');
  var nextBtn = document.getElementById('carousel-next');
  var wrapper = document.getElementById('home-videos');

  if (prevBtn) {
    prevBtn.addEventListener('click', function() {
      carouselPrev();
      resetCarouselAutoPlay();
    });
  }
  if (nextBtn) {
    nextBtn.addEventListener('click', function() {
      carouselNext();
      resetCarouselAutoPlay();
    });
  }
  if (wrapper) {
    wrapper.addEventListener('mouseenter', function() {
      carouselHovered = true;
    });
    wrapper.addEventListener('mouseleave', function() {
      carouselHovered = false;
    });
  }
}

function openVideoModal(video) {
  // Create or reuse video modal
  let modal = document.getElementById('modal-video');
  if (!modal) {
    modal = document.createElement('div');
    modal.id = 'modal-video';
    modal.className = 'modal-overlay';
    const inner = document.createElement('div');
    inner.className = 'modal';
    inner.style.maxWidth = '720px';
    inner.style.padding = '0';
    inner.style.overflow = 'hidden';
    const closeBtn = document.createElement('button');
    closeBtn.className = 'modal-close';
    closeBtn.style.cssText = 'z-index:10; top:8px; right:8px';
    closeBtn.textContent = '\u00D7';
    closeBtn.addEventListener('click', () => {
      const vid = modal.querySelector('video');
      if (vid) vid.pause();
      const iframe = modal.querySelector('iframe');
      if (iframe) iframe.src = '';
      closeModal('modal-video');
    });
    inner.appendChild(closeBtn);
    const videoBody = document.createElement('div');
    videoBody.id = 'video-body';
    inner.appendChild(videoBody);
    modal.appendChild(inner);
    document.body.appendChild(modal);
  }

  const body = document.getElementById('video-body');
  body.textContent = '';

  if (video.youtubeId) {
    // YouTube iframe embed
    const wrapper = document.createElement('div');
    wrapper.style.cssText = 'position:relative; width:100%; padding-bottom:56.25%; background:#000';
    const iframe = document.createElement('iframe');
    iframe.src = 'https://www.youtube.com/embed/' + video.youtubeId + '?autoplay=1';
    iframe.style.cssText = 'position:absolute; top:0; left:0; width:100%; height:100%; border:none';
    iframe.allow = 'autoplay; encrypted-media';
    iframe.allowFullscreen = true;
    wrapper.appendChild(iframe);
    body.appendChild(wrapper);
  } else if (video.src) {
    // HTML5 video fallback
    const videoEl = document.createElement('video');
    videoEl.controls = true;
    videoEl.autoplay = true;
    videoEl.style.cssText = 'width:100%; display:block; background:#000';
    videoEl.poster = video.poster;
    const source = document.createElement('source');
    source.src = video.src;
    source.type = 'video/mp4';
    videoEl.appendChild(source);
    body.appendChild(videoEl);
  }

  const titleBar = document.createElement('div');
  titleBar.style.cssText = 'padding:16px; font-size:16px; font-weight:700; color:var(--text-primary)';
  titleBar.textContent = video.title;
  body.appendChild(titleBar);

  modal.classList.add('active');
}

// #17 Scroll-to-top button
function initScrollToTop() {
  let btn = document.querySelector('.scroll-top');
  if (!btn) {
    btn = document.createElement('button');
    btn.className = 'scroll-top';
    btn.textContent = '\u2191';
    btn.addEventListener('click', () => {
      window.scrollTo({ top: 0, behavior: 'smooth' });
    });
    document.body.appendChild(btn);
  }

  window.addEventListener('scroll', () => {
    if (window.scrollY > 400) {
      btn.classList.add('visible');
    } else {
      btn.classList.remove('visible');
    }
  });
}

// === WIN PROBABILITY CHART (Betting Page) ===
let winProbAnimationId = null;

function generateWinProbData(numPoints) {
  // Generate realistic-looking win probability data
  const awayProb = [50];
  const homeProb = [50];
  for (let i = 1; i < numPoints; i++) {
    // Random walk with mean reversion and momentum
    const momentum = (awayProb[i - 1] - 50) * 0.02;
    const swing = (Math.random() - 0.5) * 14 - momentum;
    let val = awayProb[i - 1] + swing;
    // Clamp to 5-95 range for realism
    val = Math.max(5, Math.min(95, val));
    awayProb.push(val);
    homeProb.push(100 - val);
  }
  return { away: awayProb, home: homeProb };
}

function populateWinProbDropdown() {
  const sel = document.getElementById('win-prob-game-select');
  if (!sel) return;
  sel.textContent = '';

  // Gather all live/final games across sports
  const sportLabels = { nba: 'NBA', nfl: 'NFL', mlb: 'MLB', nhl: 'NHL', mls: 'MLS' };
  const timeLabels = {
    nba: ['Q1', 'Q2', 'Q3', 'Q4'],
    nfl: ['Q1', 'Q2', 'Q3', 'Q4'],
    mlb: ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
    nhl: ['P1', 'P2', 'P3'],
    mls: ['1H', '2H'],
  };
  let firstValue = null;

  Object.entries(SCORES).forEach(function(entry) {
    const sport = entry[0];
    const games = entry[1];
    games.forEach(function(game, idx) {
      if (game.status === 'live' || game.status === 'final') {
        const teamA = TEAMS[sport] && TEAMS[sport][game.away] ? TEAMS[sport][game.away].name : game.away;
        const teamH = TEAMS[sport] && TEAMS[sport][game.home] ? TEAMS[sport][game.home].name : game.home;
        const opt = document.createElement('option');
        const val = sport + '-' + idx;
        opt.value = val;
        opt.textContent = sportLabels[sport] + ': ' + teamA + ' @ ' + teamH + ' (' + game.status.toUpperCase() + ')';
        opt.dataset.sport = sport;
        opt.dataset.away = game.away;
        opt.dataset.home = game.home;
        sel.appendChild(opt);
        if (!firstValue) firstValue = val;
      }
    });
  });

  if (firstValue) {
    sel.value = firstValue;
  }

  sel.addEventListener('change', function() {
    drawWinProbChart();
  });
}

function drawWinProbChart() {
  const canvas = document.getElementById('win-prob-canvas');
  const sel = document.getElementById('win-prob-game-select');
  if (!canvas || !sel || !sel.value) return;

  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  const displayW = 760;
  const displayH = 340;
  canvas.width = displayW * dpr;
  canvas.height = displayH * dpr;
  canvas.style.width = displayW + 'px';
  canvas.style.height = displayH + 'px';
  ctx.scale(dpr, dpr);

  const opt = sel.options[sel.selectedIndex];
  const sport = opt.dataset.sport;
  const awayCode = opt.dataset.away;
  const homeCode = opt.dataset.home;
  const awayName = TEAMS[sport] && TEAMS[sport][awayCode] ? TEAMS[sport][awayCode].name : awayCode;
  const homeName = TEAMS[sport] && TEAMS[sport][homeCode] ? TEAMS[sport][homeCode].name : homeCode;

  const timeLabels = {
    nba: ['Q1', 'Q2', 'Q3', 'Q4'],
    nfl: ['Q1', 'Q2', 'Q3', 'Q4'],
    mlb: ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
    nhl: ['P1', 'P2', 'P3'],
    mls: ['1H', '2H'],
  };
  const labels = timeLabels[sport] || ['Q1', 'Q2', 'Q3', 'Q4'];
  const numPoints = labels.length * 5; // 5 data points per period
  const data = generateWinProbData(numPoints);

  // Chart layout
  const padL = 50, padR = 20, padT = 40, padB = 40;
  const chartW = displayW - padL - padR;
  const chartH = displayH - padT - padB;

  // Cancel previous animation
  if (winProbAnimationId) {
    cancelAnimationFrame(winProbAnimationId);
    winProbAnimationId = null;
  }

  let progress = 0;
  const totalFrames = 80;

  function drawFrame() {
    progress++;
    const drawPoints = Math.min(numPoints, Math.ceil((progress / totalFrames) * numPoints));

    ctx.clearRect(0, 0, displayW, displayH);

    // Background
    ctx.fillStyle = '#0c0c14';
    ctx.fillRect(0, 0, displayW, displayH);

    // Grid lines
    ctx.strokeStyle = '#1a1a2e';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const y = padT + (chartH / 4) * i;
      ctx.beginPath();
      ctx.moveTo(padL, y);
      ctx.lineTo(padL + chartW, y);
      ctx.stroke();
    }

    // 50% line (dashed)
    ctx.strokeStyle = '#252540';
    ctx.lineWidth = 1;
    ctx.setLineDash([6, 4]);
    const midY = padT + chartH / 2;
    ctx.beginPath();
    ctx.moveTo(padL, midY);
    ctx.lineTo(padL + chartW, midY);
    ctx.stroke();
    ctx.setLineDash([]);

    // Y-axis labels
    ctx.fillStyle = '#555';
    ctx.font = '11px Inter, sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    var yLabels = [100, 75, 50, 25, 0];
    for (var yi = 0; yi < yLabels.length; yi++) {
      var yPos = padT + (chartH / 4) * yi;
      ctx.fillText(yLabels[yi] + '%', padL - 8, yPos);
    }

    // X-axis labels (period labels)
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillStyle = '#555';
    for (var li = 0; li < labels.length; li++) {
      var xPos = padL + ((li + 0.5) / labels.length) * chartW;
      ctx.fillText(labels[li], xPos, padT + chartH + 8);
    }

    // Period dividers
    ctx.strokeStyle = '#1a1a2e';
    ctx.lineWidth = 1;
    for (var di = 1; di < labels.length; di++) {
      var dx = padL + (di / labels.length) * chartW;
      ctx.beginPath();
      ctx.moveTo(dx, padT);
      ctx.lineTo(dx, padT + chartH);
      ctx.stroke();
    }

    // Draw lines
    function drawLine(dataArr, color) {
      if (drawPoints < 2) return;
      ctx.strokeStyle = color;
      ctx.lineWidth = 2.5;
      ctx.lineJoin = 'round';
      ctx.lineCap = 'round';
      ctx.beginPath();
      for (var pi = 0; pi < drawPoints; pi++) {
        var x = padL + (pi / (numPoints - 1)) * chartW;
        var y = padT + chartH - (dataArr[pi] / 100) * chartH;
        if (pi === 0) {
          ctx.moveTo(x, y);
        } else {
          ctx.lineTo(x, y);
        }
      }
      ctx.stroke();
    }

    drawLine(data.away, '#00f0ff');
    drawLine(data.home, '#ff3e3e');

    // Legend
    ctx.font = '12px Inter, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';

    // Away legend
    ctx.fillStyle = '#00f0ff';
    ctx.fillRect(padL + 10, padT - 28, 14, 3);
    ctx.fillStyle = '#ccc';
    ctx.fillText(awayName + ' (Away)', padL + 30, padT - 26);

    // Home legend
    var awayTextW = ctx.measureText(awayName + ' (Away)').width;
    ctx.fillStyle = '#ff3e3e';
    ctx.fillRect(padL + 44 + awayTextW, padT - 28, 14, 3);
    ctx.fillStyle = '#ccc';
    ctx.fillText(homeName + ' (Home)', padL + 64 + awayTextW, padT - 26);

    if (progress < totalFrames) {
      winProbAnimationId = requestAnimationFrame(drawFrame);
    } else {
      winProbAnimationId = null;
    }
  }

  drawFrame();
}

// === PLAYER RADAR CHART (Tools Page) ===
function populateRadarDropdown() {
  const sel = document.getElementById('radar-player-select');
  if (!sel) return;
  sel.textContent = '';

  const defaultOpt = document.createElement('option');
  defaultOpt.value = '';
  defaultOpt.textContent = 'Select Player...';
  sel.appendChild(defaultOpt);

  const nbaPlayers = PLAYERS.filter(function(p) { return p.sport === 'nba'; });
  nbaPlayers.forEach(function(p) {
    const opt = document.createElement('option');
    opt.value = p.id;
    opt.textContent = p.name + ' (' + p.team + ')';
    sel.appendChild(opt);
  });

  // Default to first player
  if (nbaPlayers.length > 0) {
    sel.value = nbaPlayers[0].id;
  }

  sel.addEventListener('change', function() {
    drawRadarChart(true);
  });
}

function mapPlayerToRadar(player) {
  // Map actual stats to 0-100 scale
  // ppg: 0-40 -> Scoring
  // rpg: 0-15 -> Rebounding
  // apg: 0-12 -> Playmaking
  // fg: 0.35-0.60 -> Efficiency
  // gp: 0-82 -> Durability
  var ppgNum = typeof player.ppg === 'number' ? player.ppg : parseFloat(player.ppg) || 0;
  var rpgNum = typeof player.rpg === 'number' ? player.rpg : parseFloat(player.rpg) || 0;
  var apgNum = typeof player.apg === 'number' ? player.apg : parseFloat(player.apg) || 0;
  var fgStr = String(player.fg).replace('.', '0.');
  var fgNum = parseFloat(player.fg) || 0;
  if (fgNum > 1) fgNum = fgNum / 100;
  var gpNum = typeof player.gp === 'number' ? player.gp : parseFloat(player.gp) || 0;

  return {
    Scoring: Math.min(100, (ppgNum / 40) * 100),
    Rebounding: Math.min(100, (rpgNum / 15) * 100),
    Playmaking: Math.min(100, (apgNum / 12) * 100),
    Efficiency: Math.min(100, ((fgNum - 0.35) / 0.25) * 100),
    Durability: Math.min(100, (gpNum / 82) * 100),
  };
}

let radarAnimationId = null;
let radarCurrentValues = null;

function drawRadarChart(animate) {
  const svg = document.getElementById('radar-chart-svg');
  const sel = document.getElementById('radar-player-select');
  if (!svg || !sel) return;

  const playerId = parseInt(sel.value);
  const player = PLAYERS.find(function(p) { return p.id === playerId; });
  if (!player) return;

  const values = mapPlayerToRadar(player);
  const axes = ['Scoring', 'Rebounding', 'Playmaking', 'Efficiency', 'Durability'];
  const cx = 250, cy = 220, maxR = 160;
  const angleOffset = -Math.PI / 2; // start at top

  // If animating, interpolate from previous values
  const startValues = radarCurrentValues ? Object.assign({}, radarCurrentValues) : null;
  const targetValues = values;

  if (radarAnimationId) {
    cancelAnimationFrame(radarAnimationId);
    radarAnimationId = null;
  }

  let progress = 0;
  const totalFrames = animate && startValues ? 40 : 1;

  function getAngle(i) {
    return angleOffset + (2 * Math.PI * i) / axes.length;
  }

  function getPoint(i, val) {
    var r = (val / 100) * maxR;
    var angle = getAngle(i);
    return {
      x: cx + r * Math.cos(angle),
      y: cy + r * Math.sin(angle),
    };
  }

  function renderFrame() {
    progress++;
    var t = totalFrames === 1 ? 1 : progress / totalFrames;
    // Ease out
    t = 1 - Math.pow(1 - t, 3);

    var currentVals = {};
    for (var ai = 0; ai < axes.length; ai++) {
      var axis = axes[ai];
      if (startValues && startValues[axis] !== undefined) {
        currentVals[axis] = startValues[axis] + (targetValues[axis] - startValues[axis]) * t;
      } else {
        currentVals[axis] = targetValues[axis] * t;
      }
    }

    // Clear SVG using DOM APIs
    while (svg.firstChild) {
      svg.removeChild(svg.firstChild);
    }

    var ns = 'http://www.w3.org/2000/svg';

    // Background rect
    var bgRect = document.createElementNS(ns, 'rect');
    bgRect.setAttribute('width', '500');
    bgRect.setAttribute('height', '460');
    bgRect.setAttribute('fill', '#0c0c14');
    bgRect.setAttribute('rx', '8');
    svg.appendChild(bgRect);

    // Grid rings (5 levels)
    for (var level = 1; level <= 5; level++) {
      var ringPts = [];
      for (var ri = 0; ri < axes.length; ri++) {
        var rr = (level / 5) * maxR;
        var rAngle = getAngle(ri);
        ringPts.push((cx + rr * Math.cos(rAngle)) + ',' + (cy + rr * Math.sin(rAngle)));
      }
      var ringPoly = document.createElementNS(ns, 'polygon');
      ringPoly.setAttribute('points', ringPts.join(' '));
      ringPoly.setAttribute('fill', 'none');
      ringPoly.setAttribute('stroke', '#1a1a2e');
      ringPoly.setAttribute('stroke-width', level === 5 ? '1.5' : '1');
      svg.appendChild(ringPoly);
    }

    // Axis lines
    for (var axi = 0; axi < axes.length; axi++) {
      var endPt = getPoint(axi, 100);
      var axisLine = document.createElementNS(ns, 'line');
      axisLine.setAttribute('x1', cx);
      axisLine.setAttribute('y1', cy);
      axisLine.setAttribute('x2', endPt.x);
      axisLine.setAttribute('y2', endPt.y);
      axisLine.setAttribute('stroke', '#1a1a2e');
      axisLine.setAttribute('stroke-width', '1');
      svg.appendChild(axisLine);
    }

    // Data polygon
    var dataPts = [];
    for (var di = 0; di < axes.length; di++) {
      var dPt = getPoint(di, currentVals[axes[di]]);
      dataPts.push(dPt.x + ',' + dPt.y);
    }

    var fillPoly = document.createElementNS(ns, 'polygon');
    fillPoly.setAttribute('points', dataPts.join(' '));
    fillPoly.setAttribute('fill', 'rgba(0, 240, 255, 0.2)');
    fillPoly.setAttribute('stroke', '#00f0ff');
    fillPoly.setAttribute('stroke-width', '2');
    svg.appendChild(fillPoly);

    // Data point dots
    for (var ddi = 0; ddi < axes.length; ddi++) {
      var dotPt = getPoint(ddi, currentVals[axes[ddi]]);
      var dot = document.createElementNS(ns, 'circle');
      dot.setAttribute('cx', dotPt.x);
      dot.setAttribute('cy', dotPt.y);
      dot.setAttribute('r', '4');
      dot.setAttribute('fill', '#00f0ff');
      svg.appendChild(dot);
    }

    // Axis labels
    var labelOffsets = [
      { dx: 0, dy: -18 },   // Scoring (top)
      { dx: 18, dy: 6 },    // Rebounding (right)
      { dx: 12, dy: 16 },   // Playmaking (bottom-right)
      { dx: -12, dy: 16 },  // Efficiency (bottom-left)
      { dx: -18, dy: 6 },   // Durability (left)
    ];
    for (var li = 0; li < axes.length; li++) {
      var lPt = getPoint(li, 115);
      var label = document.createElementNS(ns, 'text');
      label.setAttribute('x', lPt.x);
      label.setAttribute('y', lPt.y);
      label.setAttribute('text-anchor', 'middle');
      label.setAttribute('dominant-baseline', 'middle');
      label.setAttribute('fill', '#888');
      label.setAttribute('font-size', '12');
      label.setAttribute('font-family', 'Inter, sans-serif');
      label.setAttribute('font-weight', '600');
      label.textContent = axes[li];
      svg.appendChild(label);

      // Value below label
      var valLabel = document.createElementNS(ns, 'text');
      valLabel.setAttribute('x', lPt.x);
      valLabel.setAttribute('y', lPt.y + 15);
      valLabel.setAttribute('text-anchor', 'middle');
      valLabel.setAttribute('dominant-baseline', 'middle');
      valLabel.setAttribute('fill', '#00f0ff');
      valLabel.setAttribute('font-size', '11');
      valLabel.setAttribute('font-family', 'JetBrains Mono, monospace');
      valLabel.setAttribute('font-weight', '700');
      valLabel.textContent = Math.round(currentVals[axes[li]]);
      svg.appendChild(valLabel);
    }

    // Player name title
    var titleText = document.createElementNS(ns, 'text');
    titleText.setAttribute('x', cx);
    titleText.setAttribute('y', 24);
    titleText.setAttribute('text-anchor', 'middle');
    titleText.setAttribute('fill', '#eee');
    titleText.setAttribute('font-size', '15');
    titleText.setAttribute('font-family', 'Inter, sans-serif');
    titleText.setAttribute('font-weight', '700');
    titleText.textContent = player.name;
    svg.appendChild(titleText);

    if (progress < totalFrames) {
      radarAnimationId = requestAnimationFrame(renderFrame);
    } else {
      radarAnimationId = null;
      radarCurrentValues = Object.assign({}, targetValues);
    }
  }

  renderFrame();
}

// === INIT ===
document.addEventListener('DOMContentLoaded', () => {
  buildTicker();

  // #9 Use IntersectionObserver instead of direct call
  observeCounters();

  renderStatTable();

  // #4 Attach query sport change listener
  const querySportSel = document.getElementById('query-sport');
  if (querySportSel) {
    querySportSel.addEventListener('change', () => updateQueryStats());
  }

  // #12 Attach signup handler to modal button
  const signupBtn = document.querySelector('#modal-signup .btn-primary');
  if (signupBtn) {
    signupBtn.removeAttribute('onclick');
    signupBtn.addEventListener('click', () => handleSignup());
  }

  // #2 Attach date switcher buttons
  const scoreHeader = document.querySelector('#page-scores .section-header');
  if (scoreHeader) {
    const btns = scoreHeader.querySelectorAll('.btn');
    if (btns[0]) btns[0].addEventListener('click', () => switchDate('yesterday'));
    if (btns[1]) btns[1].addEventListener('click', () => switchDate('today'));
    if (btns[2]) btns[2].addEventListener('click', () => switchDate('tomorrow'));
  }

  // #17 Scroll to top
  initScrollToTop();

  // Render videos on home page (3D carousel)
  initCarouselControls();
  renderVideos();

  // Start on home
  navigate('home');

  // #7 Simulate live score updates with winning class update
  setInterval(() => {
    const liveCards = document.querySelectorAll('.score-card.live');
    liveCards.forEach(card => {
      const scores = card.querySelectorAll('.score-team-score');
      scores.forEach(score => {
        if (Math.random() > 0.85) {
          const current = parseInt(score.textContent);
          if (!isNaN(current) && current > 0) {
            score.textContent = current + Math.floor(Math.random() * 3);
            score.style.color = 'var(--yellow)';
            setTimeout(() => { score.style.color = ''; }, 1000);
          }
        }
      });

      // Update winning class on both team scores after potential changes
      const allScores = card.querySelectorAll('.score-team-score');
      if (allScores.length === 2) {
        const val0 = parseInt(allScores[0].textContent);
        const val1 = parseInt(allScores[1].textContent);
        if (!isNaN(val0) && !isNaN(val1)) {
          allScores[0].classList.toggle('winning', val0 > val1);
          allScores[1].classList.toggle('winning', val1 > val0);
        }
      }
    });
  }, 5000);
});
