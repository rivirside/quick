// ==========================================
// BoxScores - Fake Data
// ==========================================

const TEAMS = {
  nba: {
    LAL: { name: 'Lakers', city: 'Los Angeles', color: '#552583', accent: '#FDB927' },
    BOS: { name: 'Celtics', city: 'Boston', color: '#007A33', accent: '#BA9653' },
    GSW: { name: 'Warriors', city: 'Golden State', color: '#1D428A', accent: '#FFC72C' },
    MIA: { name: 'Heat', city: 'Miami', color: '#98002E', accent: '#F9A01B' },
    NYK: { name: 'Knicks', city: 'New York', color: '#006BB6', accent: '#F58426' },
    CHI: { name: 'Bulls', city: 'Chicago', color: '#CE1141', accent: '#000000' },
    DAL: { name: 'Mavericks', city: 'Dallas', color: '#00538C', accent: '#B8C4CA' },
    PHX: { name: 'Suns', city: 'Phoenix', color: '#1D1160', accent: '#E56020' },
    DEN: { name: 'Nuggets', city: 'Denver', color: '#0E2240', accent: '#FEC524' },
    MIL: { name: 'Bucks', city: 'Milwaukee', color: '#00471B', accent: '#EEE1C6' },
    PHI: { name: '76ers', city: 'Philadelphia', color: '#006BB6', accent: '#ED174C' },
    MEM: { name: 'Grizzlies', city: 'Memphis', color: '#5D76A9', accent: '#12173F' },
    SAS: { name: 'Spurs', city: 'San Antonio', color: '#C4CED4', accent: '#000000' },
    OKC: { name: 'Thunder', city: 'Oklahoma City', color: '#007AC1', accent: '#EF6100' },
    MIN: { name: 'Timberwolves', city: 'Minnesota', color: '#0C2340', accent: '#236192' },
    ORL: { name: 'Magic', city: 'Orlando', color: '#0077C0', accent: '#C4CED4' },
    NOP: { name: 'Pelicans', city: 'New Orleans', color: '#0C2340', accent: '#C8102E' },
    LAC: { name: 'Clippers', city: 'LA', color: '#C8102E', accent: '#1D428A' },
  },
  nfl: {
    KC: { name: 'Chiefs', city: 'Kansas City', color: '#E31837', accent: '#FFB81C' },
    SF: { name: '49ers', city: 'San Francisco', color: '#AA0000', accent: '#B3995D' },
    BUF: { name: 'Bills', city: 'Buffalo', color: '#00338D', accent: '#C60C30' },
    DAL: { name: 'Cowboys', city: 'Dallas', color: '#003594', accent: '#869397' },
    BAL: { name: 'Ravens', city: 'Baltimore', color: '#241773', accent: '#9E7C0C' },
    DET: { name: 'Lions', city: 'Detroit', color: '#0076B6', accent: '#B0B7BC' },
    PHI: { name: 'Eagles', city: 'Philadelphia', color: '#004C54', accent: '#A5ACAF' },
    MIA: { name: 'Dolphins', city: 'Miami', color: '#008E97', accent: '#FC4C02' },
  },
  mlb: {
    NYY: { name: 'Yankees', city: 'New York', color: '#003087', accent: '#E4002C' },
    LAD: { name: 'Dodgers', city: 'Los Angeles', color: '#005A9C', accent: '#EF3E42' },
    HOU: { name: 'Astros', city: 'Houston', color: '#002D62', accent: '#EB6E1F' },
    ATL: { name: 'Braves', city: 'Atlanta', color: '#CE1141', accent: '#13274F' },
    BOS: { name: 'Red Sox', city: 'Boston', color: '#BD3039', accent: '#0C2340' },
    CHC: { name: 'Cubs', city: 'Chicago', color: '#0E3386', accent: '#CC3433' },
    SD: { name: 'Padres', city: 'San Diego', color: '#2F241D', accent: '#FFC425' },
    PHI: { name: 'Phillies', city: 'Philadelphia', color: '#E81828', accent: '#002D72' },
  },
  nhl: {
    EDM: { name: 'Oilers', city: 'Edmonton', color: '#041E42', accent: '#FF4C00' },
    FLA: { name: 'Panthers', city: 'Florida', color: '#041E42', accent: '#C8102E' },
    COL: { name: 'Avalanche', city: 'Colorado', color: '#6F263D', accent: '#236192' },
    DAL: { name: 'Stars', city: 'Dallas', color: '#006847', accent: '#8F8F8C' },
    TOR: { name: 'Maple Leafs', city: 'Toronto', color: '#00205B', accent: '#FFFFFF' },
    NYR: { name: 'Rangers', city: 'New York', color: '#0038A8', accent: '#CE1126' },
    VGK: { name: 'Golden Knights', city: 'Vegas', color: '#B4975A', accent: '#333F42' },
    CAR: { name: 'Hurricanes', city: 'Carolina', color: '#CC0000', accent: '#000000' },
  },
  mls: {
    LAFC: { name: 'LAFC', city: 'Los Angeles', color: '#C39E6D', accent: '#000000' },
    MIA: { name: 'Inter Miami', city: 'Miami', color: '#F7B5CD', accent: '#231F20' },
    ATL: { name: 'Atlanta United', city: 'Atlanta', color: '#80000A', accent: '#A19060' },
    SEA: { name: 'Sounders', city: 'Seattle', color: '#005695', accent: '#658D1B' },
    NYCFC: { name: 'NYCFC', city: 'New York', color: '#6CACE4', accent: '#041E42' },
    CIN: { name: 'FC Cincinnati', city: 'Cincinnati', color: '#003087', accent: '#FE5000' },
  },
};

const SCORES = {
  nba: [
    { away: 'LAL', home: 'BOS', awayScore: 112, homeScore: 108, quarter: 'Q4', time: '2:31', status: 'live' },
    { away: 'GSW', home: 'MIA', awayScore: 98, homeScore: 101, quarter: '', time: '', status: 'final' },
    { away: 'NYK', home: 'CHI', awayScore: 88, homeScore: 91, quarter: 'Q3', time: '5:42', status: 'live' },
    { away: 'DAL', home: 'PHX', awayScore: 105, homeScore: 99, quarter: '', time: '', status: 'final' },
    { away: 'DEN', home: 'MIL', awayScore: 0, homeScore: 0, quarter: '', time: '7:30 PM', status: 'upcoming' },
    { away: 'PHI', home: 'MEM', awayScore: 0, homeScore: 0, quarter: '', time: '8:00 PM', status: 'upcoming' },
  ],
  nfl: [
    { away: 'KC', home: 'SF', awayScore: 24, homeScore: 21, quarter: '4th', time: '4:12', status: 'live' },
    { away: 'BUF', home: 'DAL', awayScore: 31, homeScore: 17, quarter: '', time: '', status: 'final' },
    { away: 'BAL', home: 'DET', awayScore: 27, homeScore: 27, quarter: '3rd', time: '8:55', status: 'live' },
    { away: 'PHI', home: 'MIA', awayScore: 0, homeScore: 0, quarter: '', time: 'SUN 1:00 PM', status: 'upcoming' },
  ],
  mlb: [
    { away: 'NYY', home: 'BOS', awayScore: 5, homeScore: 3, quarter: '', time: '', status: 'final' },
    { away: 'LAD', home: 'HOU', awayScore: 2, homeScore: 4, quarter: 'Bot 6', time: '', status: 'live' },
    { away: 'ATL', home: 'CHC', awayScore: 7, homeScore: 1, quarter: '', time: '', status: 'final' },
    { away: 'SD', home: 'PHI', awayScore: 0, homeScore: 0, quarter: '', time: '7:10 PM', status: 'upcoming' },
  ],
  nhl: [
    { away: 'EDM', home: 'FLA', awayScore: 3, homeScore: 2, quarter: '3rd', time: '11:23', status: 'live' },
    { away: 'COL', home: 'DAL', awayScore: 4, homeScore: 5, quarter: 'OT', time: '', status: 'final' },
    { away: 'TOR', home: 'NYR', awayScore: 0, homeScore: 0, quarter: '', time: '7:00 PM', status: 'upcoming' },
    { away: 'VGK', home: 'CAR', awayScore: 1, homeScore: 1, quarter: '2nd', time: '6:44', status: 'live' },
  ],
  mls: [
    { away: 'LAFC', home: 'MIA', awayScore: 2, homeScore: 2, quarter: "72'", time: '', status: 'live' },
    { away: 'ATL', home: 'SEA', awayScore: 1, homeScore: 3, quarter: '', time: '', status: 'final' },
    { away: 'NYCFC', home: 'CIN', awayScore: 0, homeScore: 0, quarter: '', time: '8:00 PM', status: 'upcoming' },
    { away: 'SEA', home: 'LAFC', awayScore: 2, homeScore: 1, quarter: '', time: '', status: 'final' },
  ],
};

const PLAYERS = [
  // NBA
  { id: 1, name: 'LeBron James', team: 'LAL', sport: 'nba', pos: 'SF', ppg: 25.4, rpg: 7.1, apg: 8.2, fg: '.510', gp: 62 },
  { id: 2, name: 'Jayson Tatum', team: 'BOS', sport: 'nba', pos: 'SF', ppg: 27.8, rpg: 8.5, apg: 4.9, fg: '.471', gp: 68 },
  { id: 3, name: 'Stephen Curry', team: 'GSW', sport: 'nba', pos: 'PG', ppg: 26.1, rpg: 4.6, apg: 5.3, fg: '.458', gp: 58 },
  { id: 4, name: 'Luka Doncic', team: 'DAL', sport: 'nba', pos: 'PG', ppg: 33.2, rpg: 9.1, apg: 9.6, fg: '.487', gp: 55 },
  { id: 5, name: 'Nikola Jokic', team: 'DEN', sport: 'nba', pos: 'C', ppg: 26.3, rpg: 12.4, apg: 9.0, fg: '.581', gp: 70 },
  { id: 6, name: 'Giannis Antetokounmpo', team: 'MIL', sport: 'nba', pos: 'PF', ppg: 30.1, rpg: 11.5, apg: 5.8, fg: '.553', gp: 60 },
  { id: 7, name: 'Victor Wembanyama', team: 'SAS', sport: 'nba', pos: 'C', ppg: 22.8, rpg: 10.2, apg: 3.7, fg: '.461', gp: 66 },
  { id: 8, name: 'Shai Gilgeous-Alexander', team: 'OKC', sport: 'nba', pos: 'SG', ppg: 31.5, rpg: 5.5, apg: 6.2, fg: '.510', gp: 65 },
  { id: 9, name: 'Anthony Edwards', team: 'MIN', sport: 'nba', pos: 'SG', ppg: 26.7, rpg: 5.8, apg: 5.1, fg: '.462', gp: 64 },
  { id: 10, name: 'Joel Embiid', team: 'PHI', sport: 'nba', pos: 'C', ppg: 34.7, rpg: 11.0, apg: 5.6, fg: '.529', gp: 42 },
  { id: 11, name: 'Jimmy Butler', team: 'MIA', sport: 'nba', pos: 'SF', ppg: 20.8, rpg: 5.3, apg: 5.0, fg: '.497', gp: 50 },
  { id: 12, name: 'Jalen Brunson', team: 'NYK', sport: 'nba', pos: 'PG', ppg: 28.4, rpg: 3.5, apg: 6.5, fg: '.479', gp: 67 },

  // NFL
  { id: 13, name: 'Patrick Mahomes', team: 'KC', sport: 'nfl', pos: 'QB', passYds: 4839, td: 38, int: 12, qbr: 98.2, gp: 17 },
  { id: 14, name: 'Josh Allen', team: 'BUF', sport: 'nfl', pos: 'QB', passYds: 4544, td: 40, int: 15, qbr: 95.1, gp: 17 },
  { id: 15, name: 'Lamar Jackson', team: 'BAL', sport: 'nfl', pos: 'QB', passYds: 3895, td: 27, int: 7, qbr: 104.5, gp: 16 },
  { id: 16, name: 'Jalen Hurts', team: 'PHI', sport: 'nfl', pos: 'QB', passYds: 3858, td: 30, int: 10, qbr: 92.4, gp: 17 },
  { id: 17, name: 'CeeDee Lamb', team: 'DAL', sport: 'nfl', pos: 'WR', recYds: 1532, td: 12, rec: 112, targets: 158, gp: 17 },
  { id: 18, name: 'Tyreek Hill', team: 'MIA', sport: 'nfl', pos: 'WR', recYds: 1480, td: 13, rec: 105, targets: 148, gp: 16 },

  // MLB
  { id: 19, name: 'Shohei Ohtani', team: 'LAD', sport: 'mlb', pos: 'DH', avg: '.310', hr: 54, rbi: 130, ops: '1.036', gp: 159 },
  { id: 20, name: 'Aaron Judge', team: 'NYY', sport: 'mlb', pos: 'RF', avg: '.267', hr: 58, rbi: 144, ops: '.993', gp: 155 },
  { id: 21, name: 'Mookie Betts', team: 'LAD', sport: 'mlb', pos: 'SS', avg: '.295', hr: 31, rbi: 97, ops: '.905', gp: 148 },
  { id: 22, name: 'Ronald Acuna Jr', team: 'ATL', sport: 'mlb', pos: 'RF', avg: '.283', hr: 41, rbi: 106, ops: '.952', gp: 152 },
  { id: 23, name: 'Corey Seager', team: 'TEX', sport: 'mlb', pos: 'SS', avg: '.278', hr: 33, rbi: 100, ops: '.882', gp: 145 },

  // NHL
  { id: 24, name: 'Connor McDavid', team: 'EDM', sport: 'nhl', pos: 'C', goals: 52, assists: 89, points: 141, plusMinus: '+32', gp: 78 },
  { id: 25, name: 'Nathan MacKinnon', team: 'COL', sport: 'nhl', pos: 'C', goals: 45, assists: 75, points: 120, plusMinus: '+28', gp: 76 },
  { id: 26, name: 'Auston Matthews', team: 'TOR', sport: 'nhl', pos: 'C', goals: 60, assists: 39, points: 99, plusMinus: '+18', gp: 74 },
  { id: 27, name: 'Aleksander Barkov', team: 'FLA', sport: 'nhl', pos: 'C', goals: 28, assists: 52, points: 80, plusMinus: '+22', gp: 72 },
];

const ARTICLES = [
  { id: 1, title: "Wemby Drops 45 in Most Historic Rookie Performance Since LeBron", category: 'highlight', sport: 'nba', readTime: '4 min', author: 'Mike Chen', image: 'https://images.unsplash.com/photo-1546519638-68e109498ffc?w=600&h=400&fit=crop' },
  { id: 2, title: "The Chiefs Dynasty: Is This the Greatest Run in NFL History?", category: 'analysis', sport: 'nfl', readTime: '8 min', author: 'Sarah Williams', image: 'https://images.unsplash.com/photo-1566577739112-5180d4bf9de7?w=600&h=400&fit=crop' },
  { id: 3, title: "Ohtani's 50/50 Season: Breaking Down Every Single Homer", category: 'analysis', sport: 'mlb', readTime: '12 min', author: 'James Park', image: 'https://images.unsplash.com/photo-1508344928928-7165b67de128?w=600&h=400&fit=crop' },
  { id: 4, title: "TOP 10: Worst Flops of the Week (You Won't Believe #3)", category: 'humor', sport: 'nba', readTime: '3 min', author: 'BoxScores Staff', image: 'https://images.unsplash.com/photo-1504450758481-7338bbe75005?w=600&h=400&fit=crop' },
  { id: 5, title: "NFL Trade Deadline: Every Deal Graded From A+ to F", category: 'breaking', sport: 'nfl', readTime: '10 min', author: 'Tom Rivera', image: 'https://images.unsplash.com/photo-1508098682722-e99c643e7f0b?w=600&h=400&fit=crop' },
  { id: 6, title: "McDavid Is Playing a Different Sport Than Everyone Else", category: 'analysis', sport: 'nhl', readTime: '6 min', author: 'Lisa Chang', image: 'https://images.unsplash.com/photo-1580692475446-c2fabbbbf835?w=600&h=400&fit=crop' },
  { id: 7, title: "Why Your Team's Front Office Is Actually Tanking (With Receipts)", category: 'hot-take', sport: 'nba', readTime: '7 min', author: 'Derek Moore', image: 'https://images.unsplash.com/photo-1574623452334-9c0eb0e24228?w=600&h=400&fit=crop' },
  { id: 8, title: "Messi Effect: MLS Attendance Numbers Are Absolutely Insane", category: 'highlight', sport: 'mls', readTime: '5 min', author: 'Carlos Ruiz', image: 'https://images.unsplash.com/photo-1431324155629-1a6deb1dec8d?w=600&h=400&fit=crop' },
  { id: 9, title: "The Worst Calls of 2026 So Far (A Visual Breakdown)", category: 'humor', sport: 'nba', readTime: '4 min', author: 'BoxScores Staff', image: 'https://images.unsplash.com/photo-1587280501635-68a0e82cd5ff?w=600&h=400&fit=crop' },
  { id: 10, title: "Judge vs Ohtani: The Definitive Statistical Comparison", category: 'analysis', sport: 'mlb', readTime: '9 min', author: 'James Park', image: 'https://images.unsplash.com/photo-1529768167801-9173d94c2a42?w=600&h=400&fit=crop' },
  { id: 11, title: "Breaking: Blockbuster 3-Team Trade Sends Shockwaves Through the League", category: 'breaking', sport: 'nba', readTime: '2 min', author: 'Sarah Williams', image: 'https://images.unsplash.com/photo-1518063319789-7217e6706b04?w=600&h=400&fit=crop' },
  { id: 12, title: "When the Mascot Does More Than the Entire Roster (Compilation)", category: 'humor', sport: 'nba', readTime: '3 min', author: 'BoxScores Staff', image: 'https://images.unsplash.com/photo-1461896836934-bd45ba8c4e0c?w=600&h=400&fit=crop' },
  { id: 13, title: "SGA Is the MVP and It's Not Even Close Anymore", category: 'hot-take', sport: 'nba', readTime: '6 min', author: 'Derek Moore', image: 'https://images.unsplash.com/photo-1559692048-79a3f837883d?w=600&h=400&fit=crop' },
  { id: 14, title: "Every Stanley Cup Contender Ranked by Goaltending Depth", category: 'analysis', sport: 'nhl', readTime: '8 min', author: 'Lisa Chang', image: 'https://images.unsplash.com/photo-1515703407324-5f753afd8be8?w=600&h=400&fit=crop' },
  { id: 15, title: "The 5 Most Lopsided Trades of the Decade (So Far)", category: 'analysis', sport: 'nfl', readTime: '7 min', author: 'Tom Rivera', image: 'https://images.unsplash.com/photo-1560272564-c83b66b1ad12?w=600&h=400&fit=crop' },
  { id: 16, title: "MLS Power Rankings: Who's Actually Winning the Shield Race?", category: 'analysis', sport: 'mls', readTime: '5 min', author: 'Carlos Ruiz', image: 'https://images.unsplash.com/photo-1517466787929-bc90951d0974?w=600&h=400&fit=crop' },
  { id: 17, title: "Breaking: Star QB Suffers Season-Ending Injury in Practice", category: 'breaking', sport: 'nfl', readTime: '2 min', author: 'Sarah Williams', image: 'https://images.unsplash.com/photo-1521731978332-9e9e714bdd20?w=600&h=400&fit=crop' },
  { id: 18, title: "Fan Catches Foul Ball With One Hand While Holding Baby (Video)", category: 'humor', sport: 'mlb', readTime: '2 min', author: 'BoxScores Staff', image: 'https://images.unsplash.com/photo-1471295253337-3ceaaedca402?w=600&h=400&fit=crop' },
  { id: 19, title: "The Analytics Revolution: How One Team Used Data to Build a Dynasty", category: 'analysis', sport: 'nba', readTime: '11 min', author: 'Mike Chen', image: 'https://images.unsplash.com/photo-1543351611-58f69d7c1571?w=600&h=400&fit=crop' },
  { id: 20, title: "Why the Next CBA Will Change Everything About Free Agency", category: 'hot-take', sport: 'nba', readTime: '9 min', author: 'Derek Moore', image: 'https://images.unsplash.com/photo-1519861531473-9200262188bf?w=600&h=400&fit=crop' },
];

const BETTING = {
  nba: [
    { away: 'LAL', home: 'BOS', spread: 'BOS -3.5', ml_away: '+145', ml_home: '-170', ou: '224.5', time: 'LIVE' },
    { away: 'NYK', home: 'CHI', spread: 'NYK -5', ml_away: '-210', ml_home: '+175', ou: '218', time: 'LIVE' },
    { away: 'DEN', home: 'MIL', spread: 'MIL -2', ml_away: '+115', ml_home: '-135', ou: '231', time: '7:30 PM' },
    { away: 'PHI', home: 'MEM', spread: 'PHI -1.5', ml_away: '-125', ml_home: '+105', ou: '215.5', time: '8:00 PM' },
  ],
  nfl: [
    { away: 'KC', home: 'SF', spread: 'KC -1.5', ml_away: '-120', ml_home: '+100', ou: '47.5', time: 'LIVE' },
    { away: 'PHI', home: 'MIA', spread: 'PHI -3', ml_away: '-155', ml_home: '+130', ou: '51', time: 'SUN 1:00' },
  ],
  mlb: [
    { away: 'LAD', home: 'HOU', spread: 'LAD -1.5', ml_away: '-145', ml_home: '+125', ou: '8.5', time: 'LIVE' },
    { away: 'SD', home: 'PHI', spread: 'PHI -1.5', ml_away: '+130', ml_home: '-150', ou: '9', time: '7:10 PM' },
  ],
  nhl: [
    { away: 'EDM', home: 'FLA', spread: 'EDM -1.5', ml_away: '-130', ml_home: '+110', ou: '6.5', time: 'LIVE' },
    { away: 'TOR', home: 'NYR', spread: 'NYR -1.5', ml_away: '+115', ml_home: '-135', ou: '5.5', time: '7:00 PM' },
  ],
};

const INJURIES = [
  { player: 'Kawhi Leonard', team: 'LAC', sport: 'nba', injury: 'Knee (ACL)', status: 'Out', updated: 'Mar 20' },
  { player: 'Ja Morant', team: 'MEM', sport: 'nba', injury: 'Shoulder (Labrum)', status: 'Doubtful', updated: 'Mar 21' },
  { player: 'Derrick Henry', team: 'BAL', sport: 'nfl', injury: 'Hamstring', status: 'Questionable', updated: 'Mar 20' },
  { player: 'Tua Tagovailoa', team: 'MIA', sport: 'nfl', injury: 'Concussion Protocol', status: 'Out', updated: 'Mar 19' },
  { player: 'Mike Trout', team: 'LAA', sport: 'mlb', injury: 'Knee (Meniscus)', status: 'Out', updated: 'Mar 21' },
  { player: 'Fernando Tatis Jr', team: 'SD', sport: 'mlb', injury: 'Wrist', status: 'Questionable', updated: 'Mar 20' },
  { player: 'Gabriel Landeskog', team: 'COL', sport: 'nhl', injury: 'Knee (ACL)', status: 'Out', updated: 'Mar 18' },
  { player: 'Klay Thompson', team: 'DAL', sport: 'nba', injury: 'Ankle (Sprain)', status: 'Questionable', updated: 'Mar 21' },
  { player: 'Zion Williamson', team: 'NOP', sport: 'nba', injury: 'Hamstring', status: 'Out', updated: 'Mar 20' },
  { player: 'Nick Bosa', team: 'SF', sport: 'nfl', injury: 'Hip', status: 'Questionable', updated: 'Mar 21' },
  { player: 'Chris Sale', team: 'ATL', sport: 'mlb', injury: 'Elbow (UCL)', status: 'Out', updated: 'Mar 17' },
  { player: 'Sidney Crosby', team: 'PIT', sport: 'nhl', injury: 'Upper Body', status: 'Day-to-Day', updated: 'Mar 21' },
  { player: 'Paolo Banchero', team: 'ORL', sport: 'nba', injury: 'Oblique', status: 'Out', updated: 'Mar 19' },
  { player: 'Christian McCaffrey', team: 'SF', sport: 'nfl', injury: 'Calf', status: 'Doubtful', updated: 'Mar 20' },
  { player: 'Chet Holmgren', team: 'OKC', sport: 'nba', injury: 'Hip', status: 'Day-to-Day', updated: 'Mar 21' },
];

const PREDICTIONS = [
  { away: 'LAL', home: 'BOS', sport: 'nba', prediction: 'BOS', confidence: 62, awayWin: 38, homeWin: 62 },
  { away: 'DEN', home: 'MIL', sport: 'nba', prediction: 'DEN', confidence: 55, awayWin: 55, homeWin: 45 },
  { away: 'KC', home: 'SF', sport: 'nfl', prediction: 'KC', confidence: 58, awayWin: 58, homeWin: 42 },
  { away: 'PHI', home: 'MIA', sport: 'nfl', prediction: 'PHI', confidence: 71, awayWin: 71, homeWin: 29 },
  { away: 'SD', home: 'PHI', sport: 'mlb', prediction: 'PHI', confidence: 60, awayWin: 40, homeWin: 60 },
  { away: 'TOR', home: 'NYR', sport: 'nhl', prediction: 'NYR', confidence: 54, awayWin: 46, homeWin: 54 },
  { away: 'PHI', home: 'MEM', sport: 'nba', prediction: 'PHI', confidence: 57, awayWin: 57, homeWin: 43 },
  { away: 'LAD', home: 'HOU', sport: 'mlb', prediction: 'LAD', confidence: 63, awayWin: 63, homeWin: 37 },
  { away: 'EDM', home: 'FLA', sport: 'nhl', prediction: 'EDM', confidence: 51, awayWin: 51, homeWin: 49 },
  { away: 'VGK', home: 'CAR', sport: 'nhl', prediction: 'CAR', confidence: 56, awayWin: 44, homeWin: 56 },
];

const TRENDING_BETS = [
  { game: 'LAL @ BOS', type: 'Spread', public: 72, sharp: 34, side: 'BOS -3.5' },
  { game: 'KC @ SF', type: 'ML', public: 65, sharp: 58, side: 'KC' },
  { game: 'DEN @ MIL', type: 'O/U', public: 68, sharp: 42, side: 'Over 231' },
  { game: 'NYY @ BOS', type: 'ML', public: 55, sharp: 71, side: 'NYY' },
  { game: 'EDM @ FLA', type: 'Spread', public: 48, sharp: 62, side: 'EDM -1.5' },
];

const PRICING = [
  {
    name: 'Free',
    price: '$0',
    period: 'forever',
    features: ['Live scores (all sports)', 'Basic box scores', 'Stories & humor section', 'Limited player comparison (3/day)', 'Ad-supported'],
    cta: 'Get Started',
    popular: false,
  },
  {
    name: 'Pro',
    price: '$9.99',
    period: '/mo',
    features: ['Everything in Free', 'Unlimited data tools', 'Advanced stat explorer', 'Custom query builder', 'Ad-free experience', 'Custom dashboard layout', 'Export data to CSV'],
    cta: 'Start Free Trial',
    popular: true,
  },
  {
    name: 'Degen',
    price: '$19.99',
    period: '/mo',
    features: ["Everything in Pro", "Gambler's Paradise full access", 'AI predictions engine', 'Real-time injury alerts', 'Sharp money indicators', 'Betting trends & analytics', 'Priority data updates', 'Pause anytime'],
    cta: 'Go Full Degen',
    popular: false,
  },
];

// Rosters for box score modals (5 key players per team)
const ROSTERS = {
  nba: {
    LAL: ['LeBron James', 'Anthony Davis', 'D\'Angelo Russell', 'Austin Reaves', 'Rui Hachimura'],
    BOS: ['Jayson Tatum', 'Jaylen Brown', 'Derrick White', 'Jrue Holiday', 'Kristaps Porzingis'],
    GSW: ['Stephen Curry', 'Klay Thompson', 'Andrew Wiggins', 'Draymond Green', 'Kevon Looney'],
    MIA: ['Jimmy Butler', 'Bam Adebayo', 'Tyler Herro', 'Terry Rozier', 'Caleb Martin'],
    NYK: ['Jalen Brunson', 'Julius Randle', 'OG Anunoby', 'Donte DiVincenzo', 'Isaiah Hartenstein'],
    CHI: ['DeMar DeRozan', 'Zach LaVine', 'Coby White', 'Nikola Vucevic', 'Patrick Williams'],
    DAL: ['Luka Doncic', 'Kyrie Irving', 'PJ Washington', 'Daniel Gafford', 'Dereck Lively'],
    PHX: ['Kevin Durant', 'Devin Booker', 'Bradley Beal', 'Jusuf Nurkic', 'Grayson Allen'],
    DEN: ['Nikola Jokic', 'Jamal Murray', 'Michael Porter Jr', 'Aaron Gordon', 'Kentavious Caldwell-Pope'],
    MIL: ['Giannis Antetokounmpo', 'Damian Lillard', 'Khris Middleton', 'Brook Lopez', 'Bobby Portis'],
    PHI: ['Joel Embiid', 'Tyrese Maxey', 'Tobias Harris', 'Buddy Hield', 'De\'Anthony Melton'],
    MEM: ['Ja Morant', 'Desmond Bane', 'Jaren Jackson Jr', 'Marcus Smart', 'Santi Aldama'],
  },
  nfl: {
    KC: ['Patrick Mahomes', 'Travis Kelce', 'Rashee Rice', 'Isiah Pacheco', 'Chris Jones'],
    SF: ['Brock Purdy', 'Christian McCaffrey', 'George Kittle', 'Deebo Samuel', 'Nick Bosa'],
    BUF: ['Josh Allen', 'Stefon Diggs', 'James Cook', 'Matt Milano', 'Von Miller'],
    DAL: ['Dak Prescott', 'CeeDee Lamb', 'Micah Parsons', 'Zack Martin', 'DeMarcus Lawrence'],
    BAL: ['Lamar Jackson', 'Derrick Henry', 'Mark Andrews', 'Zay Flowers', 'Roquan Smith'],
    DET: ['Jared Goff', 'Amon-Ra St. Brown', 'Jahmyr Gibbs', 'Aidan Hutchinson', 'Penei Sewell'],
    PHI: ['Jalen Hurts', 'A.J. Brown', 'DeVonta Smith', 'Saquon Barkley', 'Haason Reddick'],
    MIA: ['Tua Tagovailoa', 'Tyreek Hill', 'Jaylen Waddle', 'De\'Von Achane', 'Jalen Ramsey'],
  },
  mlb: {
    NYY: ['Aaron Judge', 'Juan Soto', 'Anthony Volpe', 'Giancarlo Stanton', 'Gerrit Cole'],
    LAD: ['Shohei Ohtani', 'Mookie Betts', 'Freddie Freeman', 'Will Smith', 'Yoshinobu Yamamoto'],
    HOU: ['Jose Altuve', 'Yordan Alvarez', 'Kyle Tucker', 'Alex Bregman', 'Framber Valdez'],
    ATL: ['Ronald Acuna Jr', 'Matt Olson', 'Austin Riley', 'Ozzie Albies', 'Spencer Strider'],
    BOS: ['Rafael Devers', 'Masataka Yoshida', 'Jarren Duran', 'Trevor Story', 'Brayan Bello'],
    CHC: ['Cody Bellinger', 'Dansby Swanson', 'Nico Hoerner', 'Ian Happ', 'Justin Steele'],
    SD: ['Fernando Tatis Jr', 'Manny Machado', 'Xander Bogaerts', 'Ha-Seong Kim', 'Yu Darvish'],
    PHI: ['Bryce Harper', 'Trea Turner', 'Kyle Schwarber', 'J.T. Realmuto', 'Zack Wheeler'],
  },
  nhl: {
    EDM: ['Connor McDavid', 'Leon Draisaitl', 'Evan Bouchard', 'Ryan Nugent-Hopkins', 'Zach Hyman'],
    FLA: ['Aleksander Barkov', 'Matthew Tkachuk', 'Sam Reinhart', 'Carter Verhaeghe', 'Sergei Bobrovsky'],
    COL: ['Nathan MacKinnon', 'Cale Makar', 'Mikko Rantanen', 'Valeri Nichushkin', 'Devon Toews'],
    DAL: ['Jason Robertson', 'Roope Hintz', 'Miro Heiskanen', 'Joe Pavelski', 'Jake Oettinger'],
    TOR: ['Auston Matthews', 'Mitch Marner', 'William Nylander', 'John Tavares', 'Morgan Rielly'],
    NYR: ['Artemi Panarin', 'Adam Fox', 'Mika Zibanejad', 'Chris Kreider', 'Igor Shesterkin'],
    VGK: ['Jack Eichel', 'Mark Stone', 'Chandler Stephenson', 'Shea Theodore', 'Adin Hill'],
    CAR: ['Sebastian Aho', 'Andrei Svechnikov', 'Martin Necas', 'Brent Burns', 'Frederik Andersen'],
  },
  mls: {
    LAFC: ['Carlos Vela', 'Denis Bouanga', 'Kellyn Acosta', 'Jesus Murillo', 'Maxime Crepeau'],
    MIA: ['Lionel Messi', 'Luis Suarez', 'Jordi Alba', 'Sergio Busquets', 'Drake Callender'],
    ATL: ['Giorgos Giakoumakis', 'Thiago Almada', 'Caleb Wiley', 'Brad Guzan', 'Brooks Lennon'],
    SEA: ['Jordan Morris', 'Nicolas Lodeiro', 'Cristian Roldan', 'Stefan Frei', 'Albert Rusnak'],
  },
};

// Standings
const STANDINGS = {
  nba: [
    { team: 'BOS', w: 56, l: 14, pct: '.800', gb: '-', streak: 'W5' },
    { team: 'OKC', w: 52, l: 18, pct: '.743', gb: '4', streak: 'W3' },
    { team: 'DEN', w: 49, l: 21, pct: '.700', gb: '7', streak: 'L1' },
    { team: 'MIL', w: 46, l: 24, pct: '.657', gb: '10', streak: 'W2' },
    { team: 'NYK', w: 45, l: 25, pct: '.643', gb: '11', streak: 'W1' },
    { team: 'MIN', w: 44, l: 26, pct: '.629', gb: '12', streak: 'L2' },
    { team: 'DAL', w: 43, l: 27, pct: '.614', gb: '13', streak: 'W4' },
    { team: 'LAL', w: 40, l: 30, pct: '.571', gb: '16', streak: 'L1' },
    { team: 'PHX', w: 39, l: 31, pct: '.557', gb: '17', streak: 'W1' },
    { team: 'PHI', w: 38, l: 32, pct: '.543', gb: '18', streak: 'W2' },
  ],
  nfl: [
    { team: 'KC', w: 14, l: 3, pct: '.824', gb: '-', streak: 'W6' },
    { team: 'BAL', w: 13, l: 4, pct: '.765', gb: '1', streak: 'W2' },
    { team: 'SF', w: 12, l: 5, pct: '.706', gb: '2', streak: 'L1' },
    { team: 'BUF', w: 11, l: 6, pct: '.647', gb: '3', streak: 'W3' },
    { team: 'DET', w: 11, l: 6, pct: '.647', gb: '3', streak: 'W1' },
    { team: 'DAL', w: 10, l: 7, pct: '.588', gb: '4', streak: 'L2' },
    { team: 'PHI', w: 10, l: 7, pct: '.588', gb: '4', streak: 'W1' },
    { team: 'MIA', w: 9, l: 8, pct: '.529', gb: '5', streak: 'L1' },
  ],
  mlb: [
    { team: 'LAD', w: 95, l: 52, pct: '.646', gb: '-', streak: 'W3' },
    { team: 'ATL', w: 90, l: 57, pct: '.612', gb: '5', streak: 'W2' },
    { team: 'HOU', w: 87, l: 60, pct: '.592', gb: '8', streak: 'L1' },
    { team: 'NYY', w: 85, l: 62, pct: '.578', gb: '10', streak: 'W1' },
    { team: 'PHI', w: 84, l: 63, pct: '.571', gb: '11', streak: 'W4' },
    { team: 'BOS', w: 78, l: 69, pct: '.531', gb: '17', streak: 'L3' },
    { team: 'SD', w: 76, l: 71, pct: '.517', gb: '19', streak: 'L1' },
    { team: 'CHC', w: 74, l: 73, pct: '.503', gb: '21', streak: 'W1' },
  ],
  nhl: [
    { team: 'EDM', w: 48, l: 22, pct: '.686', gb: '-', streak: 'W4' },
    { team: 'FLA', w: 46, l: 24, pct: '.657', gb: '2', streak: 'W2' },
    { team: 'DAL', w: 44, l: 25, pct: '.638', gb: '3.5', streak: 'W1' },
    { team: 'COL', w: 43, l: 27, pct: '.614', gb: '5', streak: 'L1' },
    { team: 'NYR', w: 42, l: 28, pct: '.600', gb: '6', streak: 'W3' },
    { team: 'CAR', w: 41, l: 28, pct: '.594', gb: '6.5', streak: 'L2' },
    { team: 'TOR', w: 40, l: 30, pct: '.571', gb: '8', streak: 'W1' },
    { team: 'VGK', w: 39, l: 31, pct: '.557', gb: '9', streak: 'L1' },
  ],
};

// Videos (YouTube embeds)
const VIDEOS = [
  // NBA
  { id: 1, title: "NBA Top 10 Plays of the Night", sport: 'nba', duration: '4:32', poster: 'https://images.unsplash.com/photo-1546519638-68e109498ffc?w=800&h=450&fit=crop', youtubeId: 'aVIARfKJfFo' },
  { id: 2, title: "Most Unbelievable NBA Buzzer Beaters", sport: 'nba', duration: '12:47', poster: 'https://images.unsplash.com/photo-1504450758481-7338bbe75005?w=800&h=450&fit=crop', youtubeId: '1BRGjnBMkHU' },
  { id: 3, title: "Steph Curry's Greatest 3-Pointers Ever", sport: 'nba', duration: '8:15', poster: 'https://images.unsplash.com/photo-1574623452334-9c0eb0e24228?w=800&h=450&fit=crop', youtubeId: 'OYwhmKgNvEo' },
  // NFL
  { id: 4, title: "NFL Greatest Catches of All Time", sport: 'nfl', duration: '10:21', poster: 'https://images.unsplash.com/photo-1566577739112-5180d4bf9de7?w=800&h=450&fit=crop', youtubeId: 'XagKe14YDoo' },
  { id: 5, title: "Top 100 Plays of the NFL Season", sport: 'nfl', duration: '22:38', poster: 'https://images.unsplash.com/photo-1508098682722-e99c643e7f0b?w=800&h=450&fit=crop', youtubeId: 'H7sOHnCGn4I' },
  { id: 6, title: "Craziest Trick Plays in NFL History", sport: 'nfl', duration: '11:04', poster: 'https://images.unsplash.com/photo-1560272564-c83b66b1ad12?w=800&h=450&fit=crop', youtubeId: 'aIWofrMpcxA' },
  // MLB
  { id: 7, title: "Shohei Ohtani's Most Insane Moments", sport: 'mlb', duration: '9:33', poster: 'https://images.unsplash.com/photo-1508344928928-7165b67de128?w=800&h=450&fit=crop', youtubeId: 'P-jS0qxfSMI' },
  { id: 8, title: "Greatest Home Runs in MLB History", sport: 'mlb', duration: '15:20', poster: 'https://images.unsplash.com/photo-1529768167801-9173d94c2a42?w=800&h=450&fit=crop', youtubeId: 'hJmGafYJEto' },
  { id: 9, title: "Insane Baseball Plays You Won't Believe", sport: 'mlb', duration: '7:45', poster: 'https://images.unsplash.com/photo-1471295253337-3ceaaedca402?w=800&h=450&fit=crop', youtubeId: 'NQ6GfB5YwmM' },
  // NHL
  { id: 10, title: "Connor McDavid Highlights: Pure Speed", sport: 'nhl', duration: '6:12', poster: 'https://images.unsplash.com/photo-1580692475446-c2fabbbbf835?w=800&h=450&fit=crop', youtubeId: '6GZhDvNkJSg' },
  { id: 11, title: "Best NHL Goals of the Season", sport: 'nhl', duration: '14:55', poster: 'https://images.unsplash.com/photo-1515703407324-5f753afd8be8?w=800&h=450&fit=crop', youtubeId: 'oUAX2mFHyKY' },
  { id: 12, title: "Hardest Hits in Hockey History", sport: 'nhl', duration: '10:08', poster: 'https://images.unsplash.com/photo-1559692048-79a3f837883d?w=800&h=450&fit=crop', youtubeId: 'Yvzp5lmt2L4' },
  // MLS / Soccer
  { id: 13, title: "Best MLS Goals of the Season", sport: 'mls', duration: '5:48', poster: 'https://images.unsplash.com/photo-1431324155629-1a6deb1dec8d?w=800&h=450&fit=crop', youtubeId: 'aPnMIGd3Ivs' },
  { id: 14, title: "Messi's Greatest Inter Miami Moments", sport: 'mls', duration: '8:30', poster: 'https://images.unsplash.com/photo-1517466787929-bc90951d0974?w=800&h=450&fit=crop', youtubeId: 'JVg8ZATGjCQ' },
  { id: 15, title: "Most Dramatic Last-Minute Soccer Goals", sport: 'mls', duration: '11:22', poster: 'https://images.unsplash.com/photo-1543351611-58f69d7c1571?w=800&h=450&fit=crop', youtubeId: 'NVOamCVMfDA' },
];
