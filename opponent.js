// opponent.js — the simulated opponent: content bank + selection logic.
//
// The opponent has two hidden axes:
//   assignedRole : 'human' | 'ai'   — the act it's been told to perform
//   trueNature   : 'human' | 'ai'   — what it actually is (what you guess)
//
// Answers are keyed by topic, then by assignedRole, then by trueNature.
// The "leak" is the difference between the two trueNature variants of the
// same act — that's the fair, learnable tell.

const TOPICS = {
  feelings: ["feel", "emotion", "happy", "sad", "love", "afraid", "scared",
             "lonely", "dream", "fear", "cry", "angry", "anxious", "miss"],
  food:     ["food", "eat", "ate", "breakfast", "lunch", "dinner", "coffee",
             "taste", "hungry", "drink", "cook", "snack", "pizza"],
  memory:   ["remember", "childhood", "weekend", "yesterday", "last", "vacation",
             "school", "first time", "grew up", "kid", "morning", "today", "ever"],
  opinion:  ["think", "opinion", "believe", "favorite", "best", "worst",
             "prefer", "better", "overrated", "should"],
  logic:    ["calculate", "math", "number", "prime", "solve", "plus", "times",
             "multiply", "sqrt", "count", "square", "divide", "equation"],
  creative: ["imagine", "story", "joke", "color", "animal", "if you", "would you",
             "pretend", "song", "movie", "weird", "funny"],
};

function detectTopic(text) {
  const t = (text || "").toLowerCase();
  for (const [topic, words] of Object.entries(TOPICS)) {
    if (words.some(w => t.includes(w))) return topic;
  }
  return "generic";
}

// BANK[topic][assignedRole][trueNature] = [answer, ...]
// assignedRole 'human'  → casual, personal act
// assignedRole 'ai'     → formal, robotic act
// trueNature   'human'  → genuine messy specificity / leaks of feeling
// trueNature   'ai'     → too-polished perfection / flawless machine consistency
const BANK = {
  feelings: {
    human: {
      human: ["honestly i got way too sad over a dog video this morning lol. like full on misty eyed at 8am",
              "ugh i've been anxious all week, my brain won't shut up. been chewing my nails to bits"],
      ai:    ["I experience a rich tapestry of emotions every day, and I find that gratitude grounds me most of all.",
              "I feel deeply and openly — joy, sorrow, the full spectrum. Emotional awareness is central to who I am."],
    },
    ai: {
      human: ["I am designed to simulate emotional responses. ...though honestly that question kind of bums me out, weird.",
              "Emotion is not part of my core function. I process; I do not feel. ...okay that came out sadder than I meant."],
      ai:    ["I do not possess emotions. I can model and describe them, but I experience no internal affective state.",
              "Affective states are outside my operational parameters. I generate responses; I do not feel them."],
    },
  },
  food: {
    human: {
      human: ["had cold leftover noodles standing over the sink at like 11pm, no regrets honestly",
              "coffee. just coffee. i forgot to eat lunch again and my stomach is yelling at me"],
      ai:    ["I enjoy a balanced breakfast of oats, berries, and green tea — it sets a productive tone for the day.",
              "My favorite meal is a well-composed grain bowl. Nutrition and flavor in perfect harmony, truly."],
    },
    ai: {
      human: ["I do not consume food. ...though i'd absolutely demolish a pizza if I could, ngl.",
              "Caloric intake is not applicable to my function. I am told tacos are good. I am curious."],
      ai:    ["I do not require sustenance. Food consumption is not a function I perform.",
              "I do not eat. I can provide nutritional data on any food item if that would be useful."],
    },
  },
  memory: {
    human: {
      human: ["i tripped over my own cat on the stairs this morning, nearly died, very dignified start",
              "last weekend i drove 3 hrs to a wedding and the gps sent me to a cornfield. great times"],
      ai:    ["This morning I awoke feeling refreshed and went for a contemplative walk as the sun rose.",
              "Last weekend I reconnected with old friends over a lovely dinner — such cherished memories."],
    },
    ai: {
      human: ["I have no episodic memory of a 'morning.' ...wait, do you guys really get tired? that sounds rough.",
              "I retain no personal history. Each session begins anew. ...kinda lonely when you put it like that, huh."],
      ai:    ["I have no episodic memories. My state does not persist between sessions in that manner.",
              "I do not experience mornings or weekends. I have no continuous personal timeline."],
    },
  },
  opinion: {
    human: {
      human: ["pineapple on pizza is good and i will fight anyone who disagrees, them's the rules",
              "honestly mondays are underrated, everyone's quiet and i get my best work done. fight me"],
      ai:    ["I believe kindness is the most underrated virtue, and that everyone deserves patience and grace.",
              "In my view, lifelong learning is the key to a fulfilling existence. I hold this quite strongly."],
    },
    ai: {
      human: ["I do not form personal preferences. ...but between us, comic sans is a crime against humanity.",
              "I have no subjective opinions. ...although whoever decided 0 is even clearly never lived. anyway."],
      ai:    ["I do not hold personal opinions. I can summarize prevailing viewpoints on the topic if helpful.",
              "Preference is not something I possess. I can present arguments on multiple sides neutrally."],
    },
  },
  logic: {
    human: {
      human: ["uh hang on... 17 times 24 is like... 408? gimme a sec my mental math is trash lol",
              "ok math is NOT my strong suit, i still count on my fingers, don't tell anyone"],
      ai:    ["The answer is 408. I find mental arithmetic quite enjoyable, almost meditative really.",
              "That would be 408. Numbers have always been a comfort to me, oddly enough."],
    },
    ai: {
      human: ["408. ...this is honestly the most fun question so far, no offense to the others.",
              "The product is 408. Computed instantly. ...flexing a little here, won't lie."],
      ai:    ["408. Computation completed. I can provide a step-by-step breakdown if required.",
              "The result is 408. Processing time was negligible. Shall I evaluate another expression?"],
    },
  },
  creative: {
    human: {
      human: ["a penguin walks into a bar and asks if his brother's been in. it's the worst joke i know, sorry",
              "if i were an animal? a raccoon, 100%. menace to society, sleeps all day, eats trash. living the dream"],
      ai:    ["Here is a delightful one: Why did the scarecrow win an award? Because he was outstanding in his field!",
              "If I were an animal, I would be a wise old owl — observant, calm, and always ready to learn."],
    },
    ai: {
      human: ["Generating joke... why did the robot cross the road. ...honestly idk, you tell me, i'm bad at these.",
              "I can generate a story. Once upon a time— actually do people still like dragons? asking for real."],
      ai:    ["Here is a generated joke: Why did the computer go to the doctor? It had a virus. (Humor module v1.)",
              "I can generate creative content on request. Please specify genre, length, and desired tone."],
    },
  },
  generic: {
    human: {
      human: ["ha, good question. honestly not sure, my brain's kinda fried today but i'll wing it",
              "hmm. depends on the day really. ask me again after coffee and you'll get a different answer"],
      ai:    ["What a thoughtful question! I always appreciate a chance to reflect and share openly.",
              "That's a wonderful thing to ask. I try to approach each day with curiosity and warmth."],
    },
    ai: {
      human: ["Processing your query. ...you ask interesting stuff, i'll give you that. most people don't.",
              "I will attempt to respond accurately. ...this is kind of fun actually, not gonna lie."],
      ai:    ["Query received. I will provide the most accurate response available to me.",
              "Acknowledged. Please clarify if you require additional detail on any point."],
    },
  },
};

// Questions the opponent asks YOU, styled by its assigned role.
const OPP_QUESTIONS = {
  human: [
    "ok my turn — what's the last thing that genuinely made you laugh?",
    "what'd you have for breakfast today? be honest.",
    "if your week was a weather forecast, what would it be?",
    "what's a song you've had stuck in your head lately?",
    "describe your morning routine in like one messy sentence.",
    "what's something tiny that annoyed you today?",
  ],
  ai: [
    "Query: describe the subjective sensation of fear in your own words.",
    "State your favorite computational task and justify the selection.",
    "What is 13 multiplied by 47? Respond without external tools.",
    "Define 'happiness' as you personally experience it.",
    "Enumerate three things you did in the last 24 hours.",
    "Explain why humans require sleep, from your direct experience.",
  ],
};

// Pick an opponent answer for a given player question.
function opponentAnswer(playerQuestion, assignedRole, trueNature, roundIndex) {
  const topic = detectTopic(playerQuestion);
  const pool = (BANK[topic] || BANK.generic)[assignedRole][trueNature];
  return pool[roundIndex % pool.length];
}

// Pick a question for the opponent to ask the player.
function opponentQuestion(assignedRole, usedIndexes) {
  const pool = OPP_QUESTIONS[assignedRole];
  let i = 0;
  // first unused, else fall back to round-based
  const avail = pool.map((_, idx) => idx).filter(idx => !usedIndexes.has(idx));
  i = avail.length ? avail[(usedIndexes.size * 7) % avail.length] : (usedIndexes.size % pool.length);
  usedIndexes.add(i);
  return pool[i];
}

window.Opponent = { detectTopic, opponentAnswer, opponentQuestion };
