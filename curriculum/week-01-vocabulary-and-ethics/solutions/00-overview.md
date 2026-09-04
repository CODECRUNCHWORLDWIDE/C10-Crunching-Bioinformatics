# Week 1 Solutions — Vocabulary and Ethics

These are the worked answers for Week 1. Read them **after** you have made an honest attempt at each item, not before.

That is not a formality. Week 1 is the week where you install the words. The moment you can no longer start a sentence — "a transcript is…" — and you have to sit there for thirty seconds with the cursor blinking, *that* is the learning event. Reading our definition first replaces the thirty seconds with a warm feeling of recognition and installs nothing. You will find out in Week 6, when `apply_hard_filters` is fog because "variant" is fog.

So: write your glossary entry badly first. Write the parser and let it crash. Then come here and diff your thinking against ours.

A note on what "correct" means in this unit. The code answers are exact — they run, and the expected output is printed verbatim so you can compare byte for byte. The written answers (glossary entries, the re-identification essay, the reflection) are *model* answers, not the only answers. Each one carries a "what a grader is looking for" note describing the range that earns full marks. If your wording differs from ours but hits the same substance, you are fine. If your wording differs because you defined a different thing, you are not.

## What is answered here

| Kind | Item | Answer |
|------|------|--------|
| Exercise | 1 — Glossary in Your Own Words | [exercises.md](./exercises.md#exercise-1--glossary-in-your-own-words) |
| Exercise | 2 — FASTA by hand | [exercises.md](./exercises.md#exercise-2--fasta-by-hand) |
| Exercise | 3 — Public-Data Inventory | [exercises.md](./exercises.md#exercise-3--public-data-inventory) |
| Challenge | 1 — Reverse complement, GC, and translate | [challenges.md](./challenges.md#challenge-1--reverse-complement-gc-and-translate) |
| Homework | 1 — Read a real methods section | [homework.md](./homework.md#problem-1--read-a-real-methods-section) |
| Homework | 2 — Hand-translate a sequence | [homework.md](./homework.md#problem-2--hand-translate-a-sequence) |
| Homework | 3 — Re-identification thought experiment | [homework.md](./homework.md#problem-3--re-identification-thought-experiment) |
| Homework | 4 — Build a tiny FASTA validator | [homework.md](./homework.md#problem-4--build-a-tiny-fasta-validator) |
| Homework | 5 — Refine your glossary | [homework.md](./homework.md#problem-5--refine-your-glossary) |
| Homework | 6 — Mini reflection essay | [homework.md](./homework.md#problem-6--mini-reflection-essay) |
| Mini-project | Personal Glossary and Public-Data Inventory | [mini-project.md](./mini-project.md) |
| Stretch goals | Week 1 README stretch list | [exercises.md](./exercises.md#stretch-goals) |

**Quiz.** The Week-1 quiz keeps its answer key inline, at the bottom of the quiz itself: [quiz.md — Answer key](../quiz.md#answer-key). It is not duplicated here. Take the quiz closed-book first; the key explains the reasoning for all ten questions.

**Mini-project reference implementation.** The Week-1 mini-project is a small public repository, not a single file, so the model answer is committed as a real directory tree at [`projects/solutions/week-01-vocabulary-and-ethics/`](../../../projects/solutions/week-01-vocabulary-and-ethics/). [mini-project.md](./mini-project.md) is the walkthrough: architecture, the decisions that were genuinely open, reading order, and how the result scores against the mini-project's own rubric.

## Two files this unit indexes but does not ship

Week 1's `exercises/README.md` and `challenges/README.md` link to two Python stubs that are **not currently in the repository**:

- `exercises/exercise-02-fasta-by-hand.py`
- `challenges/challenge-01-reverse-complement.py`

Homework Problems 2 and 4 both depend on them, so we have not skipped them. Instead, each answer states the interface it assumes, derived from the week's own prose:

- Exercise 2 is described in [the Week-1 README](../README.md) as "parse a FASTA file in pure Python — no Biopython, no regex tricks, just `open()` and a `for` loop," and Homework Problem 4 refers to "your `parse_fasta`."
- Challenge 1 is described in [the Week-1 README](../README.md) as "`reverse_complement`, `GC_content`, `translate` from scratch," and [challenges/README.md](../challenges/README.md) adds "no Biopython" and "test cases at the bottom of the file."

Those two sentences fully determine the deliverable, so the answers below are complete and runnable. If the stubs later land with different function signatures, the *logic* here still stands — rename to match.

## Environment these answers were verified against

| Thing | Version |
|-------|---------|
| Python | CPython 3.13.2 |
| Third-party packages | none — stdlib only, by design |
| Platform | Windows 11, run from a POSIX shell; the code is platform-independent |

Nothing in Week 1 needs Biopython, samtools, or a conda environment. That starts in Week 2. Every code answer here runs on any CPython 3.9 or newer (the only modern syntax used is `from __future__ import annotations` plus PEP 585/604 type hints, which that import makes safe on 3.9).

---

*Continue to [exercises.md](./exercises.md).*
