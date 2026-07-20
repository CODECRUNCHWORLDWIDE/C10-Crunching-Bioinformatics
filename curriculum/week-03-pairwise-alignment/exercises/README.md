# Week 3 — Exercises

Three focused Python drills. Each one runs from the command line, finishes inside 60 minutes, and produces output you can keep. Together they take you from "I filled in the matrix on paper" to "I have a working NumPy aligner that agrees with Biopython on score."

## Index

1. **[Exercise 1 — Needleman-Wunsch in NumPy](exercise-01-needleman-wunsch-numpy.py)** — implement the global-alignment recurrence from scratch with a `+1 / -1` DNA scoring scheme and linear gap penalty. (~50 min)
2. **[Exercise 2 — Smith-Waterman in NumPy](exercise-02-smith-waterman-numpy.py)** — re-use most of Exercise 1, change the two lines that turn NW into SW. (~40 min)
3. **[Exercise 3 — Compare against Biopython](exercise-03-biopython-compare.py)** — verify your scores against `Bio.Align.PairwiseAligner` on a battery of DNA and protein inputs. (~45 min)

## How to work the exercises

- Install Biopython 1.83 and NumPy 1.26 before starting: `python -m pip install biopython==1.83 numpy==1.26.4`.
- **Type the code yourself.** Copy-paste defeats the muscle-memory point.
- Each file is runnable: `python exercise-XX.py`. All assertions must pass.
- If you get stuck for more than 10 minutes, peek at the inline hints, then walk the 5x5 example from Lecture 1 §4 by hand.
- No solutions are checked in. The exercises are self-checking — every script ends in an assertion block that confirms the outputs match the by-hand example.
