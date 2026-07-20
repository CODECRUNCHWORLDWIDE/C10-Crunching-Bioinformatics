# Week 3 — Challenges

One challenge this week. It is the algorithm-design problem that separates "I implemented Needleman-Wunsch with linear gaps" from "I can build a production-style aligner." Budget 90–120 minutes.

## Index

1. **[Challenge 1 — Affine gap penalty](challenge-01-affine-gap-penalty.md)** — extend your Needleman-Wunsch implementation to Gotoh's (1982) three-matrix recurrence so the gap penalty is `open + extend * (k - 1)` instead of `g * k`, without blowing up the time complexity past O(mn). (~110 min)

## How to work the challenge

- Read the prompt fully before writing any code. Sketch the three-matrix transitions on paper.
- **Verify against `Bio.Align.PairwiseAligner`.** A correct Gotoh implementation produces scores that match Biopython under matched gap-open / gap-extend parameters to the integer. If your score is off by 1 or 2, the bug is in your "open a gap" vs "extend a gap" transition logic.
- Pick a tie-break order and document it. Affine gap recurrences have more tie-cases than linear, and inconsistent tie-breaking is a frequent source of "my alignment looks different from Biopython but the score matches" issues.
