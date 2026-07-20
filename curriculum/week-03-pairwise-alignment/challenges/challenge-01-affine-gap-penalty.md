# Challenge 1 — Affine Gap Penalty (Gotoh 1982)

> **Estimated time:** 110 minutes.
> **Goal:** Implement Needleman-Wunsch with an affine gap penalty `open + extend * (k - 1)` using Gotoh's (1982) three-matrix recurrence in pure NumPy, then verify your scores agree with `Bio.Align.PairwiseAligner` (Biopython 1.83) to the integer across a battery of test pairs.

Lecture 2 §5 introduced affine gaps and explained why they match biology better than linear gaps. This challenge is the algorithmic substance: implementing them efficiently. A naive affine extension of Lecture 1's recurrence would blow time complexity up to O(mn(m + n)) — for each cell, you would have to look at all possible gap lengths back through the previous row/column. Gotoh's trick keeps the complexity at O(mn) by tracking three matrices simultaneously.

## Background — the three-matrix recurrence

Define three matrices, each of shape `(m+1) x (n+1)`:

- **`M[i, j]`** — the best alignment score for the prefixes `A[1..i]` and `B[1..j]` **assuming the alignment ends with a match/mismatch column** (i.e., `A[i]` aligned to `B[j]`).
- **`X[i, j]`** — the best score **assuming the alignment ends with a gap in `A`** (i.e., `B[j]` aligned to `-`). Confusingly, "gap in A" means `A` contributed `-` at the last column, so the last column is `(-, B[j])`. Different textbooks use different sign conventions for `X` and `Y`; pick one and stick with it.
- **`Y[i, j]`** — the best score **assuming the alignment ends with a gap in `B`** (i.e., `A[i]` aligned to `-`).

The recurrences (with `s` the substitution score, `o` the gap-open penalty, `e` the gap-extend penalty, both negative):

```
M[i, j] = max(  M[i-1, j-1]  ,
                X[i-1, j-1]  ,
                Y[i-1, j-1]  ) + s(A[i], B[j])

X[i, j] = max(  M[i,   j-1] + o  ,    # open a new gap in A
                X[i,   j-1] + e  )    # extend an existing gap in A

Y[i, j] = max(  M[i-1, j  ] + o  ,    # open a new gap in B
                Y[i-1, j  ] + e  )    # extend an existing gap in B
```

The interpretation: to start a new gap costs `o` (the gap-open penalty); to extend an existing gap costs `e` (the gap-extend penalty). The matrices encode "what was the state of the alignment immediately before this cell?" so that the open-vs-extend decision is made by looking at *which matrix* the predecessor came from.

The final score is `max(M[m, n], X[m, n], Y[m, n])`.

### Initialization

The initialization is where most implementations get subtle bugs. The conservative choice:

```
M[0, 0] = 0
X[0, 0] = Y[0, 0] = -infinity     # cannot end in a gap before any residue
M[i, 0] = -infinity   for i >= 1  # ending in match before B has started: impossible
M[0, j] = -infinity   for j >= 1
X[i, 0] = -infinity   for i >= 1  # cannot have a gap in A before any B residue
X[0, j] = o + e * (j - 1)  for j >= 1   # leading gap-in-A run
Y[i, 0] = o + e * (i - 1)  for i >= 1   # leading gap-in-B run
Y[0, j] = -infinity   for j >= 1
```

If you use a large negative integer (e.g. `-10**9`) for `-infinity`, your int32 NumPy matrices behave correctly — the impossible transitions never win the max.

### Traceback

Traceback is more involved than the linear case because the state at each step must remember *which of the three matrices* we are currently in:

- Start at the matrix with the maximum value in cell `(m, n)`.
- If we are in `M`, the previous step was a match/mismatch column: emit `(A[i], B[j])`, decrement both `i` and `j`, transition to whichever of `M`, `X`, `Y` achieved the max at `(i-1, j-1)`.
- If we are in `X`, the previous step was a gap-in-A column: emit `(-, B[j])`, decrement `j` only, transition to `M` (if the predecessor was the "open" branch) or stay in `X` (if the predecessor was the "extend" branch).
- If we are in `Y`, symmetric to `X`.

Continue until `i == 0` and `j == 0`. As with linear-gap NW, store traceback pointers as you fill the matrices to avoid re-deriving the transition decisions during traceback.

## Task

Write `gotoh.py` in your portfolio repo at `week-03/challenges/gotoh.py`. The script should:

1. Implement a function `affine_align(seq_a, seq_b, substitution_matrix, open_penalty, extend_penalty)` that fills the three matrices and returns the score plus aligned strings (with `-` characters for gap columns).
2. Accept either a NumPy substitution matrix + index dict (for DNA) or a `Bio.Align.substitution_matrices.Array` (for protein) — handle both gracefully.
3. Print, for each of the test pairs below, a one-line summary in the form:
   ```
   <pair-name>  score=<int>  biopython_score=<int>  match=<yes|no>
   ```
4. Exit with code 0 if all biopython comparisons match, 1 otherwise.

## Test pairs

Run your implementation on **all five** of these pairs and verify against Biopython.

### Pair 1 — DNA, no gaps possible

```python
a = "ACGT"
b = "ACGT"
```

Match score `+1`, mismatch `-1`, gap-open `-2`, gap-extend `-1`. Expected score: `4`.

### Pair 2 — DNA, one mid-sequence indel of length 1

```python
a = "ACGT"
b = "AGT"
```

Same scoring. Expected score: `1` (three matches plus a single-residue gap-open).

### Pair 3 — DNA, one mid-sequence indel of length 3

```python
a = "ACGTACGT"
b = "ACGT"
```

Same scoring. Expected score: with the affine `open=-2, extend=-1`, a single 4-residue gap costs `-2 + 3 * -1 = -5`. Plus 4 matches at `+1` each: total `-1`. Verify this on paper before trusting your code.

### Pair 4 — Protein, BLOSUM62 with BLAST-default affine gaps

```python
a = "WFQLRGFE"
b = "WFKLRFE"      # one residue deleted from the middle
```

BLOSUM62, gap-open `-11`, gap-extend `-1`. Expected: a single 1-residue gap, score on the order of `+30` (the published BLOSUM62 entries make this calculable by hand; do it).

### Pair 5 — Long-gap robustness

Build two 200 bp DNA sequences that differ by a single 50 bp insertion in the middle. Score under affine `open=-5, extend=-1`. Confirm that:

- Your affine score correctly assigns penalty `-5 + 49 * -1 = -54` to the gap (plus 150 matches at `+1` each: total `+96`).
- The same pair under a *linear* gap penalty of `-1` per residue would give `-50` for the gap (total `+100`) — confirming the affine penalty correctly discourages the long gap relative to a sparser scheme.

## Acceptance criteria

- `python gotoh.py` runs without crashing on the five test pairs.
- All five biopython comparisons report `match=yes`.
- Peak memory is bounded by O(mn) in matrix entries — for the longest test pair (200 x 250), that is ~50,000 cells per matrix, ~150,000 cells total, well under 1 MB at int32.
- Your code includes a docstring at the top citing Gotoh 1982 and the BLAST-default `(11, 1)` and EMBOSS-default `(10, 0.5)` affine settings.
- Commit the file plus a `notes/gotoh-notes.md` recording any tie-break decisions you made and any cases where your alignment differs from Biopython's at the byte level (same score, different gap placement is fine — record it explicitly).

## Hints (do not peek for at least 45 minutes)

<details>
<summary>Hint 1 — How do I avoid the -infinity overflow in int32?</summary>

Use `numpy.iinfo(np.int32).min // 2` as your stand-in for `-infinity`. Halving the minimum gives you headroom — adding `o + e * k` to it cannot underflow into a positive number even for the longest reasonable `k`. A common alternative is to use `np.full(shape, -10**9, dtype=np.int32)`; either is fine.
</details>

<details>
<summary>Hint 2 — What does my inner loop look like?</summary>

For each `(i, j)` in `1..m, 1..n`:

```python
m_diag = max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1])
M[i, j] = m_diag + s(A[i-1], B[j-1])

X[i, j] = max(M[i, j-1] + open_pen, X[i, j-1] + extend_pen)
Y[i, j] = max(M[i-1, j] + open_pen, Y[i-1, j] + extend_pen)
```

Three max-of-three (or max-of-two) calls per cell. The order matters for tie-breaking — by convention, prefer `M` over `X` over `Y` (matches over gap-in-A over gap-in-B).
</details>

<details>
<summary>Hint 3 — How do I record traceback pointers across three matrices?</summary>

Allocate three traceback matrices in parallel:
- `tM[i, j]` ∈ {0=from-M-diag, 1=from-X-diag, 2=from-Y-diag} — which matrix the diagonal predecessor came from.
- `tX[i, j]` ∈ {0=open-from-M, 1=extend-from-X} — which branch of X was taken.
- `tY[i, j]` ∈ {0=open-from-M, 1=extend-from-Y} — which branch of Y was taken.

Then the traceback walks state-by-state: "I am in matrix `M` at cell `(i, j)`; emit a column; jump to matrix `tM[i, j]` at cell `(i-1, j-1)`." The state machine is small.
</details>

<details>
<summary>Hint 4 — How do I verify against Biopython?</summary>

```python
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2     # Biopython's "open" is what you pay
aligner.extend_gap_score = -1   # for the first residue of a gap.
# WARNING: Biopython's convention is that a gap of length 1 costs
# `open_gap_score`, a gap of length 2 costs
# `open_gap_score + extend_gap_score`, etc. This matches Gotoh's
# convention as we described it above. Double-check by aligning two
# specific sequences whose gap-length you know.

print(aligner.score("ACGTACGT", "ACGT"))
```

If your `gotoh.py` reports a different score than Biopython, the bug is almost always in the initialization of `X[0, j]` and `Y[i, 0]`. Print those rows side-by-side with what Biopython's first row/column would imply.
</details>

## Stretch

If you finish under time and want more:

- Add `mode="local"` support — the local affine recurrence (Smith-Waterman + Gotoh) takes max-of-(M, X, Y, 0) and starts traceback at the global argmax of `max(M, X, Y)`. Verify against Biopython's `aligner.mode = "local"`.
- Vectorize the inner loop along the anti-diagonals. The diagonals of an `(m+1, n+1)` matrix can be filled in parallel — there are `m + n + 1` diagonals, each independent. NumPy's fancy-indexing lets you fill each diagonal with one vectorized expression. Expect ~5-10x speedup on a 500x500 alignment.
- Benchmark your vectorized implementation against Biopython's C extension on a 1000x1000 random DNA pair. Biopython will still win (it is C and SIMD-vectorized), but the gap should narrow from 100x to 10x or so. Document the timings in your `notes/`.

## What you should be able to do after this

- Implement Gotoh's three-matrix affine-gap recurrence without referring to a textbook.
- Read a real aligner's source code (BWA, BLAST, minimap2) and recognize the recurrence underneath the SIMD-vectorized inner loop.
- Choose `open` and `extend` deliberately for a given biological context, based on the expected indel-length distribution.

---

*Submit by committing `gotoh.py` and `notes/gotoh-notes.md` to your portfolio repo.*
