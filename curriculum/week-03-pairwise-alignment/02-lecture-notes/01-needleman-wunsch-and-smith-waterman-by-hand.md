# Lecture 1 — Needleman-Wunsch and Smith-Waterman, by Hand

> **Duration:** ~3 hours of reading + paper-and-pencil + Python.
> **Outcome:** You can derive the Needleman-Wunsch recurrence from first principles, fill in the dynamic-programming matrix for a 5x5 worked example, traceback the optimal alignment, and do the same for Smith-Waterman with the zero-floor rule. You can state the time and space complexity in big-O without looking it up.

If you only remember one thing from this lecture, remember this:

> **Pairwise alignment is dynamic programming over a substitution score.** You have two sequences of length `m` and `n`. You build an `(m+1) x (n+1)` matrix of partial scores. Each cell is the max of three predecessors. The optimal alignment is a path through that matrix, recovered by traceback. Needleman-Wunsch (1970) does this globally; Smith-Waterman (1981) adds a single rule — floor at zero, start traceback from the global max — and gets local alignment for free.

This lecture is paper-and-pencil-heavy. Have a notebook open. Do not skip the by-hand example in §4. It is the hour that decides whether your NumPy implementation works on Thursday.

---

## 1. The problem

Given two sequences over an alphabet (DNA letters `{A, C, G, T}` or amino acids `{A, R, N, D, ...}`), we want the **best alignment**. An alignment is a pair of equal-length strings drawn from the original sequences plus a gap character `-`:

```
A:  A C G T A C G T
B:  A C - T A C G G
```

Each column is either a **match** (the two residues are identical), a **mismatch** (the two residues differ), or a **gap** (one side is `-`).

We score an alignment by summing per-column scores:

- A match scores `s(a, a)` — by convention positive.
- A mismatch scores `s(a, b)` — by convention negative, or in protein alignment a number from a substitution matrix.
- A gap costs `-g` per residue (the **linear** gap penalty) or `-(open + extend * (k - 1))` for a gap of length `k` (the **affine** gap penalty, covered in Lecture 2).

The substitution function `s(a, b)` is the **substitution matrix**. For DNA, the simplest choice is `s(a, a) = +1`, `s(a, b) = -1` for `a ≠ b`. For protein, you reach for **BLOSUM62** (Henikoff & Henikoff 1992) or **PAM250** (Dayhoff et al. 1978). We cover those in Lecture 2.

The problem is: find the alignment with **maximum total score**. The brute-force solution enumerates all alignments and is exponential in the sequence length. Needleman-Wunsch makes it `O(mn)`.

---

## 2. The Needleman-Wunsch recurrence (global alignment)

Saul Needleman and Christian Wunsch, *Journal of Molecular Biology* 48:443, 1970. Four pages, one matrix, the foundational paper of computational sequence alignment. You can read it in an hour. You should.

Define `H[i, j]` as the **best alignment score** between the prefix `A[1..i]` and the prefix `B[1..j]`. The recurrence is:

```
                ┌  H[i-1, j-1] + s(A[i], B[j])      (align A[i] with B[j])
H[i, j] = max   │  H[i-1, j  ] + g                  (gap in B: A[i] vs -)
                └  H[i  , j-1] + g                  (gap in A:  -  vs B[j])
```

With initial conditions:

```
H[0, 0] = 0
H[i, 0] = i * g        for i = 1 .. m   (a prefix of A aligned to all gaps)
H[0, j] = j * g        for j = 1 .. n   (a prefix of B aligned to all gaps)
```

Here `g` is the **gap penalty**, conventionally negative. If `g = -2`, then `H[3, 0] = -6` — aligning the first three letters of `A` to three gaps costs six points.

The final optimal score is `H[m, n]`. The optimal alignment itself is recovered by **traceback** from `H[m, n]` to `H[0, 0]`, choosing at each step the predecessor that achieved the max.

### Why this is correct

The proof is induction on `i + j`. Any alignment of `A[1..i]` and `B[1..j]` ends in exactly one of three ways:

1. `A[i]` aligned to `B[j]` — the rest is an alignment of `A[1..i-1]` and `B[1..j-1]`.
2. `A[i]` aligned to a gap — the rest is an alignment of `A[1..i-1]` and `B[1..j]`.
3. A gap aligned to `B[j]` — the rest is an alignment of `A[1..i]` and `B[1..j-1]`.

If `H[i-1, j-1]`, `H[i-1, j]`, and `H[i, j-1]` are all optimal scores for their respective sub-problems (induction hypothesis), then the max of the three extensions is the optimal score for `H[i, j]`. The base cases `H[i, 0]` and `H[0, j]` are forced — there is only one alignment of a prefix to the empty string, namely all gaps.

This is the textbook example of **optimal substructure**: the optimum of a sub-problem is part of the optimum of the larger problem. Every dynamic-programming algorithm in this course rests on the same idea.

---

## 3. The Smith-Waterman recurrence (local alignment)

Temple Smith and Michael Waterman, *Journal of Molecular Biology* 147:195, 1981. Four pages. The single most-cited paper in bioinformatics for a long stretch. Read it.

The local-alignment problem: instead of forcing the alignment to span both sequences end to end, find the **best matching sub-sequence** of `A` against the **best matching sub-sequence** of `B`. This is what you want when you are searching a 200-residue protein against a 30,000-protein database — most of the database is irrelevant, you want only the segment that lines up.

The recurrence is identical to Needleman-Wunsch except for **two changes**:

```
                ┌  H[i-1, j-1] + s(A[i], B[j])
                │  H[i-1, j  ] + g
H[i, j] = max   │  H[i  , j-1] + g
                └  0                                ← the new "zero floor"
```

And the initial conditions are simpler:

```
H[i, 0] = 0    for all i
H[0, j] = 0    for all j
```

The two changes:

1. **Zero floor.** Any cell whose three predecessors all yield negative values is reset to zero. The interpretation: "if continuing the alignment from this point would cost more than restarting, restart." A zero in the matrix is the start of a local alignment.
2. **Starting point of traceback.** Needleman-Wunsch tracebacks from `H[m, n]`. Smith-Waterman tracebacks from the **single largest value anywhere in the matrix**, and stops the first time it hits a zero. The aligned segment lies between those two coordinates.

Everything else — the scoring function, the gap penalty, the time complexity, the matrix layout — is identical. The change is two lines. The result is one of the most-used algorithms in biology.

### Why the zero floor works

Without the floor, a long region of mismatches between two distantly related sequences would push the score negative, and the optimal *global* score would penalize the entire alignment. The floor says: "give up on this attempt, set the running score back to zero, and let a better-scoring sub-region elsewhere in the matrix win." The optimal local alignment is the sub-path that accumulates the highest local score before being undone.

The 1981 paper proves this is correct in two pages. It is the kind of insight that, in retrospect, looks obvious. Most great algorithms do.

---

## 4. A worked example by hand

We are going to align two short DNA sequences with linear gap penalty `g = -2` and the simplest substitution scheme: `s(a, a) = +1` (match), `s(a, b) = -1` (mismatch).

Let:

- `A = G A T T A` (length 5)
- `B = G C A T G C U` — actually let's use `B = G C A T C` (length 5) for symmetry.

So we will build a 6x6 Needleman-Wunsch matrix indexed `i = 0..5` (rows, `A`) and `j = 0..5` (columns, `B`).

### 4.1 Initialization

The first row and column are multiples of `g = -2`:

```
            j=0   j=1   j=2   j=3   j=4   j=5
            -     G     C     A     T     C
   i=0  -    0    -2    -4    -6    -8   -10
   i=1  G   -2     ?     ?     ?     ?     ?
   i=2  A   -4     ?     ?     ?     ?     ?
   i=3  T   -6     ?     ?     ?     ?     ?
   i=4  T   -8     ?     ?     ?     ?     ?
   i=5  A  -10     ?     ?     ?     ?     ?
```

The `0` in the corner is the score of two empty prefixes. Each subsequent row/column entry adds `g = -2`.

### 4.2 Filling in `H[1, 1]`

`A[1] = G`, `B[1] = G`. Match: `s(G, G) = +1`. The three predecessors are:

- Diagonal: `H[0, 0] + s(G, G) = 0 + 1 = 1`.
- Up: `H[0, 1] + g = -2 + (-2) = -4`.
- Left: `H[1, 0] + g = -2 + (-2) = -4`.

Max is `1`. So `H[1, 1] = 1`. The "best predecessor" is the diagonal — store this arrow for the traceback.

### 4.3 Filling in `H[1, 2]`

`A[1] = G`, `B[2] = C`. Mismatch: `s(G, C) = -1`. Predecessors:

- Diagonal: `H[0, 1] + s(G, C) = -2 + (-1) = -3`.
- Up: `H[0, 2] + g = -4 + (-2) = -6`.
- Left: `H[1, 1] + g = 1 + (-2) = -1`.

Max is `-1`. So `H[1, 2] = -1`. Best predecessor: left (a gap in `A`).

### 4.4 Continuing row 1

Same procedure for `j = 3, 4, 5`. Each step picks the max of three predecessors. The completed row 1 is:

```
   i=1  G   -2     1    -1    -3    -5    -7
```

### 4.5 The complete matrix

Filling the rest by the same rule (do this on paper — the muscle memory is the whole point):

```
            j=0   j=1   j=2   j=3   j=4   j=5
            -     G     C     A     T     C
   i=0  -    0    -2    -4    -6    -8   -10
   i=1  G   -2     1    -1    -3    -5    -7
   i=2  A   -4    -1     0     0    -2    -4
   i=3  T   -6    -3    -2    -1     1    -1
   i=4  T   -8    -5    -4    -3     0     0
   i=5  A  -10    -7    -6    -3    -2    -1
```

The optimal global alignment score is `H[5, 5] = -1`.

### 4.6 Traceback

Start at `H[5, 5] = -1`. At each step, look at the three predecessors and pick the one whose value + the appropriate score equals the current cell. Record `match`, `gap-in-A`, or `gap-in-B`.

- `H[5, 5] = -1`. Diagonal predecessor `H[4, 4] = 0`, mismatch (`A vs C`) score `-1`. So `0 + (-1) = -1`. ✓ Match column with mismatch.
- `H[4, 4] = 0`. Diagonal `H[3, 3] = -1`, mismatch (`T vs T`) — wait, `A[4] = T`, `B[4] = T` is a match. `-1 + 1 = 0`. ✓ Match column (matches).
- `H[3, 3] = -1`. Diagonal `H[2, 2] = 0`, mismatch (`T vs A`) score `-1`. `0 + (-1) = -1`. ✓ Match column with mismatch.
- `H[2, 2] = 0`. Diagonal `H[1, 1] = 1`, mismatch (`A vs C`) score `-1`. `1 + (-1) = 0`. ✓ Match column with mismatch.
- `H[1, 1] = 1`. Diagonal `H[0, 0] = 0`, match (`G vs G`) score `+1`. `0 + 1 = 1`. ✓ Match column.

We have reached `H[0, 0]`. The recovered alignment, in reverse order:

```
A:  G  A  T  T  A
B:  G  C  A  T  C
    ↑  ↑  ↑  ↑  ↑
    =  X  X  =  X        (= = match, X = mismatch)
```

Two matches (`G` at position 1, `T` at position 4), three mismatches. Total score: `2 * 1 + 3 * -1 = -1`. ✓

### 4.7 The Smith-Waterman version of the same matrix

For local alignment, all entries are floored at zero. Re-fill (negative values become zero):

```
            j=0   j=1   j=2   j=3   j=4   j=5
            -     G     C     A     T     C
   i=0  -    0     0     0     0     0     0
   i=1  G    0     1     0     0     0     0
   i=2  A    0     0     0     1     0     0
   i=3  T    0     0     0     0     2     0
   i=4  T    0     0     0     0     1     1
   i=5  A    0     0     0     1     0     0
```

The global maximum is `H[3, 4] = 2`. Traceback from `(3, 4)`:

- `H[3, 4] = 2`. Diagonal `H[2, 3] = 1`, match (`T vs T`) score `+1`. `1 + 1 = 2`. ✓
- `H[2, 3] = 1`. Diagonal `H[1, 2] = 0`, match (`A vs A`) score `+1`. `0 + 1 = 1`. ✓
- `H[1, 2] = 0`. Stop — we have hit zero.

The local alignment is the two-column segment:

```
A:  A T
B:  A T
```

Score: 2. The remaining residues are left out of the local alignment, which is the whole point of Smith-Waterman.

---

## 5. Pseudocode

For Needleman-Wunsch (global) with linear gap penalty `g`:

```
function NW(A, B, s, g):
    m = len(A); n = len(B)
    H = matrix[(m+1) x (n+1)]
    trace = matrix[(m+1) x (n+1)]      # 0=stop, 1=diag, 2=up, 3=left

    H[0][0] = 0
    for i in 1..m: H[i][0] = i * g; trace[i][0] = 2  # up
    for j in 1..n: H[0][j] = j * g; trace[0][j] = 3  # left

    for i in 1..m:
        for j in 1..n:
            diag = H[i-1][j-1] + s(A[i], B[j])
            up   = H[i-1][j  ] + g
            left = H[i  ][j-1] + g
            best = max(diag, up, left)
            H[i][j] = best
            if   best == diag: trace[i][j] = 1
            elif best == up:   trace[i][j] = 2
            else:              trace[i][j] = 3

    # traceback from (m, n) to (0, 0)
    align_A = []; align_B = []
    i = m; j = n
    while i > 0 or j > 0:
        if trace[i][j] == 1:
            align_A.append(A[i]); align_B.append(B[j]); i -= 1; j -= 1
        elif trace[i][j] == 2:
            align_A.append(A[i]); align_B.append("-"); i -= 1
        else:
            align_A.append("-"); align_B.append(B[j]); j -= 1

    return H[m][n], reverse(align_A), reverse(align_B)
```

For Smith-Waterman (local), three changes:

```
function SW(A, B, s, g):
    ...
    H[i][0] = 0  for all i      # no gap-penalty initialization
    H[0][j] = 0  for all j

    for i in 1..m:
        for j in 1..n:
            ...
            best = max(diag, up, left, 0)   # the zero floor
            H[i][j] = best
            ...

    # traceback from the GLOBAL MAX of H, not from H[m][n]
    (i_start, j_start) = argmax(H)
    ...
    while H[i][j] != 0:    # stop at first zero, not at (0, 0)
        ...

    return H[i_start][j_start], reverse(align_A), reverse(align_B), (i_start, j_start)
```

Read both. Implement both. You will do this in NumPy on Thursday.

---

## 6. Time and space complexity

- **Time.** O(mn). The two nested loops each execute `m` and `n` times, and each cell does constant work (three predecessor lookups, three additions, a max).
- **Space.** O(mn) for the score matrix `H` and the traceback matrix. For 1000 x 1000 inputs, that is one million entries per matrix — a handful of MB in NumPy. For 30,000 x 30,000 (two bacterial genomes), it is one billion entries — 8 GB for int32. **Don't do that without Hirschberg.**

The Hirschberg trick reduces space to O(min(m, n)) by computing the score in two halves and recursively reconstructing the alignment. We skip it this week. Mention it in your head whenever you see a paper align two whole chromosomes — they used Hirschberg, or a banded variant, or a totally different algorithm (BLAST-style seed-and-extend, minimap2-style chaining).

For the standard bioinformatics use cases this week — pairs of proteins, short DNA fragments, single genes — O(mn) is fine. A 500 x 500 NumPy alignment runs in milliseconds.

---

## 7. Where the bugs hide

The traceback is where everyone first gets it wrong. Three classic bugs:

**1. Off-by-one between sequence indices and matrix indices.** The matrix has rows indexed `0..m` but the sequence is indexed `1..m`. When you look up `s(A[i], B[j])`, you must use `A[i-1]` and `B[j-1]` in zero-indexed Python. Get this wrong and every score is wrong by one cell.

**2. Tie-breaking in the max.** When two predecessors are equal, you must pick one consistently. Different choices produce different-but-equivalent alignments (same score, different gap placement). Biopython's `PairwiseAligner` documents its tie-break order; your code should pick one and stick with it. By convention: prefer diagonal over up over left (matches/mismatches over gaps).

**3. Traceback start condition for Smith-Waterman.** You traceback from the **single global max of `H`**, not from `H[m, n]`. If you copy your NW traceback into SW and forget to change the start cell, your local alignment is wrong on every input. Test this against EMBOSS `water` on a known pair.

**4. Forgetting the zero floor.** A Smith-Waterman cell with all three predecessors negative becomes `0`, not the negative max. The cell's traceback pointer is `stop` (i.e., this cell starts a fresh local alignment), not `diagonal`/`up`/`left`. Encode that explicitly.

When your unit test against `Bio.Align.PairwiseAligner` disagrees on the score, the bug is one of these four. Walk the 5x5 matrix by hand — that is faster than staring at NumPy output.

---

## 8. Recap and next lecture

You should now be able to:

- State the Needleman-Wunsch and Smith-Waterman recurrences from memory.
- Fill in a 5x5 NW or SW matrix by hand, with linear gap penalty, given `s(a, b)` and `g`.
- Walk the traceback from the appropriate starting cell to the appropriate stopping condition.
- Explain the two-line change that turns NW into SW (zero floor + traceback starts at the global max).
- State the time and space complexity of both algorithms.
- Diagnose the four classic traceback bugs.

In [Lecture 2](./02-substitution-matrices-and-gap-penalties.md) we replace the toy `+1 / -1` substitution scheme with the real ones — **PAM250**, **BLOSUM62**, **BLOSUM45**, **BLOSUM80** — and we replace the linear gap penalty with the affine gap penalty (Gotoh 1982) that every production aligner uses today. By the end of Lecture 2 you will be able to *choose* a matrix and a gap-penalty pair for a given pair of sequences, instead of using whatever your tool defaulted to.

---

*Continue to [Lecture 2 — Substitution Matrices and Gap Penalties](./02-substitution-matrices-and-gap-penalties.md).*
