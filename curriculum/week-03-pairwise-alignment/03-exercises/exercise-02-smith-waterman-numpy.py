"""
Exercise 2 - Smith-Waterman in NumPy.

Goal: take your Exercise 1 Needleman-Wunsch implementation, change two
lines, and you have a local aligner. This is the smallest meaningful
change in introductory bioinformatics, and one of the most-cited
algorithms in biology (Smith & Waterman 1981).

Estimated time: 40 minutes.

Acceptance criteria:
- `python exercise-02-smith-waterman-numpy.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You implemented `build_sw_matrix` (Needleman-Wunsch with two changes:
  the zero-floor rule, and gap-free initialization of row 0/column 0)
  and `traceback_local` (traceback that starts at the global max of H
  and stops at the first zero, not at H[0, 0]).
- You returned the start and end coordinates of the local alignment in
  both sequences - this is what makes Smith-Waterman useful for
  "where in this protein does this domain hit?" queries.

Required:
    python -m pip install numpy==1.26.4

What you learn:
- The two-line change between NW and SW.
- How to find the global argmax of a NumPy array
  (`np.unravel_index(np.argmax(h), h.shape)`).
- Why the start/end coordinates matter operationally (BLAST returns
  them as `qstart`, `qend`, `sstart`, `send` - exactly these numbers).

TO COMPLETE: implement the two functions below. Run the file; all
assertions must pass.
"""

from __future__ import annotations

import numpy as np


MATCH = 1
MISMATCH = -1
GAP = -2

STOP = 0
DIAG = 1
UP = 2
LEFT = 3


def score(a: str, b: str) -> int:
    return MATCH if a == b else MISMATCH


def build_sw_matrix(seq_a: str, seq_b: str) -> tuple[np.ndarray, np.ndarray]:
    """Fill the Smith-Waterman dynamic-programming matrices.

    Differs from Needleman-Wunsch in exactly two places:
      1. Row 0 and column 0 are all zeros (no gap-penalty initialization).
      2. Each interior cell is the max of (diag, up, left, 0) - the
         "zero floor" rule. If the max is zero, the traceback pointer is
         STOP, marking the start of a fresh local alignment.

    Returns:
        (H, T) where H is the (m+1, n+1) score matrix and T is the
        (m+1, n+1) traceback-pointer matrix.
    """
    m, n = len(seq_a), len(seq_b)
    h = np.zeros((m + 1, n + 1), dtype=np.int32)
    t = np.zeros((m + 1, n + 1), dtype=np.int32)
    # Row 0 and column 0 are zero - exactly what np.zeros gave us. Done.

    # TODO: fill in the interior. For each i in 1..m and j in 1..n:
    #   diag = H[i-1, j-1] + score(seq_a[i-1], seq_b[j-1])
    #   up   = H[i-1, j  ] + GAP
    #   left = H[i  , j-1] + GAP
    #   best = max(diag, up, left, 0)        <-- the zero floor
    #   T[i, j] = DIAG / UP / LEFT / STOP    <-- STOP if best == 0
    #
    # Tie-break order: prefer DIAG > UP > LEFT > STOP.
    raise NotImplementedError("Fill the interior with the SW recurrence")


def traceback_local(
    seq_a: str,
    seq_b: str,
    h: np.ndarray,
    t: np.ndarray,
) -> tuple[str, str, tuple[int, int], tuple[int, int]]:
    """Traceback the local alignment.

    Start at the GLOBAL MAX of h (not at h[m, n]), stop at the first
    zero (not at h[0, 0]).

    Returns:
        (aligned_a, aligned_b, (a_start, a_end), (b_start, b_end))
        where a_start, a_end are 1-based, inclusive coordinates in
        seq_a (matching BLAST's `qstart`/`qend` convention) and similarly
        for b_start, b_end. If the matrix is all zeros, return four
        empty / zero values.
    """
    if h.max() == 0:
        return "", "", (0, 0), (0, 0)

    # TODO: locate the global argmax of h.
    #   Use np.unravel_index(np.argmax(h), h.shape) to get (i, j).
    #   Record (i, j) as the end coordinate of the local alignment.
    #
    # Then walk back, exactly like NW traceback, but stop when t[i, j] == STOP
    # (equivalently: when h[i, j] == 0). Record the (i, j) at which you
    # stopped as the start coordinate.
    #
    # Convert from matrix coordinates (1..m, 1..n) to sequence coordinates
    # (1..m for seq_a, 1..n for seq_b) - they happen to coincide because
    # the matrix is indexed with a leading zero row/column.

    aligned_a: list[str] = []
    aligned_b: list[str] = []

    raise NotImplementedError("Walk the local traceback")

    # When you exit the loop, return:
    # (
    #     "".join(reversed(aligned_a)),
    #     "".join(reversed(aligned_b)),
    #     (start_i, end_i),   # 1-based, inclusive
    #     (start_j, end_j),
    # )


def smith_waterman(
    seq_a: str,
    seq_b: str,
) -> tuple[int, str, str, tuple[int, int], tuple[int, int]]:
    """Convenience wrapper - score, aligned strings, coordinates."""
    h, t = build_sw_matrix(seq_a, seq_b)
    aligned_a, aligned_b, a_coords, b_coords = traceback_local(
        seq_a, seq_b, h, t
    )
    return int(h.max()), aligned_a, aligned_b, a_coords, b_coords


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-02-smith-waterman-numpy.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # The Lecture 1 §4 worked example - same input as Exercise 1.
    a = "GATTA"
    b = "GCATC"

    h, t = build_sw_matrix(a, b)

    # Shape checks.
    assert h.shape == (6, 6)
    assert t.shape == (6, 6)

    # Row 0 / column 0 are zero.
    assert h[0, :].sum() == 0
    assert h[:, 0].sum() == 0

    # A few interior cells from the by-hand Smith-Waterman matrix in
    # Lecture 1 §4.7.
    assert h[1, 1] == 1, f"H[1,1] should be 1, got {h[1, 1]}"
    assert h[2, 3] == 1, f"H[2,3] should be 1, got {h[2, 3]}"
    assert h[3, 4] == 2, f"H[3,4] should be 2, got {h[3, 4]}"

    # The global max is at (3, 4) with value 2.
    assert int(h.max()) == 2

    # Traceback the local alignment.
    score_val, aa, bb, a_coords, b_coords = smith_waterman(a, b)
    assert score_val == 2, f"local score should be 2, got {score_val}"
    assert aa == "AT" and bb == "AT", f"local alignment: {aa!r} vs {bb!r}"
    # The local alignment covers positions 2-3 of A (1-based, inclusive)
    # and positions 3-4 of B.
    assert a_coords == (2, 3), f"a_coords: {a_coords}"
    assert b_coords == (3, 4), f"b_coords: {b_coords}"

    # Identity case: aligning a sequence with itself gives a local
    # alignment that spans the whole thing with score = len * MATCH.
    a2 = "ACGTACGT"
    s2, aa2, bb2, ca2, cb2 = smith_waterman(a2, a2)
    assert s2 == len(a2) * MATCH
    assert aa2 == a2
    assert bb2 == a2
    assert ca2 == (1, len(a2))
    assert cb2 == (1, len(a2))

    # No-match case: two sequences with no shared characters should have
    # a local score of 0 (the zero floor makes the optimal local
    # alignment "the empty alignment").
    s3, aa3, bb3, ca3, cb3 = smith_waterman("AAAA", "CCCC")
    assert s3 == 0
    # An empty traceback is fine here - no positive-scoring sub-region.

    # A more interesting case: a short matching region embedded in
    # surrounding mismatches.
    a4 = "TTTACGTAAA"
    b4 = "GGGACGTCCC"
    # The shared core "ACGT" scores +4 (4 matches at +1 each); the
    # flanking regions are all mismatches/gaps and the zero floor will
    # exclude them.
    s4, aa4, bb4, ca4, cb4 = smith_waterman(a4, b4)
    assert s4 == 4, f"shared core score should be 4, got {s4}"
    assert aa4 == "ACGT", f"local aligned_a: {aa4!r}"
    assert bb4 == "ACGT", f"local aligned_b: {bb4!r}"
    assert ca4 == (4, 7), f"a coords: {ca4}"
    assert cb4 == (4, 7), f"b coords: {cb4}"

    print("All assertions passed. Your NumPy Smith-Waterman agrees with")
    print("the Lecture 1 §4.7 by-hand example (local score = 2) and four")
    print("more hand-verified cases. Move on to Exercise 3 - Biopython.")
