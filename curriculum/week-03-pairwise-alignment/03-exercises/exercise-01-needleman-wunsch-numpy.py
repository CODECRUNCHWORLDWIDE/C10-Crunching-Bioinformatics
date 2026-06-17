"""
Exercise 1 - Needleman-Wunsch in NumPy.

Goal: implement the global pairwise alignment algorithm from Lecture 1 in
pure NumPy. Use a `+1 / -1` DNA scoring scheme and a linear gap penalty
of `-2`. Verify against the 5x5 by-hand example from Lecture 1 §4.

Estimated time: 50 minutes.

Acceptance criteria:
- `python exercise-01-needleman-wunsch-numpy.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You implemented two functions: `build_dp_matrix` (fills the score and
  traceback matrices) and `traceback_alignment` (walks back from
  H[m, n] to H[0, 0] to produce the aligned strings).
- You used `numpy.zeros` to allocate the matrices and integer-typed
  arrays (`np.int32` is fine). No Python-level lists for the inner loop.

Required:
    python -m pip install numpy==1.26.4

What you learn:
- The Needleman-Wunsch recurrence translated to NumPy index arithmetic.
- How to encode the traceback as a small integer matrix
  (0 = stop, 1 = diag, 2 = up, 3 = left).
- The off-by-one between sequence indices (1-based in the recurrence)
  and Python list/array indices (0-based).
- Why tie-breaking matters and how to pick a convention you can defend.

TO COMPLETE: implement the two functions below. Run the file; all
assertions must pass.
"""

from __future__ import annotations

import numpy as np


# Toy DNA scoring scheme - the one used in Lecture 1 §4.
MATCH = 1
MISMATCH = -1
GAP = -2

# Traceback pointer codes. Use these constants throughout.
STOP = 0
DIAG = 1
UP = 2
LEFT = 3


def score(a: str, b: str) -> int:
    """Substitution score for a single residue pair under the toy scheme."""
    return MATCH if a == b else MISMATCH


def build_dp_matrix(seq_a: str, seq_b: str) -> tuple[np.ndarray, np.ndarray]:
    """Fill the Needleman-Wunsch dynamic-programming matrices.

    Returns:
        (H, T) where H is the (m+1, n+1) score matrix and T is the
        (m+1, n+1) traceback-pointer matrix with values in
        {STOP, DIAG, UP, LEFT}.
    """
    m, n = len(seq_a), len(seq_b)

    # Allocate. int32 is plenty for any sequence this week will throw at us.
    h = np.zeros((m + 1, n + 1), dtype=np.int32)
    t = np.zeros((m + 1, n + 1), dtype=np.int32)

    # TODO: initialize row 0 and column 0 with gap-penalty multiples.
    # H[0, 0] = 0 (already from np.zeros).
    # H[i, 0] = i * GAP, with T[i, 0] = UP.
    # H[0, j] = j * GAP, with T[0, j] = LEFT.
    #
    # Hint: this is a one-line vectorized fill with np.arange.
    raise NotImplementedError("Initialize row 0 and column 0 of h and t")

    # TODO: fill in the interior. For each i in 1..m and j in 1..n:
    #   diag = H[i-1, j-1] + score(seq_a[i-1], seq_b[j-1])
    #   up   = H[i-1, j  ] + GAP
    #   left = H[i  , j-1] + GAP
    #   best = max(diag, up, left)
    #   T[i, j] = whichever of DIAG / UP / LEFT was the argmax (break ties
    #             by preferring DIAG, then UP, then LEFT).
    #
    # A nested Python `for` loop is fine for this exercise - we are
    # building intuition, not chasing performance. Performance comes in
    # the mini-project, where we vectorize the anti-diagonal sweep.
    #
    # Hint: the cleanest tie-break is to compute the three values, then
    # pick DIAG if diag == best, elif up == best UP, else LEFT.

    # Implement the loop here.


def traceback_alignment(
    seq_a: str,
    seq_b: str,
    t: np.ndarray,
) -> tuple[str, str]:
    """Walk the traceback matrix from (m, n) to (0, 0).

    Returns:
        (aligned_a, aligned_b) - two equal-length strings over the
        alphabet of seq_a / seq_b plus the gap character '-'.
    """
    m, n = len(seq_a), len(seq_b)
    aligned_a: list[str] = []
    aligned_b: list[str] = []

    # TODO: walk from (m, n) back to (0, 0), appending residues/gaps to
    # aligned_a and aligned_b according to the traceback pointer.
    #
    # While i > 0 or j > 0:
    #   if t[i, j] == DIAG: emit (seq_a[i-1], seq_b[j-1]); i -= 1; j -= 1
    #   elif t[i, j] == UP:  emit (seq_a[i-1], '-');        i -= 1
    #   else (LEFT):         emit ('-', seq_b[j-1]);        j -= 1
    raise NotImplementedError("Walk the traceback")

    # The strings were built in reverse - reverse them on return.
    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def needleman_wunsch(seq_a: str, seq_b: str) -> tuple[int, str, str]:
    """Convenience wrapper - score plus aligned strings."""
    h, t = build_dp_matrix(seq_a, seq_b)
    aligned_a, aligned_b = traceback_alignment(seq_a, seq_b, t)
    return int(h[len(seq_a), len(seq_b)]), aligned_a, aligned_b


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-01-needleman-wunsch-numpy.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # The Lecture 1 §4 worked example.
    a = "GATTA"
    b = "GCATC"

    h, t = build_dp_matrix(a, b)

    # Shape checks.
    assert h.shape == (6, 6), f"H wrong shape: {h.shape}"
    assert t.shape == (6, 6), f"T wrong shape: {t.shape}"

    # Initialization (gap-penalty multiples).
    assert h[0, 0] == 0
    assert h[1, 0] == -2
    assert h[5, 0] == -10
    assert h[0, 5] == -10

    # A few interior cells from the by-hand matrix.
    assert h[1, 1] == 1, f"H[1,1] should be 1, got {h[1, 1]}"
    assert h[3, 4] == 1, f"H[3,4] should be 1, got {h[3, 4]}"
    assert h[5, 5] == -1, f"H[5,5] should be -1, got {h[5, 5]}"

    # Traceback yields the correct alignment.
    score_val, aligned_a, aligned_b = needleman_wunsch(a, b)
    assert score_val == -1, f"score should be -1, got {score_val}"
    assert len(aligned_a) == len(aligned_b)
    # The optimal alignment is one of a small set; under our tie-break
    # (prefer DIAG) we get the all-substitutions-no-gaps form:
    assert aligned_a == "GATTA", f"unexpected aligned_a: {aligned_a!r}"
    assert aligned_b == "GCATC", f"unexpected aligned_b: {aligned_b!r}"

    # A second example: two identical sequences should produce a perfect
    # all-match alignment whose score is len(seq) * MATCH.
    a2 = "ACGTACGT"
    s2, aa2, bb2 = needleman_wunsch(a2, a2)
    assert s2 == len(a2) * MATCH, f"identity case wrong score: {s2}"
    assert aa2 == a2 and bb2 == a2

    # A third example: alignment against the empty string should be all gaps.
    s3, aa3, bb3 = needleman_wunsch("ACGT", "")
    assert s3 == 4 * GAP, f"empty-vs-string score: {s3}"
    assert aa3 == "ACGT" and bb3 == "----"

    # Reverse case - empty first.
    s4, aa4, bb4 = needleman_wunsch("", "ACGT")
    assert s4 == 4 * GAP
    assert aa4 == "----" and bb4 == "ACGT"

    # Asymmetric example with an insertion: one of three valid optimal
    # alignments has a single gap in B.
    s5, aa5, bb5 = needleman_wunsch("ACGT", "AGT")
    # 3 matches + 1 gap = 3 + (-2) = 1
    assert s5 == 1, f"ACGT vs AGT score: {s5}"
    assert len(aa5) == len(bb5)
    assert aa5.replace("-", "") == "ACGT"
    assert bb5.replace("-", "") == "AGT"

    print("All assertions passed. Your NumPy Needleman-Wunsch agrees with")
    print("the Lecture 1 §4 by-hand example (score = -1) and three more")
    print("hand-verified cases. Move on to Exercise 2 - Smith-Waterman.")
