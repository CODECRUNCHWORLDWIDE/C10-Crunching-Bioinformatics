"""
Exercise 3 - Compare against Biopython.

Goal: verify the scores produced by `Bio.Align.PairwiseAligner` (the
modern Biopython 1.83 API that replaced the deprecated `Bio.pairwise2`)
match scores produced by code that uses the same scoring scheme.

We do not import your NumPy implementations from Exercises 1 and 2 -
this script stands alone. The point is to internalize Biopython's
API surface, which the mini-project's benchmark will lean on heavily.

Estimated time: 45 minutes.

Acceptance criteria:
- `python exercise-03-biopython-compare.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You configured a `PairwiseAligner` for both global and local mode.
- You loaded BLOSUM62 from `Bio.Align.substitution_matrices` and
  confirmed `s('W', 'W') == 11` and `s('L', 'I') == 2`.
- You scored a real protein pair (a short SARS-CoV-2 spike fragment
  against its closest SARS-CoV-1 homolog) under BLOSUM62 with affine
  gap-open -11 / gap-extend -1, and printed both the score and the
  alignment.

Required:
    python -m pip install biopython==1.83

What you learn:
- The `Bio.Align.PairwiseAligner` API: mode, match_score / mismatch_score,
  open_gap_score / extend_gap_score, substitution_matrix.
- The difference between `aligner.score(a, b)` (just the score, fast)
  and `aligner.align(a, b)` (iterator of `Alignment` objects, slower).
- Why BLOSUM62 entries are integers in the canonical NCBI file but
  Biopython surfaces them as floats.
- How to read an `Alignment` object's `aligned` attribute (a pair of
  NumPy-style coordinate arrays).

TO COMPLETE: implement the three functions below. Run the file; all
assertions must pass.
"""

from __future__ import annotations

import numpy as np

from Bio.Align import PairwiseAligner, substitution_matrices


def configure_dna_global() -> PairwiseAligner:
    """Build a PairwiseAligner for DNA global alignment with the toy
    +1 / -1 / gap -2 scheme used in Exercises 1 and 2.

    The aligner should:
      - have mode 'global'
      - score matches at +1 and mismatches at -1
      - use a LINEAR gap penalty of -2 per residue
        (set BOTH open_gap_score and extend_gap_score to -2 -
        this is how Biopython represents the linear case)

    Returns: a configured PairwiseAligner.
    """
    # TODO:
    # aligner = PairwiseAligner()
    # aligner.mode = "global"
    # aligner.match_score = 1
    # aligner.mismatch_score = -1
    # aligner.open_gap_score = -2
    # aligner.extend_gap_score = -2
    # return aligner
    raise NotImplementedError("Configure a DNA global aligner")


def configure_dna_local() -> PairwiseAligner:
    """Same as configure_dna_global but mode='local' for Smith-Waterman."""
    # TODO: same as above with mode = 'local'.
    raise NotImplementedError("Configure a DNA local aligner")


def configure_protein_blosum62() -> PairwiseAligner:
    """Build a PairwiseAligner for protein alignment with BLOSUM62 and
    affine gaps (open -11, extend -1) - the BLAST defaults.

    The aligner should:
      - have mode 'global' (we'll switch to local in the assertions
        when we need to)
      - use the BLOSUM62 substitution matrix from
        Bio.Align.substitution_matrices
      - set open_gap_score = -11
      - set extend_gap_score = -1
      - set end_gap_score = 0 (no penalty for unaligned residues at
        the ends - we want a 'semi-global' result for short queries)
    """
    # TODO:
    # aligner = PairwiseAligner()
    # aligner.mode = "global"
    # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # aligner.open_gap_score = -11
    # aligner.extend_gap_score = -1
    # aligner.target_end_gap_score = 0
    # aligner.query_end_gap_score = 0
    # return aligner
    raise NotImplementedError("Configure a BLOSUM62 protein aligner")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-biopython-compare.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # ------------------------------------------------------------------
    # 1) DNA global alignment - reproduce the Lecture 1 §4 by-hand score.
    # ------------------------------------------------------------------
    global_aligner = configure_dna_global()
    assert global_aligner.mode == "global"

    score_gatta_gcatc = global_aligner.score("GATTA", "GCATC")
    # By-hand answer: -1.
    assert score_gatta_gcatc == -1, (
        f"GATTA vs GCATC global score should be -1, got {score_gatta_gcatc}"
    )

    # Identity case.
    assert global_aligner.score("ACGTACGT", "ACGTACGT") == 8

    # ------------------------------------------------------------------
    # 2) DNA local alignment - reproduce the Lecture 1 §4.7 SW score.
    # ------------------------------------------------------------------
    local_aligner = configure_dna_local()
    assert local_aligner.mode == "local"

    score_local = local_aligner.score("GATTA", "GCATC")
    # By-hand answer for the local case: 2 (the "AT" shared subsequence).
    assert score_local == 2, f"local score should be 2, got {score_local}"

    # The embedded-core case from Exercise 2.
    embedded_local = local_aligner.score("TTTACGTAAA", "GGGACGTCCC")
    assert embedded_local == 4, (
        f"embedded ACGT local score should be 4, got {embedded_local}"
    )

    # ------------------------------------------------------------------
    # 3) Load BLOSUM62 and confirm canonical entries.
    # ------------------------------------------------------------------
    blosum62 = substitution_matrices.load("BLOSUM62")
    # Tryptophan vs tryptophan: the famous +11 entry.
    assert int(blosum62["W", "W"]) == 11, (
        f"BLOSUM62 W vs W should be 11, got {blosum62['W', 'W']}"
    )
    # Leucine vs isoleucine: the conservative-aliphatic substitution.
    assert int(blosum62["L", "I"]) == 2
    # Leucine vs glutamate: the harsh hydrophobic-vs-charged substitution.
    assert int(blosum62["L", "E"]) == -3
    # Alanine vs alanine: the unremarkable match.
    assert int(blosum62["A", "A"]) == 4
    # Cysteine vs cysteine: a structurally informative match.
    assert int(blosum62["C", "C"]) == 9

    # ------------------------------------------------------------------
    # 4) Protein alignment - score a small SARS-CoV-2 / SARS-CoV-1 spike
    #    fragment pair under BLOSUM62 with affine gaps.
    # ------------------------------------------------------------------
    # First 60 residues of the SARS-CoV-2 spike protein (UniProt P0DTC2).
    sars2_spike = (
        "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS"
    )
    # First 60 residues of the SARS-CoV-1 spike protein (UniProt P59594).
    # ~76% identical at this region - a "moderate divergence" case where
    # BLOSUM62 is the canonical choice.
    sars1_spike = (
        "MFIFLLFLTLTSGSDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFL"
    )

    protein_aligner = configure_protein_blosum62()
    assert protein_aligner.mode == "global"
    # The substitution_matrix attribute should be set.
    assert protein_aligner.substitution_matrix is not None

    # Both inputs are 60 residues; the alignment fits BLOSUM62 with
    # affine gaps; the score is a positive integer (Biopython returns
    # float, but the value is integer-valued for integer matrices).
    s_pair = protein_aligner.score(sars2_spike, sars1_spike)
    assert s_pair > 0, f"spike pair score should be positive, got {s_pair}"

    # Sanity check: the score is bounded by the perfect-match case for
    # the shorter sequence (60 residues * max single-residue score of
    # 11 = 660). We expect well below that.
    upper_bound = len(sars1_spike) * 11
    assert s_pair < upper_bound

    # Score a self-alignment for cross-checking. The score should be
    # sum_i BLOSUM62(a_i, a_i) for sequence a.
    self_score = protein_aligner.score(sars2_spike, sars2_spike)
    expected_self = sum(int(blosum62[r, r]) for r in sars2_spike)
    assert int(self_score) == expected_self, (
        f"self-alignment {self_score} != sum of diagonal "
        f"BLOSUM62 entries {expected_self}"
    )

    # Get the alignment itself (not just the score). Biopython's
    # PairwiseAligner.align returns an iterator of Alignment objects.
    alignments = protein_aligner.align(sars2_spike, sars1_spike)
    # Take the first one. Note: there may be multiple equally-scoring
    # alignments; Biopython hands them out lazily.
    top = next(iter(alignments))
    assert top.score == s_pair

    # The `aligned` attribute is a 2-tuple of arrays giving the
    # coordinate ranges of the aligned blocks in each sequence. For a
    # gap-free alignment this is a single block; for an alignment with
    # gaps it is multiple blocks.
    aligned_blocks = top.aligned
    assert len(aligned_blocks) == 2  # two sequences

    # Convert the alignment to a string we can print and inspect.
    # Biopython's str(Alignment) produces the three-line "track" format
    # familiar from EMBOSS needle: top sequence, match indicator, bottom.
    alignment_text = str(top)
    assert "M" in alignment_text  # both sequences start with M
    # The string format is multi-line; we print it as a sanity check.

    print("All assertions passed.")
    print()
    print("BLOSUM62 protein alignment summary:")
    print(f"  Sequences:        60-residue spike N-terminal fragments")
    print(f"  SARS-CoV-2 (P0DTC2): {sars2_spike[:40]}...")
    print(f"  SARS-CoV-1 (P59594): {sars1_spike[:40]}...")
    print(f"  Substitution:     BLOSUM62")
    print(f"  Gap-open:         -11")
    print(f"  Gap-extend:       -1")
    print(f"  Alignment score:  {int(s_pair)}")
    print()
    print("Top alignment (first 200 chars of Biopython's track format):")
    print(alignment_text[:200])
    print()
    print("You are ready for the mini-project benchmark.")

    # Touch numpy so we know it imported - the mini-project will need it.
    _ = np.zeros((2, 2), dtype=np.int32)
