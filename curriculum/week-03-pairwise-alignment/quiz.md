# Week 3 — Quiz

Ten multiple-choice questions on Needleman-Wunsch, Smith-Waterman, substitution matrices, gap penalties, and the Biopython API. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** In a Needleman-Wunsch dynamic-programming matrix `H` for two sequences of length `m` and `n`, the optimal global-alignment score is read from which cell?

- A) `H[0, 0]`.
- B) `H[m, n]`.
- C) The global maximum of `H`.
- D) `H[m, 0] + H[0, n]`.

---

**Q2.** Smith-Waterman differs from Needleman-Wunsch in exactly two places. Which two?

- A) The substitution matrix is replaced with a position-specific scoring matrix, and gaps are disallowed.
- B) The recurrence takes the max with zero (the zero floor), and traceback starts at the global maximum of `H` instead of at `H[m, n]`.
- C) The matrix is filled column-by-column instead of row-by-row, and the gap penalty is replaced with an affine penalty.
- D) The score is divided by sequence length, and a probability cutoff replaces the substitution matrix.

---

**Q3.** Which of the following statements about substitution matrices is correct?

- A) PAM250 is appropriate for closely related proteins; BLOSUM62 is appropriate for distant ones.
- B) A high BLOSUM number (e.g., BLOSUM80) is built from clusters of sequences with high percent identity, and is appropriate for **closely related** proteins.
- C) PAM and BLOSUM matrices give identical scores for any input — they are mathematically equivalent.
- D) BLOSUM62 entries are derived from physical-chemistry principles of amino acid side chains, not from observed substitution frequencies.

---

**Q4.** BLOSUM62 scores `s('W', 'W')` (tryptophan vs tryptophan) at:

- A) +1.
- B) +4.
- C) +9.
- D) +11.

---

**Q5.** Aligning a 200-residue protein against a 300,000-protein database, you want to know which database entry has the best matching sub-region (not the best end-to-end alignment). The correct algorithm choice is:

- A) Needleman-Wunsch — global alignment is always preferred for database search.
- B) Smith-Waterman — local alignment, because you do not expect the full query to match the full database entry.
- C) Edit distance — it is symmetric and easier to compute.
- D) Hamming distance — for protein, single-character substitutions are the dominant signal.

---

**Q6.** A gap of length 5 under an affine gap penalty with `open = -11` and `extend = -1` costs:

- A) -5.
- B) -11.
- C) -15.
- D) -55.

---

**Q7.** The time complexity of the standard Needleman-Wunsch algorithm on two sequences of length `m` and `n` is:

- A) O(m + n).
- B) O(m log n).
- C) O(mn).
- D) O(2^(m+n)).

---

**Q8.** In Biopython 1.83, the recommended replacement for the deprecated `Bio.pairwise2` is:

- A) `Bio.PairwiseAlign`.
- B) `Bio.Align.PairwiseAligner`.
- C) `Bio.Align.LegacyAligner`.
- D) `Bio.AlignIO.pairwise`.

---

**Q9.** Which of the following is **the correct Biopython 1.83 way** to load BLOSUM62 and look up the score for `W vs W`?

- A) `Bio.Align.substitution_matrices.load("BLOSUM62")['W', 'W']`
- B) `Bio.SubsMat.MatrixInfo.blosum62[('W', 'W')]`
- C) `Bio.pairwise2.format_alignment(blosum62, 'W', 'W')`
- D) `Bio.Align.PairwiseAligner.score('W', 'W', matrix='BLOSUM62')`

---

**Q10.** Smith-Waterman traceback should stop when:

- A) `i == 0` and `j == 0`.
- B) The traceback reaches a cell whose value is zero.
- C) The traceback reaches `H[m, n]`.
- D) The traceback has emitted `min(m, n)` columns.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — `H[m, n]` is the score of the optimal alignment that uses *all* of `A` and *all* of `B` end-to-end. That is the definition of global alignment. The global max of `H` is where Smith-Waterman tracebacks start, but for NW it is `H[m, n]`.

2. **B** — The zero floor and the global-argmax traceback start are the two changes. Smith-Waterman uses the same substitution matrix and the same gap penalty as the NW it is built from; what changes is the recurrence (max-with-zero) and the traceback start cell.

3. **B** — High BLOSUM number = clustering at high percent identity = matrix appropriate for similar sequences. PAM is the opposite: high PAM number means distant. PAM250 is for distant, BLOSUM62 is the medium-divergence default. BLOSUM matrices are derived empirically from the BLOCKS database, not from chemistry.

4. **D** — Tryptophan is the largest and rarest of the 20 standard amino acids, and a `W vs W` match is highly informative. BLOSUM62 awards it +11, the largest single diagonal entry in the matrix.

5. **B** — Smith-Waterman, every time, for "find the matching sub-region in this database." BLAST is, under the hood, a heuristic-seeded Smith-Waterman. Global alignment is appropriate when the two sequences are known a priori to span the same region.

6. **C** — `open + extend * (k - 1) = -11 + (-1) * 4 = -15`. The first residue costs the open penalty; each additional residue costs the extend penalty.

7. **C** — Two nested loops over `m` and `n` cells, constant work per cell. O(mn) is the textbook answer. Hirschberg's trick reduces *space* to O(min(m, n)) but does not change time complexity.

8. **B** — `Bio.Align.PairwiseAligner` is the modern (Biopython 1.79+) API. `Bio.pairwise2` is deprecated and slated for removal; `Bio.SubsMat` is also deprecated in favor of `Bio.Align.substitution_matrices`. Get used to the new names.

9. **A** — `Bio.Align.substitution_matrices.load("BLOSUM62")['W', 'W']`. The returned object is a NumPy-array-like subclass with letter indexing. The other options reference deprecated APIs (B), nonexistent ones (C, D), or wrong call signatures.

10. **B** — Smith-Waterman traceback stops at the first zero, not at `H[0, 0]`. The traceback start is the global argmax of `H`, and the stop is the first zero on the back-path. The aligned region is what lies between those two coordinates.

</details>

---

If you scored under 7, re-read Lecture 1 for the algorithm questions and Lecture 2 for the substitution-matrix and gap-penalty questions. If you scored 9 or 10, you are ready to start the [homework](./homework.md).
