# Week 9 — Quiz

> **Reproducibility note.** This quiz tests your knowledge of the mechanics of multiple sequence alignment, distance models, tree-building, and tree-file formats. Knowing the mechanics is the difference between a defensible tree and a tree-shaped picture; pin every parameter and the result is reproducible.

Ten multiple-choice questions on MAFFT, Clustal Omega, MUSCLE, Jukes-Cantor, Kimura two-parameter, UPGMA, neighbor-joining, maximum likelihood, bootstrap, and Newick / Nexus. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** The progressive alignment heuristic, as implemented by MAFFT, Clustal Omega, and MUSCLE, proceeds in three stages:

- A) Pick a reference sequence, BLAST every other sequence against it, take the BLAST output as the alignment.
- B) Compute pairwise distances; build a guide tree; align profiles along the guide tree from the leaves to the root.
- C) Apply the Needleman-Wunsch dynamic-programming algorithm in N dimensions for N sequences.
- D) Pick the longest sequence; align every other to it pairwise; concatenate.

---

**Q2.** The canonical MAFFT call `mafft --retree 2 --maxiterate 0 input.fasta` is preferred over `mafft --auto input.fasta` because:

- A) It is faster on every input.
- B) `--auto` is not reproducible across MAFFT versions; the algorithm choice can change between 7.490 and 7.526, producing a different alignment on the same input.
- C) `--retree 2` is more accurate than every other MAFFT mode.
- D) The exit code from `--auto` is unreliable.

---

**Q3.** The Jukes-Cantor distance correction `d_JC = -3/4 * ln(1 - 4/3 * p)` is needed because:

- A) The p-distance (proportion of differing columns) is biased downward at high `p` due to multiple substitutions per site that revert or chain.
- B) The p-distance is mathematically undefined.
- C) Biopython does not compute the p-distance directly.
- D) Newick files require corrected distances by spec.

---

**Q4.** The Kimura two-parameter (K2P) distance differs from Jukes-Cantor (JC) in that:

- A) K2P uses logarithms and JC does not.
- B) K2P allows transitions (A<->G, C<->T) and transversions (everything else) to have different rates; JC assumes all substitutions are equally likely.
- C) K2P is for proteins; JC is for nucleotides.
- D) K2P is faster to compute.

---

**Q5.** UPGMA produces a correct tree only when:

- A) The input alignment has more than 100 columns.
- B) The sequences evolved under a strict molecular clock (constant substitution rate across all lineages).
- C) The outgroup is the first sequence in the FASTA.
- D) The bootstrap is at least 90%.

---

**Q6.** Neighbor-joining, in contrast to UPGMA:

- A) Requires a substitution model.
- B) Produces an unrooted binary tree and does not assume a strict molecular clock; the tree is rooted post hoc on a named outgroup.
- C) Is slower than UPGMA.
- D) Cannot handle more than 50 sequences.

---

**Q7.** A bootstrap support value of 70% on an internal branch of a phylogenetic tree means:

- A) The branch is 70% likely to be in the true species tree.
- B) The branch appeared in 70% of the bootstrap-resampled trees; this is a measure of robustness to column resampling, not a measure of correctness in the species tree.
- C) The branch length is 70% of the alignment length.
- D) The branch is supported by 70 alignment columns.

---

**Q8.** The Felsenstein 1981 maximum-likelihood framework (free at the PMC link in resources.md) frames tree-building as:

- A) A Bayesian posterior inference under a uniform prior.
- B) The picking of the tree topology and branch lengths that maximize the probability of observing the alignment under a specified substitution model.
- C) The shortest-path problem in the alignment graph.
- D) An identical algorithm to neighbor-joining with the addition of a substitution matrix.

---

**Q9.** IQ-TREE 2's ultrafast bootstrap (UFBoot; Hoang et al. 2018) percentages are calibrated differently from the classical Felsenstein 1985 bootstrap, in the sense that:

- A) UFBoot is always higher than the classical bootstrap.
- B) A UFBoot 95% is roughly equivalent in interpretive weight to a classical bootstrap 70%; the two scales are not directly comparable without recalibration.
- C) UFBoot cannot produce values above 80%.
- D) UFBoot and classical bootstrap are identical.

---

**Q10.** When you ship a phylogenetic tree to a collaborator, the most important reproducibility metadata to include in the run-info JSON is:

- A) The Python version.
- B) The size of the alignment FASTA in bytes.
- C) The MAFFT version, the MAFFT algorithm flag, the trim threshold, the distance method, the tree-building tool and version, the substitution model (for ML), the random seed, the bootstrap replicate count, and the outgroup; plus the run date.
- D) The full Newick string copy-pasted into the email body.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — Pairwise distance, guide tree, profile-profile alignment. The progressive heuristic is the dominant family among MSA tools; all three of MAFFT, Clustal Omega, and MUSCLE are progressive aligners with different iterative-refinement strategies. Lecture 1 §2.

2. **B** — `--auto` is not reproducible across versions. MAFFT's algorithm-choice boundaries have shifted across 7.310, 7.490, and 7.526; pinning a specific algorithm (`--retree 2`, `--maxiterate 1000`, or `--localpair --maxiterate 1000`) is the only way to guarantee identical alignments across versions and machines. Lecture 1 §3, Resources style guide.

3. **A** — The p-distance is biased downward at high `p`. Multiple substitutions per site that revert or chain produce zero observed difference at a column even when several substitutions occurred. JC and K2P statistically recover the true expected distance under a substitution model. Lecture 2 §2, §3.

4. **B** — K2P separates transitions from transversions. In real DNA the transition / transversion ratio is usually 2-4 (in vertebrate mitochondrial DNA it is often 4-8). K2P captures this; JC does not. Lecture 2 §4.

5. **B** — UPGMA assumes a strict molecular clock. Under a strict clock every leaf is the same distance from the root and UPGMA recovers the true tree. Real biological sequences almost never satisfy the strict clock; UPGMA trees are systematically wrong when lineages evolve at different rates. Lecture 2 §6.

6. **B** — NJ is unrooted, does not assume a clock, and is fast `O(N^3)`. Rooting is post hoc, conventionally on a named outgroup. NJ does not need a substitution model; it works on the precomputed distance matrix. Lecture 2 §7.

7. **B** — Robustness to column resampling, not correctness in the species tree. A high bootstrap on a wrong branch is the standard failure mode of phylogenetics (long-branch attraction, alignment-induced artefacts). High support is necessary but not sufficient for correctness; report bootstrap percentages but do not over-interpret them. Lecture 2 §9.

8. **B** — Maximize the likelihood of the alignment under a substitution model, optimizing over topologies and branch lengths. The pruning algorithm computes per-column conditional probabilities recursively from the leaves to the root in time linear in the tree size. Heuristic search (NNI, SPR) explores the topology space. Lecture 2 §8, Felsenstein 1981 PMC link.

9. **B** — UFBoot 95% roughly corresponds to classical bootstrap 70%. The UFBoot calibration is by design (faster computation at the cost of a different scale). When reading other people's trees, check which bootstrap was used and interpret the numbers accordingly. Lecture 2 §9, Hoang et al. 2018.

10. **C** — All of the parameters above, plus the run date. Without them, the same FASTA can produce different trees on different days and the difference is impossible to debug. The run-info JSON is the most important artefact in any phylogenetics output directory; the trees themselves are derivative. Resources style guide; Lecture 3 §10.

</details>

---

If you scored under 7, re-read Lecture 1 §3 (MAFFT and reproducibility), Lecture 2 §3-§4 (JC and K2P), and Lecture 2 §9 (bootstrap interpretation). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
