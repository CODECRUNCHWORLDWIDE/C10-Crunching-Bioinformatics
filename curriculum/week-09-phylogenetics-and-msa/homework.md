# Week 9 Homework

> **Reproducibility note.** Every homework problem produces a file the grader can read alongside a `run-info.json`. The numbers below are illustrative; your numbers will differ slightly with different MAFFT, Biopython, or IQ-TREE 2 versions. Pin the versions, pin the seed, and commit both the result file and the run-info.

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-09/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — MAFFT vs MUSCLE 5 on the same panel

**Problem statement.** Run the Exercise 1 pipeline twice on `data/cytb_vertebrates.fasta`: once with MAFFT (`--retree 2 --maxiterate 0`) and once with MUSCLE 5 (`-align`). Compare the two alignments column-by-column. For each pair of identical taxa, count the number of columns where MAFFT and MUSCLE 5 agree on the placement of every residue. Report the agreement fraction and identify the top three regions of disagreement.

Answer in `homework/notes/p1-mafft-vs-muscle.md`:

1. How many columns are in each alignment after trimming (drop columns with > 50% gaps)?
2. What is the column-by-column agreement fraction between the two trimmed alignments? (Use the longer alignment as the reference; align the shorter to it by NW or by column-index.)
3. List the top three disagreement regions: column range, sequences involved, brief description of what each aligner does there.

**Acceptance criteria.**

- `homework/p1_compare_aligners.py` runs end to end (calls MAFFT and MUSCLE 5 via subprocess; loads both with Biopython; emits a TSV of per-column agreement).
- `homework/results/p1_alignment_diff.tsv` exists.
- `homework/notes/p1-mafft-vs-muscle.md` contains the three answers with specific numbers.
- Commit message like `p1: MAFFT vs MUSCLE5 on 10-taxon cytb, 96% column agreement, 3 disagreements in indel-rich region`.

**Hint.** A simple agreement metric: for each column index `j` in the longer alignment, find the matching column in the shorter alignment via the residue mapping in a known taxon (e.g. human). Count columns where the residue tuples match across all ten taxa.

**Estimated time.** 60 minutes.

---

## Problem 2 — JC vs K2P sensitivity on the cytochrome b panel

**Problem statement.** Take the trimmed alignment from Exercise 1. Compute JC and K2P distance matrices. For each pair of taxa, report the JC distance, the K2P distance, and the (K2P - JC) drift. Plot the drift as a function of the JC distance and confirm the relationship: drift grows roughly linearly in JC up to ~0.4 and then accelerates.

Answer in `homework/notes/p2-jc-k2p.md`:

1. What is the maximum JC distance in the matrix? Between which two taxa?
2. What is the maximum (K2P - JC) drift in the matrix? Between which two taxa?
3. Does the (K2P - JC) drift correlate with the ratio of transitions to transversions for that taxon pair? Report the correlation coefficient.
4. Build NJ trees from the JC matrix and the K2P matrix. Do the two trees have the same topology (same bipartition set)?

**Acceptance criteria.**

- `homework/p2_distance_sensitivity.py` runs end to end.
- `homework/results/p2_drift_table.tsv` exists with columns: `taxon_a`, `taxon_b`, `p_dist`, `jc_dist`, `k2p_dist`, `k2p_minus_jc`, `ratio_ti_tv`.
- `homework/results/p2_drift_plot.png` exists.
- `homework/notes/p2-jc-k2p.md` contains four numbered answers.
- Commit message like `p2: JC vs K2P drift, max 0.024 platypus-chicken, NJ topologies agree`.

**Hint.** The transition / transversion ratio for vertebrate mitochondrial DNA is typically 4-8. The K2P correction adds more distance when the observed transition fraction is far from 1/3 of total differences.

**Estimated time.** 60 minutes.

---

## Problem 3 — Bootstrap on a larger panel

**Problem statement.** Take the 20-sequence 16S rRNA panel at `data/sixteen_s.fasta`. Run MAFFT, trim, build a K2P NJ tree, and run 500 bootstrap replicates with `--seed 42`. Report the support distribution and identify the lowest-support internal branch.

Answer in `homework/notes/p3-bootstrap-16s.md`:

1. How many internal nodes does the rooted tree have? (For a binary tree with N leaves, this is N-2 for the unrooted version and N-1 for the rooted.)
2. What is the median bootstrap support across all internal nodes? The minimum? The maximum?
3. Which internal node has the lowest support? Describe the clade in one sentence.
4. If you re-run with `--seed 99` instead of `--seed 42`, do any internal node supports shift by more than 5 percentage points?

**Acceptance criteria.**

- `homework/p3_bootstrap.py` runs end to end on the 20-sequence panel in under 60 seconds.
- `homework/results/p3_tree_with_support.nwk` exists with bootstrap percentages on internal nodes.
- `homework/results/p3_tree_with_support.png` exists.
- `homework/notes/p3-bootstrap-16s.md` contains four numbered answers.
- Commit message like `p3: 500-replicate bootstrap on 20-seq 16S panel, median support 86, min 41 on Bacteroidetes/Firmicutes split`.

**Hint.** With 500 replicates the bootstrap variance is small; two runs with different seeds typically differ by less than 5 percentage points on any single node. A larger shift indicates a weakly supported branch.

**Estimated time.** 75 minutes.

---

## Problem 4 — IQ-TREE 2 ML on the cytochrome b panel

**Problem statement.** Run IQ-TREE 2 on the Exercise 1 trimmed alignment with `-m MFP -B 1000 -alrt 1000 -T 4 -seed 42`. Inspect the `.iqtree` log to find: which substitution model ModelFinder selected by AIC, which by BIC, the final log-likelihood, the AIC and BIC values for the selected model. Compare the ML tree topology (the `.treefile`) to the NJ tree from Exercise 3.

Answer in `homework/notes/p4-iqtree-cytb.md`:

1. Which model did ModelFinder pick by AIC? By BIC? Are they the same?
2. What is the final log-likelihood and the AIC value?
3. How many internal branches differ between the NJ tree (K2P, Exercise 3) and the ML tree (best model, this problem)? List each disagreement.
4. For each disagreement, state which tree's claim is more defensible and why (typically the ML tree, because it uses column-level information and a model).

**Acceptance criteria.**

- `homework/p4_iqtree.sh` is a shell script that runs the IQ-TREE 2 call.
- `homework/results/p4_iqtree.iqtree` (the log) is checked in.
- `homework/results/p4_iqtree.treefile` is checked in.
- `homework/notes/p4-iqtree-cytb.md` contains four numbered answers.
- Commit message like `p4: IQ-TREE 2 ML on cytb, ModelFinder picked GTR+F+I+G4 by both AIC and BIC, NJ-vs-ML agrees on 7/8 internal bipartitions`.

**Hint.** IQ-TREE 2's `.iqtree` log has a section titled "ModelFinder" that lists the top-10 models by AIC and BIC. The "Best-fit model" line is the one to report. The log-likelihood is on the line "BEST SCORE FOUND" near the bottom.

**Estimated time.** 75 minutes.

---

## Problem 5 — Render a publication-quality figure with ete3

**Problem statement.** Take the IQ-TREE 2 ML tree from Problem 4. Render it with ete3, with: branches coloured by UFBoot tier (>=95 green, 70-94 blue, <70 red), leaf labels showing the common name (from `data/cytb_vertebrates.tsv`), the platypus root labeled "outgroup," and a scale bar in the bottom-left corner.

Answer in `homework/notes/p5-render.md`:

1. Paste a screenshot of the rendered tree into the notes file.
2. Describe one styling decision you made (e.g. font size, leaf-label format) and the rationale.
3. Identify the deepest-supported clade in the tree (the internal branch with the highest UFBoot) and the shallowest-supported clade.

**Acceptance criteria.**

- `homework/p5_render.py` runs end to end and produces a PNG.
- `homework/results/p5_ml_tree_coloured.png` exists.
- `homework/notes/p5-render.md` contains the answers.
- Commit message like `p5: ete3 render of cytb ML tree with UFBoot-coloured branches, deepest clade (human, chimp) UFBoot 100`.

**Hint.** `ete3.NodeStyle["hz_line_color"]` and `["hz_line_width"]` are the two attributes you want. `ete3.TreeStyle.scale = 800` sets the pixels-per-substitution-per-site. Set `ete3.TextFace(..., fsize=12)` for leaf label fonts. The `tree.render(out_path, w=1200, tree_style=ts)` call writes the PNG.

**Estimated time.** 60 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 400-500 word reflection at `homework/notes/week-09-reflection.md` answering:

1. Before Week 9, what did you think a phylogenetic tree was? After Week 9, what is it actually? Pick one stage of the pipeline (alignment, distance computation, tree building, bootstrap, rendering) and say what surprised you about how arbitrary the canonical defaults are.
2. The Felsenstein 1981 ML framework is the foundation of modern phylogenetics. After reading the paper (or the lecture's summary), what is the single most important conceptual move he made over the distance-based predecessors?
3. Bootstrap support is widely reported in published trees. After Week 9, what does a bootstrap percentage actually tell you, and what does it *not* tell you? Give two concrete scenarios where a high bootstrap would mislead you.
4. The mini-project produces a tree figure plus a run-info JSON. Imagine you hand the figure to a non-bioinformatician colleague who asks "what does this tree mean?". What is the *most important* sentence you would say first? Why?

**Acceptance criteria.**

- File exists, 400-500 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The boundaries you note here are what will keep you out of trouble in any future job that touches phylogenetics.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 60 min |
| 2 | 60 min |
| 3 | 75 min |
| 4 | 75 min |
| 5 | 60 min |
| 6 | 30 min |
| **Total** | **~6 h** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
