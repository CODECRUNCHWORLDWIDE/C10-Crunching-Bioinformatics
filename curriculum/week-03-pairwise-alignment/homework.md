# Week 3 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-03/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Fill in a Needleman-Wunsch matrix by hand and confirm in Python

**Problem statement.** Pick two short DNA sequences of length 6 and 7 (your choice — make them differ by at least three residues so the alignment is non-trivial). Use `match = +1`, `mismatch = -1`, `gap = -2`. By hand, on paper, fill in the entire 7x8 NW matrix. Mark the traceback path with arrows.

Then, in `homework/p1_nw_byhand.py`, hard-code the same two sequences and use your Exercise 1 implementation to compute the matrix. Print it. Confirm:

1. Your by-hand matrix equals your Python matrix entry-by-entry.
2. Your by-hand traceback path matches one of the equally-optimal traceback paths your code returns.

**Acceptance criteria.**

- File runs as `python homework/p1_nw_byhand.py`.
- Prints the matrix and the alignment.
- Commit includes a scanned or photographed image of your by-hand work as `homework/p1_byhand.jpg` (or `.png`). Photo with a phone is fine.
- Commit message like `p1: NW by hand matches Python on AGCTAG vs AGCTTAG`.

**Hint.** If your by-hand answer disagrees with your Python answer, the bug is almost always (1) an off-by-one in the gap-penalty initialization of row 0 / column 0, or (2) a tie-break inconsistency. Walk three cells side-by-side with paper before reaching for a debugger.

**Estimated time.** 60 minutes.

---

## Problem 2 — Smith-Waterman on a real ORF pair

**Problem statement.** Download two small protein sequences from UniProt:

- Human cytochrome c (`P99999`) — 105 residues.
- Yeast (*Saccharomyces cerevisiae*) cytochrome c iso-1 (`P00044`) — 109 residues.

These are textbook orthologs separated by ~1.5 billion years of evolution. Write `homework/p2_cyt_c.py` that:

1. Fetches both sequences from UniProt (you may either pre-download the FASTA files into `homework/data/` and parse with `Bio.SeqIO`, or use `Bio.ExPASy.get_sprot_raw` to fetch live).
2. Configures a `Bio.Align.PairwiseAligner` with mode `local`, substitution matrix BLOSUM62, gap-open `-11`, gap-extend `-1`.
3. Computes the local alignment, prints the score, prints the alignment in EMBOSS-style three-line format.
4. Reports percent identity over the aligned region.

**Acceptance criteria.**

- File runs as `python homework/p2_cyt_c.py`.
- Prints score, alignment, and percent identity.
- Percent identity falls in the 50-70% range — cytochrome c is highly conserved across eukaryotes but not identical. If your result is far outside this range, your aligner is misconfigured.
- Committed.

**Hint.** Biopython's `Alignment` object has a `.aligned` attribute (a pair of coordinate arrays) and a `__str__` that produces the three-line format. Percent identity is `count(matches) / count(aligned-positions)` where you exclude gap columns from the denominator.

**Estimated time.** 60 minutes.

---

## Problem 3 — Choose the right matrix

**Problem statement.** For each of the following protein pairs, identify the appropriate BLOSUM matrix (BLOSUM45, BLOSUM62, BLOSUM80) and justify the choice in one sentence:

1. Human hemoglobin alpha (HBA1) vs human hemoglobin beta (HBB) — two paralogs from a shared duplication event in early vertebrates.
2. Human cytochrome c (P99999) vs human cytochrome c (P99999) — the same protein, sanity check.
3. Human cytochrome c (P99999) vs *Methanocaldococcus jannaschii* cytochrome b6 (Q60366) — eukaryote vs archaeon, ~3 billion years of divergence.
4. SARS-CoV-2 spike (P0DTC2) vs SARS-CoV-1 spike (P59594) — ~80% identical at the protein level.
5. *E. coli* DNA polymerase III alpha subunit vs human DNA polymerase delta catalytic subunit — distantly homologous replicative polymerases.

Write your answers in `homework/notes/p3-matrix-choice.md`. Then, in `homework/p3_verify.py`, run a BLOSUM62-default local alignment on pair (4) and confirm that the score is in the high three- or four-digit range (BLAST would report this as an "obvious homolog" hit with E-value near 0).

**Acceptance criteria.**

- `notes/p3-matrix-choice.md` has five numbered answers with one-sentence justifications each.
- `p3_verify.py` runs, prints the score and a one-line "Decision: BLOSUM62 / BLOSUM45 / BLOSUM80" tag.
- Committed.

**Hint.** Use the rule of thumb in Lecture 2 §4: BLOSUM80 for close (>50% identity), BLOSUM62 for moderate (25-50%), BLOSUM45 for distant (<25%). For pair (3), neither cytochrome c nor cytochrome b6 are universally conserved enough to expect anything alignable at this distance — your answer should note that.

**Estimated time.** 45 minutes.

---

## Problem 4 — Affine vs linear gap penalty, side by side

**Problem statement.** Build two DNA sequences:

- `a = "ACGTACGTACGTACGT"` (16 bp)
- `b = "ACGT" + "ACGTACGT"` (12 bp — `a` with a 4 bp internal deletion)

Align them with `Bio.Align.PairwiseAligner` four times:

1. Mode `global`, linear gap penalty `-1` (set `open_gap_score = -1, extend_gap_score = -1`).
2. Mode `global`, linear gap penalty `-2`.
3. Mode `global`, affine gap penalty `open = -5, extend = -1`.
4. Mode `global`, affine gap penalty `open = -10, extend = -1`.

For each, print the score and the alignment. In `homework/notes/p4-affine-vs-linear.md`, write a 200-word commentary on:

- Which configurations produced a single 4 bp gap vs multiple smaller gaps?
- How does the choice of gap-open penalty change the optimal alignment?
- Which configuration would you trust for a real DNA alignment, and why?

**Acceptance criteria.**

- `python homework/p4_affine_vs_linear.py` runs.
- Prints four alignment scores and four alignments.
- `notes/p4-affine-vs-linear.md` is committed with a substantive comparison.

**Hint.** The cost of a single 4 bp gap under each scheme: (1) `-4`; (2) `-8`; (3) `-5 - 3 = -8`; (4) `-10 - 3 = -13`. Sequences with more matches in the alignment make up for the gap cost differently in each case.

**Estimated time.** 60 minutes.

---

## Problem 5 — Benchmark your NumPy NW against Biopython

**Problem statement.** Time your Exercise 1 Needleman-Wunsch on three input sizes:

- 100 bp x 100 bp (random DNA).
- 500 bp x 500 bp.
- 1000 bp x 1000 bp.

For each size, run your implementation and `Bio.Align.PairwiseAligner` (with `mode='global'`, `match_score=1`, `mismatch_score=-1`, `open_gap_score=-2`, `extend_gap_score=-2`). Record wall-clock time with `time.perf_counter`. Take three runs and report the median.

Write `homework/p5_bench.py` and produce a Markdown table in `homework/notes/p5-benchmark.md`:

```
| Input size       | NumPy (median, s) | Biopython (median, s) | Ratio |
|------------------|------------------:|----------------------:|------:|
| 100 x 100        |                   |                       |       |
| 500 x 500        |                   |                       |       |
| 1000 x 1000      |                   |                       |       |
```

**Acceptance criteria.**

- `p5_bench.py` runs end to end.
- Table is filled with **actual numbers from your machine** (not made up).
- A one-paragraph commentary in `notes/p5-benchmark.md` interprets the ratio. Biopython will be substantially faster — by how much, and why?

**Hint.** Biopython's `Bio.Align._pairwisealigner` is a C extension. Expect a 50-200x speed advantage over your pure-Python+NumPy implementation. The mini-project on Saturday will revisit this with a tighter NumPy implementation; do not optimize here.

**Estimated time.** 45 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300-400 word reflection at `homework/notes/week-03-reflection.md` answering:

1. Which felt harder this week — filling in the matrix by hand, or implementing it in NumPy? Why?
2. Did anything you previously believed about pairwise alignment turn out to be off this week? (For example: did you think "alignment" meant a single canonical answer? Did you think BLOSUM62 was the only matrix in active use?)
3. After comparing your implementation to Biopython on a 1000 x 1000 pair, what does the speed gap tell you about the difference between teaching code and production code?
4. What is one thing you would want to learn next that this week did not cover? (Multiple sequence alignment? BWT-based aligners? GPU-accelerated alignment?)

**Acceptance criteria.**

- File exists, 300-400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. Be honest. The mistakes you note here are what you will go back and re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 1 h 0 min |
| 2 | 1 h 0 min |
| 3 | 45 min |
| 4 | 1 h 0 min |
| 5 | 45 min |
| 6 | 30 min |
| **Total** | **~5 h 0 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
