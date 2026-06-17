# Week 9 — Exercise Solutions

> **Reproducibility note.** The expected outputs below assume MAFFT 7.526, Biopython 1.84, matplotlib 3.8, and (optionally) ete3 3.1.3. Numbers will differ slightly on older or newer tool versions. The shape of the answers will not.

Each solution names the file you should write, the function bodies the reference implementation expects, and the expected numbers on the demo cytochrome b panel.

---

## Solution 1 — MAFFT via subprocess

**File:** `exercise-01-mafft-via-subprocess.py`. The starter is already runnable; the work is internalizing what the helper functions do and why each parameter is pinned.

### Expected output

Running:

```bash
python exercise-01-mafft-via-subprocess.py \
    --input data/cytb_vertebrates.fasta \
    --out-dir results/ex01
```

produces:

```
results/ex01/aligned.fasta             10 records, 1140 columns
results/ex01/alignment_summary.tsv     1141 lines (header + 1140 rows)
results/ex01/trimmed.fasta             10 records, 1098 columns
results/ex01/run-info.json             { run_date, mafft_version: 7.526, ... }
```

The MAFFT FFT-NS-2 mode aligns the ten vertebrate cytochrome b sequences into a 1,140-column alignment in ~0.4 s. The trim step (`> 50% gap` columns dropped) removes ~42 columns, leaving 1,098. The exact column counts depend on the MAFFT version and the input sequences; pin both and the numbers stay stable across runs.

### Discussion points

1. **Why `--retree 2 --maxiterate 0` and not `--auto`?** `--auto` picks the algorithm by input size and MAFFT version, which is not reproducible. Pinning the algorithm is the small price of reproducibility.
2. **Why pairwise-deletion in `column_consensus`?** Gaps are not residues; counting them as a fifth character would systematically lower the conservation fraction on indel-rich columns.
3. **Why pin the trim threshold at 0.5?** It is the conventional starting point. Tighter (0.2) drops more columns; looser (0.8) keeps more. Note the threshold in the run-info JSON and re-run with a different value to see how sensitive the alignment is to the choice.

### Sanity check: same input + same MAFFT version = byte-identical output

```bash
python exercise-01-mafft-via-subprocess.py --input data/cytb_vertebrates.fasta --out-dir results/ex01a
python exercise-01-mafft-via-subprocess.py --input data/cytb_vertebrates.fasta --out-dir results/ex01b
diff results/ex01a/aligned.fasta results/ex01b/aligned.fasta && echo "byte-identical"
```

The diff should be empty. If it is not, your MAFFT install has non-deterministic threading; drop `--thread 4` to `--thread 1` and re-test.

---

## Solution 2 — Distance matrix + NJ tree

**File:** `exercise-02-distance-and-nj-tree.py`. The reference distance functions are already implemented; the work is understanding the JC and K2P formulas and verifying the NJ output.

### Expected output

```bash
python exercise-02-distance-and-nj-tree.py \
    --input results/ex01/trimmed.fasta \
    --out-dir results/ex02 \
    --outgroup Ornithorhynchus_anatinus
```

produces:

```
results/ex02/distance_matrix_jc.tsv   11 lines (header + 10 rows; 11 cols)
results/ex02/distance_matrix_k2p.tsv  similarly
results/ex02/tree_jc_nj.nwk           Newick with branch lengths
results/ex02/tree_k2p_nj.nwk          similarly
results/ex02/tree_k2p_nj.nex          Nexus with TAXA + TREES blocks
results/ex02/run-info.json            { distance_methods: "jc,k2p", outgroup, ... }
```

Expected JC distance highlights on the demo cytochrome b panel:

| Pair | p-distance | JC distance | K2P distance |
|------|------------|-------------|--------------|
| human vs chimpanzee | 0.014 | 0.014 | 0.015 |
| human vs mouse | 0.165 | 0.184 | 0.192 |
| human vs zebrafish | 0.273 | 0.330 | 0.349 |
| platypus vs chicken | 0.299 | 0.371 | 0.395 |

JC and K2P agree to within 0.01 at low `p` and diverge by ~5-7% at the highest pairwise distances in this panel. Both stay well below saturation (`p < 0.4`); the corrections are doing useful work but not heroic work.

### Discussion points

1. **Why does JC need clamping?** As `p -> 0.75`, the log argument goes to zero and the distance diverges. Clamping at `p = 0.7499` returns `inf`; we then substitute a finite large value (10.0) in the matrix so NJ does not produce a `nan` branch length.
2. **Why is the NJ tree initially unrooted?** NJ does not place a root; it produces a star-merged binary tree. The `root_with_outgroup` step picks the root and converts the tree into a directed binary tree.
3. **Why ladderize?** Two NJ runs with different leaf-traversal orders can produce visually different but topologically identical trees. Ladderizing enforces a deterministic visual layout.

### Sanity check: same alignment + same outgroup = identical Newick

```bash
python exercise-02-distance-and-nj-tree.py --input results/ex01/trimmed.fasta --out-dir results/ex02a --outgroup Ornithorhynchus_anatinus
python exercise-02-distance-and-nj-tree.py --input results/ex01/trimmed.fasta --out-dir results/ex02b --outgroup Ornithorhynchus_anatinus
diff results/ex02a/tree_k2p_nj.nwk results/ex02b/tree_k2p_nj.nwk && echo "byte-identical"
```

---

## Solution 3 — Bootstrap + rendering

**File:** `exercise-03-render-tree-with-bootstrap.py`. The reference implements column resampling, bipartition matching for branch support, ASCII rendering, and matplotlib PNG rendering.

### Expected output

```bash
python exercise-03-render-tree-with-bootstrap.py \
    --input results/ex01/trimmed.fasta \
    --out-dir results/ex03 \
    --outgroup Ornithorhynchus_anatinus \
    --replicates 100 \
    --seed 42
```

produces:

```
results/ex03/tree_k2p_nj_with_support.nwk     internal node names are bootstrap percentages
results/ex03/tree_k2p_nj_with_support.png     matplotlib render
results/ex03/tree_k2p_nj_with_support_ete3.png  ete3 render (if ete3 available)
results/ex03/tree_ascii.txt                   ASCII rendering
results/ex03/run-info.json                    { seed: 42, n_replicates: 100, ... }
```

Expected bootstrap percentages on the cytochrome b panel (100 replicates, seed 42):

| Clade | Expected support |
|-------|-----------------:|
| (human, chimpanzee) | 100 |
| (mouse, rat) | 98-100 |
| (cow, pig) | 90-100 |
| ((human, chimp), (mouse, rat)) | 70-90 |
| ((eutherians), (cow, pig)) | 60-85 |
| (zebrafish + frog) | 50-80 |

The bootstrap is stochastic. With `--seed 42` and `--replicates 100`, two runs on the same input produce identical percentages (the seed makes the RNG deterministic). With a different seed, the numbers shift by ~5 percentage points but the topology is unchanged.

### Discussion points

1. **Why column resampling with replacement?** Felsenstein 1985 argues that alignment columns are the independent observations in phylogenetic inference. Resampling columns with replacement produces a new alignment that has the same statistical structure as the original; the variability across resamples estimates the sampling variability in the tree.
2. **Why root each bootstrap replicate before bipartition counting?** Rooting fixes the orientation of every clade so that "the clade containing human and chimp" means the same thing across replicates. Without rooting, the same topology can present different bipartitions depending on the traversal direction.
3. **Why is the support sometimes 100?** Two-leaf clades that are uncontroversial in the alignment (e.g. human-chimp) appear in every resampled tree because the signal is overwhelming. The 100% number is not "the truth"; it is "the bootstrap cannot distinguish this branch from the data."

### Sanity check: same seed = identical support

```bash
python exercise-03-render-tree-with-bootstrap.py --input results/ex01/trimmed.fasta --out-dir results/ex03a --outgroup Ornithorhynchus_anatinus --seed 42
python exercise-03-render-tree-with-bootstrap.py --input results/ex01/trimmed.fasta --out-dir results/ex03b --outgroup Ornithorhynchus_anatinus --seed 42
diff results/ex03a/tree_k2p_nj_with_support.nwk results/ex03b/tree_k2p_nj_with_support.nwk && echo "byte-identical"
```

---

## A note on running these without MAFFT or ete3 installed

The Python files compile under `python3 -m py_compile` regardless of which tools are on the PATH (Biopython, MAFFT, ete3 are imported lazily inside the functions that need them). To *run* the exercises you do need:

- MAFFT 7 on the PATH for Exercise 1.
- Biopython 1.84 (`pip install biopython`) for Exercises 1, 2, 3.
- matplotlib 3.8 (`pip install matplotlib`) for Exercise 3's PNG output.
- ete3 3.1.3 (optional; `pip install ete3` or `conda install -c bioconda ete3`) for Exercise 3's coloured PNG output.

Conda one-liner:

```bash
conda install -c bioconda -c conda-forge \
    python=3.11 mafft=7.526 biopython=1.84 matplotlib=3.8 ete3=3.1.3
```

---

## What to commit

After running all three exercises, your `week-09/exercises/` directory should contain:

```
exercises/
    exercise-01-mafft-via-subprocess.py
    exercise-02-distance-and-nj-tree.py
    exercise-03-render-tree-with-bootstrap.py
    SOLUTIONS.md
    results/
        ex01/{aligned.fasta, trimmed.fasta, alignment_summary.tsv, run-info.json}
        ex02/{distance_matrix_*.tsv, tree_*.nwk, tree_*.nex, run-info.json}
        ex03/{tree_*.nwk, tree_*.png, tree_ascii.txt, run-info.json}
```

Commit the `results/` directory. The output files are small (well under 1 MB total) and the diffs are informative when you revisit the exercises later.
