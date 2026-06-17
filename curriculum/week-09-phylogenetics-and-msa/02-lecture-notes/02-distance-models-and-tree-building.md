# Lecture 2 — Distance Models, Neighbor-Joining, and Maximum Likelihood

> **Reproducibility note.** Tree-building algorithms are randomized at several layers (bootstrap resampling, NJ tie-breaking, ML starting trees, ML branch-length optimization). Same input + same algorithm + different seed = different tree topology in the worst case. Always pin the seed and emit it in the run-info JSON. The Week 9 default is `seed = 42`.

> **Duration:** ~3 hours of reading + a brief NJ sanity check in Biopython.
> **Outcome:** You can compute a Jukes-Cantor and a Kimura two-parameter distance matrix from an alignment, build a neighbor-joining tree with `Bio.Phylo.TreeConstruction.DistanceTreeConstructor`, run IQ-TREE 2 from Python for a maximum-likelihood tree, generate bootstrap support either by hand (NJ + column resampling) or via IQ-TREE 2's ultrafast bootstrap, and explain why high bootstrap is not the same as correctness.

If you only remember one thing from this lecture, remember this:

> **A tree-building algorithm takes the alignment as the truth and emits its best estimate of the relationships. The estimate is wrong with probability that depends on the alignment quality, the model fit, the algorithm choice, and the random seed. Bootstrap support tells you how robust the tree is to column resampling, not how true it is. The defence is to pin the algorithm, pin the model, pin the seed, run the bootstrap, and report all four in the figure caption.**

Lecture 1 left off at the cleaned alignment. Lecture 2 takes that alignment, computes distances, and builds the tree.

---

## 1. Where we are in the pipeline

```
ALIGNED FASTA (cleaned; N x M' columns) ->
        distance matrix (JC, K2P, F84, or ML-derived; N x N) ->
        tree-building algorithm (UPGMA, NJ, or ML) ->
        TREE (Newick + branch lengths + optional bootstrap support).
```

The distance matrix is the bridge between the alignment and the tree under distance-based methods (UPGMA, NJ). Maximum likelihood skips the distance matrix and works directly on the columns, but the conceptual order is the same: columns -> distances or model parameters -> tree.

---

## 2. The p-distance (Hamming distance)

The simplest distance between two aligned sequences is the **proportion of differing columns**. For two rows in an alignment with `n` columns and `d` columns where the two rows differ (ignoring columns where either row has a gap, by convention), the p-distance is:

```
p = d / n
```

This is the raw Hamming distance, normalized to the alignment length.

For closely related sequences (say, two human cytochrome b paralogs at ~98% identity, p ~ 0.02), the p-distance is a good estimate of the evolutionary distance. For more divergent sequences, the p-distance is biased *downward* because of **multiple hits**: two sequences that each experienced two substitutions at the same column over evolutionary time produce an observed difference of zero at that column. The p-distance counts zero where the true count is two.

The bias gets worse as `p` increases. For DNA, the p-distance asymptotes near 0.75 (three of the four bases differ from any starting base, on average); for protein, near 0.95 (nineteen of the twenty residues differ). Beyond `p ~ 0.3` the bias is severe and a correction is required.

### Computing p-distance in Python

```python
from __future__ import annotations


def p_distance(seq_a: str, seq_b: str) -> float:
    """Return the p-distance (proportion of differing columns) between
    two aligned sequences of equal length. Columns with a gap in either
    sequence are excluded from both the numerator and the denominator
    (pairwise deletion convention).
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("p_distance requires equal-length aligned sequences.")
    differ: int = 0
    total: int = 0
    for ch_a, ch_b in zip(seq_a, seq_b):
        if ch_a == "-" or ch_b == "-":
            continue
        if ch_a == "N" or ch_b == "N":
            continue
        total += 1
        if ch_a != ch_b:
            differ += 1
    if total == 0:
        return 0.0
    return differ / total
```

The pairwise-deletion convention (skip columns where either sequence has a gap) is the conventional choice. The alternative is **complete deletion** (drop any column with a gap in any row), which is cleaner but loses more information; it is reasonable for highly-gapped alignments and the wrong choice for sparsely-gapped ones.

---

## 3. The Jukes-Cantor correction (JC; Jukes and Cantor 1969)

The Jukes-Cantor model assumes the four nucleotides are interchangeable with equal rate. Under this model, the expected number of substitutions per site `d_JC` given the observed `p` is:

```
d_JC = -3/4 * ln(1 - 4/3 * p)
```

The derivation is short:

- Let `mu` be the per-site substitution rate. Under JC, the probability that a site differs after time `t` is `p = 3/4 * (1 - exp(-4/3 * mu * t))`.
- Solve for `mu * t = d`: `d = -3/4 * ln(1 - 4/3 * p)`.

The formula has two practical caveats:

- **Saturation.** As `p -> 3/4`, the log argument goes to zero and `d_JC -> infinity`. In practice, we clamp `p` at 0.7499 and warn that the distance is saturated. Sequences with `p > 0.7` are essentially uninformative for distance-based phylogenetics; the substitution process has erased the signal.
- **The equal-rate assumption.** JC assumes every substitution (e.g. A -> G, A -> C, A -> T) happens at the same rate. In real genomes, transitions (A <-> G, C <-> T, the within-purine and within-pyrimidine substitutions) happen 2-4x faster than transversions (across the purine / pyrimidine boundary). K2P fixes this; JC does not.

### Computing JC in Python

```python
import math


def jukes_cantor_distance(p: float) -> float:
    """Return the Jukes-Cantor corrected distance from an observed p-distance.

    Clamps p at 0.7499 to avoid a log of zero. Returns float('inf') if
    the input is already at or above the saturation threshold (the caller
    should warn the user that the distance is no longer reliable).
    """
    if p < 0.0:
        raise ValueError(f"p-distance cannot be negative; got {p}.")
    if p >= 0.7499:
        return float("inf")
    return -0.75 * math.log(1.0 - 4.0 / 3.0 * p)
```

For the demo cytochrome b panel, JC distances range from 0.014 (human-chimp) to ~0.31 (platypus-fish). The 0.31 is well within the linear regime (no saturation) but high enough that the JC correction adds ~0.04 to the raw `p`.

---

## 4. The Kimura two-parameter correction (K2P; Kimura 1980)

Citation: Kimura M. "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." *Journal of Molecular Evolution* 16:111 (1980).

K2P splits the observed difference into:

- `P` — proportion of transition differences (A <-> G or C <-> T columns).
- `Q` — proportion of transversion differences (everything else).

Note `p = P + Q`. The K2P distance is:

```
d_K2P = -1/2 * ln(1 - 2*P - Q) - 1/4 * ln(1 - 2*Q)
```

The formula is more accurate than JC for inputs where the transition / transversion ratio is far from 1.0 (which is most real DNA, especially mitochondrial DNA where the ratio is often 4-8). For vertebrate cytochrome b, K2P distances are typically 5-15% larger than JC distances.

### Computing K2P in Python

```python
import math


def kimura_2p_distance(seq_a: str, seq_b: str) -> float:
    """Return the K2P-corrected distance between two aligned sequences.

    Splits observed differences into transition (P) and transversion (Q)
    proportions. Returns float('inf') if either log argument is non-positive
    (the saturation regime).
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("kimura_2p_distance requires equal-length aligned sequences.")
    transitions: int = 0
    transversions: int = 0
    total: int = 0
    purines: set[str] = {"A", "G"}
    pyrimidines: set[str] = {"C", "T"}
    for ch_a, ch_b in zip(seq_a.upper(), seq_b.upper()):
        if ch_a not in "ACGT" or ch_b not in "ACGT":
            continue
        total += 1
        if ch_a == ch_b:
            continue
        a_purine: bool = ch_a in purines
        b_purine: bool = ch_b in purines
        if a_purine == b_purine:
            transitions += 1
        else:
            transversions += 1
    if total == 0:
        return 0.0
    big_p: float = transitions / total
    big_q: float = transversions / total
    arg1: float = 1.0 - 2.0 * big_p - big_q
    arg2: float = 1.0 - 2.0 * big_q
    if arg1 <= 0.0 or arg2 <= 0.0:
        return float("inf")
    return -0.5 * math.log(arg1) - 0.25 * math.log(arg2)
```

K2P is the Week 9 default for vertebrate nucleotide alignments. For protein alignments, use a substitution matrix-based distance (PAM, BLOSUM, or LG) via Biopython's `DistanceCalculator("blosum62")` or `DistanceCalculator("identity")` followed by a matrix-specific correction; the K2P formula does not apply to protein.

---

## 5. The distance matrix

The pairwise distances between all N sequences form an N-by-N symmetric matrix with zeros on the diagonal:

```
       human  chimp  mouse  ...
human    0.0  0.014  0.165  ...
chimp  0.014    0.0  0.169  ...
mouse  0.165  0.169    0.0  ...
...
```

Biopython's `Bio.Phylo.TreeConstruction.DistanceCalculator` returns a `DistanceMatrix` object that wraps the lower-triangle as a list-of-lists; the access pattern is `matrix[i, j]` where `i, j` are taxon names. For the Week 9 exercises we build the matrix by hand (so you understand the data shape) and then convert to Biopython's `DistanceMatrix` for the NJ step.

```python
from __future__ import annotations


def build_distance_matrix(
    aligned_fasta_records: list[tuple[str, str]],
    method: str = "k2p",
) -> "DistanceMatrix":
    """Return a Biopython DistanceMatrix computed under method ('jc' or 'k2p').

    aligned_fasta_records is a list of (taxon_name, aligned_sequence) tuples.
    All sequences must be the same length.
    """
    from Bio.Phylo.TreeConstruction import DistanceMatrix

    names: list[str] = [name for name, _ in aligned_fasta_records]
    n: int = len(aligned_fasta_records)
    # DistanceMatrix expects a lower-triangle list-of-lists, including the
    # zero diagonal.
    matrix_rows: list[list[float]] = []
    for i in range(n):
        row: list[float] = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
                continue
            seq_a: str = aligned_fasta_records[i][1]
            seq_b: str = aligned_fasta_records[j][1]
            if method == "jc":
                p: float = p_distance(seq_a, seq_b)
                row.append(jukes_cantor_distance(p))
            elif method == "k2p":
                row.append(kimura_2p_distance(seq_a, seq_b))
            else:
                raise ValueError(f"Unknown method: {method!r}")
        matrix_rows.append(row)
    return DistanceMatrix(names=names, matrix=matrix_rows)
```

This function is what Exercise 2 implements.

---

## 6. UPGMA (Sokal and Michener 1958)

The Unweighted Pair-Group Method with Arithmetic Mean is the simplest tree-building algorithm:

1. Find the pair of taxa with the smallest distance.
2. Merge them into a cluster; the cluster's branch length to each child is half the pair distance.
3. Recompute the distance from the new cluster to every other taxon as the arithmetic mean.
4. Repeat until one cluster remains.

UPGMA produces a **rooted, ultrametric tree** — every leaf is the same distance from the root. This is correct *only* if the sequences evolved under a strict molecular clock (constant substitution rate across all lineages). Real biological sequences almost never satisfy the strict clock; UPGMA trees are systematically wrong when lineages evolve at different rates.

We mention UPGMA so you can recognize it in older papers. For Week 9's exercises we use neighbor-joining, which does not assume a clock.

---

## 7. Neighbor-joining (NJ; Saitou and Nei 1987)

Citation: Saitou N, Nei M. "The neighbor-joining method: a new method for reconstructing phylogenetic trees." *Molecular Biology and Evolution* 4:406 (1987). Free at <https://academic.oup.com/mbe/article/4/4/406/1029664>.

Neighbor-joining is the canonical fast distance-based tree builder. It does **not** assume a clock; the resulting tree is unrooted. The algorithm:

1. Start with a star tree (one central node, N leaves).
2. Compute a corrected distance for each pair: `Q(i, j) = (N - 2) * d(i, j) - sum_k d(i, k) - sum_k d(j, k)`.
3. Pick the pair with the smallest `Q` and merge them into a new internal node.
4. Compute the new node's distance to every remaining taxon as `d(u, k) = 1/2 * (d(i, k) + d(j, k) - d(i, j))`.
5. Repeat until two nodes remain; join with a single branch.

The output is an unrooted binary tree with branch lengths. Time complexity is `O(N^3)`; for N = 10 this is instantaneous; for N = 1,000 it is a few seconds.

NJ is implemented in `Bio.Phylo.TreeConstruction.DistanceTreeConstructor`:

```python
from __future__ import annotations


def build_nj_tree(distance_matrix: "DistanceMatrix") -> "Tree":
    """Build a neighbor-joining tree from a Biopython DistanceMatrix.

    The tree is unrooted by NJ convention; the caller may root it on
    a named outgroup with tree.root_with_outgroup(name).
    """
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

    constructor = DistanceTreeConstructor()
    return constructor.nj(distance_matrix)
```

To root on a named outgroup (the platypus is the conventional outgroup for the eutherian + bird + fish + frog cytochrome b panel):

```python
tree.root_with_outgroup("Ornithorhynchus_anatinus")
```

After rooting, the tree is binary, has explicit branch lengths, and is ready to write as Newick.

### NJ tie-breaking

When two pairs have identical `Q` values, NJ picks one by an internal tie-breaking rule. Biopython's `DistanceTreeConstructor` uses a deterministic rule (first-pair-found-wins); some other implementations randomize. For Week 9 we use Biopython's deterministic NJ, which means same input -> byte-identical Newick on a single Biopython version. Different Biopython major versions have changed the tie-breaking once (1.78 -> 1.79); pin the version.

---

## 8. Maximum likelihood (ML; Felsenstein 1981)

Citation: Felsenstein J. "Evolutionary trees from DNA sequences: a maximum likelihood approach." *Journal of Molecular Evolution* 17:368 (1981). Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>. **The most important paper of the week.** Read it.

ML phylogenetics frames tree-building as a statistical inference problem. Given an alignment and a substitution model (JC, K2P, GTR, etc.), the likelihood of a given tree topology `T` with branch lengths `b` is the probability of observing the alignment under that tree. ML picks the `(T, b)` that maximizes the likelihood.

Felsenstein's 1981 contribution was the **pruning algorithm**: for a given tree, the likelihood can be computed in time linear in the tree size by recursively combining per-column conditional probabilities from the leaves to the root. The exhaustive search over tree topologies is still intractable (the number of topologies grows super-exponentially), but heuristic search (nearest-neighbor interchanges, subtree pruning and regrafting) finds high-likelihood trees in tractable time.

Two modern fast ML implementations:

- **IQ-TREE 2** (Minh et al. 2020, *Molecular Biology and Evolution* 37:1530; free at <https://academic.oup.com/mbe/article/37/5/1530/5721363>; docs at <http://www.iqtree.org/doc/>).
- **RAxML-NG** (Kozlov et al. 2019, *Bioinformatics* 35:4453; free at <https://academic.oup.com/bioinformatics/article/35/21/4453/5487384>; docs at <https://github.com/amkozlov/raxml-ng/wiki>).

Both are free, both are conda-installable, both produce essentially the same tree on the demo input. We default to IQ-TREE 2 because the CLI is slightly cleaner and the ultrafast bootstrap (UFBoot, Hoang et al. 2018) is faster than the standard bootstrap on small alignments.

### The canonical IQ-TREE 2 call

```bash
iqtree2 \
    -s aligned.fasta \
    -m GTR+I+G \
    -B 1000 \
    -alrt 1000 \
    -T 4 \
    -seed 42 \
    -pre results/iqtree \
    -redo
```

Flag-by-flag:

- `-s aligned.fasta` — input alignment.
- `-m GTR+I+G` — the substitution model. GTR is the general-time-reversible model (six free rate parameters); `+I` adds a proportion of invariable sites; `+G` adds gamma-distributed rate heterogeneity. For nucleotide data this is the conventional default; for protein data swap to `LG+G`.
- `-B 1000` — ultrafast bootstrap with 1,000 replicates.
- `-alrt 1000` — SH-like approximate likelihood ratio test with 1,000 replicates. Provides a second, complementary support measure on every internal branch.
- `-T 4` — four threads.
- `-seed 42` — random seed for the ML search and the bootstrap.
- `-pre results/iqtree` — output file prefix.
- `-redo` — overwrite existing output files.

For the demo cytochrome b panel this runs in ~15-30 seconds. The output files:

- `results/iqtree.treefile` — the ML tree in Newick, with bootstrap percentages on internal nodes.
- `results/iqtree.iqtree` — the human-readable log: model, log-likelihood, AIC, BIC, per-branch UFBoot and aLRT support values.
- `results/iqtree.log` — the run log.
- `results/iqtree.ckp.gz`, `results/iqtree.bionj`, `results/iqtree.mldist` — intermediate files. Safe to delete; reproducible from inputs + seed.

### Calling IQ-TREE 2 from Python

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def run_iqtree2(
    aligned_fasta: Path,
    out_prefix: Path,
    model: str = "GTR+I+G",
    bootstrap: int = 1000,
    seed: int = 42,
    threads: int = 4,
) -> Path:
    """Run IQ-TREE 2 with ultrafast bootstrap. Returns the path to the .treefile."""
    cmd = [
        "iqtree2",
        "-s", str(aligned_fasta),
        "-m", model,
        "-B", str(bootstrap),
        "-alrt", str(bootstrap),
        "-T", str(threads),
        "-seed", str(seed),
        "-pre", str(out_prefix),
        "-redo",
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    treefile: Path = out_prefix.with_suffix(out_prefix.suffix + ".treefile")
    if not treefile.exists():
        raise RuntimeError(f"IQ-TREE 2 finished but {treefile} was not written.")
    return treefile
```

### Model selection: ModelFinder

If you do not know which substitution model to use, run IQ-TREE 2 once with `-m MFP` (ModelFinder Plus; Kalyaanamoorthy et al. 2017, *Nature Methods* 14:587). It tries ~80 models and picks the best one by AIC / BIC. For the demo cytochrome b panel, ModelFinder typically picks `GTR+F+I+G4`. Once you know the best model, pin it explicitly in future runs (`-m GTR+F+I+G4`) so the tree is reproducible without the model search.

---

## 9. Bootstrap support (Felsenstein 1985)

Citation: Felsenstein J. "Confidence limits on phylogenies: an approach using the bootstrap." *Evolution* 39:783 (1985). Free at <https://www.jstor.org/stable/2408678>.

The bootstrap protocol:

1. Resample the alignment columns with replacement, keeping the row count fixed. Produce a new alignment of the same shape.
2. Rebuild the tree from the resampled alignment.
3. Repeat 1-2 for `B` bootstrap replicates (usually 100-1,000).
4. For each internal branch of the original tree, count the fraction of replicate trees that contain the same bipartition. That fraction (times 100) is the bootstrap percentage.

A 70% bootstrap on a branch means: 70% of the resampled alignments produced a tree that included this same internal branch. It is a measure of **robustness to column resampling**, not a measure of how true the branch is in the underlying species tree.

### Bootstrap caveats

- **High support is not the same as correctness.** A systematic bias in the alignment or the model can produce a high-bootstrap branch that is wrong. Long-branch attraction (Felsenstein 1978) is the canonical example: two rapidly-evolving sequences cluster together because they have accumulated similar (random) substitutions, with high bootstrap support, when the true tree separates them.
- **Low support is informative.** A branch with < 50% bootstrap is essentially unsupported; the data does not distinguish it from the alternative topologies. Do not report a topology that depends on a low-bootstrap branch as a confident result.
- **The bootstrap is on the column sample, not on the model.** If your substitution model is wrong, the bootstrap will not catch it. Run with a couple of alternative models (`-m JC`, `-m GTR+I+G`) and check that the bootstrap-supported branches agree.

### Standard vs ultrafast bootstrap

The classical Felsenstein 1985 bootstrap runs the full tree search on each resampled alignment. For 1,000 replicates, this is 1,000x the cost of a single tree search.

**UFBoot** (ultrafast bootstrap, Hoang et al. 2018, *Molecular Biology and Evolution* 35:518; free at <https://academic.oup.com/mbe/article/35/2/518/4565479>) is a much faster approximation that runs the full search once and then evaluates each resampled alignment on a small set of candidate trees. UFBoot percentages are calibrated differently from standard bootstrap; a UFBoot 95% is roughly equivalent to a standard bootstrap 70%. The IQ-TREE 2 `-B 1000` flag invokes UFBoot.

For Week 9 we use UFBoot in the challenge and the mini-project. For the basic exercises (NJ tree without IQ-TREE), we implement a hand-coded bootstrap of 100 replicates, which is slow but transparent.

### Hand-coded bootstrap for NJ

```python
from __future__ import annotations

import random


def bootstrap_nj_tree(
    aligned_records: list[tuple[str, str]],
    n_replicates: int,
    distance_method: str = "k2p",
    seed: int = 42,
) -> list["Tree"]:
    """Generate n_replicates bootstrap NJ trees.

    For each replicate, resample the alignment columns with replacement
    (keeping the row count fixed), rebuild the distance matrix under
    distance_method, and build an NJ tree. Return the list of trees.
    """
    rng: random.Random = random.Random(seed)
    n_cols: int = len(aligned_records[0][1])
    trees: list["Tree"] = []
    for _ in range(n_replicates):
        col_indices: list[int] = [rng.randrange(n_cols) for _ in range(n_cols)]
        resampled: list[tuple[str, str]] = [
            (name, "".join(seq[k] for k in col_indices))
            for name, seq in aligned_records
        ]
        matrix = build_distance_matrix(resampled, method=distance_method)
        tree = build_nj_tree(matrix)
        trees.append(tree)
    return trees
```

The bootstrap-support computation (matching internal-branch bipartitions across replicates) is implemented in `Bio.Phylo.Consensus.get_support`. Exercise 3 puts the pieces together.

---

## 10. Long-branch attraction and other failure modes

- **Long-branch attraction** (Felsenstein 1978, *Systematic Biology* 27:401). Two sequences with much higher substitution rates than the rest of the panel can cluster together with high support, even when their true position in the tree is elsewhere. The fix is to use a model that captures rate heterogeneity (`+G` for gamma-distributed rates) and to include slow-evolving outgroup sequences that break up the attraction. For Week 9's demo panel this is unlikely to bite; for primate-vs-rodent vs reptile vs fish panels with large rate differences it can.
- **Saturation.** At high p-distance (> 0.5 for DNA), the JC and K2P corrections start to break down — multiple hits are common, and the inferred distance has high variance. Beyond `p > 0.7`, the corrections are useless. The fix is to use a more sophisticated model (GTR+I+G) or to switch to protein-level analysis (translate to amino acids; protein evolves more slowly).
- **Gene tree vs species tree.** A single-gene tree (e.g. cytochrome b) is not the same as the species tree. Incomplete lineage sorting, horizontal gene transfer, hybridization, and gene duplication all create gene trees that disagree with the species tree. For Week 9's vertebrate panel the disagreements are small (the cytochrome b tree closely matches the standard mammalian tree); for bacterial trees, gene-tree discordance is the rule and not the exception.
- **Alignment-induced artefacts.** A poorly-aligned column inflates the apparent substitution count and pulls the tree. The fix is to trim aggressively (Lecture 1) or to use a structure-aware alignment in the difficult regions.
- **Seed dependence.** ML searches and bootstrap replicate generation are randomized; the seed matters. Two IQ-TREE 2 runs with different seeds can produce trees that disagree on internal branches. Always pin the seed.

We name these failure modes in the mini-project write-up and ask which of them your tree could plausibly suffer from.

---

## 11. What to take to the rest of the week

By the end of Lecture 2 you should be able to:

- Compute a p-distance, a JC distance, and a K2P distance between two aligned sequences, with appropriate handling of gaps, ambiguous symbols, and saturation.
- Build an N-by-N distance matrix from an alignment and convert it to a Biopython `DistanceMatrix`.
- Run `Bio.Phylo.TreeConstruction.DistanceTreeConstructor.nj` to get an unrooted NJ tree and root it on a named outgroup.
- Call IQ-TREE 2 from Python with a pinned model, seed, and bootstrap replicate count.
- Generate bootstrap replicates by hand (for NJ) or via UFBoot (for ML) and attach the support percentages to internal branches.
- State the distinction between bootstrap support (robustness to column resampling) and topological correctness (which is unknowable from a single alignment).

Lecture 3 takes the tree and addresses the file-format and rendering layer.

## References

- Jukes TH, Cantor CR. "Evolution of protein molecules." In: Munro HN (ed.), *Mammalian Protein Metabolism*, vol. III, 21-132. Academic Press, New York (1969). The JC model. Behind a paywall as a book chapter; the JC formula is reproduced in every modern molecular-evolution textbook.
- Kimura M. "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." *Journal of Molecular Evolution* 16:111 (1980). The K2P model.
- Felsenstein J. "Cases in which parsimony or compatibility methods will be positively misleading." *Systematic Biology* 27:401 (1978). The long-branch attraction paper.
- Felsenstein J. "Evolutionary trees from DNA sequences: a maximum likelihood approach." *Journal of Molecular Evolution* 17:368 (1981). Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>. The ML paper.
- Felsenstein J. "Confidence limits on phylogenies: an approach using the bootstrap." *Evolution* 39:783 (1985). Free at <https://www.jstor.org/stable/2408678>. The bootstrap paper.
- Saitou N, Nei M. "The neighbor-joining method: a new method for reconstructing phylogenetic trees." *Molecular Biology and Evolution* 4:406 (1987). Free at <https://academic.oup.com/mbe/article/4/4/406/1029664>. The NJ paper.
- Sokal RR, Michener CD. "A statistical method for evaluating systematic relationships." *University of Kansas Science Bulletin* 38:1409 (1958). The UPGMA paper.
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. "ModelFinder: fast model selection for accurate phylogenetic estimates." *Nature Methods* 14:587 (2017). The ModelFinder paper.
- Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. "UFBoot2: improving the ultrafast bootstrap approximation." *Molecular Biology and Evolution* 35:518 (2018). Free at <https://academic.oup.com/mbe/article/35/2/518/4565479>. The ultrafast bootstrap paper.
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. "IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era." *Molecular Biology and Evolution* 37:1530 (2020). Free at <https://academic.oup.com/mbe/article/37/5/1530/5721363>.
- Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. "RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference." *Bioinformatics* 35:4453 (2019). Free at <https://academic.oup.com/bioinformatics/article/35/21/4453/5487384>. The RAxML-NG paper.
