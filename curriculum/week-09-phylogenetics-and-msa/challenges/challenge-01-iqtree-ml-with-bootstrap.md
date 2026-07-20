# Challenge 1 — Maximum-likelihood tree with IQ-TREE 2 and ultrafast bootstrap

> **Reproducibility note.** ML tree-building is randomized at the starting-tree, the topology-search, and the bootstrap layers. Two runs with different seeds can produce trees that disagree on internal branches. Pin the seed (`--seed 42`), pin the model (`-m GTR+I+G`), pin the IQ-TREE 2 version, and the result becomes reproducible. The same Newick file from a different seed is *also* defensible, but it is a *different* result.

**Estimated time:** 2.5 hours.
**Goal:** Replace the neighbor-joining tree from Exercise 3 with a maximum-likelihood tree built by IQ-TREE 2 (Minh et al. 2020) with ultrafast bootstrap (Hoang et al. 2018), under the GTR+I+G substitution model. Render the tree with branches coloured by UFBoot support. Compare the ML tree topology to the NJ tree topology; report any disagreements and which is more defensible.

This challenge is the bridge between "I have a fast NJ tree" and "I have a publishable ML tree with proper branch support." The conceptual lift is: NJ uses distances, ML uses the column-level likelihood under a substitution model, and the result is more defensible at the cost of compute time.

---

## Background — Why ML is preferred over NJ

NJ is a fast distance-based heuristic. It works on a summary statistic (the N-by-N distance matrix) and does not explore the tree space. Two consequences:

1. **Information loss.** Two alignments with the same distance matrix produce the same NJ tree, even if the column-level evolutionary patterns differ. ML retains column-level information and can distinguish them.
2. **No model.** NJ assumes the distance matrix is a faithful summary of the true tree distances. JC and K2P provide a correction for the saturation effect; ML goes further by explicitly modeling the substitution process at every column.

The ML pruning algorithm (Felsenstein 1981) computes the likelihood of a given tree topology with given branch lengths under a substitution model. The heuristic search over topologies (nearest-neighbor interchange, subtree pruning and regrafting) finds high-likelihood trees in tractable time. Modern implementations (IQ-TREE 2, RAxML-NG) make ML routine for thousands of sequences on a laptop.

---

## Task

Build a Python wrapper `iqtree_ml.py` that:

1. Validates the input alignment (a trimmed FASTA from Exercise 1).
2. Calls IQ-TREE 2 with the pinned model, seed, and bootstrap replicate count.
3. Reads the resulting `.treefile`, roots it on a named outgroup, ladderizes.
4. Renders the ML tree as a PNG with branches coloured by UFBoot tier.
5. Compares the ML topology to the NJ topology from Exercise 3 and emits a short report.
6. Writes a `run-info.json` with the IQ-TREE version, the model string, the seed, the UFBoot replicate count, and the rooted ML Newick path.

### Layout

```
crunch-bio-portfolio-<yourhandle>/
└── week-09/
    └── challenge-01/
        ├── README.md             how-to-run + write-up
        ├── env.yml               conda env file
        ├── iqtree_ml.py          the wrapper script
        ├── compare_topologies.py the NJ-vs-ML report
        ├── data/
        │   └── cytb_trimmed.fasta    trimmed alignment from Exercise 1
        └── results/
            ├── iqtree.treefile           IQ-TREE 2 output
            ├── iqtree.iqtree             IQ-TREE 2 log
            ├── tree_ml_rooted.nwk        rooted ML Newick
            ├── tree_ml_rooted.png        rendered ML tree
            ├── nj_vs_ml_report.txt       topology comparison
            └── run-info.json
```

### Required function signatures

In `iqtree_ml.py`:

```python
from __future__ import annotations

from pathlib import Path


def run_iqtree2(
    aligned_fasta: Path,
    out_prefix: Path,
    model: str = "GTR+I+G",
    bootstrap: int = 1000,
    alrt: int = 1000,
    seed: int = 42,
    threads: int = 4,
) -> Path:
    """Run IQ-TREE 2. Returns the path to the .treefile.

    Calls subprocess.run with check=True. Raises CalledProcessError if
    IQ-TREE 2 exits non-zero. Surfaces stderr on failure.
    """
    ...


def root_and_render(
    treefile: Path,
    outgroup: str,
    out_newick: Path,
    out_png: Path,
) -> None:
    """Read the IQ-TREE .treefile, root on outgroup, ladderize, write Newick
    and a PNG with branches coloured by UFBoot tier (>=95 green, 70-94
    blue, <70 red).
    """
    ...


def parse_iqtree_log(log_path: Path) -> dict[str, str]:
    """Extract the IQ-TREE 2 version, the final log-likelihood, the AIC,
    the BIC, the number of free parameters, and the wall-clock time.
    Returns a dict suitable for the run-info JSON.
    """
    ...
```

In `compare_topologies.py`:

```python
def bipartition_set(tree, all_taxa: list[str]) -> set[frozenset[str]]:
    """Return the set of non-trivial bipartitions in the tree."""
    ...


def topological_disagreements(
    nj_tree, ml_tree, all_taxa: list[str],
) -> list[frozenset[str]]:
    """Return the bipartitions present in one tree but not the other."""
    ...
```

---

## Acceptance criteria

- [ ] `iqtree_ml.py` exports `run_iqtree2(...)`, `root_and_render(...)`, and `parse_iqtree_log(...)`.
- [ ] `compare_topologies.py` exports `bipartition_set(...)` and `topological_disagreements(...)`.
- [ ] Running `python iqtree_ml.py --input data/cytb_trimmed.fasta --out-dir results/ --outgroup Ornithorhynchus_anatinus --seed 42` runs end to end in under two minutes on the demo panel.
- [ ] `results/iqtree.treefile` exists and is valid Newick with bootstrap percentages on internal nodes.
- [ ] `results/tree_ml_rooted.nwk` exists and is the rooted ladderized version.
- [ ] `results/tree_ml_rooted.png` exists and shows branches coloured by UFBoot tier.
- [ ] `results/nj_vs_ml_report.txt` exists and lists any bipartitions where the NJ tree (from Exercise 3) and the ML tree disagree.
- [ ] `results/run-info.json` records: `iqtree_version`, `iqtree_model`, `iqtree_seed`, `iqtree_bootstrap_replicates`, `iqtree_alrt_replicates`, `final_log_likelihood`, `aic`, `bic`, `wallclock_seconds`, `outgroup`, `biopython_version`, `python_version`.
- [ ] `README.md` includes a ~300-word discussion that:
  - Lists the bipartitions where NJ and ML disagree (typically 0-2 for the demo panel).
  - States which tree's claim on each disagreement is more defensible and why (ML usually wins, but think about the evidence per branch).
  - Notes the IQ-TREE 2 ModelFinder output if you ran it: which model was selected by AIC, which by BIC.
  - Names two failure modes that the demo panel is *not* susceptible to and one that it might be (the platypus is a long-branch outgroup; long-branch attraction is a plausible concern with limited outgroup diversity).
- [ ] All Python files pass `python3 -m py_compile`.

---

## Hints

1. **Pin the model explicitly.** Do not start with `-m MFP` in production; run ModelFinder once interactively to see what it picks, then pin the chosen model in your script. For the demo cytochrome b panel, ModelFinder typically picks `GTR+F+I+G4`.
2. **Pin the seed.** `--seed 42` is the curriculum default. Same input + same model + same seed = identical tree.
3. **Capture and surface stderr.** IQ-TREE 2 prints model-selection diagnostics and warnings to stderr. Print them through to your script's stderr; do not swallow.
4. **Topology comparison via bipartitions.** Two trees on the same leaf set are topologically identical iff they share the same set of non-trivial bipartitions. The Exercise 3 `get_bipartitions` function is what you want; the difference of the two sets is the disagreement.
5. **Be honest about disagreements.** On the demo cytochrome b panel, NJ and ML usually disagree on at most one internal branch and the disagreement is in the tropical part of the tree (e.g. the (cow, pig) clade's exact placement). If your disagreement count is much larger, suspect a bug.

---

## Stretch goals

- Add `-m MFP` to your IQ-TREE call once at the start, parse the AIC and BIC tables out of the `.iqtree` log, and report which model each criterion selected. AIC and BIC often disagree by one or two parameter counts; report both.
- Replace `Bio.Phylo` rendering with `ete3` and add a clade-colour annotation derived from the leaf metadata (e.g. mammals in one colour, sauropsids in another, fish in a third). The visual signal makes the tree publishable.
- Run RAxML-NG (`raxml-ng --all --msa cytb_trimmed.fasta --model GTR+G --bs-trees 1000 --threads 4 --seed 42`) on the same input and compare the RAxML-NG tree to the IQ-TREE 2 tree. They should agree on essentially every branch; if they do not, one of the two has a parameterization difference worth investigating.

---

## What to commit

Commit the `iqtree_ml.py`, `compare_topologies.py`, the `results/` directory, and the `README.md`. Gitignore the IQ-TREE intermediate files (`*.ckp.gz`, `*.bionj`, `*.mldist`, `*.uniqueseq.phy`) — they are large and reproducible from the inputs and the seed.

The commit message should be specific: `c01: IQ-TREE 2 ML on cytb panel, GTR+I+G, seed 42, 1000 UFBoot, agrees with NJ on 7/8 internal bipartitions`.
