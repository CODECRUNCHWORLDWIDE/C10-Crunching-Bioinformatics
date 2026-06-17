# Challenge 2 — Newick and Nexus round-trip with topology preservation

> **Reproducibility note.** File-format round-trips are a quiet source of phylogenetic drift. A Newick writer that drops internal node labels, a Nexus reader that strips branch lengths under 1e-6, a PhyloXML import that re-orders siblings — each is a small, defensible choice that nonetheless changes the file the next downstream tool sees. The defensive pattern is: round-trip every tree, prove topology is preserved, and document the format choices in the run-info JSON.

**Estimated time:** 2 hours.
**Goal:** Take any of the trees built in the exercises or Challenge 1 (NJ Newick, ML Newick, or a Nexus from Exercise 2), round-trip it through Newick -> Nexus -> Newick and Nexus -> Newick -> Nexus, prove that the leaf set and the internal bipartitions are preserved, and report any per-branch annotations (bootstrap support, branch length precision) that are lost in the trip.

---

## Background — Why a round-trip can lose information

The phylogenetic file formats encode the same logical object — a tree — but each format has a slightly different surface area:

- **Newick** carries: parenthetical topology, branch lengths (optional), internal node labels (optional, usually bootstrap), leaf names.
- **Nexus** carries: a TAXA block with taxon labels and (optionally) attributes, a TREES block with one or more Newick-flavored trees, optional CHARACTERS block with the alignment.
- **PhyloXML** carries: structured per-node annotations (species, sequence, accession, geographic location, custom key-value pairs) that Newick and Nexus cannot represent.

The round-trip points of failure:

1. **Branch length precision.** Newick writers usually round to 6 decimal places. Nexus readers may round further. Two round-trips can shift a branch length by 1e-7.
2. **Internal node labels (bootstrap).** Some Newick writers drop integer internal labels; others preserve them as strings. Biopython's writer drops them unless you set the `confidence_as_branch_length` flag.
3. **Sibling order.** A tree `((A,B),C);` is topologically identical to `(C,(A,B));`, but the Newick strings are not byte-identical. After a round-trip the sibling order may flip.
4. **Multi-tree files.** A Nexus file with three TREES entries should round-trip as three Newick entries; some implementations only emit the last.

Topological preservation means: same leaf set, same set of non-trivial bipartitions. Byte-identical means: same characters in the output file. Topological preservation is the strict requirement; byte-identicality is a nice-to-have that may not hold across format transitions.

---

## Task

Build a Python module `roundtrip.py` with the following structure:

```python
from __future__ import annotations

from pathlib import Path


def read_newick(path: Path):
    """Read a Newick tree. Returns a Bio.Phylo.BaseTree.Tree."""
    ...


def read_nexus(path: Path):
    """Read a Nexus tree. Returns a Bio.Phylo.BaseTree.Tree."""
    ...


def write_newick(tree, path: Path) -> None:
    """Write a tree to Newick."""
    ...


def write_nexus(tree, path: Path) -> None:
    """Write a tree to Nexus."""
    ...


def get_bipartitions(tree, all_taxa: list[str]) -> set[frozenset[str]]:
    """Return the set of non-trivial bipartitions in the tree."""
    ...


def topology_equal(tree_a, tree_b) -> bool:
    """Return True if two trees have the same leaf set and bipartitions."""
    ...


def roundtrip_newick_via_nexus(
    in_newick: Path,
    out_dir: Path,
) -> dict:
    """Round-trip: read Newick -> write Nexus -> read Nexus -> write Newick.

    Returns a dict with:
      - 'leaves_preserved': bool
      - 'bipartitions_preserved': bool
      - 'branch_lengths_max_drift': float (max abs diff after round-trip)
      - 'support_values_preserved': bool (True if all internal node confidence
        values survive the trip with absolute error < 1.0)
      - 'sibling_order_preserved': bool (True if Newick output is byte-equal
        to a canonical-ordered version of the input)
    """
    ...
```

The orchestrator function `main` should:

1. Take an input Newick path and an output directory.
2. Run `roundtrip_newick_via_nexus(input, out_dir)`.
3. Also run a Nexus -> Newick -> Nexus round-trip if the input has a corresponding Nexus.
4. Write a `roundtrip_report.json` summarizing both round-trips.
5. Print a pass/fail summary to stdout.

---

## Acceptance criteria

- [ ] `roundtrip.py` exports the function signatures listed above.
- [ ] Running `python roundtrip.py --input results/ex02/tree_k2p_nj.nwk --out-dir results/ch02/` runs end to end in under 5 seconds.
- [ ] `results/ch02/roundtrip_report.json` exists with the dict shape described above.
- [ ] On the demo cytochrome b NJ tree, the report shows: `leaves_preserved = True`, `bipartitions_preserved = True`, `branch_lengths_max_drift < 1e-5`, `support_values_preserved = True`.
- [ ] On the IQ-TREE 2 ML tree from Challenge 1 (if you completed it), the same round-trip succeeds with the same outcomes.
- [ ] `README.md` for the challenge includes a short discussion (200-300 words) that:
  - Names which fields of the tree survive a round-trip and which do not.
  - Reports the maximum branch-length drift you observed and explains its cause (float precision in the Newick writer).
  - Notes any sibling-order changes (Biopython's reader and writer normalize sibling order; this is reproducible but not byte-equal to a Newick from a different writer).
  - States the practical implication: **any pipeline that depends on byte-identical Newick across tools must pin the writer**. Use Biopython for both write and read, or use `dendropy` for both, but do not mix.
- [ ] The Python file passes `python3 -m py_compile`.

---

## Hints

1. **Use `Bio.Phylo.read` and `Bio.Phylo.write`** for both formats. The `nexus` reader requires a TAXA block; if the input Nexus does not have one, add it by parsing the Newick first.
2. **Bipartition equality** is the cleanest topology check. The Robinson-Foulds distance (Robinson and Foulds 1981) on two trees is 0 iff the bipartition sets are equal. Biopython does not implement RF; computing the set difference directly is fine.
3. **Branch-length drift** is bounded by the floating-point precision of the writer. Biopython's default Newick writer uses 5 decimal places; the drift after one round-trip is typically < 1e-5.
4. **Support values are stored as `clade.confidence` in Biopython.** Some Newick writers serialize them as the internal node *name*, others as the internal node *confidence*. Check both fields on the round-tripped tree.

---

## Stretch goals

- Add a PhyloXML leg to the round-trip: Newick -> PhyloXML -> Newick. Note any annotations that survive PhyloXML and fail back to Newick (PhyloXML carries species and sequence annotations that Newick cannot represent).
- Compare Biopython's Newick output to `dendropy`'s Newick output on the same tree. The two libraries serialize slightly differently (different precision, different sibling ordering); diff the files and report what changed.
- Round-trip an IQ-TREE 2 `.treefile` (which has both `support` and `aLRT` annotations encoded in a non-standard way) and check that both annotations survive. They usually do not without a specialized parser; document the loss.

---

## What to commit

Commit `roundtrip.py`, the `results/ch02/` directory, and the challenge `README.md`. The commit message should specify: `c02: Newick/Nexus round-trip preserves topology with 1.2e-6 max branch-length drift, support values preserved`.
