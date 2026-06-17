"""
Exercise 3 - Bootstrap a neighbor-joining tree and render it.

Goal: take the trimmed alignment from Exercise 1, generate B = 100
bootstrap-resampled alignments by column resampling with replacement,
rebuild an NJ tree per replicate under K2P, compute per-branch support
values by matching internal bipartitions across replicates, render the
final tree as ASCII to stdout, as a PNG via Bio.Phylo + matplotlib, and
(optionally) as a PNG via ete3 with branches coloured by support tier.

The exercise covers:

- Column-resample bootstrap (Felsenstein 1985).
- Bipartition matching for branch support.
- ASCII rendering (Bio.Phylo.draw_ascii).
- PNG rendering via Bio.Phylo.draw + matplotlib (Agg backend).
- Lazy ete3 import for the optional publication-quality render.

Estimated time: 90 minutes (45 min reading the bootstrap and bipartition
logic, 30 min implementing, 15 min running and inspecting).

Acceptance criteria:
- `python exercise-03-render-tree-with-bootstrap.py --input
    results/ex01/trimmed.fasta --out-dir results/ex03
    --outgroup Ornithorhynchus_anatinus --replicates 100 --seed 42`
  runs end to end without errors.
- `results/ex03/tree_k2p_nj_with_support.nwk` exists; internal node
  names are the bootstrap percentages (0-100).
- `results/ex03/tree_k2p_nj_with_support.png` exists and renders the
  tree with bootstrap labels.
- `results/ex03/tree_ascii.txt` exists and contains the ASCII tree.
- `results/ex03/run-info.json` exists with the seed, the replicate count,
  the distance method, and the Biopython version.
- A second run with the same seed produces identical bootstrap percentages.

Requirements:
    pip install biopython matplotlib
    Optional: pip install ete3 (the ete3 render is skipped if unavailable).

What you learn:
- The column-resample bootstrap pattern.
- Bipartition-based branch support.
- Headless matplotlib (matplotlib.use('Agg')).
- The lazy-import pattern for an optional dependency.

Tool versions assumed:
- Python 3.11+
- Biopython 1.84+
- matplotlib 3.8+
- ete3 3.1.3 (optional)

References:
- Bootstrap: Felsenstein 1985, Evolution 39:783
  https://www.jstor.org/stable/2408678
- ete3: Huerta-Cepas et al. 2016, Mol Biol Evol 33:1635
  https://academic.oup.com/mbe/article/33/6/1635/2579822
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Bootstrap.
# ----------------------------------------------------------------------

def resample_columns(
    records: list[tuple[str, str]],
    rng: random.Random,
) -> list[tuple[str, str]]:
    """Resample alignment columns with replacement.

    Each output record has the same length as the input. The set of
    columns is a random sample (with replacement) of the original column
    indices.
    """
    if not records:
        return []
    n_cols: int = len(records[0][1])
    indices: list[int] = [rng.randrange(n_cols) for _ in range(n_cols)]
    return [
        (name, "".join(seq[k] for k in indices))
        for name, seq in records
    ]


def build_one_bootstrap_tree(
    records: list[tuple[str, str]],
    distance_method: str,
    rng: random.Random,
) -> "Any":
    """Generate one bootstrap-resampled NJ tree."""
    # Local imports to keep this file py_compile-clean even without Biopython.
    from exercise_02_distance_and_nj_tree import (  # noqa: F401 - imported for reuse
        build_distance_matrix,
        build_nj_tree,
    )

    resampled = resample_columns(records, rng)
    matrix = build_distance_matrix(resampled, method=distance_method)
    return build_nj_tree(matrix)


# ----------------------------------------------------------------------
# Bipartition matching for branch support.
# ----------------------------------------------------------------------

def get_bipartitions(tree: "Any", all_taxa: list[str]) -> set[frozenset[str]]:
    """Return the set of non-trivial bipartitions in a tree.

    A bipartition is the set of leaf names below an internal node; the
    "non-trivial" filter drops the bipartitions that are a single leaf or
    the full taxon set.
    """
    bipartitions: set[frozenset[str]] = set()
    all_set: frozenset[str] = frozenset(all_taxa)
    for clade in tree.get_nonterminals():
        leaves_below: frozenset[str] = frozenset(
            t.name for t in clade.get_terminals()
        )
        if 1 < len(leaves_below) < len(all_set):
            bipartitions.add(leaves_below)
    return bipartitions


def compute_support(
    main_tree: "Any",
    bootstrap_trees: list["Any"],
    all_taxa: list[str],
) -> None:
    """Annotate each internal node of main_tree with bootstrap percentage.

    Modifies main_tree in place. The annotation lives on
    clade.confidence as a float in [0, 100].
    """
    n_replicates: int = len(bootstrap_trees)
    if n_replicates == 0:
        return
    rep_bipartitions: list[set[frozenset[str]]] = [
        get_bipartitions(t, all_taxa) for t in bootstrap_trees
    ]
    for clade in main_tree.get_nonterminals():
        leaves_below: frozenset[str] = frozenset(
            t.name for t in clade.get_terminals()
        )
        if 1 < len(leaves_below) < len(all_taxa):
            count: int = sum(
                1 for bset in rep_bipartitions if leaves_below in bset
            )
            clade.confidence = round(100.0 * count / n_replicates, 1)


# ----------------------------------------------------------------------
# Rendering.
# ----------------------------------------------------------------------

def render_ascii(tree: "Any", out_path: Path) -> None:
    """Render the tree as ASCII to a text file."""
    from Bio import Phylo  # lazy import

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        Phylo.draw_ascii(tree, file=fh)


def render_png_biopython(tree: "Any", out_path: Path) -> None:
    """Render the tree as a PNG via Bio.Phylo.draw + matplotlib (Agg)."""
    import matplotlib
    matplotlib.use("Agg")  # headless backend; must be set before pyplot.
    import matplotlib.pyplot as plt
    from Bio import Phylo  # lazy import

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    Phylo.draw(tree, axes=ax, do_show=False)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)


def render_png_ete3(newick_path: Path, out_path: Path) -> bool:
    """Render via ete3 with branches coloured by support tier.

    Returns True on success, False if ete3 is unavailable. Lazy import so
    this file compiles on machines without ete3.
    """
    try:
        import os
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from ete3 import Tree, TreeStyle, NodeStyle  # type: ignore[import-not-found]
    except Exception:
        print("[ex03] ete3 not available; skipping ete3 render.", file=sys.stderr)
        return False

    tree = Tree(str(newick_path), format=1)
    for node in tree.traverse():
        if node.is_leaf():
            continue
        support: float = float(getattr(node, "support", 0.0) or 0.0)
        style = NodeStyle()
        if support >= 90:
            style["hz_line_color"] = "#1a7f37"
        elif support >= 70:
            style["hz_line_color"] = "#0969da"
        else:
            style["hz_line_color"] = "#cf222e"
        style["hz_line_width"] = 2
        node.set_style(style)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.render(str(out_path), w=900, tree_style=ts)
    return True


# ----------------------------------------------------------------------
# Loaders and helpers (kept self-contained so this file is independently
# runnable, but uses the Exercise 2 helpers when available).
# ----------------------------------------------------------------------

def load_aligned_records(fasta: Path) -> list[tuple[str, str]]:
    """Mirror of the Exercise 2 helper. Local to keep this file
    independently runnable for graders.
    """
    from Bio import SeqIO  # lazy import

    records: list[tuple[str, str]] = []
    for rec in SeqIO.parse(str(fasta), "fasta"):
        records.append((rec.id, str(rec.seq)))
    if not records:
        raise ValueError(f"No records in {fasta}.")
    lengths: set[int] = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(
            f"Aligned FASTA records have inconsistent lengths: {sorted(lengths)}"
        )
    return records


def biopython_version() -> str:
    try:
        import Bio  # type: ignore[import-not-found]
        return getattr(Bio, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class Exercise3RunInfo:
    run_date: str = ""
    input_aligned_fasta: str = ""
    outgroup: str = ""
    distance_method: str = "k2p"
    n_replicates: int = 0
    seed: int = 42
    tree_builder: str = "biopython_nj"
    biopython_version: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    input_fasta: Path,
    out_dir: Path,
    outgroup: str,
    n_replicates: int,
    seed: int,
    distance_method: str = "k2p",
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    # Imports of the local Exercise 2 helpers; we do these here (not at
    # module top) so this file compiles whether or not exercise 2 is on
    # the PYTHONPATH.
    from importlib import import_module

    try:
        ex02 = import_module("exercise_02_distance_and_nj_tree")
        build_distance_matrix = ex02.build_distance_matrix
        build_nj_tree = ex02.build_nj_tree
    except Exception:
        # Inline fallback. The trade-off is a small amount of duplicated
        # code in exchange for the file being independently runnable.
        from Bio.Phylo.TreeConstruction import (  # type: ignore[import-not-found]
            DistanceMatrix,
            DistanceTreeConstructor,
        )

        def build_distance_matrix(records_, method):  # type: ignore[no-redef]
            from math import log as _log

            def _p(a, b):
                d = 0
                t = 0
                for ca, cb in zip(a.upper(), b.upper()):
                    if ca in {"-", "N"} or cb in {"-", "N"}:
                        continue
                    t += 1
                    if ca != cb:
                        d += 1
                return 0.0 if t == 0 else d / t

            def _jc(p):
                if p >= 0.7499:
                    return 10.0
                return -0.75 * _log(1.0 - 4.0 / 3.0 * p)

            def _k2p(a, b):
                purines = {"A", "G"}
                trs = 0
                trv = 0
                t = 0
                for ca, cb in zip(a.upper(), b.upper()):
                    if ca not in "ACGT" or cb not in "ACGT":
                        continue
                    t += 1
                    if ca == cb:
                        continue
                    if (ca in purines) == (cb in purines):
                        trs += 1
                    else:
                        trv += 1
                if t == 0:
                    return 0.0
                bp = trs / t
                bq = trv / t
                arg1 = 1.0 - 2.0 * bp - bq
                arg2 = 1.0 - 2.0 * bq
                if arg1 <= 0.0 or arg2 <= 0.0:
                    return 10.0
                return -0.5 * _log(arg1) - 0.25 * _log(arg2)

            names_ = [n for n, _ in records_]
            rows_ = []
            for i in range(len(records_)):
                r = []
                for j in range(i + 1):
                    if i == j:
                        r.append(0.0)
                    else:
                        a = records_[i][1]
                        b = records_[j][1]
                        if method == "jc":
                            r.append(_jc(_p(a, b)))
                        elif method == "k2p":
                            r.append(_k2p(a, b))
                        else:
                            raise ValueError(method)
                rows_.append(r)
            return DistanceMatrix(names=names_, matrix=rows_)

        def build_nj_tree(matrix):  # type: ignore[no-redef]
            return DistanceTreeConstructor().nj(matrix)

    records: list[tuple[str, str]] = load_aligned_records(input_fasta)
    all_taxa: list[str] = [name for name, _ in records]
    if outgroup not in all_taxa:
        raise ValueError(f"Outgroup {outgroup!r} not in input record IDs.")

    main_matrix = build_distance_matrix(records, method=distance_method)
    main_tree = build_nj_tree(main_matrix)
    main_tree.root_with_outgroup(outgroup)
    main_tree.ladderize()

    rng: random.Random = random.Random(seed)
    bootstrap_trees: list[Any] = []
    for _ in range(n_replicates):
        resampled = resample_columns(records, rng)
        matrix_b = build_distance_matrix(resampled, method=distance_method)
        tree_b = build_nj_tree(matrix_b)
        # Root each replicate the same way for bipartition stability.
        try:
            tree_b.root_with_outgroup(outgroup)
        except Exception:
            pass
        bootstrap_trees.append(tree_b)

    compute_support(main_tree, bootstrap_trees, all_taxa)

    # Write Newick + ASCII + PNG.
    newick_path: Path = out_dir / "tree_k2p_nj_with_support.nwk"
    ascii_path: Path = out_dir / "tree_ascii.txt"
    png_path: Path = out_dir / "tree_k2p_nj_with_support.png"
    ete3_png_path: Path = out_dir / "tree_k2p_nj_with_support_ete3.png"

    out_dir.mkdir(parents=True, exist_ok=True)
    from Bio import Phylo  # lazy import
    Phylo.write(main_tree, str(newick_path), "newick")
    render_ascii(main_tree, ascii_path)
    render_png_biopython(main_tree, png_path)
    render_png_ete3(newick_path, ete3_png_path)

    info = Exercise3RunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        input_aligned_fasta=str(input_fasta),
        outgroup=outgroup,
        distance_method=distance_method,
        n_replicates=n_replicates,
        seed=seed,
        tree_builder="biopython_nj",
        biopython_version=biopython_version(),
        python_version=python_version_string(),
    )
    run_info_path: Path = out_dir / "run-info.json"
    with run_info_path.open("w") as fh:
        json.dump(info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 3 - Column-resample bootstrap NJ + rendering.",
    )
    parser.add_argument("--input", type=Path, required=True, help="Trimmed aligned FASTA.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--outgroup", type=str, required=True, help="Outgroup taxon ID.")
    parser.add_argument("--replicates", type=int, default=100, help="Bootstrap replicates.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--distance-method",
        type=str,
        default="k2p",
        choices=["jc", "k2p"],
        help="Distance model.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        input_fasta=args.input,
        out_dir=args.out_dir,
        outgroup=args.outgroup,
        n_replicates=args.replicates,
        seed=args.seed,
        distance_method=args.distance_method,
    )
    print(f"[ex03] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
