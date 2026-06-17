"""
Exercise 2 - Compute a distance matrix and build a neighbor-joining tree.

Goal: take the trimmed alignment from Exercise 1, compute pairwise
distances under both Jukes-Cantor (Jukes and Cantor 1969) and Kimura
two-parameter (Kimura 1980), build NJ trees with Biopython, root them on
a named outgroup, and write Newick + Nexus files.

The exercise covers:

- Hand-rolling the p-distance, JC distance, and K2P distance formulas.
- Building a Biopython DistanceMatrix from a list of (name, sequence).
- Calling Bio.Phylo.TreeConstruction.DistanceTreeConstructor.nj.
- Rooting on a named outgroup.
- Writing both Newick and Nexus.
- Pinning the distance method and the outgroup in the run-info JSON.

Estimated time: 75 minutes (30 min reading the distance-model derivation,
30 min implementing the matrix builder, 15 min running and inspecting).

Acceptance criteria:
- `python exercise-02-distance-and-nj-tree.py --input
    results/ex01/trimmed.fasta --out-dir results/ex02
    --outgroup Ornithorhynchus_anatinus` runs end to end.
- `results/ex02/distance_matrix_jc.tsv` exists, square, with the right
  taxon labels along both axes and zeros on the diagonal.
- `results/ex02/distance_matrix_k2p.tsv` similarly.
- `results/ex02/tree_jc_nj.nwk` and `results/ex02/tree_k2p_nj.nwk` exist
  and contain valid Newick.
- `results/ex02/tree_k2p_nj.nex` exists and contains valid Nexus.
- `results/ex02/run-info.json` records the distance methods, the
  outgroup, and the Biopython version.
- A second run with the same inputs produces byte-identical Newick.

Requirements:
    pip install biopython
    (no other dependencies; pure Python plus Biopython.)

What you learn:
- The closed-form JC and K2P distance formulas, with saturation handling.
- The Biopython DistanceMatrix data shape (lower-triangle list-of-lists).
- The neighbor-joining algorithm via DistanceTreeConstructor.
- The Newick + Nexus write idiom.
- Why JC and K2P diverge as p-distance grows.

Tool versions assumed:
- Python 3.11+
- Biopython 1.84+

References:
- JC: Jukes and Cantor 1969 (Mammalian Protein Metabolism, vol III).
- K2P: Kimura 1980, J Mol Evol 16:111.
- NJ: Saitou and Nei 1987, Mol Biol Evol 4:406
  https://academic.oup.com/mbe/article/4/4/406/1029664
- Biopython tutorial chapter 13:
  https://biopython.org/docs/latest/Tutorial/index.html
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Distance models.
# ----------------------------------------------------------------------

SATURATION_CLAMP: float = 0.7499


def p_distance(seq_a: str, seq_b: str) -> float:
    """Return the p-distance between two aligned sequences of equal length.

    Pairwise-deletion convention: skip columns where either sequence has
    a gap or an N. Returns 0.0 if no informative columns remain.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("p_distance requires equal-length aligned sequences.")
    differ: int = 0
    total: int = 0
    for ch_a, ch_b in zip(seq_a.upper(), seq_b.upper()):
        if ch_a in {"-", "N"} or ch_b in {"-", "N"}:
            continue
        total += 1
        if ch_a != ch_b:
            differ += 1
    if total == 0:
        return 0.0
    return differ / total


def jukes_cantor_distance_from_p(p: float) -> float:
    """Return the Jukes-Cantor corrected distance for an observed p-distance.

    Clamps p at SATURATION_CLAMP. Returns float('inf') if the input is
    already saturated.
    """
    if p < 0.0:
        raise ValueError(f"p-distance cannot be negative; got {p}.")
    if p >= SATURATION_CLAMP:
        return float("inf")
    return -0.75 * math.log(1.0 - 4.0 / 3.0 * p)


def jukes_cantor_distance(seq_a: str, seq_b: str) -> float:
    """Return the Jukes-Cantor corrected distance between two aligned
    nucleotide sequences.
    """
    return jukes_cantor_distance_from_p(p_distance(seq_a, seq_b))


def kimura_2p_distance(seq_a: str, seq_b: str) -> float:
    """Return the Kimura two-parameter corrected distance between two aligned
    nucleotide sequences.

    Returns float('inf') if either log argument is non-positive (the
    saturation regime).
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("kimura_2p_distance requires equal-length aligned sequences.")
    purines: set[str] = {"A", "G"}
    transitions: int = 0
    transversions: int = 0
    total: int = 0
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


# ----------------------------------------------------------------------
# FASTA loading.
# ----------------------------------------------------------------------

def load_aligned_records(fasta: Path) -> list[tuple[str, str]]:
    """Return a list of (taxon_name, aligned_sequence) tuples.

    Verifies that all sequences are the same length.
    """
    from Bio import SeqIO  # lazy import for py_compile safety

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


# ----------------------------------------------------------------------
# Distance matrix.
# ----------------------------------------------------------------------

def build_distance_matrix(
    records: list[tuple[str, str]],
    method: str,
) -> "Any":
    """Build a Biopython DistanceMatrix under method 'jc' or 'k2p'.

    Returns a Bio.Phylo.TreeConstruction.DistanceMatrix.
    """
    from Bio.Phylo.TreeConstruction import DistanceMatrix  # lazy import

    names: list[str] = [name for name, _ in records]
    n: int = len(records)
    matrix_rows: list[list[float]] = []
    for i in range(n):
        row: list[float] = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
                continue
            seq_a: str = records[i][1]
            seq_b: str = records[j][1]
            if method == "jc":
                value: float = jukes_cantor_distance(seq_a, seq_b)
            elif method == "k2p":
                value = kimura_2p_distance(seq_a, seq_b)
            else:
                raise ValueError(f"Unknown distance method: {method!r}")
            if math.isinf(value):
                # Saturate at a finite value so NJ does not blow up.
                value = 10.0
            row.append(value)
        matrix_rows.append(row)
    return DistanceMatrix(names=names, matrix=matrix_rows)


def write_distance_tsv(
    matrix: "Any",
    out_path: Path,
) -> None:
    """Write a square TSV (lower-triangle expanded to full matrix)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    names: list[str] = matrix.names
    with out_path.open("w") as fh:
        fh.write("taxon\t" + "\t".join(names) + "\n")
        for i, name_i in enumerate(names):
            row_vals: list[str] = []
            for j, _ in enumerate(names):
                value: float = matrix[name_i, names[j]]
                row_vals.append(f"{value:.6f}")
            fh.write(name_i + "\t" + "\t".join(row_vals) + "\n")


# ----------------------------------------------------------------------
# Tree building.
# ----------------------------------------------------------------------

def build_nj_tree(matrix: "Any") -> "Any":
    """Build an NJ tree from a Biopython DistanceMatrix."""
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor  # lazy import

    constructor = DistanceTreeConstructor()
    return constructor.nj(matrix)


def root_tree(tree: "Any", outgroup: str) -> "Any":
    """Root the tree on a named outgroup taxon."""
    tree.root_with_outgroup(outgroup)
    # Ladderize for a deterministic leaf order in rendered figures.
    tree.ladderize()
    return tree


def write_newick(tree: "Any", out_path: Path) -> None:
    """Write a Biopython Tree to Newick."""
    from Bio import Phylo  # lazy import

    out_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out_path), "newick")


def write_nexus(tree: "Any", out_path: Path) -> None:
    """Write a Biopython Tree to Nexus."""
    from Bio import Phylo  # lazy import

    out_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out_path), "nexus")


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class Exercise2RunInfo:
    run_date: str = ""
    input_aligned_fasta: str = ""
    outgroup: str = ""
    distance_methods: str = "jc,k2p"
    tree_builder: str = "biopython_nj"
    n_input_records: int = 0
    biopython_version: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def biopython_version() -> str:
    try:
        import Bio  # type: ignore[import-not-found]
        return getattr(Bio, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    input_fasta: Path,
    out_dir: Path,
    outgroup: str,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    records: list[tuple[str, str]] = load_aligned_records(input_fasta)
    if outgroup not in {name for name, _ in records}:
        raise ValueError(
            f"Outgroup taxon {outgroup!r} not found in input. "
            f"Available: {sorted(name for name, _ in records)}"
        )

    matrix_jc = build_distance_matrix(records, method="jc")
    matrix_k2p = build_distance_matrix(records, method="k2p")

    write_distance_tsv(matrix_jc, out_dir / "distance_matrix_jc.tsv")
    write_distance_tsv(matrix_k2p, out_dir / "distance_matrix_k2p.tsv")

    tree_jc = build_nj_tree(matrix_jc)
    tree_k2p = build_nj_tree(matrix_k2p)
    root_tree(tree_jc, outgroup)
    root_tree(tree_k2p, outgroup)

    write_newick(tree_jc, out_dir / "tree_jc_nj.nwk")
    write_newick(tree_k2p, out_dir / "tree_k2p_nj.nwk")
    write_nexus(tree_k2p, out_dir / "tree_k2p_nj.nex")

    info = Exercise2RunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        input_aligned_fasta=str(input_fasta),
        outgroup=outgroup,
        distance_methods="jc,k2p",
        tree_builder="biopython_nj",
        n_input_records=len(records),
        biopython_version=biopython_version(),
        python_version=python_version_string(),
    )
    run_info_path: Path = out_dir / "run-info.json"
    run_info_path.parent.mkdir(parents=True, exist_ok=True)
    with run_info_path.open("w") as fh:
        json.dump(info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 2 - JC and K2P distances; NJ tree; Newick + Nexus.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input aligned FASTA (e.g. results/ex01/trimmed.fasta).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--outgroup",
        type=str,
        required=True,
        help="Outgroup taxon ID to root the tree on (must match a record ID).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        input_fasta=args.input,
        out_dir=args.out_dir,
        outgroup=args.outgroup,
    )
    print(f"[ex02] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
