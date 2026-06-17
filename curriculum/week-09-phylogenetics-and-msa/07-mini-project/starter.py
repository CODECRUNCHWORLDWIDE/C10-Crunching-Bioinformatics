"""
Mini-project starter - End-to-end phylogenetics pipeline.

This pipeline takes a FASTA of homologous sequences and emits a phylogenetic
tree with bootstrap support, plus a run-info JSON that records every
parameter needed to reproduce the run. The output looks like a publishable
tree figure; it is reproducible only if the run-info travels with it.

This file is a SKELETON. It compiles cleanly under `python3 -m py_compile`
but most functions raise NotImplementedError. Your job is to fill in the
TODOs to produce a working `build_tree(fasta_path, out_dir, outgroup, seed)`
function that runs the seven-stage pipeline described in README.md:

  1. validate_input(fasta_path)
  2. run_mafft(fasta_path, out_aligned)
  3. trim_alignment(aligned, trimmed, max_gap_fraction=0.5)
  4. build_distance_matrix_k2p(records)
  5. build_and_root_nj_tree(matrix, outgroup)
  6. bootstrap_support(main_tree, records, n_replicates, seed)
  7. write_outputs(main_tree, out_dir, run_info)

Each function has a docstring, a type signature, and a NotImplementedError
or a stub return. Replace the body with your implementation. The acceptance
criteria in README.md tell you what each function must produce.

Tool versions assumed:
- Python 3.11+
- MAFFT 7.526 (CLI tool; subprocess.run it)
- Biopython 1.84+
- matplotlib 3.8+
- IQ-TREE 2.3.6 (optional; only if you implement the stretch goal)
- ete3 3.1.3 (optional; for publication-quality rendering)
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import random
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class PipelineRunInfo:
    """Provenance metadata recorded on every run.

    Every field must be filled in before the pipeline writes the
    run-info JSON. Empty values are silently wrong; assert non-empty
    before writing.
    """
    run_date: str = ""
    input_fasta: str = ""
    input_md5: str = ""
    n_input_records: int = 0
    mafft_version: str = ""
    mafft_algorithm: str = "--retree 2 --maxiterate 0 --nuc --anysymbol --quiet"
    trim_threshold: float = 0.5
    n_aligned_columns_before_trim: int = 0
    n_aligned_columns_after_trim: int = 0
    distance_method: str = "k2p"
    tree_builder: str = "biopython_nj"
    bootstrap_replicates: int = 500
    bootstrap_method: str = "column_resample_nj"
    seed: int = 42
    outgroup: str = ""
    iqtree_version: str = ""
    iqtree_model: str = ""
    biopython_version: str = ""
    ete3_version: str = ""
    python_version: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Stage 1 - input validation.
# ----------------------------------------------------------------------

def validate_input(fasta_path: Path) -> int:
    """Confirm the input FASTA exists, has at least 5 records, no duplicate
    record IDs, and no all-N sequences. Return the record count.

    Raises FileNotFoundError or ValueError on any problem.
    """
    # TODO: implement.
    # Reference: see exercise-01-mafft-via-subprocess.py validate_input_fasta.
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")
    raise NotImplementedError("Implement validate_input.")


# ----------------------------------------------------------------------
# Stage 2 - run MAFFT.
# ----------------------------------------------------------------------

def run_mafft(fasta_path: Path, out_aligned: Path, threads: int = 4) -> None:
    """Run MAFFT in FFT-NS-2 mode. Writes to out_aligned.

    Pinned flags: --retree 2 --maxiterate 0 --nuc --anysymbol --quiet.
    """
    # TODO: implement.
    # Reference: exercise-01 run_mafft.
    raise NotImplementedError("Implement run_mafft.")


# ----------------------------------------------------------------------
# Stage 3 - trim gappy columns.
# ----------------------------------------------------------------------

def trim_alignment(
    aligned: Path,
    trimmed: Path,
    max_gap_fraction: float = 0.5,
) -> tuple[int, int]:
    """Drop columns where the gap fraction exceeds max_gap_fraction.

    Returns (n_cols_before, n_cols_after).
    """
    # TODO: implement.
    # Reference: exercise-01 trim_alignment.
    raise NotImplementedError("Implement trim_alignment.")


# ----------------------------------------------------------------------
# Stage 4 - distance matrix.
# ----------------------------------------------------------------------

def load_aligned_records(fasta: Path) -> list[tuple[str, str]]:
    """Read the trimmed alignment as a list of (taxon_name, seq) tuples.

    Verifies all sequences have the same length.
    """
    # TODO: implement.
    # Reference: exercise-02 load_aligned_records.
    raise NotImplementedError("Implement load_aligned_records.")


def build_distance_matrix_k2p(records: list[tuple[str, str]]) -> "Any":
    """Build a Biopython DistanceMatrix under the K2P model.

    Returns a Bio.Phylo.TreeConstruction.DistanceMatrix.
    """
    # TODO: implement.
    # Reference: exercise-02 build_distance_matrix (method='k2p').
    raise NotImplementedError("Implement build_distance_matrix_k2p.")


def write_distance_tsv(matrix: "Any", out_path: Path) -> None:
    """Write a square TSV of the K2P pairwise distances."""
    # TODO: implement.
    # Reference: exercise-02 write_distance_tsv.
    raise NotImplementedError("Implement write_distance_tsv.")


# ----------------------------------------------------------------------
# Stage 5 - NJ tree.
# ----------------------------------------------------------------------

def build_and_root_nj_tree(matrix: "Any", outgroup: str) -> "Any":
    """Build a Biopython NJ tree, root on outgroup, ladderize."""
    # TODO: implement.
    # Reference: exercise-02 build_nj_tree + root_tree.
    raise NotImplementedError("Implement build_and_root_nj_tree.")


# ----------------------------------------------------------------------
# Stage 6 - bootstrap.
# ----------------------------------------------------------------------

def bootstrap_support(
    main_tree: "Any",
    records: list[tuple[str, str]],
    outgroup: str,
    n_replicates: int,
    seed: int,
) -> None:
    """Resample columns, rebuild NJ trees, compute per-branch bootstrap
    support, attach to main_tree internal nodes as clade.confidence.
    """
    # TODO: implement.
    # Reference: exercise-03 resample_columns + get_bipartitions + compute_support.
    raise NotImplementedError("Implement bootstrap_support.")


# ----------------------------------------------------------------------
# Stage 7 - render and write.
# ----------------------------------------------------------------------

def render_png(tree: "Any", out_path: Path) -> None:
    """Render the tree as a PNG via Bio.Phylo + matplotlib (Agg)."""
    # TODO: implement.
    # Reference: exercise-03 render_png_biopython.
    raise NotImplementedError("Implement render_png.")


def write_tree_files(tree: "Any", out_newick: Path, out_nexus: Path) -> None:
    """Write the tree to Newick and Nexus."""
    # TODO: implement.
    raise NotImplementedError("Implement write_tree_files.")


def write_run_info(info: PipelineRunInfo, out_path: Path) -> None:
    """Write the run-info JSON. Raise on empty required fields."""
    if not info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    if not info.mafft_version:
        raise ValueError("run_info.mafft_version is empty; refusing to write.")
    if not info.outgroup:
        raise ValueError("run_info.outgroup is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

def md5_of_file(path: Path) -> str:
    """Return the hex MD5 of a file's bytes."""
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def get_mafft_version() -> str:
    """Return the MAFFT version string (from `mafft --version`)."""
    result = subprocess.run(
        ["mafft", "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    return (result.stderr or result.stdout).strip().splitlines()[0]


def biopython_version() -> str:
    try:
        import Bio  # type: ignore[import-not-found]
        return getattr(Bio, "__version__", "unknown")
    except Exception:
        return "unknown"


def ete3_version() -> str:
    try:
        import ete3  # type: ignore[import-not-found]
        return getattr(ete3, "__version__", "unknown")
    except Exception:
        return "unavailable"


def python_version_string() -> str:
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def build_tree(
    fasta_path: Path,
    out_dir: Path,
    outgroup: str,
    seed: int = 42,
    n_replicates: int = 500,
) -> Path:
    """Run the full pipeline. Returns the path to the final Newick tree.

    Orchestrates the seven stages from README.md. Each stage delegates to
    one of the functions above. The reference implementation runs in
    under three minutes on the demo COI panel.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned: Path = out_dir / "aligned.fasta"
    trimmed: Path = out_dir / "trimmed.fasta"
    matrix_tsv: Path = out_dir / "distance_matrix_k2p.tsv"
    newick: Path = out_dir / "tree_nj_bootstrap.nwk"
    nexus: Path = out_dir / "tree_final.nex"
    png: Path = out_dir / "tree_final.png"
    run_info_path: Path = out_dir / "run-info.json"

    # Stage 1.
    n_input_records: int = validate_input(fasta_path)

    # Stage 2.
    run_mafft(fasta_path, aligned)

    # Stage 3.
    n_cols_before, n_cols_after = trim_alignment(aligned, trimmed)

    # Stage 4.
    records: list[tuple[str, str]] = load_aligned_records(trimmed)
    matrix = build_distance_matrix_k2p(records)
    write_distance_tsv(matrix, matrix_tsv)

    # Stage 5.
    tree = build_and_root_nj_tree(matrix, outgroup)

    # Stage 6.
    bootstrap_support(tree, records, outgroup, n_replicates=n_replicates, seed=seed)

    # Stage 7.
    write_tree_files(tree, newick, nexus)
    render_png(tree, png)

    info = PipelineRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        input_fasta=str(fasta_path),
        input_md5=md5_of_file(fasta_path),
        n_input_records=n_input_records,
        mafft_version=get_mafft_version(),
        n_aligned_columns_before_trim=n_cols_before,
        n_aligned_columns_after_trim=n_cols_after,
        bootstrap_replicates=n_replicates,
        seed=seed,
        outgroup=outgroup,
        biopython_version=biopython_version(),
        ete3_version=ete3_version(),
        python_version=python_version_string(),
        notes="Demo COI metazoan panel; NJ tree under K2P with column-resample bootstrap.",
    )
    write_run_info(info, run_info_path)
    return newick


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="End-to-end phylogenetics pipeline (MAFFT + K2P + NJ + bootstrap).",
    )
    parser.add_argument("--input", type=Path, required=True, help="Input FASTA.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--outgroup",
        type=str,
        required=True,
        help="Outgroup taxon ID for rooting.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--replicates",
        type=int,
        default=500,
        help="Bootstrap replicate count.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    out_newick: Path = build_tree(
        fasta_path=args.input,
        out_dir=args.out_dir,
        outgroup=args.outgroup,
        seed=args.seed,
        n_replicates=args.replicates,
    )
    print(f"[mini-project] tree at {out_newick}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
