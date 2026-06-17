"""
Exercise 1 - Run MAFFT on a small FASTA via subprocess and inspect the alignment.

Goal: take a 10-sequence vertebrate cytochrome b FASTA, run MAFFT
(--retree 2 --maxiterate 0 --nuc), parse the alignment with Biopython,
report column statistics (count, per-column gap fraction, conservation),
trim columns with > 50% gap content, and save the trimmed alignment.

The exercise covers:

- Calling a CLI tool from Python with subprocess.run, check=True,
  capture_output=True.
- Parsing the aligned FASTA with Bio.AlignIO.
- Iterating alignment columns and computing per-column statistics.
- Pinning the MAFFT algorithm flag for reproducibility.
- Writing a run-info JSON alongside the alignment.

Estimated time: 60 minutes (20 min reading, 30 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-01-mafft-via-subprocess.py --input data/cytb_vertebrates.fasta
    --out-dir results/ex01` runs end to end without errors.
- `results/ex01/aligned.fasta` exists with all 10 input records, every
  record having the same length.
- `results/ex01/trimmed.fasta` exists with the same row count and a
  smaller column count.
- `results/ex01/alignment_summary.tsv` exists with at least the columns:
  column_index, gap_fraction, consensus_residue, conservation_fraction.
- `results/ex01/run-info.json` exists with the MAFFT version, the algorithm
  flag, the trim threshold, and the run date.
- A second run with the same input produces a byte-identical aligned.fasta.

Requirements:
    conda install -c bioconda mafft=7.526 biopython=1.84
    OR pip install biopython (and have MAFFT on the PATH).

What you learn:
- The subprocess.run idiom with stdout redirected to a file.
- The Bio.AlignIO API for multiple sequence alignments.
- The column-iteration pattern: alignment[:, j] is the j-th column as a string.
- The "pin the algorithm, write the run-info" reproducibility pattern.

Tool versions assumed:
- Python 3.11+
- MAFFT 7.526 (CLI tool; subprocess.run it)
- Biopython 1.84+

References:
- MAFFT: Katoh and Standley 2013, Mol Biol Evol 30:772
  https://academic.oup.com/mbe/article/30/4/772/1073398
- Biopython tutorial:
  https://biopython.org/docs/latest/Tutorial/index.html
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_MAFFT_ALGORITHM_FLAGS: list[str] = [
    "--retree", "2",
    "--maxiterate", "0",
    "--nuc",
    "--anysymbol",
    "--quiet",
]
DEFAULT_TRIM_THRESHOLD: float = 0.5


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class ExerciseRunInfo:
    """Provenance metadata for the Exercise 1 run."""
    run_date: str = ""
    input_fasta: str = ""
    input_md5: str = ""
    mafft_version: str = ""
    mafft_algorithm: str = ""
    trim_threshold: float = DEFAULT_TRIM_THRESHOLD
    n_input_records: int = 0
    n_aligned_columns_before_trim: int = 0
    n_aligned_columns_after_trim: int = 0
    biopython_version: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


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
    """Return the MAFFT version string, e.g. 'v7.526'.

    Raises subprocess.CalledProcessError if MAFFT is not on the PATH.
    """
    result = subprocess.run(
        ["mafft", "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    # MAFFT writes the version to stderr, not stdout, by convention.
    version: str = (result.stderr or result.stdout).strip().splitlines()[0]
    return version


def biopython_version() -> str:
    """Return the installed Biopython version, or 'unknown' if unimportable."""
    try:
        import Bio  # type: ignore[import-not-found]
        return getattr(Bio, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    """Return the running Python version, e.g. '3.11.9'."""
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Stage 1 - validate input.
# ----------------------------------------------------------------------

def validate_input_fasta(input_fasta: Path) -> int:
    """Confirm the input FASTA exists, is non-empty, and has at least
    two unique records. Return the record count.

    Raises FileNotFoundError or ValueError on any problem.
    """
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
    from Bio import SeqIO  # lazy import for py_compile safety

    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    if len(records) < 2:
        raise ValueError(
            f"Input FASTA needs at least two records; got {len(records)} in {input_fasta}."
        )
    names: list[str] = [r.id for r in records]
    if len(set(names)) != len(names):
        raise ValueError(
            f"Input FASTA has duplicate record IDs; dedupe before aligning."
        )
    return len(records)


# ----------------------------------------------------------------------
# Stage 2 - run MAFFT.
# ----------------------------------------------------------------------

def run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    algorithm_flags: list[str] | None = None,
    threads: int = 4,
) -> None:
    """Run MAFFT with the pinned algorithm flags. Writes to output_fasta.

    The default algorithm is FFT-NS-2 (--retree 2 --maxiterate 0), which
    is deterministic and fast. For higher accuracy on small inputs, pass
    ['--localpair', '--maxiterate', '1000', '--nuc', '--anysymbol',
    '--quiet'] (L-INS-i mode).
    """
    if algorithm_flags is None:
        algorithm_flags = list(DEFAULT_MAFFT_ALGORITHM_FLAGS)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd: list[str] = ["mafft"] + algorithm_flags + ["--thread", str(threads), str(input_fasta)]
    with output_fasta.open("w") as fh:
        result = subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
    # MAFFT prints progress to stderr; surface anything unusual.
    if result.stderr:
        for line in result.stderr.splitlines():
            if line.strip() and "done" not in line.lower():
                print(f"[mafft] {line.strip()}", file=sys.stderr)


# ----------------------------------------------------------------------
# Stage 3 - parse and summarize.
# ----------------------------------------------------------------------

@dataclass
class ColumnSummary:
    column_index: int
    gap_fraction: float
    consensus_residue: str
    conservation_fraction: float


def column_consensus(column: str) -> tuple[str, float]:
    """Return the most-common non-gap residue and its frequency."""
    non_gap: str = column.replace("-", "").upper()
    if not non_gap:
        return ("-", 0.0)
    counts: dict[str, int] = {}
    for ch in non_gap:
        counts[ch] = counts.get(ch, 0) + 1
    best: str = max(counts, key=lambda k: counts[k])
    return (best, counts[best] / len(non_gap))


def summarize_alignment(aligned_fasta: Path) -> tuple[int, int, list[ColumnSummary]]:
    """Read the alignment and return (n_rows, n_cols, per_column_summary)."""
    from Bio import AlignIO  # lazy import for py_compile safety

    alignment = AlignIO.read(str(aligned_fasta), "fasta")
    n_rows: int = len(alignment)
    n_cols: int = alignment.get_alignment_length()
    summaries: list[ColumnSummary] = []
    for j in range(n_cols):
        col: str = alignment[:, j]
        gap_fraction: float = col.count("-") / n_rows
        consensus, conservation = column_consensus(col)
        summaries.append(ColumnSummary(
            column_index=j,
            gap_fraction=gap_fraction,
            consensus_residue=consensus,
            conservation_fraction=conservation,
        ))
    return n_rows, n_cols, summaries


def write_column_summary_tsv(
    summaries: list[ColumnSummary],
    out_path: Path,
) -> None:
    """Write the per-column summary as a TSV. No pandas dependency."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("column_index\tgap_fraction\tconsensus_residue\tconservation_fraction\n")
        for s in summaries:
            fh.write(
                f"{s.column_index}\t{s.gap_fraction:.4f}\t"
                f"{s.consensus_residue}\t{s.conservation_fraction:.4f}\n"
            )


# ----------------------------------------------------------------------
# Stage 4 - trim gappy columns.
# ----------------------------------------------------------------------

def trim_alignment(
    aligned_fasta: Path,
    output_fasta: Path,
    max_gap_fraction: float = DEFAULT_TRIM_THRESHOLD,
) -> tuple[int, int]:
    """Drop columns where the gap fraction exceeds max_gap_fraction.

    Returns (n_cols_before, n_cols_after).
    """
    from Bio import AlignIO  # lazy import for py_compile safety
    from Bio.Align import MultipleSeqAlignment

    alignment = AlignIO.read(str(aligned_fasta), "fasta")
    n_rows: int = len(alignment)
    n_cols_before: int = alignment.get_alignment_length()
    keep_cols: list[int] = [
        j for j in range(n_cols_before)
        if alignment[:, j].count("-") / n_rows <= max_gap_fraction
    ]
    if not keep_cols:
        raise ValueError(
            "No columns survive the gap threshold; alignment is uniformly gappy."
        )

    # Slice each record to the kept columns.
    new_records = []
    for record in alignment:
        seq_str: str = str(record.seq)
        kept_str: str = "".join(seq_str[j] for j in keep_cols)
        new_record = record[:0]
        new_record.seq = type(record.seq)(kept_str)
        new_record.id = record.id
        new_record.name = record.name
        new_record.description = record.description
        new_records.append(new_record)

    trimmed = MultipleSeqAlignment(new_records)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    AlignIO.write(trimmed, str(output_fasta), "fasta")
    return n_cols_before, len(keep_cols)


# ----------------------------------------------------------------------
# Stage 5 - run-info JSON.
# ----------------------------------------------------------------------

def write_run_info(
    run_info: ExerciseRunInfo,
    out_path: Path,
) -> None:
    """Write the run-info JSON. Raises ValueError on empty required fields."""
    if not run_info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    if not run_info.mafft_version:
        raise ValueError("run_info.mafft_version is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    input_fasta: Path,
    out_dir: Path,
    algorithm_flags: list[str] | None = None,
    trim_threshold: float = DEFAULT_TRIM_THRESHOLD,
    threads: int = 4,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    aligned: Path = out_dir / "aligned.fasta"
    trimmed: Path = out_dir / "trimmed.fasta"
    summary_tsv: Path = out_dir / "alignment_summary.tsv"
    run_info_path: Path = out_dir / "run-info.json"

    n_input_records: int = validate_input_fasta(input_fasta)
    run_mafft(input_fasta, aligned, algorithm_flags=algorithm_flags, threads=threads)
    n_rows, n_cols_before, summaries = summarize_alignment(aligned)
    if n_rows != n_input_records:
        raise RuntimeError(
            f"MAFFT output row count {n_rows} does not match input record count "
            f"{n_input_records}; aborting."
        )
    write_column_summary_tsv(summaries, summary_tsv)
    n_cols_b, n_cols_after = trim_alignment(aligned, trimmed, max_gap_fraction=trim_threshold)
    assert n_cols_b == n_cols_before

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        input_fasta=str(input_fasta),
        input_md5=md5_of_file(input_fasta),
        mafft_version=get_mafft_version(),
        mafft_algorithm=" ".join(algorithm_flags or DEFAULT_MAFFT_ALGORITHM_FLAGS),
        trim_threshold=trim_threshold,
        n_input_records=n_input_records,
        n_aligned_columns_before_trim=n_cols_before,
        n_aligned_columns_after_trim=n_cols_after,
        biopython_version=biopython_version(),
        python_version=python_version_string(),
    )
    write_run_info(info, run_info_path)
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 1 - MAFFT via subprocess; align, summarize, trim.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input FASTA path.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--trim-threshold",
        type=float,
        default=DEFAULT_TRIM_THRESHOLD,
        help="Drop columns where the gap fraction exceeds this value. "
             f"Default: {DEFAULT_TRIM_THRESHOLD}.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="MAFFT thread count. Default: 4.",
    )
    parser.add_argument(
        "--linsi",
        action="store_true",
        help="Use L-INS-i (--localpair --maxiterate 1000) instead of FFT-NS-2.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    algorithm_flags: list[str] | None = None
    if args.linsi:
        algorithm_flags = [
            "--localpair",
            "--maxiterate", "1000",
            "--nuc",
            "--anysymbol",
            "--quiet",
        ]
    run_info_path: Path = run_exercise(
        input_fasta=args.input,
        out_dir=args.out_dir,
        algorithm_flags=algorithm_flags,
        trim_threshold=args.trim_threshold,
        threads=args.threads,
    )
    print(f"[ex01] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
