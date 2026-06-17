"""
Mini-project starter - End-to-end long-read assembly pipeline.

This pipeline takes a reference FASTA, simulates nanopore-style long reads
with badread, assembles them with Flye, polishes with Medaka, computes
N50 / L50 / GC, optionally runs BUSCO, computes QV against the reference,
and emits a run-info JSON that records every parameter needed to
reproduce the run. The output looks like a published assembly artifact;
it is reproducible only if the run-info travels with it.

This file is a SKELETON. It compiles cleanly under `python3 -m py_compile`
but most functions raise NotImplementedError. Your job is to fill in the
TODOs to produce a working `assemble_genome(reference_fasta, out_dir, seed)`
function that runs the eight-stage pipeline described in README.md:

  1. validate_reference(reference_fasta)
  2. simulate_reads(reference_fasta, out_dir, seed, coverage)
  3. filter_reads(input_fastq, output_fastq, min_length, min_mean_qv)
  4. run_flye(reads_filtered, out_dir, genome_size, input_mode, threads)
  5. run_medaka(reads_filtered, draft, out_dir, model, threads)
  6. compute_assembly_stats(polished_fasta)
  7. compute_qv_against_reference(polished_fasta, reference_fasta, out_dir)
  8. render_bandage_graph(gfa_path, out_path)

Each function has a docstring, a type signature, and a NotImplementedError
or a stub return. Replace the body with your implementation. The acceptance
criteria in README.md tell you what each function must produce.

Tool versions assumed:
- Python 3.11+
- badread 0.4.1 (CLI tool; subprocess.run it)
- Flye 2.9.5 (CLI tool)
- Medaka 1.12.0 (CLI tool; optional via --skip-medaka)
- minimap2 2.28 (CLI tool; for QV)
- BUSCO 5.7.1 (CLI tool; optional via --skip-busco)
- Bandage 0.9.0 (binary; optional)
- Biopython 1.84+

References:
- Flye: Kolmogorov et al. 2019, Nat Biotechnol 37:540
- Medaka: https://github.com/nanoporetech/medaka
- BUSCO: Manni et al. 2021, Mol Biol Evol 38:4647
- Bandage: Wick et al. 2015, Bioinformatics 31:3350
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
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
    reference_fasta: str = ""
    reference_md5: str = ""
    reference_length_bp: int = 0
    badread_version: str = ""
    badread_error_model: str = "nanopore2023"
    badread_qscore_model: str = "nanopore2023"
    badread_length: str = "15000,13000"
    badread_identity: str = "95,3,99"
    badread_coverage: str = "50x"
    seed: int = 42
    read_filter_min_length_bp: int = 1000
    read_filter_min_mean_qv: float = 10.0
    n_reads_simulated: int = 0
    n_reads_after_filter: int = 0
    flye_version: str = ""
    flye_input_mode: str = "--nano-hq"
    flye_genome_size: str = "1m"
    flye_threads: int = 4
    flye_iterations: int = 1
    medaka_version: str = ""
    medaka_model: str = "r1041_e82_400bps_sup_v4.3.0"
    medaka_skipped: bool = False
    minimap2_version: str = ""
    bandage_version: str = ""
    bandage_skipped: bool = False
    n_contigs: int = 0
    total_length_bp: int = 0
    n50_bp: int = 0
    l50: int = 0
    longest_contig_bp: int = 0
    gc_fraction: float = 0.0
    raw_qv: float = 0.0
    polished_qv: float = 0.0
    biopython_version: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Stage 1 - validate the reference.
# ----------------------------------------------------------------------

def validate_reference(reference_fasta: Path) -> int:
    """Confirm the reference FASTA exists and has total length >= 1,000 bp.

    Return the total length in bp.
    Raises FileNotFoundError or ValueError.
    """
    # TODO: implement.
    raise NotImplementedError("Stage 1 - validate_reference")


# ----------------------------------------------------------------------
# Stage 2 - simulate reads with badread.
# ----------------------------------------------------------------------

def simulate_reads(
    reference: Path,
    output_fastq: Path,
    coverage: str = "50x",
    length: str = "15000,13000",
    identity: str = "95,3,99",
    error_model: str = "nanopore2023",
    qscore_model: str = "nanopore2023",
    seed: int = 42,
) -> None:
    """Run badread simulate with stdout redirected to output_fastq.

    Raises subprocess.CalledProcessError if badread fails.
    """
    # TODO: build the subprocess.run command, redirect stdout to the output
    # file, propagate stderr selectively.
    raise NotImplementedError("Stage 2 - simulate_reads")


# ----------------------------------------------------------------------
# Stage 3 - filter reads by length and quality.
# ----------------------------------------------------------------------

def filter_reads(
    input_fastq: Path,
    output_fastq: Path,
    min_length: int = 1000,
    min_mean_qv: float = 10.0,
) -> tuple[int, int]:
    """Drop reads shorter than min_length or with mean QV below min_mean_qv.

    Returns (n_in, n_out).
    """
    # TODO: parse with Bio.SeqIO, filter, write.
    raise NotImplementedError("Stage 3 - filter_reads")


# ----------------------------------------------------------------------
# Stage 4 - run Flye.
# ----------------------------------------------------------------------

def run_flye(
    reads_fastq: Path,
    out_dir: Path,
    genome_size: str = "1m",
    input_mode: str = "--nano-hq",
    threads: int = 4,
    iterations: int = 1,
) -> Path:
    """Run Flye. Returns the path to assembly.fasta.

    Raises subprocess.CalledProcessError if Flye fails.
    """
    # TODO: build the subprocess.run command, check=True, capture_output=True.
    raise NotImplementedError("Stage 4 - run_flye")


# ----------------------------------------------------------------------
# Stage 5 - polish with Medaka.
# ----------------------------------------------------------------------

def assert_medaka_chemistry_matches(
    medaka_model: str,
    expected_chemistry_prefix: str = "r1041",
) -> None:
    """Refuse to run Medaka with a chemistry-mismatched model."""
    medaka_prefix: str = medaka_model.split("_")[0].lower()
    expected: str = expected_chemistry_prefix.lower()
    if medaka_prefix != expected:
        raise ValueError(
            f"Medaka model {medaka_model} has chemistry prefix '{medaka_prefix}' "
            f"but expected '{expected}'. Mismatching chemistry produces a worse "
            f"polish than no polish; aborting."
        )


def run_medaka(
    reads_fastq: Path,
    draft_fasta: Path,
    out_dir: Path,
    medaka_model: str = "r1041_e82_400bps_sup_v4.3.0",
    threads: int = 4,
) -> Path:
    """Run medaka_consensus. Returns the polished FASTA path.

    Raises subprocess.CalledProcessError if Medaka fails.
    """
    # TODO: assert the chemistry matches; build the subprocess.run command;
    # check=True; capture_output=True; return out_dir / 'consensus.fasta'.
    raise NotImplementedError("Stage 5 - run_medaka")


# ----------------------------------------------------------------------
# Stage 6 - compute assembly statistics.
# ----------------------------------------------------------------------

@dataclass
class AssemblyStats:
    n_contigs: int = 0
    total_length_bp: int = 0
    longest_contig_bp: int = 0
    shortest_contig_bp: int = 0
    n50_bp: int = 0
    l50: int = 0
    gc_fraction: float = 0.0


def compute_assembly_stats(fasta_path: Path) -> AssemblyStats:
    """Compute n_contigs, total length, longest, shortest, N50, L50, GC.

    For empty assemblies returns AssemblyStats with all-zero fields.
    """
    # TODO: parse with Bio.SeqIO, sort lengths descending, walk cumulative.
    raise NotImplementedError("Stage 6 - compute_assembly_stats")


# ----------------------------------------------------------------------
# Stage 7 - QV against reference via minimap2 + SAM tally.
# ----------------------------------------------------------------------

def run_minimap2_asm(
    reference: Path,
    query: Path,
    output_sam: Path,
    preset: str = "asm5",
    threads: int = 4,
) -> None:
    """Run minimap2 in assembly mode and write SAM to output_sam."""
    # TODO: build the subprocess.run command, redirect stdout to output_sam.
    raise NotImplementedError("Stage 7a - run_minimap2_asm")


@dataclass
class QvResult:
    aligned_bp: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0
    error_rate: float = 0.0
    qv: float = 0.0


def tally_sam(sam_path: Path) -> QvResult:
    """Parse SAM, sum CIGAR and NM tags, compute QV = -10 * log10(error)."""
    # TODO: walk the SAM lines, parse CIGAR with a regex, prefer the NM tag
    # if present, return the QvResult.
    raise NotImplementedError("Stage 7b - tally_sam")


# ----------------------------------------------------------------------
# Stage 8 - Bandage graph rendering (optional).
# ----------------------------------------------------------------------

def render_bandage_graph(
    gfa_path: Path,
    output_png: Path,
    height: int = 800,
) -> bool:
    """Run Bandage image. Returns True on success, False if Bandage is unavailable.

    The function should NOT raise on missing Bandage; print a warning and
    return False instead. The pipeline survives without the graph image.
    """
    # TODO: try subprocess.run; on FileNotFoundError or CalledProcessError,
    # print a warning and return False; on success, return True.
    raise NotImplementedError("Stage 8 - render_bandage_graph")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def md5_of_file(path: Path) -> str:
    """Return the hex MD5 of a file's bytes."""
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def get_tool_version(tool: str, version_flag: str = "--version") -> str:
    """Return the first non-empty version line from `<tool> <version_flag>`."""
    try:
        result = subprocess.run(
            [tool, version_flag],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return "unknown (not on PATH)"
    blob: str = (result.stdout or "") + "\n" + (result.stderr or "")
    for line in blob.splitlines():
        line = line.strip()
        if line:
            return line
    return "unknown"


def biopython_version() -> str:
    try:
        import Bio  # type: ignore[import-not-found]
        return getattr(Bio, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    return ".".join(str(part) for part in sys.version_info[:3])


def assemble_genome(
    reference_fasta: Path,
    out_dir: Path,
    seed: int = 42,
    coverage: str = "50x",
    flye_input_mode: str = "--nano-hq",
    genome_size: str = "1m",
    threads: int = 4,
    medaka_model: str = "r1041_e82_400bps_sup_v4.3.0",
    skip_medaka: bool = False,
    skip_bandage: bool = False,
    min_read_length: int = 1000,
    min_read_qv: float = 10.0,
) -> Path:
    """Run the full pipeline. Returns the path to the polished FASTA.

    Writes results/run-info.json with the full provenance.
    """
    # TODO: chain the stages together. Each stage writes its output under
    # `out_dir`. The orchestrator fills in the PipelineRunInfo and writes
    # run-info.json at the end. Use the chemistry-mismatch guard on Medaka.
    # Surface a one-line summary to stderr after each stage.
    raise NotImplementedError("Orchestrator - assemble_genome")


def write_run_info(
    run_info: PipelineRunInfo,
    out_path: Path,
) -> None:
    """Write the run-info JSON. Refuses to write if required fields are empty."""
    if not run_info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    if not run_info.flye_version:
        raise ValueError("run_info.flye_version is empty; refusing to write.")
    if not run_info.badread_version:
        raise ValueError("run_info.badread_version is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


def render_qc_report(
    run_info: PipelineRunInfo,
    out_path: Path,
) -> None:
    """Render the Markdown QC report from the PipelineRunInfo dataclass.

    The report has sections: Input, Assembly, Polish, QC, Graph, Limits.
    """
    # TODO: build a Markdown string with each section; write to out_path.
    raise NotImplementedError("QC report rendering")


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Mini-project - end-to-end long-read assembly pipeline.",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Input reference FASTA path.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="badread random seed. Default: 42.",
    )
    parser.add_argument(
        "--coverage",
        type=str,
        default="50x",
        help="Target read coverage. Default: 50x.",
    )
    parser.add_argument(
        "--flye-input-mode",
        type=str,
        default="--nano-hq",
        choices=["--nano-hq", "--nano-raw", "--pacbio-hifi", "--pacbio-raw"],
        help="Flye input mode. Default: --nano-hq.",
    )
    parser.add_argument(
        "--genome-size",
        type=str,
        default="1m",
        help="Flye --genome-size hint. Default: 1m.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Worker thread count. Default: 4.",
    )
    parser.add_argument(
        "--medaka-model",
        type=str,
        default="r1041_e82_400bps_sup_v4.3.0",
        help="Medaka model name. Default: r1041_e82_400bps_sup_v4.3.0.",
    )
    parser.add_argument(
        "--skip-medaka",
        action="store_true",
        help="Skip Medaka polish (for offline / no-Medaka environments).",
    )
    parser.add_argument(
        "--skip-bandage",
        action="store_true",
        help="Skip Bandage graph rendering.",
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=1000,
        help="Minimum read length after filter. Default: 1000.",
    )
    parser.add_argument(
        "--min-read-qv",
        type=float,
        default=10.0,
        help="Minimum mean Phred quality after filter. Default: 10.0.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    polished_path: Path = assemble_genome(
        reference_fasta=args.reference,
        out_dir=args.out_dir,
        seed=args.seed,
        coverage=args.coverage,
        flye_input_mode=args.flye_input_mode,
        genome_size=args.genome_size,
        threads=args.threads,
        medaka_model=args.medaka_model,
        skip_medaka=args.skip_medaka,
        skip_bandage=args.skip_bandage,
        min_read_length=args.min_read_length,
        min_read_qv=args.min_read_qv,
    )
    print(f"[pipeline] polished FASTA at {polished_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
