"""
Exercise 1 - Simulate long reads with badread, run Flye, parse the output.

Goal: take a small reference FASTA (1 Mb or smaller), simulate ~50x of
ONT-style long reads with badread, run Flye in --nano-hq mode, parse
the resulting assembly.fasta and assembly_info.txt, and emit a
run-info JSON recording every parameter.

The exercise covers:

- Calling badread via subprocess.run with stdout redirected to a file.
- Calling Flye via subprocess.run with --check and a pinned thread count.
- Parsing the FASTA assembly with Bio.SeqIO.
- Parsing Flye's tab-delimited assembly_info.txt by hand.
- Writing a run-info JSON alongside the assembly.

Estimated time: 75 minutes (25 min reading, 40 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-01-flye-via-subprocess.py --reference data/reference_1mb.fasta
    --out-dir results/ex01 --seed 42` runs end to end without errors.
- `results/ex01/reads.fastq` exists (badread output), > 1 MB.
- `results/ex01/flye_out/assembly.fasta` exists with at least one record.
- `results/ex01/flye_out/assembly_info.txt` exists and parses to at least one
  ContigInfo record.
- `results/ex01/run-info.json` exists with the badread version, the Flye version,
  the seed, the thread count, the input mode flag, and the run date.
- A second run with the same seed produces a byte-identical reads.fastq.
- A second run with the same reads.fastq produces an assembly.fasta whose
  contig count matches the first run.

Requirements:
    conda install -c bioconda flye=2.9.5 badread=0.4.1 biopython=1.84

What you learn:
- The subprocess.run idiom with stdout redirected to a file (badread) and
  with check=True (Flye).
- The Bio.SeqIO API for FASTA records.
- The Flye assembly_info.txt format and how to parse it.
- The "pin the seed, pin the thread count, write the run-info" reproducibility pattern.

Tool versions assumed:
- Python 3.11+
- badread 0.4.1 (CLI tool; subprocess.run it)
- Flye 2.9.5 (CLI tool; subprocess.run it)
- Biopython 1.84+

References:
- Flye: Kolmogorov et al. 2019, Nat Biotechnol 37:540
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699608/
- badread: Wick 2019, JOSS 4:1316
  https://joss.theoj.org/papers/10.21105/joss.01316
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

DEFAULT_FLYE_INPUT_MODE: str = "--nano-hq"
DEFAULT_GENOME_SIZE: str = "1m"
DEFAULT_BADREAD_LENGTH: str = "15000,13000"
DEFAULT_BADREAD_IDENTITY: str = "95,3,99"
DEFAULT_BADREAD_ERROR_MODEL: str = "nanopore2023"
DEFAULT_BADREAD_QSCORE_MODEL: str = "nanopore2023"
DEFAULT_COVERAGE: str = "50x"
DEFAULT_THREADS: int = 4
DEFAULT_FLYE_ITERATIONS: int = 1
DEFAULT_SEED: int = 42


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class ExerciseRunInfo:
    """Provenance metadata for the Exercise 1 run."""
    run_date: str = ""
    reference_fasta: str = ""
    reference_md5: str = ""
    reference_length_bp: int = 0
    badread_version: str = ""
    badread_error_model: str = DEFAULT_BADREAD_ERROR_MODEL
    badread_qscore_model: str = DEFAULT_BADREAD_QSCORE_MODEL
    badread_length: str = DEFAULT_BADREAD_LENGTH
    badread_identity: str = DEFAULT_BADREAD_IDENTITY
    badread_coverage: str = DEFAULT_COVERAGE
    seed: int = DEFAULT_SEED
    flye_version: str = ""
    flye_input_mode: str = DEFAULT_FLYE_INPUT_MODE
    flye_genome_size: str = DEFAULT_GENOME_SIZE
    flye_threads: int = DEFAULT_THREADS
    flye_iterations: int = DEFAULT_FLYE_ITERATIONS
    n_reads_simulated: int = 0
    n_contigs_assembled: int = 0
    total_assembly_length_bp: int = 0
    longest_contig_bp: int = 0
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


def get_tool_version(tool: str, version_flag: str = "--version") -> str:
    """Return the version line from `<tool> <version_flag>`.

    Many bioinformatics tools print to stderr instead of stdout; we
    capture both and return the first non-empty line.
    """
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
# Stage 1 - validate the input reference.
# ----------------------------------------------------------------------

def validate_reference(reference_fasta: Path) -> int:
    """Confirm the reference exists, is non-empty, and contains at least
    one record. Return the total length of all records.

    Raises FileNotFoundError or ValueError on any problem.
    """
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")
    from Bio import SeqIO  # lazy import for py_compile safety

    records = list(SeqIO.parse(str(reference_fasta), "fasta"))
    if not records:
        raise ValueError(f"Reference FASTA has no records: {reference_fasta}")
    total_length: int = sum(len(r.seq) for r in records)
    if total_length < 1000:
        raise ValueError(
            f"Reference FASTA total length {total_length} bp is below the "
            f"1,000 bp floor; the assembler will not produce a useful result."
        )
    return total_length


# ----------------------------------------------------------------------
# Stage 2 - simulate reads with badread.
# ----------------------------------------------------------------------

def run_badread_simulate(
    reference: Path,
    output_fastq: Path,
    coverage: str = DEFAULT_COVERAGE,
    length: str = DEFAULT_BADREAD_LENGTH,
    identity: str = DEFAULT_BADREAD_IDENTITY,
    error_model: str = DEFAULT_BADREAD_ERROR_MODEL,
    qscore_model: str = DEFAULT_BADREAD_QSCORE_MODEL,
    seed: int = DEFAULT_SEED,
) -> None:
    """Run badread simulate with stdout redirected to output_fastq.

    Raises subprocess.CalledProcessError if badread fails.
    Raises FileNotFoundError if badread is not on the PATH.
    """
    output_fastq.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "badread", "simulate",
        "--reference", str(reference),
        "--quantity", coverage,
        "--length", length,
        "--identity", identity,
        "--error_model", error_model,
        "--qscore_model", qscore_model,
        "--seed", str(seed),
    ]
    with output_fastq.open("w") as fh:
        result = subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
    if result.stderr:
        for line in result.stderr.splitlines():
            line = line.strip()
            if line.startswith("Generating") or line.startswith("Done"):
                print(f"[badread] {line}", file=sys.stderr)


def count_fastq_records(fastq_path: Path) -> int:
    """Return the number of FASTQ records in a file.

    Counts lines and divides by 4; cheaper than parsing every record.
    """
    n_lines: int = 0
    with fastq_path.open() as fh:
        for _ in fh:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(
            f"FASTQ at {fastq_path} has {n_lines} lines; not a multiple of 4."
        )
    return n_lines // 4


# ----------------------------------------------------------------------
# Stage 3 - run Flye.
# ----------------------------------------------------------------------

def run_flye(
    reads_fastq: Path,
    out_dir: Path,
    genome_size: str = DEFAULT_GENOME_SIZE,
    input_mode: str = DEFAULT_FLYE_INPUT_MODE,
    threads: int = DEFAULT_THREADS,
    iterations: int = DEFAULT_FLYE_ITERATIONS,
) -> Path:
    """Run Flye. Returns the path to assembly.fasta.

    Raises subprocess.CalledProcessError if Flye fails.
    Raises FileNotFoundError if Flye is not on the PATH.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "flye",
        input_mode, str(reads_fastq),
        "--genome-size", genome_size,
        "--out-dir", str(out_dir),
        "--threads", str(threads),
        "--iterations", str(iterations),
    ]
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    assembly_path: Path = out_dir / "assembly.fasta"
    if not assembly_path.exists():
        raise RuntimeError(
            f"Flye finished but {assembly_path} does not exist; check the log."
        )
    return assembly_path


# ----------------------------------------------------------------------
# Stage 4 - parse the assembly.
# ----------------------------------------------------------------------

@dataclass
class ContigInfo:
    """One row from Flye's assembly_info.txt."""
    seq_name: str
    length: int
    coverage: int
    circular: bool
    repeat: bool
    multiplicity: int
    alt_group: str
    graph_path: str


def parse_flye_assembly_info(info_path: Path) -> list[ContigInfo]:
    """Parse Flye's assembly_info.txt into a list of ContigInfo records.

    Skips the header (a line starting with '#') and blank lines.
    Flye writes the column order as:
        seq_name  length  cov.  circ.  repeat  mult.  alt_group  graph_path
    """
    rows: list[ContigInfo] = []
    if not info_path.exists():
        raise FileNotFoundError(f"assembly_info.txt not found: {info_path}")
    with info_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts: list[str] = line.split("\t")
            if len(parts) < 6:
                continue
            rows.append(ContigInfo(
                seq_name=parts[0],
                length=int(parts[1]),
                coverage=int(parts[2]) if parts[2].isdigit() else 0,
                circular=(parts[3] == "Y"),
                repeat=(parts[4] == "Y"),
                multiplicity=int(parts[5]) if parts[5].isdigit() else 0,
                alt_group=parts[6] if len(parts) > 6 else "",
                graph_path=parts[7] if len(parts) > 7 else "",
            ))
    return rows


def summarize_assembly(assembly_fasta: Path) -> tuple[int, int, int]:
    """Return (n_contigs, total_length_bp, longest_contig_bp) from the FASTA."""
    from Bio import SeqIO  # lazy import for py_compile safety

    records = list(SeqIO.parse(str(assembly_fasta), "fasta"))
    if not records:
        return (0, 0, 0)
    lengths: list[int] = [len(r.seq) for r in records]
    return (len(records), sum(lengths), max(lengths))


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
    if not run_info.flye_version:
        raise ValueError("run_info.flye_version is empty; refusing to write.")
    if not run_info.badread_version:
        raise ValueError("run_info.badread_version is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    reference_fasta: Path,
    out_dir: Path,
    seed: int = DEFAULT_SEED,
    coverage: str = DEFAULT_COVERAGE,
    flye_input_mode: str = DEFAULT_FLYE_INPUT_MODE,
    genome_size: str = DEFAULT_GENOME_SIZE,
    threads: int = DEFAULT_THREADS,
    iterations: int = DEFAULT_FLYE_ITERATIONS,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    reads_fastq: Path = out_dir / "reads.fastq"
    flye_out: Path = out_dir / "flye_out"
    run_info_path: Path = out_dir / "run-info.json"

    reference_length: int = validate_reference(reference_fasta)

    # Stage 2 - simulate reads.
    if not reads_fastq.exists() or reads_fastq.stat().st_size == 0:
        run_badread_simulate(
            reference=reference_fasta,
            output_fastq=reads_fastq,
            coverage=coverage,
            seed=seed,
        )
    n_reads: int = count_fastq_records(reads_fastq)

    # Stage 3 - run Flye.
    assembly_path: Path = run_flye(
        reads_fastq=reads_fastq,
        out_dir=flye_out,
        genome_size=genome_size,
        input_mode=flye_input_mode,
        threads=threads,
        iterations=iterations,
    )
    n_contigs, total_bp, longest_bp = summarize_assembly(assembly_path)

    # Stage 4 - parse assembly_info.txt and print a one-line summary.
    info_path: Path = flye_out / "assembly_info.txt"
    contigs: list[ContigInfo] = parse_flye_assembly_info(info_path)
    print(f"[ex01] Flye produced {len(contigs)} contigs:", file=sys.stderr)
    for c in contigs:
        circular: str = "circular" if c.circular else "linear"
        print(
            f"[ex01]   {c.seq_name}  length={c.length:>10d}  "
            f"coverage={c.coverage:>4d}  {circular}",
            file=sys.stderr,
        )

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        reference_fasta=str(reference_fasta),
        reference_md5=md5_of_file(reference_fasta),
        reference_length_bp=reference_length,
        badread_version=get_tool_version("badread"),
        seed=seed,
        flye_version=get_tool_version("flye"),
        flye_input_mode=flye_input_mode,
        flye_genome_size=genome_size,
        flye_threads=threads,
        flye_iterations=iterations,
        badread_coverage=coverage,
        n_reads_simulated=n_reads,
        n_contigs_assembled=n_contigs,
        total_assembly_length_bp=total_bp,
        longest_contig_bp=longest_bp,
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
        description="Exercise 1 - badread + Flye end to end.",
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
        default=DEFAULT_SEED,
        help=f"badread random seed. Default: {DEFAULT_SEED}.",
    )
    parser.add_argument(
        "--coverage",
        type=str,
        default=DEFAULT_COVERAGE,
        help=f"badread target coverage. Default: {DEFAULT_COVERAGE}.",
    )
    parser.add_argument(
        "--flye-input-mode",
        type=str,
        default=DEFAULT_FLYE_INPUT_MODE,
        choices=["--nano-hq", "--nano-raw", "--pacbio-hifi", "--pacbio-raw"],
        help=f"Flye input mode. Default: {DEFAULT_FLYE_INPUT_MODE}.",
    )
    parser.add_argument(
        "--genome-size",
        type=str,
        default=DEFAULT_GENOME_SIZE,
        help=f"Flye --genome-size hint. Default: {DEFAULT_GENOME_SIZE}.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help=f"Flye thread count. Default: {DEFAULT_THREADS}.",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=DEFAULT_FLYE_ITERATIONS,
        help=f"Flye internal polish iterations. Default: {DEFAULT_FLYE_ITERATIONS}.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        reference_fasta=args.reference,
        out_dir=args.out_dir,
        seed=args.seed,
        coverage=args.coverage,
        flye_input_mode=args.flye_input_mode,
        genome_size=args.genome_size,
        threads=args.threads,
        iterations=args.iterations,
    )
    print(f"[ex01] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
