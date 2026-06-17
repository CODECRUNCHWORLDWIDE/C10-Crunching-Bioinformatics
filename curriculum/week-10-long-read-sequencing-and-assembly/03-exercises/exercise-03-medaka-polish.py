"""
Exercise 3 - Polish a Flye assembly with Medaka and quantify the improvement.

Goal: take the Flye assembly from Exercise 1, run Medaka with a pinned
basecaller-matched model, and compute the QV improvement against the
original reference. Emit a Markdown polish report that names the model
used, the raw-vs-polished QV, and the per-error-class breakdown.

The exercise covers:

- Calling medaka_consensus via subprocess.run with a pinned model name.
- Aligning the assembly to the reference with minimap2 (-x asm5).
- Parsing the resulting SAM to tally matches, mismatches, and indels.
- Computing QV = -10 * log10(error_rate) by hand.
- Writing a Markdown comparison report with raw and polished QVs side-by-side.

Estimated time: 60 minutes (15 min reading, 35 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-03-medaka-polish.py --reads results/ex01/reads.fastq
    --draft results/ex01/flye_out/assembly.fasta --reference data/reference_1mb.fasta
    --out-dir results/ex03 --skip-medaka` runs end to end (even without Medaka
  installed; the --skip-medaka flag uses the draft as the "polished" output).
- `results/ex03/polished.fasta` exists (either the medaka_consensus output or
  the draft as a fallback).
- `results/ex03/qv_report.md` exists with raw_qv, polished_qv, and the
  per-error-class table.
- `results/ex03/run-info.json` exists with the Medaka model, the minimap2 version,
  and the QV numbers.

Requirements:
    conda install -c bioconda medaka=1.12.0 minimap2=2.28 biopython=1.84

What you learn:
- The medaka_consensus CLI and the model-name pinning pattern.
- The minimap2 -x asm5 preset for assembly-to-reference alignment.
- Parsing CIGAR strings and the NM tag to tally edit distance.
- The QV calculation: QV = -10 * log10(error_rate).
- The graceful-skip pattern for optional tools.

Tool versions assumed:
- Python 3.11+
- Medaka 1.12.0 (CLI tool; optional)
- minimap2 2.28 (CLI tool; required if --compare-to-reference is set)
- Biopython 1.84+

References:
- Medaka: https://github.com/nanoporetech/medaka
- minimap2: Li 2018, Bioinformatics 34:3094
  https://academic.oup.com/bioinformatics/article/34/18/3094/4994778
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_MEDAKA_MODEL: str = "r1041_e82_400bps_sup_v4.3.0"
DEFAULT_THREADS: int = 4
DEFAULT_MINIMAP2_PRESET: str = "asm5"


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class QvResult:
    """Result of comparing an assembly to a reference."""
    aligned_bp: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0
    error_rate: float = 0.0
    qv: float = 0.0


@dataclass
class ExerciseRunInfo:
    """Provenance metadata for Exercise 3."""
    run_date: str = ""
    reads_fastq: str = ""
    draft_fasta: str = ""
    polished_fasta: str = ""
    reference_fasta: str = ""
    medaka_version: str = ""
    medaka_model: str = DEFAULT_MEDAKA_MODEL
    medaka_skipped: bool = False
    minimap2_version: str = ""
    minimap2_preset: str = DEFAULT_MINIMAP2_PRESET
    threads: int = DEFAULT_THREADS
    raw_qv: QvResult = None  # type: ignore[assignment]
    polished_qv: QvResult = None  # type: ignore[assignment]
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_date": self.run_date,
            "reads_fastq": self.reads_fastq,
            "draft_fasta": self.draft_fasta,
            "polished_fasta": self.polished_fasta,
            "reference_fasta": self.reference_fasta,
            "medaka_version": self.medaka_version,
            "medaka_model": self.medaka_model,
            "medaka_skipped": self.medaka_skipped,
            "minimap2_version": self.minimap2_version,
            "minimap2_preset": self.minimap2_preset,
            "threads": self.threads,
            "raw_qv": asdict(self.raw_qv) if self.raw_qv else {},
            "polished_qv": asdict(self.polished_qv) if self.polished_qv else {},
            "python_version": self.python_version,
        }


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

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


def python_version_string() -> str:
    """Return the running Python version, e.g. '3.11.9'."""
    return ".".join(str(part) for part in sys.version_info[:3])


def assert_medaka_chemistry_matches(
    medaka_model: str,
    expected_chemistry_prefix: str = "r1041",
) -> None:
    """Refuse to run Medaka with a chemistry-mismatched model.

    Compares the leading 'r10' or 'r9' prefix and raises if they disagree.
    The default expects R10.4.1; override the prefix for older chemistries.
    """
    medaka_prefix: str = medaka_model.split("_")[0].lower()
    expected_prefix: str = expected_chemistry_prefix.lower()
    if medaka_prefix != expected_prefix:
        raise ValueError(
            f"Medaka model {medaka_model} has chemistry prefix '{medaka_prefix}' "
            f"but expected '{expected_prefix}'. Mismatching chemistry produces "
            f"a worse polish than no polish; aborting."
        )


# ----------------------------------------------------------------------
# Stage 1 - run Medaka.
# ----------------------------------------------------------------------

def run_medaka(
    reads_fastq: Path,
    draft_fasta: Path,
    out_dir: Path,
    medaka_model: str = DEFAULT_MEDAKA_MODEL,
    threads: int = DEFAULT_THREADS,
) -> Path:
    """Run medaka_consensus. Returns the polished FASTA path.

    Raises subprocess.CalledProcessError if Medaka fails.
    Raises FileNotFoundError if Medaka is not on the PATH.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "medaka_consensus",
        "-i", str(reads_fastq),
        "-d", str(draft_fasta),
        "-o", str(out_dir),
        "-m", medaka_model,
        "-t", str(threads),
        "-f",  # force-overwrite the output directory
    ]
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    polished_path: Path = out_dir / "consensus.fasta"
    if not polished_path.exists():
        raise RuntimeError(
            f"Medaka finished but {polished_path} does not exist; check the log."
        )
    return polished_path


# ----------------------------------------------------------------------
# Stage 2 - align assembly to reference with minimap2.
# ----------------------------------------------------------------------

def run_minimap2_asm(
    reference: Path,
    query: Path,
    output_sam: Path,
    preset: str = DEFAULT_MINIMAP2_PRESET,
    threads: int = DEFAULT_THREADS,
) -> None:
    """Run minimap2 in assembly mode. Writes SAM to output_sam.

    Raises subprocess.CalledProcessError if minimap2 fails.
    Raises FileNotFoundError if minimap2 is not on the PATH.
    """
    output_sam.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "minimap2",
        "-a",
        "-x", preset,
        "-t", str(threads),
        str(reference),
        str(query),
    ]
    with output_sam.open("w") as fh:
        subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )


# ----------------------------------------------------------------------
# Stage 3 - parse SAM, tally edit distance, compute QV.
# ----------------------------------------------------------------------

_CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")


def tally_sam_alignments(sam_path: Path) -> QvResult:
    """Parse a SAM file, sum CIGAR / NM tags, compute QV.

    Skips unmapped (flag 0x4) and supplementary (flag 0x800) records.
    Prefers the NM tag for total edit distance if the CIGAR uses 'M'
    (match-or-mismatch) rather than '=' / 'X' (explicit match / mismatch).
    """
    aligned_bp: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0

    if not sam_path.exists():
        return QvResult()

    with sam_path.open() as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            parts: list[str] = line.split("\t")
            if len(parts) < 11:
                continue
            flag: int = int(parts[1])
            if flag & 0x4 or flag & 0x800:
                continue
            cigar: str = parts[5]
            if cigar == "*":
                continue

            cigar_aligned_bp: int = 0
            cigar_explicit_mismatch: int = 0
            cigar_ins: int = 0
            cigar_del: int = 0
            uses_eq_x: bool = ("=" in cigar) or ("X" in cigar)
            for length_str, op in _CIGAR_TOKEN.findall(cigar):
                length: int = int(length_str)
                if op in ("M", "=", "X"):
                    cigar_aligned_bp += length
                if op == "X":
                    cigar_explicit_mismatch += length
                if op == "I":
                    cigar_ins += length
                if op == "D":
                    cigar_del += length

            nm: int | None = None
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag.split(":")[2])
                    break

            aligned_bp += cigar_aligned_bp
            insertions += cigar_ins
            deletions += cigar_del
            if uses_eq_x:
                mismatches += cigar_explicit_mismatch
            elif nm is not None:
                # NM = mismatches + insertions + deletions; subtract indels.
                m: int = max(0, nm - cigar_ins - cigar_del)
                mismatches += m
            else:
                # Conservative fallback: assume zero mismatches if no NM and no X.
                pass

    if aligned_bp == 0:
        return QvResult()
    total_errors: int = mismatches + insertions + deletions
    error_rate: float = total_errors / max(1, aligned_bp + insertions + deletions)
    if error_rate <= 0:
        qv: float = 60.0
    else:
        qv = -10.0 * math.log10(error_rate)
    return QvResult(
        aligned_bp=aligned_bp,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        error_rate=error_rate,
        qv=qv,
    )


# ----------------------------------------------------------------------
# Stage 4 - render the polish report.
# ----------------------------------------------------------------------

def render_qv_report(
    raw_qv: QvResult,
    polished_qv: QvResult,
    medaka_model: str,
    medaka_skipped: bool,
    out_path: Path,
) -> None:
    """Write a one-page Markdown polish report."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Polish QV Report")
    lines.append("")
    lines.append(
        f"Generated: {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}"
    )
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    if medaka_skipped:
        lines.append("- Medaka: SKIPPED (no medaka_consensus binary available).")
        lines.append("- 'Polished' QV below equals raw QV.")
    else:
        lines.append(f"- Medaka model: `{medaka_model}`")
    lines.append("")
    lines.append("## QV comparison")
    lines.append("")
    lines.append("| Metric         |       Raw |  Polished |")
    lines.append("|----------------|----------:|----------:|")
    lines.append(
        f"| Aligned bp     | {raw_qv.aligned_bp:>9,} | {polished_qv.aligned_bp:>9,} |"
    )
    lines.append(
        f"| Mismatches     | {raw_qv.mismatches:>9,} | {polished_qv.mismatches:>9,} |"
    )
    lines.append(
        f"| Insertions     | {raw_qv.insertions:>9,} | {polished_qv.insertions:>9,} |"
    )
    lines.append(
        f"| Deletions      | {raw_qv.deletions:>9,} | {polished_qv.deletions:>9,} |"
    )
    lines.append(
        f"| Error rate     | {raw_qv.error_rate:>9.6f} | {polished_qv.error_rate:>9.6f} |"
    )
    lines.append(
        f"| QV (Phred)     | {raw_qv.qv:>9.2f} | {polished_qv.qv:>9.2f} |"
    )
    lines.append("")
    if not medaka_skipped:
        if polished_qv.qv > raw_qv.qv:
            delta: float = polished_qv.qv - raw_qv.qv
            lines.append(
                f"Polish improved QV by +{delta:.2f} (lower error rate by "
                f"~{10 ** (delta / 10):.1f}x)."
            )
        else:
            lines.append(
                "Polish did NOT improve QV. Likely causes: model mismatch with "
                "the basecaller, very low coverage, or a polish step that "
                "introduced indels at homopolymers."
            )
    lines.append("")
    lines.append("## Limits and caveats")
    lines.append("")
    lines.append("- QV is computed against a reference; if the 'reference' is not")
    lines.append("  truly the source genome the QV is biased downward by real biology.")
    lines.append("- Medaka model must match the basecaller model that produced the")
    lines.append("  reads. The CLI enforces a chemistry-prefix check before running.")
    lines.append("- Homopolymer-region indels remain the dominant residual error in")
    lines.append("  ONT polishes; expect ~0.01% even after a correct Medaka run.")
    lines.append("")
    out_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------
# Stage 5 - run-info JSON.
# ----------------------------------------------------------------------

def write_run_info(
    run_info: ExerciseRunInfo,
    out_path: Path,
) -> None:
    """Write the run-info JSON. Refuses to write if run_date is empty."""
    if not run_info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    reads_fastq: Path,
    draft_fasta: Path,
    reference_fasta: Path,
    out_dir: Path,
    medaka_model: str = DEFAULT_MEDAKA_MODEL,
    threads: int = DEFAULT_THREADS,
    skip_medaka: bool = False,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    medaka_out_dir: Path = out_dir / "medaka_out"
    polished_fasta: Path = out_dir / "polished.fasta"
    raw_sam: Path = out_dir / "raw_vs_ref.sam"
    polished_sam: Path = out_dir / "polished_vs_ref.sam"
    qv_report_path: Path = out_dir / "qv_report.md"
    run_info_path: Path = out_dir / "run-info.json"

    if not draft_fasta.exists():
        raise FileNotFoundError(f"Draft FASTA not found: {draft_fasta}")
    if not reads_fastq.exists():
        raise FileNotFoundError(f"Reads FASTQ not found: {reads_fastq}")
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")

    # Stage 1 - Medaka polish (or skip).
    medaka_skipped: bool = skip_medaka
    if skip_medaka:
        # Copy the draft to polished.fasta so the rest of the pipeline can proceed.
        polished_fasta.write_bytes(draft_fasta.read_bytes())
        print("[ex03] Medaka skipped by --skip-medaka flag.", file=sys.stderr)
    else:
        try:
            assert_medaka_chemistry_matches(medaka_model)
            produced: Path = run_medaka(
                reads_fastq=reads_fastq,
                draft_fasta=draft_fasta,
                out_dir=medaka_out_dir,
                medaka_model=medaka_model,
                threads=threads,
            )
            polished_fasta.write_bytes(produced.read_bytes())
            print(f"[ex03] Medaka wrote {produced}.", file=sys.stderr)
        except FileNotFoundError:
            medaka_skipped = True
            polished_fasta.write_bytes(draft_fasta.read_bytes())
            print("[ex03] medaka_consensus not on PATH; skipping.", file=sys.stderr)
        except subprocess.CalledProcessError as exc:
            medaka_skipped = True
            polished_fasta.write_bytes(draft_fasta.read_bytes())
            print(f"[ex03] Medaka failed; skipping: {exc.stderr}", file=sys.stderr)

    # Stage 2 - align both drafts to the reference and tally QV.
    raw_qv: QvResult = QvResult()
    polished_qv: QvResult = QvResult()
    minimap2_v: str = "unknown (not on PATH)"
    try:
        minimap2_v = get_tool_version("minimap2")
        run_minimap2_asm(
            reference=reference_fasta,
            query=draft_fasta,
            output_sam=raw_sam,
            threads=threads,
        )
        run_minimap2_asm(
            reference=reference_fasta,
            query=polished_fasta,
            output_sam=polished_sam,
            threads=threads,
        )
        raw_qv = tally_sam_alignments(raw_sam)
        polished_qv = tally_sam_alignments(polished_sam)
        print(
            f"[ex03] raw QV={raw_qv.qv:.2f}; polished QV={polished_qv.qv:.2f}",
            file=sys.stderr,
        )
    except FileNotFoundError:
        print(
            "[ex03] minimap2 not on PATH; QV report will have zero values.",
            file=sys.stderr,
        )

    render_qv_report(
        raw_qv=raw_qv,
        polished_qv=polished_qv,
        medaka_model=medaka_model,
        medaka_skipped=medaka_skipped,
        out_path=qv_report_path,
    )

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        reads_fastq=str(reads_fastq),
        draft_fasta=str(draft_fasta),
        polished_fasta=str(polished_fasta),
        reference_fasta=str(reference_fasta),
        medaka_version=get_tool_version("medaka", "--version"),
        medaka_model=medaka_model,
        medaka_skipped=medaka_skipped,
        minimap2_version=minimap2_v,
        threads=threads,
        raw_qv=raw_qv,
        polished_qv=polished_qv,
        python_version=python_version_string(),
    )
    write_run_info(info, run_info_path)
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 3 - Medaka polish + QV comparison against a reference.",
    )
    parser.add_argument(
        "--reads",
        type=Path,
        required=True,
        help="Input nanopore reads FASTQ (the reads that produced the draft).",
    )
    parser.add_argument(
        "--draft",
        type=Path,
        required=True,
        help="Draft assembly FASTA (e.g. flye_out/assembly.fasta).",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Reference FASTA for QV calculation.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--medaka-model",
        type=str,
        default=DEFAULT_MEDAKA_MODEL,
        help=f"Medaka model name. Default: {DEFAULT_MEDAKA_MODEL}.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help=f"Worker thread count. Default: {DEFAULT_THREADS}.",
    )
    parser.add_argument(
        "--skip-medaka",
        action="store_true",
        help="Skip the Medaka call (for offline / no-Medaka environments).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        reads_fastq=args.reads,
        draft_fasta=args.draft,
        reference_fasta=args.reference,
        out_dir=args.out_dir,
        medaka_model=args.medaka_model,
        threads=args.threads,
        skip_medaka=args.skip_medaka,
    )
    print(f"[ex03] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
