"""
Exercise 2 - Compute N50, L50, contig stats, parse a BUSCO summary.

Goal: take the Flye assembly from Exercise 1, compute the assembly
statistics by hand (n_contigs, total length, longest contig, N50, L50,
GC fraction), optionally run BUSCO against bacteria_odb10, parse the
BUSCO short summary text file, and emit a Markdown QC report.

The exercise covers:

- Implementing N50 / L50 in plain Python without relying on third-party
  assembly-stats tools.
- Calling BUSCO via subprocess.run with --offline and a pinned lineage.
- Parsing BUSCO's plain-text short_summary file with regular expressions.
- Writing a Markdown QC report that names every parameter and every
  metric in the C10 lab-notebook voice.

Estimated time: 60 minutes (20 min reading, 30 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-02-asmstats-n50-busco.py --assembly results/ex01/flye_out/assembly.fasta
    --out-dir results/ex02 --lineage bacteria_odb10 --skip-busco` runs end to end.
- `results/ex02/asmstats.tsv` exists with one row per contig plus a summary row.
- `results/ex02/qc_report.md` exists with the asmstats table and the optional
  BUSCO block.
- `results/ex02/run-info.json` exists with the BUSCO version, the lineage,
  the dataset date, the n50 and l50 numbers, and the run date.

The --skip-busco flag is for offline use (or for the synthetic demo
reference which would score BUSCO C ~= 0). The compile-and-run does
not require BUSCO to be installed.

Requirements:
    conda install -c bioconda busco=5.7.1 biopython=1.84

What you learn:
- The N50 / L50 calculation by walking a sorted-descending length list.
- The BUSCO short_summary parser with regex.
- The "Markdown report driven by a dataclass" pattern.
- The graceful-skip pattern for optional external tools.

Tool versions assumed:
- Python 3.11+
- BUSCO 5.7.1 (CLI tool; optional)
- Biopython 1.84+

References:
- BUSCO: Manni et al. 2021, Mol Biol Evol 38:4647
  https://academic.oup.com/mbe/article/38/10/4647/6329644
- BUSCO lineage list:
  https://busco-data.ezlab.org/v5/data/lineages/
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_BUSCO_LINEAGE: str = "bacteria_odb10"
DEFAULT_BUSCO_MODE: str = "genome"
DEFAULT_THREADS: int = 4


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class AssemblyStats:
    """N50 / L50 plus the basic contig statistics."""
    n_contigs: int = 0
    total_length_bp: int = 0
    longest_contig_bp: int = 0
    shortest_contig_bp: int = 0
    n50_bp: int = 0
    l50: int = 0
    gc_fraction: float = 0.0


@dataclass
class BuscoSummary:
    """Parsed BUSCO short_summary fields."""
    busco_version: str = ""
    lineage_dataset: str = ""
    lineage_creation_date: str = ""
    n_total: int = 0
    n_complete: int = 0
    n_single_copy: int = 0
    n_duplicated: int = 0
    n_fragmented: int = 0
    n_missing: int = 0
    pct_complete: float = 0.0
    pct_single_copy: float = 0.0
    pct_duplicated: float = 0.0
    pct_fragmented: float = 0.0
    pct_missing: float = 0.0


@dataclass
class ExerciseRunInfo:
    """Provenance metadata for Exercise 2."""
    run_date: str = ""
    assembly_fasta: str = ""
    asmstats: AssemblyStats = None  # type: ignore[assignment]
    busco_summary: BuscoSummary = None  # type: ignore[assignment]
    busco_skipped: bool = False
    busco_lineage_requested: str = DEFAULT_BUSCO_LINEAGE
    biopython_version: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_date": self.run_date,
            "assembly_fasta": self.assembly_fasta,
            "asmstats": asdict(self.asmstats) if self.asmstats else {},
            "busco_summary": (
                asdict(self.busco_summary) if self.busco_summary else {}
            ),
            "busco_skipped": self.busco_skipped,
            "busco_lineage_requested": self.busco_lineage_requested,
            "biopython_version": self.biopython_version,
            "python_version": self.python_version,
        }


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

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
# Stage 1 - compute the assembly statistics.
# ----------------------------------------------------------------------

def compute_assembly_stats(fasta_path: Path) -> AssemblyStats:
    """Compute n_contigs, total length, longest, shortest, N50, L50, GC.

    For empty assemblies returns AssemblyStats with all-zero fields.
    """
    from Bio import SeqIO  # lazy import for py_compile safety

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        return AssemblyStats()

    lengths: list[int] = sorted((len(r.seq) for r in records), reverse=True)
    total: int = sum(lengths)

    cumulative: int = 0
    n50_bp: int = 0
    l50: int = 0
    for rank, length in enumerate(lengths, start=1):
        cumulative += length
        if cumulative >= total / 2:
            n50_bp = length
            l50 = rank
            break

    gc_bases: int = 0
    for r in records:
        seq_upper: str = str(r.seq).upper()
        gc_bases += seq_upper.count("G") + seq_upper.count("C")
    gc_fraction: float = gc_bases / total if total else 0.0

    return AssemblyStats(
        n_contigs=len(records),
        total_length_bp=total,
        longest_contig_bp=lengths[0],
        shortest_contig_bp=lengths[-1],
        n50_bp=n50_bp,
        l50=l50,
        gc_fraction=gc_fraction,
    )


def write_asmstats_tsv(
    fasta_path: Path,
    stats: AssemblyStats,
    out_path: Path,
) -> None:
    """Write one row per contig plus a summary row to a TSV.

    Columns: contig_name, length_bp, gc_fraction.
    """
    from Bio import SeqIO  # lazy import for py_compile safety

    out_path.parent.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    with out_path.open("w") as fh:
        fh.write("contig_name\tlength_bp\tgc_fraction\n")
        for r in records:
            seq_upper: str = str(r.seq).upper()
            gc: int = seq_upper.count("G") + seq_upper.count("C")
            gc_frac: float = gc / len(seq_upper) if len(seq_upper) else 0.0
            fh.write(f"{r.id}\t{len(r.seq)}\t{gc_frac:.4f}\n")
        # Summary line.
        fh.write(
            f"# n_contigs={stats.n_contigs}\ttotal_bp={stats.total_length_bp}\t"
            f"n50_bp={stats.n50_bp}\tl50={stats.l50}\n"
        )


# ----------------------------------------------------------------------
# Stage 2 - run BUSCO (optional).
# ----------------------------------------------------------------------

def run_busco(
    assembly_fasta: Path,
    out_dir: Path,
    lineage: str = DEFAULT_BUSCO_LINEAGE,
    mode: str = DEFAULT_BUSCO_MODE,
    threads: int = DEFAULT_THREADS,
    download_path: Path | None = None,
    offline: bool = True,
) -> Path:
    """Run BUSCO. Returns the path to the short_summary text file.

    Raises subprocess.CalledProcessError if BUSCO fails.
    Raises FileNotFoundError if BUSCO is not on the PATH.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    run_name: str = "busco_run"
    cmd: list[str] = [
        "busco",
        "-i", str(assembly_fasta),
        "-l", lineage,
        "-o", run_name,
        "--out_path", str(out_dir),
        "-m", mode,
        "-c", str(threads),
        "-f",  # force overwrite of any prior run
    ]
    if offline:
        cmd.append("--offline")
    if download_path is not None:
        cmd.extend(["--download_path", str(download_path)])
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )

    summary_path: Path = (
        out_dir / run_name / f"short_summary.specific.{lineage}.{run_name}.txt"
    )
    if not summary_path.exists():
        # Fall back to a glob if the filename convention drifts across BUSCO versions.
        candidates: list[Path] = list((out_dir / run_name).glob("short_summary*.txt"))
        if candidates:
            summary_path = candidates[0]
        else:
            raise RuntimeError(
                f"BUSCO finished but no short_summary found in {out_dir / run_name}."
            )
    return summary_path


def parse_busco_summary(summary_path: Path) -> BuscoSummary:
    """Parse the BUSCO short_summary text file."""
    text: str = summary_path.read_text()
    version_match = re.search(r"BUSCO version is:\s+(\S+)", text)
    lineage_match = re.search(
        r"lineage dataset is:\s+(\S+).*?Creation date:\s+([\d-]+)", text, re.S,
    )
    pct_match = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    n_complete_match = re.search(r"(\d+)\s+Complete BUSCOs \(C\)", text)
    n_single_match = re.search(
        r"(\d+)\s+Complete and single-copy BUSCOs \(S\)", text,
    )
    n_dup_match = re.search(
        r"(\d+)\s+Complete and duplicated BUSCOs \(D\)", text,
    )
    n_frag_match = re.search(r"(\d+)\s+Fragmented BUSCOs \(F\)", text)
    n_missing_match = re.search(r"(\d+)\s+Missing BUSCOs \(M\)", text)

    if not (
        version_match
        and lineage_match
        and pct_match
        and n_complete_match
        and n_single_match
        and n_dup_match
        and n_frag_match
        and n_missing_match
    ):
        raise ValueError(f"Could not parse BUSCO summary at {summary_path}")

    return BuscoSummary(
        busco_version=version_match.group(1),
        lineage_dataset=lineage_match.group(1),
        lineage_creation_date=lineage_match.group(2),
        n_total=int(pct_match.group(6)),
        n_complete=int(n_complete_match.group(1)),
        n_single_copy=int(n_single_match.group(1)),
        n_duplicated=int(n_dup_match.group(1)),
        n_fragmented=int(n_frag_match.group(1)),
        n_missing=int(n_missing_match.group(1)),
        pct_complete=float(pct_match.group(1)),
        pct_single_copy=float(pct_match.group(2)),
        pct_duplicated=float(pct_match.group(3)),
        pct_fragmented=float(pct_match.group(4)),
        pct_missing=float(pct_match.group(5)),
    )


# ----------------------------------------------------------------------
# Stage 3 - render the Markdown QC report.
# ----------------------------------------------------------------------

def render_qc_report(
    assembly_fasta: Path,
    stats: AssemblyStats,
    busco: BuscoSummary | None,
    busco_skipped: bool,
    out_path: Path,
) -> None:
    """Write a one-page Markdown QC report.

    The report names every metric and the parameters that produced it.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Assembly QC Report")
    lines.append("")
    lines.append(
        f"Generated: {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}"
    )
    lines.append("")
    lines.append("## Assembly statistics")
    lines.append("")
    lines.append(f"- Source: `{assembly_fasta}`")
    lines.append(f"- Contigs: {stats.n_contigs:,}")
    lines.append(f"- Total length: {stats.total_length_bp:,} bp")
    lines.append(f"- Longest contig: {stats.longest_contig_bp:,} bp")
    lines.append(f"- Shortest contig: {stats.shortest_contig_bp:,} bp")
    lines.append(f"- N50: {stats.n50_bp:,} bp")
    lines.append(f"- L50: {stats.l50}")
    lines.append(f"- GC fraction: {stats.gc_fraction:.4f}")
    lines.append("")

    if busco_skipped:
        lines.append("## BUSCO")
        lines.append("")
        lines.append("- Skipped (`--skip-busco` flag set or BUSCO not on PATH).")
        lines.append("")
    elif busco is not None:
        lines.append("## BUSCO")
        lines.append("")
        lines.append(f"- BUSCO version: {busco.busco_version}")
        lines.append(
            f"- Lineage: {busco.lineage_dataset} (created {busco.lineage_creation_date})"
        )
        lines.append(
            f"- Complete (C): {busco.pct_complete:.1f}% "
            f"(single-copy {busco.pct_single_copy:.1f}%, "
            f"duplicated {busco.pct_duplicated:.1f}%)"
        )
        lines.append(f"- Fragmented (F): {busco.pct_fragmented:.1f}%")
        lines.append(f"- Missing (M): {busco.pct_missing:.1f}%")
        lines.append(f"- Total BUSCO groups: {busco.n_total}")
        lines.append("")
    lines.append("## Limits and caveats")
    lines.append("")
    lines.append("- N50 alone can be inflated by collapsed repeats; cross-check")
    lines.append("  the per-contig coverage in `assembly_info.txt` against the")
    lines.append("  genome median. Coverage > 1.5x median is a flag.")
    lines.append("- BUSCO scores depend on the lineage dataset; pin the lineage")
    lines.append("  name AND the dataset creation date in the run-info JSON.")
    lines.append("- For real-data assemblies, validate against a reference where")
    lines.append("  possible; for novel genomes, use Merqury for reference-free QV.")
    lines.append("")
    lines.append("See `run-info.json` for full version and parameter pinning.")
    lines.append("")
    out_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------
# Stage 4 - run-info JSON.
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
    assembly_fasta: Path,
    out_dir: Path,
    lineage: str = DEFAULT_BUSCO_LINEAGE,
    skip_busco: bool = False,
    threads: int = DEFAULT_THREADS,
    download_path: Path | None = None,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    asmstats_tsv: Path = out_dir / "asmstats.tsv"
    qc_report_path: Path = out_dir / "qc_report.md"
    run_info_path: Path = out_dir / "run-info.json"

    if not assembly_fasta.exists():
        raise FileNotFoundError(f"Assembly FASTA not found: {assembly_fasta}")

    stats: AssemblyStats = compute_assembly_stats(assembly_fasta)
    write_asmstats_tsv(assembly_fasta, stats, asmstats_tsv)
    print(
        f"[ex02] n_contigs={stats.n_contigs} total_bp={stats.total_length_bp} "
        f"N50={stats.n50_bp} L50={stats.l50}",
        file=sys.stderr,
    )

    busco: BuscoSummary | None = None
    busco_skipped: bool = False
    if skip_busco:
        busco_skipped = True
        print("[ex02] BUSCO skipped by --skip-busco flag.", file=sys.stderr)
    else:
        try:
            summary_path: Path = run_busco(
                assembly_fasta=assembly_fasta,
                out_dir=out_dir / "busco_out",
                lineage=lineage,
                threads=threads,
                download_path=download_path,
            )
            busco = parse_busco_summary(summary_path)
            print(
                f"[ex02] BUSCO C={busco.pct_complete:.1f}% "
                f"F={busco.pct_fragmented:.1f}% M={busco.pct_missing:.1f}%",
                file=sys.stderr,
            )
        except FileNotFoundError:
            busco_skipped = True
            print("[ex02] BUSCO not on PATH; skipping.", file=sys.stderr)
        except subprocess.CalledProcessError as exc:
            busco_skipped = True
            print(
                f"[ex02] BUSCO call failed; skipping: {exc.stderr}",
                file=sys.stderr,
            )

    render_qc_report(
        assembly_fasta=assembly_fasta,
        stats=stats,
        busco=busco,
        busco_skipped=busco_skipped,
        out_path=qc_report_path,
    )

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        assembly_fasta=str(assembly_fasta),
        asmstats=stats,
        busco_summary=busco if busco is not None else BuscoSummary(),
        busco_skipped=busco_skipped,
        busco_lineage_requested=lineage,
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
        description="Exercise 2 - assembly statistics + optional BUSCO.",
    )
    parser.add_argument(
        "--assembly",
        type=Path,
        required=True,
        help="Path to the assembly FASTA (e.g. flye_out/assembly.fasta).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--lineage",
        type=str,
        default=DEFAULT_BUSCO_LINEAGE,
        help=f"BUSCO lineage dataset. Default: {DEFAULT_BUSCO_LINEAGE}.",
    )
    parser.add_argument(
        "--skip-busco",
        action="store_true",
        help="Skip the BUSCO call (for offline / no-BUSCO environments).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help=f"BUSCO thread count. Default: {DEFAULT_THREADS}.",
    )
    parser.add_argument(
        "--download-path",
        type=Path,
        default=None,
        help="Path to a pre-downloaded BUSCO lineage cache.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        assembly_fasta=args.assembly,
        out_dir=args.out_dir,
        lineage=args.lineage,
        skip_busco=args.skip_busco,
        threads=args.threads,
        download_path=args.download_path,
    )
    print(f"[ex02] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
