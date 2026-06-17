"""
Exercise 2 - Run FilterMutectCalls and CalculateContamination, then report
the filter-reason tally.

Educational and research use only. Do not use any output of this pipeline
to guide patient care.

Goal: take the unfiltered Mutect2 VCF from Exercise 1, run
GetPileupSummaries on both BAMs against a common-biallelic VCF, run
CalculateContamination to produce a contamination.table, run
FilterMutectCalls with --contamination-table, parse the filtered VCF,
and emit a tally of FILTER reasons with a Markdown summary.

The exercise covers:

- Calling GATK GetPileupSummaries via subprocess for each sample.
- Calling GATK CalculateContamination via subprocess (matched-normal mode).
- Calling GATK FilterMutectCalls via subprocess with --contamination-table
  and --stats.
- Reading the contamination.table TSV by hand.
- Tallying the FILTER column of the filtered VCF via pysam.
- Emitting a Markdown summary alongside the run-info JSON.

Estimated time: 90 minutes (30 min reading, 50 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-02-filter-and-contamination.py
    --reference data/chr22_GRCh38.fasta
    --unfiltered-vcf results/ex01/unfiltered.vcf.gz
    --tumor-bam data/tumor_chr22.bam
    --normal-bam data/normal_chr22.bam
    --common-biallelic data/chr22_common_biallelic.vcf.gz
    --out-dir results/ex02` runs end to end without errors when
  GATK is on the PATH; gracefully skips otherwise.
- `results/ex02/contamination.table` exists.
- `results/ex02/filtered.vcf.gz` exists.
- `results/ex02/filter_tally.md` shows the count for each FILTER reason.
- `results/ex02/run-info.json` records the contamination estimate, the
  PASS variant count, and the per-filter counts.

Requirements:
    conda install -c bioconda gatk4=4.5.0.0 samtools=1.20 pysam=0.22.1

What you learn:
- The two-step pileup + contamination flow.
- The FilterMutectCalls filter set and how to tally it.
- The Markdown-report-alongside-JSON output convention.

Tool versions assumed:
- Python 3.11+
- GATK 4.5.0.0
- samtools 1.20
- pysam 0.22.1

References:
- Mutect2: Cibulskis et al. 2013, Nat Biotechnol 31:213
- FilterMutectCalls:
  https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-FilterMutectCalls
- GetPileupSummaries:
  https://gatk.broadinstitute.org/hc/en-us/articles/360037593871-GetPileupSummaries
- CalculateContamination:
  https://gatk.broadinstitute.org/hc/en-us/articles/360037593751-CalculateContamination
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_JAVA_HEAP: str = "4g"


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class ExerciseRunInfo:
    """Provenance metadata for the Exercise 2 run."""
    run_date: str = ""
    unfiltered_vcf: str = ""
    filtered_vcf: str = ""
    tumor_bam: str = ""
    normal_bam: str = ""
    reference_fasta: str = ""
    common_biallelic_vcf: str = ""
    contamination_table: str = ""
    contamination_fraction: float = 0.0
    contamination_error: float = 0.0
    gatk_version: str = ""
    pysam_version: str = ""
    java_heap: str = DEFAULT_JAVA_HEAP
    n_total: int = 0
    n_pass: int = 0
    filter_counts: dict[str, int] = field(default_factory=dict)
    skipped: bool = False
    skip_reason: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

def md5_of_file(path: Path) -> str:
    """Return the hex MD5 of a file, or '' if missing or unreadable."""
    if not path.exists():
        return ""
    h = hashlib.md5()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 16), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return ""


def get_tool_version(tool: str, version_flag: str = "--version") -> str:
    """Return the first non-empty line from `<tool> <version_flag>`."""
    if not shutil.which(tool):
        return "unknown (not on PATH)"
    try:
        result = subprocess.run(
            [tool, version_flag],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return "unknown (call failed)"
    blob: str = (result.stdout or "") + "\n" + (result.stderr or "")
    for line in blob.splitlines():
        line = line.strip()
        if line:
            return line
    return "unknown"


def pysam_version_string() -> str:
    """Return the installed pysam version, or 'unknown'."""
    try:
        import pysam  # type: ignore[import-not-found]
        return getattr(pysam, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    """Return the running Python version."""
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Stage 1 - GetPileupSummaries.
# ----------------------------------------------------------------------

def run_get_pileup_summaries(
    bam: Path,
    common_biallelic_vcf: Path,
    out_table: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> None:
    """Run gatk GetPileupSummaries on a single BAM.

    Raises subprocess.CalledProcessError on failure.
    Raises FileNotFoundError if gatk is not on the PATH.
    """
    if not shutil.which("gatk"):
        raise FileNotFoundError("gatk is not on the PATH")
    out_table.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "gatk", "--java-options", f"-Xmx{java_heap}",
        "GetPileupSummaries",
        "-I", str(bam),
        "-V", str(common_biallelic_vcf),
        "-L", str(common_biallelic_vcf),
        "-O", str(out_table),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


# ----------------------------------------------------------------------
# Stage 2 - CalculateContamination.
# ----------------------------------------------------------------------

def run_calculate_contamination(
    tumor_pileups: Path,
    normal_pileups: Path,
    out_table: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> None:
    """Run gatk CalculateContamination with the matched normal pileups."""
    if not shutil.which("gatk"):
        raise FileNotFoundError("gatk is not on the PATH")
    out_table.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "gatk", "--java-options", f"-Xmx{java_heap}",
        "CalculateContamination",
        "-I", str(tumor_pileups),
        "-matched", str(normal_pileups),
        "-O", str(out_table),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_contamination_table(table_path: Path) -> tuple[float, float]:
    """Read the contamination.table TSV and return (fraction, error).

    The CalculateContamination output has three columns:
        sample  contamination  error
    Returns (0.0, 0.0) on parse failure.
    """
    if not table_path.exists():
        return (0.0, 0.0)
    with table_path.open() as fh:
        header_seen: bool = False
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts: list[str] = line.split("\t")
            if not header_seen and parts and not _is_float(parts[1] if len(parts) > 1 else ""):
                header_seen = True
                continue
            if len(parts) >= 3:
                try:
                    return (float(parts[1]), float(parts[2]))
                except ValueError:
                    continue
    return (0.0, 0.0)


def _is_float(s: str) -> bool:
    """Return True if s parses as a float."""
    try:
        float(s)
        return True
    except ValueError:
        return False


# ----------------------------------------------------------------------
# Stage 3 - FilterMutectCalls.
# ----------------------------------------------------------------------

def run_filter_mutect_calls(
    reference: Path,
    unfiltered_vcf: Path,
    contamination_table: Path,
    out_vcf: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> None:
    """Run gatk FilterMutectCalls with the contamination table.

    Mutect2 writes a `.stats` file alongside the unfiltered VCF; we pass
    it explicitly so FilterMutectCalls knows where to find it.
    """
    if not shutil.which("gatk"):
        raise FileNotFoundError("gatk is not on the PATH")
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    stats_path: Path = Path(str(unfiltered_vcf) + ".stats")
    cmd: list[str] = [
        "gatk", "--java-options", f"-Xmx{java_heap}",
        "FilterMutectCalls",
        "-R", str(reference),
        "-V", str(unfiltered_vcf),
        "--contamination-table", str(contamination_table),
        "-O", str(out_vcf),
    ]
    if stats_path.exists():
        cmd.extend(["--stats", str(stats_path)])
    subprocess.run(cmd, check=True, capture_output=True, text=True)


# ----------------------------------------------------------------------
# Stage 4 - tally filter reasons.
# ----------------------------------------------------------------------

def tally_filters(vcf_path: Path) -> tuple[int, int, dict[str, int]]:
    """Return (n_total, n_pass, per_filter_counts) from a filtered VCF.

    A record's FILTER may be:
    - empty (we count as PASS)
    - ['PASS']
    - a comma-separated list of filter names (each contributes one count)

    Multi-filter variants increment each named filter's count.
    """
    try:
        import pysam  # type: ignore[import-not-found]
    except ImportError as exc:
        raise FileNotFoundError("pysam is not importable") from exc
    counts: Counter = Counter()
    n_total: int = 0
    n_pass: int = 0
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            n_total += 1
            filt: list[str] = list(rec.filter.keys())
            if not filt or filt == ["PASS"]:
                n_pass += 1
                counts["PASS"] += 1
                continue
            for f in filt:
                counts[f] += 1
    return n_total, n_pass, dict(counts)


# ----------------------------------------------------------------------
# Stage 5 - render Markdown.
# ----------------------------------------------------------------------

def render_filter_markdown(
    out_path: Path,
    n_total: int,
    n_pass: int,
    filter_counts: dict[str, int],
    contamination_fraction: float,
    contamination_error: float,
) -> None:
    """Write a Markdown report of the filter tally and contamination."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_filters: list[tuple[str, int]] = sorted(
        filter_counts.items(),
        key=lambda kv: kv[1],
        reverse=True,
    )
    pass_pct: float = (n_pass / n_total * 100.0) if n_total > 0 else 0.0
    lines: list[str] = [
        "# Exercise 2 - FilterMutectCalls Tally",
        "",
        "**Educational and research use only.** Not for clinical decisions.",
        "",
        f"- Total candidate variants in filtered VCF: {n_total}",
        f"- PASS variants: {n_pass} ({pass_pct:.1f}% of total)",
        f"- Estimated tumor contamination: {contamination_fraction:.4f} "
        f"(error {contamination_error:.4f})",
        "",
        "## Per-filter occurrences",
        "",
        "| FILTER reason          | Count |",
        "|------------------------|------:|",
    ]
    for name, count in sorted_filters:
        lines.append(f"| {name:<22s} | {count:>5d} |")
    lines.extend([
        "",
        "Note: a variant can match multiple filter reasons; the counts are not exclusive.",
        "PASS is mutually exclusive with all named filters.",
        "",
    ])
    out_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------
# Stage 6 - run-info JSON.
# ----------------------------------------------------------------------

def write_run_info(run_info: ExerciseRunInfo, out_path: Path) -> None:
    """Write the run-info JSON. Raises ValueError on empty run_date."""
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
    reference: Path,
    unfiltered_vcf: Path,
    tumor_bam: Path,
    normal_bam: Path,
    common_biallelic: Path,
    out_dir: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
    skip_if_missing: bool = True,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tumor_pileups: Path = out_dir / "tumor.pileups.table"
    normal_pileups: Path = out_dir / "normal.pileups.table"
    contam_table: Path = out_dir / "contamination.table"
    filtered_vcf: Path = out_dir / "filtered.vcf.gz"
    filter_md: Path = out_dir / "filter_tally.md"
    run_info_path: Path = out_dir / "run-info.json"

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        unfiltered_vcf=str(unfiltered_vcf),
        filtered_vcf=str(filtered_vcf),
        tumor_bam=str(tumor_bam),
        normal_bam=str(normal_bam),
        reference_fasta=str(reference),
        common_biallelic_vcf=str(common_biallelic),
        contamination_table=str(contam_table),
        gatk_version=get_tool_version("gatk", "--version"),
        pysam_version=pysam_version_string(),
        java_heap=java_heap,
        python_version=python_version_string(),
    )

    inputs: list[tuple[str, Path]] = [
        ("reference", reference),
        ("unfiltered_vcf", unfiltered_vcf),
        ("tumor_bam", tumor_bam),
        ("normal_bam", normal_bam),
        ("common_biallelic", common_biallelic),
    ]
    missing: list[str] = [label for label, path in inputs if not path.exists()]
    if missing:
        if skip_if_missing:
            info.skipped = True
            info.skip_reason = f"missing inputs: {', '.join(missing)}"
            write_run_info(info, run_info_path)
            print(f"[ex02] skipped: missing {missing}", file=sys.stderr)
            return run_info_path
        raise FileNotFoundError(f"missing inputs: {missing}")

    if not shutil.which("gatk"):
        if skip_if_missing:
            info.skipped = True
            info.skip_reason = "gatk is not on the PATH"
            write_run_info(info, run_info_path)
            print("[ex02] skipped: gatk not on PATH", file=sys.stderr)
            return run_info_path
        raise FileNotFoundError("gatk is not on the PATH")

    # Stage 1 - GetPileupSummaries on each BAM.
    try:
        run_get_pileup_summaries(tumor_bam, common_biallelic, tumor_pileups, java_heap=java_heap)
        run_get_pileup_summaries(normal_bam, common_biallelic, normal_pileups, java_heap=java_heap)
    except subprocess.CalledProcessError as exc:
        info.skipped = True
        info.skip_reason = f"GetPileupSummaries failed: {exc.returncode}"
        write_run_info(info, run_info_path)
        return run_info_path

    # Stage 2 - CalculateContamination.
    try:
        run_calculate_contamination(
            tumor_pileups=tumor_pileups,
            normal_pileups=normal_pileups,
            out_table=contam_table,
            java_heap=java_heap,
        )
    except subprocess.CalledProcessError as exc:
        info.skipped = True
        info.skip_reason = f"CalculateContamination failed: {exc.returncode}"
        write_run_info(info, run_info_path)
        return run_info_path

    contam_fraction, contam_error = parse_contamination_table(contam_table)
    info.contamination_fraction = contam_fraction
    info.contamination_error = contam_error

    # Stage 3 - FilterMutectCalls.
    try:
        run_filter_mutect_calls(
            reference=reference,
            unfiltered_vcf=unfiltered_vcf,
            contamination_table=contam_table,
            out_vcf=filtered_vcf,
            java_heap=java_heap,
        )
    except subprocess.CalledProcessError as exc:
        info.skipped = True
        info.skip_reason = f"FilterMutectCalls failed: {exc.returncode}"
        write_run_info(info, run_info_path)
        return run_info_path

    # Stage 4 - tally filters and render Markdown.
    try:
        n_total, n_pass, filter_counts = tally_filters(filtered_vcf)
        info.n_total = n_total
        info.n_pass = n_pass
        info.filter_counts = filter_counts
        render_filter_markdown(
            out_path=filter_md,
            n_total=n_total,
            n_pass=n_pass,
            filter_counts=filter_counts,
            contamination_fraction=contam_fraction,
            contamination_error=contam_error,
        )
    except Exception as exc:
        print(f"[ex02] tally / render failed: {exc}", file=sys.stderr)

    write_run_info(info, run_info_path)
    print(
        f"[ex02] filtered.vcf.gz: {info.n_pass} PASS of {info.n_total} total; "
        f"contamination {info.contamination_fraction:.4f}",
        file=sys.stderr,
    )
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 2 - FilterMutectCalls + CalculateContamination.",
    )
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--unfiltered-vcf", type=Path, required=True)
    parser.add_argument("--tumor-bam", type=Path, required=True)
    parser.add_argument("--normal-bam", type=Path, required=True)
    parser.add_argument("--common-biallelic", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--java-heap", type=str, default=DEFAULT_JAVA_HEAP)
    parser.add_argument("--no-skip-if-missing", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        reference=args.reference,
        unfiltered_vcf=args.unfiltered_vcf,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        common_biallelic=args.common_biallelic,
        out_dir=args.out_dir,
        java_heap=args.java_heap,
        skip_if_missing=not args.no_skip_if_missing,
    )
    print(f"[ex02] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
