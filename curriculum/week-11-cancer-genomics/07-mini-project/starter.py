"""
Mini-Project Starter - Somatic variant calling pipeline (Mutect2 +
FilterMutectCalls + SigProfilerAssignment).

Educational and research use only. Do not use any output of this pipeline
to guide patient care.

This is a skeleton you complete in the mini-project. The stages are
sketched out as functions; each has a TODO marking where you fill in
the call to the Exercise 1 / 2 / 3 helpers (or reimplement them here).

The expected end state is a single function:

    run_cancer_pipeline(
        tumor_bam: Path,
        normal_bam: Path,
        reference: Path,
        pon: Path,
        germline: Path,
        common_biallelic: Path,
        intervals: Path,
        out_dir: Path,
    ) -> Path

that returns the path to the final run-info.json after running all seven
stages.

The pipeline stages:

    Stage 1: Validate inputs and verify BAM @RG SM: sample names.
    Stage 2: Run Mutect2 in tumor-normal mode.
    Stage 3: Run GetPileupSummaries on both BAMs.
             Run CalculateContamination with the matched normal pileups.
    Stage 4: Run FilterMutectCalls with --contamination-table.
    Stage 5: Build the 96-class trinucleotide spectrum from PASS SNVs.
    Stage 6: Run SigProfilerAssignment against COSMIC v3.3 SBS catalog.
    Stage 7: Render filter_tally.md, signature_summary.md, qc_report.md,
             and run-info.json.

Pin every tool version in the run-info; the JSON is the deliverable as
much as the VCFs are.

References:
    - Mutect2: Cibulskis et al. 2013, Nat Biotechnol 31:213
    - Strelka2: Kim et al. 2018, Nat Methods 15:591
    - COSMIC signatures: Alexandrov et al. 2020, Nature 578:94
    - GATK Best Practices for somatic short-variant discovery:
      https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_THREADS: int = 4
DEFAULT_JAVA_HEAP: str = "4g"
DEFAULT_GENOME_BUILD: str = "GRCh38"
DEFAULT_COSMIC_VERSION: str = "3.3"
DEFAULT_TUMOR_SAMPLE_NAME: str = "TUMOR"


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class PipelineRunInfo:
    """Provenance for the full mini-project pipeline run."""
    run_date: str = ""
    tumor_bam: str = ""
    tumor_sample_name: str = ""
    normal_bam: str = ""
    normal_sample_name: str = ""
    reference_fasta: str = ""
    reference_build: str = DEFAULT_GENOME_BUILD
    pon_vcf: str = ""
    germline_vcf: str = ""
    common_biallelic_vcf: str = ""
    intervals: str = ""
    threads: int = DEFAULT_THREADS
    java_heap: str = DEFAULT_JAVA_HEAP
    gatk_version: str = ""
    samtools_version: str = ""
    bcftools_version: str = ""
    pysam_version: str = ""
    sigprofiler_version: str = ""
    cosmic_version: str = DEFAULT_COSMIC_VERSION
    n_candidate_variants: int = 0
    n_candidate_snvs: int = 0
    n_candidate_indels: int = 0
    n_pass_variants: int = 0
    n_pass_snvs_used_for_spectrum: int = 0
    filter_counts: dict[str, int] = field(default_factory=dict)
    contamination_fraction: float = 0.0
    contamination_error: float = 0.0
    top_signatures: list[dict[str, Any]] = field(default_factory=list)
    cosine_similarity: float = 0.0
    educational_use_only: bool = True
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

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


def sigprofiler_version_string() -> str:
    """Return the installed SigProfilerAssignment version, or 'unknown'."""
    try:
        import SigProfilerAssignment  # type: ignore[import-not-found]
        return getattr(SigProfilerAssignment, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    """Return the running Python version."""
    return ".".join(str(part) for part in sys.version_info[:3])


def get_sample_name_via_samtools(bam_path: Path) -> str:
    """Return the SM: value from the first @RG header line of a BAM."""
    if not shutil.which("samtools"):
        raise FileNotFoundError("samtools is not on the PATH")
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        check=True, capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        if not line.startswith("@RG"):
            continue
        for field_str in line.split("\t"):
            if field_str.startswith("SM:"):
                return field_str[len("SM:"):].strip()
    raise ValueError(f"No @RG SM: tag found in {bam_path}")


# ----------------------------------------------------------------------
# Stage 1 - validate inputs.
# ----------------------------------------------------------------------

def stage_1_validate(
    tumor_bam: Path,
    normal_bam: Path,
    reference: Path,
    pon: Path,
    germline: Path,
    common_biallelic: Path,
    intervals: Path,
) -> tuple[str, str]:
    """Validate every input and return (tumor_sample, normal_sample)."""
    for label, path in [
        ("tumor_bam", tumor_bam), ("normal_bam", normal_bam),
        ("reference", reference), ("pon", pon),
        ("germline", germline), ("common_biallelic", common_biallelic),
        ("intervals", intervals),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
    fai: Path = Path(str(reference) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(f"Reference index not found: {fai}")
    tumor_sample: str = get_sample_name_via_samtools(tumor_bam)
    normal_sample: str = get_sample_name_via_samtools(normal_bam)
    if tumor_sample == normal_sample:
        raise ValueError(
            f"Tumor and normal share sample name '{tumor_sample}'. Aborting."
        )
    return tumor_sample, normal_sample


# ----------------------------------------------------------------------
# Stage 2 - Mutect2.
# ----------------------------------------------------------------------

def stage_2_mutect2(
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    tumor_sample: str,
    normal_sample: str,
    pon: Path,
    germline: Path,
    intervals: Path,
    out_vcf: Path,
    threads: int = DEFAULT_THREADS,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> None:
    """Run GATK Mutect2 in tumor-normal mode.

    TODO (mini-project): reuse the Exercise 1 `run_mutect2` helper or
    reimplement the subprocess call here.
    """
    # TODO: implement the subprocess.run() call to `gatk Mutect2`.
    # See exercises/exercise-01-mutect2-via-subprocess.py for the canonical pattern.
    raise NotImplementedError("stage_2_mutect2: complete in the mini-project")


# ----------------------------------------------------------------------
# Stage 3 - contamination.
# ----------------------------------------------------------------------

def stage_3_contamination(
    tumor_bam: Path,
    normal_bam: Path,
    common_biallelic: Path,
    out_dir: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> tuple[Path, float, float]:
    """Run GetPileupSummaries on each BAM, then CalculateContamination.

    TODO (mini-project): reuse the Exercise 2 helpers `run_get_pileup_summaries`,
    `run_calculate_contamination`, and `parse_contamination_table`.

    Returns (contamination_table_path, contamination_fraction, contamination_error).
    """
    # TODO: implement the two GetPileupSummaries calls + the matched-normal
    # CalculateContamination call.
    raise NotImplementedError("stage_3_contamination: complete in the mini-project")


# ----------------------------------------------------------------------
# Stage 4 - FilterMutectCalls.
# ----------------------------------------------------------------------

def stage_4_filter(
    reference: Path,
    unfiltered_vcf: Path,
    contamination_table: Path,
    out_vcf: Path,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> dict[str, int]:
    """Run FilterMutectCalls and return the FILTER tally.

    TODO (mini-project): reuse the Exercise 2 helpers `run_filter_mutect_calls`
    and `tally_filters`.

    Returns a dict keyed on FILTER reason (including 'PASS').
    """
    # TODO: implement the FilterMutectCalls call + the FILTER tally read.
    raise NotImplementedError("stage_4_filter: complete in the mini-project")


# ----------------------------------------------------------------------
# Stage 5 - 96-class spectrum.
# ----------------------------------------------------------------------

def stage_5_spectrum(
    filtered_vcf: Path,
    reference: Path,
    sample_name: str,
    out_tsv: Path,
) -> tuple[int, int]:
    """Build the 96-class spectrum from PASS SNVs.

    TODO (mini-project): reuse the Exercise 3 helpers `build_96_class_spectrum`
    and `write_spectrum_tsv`.

    Returns (n_pass_records, n_snvs_used).
    """
    # TODO: implement the spectrum build + TSV write.
    raise NotImplementedError("stage_5_spectrum: complete in the mini-project")


# ----------------------------------------------------------------------
# Stage 6 - SigProfilerAssignment.
# ----------------------------------------------------------------------

def stage_6_signatures(
    spectrum_tsv: Path,
    out_dir: Path,
    genome_build: str = DEFAULT_GENOME_BUILD,
    cosmic_version: str = DEFAULT_COSMIC_VERSION,
) -> tuple[list[dict[str, Any]], float]:
    """Run SigProfilerAssignment.cosmic_fit() and parse the result.

    TODO (mini-project): reuse the Exercise 3 helpers `run_sigprofiler_assignment`,
    `parse_activities_tsv`, `find_stats_file`, and `parse_sample_stats`.

    Returns (top_signatures_list, cosine_similarity).
    """
    # TODO: implement the SigProfilerAssignment call + the activities + stats parse.
    raise NotImplementedError("stage_6_signatures: complete in the mini-project")


# ----------------------------------------------------------------------
# Stage 7 - render reports.
# ----------------------------------------------------------------------

def stage_7_render_reports(
    info: PipelineRunInfo,
    filter_tally_md: Path,
    signature_summary_md: Path,
    qc_report_md: Path,
    run_info_json: Path,
) -> None:
    """Render the Markdown reports and the run-info JSON.

    TODO (mini-project): reuse the Exercise 2 `render_filter_markdown` and
    the Exercise 3 `render_signature_markdown`; add a top-level
    qc_report.md that combines both.
    """
    # TODO: implement the three Markdown renders + the JSON write.
    # The educational-use disclaimer is mandatory at the top of each Markdown.
    raise NotImplementedError("stage_7_render_reports: complete in the mini-project")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_cancer_pipeline(
    tumor_bam: Path,
    normal_bam: Path,
    reference: Path,
    pon: Path,
    germline: Path,
    common_biallelic: Path,
    intervals: Path,
    out_dir: Path,
    threads: int = DEFAULT_THREADS,
    java_heap: str = DEFAULT_JAVA_HEAP,
    cosmic_version: str = DEFAULT_COSMIC_VERSION,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> Path:
    """Run the full mini-project pipeline. Returns the path to run-info.json.

    The orchestrator is intentionally skeletal; each stage is filled in
    using the Exercise 1 / 2 / 3 helpers.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    unfiltered_vcf: Path = out_dir / "unfiltered.vcf.gz"
    filtered_vcf: Path = out_dir / "filtered.vcf.gz"
    spectrum_tsv: Path = out_dir / "spectrum_96.tsv"
    sigprofiler_out: Path = out_dir / "sigprofiler_out"
    filter_tally_md: Path = out_dir / "filter_tally.md"
    signature_summary_md: Path = out_dir / "signature_summary.md"
    qc_report_md: Path = out_dir / "qc_report.md"
    run_info_json: Path = out_dir / "run-info.json"

    info = PipelineRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        tumor_bam=str(tumor_bam),
        normal_bam=str(normal_bam),
        reference_fasta=str(reference),
        reference_build=genome_build,
        pon_vcf=str(pon),
        germline_vcf=str(germline),
        common_biallelic_vcf=str(common_biallelic),
        intervals=str(intervals),
        threads=threads,
        java_heap=java_heap,
        gatk_version=get_tool_version("gatk", "--version"),
        samtools_version=get_tool_version("samtools", "--version"),
        bcftools_version=get_tool_version("bcftools", "--version"),
        pysam_version=pysam_version_string(),
        sigprofiler_version=sigprofiler_version_string(),
        cosmic_version=cosmic_version,
        python_version=python_version_string(),
    )

    # Stage 1.
    tumor_sample, normal_sample = stage_1_validate(
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        reference=reference,
        pon=pon,
        germline=germline,
        common_biallelic=common_biallelic,
        intervals=intervals,
    )
    info.tumor_sample_name = tumor_sample
    info.normal_sample_name = normal_sample

    # Stages 2-6: each TODO is the learner's work.
    stage_2_mutect2(
        reference=reference,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        tumor_sample=tumor_sample,
        normal_sample=normal_sample,
        pon=pon,
        germline=germline,
        intervals=intervals,
        out_vcf=unfiltered_vcf,
        threads=threads,
        java_heap=java_heap,
    )
    # info.n_candidate_variants = ...  fill from VCF parse

    contam_table, contam_fraction, contam_error = stage_3_contamination(
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        common_biallelic=common_biallelic,
        out_dir=out_dir,
        java_heap=java_heap,
    )
    info.contamination_fraction = contam_fraction
    info.contamination_error = contam_error

    filter_counts: dict[str, int] = stage_4_filter(
        reference=reference,
        unfiltered_vcf=unfiltered_vcf,
        contamination_table=contam_table,
        out_vcf=filtered_vcf,
        java_heap=java_heap,
    )
    info.filter_counts = filter_counts
    info.n_pass_variants = filter_counts.get("PASS", 0)

    n_pass_records, n_snvs_used = stage_5_spectrum(
        filtered_vcf=filtered_vcf,
        reference=reference,
        sample_name=tumor_sample,
        out_tsv=spectrum_tsv,
    )
    info.n_pass_snvs_used_for_spectrum = n_snvs_used

    top_signatures, cosine = stage_6_signatures(
        spectrum_tsv=spectrum_tsv,
        out_dir=sigprofiler_out,
        genome_build=genome_build,
        cosmic_version=cosmic_version,
    )
    info.top_signatures = top_signatures
    info.cosine_similarity = cosine

    # Stage 7.
    stage_7_render_reports(
        info=info,
        filter_tally_md=filter_tally_md,
        signature_summary_md=signature_summary_md,
        qc_report_md=qc_report_md,
        run_info_json=run_info_json,
    )
    return run_info_json


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Mini-Project Starter - somatic variant pipeline.",
    )
    parser.add_argument("--tumor-bam", type=Path, required=True)
    parser.add_argument("--normal-bam", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--pon", type=Path, required=True)
    parser.add_argument("--germline-resource", type=Path, required=True)
    parser.add_argument("--common-biallelic", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--java-heap", type=str, default=DEFAULT_JAVA_HEAP)
    parser.add_argument("--cosmic-version", type=str, default=DEFAULT_COSMIC_VERSION)
    parser.add_argument("--genome-build", type=str, default=DEFAULT_GENOME_BUILD)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    out_path: Path = run_cancer_pipeline(
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        reference=args.reference,
        pon=args.pon,
        germline=args.germline_resource,
        common_biallelic=args.common_biallelic,
        intervals=args.intervals,
        out_dir=args.out_dir,
        threads=args.threads,
        java_heap=args.java_heap,
        cosmic_version=args.cosmic_version,
        genome_build=args.genome_build,
    )
    print(f"[mini-project] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
