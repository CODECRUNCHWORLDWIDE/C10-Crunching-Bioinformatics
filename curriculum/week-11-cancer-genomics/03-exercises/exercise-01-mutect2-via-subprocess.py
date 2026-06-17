"""
Exercise 1 - Run Mutect2 in tumor-normal mode via subprocess and parse the output.

Educational and research use only. Do not use any output of this pipeline
to guide patient care.

Goal: take a tumor BAM, a matched normal BAM, a reference FASTA, a
panel-of-normals VCF, a germline-resource VCF, and an intervals file;
call Mutect2 in tumor-normal mode via subprocess; parse the resulting
unfiltered VCF; emit a run-info JSON recording every parameter.

The exercise covers:

- Reading BAM @RG SM: headers to verify sample names before calling.
- Calling Mutect2 via subprocess.run with check=True and capture_output=True.
- Parsing the Mutect2 VCF with pysam.
- Counting candidate variants and the depth distribution.
- Writing a run-info JSON alongside the unfiltered VCF.

Estimated time: 90 minutes (30 min reading, 50 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-01-mutect2-via-subprocess.py
    --reference data/chr22_GRCh38.fasta
    --tumor-bam data/tumor_chr22.bam
    --normal-bam data/normal_chr22.bam
    --pon data/chr22_pon.vcf.gz
    --germline-resource data/chr22_gnomad.vcf.gz
    --intervals data/chr22_intervals.bed
    --out-dir results/ex01` runs end to end without errors when
  GATK is on the PATH; gracefully skips with a message otherwise.
- `results/ex01/unfiltered.vcf.gz` exists with at least one record (in non-skip mode).
- `results/ex01/run-info.json` exists with the GATK version, the
  sample names from the BAM headers, the input mode, the threads, the
  reference build, the PON / germline-resource paths and MD5s, and the run date.
- The sample-name verification step prevents a Mutect2 call with
  reversed tumor / normal flags from running.

Requirements:
    conda install -c bioconda gatk4=4.5.0.0 samtools=1.20 pysam=0.22.1 biopython=1.84

What you learn:
- The subprocess.run idiom with check=True and capture_output=True for GATK.
- Reading the BAM @RG SM: header via samtools and via pysam.
- The pysam VariantFile API for reading VCFs.
- The "verify sample names before you call" reproducibility pattern.

Tool versions assumed:
- Python 3.11+
- GATK 4.5.0.0 (CLI tool; subprocess.run it)
- samtools 1.20 (CLI; for the @RG SM verification step)
- pysam 0.22.1 (Python library)

References:
- Mutect2: Cibulskis et al. 2013, Nat Biotechnol 31:213
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/
- GATK Best Practices for somatic short-variant discovery:
  https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
- gnomAD: Karczewski et al. 2020, Nature 581:434
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_THREADS: int = 4
DEFAULT_JAVA_HEAP: str = "4g"
DEFAULT_REFERENCE_BUILD: str = "GRCh38"
DEFAULT_TUMOR_SAMPLE_NAME_HINT: str = "TUMOR"
DEFAULT_NORMAL_SAMPLE_NAME_HINT: str = "NORMAL"


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class ExerciseRunInfo:
    """Provenance metadata for the Exercise 1 run."""
    run_date: str = ""
    reference_fasta: str = ""
    reference_md5: str = ""
    reference_build: str = DEFAULT_REFERENCE_BUILD
    tumor_bam: str = ""
    tumor_bam_md5: str = ""
    tumor_sample_name: str = ""
    normal_bam: str = ""
    normal_bam_md5: str = ""
    normal_sample_name: str = ""
    pon_vcf: str = ""
    pon_vcf_md5: str = ""
    germline_vcf: str = ""
    germline_vcf_md5: str = ""
    intervals: str = ""
    out_vcf: str = ""
    gatk_version: str = ""
    samtools_version: str = ""
    pysam_version: str = ""
    threads: int = DEFAULT_THREADS
    java_heap: str = DEFAULT_JAVA_HEAP
    candidate_variants: int = 0
    candidate_snvs: int = 0
    candidate_indels: int = 0
    skipped: bool = False
    skip_reason: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

def md5_of_file(path: Path) -> str:
    """Return the hex MD5 of a file's bytes, or empty string on read failure."""
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
    """Return the first non-empty line from `<tool> <version_flag>`.

    GATK and samtools both write the version line to one of stdout / stderr;
    we capture both and pick the first informative line.
    """
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
    """Return the installed pysam version, or 'unknown' if unimportable."""
    try:
        import pysam  # type: ignore[import-not-found]
        return getattr(pysam, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    """Return the running Python version, e.g. '3.11.9'."""
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Stage 1 - read the BAM @RG SM: header.
# ----------------------------------------------------------------------

def get_sample_name_via_samtools(bam_path: Path) -> str:
    """Return the SM: value from the first @RG header line of a BAM.

    Uses `samtools view -H`. Raises ValueError if no @RG SM: is found.
    Raises FileNotFoundError if samtools is not on the PATH.
    """
    if not shutil.which("samtools"):
        raise FileNotFoundError("samtools is not on the PATH")
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        if not line.startswith("@RG"):
            continue
        for field in line.split("\t"):
            if field.startswith("SM:"):
                return field[len("SM:"):].strip()
    raise ValueError(f"No @RG SM: tag found in {bam_path}")


def get_sample_name_via_pysam(bam_path: Path) -> str:
    """Return the SM: value via pysam.AlignmentFile.header.

    Fallback for when samtools is not on the PATH.
    """
    try:
        import pysam  # type: ignore[import-not-found]
    except ImportError as exc:
        raise FileNotFoundError("pysam is not importable") from exc
    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        header_dict: dict[str, Any] = dict(af.header)
        for rg in header_dict.get("RG", []):
            sm: str = rg.get("SM", "")
            if sm:
                return sm
    raise ValueError(f"No @RG SM: tag found in {bam_path}")


def get_sample_name(bam_path: Path) -> str:
    """Return the SM: value via samtools if available, else pysam."""
    if shutil.which("samtools"):
        return get_sample_name_via_samtools(bam_path)
    return get_sample_name_via_pysam(bam_path)


# ----------------------------------------------------------------------
# Stage 2 - sanity-check the input set.
# ----------------------------------------------------------------------

def verify_inputs(
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    pon: Path,
    germline: Path,
    intervals: Path,
) -> tuple[str, str]:
    """Verify every input exists and return the tumor / normal sample names.

    Raises FileNotFoundError or ValueError on any problem.
    """
    for label, path in [
        ("reference", reference),
        ("tumor_bam", tumor_bam),
        ("normal_bam", normal_bam),
        ("pon", pon),
        ("germline", germline),
        ("intervals", intervals),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} input not found: {path}")
    fai: Path = Path(str(reference) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(
            f"Reference index not found: {fai}. Run `samtools faidx {reference}`."
        )
    tumor_sample: str = get_sample_name(tumor_bam)
    normal_sample: str = get_sample_name(normal_bam)
    if not tumor_sample:
        raise ValueError(f"Empty tumor sample name from {tumor_bam}")
    if not normal_sample:
        raise ValueError(f"Empty normal sample name from {normal_bam}")
    if tumor_sample == normal_sample:
        raise ValueError(
            f"Tumor and normal samples share the same name '{tumor_sample}'. "
            f"This usually means the BAMs were mis-tagged at library prep. Aborting."
        )
    return tumor_sample, normal_sample


# ----------------------------------------------------------------------
# Stage 3 - run Mutect2.
# ----------------------------------------------------------------------

def run_mutect2(
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    tumor_sample: str,
    normal_sample: str,
    pon_vcf: Path,
    germline_vcf: Path,
    intervals: Path,
    out_vcf: Path,
    threads: int = DEFAULT_THREADS,
    java_heap: str = DEFAULT_JAVA_HEAP,
) -> None:
    """Run GATK Mutect2 in tumor-normal mode.

    Raises subprocess.CalledProcessError if Mutect2 fails.
    Raises FileNotFoundError if gatk is not on the PATH.
    """
    if not shutil.which("gatk"):
        raise FileNotFoundError("gatk is not on the PATH")
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "gatk", "--java-options", f"-Xmx{java_heap}",
        "Mutect2",
        "-R", str(reference),
        "-I", str(tumor_bam),
        "-I", str(normal_bam),
        "-tumor", tumor_sample,
        "-normal", normal_sample,
        "--panel-of-normals", str(pon_vcf),
        "--germline-resource", str(germline_vcf),
        "-L", str(intervals),
        "-O", str(out_vcf),
        "--native-pair-hmm-threads", str(threads),
    ]
    print(f"[ex01] Mutect2 cmd: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    log_path: Path = out_vcf.parent / "mutect2.log"
    with log_path.open("w") as fh:
        fh.write("# stdout\n")
        fh.write(result.stdout or "")
        fh.write("\n# stderr\n")
        fh.write(result.stderr or "")


# ----------------------------------------------------------------------
# Stage 4 - parse the unfiltered VCF.
# ----------------------------------------------------------------------

@dataclass
class VariantCounts:
    """Summary counts from an unfiltered VCF."""
    total: int = 0
    snvs: int = 0
    indels: int = 0
    by_filter: dict[str, int] | None = None


def parse_unfiltered_vcf(vcf_path: Path) -> VariantCounts:
    """Parse the Mutect2 unfiltered VCF and return summary counts.

    Counts total records, SNVs (REF len == 1, all ALT len == 1), indels (otherwise),
    and per-FILTER occurrences (the unfiltered VCF should have all PASS or empty).
    """
    try:
        import pysam  # type: ignore[import-not-found]
    except ImportError as exc:
        raise FileNotFoundError("pysam is not importable") from exc
    counts = VariantCounts(by_filter={})
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            counts.total += 1
            alts: list[str] = list(rec.alts or [])
            is_snv: bool = (
                len(rec.ref) == 1
                and bool(alts)
                and all(len(a) == 1 for a in alts)
            )
            if is_snv:
                counts.snvs += 1
            else:
                counts.indels += 1
            filt: list[str] = list(rec.filter.keys())
            if not filt:
                counts.by_filter["PASS"] = counts.by_filter.get("PASS", 0) + 1
            else:
                for f in filt:
                    counts.by_filter[f] = counts.by_filter.get(f, 0) + 1
    return counts


# ----------------------------------------------------------------------
# Stage 5 - run-info JSON.
# ----------------------------------------------------------------------

def write_run_info(run_info: ExerciseRunInfo, out_path: Path) -> None:
    """Write the run-info JSON. Raises ValueError on empty required fields."""
    if not run_info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    if not run_info.tumor_sample_name and not run_info.skipped:
        raise ValueError("run_info.tumor_sample_name is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    pon: Path,
    germline: Path,
    intervals: Path,
    out_dir: Path,
    threads: int = DEFAULT_THREADS,
    java_heap: str = DEFAULT_JAVA_HEAP,
    skip_if_no_gatk: bool = True,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    unfiltered_vcf: Path = out_dir / "unfiltered.vcf.gz"
    run_info_path: Path = out_dir / "run-info.json"
    out_dir.mkdir(parents=True, exist_ok=True)

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        reference_fasta=str(reference),
        reference_md5=md5_of_file(reference),
        tumor_bam=str(tumor_bam),
        tumor_bam_md5=md5_of_file(tumor_bam),
        normal_bam=str(normal_bam),
        normal_bam_md5=md5_of_file(normal_bam),
        pon_vcf=str(pon),
        pon_vcf_md5=md5_of_file(pon),
        germline_vcf=str(germline),
        germline_vcf_md5=md5_of_file(germline),
        intervals=str(intervals),
        out_vcf=str(unfiltered_vcf),
        threads=threads,
        java_heap=java_heap,
        gatk_version=get_tool_version("gatk", "--version"),
        samtools_version=get_tool_version("samtools", "--version"),
        pysam_version=pysam_version_string(),
        python_version=python_version_string(),
    )

    # Stage 1+2 - verify inputs and read sample names.
    try:
        tumor_sample, normal_sample = verify_inputs(
            reference=reference,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            pon=pon,
            germline=germline,
            intervals=intervals,
        )
        info.tumor_sample_name = tumor_sample
        info.normal_sample_name = normal_sample
    except (FileNotFoundError, ValueError) as exc:
        if skip_if_no_gatk:
            info.skipped = True
            info.skip_reason = f"input verification failed: {exc}"
            write_run_info(info, run_info_path)
            print(f"[ex01] skipped: {exc}", file=sys.stderr)
            return run_info_path
        raise

    # Stage 3 - run Mutect2.
    if not shutil.which("gatk"):
        if skip_if_no_gatk:
            info.skipped = True
            info.skip_reason = "gatk is not on the PATH"
            write_run_info(info, run_info_path)
            print("[ex01] skipped: gatk not on PATH", file=sys.stderr)
            return run_info_path
        raise FileNotFoundError("gatk is not on the PATH")

    try:
        run_mutect2(
            reference=reference,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            tumor_sample=info.tumor_sample_name,
            normal_sample=info.normal_sample_name,
            pon_vcf=pon,
            germline_vcf=germline,
            intervals=intervals,
            out_vcf=unfiltered_vcf,
            threads=threads,
            java_heap=java_heap,
        )
    except subprocess.CalledProcessError as exc:
        info.skipped = True
        info.skip_reason = f"Mutect2 failed: returncode={exc.returncode}"
        write_run_info(info, run_info_path)
        print(f"[ex01] Mutect2 failed: {exc}", file=sys.stderr)
        return run_info_path

    # Stage 4 - parse the VCF.
    if unfiltered_vcf.exists():
        try:
            counts = parse_unfiltered_vcf(unfiltered_vcf)
            info.candidate_variants = counts.total
            info.candidate_snvs = counts.snvs
            info.candidate_indels = counts.indels
        except Exception as exc:
            print(f"[ex01] VCF parse failed: {exc}", file=sys.stderr)

    write_run_info(info, run_info_path)
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 1 - Mutect2 via subprocess (educational/research use only).",
    )
    parser.add_argument(
        "--reference", type=Path, required=True,
        help="Reference FASTA path. Must have a .fai index.",
    )
    parser.add_argument(
        "--tumor-bam", type=Path, required=True,
        help="Tumor BAM path. Must have a .bai index.",
    )
    parser.add_argument(
        "--normal-bam", type=Path, required=True,
        help="Matched-normal BAM path. Must have a .bai index.",
    )
    parser.add_argument(
        "--pon", type=Path, required=True,
        help="Panel-of-normals VCF path.",
    )
    parser.add_argument(
        "--germline-resource", type=Path, required=True,
        help="gnomAD allele-frequency VCF path.",
    )
    parser.add_argument(
        "--intervals", type=Path, required=True,
        help="Intervals BED or list path.",
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True,
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--threads", type=int, default=DEFAULT_THREADS,
        help=f"Mutect2 --native-pair-hmm-threads. Default: {DEFAULT_THREADS}.",
    )
    parser.add_argument(
        "--java-heap", type=str, default=DEFAULT_JAVA_HEAP,
        help=f"JVM -Xmx value. Default: {DEFAULT_JAVA_HEAP}.",
    )
    parser.add_argument(
        "--no-skip-if-missing",
        action="store_true",
        help="Raise an error if GATK or any input is missing (default: graceful skip).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        reference=args.reference,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        pon=args.pon,
        germline=args.germline_resource,
        intervals=args.intervals,
        out_dir=args.out_dir,
        threads=args.threads,
        java_heap=args.java_heap,
        skip_if_no_gatk=not args.no_skip_if_missing,
    )
    print(f"[ex01] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
