"""
Exercise 1 - Call variants on a sorted, indexed BAM end to end with bcftools.

Goal: take a sorted, indexed, duplicate-marked BAM file (the output of
Week 5 Exercise 1 or the Week 5 mini-project), run the canonical bcftools
mpileup + bcftools call pipeline, hard-filter the result, normalize indel
representation, and verify the output VCF with pysam.

Estimated time: 50 minutes (most spent on toolchain debugging the first
time you run bcftools on your machine).

Acceptance criteria:
- `python exercise-01-call-with-bcftools.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The BAM at exercises/aln/lambda.sorted.bam from Week 5 Exercise 1 is
  read, NOT re-built; this exercise assumes Week 5 is done.
- After the script runs, calls/lambda.raw.vcf.gz, calls/lambda.filtered.vcf.gz,
  and calls/lambda.norm.vcf.gz all exist on disk, each with a .tbi index.
- You implemented five functions: `index_reference`, `call_variants`,
  `filter_variants`, `normalize_variants`, and `summarize_vcf`.

Requirements:
    conda install -c bioconda bcftools=1.19 samtools=1.19 pysam=0.22

What you learn:
- The full bcftools mpileup -> bcftools call -> bcftools filter
  -> bcftools norm pipeline, invoked from Python via subprocess.
- How pysam.VariantFile reads a sorted+indexed VCF and lets you iterate
  over records.
- The difference between QUAL, INFO/DP, INFO/MQ, and INFO/SP and why
  each one matters as a hard-filter axis.
- Why bcftools norm is required before any cross-VCF comparison.

TO COMPLETE: implement the five functions below. Run the file; all
assertions must pass.

Reference: bacteriophage lambda is the textbook tiny genome - same
reference as Week 5 Exercise 1. The full bcftools pipeline on lambda
takes ~5 seconds end to end, so you can iterate fast.

Tool versions assumed:
- bcftools 1.19
- samtools 1.19
- pysam 0.22
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any


# Reference accession. NC_001416.1 is the bacteriophage lambda complete
# genome, 48,502 bp - small enough that bcftools runs in seconds.
EXPECTED_CONTIG = "NC_001416.1"
EXPECTED_REF_LENGTH = 48502


# The input BAM. We assume Week 5 Exercise 1 has been completed and the
# sorted, indexed BAM is at this path. If not, the script errors out
# with a helpful message.
DEFAULT_INPUT_BAM = (
    Path(__file__).parent.parent.parent
    / "week-05-read-alignment"
    / "exercises"
    / "aln"
    / "lambda.sorted.bam"
)

# Lambda is a haploid virus (essentially a single linear DNA molecule),
# so we set --ploidy 1 just like for bacteria.
PLOIDY = 1


def index_reference(ref_path: Path) -> None:
    """Run `samtools faidx` on `ref_path` if not already indexed.

    bcftools mpileup requires the .fai index to do random-access lookups
    on the reference FASTA. Without it, mpileup errors out with an
    unhelpful message.

    Hint: check for `{ref_path}.fai`; if missing, run
    `subprocess.run(["samtools", "faidx", str(ref_path)], check=True)`.
    """
    fai_marker = ref_path.parent / f"{ref_path.name}.fai"

    # TODO: if fai_marker does not exist, run `samtools faidx ref_path`.
    raise NotImplementedError("Run samtools faidx on the reference")


def call_variants(
    ref_path: Path, bam_path: Path, out_vcf: Path, ploidy: int = 1
) -> Path:
    """Run `bcftools mpileup | bcftools call` to emit a raw VCF.

    This invokes the canonical two-step pipeline:

        bcftools mpileup -Ou -f ref.fa --max-depth 250 -a 'AD,DP,SP' \\
            aln.bam \\
        | bcftools call -m -v --ploidy <ploidy> -Oz -o out.vcf.gz -

    Then index the resulting VCF with `bcftools index -t`.

    Args:
        ref_path: path to the reference FASTA (with .fai index alongside).
        bam_path: path to a sorted, indexed BAM.
        out_vcf: where to write the bgzipped VCF.
        ploidy: 1 for haploid (lambda, bacteria), 2 for diploid (human germline).

    Returns the path to the indexed VCF.

    Hint: use `subprocess.run(cmd, shell=True, check=True)` with the
    full pipe string. shell=True is appropriate here because the pipe
    is part of the command.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    # TODO: build the shell command. Use bcftools mpileup -Ou -f ref.fa
    #       --max-depth 250 -a 'AD,DP,SP' bam | bcftools call -m -v
    #       --ploidy <ploidy> -Oz -o out.vcf.gz -
    # TODO: run it with subprocess.run(cmd, shell=True, check=True).
    # TODO: run bcftools index -t out_vcf to produce the .tbi.
    raise NotImplementedError("Run bcftools mpileup | bcftools call and index")


def filter_variants(in_vcf: Path, out_vcf: Path) -> Path:
    """Apply a simplified hard-filter expression to a raw VCF.

    Uses the bacterial-friendly recipe from Lecture 2 Section 4.3:

        bcftools filter -Oz -o out.vcf.gz \\
            -s LowQual \\
            -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' \\
            in.vcf.gz

    Variants matching the expression are tagged 'LowQual' in the FILTER
    column; the rest are tagged PASS.

    Args:
        in_vcf: path to the raw bgzipped VCF (output of call_variants).
        out_vcf: where to write the filtered bgzipped VCF.

    Returns the path to the indexed filtered VCF.

    Hint: use subprocess.run with a list (not shell=True). The expression
    must be a single argument; quoting is handled by subprocess.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    expr = "QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60"

    # TODO: run bcftools filter -Oz -o out_vcf -s LowQual -e <expr> in_vcf.
    # TODO: run bcftools index -t out_vcf to produce the .tbi.
    raise NotImplementedError("Run bcftools filter and index")


def normalize_variants(
    ref_path: Path, in_vcf: Path, out_vcf: Path
) -> Path:
    """Left-align indels and split multiallelic records.

    Uses the canonical normalization recipe:

        bcftools norm -f ref.fa -m -any -Oz -o out.vcf.gz in.vcf.gz

    -f ref.fa is required for left-alignment.
    -m -any splits multiallelic records into one per ALT allele.

    Args:
        ref_path: path to the reference FASTA (with .fai index).
        in_vcf: path to the filtered VCF.
        out_vcf: where to write the normalized VCF.

    Returns the path to the indexed normalized VCF.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    # TODO: run bcftools norm -f ref -m -any -Oz -o out_vcf in_vcf.
    # TODO: run bcftools index -t out_vcf to produce the .tbi.
    raise NotImplementedError("Run bcftools norm and index")


def summarize_vcf(vcf_path: Path) -> dict:
    """Open a VCF with pysam.VariantFile and return summary statistics.

    Returns a dict with keys:
        "n_total":      int, total variant records
        "n_snp":        int, records where REF and all ALTs are length 1
        "n_indel":      int, records where REF or any ALT length is > 1
        "n_pass":       int, records with FILTER == PASS (or empty FILTER)
        "n_low_qual":   int, records with FILTER containing 'LowQual'
        "mean_qual":    float, mean of QUAL across all records
        "mean_dp":      float, mean of INFO/DP across all records
        "samples":      list of sample names in the VCF header

    Hint: use `pysam.VariantFile(str(vcf_path))`. Iterate with
    `for rec in vf:`. Use `rec.ref`, `rec.alts`, `rec.qual`,
    `rec.filter.keys()`, `rec.info["DP"]`, `vf.header.samples`.

    A 'PASS' record has either FILTER == [] (no filters applied yet)
    or FILTER == ['PASS'] (filtered and passed).
    """
    import pysam

    # TODO: open the VCF with pysam.VariantFile.
    # TODO: iterate over records, tallying the counts above.
    # TODO: compute mean QUAL and mean INFO/DP.
    # TODO: return the dict.
    raise NotImplementedError("Open VCF with pysam and tally fields")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-01-call-with-bcftools.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    here = Path(__file__).parent
    ref_dir = here / "ref"
    aln_dir = here / "aln"
    calls_dir = here / "calls"

    # The reference and BAM are expected to live at the Week 5 exercise
    # paths. If they have been moved, edit DEFAULT_INPUT_BAM and the
    # ref_path below.
    ref_path = (
        here.parent.parent
        / "week-05-read-alignment"
        / "exercises"
        / "ref"
        / "lambda.fa"
    )
    bam_path = DEFAULT_INPUT_BAM

    raw_vcf = calls_dir / "lambda.raw.vcf.gz"
    filt_vcf = calls_dir / "lambda.filtered.vcf.gz"
    norm_vcf = calls_dir / "lambda.norm.vcf.gz"

    # Toolchain sanity check.
    for tool in ("bcftools", "samtools"):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"[exercise-01] {tool!r} is not on PATH. Install it:\n"
                f"    conda install -c bioconda {tool}=1.19\n"
                f"and re-run."
            )

    # The reference and BAM must already exist (from Week 5 Exercise 1).
    if not ref_path.exists():
        raise SystemExit(
            f"[exercise-01] Reference FASTA not found at:\n"
            f"    {ref_path}\n"
            f"Run Week 5 Exercise 1 first to fetch the lambda reference."
        )
    if not bam_path.exists():
        raise SystemExit(
            f"[exercise-01] Input BAM not found at:\n"
            f"    {bam_path}\n"
            f"Run Week 5 Exercise 1 first to produce the sorted, indexed BAM."
        )

    print("[exercise-01] Step 1: indexing reference (samtools faidx) ...")
    index_reference(ref_path)
    assert (ref_path.parent / f"{ref_path.name}.fai").exists(), (
        "samtools faidx did not produce .fai"
    )

    print("[exercise-01] Step 2: calling variants (bcftools mpileup + call) ...")
    call_variants(ref_path, bam_path, raw_vcf, ploidy=PLOIDY)
    assert raw_vcf.exists(), f"raw VCF not at {raw_vcf}"
    assert raw_vcf.with_suffix(".gz.tbi").exists(), (
        f"VCF index .tbi not at {raw_vcf}.tbi - did bcftools index run?"
    )

    print("[exercise-01] Step 3: hard-filtering (bcftools filter) ...")
    filter_variants(raw_vcf, filt_vcf)
    assert filt_vcf.exists(), f"filtered VCF not at {filt_vcf}"
    assert filt_vcf.with_suffix(".gz.tbi").exists(), (
        f"filtered VCF index .tbi not at {filt_vcf}.tbi"
    )

    print("[exercise-01] Step 4: normalizing indel representation (bcftools norm) ...")
    normalize_variants(ref_path, filt_vcf, norm_vcf)
    assert norm_vcf.exists(), f"normalized VCF not at {norm_vcf}"
    assert norm_vcf.with_suffix(".gz.tbi").exists(), (
        f"normalized VCF index .tbi not at {norm_vcf}.tbi"
    )

    print("[exercise-01] Step 5: summarizing the normalized VCF ...")
    summary = summarize_vcf(norm_vcf)

    # Field presence.
    for key in (
        "n_total", "n_snp", "n_indel", "n_pass", "n_low_qual",
        "mean_qual", "mean_dp", "samples",
    ):
        assert key in summary, f"summary missing field {key!r}"

    # Sanity. For 1000 simulated paired-end reads against the lambda
    # reference with 0.1% error, expect:
    # - n_total: most likely 0-5 (wgsim does not introduce real variants;
    #   any non-zero count is from sequencing errors stacking up).
    # - n_snp + n_indel == n_total.
    # - samples has exactly one entry.
    # - mean_dp is finite (not NaN).
    assert summary["n_total"] >= 0, "n_total cannot be negative"
    assert summary["n_snp"] + summary["n_indel"] == summary["n_total"], (
        f"snp+indel = {summary['n_snp']+summary['n_indel']} != "
        f"n_total = {summary['n_total']}"
    )
    assert len(summary["samples"]) == 1, (
        f"expected 1 sample; got {len(summary['samples'])}"
    )

    print()
    print("[exercise-01] VCF summary:")
    print(f"  Total variants:       {summary['n_total']}")
    print(f"  SNPs:                 {summary['n_snp']}")
    print(f"  Indels:               {summary['n_indel']}")
    print(f"  PASS variants:        {summary['n_pass']}")
    print(f"  LowQual variants:     {summary['n_low_qual']}")
    print(f"  Mean QUAL:            {summary['mean_qual']:.2f}")
    print(f"  Mean depth:           {summary['mean_dp']:.2f}")
    print(f"  Samples:              {summary['samples']}")
    print()
    print("[exercise-01] All assertions passed. The VCF at")
    print(f"[exercise-01]   {norm_vcf}")
    print("[exercise-01] is ready for Exercise 2 (VCF parsing by hand).")
    print("[exercise-01] Continue to exercise-02.")
