"""
Exercise 1 - Align a small genome end to end with bwa mem.

Goal: fetch the bacteriophage lambda reference (NC_001416.1, 48,502 bp),
simulate ~1000 paired-end Illumina-like reads against it (or use the
pre-generated reads at reads/lambda_R1.fq.gz and reads/lambda_R2.fq.gz
if present), index the reference with `bwa index`, align with `bwa mem`,
sort and index the BAM with samtools, and verify the resulting BAM with
pysam.

Estimated time: 45 minutes (most spent on toolchain debugging the first
time you run bwa/samtools on your machine).

Acceptance criteria:
- `python exercise-01-align-small-genome.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The reference FASTA is cached to `ref/lambda.fa` on first run; on
  subsequent runs the script reads from cache and does NOT hit NCBI.
- After the script runs, `aln/lambda.sorted.bam` and
  `aln/lambda.sorted.bam.bai` exist on disk.
- You implemented four functions: `fetch_reference`, `index_reference`,
  `align_reads`, and `summarize_bam`.

Requirements:
    conda install -c bioconda bwa=0.7.17 samtools=1.19 pysam=0.22 biopython=1.83
    (and a real network connection to NCBI on first run for the reference)

What you learn:
- The full bwa index -> bwa mem -> samtools sort -> samtools index pipeline,
  invoked from Python via subprocess.
- How pysam reads a sorted+indexed BAM and lets you iterate over reads.
- How to read a SAM flag programmatically (read.is_paired, read.is_unmapped,
  read.is_reverse, etc.) instead of decoding bits by hand.
- Why every downstream tool requires the .bai index.

TO COMPLETE: implement the four functions below. Run the file; all
assertions must pass.

Reference: bacteriophage lambda is the textbook tiny genome - it shows
up in every bioinformatics tutorial for a reason. The full BWA pipeline
on lambda takes ~5 seconds end to end, so you can iterate fast.

Tool versions assumed:
- BWA 0.7.17
- samtools 1.19
- pysam 0.22
- Biopython 1.83
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

from Bio import Entrez, SeqIO


# NCBI requires a contact email on every E-utilities call.
# REPLACE THIS with your real address before running.
Entrez.email = "you@example.com"


# Reference accession. NC_001416.1 is the bacteriophage lambda complete
# genome, 48,502 bp - small enough that bwa index runs in ~1 second.
REFERENCE_ACCESSION = "NC_001416.1"
EXPECTED_REF_LENGTH = 48502


# Read files. The exercise expects these to exist in reads/ already
# (generated upstream with wgsim or similar). If they do not exist,
# the script will produce a helpful error message.
# To generate them yourself once with wgsim:
#   wgsim -N 1000 -1 150 -2 150 -e 0.001 ref/lambda.fa \
#       reads/lambda_R1.fq reads/lambda_R2.fq
#   gzip reads/lambda_R1.fq reads/lambda_R2.fq
EXPECTED_R1 = "reads/lambda_R1.fq.gz"
EXPECTED_R2 = "reads/lambda_R2.fq.gz"


def fetch_reference(accession: str, out_path: Path) -> Path:
    """Fetch a single nucleotide sequence from NCBI nuccore.

    If `out_path` already exists, do nothing (idempotent). Otherwise,
    fetch via `Bio.Entrez.efetch(db="nuccore", id=accession, rettype="fasta")`,
    write the FASTA text to `out_path`, and return the path.

    Hint: `Entrez.efetch(...)` returns a handle; `.read()` gives the
    text. Write it to disk with `out_path.write_text(...)`.

    Returns:
        The path to the FASTA file on disk.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        return out_path

    print(f"[exercise-01] Fetching {accession} from NCBI ...")
    # TODO: use Entrez.efetch to download the FASTA for `accession`.
    # Write the text response to out_path.
    raise NotImplementedError("Fetch the reference via Entrez.efetch")


def index_reference(ref_path: Path) -> None:
    """Run `bwa index` on `ref_path` if not already indexed.

    BWA produces five auxiliary files: ref.fa.amb, ref.fa.ann,
    ref.fa.bwt, ref.fa.pac, ref.fa.sa. Their presence is the indicator
    that the index already exists; do not re-index.

    Hint: `subprocess.run(["bwa", "index", str(ref_path)], check=True)`
    is the canonical call. Check for one of the five output files
    (e.g. `{ref_path}.bwt`) to decide whether to skip.

    Also run `samtools faidx ref_path` to produce ref.fa.fai for
    random access by pysam and bcftools.
    """
    bwt_marker = ref_path.parent / f"{ref_path.name}.bwt"
    fai_marker = ref_path.parent / f"{ref_path.name}.fai"

    # TODO: if bwt_marker does not exist, run `bwa index ref_path`.
    # TODO: if fai_marker does not exist, run `samtools faidx ref_path`.
    # Use `subprocess.run([...], check=True)` so a failure raises
    # CalledProcessError loud and clear.
    raise NotImplementedError("Run bwa index and samtools faidx")


def align_reads(ref_path: Path, r1: Path, r2: Path, out_bam: Path) -> Path:
    """Align paired-end reads, pipe to samtools sort, write sorted BAM.

    This invokes the canonical bwa mem + samtools sort pipeline:

        bwa mem -t 2 \\
            -R '@RG\\tID:lambda\\tSM:lambda\\tLB:lib1\\tPL:ILLUMINA' \\
            ref.fa R1.fq.gz R2.fq.gz \\
        | samtools sort -@ 2 -o out.bam -

    Then index the resulting BAM with `samtools index`.

    Hint: use `subprocess.run(cmd, shell=True, check=True)` with the
    full pipe string. shell=True is appropriate here because the pipe
    is part of the command. Use proper shell-quoting for the read group
    (the @RG string contains tabs - in shell, "@RG\\tID:..." gives the
    literal \\t sequence that bwa parses into a tab).

    Returns the path to the sorted+indexed BAM.
    """
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    # The @RG string. Note the literal backslash-t for the tab separator
    # that bwa mem requires.
    read_group = r"@RG\tID:lambda_test\tSM:lambda\tLB:lib1\tPL:ILLUMINA"

    # TODO: build the shell command as a string, in the canonical
    # bwa mem | samtools sort - pipe form. Use the read_group above.
    # TODO: run it with subprocess.run(cmd, shell=True, check=True).
    # TODO: run samtools index on the resulting BAM.
    raise NotImplementedError("Run bwa mem | samtools sort - and index")


def summarize_bam(bam_path: Path) -> dict:
    """Open a sorted+indexed BAM with pysam and return summary stats.

    Returns a dict with keys:
        "n_total":         int, total reads in BAM (including unmapped)
        "n_mapped":        int, reads with the 0x4 flag bit unset
        "n_paired":        int, reads with 0x1 set
        "n_proper_pair":   int, reads with 0x2 set
        "n_reverse":       int, reads with 0x10 set (reverse-strand)
        "mapping_rate":    float, n_mapped / n_total, 0.0-1.0
        "contigs":         list of (name, length) tuples from the header

    Hint: use `pysam.AlignmentFile(bam_path, "rb")`. Iterate over all
    reads with `for read in af:`. Use the boolean accessors
    `read.is_unmapped`, `read.is_paired`, etc., rather than decoding the
    flag by hand (Exercise 2 will exercise the hand-decoding).

    The `.references` and `.lengths` attributes give the contig list.
    """
    import pysam  # imported here so the rest of the file can be parsed
    # without pysam present, for the assertion-free dry run.

    # TODO: open the BAM in read-binary mode.
    # TODO: iterate over reads, tallying the flag bits above.
    # TODO: build and return the result dict.
    raise NotImplementedError("Open BAM with pysam and tally flags")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-01-align-small-genome.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    ref_dir = Path(__file__).parent / "ref"
    reads_dir = Path(__file__).parent / "reads"
    aln_dir = Path(__file__).parent / "aln"
    ref_path = ref_dir / "lambda.fa"
    r1 = reads_dir / "lambda_R1.fq.gz"
    r2 = reads_dir / "lambda_R2.fq.gz"
    bam_path = aln_dir / "lambda.sorted.bam"

    # Toolchain sanity check.
    for tool in ("bwa", "samtools"):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"[exercise-01] {tool!r} is not on PATH. Install it:\n"
                f"    conda install -c bioconda {tool}=...\n"
                f"and re-run."
            )

    print("[exercise-01] Step 1: fetching reference ...")
    fetch_reference(REFERENCE_ACCESSION, ref_path)
    assert ref_path.exists(), f"reference FASTA not at {ref_path}"
    record = SeqIO.read(ref_path, "fasta")
    assert len(record.seq) == EXPECTED_REF_LENGTH, (
        f"lambda reference length {len(record.seq)} != expected "
        f"{EXPECTED_REF_LENGTH}; did NCBI swap the version?"
    )
    print(f"[exercise-01] Reference: {record.id}, {len(record.seq)} bp.")

    print("[exercise-01] Step 2: indexing reference ...")
    index_reference(ref_path)
    assert (ref_dir / "lambda.fa.bwt").exists(), "bwa index did not produce .bwt"
    assert (ref_dir / "lambda.fa.fai").exists(), "samtools faidx did not produce .fai"

    # Confirm reads exist before alignment.
    if not r1.exists() or not r2.exists():
        raise SystemExit(
            f"[exercise-01] Expected paired-end reads at:\n"
            f"    {r1}\n"
            f"    {r2}\n"
            f"Generate them with wgsim before running this exercise:\n"
            f"    wgsim -N 1000 -1 150 -2 150 -e 0.001 {ref_path} \\\n"
            f"        {reads_dir}/lambda_R1.fq {reads_dir}/lambda_R2.fq\n"
            f"    gzip {reads_dir}/lambda_R1.fq {reads_dir}/lambda_R2.fq\n"
            f"(wgsim ships with samtools-utils or its own bioconda package.)"
        )

    print("[exercise-01] Step 3: aligning reads ...")
    align_reads(ref_path, r1, r2, bam_path)
    assert bam_path.exists(), f"sorted BAM not at {bam_path}"
    assert bam_path.with_suffix(".bam.bai").exists(), (
        f"BAM index .bai not at {bam_path}.bai - did samtools index run?"
    )

    print("[exercise-01] Step 4: summarizing BAM ...")
    summary = summarize_bam(bam_path)

    # Field presence.
    for key in (
        "n_total", "n_mapped", "n_paired",
        "n_proper_pair", "n_reverse", "mapping_rate", "contigs",
    ):
        assert key in summary, f"summary missing field {key!r}"

    # Sanity. For 1000 simulated paired-end reads against the lambda
    # reference with 0.1% error, expect:
    # - n_total = 2000 (paired reads, both ends)
    # - n_mapped > 1990 (mapping rate > 99.5%)
    # - n_paired = 2000 (every read is paired)
    # - n_proper_pair > 1960 (proper-pair rate > 98%)
    # - contigs has exactly one entry, ('NC_001416.1', 48502)
    assert summary["n_total"] >= 1000, (
        f"n_total = {summary['n_total']}; expected at least 1000"
    )
    assert summary["mapping_rate"] > 0.95, (
        f"mapping rate = {summary['mapping_rate']:.3f}; expected > 0.95"
    )
    assert len(summary["contigs"]) == 1, (
        f"expected 1 contig (lambda only); got {len(summary['contigs'])}"
    )
    assert summary["contigs"][0][1] == EXPECTED_REF_LENGTH, (
        f"contig length {summary['contigs'][0][1]} != "
        f"{EXPECTED_REF_LENGTH}"
    )

    print()
    print("[exercise-01] Alignment summary:")
    print(f"  Total reads:    {summary['n_total']}")
    print(f"  Mapped:         {summary['n_mapped']} "
          f"({100*summary['mapping_rate']:.2f}%)")
    print(f"  Paired:         {summary['n_paired']}")
    print(f"  Proper pair:    {summary['n_proper_pair']}")
    print(f"  Reverse strand: {summary['n_reverse']}")
    print(f"  Contigs:        {summary['contigs']}")
    print()
    print("[exercise-01] All assertions passed. The BAM at")
    print(f"[exercise-01]   {bam_path}")
    print("[exercise-01] is ready for Exercise 3 (coverage plot).")
    print("[exercise-01] Continue to exercise-02 (SAM by hand).")
