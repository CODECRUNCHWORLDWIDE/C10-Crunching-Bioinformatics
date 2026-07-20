# Week 6 — Exercises

Three focused drills. Each one runs from the command line, finishes inside 60 minutes, and produces output you can keep. Together they take you from "I read about VCF" to "I can call variants, parse a VCF by hand, and annotate the result with VEP." The mini-project then composes these three into a full BAM-to-annotated-VCF pipeline on the SRR1770413 read set you aligned in Week 5.

## Index

1. **[Exercise 1 — Call variants with bcftools](exercise-01-call-with-bcftools.py)** — index the lambda reference with `samtools faidx`, run the canonical `bcftools mpileup | bcftools call` pipeline against the Week 5 Exercise 1 BAM, apply hard filters with `bcftools filter`, normalize indel representation with `bcftools norm`, and verify the output with `pysam.VariantFile`. (~50 min, mostly local — no internet needed.)
2. **[Exercise 2 — Parse VCF by hand](exercise-02-vcf-parse.py)** — read a small VCF file as plain text, split each record into the eight mandatory columns plus FORMAT and per-sample columns, decode INFO and FORMAT fields, and classify each variant as SNP or indel. (~45 min, pure Python — no bcftools or VEP.)
3. **[Exercise 3 — VEP annotation](exercise-03-vep-annotate.py)** — send a small set of variants to the Ensembl VEP REST API, parse the JSON response, and tabulate per-transcript consequence annotations. (~60 min, requires an internet connection for the REST API.)

## How to work the exercises

- Install the toolchain *first*: `conda install -c bioconda bcftools=1.19 samtools=1.19 pysam=0.22 biopython=1.83 pandas matplotlib` and `pip install requests==2.31`. The exercises assume `bcftools` and `samtools` are on your `PATH`.
- Exercise 1 reads from the Week 5 Exercise 1 BAM (`exercises/aln/lambda.sorted.bam` in the week-05 directory). Run Week 5 Exercise 1 first; this exercise *does not* re-align reads.
- **Type the code yourself.** Copy-paste defeats the muscle-memory point.
- Each file is runnable: `python exercise-XX.py`. All assertions must pass.
- Exercises 1 and 3 produce on-disk artifacts (`.vcf.gz`, `.tbi`, `.tsv`). Commit the TSVs and small VCFs to your portfolio repo; large VCFs (the mini-project's) should be rebuilt by `bash run.sh` when needed.
- If you get stuck for more than 10 minutes, peek at the inline hints, then re-read Lecture 1 (for Exercises 1 and 2) or Lecture 2 (for Exercise 3).
- No solutions are checked in. The exercises are self-checking — every script ends in an assertion block.

## Accession IDs used this week

Cited by NCBI / SRA accession so you can verify your data is the same as the curriculum's:

- **`NC_001416.1`** — Bacteriophage lambda complete genome (48,502 bp). Used in Exercise 1 (called against the Week 5 Exercise 1 BAM).
- **`NC_000913.3`** — *Escherichia coli* str. K-12 substr. MG1655 complete genome (4,641,652 bp). Used in the mini-project.
- **`SRR1770413`** — Illumina HiSeq paired-end resequencing of *E. coli* K-12 MG1655 (~5 GB compressed). The mini-project read set; not used in the exercises themselves.

If any of these accessions has been retired or updated by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt.

## Tool versions used in the exercises

- bcftools 1.19
- samtools 1.19
- pysam 0.22
- requests 2.31 (Exercise 3 only)
- pandas 2+
- ensembl-vep 110 (optional offline alternative to the REST API in Exercise 3)
- Python 3.11+
