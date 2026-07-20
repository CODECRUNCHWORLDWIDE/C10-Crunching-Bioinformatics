# Week 5 — Exercises

Three focused drills. Each one runs from the command line, finishes inside 60 minutes, and produces output you can keep. Together they take you from "I read about FM-indexes" to "I can align reads, parse SAM by hand, and produce a coverage plot." The mini-project then composes these three into a full FASTQ-to-coverage-plot pipeline on a real SRA dataset.

## Index

1. **[Exercise 1 — Align a small genome](exercise-01-align-small-genome.py)** — index the bacteriophage lambda reference (`NC_001416.1`, 48.5 kb), align a tiny simulated read set against it with `bwa mem` via `subprocess`, sort + index the resulting BAM, and verify the output with `pysam`. (~45 min, mostly local — no NCBI calls beyond the one-time reference fetch.)
2. **[Exercise 2 — Parse SAM by hand](exercise-02-parse-sam-by-hand.py)** — read a small SAM file as plain text, decode the FLAG field bit by bit, parse the CIGAR string into operations, and compute the query/reference spans without any library help. (~40 min, pure Python — no BWA or samtools.)
3. **[Exercise 3 — Coverage plot](exercise-03-coverage-plot.py)** — use `pysam` to walk a sorted BAM, compute per-position coverage, and render a matplotlib plot with mean/median lines. (~50 min, pure Python after the BAM is built.)

## How to work the exercises

- Install the toolchain *first*: `conda install -c bioconda bwa=0.7.17 minimap2=2.26 samtools=1.19 pysam=0.22 biopython=1.83 pandas matplotlib`. The exercises assume `bwa`, `samtools`, and `minimap2` are on your `PATH`.
- **Set `Bio.Entrez.email`** to a real address (yours) at the top of every script that touches NCBI. The reference-fetch step in Exercise 1 calls Entrez once.
- **Type the code yourself.** Copy-paste defeats the muscle-memory point.
- Each file is runnable: `python exercise-XX.py`. All assertions must pass.
- Exercises 1 and 3 produce on-disk artifacts (BAM, BAI, PNG). Commit the PNGs to your portfolio repo; the BAMs are usually too large to commit and should be re-built by `bash run.sh` when needed.
- If you get stuck for more than 10 minutes, peek at the inline hints, then re-read Lecture 1 (for Exercise 1) or Lecture 2 (for Exercises 2 and 3).
- No solutions are checked in. The exercises are self-checking — every script ends in an assertion block.

## Accession IDs used this week

Cited by NCBI / SRA accession so you can verify your data is the same as the curriculum's:

- **`NC_001416.1`** — Bacteriophage lambda complete genome (48,502 bp). The toy reference used in Exercise 1 so a full pipeline runs in seconds.
- **`NC_000913.3`** — *Escherichia coli* str. K-12 substr. MG1655 complete genome (4,641,652 bp). The canonical bacterial reference. Used in the mini-project.
- **`SRR1770413`** — Illumina HiSeq paired-end resequencing of *E. coli* K-12 MG1655 (~5 GB compressed). The mini-project read set; not used in the exercises themselves.

If any of these accessions has been retired or updated by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt.

## Tool versions used in the exercises

- BWA 0.7.17
- minimap2 2.26
- samtools 1.19
- pysam 0.22
- Biopython 1.83
- matplotlib 3.8+
- numpy 1.26.4
- pandas 2+
