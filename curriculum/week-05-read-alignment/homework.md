# Week 5 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-05/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Align E. coli reads against the K-12 reference

**Problem statement.** Fetch the *E. coli* K-12 MG1655 reference (`NC_000913.3`, ~4.6 Mb) and a small read subset (e.g. 100,000 paired-end reads from `SRR1770413` via `fasterq-dump --split-files -X 100000 SRR1770413`). Index the reference with `bwa index`, align with `bwa mem -t 4`, sort to a BAM with `samtools sort`, and index with `samtools index`. Time the alignment with `time` and record the wall-clock duration in `homework/notes/p1-ecoli-align.md`.

Answer:

1. How long did `bwa mem` take?
2. What was the mapped-read percentage from `samtools flagstat`?
3. What was the properly-paired percentage?
4. What was the mean coverage from `samtools coverage`?

**Acceptance criteria.**

- `homework/p1_align.sh` runs end to end on the read subset.
- BAM at `homework/aln/p1.sorted.bam` is committable (or in `.gitignore` if too large; commit `samtools flagstat` output instead).
- `notes/p1-ecoli-align.md` contains four numbered, numeric answers.
- Commit message like `p1: SRR1770413 subset against NC_000913.3, X.X% mapped`.

**Hint.** Use `fasterq-dump --split-files -X 100000 SRR1770413 -O reads/` to get the first 100,000 read pairs. The resulting two FASTQ files are ~80 MB compressed each, manageable on a laptop.

**Estimated time.** 60 minutes.

---

## Problem 2 — minimap2 vs bwa mem on the same dataset

**Problem statement.** Repeat Problem 1's alignment with `minimap2 -ax sr -t 4` instead of `bwa mem`. Sort and index the resulting BAM. In `homework/notes/p2-minimap2-vs-bwa.md`, compare:

1. Wall-clock time for each aligner.
2. Number of mapped reads (`samtools view -F 4 -c`).
3. Number of reads with `MAPQ ≥ 30`.
4. Sample 1000 read names that are mapped in both BAMs. For each pair of records (one from each BAM), check whether the `RNAME` and `POS` agree. What fraction agree exactly? What fraction disagree by ≤ 5 bp?
5. Pick one read where bwa mem and minimap2 disagree on position by > 100 bp. Examine the CIGAR strings and explain the disagreement.

**Acceptance criteria.**

- `homework/p2_minimap2_align.sh` runs.
- BAM at `homework/aln/p2.sorted.bam` is committable or has flagstat output committed.
- `notes/p2-minimap2-vs-bwa.md` contains five numbered answers with specific numbers.
- Commit message like `p2: bwa vs minimap2 on SRR1770413 subset`.

**Hint.** For the per-read comparison, use `pysam` to load both BAMs, build a dict keyed on QNAME from one, and look up each read of the other. Read names in paired-end BAMs appear twice (once per mate); decide whether to compare R1-only, R2-only, or both.

**Estimated time.** 75 minutes.

---

## Problem 3 — Decode SAM flags for 20 reads by hand

**Problem statement.** Take the first 20 alignment records from your Problem 1 BAM: `samtools view homework/aln/p1.sorted.bam | head -20`. For each record, manually decode the FLAG field bit by bit and write a Markdown table at `homework/notes/p3-flag-decoding.md` with columns:

| QNAME | FLAG | Decoded bits | Interpretation in plain English |
|-------|-----:|--------------|---------------------------------|

E.g. for FLAG 99 the row would be:

| SRR1770413.1 | 99 | PAIRED + PROPER_PAIR + MREVERSE + READ1 | Properly-paired forward read whose mate is on the reverse strand. |

Do not use `samtools view -h` or any decoder — work from the binary representation. Verify your decoding at the end with the Broad's "explain flags" tool at <https://broadinstitute.github.io/picard/explain-flags.html>.

**Acceptance criteria.**

- `notes/p3-flag-decoding.md` has 20 rows, each correctly decoded.
- Each row's interpretation is one English sentence.
- At least one row is not a "properly paired" record (you should encounter at least one secondary, supplementary, unmapped, or duplicate read in 20 records — the dataset has a few percent of these).

**Hint.** Treat FLAG as a 12-bit field. `flag & 0x40` is nonzero ↔ READ1 bit is set. Walk through all 12 bits once and you have the decoding.

**Estimated time.** 45 minutes.

---

## Problem 4 — Long-read alignment with minimap2

**Problem statement.** Fetch a small Oxford Nanopore *E. coli* dataset, `SRR2014925` (or a subset via `fasterq-dump --split-files -X 5000 SRR2014925`). Align it against `NC_000913.3` with `minimap2 -ax map-ont -t 4`. Sort and index. Compare the resulting BAM to the short-read BAM from Problem 1.

In `homework/notes/p4-nanopore.md`:

1. What is the median read length in the nanopore FASTQ? (Compare to the 150 bp Illumina median.) Use `awk` on the FASTQ: `awk 'NR%4==2 {print length($0)}' reads.fq | sort -n | awk 'NR==n/2 {print}'`.
2. What is the mapped-read fraction for the long-read BAM?
3. What is the distribution of MAPQ for long reads vs short reads? (Compute the fraction at MAPQ ≥ 30 for each.)
4. Pick one long read with a CIGAR containing more than 5 indel operations. Print the CIGAR and the read length. Interpret: are the indels real biology (a structural variant) or noise (Nanopore's ~5–10% error rate)?

**Acceptance criteria.**

- `homework/p4_nanopore_align.sh` runs.
- `notes/p4-nanopore.md` contains four numbered answers.
- Long-read mapping rate is comparable to short-read (within ±5%); if it is wildly different, something is wrong with the preset choice.

**Hint.** Long reads have ~5–15% error in Oxford Nanopore. Most "indels" in the CIGAR are sequencing errors, not biology. A real structural variant would show as a single long `I` or `D` operation (> 50 bp), not as many small ones scattered through the read.

**Estimated time.** 60 minutes.

---

## Problem 5 — Coverage QC on a low-coverage region

**Problem statement.** Using your Problem 1 BAM, run `samtools depth -a` and find the lowest-coverage 1 kb window in the reference. Investigate what the underlying biology is.

In `homework/p5_coverage_qc.py`:

1. Read the `samtools depth -a` output into pandas.
2. Compute mean coverage in non-overlapping 1 kb windows.
3. Sort by mean coverage ascending.
4. Print the bottom 5 windows with their coordinates.
5. Pick the bottom window. Look it up in NCBI's genome browser or `samtools faidx ecoli.fa NC_000913.3:<start>-<end>` and ask: what gene or region is here? Is it a known repetitive or hard-to-align region (rRNA operon, transposon, low-complexity)?

Write the answer in `homework/notes/p5-low-coverage.md` (150–250 words).

**Acceptance criteria.**

- `p5_coverage_qc.py` runs and prints the bottom 5 windows.
- `notes/p5-low-coverage.md` identifies the biological feature at the lowest-coverage window.
- Speculative claims are flagged ("likely a paralogous region") rather than stated as fact.

**Hint.** *E. coli* MG1655 has 7 rRNA operons (`rrnA` through `rrnH`) that each carry essentially identical 16S, 23S, and 5S rRNA genes. Reads from these regions multimap and many tools (including BWA-MEM with default `-c 500`) drop them, producing characteristic coverage dips. If your low-coverage region matches an `rrn*` operon, that is the textbook answer.

**Estimated time.** 60 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300–400 word reflection at `homework/notes/week-05-reflection.md` answering:

1. Before Week 5, what did you think SAM/BAM was? What is it actually? Pick one column of the eleven you found most surprising and say why.
2. The first time you saw a `MAPQ = 0` read in a real BAM, what did you assume it meant? After Week 5, what does it actually mean? In what way does the distinction matter for variant calling?
3. The FM-index lets BWA-MEM align 150 bp reads against a 3 Gb human reference in `O(m)` time per read instead of `O(m·n)`. In your own words, what is the data structure doing — and why is it stunning that a 3 GB index can answer "exact occurrences of this 150-character string" in microseconds?
4. The mini-project asks you to produce a coverage plot from a real SRA dataset. What QC signal in that plot — coverage uniformity, mean depth, fraction of zero-coverage bases, distribution of MAPQ — would convince you the alignment is healthy? Name one specific number you would expect to see and one number that would worry you.

**Acceptance criteria.**

- File exists, 300–400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The mistakes you note here are what you will re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 1 h 0 min |
| 2 | 1 h 15 min |
| 3 | 45 min |
| 4 | 1 h 0 min |
| 5 | 1 h 0 min |
| 6 | 30 min |
| **Total** | **~5 h 30 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
