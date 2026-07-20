# Week 5 — Read Alignment

In Week 4 you built a BLAST-driven taxonomy classifier that takes one ~1.5 kb 16S rRNA sequence at a time and pushes it through NCBI's heuristic seed-and-extend pipeline against a curated 25,000-entry database. That worked at the "one query, big database" end of the search problem. In Week 5 we rotate the problem 90 degrees: **one reference, millions of short reads**. The typical short-read sequencing run on an Illumina NovaSeq produces 1–10 billion 150 bp reads. Mapping every one of those against a 3 Gb human genome with Smith-Waterman is `O(read_count · m · n) ≈ 10^9 · 150 · 3·10^9 ≈ 4.5·10^20` cell updates — about 14,000 years on a single core at 10^9 cells/second. The wall is the same wall BLAST hit in 1990, but in a regime where the queries are tiny and uniform and the reference is huge and fixed. The algorithmic answer is the **FM-index** (Ferragina & Manzini 2000) layered on the **Burrows-Wheeler Transform** (Burrows & Wheeler 1994), and the production tools are `bwa mem` (Li 2013) for short Illumina reads and `minimap2` (Li 2018) for long Oxford Nanopore / PacBio reads. Both ship with all major Linux distributions, both are < 2 MB binaries, both publish papers worth reading, and both are what you will be running on real data from Monday onward.

By Friday of Week 5 you will be able to index a reference genome with `bwa index` and `minimap2 -d`, align paired-end Illumina reads to it with `bwa mem` producing valid SAM output, convert SAM to sorted-and-indexed BAM with `samtools sort` + `samtools index`, read a BAM file with `pysam` and pull out per-read information (the reference position, the CIGAR string, the mapping quality, the SAM flag), compute and plot per-position coverage from a BAM file, mark PCR/optical duplicates with `samtools markdup` (or `picard MarkDuplicates`), and explain in one sentence each what a SAM flag bit means, what a CIGAR operation does, and why mapping quality is on a `-10·log10(P_wrong)` scale. The mini-project takes a real public small-genome dataset (the E. coli K-12 SRA run SRR1770413, ~5 GB of paired-end Illumina reads against the 4.6 Mb MG1655 reference `NC_000913.3`) and walks it through the full FASTQ-to-coverage-plot pipeline, ending in a duplicate-marked sorted BAM, a per-position coverage plot, and a short written interpretation of the coverage distribution.

The other half of the week is **what the FM-index gives you that hash-table-based aligners did not**. BLAST's `blastn` builds a hash table of k-mer occurrences in the database — fine for one query against a 25 kb database, but a hash table of 28-mers in a 3 Gb human genome is itself ~12 GB. BWA's FM-index of the human genome is ~3 GB on disk and queryable in O(read_length) time per read, independent of the genome's size. The price is conceptual complexity (the FM-index is a non-trivial data structure that took two decades of indexing-theory work to mature) and a one-time indexing cost (`bwa index hg38.fa` takes ~90 minutes on a laptop). After indexing, alignment is roughly **a million reads per minute per core**. That is the speed that makes whole-genome sequencing a routine clinical workflow rather than a heroic research project. Once you have internalized the FM-index trick you will see it everywhere: in `bowtie2` (Langmead & Salzberg 2012), in `minimap2`'s minimizer-and-chain seeding, in HISAT2's hierarchical FM-index for spliced RNA alignment, in compressed text indexes for any large reference.

## Learning objectives

By the end of this week, you will be able to:

- **Describe** the FM-index + BWT data structure in two paragraphs, naming the rank/select operations and explaining why backward search runs in `O(m)` time independent of reference length `n`.
- **Choose** between `bwa mem` and `minimap2` for a given read profile (read length, expected error rate, paired vs single-end, RNA vs DNA) and defend the choice in one sentence.
- **Run** `bwa index` on a reference FASTA to produce the five-file index, then `bwa mem -t 4 ref.fa R1.fq.gz R2.fq.gz > out.sam` to align paired-end reads.
- **Run** `minimap2 -ax map-ont ref.fa nanopore.fq.gz > out.sam` for long-read alignment with the appropriate preset (`map-ont`, `map-pb`, `sr`, `splice`).
- **Convert** SAM to sorted-and-indexed BAM via `samtools view -bS in.sam | samtools sort -o out.bam` then `samtools index out.bam`, and explain why every downstream tool requires the index.
- **Read** the eleven mandatory SAM fields by name (`QNAME`, `FLAG`, `RNAME`, `POS`, `MAPQ`, `CIGAR`, `RNEXT`, `PNEXT`, `TLEN`, `SEQ`, `QUAL`) and decode a SAM `FLAG` field bit by bit (`0x1` = paired, `0x4` = unmapped, `0x10` = reverse-strand, etc.).
- **Parse** a CIGAR string (e.g. `36M2I12M1D80M`) into a list of `(op, length)` tuples and compute the implied reference span and query span.
- **Interpret** a mapping quality value `MAPQ = 60` (probability of misalignment `≈ 10^-6`) vs `MAPQ = 0` (multimapper, no unique placement).
- **Compute** per-position coverage from a sorted BAM via `samtools depth` and visualize it as a coverage plot with matplotlib.
- **Mark** PCR/optical duplicates with `samtools markdup` (after `samtools fixmate -m`), and explain the difference between *removing* duplicates and *marking* them with the `0x400` flag bit.
- **Identify** at least three failure modes of short-read alignment (low-complexity regions giving high multimapper rates, large structural variants exceeding the aligner's gap-open threshold, contamination from adapters or other organisms) and the standard QC signal each one produces.

## Prerequisites

This week assumes Weeks 1, 2, 3, and 4 are **done and committed**. Specifically:

- You can parse a FASTQ file with `Bio.SeqIO.parse` and pull a list of reads with quality scores into Python (Week 2 Exercises 1 and 2).
- You can read a BLAST tabular output into pandas, filter by E-value, and merge with a metadata table (Week 4 Exercise 3).
- You have a working `crunch-bio-portfolio-<yourhandle>/` repo with `week-04/` committed. The Week 5 mini-project lives in `week-05/` alongside it.
- You have Python 3.11+, Biopython 1.83, and pandas installed from Week 4. You will need to install BWA 0.7.17, minimap2 2.26, samtools 1.19, and pysam 0.22 this week — `conda install -c bioconda bwa=0.7.17 minimap2=2.26 samtools=1.19 pysam=0.22` is the canonical path.

You do not need biology beyond "DNA is double-stranded so a read can map to either strand, and Illumina paired-end reads come from the two ends of the same ~300 bp DNA fragment." You do need disk space — even our small mini-project dataset (SRR1770413) is ~5 GB compressed and ~15 GB uncompressed, and the resulting BAM file is ~2 GB. Plan ~25 GB of free disk space for the week.

## Topics covered

- The short-read mapping problem at scale: 1 billion 150 bp reads against a 3 Gb reference is `O(10^20)` Smith-Waterman cell updates per run. We need both an index (sub-linear seed lookup) and a heuristic (skip the obvious non-matches).
- The Burrows-Wheeler Transform: a reversible permutation of the reference that brings similar contexts together, enabling efficient compression and indexing. Reference: Burrows & Wheeler, *DEC SRC Research Report* 124, 1994.
- The FM-index: a compressed self-index over the BWT supporting **rank** and **select** queries in O(1) time. Reference: Ferragina & Manzini, *FOCS* 2000. Practical implementation in `bwa` (Li & Durbin, *Bioinformatics* 25:1754, 2009; revised in Li, *arXiv:1303.3997*, 2013 for `bwa mem`).
- Backward search: query a length-`m` read against an FM-indexed reference in `O(m)` time, returning the list of all exact occurrences as an interval in the suffix array. The seeding pattern for BWA.
- Seed-and-chain (minimap2): split the read into minimizers, look up each in the indexed reference, chain colinear minimizer matches, run a banded dynamic-programming extension on the chains. Reference: Li, *Bioinformatics* 34:3094, 2018.
- BWA-MEM's seeding (Maximal Exact Matches): for each position in the read, find the longest exact match against the reference such that extending by one more base would mismatch. MEMs are the seeds; the chain-and-extend phase fills the gaps. Reference: Li 2013 (`arXiv:1303.3997`).
- The SAM specification (Li et al., *Bioinformatics* 25:2078, 2009): the eleven mandatory tab-separated columns, the optional `tag:type:value` triples, the header lines (`@HD`, `@SQ`, `@RG`, `@PG`). The full spec is at <https://samtools.github.io/hts-specs/SAMv1.pdf> — ~25 pages, surprisingly readable.
- SAM flags: a 12-bit field encoded as an integer. Bit `0x1` = read is paired, `0x2` = both ends mapped in proper pair, `0x4` = read unmapped, `0x10` = read reverse-complemented, `0x40` = first in pair, `0x80` = second in pair, `0x100` = secondary alignment, `0x400` = PCR/optical duplicate, `0x800` = supplementary alignment. Memorize the common ones; look the rest up.
- CIGAR operations: `M` (alignment match — could be sequence match or mismatch), `I` (insertion to reference), `D` (deletion from reference), `S` (soft clip — bases in read but not aligned), `H` (hard clip — bases not in the SEQ at all), `N` (skipped region, e.g. an intron in RNA-seq), `=` (exact match), `X` (mismatch). Modern aligners emit `=` and `X` if asked (`bwa mem -M`); older convention collapses both into `M`.
- Mapping quality (`MAPQ`): the `-10·log10(P[read is misaligned])` Phred-scaled probability the placement is wrong. `MAPQ = 60` is the BWA-MEM ceiling (P_wrong ≤ 10^-6). `MAPQ = 0` is "I have no idea which of N equally-good positions is right" (multimapper). Anything in between is interpolated.
- Coverage: the number of reads spanning each reference position. `samtools depth -a ref.bam` gives per-position coverage; `samtools coverage` gives a per-contig summary. Mean coverage = (total aligned bases) / (reference length); typical whole-genome runs target 30x mean coverage for variant calling.
- PCR and optical duplicates: reads that arose from the same original DNA fragment via PCR amplification (PCR duplicates) or from adjacent optical features on the flow cell (optical duplicates). Both inflate apparent coverage at a single position. `samtools markdup` (or `picard MarkDuplicates`) detects them by 5' alignment position and read orientation, marks them with flag `0x400`, and either leaves them in the BAM (default — preserves data) or removes them (`samtools markdup -r`).
- Common failure modes: **low-complexity regions** (`(AT)_n` repeats, microsatellites) produce high multimapper rates and `MAPQ = 0` reads; **large structural variants** (insertions or deletions > 50 bp) typically exceed BWA-MEM's affine-gap penalty and produce soft-clipped reads instead of aligned indels (need a long-read aligner or a structural-variant caller); **contamination** (adapter sequences, host DNA in a microbial sample) produces reads that either fail to align or align to unexpected references (always run `fastqc` or `kraken` upstream).

## Weekly schedule

The schedule below adds up to approximately **36 hours**. Treat it as a target. Monday's lecture on the FM-index is the hour that decides whether the rest of the week makes algorithmic sense — read it twice if needed.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | FM-index, BWT, backward search, BWA-MEM seeding    |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | Running bwa mem and minimap2 end to end            |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Wednesday | SAM/BAM, samtools, flags, CIGAR, MAPQ              |    1h    |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     6.5h    |
| Thursday  | Coverage, duplicates, pysam parsing                |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Mini-project deep work + plot polish               |    0h    |    1h     |     0h     |    0.5h   |   1h     |     2h       |    0h      |     4.5h    |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **7.5h**  | **2h**     | **3h**    | **6h**   | **8h**       | **2h**     | **34.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | BWA, minimap2, samtools, htslib docs + reference papers |
| [lecture-notes/01-from-fastq-to-bam.md](./lecture-notes/01-from-fastq-to-bam.md) | `bwa index`, `bwa mem`, `minimap2 -ax`, the FASTQ→SAM→BAM flow, choice criteria between BWA and minimap2 |
| [lecture-notes/02-sam-bam-format-and-samtools.md](./lecture-notes/02-sam-bam-format-and-samtools.md) | The SAM spec column by column, decoding FLAG, MAPQ, CIGAR, sorting and indexing, the samtools toolchain |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-align-small-genome.py](./exercises/exercise-01-align-small-genome.py) | Align a tiny E. coli read subset against `NC_000913.3` end to end with `bwa mem` via `subprocess` |
| [exercises/exercise-02-parse-sam-by-hand.py](./exercises/exercise-02-parse-sam-by-hand.py) | Parse a small SAM file by reading text lines, decode FLAG and CIGAR with no library help |
| [exercises/exercise-03-coverage-plot.py](./exercises/exercise-03-coverage-plot.py) | Use `pysam` to walk a sorted BAM, compute per-position coverage, and render a matplotlib plot |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-detect-duplicates.md](./challenges/challenge-01-detect-duplicates.md) | Compute the duplication rate from a BAM file by hand and compare to `samtools markdup`'s answer |
| [quiz.md](./quiz.md) | 10 multiple-choice questions on FM-index, SAM, MAPQ, and the alignment toolchain |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | Align SRR1770413 against `NC_000913.3`, mark duplicates, produce coverage plots, write up |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("BWA 0.7.17," "samtools 1.19," "minimap2 2.26"). We cite tools by their paper ("BWA-MEM, Li *arXiv:1303.3997* 2013"). We say "mean coverage 47.3x, median 46x, coefficient of variation 0.18 across 4.6 Mb of *E. coli* MG1655 reference, with two regions of < 5x coverage at positions 1.06–1.07 Mb and 3.94–3.95 Mb corresponding to the rRNA operons" not "coverage looks even." A coverage plot is a number on a known scale. If your methods section uses the words "high coverage" or "good alignment" without a number, you have not written one yet.

## A note on the data size

Real-world short-read datasets are large. The mini-project SRA run SRR1770413 is:

- **5.0 GB** compressed (`.fastq.gz` paired-end, R1 and R2 combined).
- **15 GB** uncompressed if you decompress (do not, unless you have to — most tools accept `.gz` directly).
- **~2 GB** sorted BAM after alignment.
- **~50 MB** BAM index (`.bai`).

Add the *E. coli* reference (4.6 Mb FASTA, ~5 MB) and the BWA index (5 files totaling ~12 MB), and you are at ~7 GB of disk for the mini-project. That fits on a laptop. The full human-genome equivalent would be ~30 GB of FASTQ, ~80 GB of sorted BAM, and a ~3 GB FM-index. Plan accordingly when you move to human data in Week 6.

If your laptop is short on disk, the SRA toolkit lets you stream the first N reads only:

```bash
fasterq-dump --split-files -p --threads 4 -X 1000000 SRR1770413
```

(takes the first 1,000,000 read pairs — ~200 MB compressed instead of 5 GB, enough to verify the pipeline end to end without committing the full run).

## Stretch goals

If you finish early and want to push further, try any of the following:

- Read the BWA-MEM paper (Li, *arXiv:1303.3997* 2013) and the minimap2 paper (Li, *Bioinformatics* 34:3094, 2018) end to end. The minimap2 paper is particularly readable as a "what changed since BWA" tour of the field.
- Reproduce the FM-index by hand on a small reference. Take the 12-base reference `ACGTACGTACGT$`, build its BWT, build the count array `C[c]` and the occurrence array `Occ[c, i]`, and walk a backward search for the query `CGT` step by step. The Ferragina-Manzini 2000 paper is the rigorous source; section 4 of the BWA 2009 paper gives a worked example in three pages.
- Install [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (Langmead & Salzberg 2012) and benchmark it against `bwa mem` on the mini-project FASTQ. Bowtie2 is the other production short-read aligner; the comparison is informative.
- Run `samtools flagstat` and `samtools idxstats` on your mini-project BAM and learn to read both. `flagstat` summarizes the SAM flag distribution (total reads, paired, properly paired, duplicates, supplementary, mapped %). `idxstats` reports per-reference-contig read counts. Both are one-liners and both appear in every methods section of every paper using short-read data.
- Mark duplicates with **both** `samtools markdup` and `picard MarkDuplicates`, and compare. Picard is the GATK ecosystem's standard; `samtools markdup` is newer (1.10+) and faster but slightly different in optical-duplicate handling. The difference is the kind of detail that comes up in interviews.

## Up next

Continue to [Week 6 — Variant calling](../week-06/) once you have pushed your mini-project to GitHub. Week 6 takes the sorted-and-indexed BAM you produce this week and runs `bcftools mpileup` + `bcftools call` (or GATK `HaplotypeCaller`) to produce a VCF of single-nucleotide variants and short indels, with hard-filter quality control. The coverage plot you build this week is the QC signal that decides whether variant calling will work at all.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
