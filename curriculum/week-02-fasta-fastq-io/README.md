# Week 2 — FASTA, FASTQ, and Biopython IO

In Week 1 you parsed a FASTA file in pure Python so you would understand it from the inside out. In Week 2 we **graduate to Biopython** and we add the file format that every short-read pipeline begins with: **FASTQ**. By Friday of Week 2 you will be able to parse, inspect, filter, and re-emit FASTA and FASTQ files at production scale — and you will know what a Phred quality score is, what FastQC actually computes, and why a single tab character in the wrong column can destroy a downstream pipeline.

The other half of the week is **quality**. Sequencers do not produce perfect reads. Every base in a FASTQ file carries a probability that the machine is wrong about it. We will learn the Phred math (Q = -10 log10(P_error)), the encoding zoo (Sanger / Solexa / Illumina 1.3 / 1.5 / 1.8 — the historical disaster you must know enough about to recognize and route around), and the trimming heuristics that decide which bases survive into your downstream analysis. The mini-project is a real QC report on a 1000 Genomes subset.

## Learning objectives

By the end of this week, you will be able to:

- **Read** a FASTA file with `Bio.SeqIO.parse` and explain when you would use `parse` vs `read` vs `to_dict` vs the lower-level `Bio.SeqIO.FastaIO` parser.
- **Write** a FASTA file with `Bio.SeqIO.write` using a list of `SeqRecord` objects you constructed yourself.
- **Parse** a FASTQ file, extract per-base quality scores, and compute the mean quality and length distribution.
- **Convert** a Phred quality score to an error probability and back — by hand, on paper, and in Python.
- **Recognize** Sanger / Phred+33 vs Illumina 1.3 / Phred+64 encodings from the ASCII histogram of a file, and articulate why this used to be a real production problem.
- **Trim** reads on quality using a sliding-window or 3'-end heuristic, and quantify how many bases survived.
- **Run** FastQC on a real FASTQ file and read its HTML output critically — what's a real problem, what's a cosmetic warning.
- **Distinguish** between a `SeqRecord` and a plain `str`, and choose the right one for the job (hint: `str` for tight inner loops, `SeqRecord` whenever you need to round-trip metadata).

## Prerequisites

This week assumes Week 1 is **done and committed**. Specifically:

- You can parse a FASTA file in pure Python (Exercise 1.2). If `parse_fasta(text)` is not in your muscle memory, go back and finish it.
- You can compute reverse complement and GC content from scratch (Challenge 1.1).
- You have a public GitHub portfolio repo (`crunch-bio-portfolio-<yourhandle>`) with `week-01/` already committed.
- You have **Python 3.11+** and can `pip install`. We will pin Biopython 1.83.

You also need ~2 GB of free disk for the mini-project FASTQ files. The 1000 Genomes subset we use is small (a handful of samples, chr22 only) but FASTQ compresses poorly compared to BAM, and you will want headroom.

## Topics covered

- The FASTA line format — the `>header` convention, the multi-FASTA file, the "all lines until next `>`" rule
- Description-line conventions across NCBI (`>NC_045512.2 description`), GenBank, UniProt (`>sp|P0DTC2|SPIKE_SARS2 ...`)
- `Bio.SeqIO.parse`, `Bio.SeqIO.read`, `Bio.SeqIO.to_dict`, `Bio.SeqIO.write` — the four functions that cover 95% of daily work
- `SeqRecord`, `Seq`, `SeqFeature` — what each is, when to use which
- The FASTQ 4-line-per-record format — `@header` / sequence / `+` / quality
- Phred quality scores: the math (Q = -10 log10(P_error)), the ASCII offset, the dynamic range
- The encoding zoo: Sanger Phred+33, Solexa Phred+64, Illumina 1.3 Phred+64, Illumina 1.5 Phred+64-with-quirk, Illumina 1.8+ Phred+33 (now the de facto standard)
- Read-trimming heuristics: hard 3'-end cut, sliding window, BWA-style soft-clip-ready trim
- FastQC: what each of its eleven modules tests, which warnings are real, which are cosmetic
- Streaming vs in-memory parsing — when a 50 GB FASTQ file forces you to write a generator

## Weekly schedule

The schedule below adds up to approximately **36 hours**. Treat it as a target. If you have done the Week 1 exercises and committed them, this week's IO work goes fast; the QC mini-project is where the time accumulates.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | FASTA with Biopython, SeqRecord                    |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | FASTQ format, Phred quality scores                 |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Wednesday | Encoding zoo, FastQC walk-through                  |    1h    |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     6.5h    |
| Thursday  | Quality trimming, filtering                        |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Streaming large files, QC report writing           |    0h    |    1h     |     1h     |    0.5h   |   1h     |     2h       |    0h      |     5.5h    |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **7.5h**  | **3h**     | **3h**    | **6h**   | **8h**       | **2h**     | **35.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | Curated Biopython tutorial chapters, FASTQ specs, FastQC docs, 1000 Genomes data portal |
| [lecture-notes/01-fasta-the-line-format.md](./lecture-notes/01-fasta-the-line-format.md) | FASTA in depth, Biopython `SeqIO`, `SeqRecord`, write-out |
| [lecture-notes/02-fastq-and-quality-scores.md](./lecture-notes/02-fastq-and-quality-scores.md) | FASTQ format, Phred scores, encoding zoo, FastQC |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-parse-fasta-biopython.py](./exercises/exercise-01-parse-fasta-biopython.py) | Re-parse last week's FASTA with Biopython |
| [exercises/exercise-02-fastq-quality-plot.py](./exercises/exercise-02-fastq-quality-plot.py) | Plot per-base quality across a FASTQ file |
| [exercises/exercise-03-filter-by-quality.py](./exercises/exercise-03-filter-by-quality.py) | Filter reads on mean quality and length |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-streaming-large-fasta.md](./challenges/challenge-01-streaming-large-fasta.md) | Stream a multi-GB FASTA without loading it into RAM |
| [quiz.md](./quiz.md) | 10 multiple-choice questions, format and quality |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | QC report on a 1000 Genomes chr22 FASTQ subset |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("Biopython 1.83," "FastQC 0.12.1"). We cite accessions ("ERR1019034," "NC_045512.2," "GRCh38.p14"). We say "the mean per-base quality across read positions 1–35 is Q30," not "the quality looks fine." A Phred score is a probability statement, not a vibe. If your QC report uses the word "good" without a number attached, you have not written a QC report yet.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Open `Bio/SeqIO/FastaIO.py` and `Bio/SeqIO/QualityIO.py` in the Biopython source. Read each in full. Compare the FASTA parser to the one you wrote in Week 1 Exercise 2.
- Download a recent (2024 or later) SARS-CoV-2 Nanopore FASTQ from ENA — note that long-read quality profiles look very different from short-read Illumina, and that FastQC's interpretation of "good" was calibrated for Illumina.
- Re-implement `Bio.SeqIO.parse` as a generator that takes a file path and yields `SeqRecord` objects. You will end up with ~40 lines and a much better intuition for what the library is hiding from you.
- Read the FASTQ Wikipedia article end-to-end. The history section is genuinely fascinating — Illumina changed their quality encoding three times in five years and the field is still scarred.

## Up next

Continue to [Week 3 — Pairwise alignment](../week-03/) once you have pushed your mini-project QC report to GitHub.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
