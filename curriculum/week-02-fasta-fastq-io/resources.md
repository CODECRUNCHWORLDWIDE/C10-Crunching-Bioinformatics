# Week 2 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (Biopython 1.83, FastQC 0.12.1, samtools 1.20), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Biopython Tutorial and Cookbook (1.83)** — full PDF, free. Read chapters 2 (Quick Start), 5 (Sequence Input/Output), 20 (Cookbook):
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **FASTQ format specification** — the Cock et al. 2010 *NAR* paper that finally documented FASTQ formally:
  <https://academic.oup.com/nar/article/38/6/1767/3112533>
- **FastQC documentation** — what each module actually tests, and how to read the report:
  <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/>
- **1000 Genomes Project — data portal**, where the mini-project FASTQ files live:
  <https://www.internationalgenome.org/data-portal/sample>
- **NCBI Sequence Read Archive (SRA)** — the public archive that hosts most published FASTQ:
  <https://www.ncbi.nlm.nih.gov/sra>

## Format references (have these open in tabs)

- **FASTA format — NCBI quick reference**:
  <https://blast.ncbi.nlm.nih.gov/doc/blast-topics/>
- **FASTQ format — Wikipedia** (yes, really — the history section is the best one online):
  <https://en.wikipedia.org/wiki/FASTQ_format>
- **Phred quality scores — Ewing & Green 1998** (foundational, free PDF via *Genome Research*):
  <https://genome.cshlp.org/content/8/3/186.full>
- **NCBI sequence-format conventions** (description-line conventions across NCBI, GenBank, UniProt):
  <https://www.ncbi.nlm.nih.gov/genbank/fastaformat/>

## Tools you will install this week

- **Biopython 1.83** — `pip install biopython==1.83`. Tutorial linked above.
- **FastQC 0.12.1** — `conda install -c bioconda fastqc=0.12.1` (or download from Babraham).
- **seqkit 2.8** — a fast FASTA/FASTQ command-line utility we use for sanity checks: <https://bioinf.shenwei.me/seqkit/>
- **samtools 1.20** (only used in week 5 onwards, but installing it now smooths setup): <https://www.htslib.org/>

## Free books (chapter-level)

- **Bioinformatics Data Skills**, Vince Buffalo — chapters 3, 8, 10 (free):
  <https://vincebuffalo.com/bds/>
- **Biopython Tutorial and Cookbook** — full text linked above; chapter 5 is the canonical FASTQ-in-Python read:
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **Computational Genomics with R**, Akalin — chapter 7 on read QC (free):
  <https://compgenomr.github.io/book/>

## Public FASTQ datasets you will meet this week

- **1000 Genomes Project (phase 3)** — high-coverage Illumina FASTQ for population samples:
  <https://www.internationalgenome.org/data-portal/sample>
- **NA12878** — the most-sequenced human genome on Earth, a great QC sandbox:
  <https://www.internationalgenome.org/data-portal/sample/NA12878>
- **SARS-CoV-2 reference** (NC_045512.2) — small FASTA, great for FASTA-format drills:
  <https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2>
- **ENA / SRA** — search by accession (e.g. `ERR1019034`) and download FASTQ directly:
  <https://www.ebi.ac.uk/ena/browser/home>

## Videos (free, no signup)

- **Babraham Bioinformatics — FastQC walk-through** (45 min):
  <https://www.youtube.com/@BabrahamBioinf>
- **Galaxy Training Network — FASTQ quality control** (text + video):
  <https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html>
- **iBiology — Next-generation sequencing primer** (15 min, free):
  <https://www.ibiology.org/genetics-and-gene-regulation/next-generation-sequencing/>

## Open-source code to read this week

You can learn more from one hour reading other people's code than from three hours of tutorials. Pick one:

- **Biopython** — open `Bio/SeqIO/FastaIO.py` and `Bio/SeqIO/QualityIO.py` side by side:
  <https://github.com/biopython/biopython/tree/master/Bio/SeqIO>
- **seqkit** — Go, but the source is readable; the `seq` and `stats` subcommands are a tour of FASTQ in 300 lines:
  <https://github.com/shenwei356/seqkit>
- **FastQC** — Java; look at `uk.ac.babraham.FastQC.Modules` to see exactly what each report panel measures:
  <https://github.com/s-andrews/FastQC>

## FASTQ cheat sheet

Keep this open while you work the exercises.

| Line | Starts with | Contents |
|------|-------------|----------|
| 1 | `@` | Read identifier (and optional description) |
| 2 |     | The sequence — `{A, C, G, T, N}` |
| 3 | `+` | Separator (sometimes repeats the identifier; usually empty) |
| 4 |     | Quality string — one ASCII character per base, same length as line 2 |

## Phred quality cheat sheet

| Phred Q | P(error) | Plain English |
|--------:|---------:|---------------|
| 10 | 1 in 10 | 90% confident in the base |
| 20 | 1 in 100 | 99% confident |
| 30 | 1 in 1,000 | 99.9% confident — "Q30" is the common Illumina spec |
| 40 | 1 in 10,000 | 99.99% — typical for the best 25 bases of an Illumina run |
| 50 | 1 in 100,000 | PacBio HiFi / Sanger territory |
| 60 | 1 in 1,000,000 | Mostly only on consensus calls |

## Encoding cheat sheet

| Encoding | ASCII offset | Phred range | Era |
|----------|-------------:|------------:|-----|
| Sanger / Phred+33 | 33 | 0–93 | 1990s onward; the modern default |
| Solexa+64 | 64 | -5–62 | 2004–2006; pre-Illumina-acquisition |
| Illumina 1.3+ | 64 | 0–62 | 2009–2010 |
| Illumina 1.5+ | 64 | 3–62 (with quirky 'B' tag) | 2009–2011 |
| Illumina 1.8+ | 33 | 0–41 | 2011 onward — now the de facto standard |

If you see ASCII characters below `;` (semicolon, ASCII 59) in your quality string, the file is Phred+33. If the lowest character you see is in the range `@`–`B`, the file is most likely Phred+64. The exercises walk you through detecting this programmatically.

---

*If a link 404s, please open an issue so we can replace it.*
