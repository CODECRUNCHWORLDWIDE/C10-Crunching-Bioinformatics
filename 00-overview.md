# C10 · Crunching Bioinformatics

> A free, open-source **12-week bioinformatics track** for engineers (with biology curiosity) and biologists (with Python curiosity). From your first FASTA file to a published-quality variant-calling pipeline and a phylogenetic tree of your own SARS-CoV-2 samples. Built around public open-data sources and reproducibility.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Python · Biopython · R](https://img.shields.io/badge/stack-Python_·_Biopython_·_R-15803D.svg)](#stack)
[![Built in the open](https://img.shields.io/badge/built-in%20the%20open-15803D.svg)](https://github.com/CODECRUNCHWORLDWIDE)

C10 sits at the rare intersection where Python, statistics, and biology meet. The bioinformatics field is structurally **open-source-first** (NCBI, Ensembl, EBI, Bioconductor are all free), which is why a Code Crunch track on it is feasible and why we treat reproducibility — pinned environments, public datasets, deterministic pipelines — as non-negotiable.

---

## Pathway summary

- **Full-time:** 12 weeks · ~36 hrs/week · ~432 hours
- **Working-engineer pace:** 6 months · ~18 hrs/week
- **Evening / undergraduate pace:** 1 year · ~9 hrs/week

See [`SYLLABUS.md`](SYLLABUS.md) for the full 12-week breakdown.

---

## What you will be able to do at the end of 12 weeks

- **Parse and manipulate** biological sequence files: FASTA, FASTQ, SAM/BAM, VCF, GFF/GTF.
- **Use Biopython** fluently — `SeqIO`, `Entrez`, alignment, BLAST automation.
- **Run a small variant-calling pipeline** end-to-end: read → align → variant → annotate, with quality-control at every step.
- **Build a phylogenetic tree** from your own data and interpret its branches with statistical support.
- **Perform a basic transcriptomic analysis** — differential expression, basic clustering, volcano plots.
- **Use R / Bioconductor** for the parts where it dominates (DESeq2, edgeR, complex visualization).
- **Set up a reproducible pipeline** with `snakemake` or `nextflow` so a reviewer can re-run your analysis with one command.
- **Read a paper from Nature Methods / Bioinformatics** and reproduce one of its supplementary figures from public data.

---

## Who this is for

- **Biology / pre-med / pharma student** with C1-equivalent Python comfort who wants to work with -omics data.
- **Software engineer at a biotech, pharma, or academic lab** preparing for a bioinformatics-adjacent role.
- **Quantitative-biology / computational-biology grad student** seeking a structured open-source curriculum.
- **Undergraduate researcher** who wants to add bioinformatics to a CV honestly.

Not for: people who want a deep biochemistry course (this is the *computational* side; we explain biology only as needed) or pure ML researchers (see [C5](../C5-CRUNCH-AI-DATA-SCIENCE/) and the upcoming C23 Crunch Agents for that).

---

## Prerequisites

- **C1 Weeks 1–11** completed (Python, file IO, basic data structures, pandas).
- **High-school biology** — DNA, RNA, proteins, the central dogma. We'll re-teach as needed but assume the vocabulary.
- A computer with ≥16 GB RAM (or willingness to use a free cloud VM for the genome-scale weeks).
- Time. Bioinformatics analyses are *slow*. Plan for waits.

---

## What you ship

By the end of the 12 weeks, your `crunch-bio-portfolio-<yourhandle>` GitHub repo contains:

1. A **FASTA / FASTQ exploration notebook** with quality plots (Week 2).
2. A **pairwise sequence alignment tool** built from scratch in pure NumPy (Week 3).
3. A **BLAST-driven taxonomy classifier** for unknown sequences (Week 4).
4. A **variant-calling mini-pipeline** on a small public dataset (Week 6).
5. A **gene-expression analysis** from public RNA-seq data with volcano plots (Week 8).
6. A **phylogenetic tree** built from public viral sequences with bootstrap support (Week 9).
7. A **reproducible Snakemake pipeline** for the Week-6 variant calling (Week 11).
8. **Capstone:** a 4-page write-up on a real biological question you investigated end-to-end using public data, with a reproducible repo and a one-figure visualization (Week 12).

---

## Tools (all free, all open-source)

| Tool | Role |
|------|------|
| **Python 3.11+** | Primary language |
| **Biopython** | Sequence parsing, BLAST, Entrez |
| **pandas · NumPy · matplotlib · seaborn** | Data |
| **R + Bioconductor (DESeq2, edgeR, ape)** | Specific analyses |
| **samtools · bcftools · bwa · minimap2** | Variant-calling chain |
| **IGV (free desktop)** | Genome browsing |
| **MAFFT / MUSCLE / RAxML / IQ-TREE** | Alignment & phylogenetics |
| **Snakemake / Nextflow** | Reproducible pipelines |
| **Conda / Mamba** | Reproducible environments |
| **Public data: NCBI, Ensembl, UniProt, EBI, GISAID, GTEx, TCGA** | Datasets |

---

## Notes on data ethics

Many bioinformatics datasets contain human genotypes. Even "de-identified" data can sometimes be re-identified. We address this in Week 1 with a small ethics module. Two rules across the whole track:

1. **Use public, consent-cleared datasets** for all exercises (1000 Genomes, GTEx-public, GISAID, etc.).
2. **Do not analyze a friend or family member's DNA** for the course, even if they "say it's fine." There are consent and harm dynamics in genetic information that this curriculum does not equip you to handle. If you have a research interest in your family's genetics, do it under IRB-supervised protocol at a university.

---

## Next track after C10

- **[C5 · Crunch AI / Data Science](../C5-CRUNCH-AI-DATA-SCIENCE/)** — for deeper ML on biological data.
- **[C17 · Crunch Pro Python Advanced](../C17-CRUNCH-PRO-PYTHON-ADVANCED/)** — for the performance-of-large-pipelines side.
- **[C15 · Crunch DevOps](../C15-CRUNCH-DEVOPS/)** — for cluster-scale pipeline execution.

---

## License

GPL-3.0.

---

*C10 is part of the Code Crunch open-source curriculum.* [Master catalog ↗](../MASTER-CURRICULUM.md) · [Brand family ↗](../../assets/brand/BRAND-FAMILY.md)


---

<!-- CCWW:AUTO-INDEX:START — generated by scripts/restructure_course_repos.py; edit ABOVE this marker -->

## Course at a glance

| Section | Count |
| --- | --- |
| Curriculum entries | 13 |
| Projects | 0 |
| Past sessions | 1 |

## Curriculum

- [SYLLABUS](curriculum/SYLLABUS.md)
- [week 01 vocabulary and ethics](./curriculum/week-01-vocabulary-and-ethics/00-overview.md)
- [week 02 fasta fastq io](./curriculum/week-02-fasta-fastq-io/00-overview.md)
- [week 03 pairwise alignment](./curriculum/week-03-pairwise-alignment/00-overview.md)
- [week 04 blast and taxonomy](./curriculum/week-04-blast-and-taxonomy/00-overview.md)
- [week 05 read alignment](./curriculum/week-05-read-alignment/00-overview.md)
- [week 06 variant calling](./curriculum/week-06-variant-calling/00-overview.md)
- [week 07 transcriptomics and rna seq](./curriculum/week-07-transcriptomics-and-rna-seq/00-overview.md)
- [week 08 variant annotation and interpretation](./curriculum/week-08-variant-annotation-and-interpretation/00-overview.md)
- [week 09 phylogenetics and msa](./curriculum/week-09-phylogenetics-and-msa/00-overview.md)
- [week 10 long read sequencing and assembly](./curriculum/week-10-long-read-sequencing-and-assembly/00-overview.md)
- [week 11 cancer genomics](./curriculum/week-11-cancer-genomics/00-overview.md)
- [week 12 capstone end to end pipeline](./curriculum/week-12-capstone-end-to-end-pipeline/00-overview.md)

## In this course

- **Community** — [community/](community/)
- **Curriculum** — [curriculum/](curriculum/)
- **Projects** — [projects/](projects/)
- **Resources** — [resources/](resources/)
- **Past sessions** — [past-sessions/](past-sessions/)

<!-- CCWW:AUTO-INDEX:END -->
