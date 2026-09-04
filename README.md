# C10 · Crunching Bioinformatics

> A free, open-source **12-week bioinformatics track** for engineers (with biology curiosity) and biologists (with Python curiosity). From your first FASTA file to a published-quality variant-calling pipeline and a phylogenetic tree of your own SARS-CoV-2 samples. Built around public open-data sources and reproducibility.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Python · Biopython · R](https://img.shields.io/badge/stack-Python_·_Biopython_·_R-15803D.svg)](#stack)
[![Built in the open](https://img.shields.io/badge/built-in%20the%20open-15803D.svg)](https://github.com/CODE-CRUNCH-CLUB)

C10 sits at the rare intersection where Python, statistics, and biology meet. The bioinformatics field is structurally **open-source-first** (NCBI, Ensembl, EBI, Bioconductor are all free), which is why a Code Crunch track on it is feasible and why we treat reproducibility — pinned environments, public datasets, deterministic pipelines — as non-negotiable.

---

## Standards & equivalency

> C10 stands in for a university's undergraduate bioinformatics course.

**University equivalent.** Bioinformatics — `BSC 4934`, `CAP 4553`, `CS 4775`. Coverage: full.

C10 carries no credit, no transcript entry, no accreditation and no proctored exam. The equivalence is one of **content and skill**: everything an accredited section of that course teaches, taught here at the same depth or deeper, and assessed. What a registrar records is not something an open repository can give you.

| University outcome | Where this course teaches it | Depth |
| --- | --- | --- |
| Describe the molecular biology, the canonical sequence data formats, and the data-governance constraints that computational analysis of biological data operates under | [Week 01](curriculum/week-01-vocabulary-and-ethics/) | deeper |
| Read, write and manipulate biological sequence data programmatically, including quality-score encodings and the arithmetic behind them | [Week 02](curriculum/week-02-fasta-fastq-io/) | same |
| Derive, implement and analyse the dynamic-programming algorithms for global and local pairwise alignment, with substitution matrices and gap models | [Week 03](curriculum/week-03-pairwise-alignment/) | deeper |
| Explain heuristic sequence-database search and its scoring statistics, and interpret a result set by E-value and bit score | [Week 04](curriculum/week-04-blast-and-taxonomy/) | deeper |
| Map high-throughput sequencing reads to a reference genome and work with the resulting alignment formats | [Week 05](curriculum/week-05-read-alignment/) | deeper |
| Call and filter genetic variants from aligned reads, and evaluate the quality of a variant call set | [Week 06](curriculum/week-06-variant-calling/) | same |
| Quantify gene expression from RNA sequencing data, normalize it across samples and transcript lengths, and test for differential expression | quantification and normalization in [Week 07](curriculum/week-07-transcriptomics-and-rna-seq/); the differential-expression test is capstone Track 2 in [Week 12](curriculum/week-12-capstone-end-to-end-pipeline/) | same |
| Annotate variants for functional consequence and interpret them against population and clinical evidence | [Week 08](curriculum/week-08-variant-annotation-and-interpretation/) | deeper |
| Construct multiple sequence alignments and infer phylogenetic trees, with an assessment of branch support | [Week 09](curriculum/week-09-phylogenetics-and-msa/) | deeper |
| Assemble a genome de novo from sequencing reads and assess the quality of the assembly | [Week 10](curriculum/week-10-long-read-sequencing-and-assembly/) | deeper |
| Apply the methods of the course to a specialised genomics domain and interpret the results in their biological context | [Week 11](curriculum/week-11-cancer-genomics/) | deeper |
| Read the primary computational-biology literature and evaluate a published method's claims | [Week 01](curriculum/week-01-vocabulary-and-ethics/homework.md), problem 1 — read a real methods section — and again in the [Week 12](curriculum/week-12-capstone-end-to-end-pipeline/) write-up | same |
| Complete a substantial independent analysis on real data and report it | [Week 12](curriculum/week-12-capstone-end-to-end-pipeline/) | deeper |

Every row above points at a week that **assigns work** on that outcome — exercises, a challenge, six homework problems, a ten-question quiz and a mini-project — not merely a week that mentions it.

**The industry bar.** What an employer expects of somebody paid to do bioinformatics, and where this course makes the learner do it. Two rows say plainly what C10 does *not* ship.

| What the job expects | Where this course does it |
| --- | --- |
| Work lands as a commit in a repository you own, not a file on your desktop | The capstone deliverable is a tagged Git repository, not a notebook — [`curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/03-project-tracks-and-deposit.md`](curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/03-project-tracks-and-deposit.md) |
| You examine work you did not produce and form a judgement on it | C10 ships **no code-review assignment**. What it does instead, from Week 5 onward, is make the learner run a second, independently written implementation of the same step and account for every disagreement — [`curriculum/week-06-variant-calling/challenges/challenge-01-compare-callers.md`](curriculum/week-06-variant-calling/challenges/challenge-01-compare-callers.md), [`curriculum/week-11-cancer-genomics/challenges/challenge-01-strelka2-cross-check.md`](curriculum/week-11-cancer-genomics/challenges/challenge-01-strelka2-cross-check.md) |
| Tests exist, and the command to run them is written down | **One** unit requires a test suite — [`curriculum/week-08-variant-annotation-and-interpretation/challenges/challenge-01-implement-acmg-criteria.md`](curriculum/week-08-variant-annotation-and-interpretation/challenges/challenge-01-implement-acmg-criteria.md) asks for `test_acmg_classifier.py` with a `pytest` case per criterion. Everywhere else the check is an assertion block at the end of each exercise script, described in the exercise index — [`curriculum/week-02-fasta-fastq-io/exercises/README.md`](curriculum/week-02-fasta-fastq-io/exercises/README.md) |
| Failure is taught, not discovered | Each week closes its lecture notes on the failure modes of that step and the diagnostic signal each one produces, and each mini-project write-up carries a `Failure modes observed` section — [`curriculum/week-10-long-read-sequencing-and-assembly/lecture-notes/03-polishing-and-qc.md`](curriculum/week-10-long-read-sequencing-and-assembly/lecture-notes/03-polishing-and-qc.md), [`curriculum/week-05-read-alignment/mini-project/README.md`](curriculum/week-05-read-alignment/mini-project/README.md) |
| Dependencies are isolated and every version is pinned | Conda environments from Week 2 onward, and an `environment.yml` plus an explicit lockfile at the capstone — [`curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/02-conda-bioconda-and-containers.md`](curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/02-conda-bioconda-and-containers.md) |
| The analysis runs unattended as a pipeline, not as commands typed in order | [`curriculum/week-07-transcriptomics-and-rna-seq/challenges/challenge-01-build-a-mini-pipeline-with-snakemake.md`](curriculum/week-07-transcriptomics-and-rna-seq/challenges/challenge-01-build-a-mini-pipeline-with-snakemake.md), then the whole capstone |
| Every run records what produced it | A `run-info.json` — tool versions, seeds, accessions, commit — is required with every deliverable from Week 8 on — [`curriculum/week-09-phylogenetics-and-msa/mini-project/README.md`](curriculum/week-09-phylogenetics-and-msa/mini-project/README.md) |
| It runs from a clean clone by following the README | The capstone is graded on its reproducibility profile, and a re-run on a clean machine is part of the deliverable — [`curriculum/week-12-capstone-end-to-end-pipeline/README.md`](curriculum/week-12-capstone-end-to-end-pipeline/README.md) |
| A pipeline that re-runs itself on every push | **C10 ships no continuous-integration unit.** A CI workflow appears once, as a stretch goal at the end of Week 12 — [`curriculum/week-12-capstone-end-to-end-pipeline/README.md`](curriculum/week-12-capstone-end-to-end-pipeline/README.md) |

**Beyond both bars.** Clearing the two floors is entry, not success. Open any of these and check in under a minute.

| What we add | Which bar it beats | Where it lives |
| --- | --- | --- |
| Worked answers ship with the course rather than behind a deadline — a full solution set for Week 1 and a `SOLUTIONS.md` for Weeks 7 through 11. Weeks 2 to 6 are self-checking by design and say so in their exercise index | both | [`curriculum/week-09-phylogenetics-and-msa/exercises/SOLUTIONS.md`](curriculum/week-09-phylogenetics-and-msa/exercises/SOLUTIONS.md) |
| Every weekly quiz carries its own answer key at the bottom of the same file, so a learner can mark themselves the hour they sit it | both | [`curriculum/week-03-pairwise-alignment/quiz.md`](curriculum/week-03-pairwise-alignment/quiz.md) |
| Every method is taught by running two independently written implementations over the same input and requiring an account of where they disagree — bcftools against GATK, kallisto against Salmon against HISAT2, Mutect2 against Strelka2, Flye against Hifiasm | both | [`curriculum/week-06-variant-calling/challenges/challenge-01-compare-callers.md`](curriculum/week-06-variant-calling/challenges/challenge-01-compare-callers.md) |
| Data ethics is a taught and assessed unit in Week 1 — consent, re-identification, the IRB route, and the line between "associated with" and genetic determinism — not a paragraph in a course policy | both | [`curriculum/week-01-vocabulary-and-ethics/lecture-notes/02-data-ethics-and-public-data-sources.md`](curriculum/week-01-vocabulary-and-ethics/lecture-notes/02-data-ethics-and-public-data-sources.md) |
| Two weeks sit past the standard undergraduate outcome set: matched tumour–normal somatic calling with COSMIC signature decomposition, and long-read de novo assembly with polishing and graph inspection | university | [`curriculum/week-11-cancer-genomics/`](curriculum/week-11-cancer-genomics/) |
| The learner ends holding a deposited, DOI-citeable pipeline with a provenance record, not a grade only a registrar can see | both | [`curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/03-project-tracks-and-deposit.md`](curriculum/week-12-capstone-end-to-end-pipeline/lecture-notes/03-project-tracks-and-deposit.md) |

**Gaps we declare.** None against the bioinformatics outcome set. Two things C10 does not claim, and does not pretend to: it is not a statistics course — the negative-binomial model behind a differential-expression test is used and cited, not derived — and it is not a structural-biology course, so protein structure prediction and molecular docking are absent. Those belong elsewhere, and C10 does not claim them.

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

- **Biology / pre-med / pharma learner** with C1-equivalent Python comfort who wants to work with -omics data.
- **Software engineer at a biotech, pharma, or academic lab** preparing for a bioinformatics-adjacent role.
- **Quantitative-biology / computational-biology grad learner** seeking a structured open-source curriculum.
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
