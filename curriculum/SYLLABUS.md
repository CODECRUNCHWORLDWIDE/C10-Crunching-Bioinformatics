# C10 · Crunching Bioinformatics — Syllabus

**12 weeks · ~36 hrs/week intensive (or scaled) · C1 + biology vocabulary → bioinformatics analyst**

---

## Program at a glance

| Phase | Weeks | Outcome |
|-------|-------|---------|
| **Phase 1 — Sequences** | 01 – 03 | FASTA, FASTQ, alignment, BLAST |
| **Phase 2 — Genomes** | 04 – 06 | Read alignment, variants, annotation |
| **Phase 3 — Expression & Phylogenetics** | 07 – 09 | RNA-seq, differential expression, trees |
| **Phase 4 — Reproducibility & Capstone** | 10 – 12 | Snakemake, real research question |

---

## Weekly breakdown

**Week 1 — Vocabulary & ethics.** What bioinformatics actually is. The central dogma in 90 minutes. Data-ethics module (consent, re-identification, family DNA). Open data sources.

- *Mini-project:* A 1-page glossary in your own words + a public-data inventory you'll use through the course.

**Week 2 — FASTA, FASTQ, IO.** Parsing with Biopython. Quality scores. FastQC. Common gotchas.

- *Mini-project:* A QC report on a real public dataset (1000 Genomes subset).

**Week 3 — Pairwise alignment.** Needleman-Wunsch, Smith-Waterman, by hand and in NumPy. Substitution matrices.

- *Mini-project:* A pure-NumPy Smith-Waterman implementation, benchmarked vs Biopython.

**Week 4 — BLAST and taxonomy.** Running BLAST locally and via NCBI. Interpreting E-values. Building a classifier.

- *Mini-project:* A BLAST-driven taxonomy classifier for unknown sequences.

**Week 5 — Read alignment.** `bwa`, `minimap2`, SAM/BAM. Coverage. Duplicate marking.

- *Mini-project:* Align a public small-genome dataset; produce coverage plots.

**Week 6 — Variant calling.** `bcftools`, GATK basics, VCF parsing, hard filters, annotation with VEP.

- *Mini-project:* End-to-end variant call on a small public dataset with QC.

**Week 7 — Transcriptomics intro.** RNA-seq design. Counting reads. Normalization (TPM, FPKM, DESeq2's median-of-ratios).

- *Mini-project:* From counts to a clean expression matrix.

**Week 8 — Differential expression.** DESeq2 (R), edgeR. Volcano plots. Multiple-testing correction.

- *Mini-project:* A differential-expression analysis with FDR control and a volcano plot.

**Week 9 — Phylogenetics.** Multiple sequence alignment (MAFFT). Tree inference (IQ-TREE, RAxML). Bootstrap. Interpreting branches.

- *Mini-project:* A phylogenetic tree of public viral sequences with bootstrap support.

**Week 10 — Visualization & communication.** `ggplot2` (R) and `seaborn`. Genome browsers (IGV). What a publication-quality figure looks like.

- *Mini-project:* Re-create one figure from a published paper using public data.

**Week 11 — Reproducible pipelines.** Snakemake, Nextflow. Conda environment files. Containerization. Versioning data.

- *Mini-project:* Convert the Week-6 variant-calling pipeline to Snakemake, with full reproducibility.

**Week 12 — Capstone.** A self-directed research question on public data, end-to-end, with a 4-page write-up.

- *Capstone:* Public repo with a one-command pipeline + a figure + a written interpretation.

---

## Weekly load

| Component | hrs/wk |
|-----------|------:|
| Lectures / readings | 6 |
| Hands-on exercises | 8 |
| Pipeline / analysis runs | 4 |
| Quiz | 3 |
| Homework | 6 |
| Mini-project | 7 |
| Self-study | 2 |
| **Total** | **36** |

---

## License

GPL-3.0.
