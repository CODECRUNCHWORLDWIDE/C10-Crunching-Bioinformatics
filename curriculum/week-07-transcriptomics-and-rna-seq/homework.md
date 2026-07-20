# Week 7 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-07/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Trim two more SRA samples with fastp

**Problem statement.** Run `fastp` on the other two mini-project samples (`SRR453566` glucose-rep-1 and `SRR453567` glucose-rep-2) with the same flags as Exercise 1. For each sample, parse the JSON report and record the seven QC fields (`reads_in`, `reads_out`, `pct_retained`, `pct_q30`, `adapter_pct`, `duplication_rate`, `insert_size_peak`). Build a 3-row table comparing the two new glucose samples to the galactose sample from Exercise 1.

Answer:

1. Which of the three samples has the highest duplication rate? Does that align with a known property of the sample (e.g. PCR amplification, condition-specific complexity)?
2. Which has the lowest `pct_q30`? Does it correlate with sequencing date / position in the run / batch (if you have metadata)?
3. Are all three samples in the "Healthy" band? If any sample is in "Marginal" or "Bad," note it and decide whether to include or drop.

**Acceptance criteria.**

- `homework/p1_trim.sh` runs the two fastp jobs end to end.
- `qc/SRR453566.fastp.html`, `qc/SRR453566.fastp.json`, `qc/SRR453567.fastp.html`, `qc/SRR453567.fastp.json` exist.
- `homework/notes/p1-fastp-comparison.md` contains a 3-row table and three numbered answers.
- Commit message like `p1: fastp on glucose replicates, both healthy`.

**Hint.** The fastp command is identical to Exercise 1 except for the sample name. Use a loop: `for s in SRR453566 SRR453567; do fastp ... ; done`. Or write the loop into a small bash script.

**Estimated time.** 45 minutes.

---

## Problem 2 — kallisto quant on all three samples

**Problem statement.** With the same kallisto index from Exercise 2, run `kallisto quant -b 100 -t 4` on all three samples. For each sample, record:

1. `n_processed` (total trimmed reads).
2. `n_pseudoaligned` (reads that matched the transcriptome).
3. `p_pseudoaligned` (alignment rate, %).
4. `p_unique` (% of pseudoaligned reads with a single-transcript compatibility class).
5. Mean fragment length (from `frag_length_mean` in `run_info.json`).
6. EM rounds to convergence.

Build a per-sample comparison table.

Answer:

1. Are the alignment rates similar across samples (within 2-3 percentage points)?
2. Is one sample much higher in `p_unique` than the others? Why might that be?
3. Are the fragment-length estimates similar (within 20 bp of each other)?

**Acceptance criteria.**

- `homework/p2_quant.sh` runs the three kallisto jobs.
- `quant/SRR453566/`, `quant/SRR453567/`, `quant/SRR453568/` each contain `abundance.tsv`, `abundance.h5`, `run_info.json`.
- `homework/notes/p2-kallisto-comparison.md` contains a 3-row table and three numbered answers.
- Commit message like `p2: kallisto quant on three samples, alignment rates 94-95%`.

**Hint.** kallisto's command line is identical across the three samples; only the input FASTQs and output directory change. Use a loop or a Snakemake rule.

**Estimated time.** 30 minutes (most of which is the EM step waiting).

---

## Problem 3 — Build a 3-sample gene-level counts matrix

**Problem statement.** Write `homework/p3_counts_matrix.py` that:

1. Loads each sample's `abundance.tsv`.
2. Aggregates per-transcript counts to per-gene counts using the Ensembl release 110 yeast GTF (or just rename target_id to gene_id; for yeast they are nearly identical).
3. Joins the three per-sample frames on `gene_id` and produces a single TSV with columns `gene_id, S1, S2, S3` (the three sample names of your choice — suggest `glucose_rep1, glucose_rep2, galactose_rep1`).
4. Rounds the per-cell counts to integers (DESeq2 wants ints).

Answer in `homework/notes/p3-counts-matrix.md`:

1. How many rows (genes) does the matrix have?
2. How many genes have a count ≥ 10 in all three samples? (Filtering threshold for downstream DE.)
3. What is the per-sample column sum? Should it match the kallisto `n_pseudoaligned` for each sample?

**Acceptance criteria.**

- `homework/p3_counts_matrix.py` runs and writes `homework/counts/all_samples.counts.tsv`.
- The TSV has ~6,975 rows + 1 header line.
- `homework/notes/p3-counts-matrix.md` contains three numbered answers.
- Commit message like `p3: 3-sample counts matrix, 6975 genes, 3812 expressed in all`.

**Hint.** For yeast, you can skip the GTF parsing if you treat `target_id` as `gene_id` directly (most yeast transcripts have IDs identical to their genes). For human/mouse, you would need `tximport`-style aggregation.

**Estimated time.** 45 minutes.

---

## Problem 4 — Compute TPM by hand for the 3-sample matrix

**Problem statement.** Write `homework/p4_tpm_matrix.py` that takes the 3-sample counts matrix from Problem 3 and:

1. Loads the per-gene effective length from any one of the `abundance.tsv` files (effective length is per-transcript; for yeast where transcript ≈ gene, this works directly).
2. Computes TPM per sample using the Lecture 3 formula.
3. Verifies that each sample column of the TPM matrix sums to 10^6 (within float tolerance).
4. Saves `homework/counts/all_samples.tpm.tsv`.

Answer in `homework/notes/p4-tpm-matrix.md`:

1. What is the TPM of GAL1 (`YBR020W`) in each of the three samples?
2. What is the fold change of GAL1 (galactose vs mean of two glucose samples)?
3. Pick two other genes from the GAL regulon (GAL7, GAL10, GAL2, GAL3, GCY1) and report their TPMs.

**Acceptance criteria.**

- `homework/p4_tpm_matrix.py` runs, writes the TPM TSV, and prints the per-sample sums (each ~10^6).
- `homework/notes/p4-tpm-matrix.md` contains three numbered answers with specific TPM values.
- The galactose-vs-glucose fold change for GAL1 should be ≥ 100x.
- Commit message like `p4: TPM matrix; GAL1 fold change galactose/glucose = 387x`.

**Hint.** The TPM denominator is the per-sample sum of `count / eff_length`. Compute it once per sample, then divide.

**Estimated time.** 45 minutes.

---

## Problem 5 — Compare kallisto and featureCounts on one sample

**Problem statement.** Pick one of your three samples (`SRR453568` is fine; you already have the kallisto output). Run the alignment-based pipeline on it: `hisat2-build`, `hisat2 | samtools sort`, `samtools index`, `featureCounts -t exon -g gene_id -p --countReadPairs -s 2`. Build a per-gene comparison table with columns `gene_id, counts_kallisto, counts_fc, tpm_kallisto, tpm_fc, log2_ratio`. Compute the Pearson correlation of `log2(TPM+1)` on the well-expressed (TPM > 10) subset.

Answer in `homework/notes/p5-kallisto-vs-fc.md`:

1. What is the Pearson r? Is it ≥ 0.97? (If much less, debug the strandedness flag.)
2. What are the top 5 most discordant genes (by `|log2_ratio|`)? Are they rRNA or paralog-family genes?
3. Pick one discordant gene and explain in 50 words why the two methods disagree on it.

**Acceptance criteria.**

- `homework/p5_compare.py` runs and prints the Pearson r.
- `homework/notes/p5-kallisto-vs-fc.md` contains three numbered answers with specific numbers and one named case study.
- Pearson r ≥ 0.97.
- Commit message like `p5: kallisto vs featureCounts on SRR453568, r = 0.987`.

**Hint.** featureCounts requires a sorted BAM. HISAT2 → `samtools sort` → `samtools index` is the standard chain. The `-s 2` flag (reverse-stranded) is critical for Illumina dUTP libraries; the wrong setting will halve your counts. If you are unsure, run `infer_experiment.py` from RSeQC (Wang et al. 2012) or check the SRA metadata.

**Estimated time.** 90 minutes (most in HISAT2 alignment).

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300-400 word reflection at `homework/notes/week-07-reflection.md` answering:

1. Before Week 7, what did you think "RNA-seq" was? What is it actually? Pick one step of the molecular workflow you found most surprising and say why.
2. The first time you ran kallisto, how did the speed compare to what you expected? After Week 7, when does kallisto's speed matter vs when does the BAM-producing classical pipeline matter? Give one concrete example of each.
3. The TPM/CPM/RPKM/FPKM normalizations look interchangeable at first glance. After Week 7, what is each one actually for? Pick the one you find most counterintuitive (probably "TPM but not CPM divides by effective length") and explain in your own words what the division does.
4. The mini-project produces a 3-sample counts matrix and asks for a biological interpretation. What is the difference between a counts matrix and a differential-expression table? Why does the second exist? What kind of question can you answer with the DE table that you cannot answer with the counts matrix?

**Acceptance criteria.**

- File exists, 300-400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The mistakes you note here are what you will re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 45 min |
| 2 | 30 min |
| 3 | 45 min |
| 4 | 45 min |
| 5 | 1 h 30 min |
| 6 | 30 min |
| **Total** | **~4 h 45 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
