# Mini-Project — Yeast RNA-seq counts matrix and PCA

> Build a reproducible RNA-seq pipeline that takes three yeast paired-end FASTQ samples (two glucose replicates and one galactose replicate), trims them with fastp, pseudoaligns each with kallisto, aggregates per-transcript to per-gene counts, computes log-TPM, runs PCA, and writes up a one-page biological interpretation. End with: a 3-sample counts matrix, a 3-sample TPM matrix, a PCA scatter, and a 600-800 word write-up that defends every pipeline parameter and names the biology of the result.

This is the C10 mini-project that produces a **methods-section-quality counts matrix with measured QC**, not just a single-sample demonstration. By the end of it you will have a `quantify_three_samples.py` script and a `run.sh` wrapper you can point a recruiter at, a results directory with the matrix, TPM, PCA plot, and per-sample QC, and a write-up that defends each pipeline parameter and names the biological signal you observed.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-07/mini-project/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 7 section)
└── week-07/
    └── mini-project/
        ├── README.md               one-page report (~600-800 words)
        ├── run.sh                  one-command reproduction script
        ├── env.yml                 conda environment file pinning all tool versions
        ├── data/
        │   ├── sce.cdna.fa.gz      Ensembl release 110 yeast cDNA
        │   ├── sce.gtf             Ensembl release 110 yeast annotation
        │   └── samples.tsv         per-sample metadata (SRA, condition, replicate)
        ├── quantify_three_samples.py   the orchestration script
        ├── build_matrix.py         per-transcript -> per-gene aggregation
        ├── pca.py                  log-TPM PCA + plot
        ├── raw/                    raw FASTQs (gitignored)
        ├── trim/                   fastp-trimmed FASTQs (gitignored)
        ├── qc/                     fastp HTML/JSON (committed)
        ├── index/                  kallisto index (gitignored, fast to rebuild)
        ├── quant/                  per-sample kallisto output (abundance.tsv committed)
        └── results/
            ├── counts_matrix.tsv      3-sample integer counts
            ├── tpm_matrix.tsv         3-sample TPM (each column sums to 1e6)
            ├── log_tpm_matrix.tsv     log2(TPM+1) for PCA / heatmap input
            ├── pca_scores.tsv         per-sample PC1, PC2 scores
            ├── pca_loadings.tsv       per-gene PC1, PC2 loadings (top 50 each)
            ├── pca_scatter.png        2D PCA scatter, samples colored by condition
            ├── per_sample_qc.tsv      per-sample fastp + kallisto QC summary
            └── gal_regulon.tsv        GAL1/2/3/7/10/GCY1 expression in each sample
```

By the end you will have a clean, reproducible Week 7 directory you can point a recruiter at — and `quantify_three_samples.py` is the kind of pipeline that opens conversations with working bioinformaticians and biotech shops.

---

## The dataset

You will work with three publicly available *Saccharomyces cerevisiae* RNA-seq samples from Gierlinski et al. 2015 (*Bioinformatics* 31:3625), accessible through the NCBI SRA:

| Sample handle      | SRA accession | Carbon source | Replicate |
|--------------------|---------------|---------------|-----------|
| `glucose_rep1`     | `SRR453566`   | Glucose       | 1         |
| `glucose_rep2`     | `SRR453567`   | Glucose       | 2         |
| `galactose_rep1`   | `SRR453568`   | Galactose     | 1         |

Each sample is ~3 M paired-end 75 bp reads, polyA-selected, sequenced on Illumina HiSeq 2000. Total download: ~750 MB across the three samples.

Expected biology:

- The **GAL regulon** (`GAL1, GAL2, GAL3, GAL7, GAL10, GCY1`) is induced ~100-500-fold on galactose vs glucose. This is the textbook positive control for yeast RNA-seq.
- The two glucose replicates should cluster tightly in PCA (PC1 distance < PC2 distance to the galactose sample).
- The galactose sample should be well separated from the two glucose samples on PC1.
- A handful of other carbon-source-related genes will move (HXT family hexose transporters, MIG1-target genes), but the GAL regulon is the loudest signal.

If your PCA does not show this pattern, debug *before* writing the report. The "Common pitfalls" section below is the first place to look.

### Why yeast and not human

*Saccharomyces cerevisiae* is the canonical small-eukaryote RNA-seq benchmark for the same three reasons as Week 5/6's *E. coli*: small genome and transcriptome (~12 Mb, ~7,000 transcripts) so the full pipeline fits on a laptop and runs in under 10 minutes per sample, well-curated reference annotations from Ensembl and SGD, characteristic experimentally validated regulons (GAL, ZAP1, GCN4) that let you verify the pipeline works. Week 8 moves to differential expression (the formal DESeq2 / edgeR analysis on this same matrix); Week 7's mini-project produces the input.

---

## Rules

- **You may** use fastp 0.23.4, kallisto 0.50.1, salmon 1.10.2 (optional for comparison), hisat2 2.2.1 (optional for comparison), subread 2.0.6 (optional for comparison), samtools 1.19, pysam 0.22, Biopython 1.83, pandas 2.2, numpy 1.26, scikit-learn 1.4, matplotlib 3.8, and the standard library.
- **You may** consult Lectures 1, 2, 3, the kallisto paper (Bray et al. 2016), the Salmon paper (Patro et al. 2017), the fastp paper (Chen et al. 2018), the "Hitchhiker's Guide" (Conesa et al. 2016), the GFF3 spec, and your Week 7 exercises and challenges.
- **You may NOT** copy a pre-written RNA-seq pipeline from the internet. The point is to *build* the pipeline. Reading the `nf-core/rnaseq` docs for inspiration is fine; copy-pasting a `nf-core/rnaseq` `main.nf` is not.
- **You must** cache the kallisto index to disk. A second run of `bash run.sh` should not re-download the cDNA FASTA or re-build the index if they are still on disk.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file. The only network access on a fresh run should be the initial Ensembl cDNA + GTF + SRA downloads.

---

## Acceptance criteria

- [ ] `mini-project/quantify_three_samples.py` exports a function `quantify(samples_tsv: Path, ref_dir: Path, out_dir: Path) -> Path` that runs the full FASTQ → counts matrix pipeline and returns the path to the counts matrix.
- [ ] The pipeline implements **six stages**:
  1. Download/locate raw FASTQ pair per sample (or use existing).
  2. `fastp -i R1.fq.gz -I R2.fq.gz -o R1.trim.fq.gz -O R2.trim.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 36 --trim_poly_g -h qc/<sample>.fastp.html -j qc/<sample>.fastp.json -w 4`.
  3. `kallisto index -i index/sce.idx -k 31 ref/sce.cdna.fa.gz` (skip if `sce.idx` exists).
  4. `kallisto quant -i index/sce.idx -o quant/<sample>/ -t 4 -b 100 trim/<sample>_1.trim.fq.gz trim/<sample>_2.trim.fq.gz`.
  5. Aggregate per-transcript counts to per-gene counts (for yeast, mostly a rename); join across samples; write `counts_matrix.tsv`.
  6. Compute TPM and log2(TPM+1); run PCA; save scores, loadings, and a scatter PNG.
- [ ] `build_matrix.py` produces three TSVs:
  - `results/counts_matrix.tsv` — integer counts, one row per gene, three columns per sample.
  - `results/tpm_matrix.tsv` — TPM per sample (each column sums to ~10^6, verify in code).
  - `results/log_tpm_matrix.tsv` — log2(TPM+1) for visualization.
- [ ] `pca.py` produces:
  - `results/pca_scores.tsv` — per-sample PC1, PC2 scores.
  - `results/pca_loadings.tsv` — top 50 PC1 loadings and top 50 PC2 loadings (which genes drive each PC).
  - `results/pca_scatter.png` — 2D scatter, sample names labeled, colored by condition.
- [ ] `results/per_sample_qc.tsv` aggregates per-sample QC across the three samples (fastp + kallisto fields).
- [ ] `results/gal_regulon.tsv` reports the GAL1/2/3/7/10/GCY1 TPM in each sample.
- [ ] `mini-project/README.md` is a one-page (~600-800 word) report containing:
  - One-sentence description of the dataset, the reference, and the pipeline's six stages.
  - Methods section in C10 voice: every tool pinned ("fastp 0.23.4", "kallisto 0.50.1", "Ensembl release 110"), every parameter explicit (`fastp --qualified_quality_phred 20`, `kallisto quant -b 100`, `--length_required 36`).
  - Results section in C10 voice: per-sample numbers (`n_processed`, `n_pseudoaligned`, `p_pseudoaligned`, `est_counts` of the top 5 genes), the PCA result ("PC1 separates carbon source, explains 78% variance"), the top 10 PC1 loadings ("GAL1, GAL7, GAL10, MIG1, ..."), and the GAL regulon table.
  - Discussion section: 100-200 words on biological interpretation. What does the PCA tell you? What does the GAL regulon induction confirm? Are there any unexpected loadings worth flagging for follow-up?
- [ ] `run.sh` is a single bash script that, given a fresh checkout + `conda env create -f env.yml`, reproduces the entire pipeline from scratch in under 20 minutes (with SRA downloads being the slowest step).
- [ ] The repo is **public** and at least one classmate or instructor has been added as a collaborator.

---

## Suggested approach (rough timeline)

### Thursday (3 hours)

1. (15 min) `git clone`, set up `mini-project/` directory.
2. (30 min) Write `env.yml` with pinned tool versions; create conda env; verify each tool is on the PATH.
3. (60 min) Download the three SRA samples via `prefetch` + `fasterq-dump`. This is the slow step. Run in background while you work on the next item.
4. (45 min) Download Ensembl release 110 yeast cDNA FASTA and GTF; build the kallisto index. Verify with `kallisto inspect`.
5. (30 min) Write the fastp wrapper in `quantify_three_samples.py`. Run on one sample to verify the flag set works.

### Friday (3 hours)

1. (60 min) Run fastp on all three samples; verify QC reports look healthy.
2. (60 min) Run kallisto quant on all three samples (~30 seconds each + 100 bootstrap each ≈ 4 minutes total); verify `p_pseudoaligned ≥ 90%` for each.
3. (60 min) Write `build_matrix.py` to aggregate the three `abundance.tsv` into the counts matrix, TPM matrix, and log-TPM matrix. Verify TPM column sums = 10^6 ± 1.

### Saturday (2 hours)

1. (45 min) Write `pca.py` using scikit-learn `PCA(n_components=2)`. Filter to genes with `CPM > 1 in ≥ 2 samples`. Compute scores and loadings.
2. (30 min) Make the matplotlib scatter plot. Color points by condition, label each with the sample name.
3. (45 min) Write the README. Lead with the dataset and methods; then the results (with specific numbers); then the discussion paragraph.

---

## Methods recipe (the pipeline you must reproduce in `run.sh`)

```bash
#!/usr/bin/env bash
set -euo pipefail

# Stage 0: data sources.
mkdir -p data trim qc index quant results
curl -sLo data/sce.cdna.fa.gz \
    http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
curl -sLo data/sce.gtf.gz \
    http://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
gunzip -kf data/sce.gtf.gz

# Stage 1: per-sample fastp + kallisto quant.
for sra in SRR453566 SRR453567 SRR453568; do
    if [ ! -f raw/${sra}_1.fastq.gz ]; then
        prefetch $sra -O raw/${sra}_sra/
        fasterq-dump raw/${sra}_sra/${sra}/${sra}.sra --split-files -O raw/
        gzip raw/${sra}_1.fastq raw/${sra}_2.fastq
        rm -rf raw/${sra}_sra/
    fi
    fastp \
        -i raw/${sra}_1.fastq.gz -I raw/${sra}_2.fastq.gz \
        -o trim/${sra}_1.trim.fq.gz -O trim/${sra}_2.trim.fq.gz \
        --detect_adapter_for_pe --qualified_quality_phred 20 \
        --length_required 36 --trim_poly_g \
        -h qc/${sra}.fastp.html -j qc/${sra}.fastp.json -w 4
done

# Stage 2: kallisto index.
if [ ! -f index/sce.idx ]; then
    kallisto index -i index/sce.idx -k 31 data/sce.cdna.fa.gz
fi

# Stage 3: per-sample kallisto quant.
for sra in SRR453566 SRR453567 SRR453568; do
    if [ ! -f quant/${sra}/abundance.tsv ]; then
        kallisto quant -i index/sce.idx -o quant/${sra}/ -t 4 -b 100 \
            trim/${sra}_1.trim.fq.gz trim/${sra}_2.trim.fq.gz
    fi
done

# Stage 4: build matrices.
python build_matrix.py --quant-dir quant --samples data/samples.tsv --out-dir results

# Stage 5: PCA.
python pca.py --tpm results/log_tpm_matrix.tsv --samples data/samples.tsv --out-dir results

echo "Done. Counts matrix at results/counts_matrix.tsv. PCA scatter at results/pca_scatter.png."
```

---

## Expected results (your numbers should be within ~5% of these)

### Per-sample QC

| Sample          | reads_in    | reads_out   | pct_retained | pct_q30 | n_pseudoaligned | p_pseudoaligned |
|-----------------|------------:|------------:|-------------:|--------:|----------------:|----------------:|
| glucose_rep1    | 6,400,000   | 6,144,824   | 96.0%        | 93.5%   | 2,914,012       | 94.84%          |
| glucose_rep2    | 6,200,000   | 5,970,123   | 96.3%        | 93.1%   | 2,832,401       | 94.92%          |
| galactose_rep1  | 6,400,000   | 6,144,824   | 96.0%        | 93.5%   | 2,914,012       | 94.84%          |

### Counts matrix

- Shape: ~6,975 rows (one per yeast transcript), 3 columns (one per sample).
- Per-sample column sum: ~2.8-2.9 M (matches kallisto's `n_pseudoaligned` ± small rounding).
- Genes with count ≥ 10 in all three: ~3,800.

### PCA (on log2(TPM+1), filtered to CPM > 1 in ≥ 2 samples)

- PC1 explains ~75-80% of variance. PC1 separates galactose from glucose.
- PC2 explains ~10-12% of variance. PC2 captures the glucose-rep1-vs-rep2 batch difference.
- Top 10 PC1 loadings (positive end, galactose): `GAL1, GAL7, GAL10, GAL2, GAL3, GCY1, MTH1, IMA1, IMA5, FBP1`. The first six are the GAL regulon; the rest are carbon-source-responsive non-GAL genes.
- Top 10 PC1 loadings (negative end, glucose): `HXT1, HXT3, HXT4, MIG2, RGT1`. Hexose transporters and glucose-repressed regulators.

### GAL regulon expression

| Gene         | glucose_rep1 (TPM) | glucose_rep2 (TPM) | galactose_rep1 (TPM) | log2 FC (gal/glu_mean) |
|--------------|-------------------:|-------------------:|---------------------:|----------------------:|
| GAL1 (YBR020W)| 5.2               | 6.1                | 2,894.4              | 8.97                  |
| GAL2 (YLR081W)| 12.4              | 11.8               | 1,542.3              | 7.04                  |
| GAL3 (YDR009W)| 28.5              | 31.2               | 872.2                | 4.79                  |
| GAL7 (YBR018C)| 8.3               | 9.1                | 2,547.8              | 8.21                  |
| GAL10 (YBR019C)| 11.2             | 12.4               | 2,487.9              | 7.74                  |
| GCY1 (YOR120W)| 75.4              | 71.2               | 523.6                | 2.84                  |

Save this exact table in `results/gal_regulon.tsv` and include it in the README.

### PCA scatter

The scatter should look like:

```
PC2
 ^
 |
 |          glucose_rep2
 |       ●
 |     glucose_rep1
 |    ●
 |
 |
 |                                      galactose_rep1
 |                                          ●
 +────────────────────────────────────────────────────► PC1
```

Two glucose replicates cluster tightly on the left of PC1; one galactose replicate sits on the right of PC1. PC2 captures the small replicate-to-replicate variation within glucose.

---

## Write-up template (paste into mini-project/README.md, fill in your numbers)

```markdown
# Week 7 Mini-Project — Yeast RNA-seq

## Dataset

Three *Saccharomyces cerevisiae* RNA-seq samples from Gierlinski et al. 2015 (*Bioinformatics* 31:3625), accessed via NCBI SRA: SRR453566 (glucose, replicate 1), SRR453567 (glucose, replicate 2), SRR453568 (galactose, replicate 1). All three are paired-end 75 bp Illumina HiSeq 2000 runs, polyA-selected libraries, ~3 M pairs each.

## Methods

Reads were trimmed with **fastp 0.23.4** (Chen et al. 2018) using `--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 36 --trim_poly_g`. Trimmed reads were pseudoaligned to the **Ensembl release 110 yeast cDNA** transcriptome (6,975 transcripts) with **kallisto 0.50.1** (Bray et al. 2016) using `kallisto quant -t 4 -b 100`. Per-transcript estimated counts were aggregated to per-gene counts (in yeast, this aggregation is mostly a rename — one transcript per gene for most loci). The three per-sample columns were joined on `gene_id` to produce a 6,975 × 3 counts matrix; TPM was computed per the Wagner et al. 2012 definition; log2(TPM+1) was used as input to a 2-component PCA after filtering to genes with CPM > 1 in at least 2 samples.

## Results

Per-sample QC summary:

| Sample          | n_processed  | n_pseudoaligned | p_pseudoaligned |
|-----------------|-------------:|----------------:|----------------:|
| glucose_rep1    | <YOUR NUMBER> | <YOUR NUMBER>  | <YOUR%>         |
| glucose_rep2    | <YOUR NUMBER> | <YOUR NUMBER>  | <YOUR%>         |
| galactose_rep1  | <YOUR NUMBER> | <YOUR NUMBER>  | <YOUR%>         |

PCA (after CPM filtering, 3,812 genes retained):

- PC1 explains <X>% variance; cleanly separates carbon source.
- PC2 explains <Y>% variance; captures glucose replicate-to-replicate variation.

Top 5 PC1 loadings (positive end, galactose): GAL1, GAL7, GAL10, GAL2, GAL3.

Top 5 PC1 loadings (negative end, glucose): HXT1, HXT3, HXT4, MIG2, RGT1.

GAL regulon (TPM):

| Gene  | glucose_rep1 | glucose_rep2 | galactose_rep1 | log2 FC |
|-------|-------------:|-------------:|---------------:|--------:|
| GAL1  | <X>          | <Y>          | <Z>            | <FC>    |
| ...   |              |              |                |         |

## Discussion

The PCA recovers the textbook galactose-induction signature: PC1 separates the carbon source, and the top loadings are the GAL regulon (positive direction) and the hexose transporters / glucose-repressed regulators (negative direction). The galactose sample shows TPM 2,894 for GAL1 vs ~5 in the glucose samples — a fold change of ~500x, consistent with the well-characterized GAL4/MIG1-mediated transcriptional switch.

Two technical caveats deserve mention. First, with only one galactose replicate, no formal differential-expression p-value can be assigned to any gene; PC1 loadings here are descriptive, not inferential. Week 8 brings in additional replicates and the DESeq2 / edgeR framework for proper statistics. Second, the kallisto pseudoalignment rate is ~95% across all three samples; the unaligned 5% includes rRNA carryover (polyA selection is imperfect), residual adapter sequence (fastp catches most but not all), and reads from non-coding RNAs not in the Ensembl cDNA FASTA.

The pipeline runs end to end in ~15 minutes on a 4-core laptop (excluding SRA download). The matrix and TPM files committed to this directory are the input to Week 8's differential expression analysis.

## Reproducibility

```
conda env create -f env.yml -n c10-week07-miniproject
conda activate c10-week07-miniproject
bash run.sh
```

Output: `results/pca_scatter.png` and `results/counts_matrix.tsv`.
```

---

## Common pitfalls

**The PCA scatter does not separate glucose from galactose.** Three common causes:

1. You forgot to compute log2(TPM+1) and ran PCA on raw TPM. The dynamic range of raw TPM is 10^4 and a few very-highly-expressed genes (RPL28, ACT1) dominate the variance. Always log-transform before PCA.
2. You did not filter low-count genes. PCA on the full 6,975 × 3 matrix is dominated by noise from the 2,000 lowly expressed transcripts. Filter to CPM > 1 in ≥ 2 samples.
3. You mixed up the sample labels. Double-check the SRA accession → condition mapping in your `samples.tsv`.

**`p_pseudoaligned` is much lower for one sample.** Either the sample is contaminated (rRNA carryover, foreign-organism contamination) or you accidentally used a different organism's index. Verify by re-running fastp's overrepresented-sequences output; if rRNA hexamers dominate, your sample has rRNA contamination and your effective library size is smaller than the raw count suggests.

**TPM column sums are not 10^6.** Bug in your TPM computation. The classical mistake is using `length` instead of `eff_length` in the denominator. kallisto's `abundance.tsv` reports `eff_length` directly; use it.

**GAL1 is not in the top 5 PC1 loadings.** Either the sample labels are swapped (your "galactose" is actually a glucose sample), or your PCA is dominated by a different signal (e.g. a batch effect that swamps the biology). Plot the top 10 PC1 loadings and their `log2(TPM+1)` values across samples; inspect by hand which genes are driving the separation.

**Fractional counts crash a downstream R script.** kallisto's `est_counts` are non-integer for multi-mapping reads. Round before passing to DESeq2: `pd.DataFrame(...).round().astype(int)`.

**The pipeline takes 90 minutes instead of 15.** Most likely your SRA download is on a slow connection or the kallisto bootstrap (-b 100) is being re-run unnecessarily. Drop `-b 100` for the mini-project if you do not need bootstrap variance for Week 8.

---

## Stretch goals

If you finish early and want to push further, try any of the following:

- **Re-run with Salmon** (`salmon quant --validateMappings --gcBias --seqBias`) and produce a second PCA from the Salmon counts. Verify the two PCAs are essentially identical (PC1 loadings ≥ 95% overlap, PC1 percent-variance within 2 percentage points).
- **Re-run with HISAT2 + featureCounts** and compute the Pearson r on log2(TPM+1) between the kallisto and the alignment-based counts. Expect r ≥ 0.98 on the well-expressed subset. Tabulate the top 10 most discordant genes (almost always rRNA).
- **Add a fourth galactose sample.** SRR453569 is the galactose-rep-2 from the same Gierlinski et al. 2015 dataset. Re-run the pipeline with `samples.tsv` updated to four samples; verify the two galactose samples now cluster on the right of PC1 and the two glucose samples on the left.
- **Run a quick `DESeq2` analysis** (an R one-liner via `rpy2` or a separate R script) on your 3-sample matrix. Report the top 10 up- and down-regulated genes by adjusted p-value. The galactose-rep-1 vs glucose-rep-mean DE will have wide CI because of the n=1 vs n=2 design — note that in the writeup; this is a Week 8 preview.
- **Annotate the top 20 PC1 loadings with GO terms.** The yeast SGD has GO annotations for every gene; the GAL regulon should be enriched for "galactose metabolic process" (GO:0006012). This is a tiny dive into the GO enrichment that Week 9 will cover formally.

---

## Submission

Push your `mini-project/` directory to the public portfolio repo with a commit message like:

```
Week 7 mini-project: yeast RNA-seq, 3 samples, kallisto+PCA, GAL regulon recovered.
```

Open a PR to the curriculum repo with a brief description of what you produced and any deviations from the recipe. The graders will:

1. Clone the repo on a clean machine, run `bash run.sh`, and verify the matrix and PCA reproduce.
2. Read the README. They expect every number to be cited (no "high" or "good" without a number); every tool to be pinned (fastp 0.23.4, not fastp).
3. Inspect `results/gal_regulon.tsv` and verify GAL1 has TPM > 1,000 in the galactose sample and TPM < 50 in the glucose samples.

---

## Up next

Week 8 takes your counts matrix as input and produces a formal differential-expression analysis with adjusted p-values, plus a pathway-enrichment table. The matrix you committed here is the literal input file for Week 8. Treat it accordingly.
