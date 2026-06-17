# Lecture 3 — Project Tracks, run-info.json, and the Deposit

> **Educational and research use only.** The four capstone project tracks below use real public data. Outputs are educational artefacts; no variant call, expression difference, phylogeny, or assembly produced in the capstone is clinical software. Where the upstream data is human-derived (Track 1: GIAB; Track 2: GEO RNA-seq), the original donors consented to public deposit and the data is fully de-identified per the deposit's terms of use; respect those terms.

Lectures 1 and 2 set up the workflow-manager and container layers. Lecture 3 covers what you put through them: the four project tracks, the canonical provenance record (`run-info.json`), the Zenodo / DOI deposit step, and the wrap-up sidebar that closes C10.

## 3.1 — Picking a track

The capstone has four pre-designed tracks. Pick **one**. The tracks are sized to fit a 16 GB laptop and a 5-day work week.

| Track | Question | Data source | Pipeline | Output | C-week scaffold |
|-------|----------|-------------|----------|--------|-----------------|
| 1: Variant discovery | "How accurate is my variant caller on a known truth-set sample?" | GIAB HG002 chr20 (NIST) | align → call → benchmark | Precision / recall / F1 vs GIAB truth | W5, W6, W8 |
| 2: RNA-seq DE | "Which genes are differentially expressed in condition A vs B?" | GEO series (e.g. GSE52778) | trim → quantify → aggregate → test | Volcano plot, top-30 gene table | W7 |
| 3: MSA + phylogeny | "What does the phylogeny of N viral genomes look like?" | NCBI Virus (e.g. SARS-CoV-2 lineage) | align → trim → infer tree → visualize | Tree image, bootstrap-support table | W9 |
| 4: Long-read assembly | "Can I assemble a small bacterial genome from Nanopore reads?" | SRA (e.g. *E. coli* K-12 R10.4) | assemble → polish → QC | Polished FASTA, BUSCO / Quast report | W10 |

The track choice is one of taste and prior experience. Track 1 is the most "wet-lab adjacent" — the GIAB truth set lets you measure your pipeline's accuracy quantitatively. Track 2 produces the most-publishable artefact for a portfolio (volcano plots are visually striking). Track 3 has the shortest runtime (5-20 minutes) and is the gentlest on hardware. Track 4 is the most algorithmically interesting and the heaviest on RAM (Flye on a 4-Mb genome wants 8-12 GB).

You do not need to extend the analysis beyond the canonical pipeline. A clean Track 2 with eight samples, a `Snakefile` of 6-10 rules, a pinned environment, and a Zenodo-deposited tagged release is a passing capstone. A bespoke 15-rule pipeline with 100 samples but no version pinning is not.

## 3.2 — Track 1 deep dive: variant discovery

**Question:** I align Illumina short reads from HG002 (a Genome in a Bottle sample) to GRCh38 chr20, call variants with DeepVariant or GATK HaplotypeCaller, and compare against the GIAB chr20 truth VCF. What are the precision, recall, and F1?

**Input data:**

- HG002 30x WGS reads, chr20 only. Free at <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/>. The full WGS is 60-100 GB; the chr20 subset is ~2 GB.
- GRCh38 chr20 reference FASTA. Free from the GATK resource bundle at <https://console.cloud.google.com/storage/browser/genomics-public-data/references/GRCh38_Verily/>.
- GIAB HG002 truth VCF v4.2.1. Free at <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/>.

**Pipeline (Snakemake rule sketch):**

```
rule all:
    input:
        "results/happy_summary.tsv",
        "figures/precision_recall.png"

rule trim_fastq:
    input:  "data/HG002_chr20_R1.fastq.gz", "data/HG002_chr20_R2.fastq.gz"
    output: "qc/HG002_chr20_R1.trim.fastq.gz", "qc/HG002_chr20_R2.trim.fastq.gz"
    shell:  "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} ..."

rule align:
    input:  "qc/HG002_chr20_R1.trim.fastq.gz", "qc/HG002_chr20_R2.trim.fastq.gz", "ref/GRCh38.chr20.fa"
    output: "bam/HG002_chr20.bam", "bam/HG002_chr20.bam.bai"
    shell:  "bwa mem -t {threads} {input[2]} {input[0]} {input[1]} | samtools sort -o {output[0]} && samtools index {output[0]}"

rule call_variants:
    input:  "bam/HG002_chr20.bam", "ref/GRCh38.chr20.fa"
    output: "vcf/HG002_chr20.vcf.gz", "vcf/HG002_chr20.vcf.gz.tbi"
    shell:  "deepvariant --model_type=WGS --ref={input[1]} --reads={input[0]} --output_vcf={output[0]}"

rule benchmark:
    input:  "vcf/HG002_chr20.vcf.gz", "truth/HG002_chr20_truth.vcf.gz", "truth/HG002_chr20_confident.bed"
    output: "results/happy_summary.tsv"
    shell:  "hap.py {input[1]} {input[0]} -f {input[2]} -r ref/GRCh38.chr20.fa -o results/happy"

rule plot:
    input:  "results/happy_summary.tsv"
    output: "figures/precision_recall.png"
    shell:  "python scripts/plot_precision_recall.py {input} {output}"
```

**Expected output:** A `happy_summary.tsv` with one row per variant type (SNVs vs indels) and columns for precision, recall, F1. Typical numbers on chr20 HG002 with DeepVariant: precision ~0.998, recall ~0.997, F1 ~0.998 on SNVs; precision ~0.990, recall ~0.990, F1 ~0.990 on indels. Numbers vary by caller; the capstone reports its own.

**Runtime:** 30-90 minutes on a 16 GB laptop with 8 cores.

## 3.3 — Track 2 deep dive: RNA-seq differential expression

**Question:** Given 8 RNA-seq samples from GSE52778 (4 control + 4 dexamethasone-treated airway smooth muscle cells), which genes are differentially expressed?

**Input data:**

- GSE52778 metadata at <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778>. 8 samples (4 control, 4 dexamethasone). FASTQ files available via SRA at SRP033351.
- GENCODE v44 transcriptome (chr1+chr19+chr21 subset for the teaching version; full transcriptome for the production version). Free at <https://www.gencodegenes.org/human/release_44.html>.

**Pipeline (Snakemake rule sketch):**

```
SAMPLES = ["SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
           "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521"]

rule all:
    input:
        "results/dge_top30.tsv",
        "figures/volcano.png",
        "figures/pca.png"

rule download_fastq:
    output: "data/{sample}_R1.fastq.gz", "data/{sample}_R2.fastq.gz"
    shell:  "fasterq-dump --split-files {wildcards.sample} -O data/ && gzip data/{wildcards.sample}_1.fastq && gzip data/{wildcards.sample}_2.fastq && mv data/{wildcards.sample}_1.fastq.gz {output[0]} && mv data/{wildcards.sample}_2.fastq.gz {output[1]}"

rule trim_fastq:
    input:  "data/{sample}_R1.fastq.gz", "data/{sample}_R2.fastq.gz"
    output: "qc/{sample}_R1.trim.fastq.gz", "qc/{sample}_R2.trim.fastq.gz"
    shell:  "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} ..."

rule salmon_index:
    input:  "ref/gencode.v44.transcripts.fa.gz"
    output: directory("ref/salmon_index")
    shell:  "salmon index -t {input} -i {output} -k 31"

rule salmon_quant:
    input:  "qc/{sample}_R1.trim.fastq.gz", "qc/{sample}_R2.trim.fastq.gz", "ref/salmon_index"
    output: directory("quants/{sample}")
    shell:  "salmon quant -i {input[2]} -l A -1 {input[0]} -2 {input[1]} -o {output} -p {threads}"

rule deseq2:
    input:  expand("quants/{sample}", sample=SAMPLES), "metadata.tsv", "ref/tx2gene.tsv"
    output: "results/deseq2_results.rds", "results/dge_top30.tsv", "figures/volcano.png", "figures/pca.png"
    shell:  "Rscript scripts/run_deseq2.R"
```

**Expected output:** A `dge_top30.tsv` ranked by p-adjusted, with columns `gene_id`, `gene_name`, `log2FC`, `lfcSE`, `pvalue`, `padj`, `baseMean`. A `volcano.png` plotting `-log10(padj)` against `log2FC`. A `pca.png` showing the eight samples colored by condition.

**Expected genes for GSE52778:** the canonical hits are DUSP1, KLF15, PER1, FKBP5 (all glucocorticoid-responsive). The capstone should recover at least three of the four in the top 30.

**Runtime:** 1-3 hours on a 16 GB laptop with 8 cores. The Salmon index step is the slow one (5-10 minutes); the actual quant is 1-5 minutes per sample.

## 3.4 — Track 3 deep dive: MSA + phylogeny

**Question:** I download 50 SARS-CoV-2 genomes (lineage B.1.617.2 / Delta, one country, one month), align them with MAFFT, trim with trimAl, infer a tree with IQ-TREE 2 under model `GTR+G+I` with 1000 ultrafast bootstraps, and visualize the tree annotated by sampling date.

**Input data:**

- 50 SARS-CoV-2 genomes from NCBI Virus. Free at <https://www.ncbi.nlm.nih.gov/labs/virus/>. Filter by lineage = `B.1.617.2`, country, collection date range, and download FASTA + CSV metadata.
- Reference genome: Wuhan-Hu-1 (NC_045512.2). Free.

**Pipeline (Snakemake rule sketch):**

```
rule all:
    input:
        "figures/tree.png",
        "results/bootstrap_table.tsv"

rule combine_sequences:
    input:  "data/lineage_B.1.617.2_50.fasta", "data/wuhan_hu_1.fasta"
    output: "work/combined.fasta"
    shell:  "cat {input} > {output}"

rule align:
    input:  "work/combined.fasta"
    output: "work/aligned.fasta"
    shell:  "mafft --auto --thread {threads} {input} > {output}"

rule trim:
    input:  "work/aligned.fasta"
    output: "work/trimmed.fasta"
    shell:  "trimal -in {input} -out {output} -gappyout"

rule infer_tree:
    input:  "work/trimmed.fasta"
    output: "work/iqtree.treefile", "work/iqtree.iqtree"
    shell:  "iqtree2 -s {input} -m GTR+G+I -B 1000 -T {threads} -seed 42 --prefix work/iqtree"

rule visualize:
    input:  "work/iqtree.treefile", "data/metadata.tsv"
    output: "figures/tree.png", "results/bootstrap_table.tsv"
    shell:  "python scripts/plot_tree.py {input[0]} {input[1]} {output[0]} {output[1]}"
```

**Expected output:** A `tree.png` showing the 50-taxon tree with tip labels color-coded by sampling date. A `bootstrap_table.tsv` listing every internal node and its bootstrap support; a properly-supported phylogeny has >75% of internal nodes with >70% bootstrap support.

**Runtime:** 5-20 minutes for the whole pipeline on 8 cores. MAFFT on 50 genomes takes 30-60 seconds; IQ-TREE 2 with 1000 ultrafast bootstraps on 50 taxa of ~30 kb each takes 2-10 minutes.

## 3.5 — Track 4 deep dive: long-read assembly

**Question:** I download Nanopore reads for *E. coli* K-12 MG1655 at 60x coverage, assemble with Flye, polish with Medaka, and quality-check with BUSCO and Quast. How close is my assembly to the reference?

**Input data:**

- Nanopore R10.4 reads for *E. coli* K-12 MG1655 from SRA (e.g. SRR21070078 or any comparable accession). Free. Trim to 60x for the teaching subset (~3 GB).
- Reference genome U00096.3 (E. coli K-12 MG1655). Free at <https://www.ncbi.nlm.nih.gov/nuccore/U00096.3>.

**Pipeline (Snakemake rule sketch):**

```
rule all:
    input:
        "results/busco_summary.txt",
        "results/quast_report.tsv",
        "results/polished.fasta"

rule filter_reads:
    input:  "data/ecoli_nanopore.fastq.gz"
    output: "qc/filtered.fastq.gz"
    shell:  "filtlong --min_length 1000 --min_mean_q 7 {input} | gzip > {output}"

rule assemble:
    input:  "qc/filtered.fastq.gz"
    output: directory("work/flye_out"), "work/flye_out/assembly.fasta"
    shell:  "flye --nano-raw {input} --out-dir work/flye_out --genome-size 5m --threads {threads} --seed 42"

rule polish:
    input:  "work/flye_out/assembly.fasta", "qc/filtered.fastq.gz"
    output: "results/polished.fasta"
    shell:  "medaka_consensus -i {input[1]} -d {input[0]} -o work/medaka_out -t {threads} -m r1041_e82_400bps_sup_v5.0.0 && cp work/medaka_out/consensus.fasta {output}"

rule busco:
    input:  "results/polished.fasta"
    output: "results/busco_summary.txt"
    shell:  "busco -i {input} -o work/busco_out --lineage_dataset bacteria_odb10 --mode genome --cpu {threads} && cp work/busco_out/short_summary*.txt {output}"

rule quast:
    input:  "results/polished.fasta", "ref/ecoli_reference.fasta"
    output: "results/quast_report.tsv"
    shell:  "quast.py {input[0]} -r {input[1]} -o work/quast_out -t {threads} && cp work/quast_out/report.tsv {output}"
```

**Expected output:** A polished `assembly.fasta` of one to three contigs spanning ~4.6 Mb. A BUSCO report showing ~99% complete-single-copy genes on `bacteria_odb10`. A Quast report comparing assembly size, NG50, and number of mismatches against the reference. Typical numbers for E. coli at 60x with R10.4: one circular contig of ~4.64 Mb, BUSCO C:99.2% [S:99.0%, D:0.2%], Quast NG50 ~4.6 Mb, mismatches per 100 kb ~0.5.

**Runtime:** 45-90 minutes on a 16 GB laptop with 8 cores. Flye on 60x is the slow rule (20-40 minutes); Medaka polish is 10-20 minutes; BUSCO and Quast are 5-15 minutes.

## 3.6 — run-info.json: the inline provenance record

Every capstone run writes a `run-info.json` next to the output. The pattern from Weeks 8-11, extended for Week 12's workflow-manager context.

A minimal `run-info.json` for a Snakemake run:

```json
{
  "project": "C10-Capstone-RNA-seq-GSE52778",
  "version": "v1.0",
  "run_date_utc": "2026-05-14T19:42:03Z",
  "git_commit": "def4567890abcdef1234567890abcdef12345678",
  "git_tag": "v1.0",
  "workflow_manager": "snakemake",
  "workflow_manager_version": "7.32.4",
  "workflow_file": "Snakefile",
  "workflow_file_sha256": "abc...",
  "conda_environment": "environment.yml",
  "conda_lockfile": "environment.lock.linux-64.txt",
  "conda_lockfile_sha256": "ghi...",
  "container_image": "c10-capstone-rnaseq-v1.0.sif",
  "container_image_sha256": "jkl...",
  "container_build_date_utc": "2026-05-13T11:20:14Z",
  "host": {
    "hostname": "lab-laptop-04",
    "os": "Ubuntu 22.04.4 LTS",
    "kernel": "5.15.0-105-generic",
    "cpu_model": "Intel(R) Core(TM) i7-1165G7",
    "cpu_count_logical": 8,
    "ram_gb": 16
  },
  "input_data": {
    "study_accession": "GSE52778",
    "sra_run_accessions": ["SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
                           "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521"],
    "metadata_sha256": "mno...",
    "transcriptome_release": "GENCODE v44",
    "transcriptome_sha256": "pqr..."
  },
  "pinned_parameters": {
    "fastp_qualified_quality_phred": 20,
    "fastp_length_required": 25,
    "salmon_index_kmer_size": 31,
    "salmon_quant_libtype": "A",
    "deseq2_alpha": 0.05,
    "deseq2_lfcThreshold": 0.0,
    "rng_seed": 42
  },
  "outputs": {
    "results/dge_top30.tsv": {"sha256": "stu...", "rows": 30, "cols": 7},
    "figures/volcano.png": {"sha256": "vwx..."},
    "figures/pca.png": {"sha256": "yza..."}
  },
  "runtime_seconds": 4523,
  "license": "CC-BY-4.0 for data, MIT for code",
  "deposit_doi": "10.5281/zenodo.99999999",
  "disclaimer": "Educational and research use only. Not validated for clinical use."
}
```

Read it section by section:

- **project / version / run_date_utc / git_commit / git_tag** — what was run.
- **workflow_manager, workflow_manager_version, workflow_file_sha256** — the workflow layer.
- **conda_environment, conda_lockfile_sha256** — the package-manager layer.
- **container_image, container_image_sha256, container_build_date_utc** — the container layer.
- **host** — where the run happened.
- **input_data** — the source data, accessions, and hashes.
- **pinned_parameters** — every non-default tool parameter.
- **outputs** — every output file with its hash.
- **runtime_seconds** — wall-clock.
- **license / deposit_doi / disclaimer** — the deposit metadata.

Write this from inside the workflow. Snakemake's `onsuccess:` handler is a good place:

```python
onsuccess:
    import json, hashlib, subprocess, datetime, platform, os, time
    info = {...}  # populate as above
    with open("results/run-info.json", "w") as fh:
        json.dump(info, fh, indent=2)
```

Nextflow's `workflow.onComplete` handler is the equivalent.

## 3.7 — The Zenodo deposit

Zenodo (<https://zenodo.org/>) is the CERN-hosted free repository for academic data deposits. It issues a DOI for every uploaded artefact. Files up to 50 GB. The GitHub-Zenodo integration archives every tagged release automatically.

**Steps to deposit the capstone:**

1. Tag a release in your GitHub repo: `git tag -a v1.0 -m "C10 capstone v1.0" && git push origin v1.0`.
2. Sign in to Zenodo with GitHub OAuth (free).
3. Visit <https://zenodo.org/account/settings/github/> and flip the toggle for your capstone repo to "On."
4. Push another tag (or re-tag). Zenodo automatically archives the release and issues a DOI within a few minutes.
5. Copy the DOI badge markdown from Zenodo's repo page and paste it into your `README.md`.
6. Commit the badge.

The result: every tagged release of your capstone is archived under a DOI. The DOI is permanent (Zenodo guarantees long-term hosting). The artefact is citeable.

**What to include in the deposit:**

- The source tree (the workflow files, the scripts, the `environment.yml`, the `Singularity.def`).
- The `run-info.json` (the canonical run record).
- The output figures and tables (small; <10 MB total).
- The `report.md` and the `wrap-up.md`.

**What to exclude:**

- Large input data files (link instead via `download.sh` with checksums).
- Container images (link to Docker Hub / Zenodo as separate deposits if needed; `.sif` files can be uploaded separately).
- Conda environment caches.
- The workflow manager's working directory (`work/`, `.snakemake/`).

## 3.8 — CITATION.cff

GitHub honours a `CITATION.cff` file in the repo root and surfaces a "Cite this repository" button. Add one:

```yaml
cff-version: 1.2.0
title: "C10 Capstone: GSE52778 RNA-seq differential expression"
message: "If you use this software, please cite it as below."
type: software
authors:
  - family-names: "YourLastName"
    given-names: "YourFirstName"
    orcid: "https://orcid.org/0000-0000-0000-0000"
date-released: "2026-05-14"
version: "1.0"
doi: "10.5281/zenodo.99999999"
url: "https://github.com/yourname/c10-capstone"
license: "MIT"
keywords:
  - bioinformatics
  - rna-seq
  - snakemake
  - differential-expression
references:
  - type: article
    authors:
      - family-names: "Mölder"
        given-names: "Felix"
    title: "Sustainable data analysis with Snakemake"
    journal: "F1000Research"
    year: 2021
    volume: 10
    start: 33
```

This is a short file but it is the canonical metadata for software-citation tooling.

## 3.9 — The report.md

A one-page write-up of the capstone. Sections:

- **Question.** One sentence: "I investigated *which genes are differentially expressed in airway smooth muscle cells treated with dexamethasone* using public data from GSE52778."
- **Data.** Accession, sample count, conditions, total bytes.
- **Pipeline.** One paragraph: "FASTQ files were trimmed with fastp 0.23.4, quantified with Salmon 1.10.2 against GENCODE v44 transcripts, aggregated to gene level with tximport 1.30.0, and tested for differential expression with DESeq2 1.42.0. The pipeline is implemented in Snakemake 7.32.4 with the Conda environment captured to `environment.lock.linux-64.txt` (47 packages) and built into an Apptainer 1.2.5 container `c10-capstone-rnaseq-v1.0.sif`."
- **Results.** Two paragraphs and one figure. Numbers, gene names, what was recovered.
- **Limitations.** Which Weeks 1-12 caveats apply. The capstone is not clinical; the n is small; the technical replicates are not biological replicates; the gene-level aggregation discards transcript isoform information; etc.
- **Reproducibility.** One paragraph: "The pipeline runs end-to-end on a Linux x86_64 16 GB laptop in 1 hour 14 minutes with `apptainer run c10-capstone-rnaseq-v1.0.sif snakemake --cores 8 --use-conda`. The byte-identical reproduction check passes on a clean Ubuntu 22.04 VM. The run is deposited at DOI 10.5281/zenodo.99999999."
- **Citations.** The papers behind every tool: Snakemake (Mölder et al. 2021), bioconda (Grüning et al. 2018), Singularity (Kurtzer et al. 2017), fastp (Chen et al. 2018), Salmon (Patro et al. 2017), tximport (Soneson et al. 2015), DESeq2 (Love et al. 2014), GENCODE (Frankish et al. 2023).

One page. No more. The deposit's `report.md` is the executive summary; the `Snakefile` and the code are the technical depth.

## 3.10 — The wrap-up sidebar

`wrap-up.md` lives next to `report.md` and is the C10 retrospective. One paragraph per Week 1-12. The template:

```markdown
# C10 Wrap-Up — Crunching Bioinformatics

A retrospective on the 12 weeks of C10. Author: <your name>. Date: <YYYY-MM-DD>.

## Week 1 — Vocabulary and Ethics
What I learned. What stuck. What I would do differently.

## Week 2 — FASTA / FASTQ I/O
...

## Week 12 — Capstone
...

## Overall
One paragraph summarizing the arc of C10 from your perspective.
```

The wrap-up is not graded for length; it is graded for honesty. Name the weeks you struggled with. Name the weeks that clicked. Name the tools you would pick differently if you redid the curriculum. The wrap-up is the only piece of the capstone where you are encouraged to be opinionated.

## 3.11 — The reproducibility check

The capstone's pass criterion is reproducibility. The check:

1. Push the final commit to GitHub.
2. Spin up a clean Linux VM (a freshly-launched Codespace, a clean Vagrant box, a 50-cent EC2 instance — anything that has never seen your environment).
3. Install Apptainer (one `apt install apptainer` on Ubuntu 22.04).
4. Clone the repo. `cd` into it.
5. Build the container: `apptainer build c10-capstone.sif Singularity.def`.
6. Run the pipeline: `apptainer run c10-capstone.sif snakemake --cores 4 --use-conda`.
7. Compare the output hashes: `sha256sum -c results/hashes.txt`.

If step 7 passes, the capstone is byte-identically reproducible. If it does not, you have a non-determinism somewhere (un-pinned seed, un-pinned thread count, time-of-day in a header, filesystem-order dependency). The Week 12 challenges and homework walk through diagnosing these.

A passing capstone is byte-identical on the deterministic rules and within-tolerance on the documented stochastic rules. A failing capstone is byte-identical-by-accident only (passes on one machine, fails on another).

## 3.12 — Where the capstone goes next

After the deposit:

- **Add a CI workflow** (GitHub Actions): re-run the pipeline on a tiny test subset on every push. Pass / fail badge in the README.
- **Write a blog post.** The capstone is portfolio material. A 1500-word blog post on your personal site (free hosting on GitHub Pages) describing the question, the pipeline, and the deposit is a strong portfolio item.
- **Submit a nf-core PR.** If your capstone is Nextflow, the nf-core community accepts community contributions for new pipelines (<https://nf-co.re/contributing/>).
- **Cite it in your thesis / honours project.** The capstone has a DOI; it is a citeable artefact. Use it.

## 3.13 — Recap

- The four capstone tracks are sized to fit a 16 GB laptop and a 5-day work week. Pick one.
- Every track has a worked Snakefile sketch (Section 3.2-3.5). Start from the sketch, scale to your dataset, polish.
- The `run-info.json` (Section 3.6) is the inline provenance record. Write it at the end of every run.
- Zenodo (Section 3.7) is the free DOI deposit. Tag a GitHub release; Zenodo archives automatically.
- `CITATION.cff` (Section 3.8) makes the deposit citeable from GitHub's UI.
- `report.md` (Section 3.9) is the one-page executive summary.
- `wrap-up.md` (Section 3.10) is the C10 retrospective. One paragraph per week.
- The byte-identical reproduction check (Section 3.11) is the capstone's pass criterion.

You now have everything you need to build the capstone. The exercises walk you through the workflow-manager basics on a toy example; the challenges push you on DAG rendering and reproduction; the mini-project is the capstone itself.
