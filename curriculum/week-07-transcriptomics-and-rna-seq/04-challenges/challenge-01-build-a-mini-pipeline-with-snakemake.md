# Challenge 1 — Build a mini-pipeline with Snakemake

> **Estimated time:** 2 hours.
> **Goal:** Wire `fastp` and `kallisto` into a 3-sample Snakemake pipeline that takes a directory of paired-end FASTQs and produces a gene-level counts matrix on demand. Demonstrate that the pipeline is *reproducible* (a fresh checkout + `snakemake --cores 4` recovers the same matrix), *idempotent* (re-running does no extra work), and *parallelizable* (3 samples run concurrently on a multi-core laptop).

This challenge is the bridge between "I ran fastp and kallisto by hand on one sample" (Exercises 1 and 2) and "I have a 10-sample RNA-seq study that re-runs from a single command when the GTF version updates." Snakemake (Köster and Rahmann 2012, *Bioinformatics* 28:2520) is the de-facto standard workflow manager for bioinformatics pipelines and is what `nf-core/rnaseq` and most large-scale RNA-seq efforts use under the hood.

---

## Background — Snakemake in 5 minutes

Snakemake is a Python-based workflow manager where each step of a pipeline is declared as a **rule**:

```python
rule trim:
    input:
        r1 = "raw/{sample}_1.fastq.gz",
        r2 = "raw/{sample}_2.fastq.gz",
    output:
        r1 = "trim/{sample}_1.trim.fq.gz",
        r2 = "trim/{sample}_2.trim.fq.gz",
        html = "qc/{sample}.fastp.html",
        json = "qc/{sample}.fastp.json",
    threads: 4
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --trim_poly_g \
            -h {output.html} -j {output.json} \
            -w {threads}
        """
```

Snakemake reads the rules, builds a dependency DAG (directed acyclic graph) between outputs and inputs, and schedules rules in topological order on the available cores. Rules that are independent run in parallel. Rules whose outputs already exist (and are newer than the inputs) are skipped — this is the **idempotency** property.

Reference: Köster J, Rahmann S 2012, *Bioinformatics* 28:2520, "Snakemake — a scalable bioinformatics workflow engine," <https://academic.oup.com/bioinformatics/article/28/19/2520/290322>.

---

## Task

Build a Snakefile that orchestrates a 3-sample yeast RNA-seq pipeline.

### Layout

```
crunch-bio-portfolio-<yourhandle>/
└── week-07/
    └── challenge-01/
        ├── Snakefile
        ├── config.yaml         (sample IDs, reference URL, threads)
        ├── env.yml             (conda env file)
        ├── README.md           (how-to-run + DAG image)
        ├── raw/                (downloaded FASTQs; gitignored)
        ├── ref/                (Ensembl transcriptome FASTA; gitignored)
        ├── trim/               (fastp output; gitignored)
        ├── qc/                 (fastp HTML/JSON; committed)
        ├── index/              (kallisto index; gitignored)
        ├── quant/              (kallisto quant output; abundance.tsv committed)
        ├── counts/             (final 3-sample matrix; committed)
        └── logs/               (per-rule logs)
```

### Sample IDs

Use the same three yeast samples as the mini-project:

| Sample ID | SRA accession | Condition           |
|-----------|---------------|---------------------|
| `S1`      | `SRR453566`   | glucose, replicate 1 |
| `S2`      | `SRR453567`   | glucose, replicate 2 |
| `S3`      | `SRR453568`   | galactose, replicate 1 |

### Required rules

The Snakefile must contain rules for:

1. **`download`** — fetch a sample's raw FASTQs from SRA (or a no-cost mirror). Inputs: SRA accession (from `config.yaml`). Outputs: `raw/{sample}_1.fastq.gz`, `raw/{sample}_2.fastq.gz`.
2. **`fetch_transcriptome`** — download the yeast Ensembl cDNA FASTA. Inputs: URL (config). Output: `ref/sce.cdna.fa.gz`.
3. **`build_index`** — `kallisto index`. Input: `ref/sce.cdna.fa.gz`. Output: `index/sce.idx`.
4. **`trim`** — `fastp` on a sample's raw FASTQs. Inputs: `raw/{sample}_{1,2}.fastq.gz`. Outputs: `trim/{sample}_{1,2}.trim.fq.gz`, `qc/{sample}.fastp.{html,json}`.
5. **`quant`** — `kallisto quant` on a sample's trimmed FASTQs. Inputs: `index/sce.idx`, `trim/{sample}_{1,2}.trim.fq.gz`. Outputs: `quant/{sample}/abundance.tsv`, `quant/{sample}/abundance.h5`, `quant/{sample}/run_info.json`.
6. **`build_matrix`** — Python script that loads each sample's `abundance.tsv`, aggregates per-transcript counts to per-gene counts using a transcript-to-gene map, and writes `counts/all_samples.counts.tsv`. (For yeast, transcripts and genes are nearly identical, so the aggregation is mostly a renaming.) Inputs: all three `quant/{sample}/abundance.tsv`. Output: `counts/all_samples.counts.tsv`.
7. **`all`** — a target rule with `input: counts/all_samples.counts.tsv` that triggers the entire pipeline.

### The Snakefile (reference outline)

```python
configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())   # ["S1", "S2", "S3"]


rule all:
    input:
        "counts/all_samples.counts.tsv",
        expand("qc/{sample}.fastp.html", sample=SAMPLES),


rule download:
    output:
        r1 = "raw/{sample}_1.fastq.gz",
        r2 = "raw/{sample}_2.fastq.gz",
    params:
        sra = lambda wc: config["samples"][wc.sample]["sra"],
    log:
        "logs/download.{sample}.log",
    shell:
        """
        prefetch {params.sra} -O raw/{wildcards.sample}_sra/ > {log} 2>&1
        fasterq-dump raw/{wildcards.sample}_sra/{params.sra}/{params.sra}.sra \
            --split-files -O raw/ >> {log} 2>&1
        gzip raw/{params.sra}_1.fastq && mv raw/{params.sra}_1.fastq.gz {output.r1}
        gzip raw/{params.sra}_2.fastq && mv raw/{params.sra}_2.fastq.gz {output.r2}
        rm -rf raw/{wildcards.sample}_sra/
        """


rule fetch_transcriptome:
    output:
        "ref/sce.cdna.fa.gz",
    params:
        url = config["transcriptome_url"],
    shell:
        "curl -sLo {output} {params.url}"


rule build_index:
    input:
        "ref/sce.cdna.fa.gz",
    output:
        "index/sce.idx",
    log:
        "logs/build_index.log",
    shell:
        "kallisto index -i {output} -k 31 {input} > {log} 2>&1"


rule trim:
    input:
        r1 = "raw/{sample}_1.fastq.gz",
        r2 = "raw/{sample}_2.fastq.gz",
    output:
        r1 = "trim/{sample}_1.trim.fq.gz",
        r2 = "trim/{sample}_2.trim.fq.gz",
        html = "qc/{sample}.fastp.html",
        json = "qc/{sample}.fastp.json",
    threads: 4
    log:
        "logs/trim.{sample}.log",
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 36 \
            --trim_poly_g \
            -h {output.html} -j {output.json} \
            -w {threads} \
            > {log} 2>&1
        """


rule quant:
    input:
        idx = "index/sce.idx",
        r1 = "trim/{sample}_1.trim.fq.gz",
        r2 = "trim/{sample}_2.trim.fq.gz",
    output:
        abundance = "quant/{sample}/abundance.tsv",
        h5 = "quant/{sample}/abundance.h5",
        info = "quant/{sample}/run_info.json",
    threads: 4
    log:
        "logs/quant.{sample}.log",
    shell:
        """
        kallisto quant \
            -i {input.idx} \
            -o quant/{wildcards.sample}/ \
            -t {threads} \
            -b 100 \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule build_matrix:
    input:
        expand("quant/{sample}/abundance.tsv", sample=SAMPLES),
    output:
        "counts/all_samples.counts.tsv",
    script:
        "scripts/build_matrix.py"
```

### The config.yaml

```yaml
transcriptome_url: "http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"

samples:
  S1:
    sra: SRR453566
    condition: glucose
    replicate: 1
  S2:
    sra: SRR453567
    condition: glucose
    replicate: 2
  S3:
    sra: SRR453568
    condition: galactose
    replicate: 1
```

### The scripts/build_matrix.py

```python
"""Aggregate per-sample kallisto abundance.tsv into a 3-sample counts matrix."""

from __future__ import annotations
import pandas as pd
from pathlib import Path


def build_matrix(abundance_paths: list[str], output_path: str) -> None:
    columns: dict[str, pd.Series] = {}
    for path in abundance_paths:
        sample: str = Path(path).parent.name
        df = pd.read_csv(path, sep="\t")
        columns[sample] = df.set_index("target_id")["est_counts"]
    matrix = pd.DataFrame(columns).fillna(0.0)
    matrix = matrix.round().astype(int)
    matrix.to_csv(output_path, sep="\t", index_label="gene_id")


if __name__ == "__main__":
    # Snakemake injects `snakemake.input` and `snakemake.output`.
    build_matrix(list(snakemake.input), str(snakemake.output))
```

### How to run

```bash
# 1. Create the conda env.
conda env create -f env.yml -n c10-week07
conda activate c10-week07

# 2. Sanity check.
snakemake --cores 4 --dry-run

# 3. Run.
snakemake --cores 4

# 4. Visualize the DAG.
snakemake --dag | dot -Tpng > dag.png
```

### Verification

After a successful run:

```bash
ls counts/all_samples.counts.tsv
wc -l counts/all_samples.counts.tsv
# Expected: ~6,975 rows + 1 header = 6,976 lines.

head -5 counts/all_samples.counts.tsv
# Expected:
# gene_id    S1     S2     S3
# YAL001C    243    251    198
# YAL002W    189    195    156
# YAL003W    4521   4612   3987
# YAL005C    1212   1198   1087
```

Run idempotency check:

```bash
snakemake --cores 4
# Expected: "Nothing to be done (all requested files are present and up to date)."
```

Touch a single input and re-run:

```bash
touch raw/S1_1.fastq.gz
snakemake --cores 4 --dry-run
# Expected: shows that `trim S1`, `quant S1`, `build_matrix` would be re-run.
# Other samples (S2, S3) are not affected.
```

This is the value of Snakemake: only the work that needs to be redone is redone.

---

## Acceptance criteria

- [ ] `Snakefile`, `config.yaml`, `env.yml`, and `scripts/build_matrix.py` exist and pass `snakemake --dry-run`.
- [ ] `snakemake --cores 4` builds `counts/all_samples.counts.tsv` from a fresh checkout (or from a pre-populated `raw/` directory if SRA download is impractical).
- [ ] `counts/all_samples.counts.tsv` has ~6,975 rows and 4 columns (`gene_id`, `S1`, `S2`, `S3`).
- [ ] Re-running `snakemake --cores 4` does no extra work (idempotent).
- [ ] A DAG image (`dag.png`) is committed.
- [ ] All `qc/*.fastp.html` reports are committed.
- [ ] `README.md` documents the pipeline, the config, and the expected runtime.
- [ ] Commit message like `c1: 3-sample kallisto pipeline with Snakemake`.

---

## Stretch goals

- **Add a Salmon quant rule** alongside the kallisto rule and produce two parallel counts matrices. Compare them in `scripts/compare_quantifiers.py`.
- **Add a STAR or HISAT2 alignment rule** and a featureCounts rule to produce a third matrix from alignment-based counting. Compare to the pseudoalignment matrices.
- **Add a `multiqc` rule** that aggregates the per-sample fastp + kallisto reports into a single HTML dashboard. MultiQC (Ewels et al. 2016, *Bioinformatics* 32:3047) is the standard QC aggregator.
- **Add a snakemake-conda mode** with per-rule conda environments specified by `conda:` directives. This is how production pipelines pin tool versions per rule.
- **Parameterize the species** so the pipeline can be re-run on a human dataset by changing one line in `config.yaml`. The yeast-vs-human difference is essentially just the URL and the index-build runtime.

---

## What you learned

- Snakemake's rule-based model is the standard way to express RNA-seq pipelines. The `input → output` declarations are the dependency graph; everything else is bookkeeping.
- Idempotency comes for free: outputs newer than inputs are skipped. This is the property that makes a 10-sample pipeline feel like a 1-sample pipeline after the first run.
- Parallelism comes for free: independent rules run concurrently up to the `--cores` limit. A 3-sample pipeline with `--cores 4` runs three `trim` and three `quant` calls in parallel.
- The DAG is the documentation. `snakemake --dag` produces a graph you can paste into a methods section to show the pipeline structure.
- Real RNA-seq projects use `nf-core/rnaseq` (a Nextflow port of essentially this pipeline at 10x the complexity), but the underlying idea is identical: rules + DAG + parallel execution. Snakemake at this scale is a tractable substitute and is easier to read.

Continue to [Challenge 2 — Compare tools on the same sample](./challenge-02-compare-tools-on-same-sample.md).
