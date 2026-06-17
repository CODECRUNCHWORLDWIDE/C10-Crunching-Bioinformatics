# Lecture 1 — The Reproducibility Crisis and Pipeline-as-Code

> **Educational and research use only.** This lecture introduces the reproducibility crisis in computational biology and the workflow-manager response. Snakemake (Mölder et al. 2021, *F1000Research* 10:33) and Nextflow (Di Tommaso et al. 2017, *Nature Biotechnology* 35:316) are both free open-source tools. The capstone pipeline you ship is an educational artefact only; outputs are not clinical software.

## 1.1 — The 26% number

In 2009, a group of statisticians at the National Institute of Statistical Sciences (Ioannidis et al. 2009, *Nature Genetics* 41:149; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2778103/>) tried to reproduce 18 published microarray-analysis papers. The papers were a representative sample of *Nature Genetics* publications from 2005-2006; the team had access to the raw data deposits and to the methods sections. They reproduced **two of the 18**.

The reasons for the 16 failures were not exotic. The most common: the methods section was insufficient to reconstruct the analysis. Specific parameter values were missing. The software was no longer available at the cited URL. The R version had drifted. The package versions had drifted. A normalization step was described as "quantile normalization" without specifying which of several quantile-normalization variants was used. The team's conclusion was that the field had a reproducibility problem that began at the methods-writing step and compounded through every downstream non-pinning decision.

The 2009 study was not an outlier. Beaulieu-Jones and Greene (2017, *Nature Biotechnology* 35:342; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6103790/>) ran a similar exercise and proposed *continuous analysis*: a CI workflow that re-runs the pipeline on every code commit and verifies the output hash. The 2019 Mangul et al. paper (*PLOS Biology* 17:e3000333; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6592561/>) surveyed 50 bioinformatics tools and found that 26 of them could not be installed from the published instructions without modification. The 26-of-50 number became shorthand for "the canonical bioinformatics tool, picked at random, is 50/50 to bootstrap on a clean machine."

These numbers are the reason the field invented workflow managers. The response is twofold: (a) **pin everything** (the tool version, the reference data, the random seed, the thread count) and (b) **express the pipeline as code** (a `Snakefile` or a `main.nf`, not a shell script). The C10 capstone graduate ships a pipeline whose reproducibility is provable in the byte-identical sense.

## 1.2 — What "reproducible" means

The field uses three terms that are easy to confuse: **reproducible**, **replicable**, **generalizable**.

- **Reproducible:** the same code on the same input data produces the same output. The bar is byte-identical (for deterministic tools) or within-tolerance (for tools with documented randomness).
- **Replicable:** the same scientific question, asked of fresh data with the same protocol, produces the same answer. This is a *scientific* property, not a *computational* property. The capstone is graded on reproducibility, not replicability.
- **Generalizable:** the answer extends beyond the sampled population. The capstone is not graded on generalizability.

The 2018 NASEM consensus report (*Reproducibility and Replicability in Science*, National Academies Press, free at <https://nap.nationalacademies.org/catalog/25303>) is the canonical reference for the three-way distinction. The C10 capstone targets reproducibility — the lowest of the three bars, the one that is achievable by one person on one laptop in one week, and the bar at which the published bioinformatics literature is empirically failing.

## 1.3 — The pipeline-as-code argument

A shell script is an executable description of a pipeline. A `Snakefile` is also an executable description of a pipeline. The difference is what each captures:

A shell script captures the **sequence of commands**. Run them top to bottom; if one fails, you re-start from the top (or from the last `if` checkpoint you remembered to add). The shell script does not know that step 3 depends on step 1's output; it just runs in order. If you change step 1, you re-run everything downstream by hand; if you change step 3, the script does not warn you that the rest is now stale.

A `Snakefile` (or a Nextflow `main.nf`) captures the **dependency graph** of inputs and outputs. Each rule declares its `input`, `output`, `params`, and `shell` command. The workflow manager computes the DAG (directed acyclic graph) at runtime, figures out which rules need to run for the requested target file, and runs only those. If a rule's input is newer than its output, the rule re-runs; if a rule's output is missing, the rule runs; if both are present and the output is up to date, the rule is skipped. The DAG is the engineering payoff: re-running on a small input change costs only the affected rules, not the whole pipeline.

The 2021 Mölder et al. paper formalizes this in §2. The DAG model is the technical core of Snakemake. Nextflow does the same with a channel-based dataflow model (Di Tommaso et al. 2017 §2-§3): processes consume from channels, channels carry tuples, and Nextflow runs each process when its input channel has a tuple waiting. Either model — file-pattern DAG (Snakemake) or channel dataflow (Nextflow) — solves the same problem the shell script does not.

## 1.4 — Snakemake at the rule level

A Snakemake rule looks like this:

```python
rule align_reads:
    input:
        fastq_r1 = "data/{sample}_R1.fastq.gz",
        fastq_r2 = "data/{sample}_R2.fastq.gz",
        index = "ref/GRCh38.chr20.bwa.amb"
    output:
        bam = "results/{sample}.bam",
        bai = "results/{sample}.bam.bai"
    params:
        rg = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime_min = 60
    conda:
        "envs/bwa.yaml"
    log:
        "logs/align_{sample}.log"
    shell:
        "bwa mem -t {threads} -R '{params.rg}' "
        "ref/GRCh38.chr20 {input.fastq_r1} {input.fastq_r2} "
        "| samtools sort -@ {threads} -o {output.bam} - 2> {log} && "
        "samtools index {output.bam}"
```

Read the rule top to bottom:

- **`rule align_reads:`** — the rule's name. Used in error messages and DAG rendering.
- **`input:`** — the files the rule consumes. The `{sample}` is a *wildcard*: it matches any string and Snakemake resolves it from the target filename.
- **`output:`** — the files the rule produces. Snakemake derives the dependency graph from `input` and `output`.
- **`params:`** — variables that are not files. Here the read-group string. `params` flow into the `shell` command via `{params.rg}`.
- **`threads:`** — the parallelism level. Snakemake honours this when running multi-threaded tools and respects a global `--cores` budget across rules.
- **`resources:`** — generic resource declarations. `mem_mb` and `runtime_min` are conventional but the names are user-defined. The cluster scheduler (SLURM, etc.) reads these to size each job.
- **`conda:`** — the Conda environment file for this rule. Snakemake creates and activates the environment automatically when run with `--use-conda`.
- **`log:`** — the path for stderr / log output. Convention: one log per rule per wildcard combination.
- **`shell:`** — the command to run. Triple-quoted strings are common for multi-line commands.

This rule expresses one node of the DAG. A pipeline of 10-30 rules forms the whole DAG. The DAG renders to SVG with one command. The DAG image goes in the capstone repo.

## 1.5 — Nextflow at the process level

The Nextflow equivalent of the same rule:

```groovy
process ALIGN_READS {
    tag "${sample_id}"
    cpus 8
    memory '16 GB'
    time '60.min'

    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.20'

    input:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \\
        GRCh38.chr20 ${fastq_r1} ${fastq_r2} \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    """
}
```

Same logical content, different syntax. The two-axis comparison:

- **Snakemake** is Python with a domain-specific syntax for rules. Wildcards are inline in filenames; Snakemake matches them against the target.
- **Nextflow** is Groovy with a domain-specific syntax for processes. Inputs and outputs are tuples flowing through channels; Nextflow runs the process for each tuple.

Pick whichever syntax you prefer. The Snakemake DAG and the Nextflow DAG of a four-rule pipeline are isomorphic; both run the same shell commands in the same dependency order.

## 1.6 — Running Snakemake

The basic invocation:

```bash
snakemake --cores 4 --use-conda
```

What this does:

- Looks for a `Snakefile` in the current directory.
- Computes the DAG, picks the default target (the first rule's first output, conventionally `rule all:` listing all final outputs).
- Creates / activates Conda environments declared in `conda:` directives.
- Runs the rules in dependency order, up to 4 cores in parallel.

For a target file:

```bash
snakemake --cores 4 --use-conda results/sample01.vcf.gz
```

For the DAG:

```bash
snakemake --dag | dot -Tsvg > dag.svg
```

For a dry run (compute the DAG but do not execute):

```bash
snakemake --cores 4 -n
```

For containerized execution:

```bash
snakemake --cores 4 --use-conda --use-singularity --singularity-args "-B /scratch"
```

The `-B /scratch` flag binds the host's `/scratch` into the container; Singularity / Apptainer does not by default bind paths outside `$HOME` and `/tmp`.

## 1.7 — Running Nextflow

The basic invocation:

```bash
nextflow run main.nf
```

For a profile (a Nextflow conventional way to switch between local / SLURM / cloud):

```bash
nextflow run main.nf -profile slurm
```

For containerized execution:

```bash
nextflow run main.nf -with-singularity 'apptainer-image.sif'
```

For the DAG:

```bash
nextflow run main.nf -with-dag dag.svg
```

Nextflow caches every process invocation in `work/`. Re-running with the same inputs picks up from the cache — equivalent to Snakemake's incremental re-build behaviour.

## 1.8 — The DAG and the target file

The mental model for both workflow managers: you ask for an output file (the *target*), the workflow manager computes the set of rules whose outputs are upstream of the target, it builds the DAG of input-output dependencies among those rules, and it runs each rule in dependency order. The target is the entry point.

For a typical RNA-seq pipeline, the target is something like `results/dge_top30.tsv` (the top 30 differentially-expressed genes). The DAG that builds it:

```
fastq files → trim → quantify → tximport → DESeq2 → top30.tsv
                                                  ↘ volcano.png
```

Snakemake sees `results/dge_top30.tsv`, finds the rule whose output matches, identifies that rule's inputs (`results/deseq2_results.rds`), finds the rule whose output is `results/deseq2_results.rds`, identifies *its* inputs, and so on up the DAG until it reaches files that exist on disk (the input FASTQs). Then it runs the DAG bottom-up.

The DAG is the artefact you ship. It tells a reviewer the pipeline's structure in one image. It is also the thing that makes the pipeline maintainable: changing one rule changes one node of the DAG, not the whole script.

## 1.9 — Wildcards and parallelism

Snakemake's wildcards (and Nextflow's channels) generalize the DAG over multiple samples. The same rule runs for each value of the wildcard:

```python
SAMPLES = ["S01", "S02", "S03", "S04", "S05"]

rule all:
    input:
        expand("results/{sample}.bam", sample=SAMPLES)

rule align_reads:
    input:
        fastq = "data/{sample}.fastq.gz"
    output:
        bam = "results/{sample}.bam"
    shell:
        "bwa mem ref/genome {input.fastq} | samtools sort -o {output.bam}"
```

`expand("results/{sample}.bam", sample=SAMPLES)` produces `["results/S01.bam", "results/S02.bam", ...]`. Snakemake matches each filename against the `align_reads` rule's `output` pattern, sets `{wildcards.sample}` to the matching value, and runs the rule. With `--cores 4`, four samples run in parallel.

The Nextflow equivalent uses a channel:

```groovy
Channel.fromFilePairs("data/*.{R1,R2}.fastq.gz")
       .set { reads_ch }

workflow {
    ALIGN_READS(reads_ch)
}
```

`Channel.fromFilePairs` enumerates the input files and emits one tuple per sample; the `workflow` block wires the channel to the process. The process runs once per tuple.

## 1.10 — Determinism and seeds

A rule is **deterministic** if the same input produces the same output every time. Most bioinformatics tools are deterministic given fixed inputs, with the following common exceptions:

- **Multi-threaded reductions** are not associative in floating point. A 4-thread MD5 over the same byte stream produces the same hash because MD5 is sequential; a 4-thread reduction of read-quality scores does not necessarily produce the same sum because floating-point addition is not associative across thread orderings. Tools like Flye, MaSuRCA, and SPAdes warn about this in their docs. Pin `--threads` to a fixed value if byte-identical output matters.
- **Random seeds.** Tools that sample internally (Flye's read-clustering, IQ-TREE's tree-search, DESeq2's `set.seed`) need their seed pinned. Most tools have a `--seed N` flag; some default to a deterministic seed (good), some default to `$(date +%s)` (bad). Read the docs for each tool you use.
- **File-system order.** `os.listdir` and the shell glob both return files in filesystem-dependent order. On Linux ext4 this is hash-order; on macOS APFS this is creation-order. Always sort your inputs before iterating. Snakemake's `expand` is sorted; the shell's `*.fastq` is not.
- **Time-of-day in output headers.** Some tools embed `date +%s` in their output (e.g. VCF header `##fileDate=2026-05-14`). Two runs on consecutive days produce non-byte-identical VCFs even though the data is identical. Pre-process the output (strip the `##fileDate` line) before hashing.

The capstone's deterministic rules should produce byte-identical output across runs. The capstone's non-deterministic rules should document the source of non-determinism in `report.md`.

## 1.11 — Snakemake and Nextflow side by side

Both workflow managers solve the same problem with different ergonomics. The day-2 lecture (in resources.md) goes deeper into bioconda; the day-3 lecture goes into containers. For now, the takeaways are:

- **Pick one and commit to it for the capstone.** Mixing both inside one pipeline is technically possible (Snakemake can call Nextflow as a shell step and vice-versa) but it doubles the cognitive load and is not graded credit.
- **Both honour the same Conda environment files.** The `environment.yml` you write for Snakemake works for Nextflow with one syntactic adjustment (Snakemake: `conda: "envs/foo.yaml"`; Nextflow: `conda 'bioconda::foo=1.0'` inline or in `nextflow.config`).
- **Both produce a renderable DAG.** Commit the SVG.
- **Both have community libraries** (Snakemake-wrappers; nf-core). Borrow shamelessly.

The choice between them is taste. C10 has used CLI-script + run-info.json patterns through Week 11; the capstone is your first end-to-end workflow-manager pipeline. Snakemake is the recommended starting point for Python-centric students; Nextflow is the recommended starting point for students who want the channel-dataflow mental model.

## 1.12 — Failure modes the workflow manager fixes

A list of reproducibility failures that a workflow manager structurally prevents (and a list of failures it does not):

**Fixed by workflow manager:**

- "I forgot to re-run step 5 after I changed step 2." The DAG figures out the staleness.
- "The pipeline ran on Sam's laptop but not on Lin's." The `environment.yml` and the container fix this.
- "I cannot remember which version of GATK produced the published VCF." The `conda:` directive pins it.
- "I lost the log files for the third run." The `log:` directive captures stderr to a fixed path.
- "Step 3 failed and now I do not know where it died." The workflow manager reports the failing rule, the failing wildcard, the log path, and the command line.

**Not fixed by workflow manager:**

- "I picked the wrong reference build." The workflow manager runs whatever you tell it; if you point it at GRCh37 when the truth set is GRCh38, you get wrong answers in a perfectly-reproducible way.
- "I copied the wrong sample sheet." The workflow manager runs the rules; it does not verify that `data/sample01_R1.fastq.gz` is in fact sample 1.
- "I picked a noise floor that does not match the data." The workflow manager runs FilterMutectCalls with whatever thresholds you specify; it does not tell you the thresholds are wrong.

The workflow manager is a reproducibility tool, not a correctness tool. The capstone needs both: a reproducible pipeline (the workflow manager's job) and a correct pipeline (your job, with the help of orthogonal cross-checks like GIAB truth sets).

## 1.13 — What ships with the capstone

The deposited capstone artefact contains:

- The `Snakefile` (or `main.nf`).
- The `environment.yml` (human-readable) and `environment.lock.txt` (byte-identical).
- The `Singularity.def` (or `Dockerfile`).
- The `.sif` image (or instructions to build it).
- The `run-info.json` (the canonical provenance record).
- The DAG SVG (`dag.svg`).
- The `README.md` with the one-line run command.
- The `report.md` (the one-page write-up).
- The `wrap-up.md` (the C10 retrospective).
- The output figures and tables.
- A `CITATION.cff` (with the deposit DOI).

That is the deposit. The next lecture goes into Conda / bioconda and the container layer; the third lecture into the deposit mechanics and the capstone tracks.

## 1.14 — Recap

- The reproducibility crisis is empirically real: the 2009 *Nat Genet* study reproduced 2 of 18 microarray papers; the 2019 *PLOS Biology* survey found 26 of 50 bioinformatics tools un-installable from published instructions.
- The response is twofold: pin everything and express the pipeline as code.
- Snakemake (Mölder et al. 2021) and Nextflow (Di Tommaso et al. 2017) are the two canonical workflow managers. Both free, both open source, both ship the DAG-rendering tooling.
- A workflow manager structurally prevents staleness, environment drift, and lost logs. It does not prevent wrong reference builds, wrong sample sheets, or wrong noise floors.
- The capstone deposit is the trio of code + data + run instructions, packaged for re-execution by a reviewer.

Lecture 2 covers Conda and bioconda for environment pinning, and Singularity / Apptainer for the container. Lecture 3 covers the four capstone project tracks, the `run-info.json` pattern, the Zenodo DOI deposit, and the wrap-up sidebar.
