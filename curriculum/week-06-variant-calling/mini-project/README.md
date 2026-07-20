# Mini-Project — Variant Call on a Public Small-Genome Dataset

> Build a reproducible BAM → raw VCF → filtered VCF → VEP-annotated VCF pipeline on the Week 5 mini-project BAM (Illumina paired-end resequencing of *E. coli* K-12 MG1655, SRA accession `SRR1770413` aligned against `NC_000913.3`). End with an annotated, filtered VCF, a `bcftools stats` summary, a consequence-distribution table, and a short written interpretation of the variant set.

This is the C10 mini-project that produces a **methods-section-quality variant call set with measured QC**, not just a single-variant demonstration. By the end of it you will have a `call_and_annotate.py` script and a `run.sh` wrapper you can point a recruiter at, a results directory with summary stats and an annotation table, and a write-up that defends each pipeline parameter and names the failure modes you observed.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-06/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 6 section)
└── week-06/
    ├── README.md                   one-page report (~800-1000 words)
    ├── run.sh                      one-command reproduction script
    ├── env.yml                     conda environment file pinning all tool versions
    ├── data/
    │   ├── ecoli.fa                NC_000913.3 reference (symlink to Week 5)
    │   ├── ecoli.fa.fai            samtools faidx index
    │   └── vep_cache/              VEP offline cache for E. coli K-12
    ├── call_and_annotate.py        the variant-calling + annotation pipeline
    ├── summarize_vcf.py            the per-variant summary script
    ├── calls/
    │   ├── SRR1770413.raw.vcf.gz       raw bcftools call output
    │   ├── SRR1770413.filtered.vcf.gz  after bcftools filter
    │   ├── SRR1770413.norm.vcf.gz      after bcftools norm
    │   ├── SRR1770413.vep.vcf          after VEP annotation
    │   └── SRR1770413.pass.vcf.gz      PASS-only subset for downstream use
    └── results/
        ├── bcftools_stats.txt      bcftools stats summary
        ├── consequence_counts.tsv  per-consequence counts from VEP CSQ
        ├── impact_counts.tsv       per-IMPACT category counts
        ├── filter_breakdown.tsv    per-filter-tag counts
        ├── qual_histogram.png      QUAL value distribution
        ├── dp_histogram.png        per-variant DP distribution
        └── variants_per_kb.png     variants-per-kb across the reference
```

By the end you will have a clean, reproducible Week 6 directory you can point a recruiter at — and `call_and_annotate.py` is the kind of pipeline that opens conversations with working bioinformaticians and clinical-sequencing shops.

---

## The dataset

You will work with the **`SRR1770413`** SRA run, *already aligned* against **`NC_000913.3`** in the Week 5 mini-project. The Week 5 output is the input to Week 6: the sorted, indexed, duplicate-marked BAM at `<week-05-dir>/aln/SRR1770413.markdup.bam`.

If you did not complete Week 5's mini-project, you can substitute the smaller Week 5 Exercise 1 BAM (`<week-05-dir>/exercises/aln/lambda.sorted.bam`) — the pipeline is the same in shape, just with a smaller reference and a smaller variant set. Note this in your reproducibility receipt.

Expected variant counts on `SRR1770413` against `NC_000913.3`:

- **Raw**: ~180-220 variants (varies by `bcftools` version and parameters).
- **After hard filters**: ~150-180 (typically 80-90% of raw).
- **After normalization**: same count, but with multiallelic records split and indels left-aligned.
- **By type**: ~85-90% SNPs, ~10-15% indels.
- **Ts/Tv ratio**: ~0.9-1.1 (consistent with random bacterial mutation).
- **By VEP impact**: ~50-60% missense, ~30-40% synonymous, ~5-10% intergenic, a handful of HIGH-impact (stop-gained, frameshift).

If your numbers are wildly different, debug *before* writing the report. The "Common pitfalls" section at the bottom of this file is the first place to look.

### Why E. coli and not human

*E. coli* MG1655 is the standard small-bacterial-reference benchmark for the same three reasons as Week 5: small genome (4.6 Mb) so the full pipeline fits on a laptop and runs in under 5 minutes, well-curated reference, characteristic rRNA-operon multimapper signature that lets you verify the pipeline works. Week 7 moves to RNA-seq (transcriptomics on a different organism); Week 6's *E. coli* is the variant-calling closer of the DNA-sequencing arc.

---

## Rules

- **You may** use bcftools 1.19, GATK 4.5.0.0, ensembl-vep 110, samtools 1.19, pysam 0.22, Biopython 1.83, pandas, matplotlib, numpy, and the standard library.
- **You may** consult Lectures 1 and 2, the bcftools paper (Li 2011, Danecek 2021), the GATK Best Practices, the VCF spec, and your Week 6 exercises and challenge.
- **You may NOT** copy a pre-written variant-calling pipeline from the internet. The point is to *build* the pipeline. Reading the GATK Best Practices for inspiration is fine; copy-pasting a `nf-core/sarek` config is not.
- **You must** cache the VEP database to disk. A second run of `bash run.sh` should not re-download the cache or re-call the variants if the BAM has not changed.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file. The only network access on a fresh run should be the initial VEP cache download (~25 MB for *E. coli*).

---

## Acceptance criteria

- [ ] `week-06/call_and_annotate.py` exports a function `call_and_annotate(ref_fa, bam, out_dir, *, sample, vep_cache_dir, ploidy=1) -> Path` that runs the full BAM → annotated VCF pipeline and returns the path to the PASS-only annotated VCF.
- [ ] The pipeline implements **six stages**:
  1. `samtools faidx` (skip if `.fai` exists).
  2. `bcftools mpileup -Ou -f ref --max-depth 250 -a 'AD,DP,SP' aln.bam | bcftools call -m -v --ploidy 1 -Oz -o raw.vcf.gz -`.
  3. `bcftools filter -s LowQual -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' raw.vcf.gz -Oz -o filtered.vcf.gz`.
  4. `bcftools norm -f ref -m -any filtered.vcf.gz -Oz -o norm.vcf.gz`.
  5. `vep --cache --dir_cache vep_cache/ --species escherichia_coli_str_k_12_substr_mg1655 --vcf --symbol --canonical -i norm.vcf.gz -o vep.vcf`.
  6. `bcftools view -f PASS vep.vcf -Oz -o pass.vcf.gz` + `bcftools index -t`.
- [ ] `week-06/summarize_vcf.py` produces three PNGs and three TSV summary tables:
  - `results/qual_histogram.png` — distribution of QUAL across PASS variants (expect a right-skewed unimodal peak).
  - `results/dp_histogram.png` — distribution of per-variant DP across PASS variants (expect a peak near mean coverage).
  - `results/variants_per_kb.png` — variants per 10 kb across the reference (expect roughly uniform, with notable peaks/dips at specific loci).
  - `results/consequence_counts.tsv` — count per VEP consequence term (`missense_variant`, `synonymous_variant`, etc.).
  - `results/impact_counts.tsv` — count per VEP IMPACT category (`LOW`, `MODERATE`, `HIGH`, `MODIFIER`).
  - `results/filter_breakdown.tsv` — count per FILTER tag (PASS, LowQual, etc.).
- [ ] `week-06/README.md` is a one-page (≤ 1,000 word) report containing:
  - One-sentence description of the dataset, the reference, and the pipeline's six stages.
  - Methods section in C10 voice: every tool pinned ("bcftools 1.19", "ensembl-vep 110", "GATK 4.5.0.0"), every parameter explicit (`bcftools call -m -v --ploidy 1`, `bcftools filter -s LowQual -e '...'`, `vep --cache --dir_cache vep_cache/`).
  - Quantitative findings: raw variant count, PASS count, SNP/indel breakdown, Ts/Tv ratio, depth and QUAL distributions (mean, median, IQR), consequence breakdown by VEP CSQ field. Each value to one decimal place.
  - An interpretation of the consequence distribution: how many missense vs synonymous, ratio of the two, fraction in coding regions, any HIGH-impact variants and their genes.
  - A failure-modes section listing every QC signal that surprised you, with one-sentence diagnoses.
  - A reproducibility receipt block.
- [ ] `week-06/run.sh` runs end-to-end on a fresh clone:
  - Activates the conda env from `env.yml`.
  - Confirms the Week 5 BAM exists at the expected path; errors out with a clear message if not.
  - Builds the VEP cache if not already on disk.
  - Runs `call_and_annotate.py` to call, filter, normalize, and annotate.
  - Runs `summarize_vcf.py` to produce the three PNGs and three TSVs.
- [ ] `week-06/env.yml` pins every tool to an exact version:
  ```yaml
  name: c10-week-06
  channels:
    - conda-forge
    - bioconda
  dependencies:
    - python=3.11
    - numpy=1.26.4
    - matplotlib
    - pandas
    - biopython=1.83
    - samtools=1.19
    - bcftools=1.19
    - pysam=0.22
    - gatk4=4.5.0.0
    - ensembl-vep=110
    - pip
    - pip:
        - cyvcf2==0.30
  ```
- [ ] `week-06/results/bcftools_stats.txt` contains a `bcftools stats` summary with the variant counts and Ts/Tv ratio.
- [ ] `week-06/results/consequence_counts.tsv` is a sorted-by-count TSV with at least 5 distinct consequence terms.
- [ ] `week-06/calls/SRR1770413.pass.vcf.gz` is a bgzipped, tabix-indexed VCF with FILTER=PASS variants only.
- [ ] The repo passes a fresh-clone test: `git clone`, `cd week-06`, `bash run.sh` reproduces everything in `results/` (modulo a ±5% difference in counts if `bcftools`/`GATK`/`VEP` versions have moved).

---

## Suggested order of operations

### Phase 1 — Environment setup (~30 min)

1. Create `week-06/env.yml` (see acceptance criteria).
2. `conda env create -f week-06/env.yml`. Activate it. Confirm: `bcftools --version | head -1` → `bcftools 1.19`, `vep --help | head -2` → `Versions: ensembl ... ensembl-funcgen ... ensembl-io ... ensembl-variation`, `python -c "import pysam; print(pysam.__version__)"` → `0.22`.
3. Commit: `Week 6 env.yml`.

### Phase 2 — Confirm the Week 5 BAM is reachable (~10 min)

1. Symlink (or copy) the Week 5 reference and BAM into the Week 6 directory:
   ```bash
   ln -sf ../../week-05/data/ecoli.fa data/ecoli.fa
   ln -sf ../../week-05/aln/SRR1770413.markdup.bam aln/
   ln -sf ../../week-05/aln/SRR1770413.markdup.bam.bai aln/
   ```
2. Verify: `samtools quickcheck aln/SRR1770413.markdup.bam && echo OK` should print `OK`.
3. Commit: `Link Week 5 reference and BAM`.

### Phase 3 — Build the variant-calling pipeline (~2 hours)

1. Write `call_and_annotate.py` with the `call_and_annotate(...)` function. Outline:
   ```python
   def call_and_annotate(
       ref_fa: Path, bam: Path, out_dir: Path, *,
       sample: str = "ecoli_K12",
       vep_cache_dir: Path,
       vep_species: str = "escherichia_coli_str_k_12_substr_mg1655",
       ploidy: int = 1,
   ) -> Path:
       out_dir.mkdir(parents=True, exist_ok=True)
       # Stage 1: index reference (idempotent).
       if not Path(f"{ref_fa}.fai").exists():
           subprocess.run(["samtools", "faidx", str(ref_fa)], check=True)
       # Stage 2: mpileup | call.
       raw_vcf = out_dir / "raw.vcf.gz"
       cmd = (
           f"bcftools mpileup -Ou -f {ref_fa} --max-depth 250 "
           f"-a 'AD,DP,SP' {bam} "
           f"| bcftools call -m -v --ploidy {ploidy} "
           f"-Oz -o {raw_vcf} -"
       )
       subprocess.run(cmd, shell=True, check=True)
       subprocess.run(["bcftools", "index", "-t", str(raw_vcf)], check=True)
       # Stage 3: filter.
       filt_vcf = out_dir / "filtered.vcf.gz"
       expr = "QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60"
       subprocess.run([
           "bcftools", "filter", "-Oz", "-o", str(filt_vcf),
           "-s", "LowQual", "-e", expr, str(raw_vcf),
       ], check=True)
       subprocess.run(["bcftools", "index", "-t", str(filt_vcf)], check=True)
       # Stage 4: norm.
       norm_vcf = out_dir / "norm.vcf.gz"
       subprocess.run([
           "bcftools", "norm", "-f", str(ref_fa), "-m", "-any",
           "-Oz", "-o", str(norm_vcf), str(filt_vcf),
       ], check=True)
       subprocess.run(["bcftools", "index", "-t", str(norm_vcf)], check=True)
       # Stage 5: VEP.
       vep_vcf = out_dir / "vep.vcf"
       subprocess.run([
           "vep",
           "--input_file", str(norm_vcf),
           "--output_file", str(vep_vcf),
           "--species", vep_species,
           "--cache", "--dir_cache", str(vep_cache_dir),
           "--vcf", "--symbol", "--canonical", "--force_overwrite",
       ], check=True)
       # Stage 6: PASS-only.
       pass_vcf = out_dir / "pass.vcf.gz"
       subprocess.run([
           "bcftools", "view", "-f", "PASS", "-Oz",
           "-o", str(pass_vcf), str(vep_vcf),
       ], check=True)
       subprocess.run(["bcftools", "index", "-t", str(pass_vcf)], check=True)
       return pass_vcf
   ```
2. Add input validation: every path must exist; ploidy must be ≥ 1; the sample name must be non-empty.
3. Test on a tiny subset first (e.g., restrict mpileup to a small region with `-r NC_000913.3:1-100000`) to verify the pipeline works before running on the full BAM.
4. Commit: `call_and_annotate.py end to end on test region`.

### Phase 4 — Install the VEP cache (~10 min)

1. `vep_install --AUTO c --SPECIES escherichia_coli_str_k_12_substr_mg1655 --CACHEDIR data/vep_cache/`.
2. Confirm: `ls data/vep_cache/escherichia_coli_str_k_12_substr_mg1655_*` should show the cache files (~25 MB).
3. Commit: `VEP cache installed for E. coli K-12 MG1655`.

### Phase 5 — Run on the full dataset (~5 min wall-clock)

1. `python call_and_annotate.py data/ecoli.fa aln/SRR1770413.markdup.bam calls/`.
2. Confirm `calls/SRR1770413.raw.vcf.gz`, `calls/SRR1770413.filtered.vcf.gz`, `calls/SRR1770413.norm.vcf.gz`, `calls/SRR1770413.vep.vcf`, and `calls/SRR1770413.pass.vcf.gz` all exist.
3. Confirm with `bcftools stats calls/SRR1770413.pass.vcf.gz | head -50`:
   - Total variants ~150-180.
   - SNPs ~85-90% of total.
   - Indels ~10-15% of total.
   - Ts/Tv ~0.9-1.1.
4. Commit: `Called and annotated SRR1770413 variants`.

### Phase 6 — Summary statistics and plots (~1 hour)

1. Write `summarize_vcf.py`:
   - Read `calls/SRR1770413.vep.vcf` (the full VEP output, before PASS-only subset).
   - Parse the `CSQ` field per variant (the CSQ Format spec is in the VCF header `##INFO=<ID=CSQ,...>`).
   - Tabulate per-consequence and per-IMPACT counts; write as TSV.
   - Tabulate per-FILTER-tag counts; write as TSV.
   - Plot QUAL histogram, DP histogram, and per-10-kb-window variant density.
2. Run `bcftools stats calls/SRR1770413.pass.vcf.gz > results/bcftools_stats.txt`.
3. Commit: `Variant set QC summaries and plots`.

### Phase 7 — Write the report (~2 hours)

Open `week-06/README.md`. Structure:

```
# Week 6 — Variant calling of SRR1770413 against NC_000913.3

## Dataset and biological question
- Input BAM: SRR1770413 aligned against NC_000913.3 (E. coli K-12
  MG1655) by the Week 5 mini-project. ~750,000 read pairs at ~49x
  mean coverage, sorted, indexed, duplicate-marked.
- Reference: NC_000913.3, 4,641,652 bp.
- Question: what variants distinguish this re-sequenced isolate from
  the canonical MG1655 reference? What is their functional impact?

## Methods
[200 words. bcftools 1.19 with the canonical mpileup | call pipeline,
--ploidy 1 for haploid bacteria, hard-filter recipe QUAL<30 ||
INFO/DP<10 || INFO/MQ<40 || INFO/SP>60. bcftools norm -m -any for
indel left-alignment and multiallelic split. VEP 110 with the offline
cache for E. coli K-12. Implementation in Python via subprocess;
no nf-core or snakemake.]

## Findings
- Raw variants: NNN
- PASS variants: NNN (NN.N% of raw)
- LowQual variants: NN
- SNPs (PASS): NNN (NN.N% of PASS)
- Indels (PASS): NN (NN.N% of PASS)
- Ts/Tv (PASS SNPs): 0.NN (consistent with bacterial random mutation)
- Mean QUAL (PASS): NNN.N
- Mean DP (PASS): NN.N
- Variants per Mb: NN.N

## Functional annotation
- VEP consequence distribution (PASS variants):
  - missense_variant: NN (NN.N%)
  - synonymous_variant: NN (NN.N%)
  - intergenic_variant: NN (NN.N%)
  - [other terms ...]
- IMPACT distribution:
  - LOW: NN
  - MODERATE: NN
  - HIGH: NN (list specific variants if any)
  - MODIFIER: NN
- Notable HIGH-impact variants: [gene name, type, biological note]

## Failure modes observed
- [Any low-confidence variants in repetitive regions; reasoning.]
- [Any indel clusters near homopolymer runs; reasoning.]

## Reproducibility
[the receipt block, below]
```

Aim for **clear, quantitative, undramatic**. C10 voice. No "many variants found"; give numbers.

Commit: `Week 6 report v1`.

### Phase 8 — Reproducibility receipt and polish (~30 min)

Add the receipt block at the bottom of `week-06/README.md`:

```
+--------------------------------------------------------------------+
|  REPRODUCIBILITY                                                   |
|                                                                    |
|  Data source:   Week 5 mini-project BAM                            |
|                  SRR1770413.markdup.bam (built YYYY-MM-DD)         |
|                  + RefSeq NC_000913.3 (E. coli K-12 MG1655)        |
|  Pipeline:      bash run.sh                                        |
|  Container:     conda env (env.yml pinned to bcftools 1.19,        |
|                  samtools 1.19, gatk 4.5.0.0, ensembl-vep 110,     |
|                  pysam 0.22, Python 3.11)                          |
|  VEP cache:     escherichia_coli_str_k_12_substr_mg1655 (release   |
|                  110, ~25 MB, installed YYYY-MM-DD)                |
|  Command:       bash run.sh                                        |
|  Wall time:     ~3 min (M1, 16 GB RAM, 4 threads, network only for |
|                  first run for VEP cache)                          |
+--------------------------------------------------------------------+
```

Sanity-check: every file path in the README is correct, every URL works, `bash run.sh` succeeds on a fresh checkout. Commit: `Mini-project reproducibility receipt`.

### Phase 9 — Update the repo root README (~15 min)

Add a Week 6 section to your portfolio repo's top-level `README.md` linking to `week-06/README.md` with a one-paragraph summary. Commit: `Week 6 entry in portfolio README`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Correctness | 25% | PASS rate > 80%, Ts/Tv 0.9-1.1, consequence distribution has at least 5 distinct terms. Pipeline runs end to end without manual intervention. |
| Reproducibility | 20% | `bash run.sh` works on a fresh clone. Versions pinned. VEP cache cached. |
| Code quality | 15% | `call_and_annotate.py` has module docstrings, structured return types, input validation, subprocess errors propagated. Idempotent on re-run. |
| Quantitative report | 15% | Every claim has a number. C10 voice throughout. Consequence interpretation names specific genes for any HIGH-impact variants. |
| Plot quality | 15% | All three plots are labelled, captioned, biologically interpretable. QUAL and DP distributions show clear unimodal peaks. Variant density plot identifies any hot/cold regions. |
| Voice and precision | 10% | Reads like a methods section. No "many" without numbers, no "high quality" without thresholds. |

---

## What this prepares you for

- **Week 7** leaves DNA behind for RNA-seq. The pysam patterns from Weeks 5 and 6 and the pandas patterns from Week 4 both come back; the variant-calling machinery is parked.
- **Week 8** runs DESeq2 on the read counts from Week 7. The hard-filter discipline from Week 6 ("every threshold is a trade-off; document the threshold and the alternatives") applies directly to FDR control in differential expression.
- **Week 11** asks you to convert this pipeline to Snakemake. The current `run.sh` is the rehearsal: every step is a `subprocess.run(...)` invocation that can become a Snakemake rule.
- **Week 12 capstone** will likely include some form of variant calling or population analysis. The reproducibility template you build here is the model.

---

## Common pitfalls

A short list, from instructor experience:

- **Forgetting `--ploidy 1` for *E. coli*.** Every variant gets called as heterozygous; QUAL values are half their expected magnitude; the PASS set is unexpectedly small. Always set ploidy explicitly.
- **Skipping `bcftools norm`.** Indels look different in the output than in the input VCF; cross-VCF comparison (Problem 2, Challenge 1) reports false disagreements. Always normalize before any cross-tool comparison.
- **Running VEP without `--canonical`.** The output CSQ field has one entry per overlapping transcript; without `--canonical`, you do not know which entry corresponds to the gene's canonical transcript, so the consequence-distribution summary is ambiguous. Always pass `--canonical`.
- **Hardcoding paths.** `run.sh` should compute paths relative to its own location (`SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)`) so a fresh clone in any directory works.
- **Committing the VEP cache.** The cache for *E. coli* is small (~25 MB) but committing it inflates the repo unnecessarily. Add `data/vep_cache/` to `.gitignore` and rebuild on `run.sh`. Same for the BAM and BAI.

---

## Submission

When done:

1. Confirm `bash run.sh` works on a fresh clone (delete `calls/` and re-run; it should rebuild from the linked Week 5 BAM and cached reference + VEP database).
2. Confirm `week-06/README.md` renders cleanly on GitHub.
3. Confirm `results/qual_histogram.png`, `results/dp_histogram.png`, `results/variants_per_kb.png` are committed and meaningful.
4. Confirm `results/bcftools_stats.txt`, `results/consequence_counts.tsv`, `results/impact_counts.tsv`, `results/filter_breakdown.tsv` are committed.
5. Confirm `calls/SRR1770413.pass.vcf.gz` and its `.tbi` index are committed.
6. Push.
7. Open Week 7 — only after the report is committed and the failure-modes section explains every QC surprise.

---

## Resources

- [bcftools howtos](https://samtools.github.io/bcftools/howtos/) — the command-line reference for the toolchain.
- [bcftools paper (Danecek 2021)](https://academic.oup.com/gigascience/article/10/2/giab008/6137722) — the modern reference.
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471) — the canonical hard-filter recipe.
- [Ensembl VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html) — every command with examples.
- [VCFv4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) — the ~40-page canonical reference.
- [pysam documentation](https://pysam.readthedocs.io/) — the Python BAM/VCF reader API.
- [Lecture 1 — From BAM to VCF with bcftools](../lecture-notes/01-from-bam-to-vcf-bcftools.md)
- [Lecture 2 — GATK Best Practices and hard filters](../lecture-notes/02-gatk-best-practices-and-hard-filters.md)
- [Challenge 1 — Compare bcftools and GATK on the same BAM](../challenges/challenge-01-compare-callers.md)
- [resources.md](../resources.md) — full week resource list.
