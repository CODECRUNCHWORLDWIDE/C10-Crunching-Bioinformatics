# Mini-Project — End-to-end variant annotation and interpretation pipeline

> **Educational and research use only.** This mini-project builds an end-to-end variant annotation pipeline. Its output looks like the output of a clinical pipeline, but it is not one. Clinical interpretation requires a CLIA-certified laboratory and a board-certified clinical geneticist. The disclaimer that opens every Week 8 file applies. **You must not act on the report this pipeline produces, nor share it with anyone who might.**

Build a reproducible variant annotation pipeline that takes a VCF (the bundled `data/demo.vcf` of ~50 hand-picked variants), annotates each variant with VEP and SnpEff, queries gnomAD and ClinVar for the population and clinical-knowledge axes, computes SIFT/PolyPhen-2 scores via VEP plugins, applies the mechanically computable subset of the ACMG/AMP 2015 criteria, and emits a per-variant interpretation report as CSV plus a standalone HTML page. End with: an annotated VCF, a CSV report, an HTML report, a run-info JSON, and a 600-800 word write-up that defends every pipeline parameter and names every limit.

This is the C10 mini-project that produces a **methods-section-quality variant interpretation report with measured provenance**, not just a one-shot demonstration. By the end of it you will have an `annotate_pipeline.py` script and a `run.sh` wrapper, a results directory with the report and the run-info, and a write-up that defends every parameter and explicitly names what the report cannot do.

**Estimated time:** 7 hours (split across Wednesday, Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-08/mini-project/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 8 section)
└── week-08/
    └── mini-project/
        ├── README.md               one-page report (~600-800 words)
        ├── run.sh                  one-command reproduction script
        ├── env.yml                 conda environment file pinning all tool versions
        ├── data/
        │   ├── demo.vcf            input VCF (~50 variants)
        │   ├── lof_gene_set.tsv    curated LOF-disease-mechanism gene list
        │   └── samples.tsv         per-variant metadata (gene hint, expected category)
        ├── annotate_pipeline.py    the orchestration script
        ├── acmg_classifier.py      the ACMG classifier (from Challenge 1)
        ├── starter.py              skeleton implementations with TODOs
        ├── caches/                 gnomAD + ClinVar response caches (gitignored)
        └── results/
            ├── demo.vep.vcf           VEP-annotated VCF
            ├── demo.snpeff.vcf        SnpEff-annotated VCF
            ├── variants.csv           per-variant report
            ├── variants.html          standalone HTML report
            ├── run-info.json          versions and run date
            └── disagreement_report.txt    VEP-vs-SnpEff impact-tier disagreements
```

By the end you will have a clean, reproducible Week 8 directory you can point a recruiter at — and `annotate_pipeline.py` is the kind of pipeline that opens conversations with working bioinformaticians and biotech / clinical-diagnostics shops, *as long as* you discuss it with the limits front and centre.

---

## The dataset

You will work with the bundled `data/demo.vcf` — a **synthetic teaching VCF** of ~50 variants. The composition:

- **~30 variants** drawn from the GIAB NA12878 truth set (Zook et al. 2014, *Nature Biotechnology* 32:246, free at <https://www.nature.com/articles/nbt.2835>) on chromosomes 1, 3, 7, 13, 17, 19. These are real variants confirmed by orthogonal methods (Sanger sequencing, array, manual review). The GIAB truth set is the canonical benchmark for variant calling.
- **~20 hand-picked variants** from ClinVar with known Pathogenic or Likely Pathogenic classifications across well-characterized disease genes (BRCA1, BRCA2, CFTR, MLH1, MSH2, APC, TP53, FBN1, LDLR, HEXA, ATM, RB1, NF1, PTEN, VHL, COL3A1, RYR1, SCN1A, MEN1). The hand-picked variants are *for educational examples only* — they are real published pathogenic variants but they are not from a patient sample and the inclusion of a variant in this teaching set does not constitute clinical information about any individual.

The file is ~10 KB. Every variant is bi-allelic and has a populated FILTER=PASS column. The QUAL column is set to 500 across the board (the variants are pre-validated; QUAL is informational).

Expected coverage by the annotation pipeline:

- All ~50 variants should produce at least one VEP CSQ record and at least one SnpEff ANN record.
- ~30-35 variants should have a gnomAD entry (the GIAB-derived subset; the hand-picked Pathogenic subset is mostly absent from gnomAD by selection).
- ~20-25 variants should have a ClinVar entry (the hand-picked subset plus a few of the GIAB-derived ones that overlap known SNPs).
- After the ACMG classifier:
  - ~4-6 variants Likely_pathogenic (PVS1 + PM2 hits in known LOF genes).
  - ~12-15 variants Likely_benign or Benign (BA1 stand-alone or BS+BP combinations).
  - ~25-30 variants Uncertain_significance (insufficient evidence; expected for most teaching variants).

### Why this composition

The didactic point is that *most variants in a typical exome are VUS*, even when you have VEP + SnpEff + gnomAD + ClinVar + SIFT + PolyPhen all firing. The pipeline produces a defensible per-variant row for every variant. The ACMG classifier produces a defensible classification for ~30% of them. The remaining ~70% need manual review by a clinical geneticist — which is exactly the point. The pipeline is the input to the geneticist's review, not a replacement for it.

---

## Rules

- **You may** use VEP 110.1, SnpEff 5.2, cyvcf2 0.30.28, pysam 0.22.1, requests 2.32, pandas 2.2, numpy 1.26, jinja2 3.1, and the standard library.
- **You may** consult Lectures 1, 2, 3, the VEP paper (McLaren et al. 2016), the SnpEff paper (Cingolani et al. 2012), the ACMG/AMP 2015 paper (Richards et al. 2015), the ClinVar paper (Landrum et al. 2018), the gnomAD paper (Karczewski et al. 2020), and the Week 8 exercises and challenges.
- **You may NOT** copy a pre-written variant interpretation pipeline from the internet. The point is to *build* the pipeline. Reading the `Ensembl/VEP_plugins` README for inspiration is fine; copy-pasting a complete interpretation pipeline is not.
- **You must** cache the VEP cache to disk once. A second run of `bash run.sh` should not re-download the cache.
- **You must** cache the gnomAD and ClinVar query responses in `caches/`. A second run should not re-hit either API.
- **You must** record the run-info JSON with database versions and run date. Without these, the report is not reproducible.
- **The report must contain the educational-use-only disclaimer prominently.** If you are tempted to remove it because "the pipeline is good," that is the moment to leave it in.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file.

---

## Acceptance criteria

- [ ] `mini-project/annotate_pipeline.py` exports a function `annotate(vcf_path: Path, ref_dir: Path, out_dir: Path) -> Path` that runs the full pipeline and returns the path to the HTML report.
- [ ] The pipeline implements **seven stages**:
  1. **Validate input.** Confirm the input VCF exists, has a `##reference=GRCh38` header, and at least one variant.
  2. **Annotate with VEP.** Run `vep ... --canonical --mane --sift b --polyphen b --af --af_gnomadg --af_gnomade --check_existing` against the cached VEP database. Skip if the output VCF is newer than the input.
  3. **Annotate with SnpEff.** Run `snpEff -v GRCh38.105 ...` for the second-opinion impact tier. Compare to VEP and emit `disagreement_report.txt` for variants where the impact tier disagrees.
  4. **Query gnomAD.** For each variant, query the gnomAD GraphQL API for population frequency. Use the shelve-based cache from Lecture 2 to avoid re-querying.
  5. **Query ClinVar.** For each variant, query the local ClinVar VCF release via tabix (or NCBI E-utilities as a fallback) for the clinical assertion.
  6. **Apply ACMG.** Run the classifier from Challenge 1 (`acmg_classifier.py`) on the integrated frame. Compute the 8 mechanically computable criteria; mark the rest as skipped.
  7. **Emit reports.** Write `variants.csv`, `variants.html`, and `run-info.json`.
- [ ] `annotate_pipeline.py` produces:
  - `results/variants.csv` — per-variant frame, one row per input variant, at least 18 columns including the ACMG fields.
  - `results/variants.html` — standalone HTML report with the disclaimer, the metadata box, the per-variant table colour-coded by impact tier, and the skipped-criteria section.
  - `results/run-info.json` — versions, dates, and the disclaimer text.
  - `results/disagreement_report.txt` — text report of variants where VEP and SnpEff disagree on impact tier (typically ~5% of variants).
- [ ] `mini-project/README.md` is a one-page (~600-800 word) report containing:
  - One-sentence description of the dataset, the reference, and the pipeline's seven stages.
  - Methods section in C10 voice: every tool pinned ("VEP 110.1", "SnpEff 5.2", "ClinVar release 2024-09-01", "gnomAD v4.1.0"), every parameter explicit.
  - Results section in C10 voice: per-axis counts (`n_variants` with a gnomAD record, `n_variants` with a ClinVar record, `n_variants` with each impact tier), the ACMG classification counts (Pathogenic / Likely_pathogenic / VUS / Likely_benign / Benign), the VEP-vs-SnpEff agreement rate, and 3 hand-picked variants with their full row of evidence.
  - Discussion section: 150-250 words on the limits of the pipeline. What is *not* in the report? What would a clinical pipeline add? Why does the disclaimer exist?
- [ ] `run.sh` is a single bash script that, given a fresh checkout + `conda env create -f env.yml`, reproduces the entire pipeline from scratch in under 5 minutes (after the VEP cache and ClinVar VCF are downloaded).
- [ ] The repo is **public** and at least one classmate or instructor has been added as a collaborator.
- [ ] **Most importantly**: the report's disclaimer is prominent and the discussion section is honest about the limits.

---

## Suggested approach (rough timeline)

### Wednesday (1 hour)

1. (15 min) `git clone`, set up `mini-project/` directory.
2. (30 min) Write `env.yml` with pinned tool versions; create conda env; verify each tool is on the PATH (VEP, SnpEff, cyvcf2).
3. (15 min) Read this README end to end and the Challenge 1 README. Sketch the seven-stage flow on paper.

### Thursday (2 hours)

1. (15 min) Download (or symlink) the demo VCF, the LOF gene set, and the ClinVar VCF release.
2. (30 min) Run VEP and SnpEff by hand on the demo VCF; verify the CSQ and ANN fields are present.
3. (60 min) Write `annotate_pipeline.py` stages 1-3 (validate, VEP, SnpEff). Save the disagreement report.
4. (15 min) Commit and push.

### Friday (2 hours)

1. (45 min) Write stage 4 (gnomAD query with cache).
2. (45 min) Write stage 5 (ClinVar query via local VCF tabix or E-utilities).
3. (30 min) Verify the per-variant frame has all the columns; sanity-check a few rows by hand.

### Saturday (2 hours)

1. (60 min) Integrate `acmg_classifier.py` (from Challenge 1) into stage 6. Run on the full frame and tabulate the classifications.
2. (30 min) Write stage 7 (CSV + HTML output). Build the HTML template, embed the disclaimer.
3. (30 min) Write the README. Lead with the methods; then the results (with specific numbers); then the discussion with the limits.

---

## Methods recipe (the pipeline you must reproduce in `run.sh`)

```bash
#!/usr/bin/env bash
set -euo pipefail

# Stage 0: data sources and environment.
mkdir -p data results caches

# Stage 0a: VEP cache (~25 GB, one-time install).
if [ ! -d "${HOME}/.vep/homo_sapiens_vep_110_GRCh38" ]; then
    vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 \
        --CACHEDIR ${HOME}/.vep --NO_HTSLIB --CACHE_VERSION 110
fi

# Stage 0b: SnpEff database (~2.5 GB, one-time install).
if [ ! -d "${HOME}/snpEff/data/GRCh38.105" ]; then
    snpEff download -v GRCh38.105
fi

# Stage 0c: ClinVar VCF release (~200 MB, refresh every 2 weeks).
if [ ! -f data/clinvar.vcf.gz ]; then
    curl -sLo data/clinvar.vcf.gz \
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    curl -sLo data/clinvar.vcf.gz.tbi \
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi

# Stage 1: VEP annotation.
if [ ! -f results/demo.vep.vcf ]; then
    vep \
        -i data/demo.vcf \
        -o results/demo.vep.vcf \
        --vcf \
        --cache --dir_cache ${HOME}/.vep --offline \
        --species homo_sapiens --assembly GRCh38 --cache_version 110 \
        --sift b --polyphen b \
        --canonical --mane --symbol --biotype \
        --af --af_gnomadg --af_gnomade \
        --check_existing \
        --fork 4 --force_overwrite --no_stats
fi

# Stage 2: SnpEff annotation.
if [ ! -f results/demo.snpeff.vcf ]; then
    snpEff -v -csvStats results/snpeff_stats.csv \
           GRCh38.105 data/demo.vcf > results/demo.snpeff.vcf
fi

# Stages 3-7: gnomAD + ClinVar + ACMG + reports.
python annotate_pipeline.py \
    --vep-vcf results/demo.vep.vcf \
    --snpeff-vcf results/demo.snpeff.vcf \
    --clinvar-vcf data/clinvar.vcf.gz \
    --lof-gene-set data/lof_gene_set.tsv \
    --cache-dir caches \
    --out-csv results/variants.csv \
    --out-html results/variants.html \
    --out-run-info results/run-info.json \
    --out-disagreement results/disagreement_report.txt

echo "Done. HTML report at results/variants.html."
```

---

## Expected results (your numbers should be within ~10% of these)

### Per-axis variant counts

| Axis                                 | Count       |
|--------------------------------------|------------:|
| Total variants in `demo.vcf`        | 50          |
| Variants with a VEP CSQ record       | 50 (100%)   |
| Variants with a SnpEff ANN record    | 50 (100%)   |
| Variants with a gnomAD v4.1 entry    | ~33 (66%)   |
| Variants with a ClinVar record       | ~22 (44%)   |
| Variants with SIFT score (missense)  | ~25 (50%)   |
| Variants with PolyPhen score         | ~25 (50%)   |
| VEP-vs-SnpEff impact-tier agreement  | ~47 (94%)   |

### Impact-tier counts

| Impact   | Count |
|----------|------:|
| HIGH     | ~5    |
| MODERATE | ~28   |
| LOW      | ~9    |
| MODIFIER | ~8    |

### ACMG classification counts

| Classification              | Count |
|-----------------------------|------:|
| Pathogenic                  | 0-1   |
| Likely_pathogenic           | 4-6   |
| Uncertain_significance      | 28-32 |
| Likely_benign               | 4-6   |
| Benign                      | 8-12  |

### Three hand-picked variants you should be able to defend

| chrom | pos       | gene  | consequence            | impact   | clinvar_clnsig          | gnomad_popmax_af | acmg_classification |
|-------|----------:|-------|------------------------|----------|--------------------------|-----------------:|---------------------|
| chr17 | 43094077  | BRCA1 | missense_variant       | MODERATE | Conflicting              | 0.00018          | VUS                 |
| chr5  | 112815473 | APC   | frameshift_variant     | HIGH     | Pathogenic (4-star)      | absent           | Likely_pathogenic   |
| chr7  | 117559590 | CFTR  | synonymous_variant     | LOW      | Benign                   | 0.18             | Benign              |

For each of these three, your `variants.html` should show every column of evidence: the gene, the consequence, the impact, the rsID, the SIFT/PolyPhen scores where applicable, the gnomAD AF (overall and popmax), the ClinVar CLNSIG and CLNREVSTAT, and the ACMG classification with the evidence set.

---

## Write-up template (paste into mini-project/README.md, fill in your numbers)

```markdown
# Week 8 Mini-Project — Variant Annotation and Interpretation

> Educational and research use only. This report is not a clinical interpretation.
> Variant interpretation in a clinical context requires CLIA-certified laboratory
> review and board-certified clinical geneticist sign-off. Do not act on this report.

## Dataset

A 50-variant synthetic teaching VCF combining ~30 GIAB NA12878 truth-set variants
(Zook et al. 2014, *Nat Biotechnol* 32:246) with ~20 hand-picked ClinVar-curated
Pathogenic/Likely-Pathogenic variants across well-characterized disease genes
(BRCA1, BRCA2, CFTR, MLH1, MSH2, APC, TP53, FBN1, LDLR, HEXA). All variants are
bi-allelic, PASS-quality, on GRCh38.

## Methods

Variants were annotated against the **Ensembl release 110** gene model using
**VEP 110.1** (McLaren et al. 2016, *Genome Biol* 17:122) with `--canonical
--mane --sift b --polyphen b --af --af_gnomadg --af_gnomade --check_existing`,
and independently annotated with **SnpEff 5.2** (Cingolani et al. 2012, *Fly*
6:80) against the `GRCh38.105` database for a second-opinion impact tier.
Population frequencies were retrieved from **gnomAD v4.1.0** (Karczewski et al.
2020, *Nature* 581:434) via the public GraphQL API at
`https://gnomad.broadinstitute.org/api`, using a local shelve cache for
idempotent re-runs. Clinical assertions were retrieved from the **ClinVar VCF
release 2024-09-01** (Landrum et al. 2018, *Nucleic Acids Res* 46:D1062) via
tabix lookups on the local indexed VCF. The mechanically computable subset of
the **ACMG/AMP 2015 criteria** (Richards et al. 2015, *Genet Med* 17:405, free
at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/>) was applied:
PVS1 (with the ClinGen Haploinsufficiency Tier 3 gene list), PM2 (popmax <
0.0001), PM4, PP3 (SIFT < 0.05 AND PolyPhen > 0.85), BA1 (popmax > 0.05), BS2,
BP4, BP7. The remaining 20 criteria — PS2, PS3, PS4, PM3, PM5, PM6, PP1, PP2,
PP4, PP5, BS1, BS3, BS4, BP1, BP2, BP5, BP6 — are explicitly flagged as
not evaluated in the per-variant `skipped_criteria` column.

## Results

Per-axis variant counts:

| Axis | Count |
|------|------:|
| Total variants | <YOUR N> |
| Variants with gnomAD entry | <YOUR N> |
| Variants with ClinVar entry | <YOUR N> |
| VEP-vs-SnpEff impact-tier agreement | <YOUR %> |

Impact-tier counts:

| Impact | Count |
|--------|------:|
| HIGH | <X> |
| MODERATE | <Y> |
| LOW | <Z> |
| MODIFIER | <W> |

ACMG classification counts:

| Classification | Count |
|----------------|------:|
| Pathogenic | <X> |
| Likely_pathogenic | <Y> |
| Uncertain_significance | <Z> |
| Likely_benign | <W> |
| Benign | <V> |

Three hand-picked variants (full per-variant rows in `variants.html`):

| Gene | Consequence | gnomAD popmax | ClinVar CLNSIG | ACMG |
|------|-------------|---------------:|----------------|------|
| BRCA1 | missense   | <X>            | Conflicting     | VUS  |
| APC   | frameshift | absent         | Pathogenic 4*   | LP   |
| CFTR  | synonymous | 0.18           | Benign          | Benign |

## Discussion

The pipeline produces a defensible per-variant row for every input variant.
The ACMG classifier produces a defensible automated classification for ~30% of
variants (those with strong PVS1+PM2 or BA1 evidence). The remaining ~70% are
Uncertain Significance — not because the pipeline failed, but because the
mechanically computable subset of the ACMG criteria is insufficient for
classifying most variants. The 20 criteria that require functional studies,
pedigree segregation, or patient phenotype are the missing pieces; they cannot
be supplied by a VCF + databases alone.

Three concrete limits this pipeline has (and a clinical pipeline does not):

1. **No phasing.** Multi-variant haplotypes (e.g. compound heterozygosity for a
   recessive disorder) require phasing across the gene. Short-read data does
   not phase reliably; long-read sequencing or family trio data would.
2. **No patient phenotype.** The PP4 criterion ("phenotype highly specific
   for the disease") requires HPO terms from the patient's clinical record,
   which we do not have. The PVS1 + PM2 + PP3 combination might support a
   classification, but without phenotype confirmation it cannot make the
   transition to a clinical action.
3. **No clinician judgment.** Conflicting ClinVar interpretations (~3% of
   records, several in this dataset) require a human reviewer to weigh
   evidence — a 4-star Pathogenic from an expert panel vs a 1-star Likely
   Benign from a single submitter, both for the same variant, is a case the
   pipeline correctly flags but cannot resolve.

The disclaimer at the top of this report is not boilerplate. It marks the
exact boundary between what a pipeline can produce and what a clinical team
must produce. Reading the report, a clinical team has ~30% of the variants
classified for them and ~70% prepared for them to classify. That is the
correct division of labour.

## Reproducibility

```
conda env create -f env.yml -n c10-week08-miniproject
conda activate c10-week08-miniproject
bash run.sh
```

Output: `results/variants.html` and `results/variants.csv`. Run-info JSON at
`results/run-info.json` records database versions and run date.

Run-info excerpt:

```json
{
  "run_date": "<YOUR DATE>",
  "vep_version": "110.1",
  "vep_cache_version": "110",
  "snpeff_version": "5.2",
  "snpeff_db": "GRCh38.105",
  "clinvar_release": "2024-09-01",
  "gnomad_version": "v4.1.0",
  "assembly": "GRCh38"
}
```
```

---

## Common pitfalls

**The HTML report lacks the disclaimer.** This is the most important pitfall. The disclaimer is not optional. A report without it is not a teaching artifact; it is a problem waiting to be misinterpreted.

**The ACMG classifier produces PVS1 + BA1 on the same variant.** This is an evidence inconsistency. Either popmax > 0.05 (BA1) or popmax < 0.0001 (PM2 + thus possibly PVS1 + PM2) — both should not fire on the same row. Add a warning in `acmg_warnings` and inspect the input data.

**Variants are absent from gnomAD because the variant ID is malformed.** gnomAD's variant ID format is `chrom-pos-ref-alt` with no "chr" prefix. If your input VCF uses `chr17` and you query `chr17-43094077-G-A`, gnomAD returns nothing. Strip the prefix before constructing the query.

**The VEP cache is the wrong version.** `--cache_version 110` requires `${HOME}/.vep/homo_sapiens_vep_110_GRCh38/` on disk. Mismatched versions error out. Confirm with `ls ${HOME}/.vep/` before running.

**ClinVar tabix returns nothing.** The ClinVar VCF uses chromosome names with `chr` prefix; some other tools do not. If `tabix clinvar.vcf.gz 17:43094077-43094077` returns nothing, try `tabix clinvar.vcf.gz chr17:43094077-43094077`.

**The run-info JSON is missing fields.** A missing field is silently wrong; the report renders fine but cannot be reproduced. Add an `assert` that every required field is non-empty before writing the JSON.

**The SIFT/PolyPhen scores are not parsing.** VEP's `--sift b --polyphen b` emits strings like `deleterious(0.02)`. If you only pass `--sift p`, you get `deleterious` (no score), and the PP3 criterion cannot fire. Always use `b` for "both."

---

## Stretch goals

If you finish early and want to push further:

- **Add SpliceAI as a VEP plugin.** SpliceAI (Jaganathan et al. 2019, *Cell* 176:535) scores any variant for its likelihood of disrupting a splice site, including non-canonical ones. Add it via `--plugin SpliceAI,snv=spliceai_scores.raw.snv.hg38.vcf.gz` and add a `spliceai_score` column to the report. Use the score in BP7 (if synonymous AND SpliceAI < 0.1, no splice impact, BP7 fires) and PP3 (if SpliceAI > 0.5, splice impact, additional PP3 evidence).
- **Add CADD as a third predictor.** CADD (Rentzsch et al. 2019, *Nucleic Acids Res* 47:D886) combines ~60 features into a single deleteriousness score, scaled so >= 20 is "top 1% of likely-deleterious variants." Adds to PP3/BP4 in the 2024 refined framework.
- **Add a multi-sample mode.** Extend the pipeline to handle a multi-sample VCF (e.g. the GIAB NA12878 + NA24385 + NA24149 trio). Add a per-trio inheritance pattern column (de novo, autosomal recessive, autosomal dominant, X-linked).
- **Implement PS1 and PM5 by joining against ClinVar.** Both criteria require a ClinVar-derived table of known-pathogenic missense at specific amino acid positions. Build the table once at pipeline start; use it in `acmg_classifier.py`.
- **Build a per-gene summary page.** Group the per-variant rows by gene; for each gene with at least one variant, render a separate HTML section with the variants in that gene, their classifications, and a per-gene "gene-burden" summary (e.g. "this patient has 1 Likely Pathogenic, 2 VUS, and 3 Benign variants in BRCA1").
- **Reproduce a published variant interpretation.** Pick a variant from a published clinical case report (the OMIM database has thousands). Run your pipeline on a VCF containing that variant, compare the classification to the published one, write up the gap. The gap is where the clinical judgment lives.

---

## Submission

Push your `mini-project/` directory to the public portfolio repo with a commit message like:

```
Week 8 mini-project: end-to-end variant annotation, VEP+SnpEff+gnomAD+ClinVar+ACMG,
educational use only.
```

Open a PR to the curriculum repo with a brief description of what you produced and any deviations from the recipe. The graders will:

1. Clone the repo on a clean machine, run `bash run.sh`, and verify the pipeline reproduces the expected counts within ~10%.
2. Read the README. They expect every number to be cited (no "high" or "good" without a number); every tool to be pinned (VEP 110.1, not VEP); every disclaimer to be present.
3. Inspect `results/variants.html` in a browser and confirm the disclaimer is visible at the top, the metadata box is filled in, and the skipped-criteria section names ≥ 12 criteria.

---

## Up next

Week 9 leaves single-genome variant interpretation behind and moves to **phylogenetics and comparative genomics**: multiple-sequence alignment, tree building with RAxML and IQ-TREE, the dN/dS ratio. The cyvcf2 patterns from Week 8 are parked; the pandas reporting patterns come back. The disclaimer-on-every-output discipline you internalized this week applies to any future bioinformatics work that touches human data, in any week, in any year.
