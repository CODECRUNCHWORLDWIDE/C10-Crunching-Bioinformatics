# Week 11 — Resources

> **Educational and research use only.** None of the tools, papers, or knowledge bases referenced here are clinical software. Their outputs are valid as research artefacts and as teaching material. They are not substitutes for a clinical assay. Real patient-care decisions require an accredited laboratory operating under CAP / CLIA or equivalent oversight.

Use this file as the index of citations, papers, software pages, knowledge bases, and reference data. Open the **primary citation** for each tool when in doubt about a default; the GitHub README is a quick-start, the paper is the model.

---

## Primary papers (read these)

### Somatic variant calling

- **Mutect2 / Mutect.** Cibulskis K, Lawrence MS, Carter SL, et al. *Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples.* **Nature Biotechnology** 31:213-219 (2013). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/>. The original Mutect paper. Establishes the matched-tumor-normal Bayesian model: the prior on somatic mutation rate, the likelihood under the tumor allele-fraction model, the posterior probability that a candidate variant is somatic versus a germline / artifact. Mutect2 inside GATK is the modern rewrite; it adds joint indel calling, the panel-of-normals concept, and the active-region haplotype assembler. Read §2 (the model) carefully.

- **Strelka2.** Kim S, Scheffler K, Halpern AL, et al. *Strelka2: fast and accurate calling of germline and somatic variants.* **Nature Methods** 15:591-594 (2018). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>. Illumina's matched tumor-normal caller. Different statistical model from Mutect2 (the empirical variant scoring component is a Random Forest trained on truth sets); similar overall performance; commonly run alongside Mutect2 as a cross-check. The PASS-PASS overlap between Strelka2 and Mutect2 is typically 80-95% at AF >= 10%, lower at AF < 5%.

- **GATK4.** Van der Auwera GA, O'Connor BD. *Genomics in the Cloud: Using Docker, GATK, and WDL in Terra.* O'Reilly Media (2020). Not a paper but the canonical GATK4 reference. The GATK Best Practices for somatic variant calling are at <https://gatk.broadinstitute.org/hc/en-us/articles/360035531132>. The Mutect2 tool documentation is at <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2>.

### Population resources

- **gnomAD.** Karczewski KJ, Francioli LC, Tiao G, et al. *The mutational constraint spectrum quantified from variation in 141,456 humans.* **Nature** 581:434-443 (2020). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/>. The Genome Aggregation Database. The canonical source of population allele-frequency annotations used as the `--germline-resource` for Mutect2. Browser at <https://gnomad.broadinstitute.org/>. The VCFs are free downloads; we use the chr22 subset for the didactic dataset.

### Mutational signatures

- **COSMIC v3 signatures.** Alexandrov LB, Kim J, Haradhvala NJ, et al. *The repertoire of mutational signatures in human cancer.* **Nature** 578:94-101 (2020). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7054213/>. The Alexandrov lab's pan-cancer signature catalog derived from non-negative matrix factorization on the PCAWG cohort. 78 SBS signatures, plus DBS and ID signatures. The canonical reference; every modern signature tool maps to this catalog. Read §1 (the 96-class formulation) and §2 (the NMF model) for the mechanics.

- **deconstructSigs.** Rosenthal R, McGranahan N, Herrero J, et al. *deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution.* **Genome Biology** 17:31 (2016). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4762164/>. An R package for assigning a single tumor's spectrum to a reference catalog. Iterative least-squares with optional cosine-similarity threshold. Older than SigProfilerAssignment but still widely cited; commonly used as a sanity check.

- **SigProfiler.** SigProfilerAssignment, SigProfilerExtractor, SigProfilerMatrixGenerator, SigProfilerSimulator — the canonical Alexandrov-lab Python toolchain. Documentation at <https://github.com/AlexandrovLab>. The 2020 *Nature* paper describes the methodology; the tools themselves are free and open source under a BSD-3 license.

### Clinical-interpretation knowledge bases

- **COSMIC.** Sondka Z, Dhir NB, Carvalho-Silva D, et al. *COSMIC: a curated database of somatic variants and clinical data for cancer.* **Nucleic Acids Research** 52:D1210-D1217 (2024). Free at the journal: <https://academic.oup.com/nar/article/52/D1/D1210/7416441>. The Catalogue Of Somatic Mutations In Cancer at the Sanger Institute. Free for academic use; registration required. Browser at <https://cancer.sanger.ac.uk/cosmic>. The Census of Cancer Genes (the curated list of cancer driver genes) is at <https://cancer.sanger.ac.uk/census>.

- **OncoKB.** Chakravarty D, Gao J, Phillips SM, et al. *OncoKB: A Precision Oncology Knowledge Base.* **JCO Precision Oncology** 1:1-16 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5586540/>. The Memorial Sloan Kettering knowledge base of oncogenic alterations mapped to FDA evidence levels. The public-tier annotation API is at <https://www.oncokb.org/apiAccess>; the full set requires institutional registration but the public tier covers the standard panel of clinically actionable alterations. Evidence levels: 1 (FDA-approved in this tumor type), 2 (standard care), 3A / 3B (clinical evidence), 4 (biological evidence).

- **CIViC.** Griffith M, Spies NC, Krysiak K, et al. *CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer.* **Nature Genetics** 49:170-174 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5367263/>. The Washington University community-curated database. Free, open, no registration required. Browser at <https://civicdb.org/>. The TSV data dumps are at <https://civicdb.org/releases/main>; the API is at <https://docs.civicdb.org/en/latest/model/data_releases.html>.

### Variant annotation (referenced but not deeply used)

- **Funcotator.** GATK's functional annotation tool. Documentation at <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Funcotator>. Used to annotate a VCF with gene names, transcript IDs, amino-acid changes, dbSNP IDs, ClinVar significance, and COSMIC hits.

- **VEP (Ensembl Variant Effect Predictor).** McLaren W, Gil L, Hunt SE, et al. *The Ensembl Variant Effect Predictor.* **Genome Biology** 17:122 (2016). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4893825/>. The Ensembl annotation tool. Functionally similar to Funcotator; differs in default transcript sets and in plugin ecosystem.

- **SnpEff.** Cingolani P, Platts A, Wang LL, et al. *A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.* **Fly** 6:80-92 (2012). The third widely-used annotator; faster than VEP but with fewer database integrations.

### Background and method comparisons

- **PCAWG / Pan-Cancer Analysis of Whole Genomes.** The ICGC / TCGA Pan-Cancer Analysis of Whole Genomes Consortium. *Pan-cancer analysis of whole genomes.* **Nature** 578:82-93 (2020). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7025898/>. The umbrella paper for the cohort that the COSMIC v3 signatures were derived from. 2,658 whole-cancer-genome samples across 38 tumor types.

- **TCGA.** The Cancer Genome Atlas. <https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga>. The NCI's public cancer-genomics cohort. The MAF (Mutation Annotation Format) files distributed by TCGA are the canonical training set for cancer-bioinformatics pedagogy.

- **DREAM Somatic Mutation Calling Challenge.** Ewing AD, Houlahan KE, Hu Y, et al. *Combining tumor genome simulation with crowdsourcing to benchmark somatic single-nucleotide-variant detection.* **Nature Methods** 12:623-630 (2015). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509593/>. The community benchmark for somatic callers; the simulated tumor-normal pairs distributed by the challenge are the standard test set for new callers.

- **Mutect2 / Strelka2 head-to-head benchmarks.** Multiple. A useful one: Krøigård AB, Thomassen M, Lænkholm AV, et al. *Evaluation of Nine Somatic Variant Callers for Detection of Somatic Mutations in Exome and Targeted Deep Sequencing Data.* **PLoS ONE** 11:e0151664 (2016). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4801202/>. Pre-Mutect2 but methodologically representative; subsequent benchmarks (Xu 2018, Nat Sci Rep 8:13830) include Mutect2 and Strelka2 and agree that they are the top two callers on most input types.

---

## Software pages

| Tool | URL | Conda channel | Pinned version this week |
|------|-----|---------------|---------------------------|
| GATK4 (Mutect2 + FilterMutectCalls + CalculateContamination) | <https://gatk.broadinstitute.org/> | `bioconda::gatk4=4.5.0.0` | 4.5.0.0 |
| Strelka2 | <https://github.com/Illumina/strelka> | `bioconda::strelka=2.9.10` | 2.9.10 |
| samtools | <https://www.htslib.org/> | `bioconda::samtools=1.20` | 1.20 |
| bcftools | <https://www.htslib.org/> | `bioconda::bcftools=1.20` | 1.20 |
| pysam | <https://pysam.readthedocs.io/> | `bioconda::pysam=0.22.1` | 0.22.1 |
| SigProfilerAssignment | <https://github.com/AlexandrovLab/SigProfilerAssignment> | `bioconda::sigprofilerassignment=0.1.4` | 0.1.4 |
| SigProfilerMatrixGenerator | <https://github.com/AlexandrovLab/SigProfilerMatrixGenerator> | `bioconda::sigprofilermatrixgenerator=1.2.26` | 1.2.26 |
| deconstructSigs (R; via rpy2 if used) | <https://github.com/raerose01/deconstructSigs> | CRAN; install via R | 1.9.0 |
| Funcotator data sources (optional) | <https://gatk.broadinstitute.org/hc/en-us/articles/360035889931> | (bundle download) | funcotator_dataSources.v1.7.20200521s |
| Biopython | <https://biopython.org/> | `conda-forge::biopython=1.84` | 1.84 |
| pandas | <https://pandas.pydata.org/> | `conda-forge::pandas=2.2.2` | 2.2.2 |
| numpy | <https://numpy.org/> | `conda-forge::numpy=1.26.4` | 1.26.4 |
| scipy | <https://scipy.org/> | `conda-forge::scipy=1.13.1` | 1.13.1 |
| matplotlib (for plots) | <https://matplotlib.org/> | `conda-forge::matplotlib=3.8.4` | 3.8.4 |

### Canonical install command

```bash
conda create -n cancer-w11 -c bioconda -c conda-forge \
  gatk4=4.5.0.0 strelka=2.9.10 samtools=1.20 bcftools=1.20 \
  pysam=0.22.1 sigprofilerassignment=0.1.4 sigprofilermatrixgenerator=1.2.26 \
  biopython=1.84 pandas=2.2.2 numpy=1.26.4 scipy=1.13.1 matplotlib=3.8.4 \
  python=3.11
conda activate cancer-w11
```

The SigProfiler reference genome data (GRCh38) is downloaded on first use by `SigProfilerMatrixGenerator install GRCh38`; this is a ~3 GB one-time download.

---

## Reference data

The Week 11 didactic dataset ships chr22-only subsets of the standard references in `data/` (see the mini-project README). The full datasets are downloadable from the Broad Institute resource bundle at <https://console.cloud.google.com/storage/browser/genomics-public-data/references>.

| Resource | Full source | Didactic subset shipped |
|----------|-------------|--------------------------|
| GRCh38 / hg38 reference FASTA | `Homo_sapiens_assembly38.fasta` (Broad) | `data/chr22_GRCh38.fasta` (~50 MB) |
| GRCh38 reference index `.fai` | (Broad) | `data/chr22_GRCh38.fasta.fai` |
| GRCh38 reference dictionary `.dict` | (Broad) | `data/chr22_GRCh38.dict` |
| Panel of Normals (PON) | `1000g_pon.hg38.vcf.gz` (Broad somatic-hg38 bundle) | `data/chr22_pon.vcf.gz` (~5 MB) |
| Germline resource (gnomAD) | `af-only-gnomad.hg38.vcf.gz` (Broad somatic-hg38 bundle) | `data/chr22_gnomad.vcf.gz` (~30 MB) |
| Tumor BAM (didactic) | (simulated; see mini-project README) | `data/tumor_chr22.bam` (~80 MB) |
| Normal BAM (didactic) | (simulated; see mini-project README) | `data/normal_chr22.bam` (~80 MB) |
| Tumor BAM index | (samtools index) | `data/tumor_chr22.bam.bai` |
| Normal BAM index | (samtools index) | `data/normal_chr22.bam.bai` |

The full Broad somatic-hg38 bundle is documented at <https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle>.

### SigProfiler reference genome

`SigProfilerMatrixGenerator install GRCh38` writes ~3 GB of trinucleotide-context lookup tables under `~/.SigProfilerMatrixGenerator/`. This is a one-time download; subsequent runs use the cache. The same step is required for `GRCh37` if you work with hg19 data.

### COSMIC mutational-signature catalog

The v3.3 catalog is bundled with SigProfilerAssignment 0.1.4 (`COSMIC_v3.3.1_SBS_GRCh38.txt` and similar for ID and DBS). The standalone download is at <https://cancer.sanger.ac.uk/signatures/sbs/> (free; no registration). The 96-class trinucleotide ordering follows the Alexandrov-lab convention: substitutions are normalized to the pyrimidine reference (C>A, C>G, C>T, T>A, T>C, T>G) and within each substitution the 16 trinucleotide flanks are ordered ACA, ACC, ACG, ACT, CCA, CCC, ..., TCT (the central base is the substituted base; the flanking bases iterate A < C < G < T).

---

## Knowledge-base look-up

### COSMIC

- Browser: <https://cancer.sanger.ac.uk/cosmic>
- Cancer Gene Census: <https://cancer.sanger.ac.uk/census>
- API documentation: <https://cancer.sanger.ac.uk/cosmic/help/api>
- Academic use: free with registration; commercial use requires a licence.

To look up a variant: <https://cancer.sanger.ac.uk/cosmic/search?q=TP53+R175H>. The page reports the cancer-type frequency, the most-cited papers, and the linked tissue distribution.

### OncoKB

- Browser: <https://www.oncokb.org/>
- API: <https://api.oncokb.org/oncokb-website/api>
- Public-tier annotation: <https://www.oncokb.org/apiAccess> (rate-limited; sign in for a free token).
- Evidence levels: <https://www.oncokb.org/levels>

Example variant lookup: <https://www.oncokb.org/gene/TP53/R175H>. The page reports the OncoKB Mutation Effect ("Loss-of-function"), the FDA evidence level ("Level 3B"), and the linked clinical trials.

### CIViC

- Browser: <https://civicdb.org/>
- TSV release: <https://civicdb.org/releases/main>
- API: <https://docs.civicdb.org/en/latest/model/data_releases.html>
- Documentation: <https://docs.civicdb.org/>

Example variant lookup: <https://civicdb.org/variants/12/summary>. The page reports the evidence items, the supporting publications, and the assertion-level summary.

---

## Datasets for practice

### Didactic (shipped with this week)

- `data/tumor_chr22.bam` + `data/normal_chr22.bam`: a small simulated chr22 tumor-normal pair derived from a public GRCh38 reference plus synthetic somatic variants seeded at known positions. ~80 MB each. The known-truth VCF is at `data/truth_chr22.vcf.gz`; you can use it to compute precision and recall after FilterMutectCalls.

### Public (download separately)

- **DREAM Somatic Mutation Calling Challenge** simulated pairs: ICGC-DREAM Mutation Calling Challenge synthetic datasets, free at <https://www.synapse.org/#!Synapse:syn312572>. ~10 GB each pair.
- **GIAB (Genome in a Bottle)** germline benchmark: HG001 / HG002 / HG003 etc. <https://www.nist.gov/programs-projects/genome-bottle>. Useful for the germline side of the equation; the GIAB family also has a Mendelian-trio version helpful for germline filter tuning.
- **TCGA MAF dumps**: <https://portal.gdc.cancer.gov/>. The TCGA MAFs are the canonical reference for "what does a real tumor look like at the variant level?" Use the PCAWG-derived Pan-Cancer file <https://gdc.cancer.gov/about-data/publications/pancanatlas> if you want the harmonized version.

### Mutational-signature reference profiles

- **COSMIC v3.3 SBS / DBS / ID profiles**: <https://cancer.sanger.ac.uk/signatures/sbs/>. Free, no registration.
- **PCAWG signature attributions per sample**: <https://www.synapse.org/#!Synapse:syn11801791>. Free with Synapse registration.

---

## Helpful tutorials and walkthroughs

- **GATK Best Practices somatic SNV/indel pipeline.** <https://gatk.broadinstitute.org/hc/en-us/articles/360035531132>. The canonical end-to-end walkthrough.
- **Mutect2 Tumor-Normal somatic short variant discovery** (Broad article). <https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-Mutect2-Somatic-Short-Variant-Discovery-SNVs-Indels-->. The step-by-step subset that we are translating to Python.
- **SigProfilerAssignment tutorial.** <https://github.com/AlexandrovLab/SigProfilerAssignment#tutorials>. The example notebook covers the canonical inputs and outputs.
- **The Galaxy Training Network "Somatic variant calling" hands-on.** <https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/somatic-variants/tutorial.html>. A GUI-based walkthrough using public test data; the conceptual steps are identical to the CLI / Python pipeline this week.

---

## Citation snippets (drop-in for your write-ups)

- Mutect2 / Mutect: "We called somatic SNVs and indels with Mutect2 inside GATK 4.5.0.0 in tumor-normal mode (Cibulskis et al. 2013, *Nature Biotechnology* 31:213) against GRCh38, using the Broad-published 1000G Panel of Normals (`1000g_pon.hg38.vcf.gz`) and the gnomAD v2 allele-frequency annotation as the germline resource (`af-only-gnomad.hg38.vcf.gz`)."
- FilterMutectCalls: "Raw variants were filtered with FilterMutectCalls using the default thresholds, with `--contamination-table` set to the per-sample contamination estimate from GATK CalculateContamination on pileups derived from common biallelic GRCh38 sites."
- Strelka2: "We cross-checked the Mutect2 PASS set against an independent Strelka2 2.9.10 run on the same BAM pair (Kim et al. 2018, *Nature Methods* 15:591); the two callers agreed on 87.6% of PASS SNVs at allele frequency >= 10%."
- COSMIC signatures: "We decomposed the somatic SNV trinucleotide spectrum against the COSMIC v3.3 SBS reference (Alexandrov et al. 2020, *Nature* 578:94) using SigProfilerAssignment 0.1.4 (Alexandrov-lab); the top contributing signatures were SBS1, SBS5, and SBS18 with a reconstructed-spectrum cosine similarity of 0.93."
- OncoKB: "Clinically actionable somatic variants were annotated against OncoKB (Chakravarty et al. 2017, *JCO Precision Oncology* 1:1) using the public-tier annotation API."
- CIViC: "Variant-level clinical evidence was cross-referenced against the CIViC community knowledge base (Griffith et al. 2017, *Nature Genetics* 49:170; data release 2025-03)."

---

## A note on what is NOT in this resource list

We do **not** link to any tool that gates its core functionality behind a commercial licence (Strelka1's old single-sample license is irrelevant — Strelka2 is open source; we are not using SomaticSeq's commercial bundle or the proprietary callers). We do not link to any clinical-grade pipeline (Foundation Medicine, Guardant, Caris, Tempus) because those are subscription clinical-laboratory services that you cannot reproduce from a paper alone, and because this week's deliverables are not clinical. Every tool we reference is free for academic and research use; every paper we cite is free at PMC or at the publisher's open-access link. If you need a tool that is not in the table above and is not free, mention it in your write-up and explain why; do not silently swap it in.

---

## Quick-start commands (copy-paste reference)

### Set up the conda environment

```bash
conda create -n cancer-w11 -c bioconda -c conda-forge \
  gatk4=4.5.0.0 strelka=2.9.10 samtools=1.20 bcftools=1.20 \
  pysam=0.22.1 sigprofilerassignment=0.1.4 sigprofilermatrixgenerator=1.2.26 \
  biopython=1.84 pandas=2.2.2 numpy=1.26.4 scipy=1.13.1 matplotlib=3.8.4 \
  python=3.11
conda activate cancer-w11
python -c "from SigProfilerMatrixGenerator import install; install.install('GRCh38')"
```

### Index the reference and BAMs

```bash
samtools faidx data/chr22_GRCh38.fasta
gatk CreateSequenceDictionary -R data/chr22_GRCh38.fasta
samtools index data/tumor_chr22.bam
samtools index data/normal_chr22.bam
```

### Verify BAM @RG SM: headers

```bash
samtools view -H data/tumor_chr22.bam  | grep '^@RG'
samtools view -H data/normal_chr22.bam | grep '^@RG'
```

### Mutect2 tumor-normal call

```bash
gatk --java-options "-Xmx4g" Mutect2 \
  -R data/chr22_GRCh38.fasta \
  -I data/tumor_chr22.bam \
  -I data/normal_chr22.bam \
  -tumor TUMOR_SAMPLE \
  -normal NORMAL_SAMPLE \
  --panel-of-normals data/chr22_pon.vcf.gz \
  --germline-resource data/chr22_gnomad.vcf.gz \
  -L data/chr22_intervals.bed \
  -O results/unfiltered.vcf.gz \
  --native-pair-hmm-threads 4
```

### GetPileupSummaries + CalculateContamination

```bash
gatk --java-options "-Xmx4g" GetPileupSummaries \
  -I data/tumor_chr22.bam \
  -V data/chr22_common_biallelic.vcf.gz \
  -L data/chr22_common_biallelic.vcf.gz \
  -O results/tumor.pileups.table

gatk --java-options "-Xmx4g" GetPileupSummaries \
  -I data/normal_chr22.bam \
  -V data/chr22_common_biallelic.vcf.gz \
  -L data/chr22_common_biallelic.vcf.gz \
  -O results/normal.pileups.table

gatk --java-options "-Xmx4g" CalculateContamination \
  -I results/tumor.pileups.table \
  -matched results/normal.pileups.table \
  -O results/contamination.table
```

### FilterMutectCalls

```bash
gatk --java-options "-Xmx4g" FilterMutectCalls \
  -R data/chr22_GRCh38.fasta \
  -V results/unfiltered.vcf.gz \
  --contamination-table results/contamination.table \
  --stats results/unfiltered.vcf.gz.stats \
  -O results/filtered.vcf.gz
```

### PASS-only subset

```bash
bcftools view -f PASS results/filtered.vcf.gz -O z -o results/pass.vcf.gz
tabix -p vcf results/pass.vcf.gz
```

### 96-class spectrum + SigProfilerAssignment (Python)

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("../exercises").resolve()))
# Import helpers from exercise-03 (rename or import the file as a module).
```

### Strelka2 (optional cross-check)

```bash
configureStrelkaSomaticWorkflow.py \
  --normalBam data/normal_chr22.bam \
  --tumorBam data/tumor_chr22.bam \
  --referenceFasta data/chr22_GRCh38.fasta \
  --runDir results/strelka_run \
  --callRegions data/chr22_intervals.bed.gz

results/strelka_run/runWorkflow.py -m local -j 4
```

### Normalize VCFs for cross-caller comparison

```bash
bcftools norm -f data/chr22_GRCh38.fasta -m -any results/filtered.vcf.gz \
  | bcftools sort -O z -o results/mutect2_norm.vcf.gz
tabix -p vcf results/mutect2_norm.vcf.gz
```

---

## OncoTree codes (selected)

OncoKB queries are tumor-type-specific via OncoTree codes. The most-common codes:

| Code  | Tumor type |
|-------|------------|
| BRCA  | Breast Cancer |
| LUAD  | Lung Adenocarcinoma |
| LUSC  | Lung Squamous Cell Carcinoma |
| COAD  | Colon Adenocarcinoma |
| READ  | Rectal Adenocarcinoma |
| PRAD  | Prostate Adenocarcinoma |
| OV    | Ovarian Cancer |
| GBM   | Glioblastoma Multiforme |
| LGG   | Low-Grade Glioma |
| SKCM  | Skin Cutaneous Melanoma |
| AML   | Acute Myeloid Leukemia |
| ALL   | Acute Lymphoblastic Leukemia |
| CML   | Chronic Myeloid Leukemia |
| CLL   | Chronic Lymphocytic Leukemia |
| MM    | Multiple Myeloma |
| PAAD  | Pancreatic Adenocarcinoma |
| HCC   | Hepatocellular Carcinoma |
| STAD  | Stomach Adenocarcinoma |
| ESCA  | Esophageal Carcinoma |
| BLCA  | Bladder Urothelial Carcinoma |
| KIRC  | Kidney Renal Clear Cell Carcinoma |
| THCA  | Thyroid Carcinoma |
| HNSC  | Head and Neck Squamous Cell Carcinoma |

The full OncoTree is at <http://oncotree.mskcc.org/>. The code matters: a variant might be Level 1 in BRCA and Level 3B in LUAD.

---

## COSMIC mutational-signature aetiology cheat sheet

Selected SBS signatures from COSMIC v3.3 (full catalog at <https://cancer.sanger.ac.uk/signatures/sbs/>):

| Signature | Aetiology | Tumor-type association | Notes |
|-----------|-----------|------------------------|-------|
| SBS1 | 5-methylcytosine deamination | All, age-correlated | Clock-like; correlates with age at diagnosis |
| SBS2 | APOBEC cytidine deaminase | Bladder, cervix, breast | TpC trinucleotide context |
| SBS3 | Homologous-recombination deficiency | Breast, ovarian, pancreatic | BRCA1/2 carriers; PARPi-sensitivity hint |
| SBS4 | Tobacco smoking | Lung, head & neck | CpC > TpA bias |
| SBS5 | Clock-like; unknown | All | Often co-dominates with SBS1 |
| SBS6 | MMR deficiency | Colorectal, endometrial | MSI tumors |
| SBS7a/b/c/d | UV light | Skin (melanoma, SCC) | C>T at TpC |
| SBS9 | POLH activity | Lymphoid; SHM-associated | Hypermutation at lymphoid genes |
| SBS10a/b | POLE exonuclease deficiency | Hypermutated tumors | TCT > TAT or TCG > TTG signatures |
| SBS11 | Temozolomide | Treated tumors | C>T at NpC |
| SBS13 | APOBEC | Bladder, cervix, breast | TpC > GpC bias |
| SBS15 | MMR deficiency | Colorectal | Co-occurs with MSI |
| SBS18 | Reactive oxygen damage | Various | G>T bias |
| SBS20 | MMR + POLD1 mutation | Colorectal | Specific MMR variant |
| SBS22 | Aristolochic acid | Liver, urothelial | A>T at TpA |
| SBS24 | Aflatoxin | Liver | C>A bias |
| SBS26 | MMR deficiency | Colorectal, endometrial | MSI tumors |
| SBS29 | Tobacco chewing | Oral | Different from smoking |
| SBS31 | Platinum chemotherapy | Treated | C>T at CpC |
| SBS35 | Platinum chemotherapy | Treated | Similar to SBS31 |
| SBS36 | MUTYH BER deficiency | Colorectal | Adenomatous polyposis |
| SBS39 | Unknown | (Degenerate with SBS3) | Report both as ambiguous |
| SBS40 | Clock-like (kidney) | Kidney; aging | Often co-linear with SBS5 |

The aetiology column is well-supported for many signatures, suggestive for others, and unknown for some. Read the COSMIC documentation for the citations.

---

## Stretch reading

- **Clonal evolution.** Nowell PC. *The clonal evolution of tumor cell populations.* **Science** 194:23-28 (1976). The original clonal-evolution framework that still organizes how we think about why somatic variants partition into drivers and passengers. Closed-access at the journal but the concept is well-summarized in any modern cancer-biology textbook.
- **Driver and passenger mutations.** Stratton MR, Campbell PJ, Futreal PA. *The cancer genome.* **Nature** 458:719-724 (2009). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2821689/>. The canonical "what is the cancer genome?" review; defines drivers vs passengers, the classes of cancer genes (oncogenes, tumor suppressors, stability genes), and the genomic landscape concept.
- **Tumor heterogeneity.** McGranahan N, Swanton C. *Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future.* **Cell** 168:613-628 (2017). Free at the journal: <https://www.cell.com/cell/fulltext/S0092-8674(17)30100-7>. The modern review of intratumor heterogeneity, sub-clonal variants, and why a single biopsy is an incomplete sampling.
- **TP53.** Kastenhuber ER, Lowe SW. *Putting p53 in Context.* **Cell** 170:1062-1078 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5743327/>. The canonical TP53 review; useful for grounding "TP53 R175H" into the broader picture of the most-mutated gene in cancer.

---

## A final note on databases that require institutional access

OncoKB and COSMIC have public tiers and licensed tiers. The public tiers cover the standard panel of variants that this week needs; the licensed tiers add up-to-the-week curated updates and the full annotation API. If you are at a university or non-profit, you can usually obtain the academic licence at no cost via your institution's grants office. CIViC is fully free and open. None of the public tiers gates the variants needed for this week.

The Week 11 deliverables are reproducible from free resources only. If you find a tool we point to has become non-free, please open an issue.
