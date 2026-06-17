# Lecture 1 — Somatic vs Germline: The Matched-Pair Model, Purity, Contamination, and CNV at the Concept Level

> **Educational and research use only.** The model in this lecture is the same model used by accredited clinical laboratories, but the implementation we build does not have the validation, the version locking, the quality-management system, or the regulatory oversight that turns a research pipeline into a clinical assay. Treat every output of this week as a learning artefact. Do not apply it to patient care.

## 1. The question and the conceptual model

The central question of cancer genomics: given DNA sequencing data from a patient with cancer, **which of the variants observed in the tumor are causally relevant and which are background noise?** Two kinds of answer are wrong before we even start. Treating every variant in the tumor as a candidate causal variant ignores the fact that the patient was born with ~4-5 million germline variants (mostly common, mostly inherited) that have nothing to do with their cancer. Treating only the rare, novel variants as candidates filters out the recurrent driver mutations (TP53 R175H, KRAS G12D, BRAF V600E) that appear in thousands of tumors and are present in dbSNP at non-trivial frequency precisely because they are recurrent in disease.

The standard resolution: **sequence two samples from the same patient — a tumor sample and a matched normal sample — and call as "somatic" the variants present in the tumor that are not also present in the normal**. The matched normal is the patient's germline baseline; any variant in the tumor that is not in the normal must have arisen in the tumor's clonal lineage. This is the **matched tumor-normal model** and it is the foundation of every modern somatic variant-calling pipeline.

The model assumes:

- The normal sample is genuinely *not tumor*. In practice this is blood (most common), buccal swab, or histologically-confirmed normal tissue adjacent to the tumor.
- The tumor sample contains tumor cells. Most clinical specimens are mixtures of tumor cells, immune infiltrate, stroma, blood vessels, and necrotic regions; the fraction of tumor cells is the **purity**.
- The two samples come from the same patient and have not been swapped or cross-contaminated in the lab.
- The two samples have been processed equivalently: same library prep chemistry, same sequencing run if possible, same alignment, same duplicate-marking, same base-quality recalibration. Differences in any of those introduce systematic biases between the samples that the caller will interpret as "tumor-specific" signal.

When any of those assumptions fails, the caller produces defensible-looking results that are wrong. The rest of this lecture is the catalog of how the assumptions fail and how to detect each failure.

## 2. Clonal evolution and what "somatic" actually means

Cancer is a **clonal disease**. The conventional model (Nowell 1976, *Science* 194:23; the original framing) goes like this: in a normal tissue, a single cell acquires a driver mutation (TP53 inactivation, an APC truncation, a KRAS gain-of-function) that confers a fitness advantage in that tissue micro-environment. That cell divides; its descendants share the founding driver and accumulate additional mutations as they divide. The descendants compete with each other; sub-clones with additional drivers expand at the expense of clones without those additional drivers; over months to years the population evolves into a clinically detectable tumor with a polyclonal structure. A typical solid tumor at biopsy contains:

- **Truncal mutations.** Variants present in every tumor cell because they were present in the founding clone. These are the earliest drivers and the easiest to detect because they appear at allele frequency ~50% (heterozygous) or ~100% (homozygous after loss-of-heterozygosity) in every cell of the tumor.
- **Sub-clonal mutations.** Variants present in only a fraction of the tumor cells because they arose in a sub-clone after the founding event. These appear at lower allele frequencies — proportional to the sub-clone's fraction of the tumor cell population.
- **Passenger mutations.** Variants that arose at the same time as drivers but had no fitness consequence. The vast majority of mutations in any tumor are passengers; the few drivers do the causal work but the passengers vastly outnumber them.

The somatic-variant caller's job is to find all of the above and report each with an allele frequency. Distinguishing drivers from passengers is the downstream interpretation layer (OncoKB, CIViC, COSMIC Cancer Gene Census). The caller does not do that; it just produces a list with frequencies and the interpretation tools annotate the list.

The matched-normal subtraction is what removes the **germline** variants — the variants the patient was born with — from the candidate list. Germline variants are by definition present in every cell of the body, including the normal sample, at approximately the same allele frequency (~50% heterozygous, ~100% homozygous; modulo small fluctuations from sequencing depth). When the caller sees a variant at 50% AF in both tumor and normal, it concludes "germline" and removes it from the somatic candidate set.

## 3. The matched normal — what to sequence and what can go wrong

The matched normal is conventionally **peripheral blood**: drawn at biopsy or at diagnosis, processed into a buffy-coat sample, DNA-extracted, library-prepped alongside the tumor, sequenced on the same instrument. Blood is the standard because it is non-invasive to collect, because the DNA is uncontaminated by adjacent tumor tissue, and because the cell turnover is rapid enough that the somatic mutation burden of the blood compartment is much lower than the tumor's.

Three exceptions where blood is not the right normal:

- **Hematologic malignancies.** A leukemia patient's blood is, by definition, leukemic. The "normal" blood sample contains the tumor; calling against it will subtract the leukemic variants out as if they were germline. The fix is to use a non-hematopoietic tissue: a buccal swab, a skin biopsy, or hair follicles. Each has its own limitations (buccal swabs are contaminated with oral microbiome and food DNA; skin biopsies have UV-induced field effects; hair follicles have low yield).
- **Field-effect cancers.** Smokers' lung tissue and Barrett's-esophagus patients' esophageal tissue carry pre-cancerous mutations across the entire affected tissue. A "normal" adjacent-tissue sample will share field-effect variants with the tumor and the caller will subtract them out. Blood is the better normal in these cases.
- **Constitutional mosaicism.** Rare patients have somatic variants that arose during embryonic development and are present in some tissues but not others. A normal sample from one tissue may miss variants present in another tissue. Detection is rare and difficult; usually requires WGS of multiple tissues.

The standard workflow assumes blood is fine. For solid tumors that have not metastasized into bone marrow and for non-hematologic malignancies, that is correct.

## 4. Purity: the fraction of the tumor sample that is tumor cells

A clinical tumor specimen is never 100% tumor cells. A typical sample (a punch biopsy or a resection wedge sent for sequencing) is a mixture:

- Tumor cells (the clone of interest).
- Stromal cells (fibroblasts, endothelial cells, smooth muscle, adipocytes — host tissue that the tumor is growing within).
- Immune infiltrate (tumor-infiltrating lymphocytes, macrophages, neutrophils — host immune cells responding to the tumor).
- Blood vessels (red and white blood cells in the vasculature).
- Necrotic tumor (dead tumor cells that have lost membrane integrity; their DNA is often degraded and contributes broken short reads).

The **purity** of the sample is the fraction of cells that are tumor cells. A common pathology-estimated value for a solid tumor is 40-70%. Some indolent tumors (low-grade prostate cancers, low-grade hematopoietic malignancies infiltrating a sampled tissue) can be much lower; some aggressive tumors (high-grade osteosarcoma, treatment-naive triple-negative breast cancers) can exceed 90%.

Purity matters for the allele frequencies the caller sees. A truncal heterozygous somatic variant in a pure tumor has allele frequency 50% (one mutant + one wild-type chromosome, half the reads cover the mutant). In a 60%-pure sample, the apparent AF is 60% * 50% + 40% * 0% = 30% (the 60% tumor cells contribute their 50% AF; the 40% non-tumor cells contribute 0% AF since they do not carry the variant). At 30% purity, the same variant appears at 15%; at 10% purity, at 5%; at 5% purity, at 2.5% — which is at or below the noise floor of most callers' default thresholds.

The relationship between purity, the true allele fraction in the tumor cells, and the observed allele fraction in the mixture is:

```text
observed_AF = purity * tumor_cell_AF + (1 - purity) * normal_AF
```

For a tumor cell variant not present in the normal (`normal_AF = 0`), this simplifies to `observed_AF = purity * tumor_cell_AF`. For a heterozygous truncal variant (`tumor_cell_AF = 0.5`), the relationship is `observed_AF = 0.5 * purity`. The caller does not know the purity; it infers it from the distribution of observed AFs across many variants.

**Purity estimation tools.** Several free tools estimate purity from sequencing data: FACETS, Sequenza, PureCN, ABSOLUTE. They typically use both the allele-frequency distribution and the copy-number signal jointly; we describe CNV in §6. The Week 11 pipeline does not include a dedicated purity estimator (each has its own setup cost); the mini-project write-up names the limitation.

## 5. Contamination: a different patient's DNA in the sample

Contamination is *not* purity. Purity is "what fraction of the cells in this sample are tumor cells from this patient?" Contamination is "what fraction of the DNA in this sample is from a *different* individual?"

Cross-sample contamination in a sequencing facility is unavoidable at low levels. A 1% contamination rate is realistic in a high-throughput lab; a 0.1% rate is excellent; a 5% rate is a serious quality flag. The contaminant is usually another patient's sample that was processed nearby (the index sequence partially overlapped; a robotic pipette aerosolized between wells; the sequencer's flowcell carried over from a previous run).

A 1% contamination rate from a high-AF germline variant in the contaminant manifests as a 1% AF "somatic" call in the patient. From the caller's perspective:

- The variant is absent from the matched normal (which is also at ~1% contamination, but with a *different* contaminant typically, or with too little signal to detect).
- The variant is present at 1% in the tumor.
- The caller flags it as somatic at 1% AF.

This is a false positive. In a 200 Mb genome at 50x coverage, a 1% contamination rate produces tens of thousands of these false positives — enough to drown out the genuine low-AF somatic signal.

**Contamination estimation.** GATK's `GetPileupSummaries` walks a list of known biallelic germline sites (common SNPs from a public source) and counts the allele depth at each. `CalculateContamination` then takes these pileups and computes the maximum-likelihood contamination fraction under the assumption that contamination is by an unrelated individual sampled from the same population. The output is a small `contamination.table` file with one row per sample and one column for the estimated contamination fraction.

The standard call:

```bash
gatk GetPileupSummaries \
  -I tumor.bam \
  -V common_biallelic.vcf.gz \
  -L common_biallelic.vcf.gz \
  -O tumor.pileups.table

gatk GetPileupSummaries \
  -I normal.bam \
  -V common_biallelic.vcf.gz \
  -L common_biallelic.vcf.gz \
  -O normal.pileups.table

gatk CalculateContamination \
  -I tumor.pileups.table \
  -matched normal.pileups.table \
  -O contamination.table
```

The result is passed to `FilterMutectCalls --contamination-table contamination.table`. Variants whose AF is consistent with the estimated contamination fraction get the `contamination` filter reason.

A clean lab will produce contamination estimates in the 0-2% range; a 5%+ estimate is a flag for the lab. The Week 11 didactic dataset is simulated and noise-free, so the contamination is 0; on real data, expect a small positive value.

## 6. Copy-number variation (CNV) at the conceptual level

Most somatic variant calling is about *single-nucleotide variants and small indels*: a substitution, an insertion of one base, a deletion of a few bases. Cancer also produces large-scale genomic changes: deletions of millions of bases, duplications, focal amplifications, whole-chromosome gains and losses, complex rearrangements (chromothripsis), and structural variants between chromosomes (translocations). These are the **copy-number variations** and **structural variants** of cancer.

We do not call CNVs in the Week 11 pipeline. Each CNV caller has substantial setup cost (a calibrated normal-coverage profile, an interval list, a segmentation model), and CNV interpretation is its own week-long topic. But we cannot teach somatic SNV interpretation without naming CNVs because **CNVs change the interpretation of SNVs**.

Two examples:

- **Loss of heterozygosity (LOH) at TP53.** A patient inherits two TP53 alleles. A somatic mutation in one allele produces a heterozygous TP53 mutation; the wild-type allele is still present and is still producing functional p53 protein. If the tumor *also* loses the wild-type allele (through a focal deletion that removes the region containing the wild-type copy), the tumor has *biallelic* TP53 loss — no functional p53. This is the disease-causing event. The SNV alone does not tell you about the LOH; the CNV call is what reveals it. A research-grade SNV pipeline that does not also call CNVs will report the TP53 mutation but miss that it is biallelic, and will likely under-estimate its causal significance.
- **Focal amplification of MYC or ERBB2.** Some tumors amplify a small region containing a single oncogene (HER2 in breast cancer; MYC in lymphoma). The oncogene is wild-type at the sequence level — no SNVs — but is present at 10x to 50x normal copy number. The SNV caller sees nothing; the CNV caller sees the amplification. Without the CNV call, the tumor's primary driver is invisible.

A research-grade SNV pipeline should be honest about this: the somatic SNV report covers only the SNV layer of cancer genomics; CNVs, structural variants, and gene fusions are separate layers that require separate tools.

**The free CNV tools.** GATK CNV (gatk CallCopyRatioSegments and friends), FACETS, Sequenza, CNVkit, PureCN. All free and open source. Most use the log2 read-depth ratio of tumor to normal, often co-modeled with the B-allele frequency of common heterozygous SNPs. We name them so you know the landscape; we do not call them in the Week 11 mini-project.

## 7. The data flow of a somatic SNV pipeline (preview)

The full flow of the Week 11 mini-project, from raw input to final report:

```text
                  ┌──────────────────┐
                  │ Tumor BAM        │
                  │ Normal BAM       │
                  │ Reference FASTA  │
                  │ PON VCF          │
                  │ gnomAD VCF       │
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ Mutect2          │   tumor-normal mode
                  │ (GATK 4.5.0.0)   │   produces unfiltered VCF
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ GetPileupSummaries│  one per sample
                  │ + CalculateContam│   produces contamination.table
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ FilterMutectCalls│   applies standard filter set
                  │                  │   produces filtered VCF
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ Parse PASS VCF   │   pysam / cyvcf2
                  │ Compute 96-class │   trinucleotide context
                  │ spectrum         │
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ SigProfiler      │   COSMIC v3.3 catalog
                  │ Assignment       │   per-signature contribution
                  │                  │   reconstructed-spectrum cosine
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ Optional OncoKB /│   per-variant evidence levels
                  │ CIViC look-up    │   no auto-calls
                  └─────────┬────────┘
                            │
                            ▼
                  ┌──────────────────┐
                  │ Markdown report  │   methods + results + limits
                  │ run-info.json    │   versions + parameters
                  └──────────────────┘
```

Each arrow is a tool call; each tool has a version, an input mode, and a set of parameters that must be pinned in the `run-info.json`. The mini-project's job is to make this flow reproducible from a single `bash run.sh`.

## 8. Why the matched normal makes the call sensitive

Without a matched normal, a somatic caller has to distinguish "rare tumor variant" from "rare germline variant" from "sequencing error" purely from the AF distribution in the tumor. The error floor of Illumina sequencing is around 0.1-0.5% per base; rare germline variants in a tumor BAM look exactly like 50%-AF heterozygous sites; ultra-rare germline variants (population frequency < 0.001) look like rare somatic variants. The tumor-only caller has to discriminate with no help except a public allele-frequency database (gnomAD); the discrimination is inherently noisier.

The matched normal collapses two of these confusions. If the variant is present in the normal at any non-trivial AF, it is germline (it was in the patient's genome before the cancer). If it is absent from the normal but present in the tumor, it is somatic (the tumor's clonal lineage acquired it). The remaining ambiguity is between "true low-AF somatic" and "sequencing error in the tumor + sequencing error in the normal happens to leave the variant only in tumor"; the caller addresses that with the joint likelihood model.

Mutect2's matched-normal model (Cibulskis et al. 2013, *Nat Biotechnol* 31:213) computes the posterior probability that a candidate variant is somatic versus a germline or artefact, given the read pileups at the variant position in both tumor and normal, the population allele frequency from the germline resource, and the panel-of-normals annotation. The math is in §2 of the paper. The intuition: the caller's prior on "this position has a somatic variant" is small (~1e-6 per base per generation in normal tissue, somewhat higher in tumor); the prior is updated by the tumor and normal read counts. A variant at 30% AF in tumor and 0% AF in normal at 50x coverage will receive an extreme posterior in favour of "somatic"; a variant at 50% AF in both will receive an extreme posterior in favour of "germline".

## 9. The panel of normals (PON)

The matched normal handles the germline confusion. The panel of normals handles the **technical-artefact** confusion. Some variants appear as candidate somatic calls across many samples not because they are recurrent disease drivers but because the alignment is systematically wrong at that position (mapping into a low-complexity region; near a known repeat boundary; in a region of structural copy-number-variation in the reference). These recurrent technical artefacts are easy to identify in a large cohort: a variant that appears as "somatic" in 5% of normal samples in your sequencing core's archive is not somatic; it is a recurrent artefact of the library prep or the alignment.

A **panel of normals** is a multi-sample VCF distilled from many technically-similar normal samples. The standard Broad Institute PON for GRCh38 (`1000g_pon.hg38.vcf.gz`) is derived from a thousand normal samples sequenced at the Broad; variants present in more than two of them are recorded. Mutect2 takes the PON via `--panel-of-normals` and down-weights or filters variants that match a PON entry.

A PON is **platform-specific**. A PON built from HiSeq X reads will not catch NovaSeq-specific artefacts; an exome PON will not handle the whole-genome artefacts; a PON from one library prep chemistry will not match another. For a research pipeline working with a published dataset, the public Broad PON is usually adequate. For a pipeline being deployed in a new sequencing facility, build your own PON from internal controls (`gatk Mutect2 --tumor` on each control, then `CreateSomaticPanelOfNormals` to combine).

## 10. The germline resource

The germline resource is a population-level allele-frequency VCF. gnomAD (Karczewski et al. 2020, *Nature* 581:434) is the canonical source; the Broad publishes a Mutect2-compatible `af-only-gnomad.hg38.vcf.gz` that contains every gnomAD variant with its population allele frequency. Mutect2 uses the germline resource to compute the probability that a variant in the tumor (and possibly in the normal at low AF, due to contamination or sub-clonality) is in fact a known germline variant; FilterMutectCalls uses the same annotation to filter likely germline contamination.

The germline resource is also reference-build-specific. The hg38 gnomAD VCF works with GRCh38 alignments; the hg19 gnomAD VCF works with GRCh37 alignments. **Mixing builds silently produces wrong calls**: the coordinates in the germline resource will not match the coordinates in the BAM, and Mutect2 will fail to annotate variants that should be flagged as germline, and will annotate variants that should not be. We pin GRCh38 throughout this week.

## 11. The FilterMutectCalls filter set

`FilterMutectCalls` is GATK's post-Mutect2 filter applied to the raw VCF. It annotates each variant's FILTER column with one of `PASS` (no filter triggered) or a comma-separated list of filter reasons. The standard filter set:

- `germline` — the variant is likely a germline variant based on the gnomAD allele frequency and the tumor-normal read-count distribution.
- `panel_of_normals` — the variant matches a PON entry.
- `clustered_events` — the variant is one of many SNVs within a short window; usually a sign of a misaligned region or a complex variant.
- `multiallelic` — the variant has more than two alleles at the site; complex calls require manual review.
- `weak_evidence` — the variant's likelihood does not meet the caller's threshold.
- `strand_bias` — the variant's supporting reads are predominantly on one strand; suggests a strand-specific sequencing artefact.
- `base_qual` — the variant's supporting reads have low base qualities at the variant position.
- `map_qual` — the variant's supporting reads have low mapping qualities; suggests the region is ambiguous.
- `position` — the variant occurs predominantly at the start or end of reads; suggests a read-end artefact.
- `contamination` — the variant's AF is consistent with the contamination estimate.
- `n_ratio` — too many reads at the position are softclipped or low-quality.
- `read_position` — similar to `position`, with a different test statistic.
- `fragment_length` — the variant's supporting reads have an unusual fragment-length distribution.
- `slippage` — the variant is in a homopolymer / short-tandem-repeat context and the supporting reads may be slippage artefacts.

A variant can have multiple filter reasons; a clean PASS variant has none. The mini-project counts the variants in each filter category and reports the counts in the Markdown report.

## 12. What we have not done yet

Lecture 1 has built the conceptual scaffolding: matched pairs, clonal evolution, purity, contamination, CNV-at-a-distance, PON, germline resource, filter set. We have not yet:

- Implemented the Mutect2 call in Python via subprocess — that is Lecture 2 and Exercise 1.
- Run FilterMutectCalls with the contamination table — Lecture 2 and Exercise 2.
- Run Strelka2 as a cross-check — Challenge 1.
- Computed the 96-class trinucleotide spectrum — Lecture 3 and Exercise 3.
- Decomposed the spectrum into COSMIC v3.3 signatures — Lecture 3 and Exercise 3.
- Looked up variants in OncoKB and CIViC — Challenge 2.
- Wrapped the full pipeline as a CLI script with a `run-info.json` — the mini-project.

Each of those is a focused chunk of work. The matched-pair model and the assumptions of this lecture organize them all.

## 13. A worked example — interpreting an AF distribution

Suppose your tumor BAM has 50x mean coverage and Mutect2 emits 240 candidate variants. You plot the distribution of tumor allele frequencies:

```text
AF bin        count
0.00-0.05      82
0.05-0.10      37
0.10-0.20      28
0.20-0.30      24
0.30-0.40      31
0.40-0.50      19
0.50-0.60      14
0.60-1.00       5
```

Three observations:

- The low-AF bins (0.00-0.10) dominate. These are not all true low-AF somatic variants; many are sequencing errors that survived the TLOD threshold and many are likely to fail the `weak_evidence` filter. Without a matched normal and a PON, distinguishing real low-AF somatic from noise here is hard.
- The 0.30-0.40 bin is the peak of "truncal heterozygous somatic at this purity". If you computed `purity * 0.5` and got 0.35, this peak is where you expected it.
- The 0.40-0.50 bin is suspicious. A heterozygous truncal somatic at 80% purity would peak at 0.40; a homozygous germline variant *also* peaks at 0.50. The `germline` filter should catch the germline-derived peak; if it does not, the matched normal is not doing its job (low coverage in the normal, or the variant is at low population frequency in gnomAD).
- The 0.50-1.00 tail is informative. Most "purely tumor" variants do not reach AF > 0.5 unless there is copy-number loss of the wild-type allele. A spike of high-AF variants in a particular chromosomal region is a flag for LOH; a CNV caller would confirm.

The AF distribution is one of the first diagnostics to run after a Mutect2 / FilterMutectCalls pass. The shape tells you about purity (the peak position), about contamination (a sub-peak at the contaminant's germline AF), about LOH (high-AF tails), and about caller noise (low-AF excess).

## 14. The sample-naming and provenance contract

The single most damaging silent failure in a Mutect2 pipeline is mis-tagging the tumor and normal BAMs. The contract:

- Each BAM has an `@RG SM:` value in its header. This is the sample name Mutect2 will look up.
- The Mutect2 command line passes `-tumor SAMPLE_NAME` and `-normal SAMPLE_NAME`. The names must match the BAM headers.
- The pipeline's `run-info.json` records both sample names and confirms they came from the BAM headers (not from the command line). This way, a downstream consumer can verify the sample-name attribution is correct.

If the same BAM is used as both tumor and normal (a self-self test), the verification step rejects it. If two BAMs share the same `SM:` value (a library-prep mis-tag), the verification step rejects it. The wrapper script for Exercise 1 implements this and the mini-project's pipeline calls the same wrapper.

The reproducibility chain is: BAM headers → wrapper verification → Mutect2 → FilterMutectCalls → run-info JSON. Each step traces the sample identity forward; a bad input cannot silently produce a defensible-looking output.

## 15. References

- Cibulskis K, Lawrence MS, Carter SL, et al. *Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples.* **Nature Biotechnology** 31:213-219 (2013). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/>.
- Karczewski KJ, Francioli LC, Tiao G, et al. *The mutational constraint spectrum quantified from variation in 141,456 humans.* **Nature** 581:434-443 (2020). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/>.
- Nowell PC. *The clonal evolution of tumor cell populations.* **Science** 194:23-28 (1976).
- Stratton MR, Campbell PJ, Futreal PA. *The cancer genome.* **Nature** 458:719-724 (2009). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2821689/>.
- McGranahan N, Swanton C. *Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future.* **Cell** 168:613-628 (2017). Free at <https://www.cell.com/cell/fulltext/S0092-8674(17)30100-7>.
- Kastenhuber ER, Lowe SW. *Putting p53 in Context.* **Cell** 170:1062-1078 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5743327/>.
- GATK Best Practices for somatic short-variant discovery: <https://gatk.broadinstitute.org/hc/en-us/articles/360035531132>.

## 16. Self-check

You should be able to answer the following without looking back:

- Why does the matched-normal model work and what does it not catch?
- Define purity and contamination. Which one is corrected by `--contamination-table`?
- A heterozygous somatic variant in a 40%-pure tumor sample at 50x coverage. What is the expected observed AF?
- The patient has chronic lymphocytic leukemia. You sequence blood as the normal. What goes wrong and how do you fix it?
- A tumor sample produces a PASS variant in TP53 at 30% AF and a CNV-deletion call covering the wild-type TP53 allele. Is the TP53 lesion monoallelic or biallelic? Why does this matter for interpretation?
- What is a panel of normals and why is it platform-specific?
- A pipeline pinned to GRCh38 is run with a gnomAD VCF on hg19 coordinates. What goes wrong and how would you notice?
- Why is the `germline` filter the most common filter reason in a typical Mutect2 PASS report?

If any of these are not crisp, re-read §3-§5 (matched normal, purity, contamination), §6 (CNV-at-distance), §9-§10 (PON, germline resource), and §11 (filter set) before moving on.

---

Continue to [Lecture 2 — Mutect2 and Strelka2: callers, filters, and the pipeline](./02-mutect2-and-strelka2.md).
