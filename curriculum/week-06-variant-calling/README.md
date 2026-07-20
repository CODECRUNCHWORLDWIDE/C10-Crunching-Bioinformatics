# Week 6 — Variant Calling

In Week 5 you took a billion-Smith-Waterman-cell problem and reduced it to `O(m)` per read using an FM-index, then sorted, indexed, duplicate-marked, and coverage-plotted the result. The output was a methods-section-quality BAM: 99.6% mapped, 98% properly paired, mean coverage 47.3x, CV 0.18 across 4.6 Mb of *E. coli* MG1655. Week 6 takes that BAM and asks the next question every short-read project asks: **where does this sample differ from the reference, and how much do I trust each difference?** The output is a **Variant Call Format (VCF)** file: a tab-separated record per polymorphic position, with an allele, a quality score, a depth, and enough metadata to defend the call to a reviewer. The canonical tools are **`bcftools mpileup` + `bcftools call`** (Li 2011, *Bioinformatics* 27:2987; Danecek et al. 2021, *GigaScience* 10:giab008) for the everyday Bayesian SNP-and-indel caller, **GATK `HaplotypeCaller`** (McKenna et al. 2010, *Genome Research* 20:1297; Poplin et al. 2018, *bioRxiv*) for the Broad Institute's local-haplotype-reassembly caller used in clinical pipelines, and **Ensembl VEP** (McLaren et al. 2016, *Genome Biology* 17:122) for annotating each variant with its predicted effect on every overlapping transcript. All three are free, all three are documented, and all three are the tools you will be running on real data from Monday onward.

By Friday of Week 6 you will be able to call SNPs and short indels from a sorted, duplicate-marked BAM with `bcftools mpileup -f ref.fa aln.bam | bcftools call -mv -Oz -o calls.vcf.gz`, read a VCF record column by column and explain every field (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, `<sample1>`), parse a VCF with `pysam.VariantFile` or `cyvcf2` and extract per-variant depth (`DP`), allele depth (`AD`), genotype quality (`GQ`), and strand bias (`SP`/`FS`), apply hard filters from the GATK Best Practices recipe (`QD < 2`, `FS > 60`, `MQ < 40`, `MQRankSum < -12.5`, `ReadPosRankSum < -8`, `SOR > 3` for SNPs; a different set for indels) and explain what each one is filtering out, run **Variant Effect Predictor** (VEP, online via the REST API or offline via the Docker image) against your filtered VCF and read the resulting consequence annotations (`missense_variant`, `synonymous_variant`, `stop_gained`, etc.), and compare two callers on the same BAM to quantify their disagreement at the per-position level. The mini-project takes the sorted+markdup BAM you built in Week 5 (SRR1770413 aligned against `NC_000913.3`), runs the full `bcftools` pipeline to produce a filtered, annotated VCF, and writes up the findings: how many variants, in which functional categories, with which transition/transversion ratio, and with what fraction filtered by each hard-filter rule.

The other half of the week is **what makes a variant real**. A naive caller would emit a variant at every position where at least one read disagrees with the reference, which on a 30x-coverage human BAM is tens of millions of "variants" — almost all of them sequencing errors, alignment artifacts, or duplicate-driven false positives. The Bayesian step in `bcftools call -m` and the local-reassembly step in `HaplotypeCaller` are filters that ask "given the reads here, what is the posterior probability that the reference allele is wrong?" — a question that integrates depth, base quality, mapping quality, allele frequency, strand bias, and the alignment-error model of the platform. The hard filters layered on top discard variants whose annotations look pathological even before genotype inference: variants called from one strand only, variants in alignments with anomalous mapping quality, variants at the very edge of read positions. The GATK Best Practices recipe codifies the thresholds that fifteen years of clinical-genomics experience have converged on; you will memorize the six SNP filters this week and never again ask "should I trust a `QD = 1.2` variant" (no — `QD < 2` is the standard exclusion). Variant effect prediction with VEP closes the loop: a position is now a *call*, a call is now a *consequence* (missense, synonymous, splice-region, intergenic), and a consequence is what a biologist actually wants to read.

## Learning objectives

By the end of this week, you will be able to:

- **Describe** the `bcftools mpileup` + `bcftools call -m` pipeline in two paragraphs, naming the per-position genotype-likelihood model (binomial sampling of reads from the two alleles of a diploid genotype, or the haploid equivalent for bacteria), and explain why this is faster than but lower-accuracy-on-indels than GATK `HaplotypeCaller`.
- **Choose** between `bcftools call` and GATK `HaplotypeCaller` for a given dataset (small bacterial genome vs human germline vs somatic tumor) and defend the choice in one sentence.
- **Run** `bcftools mpileup -f ref.fa aln.bam | bcftools call -mv -Oz -o calls.vcf.gz` end to end on a sorted, indexed, markdup BAM and produce a well-formed VCF.
- **Read** the eight mandatory VCF columns by name (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`) plus the per-sample columns (`FORMAT` + one column per sample), and decode every `INFO` and `FORMAT` field used by `bcftools` (`DP`, `AD`, `MQ`, `SP`, `GQ`, `GT`, `PL`).
- **Parse** a VCF with `pysam.VariantFile` or `cyvcf2`, iterate over records, and extract per-variant depth, allele depth, genotype quality, and predicted consequence.
- **Apply** the GATK Best Practices hard filters for SNPs (`QD < 2 || FS > 60 || MQ < 40 || MQRankSum < -12.5 || ReadPosRankSum < -8 || SOR > 3`) and indels (`QD < 2 || FS > 200 || ReadPosRankSum < -20 || SOR > 10`) using `bcftools filter -s 'LowQual' -e '<expression>'`, and explain what each filter is removing.
- **Annotate** a filtered VCF with Ensembl VEP (offline cache mode for production, or REST API for one-off lookups), and read the `CSQ` field to identify the consequence of each variant on every overlapping transcript.
- **Compute** standard variant-set QC metrics: transition/transversion ratio (`Ts/Tv`, ~2.0 for human germline, ~1.0 for bacterial random mutation), variants per Mb, fraction passing each filter, fraction by consequence category, and explain what an anomalous value of each metric is diagnostic of.
- **Compare** two variant callers (`bcftools call` vs GATK `HaplotypeCaller`) on the same BAM with `bcftools isec`, and interpret the intersection / private-to-A / private-to-B counts.
- **Identify** at least three failure modes of variant calling (low-coverage positions producing false negatives, PCR-duplicate stacks producing false positives if dedup was skipped, mapping-quality-zero reads in repetitive regions producing both) and the standard QC signal each one produces in a VCF.

## Prerequisites

This week assumes Weeks 1-5 are **done and committed**. Specifically:

- You have a sorted, indexed, duplicate-marked BAM from Week 5's mini-project (`SRR1770413.markdup.bam` aligned against `NC_000913.3`). Week 6 reuses it as input — no re-alignment.
- You can read SAM/BAM records with `pysam.AlignmentFile` and decode FLAG/CIGAR/MAPQ from Week 5 Exercise 2 and 3.
- You can call shell tools from Python with `subprocess.run(..., check=True)` from Week 5 Exercise 1.
- You have Python 3.11+, pandas, matplotlib, and the Week 5 conda env. You will need to install bcftools 1.19 and ensembl-vep 110 this week — `conda install -c bioconda bcftools=1.19 ensembl-vep=110` is the canonical path.

You do not need biology beyond "a variant is a position where the sample's sequence differs from the reference; SNPs are single-base substitutions, indels are small insertions or deletions, and structural variants are big rearrangements out of scope for this week." You do need familiarity with the SAM/BAM format from Week 5 — most VCF debugging eventually goes back to looking at the underlying BAM with `samtools tview` or `pysam.pileup`.

## Topics covered

- The variant-calling problem at scale: a 30x human genome BAM has ~3 Gb of reference, each position covered by ~30 reads, each read with ~0.1% base-call error. The naive "every disagreement is a variant" approach yields ~3 · 10^9 · 30 · 10^-3 ≈ 10^8 candidate variants per genome, of which ~3 · 10^6 (less than 5%) are real. Variant calling is fundamentally a *filtering* problem.
- The genotype-likelihood model: for a diploid position with two possible alleles `A` (reference) and `B` (alternative), the genotype is one of `AA`, `AB`, `BB`. Given `n` reads observed at the position with `k` carrying the `B` allele, `P(reads | genotype)` is binomial(`k`, `n`, allele-frequency-in-genotype) (`0` for AA, `0.5` for AB, `1` for BB). Apply Bayes' theorem with a prior on allele frequency (Hardy-Weinberg) to get a posterior over genotypes. Reference: Li 2011 (`bcftools`'s model), DePristo et al. 2011 (GATK's model with a more sophisticated likelihood).
- The `bcftools mpileup` + `bcftools call` pipeline: `mpileup` walks the BAM column by column, builds per-position read-allele matrices, and emits a BCF (binary VCF) with per-genotype likelihoods. `call -m` applies the multiallelic-and-rare-variant model (Li 2011 + Danecek 2014, *Bioinformatics* 30:2837) to assign a most-likely genotype per sample.
- The GATK `HaplotypeCaller` pipeline: in regions where the reads disagree with the reference, locally re-assemble candidate haplotypes from the reads using a De Bruijn graph (the same data structure as a short-read assembler — covered briefly here, in detail in a hypothetical Week 14). Re-align reads against the candidate haplotypes. Compute likelihoods of each haplotype. Output the variant calls implied by the most-likely haplotypes. Reference: Poplin et al. 2018.
- The VCF specification (Danecek et al. 2011, *Bioinformatics* 27:2156): eight mandatory tab-separated columns (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`) plus optional `FORMAT` and per-sample columns. Header lines start with `##` for metadata and `#` for the column header. The full spec is at <https://samtools.github.io/hts-specs/VCFv4.3.pdf> — ~40 pages, denser than SAM but worth reading once.
- The `INFO` field: a semicolon-separated list of `KEY=VALUE` annotations applied per variant (not per sample). Common ones from `bcftools`: `DP` (total depth), `AF` (allele frequency), `MQ` (RMS mapping quality), `SP` (Phred-scaled strand-bias p-value), `INDEL` (set if the variant is an indel).
- The `FORMAT` field and per-sample columns: a colon-separated list of keys (typical: `GT:PL:DP:SP:AD`) whose values are colon-separated in each sample column. `GT` is the genotype (`0/0`, `0/1`, `1/1` for diploid; `0`, `1` for haploid; `.` for missing). `PL` is the Phred-scaled genotype likelihoods.
- Hard filters from the GATK Best Practices (<https://gatk.broadinstitute.org/hc/en-us/articles/360035890471>). For SNPs: `QD < 2` (low quality per depth), `FS > 60` (Fisher's strand bias), `MQ < 40` (mean mapping quality), `MQRankSum < -12.5`, `ReadPosRankSum < -8`, `SOR > 3`. For indels: a different set with looser `FS` and `SOR` thresholds. Apply with `bcftools filter -s LowQual -e <expression>`.
- VQSR (Variant Quality Score Recalibration): GATK's machine-learning replacement for hard filters, trained on known high-confidence variant sets (dbSNP, HapMap, 1000 Genomes). VQSR needs ≥ 30 samples and a curated truth set, so it is out of scope for the Week 6 *E. coli* mini-project (one sample, no truth set) — hard filters are the right tool here. We mention VQSR for completeness.
- Variant Effect Predictor (VEP): given a VCF and a species transcript database, for each variant find all overlapping transcripts and predict the consequence on each (`missense_variant`, `synonymous_variant`, `splice_region_variant`, `stop_gained`, `frameshift_variant`, ~30 consequence terms total; the full list is at <https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html>). Annotations land in the `CSQ` field of `INFO`.
- Quality metrics for a variant set: **Ts/Tv ratio** (transitions: A↔G, C↔T; transversions: A↔T, A↔C, G↔T, G↔C; the ratio is ~2.0–2.1 for human whole-genome data, much lower for bacteria and very high for coding-region-only subsets). **Variants per Mb**: ~3,000 per Mb for human heterozygous SNPs, ~50–500 per Mb for bacterial sample-vs-reference, depending on phylogenetic distance. **% PASS**: fraction of called variants that survive hard filters; expect 70-95% for well-behaved data.
- Common failure modes: **low-coverage positions** (< 5x) produce false negatives because the likelihood model has too few observations; **PCR-duplicate stacks** (if duplicates were not marked) produce false-positive variant calls at the duplicate's position because all the same fragment is counted as independent observations; **multimappers** in repetitive regions produce false positives at every paralogous copy because the reads were placed by a coin-flip in `bwa mem`'s MAPQ-0 case; **soft-clipped reads at structural-variant breakpoints** produce false-positive indels because the caller sees the local mismatches but not the underlying large rearrangement.

## Weekly schedule

The schedule below adds up to approximately **35 hours**. Treat it as a target. Monday's lecture on the genotype-likelihood model is the hour that decides whether the rest of the week makes statistical sense — read it twice if needed.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | Genotype likelihoods, `bcftools mpileup`/`call`    |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | Running bcftools end to end, VCF format            |    1.5h  |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5h      |
| Wednesday | GATK Best Practices, hard filters                  |    1.5h  |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     7h      |
| Thursday  | VEP annotation, consequence categories             |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Mini-project deep work + write-up                  |    0h    |    1h     |     0h     |    0.5h   |   1h     |     2h       |    0h      |     4.5h    |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **7.5h**  | **2h**     | **3h**    | **6h**   | **8h**       | **2h**     | **34.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | bcftools, GATK Best Practices, VEP docs + reference papers |
| [lecture-notes/01-from-bam-to-vcf-bcftools.md](./lecture-notes/01-from-bam-to-vcf-bcftools.md) | The genotype-likelihood model, `bcftools mpileup` + `call`, the VCF format column by column |
| [lecture-notes/02-gatk-best-practices-and-hard-filters.md](./lecture-notes/02-gatk-best-practices-and-hard-filters.md) | GATK `HaplotypeCaller` overview, the six SNP hard filters and four indel hard filters, VEP annotation |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-call-with-bcftools.py](./exercises/exercise-01-call-with-bcftools.py) | Call variants on a tiny BAM end to end with `bcftools mpileup` + `call` via `subprocess` |
| [exercises/exercise-02-vcf-parse.py](./exercises/exercise-02-vcf-parse.py) | Parse a VCF file by reading text lines, decode `INFO` and `FORMAT` fields, build a per-variant pandas frame |
| [exercises/exercise-03-vep-annotate.py](./exercises/exercise-03-vep-annotate.py) | Annotate a filtered VCF with VEP (REST API for the exercise, offline cache documented for production), parse the `CSQ` field |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-compare-callers.md](./challenges/challenge-01-compare-callers.md) | Run `bcftools call` and GATK `HaplotypeCaller` on the same BAM, compute the intersection and private sets, explain the disagreements |
| [quiz.md](./quiz.md) | 10 multiple-choice questions on VCF, hard filters, VEP, and the variant-calling toolchain |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | Variant-call SRR1770413 against `NC_000913.3`, apply hard filters, annotate with VEP, write up findings |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("bcftools 1.19," "ensembl-vep 110," "GATK 4.5.0.0"). We cite tools by their paper ("bcftools, Li 2011 *Bioinformatics* 27:2987; Danecek 2021 *GigaScience* 10:giab008"). We say "168 SNPs and 22 indels called in 4.6 Mb of *E. coli* MG1655 reference at mean coverage 47.3x, of which 152 SNPs (90.5%) and 19 indels (86.4%) PASS the GATK hard filters, Ts/Tv = 0.94 (consistent with random bacterial mutation), and 87 of 152 PASS SNPs (57.2%) are coding-region missense" not "we found some SNPs." A variant count is a number. If your methods section uses the words "many variants" or "high quality" without a number, you have not written one yet.

## A note on the data size

VCF files are small compared to BAM. The mini-project VCF is:

- **~50 KB** as bgzipped + indexed VCF (`.vcf.gz` + `.tbi`).
- **~2 MB** uncompressed.
- **~3 MB** annotated VCF after VEP adds `CSQ` fields (one per overlapping transcript).

Add the *E. coli* GFF3 annotation file (~2 MB) for the VEP offline cache, and you are at ~10 MB of disk for the entire variant-calling step on top of the Week 5 BAM. This fits in a git commit; the BAM does not. Commit the VCFs and the QC text outputs; gitignore the BAMs.

For a human germline genome, the equivalent VCF is ~150 MB compressed (~3 million SNPs and ~500,000 indels). VEP-annotated, it grows to ~1 GB because every variant has 1-10 consequence annotations (one per overlapping transcript). Plan accordingly when you move to human data.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Read the `bcftools` paper (Li 2011) and the GATK `HaplotypeCaller` paper (Poplin 2018) end to end. The Poplin paper is particularly readable as a "what we changed from the 2010 paper" tour.
- Run GATK `HaplotypeCaller` on your mini-project BAM and run `bcftools isec` to compute the intersection. The expected agreement on a well-behaved bacterial BAM is ≥ 95% at the SNP level, ≥ 85% at the indel level. Indel disagreements are almost always near homopolymer runs (both callers struggle there in different ways).
- Apply VQSR (Variant Quality Score Recalibration) to the mini-project VCF using the GATK 4 documentation. You will fail — VQSR needs ≥ 30 samples and a curated truth set. The point is to understand *why* you failed: the technique is data-hungry by design.
- Install **`SnpEff`** (Cingolani et al. 2012, *Fly* 6:80) and re-annotate your VCF. Compare to VEP. The two annotators agree on ~95% of consequence calls; the disagreements are at splice-region edges and at variants spanning multiple transcripts where they choose different "canonical" transcripts.
- Reproduce the GATK hard-filter thresholds yourself: take a high-confidence truth set (the GIAB NA12878 small-variant truth VCF), plot each metric (`QD`, `FS`, `MQ`, etc.) for the truth-positive vs the truth-negative variants, and pick the threshold that gives 99% sensitivity. You will recover something close to (but not identical to) the Best Practices values — the Best Practices are a compromise across many datasets, and dataset-specific filters can do slightly better.

## Up next

Continue to [Week 7 — Transcriptomics intro](../week-07/) once you have pushed your mini-project to GitHub. Week 7 leaves DNA behind: instead of aligning DNA reads and calling variants, you align RNA reads (often with `STAR` or `salmon` instead of `bwa mem`) and count them against gene features to estimate expression. The pysam patterns from Week 5 and the pandas patterns from Week 4 will both come back; the variant-calling machinery from Week 6 is parked.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
