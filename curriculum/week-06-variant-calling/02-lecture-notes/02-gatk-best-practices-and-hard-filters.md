# Lecture 2 — GATK Best Practices and Hard Filters

> **Duration:** ~2.5 hours of reading + hands-on command line + a short VEP annotation walk-through.
> **Outcome:** You can run GATK `HaplotypeCaller` on a sorted, indexed BAM, apply the six SNP hard filters and four indel hard filters from the Best Practices recipe with `bcftools filter`, annotate the surviving variants with Ensembl VEP, and read the resulting `CSQ` field to identify the consequence of each variant on every overlapping transcript.

Lecture 1 produced a raw VCF from `bcftools call`. Lecture 2 turns that raw VCF into the methods-section-quality output you can hand to a biologist: a filtered, annotated VCF with per-variant consequences and a transparent provenance trail.

If you only remember one thing from this lecture, remember this:

> **The GATK Best Practices hard filters are not a magic recipe; they are six (for SNPs) or four (for indels) thresholds on metrics that capture systematic alignment and base-call errors. Each threshold was calibrated against high-confidence truth sets (the original NA12878 1000 Genomes data, the GIAB benchmark set) and tuned to remove ≥ 95% of false positives while preserving ≥ 99% of true positives. You can run a stricter recipe, you can run a looser recipe, but if you cannot say what each metric measures and why its threshold is where it is, you are running the recipe on faith. This lecture removes the faith.**

---

## 1. GATK `HaplotypeCaller` — the local-reassembly caller

`bcftools call` is a column-by-column caller: at each reference position, look at the reads, run the likelihood model, emit a call. This works well when the reference alignment is correct. It struggles in regions where the reads disagree with the reference because of *indels* — an unaligned indel in the read produces mismatches at the surrounding positions, and the per-column likelihood model interprets those mismatches as multiple independent SNPs rather than one shifted alignment.

GATK `HaplotypeCaller` (Poplin et al. 2018) solves this with **local haplotype reassembly**: in regions where the reads disagree with the reference, it re-assembles candidate haplotypes from the reads using a De Bruijn graph (the same data structure as a short-read assembler), then re-aligns the reads against each candidate haplotype and picks the haplotype combination that best explains the data. The variants emitted are the differences between the chosen haplotypes and the reference. This is more accurate than `bcftools` at indels and at clusters of nearby SNPs, but more expensive (~10x slower).

### 1.1 The four steps of HaplotypeCaller

1. **Active region detection.** Walk the BAM. At each position, compute a per-position activity score from the reads (mismatch count, soft-clip count, indel count). Stretches of high-activity positions become "active regions." Inactive regions are skipped — every variant call comes from an active region.
2. **Local assembly.** For each active region (typically 100-300 bp), build a De Bruijn graph from the reads. The reference is also added as a path through the graph. Each non-reference path through the graph is a candidate haplotype.
3. **Likelihood computation.** For each read, compute the likelihood of the read under each candidate haplotype using a pair-Hidden-Markov-Model that accounts for base-call errors, indel-call errors, and a position-dependent mismatch model. This is the heavy step computationally.
4. **Genotype assignment.** Pick the pair of haplotypes (for diploid; one haplotype for haploid) that maximizes the posterior over the reads. Emit the variants implied by the chosen haplotypes.

### 1.2 A minimal HaplotypeCaller command

```bash
gatk HaplotypeCaller \
    -R ref/ecoli.fa \
    -I aln/SRR1770413.markdup.bam \
    -O calls/SRR1770413.gatk.vcf.gz \
    --sample-ploidy 1
```

For *E. coli* at 50x coverage, this takes ~5 minutes (vs ~30 seconds for `bcftools`). For human germline at 30x, ~12-24 hours single-threaded; in production it is parallelized over genomic intervals (e.g., per chromosome) with `--intervals chr1.bed`.

The output is a single-sample VCF, structurally identical to `bcftools call`'s output. The `INFO` field has GATK-specific keys (`QD`, `FS`, `MQ`, `MQRankSum`, `ReadPosRankSum`, `SOR`) instead of `bcftools`'s (`SP`, `INDEL`, `MQ0`), and the genotype likelihoods are computed by the pair-HMM model rather than by the bcftools binomial-plus-error model.

### 1.3 Joint genotyping with gVCFs

For multi-sample studies (e.g., the 1000 Genomes Project, or a cancer trio of tumor + normal + matched-normal), HaplotypeCaller can emit a **gVCF** (`-ERC GVCF` flag), which contains a record at every position including non-variant ones, with a "could be a variant if you knew about other samples" block representation. Multiple gVCFs are then joint-genotyped with `gatk GenotypeGVCFs` to produce a single multi-sample VCF where every variant is genotyped in every sample. This is the production pipeline for large cohorts; for single-sample analyses (the mini-project), the simpler `-O calls/x.vcf.gz` direct VCF output is enough.

---

## 2. `bcftools call` vs GATK `HaplotypeCaller` — when to use each

| Question | Use |
|----------|-----|
| Single bacterial sample, ~50x coverage, fast iteration | `bcftools call -m --ploidy 1` |
| Single human germline sample, 30x coverage | Either. GATK is the clinical-grade default. `bcftools` is faster and ~95% equivalent. |
| Multi-sample cohort (≥ 10 samples) for joint analysis | GATK `HaplotypeCaller -ERC GVCF` + `GenotypeGVCFs` |
| Cancer somatic variant calling | Neither; use a somatic caller (Mutect2 from GATK, or `bcftools call -C alleles` for known sites) |
| Population-genetics study (allele frequencies) | `bcftools` with `--ploidy-file` for sex chromosomes |
| Indel-heavy region (homopolymers, microsatellites) | GATK (the local reassembly is more accurate at indels) |
| Production clinical pipeline | GATK (the validated default in CAP/CLIA-certified labs) |

**Default to `bcftools` for single-sample analyses outside clinical contexts.** It is faster, the binary is smaller (4 MB vs 1 GB), and the output is ~95% equivalent. Switch to GATK when (a) the downstream consumer requires GATK output (clinical pipelines, the GATK ecosystem), (b) you have ≥ 10 samples to joint-genotype, or (c) your indel calls matter enough to pay the 10x runtime cost.

The mini-project this week uses `bcftools` (matching the bacterial single-sample profile). Challenge 1 asks you to run both and compare.

---

## 3. Indel normalization — why `bcftools norm` is non-optional

A single biological indel can be represented in multiple ways in a VCF. Consider a 2-bp deletion in the reference sequence `AAATAAA`:

```
Position:  1 2 3 4 5 6 7
Reference: A A A T A A A
Sample:    A A T A A A         (deletion of two A's)
```

The deletion can equivalently be encoded as:

- `POS=1, REF=AAA, ALT=A`  (left-most: anchor at position 1)
- `POS=2, REF=AAT, ALT=A`  (shifted one to the right)
- `POS=3, REF=ATA, ALT=A`  (shifted two to the right)

All three encodings describe the *same* biological event but record different `POS` and `REF` strings. When you compare two VCFs (e.g., `bcftools` vs GATK; or a sample VCF vs the GIAB truth set), records that describe the same variant in different encodings will look different and be miscounted as disagreements.

**The convention is left-alignment.** `bcftools norm -f ref.fa` shifts every indel as far left as the reference allows. After normalization, the deletion above is always encoded as `POS=1, REF=AAA, ALT=A`. Always normalize before comparing, before annotating with VEP, and before intersecting with public variant databases.

`bcftools norm` also **splits multiallelic records** into one record per ALT allele (when given `-m -any`), and **joins** them back (`-m +any`). The split form is required by most downstream tools (VEP, the comparison harness); the joined form is more compact and is the bcftools default emission.

```bash
bcftools norm -f ref/ecoli.fa -m -any -Oz -o calls/sample.norm.vcf.gz calls/sample.vcf.gz
bcftools index -t calls/sample.norm.vcf.gz
```

After this step, every indel is left-aligned and every multiallelic site is split. The VCF is now in **canonical form**, the form every downstream tool expects.

---

## 4. The GATK Best Practices hard filters

The Best Practices recipe (<https://gatk.broadinstitute.org/hc/en-us/articles/360035890471>) defines six hard filters for SNPs and four for indels. Each filter is a threshold on an `INFO` field computed by GATK `HaplotypeCaller` (`bcftools` computes some of these but not all; see §4.2 for the equivalents).

### 4.1 The six SNP filters

| Filter | Expression | What it removes | Failure mode |
|--------|------------|-----------------|--------------|
| `QD2` | `QD < 2.0` | Low quality per unit of depth. | Variants supported by many low-quality reads (artifact at high coverage). |
| `FS60` | `FS > 60.0` | High Fisher's strand bias. | Variants seen on one strand only (PCR artifact or strand-specific sequencing error). |
| `MQ40` | `MQ < 40.0` | Low average mapping quality. | Variants in multimapper / repetitive regions where reads were placed by coin-flip. |
| `MQRankSum-12.5` | `MQRankSum < -12.5` | Mapping-quality difference between ref and alt reads. | Variants where the alt-supporting reads are systematically worse-mapped than the ref-supporting reads. |
| `ReadPosRankSum-8` | `ReadPosRankSum < -8.0` | Read-position bias. | Variants where the alt allele is preferentially at the ends of reads (alignment-end soft-clip artifact). |
| `SOR3` | `SOR > 3.0` | Symmetric strand-odds-ratio. | An alternate strand-bias metric, more robust to large allele counts than FS. |

A SNP is excluded if **any** of these filters fires. The expression to apply all six at once with `bcftools filter`:

```bash
bcftools filter -Oz -o calls/sample.snp.filt.vcf.gz \
    -s LowQual \
    -e 'TYPE="snp" && (QD<2.0 || FS>60.0 || MQ<40.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0)' \
    calls/sample.gatk.vcf.gz
```

### 4.2 The four indel filters

| Filter | Expression | What it removes |
|--------|------------|-----------------|
| `QD2` | `QD < 2.0` | Same as SNPs. |
| `FS200` | `FS > 200.0` | High strand bias. The threshold is **much looser** than for SNPs (200 vs 60) because indels are intrinsically more strand-biased: a real indel is preferentially detected on one strand because the alignment algorithm soft-clips the other strand's reads. |
| `ReadPosRankSum-20` | `ReadPosRankSum < -20.0` | Same idea as SNPs with a looser threshold. |
| `SOR10` | `SOR > 10.0` | Loose SOR for indels. |

The indel-specific expression:

```bash
bcftools filter -Oz -o calls/sample.indel.filt.vcf.gz \
    -s LowQual \
    -e 'TYPE="indel" && (QD<2.0 || FS>200.0 || ReadPosRankSum<-20.0 || SOR>10.0)' \
    calls/sample.gatk.vcf.gz
```

For a complete pipeline, filter SNPs and indels separately, then merge:

```bash
# Filter SNPs.
bcftools view -v snps calls/sample.gatk.vcf.gz \
| bcftools filter -e '<SNP-expression>' -s LowQualSNP -Oz -o calls/snp.vcf.gz

# Filter indels.
bcftools view -v indels calls/sample.gatk.vcf.gz \
| bcftools filter -e '<INDEL-expression>' -s LowQualIndel -Oz -o calls/indel.vcf.gz

# Merge.
bcftools concat -a -Oz -o calls/sample.filtered.vcf.gz calls/snp.vcf.gz calls/indel.vcf.gz
bcftools index -t calls/sample.filtered.vcf.gz
```

### 4.3 The bcftools equivalents

If you called variants with `bcftools` rather than GATK, you do not have all the GATK metrics (`QD`, `MQRankSum`, `ReadPosRankSum`, `SOR`). `bcftools` computes a similar set under different names:

| GATK metric | bcftools equivalent | Notes |
|-------------|---------------------|-------|
| `QD` | `INFO/QUAL / INFO/DP` | Compute as `(QUAL/DP)`; no direct field. |
| `FS` | `INFO/SP` | Phred-scaled instead of GATK's natural-log scale, but same idea. |
| `MQ` | `INFO/MQ` | Identical definition. |
| `MQRankSum` | (not computed) | bcftools does not emit this. |
| `ReadPosRankSum` | (not computed) | bcftools does not emit this. |
| `SOR` | (not computed) | bcftools does not emit this. |

For the `bcftools`-only bacterial mini-project, the simplified recipe is:

```bash
bcftools filter -Oz -o calls/sample.filtered.vcf.gz \
    -s LowQual \
    -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' \
    calls/sample.raw.vcf.gz
```

This catches the most obvious failure modes (low QUAL, low depth, low mapping quality, high strand bias) without requiring the GATK-specific metrics. For human germline data, run GATK `HaplotypeCaller` and apply the full six-filter recipe.

### 4.4 What "passes filter" actually means

After `bcftools filter`, every variant has a `FILTER` value of either `PASS` or one of your filter tag names. To extract only the PASS variants:

```bash
bcftools view -f PASS calls/sample.filtered.vcf.gz -Oz -o calls/sample.pass.vcf.gz
```

`PASS` does **not** mean "this variant is biologically real." It means "this variant survived the filters I configured." If you configured loose filters, your PASS set is large and includes false positives. If you configured strict filters, your PASS set is small and excludes true positives. The job of QC is to pick filters whose `PASS` set has the right precision/recall trade-off for your downstream consumer.

A common confusion: variants emitted by `bcftools call` initially have `FILTER=.` (missing, "not yet filtered"). They are *not* `PASS` until `bcftools filter` runs. Some downstream tools treat `.` as equivalent to `PASS`, but it is best practice to run the filter step and produce explicit `PASS` / `LowQual` values.

---

## 5. Ensembl VEP — turning positions into consequences

A VCF is a list of positions where the sample differs from the reference. What a biologist actually wants is a list of *consequences*: this variant breaks gene X, this variant has no functional impact, this variant hits a known disease locus. The conversion is done by a **variant effect predictor**, of which the two production options are:

- **Ensembl VEP** (McLaren et al. 2016) — the Ensembl Project's annotator. Handles every species with an Ensembl genome assembly. Supports REST API for one-off lookups and offline cache for production.
- **SnpEff** (Cingolani et al. 2012) — a Java-based annotator with broadly similar functionality. Slightly different consequence ranking; agrees with VEP on ~95% of calls.

VEP is the C10 default. The output is a VCF identical to the input plus a `CSQ` field in `INFO` containing a comma-separated list of consequence annotations, one per overlapping transcript.

### 5.1 The VEP offline cache

VEP runs in three modes: (a) online, querying the Ensembl REST API at <https://rest.ensembl.org/>; (b) offline cache, querying a pre-downloaded species cache on disk; (c) custom, querying a local GFF/GTF file. The cache is the production choice — fast (10,000 variants/sec for human), reproducible, no network dependency.

```bash
# One-time cache install for E. coli K-12.
vep_install --AUTO c --SPECIES escherichia_coli_str_k_12_substr_mg1655 --CACHEDIR vep_cache/

# Cache install for human (~25 GB).
vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR vep_cache/
```

The first run downloads ~25 MB (*E. coli*) or ~25 GB (human) of pre-computed transcript data. Subsequent runs use the on-disk cache.

### 5.2 Running VEP on a filtered VCF

```bash
vep --input_file calls/sample.filtered.vcf.gz \
    --output_file calls/sample.vep.vcf \
    --species escherichia_coli_str_k_12_substr_mg1655 \
    --cache --dir_cache vep_cache/ \
    --vcf --symbol --canonical --force_overwrite
```

Flag-by-flag:

- **`--input_file`** — the filtered VCF (bgzipped or plain).
- **`--output_file`** — output path. If `--vcf` is set, the output is a VCF; otherwise, a tab-separated text table.
- **`--species`** — species identifier. The full list is at <https://www.ensembl.org/info/about/species.html>.
- **`--cache --dir_cache`** — use the offline cache rather than the REST API.
- **`--vcf`** — emit output as VCF with `CSQ` field, not as tab-separated text. This is what you want when feeding downstream tools.
- **`--symbol`** — add HUGO gene symbols (e.g. `rpoB` rather than `Ensembl ID b3987`).
- **`--canonical`** — mark the canonical transcript per gene (the one Ensembl considers most representative).
- **`--force_overwrite`** — overwrite the output file if it exists. Default is to error.

For an *E. coli* VCF with 200 variants, VEP runs in ~2 seconds. For a human VCF with 3 million variants, ~30 minutes offline.

### 5.3 The `CSQ` field — what comes back

After VEP, every variant has a `CSQ` field in `INFO` containing a pipe-delimited string with one entry per overlapping transcript. The fields and their order are documented in the VCF header line `##INFO=<ID=CSQ,...>`:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|...">
```

A real annotated record:

```
NC_000913.3 3947830 . C T 198.0 PASS DP=42;MQ=60;AF=1;CSQ=T|missense_variant|MODERATE|rpoB|b3987|Transcript|AAC76964.1|protein_coding|1/1||c.1547C>T|p.Ala516Val|1547|1547|516|A/V|gCa/gTa|||  GT:DP:AD 1:42:0,42
```

The `CSQ` field decodes as:

- `T` — the alt allele.
- `missense_variant` — the consequence.
- `MODERATE` — the IMPACT category (LOW, MODERATE, HIGH, MODIFIER).
- `rpoB` — the gene symbol (RNA polymerase β subunit, a famous antibiotic-resistance locus).
- `b3987` — the Ensembl gene ID.
- `Transcript` — the feature type.
- `AAC76964.1` — the transcript ID.
- `protein_coding` — the biotype.
- `1/1` — exon 1 of 1.
- `c.1547C>T` — the HGVS coding-sequence notation.
- `p.Ala516Val` — the HGVS protein notation (alanine to valine at codon 516).
- `516` — the protein position.
- `A/V` — the amino acid change.
- `gCa/gTa` — the codon change (lowercase = unchanged bases, uppercase = the variant base).

If the variant overlaps multiple transcripts, `CSQ` has multiple comma-separated entries. The `--canonical` flag adds a `CANONICAL` field to mark which entry corresponds to the canonical transcript; use this to pick one consequence per variant for downstream summary tables.

### 5.4 Parsing CSQ in Python

```python
import pysam

vcf = pysam.VariantFile("calls/sample.vep.vcf")
csq_format = vcf.header.info["CSQ"].description.split("Format: ")[1].split("|")

for rec in vcf:
    for csq in rec.info.get("CSQ", []):
        fields = dict(zip(csq_format, csq.split("|")))
        if fields.get("CANONICAL") == "YES":
            print(rec.chrom, rec.pos, fields["SYMBOL"],
                  fields["Consequence"], fields.get("Amino_acids", ""))
```

Output for the *E. coli* mini-project:

```
NC_000913.3 3947830 rpoB missense_variant A/V
NC_000913.3 4148500 gyrA missense_variant S/L
...
```

This is the format a biologist actually wants to read. The variant-calling pipeline ends here.

---

## 6. The complete pipeline end to end

Putting it all together, the canonical Week 6 pipeline from sorted-markdup BAM to annotated VCF:

```bash
# Step 1: variant calling.
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/SRR1770413.markdup.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/SRR1770413.raw.vcf.gz - \
&& bcftools index -t calls/SRR1770413.raw.vcf.gz

# Step 2: hard filtering.
bcftools filter -Oz -o calls/SRR1770413.filtered.vcf.gz \
    -s LowQual \
    -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' \
    calls/SRR1770413.raw.vcf.gz \
&& bcftools index -t calls/SRR1770413.filtered.vcf.gz

# Step 3: normalization.
bcftools norm -f ref/ecoli.fa -m -any -Oz \
    -o calls/SRR1770413.norm.vcf.gz \
    calls/SRR1770413.filtered.vcf.gz \
&& bcftools index -t calls/SRR1770413.norm.vcf.gz

# Step 4: VEP annotation.
vep --input_file calls/SRR1770413.norm.vcf.gz \
    --output_file calls/SRR1770413.vep.vcf \
    --species escherichia_coli_str_k_12_substr_mg1655 \
    --cache --dir_cache vep_cache/ \
    --vcf --symbol --canonical --force_overwrite

# Step 5: PASS-only set for downstream consumers.
bcftools view -f PASS calls/SRR1770413.vep.vcf -Oz \
    -o calls/SRR1770413.pass.vcf.gz \
&& bcftools index -t calls/SRR1770413.pass.vcf.gz

# Step 6: stats summary.
bcftools stats calls/SRR1770413.pass.vcf.gz > results/vcf_stats.txt
```

End-to-end runtime on *E. coli*: ~1 minute. End-to-end runtime on human germline: ~3-5 hours. This is the pipeline the mini-project asks you to build.

---

## 7. Comparing two callers — `bcftools isec`

Challenge 1 this week asks you to run `bcftools call` and GATK `HaplotypeCaller` on the same BAM and compare. The comparison tool is `bcftools isec`:

```bash
bcftools isec -p compare/ -Oz \
    calls/SRR1770413.bcftools.norm.vcf.gz \
    calls/SRR1770413.gatk.norm.vcf.gz
```

The `-p compare/` flag writes four output VCFs:

- `compare/0000.vcf.gz` — variants private to the first caller (in bcftools only).
- `compare/0001.vcf.gz` — variants private to the second caller (in GATK only).
- `compare/0002.vcf.gz` — variants in both callers (intersection, with bcftools annotations).
- `compare/0003.vcf.gz` — variants in both callers (intersection, with GATK annotations).

The expected agreement on a well-behaved bacterial BAM is ~95% at the SNP level, ~85% at the indel level. Indel disagreements cluster at homopolymer runs; SNP disagreements cluster at multiallelic sites and at low-coverage positions. Challenge 1 walks you through the analysis in detail.

---

## 8. Common misconceptions

- **"The hard filters always reject false positives and accept true positives."** No. Hard filters are thresholds on noisy metrics. At any threshold, some true positives fall below the threshold (false negatives) and some false positives slip above (false positives). The thresholds in the Best Practices recipe trade off ~1% sensitivity for ~95% specificity. Stricter thresholds buy more specificity at the cost of sensitivity; looser thresholds the reverse.
- **"VEP picks the right consequence per variant."** It picks *every* consequence — one per overlapping transcript. If a variant overlaps three transcripts, you get three CSQ entries. The `--pick`/`--canonical`/`--most_severe` flags choose one for summary purposes, but the underlying data is one entry per transcript, and which one to keep is a downstream-analysis choice.
- **"`bcftools norm` is optional."** No. Without normalization, two VCFs representing the same biological variants will appear to disagree. Always run `bcftools norm -f ref.fa -m -any` before any cross-VCF operation.
- **"GATK is strictly better than bcftools."** No. GATK is more accurate at indels and at clusters; bcftools is faster, simpler, and equally accurate at SNPs in well-behaved data. The "right" caller depends on the question. For single-sample bacterial WGS, bcftools is the better choice (faster, simpler, the indel advantage of GATK is moot at typical bacterial sequencing depths).
- **"PASS variants are real, non-PASS are not."** PASS is a *traceability* tag, not a *truth* tag. The exact filter expression determines what falls into PASS; change the expression and the set changes. Always document the filter expression in your methods section.

---

## 9. Where this lecture lands you for the mini-project

After this lecture you should be able to:

- Run `bcftools call` and GATK `HaplotypeCaller` end to end on a sorted-markdup BAM.
- Normalize indel representation with `bcftools norm -f ref.fa -m -any`.
- Apply the six SNP and four indel hard filters from the Best Practices recipe.
- Annotate a filtered VCF with VEP, using the offline cache for bacterial species.
- Read a `CSQ` field and identify the consequence on the canonical transcript per gene.
- Compare two callers with `bcftools isec` and interpret the intersection / private sets.

The mini-project (`mini-project/README.md`) builds this pipeline end to end on the SRR1770413 BAM from Week 5 and writes up the findings. Every step in the above list is a numbered phase of the mini-project.

---

## Self-check questions

Before you move on, answer these without looking.

1. State the four steps of GATK `HaplotypeCaller`'s algorithm. (§1.1)
2. Why is local reassembly more accurate than column-by-column calling at indels? (§1)
3. What does `bcftools norm -f ref.fa` do, and why is it required before comparing two VCFs? (§3)
4. Name the six SNP hard filters and what each one removes. (§4.1)
5. Why is the FS threshold for indels (200) so much looser than for SNPs (60)? (§4.2)
6. What does the IMPACT category in VEP mean? Name the four values. (§5.3)
7. Decode this CSQ field: `T|missense_variant|MODERATE|rpoB|b3987|Transcript|AAC76964.1|protein_coding|1/1||c.1547C>T|p.Ala516Val`. (§5.3)
8. Why does `bcftools view -f PASS` not guarantee biological correctness? (§4.4)
9. Under what scenario would you choose `bcftools call` over GATK `HaplotypeCaller`? Name two. (§2)
10. What does `bcftools isec -p out/ A.vcf.gz B.vcf.gz` produce, and what are the four output files? (§7)

Answers are not provided. The answers are in the section references; do the work.

---

## Further reading

- Poplin, R. et al. (2018). Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv* 201178.
- McKenna, A. et al. (2010). The Genome Analysis Toolkit. *Genome Research* 20(9):1297-1303.
- McLaren, W. et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology* 17:122.
- Cingolani, P. et al. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly* 6(2):80-92.
- The GATK Best Practices: <https://gatk.broadinstitute.org/hc/en-us/articles/360035890471>.
- The VEP documentation: <https://www.ensembl.org/info/docs/tools/vep/index.html>.

---

*Continue to the [exercises](../03-exercises/00-overview.md) and the [mini-project](../07-mini-project/00-overview.md) once you have answered the self-check questions.*
