# Challenge 1 — Strelka2 Cross-Check on the Same BAM Pair

> **Educational and research use only.** Cross-caller comparisons are a research-quality validation step. They do not replace clinical-pipeline validation. Do not use Week 11 output to guide patient care.

## What you will do

Run **Strelka2** (Kim et al. 2018, *Nature Methods* 15:591; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>) on the same tumor-normal BAM pair you used for Exercise 1. Compute the PASS-PASS overlap between Strelka2's somatic SNV VCF and Mutect2's PASS variants from Exercise 2. Report the Jaccard index, the union, the intersection, and the per-caller-only sets. Identify three variants that disagree (PASS in one, filtered in the other) and explain *why* the callers disagreed.

This challenge is independent of Exercise 1 / 2 / 3; it does not require SigProfilerAssignment. It is a self-contained ~3 hour exercise that exercises both Strelka2's run-configuration model and the cross-caller comparison pattern.

## Why

Every variant caller has its own statistical model, its own filter heuristics, and its own systematic blind spots. A variant that PASSes in both Mutect2 and Strelka2 carries higher confidence than a variant that PASSes in only one. The Mutect2 + Strelka2 cross-check pattern is widely used in the cancer-genomics literature; the literature reports PASS-PASS overlap at AF >= 10% in the 80-95% range on real tumor-normal data.

For a clinical pipeline the cross-check is usually a sensitivity floor: "we report variants that PASS in either caller" with a per-variant flag for "PASS in both" vs "PASS in one only". The single-caller mode is faster but loses ~5-10% sensitivity on low-AF variants where the disagreement is largest.

## Prerequisites

- Exercises 1 and 2 done and `results/ex02/filtered.vcf.gz` written.
- Strelka2 2.9.10 installed (`conda install -c bioconda strelka=2.9.10`).
- samtools 1.20 and pysam 0.22.1.
- The same tumor and normal BAMs from Exercise 1.

If Strelka2 is not installed and you cannot install it, the challenge still has a value: read the Strelka2 paper (Kim et al. 2018), simulate the comparison using the same dataset's known-truth VCF (`data/truth_chr22.vcf.gz`), and report the expected agreement based on the published Strelka2 sensitivity / specificity.

## Suggested approach

### Step 1 — read the Strelka2 paper (20 min)

Skim Kim et al. 2018 (free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>). Pay attention to:

- The **EVS (Empirical Variant Score)**: a Random Forest trained on truth sets. Strelka2's filter strategy is a single EVS threshold (~5-10 depending on variant type), not the multi-filter approach of Mutect2.
- The **separate SNV and indel pipelines**: Strelka2 emits two VCFs (`somatic.snvs.vcf.gz`, `somatic.indels.vcf.gz`); Mutect2 emits one combined VCF.
- The **per-position Bayesian model**: no active-region haplotype assembly; the SNV caller is position-by-position. This makes Strelka2 faster than Mutect2 on simple variants and weaker on indels in complex regions.

### Step 2 — run Strelka2 (30 min)

Strelka2 has a two-phase run model. Phase 1 builds a working directory and a runtime configuration; phase 2 runs the actual calling.

```bash
mkdir -p results/ch01
configureStrelkaSomaticWorkflow.py \
  --normalBam data/normal_chr22.bam \
  --tumorBam data/tumor_chr22.bam \
  --referenceFasta data/chr22_GRCh38.fasta \
  --runDir results/ch01/strelka_run \
  --callRegions data/chr22_intervals.bed.gz

# Phase 2:
results/ch01/strelka_run/runWorkflow.py -m local -j 4
```

Output VCFs land in `results/ch01/strelka_run/results/variants/`:

- `somatic.snvs.vcf.gz` — SNV calls with EVS.
- `somatic.indels.vcf.gz` — indel calls with EVS.

PASS variants are those with FILTER = `PASS`.

Note: `--callRegions` expects a *bgzipped* BED with a tabix index. Run `bgzip data/chr22_intervals.bed && tabix -p bed data/chr22_intervals.bed.gz` once.

### Step 3 — normalize both VCFs (15 min)

Cross-caller comparison requires that both VCFs are *normalized* to the same canonical form. This means splitting multi-allelic sites and left-aligning indels.

```bash
bcftools norm -f data/chr22_GRCh38.fasta -m -any results/ex02/filtered.vcf.gz \
  | bcftools sort -O z -o results/ch01/mutect2_norm.vcf.gz
tabix -p vcf results/ch01/mutect2_norm.vcf.gz

bcftools norm -f data/chr22_GRCh38.fasta -m -any \
  results/ch01/strelka_run/results/variants/somatic.snvs.vcf.gz \
  | bcftools sort -O z -o results/ch01/strelka_snvs_norm.vcf.gz
tabix -p vcf results/ch01/strelka_snvs_norm.vcf.gz

bcftools norm -f data/chr22_GRCh38.fasta -m -any \
  results/ch01/strelka_run/results/variants/somatic.indels.vcf.gz \
  | bcftools sort -O z -o results/ch01/strelka_indels_norm.vcf.gz
tabix -p vcf results/ch01/strelka_indels_norm.vcf.gz
```

### Step 4 — compute the PASS-PASS overlap (45 min)

Write a Python script `cross_check.py` that:

- Reads each normalized VCF.
- Extracts PASS variants as `(chrom, pos, ref, alt)` tuples.
- Computes the union, intersection, Mutect2-only, and Strelka2-only sets.
- Computes the Jaccard index `J = |intersection| / |union|`.
- Writes a Markdown report `cross_check.md` with the counts and the per-set first 10 variants.

A reference skeleton:

```python
from pathlib import Path
import pysam
from typing import Set

VariantKey = tuple[str, int, str, str]

def pass_variants(vcf_path: Path) -> set[VariantKey]:
    """Return the set of (chrom, pos, ref, alt) tuples for PASS variants."""
    out: set[VariantKey] = set()
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            filt: list[str] = list(rec.filter.keys())
            if filt and filt != ["PASS"]:
                continue
            for alt in (rec.alts or []):
                out.add((rec.chrom, rec.pos, rec.ref, alt))
    return out


def jaccard(a: set[VariantKey], b: set[VariantKey]) -> float:
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def main() -> None:
    m2 = pass_variants(Path("results/ch01/mutect2_norm.vcf.gz"))
    s2_snv = pass_variants(Path("results/ch01/strelka_snvs_norm.vcf.gz"))
    s2_indel = pass_variants(Path("results/ch01/strelka_indels_norm.vcf.gz"))
    s2 = s2_snv | s2_indel
    print(f"Mutect2 PASS:    {len(m2)}")
    print(f"Strelka2 PASS:   {len(s2)}")
    print(f"Intersection:    {len(m2 & s2)}")
    print(f"Mutect2-only:    {len(m2 - s2)}")
    print(f"Strelka2-only:   {len(s2 - m2)}")
    print(f"Jaccard:         {jaccard(m2, s2):.3f}")
```

### Step 5 — pick three disagreements and explain (30 min)

For each disagreement (PASS in one caller, filtered or absent in the other), open both VCFs at the variant position and read:

- The allele frequency reported by each caller.
- The depth at the position.
- The supporting read count.
- The filter reason in the caller that did not PASS.

Three common reasons for disagreement:

- **Low AF (3-7%) in tumor.** Mutect2 may PASS it (the haplotype assembler builds confidence from the local read pattern); Strelka2 may filter it (the EVS threshold is calibrated for higher AF). The variant is likely a true low-AF somatic, harder to detect.
- **Indel in a homopolymer.** Strelka2's indel caller is conservative around homopolymers; Mutect2 may PASS it via the haplotype assembler. The "true" answer is hard to determine without orthogonal validation.
- **Variant in a repeat region.** Mutect2 may filter as `clustered_events` (the local pileup looks noisy); Strelka2 may PASS it. The "true" answer requires inspection of the repeat structure.

For each of your three disagreements, write a one-paragraph note explaining the most likely cause.

### Step 6 — write the Markdown report (20 min)

Combine the counts and the three disagreements into `cross_check.md`:

```markdown
# Challenge 1 — Mutect2 vs Strelka2 PASS-PASS Overlap

## Counts

| Set                 | Count |
|---------------------|------:|
| Mutect2 PASS        |   158 |
| Strelka2 PASS       |   149 |
| Intersection (both) |   137 |
| Mutect2-only        |    21 |
| Strelka2-only       |    12 |
| Union (either)      |   170 |
| Jaccard index       |  0.806|

## Disagreement examples

### Example 1: chr22:23456789 C>T

Mutect2: PASS, AF=0.06, DP=84, alt_count=5
Strelka2: filtered as LowEVS (EVS=2.3), AF=0.06, DP=82
Most likely cause: low AF puts it below Strelka2's EVS threshold; Mutect2's
haplotype assembler is more permissive on low-AF SNVs.

### Example 2: ...

### Example 3: ...
```

## Acceptance criteria

- [ ] `results/ch01/strelka_run/results/variants/somatic.snvs.vcf.gz` exists.
- [ ] `results/ch01/mutect2_norm.vcf.gz`, `results/ch01/strelka_snvs_norm.vcf.gz`, and `results/ch01/strelka_indels_norm.vcf.gz` exist.
- [ ] `results/ch01/cross_check.md` reports the Jaccard index, the four counts (Mutect2-PASS, Strelka2-PASS, intersection, union), and three worked disagreement examples.
- [ ] `results/ch01/run-info.json` records the Strelka2 version, the Mutect2 version, the normalization steps (bcftools version), and the Jaccard index.
- [ ] You can articulate, in two sentences, why one caller PASSed a variant the other filtered, for each of the three examples.

## Optional extensions

- **Run both callers on a real-data set.** The DREAM Somatic Mutation Calling Challenge simulated pairs (<https://www.synapse.org/#!Synapse:syn312572>) have ground truth; you can compute precision and recall for each caller separately, then for the union and intersection.
- **Test a third caller.** VarDict (<https://github.com/AstraZeneca-NGS/VarDict>) is free and open source; a three-caller union typically gains another 2-5% sensitivity for ~10% false-positive cost.
- **Compute the per-AF-bin agreement.** Bin variants by tumor AF (0-5%, 5-10%, 10-20%, 20-30%, 30%+); compute the Jaccard within each bin. Expect Jaccard to be ~0.5 at low AF and ~0.95 at high AF.

## A note on what this challenge does not show

A higher Jaccard index does not mean a more correct call set. Two callers that share a systematic bias (e.g. both trained on the same truth sets) will agree even when both are wrong. The cross-check is informative about *disagreement*, not about absolute correctness. To validate the calls, you need orthogonal evidence: a second sequencing platform, a different library prep, a clinical-grade panel, or direct experimental confirmation. The Week 11 didactic dataset has a known-truth VCF (`data/truth_chr22.vcf.gz`) for this purpose; the cross-check tells you "where the callers disagreed", and the truth-set comparison tells you "which of them was right".

## References

- Kim S, Scheffler K, Halpern AL, et al. *Strelka2: fast and accurate calling of germline and somatic variants.* **Nature Methods** 15:591-594 (2018). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>.
- Cibulskis K, Lawrence MS, Carter SL, et al. *Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples.* **Nature Biotechnology** 31:213-219 (2013). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/>.
- Strelka2 documentation: <https://github.com/Illumina/strelka>.
- Ewing AD, Houlahan KE, Hu Y, et al. *Combining tumor genome simulation with crowdsourcing to benchmark somatic single-nucleotide-variant detection.* **Nature Methods** 12:623-630 (2015). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509593/>.
