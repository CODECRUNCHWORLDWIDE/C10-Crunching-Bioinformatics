# Lecture 2 — Mutect2 and Strelka2: The Callers, the Filters, and the Pipeline

> **Educational and research use only.** This lecture describes the open-source somatic variant-calling pipelines we use in the exercises and the mini-project. The same tools, run in clinical laboratories with full validation, version locking, and quality oversight, produce clinical reports. Our pipeline does not have any of that scaffolding. Do not use Week 11 outputs to guide patient care.

## 1. Why two callers and not one

The pipeline we build runs **Mutect2** (Cibulskis et al. 2013, *Nature Biotechnology* 31:213; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/>) as the primary somatic SNV / indel caller, and optionally **Strelka2** (Kim et al. 2018, *Nature Methods* 15:591; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>) as an independent cross-check. The reason for two: every caller has its own model assumptions, its own filter heuristics, and its own systematic blind spots. A variant that PASSes in both callers carries higher confidence than a variant that PASSes in only one. The literature reports the Mutect2 / Strelka2 PASS-PASS overlap at AF >= 10% in the 80-95% range on real tumor-normal data, with the disagreement concentrated at low allele frequencies, in repetitive regions, and at small indels.

The pattern is widely used in cancer-genomics analysis pipelines: report variants that PASS in either caller for sensitivity, report variants that PASS in both for specificity, and quote both numbers. The pipeline this week uses Mutect2 as the default and treats Strelka2 as a Challenge 1 add-on so the core workflow stays small.

## 2. Inside Mutect2: the active-region haplotype caller

Mutect2 is not a simple position-by-position variant caller. It is an **active-region haplotype assembler** with a matched-normal Bayesian model on top. The high-level flow:

1. **Scan the BAM for active regions.** An active region is a stretch of the reference where the read pileup suggests at least one read disagrees with the reference at non-trivial frequency. Mutect2 walks the BAMs in chunks and identifies stretches that are worth assembling.
2. **Assemble local haplotypes.** For each active region, Mutect2 builds a small De Bruijn graph from the reads' k-mers and enumerates the haplotypes (sequences of bases that traverse the graph) that the reads support. This step is conceptually similar to short-read OLC assembly but local: each active region is ~100-300 bp.
3. **Realign reads to candidate haplotypes.** For each read overlapping the active region, Mutect2 computes the likelihood under each candidate haplotype using a pair-HMM (Hidden Markov Model) that accounts for sequencing-error rates per base.
4. **Compute the variant likelihood.** For each candidate variant supported by a haplotype, Mutect2 computes the likelihood under "somatic variant exists" versus "no variant exists" and emits a TLOD (tumor log-odds) score. With a matched normal, it also computes a NLOD (normal log-odds, the log-odds against the variant being in the normal).
5. **Emit unfiltered VCF.** The output is a VCF with one record per candidate variant, annotated with TLOD, NLOD, depth, allele depth, mapping quality, and the supporting haplotype information.

The active-region + haplotype-assembly approach lets Mutect2 call indels and complex variants reliably; a naive position-by-position caller would emit each base of a multi-base substitution as a separate (and inconsistent) call.

The full canonical call:

```bash
gatk Mutect2 \
  -R reference.fasta \
  -I tumor.bam \
  -I normal.bam \
  -tumor TUMOR_SAMPLE_NAME \
  -normal NORMAL_SAMPLE_NAME \
  --panel-of-normals pon.vcf.gz \
  --germline-resource gnomad.vcf.gz \
  -L intervals.list \
  -O unfiltered.vcf.gz
```

Every flag matters:

- `-R reference.fasta` — the same reference both BAMs were aligned to. Coordinates must agree.
- `-I tumor.bam` and `-I normal.bam` — passed as separate `-I` flags. Both have indexes (`.bai`).
- `-tumor` and `-normal` — the *sample names* (not file paths) as listed in each BAM's `@RG SM:` header. If these are reversed, Mutect2 produces inverted somatic calls; we cover the verification check below.
- `--panel-of-normals` — the PON VCF. Optional but strongly recommended.
- `--germline-resource` — the gnomAD allele-frequency VCF. Required for the FilterMutectCalls `germline` filter to work.
- `-L intervals.list` — restrict to specified intervals. For exome calling, pass the exome capture intervals; for whole-genome, omit or pass a chromosome list.
- `-O unfiltered.vcf.gz` — output VCF, BGZF-compressed.

Mutect2 also emits a `.tbi` index alongside the VCF and a `.stats` file that FilterMutectCalls consumes downstream.

## 3. Verifying the sample names before you call

The single most damaging bug in a Mutect2 pipeline is mis-tagging the tumor and normal BAMs — if Mutect2 gets the sample names backward, it produces calls against the inverse comparison: "variants in the normal that are not in the tumor", which is a list of variants where the tumor has *reverted* a germline variant to the reference (a rare event; the call set will be tiny). The pipeline silently produces a defensible-looking result that is wrong.

The protection: every Mutect2 wrapper script we write **reads the BAM @RG SM: header** with `samtools view -H` or `pysam.AlignmentFile.header`, extracts the sample name, and verifies it matches what the user passed on the command line.

```python
import subprocess
def get_sample_name(bam_path: Path) -> str:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam_path)],
        check=True, capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        if line.startswith("@RG") and "SM:" in line:
            for field in line.split("\t"):
                if field.startswith("SM:"):
                    return field[len("SM:"):]
    raise ValueError(f"No SM: tag found in {bam_path}")
```

Then the wrapper asserts that `get_sample_name(tumor_bam) == args.tumor_sample` and similarly for the normal. This catches mis-tags before Mutect2 runs.

## 4. Calling Mutect2 from Python with subprocess

The exercises use a thin wrapper around `subprocess.run(..., check=True, capture_output=True)`:

```python
def run_mutect2(
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    tumor_sample: str,
    normal_sample: str,
    pon_vcf: Path,
    germline_vcf: Path,
    intervals: Path,
    out_vcf: Path,
    threads: int = 4,
) -> None:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "gatk", "Mutect2",
        "-R", str(reference),
        "-I", str(tumor_bam),
        "-I", str(normal_bam),
        "-tumor", tumor_sample,
        "-normal", normal_sample,
        "--panel-of-normals", str(pon_vcf),
        "--germline-resource", str(germline_vcf),
        "-L", str(intervals),
        "-O", str(out_vcf),
        "--native-pair-hmm-threads", str(threads),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
```

Three patterns worth naming:

- `check=True` raises `CalledProcessError` if GATK exits non-zero; we catch and re-raise with a more informative message.
- `capture_output=True, text=True` collects stdout and stderr as strings (Mutect2 writes long progress logs to stderr). We log them to a file rather than printing to console.
- `--native-pair-hmm-threads` controls the pair-HMM parallelism inside Mutect2. Distinct from process-level threading.

The wrapper does not pass `--af-of-alleles-not-in-resource`, `--genotype-germline-sites`, `--mitochondria-mode`, or `--linked-de-bruijn-graph`; each is a specialized option and the defaults are appropriate for the didactic dataset. Always check the canonical GATK documentation for the production calls.

## 5. FilterMutectCalls and the filter set

The raw Mutect2 VCF contains every candidate variant whose TLOD exceeds a threshold (default: 6.3). Most of those will fail one or more quality filters. `FilterMutectCalls` is the GATK tool that applies the filter set and annotates each variant's FILTER column.

```bash
gatk FilterMutectCalls \
  -R reference.fasta \
  -V unfiltered.vcf.gz \
  --contamination-table contamination.table \
  --stats unfiltered.vcf.gz.stats \
  -O filtered.vcf.gz
```

The filter set is documented at <https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-FilterMutectCalls>. We summarized them in Lecture 1 §11; here is the practical reading:

- A clean PASS variant has FILTER `PASS`.
- A germline variant has FILTER `germline`; FilterMutectCalls believed it was inherited.
- A low-AF variant in a contaminated sample has FILTER `contamination`; the AF is consistent with the contamination estimate.
- A variant in a repeat region has FILTER `slippage` or `clustered_events`; the supporting reads may be slippage artefacts.
- A multi-filter variant has FILTER `germline;weak_evidence` or `slippage;strand_bias;base_qual`; each filter independently flagged it.

In our pipeline we tally each filter reason and report the counts. The mini-project's Markdown report includes a table:

```text
FILTER reason             count
PASS                      158
germline                   38
weak_evidence              12
panel_of_normals            9
clustered_events            6
strand_bias                 4
slippage                    3
contamination               1
(others)                    2
```

The pattern of filter reasons is informative. A run that flags 50% of candidate variants as `germline` means the matched normal is doing its job. A run that flags 50% as `weak_evidence` means the coverage was inadequate. A run that flags 20% as `contamination` is a major lab quality issue.

## 6. Running CalculateContamination

The contamination table is computed in two steps: `GetPileupSummaries` walks a list of common biallelic germline sites and counts the allele depth at each; `CalculateContamination` takes the pileups and computes the maximum-likelihood contamination fraction.

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

The `common_biallelic.vcf.gz` is a subset of gnomAD restricted to biallelic SNVs at population frequency 0.01-0.2 — common enough that most samples carry many, rare enough that the AF in the tumor is unaffected by tumor evolution. The Broad publishes a chr-distributed version (`small_exac_common_3.hg38.vcf.gz`) in the somatic-hg38 resource bundle.

The output `contamination.table` is a TSV:

```text
sample              contamination       error
TUMOR_SAMPLE        0.012               0.003
```

Pass this to FilterMutectCalls via `--contamination-table contamination.table`. Mutect2 will flag low-AF variants whose AF is consistent with the contamination as `contamination`.

The same flow works for tumor-only mode (omit `-matched`) but the estimate is less precise.

## 7. Inside Strelka2

Strelka2 (Kim et al. 2018) is Illumina's matched tumor-normal caller. It is a different statistical model from Mutect2:

- Strelka2 separates the SNV and indel calling pipelines; they emit separate VCFs.
- The Empirical Variant Score (EVS) is a Random Forest trained on labeled truth sets; Strelka2's filter set is a single EVS threshold (~5-10 depending on the variant type) rather than the multi-filter approach of FilterMutectCalls.
- Strelka2's somatic SNV model is a position-by-position Bayesian model; it does not do active-region haplotype assembly.
- The indel caller does small-window haplotype enumeration but is more conservative than Mutect2.

Strelka2 runs in two phases: a `configureStrelkaSomaticWorkflow.py` step that creates a working directory and a runtime configuration, then a `runWorkflow.py` step that does the actual calling.

```bash
configureStrelkaSomaticWorkflow.py \
  --normalBam normal.bam \
  --tumorBam tumor.bam \
  --referenceFasta reference.fasta \
  --runDir strelka_run

strelka_run/runWorkflow.py -m local -j 4
```

Output VCFs land in `strelka_run/results/variants/`:

- `somatic.snvs.vcf.gz` — SNV calls with EVS.
- `somatic.indels.vcf.gz` — indel calls with EVS.

PASS variants are those with FILTER = `PASS`; the FILTER column in Strelka2 is `PASS` if the EVS score exceeds the threshold and there are no Strelka-specific quality flags.

Challenge 1 walks through Strelka2 end-to-end and compares the PASS sets to Mutect2.

## 8. Comparing Mutect2 and Strelka2 PASS sets

The standard comparison is the **PASS-PASS overlap** at matched chromosome / position / REF / ALT. Both VCFs are normalized first (`bcftools norm -f reference.fasta -m -any input.vcf.gz`) to split multi-allelic sites and left-align indels. The comparison reduces to set arithmetic on `(chrom, pos, ref, alt)` tuples.

```python
from typing import Set
from pathlib import Path
import pysam

VariantKey = tuple[str, int, str, str]

def pass_variants(vcf_path: Path) -> set[VariantKey]:
    out: set[VariantKey] = set()
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            if rec.filter.keys() != ["PASS"] and "PASS" not in rec.filter.keys():
                continue
            for alt in (rec.alts or []):
                out.add((rec.chrom, rec.pos, rec.ref, alt))
    return out

m2 = pass_variants(Path("filtered.mutect2.vcf.gz"))
s2_snv = pass_variants(Path("strelka_run/results/variants/somatic.snvs.vcf.gz"))
s2_indel = pass_variants(Path("strelka_run/results/variants/somatic.indels.vcf.gz"))
s2 = s2_snv | s2_indel

union = m2 | s2
both = m2 & s2
mutect2_only = m2 - s2
strelka2_only = s2 - m2

print(f"Mutect2 PASS: {len(m2)}")
print(f"Strelka2 PASS: {len(s2)}")
print(f"Intersection (both): {len(both)}")
print(f"Union (either): {len(union)}")
print(f"Mutect2 only: {len(mutect2_only)}")
print(f"Strelka2 only: {len(strelka2_only)}")
print(f"Concordance (intersection / union): {len(both) / len(union):.3f}")
```

The standard reporting is the **Jaccard index** (intersection over union). On the simulated didactic pair we ship, the Jaccard is typically 0.85-0.92 for SNVs and 0.70-0.85 for indels (indels are harder to agree on). On real tumor-normal data the numbers are similar at high AF (>= 10%) and lower at low AF.

## 9. Reading a Mutect2 VCF in Python

The PASS VCF is the input to the mutational-signature decomposition (Lecture 3, Exercise 3) and to the OncoKB / CIViC lookup (Challenge 2). The standard parsing library is **pysam** (it wraps htslib) or **cyvcf2** (a slightly faster wrapper); we use pysam for consistency with Weeks 5 and 8.

```python
import pysam
from pathlib import Path

def parse_pass_variants(vcf_path: Path) -> list[dict]:
    rows: list[dict] = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            filters: list[str] = list(rec.filter.keys())
            if filters and filters != ["PASS"]:
                continue
            tumor_sample: str = list(rec.samples.keys())[0]
            sample_data = rec.samples[tumor_sample]
            ad = sample_data.get("AD")
            af = sample_data.get("AF")
            dp = sample_data.get("DP")
            for alt in (rec.alts or []):
                rows.append({
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": alt,
                    "qual": rec.qual,
                    "filter": "PASS",
                    "tumor_ad_ref": ad[0] if ad else None,
                    "tumor_ad_alt": ad[1] if ad and len(ad) > 1 else None,
                    "tumor_af": af[0] if isinstance(af, tuple) else af,
                    "tumor_dp": dp,
                })
    return rows
```

The `tumor_sample` name is the first sample in the VCF; Mutect2 orders the FORMAT samples as tumor-first by default. Allele depth (`AD`) is a tuple of (ref_count, alt_count); allele fraction (`AF`) is the alt fraction; depth (`DP`) is total reads at the site.

For each PASS SNV we will need (chrom, pos, ref, alt, trinucleotide context); see Lecture 3 §3.

## 10. Reading the FILTER tally

A useful smoke check after FilterMutectCalls is the histogram of filter reasons:

```python
from collections import Counter
import pysam
from pathlib import Path

def filter_histogram(vcf_path: Path) -> Counter:
    counts: Counter = Counter()
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            filters: list[str] = list(rec.filter.keys())
            if not filters:
                counts["PASS"] += 1
            elif filters == ["PASS"]:
                counts["PASS"] += 1
            else:
                for f in filters:
                    counts[f] += 1
    return counts
```

A healthy report has:

- PASS as the largest category.
- `germline` as the second-largest (the matched normal is doing its job).
- `weak_evidence`, `clustered_events`, `panel_of_normals` in the moderate-frequency range.
- `contamination` rare (less than 1-2% of total candidates).
- `strand_bias`, `base_qual`, `map_qual` very rare.

A skewed distribution is a flag for a quality issue: very high `weak_evidence` means low coverage; very high `contamination` means lab cross-contamination; very high `clustered_events` means the alignment was poor and Mutect2 is calling correlated noise.

## 11. Calling indels with Mutect2

Mutect2 calls SNVs and indels in the same pass; both end up in the same VCF. Indels show up as variants whose REF length and ALT length differ. The 96-class mutational-signature decomposition (Lecture 3) is SNV-only; indels are decomposed into a separate ID-class catalog by SigProfilerAssignment, with its own COSMIC v3.3 reference set.

A quick filter for SNV-only:

```python
def is_snv(rec) -> bool:
    return len(rec.ref) == 1 and all(len(a) == 1 for a in (rec.alts or []))
```

The mini-project pipeline separates SNVs and indels after FilterMutectCalls; the SNVs feed into SigProfilerAssignment SBS mode, the indels into SigProfilerAssignment ID mode (optional).

## 12. The `--native-pair-hmm-threads` knob

Mutect2's compute is dominated by the pair-HMM step (computing read likelihoods under each candidate haplotype). The `--native-pair-hmm-threads` flag controls the per-process thread count for that step; 4 is a reasonable default for a laptop, 16 for a workstation. Higher counts give diminishing returns past 8 because the active-region scanner becomes the bottleneck.

Mutect2 does not parallelize across chromosomes by default; for large input you split by interval (one Mutect2 call per chromosome) and merge the resulting VCFs with `MergeVcfs`. The didactic dataset is chr22 only, so this is unnecessary.

## 13. The GATK Java options

GATK4 is a Java application. The default JVM heap may be insufficient on real-data BAMs; pass `--java-options "-Xmx8g"` (or larger) for whole-genome runs. The Week 11 didactic pipeline runs in the default heap, but the `run_mutect2` wrapper supports an optional `java_options` parameter for production runs.

## 14. Quick troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Mutect2 found no variants` on a known-positive dataset | Sample names wrong (tumor and normal flipped) | Re-check `-tumor` / `-normal` against the BAM @RG SM tags |
| `Mutect2 found no variants` on an interval | Interval falls outside the BAM coverage | Verify `-L` against `samtools coverage` |
| Hundreds of `clustered_events` filters | Repetitive region; misaligned reads | Restrict to high-confidence intervals; consider a stricter PON |
| Very high `contamination` filter count | Lab cross-contamination | Re-run with `--contamination-table` from CalculateContamination |
| Mutect2 OOM (out of memory) | JVM heap too small | `--java-options "-Xmx16g"` |
| Strelka2 fails at configure step | BAM not indexed | `samtools index tumor.bam normal.bam` |
| Strelka2 fails at runWorkflow step | Reference not indexed | `samtools faidx reference.fasta && samtools dict reference.fasta > reference.dict` |
| FilterMutectCalls fails with "stats file not found" | The `.stats` Mutect2 emits was deleted | Re-run Mutect2 with `--stats` explicit and use it for FilterMutectCalls |

## 15. The reproducibility checklist for a Mutect2 + FilterMutectCalls pipeline

Before you publish a somatic variant call, the `run-info.json` records:

- GATK version (`gatk --version`).
- Strelka2 version (`configureStrelkaSomaticWorkflow.py --version`) if used.
- samtools and bcftools versions.
- Reference build (GRCh38; the genome FASTA URL or MD5).
- PON source and date (URL or local path + MD5).
- Germline-resource VCF source and date.
- Tumor and normal BAM paths and their `@RG SM:` sample names.
- Intervals file (or "whole genome" if none).
- Contamination estimate (the value from CalculateContamination, or "not run").
- FilterMutectCalls version and parameter set.
- Run date and the host.

Without this, a published call cannot be reproduced.

## 16. References

- Cibulskis K, Lawrence MS, Carter SL, et al. *Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples.* **Nature Biotechnology** 31:213-219 (2013). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3833702/>.
- Kim S, Scheffler K, Halpern AL, et al. *Strelka2: fast and accurate calling of germline and somatic variants.* **Nature Methods** 15:591-594 (2018). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6314977/>.
- Van der Auwera GA, O'Connor BD. *Genomics in the Cloud.* O'Reilly Media (2020). The canonical GATK4 reference.
- GATK Best Practices for somatic short-variant discovery: <https://gatk.broadinstitute.org/hc/en-us/articles/360035531132>.
- Mutect2 documentation: <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2>.
- FilterMutectCalls documentation: <https://gatk.broadinstitute.org/hc/en-us/articles/360036726891-FilterMutectCalls>.
- Strelka2 documentation: <https://github.com/Illumina/strelka>.

## 17. Self-check

You should be able to answer the following without looking back:

- What is an active region in Mutect2 and why does the haplotype-assembly approach matter for indel calling?
- What does the `-tumor` / `-normal` flag refer to and what is the canonical failure mode if it is misset?
- What goes wrong if `--germline-resource` and the reference BAM are on different builds (hg19 vs GRCh38)?
- What does `FilterMutectCalls --contamination-table` do that you cannot do with the default contamination filter?
- How is Strelka2's filter strategy different from Mutect2's?
- Define the PASS-PASS overlap and the Jaccard index. What is the typical value for AF >= 10% SNVs on real data?
- What is the difference between `--native-pair-hmm-threads` and JVM heap size, and which one matters for memory pressure?

If any of these are not crisp, re-read §2 (Mutect2 active-region model), §3 (sample-name verification), §5-§6 (filters and contamination), §7-§8 (Strelka2 and the PASS-PASS overlap) before moving on.

---

Continue to [Lecture 3 — Mutational signatures and the COSMIC / OncoKB / CIViC interpretation layer](./03-mutational-signatures-and-interpretation.md).
