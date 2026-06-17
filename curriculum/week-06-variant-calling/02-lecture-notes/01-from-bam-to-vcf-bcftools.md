# Lecture 1 — From BAM to VCF with bcftools

> **Duration:** ~3 hours of reading + paper-and-pencil probability + a brief Python sanity check.
> **Outcome:** You can describe the genotype-likelihood model in two paragraphs, write down the probability `P(reads | genotype)` for a diploid site by hand, run `bcftools mpileup | bcftools call -m` end to end on a sorted, indexed BAM, and read the resulting VCF column by column.

If you only remember one thing from this lecture, remember this:

> **Variant calling is a hypothesis test, not a counting exercise. For each reference position, the caller compares the likelihood of the observed reads under "this sample is homozygous reference" vs "this sample carries an alternative allele," and emits a variant only when the latter is favored by Bayes' rule. The likelihood model integrates depth (more reads is more evidence), base quality (a Q30 base is 1000x more informative than a Q3 base), and mapping quality (a MAPQ-60 read is 10^6 times more trustworthy than a MAPQ-0 read). `bcftools mpileup` builds the per-position read-allele matrix; `bcftools call -m` applies the likelihood model and emits a VCF. Everything else this week — filtering, annotation, comparison across callers — is downstream of that VCF.**

Week 5's BAM is the input. This week's VCF is the output. The transformation is the probability calculation in between.

---

## 1. The problem `bcftools` solves

A typical 30x whole-genome human BAM has ~3 · 10^9 reference positions, each covered by ~30 reads, each read with ~0.1% base-call error (Q30). The naive approach — emit a variant wherever any read disagrees with the reference — produces:

```
N_naive ≈ N_positions × coverage × error_rate
       ≈ 3·10^9 × 30 × 10^-3
       ≈ 10^8 candidate variants.
```

Of these, the real number is ~3 · 10^6 (3 million SNPs per human, the textbook estimate from the 1000 Genomes project). So the naive approach has a **false discovery rate of ~97%**. Variant calling is fundamentally a *filtering* problem: how do you tell the 3% of disagreements that reflect a true biological allele from the 97% that reflect sequencing or alignment error?

The answer is the **genotype-likelihood model**: for each position, ask the probabilistic question "given the reads I observe here, with their base qualities and mapping qualities, what is the posterior probability of each possible genotype?" The genotype with the highest posterior is the called genotype. The variant is emitted if the called genotype differs from homozygous reference. The QUAL value in the VCF is `-10 · log10(P[wrong call])`, on the same Phred scale as MAPQ in the BAM.

This is the same idea as MAPQ in Lecture 5 — convert a probability to a Phred score for easy comparison — but applied at the per-position level instead of the per-read level. And like MAPQ, the conversion is exact when the model assumptions hold (Bayesian inference, conditional independence of reads) and approximate when they do not (correlated errors at homopolymer runs, mapping-quality-zero clusters in repetitive regions). Hard filters in §6 are the patches we apply where the model breaks down.

---

## 2. The genotype-likelihood model in three concepts

The probability calculation behind `bcftools call -m` is the simplest non-trivial Bayesian inference you will see this semester. You do not need to implement it to use `bcftools` — but a working mental model of the three components is the difference between "the caller is magic" and "I can predict what the caller will and will not get right."

### 2.1 The likelihood: `P(reads | genotype)`

For a diploid sample at a position with two possible alleles `A` (reference) and `B` (alternative), the three possible genotypes are:

- `AA` — homozygous reference.
- `AB` — heterozygous.
- `BB` — homozygous alternative.

For each genotype, the expected allele frequency of `B` in the reads is:

```
f(AA) = 0.0    (every read should be A)
f(AB) = 0.5    (each read is A or B with 50/50 odds)
f(BB) = 1.0    (every read should be B)
```

Given `n` reads observed at the position with `k` of them carrying the `B` allele, the binomial probability of that observation under each genotype is:

```
P(k of n | AA) = C(n, k) · 0.0^k · 1.0^(n-k)           (= 0 if k > 0, else 1)
P(k of n | AB) = C(n, k) · 0.5^k · 0.5^(n-k)           = C(n, k) · 0.5^n
P(k of n | BB) = C(n, k) · 1.0^k · 0.0^(n-k)           (= 0 if k < n, else 1)
```

The pure-binomial form above breaks if you ever observe a single `B` read in an `AA` site or a single `A` read in a `BB` site (the likelihood becomes 0, which is nonsense — a Q30 sequencing error happens about 1 in 1000 bases). To fix this, the model substitutes each read's likelihood with one that accounts for the base-call error probability `ε = 10^(-Q/10)`:

```
P(read = A | genotype with allele frequency f) = (1 - f) · (1 - ε) + f · ε
P(read = B | genotype with allele frequency f) = f · (1 - ε) + (1 - f) · ε
```

So a Q30 base (ε = 0.001) seen as `B` at an `AA` site contributes a likelihood of `0.001`, not `0` — and accumulating 30 such reads multiplies down to `10^-90`, low but finite. The model accommodates real-world error rates without collapsing to zero on the first sequencing miscall.

For haploid samples (bacteria, mitochondria), there are only two genotypes: `A` (reference) and `B` (alternative), with `f(A) = 0` and `f(B) = 1`. The `--ploidy 1` flag in `bcftools call` switches to this two-genotype model.

### 2.2 The prior: `P(genotype)`

The genotype-likelihood model converts to a posterior via Bayes:

```
P(genotype | reads) ∝ P(reads | genotype) · P(genotype)
```

The prior `P(genotype)` comes from population genetics. For a known variant in a population, you might use Hardy-Weinberg equilibrium with the allele frequency `p` from the population:

```
P(AA) = (1 - p)^2
P(AB) = 2 · p · (1 - p)
P(BB) = p^2
```

For a *new* variant (most variants in a single-sample WGS run are not in dbSNP), the prior is harder to pin. `bcftools call -m` uses a default per-site mutation rate of `θ = 1.1 · 10^-3` (the human nucleotide diversity from the literature) and computes the prior from a coalescent model. You can override with `-P θ` if you have a better estimate for your species (much higher for many bacteria, much lower for many plants).

The prior matters less than you might think. At positions with ≥ 10 reads, the likelihood is sharp enough that the prior is essentially decorative — the posterior is dominated by the data. The prior matters most at low-coverage positions where the data is ambiguous.

### 2.3 The posterior: pick the maximum

After Bayes, you have:

```
P(AA | reads), P(AB | reads), P(BB | reads)
```

The caller emits the genotype with the maximum posterior. The `QUAL` value in the VCF is `-10 · log10(1 - max_posterior)` — Phred-scaled probability that the call is wrong. The `PL` field in `FORMAT` is the full vector of Phred-scaled likelihoods for each genotype (with the maximum-likelihood genotype set to 0 and the others as differences from that).

A concrete example. Suppose at position 100 you observe 30 reads, 14 carrying allele A and 16 carrying allele B, all at Q30 (ε = 0.001). The likelihoods are:

```
log10 P(reads | AA) ≈ 14 · log10(0.999) + 16 · log10(0.001)
                  ≈ 14 · (-0.0004) + 16 · (-3)
                  ≈ -48     (essentially zero probability)
log10 P(reads | AB) ≈ 30 · log10(0.5)
                  ≈ -9.0
log10 P(reads | BB) ≈ 14 · log10(0.001) + 16 · log10(0.999)
                  ≈ 14 · (-3) + 16 · (-0.0004)
                  ≈ -42     (essentially zero probability)
```

The heterozygous genotype `AB` is overwhelmingly favored. Apply a flat prior (1/3 each), and the posterior `P(AB | reads)` is essentially 1.0, so `QUAL` is at the cap (~99). The `PL` field is `(390, 0, 330)` — `AA` is 39 log units below `AB`, `BB` is 33 log units below.

This is the calculation `bcftools call -m` does at every variant position. It is straightforward; the engineering question is making it fast enough to run on 3 billion positions in a couple of hours.

---

## 3. `bcftools mpileup` — building the per-position read-allele matrix

`bcftools mpileup` is the data-collection step. It walks the BAM column by column (the same `pileup` iteration pattern as `pysam.AlignmentFile.pileup`), and for each position emits a record containing:

- The reference base at this position (from the FASTA).
- The list of (read, base, base-quality, mapping-quality) tuples observed at this position.
- Various per-position summary stats (total depth, allele depths, strand counts).

The output is binary VCF (BCF) with one record per *covered* position (positions with zero reads are skipped). The records do not yet carry variant calls — they carry the *data* the caller will use.

```bash
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/SRR1770413.markdup.bam \
    -o calls/SRR1770413.mpileup.bcf
```

Key flags:

- **`-f ref/ecoli.fa`** — the reference FASTA. Required; `mpileup` cannot proceed without it (it needs to know the reference base at every position).
- **`-Ou`** — output uncompressed BCF. Use this when piping to `bcftools call`; the next step will re-compress. Alternatives: `-Oz` (bgzipped VCF), `-Ov` (uncompressed VCF, human-readable), `-Ob` (compressed BCF).
- **`--max-depth 250`** — cap the number of reads considered at any position. The default is 250; raise for high-coverage targeted sequencing. The cap exists because very-high-depth positions (1000x+) are usually in repetitive regions where the calls are unreliable anyway.
- **`-a 'AD,DP,SP'`** — extra per-record annotations. `AD` = allele depths (the (k, n-k) numbers from §2.1), `DP` = total depth, `SP` = Phred-scaled strand-bias p-value. All three are essential for downstream filtering; do not omit.
- **`-r chr:start-end`** — restrict to a region. Useful for testing on a slice before running on the whole genome.

The output BCF is the input to `bcftools call`.

### 3.1 What `mpileup` does at each position

For each position covered by at least one read, `mpileup` does:

1. Look up the reference base from `ref/ecoli.fa` (requires `samtools faidx` to have been run on the FASTA — see Week 5).
2. Iterate over reads spanning this position (via the pysam pileup machinery in `htslib`).
3. For each read, extract the base at this position (using the CIGAR to map the query offset to the reference offset), the base quality, and the read's MAPQ.
4. Filter reads by MAPQ threshold (default 0; `-q 20` is a common stricter setting) and base quality threshold (default 13; `-Q 20` is stricter).
5. Filter reads by FLAG: drop secondary alignments (`0x100`), supplementary alignments (`0x800`), QC-failed (`0x200`), and duplicates (`0x400`) by default. (This is why marking duplicates in Week 5 matters: an unmarked-duplicate BAM passed to `mpileup` will count the same physical fragment multiple times, inflating allele depths and producing false-positive variants.)
6. Tally per-allele depths, strand counts, and quality sums.
7. Compute the per-allele Phred-scaled genotype likelihoods (the `PL` field).
8. Emit one BCF record with all of the above.

The output is several gigabytes for a human genome, ~30 MB for the *E. coli* mini-project.

---

## 4. `bcftools call` — applying the genotype model

`bcftools call` is the modeling step. It takes the BCF from `mpileup` and applies the genotype-likelihood model from §2 to produce a VCF with one record per *variant* position (positions with `GT = 0/0` for all samples are dropped if `-v` is set).

```bash
bcftools call -m -v --ploidy 1 -Oz -o calls/SRR1770413.vcf.gz \
    calls/SRR1770413.mpileup.bcf \
&& bcftools index -t calls/SRR1770413.vcf.gz
```

Key flags:

- **`-m`** — the multiallelic-and-rare-variant calling model (Li 2011 + Danecek 2014). This is the modern default; the legacy `-c` (consensus) model is from 2011 and should not be used for new work.
- **`-v`** — output variant sites only. Without `-v`, every position covered by a read is emitted (useful for some downstream tools like population genotyping; useless for routine variant discovery).
- **`--ploidy 1`** — haploid model for *E. coli*. **Always set this for bacteria.** The default is 2 (human diploid).
- **`-Oz`** — output bgzipped VCF, suitable for `tabix` indexing. Alternatives: `-Ov` (uncompressed VCF), `-Ob` (compressed BCF, smaller and faster to read with pysam).
- **`-o calls/SRR1770413.vcf.gz`** — output file. Defaults to stdout if omitted.

The two-step `mpileup | call` is the canonical pattern; the steps are commonly piped together:

```bash
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/SRR1770413.markdup.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/SRR1770413.vcf.gz - \
&& bcftools index -t calls/SRR1770413.vcf.gz
```

For an *E. coli* mini-project BAM (4.6 Mb, ~50x coverage), this takes ~30 seconds and emits ~150-200 variants. For a human germline BAM (3 Gb, 30x), ~2-4 hours and ~3 million variants. The two-step pipeline can be parallelized over chromosomes with GNU `parallel`; each chromosome takes ~10 minutes on a human genome.

---

## 5. The VCF format, column by column

The output of `bcftools call` is a Variant Call Format (VCF) file. Like SAM, VCF is plain text with tab-separated columns and a header. The spec is Danecek et al. 2011 + the VCFv4.3 specification at <https://samtools.github.io/hts-specs/VCFv4.3.pdf>.

### 5.1 Header lines

VCF headers start with `##` for metadata and a single `#` for the column header. A typical bcftools-produced header:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=NC_000913.3,length=4641652>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##INFO=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##bcftools_callCommand=call -m -v --ploidy 1 -Oz -o calls/SRR1770413.vcf.gz; Date=2026-05-13
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SRR1770413
```

The metadata lines (`##`) define every `INFO` and `FORMAT` field that may appear in the records. The column header line (`#`) is the start of the table. Sample names are columns 10+; in this example there is one sample, `SRR1770413`.

### 5.2 Variant records — the eight mandatory columns

A typical bcftools-emitted SNP record:

```
NC_000913.3  150123  .  A  G  227.0  .  DP=42;MQ=60;SP=2;AF=1  GT:PL:DP:AD  1:255,0:42:0,42
```

Column by column:

- **`CHROM = NC_000913.3`** — the contig from the BAM `@SQ` headers.
- **`POS = 150123`** — 1-based position of the variant on the reference.
- **`ID = .`** — no known variant ID (we did not annotate from dbSNP; `.` means "missing").
- **`REF = A`** — the reference base at this position.
- **`ALT = G`** — the alternative allele. Multiple alts are comma-separated (`G,T`).
- **`QUAL = 227.0`** — Phred-scaled probability the variant call is wrong. 227 ≈ 10^-22.7, essentially certain.
- **`FILTER = .`** — not yet filtered (`bcftools call` does not apply hard filters; that is the next step). After `bcftools filter`, this will be `PASS` or a filter name.
- **`INFO = DP=42;MQ=60;SP=2;AF=1`** — semicolon-separated `KEY=VALUE` pairs.
  - `DP=42` — 42 reads support this position.
  - `MQ=60` — average mapping quality of supporting reads is 60 (BWA-MEM cap; high confidence).
  - `SP=2` — Phred-scaled strand-bias p-value is 2 (a *low* number means *low* bias; the alt allele is seen on both strands).
  - `AF=1` — allele frequency is 1.0 (haploid sample, the call is homozygous-alt).
- **`FORMAT = GT:PL:DP:AD`** — the colon-separated keys describing the per-sample columns.
- **`SRR1770413 = 1:255,0:42:0,42`** — the per-sample values:
  - `GT=1` — genotype is `1` (the alt allele; haploid).
  - `PL=255,0` — Phred-scaled likelihoods: 255 for `0` (homozygous ref), 0 for `1` (homozygous alt). The alt allele is 255 log units more likely than the ref allele.
  - `DP=42` — 42 reads at this position in this sample.
  - `AD=0,42` — 0 reads carrying the ref allele, 42 reads carrying the alt allele.

For a diploid sample, the `GT` would be `0/0`, `0/1`, or `1/1`, and the `PL` would have three values (one per genotype).

### 5.3 An indel record

A 1-base insertion:

```
NC_000913.3  450200  .  T  TC  198.0  .  INDEL;DP=38;MQ=60;SP=5  GT:PL:DP:AD  1:255,0:38:0,38
```

Note that the `REF` is `T` (the anchor base, one bp before the inserted base, included for unambiguity) and the `ALT` is `TC` (the anchor plus the inserted base). The `INFO` field carries `INDEL` to flag that this is not a SNP.

A 2-base deletion:

```
NC_000913.3  600789  .  TCG  T  256.0  .  INDEL;DP=45;MQ=60;SP=1  GT:PL:DP:AD  1:255,0:45:0,45
```

`REF = TCG`, `ALT = T` (the anchor plus what is *kept*). The deleted bases are the 2 in `REF` after the anchor.

Indel representation has edge cases (left-alignment, normalization, multiallelic splits) that `bcftools norm` standardizes; see Lecture 2 §3 for the details.

---

## 6. A worked end-to-end pipeline

Here is the canonical short-variant calling pipeline from a sorted, indexed, markdup BAM to a filtered VCF. You will run it dozens of times this semester; memorize the steps.

### 6.1 Step 1 — Index the reference (one-time)

```bash
samtools faidx ref/ecoli.fa
# Produces ref/ecoli.fa.fai.
```

If you ran Week 5's pipeline, this is already done. `bcftools mpileup` requires the `.fai`; it cannot proceed without random-access into the reference FASTA.

### 6.2 Step 2 — Call variants

```bash
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/SRR1770413.markdup.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/SRR1770413.raw.vcf.gz - \
&& bcftools index -t calls/SRR1770413.raw.vcf.gz
```

For *E. coli* at 50x coverage, ~30 seconds. The output `.vcf.gz` is ~30 KB; the `.tbi` index is ~4 KB.

### 6.3 Step 3 — Hard-filter the raw calls

```bash
bcftools filter -Oz -o calls/SRR1770413.filtered.vcf.gz \
    -s LowQual \
    -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' \
    calls/SRR1770413.raw.vcf.gz \
&& bcftools index -t calls/SRR1770413.filtered.vcf.gz
```

Variants matching the expression get `FILTER=LowQual`; the rest get `FILTER=PASS`. The expression above is the simplified bacterial recipe from `resources.md`; the GATK Best Practices use a more elaborate expression for human germline data (see Lecture 2 §4).

The `-s` (filter tag name) and `-e` (exclude expression) flags together apply the named filter. The expression syntax is `INFO/FIELD` for `INFO` fields, `FORMAT/FIELD` for per-sample fields, and standard arithmetic + boolean operators.

### 6.4 Step 4 — Normalize indel representation

```bash
bcftools norm -f ref/ecoli.fa -Oz -o calls/SRR1770413.norm.vcf.gz \
    calls/SRR1770413.filtered.vcf.gz \
&& bcftools index -t calls/SRR1770413.norm.vcf.gz
```

`bcftools norm` does two things: **left-align** indels (an indel `AAA → AA` can equivalently be `AA → A` shifted right by one; the convention is to shift left) and **split multiallelic** variants (`A → G,T` becomes two records, `A → G` and `A → T`). Both are required by most downstream tools (VEP, the comparison harness, gnomAD intersection queries).

### 6.5 Step 5 — Sanity check with `bcftools stats`

```bash
bcftools stats calls/SRR1770413.norm.vcf.gz > results/vcf_stats.txt
```

The output includes the total count by variant type (SNP, indel), the Ts/Tv ratio, the per-FILTER counts, and the allele-frequency distribution. For an *E. coli* mini-project the expected output:

```
SN  0  number of SNPs:    168
SN  0  number of indels:   22
SN  0  number of others:    0
SN  0  ts/tv:            0.94
```

Ts/Tv around 1.0 is the expected ratio for random bacterial mutation. Human germline is ~2.0 (transitions are biochemically more common); deviation from the species-typical value is a QC signal that something is wrong.

---

## 7. A worked sanity check in Python

The full pipeline runs from Python via `subprocess`, the same pattern as Week 5 Exercise 1. Here is a minimal sanity check:

```python
from __future__ import annotations
import subprocess
from pathlib import Path

import pysam

ref = Path("ref/ecoli.fa")
bam = Path("aln/SRR1770413.markdup.bam")
vcf = Path("calls/SRR1770413.raw.vcf.gz")
vcf.parent.mkdir(parents=True, exist_ok=True)

# Step 1: index reference (if not already done).
if not Path(f"{ref}.fai").exists():
    subprocess.run(["samtools", "faidx", str(ref)], check=True)

# Step 2: call variants. Use shell=True for the pipe.
cmd = (
    f"bcftools mpileup -Ou -f {ref} --max-depth 250 -a 'AD,DP,SP' {bam} "
    f"| bcftools call -m -v --ploidy 1 -Oz -o {vcf} -"
)
subprocess.run(cmd, shell=True, check=True)

# Step 3: index the VCF.
subprocess.run(["bcftools", "index", "-t", str(vcf)], check=True)

# Step 4: sanity-check with pysam.VariantFile.
vf = pysam.VariantFile(str(vcf))
n_total = 0
n_snp = 0
n_indel = 0
for rec in vf:
    n_total += 1
    if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts):
        n_snp += 1
    else:
        n_indel += 1

print(f"Total variants: {n_total}")
print(f"SNPs:           {n_snp}")
print(f"Indels:         {n_indel}")
```

Expected output for SRR1770413 against NC_000913.3:

```
Total variants: 190
SNPs:           168
Indels:         22
```

If your numbers are wildly different (total < 50 or > 1000), the typical causes are: (a) ploidy not set (every variant becomes heterozygous and many false positives slip through), (b) duplicates not marked (PCR-duplicate stacks produce false-positive variants), (c) reference mismatch (the BAM was aligned to a different reference version than the FASTA passed to `mpileup`).

---

## 8. Common misconceptions

A short list of "things that seem right but are not":

- **"More coverage is always better for variant calling."** Up to a point. Below 10x, the likelihood model has too few observations and the variant calls are unreliable (false negatives at low coverage). Above 100x, additional coverage produces diminishing returns and can mask PCR-duplicate-driven artifacts. For 1000x panel sequencing, you genuinely need that depth; for routine WGS, 30x is the sweet spot.
- **"A variant with QUAL = 30 is good enough."** Sometimes. For a well-behaved human germline call set, `QUAL >= 30` is the standard threshold. For a bacterial single-sample call set, `QUAL >= 50` is more conservative. For somatic variant calling (tumor), `QUAL >= 100` is the typical floor — and the per-variant QUAL is not the right metric anyway (somatic callers use allele-fraction-aware metrics).
- **"`bcftools call` and GATK `HaplotypeCaller` produce the same calls."** They agree on ~95% of SNPs and ~85% of indels. The disagreements are at near-homopolymer indels (different alignment-error models), in low-coverage regions (different priors), and at multi-allelic sites (different ploidy assumptions). Challenge 1 quantifies this on a real BAM.
- **"The `FILTER=PASS` variants are the real ones; everything else is junk."** Mostly true, but `PASS` is *only the answer for the filters you applied*. Apply a stricter filter set and your `PASS` set shrinks. Apply a looser one and it grows. The hard-filter thresholds in the Best Practices are a compromise; for clinical sequencing you may want stricter, for population genetics you may want looser. The point of a `FILTER` column is that the filtering is *traceable*, not that any specific threshold is universally right.
- **"`bcftools call --ploidy 1` only matters for bacteria."** It also matters for human mitochondrial DNA (ploidy is 1 — there are many mitochondria per cell but they replicate clonally), for plant chloroplasts (same reason), and for cancer cell lines with known copy-number changes. The flag accepts a per-region file via `--ploidy-file` for these cases.

If any of these surprised you, re-read sections 2, 4, and 6.

---

## 9. Where this lecture lands you for Lecture 2

After this lecture you should be able to:

- Describe the genotype-likelihood model in two paragraphs (binomial likelihood, Bayes' rule, MAP estimate).
- Run `bcftools mpileup -f ref.fa aln.bam | bcftools call -mv --ploidy 1 -Oz -o calls.vcf.gz -` end to end.
- Read a VCF record column by column and decode every `INFO` and `FORMAT` field.
- Distinguish the *likelihood* (purely data-driven) from the *posterior* (likelihood × prior).
- Name the three reasons a variant might be wrong: low coverage, mapping artifacts, sequencing errors.

Lecture 2 takes the raw VCF and applies the GATK Best Practices hard filters, then annotates the survivors with VEP to convert "position 150,123 has a `G` allele" into "position 150,123 is a missense variant in gene `rpoB` causing `D516G`." By the end of Lecture 2 you will be able to read a fully-annotated VCF as "here are the 152 PASS SNPs in this sample, 87 of which are coding-region missense, of which 12 hit known antibiotic-resistance loci."

---

## Self-check questions

Before you move on, answer these without looking. If you cannot answer one, re-read the relevant section.

1. State the three components of the Bayesian inference behind `bcftools call -m`. (§2)
2. For a haploid sample with 20 reads of which 18 carry the alt allele at Q30, write down (or sketch) the per-genotype log-likelihoods. Which genotype wins, and by how many log units? (§2.1)
3. What does `--max-depth 250` do in `bcftools mpileup`? Why is the cap there? (§3)
4. What does the `INFO/SP` field mean? At what value would you start to worry? (§5.2)
5. Decode the FORMAT/sample column `GT:PL:DP:AD = 1:255,0:42:0,42`. What is the genotype, the depth, and the per-allele depth? (§5.2)
6. Why does an indel record have `REF=T` and `ALT=TC` rather than `REF=` and `ALT=C`? (§5.3)
7. Why is `--ploidy 1` required for *E. coli* and not for human germline data? (§4, §8)
8. Under what biological scenario would the likelihood model in `bcftools call -m` produce a false-positive variant call? Name two. (§8)
9. The QUAL field is on what scale? What does QUAL = 60 mean numerically? (§5.2)
10. What is the difference between `bcftools call -v` and `bcftools call` without `-v`? When would you skip `-v`? (§4)

Answers are not provided. If you struggle, the answers are in the section references; do the work.

---

## Further reading

- Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics* 27(21):2987–2993.
- Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- DePristo, M. et al. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics* 43:491–498.
- Danecek, P. et al. (2011). The variant call format and VCFtools. *Bioinformatics* 27(15):2156–2158.
- The VCFv4.3 specification: <https://samtools.github.io/hts-specs/VCFv4.3.pdf>.
- The bcftools source code: <https://github.com/samtools/bcftools>.
- The bcftools howtos: <https://samtools.github.io/bcftools/howtos/>.

---

*Continue to [Lecture 2 — GATK Best Practices and hard filters](./02-gatk-best-practices-and-hard-filters.md) once you have answered the self-check questions.*
