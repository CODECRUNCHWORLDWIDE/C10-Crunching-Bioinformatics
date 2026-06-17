# Challenge 1 — Compare bcftools and GATK on the same BAM

> **Estimated time:** 90 minutes.
> **Goal:** Run `bcftools call -m` and GATK `HaplotypeCaller` on the same sorted, duplicate-marked BAM, normalize both VCFs, compute the per-variant intersection and private sets with `bcftools isec`, and explain at least two specific variants where the callers disagree. Quantify the agreement at the SNP and indel level. The point is not to declare one caller the winner. The point is to map the corners where the two production callers diverge — and to defend that map in writing.

Lecture 2 §2 introduced the bcftools-vs-GATK trade-off as "bcftools is faster and simpler; GATK is more accurate at indels because of local reassembly." This challenge has you measure that trade-off on a real BAM and produce a per-variant comparison you can hand to a reviewer.

## Background — the algorithms in one paragraph each

**`bcftools call -m`** walks the BAM column by column. At each position, it builds the per-allele depth (`AD`) and per-base-quality matrix, applies the Bayesian likelihood model from Lecture 1 §2 (binomial-with-error-rate, per-genotype posterior, MAP estimate), and emits a variant if the MAP genotype differs from homozygous-reference. It does not re-align reads. It does not consider haplotype context. It is fast (~30 sec on the *E. coli* mini-project BAM) and accurate at SNPs in well-aligned regions.

**GATK `HaplotypeCaller`** identifies *active regions* of high mismatch/indel/soft-clip density, builds a De Bruijn graph from the reads in each active region, enumerates candidate haplotypes as paths through the graph, re-aligns the reads against each candidate haplotype using a pair-Hidden-Markov-Model, and emits the variants implied by the most likely haplotype pair (for diploid; one haplotype for haploid). The local-reassembly step means that **a true indel that was scattered across the reads as soft-clips is correctly called as an indel**, whereas `bcftools` may emit a noisy cluster of SNPs at the indel breakpoint or miss the indel entirely.

## Task

Run both callers on the Week 5 mini-project BAM (or, if that is too large, on a tractable slice), normalize, intersect, and write up.

### 1. The input BAM

Use the sorted+indexed+markdup BAM from the **Week 5 mini-project** (`SRR1770413.markdup.bam`). This is ~750k read pairs against `NC_000913.3` at ~49x mean coverage. If you did not complete the Week 5 mini-project, you can substitute the Exercise 1 BAM from Week 5 (`lambda.sorted.bam`) — the comparison is the same in shape, just smaller.

### 2. Run both callers

#### bcftools

```bash
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/SRR1770413.markdup.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/SRR1770413.bcftools.vcf.gz - \
&& bcftools index -t calls/SRR1770413.bcftools.vcf.gz
```

Expected: ~150-200 variants total, ~30 seconds wall-clock.

#### GATK

```bash
gatk HaplotypeCaller \
    -R ref/ecoli.fa \
    -I aln/SRR1770413.markdup.bam \
    -O calls/SRR1770413.gatk.vcf.gz \
    --sample-ploidy 1
bcftools index -t calls/SRR1770413.gatk.vcf.gz
```

Expected: ~150-200 variants total, ~5 minutes wall-clock.

If GATK is not installed and you do not want to install ~1 GB of Java, you can substitute a second `bcftools` run with different parameters (e.g. stricter `-q` and `-Q` thresholds in mpileup) — the comparison machinery is the same. The point is to compare two callers, not specifically to compare bcftools to GATK.

### 3. Normalize both VCFs

```bash
bcftools norm -f ref/ecoli.fa -m -any -Oz \
    -o calls/SRR1770413.bcftools.norm.vcf.gz \
    calls/SRR1770413.bcftools.vcf.gz
bcftools index -t calls/SRR1770413.bcftools.norm.vcf.gz

bcftools norm -f ref/ecoli.fa -m -any -Oz \
    -o calls/SRR1770413.gatk.norm.vcf.gz \
    calls/SRR1770413.gatk.vcf.gz
bcftools index -t calls/SRR1770413.gatk.norm.vcf.gz
```

This is **non-optional**. Without `-m -any` and `-f ref.fa`, indels in different equivalent representations will look like distinct variants in the intersection step.

### 4. Intersect and tabulate

```bash
bcftools isec -p compare/ -Oz \
    calls/SRR1770413.bcftools.norm.vcf.gz \
    calls/SRR1770413.gatk.norm.vcf.gz
```

This writes four VCFs to `compare/`:

- `compare/0000.vcf.gz` — variants private to bcftools (only in bcftools).
- `compare/0001.vcf.gz` — variants private to GATK (only in GATK).
- `compare/0002.vcf.gz` — variants in both (from the bcftools side, with bcftools annotations).
- `compare/0003.vcf.gz` — variants in both (from the GATK side, with GATK annotations).

Count each:

```bash
for f in compare/000{0,1,2,3}.vcf.gz; do
    echo "$f: $(bcftools view -H "$f" | wc -l)"
done
```

Split each set by variant type (SNP vs indel) using `bcftools view -v snps` / `-v indels`. Report a 2x2 table:

```
              | Both callers | bcftools only | GATK only | Total
--------------|-------------:|--------------:|----------:|------:
SNPs          |              |               |           |
Indels        |              |               |           |
Total         |              |               |           |
```

### 5. Inspect specific disagreements

Pick at least **two** specific variants where the callers disagree. For each, run:

```bash
# Look at the underlying alignment.
samtools tview -d T aln/SRR1770413.markdup.bam ref/ecoli.fa \
    -p NC_000913.3:<position>

# Look at the raw read pileup.
samtools mpileup -f ref/ecoli.fa -r NC_000913.3:<position>-<position> \
    aln/SRR1770413.markdup.bam
```

For each disagreement, in `challenges/notes/caller-comparison.md`:

- Print the position, REF, ALT, and the bcftools and GATK QUAL/DP/MQ/whatever-they-report.
- Describe the underlying alignment in one sentence (clean? soft-clipped? near a homopolymer? near a repetitive region?).
- Hypothesize which caller is right (or whether the position is ambiguous). Cite the algorithmic reason — "GATK reassembled and got the indel; bcftools saw scattered mismatches and called nothing" is a typical bcftools-private-indel diagnosis; "bcftools called a SNP at a position where GATK reassembled and saw a 1-bp deletion, both consistent with the data" is a typical SNP/indel-disagreement at near-homopolymer positions.

### 6. Write up

In `challenges/notes/caller-comparison.md` (300-500 words), answer:

- What was the total variant count from each caller? What was the intersection?
- What was the SNP agreement rate? The indel agreement rate? Compare to the textbook expectation (~95% SNPs, ~85% indels).
- For at least **two** specific variants in the disagreement set, paste the bcftools and GATK records side by side. Explain the alignment context and the algorithmic reason for the disagreement.
- For your dataset, which caller would you trust more if you had to pick one? Why? Name one downstream use case (clinical reporting, population genetics, antibiotic-resistance screening) and say which caller you would use and why.
- In one sentence, what would change in your answer if the dataset were a 30x human germline WGS instead of a 50x *E. coli* WGS?

## Acceptance criteria

- `python challenges/compare_callers.py <bam>` (or a `bash compare.sh` equivalent) runs without crashing.
- Both callers run end to end on the input BAM and produce indexed VCFs.
- Both VCFs are normalized with `bcftools norm -f ref.fa -m -any`.
- `compare/000{0,1,2,3}.vcf.gz` files exist with non-empty contents.
- A comparison table (3 rows: SNPs, Indels, Total; 4 columns: Both, bcftools-only, GATK-only, Total) is in `notes/caller-comparison.md`.
- At least two specific disagreement examples are documented with bcftools and GATK records and an algorithmic interpretation.
- The writeup is 300-500 words and addresses all five prompts.
- The script is < 250 lines (it should orchestrate, not re-implement the callers).

## Hints (do not peek for at least 20 minutes)

<details>
<summary>Hint 1 — How do I run GATK without installing 1 GB of Java?</summary>

If GATK is not available, do the comparison between two bcftools runs at different stringency:

```bash
# Looser: defaults.
bcftools mpileup -Ou -f ref.fa --max-depth 250 aln.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/a.vcf.gz -

# Stricter: -q 20 (MAPQ filter), -Q 20 (base-quality filter).
bcftools mpileup -Ou -f ref.fa --max-depth 250 -q 20 -Q 20 aln.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/b.vcf.gz -
```

The comparison is the same in shape; both call sets will have different counts and the disagreements will localize at low-quality reads and low-mapping-quality positions. Not as instructive as a true bcftools-vs-GATK comparison, but a valid substitute.

</details>

<details>
<summary>Hint 2 — How do I pick "interesting" disagreements?</summary>

After running `bcftools isec -p compare/`, look at:

```bash
# bcftools-private indels (most common GATK-bcftools disagreement category).
bcftools view -v indels compare/0000.vcf.gz | head

# GATK-private SNPs near indels (the second-most-common disagreement).
bcftools view -v snps compare/0001.vcf.gz | head
```

Pick a position from each, run `samtools tview` to see the alignment, and you have your two case studies. Avoid positions in the rRNA operons (4.0-4.2 Mb on the *E. coli* reference) — those are multimapper-heavy and both callers will produce noisy calls there; the disagreement is not informative.

</details>

<details>
<summary>Hint 3 — Why is normalization required before isec?</summary>

Consider a 2-bp deletion. bcftools might emit it as `POS=100 REF=AGG ALT=A`. GATK might emit it as `POS=101 REF=GG ALT=` (a different anchor base). Both encodings describe the same biological event, but at the byte-level they have different `POS` and `REF` values, so `bcftools isec` will mark them as different variants.

`bcftools norm -f ref.fa -m -any` does two things: (a) left-aligns indels by shifting them as far left as the reference allows, producing a canonical `POS` for every indel; (b) splits multiallelic records into one record per ALT allele, so a `POS=100 REF=A ALT=G,T` record in one VCF and two `REF=A ALT=G` + `REF=A ALT=T` records in the other become directly comparable.

Without `-f ref.fa`, left-alignment cannot be done; without `-m -any`, the multiallelic mismatch persists. Always run both.

</details>

<details>
<summary>Hint 4 — How do I parse `bcftools isec` output into a counts table?</summary>

```python
import subprocess
from pathlib import Path

compare_dir = Path("compare")
labels = {
    "0000.vcf.gz": "bcftools_only",
    "0001.vcf.gz": "gatk_only",
    "0002.vcf.gz": "both",
}

counts = {}
for fname, label in labels.items():
    path = compare_dir / fname
    # Count SNPs.
    snp = subprocess.check_output(
        ["bcftools", "view", "-H", "-v", "snps", str(path)],
        text=True,
    ).count("\n")
    # Count indels.
    indel = subprocess.check_output(
        ["bcftools", "view", "-H", "-v", "indels", str(path)],
        text=True,
    ).count("\n")
    counts[label] = {"snps": snp, "indels": indel, "total": snp + indel}

print(counts)
```

The intersection-from-bcftools (`0002`) and intersection-from-GATK (`0003`) should have the same *count* but may differ in annotations; use `0002` for the "both" cell of the table.

</details>

<details>
<summary>Hint 5 — What if my caller agreement rate is suspiciously high (or low)?</summary>

If the agreement rate is > 98%, double-check: did you actually run two different callers, or did you run the same caller twice? Did you normalize both VCFs? Did you set `--ploidy 1` for both? A common bug: running `bcftools call` without `--ploidy 1` against *E. coli* produces a different call set than GATK with `--sample-ploidy 1`, because the ploidy mismatch changes the priors.

If the agreement rate is < 80%, double-check: are both VCFs over the same coordinate space? Are both run on the same BAM? Did you accidentally compare a filtered VCF to an unfiltered one (run `bcftools isec` *before* `bcftools filter` so you are comparing the raw calls)?

</details>

## Stretch

If you finish under time and want more:

- Add a third caller: **DeepVariant** (Poplin et al. 2018 — same year as the HaplotypeCaller paper, different team, very different algorithm — deep learning rather than HMM). Compare its calls to bcftools and GATK on the same BAM. DeepVariant is more accurate than both at SNPs in well-trained genomes (human) but less accurate at bacteria.
- Compute the **Ts/Tv ratio** for each caller's PASS variants. Bacterial *E. coli* against `NC_000913.3` should have Ts/Tv ≈ 1.0 (random mutation has no transition/transversion bias). Significant deviation from 1.0 is a caller-bias signal: bcftools historically had slightly elevated Ts/Tv from its base-quality model, while GATK was more neutral.
- Compute a **precision/recall** estimate by treating one caller's PASS variants as the "truth" and computing precision and recall of the other. For a one-sample bacterial dataset there is no real truth set, so this is a *relative* metric (precision-relative-to-X), but it is instructive.
- Read the relevant pages of `mcall.c` (the bcftools multiallelic caller, ~2000 lines of C) and one of the GATK HaplotypeCaller Java source files. The two implementations are an interesting study in how the same statistical problem is approached by two engineering teams.

## What you should be able to do after this

- Run `bcftools call` and GATK `HaplotypeCaller` on the same BAM.
- Normalize VCFs with `bcftools norm -f ref.fa -m -any` and explain why.
- Compute the per-variant intersection and private sets with `bcftools isec`.
- Identify specific positions where the two callers disagree and explain the algorithmic reason in one sentence each.
- Defend a choice of caller in writing, with reference to the dataset's profile (bacterial vs human, single vs cohort, germline vs somatic).

---

*Submit by committing `challenges/compare_callers.py` (or `compare.sh`) and `challenges/notes/caller-comparison.md` to your portfolio repo.*
