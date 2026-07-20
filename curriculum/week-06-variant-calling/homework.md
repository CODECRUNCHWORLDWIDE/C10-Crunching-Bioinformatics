# Week 6 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-06/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Call variants on the Week 5 mini-project BAM

**Problem statement.** Take the sorted, indexed, duplicate-marked BAM you produced in Week 5's mini-project (`SRR1770413.markdup.bam` aligned against `NC_000913.3`). Run the canonical `bcftools mpileup | bcftools call -m --ploidy 1` pipeline against it. Save the raw VCF, hard-filter it with `bcftools filter` using the bacterial-friendly recipe (`QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60`), normalize with `bcftools norm`, and produce a `bcftools stats` summary. Time the call with `time` and record the wall-clock duration in `homework/notes/p1-call-variants.md`.

Answer:

1. How long did `bcftools mpileup | bcftools call` take?
2. How many variants did the raw VCF contain?
3. How many variants PASS the hard filters?
4. What is the Ts/Tv ratio reported by `bcftools stats`?

**Acceptance criteria.**

- `homework/p1_call.sh` runs end to end on the markdup BAM.
- VCFs at `homework/calls/p1.raw.vcf.gz`, `homework/calls/p1.filtered.vcf.gz`, `homework/calls/p1.norm.vcf.gz` are committable.
- `notes/p1-call-variants.md` contains four numbered, numeric answers.
- Commit message like `p1: SRR1770413 variants against NC_000913.3, NNN raw, NNN PASS`.

**Hint.** The full pipeline runs in ~1 minute on a laptop for the *E. coli* dataset. If yours takes much longer, you may be missing the `.fai` index on the reference. Run `samtools faidx ref/ecoli.fa` first.

**Estimated time.** 45 minutes.

---

## Problem 2 — GATK HaplotypeCaller on the same BAM

**Problem statement.** Repeat Problem 1's variant call with GATK `HaplotypeCaller -R ref/ecoli.fa -I aln/SRR1770413.markdup.bam -O calls/p2.gatk.vcf.gz --sample-ploidy 1` instead of `bcftools`. Normalize with `bcftools norm`. In `homework/notes/p2-gatk-vs-bcftools.md`, compare:

1. Wall-clock time for each caller.
2. Total variants called by each (raw, before any filtering).
3. Number of SNPs and indels for each.
4. Run `bcftools isec -p compare/` between the two normalized VCFs. Report the counts in the four output files (private to bcftools, private to GATK, intersection-from-bcftools, intersection-from-GATK).
5. Pick one variant where bcftools called something but GATK did not, and one where GATK called something but bcftools did not. Paste both records side by side and explain the disagreement.

**Acceptance criteria.**

- `homework/p2_gatk.sh` runs.
- VCF at `homework/calls/p2.gatk.norm.vcf.gz` is committable.
- `notes/p2-gatk-vs-bcftools.md` contains five numbered answers with specific numbers and two named-position case studies.
- Commit message like `p2: bcftools vs GATK on SRR1770413`.

**Hint.** If GATK is unavailable, substitute a second `bcftools` run with stricter mpileup parameters (`-q 20 -Q 20`). The comparison machinery is identical; you just lose the most interesting (algorithmic) source of disagreement.

**Estimated time.** 75 minutes.

---

## Problem 3 — Decode INFO and FORMAT for 20 variants by hand

**Problem statement.** Take the first 20 PASS variants from your Problem 1 normalized VCF:

```bash
bcftools view -f PASS homework/calls/p1.norm.vcf.gz \
| grep -v '^#' | head -20 > homework/notes/p3-variants.txt
```

For each record, write a Markdown table at `homework/notes/p3-vcf-decoding.md` with columns:

| CHROM | POS | REF | ALT | QUAL | DP | MQ | SP | GT | AD | One-sentence interpretation |
|-------|----:|-----|-----|-----:|---:|---:|---:|-----|-----|-----------------------------|

E.g.:

| NC_000913.3 | 150123 | A | G | 227 | 42 | 60 | 2 | 1 | 0,42 | Homozygous-alt SNP, depth 42, all alt reads, mapping quality 60 — high-confidence call. |

Do not use any decoding library — work from the text VCF and Lecture 1 §5. Verify your work by spot-checking with `bcftools view -v snps -f PASS` and comparing field-by-field.

**Acceptance criteria.**

- `notes/p3-vcf-decoding.md` has 20 rows, each correctly decoded.
- Each row's interpretation is one English sentence.
- At least one row is an indel (REF or ALT longer than 1 base).

**Hint.** The INFO column is semicolon-separated; the FORMAT and per-sample columns are colon-separated and parallel. Walk through one record with Lecture 1 §5.2 open, then the next 19 are mechanical.

**Estimated time.** 45 minutes.

---

## Problem 4 — VEP annotation on the filtered VCF

**Problem statement.** Install a VEP cache for *E. coli* (or use the REST API if installation is impractical):

```bash
vep_install --AUTO c --SPECIES escherichia_coli_str_k_12_substr_mg1655 --CACHEDIR vep_cache/
```

Annotate your Problem 1 normalized VCF:

```bash
vep --input_file homework/calls/p1.norm.vcf.gz \
    --output_file homework/calls/p4.vep.vcf \
    --species escherichia_coli_str_k_12_substr_mg1655 \
    --cache --dir_cache vep_cache/ \
    --vcf --symbol --canonical --force_overwrite
```

In `homework/p4_vep_summary.py`:

1. Read `p4.vep.vcf` with `pysam.VariantFile`.
2. For each PASS variant, parse the `CSQ` field and extract the canonical-transcript consequence.
3. Build a counts table by consequence term (`missense_variant`, `synonymous_variant`, `intergenic_variant`, etc.).
4. Build a counts table by IMPACT category (`HIGH`, `MODERATE`, `LOW`, `MODIFIER`).
5. Save both tables to `notes/p4-vep-summary.md` and write 50-100 words interpreting them.

Answer in the writeup: are there any HIGH-impact variants? Which gene(s)? Does the consequence-distribution look biologically reasonable for an *E. coli* sample-vs-reference comparison?

**Acceptance criteria.**

- `homework/p4_vep_summary.py` runs and produces the two counts tables.
- `notes/p4-vep-summary.md` contains both tables and the 50-100 word interpretation.
- `p4.vep.vcf` is committable.

**Hint.** If `vep_install` is too painful, use the REST API approach from Exercise 3 for a small number of variants. The CSQ field has the same format whichever route you take.

**Estimated time.** 60 minutes.

---

## Problem 5 — Compare hard-filter recipes

**Problem statement.** Take your Problem 1 raw VCF (before filtering) and apply three different hard-filter recipes:

1. **Strict**: `QUAL<50 || INFO/DP<15 || INFO/MQ<50 || INFO/SP>40`.
2. **Standard** (the one from Lecture 2 §4.3): `QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60`.
3. **Loose**: `QUAL<10 || INFO/DP<5 || INFO/MQ<20 || INFO/SP>100`.

For each recipe, count: total variants, PASS, LowQual, and (using `bcftools view -v snps/indels`) PASS-SNPs and PASS-indels.

In `homework/notes/p5-filter-comparison.md`:

1. Tabulate the counts (3 rows for the 3 recipes, columns for PASS, LowQual, PASS-SNPs, PASS-indels).
2. Interpret: how many variants does each recipe accept? Which recipe is most conservative? Which is most liberal?
3. Speculate: if your downstream use case were *antibiotic-resistance screening* (where false negatives are dangerous), which recipe would you pick? If it were *population-level variant cataloguing* (where false positives are dangerous), which would you pick?

**Acceptance criteria.**

- `homework/p5_filter_comparison.py` runs all three filters and prints the table.
- `notes/p5-filter-comparison.md` contains the table and the two interpretation paragraphs.

**Hint.** `bcftools filter -e <expression>` does not modify the variant set; it adds a `FILTER` tag. To get the actual PASS set, follow up with `bcftools view -f PASS`. Or use `bcftools filter -i <expression>` (include instead of exclude) for a different idiom that produces a PASS-only set directly.

**Estimated time.** 60 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300-400 word reflection at `homework/notes/week-06-reflection.md` answering:

1. Before Week 6, what did you think a VCF was? What is it actually? Pick one column of the eight you found most surprising and say why.
2. The first time you saw `bcftools call` emit a variant in a real BAM, what did you assume the QUAL value meant? After Week 6, what does it actually mean? In what way does the distinction matter for downstream filtering?
3. The GATK Best Practices hard filters look like an arbitrary set of magic numbers. After Week 6, what is each one actually filtering for? Pick the one you find most counterintuitive (probably `MQRankSum` or `ReadPosRankSum`) and explain in your own words what it measures.
4. The mini-project asks you to produce a VEP-annotated VCF from a real SRA dataset. What is the difference between a "called variant" and an "annotated variant"? Why does the second exist? What kind of question can you answer with the annotated VCF that you cannot answer with the called VCF?

**Acceptance criteria.**

- File exists, 300-400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The mistakes you note here are what you will re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 45 min |
| 2 | 1 h 15 min |
| 3 | 45 min |
| 4 | 1 h 0 min |
| 5 | 1 h 0 min |
| 6 | 30 min |
| **Total** | **~5 h 15 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
