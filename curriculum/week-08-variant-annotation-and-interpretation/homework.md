# Week 8 Homework

> **Educational and research use only.** Every problem in this homework produces output that looks like clinical variant-interpretation output. None of it is clinical interpretation. The disclaimer that opens every Week 8 file applies here.

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-08/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Annotate three new sample VCFs with VEP and SnpEff

**Problem statement.** Take the three SRA-derived demo VCFs bundled at `data/p1_sample_{A,B,C}.vcf.gz` (each ~30 variants from the GIAB NA12878 truth set, restricted to chromosome 17). For each sample, run `vep ... --canonical --mane --sift b --polyphen b --af --af_gnomadg --af_gnomade --check_existing` AND `snpEff -v GRCh38.105 ...`. Confirm that VEP and SnpEff agree on the impact tier (`HIGH` / `MODERATE` / `LOW` / `MODIFIER`) for ≥ 95% of variants. For the disagreement subset, write a short note per variant explaining why the two disagree (almost always: missense-vs-inframe-deletion on edge cases, or splice-region vs intron-near-splice).

Answer in `homework/notes/p1-vep-vs-snpeff.md`:

1. Across all three samples, how many variants total? How many had the same impact tier under both VEP and SnpEff?
2. For the disagreement subset (~5% of variants), list each variant and the explanation in one sentence.
3. Which tool's choice would you trust for the disagreements? Why?

**Acceptance criteria.**

- `homework/p1_annotate.sh` runs the six annotation jobs end to end (three samples x two tools).
- `homework/results/p1_sample_{A,B,C}.{vep,snpeff}.vcf` exist.
- `homework/notes/p1-vep-vs-snpeff.md` contains the three answers with specific numbers.
- Commit message like `p1: VEP+SnpEff on 3 samples, 94% impact-tier agreement`.

**Hint.** Use a loop: `for s in A B C; do vep ... -i p1_sample_${s}.vcf.gz -o p1_sample_${s}.vep.vcf; snpEff ... > p1_sample_${s}.snpeff.vcf; done`. The impact-tier comparison can be a one-liner: parse both annotated VCFs with cyvcf2, pull the IMPACT field, join on coordinate, compare.

**Estimated time.** 60 minutes.

---

## Problem 2 — Programmatic gnomAD queries with a cache

**Problem statement.** Extend the Exercise 2 script to handle a 100-variant VCF efficiently. Read the input VCF with cyvcf2, build the gnomAD variant IDs (`chrom-pos-ref-alt`), query gnomAD for each, save to a TSV. Use the shelve-based cache from Exercise 2. Run twice; verify that the second run is at least 100x faster than the first.

Answer in `homework/notes/p2-gnomad-bulk.md`:

1. How many variants are present in gnomAD v4.1? How many absent? (The absent ones are the PM2 evidence subset.)
2. What is the median popmax AF in the "present" subset?
3. How many variants have popmax > 0.05 (BA1 criterion)?
4. First-run time vs second-run time (you should see ~30-60 seconds first, < 1 second second).

**Acceptance criteria.**

- `homework/p2_gnomad_bulk.py` runs end to end on the bundled 100-variant VCF.
- `homework/results/p2_gnomad.tsv` exists with at least the columns: chrom, pos, ref, alt, gnomad_af, gnomad_popmax_af, gnomad_filters, gnomad_present.
- `homework/notes/p2-gnomad-bulk.md` contains four answers with specific numbers.
- Commit message like `p2: gnomAD bulk query, 73/100 present, popmax median 0.002`.

**Hint.** The cache key should include the gnomAD dataset (`gnomad_r4`) so that switching datasets does not silently return stale cached results. The "absent in gnomAD" case is informative for PM2; do not treat it as an error.

**Estimated time.** 60 minutes.

---

## Problem 3 — Build the ACMG classifier on 50 variants

**Problem statement.** Extend Challenge 1's ACMG classifier to read your Problem 1 + Problem 2 outputs and produce a classification for each variant. Implement the eight cleanly mechanically computable criteria (PVS1, PM2, PM4, PP3, BA1, BS2, BP4, BP7). For each variant, report the classification and the evidence set, and explicitly list the skipped criteria.

Answer in `homework/notes/p3-acmg-counts.md`:

1. Across the 100 variants, how many classify as Pathogenic / Likely_pathogenic / VUS / Likely_benign / Benign?
2. How many fire BA1 (stand-alone benign)?
3. How many fire PVS1?
4. How many produce inconsistent evidence (e.g. PVS1 + BA1, which should be impossible if the data is consistent)?

**Acceptance criteria.**

- `homework/p3_acmg.py` runs end to end.
- `homework/results/p3_classified.csv` exists with columns including `acmg_classification`, `acmg_evidence`, `skipped_criteria`, `acmg_warnings`.
- `homework/notes/p3-acmg-counts.md` contains four numbered answers.
- Commit message like `p3: ACMG classifier on 100 variants, 3 LP, 12 LB, 47 VUS, 38 Benign`.

**Hint.** Use the LOF gene set from Challenge 1. The "inconsistent evidence" check is the most informative part — it catches the cases where your input data has a bug (e.g. popmax 0.07 reported for a variant that should be rare).

**Estimated time.** 75 minutes.

---

## Problem 4 — Render the report as HTML

**Problem statement.** Take the Problem 3 classification table and render it as a standalone HTML report. The report must include: the disclaimer at the top, a metadata box (run date, VEP version, ClinVar release, gnomAD version), a per-variant table colour-coded by impact tier, and a separate "skipped criteria" section that lists which ACMG criteria were *not* evaluated.

The HTML must render correctly offline (no external CSS, no external fonts).

**Acceptance criteria.**

- `homework/p4_render.py` runs end to end.
- `homework/results/p4_report.html` exists and renders in a browser.
- The disclaimer is visible at the top of the page.
- The metadata box lists the run date and database versions.
- The per-variant table has at least 12 columns including ACMG classification.
- The "skipped criteria" section lists at least 12 criteria as not-evaluated.
- Commit message like `p4: HTML report with disclaimer, metadata, and skipped-criteria list`.

**Hint.** The Exercise 3 reference template is a starting point; you can extend it. If you prefer Jinja2, that is fine; the standard-library `string.Template` works too.

**Estimated time.** 60 minutes.

---

## Problem 5 — Pharmacogenomics on a small VCF

**Problem statement.** Take the bundled `data/p5_pgx_demo.vcf.gz` (a small VCF of ~12 pharmacogene variants on CYP2D6, CYP2C19, CYP2C9, TPMT, SLCO1B1, VKORC1, DPYD). Query PharmGKB by rsID for each variant. For matched variants, look up the CPIC tier-1 recommendation for the relevant drug(s). Report a per-drug recommendation table.

Answer in `homework/notes/p5-pgx.md`:

1. How many of the 12 variants are PharmGKB-indexed (have a known star-allele mapping)?
2. What is the CPIC recommendation for codeine in a patient heterozygous for CYP2D6*4?
3. What is the recommendation for clopidogrel in a CYP2C19*2/*2 (poor metabolizer) patient?
4. What is the recommendation for simvastatin in a patient heterozygous for SLCO1B1 rs4149056 (T/C)?

**Acceptance criteria.**

- `homework/p5_pgx.py` runs end to end.
- `homework/results/p5_recommendations.csv` exists with at least the columns: gene, rsid, star_allele, genotype, phenotype, drug, recommendation, evidence_level.
- `homework/notes/p5-pgx.md` contains four numbered answers.
- The HTML version of the report prominently includes the "this is not a prescription" disclaimer.
- Commit message like `p5: PharmGKB lookups on 12 PGx variants, 9 indexed, 14 drug recs`.

**Hint.** The simplified phenotype-assignment rule (heterozygous LOF = intermediate, homozygous LOF = poor, heterozygous GOF = ultrarapid, wild-type = normal) is fine for the homework. Real star-allele calling needs Aldy or Stargazer; note this limit in your write-up.

**Estimated time.** 75 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 400-500 word reflection at `homework/notes/week-08-reflection.md` answering:

1. Before Week 8, what did you think "variant interpretation" was? After Week 8, what is it actually? Pick one of the three axes (functional consequence / population frequency / clinical knowledge) and say what surprised you about how the canonical free database for that axis is built.
2. The first time you ran VEP, how did the output compare to what you expected? After Week 8, what does the CSQ field tell you, and what is *not* in it? Give two concrete examples of clinical information that VEP alone cannot provide.
3. The ACMG/AMP 2015 framework has 28 criteria. Your Week 8 implementation covers ~8 of them mechanically. Pick one criterion from the *uncovered* 20 (e.g. PS3, PM3, PP1) and explain in your own words why it cannot be evaluated from a VCF alone. What would you need to evaluate it?
4. The mini-project produces an HTML report. Imagine you hand the report to a non-bioinformatician colleague. What is the *most important* thing the colleague should understand about the report before reading any of the per-variant rows? Why is that the most important thing? (The expected answer involves the disclaimer, but say *why* the disclaimer matters in two sentences.)

**Acceptance criteria.**

- File exists, 400-500 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The boundaries you note here are what will keep you out of trouble in any future job that touches clinical bioinformatics.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 60 min |
| 2 | 60 min |
| 3 | 75 min |
| 4 | 60 min |
| 5 | 75 min |
| 6 | 30 min |
| **Total** | **~6 h** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
