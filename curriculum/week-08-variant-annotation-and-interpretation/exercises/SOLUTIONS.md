# Exercise Solutions — Week 8

> **Educational and research use only.** The reference implementations below produce the same kind of output you would see from a clinical pipeline, but they are not clinical pipelines and their output is not a clinical interpretation. Read these only after attempting your own solution.

Reference implementations and expected numerical outputs for Exercises 1, 2, and 3. The point of an exercise is the friction of doing it the first time; the solution exists to clarify after the fact.

---

## Exercise 1 — VEP annotation of a small VCF

### Reference command

```bash
vep \
    -i demo.vcf \
    -o demo.vep.vcf \
    --vcf \
    --cache \
    --dir_cache ${HOME}/.vep \
    --offline \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache_version 110 \
    --sift b --polyphen b \
    --canonical --mane --symbol --biotype \
    --af --af_gnomadg --af_gnomade \
    --check_existing \
    --fork 4 \
    --force_overwrite \
    --no_stats
```

### Expected output shape

| Field | Value (typical) |
|-------|----------------|
| input variants     | ~50 |
| output variants    | ~50 (no variant should be dropped) |
| Variants with CSQ records | 100% (every variant must overlap >= 1 transcript) |
| Mean transcripts annotated per variant | 5-15 (multi-isoform genes inflate this; intergenic variants give 1 or 2) |
| `IMPACT=HIGH` count | 3-5 |
| `IMPACT=MODERATE` count | 25-30 |
| `IMPACT=LOW` count | 8-12 |
| `IMPACT=MODIFIER` count | 5-12 |

### Reference summary table (paste into `notes/exercise-01-vep-summary.md`)

```markdown
# Exercise 1 — VEP summary

- VEP version: 110.1
- VEP cache version: 110
- Reference assembly: GRCh38
- Date of run: <YYYY-MM-DD>

## Impact-tier counts

| Impact   | Count |
|----------|------:|
| HIGH     | 3     |
| MODERATE | 27    |
| LOW      | 8     |
| MODIFIER | 12    |

## Three hand-curated examples

| chrom | pos       | ref | alt | gene  | consequence            | impact   | rsid       | sift                 | polyphen              | gnomADe_AF |
|-------|----------:|-----|-----|-------|------------------------|----------|------------|----------------------|-----------------------|-----------:|
| chr17 |  43094077 | G   | A   | BRCA1 | missense_variant       | MODERATE | rs28897696 | deleterious(0.02)    | probably_damaging(0.94) | 0.00012  |
| chr5  | 112815473 | A   | T   | APC   | frameshift_variant     | HIGH     | -          | -                    | -                     | -          |
| chr7  | 117559590 | G   | A   | CFTR  | synonymous_variant     | LOW      | rs1042077  | -                    | -                     | 0.07       |
```

### Common mistakes

- **Cache version mismatch.** `--cache_version 110` does not match `${HOME}/.vep/homo_sapiens_vep_111_GRCh38/`. Either reinstall the matching cache or pass the version that is on disk.
- **Forgetting `--mane`.** Without `--mane`, the `MANE_SELECT` field is absent from the CSQ, and Exercise 3's "pick canonical, fall back to MANE" logic falls back to the first record. The annotation still works; it just picks a less-good transcript.
- **Reading the wrong column for the impact tier.** The CSQ field order matters. Field index 2 (0-indexed) is `Consequence`; field index 3 (0-indexed) is `IMPACT`. Confirm with the `##INFO=<ID=CSQ` header line before grepping.

---

## Exercise 2 — Programmatic ClinVar and gnomAD queries

### Reference gnomAD GraphQL query

```graphql
query VariantQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    rsid
    exome {
      ac
      an
      af
      filters
      populations {
        id
        ac
        an
      }
    }
    genome { ... }
  }
}
```

### Reference rate-limit pattern

```python
def fetch_with_retry(method, url, *, json_payload=None, params=None,
                     max_attempts=4, base_sleep=1.0,
                     rate_limit_sleep=0.2, timeout=30.0):
    for attempt in range(max_attempts):
        try:
            if method == "POST":
                r = requests.post(url, json=json_payload, timeout=timeout)
            else:
                r = requests.get(url, params=params, timeout=timeout)
        except requests.RequestException:
            if attempt == max_attempts - 1:
                raise
            time.sleep(base_sleep * (2 ** attempt))
            continue
        if r.status_code in (429, 503):
            time.sleep(base_sleep * (2 ** attempt))
            continue
        r.raise_for_status()
        time.sleep(rate_limit_sleep)
        return r.json()
    raise RuntimeError("all attempts failed")
```

The key elements:

- **Exponential backoff** (`base_sleep * 2 ** attempt`) on 429/503. Doubles per retry to give the server time to recover.
- **Rate-limit sleep** (`time.sleep(0.2)`) after every successful call. Keeps you under 5 req/sec; well below gnomAD's ~10 req/sec cap and NCBI's 3 req/sec cap without an API key.
- **Fail-fast on other errors** (`response.raise_for_status()`). Do not retry 400-level errors that are not 429 — those indicate a malformed request that retrying will not fix.

### Expected output for the 10 demo variants

| chrom | pos       | gene_hint | gnomad_af | gnomad_popmax_af | clinvar_clnsig                          | clinvar_clnrevstat                                |
|-------|----------:|-----------|----------:|-----------------:|------------------------------------------|----------------------------------------------------|
| 17    |  43094077 | BRCA1     | 0.00012   | 0.00018          | Conflicting_interpretations_of_pathogenicity | criteria_provided,_conflicting_classifications  |
| 13    |  32398489 | BRCA2     | 0.005     | 0.012            | Benign/Likely_benign                     | criteria_provided,_multiple_submitters,_no_conflicts |
| 7     | 117559590 | CFTR      | 0.07      | 0.18             | Benign                                   | criteria_provided,_multiple_submitters,_no_conflicts |
| 2     |  47403534 | MSH2      | 0.0008    | 0.002            | Likely_benign                            | criteria_provided,_single_submitter                  |
| 3     |  37001008 | MLH1      | 0.0001    | 0.00015          | Uncertain_significance                   | criteria_provided,_single_submitter                  |
| 5     | 112815473 | APC       | -         | -                | Pathogenic                               | reviewed_by_expert_panel                              |
| 17    |   7674220 | TP53      | -         | -                | Pathogenic                               | reviewed_by_expert_panel                              |
| 11    | 108235829 | ATM       | -         | -                | Likely_pathogenic                        | criteria_provided,_single_submitter                  |
| 15    |  48700503 | FBN1      | -         | -                | Likely_pathogenic                        | criteria_provided,_single_submitter                  |
| 19    |  11116107 | LDLR      | -         | -                | Pathogenic                               | criteria_provided,_multiple_submitters,_no_conflicts |

(A "-" means absent from gnomAD; this is informative for the ACMG PM2 criterion, not an error.)

### Why some variants are absent from gnomAD

The expected pattern is: well-characterized severe-pediatric pathogenic variants (the kind in `APC`, `TP53`, `LDLR`) are usually *absent* or extremely rare in gnomAD because gnomAD's cohort selection excluded severe pediatric disease cases. The absence is the PM2 signal. Common polymorphisms (the kind in `CFTR` at the demo position) are well-represented in gnomAD and have AFs in the 0.01-0.10 range.

### Common mistakes

- **No rate-limit sleep.** Hammers the API at full speed, gets 429, retries with exponential backoff, eventually succeeds, but is rude. Always sleep 0.2-0.4 seconds between requests as a polite-citizen baseline.
- **Caching the empty response on a transient error.** If a request times out and you cache the empty result, the next run will see "absent" forever even after the API recovers. Cache only successful responses; let errors fall through to a retry.
- **Forgetting the variant-ID format.** gnomAD's variant ID is `chrom-pos-ref-alt` with no "chr" prefix. ClinVar's E-utilities accepts the same coordinates but in different query syntax. Read the API docs once, write a constructor function, and reuse it.

---

## Exercise 3 — Build the per-variant report

### Reference CSQ-parsing pattern

```python
def parse_csq_field(csq_field, csq_header):
    if not csq_field:
        return []
    records = []
    for raw in csq_field.split(","):
        fields = raw.split("|")
        if len(fields) != len(csq_header):
            continue
        records.append(dict(zip(csq_header, fields)))
    return records


def pick_canonical_record(records):
    canonical = next((r for r in records if r.get("CANONICAL") == "YES"), None)
    if canonical is not None:
        return canonical
    mane = next((r for r in records if r.get("MANE_SELECT")), None)
    if mane is not None:
        return mane
    return records[0] if records else {}
```

### Reference SIFT/PolyPhen score parser

```python
import re

_SCORE_PATTERN = re.compile(r"([a-zA-Z_]+)\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)")


def parse_prediction_and_score(field: str) -> tuple[str, float | None]:
    if not field:
        return "", None
    m = _SCORE_PATTERN.match(field)
    if m is None:
        return field, None
    return m.group(1), float(m.group(2))
```

### Expected CSV columns (in order)

```
chrom, pos, ref, alt, gene, consequence, impact, hgvsc, hgvsp,
sift_pred, sift_score, polyphen_pred, polyphen_score,
rsid, gnomad_af, gnomad_popmax_af,
clinvar_clnsig, clinvar_clnrevstat, clinvar_clndn
```

### Expected HTML report structure

- An H1 title: "Variant Annotation Report"
- A yellow-background disclaimer box: "Educational and research use only..." 
- A grey metadata box: run date, input VCF, VEP cache version, ClinVar release, gnomAD version
- A coloured per-variant table (HIGH red-tinted, MODERATE yellow-tinted, LOW green-tinted, MODIFIER grey)
- A footer: "Generated by exercise-03-build-variant-report.py..."

The HTML must include the disclaimer, and it must include the database version metadata. Without these, the report is not reproducible.

### Reference run-info JSON

```json
{
  "run_date": "2024-10-15",
  "input_vcf": "demo.vep.vcf",
  "vep_cache_version": "110",
  "clinvar_release": "2024-09-01",
  "gnomad_version": "v4.1.0",
  "n_variants": 10,
  "disclaimer": "EDUCATIONAL AND RESEARCH USE ONLY. ..."
}
```

### Why "pick canonical, fall back to MANE, fall back to first"

A single variant can overlap 5-20 transcripts in the multi-isoform genes. VEP emits one CSQ record per transcript, with no guarantee on the order. The choice of which record to report as "the" annotation matters for clinical reporting.

The hierarchy:

1. **Ensembl canonical** (the longest-protein-coding-transcript-by-default in each gene). VEP's `--canonical` flag adds a `CANONICAL` field that is "YES" or empty. Pick the YES record.
2. **MANE Select** (the joint Ensembl-NCBI canonical, available since 2022 for ~95% of human protein-coding genes). VEP's `--mane` flag adds the `MANE_SELECT` field. Use this as the fallback when CANONICAL is not defined.
3. **First record** (whatever VEP put first). The last-resort fallback. For intergenic and intron-only variants, this is often the only record.

The Exercise 3 reference implementation chains these three with `next((... if ...), None) or next((... if ...), None) or records[0]`. The chain is short, explicit, and easy to debug.

### Common mistakes

- **Picking the wrong transcript.** If you skip the `pick_canonical_record` step and just take the first CSQ record, you may report a non-canonical transcript's annotation. The gene and consequence will be correct (both VEP and the underlying gene model agree on the gene), but the HGVS coordinates will refer to a non-canonical transcript that the clinical literature does not use.
- **Dropping the disclaimer from the HTML output.** The disclaimer must be present on every report. Removing it because "the pipeline is good" is the moment to put it back in.
- **Forgetting the run-info JSON.** A report without database versions is not reproducible. Re-running the pipeline a year later with newer databases will produce different annotations, and you will have no record of what the old run actually queried.
- **HTML-escaping issues.** If you serialize CSQ values with embedded `<` or `&` characters (rare but possible), the HTML output breaks. Use pandas' `to_html(escape=True)` to handle this safely.
- **Joining on the wrong key.** The Exercise 2 TSV has chrom without "chr" prefix; the VEP VCF has "chr". Normalize before joining or you will get an empty merge.

---

## Putting it together

Exercises 1, 2, and 3 form a complete annotation pipeline:

```
demo.vcf
    |
    | (Exercise 1: vep ... --canonical --mane --sift b --polyphen b --check_existing)
    v
demo.vep.vcf
    |
    | (Exercise 2: python exercise-02-query-clinvar-gnomad.py)
    v
exercise-02-clinvar-gnomad.tsv
    |
    | (Exercise 3: python exercise-03-build-variant-report.py)
    v
exercise-03-variants.csv
exercise-03-variants.html
exercise-03-run-info.json
```

The mini-project (`mini-project/README.md`) wraps these three steps into a single `run.sh` and adds an ACMG classifier on top.

---

## Where to go next

- The mini-project takes a fresh VCF (not the demo VCF) and runs the same pipeline end to end, producing the report that a non-bioinformatician colleague can read.
- Challenge 1 implements the mechanically computable subset of the ACMG criteria from Lecture 3 and emits the five-tier classification per variant.
- Challenge 2 adds a pharmacogenomics axis: take a VCF, query PharmGKB, produce a per-drug recommendation report.
