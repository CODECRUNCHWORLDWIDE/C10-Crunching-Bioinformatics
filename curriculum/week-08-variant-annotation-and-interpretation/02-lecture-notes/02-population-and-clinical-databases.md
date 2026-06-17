# Lecture 2 — Population and Clinical Databases

> **Educational and research use only.** This lecture explains how to query the public population (gnomAD) and clinical (ClinVar, dbSNP) databases. The queries are mechanical lookups against public APIs and downloads. They produce evidence rows, not clinical interpretations. The same disclaimer that opens every Week 8 file applies here.

> **Duration:** ~2.5 hours of reading + a small Python query session.
> **Outcome:** You can query the gnomAD GraphQL API and the NCBI E-utilities programmatically, parse the JSON responses into structured per-variant rows, handle the rate-limits and the retry logic correctly, and join the database evidence with the VEP annotation from Lecture 1.

If you only remember one thing from this lecture, remember this:

> **gnomAD answers "how OFTEN" in a sub-population-stratified way. ClinVar answers "what is KNOWN" with a review-status star rating. Both are queried programmatically over free public APIs with rate limits in the 3-10 requests/second range, so for any non-trivial VCF the right pattern is "batch the variants, sleep between requests, retry on 429/503, and cache the responses locally." Without caching, your annotation re-fetches the same data on every run and your pipeline takes ten times as long.**

Lecture 1 closed at a 13-column per-variant table built from VEP's CSQ field. Lecture 2 enriches each row with the population frequency and clinical-knowledge axes.

---

## 1. gnomAD — the population frequency database

gnomAD (the Genome Aggregation Database; Karczewski et al. 2020, *Nature* 581:434) is the canonical answer to "how common is this variant in healthy populations." The v4.1 release (April 2024) aggregates **807,162 samples** from biobanks, case studies, and consented research cohorts. The cohort intentionally *excludes* severe pediatric disease cases and the immediate relatives of any case in a database used for case-control studies, so the aggregate is closer to "the apparently healthy population" than to "the general population."

The v4.1 sample composition:

| Subpopulation              | Code | Exome samples | Genome samples |
|----------------------------|------|---------------|----------------|
| African / African-American | afr  | 21,141        | 22,484         |
| Admixed American / Latino  | amr  | 30,328        | 8,254          |
| Ashkenazi Jewish           | asj  | 5,040         | 1,977          |
| East Asian                 | eas  | 19,872        | 2,612          |
| Finnish                    | fin  | 25,158        | 5,316          |
| Middle Eastern             | mid  | 3,031         | 158            |
| Non-Finnish European       | nfe  | 414,754       | 34,029         |
| South Asian                | sas  | 102,420       | 2,419          |
| Other                      | oth  | 8,066         | 1,047          |
| **Total**                  |      | **629,810**   | **78,296**     |

Two things to notice. First, the **non-Finnish European** bucket is by far the largest, which means the per-variant allele count for any rare variant is most precise in that sub-population. The Middle Eastern and Other buckets are tiny, so per-variant frequencies in those buckets have large confidence intervals (a variant seen once in `mid` has AF 1/3031 ≈ 0.0003, but the 95% CI is 0.00002 to 0.002). Second, the **popmax** (popmax allele frequency) is the maximum AF across non-bottlenecked sub-populations (Finnish is excluded from popmax because it is a known bottleneck, so a single carrier in a Finnish founder family does not push popmax). Most ACMG-flavored filters pivot on popmax, not on global AF.

### Querying gnomAD via the GraphQL API

gnomAD exposes a GraphQL endpoint at `https://gnomad.broadinstitute.org/api`. The query format is GraphQL — you declare the fields you want, the API returns exactly those fields, and rate limiting is per-IP at ~10 requests/second. A minimal query for one variant:

```python
from __future__ import annotations
import time
from typing import Any

import requests


GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_RATE_LIMIT_SLEEP_SEC = 0.2  # 5 req/sec; well below the ~10 req/sec cap.

GNOMAD_QUERY = """
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
    genome {
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
  }
}
"""


def query_gnomad(variant_id: str, dataset: str = "gnomad_r4", timeout: float = 30.0) -> dict[str, Any]:
    """Look up a single variant in gnomAD.

    Args:
        variant_id: gnomAD variant identifier in the form
                    "<chrom>-<pos>-<ref>-<alt>" with no "chr" prefix,
                    e.g. "7-117559590-G-A".
        dataset:    gnomAD dataset ID; one of:
                      gnomad_r2_1  (v2 exomes, GRCh37)
                      gnomad_r3    (v3 genomes, GRCh38)
                      gnomad_r4    (v4 combined, GRCh38; the default)
        timeout:    request timeout in seconds.

    Returns:
        the "variant" subtree of the JSON response, or {} if the
        variant is not in gnomAD (which is informative — absence
        from gnomAD is the PM2 criterion of the ACMG framework).
    """
    payload = {
        "query": GNOMAD_QUERY,
        "variables": {"variantId": variant_id, "dataset": dataset},
    }
    response = requests.post(GNOMAD_API, json=payload, timeout=timeout)
    response.raise_for_status()
    data = response.json()
    variant: dict[str, Any] | None = data.get("data", {}).get("variant")
    return variant or {}
```

The response, for a variant present in gnomAD:

```json
{
  "variant_id": "7-117559590-G-A",
  "rsid": "rs113993960",
  "exome": {
    "ac": 14,
    "an": 19532,
    "af": 0.000717,
    "filters": [],
    "populations": [
      {"id": "afr", "ac": 0, "an": 1024},
      {"id": "amr", "ac": 1, "an": 2148},
      {"id": "asj", "ac": 0, "an": 312},
      {"id": "eas", "ac": 0, "an": 1572},
      {"id": "fin", "ac": 2, "an": 1280},
      {"id": "nfe", "ac": 11, "an": 11512},
      {"id": "sas", "ac": 0, "an": 1684}
    ]
  },
  "genome": { ... }
}
```

The `af` field is the global allele frequency (allele count / allele number). The `populations` array gives per-subpopulation counts. The `filters` array is empty for PASS variants and non-empty (e.g. `["RF"]`) for variants flagged as low quality.

To compute popmax from the response:

```python
def compute_popmax_af(populations: list[dict[str, Any]],
                     exclude: tuple[str, ...] = ("fin", "asj", "oth", "mid")) -> float:
    """Return the maximum AF across non-bottlenecked sub-populations.

    Following gnomAD's convention, Finnish, Ashkenazi Jewish, Other,
    and Middle Eastern are excluded from popmax because they are
    bottlenecked or small.

    Args:
        populations: gnomAD populations array (one dict per subpop).
        exclude:     subpop IDs to exclude from popmax.

    Returns:
        max AF over the included subpops, 0.0 if no subpop has data.
    """
    best: float = 0.0
    for entry in populations:
        if entry["id"] in exclude:
            continue
        ac: int = int(entry["ac"])
        an: int = int(entry["an"])
        if an == 0:
            continue
        af = ac / an
        if af > best:
            best = af
    return best
```

### Rate limits, retries, and caching

The gnomAD API rate limit is ~10 requests/second per IP, with bursts allowed. Hit the limit and you get HTTP 429. Sustained high load returns 503. Production code must handle both. The polite pattern, in 30 lines:

```python
import time
from typing import Any
import requests


def fetch_with_retry(url: str, payload: dict[str, Any],
                     max_attempts: int = 4, base_sleep: float = 1.0,
                     rate_limit_sleep: float = 0.2,
                     timeout: float = 30.0) -> dict[str, Any]:
    """Fetch with exponential backoff on 429 / 503.

    Sleeps `rate_limit_sleep` after every successful call to space out
    requests. Retries up to `max_attempts` times on 429 or 503 with
    exponential backoff (base_sleep * 2 ** attempt).

    Raises requests.HTTPError on any other non-2xx response.
    """
    for attempt in range(max_attempts):
        try:
            response = requests.post(url, json=payload, timeout=timeout)
        except requests.RequestException:
            if attempt == max_attempts - 1:
                raise
            time.sleep(base_sleep * (2 ** attempt))
            continue
        if response.status_code in (429, 503):
            time.sleep(base_sleep * (2 ** attempt))
            continue
        response.raise_for_status()
        time.sleep(rate_limit_sleep)
        return response.json()
    raise RuntimeError(f"All {max_attempts} attempts failed for {url}")
```

Always **cache** the responses locally. The simplest cache is a single SQLite table keyed by `(variant_id, dataset, gnomad_version)`. The Week 8 exercises use a small `shelve`-based cache because it is in the standard library; for production you would use SQLite or a key-value store.

### When the API is the wrong choice

For a VCF of 50-500 variants, the API is fine. For a VCF of 50,000 variants (typical whole exome), per-variant API calls take 50,000 * 0.2 sec = ~3 hours. Download the bundled gnomAD VCF instead. For chr17:

```bash
curl -sLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz
curl -sLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz.tbi
```

Each chromosome is 1-5 GB. Tabix-extract the per-variant records you need:

```bash
tabix gnomad.exomes.v4.1.sites.chr17.vcf.bgz chr17:43094077-43094077
```

Then parse with cyvcf2 the same way you parsed the ClinVar VCF.

Use the API when the input VCF is small (< 1,000 variants) and the latency of a download is not justified. Use the bundled VCF when the input VCF is large or you need to run repeatedly.

---

## 2. ClinVar — the clinical-knowledge database

ClinVar (Landrum et al. 2018, *Nucleic Acids Research* 46:D1062) is the NCBI clinical-variant database. As of 2024 it contains ~3 million records contributed by clinical labs, expert panels, and curated literature reviews. Each variant can have multiple submissions; ClinVar resolves the per-variant *aggregate* classification by counting submitters and review status.

### The CLNSIG field

The per-variant clinical significance is one of:

- **Pathogenic** — the variant causes the disease.
- **Likely_pathogenic** — strong evidence, not quite over the threshold.
- **Uncertain_significance** (VUS) — insufficient evidence.
- **Likely_benign** — strong evidence against pathogenicity.
- **Benign** — no role in the disease.
- **Conflicting_interpretations_of_pathogenicity** — submitters disagree, ~3% of records.
- **drug_response** — variant affects a drug response (overlap with pharmacogenomics).
- **risk_factor** — increases disease risk but is not causal.
- **association** — statistically associated, mechanism unclear.
- **protective** — decreases risk.
- **other** / **not_provided**.

The clean five-tier ACMG-style classification is the first five above; the others are "outside the ACMG framework" or "no consensus."

### The CLNREVSTAT field (review status)

ClinVar grades each record on a 4-star scale:

| Stars | CLNREVSTAT string                                                  | Meaning |
|-------|--------------------------------------------------------------------|---------|
| 4     | `practice_guideline`                                               | In a published clinical guideline (e.g. AAP, ACMG-recommended screening) |
| 3     | `reviewed_by_expert_panel`                                         | ClinGen, ENIGMA, InSiGHT, etc. expert curation |
| 2     | `criteria_provided,_multiple_submitters,_no_conflicts`             | Multiple labs assert the same classification with criteria |
| 2     | `criteria_provided,_conflicting_classifications`                   | Multiple labs, conflicting; the aggregate classification is "Conflicting" |
| 1     | `criteria_provided,_single_submitter`                              | One lab with documented criteria |
| 0     | `no_assertion_criteria_provided`                                   | A submission without criteria |
| 0     | `no_assertion_provided`                                            | A submission without a classification |

For pipeline filtering, 3-star and 4-star records are the most authoritative. 2-star is good. 1-star is informative but should be checked. 0-star records exist for historical reasons and should be down-weighted.

### Querying ClinVar

There are three ways to query ClinVar programmatically, each with its own use case:

#### a. Download the ClinVar VCF release (the production path)

```bash
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

The VCF is ~200 MB compressed, indexed, updated every two weeks. Each variant record has the INFO fields documented in the header: `CLNSIG`, `CLNREVSTAT`, `CLNDN`, `CLNDISDB`, `CLNHGVS`, `MC` (molecular consequence), `RS` (dbSNP rsID), and `GENEINFO`. To look up a variant by coordinate:

```bash
tabix clinvar.vcf.gz chr17:43094077-43094077
```

Or with cyvcf2:

```python
from cyvcf2 import VCF


def query_clinvar_by_position(clinvar_vcf: VCF, chrom: str, pos: int,
                              ref: str, alt: str) -> dict[str, str]:
    """Look up the per-variant ClinVar record by exact-match coordinate + alleles."""
    region: str = f"{chrom}:{pos}-{pos}"
    for record in clinvar_vcf(region):
        if record.REF != ref:
            continue
        if alt not in record.ALT:
            continue
        return {
            "clnsig": record.INFO.get("CLNSIG", ""),
            "clnrevstat": record.INFO.get("CLNREVSTAT", ""),
            "clndn": record.INFO.get("CLNDN", ""),
            "clndisdb": record.INFO.get("CLNDISDB", ""),
            "rs": record.INFO.get("RS", ""),
            "geneinfo": record.INFO.get("GENEINFO", ""),
            "mc": record.INFO.get("MC", ""),
        }
    return {}
```

This is the path the mini-project uses. Tabix is fast (~50 ms per lookup), the data is local, and the database version is pinned to the file.

#### b. Query the NCBI E-utilities (the small-VCF path)

```python
import requests


def query_clinvar_eutils(rsid: str, timeout: float = 30.0) -> dict[str, Any]:
    """Look up a ClinVar record by dbSNP rsID via E-utilities.

    Slower than tabix (network round-trip), but does not require the
    local download.
    """
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "clinvar", "term": f"{rsid}[rs]", "retmode": "json"},
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()
```

E-utilities are rate-limited to 3 req/sec without an API key (10 req/sec with). For VCFs under ~500 variants this works; above that, download the VCF.

#### c. Use the VEP plugin (the simplest path)

VEP can include ClinVar fields via the `--custom` flag pointing at the ClinVar VCF:

```bash
vep -i input.vcf -o output.vep.vcf \
    --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --cache --offline --species homo_sapiens --assembly GRCh38 ...
```

This embeds the ClinVar fields directly in the CSQ output and is the most ergonomic for the per-variant report. The mini-project uses this.

---

## 3. dbSNP — the variant naming authority

dbSNP (Sherry et al. 2001) is the NCBI repository of named variants. The rsID is the canonical handle for a variant in the literature: when a paper says "carriers of rs121908755 had increased risk of X," dbSNP is where you confirm the variant coordinates, ALT allele, and surrounding context.

For Week 8 we use dbSNP for two things:

- **Confirming the rsID** of a variant before reporting it.
- **Cross-linking** an rsID to ClinVar (most ClinVar records carry the RS field).

Queries:

```python
import requests


def query_dbsnp(rsid_int: int, timeout: float = 30.0) -> dict[str, Any]:
    """Look up a dbSNP record by integer rsID (no rs prefix).

    Example: rsid_int=121908755 for rs121908755.
    """
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "snp", "id": str(rsid_int), "retmode": "json"},
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()
```

Most of what you would want from dbSNP (the canonical rsID, the variant class, the allele frequencies in the global cohort) is already in the VEP CSQ output via the `Existing_variation` field. dbSNP-as-a-separate-query is mostly useful when you have an rsID from the literature and want to translate it to a coordinate, which is the inverse of the VCF-annotation flow.

---

## 4. Joining the three axes in Python

The Lecture 1 code emits one row per variant from the VEP CSQ field. The Lecture 2 code adds per-variant gnomAD and ClinVar columns. Putting them together:

```python
from __future__ import annotations
from pathlib import Path
from typing import Any

import pandas as pd
from cyvcf2 import VCF


def annotate_variants(vep_vcf: Path,
                      clinvar_vcf: Path,
                      gnomad_dataset: str = "gnomad_r4") -> pd.DataFrame:
    """Build the per-variant interpretation frame.

    For each variant in the VEP-annotated VCF, emit a row with:
      chrom, pos, ref, alt, rsid, gene, consequence, impact,
      hgvsc, hgvsp, sift, polyphen,
      clinvar_clnsig, clinvar_clnrevstat, clinvar_clndn,
      gnomad_af, gnomad_popmax_af, gnomad_filters

    Args:
        vep_vcf:        path to a VEP-annotated VCF (Lecture 1 output).
        clinvar_vcf:    path to the local ClinVar VCF release (.vcf.gz with .tbi).
        gnomad_dataset: gnomAD dataset ID (default gnomad_r4).

    Returns:
        pandas DataFrame, one row per variant.
    """
    from cyvcf2 import VCF

    csq_header = extract_csq_header(vep_vcf)
    clinvar = VCF(str(clinvar_vcf))
    rows: list[dict[str, Any]] = []

    for variant in VCF(str(vep_vcf)):
        base = summarize_variant(variant, csq_header)
        # ClinVar local-VCF join.
        cv = query_clinvar_by_position(clinvar,
                                        variant.CHROM, int(variant.POS),
                                        variant.REF, variant.ALT[0])
        base["clinvar_clnsig"] = cv.get("clnsig", "")
        base["clinvar_clnrevstat"] = cv.get("clnrevstat", "")
        base["clinvar_clndn"] = cv.get("clndn", "")
        # gnomAD API call.
        variant_id = f"{variant.CHROM.lstrip('chr')}-{variant.POS}-{variant.REF}-{variant.ALT[0]}"
        try:
            g = query_gnomad(variant_id, dataset=gnomad_dataset)
        except Exception:  # noqa: BLE001
            g = {}
        exome = (g or {}).get("exome") or {}
        populations = exome.get("populations", [])
        base["gnomad_af"] = exome.get("af", "")
        base["gnomad_popmax_af"] = compute_popmax_af(populations) if populations else ""
        base["gnomad_filters"] = ";".join(exome.get("filters", []) or [])
        rows.append(base)

    return pd.DataFrame(rows)
```

The output is a single pandas frame that can be saved as CSV, rendered as HTML, or passed into the ACMG classifier from Lecture 3. The frame is the input to the mini-project.

---

## 5. Caveats, sub-population stratification, and the small-cohort problem

A few things to be careful about when reading the gnomAD-and-ClinVar columns:

### Sub-population stratification

The popmax is the maximum AF across non-bottlenecked sub-populations. A variant with global AF 0.0001 but Finnish AF 0.07 has popmax 0.0001 (Finnish is excluded from popmax by gnomAD convention) — but for a Finnish patient that variant is *common*. The ACMG BA1 criterion uses popmax, but for an individual interpretation you should also check the patient's ancestry and the per-sub-population AF.

For the per-variant report, expose both the popmax and the matched-ancestry AF where the patient's ancestry is known. The Week 8 demo VCF does not have ancestry metadata, so we report popmax and the full sub-population breakdown.

### The "absent in gnomAD" case

A variant absent from gnomAD (~30-50% of the variants in a typical exome) is a *strong* PM2 ACMG criterion, *if* the gene was well-covered in the gnomAD samples. Some genes (e.g. the highly homologous gene families, the GC-rich genes, the centromere-proximal genes) have poor gnomAD coverage and "absent" is uninformative. gnomAD reports per-gene coverage tracks; the `coverage` field on a variant indicates whether the position was well-covered in the gnomAD callset.

### Conflicting ClinVar interpretations

When a variant has multiple ClinVar submissions with conflicting classifications, ClinVar resolves the per-variant `CLNSIG` to `Conflicting_interpretations_of_pathogenicity` and the `CLNREVSTAT` to `criteria_provided,_conflicting_classifications`. The right action is to flag the variant for manual review and inspect the individual submissions on the ClinVar web UI. Do not paper over the conflict by picking the most recent or the most-stars submission.

### Cross-version drift

gnomAD v2 (~125K samples, GRCh37) and gnomAD v4 (~807K samples, GRCh38) report different AFs for the same variant — sometimes by a factor of 2 or more — because the sample composition changed substantially between releases. Always pin the gnomAD version in the run-info JSON. Same for ClinVar: the same variant's CLNSIG can change between ClinVar releases as new submissions arrive and old ones are reclassified.

The run-info JSON should record, at minimum:

```json
{
  "run_date": "2024-10-15",
  "vep_version": "110.1",
  "vep_cache_version": "110",
  "snpeff_version": "5.2",
  "snpeff_db": "GRCh38.105",
  "clinvar_release": "2024-09-01",
  "gnomad_version": "v4.1.0",
  "assembly": "GRCh38"
}
```

Without these, a repeat run a year later will produce different answers and you will be unable to debug the difference.

---

## 6. The `IMPACT=HIGH but ClinVar=Benign` case

A variant can have all the signs of pathogenicity (frameshift, in a disease gene, predicted HIGH impact) and still be ClinVar Benign. This happens for several reasons:

- The variant is in a gene where loss-of-function is *not* a disease mechanism. Many genes are haplosufficient — losing one copy is well-tolerated. The ACMG PVS1 criterion explicitly conditions on "LOF is a known disease mechanism."
- The frameshift is in the last exon, past the last NMD-rule boundary. Frameshifts in the last ~50 bp of the penultimate exon or anywhere in the last exon escape nonsense-mediated decay, so the truncated transcript is translated and may be functional.
- The variant is in a non-canonical transcript whose loss is tolerated; the canonical transcript is unaffected.
- The variant has been *reclassified* over time. ClinVar records reflect the current consensus, and consensus changes.

The right reading of "IMPACT=HIGH and ClinVar=Benign with 3+ stars" is: trust ClinVar. The clinical evidence has overruled the in-silico prediction. The in-silico prediction is a hypothesis; the clinical evidence is the test.

The opposite case — "IMPACT=LOW and ClinVar=Pathogenic" — is rarer but exists, usually for splice-region variants that VEP buckets as LOW but that disrupt splicing in vivo. SpliceAI (Jaganathan et al. 2019) catches most of these.

---

## 7. A Python sanity check: caching gnomAD responses

The Week 8 exercises use a small shelve-based cache to avoid hitting the gnomAD API on every run:

```python
from __future__ import annotations
import shelve
import time
from pathlib import Path
from typing import Any

import requests


class GnomadCache:
    """A simple disk-backed cache for gnomAD API responses.

    Keyed by (variant_id, dataset). Use as a context manager.
    """

    def __init__(self, cache_path: Path,
                 rate_limit_sleep: float = 0.2) -> None:
        self.cache_path: Path = cache_path
        self.rate_limit_sleep: float = rate_limit_sleep
        self._db: shelve.Shelf[Any] | None = None

    def __enter__(self) -> "GnomadCache":
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self._db = shelve.open(str(self.cache_path))
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        if self._db is not None:
            self._db.close()
            self._db = None

    def get(self, variant_id: str, dataset: str = "gnomad_r4") -> dict[str, Any]:
        if self._db is None:
            raise RuntimeError("Use GnomadCache as a context manager.")
        key = f"{dataset}:{variant_id}"
        if key in self._db:
            return self._db[key]
        result = query_gnomad(variant_id, dataset=dataset)
        self._db[key] = result
        time.sleep(self.rate_limit_sleep)
        return result
```

The cache file lives at `~/.cache/c10-week08/gnomad.shelve` and persists across runs. The first run of a new VCF makes ~N API calls; subsequent runs read from disk in a few milliseconds total. For the mini-project's 50-variant VCF, the first run takes ~10 seconds (50 * 0.2 sec rate-limit sleep); cached re-runs take < 1 second.

---

## 8. What to remember

- **gnomAD answers "how common." popmax is the AF you compare against the ACMG thresholds (BA1: > 0.05; PM2: < 0.0001).**
- **ClinVar answers "what is known." Trust 3- and 4-star records; flag conflicts; check sub-version drift.**
- **dbSNP gives you the rsID** for cross-reference with the literature. Most of the time you do not need a separate dbSNP query; the rsID is in the VEP CSQ output.
- **Sub-population stratification matters**. A variant common in one sub-population is still common in patients from that sub-population, regardless of the global AF.
- **Pin everything**. Run-info JSON with gnomAD version, ClinVar release, assembly, run date.
- **Cache responses**. The first run of a pipeline pays the API cost; subsequent runs read from disk.
- **Educational use only**. Annotation is mechanical; interpretation requires a clinician.

Continue to **Lecture 3 — ACMG Classification and Pharmacogenomics** to see how the population and clinical evidence combine into the ACMG five-tier classification, which criteria are mechanically computable, and how PharmGKB adds the pharmacogenomics axis on top.
