"""
Exercise 2 - Query ClinVar and gnomAD programmatically.

EDUCATIONAL AND RESEARCH USE ONLY. This file demonstrates the API call
patterns for the free public ClinVar (via NCBI E-utilities) and gnomAD
(via the public GraphQL endpoint) databases. The output is not a
clinical interpretation. Do not use this code to make clinical
decisions. The same disclaimer that opens every Week 8 file applies.

Goal: for each of ~10 hand-picked variants, query ClinVar for the
clinical assertion (CLNSIG, CLNREVSTAT, CLNDN) and gnomAD for the
allele frequency (AF, popmax AF, sub-population AFs). Save the joined
table as a TSV. The exercise covers:

- The polite-retry pattern (exponential backoff on 429/503).
- The per-request rate-limit sleep.
- A local shelve-based response cache.
- Parsing the gnomAD GraphQL JSON response.
- Parsing the NCBI E-utilities JSON response for ClinVar.
- Handling "variant absent from gnomAD" cleanly (this is informative
  for the ACMG PM2 criterion, not an error).

Estimated time: 75 minutes (30 min reading, 30 min implementing
the wrappers, 15 min running and inspecting results).

Acceptance criteria:
- `python exercise-02-query-clinvar-gnomad.py` runs end to end without
  network errors (it will print warnings for variants that 404; this
  is expected).
- `results/exercise-02-clinvar-gnomad.tsv` exists with columns:
  chrom, pos, ref, alt, gnomad_af, gnomad_popmax_af, gnomad_subpop_max,
  clinvar_clnsig, clinvar_clnrevstat, clinvar_clndn, query_status.
- The output table has at least 8 of 10 rows with a non-empty
  gnomad_af and at least 5 of 10 rows with a non-empty clinvar_clnsig.
  (Some demo variants are intentionally rare or absent.)
- The script's `--cache-dir` flag works; running twice with the same
  cache produces identical output and the second run is much faster.

Requirements:
    pip install requests pandas
    (no other dependencies; this is pure HTTP plus pandas.)

What you learn:
- The polite-citizen API pattern: rate-limit sleep, exponential
  backoff, retry on 429/503, fail fast on other errors.
- The shelve-based local cache for idempotent re-runs.
- The gnomAD GraphQL query shape and the response JSON.
- The NCBI E-utilities query shape and the response JSON.
- The difference between "variant not in gnomAD" (informative for
  PM2) and "API error" (need to retry).

Tool versions assumed:
- Python 3.11+
- requests 2.32+
- pandas 2.2+
- gnomAD v4.1 (the default dataset)
- NCBI ClinVar via E-utilities (the public release)

References:
- gnomAD: Karczewski et al. 2020, Nature 581:434
  https://www.nature.com/articles/s41586-020-2308-7
- ClinVar: Landrum et al. 2018, NAR 46:D1062
  https://academic.oup.com/nar/article/46/D1/D1062/4641904
- NCBI E-utilities documentation:
  https://www.ncbi.nlm.nih.gov/books/NBK25500/
- gnomAD API:
  https://gnomad.broadinstitute.org/api
"""

from __future__ import annotations

import argparse
import json
import shelve
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import requests


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Polite-citizen rate limits. gnomAD is ~10 req/sec; NCBI without
# an API key is 3 req/sec. We stay well below both.
GNOMAD_SLEEP_SEC: float = 0.25
NCBI_SLEEP_SEC: float = 0.40
MAX_ATTEMPTS: int = 4
RETRY_BASE_SLEEP_SEC: float = 1.0
REQUEST_TIMEOUT_SEC: float = 30.0

# Default cache location.
DEFAULT_CACHE_DIR: Path = Path.home() / ".cache" / "c10-week08"

# Default gnomAD dataset (v4 combined; GRCh38).
DEFAULT_GNOMAD_DATASET: str = "gnomad_r4"

# The ten demo variants for this exercise. They are intentionally
# chosen to span a range of database presence: well-characterized
# (CFTR delta-F508, BRCA1 hotspots), moderate (a few missense in
# DNA-repair genes), and absent (synthetic coordinates that should
# return empty responses).
# Each entry: (chrom_no_chr_prefix, pos, ref, alt, gene_hint).
DEMO_VARIANTS: list[tuple[str, int, str, str, str]] = [
    ("7",  117559590, "G", "A", "CFTR"),       # near delta-F508 region
    ("17", 43094077,  "G", "A", "BRCA1"),
    ("13", 32398489,  "G", "A", "BRCA2"),
    ("2",  47403534,  "C", "A", "MSH2"),
    ("3",  37001008,  "G", "A", "MLH1"),
    ("5",  112815473, "A", "T", "APC"),
    ("17", 7674220,   "C", "T", "TP53"),
    ("11", 108235829, "G", "A", "ATM"),
    ("15", 48700503,  "G", "A", "FBN1"),
    ("19", 11116107,  "G", "A", "LDLR"),
]


# ----------------------------------------------------------------------
# gnomAD client.
# ----------------------------------------------------------------------

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


def fetch_with_retry(method: str, url: str, *,
                     params: dict[str, Any] | None = None,
                     json_payload: dict[str, Any] | None = None,
                     headers: dict[str, str] | None = None,
                     max_attempts: int = MAX_ATTEMPTS,
                     base_sleep: float = RETRY_BASE_SLEEP_SEC,
                     timeout: float = REQUEST_TIMEOUT_SEC) -> dict[str, Any]:
    """HTTP fetch with exponential backoff on 429 / 503.

    Args:
        method:       "GET" or "POST".
        url:          target URL.
        params:       query-string parameters (GET).
        json_payload: JSON body (POST).
        headers:      optional headers.
        max_attempts: total attempts including the first try.
        base_sleep:   sleep before retry; doubles on each attempt.
        timeout:      per-request timeout.

    Returns:
        Parsed JSON body of the successful response.

    Raises:
        requests.HTTPError on any final non-2xx response.
        requests.RequestException on network failure after retries.
        RuntimeError if all attempts return 429/503.
    """
    last_status: int = 0
    for attempt in range(max_attempts):
        try:
            if method.upper() == "GET":
                response = requests.get(url, params=params, headers=headers, timeout=timeout)
            elif method.upper() == "POST":
                response = requests.post(url, json=json_payload, headers=headers, timeout=timeout)
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")
        except requests.RequestException as exc:
            if attempt == max_attempts - 1:
                raise
            time.sleep(base_sleep * (2 ** attempt))
            continue

        last_status = response.status_code
        if response.status_code in (429, 503):
            sleep_for: float = base_sleep * (2 ** attempt)
            time.sleep(sleep_for)
            continue
        response.raise_for_status()
        return response.json()
    raise RuntimeError(
        f"All {max_attempts} attempts to {url} returned {last_status}. "
        "Consider widening the retry window or checking service status."
    )


def query_gnomad(variant_id: str,
                 dataset: str = DEFAULT_GNOMAD_DATASET) -> dict[str, Any]:
    """Look up a single variant in gnomAD.

    Args:
        variant_id: gnomAD ID in the form "<chrom>-<pos>-<ref>-<alt>"
                    with no "chr" prefix; e.g. "7-117559590-G-A".
        dataset:    gnomAD dataset ID; the default is gnomad_r4 (v4 combined).

    Returns:
        The "variant" subtree of the JSON response (a dict), or an
        empty dict {} if the variant is not in gnomAD. An empty
        response is INFORMATIVE: it supports the ACMG PM2 criterion.
    """
    payload: dict[str, Any] = {
        "query": GNOMAD_QUERY,
        "variables": {"variantId": variant_id, "dataset": dataset},
    }
    data = fetch_with_retry("POST", GNOMAD_API_URL, json_payload=payload)
    variant: dict[str, Any] | None = data.get("data", {}).get("variant")
    return variant or {}


def compute_popmax_af(populations: list[dict[str, Any]],
                      exclude: tuple[str, ...] = ("fin", "asj", "oth", "mid")) -> float:
    """Compute the popmax allele frequency.

    Following gnomAD convention, Finnish (fin), Ashkenazi Jewish (asj),
    Other (oth), and Middle Eastern (mid) sub-populations are excluded
    from popmax because they are bottlenecked or small.

    Args:
        populations: list of {"id": <code>, "ac": <int>, "an": <int>} dicts.
        exclude:     subpop codes to exclude from popmax.

    Returns:
        The max AF over the included subpops, or 0.0 if no eligible
        subpop has data.
    """
    best: float = 0.0
    for entry in populations:
        if entry.get("id") in exclude:
            continue
        ac: int = int(entry.get("ac", 0))
        an: int = int(entry.get("an", 0))
        if an == 0:
            continue
        af: float = ac / an
        if af > best:
            best = af
    return best


def compute_subpop_max_string(populations: list[dict[str, Any]]) -> str:
    """Build a "max-subpop AF; max-subpop-ID" string for the report.

    Returns "" if no population has data.
    """
    best_id: str = ""
    best_af: float = 0.0
    for entry in populations:
        an: int = int(entry.get("an", 0))
        if an == 0:
            continue
        af: float = int(entry.get("ac", 0)) / an
        if af > best_af:
            best_af = af
            best_id = str(entry.get("id", ""))
    if not best_id:
        return ""
    return f"{best_af:.6f} ({best_id})"


# ----------------------------------------------------------------------
# ClinVar client (via NCBI E-utilities).
# ----------------------------------------------------------------------

def query_clinvar_by_position(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    """Look up ClinVar records overlapping a coordinate.

    Uses E-utilities esearch + esummary against the ClinVar database.
    Returns the first matching record's metadata, or {} if no match.

    Args:
        chrom: chromosome (no "chr" prefix), e.g. "17".
        pos:   1-based position.
        ref:   reference allele.
        alt:   alternate allele.

    Returns:
        dict with clinvar_clnsig, clinvar_clnrevstat, clinvar_clndn,
        clinvar_variation_id, clinvar_title.
        Empty dict {} if no match.
    """
    # 1. esearch: find ClinVar records overlapping the coordinate.
    search_params: dict[str, Any] = {
        "db": "clinvar",
        "term": f'{chrom}[Chromosome] AND {pos}:{pos}[Base Position for Assembly GRCh38]',
        "retmode": "json",
        "retmax": 20,
    }
    search_url: str = f"{NCBI_EUTILS_URL}/esearch.fcgi"
    try:
        search_resp: dict[str, Any] = fetch_with_retry("GET", search_url, params=search_params)
    except requests.HTTPError:
        return {}
    id_list: list[str] = search_resp.get("esearchresult", {}).get("idlist", []) or []
    if not id_list:
        return {}

    # 2. esummary: pull metadata for the matched IDs.
    summary_params: dict[str, Any] = {
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json",
    }
    summary_url: str = f"{NCBI_EUTILS_URL}/esummary.fcgi"
    try:
        summary_resp: dict[str, Any] = fetch_with_retry("GET", summary_url, params=summary_params)
    except requests.HTTPError:
        return {}

    result: dict[str, Any] = summary_resp.get("result", {})
    for vid in id_list:
        record: dict[str, Any] = result.get(vid, {})
        if not record:
            continue
        # Filter to records that match the ref/alt alleles. E-utilities
        # returns coordinate-overlapping records but not always exact
        # allele matches.
        variation_set: list[dict[str, Any]] = record.get("variation_set", []) or []
        match: dict[str, Any] | None = None
        for entry in variation_set:
            allele_info: list[dict[str, Any]] = entry.get("variation_loc", []) or []
            for loc in allele_info:
                if loc.get("assembly_name", "").upper().startswith("GRCH38"):
                    if str(loc.get("chr", "")).lstrip("chr") == chrom:
                        match = entry
                        break
            if match is not None:
                break
        if match is None:
            continue
        clinical_significance: dict[str, Any] = record.get("germline_classification", {}) or {}
        return {
            "clinvar_variation_id": vid,
            "clinvar_clnsig": clinical_significance.get("description", ""),
            "clinvar_clnrevstat": clinical_significance.get("review_status", ""),
            "clinvar_clndn": ";".join(t.get("trait_name", "")
                                       for t in clinical_significance.get("trait_set", []) or []),
            "clinvar_title": record.get("title", ""),
        }
    return {}


# ----------------------------------------------------------------------
# Cache wrapper.
# ----------------------------------------------------------------------

@dataclass
class QueryCaches:
    """Bundle the gnomAD and ClinVar disk caches."""
    gnomad_path: Path
    clinvar_path: Path

    def open(self) -> tuple[shelve.Shelf[Any], shelve.Shelf[Any]]:
        self.gnomad_path.parent.mkdir(parents=True, exist_ok=True)
        gnomad = shelve.open(str(self.gnomad_path))
        clinvar = shelve.open(str(self.clinvar_path))
        return gnomad, clinvar


def get_gnomad_cached(cache: shelve.Shelf[Any], variant_id: str,
                      dataset: str = DEFAULT_GNOMAD_DATASET) -> dict[str, Any]:
    """Cache-aware gnomAD lookup."""
    key: str = f"{dataset}:{variant_id}"
    if key in cache:
        return cache[key]
    try:
        data = query_gnomad(variant_id, dataset=dataset)
    except Exception as exc:  # noqa: BLE001
        # Cache nothing on hard error; allow retry next run.
        print(f"  [gnomad] error for {variant_id}: {exc}", file=sys.stderr)
        return {}
    cache[key] = data
    time.sleep(GNOMAD_SLEEP_SEC)
    return data


def get_clinvar_cached(cache: shelve.Shelf[Any], chrom: str, pos: int,
                       ref: str, alt: str) -> dict[str, Any]:
    """Cache-aware ClinVar lookup."""
    key: str = f"GRCh38:{chrom}:{pos}:{ref}:{alt}"
    if key in cache:
        return cache[key]
    try:
        data = query_clinvar_by_position(chrom, pos, ref, alt)
    except Exception as exc:  # noqa: BLE001
        print(f"  [clinvar] error for {chrom}:{pos} {ref}>{alt}: {exc}", file=sys.stderr)
        return {}
    cache[key] = data
    time.sleep(NCBI_SLEEP_SEC)
    return data


# ----------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------

def build_query_row(chrom: str, pos: int, ref: str, alt: str, gene_hint: str,
                     gnomad_cache: shelve.Shelf[Any],
                     clinvar_cache: shelve.Shelf[Any],
                     dataset: str = DEFAULT_GNOMAD_DATASET) -> dict[str, Any]:
    """Query gnomAD + ClinVar for one variant; return a flat dict.

    Output columns:
        chrom, pos, ref, alt, gene_hint,
        gnomad_variant_id, gnomad_rsid, gnomad_af, gnomad_popmax_af,
        gnomad_subpop_max, gnomad_filters,
        clinvar_variation_id, clinvar_clnsig, clinvar_clnrevstat,
        clinvar_clndn, query_status.
    """
    variant_id: str = f"{chrom}-{pos}-{ref}-{alt}"
    g: dict[str, Any] = get_gnomad_cached(gnomad_cache, variant_id, dataset=dataset)
    exome: dict[str, Any] = g.get("exome") or {}
    populations: list[dict[str, Any]] = exome.get("populations", []) or []
    popmax_af: float = compute_popmax_af(populations) if populations else 0.0
    subpop_max: str = compute_subpop_max_string(populations) if populations else ""
    filters: list[str] = exome.get("filters", []) or []

    cv: dict[str, Any] = get_clinvar_cached(clinvar_cache, chrom, pos, ref, alt)

    status_parts: list[str] = []
    if g:
        status_parts.append("gnomad_found")
    else:
        status_parts.append("gnomad_absent")
    if cv:
        status_parts.append("clinvar_found")
    else:
        status_parts.append("clinvar_absent")
    status: str = ";".join(status_parts)

    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene_hint": gene_hint,
        "gnomad_variant_id": g.get("variant_id", "") if g else "",
        "gnomad_rsid": g.get("rsid", "") if g else "",
        "gnomad_af": exome.get("af", "") if g else "",
        "gnomad_popmax_af": popmax_af if g else "",
        "gnomad_subpop_max": subpop_max,
        "gnomad_filters": ";".join(filters),
        "clinvar_variation_id": cv.get("clinvar_variation_id", ""),
        "clinvar_clnsig": cv.get("clinvar_clnsig", ""),
        "clinvar_clnrevstat": cv.get("clinvar_clnrevstat", ""),
        "clinvar_clndn": cv.get("clinvar_clndn", ""),
        "query_status": status,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query ClinVar and gnomAD for a list of variants. "
                    "Educational and research use only. Not for clinical decision-making."
    )
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).parent / "results" / "exercise-02-clinvar-gnomad.tsv",
                        help="Output TSV path.")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR,
                        help="Local cache directory (gnomAD + ClinVar shelves).")
    parser.add_argument("--dataset", type=str, default=DEFAULT_GNOMAD_DATASET,
                        choices=["gnomad_r2_1", "gnomad_r3", "gnomad_r4"],
                        help="gnomAD dataset to query.")
    parser.add_argument("--offline-demo", action="store_true",
                        help="Skip network calls; use bundled demo responses. Useful for CI.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_path: Path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.offline_demo:
        return _run_offline_demo(out_path)

    caches = QueryCaches(
        gnomad_path=args.cache_dir / "gnomad.shelve",
        clinvar_path=args.cache_dir / "clinvar.shelve",
    )
    gnomad_cache, clinvar_cache = caches.open()
    try:
        rows: list[dict[str, Any]] = []
        for chrom, pos, ref, alt, gene_hint in DEMO_VARIANTS:
            print(f"[exercise-02] querying {chrom}-{pos}-{ref}-{alt} ({gene_hint}) ...")
            row = build_query_row(chrom, pos, ref, alt, gene_hint,
                                  gnomad_cache, clinvar_cache, dataset=args.dataset)
            rows.append(row)
    finally:
        gnomad_cache.close()
        clinvar_cache.close()

    df: pd.DataFrame = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)

    print()
    print(f"[exercise-02] wrote {out_path} with {len(df)} rows.")
    print()
    print("[exercise-02] per-row query status:")
    for _, row in df.iterrows():
        print(f"  {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}  status={row['query_status']}")
    return 0


def _run_offline_demo(out_path: Path) -> int:
    """Produce a synthetic table for CI / disconnected environments.

    Returns 0 on success. Does no network I/O.
    """
    rows: list[dict[str, Any]] = []
    for chrom, pos, ref, alt, gene_hint in DEMO_VARIANTS:
        rows.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "gene_hint": gene_hint,
            "gnomad_variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "gnomad_rsid": "",
            "gnomad_af": 0.0001 if gene_hint in ("BRCA1", "BRCA2", "MSH2") else "",
            "gnomad_popmax_af": 0.0002 if gene_hint in ("BRCA1", "BRCA2", "MSH2") else "",
            "gnomad_subpop_max": "0.000200 (nfe)" if gene_hint in ("BRCA1", "BRCA2", "MSH2") else "",
            "gnomad_filters": "",
            "clinvar_variation_id": "",
            "clinvar_clnsig": "Likely_pathogenic" if gene_hint in ("BRCA1", "BRCA2") else "",
            "clinvar_clnrevstat": "criteria_provided,_single_submitter" if gene_hint in ("BRCA1", "BRCA2") else "",
            "clinvar_clndn": "Hereditary breast and ovarian cancer syndrome" if gene_hint in ("BRCA1", "BRCA2") else "",
            "query_status": "offline_demo",
        })
    df: pd.DataFrame = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[exercise-02] offline-demo wrote {out_path} with {len(df)} rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
