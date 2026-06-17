"""
Exercise 3 - Annotate a VCF with VEP via the Ensembl REST API.

Goal: take a small VCF (the output of Exercise 1, or the bundled example),
send each variant to the Ensembl VEP REST API at
https://rest.ensembl.org/vep/<species>/region/, parse the JSON response
to extract the consequence(s) for each variant on each overlapping
transcript, and tabulate the per-variant consequence summary.

The REST API is appropriate for exercises (few variants, no offline
cache to install). For production use against a full-genome VCF, use
the offline cache: `vep --cache --dir_cache vep_cache/`. The cache
delivers ~10,000x the throughput.

Estimated time: 60 minutes (10 minutes setup, 30 minutes implementation,
20 minutes interpretation).

Acceptance criteria:
- `python exercise-03-vep-annotate.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The output table at results/vep_annotations.tsv exists with at least
  one row per input variant.
- You implemented four functions: `read_variants_from_vcf`,
  `query_vep_rest`, `parse_vep_response`, and `tabulate_annotations`.

Requirements:
    pip install requests pandas
    (and an internet connection to reach https://rest.ensembl.org/)

What you learn:
- How the Ensembl VEP REST API works (HTTP GET against
  https://rest.ensembl.org/vep/<species>/region/<chr>:<pos>-<pos>/<alt>?...).
- How VEP encodes a variant's effect on a transcript: consequence terms
  (missense_variant, synonymous_variant, etc.), IMPACT category (LOW,
  MODERATE, HIGH, MODIFIER), HGVS notation (c.1547C>T, p.Ala516Val).
- The full list of ~30 consequence terms in the Sequence Ontology.
- Why production VEP uses the offline cache and not the REST API for
  whole-genome work (REST is ~15 req/sec rate-limited; offline is
  ~10000 req/sec).

TO COMPLETE: implement the four functions below. Run the file; all
assertions must pass.

Tool versions assumed:
- Python 3.11+
- requests 2.31+
- pandas 2+

Note on networking: the REST API may be slow or unavailable. If you get
HTTP 503 or timeouts, retry. If consistently failing, fall back to the
offline VEP cache: `vep_install --AUTO c --SPECIES homo_sapiens --CACHEDIR
vep_cache/` then `vep --cache --dir_cache vep_cache/ -i input.vcf -o output.vcf`.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, NamedTuple


# A small bundled example: human variants for which we know the
# canonical consequence. These come from the published literature
# (BRCA1 p.E1346K, p53 R175H, CFTR F508del - well-known disease loci).
# We bundle them so this exercise runs even if Exercise 1's VCF is empty.
EXAMPLE_VARIANTS_TSV = """\
chrom\tpos\tref\talt\texpected_gene
17\t43093843\tG\tA\tBRCA1
17\t7674220\tG\tA\tTP53
7\t117559590\tCTT\tC\tCFTR
"""

# REST API base URL.
VEP_REST_URL = "https://rest.ensembl.org"

# The Ensembl species identifier for the variants above.
VEP_SPECIES = "homo_sapiens"

# Output location for the annotation table.
OUTPUT_TSV = Path(__file__).parent / "results" / "vep_annotations.tsv"


class VariantRow(NamedTuple):
    """A minimal variant record for VEP input."""
    chrom: str
    pos: int
    ref: str
    alt: str


class VEPAnnotation(NamedTuple):
    """A single per-transcript VEP consequence annotation."""
    chrom: str
    pos: int
    ref: str
    alt: str
    transcript_id: str
    gene_symbol: str
    consequence: str  # e.g. 'missense_variant'
    impact: str       # 'LOW' / 'MODERATE' / 'HIGH' / 'MODIFIER'
    hgvs_c: str       # 'c.<position><ref>><alt>', if available
    hgvs_p: str       # 'p.<aa><pos><aa>', if available
    canonical: bool   # True if VEP marks this transcript as canonical


def read_variants_from_vcf(vcf_path: Path) -> list[VariantRow]:
    """Read variants from a VCF (plain text or bgzipped).

    Returns a list of VariantRow with one entry per VCF variant line.
    Skips header lines (lines starting with '#').

    If the VCF is bgzipped (.vcf.gz), open with gzip.open. Otherwise,
    open as plain text.

    Hint: the VCF columns are tab-separated. CHROM is column 0, POS is
    column 1 (1-based, int), REF is column 3, ALT is column 4 (split on
    commas if multi-allelic; for this exercise we only handle the first
    ALT per record).
    """
    import gzip

    variants = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    mode = "rt" if str(vcf_path).endswith(".gz") else "r"

    with opener(vcf_path, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            # TODO: split the line on tabs.
            # TODO: extract chrom, pos (int), ref, alts (split on comma).
            # TODO: append a VariantRow per alt allele.
            # TODO: implement and remove the placeholder below.
            raise NotImplementedError("Parse the VCF variant line")

    return variants


def query_vep_rest(
    species: str, chrom: str, pos: int, ref: str, alt: str,
    timeout: float = 10.0, retries: int = 3,
) -> dict[str, Any] | None:
    """Query the Ensembl VEP REST API for a single variant.

    The REST endpoint is:
        GET https://rest.ensembl.org/vep/<species>/region/<chr>:<start>-<end>/<alt>

    For a SNP, start == end == pos. For an insertion, start = pos+1,
    end = pos. For a deletion, start = pos+1, end = pos+len(ref)-1
    and alt is encoded as '-'.

    For this exercise, we only handle SNPs and small indels in the
    simple form. The function returns the parsed JSON response (a list
    of one dict per variant; each dict has a 'transcript_consequences'
    field with a list of per-transcript annotations).

    Returns None on persistent failure (after retries).

    Hint: use `requests.get(url, headers={'Content-Type': 'application/json'})`.
    Sleep briefly between retries (the API has a rate limit of ~15 req/sec).
    """
    import requests

    # Build the URL for a SNP. For indels, the encoding is different;
    # see https://rest.ensembl.org/documentation/info/vep_region_get.
    # For this exercise we restrict to SNPs and 1-bp insertions/deletions
    # in their simplest form.
    if len(ref) == 1 and len(alt) == 1:
        # SNP
        url = f"{VEP_REST_URL}/vep/{species}/region/{chrom}:{pos}-{pos}/{alt}"
    elif len(ref) > len(alt) and ref.startswith(alt):
        # Deletion. REF=CTT, ALT=C means delete the trailing 'TT' starting
        # at pos+1.
        start = pos + len(alt)
        end = pos + len(ref) - 1
        url = f"{VEP_REST_URL}/vep/{species}/region/{chrom}:{start}-{end}/-"
    elif len(alt) > len(ref) and alt.startswith(ref):
        # Insertion. REF=A, ALT=AGT means insert 'GT' after pos.
        start = pos + len(ref)
        end = pos + len(ref) - 1  # zero-length range -> insertion in Ensembl
        ins_seq = alt[len(ref):]
        url = f"{VEP_REST_URL}/vep/{species}/region/{chrom}:{start}-{end}/{ins_seq}"
    else:
        # MNP or complex indel; not handled in this exercise.
        return None

    headers = {"Content-Type": "application/json"}

    for attempt in range(retries):
        # TODO: try requests.get(url, headers=headers, timeout=timeout).
        # TODO: if status_code == 200, return response.json().
        # TODO: if status_code in (429, 503), sleep then retry.
        # TODO: on persistent failure or other status, return None.
        raise NotImplementedError("Make the HTTP request and parse JSON")

    return None


def parse_vep_response(
    response: dict[str, Any] | list[dict[str, Any]] | None,
    chrom: str, pos: int, ref: str, alt: str,
) -> list[VEPAnnotation]:
    """Parse a VEP REST response into a list of VEPAnnotation rows.

    The REST response shape is a list of one dict per input variant. We
    pass one variant at a time, so we expect a 1-element list. Each dict
    has a 'transcript_consequences' key, which is a list of per-transcript
    consequence dicts. Each consequence dict has keys like:

        'transcript_id': 'ENST00000357654'
        'gene_symbol':   'BRCA1'
        'consequence_terms': ['missense_variant']  (a list; can have multiple)
        'impact':        'MODERATE'
        'hgvsc':         'ENST00000357654.7:c.4034G>A'
        'hgvsp':         'ENSP00000350283.3:p.Glu1346Lys'
        'canonical':     1  (or absent if not canonical)

    For each consequence in the response, build a VEPAnnotation row.
    If the response is None or empty, return an empty list.

    Hint: handle the missing-key case gracefully — REST responses
    sometimes omit fields. Use `consequence.get('hgvsp', '')` rather
    than `consequence['hgvsp']`.
    """
    annotations: list[VEPAnnotation] = []
    if response is None:
        return annotations
    if isinstance(response, list):
        if not response:
            return annotations
        records = response
    else:
        records = [response]

    for rec in records:
        # TODO: iterate over rec.get("transcript_consequences", []).
        # TODO: for each consequence, build a VEPAnnotation with the
        #       fields named in the docstring above. consequence_terms
        #       is a list; join with ',' or just take the first.
        raise NotImplementedError("Parse the VEP response into annotations")

    return annotations


def tabulate_annotations(annotations: list[VEPAnnotation], out_path: Path) -> None:
    """Write annotations to a tab-separated table.

    The columns: chrom, pos, ref, alt, transcript_id, gene_symbol,
    consequence, impact, hgvs_c, hgvs_p, canonical.

    Use pandas.DataFrame.to_csv(sep='\\t', index=False) for the write.

    If annotations is empty, write a header-only file.
    """
    import pandas as pd

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # TODO: convert annotations (a list of VEPAnnotation namedtuples)
    #       to a pandas DataFrame and write as TSV.
    raise NotImplementedError("Write the annotation table as TSV")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-vep-annotate.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Try to import requests; if missing, fail with a helpful message.
    try:
        import requests  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "[exercise-03] requests is not installed. Install with:\n"
            "    pip install requests==2.31\n"
            "and re-run."
        ) from exc

    try:
        import pandas as pd  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "[exercise-03] pandas is not installed. Install with:\n"
            "    pip install pandas\n"
            "and re-run."
        ) from exc

    # Build the input variants from the bundled example.
    variants = []
    for line in EXAMPLE_VARIANTS_TSV.splitlines()[1:]:  # skip header
        chrom, pos_str, ref, alt, expected_gene = line.split("\t")
        variants.append((VariantRow(chrom, int(pos_str), ref, alt), expected_gene))

    # If Exercise 1 produced a VCF, also try reading that. (Not assertion-
    # checked; just demonstrates the read path.)
    here = Path(__file__).parent
    exercise_1_vcf = here / "calls" / "lambda.norm.vcf.gz"
    if exercise_1_vcf.exists():
        print(f"[exercise-03] Found Exercise 1 VCF at {exercise_1_vcf}.")
        try:
            ex1_variants = read_variants_from_vcf(exercise_1_vcf)
            print(f"[exercise-03] Read {len(ex1_variants)} variants from Exercise 1.")
        except NotImplementedError:
            # The student has not implemented read_variants_from_vcf yet.
            # This is non-fatal; we proceed with the bundled examples.
            pass

    # Annotate each example variant.
    all_annotations: list[VEPAnnotation] = []

    print()
    print("[exercise-03] Annotating example variants via Ensembl VEP REST API ...")
    print("[exercise-03] (this requires an internet connection)")
    print()

    for variant, expected_gene in variants:
        print(f"  {variant.chrom}:{variant.pos} {variant.ref}>{variant.alt} "
              f"(expected gene: {expected_gene}) ...")
        # Sleep briefly between requests to respect the rate limit.
        time.sleep(0.2)
        response = query_vep_rest(
            VEP_SPECIES, variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        annotations = parse_vep_response(
            response, variant.chrom, variant.pos, variant.ref, variant.alt,
        )

        # Assertion: at least one annotation should come back per variant
        # (since each variant is in a well-known gene with multiple
        # transcripts).
        assert len(annotations) >= 1, (
            f"VEP returned no annotations for {variant.chrom}:{variant.pos} "
            f"{variant.ref}>{variant.alt} — is the REST API reachable?"
        )

        # Assertion: at least one annotation should match the expected gene.
        gene_symbols = {a.gene_symbol for a in annotations}
        assert expected_gene in gene_symbols, (
            f"Expected gene {expected_gene!r} not in VEP annotations: "
            f"{gene_symbols}"
        )

        all_annotations.extend(annotations)

    # Tabulate.
    print()
    print(f"[exercise-03] Total annotations across all variants: "
          f"{len(all_annotations)}")
    tabulate_annotations(all_annotations, OUTPUT_TSV)
    assert OUTPUT_TSV.exists(), f"output TSV not at {OUTPUT_TSV}"
    assert OUTPUT_TSV.stat().st_size > 0, (
        f"output TSV at {OUTPUT_TSV} is empty"
    )

    # Field summary.
    import pandas as pd
    df = pd.read_csv(OUTPUT_TSV, sep="\t")
    print()
    print("[exercise-03] Annotation summary:")
    print(f"  Total rows:       {len(df)}")
    print(f"  Unique genes:     {df['gene_symbol'].nunique()}")
    print(f"  Unique consequences: {df['consequence'].nunique()}")
    print(f"  Canonical rows:   {df['canonical'].sum()}")
    print()
    print("  By consequence:")
    print(df["consequence"].value_counts().to_string())
    print()
    print("  By IMPACT:")
    print(df["impact"].value_counts().to_string())
    print()

    # Assert that at least one consequence is HIGH or MODERATE impact
    # (the three example variants are all known damaging variants).
    high_or_moderate = df[df["impact"].isin(["HIGH", "MODERATE"])]
    assert len(high_or_moderate) >= 1, (
        "expected at least one HIGH/MODERATE-impact annotation among the "
        "example variants (BRCA1, TP53, CFTR are damaging)"
    )

    print("[exercise-03] All assertions passed.")
    print(f"[exercise-03] Annotation table saved to {OUTPUT_TSV}")
    print("[exercise-03] You now have the end-to-end VCF -> VEP annotation")
    print("[exercise-03] pipeline running. Continue to the challenge or")
    print("[exercise-03] mini-project.")
