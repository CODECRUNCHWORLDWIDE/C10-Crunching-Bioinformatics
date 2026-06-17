"""
Exercise 3 - Parse BLAST output (tabular and XML) into pandas.

Goal: take a pair of pre-recorded BLAST output files (one in `-outfmt 6`
tabular form, one in `-outfmt 5` XML form) and write the parsing layer
that the mini-project's classifier will sit on top of. Two functions
that read each format into the same canonical DataFrame; a third that
joins them with a taxonomy lineage; a fourth that picks the top hit
per query under an E-value cutoff.

Estimated time: 50 minutes.

Acceptance criteria:
- `python exercise-03-parse-blast-output.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You implemented four functions: `read_tabular`, `read_xml`,
  `filter_top_hit_per_query`, and `attach_lineage_stub`.
- The script does NOT make any network calls. It runs entirely on
  the synthetic test data embedded below.
- The two readers produce DataFrames with the same columns in the
  same order: ("qseqid", "sseqid", "pident", "length", "evalue",
  "bitscore", "identities", "align_length").

Requirements:
    python -m pip install biopython==1.83 pandas

What you learn:
- How to handle the two most-common BLAST output formats from one
  consistent pandas API.
- How to compute percent identity and HSP fields from the XML object
  graph (vs the tabular shortcut where it is a column).
- How to write "top hit per query" reductions in pandas without
  shooting yourself in the foot on the SettingWithCopyWarning.
- The shape of a taxonomy-lineage join, in stub form, ready for
  Thursday's real Entrez call.

TO COMPLETE: implement the four functions below.

Reference accessions in the synthetic test data:
    NR_117741.1  Staphylococcus aureus 16S rRNA  (taxid 1280)
    NR_074549.1  Escherichia coli K-12 MG1655 16S rRNA  (taxid 511145)
    NR_119213.1  Bacillus subtilis 168 16S rRNA  (taxid 224308)
"""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
from Bio.Blast import NCBIXML


# Canonical column order. Both readers must return these columns
# in this order.
CANONICAL_COLUMNS: tuple[str, ...] = (
    "qseqid", "sseqid", "pident", "length",
    "evalue", "bitscore", "identities", "align_length",
)


# A small, hand-crafted `-outfmt 6` tabular output. Three queries, three
# subjects, six rows. Real BLAST runs may have many more HSPs per
# (query, subject) but this is enough to exercise the parser.
SYNTHETIC_TABULAR = """\
query_001\tNR_117741.1\t99.10\t1421\t12\t1\t1\t1421\t1\t1421\t0.0\t2598
query_001\tNR_074549.1\t87.30\t1390\t160\t12\t10\t1395\t5\t1390\t1.2e-67\t245
query_002\tNR_074549.1\t99.50\t1535\t7\t1\t1\t1535\t1\t1535\t0.0\t2780
query_002\tNR_117741.1\t87.70\t1380\t150\t13\t12\t1390\t8\t1385\t3.4e-65\t235
query_003\tNR_119213.1\t98.75\t1432\t16\t2\t1\t1432\t1\t1432\t0.0\t2510
query_003\tNR_074549.1\t85.60\t1340\t180\t12\t8\t1345\t5\t1340\t4.2e-50\t200
"""


# A small synthetic BLAST XML. Two queries, two alignments per query.
# Hand-written to match the structure that Bio.Blast.NCBIXML emits.
SYNTHETIC_XML = """\
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"
  "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.15.0+</BlastOutput_version>
  <BlastOutput_reference>Altschul et al. 1997</BlastOutput_reference>
  <BlastOutput_db>16S_ribosomal_RNA</BlastOutput_db>
  <BlastOutput_query-ID>query_001</BlastOutput_query-ID>
  <BlastOutput_query-def>query_001 16S environmental</BlastOutput_query-def>
  <BlastOutput_query-len>1421</BlastOutput_query-len>
  <BlastOutput_param><Parameters>
    <Parameters_expect>1e-50</Parameters_expect>
    <Parameters_sc-match>1</Parameters_sc-match>
    <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
    <Parameters_gap-open>5</Parameters_gap-open>
    <Parameters_gap-extend>2</Parameters_gap-extend>
    <Parameters_filter>L</Parameters_filter>
  </Parameters></BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>query_001</Iteration_query-ID>
      <Iteration_query-def>query_001 16S environmental</Iteration_query-def>
      <Iteration_query-len>1421</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>NR_117741.1 Staphylococcus aureus subsp. aureus 16S rRNA</Hit_def>
          <Hit_accession>NR_117741.1</Hit_accession>
          <Hit_len>1547</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>2598</Hsp_bit-score>
              <Hsp_score>1407</Hsp_score>
              <Hsp_evalue>0</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>1421</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>1421</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>1408</Hsp_identity>
              <Hsp_positive>1408</Hsp_positive>
              <Hsp_gaps>1</Hsp_gaps>
              <Hsp_align-len>1421</Hsp_align-len>
              <Hsp_qseq>ACGT</Hsp_qseq>
              <Hsp_hseq>ACGT</Hsp_hseq>
              <Hsp_midline>||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
          <Hit_def>NR_074549.1 Escherichia coli K-12 MG1655 16S rRNA</Hit_def>
          <Hit_accession>NR_074549.1</Hit_accession>
          <Hit_len>1541</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>245</Hsp_bit-score>
              <Hsp_score>132</Hsp_score>
              <Hsp_evalue>1.2e-67</Hsp_evalue>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>1395</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>1390</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>1213</Hsp_identity>
              <Hsp_positive>1213</Hsp_positive>
              <Hsp_gaps>12</Hsp_gaps>
              <Hsp_align-len>1390</Hsp_align-len>
              <Hsp_qseq>ACGT</Hsp_qseq>
              <Hsp_hseq>ACGT</Hsp_hseq>
              <Hsp_midline>||..</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


# A toy taxonomy map: accession -> (taxid, full lineage string).
# In the real pipeline (mini-project) this comes from Bio.Entrez.efetch
# against db="taxonomy". For the exercise we hard-code it.
TAXONOMY_STUB: dict[str, tuple[int, str]] = {
    "NR_117741.1": (
        1280,
        "Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcaceae; "
        "Staphylococcus; Staphylococcus aureus",
    ),
    "NR_074549.1": (
        511145,
        "Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; "
        "Enterobacteriaceae; Escherichia; Escherichia coli",
    ),
    "NR_119213.1": (
        224308,
        "Bacteria; Firmicutes; Bacilli; Bacillales; Bacillaceae; "
        "Bacillus; Bacillus subtilis",
    ),
}


def read_tabular(tsv_text: str) -> pd.DataFrame:
    """Parse a `-outfmt 6` tabular BLAST output into the canonical DataFrame.

    Input is the TSV text (e.g., the contents of a `blastn -out` file
    or the SYNTHETIC_TABULAR constant above).

    Returns a DataFrame with CANONICAL_COLUMNS, in that order.

    For tabular output, `identities` and `align_length` are derived:
        align_length = `length` column from the tabular output
        identities   = round(pident / 100 * align_length)
    (this is the canonical reconstruction; the tabular format does not
    emit identities as a separate column unless you pass a custom
    -outfmt template.)

    Hint: parse with pandas first into the 12-column standard layout,
    then derive `identities` and reorder/select the eight columns
    listed in CANONICAL_COLUMNS.
    """
    twelve_cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore",
    ]
    # TODO: pd.read_csv on io.StringIO(tsv_text). Derive identities.
    # Re-select to CANONICAL_COLUMNS and return.
    raise NotImplementedError("Parse the tabular BLAST output")


def read_xml(xml_text: str) -> pd.DataFrame:
    """Parse a `-outfmt 5` XML BLAST output into the canonical DataFrame.

    Input is the XML text. Returns a DataFrame with CANONICAL_COLUMNS,
    one row per HSP per alignment (so a single query with two
    alignments each with one HSP yields two rows).

    Hint: use `Bio.Blast.NCBIXML.parse(io.StringIO(xml_text))` to get
    an iterator of BlastRecord objects. For each record, iterate
    `record.alignments`, then each alignment's `.hsps`. Compute
    `pident = identities / align_length * 100`.

    The `qseqid` should be the record's `.query` attribute, split at
    the first whitespace to drop the optional description.
    """
    rows: list[dict] = []
    # TODO: parse the XML via NCBIXML.parse and accumulate rows.
    # Each row should have keys matching CANONICAL_COLUMNS.
    raise NotImplementedError("Parse the XML BLAST output")

    return pd.DataFrame(rows, columns=list(CANONICAL_COLUMNS))


def filter_top_hit_per_query(
    df: pd.DataFrame,
    *,
    evalue_cutoff: float = 1e-10,
) -> pd.DataFrame:
    """Reduce to one row per query: the hit with the lowest E-value.

    Drops queries with no hit below `evalue_cutoff`. Returns a DataFrame
    indexed by `qseqid` with the remaining canonical columns.

    Hint: filter, then `.sort_values(["qseqid", "evalue"])`, then
    `.groupby("qseqid").first()`. Use `.copy()` to avoid pandas's
    SettingWithCopyWarning in downstream code.
    """
    # TODO: filter by evalue_cutoff, sort, groupby first, return.
    raise NotImplementedError("Reduce to top hit per query")


def attach_lineage_stub(
    df: pd.DataFrame,
    tax_map: dict[str, tuple[int, str]] = TAXONOMY_STUB,
) -> pd.DataFrame:
    """Join the top-hit DataFrame with a taxonomy lineage (stub version).

    For each row, look up `sseqid` in `tax_map` and add two columns:
    `taxid` (int) and `lineage` (str). If the accession is missing
    from `tax_map`, fill with NaN / empty string.

    This is the stub the mini-project replaces with a real
    `Bio.Entrez.efetch(db="taxonomy", ...)` call.
    """
    # TODO: map sseqid -> taxid, sseqid -> lineage. Assign as columns.
    raise NotImplementedError("Attach the taxonomy lineage stub")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-parse-blast-output.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("[exercise-03] Step 1 - read the synthetic tabular output ...")
    tab_df = read_tabular(SYNTHETIC_TABULAR)
    assert list(tab_df.columns) == list(CANONICAL_COLUMNS), (
        f"tabular reader produced wrong columns: {list(tab_df.columns)}"
    )
    # Six rows: three queries x two hits each in the synthetic data.
    assert len(tab_df) == 6, f"tabular row count {len(tab_df)} != 6"

    # Sanity-check derived `identities`. For query_001 vs NR_117741.1:
    # pident=99.10, length=1421 -> identities should be round(1408.31) = 1408.
    row = tab_df[
        (tab_df["qseqid"] == "query_001") & (tab_df["sseqid"] == "NR_117741.1")
    ].iloc[0]
    assert row["identities"] == 1408, (
        f"derived identities = {row['identities']}; expected 1408"
    )

    print("[exercise-03] Step 2 - read the synthetic XML output ...")
    xml_df = read_xml(SYNTHETIC_XML)
    assert list(xml_df.columns) == list(CANONICAL_COLUMNS), (
        f"XML reader produced wrong columns: {list(xml_df.columns)}"
    )
    # Two alignments x one HSP each = 2 rows.
    assert len(xml_df) == 2, f"XML row count {len(xml_df)} != 2"

    # The top XML hit is NR_117741.1 with identity 1408/1421 = 99.08%
    # (allowing some float wiggle).
    top_xml = xml_df.sort_values("evalue").iloc[0]
    assert top_xml["sseqid"] == "NR_117741.1"
    assert abs(top_xml["pident"] - 99.08) < 0.1, (
        f"top XML pident = {top_xml['pident']}; expected ~99.08"
    )

    print("[exercise-03] Step 3 - reduce tabular to top hit per query ...")
    top_per_query = filter_top_hit_per_query(tab_df, evalue_cutoff=1e-10)
    assert len(top_per_query) == 3, (
        f"expected 3 queries after top-hit reduction; got {len(top_per_query)}"
    )
    # query_001's top hit must be NR_117741.1 (S. aureus).
    assert top_per_query.loc["query_001", "sseqid"] == "NR_117741.1"
    # query_002's top hit must be NR_074549.1 (E. coli).
    assert top_per_query.loc["query_002", "sseqid"] == "NR_074549.1"
    # query_003's top hit must be NR_119213.1 (B. subtilis).
    assert top_per_query.loc["query_003", "sseqid"] == "NR_119213.1"

    print("[exercise-03] Step 4 - attach the taxonomy lineage stub ...")
    classified = attach_lineage_stub(top_per_query)
    assert "taxid" in classified.columns
    assert "lineage" in classified.columns

    # Spot-check one row.
    s_aureus_row = classified.loc["query_001"]
    assert s_aureus_row["taxid"] == 1280
    assert "Staphylococcus" in s_aureus_row["lineage"]
    assert "Bacteria" in s_aureus_row["lineage"]

    print()
    print("[exercise-03] Final classified table:")
    print(classified[["sseqid", "pident", "evalue", "taxid", "lineage"]])

    # Filter cutoff behavior: a strict cutoff should drop the weak
    # second hits, leaving the three strong top hits intact.
    print()
    print("[exercise-03] Step 5 - check filter cutoff behavior ...")
    strict = filter_top_hit_per_query(tab_df, evalue_cutoff=1e-200)
    # The three "0.0" top hits should still pass; the 1.2e-67 etc.
    # would have been excluded anyway by the top-hit reduction.
    assert len(strict) == 3, (
        f"strict cutoff retained {len(strict)} queries; expected 3"
    )

    # An impossibly-strict cutoff drops everything.
    impossible = filter_top_hit_per_query(tab_df, evalue_cutoff=1e-1000)
    # Behavior: float('1e-1000') underflows to 0.0 in Python. So
    # nothing should be < 0. Be explicit about this:
    assert len(impossible) == 0, (
        "impossibly-strict cutoff should retain zero queries"
    )

    print()
    print("[exercise-03] All assertions passed. You can now read both")
    print("[exercise-03] tabular and XML BLAST output, reduce to top")
    print("[exercise-03] hits, and join with taxonomy. Time to start")
    print("[exercise-03] the challenge and mini-project.")
