"""
Exercise 1 - BLAST an unknown sequence against NCBI nt via Biopython.

Goal: take a single unknown 16S rRNA sequence, submit it to NCBI BLAST
using `Bio.Blast.NCBIWWW.qblast` against the curated `16S_ribosomal_RNA`
database, parse the top hit out of the returned XML, and identify the
organism. This is the one-query happy path that the mini-project will
generalize to ~20 queries.

Estimated time: 50 minutes (mostly waiting on NCBI).

Acceptance criteria:
- `python exercise-01-blast-an-unknown-sequence.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The XML response is cached to `cache/exercise-01.xml` on first run; on
  subsequent runs the script reads from cache and does NOT hit NCBI.
- You implemented three functions: `load_query` (reads the FASTA),
  `run_or_load_blast` (submits to NCBI or reads cache), and
  `extract_top_hit` (parses the XML and returns the top hit as a
  dict with keys "accession", "title", "evalue", "bitscore",
  "pident", "align_length").

Requirements:
    python -m pip install biopython==1.83
    (and a real network connection to NCBI on first run)

What you learn:
- The `Bio.Blast.NCBIWWW.qblast` call surface and its key arguments.
- How to cache an XML response to disk to avoid re-querying NCBI.
- How to navigate a `Bio.Blast.NCBIXML.BlastRecord` and pull out the
  fields you care about.
- Why setting `Bio.Entrez.email` is non-negotiable for any script that
  touches the NCBI servers.

TO COMPLETE: implement the three functions below. Run the file; all
assertions must pass. The first run takes 30-120 seconds (network);
subsequent runs are instant (cache).

Reference accession: the unknown sequence used for this exercise is
the published 16S rRNA gene of an environmental isolate that BLAST
will identify unambiguously to the genus *Staphylococcus*. The closest
reference in the NCBI 16S database is NR_117740.1 (S. aureus type
strain).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


# NCBI requires a contact email on every E-utilities and BLAST call.
# REPLACE THIS with your real address before running.
Entrez.email = "you@example.com"


# A ~1400 bp 16S rRNA sequence from an unknown environmental isolate.
# The published top BLAST hit against `16S_ribosomal_RNA` is
# Staphylococcus aureus, NR_117740.1, > 99% identity.
UNKNOWN_FASTA = """>unknown_16S_001 16S rRNA, environmental isolate
TGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAACAGATAAGGAGCTTGCTCCTTTGACGTTAG
CGGCGGACGGGTGAGTAACACGTGGATAACCTACCTATAAGACTGGGATAACTCCGGGAAACCGGAGCTA
ATACCGGATAATATTTCGAACCGCATGGTTCGATAGTGAAAGATGGCTTTGCTATCACTTATAGATGGAC
CCGCGCCGTATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAACGATACGTAGCCGACCTGAGAGGG
TGATCGGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCA
ATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATT
AGGGAAGAACAAACGTGTAAGTAACTGTGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTA
CGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTA
GGCGGTTTCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGAAACTT
GAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGT
GGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGAT
ACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGC
TAACGCATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCC
GCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAAATCTTGACATCCTTT
GACAACTCTAGAGATAGAGCTTTCCCTTCGGGGACAAAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCG
TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTAAGCTTAGTTGCCATCATTCAGTTGG
GCACTCTAAGTTGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTA
TGACCTGGGCTACACACGTGCTACAATGGACAATACAAAGGGCAGCGAAACCGCGAGGTCAAGCAAATCC
CATAAAGTTGTTCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGAAGCTGGAATCGCTAGTAATCG
TAGATCAGCATGCTACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTG
TAACACCCGAAGCCGGTGGAGTAACCTTTTAGGAGCTAGCCGTCGAAGGTGGGACAAATGATTGGGGTGA
AGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAA
"""


def load_query(fasta_text: str = UNKNOWN_FASTA) -> SeqIO.SeqRecord:
    """Parse the FASTA text into a Biopython SeqRecord.

    Hint: `Bio.SeqIO.read` reads a single record from a file-like
    object. `io.StringIO(fasta_text)` gives you a file-like view on
    the multi-line string above.
    """
    # TODO: parse `fasta_text` and return the single SeqRecord.
    raise NotImplementedError("Parse the FASTA text into a SeqRecord")


def run_or_load_blast(
    query_seq: str,
    cache_path: Path,
    *,
    database: str = "16S_ribosomal_RNA",
    expect: float = 1e-50,
    hitlist_size: int = 10,
) -> Any:
    """Submit query to NCBI BLAST or read from the on-disk cache.

    If `cache_path` exists, open it and return the parsed BlastRecord.
    Otherwise, submit to NCBI via `NCBIWWW.qblast`, write the XML to
    `cache_path`, then read it back and return the parsed record.

    Returns:
        A `Bio.Blast.NCBIXML.BlastRecord` (one record per single query).

    Hint: `NCBIWWW.qblast(...)` returns a handle whose `.read()` gives
    you the XML as a string. Write it to disk; then re-open the file
    and pass it to `NCBIXML.read`.
    """
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists():
        # Cache hit: read XML from disk.
        # TODO: open cache_path, parse with NCBIXML.read, return record.
        raise NotImplementedError("Read the cached XML from disk")

    # Cache miss: query NCBI and save.
    print(
        f"[exercise-01] No cache at {cache_path}; querying NCBI BLAST "
        f"against {database!r}. This takes 30-120 seconds."
    )
    # TODO: call NCBIWWW.qblast with the arguments above.
    # `program="blastn"` for nucleotide queries.
    # Write the response XML to `cache_path`.
    # Re-open and parse with NCBIXML.read; return the record.
    raise NotImplementedError("Submit the qblast call and cache the XML")


def extract_top_hit(blast_record: Any) -> dict:
    """Pull the top hit out of a BlastRecord as a flat dict.

    The top hit is the alignment with the lowest E-value (the first
    one in `blast_record.alignments` when the record has been parsed
    from XML by `NCBIXML`, which sorts by E-value ascending).

    Returns a dict with keys:
        "accession":    e.g. "NR_117740.1"
        "title":        full subject title (organism + reference info)
        "evalue":       float (the HSP's E-value)
        "bitscore":     float (the HSP's bit score)
        "pident":       float (percent identity, 0-100)
        "align_length": int (number of aligned columns, including gaps)
        "identities":   int (number of identical positions)

    Raise `ValueError` if the record has no alignments.

    Hint: `pident = identities / align_length * 100`. Use the FIRST
    HSP of the top alignment (BLAST sorts HSPs within an alignment
    by score descending).
    """
    # TODO: check that blast_record.alignments is non-empty; raise
    # ValueError if not.
    # Take the first alignment; take its first HSP. Build and return
    # the dict.
    raise NotImplementedError("Extract the top hit fields")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-01-blast-an-unknown-sequence.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    cache_dir = Path(__file__).parent / "cache"
    cache_path = cache_dir / "exercise-01.xml"

    print("[exercise-01] Loading query sequence ...")
    query = load_query()
    assert query is not None, "load_query returned None"
    assert len(query.seq) > 1300, (
        f"unexpected query length {len(query.seq)}; expected ~1400 bp"
    )
    print(f"[exercise-01] Query: {query.id}, length: {len(query.seq)} bp")

    print("[exercise-01] Running or loading BLAST ...")
    blast_record = run_or_load_blast(str(query.seq), cache_path)
    assert blast_record is not None, "run_or_load_blast returned None"
    assert len(blast_record.alignments) > 0, (
        "no alignments in BLAST record; check query and database"
    )

    print("[exercise-01] Extracting top hit ...")
    top = extract_top_hit(blast_record)

    # Field presence.
    for key in (
        "accession", "title", "evalue", "bitscore",
        "pident", "align_length", "identities",
    ):
        assert key in top, f"top hit missing field {key!r}"

    # Sanity. Against 16S_ribosomal_RNA, this query should be:
    # - Very high identity (> 95%).
    # - Highly significant (E < 1e-100).
    # - Aligned over essentially the full ~1400 bp query.
    assert top["pident"] > 95.0, (
        f"top hit pident = {top['pident']}; expected > 95"
    )
    assert top["evalue"] < 1e-100, (
        f"top hit E-value = {top['evalue']}; expected < 1e-100"
    )
    assert top["align_length"] > 1300, (
        f"top hit align_length = {top['align_length']}; expected > 1300"
    )

    # The genus should be Staphylococcus per the curriculum's expected
    # answer. Check that the title mentions it.
    assert "Staphylococcus" in top["title"], (
        f"top hit title did not contain 'Staphylococcus':\n  {top['title']!r}"
    )

    print()
    print("[exercise-01] Top hit:")
    print(f"  Accession:    {top['accession']}")
    print(f"  Title:        {top['title'][:80]}")
    print(f"  E-value:      {top['evalue']:.2e}")
    print(f"  Bit score:    {top['bitscore']:.1f}")
    print(f"  Identity:     {top['identities']}/{top['align_length']}  "
          f"({top['pident']:.2f}%)")

    print()
    print("[exercise-01] All assertions passed. The unknown isolate is")
    print("[exercise-01] confidently identified as genus Staphylococcus.")
    print("[exercise-01] Continue to exercise-02 (local database).")
