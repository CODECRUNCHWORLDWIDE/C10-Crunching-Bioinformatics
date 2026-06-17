"""
Exercise 2 - Build a local BLAST database and query it.

Goal: download a small reference FASTA (the published 16S rRNA gene of
*Escherichia coli* K-12 MG1655, NCBI accession NR_074549.1, ~1,541 bp)
plus a handful of related 16S sequences, build a local nucleotide BLAST
database with `makeblastdb`, query an unknown 16S sequence against it
with `blastn`, and parse the tabular output.

Estimated time: 45 minutes.

Acceptance criteria:
- `python exercise-02-build-local-db.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The script creates a `data/` directory next to itself containing the
  downloaded FASTA, the BLAST database files, and the query results.
- On a second run, the script SKIPS the download and the makeblastdb
  step (idempotent caching).
- You implemented four functions: `fetch_references`, `build_database`,
  `run_blastn`, and `parse_tabular`.

Requirements:
    python -m pip install biopython==1.83 pandas
    conda install -c bioconda blast=2.15
    (the `makeblastdb`, `blastn`, and `blastdbcmd` executables must be
     on your PATH)

What you learn:
- How to call `makeblastdb` from Python via `subprocess.run`.
- How to call `blastn` and read its tabular output.
- How to wire up an idempotent download/build/query pipeline that does
  not re-do work on a second run.
- How to parse `-outfmt 6` output into a pandas DataFrame with the
  twelve canonical columns.

TO COMPLETE: implement the four functions below.

Reference accessions used in this exercise (the small local database):
    NR_074549.1  Escherichia coli K-12 MG1655 16S rRNA
    NR_117741.1  Staphylococcus aureus subsp. aureus 16S rRNA
    NR_119213.1  Bacillus subtilis 168 16S rRNA
    NR_117150.1  Pseudomonas aeruginosa PAO1 16S rRNA
    NR_118889.1  Mycobacterium tuberculosis H37Rv 16S rRNA
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Sequence

import pandas as pd
from Bio import Entrez, SeqIO


# REPLACE THIS with your real address before running.
Entrez.email = "you@example.com"


# The small reference set: five 16S rRNA reference sequences, one per
# genus, chosen to give the unknown query an unambiguous best hit.
REFERENCE_ACCESSIONS: tuple[str, ...] = (
    "NR_074549.1",   # Escherichia coli K-12 MG1655
    "NR_117741.1",   # Staphylococcus aureus subsp. aureus
    "NR_119213.1",   # Bacillus subtilis 168
    "NR_117150.1",   # Pseudomonas aeruginosa PAO1
    "NR_118889.1",   # Mycobacterium tuberculosis H37Rv
)


# An unknown 16S sequence. This is the same Staphylococcus aureus-like
# environmental isolate sequence used in exercise-01. Against the small
# local database below, the top hit must be NR_117741.1 (S. aureus).
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
"""


TABULAR_COLUMNS: tuple[str, ...] = (
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
)


def fetch_references(
    accessions: Sequence[str],
    out_fasta: Path,
) -> Path:
    """Download the listed accessions from NCBI and concatenate into one FASTA.

    Idempotent: if `out_fasta` already exists, skip the download and
    return the existing path.

    Hint: `Bio.Entrez.efetch(db="nuccore", id=",".join(accessions),
    rettype="fasta", retmode="text")` returns a handle that yields a
    multi-FASTA. Write its contents to `out_fasta` and return the path.
    """
    if out_fasta.exists():
        print(f"[exercise-02] References already at {out_fasta}; skipping fetch.")
        return out_fasta

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    print(
        f"[exercise-02] Fetching {len(accessions)} accessions from NCBI: "
        f"{', '.join(accessions)}"
    )
    # TODO: call Entrez.efetch, read the handle, write to out_fasta.
    raise NotImplementedError("Fetch references and write to out_fasta")


def build_database(in_fasta: Path, db_prefix: Path) -> Path:
    """Run `makeblastdb` to build a nucleotide BLAST database.

    Idempotent: if `<db_prefix>.nhr` already exists, skip the build.

    Returns the `db_prefix` path (for use as the `-db` argument to
    `blastn`).

    Hint: `subprocess.run(["makeblastdb", "-in", str(in_fasta), ...],
    check=True)`. The mandatory flags are `-in`, `-dbtype nucl`,
    `-out`, and `-parse_seqids`. Pass `capture_output=True` so the
    `makeblastdb` chatter does not pollute your script's stdout.
    """
    nhr_path = db_prefix.with_suffix(".nhr")
    if nhr_path.exists():
        print(f"[exercise-02] Database already at {db_prefix}; skipping build.")
        return db_prefix

    db_prefix.parent.mkdir(parents=True, exist_ok=True)

    print(f"[exercise-02] Building BLAST database at {db_prefix} ...")
    # TODO: subprocess.run(...) the makeblastdb command.
    # Check the return code (or use check=True).
    raise NotImplementedError("Run makeblastdb")


def run_blastn(
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    *,
    evalue: float = 1e-5,
    task: str = "blastn",
) -> Path:
    """Run `blastn` with `-outfmt 6` and the default 12 columns.

    Returns the path to the written results TSV.

    Hint: subprocess.run with `["blastn", "-query", ..., "-db", ...,
    "-out", ..., "-outfmt", "6", "-evalue", str(evalue), "-task", task]`.
    Use `check=True` so non-zero exit codes raise.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    print(f"[exercise-02] Running blastn (task={task}) on {query_fasta} ...")
    # TODO: build the argument list and run subprocess.run with check=True.
    raise NotImplementedError("Run blastn")


def parse_tabular(tsv_path: Path) -> pd.DataFrame:
    """Read a `-outfmt 6` BLAST output into a pandas DataFrame.

    Returns a DataFrame with the twelve canonical columns
    (TABULAR_COLUMNS).

    Hint: `pd.read_csv(tsv_path, sep="\\t", names=list(TABULAR_COLUMNS))`.
    Empty file (zero hits) should produce an empty DataFrame with the
    twelve columns, not crash.
    """
    # TODO: implement, including the empty-file case.
    raise NotImplementedError("Parse the tabular output")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-02-build-local-db.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    here = Path(__file__).parent
    data_dir = here / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = data_dir / "references_16S.fasta"
    db_prefix = data_dir / "ref_16S_db"
    query_fasta = data_dir / "unknown.fasta"
    results_tsv = data_dir / "results.tsv"

    # Stage the query FASTA on disk so blastn can read it.
    query_fasta.write_text(UNKNOWN_FASTA)

    print("[exercise-02] Step 1 - fetch reference FASTA from NCBI ...")
    fetch_references(REFERENCE_ACCESSIONS, reference_fasta)
    assert reference_fasta.exists(), "reference FASTA was not written"

    # Confirm we got the expected number of records.
    refs = list(SeqIO.parse(reference_fasta, "fasta"))
    assert len(refs) == len(REFERENCE_ACCESSIONS), (
        f"expected {len(REFERENCE_ACCESSIONS)} reference records, "
        f"got {len(refs)}"
    )

    print("[exercise-02] Step 2 - build the local BLAST database ...")
    build_database(reference_fasta, db_prefix)
    assert db_prefix.with_suffix(".nhr").exists(), "BLAST DB not built"

    print("[exercise-02] Step 3 - blastn the unknown query ...")
    run_blastn(query_fasta, db_prefix, results_tsv, evalue=1e-10, task="blastn")
    assert results_tsv.exists(), "blastn did not write a results file"

    print("[exercise-02] Step 4 - parse the tabular output ...")
    df = parse_tabular(results_tsv)

    # Shape and column checks.
    assert list(df.columns) == list(TABULAR_COLUMNS), (
        f"DataFrame columns wrong: {list(df.columns)}"
    )
    assert len(df) > 0, "expected at least one hit; got zero"

    # Top hit (lowest E-value, highest bit score) must be the
    # S. aureus reference NR_117741.1.
    top = df.sort_values("evalue").iloc[0]
    assert top["sseqid"].startswith("NR_117741"), (
        f"expected top hit NR_117741.1, got {top['sseqid']!r}"
    )
    # Identity should be very high (the query is essentially S. aureus).
    assert top["pident"] > 95.0, (
        f"top hit pident = {top['pident']}; expected > 95"
    )
    # E-value should be tiny.
    assert top["evalue"] < 1e-100, (
        f"top hit E-value = {top['evalue']}; expected < 1e-100"
    )

    print()
    print("[exercise-02] Top hit:")
    print(f"  Subject:    {top['sseqid']}")
    print(f"  pident:     {top['pident']:.2f}%")
    print(f"  length:     {top['length']}")
    print(f"  E-value:    {top['evalue']:.2e}")
    print(f"  bit score:  {top['bitscore']:.1f}")

    # Confirm idempotency: a second call to fetch_references should
    # NOT re-fetch (we cannot trivially observe this without a network
    # mock; instead, check that the file's mtime is unchanged on a
    # repeat call).
    print()
    print("[exercise-02] Step 5 - check idempotency on a repeat call ...")
    mtime_before = reference_fasta.stat().st_mtime
    fetch_references(REFERENCE_ACCESSIONS, reference_fasta)
    mtime_after = reference_fasta.stat().st_mtime
    assert mtime_before == mtime_after, (
        "fetch_references re-fetched on a second call; "
        "your idempotency check is missing"
    )

    print()
    print("[exercise-02] All assertions passed. Local database built,")
    print("[exercise-02] queried, and parsed. Continue to exercise-03.")
