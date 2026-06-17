"""
Exercise 2 — FASTA by hand (no Biopython).

Goal: Write a working FASTA parser using only the Python standard library.
No Biopython, no regex tricks, no third-party packages.

Estimated time: 35 minutes.

Why we do this: in Week 2 we will use `Bio.SeqIO.parse` and you will be
thankful for the abstraction. But if your first experience of FASTA is
through a library, you will never debug a malformed file. FASTA is one
of the simplest formats in computational biology. Parsing it by hand
takes 30 lines. Do it once. You will read other people's parsers more
confidently for the rest of your career.

Acceptance criteria:
- `python exercise-02-fasta-by-hand.py` runs without crashing.
- The two functions below pass all `assert` checks at the bottom.
- Your implementation does NOT import Biopython, regex, or any non-stdlib package.
- Your implementation handles:
    * Multiple records per file.
    * Sequences spanning multiple lines.
    * Optional trailing whitespace / blank lines.
    * An empty file (returns an empty iterator).

What you learn:
- The shape of a FASTA record (header `>...` then sequence lines).
- The difference between iterating and accumulating in a parser.
- Why writing your own parser before reading the Biopython source makes
  reading the Biopython source much easier.

TO COMPLETE: implement `parse_fasta` and `fasta_stats`. Do not look at
the hint until you have tried at least 15 minutes.
"""

from __future__ import annotations

from typing import Iterator


# A small inline test FASTA. In Week 2 we will read real files from disk.
TEST_FASTA = """>seq1 a tiny test record
ATGCGT
ACGT
>seq2 another record, lowercase mixed in
atgc
GCGC
>seq3 single-line record
AAAAAA
"""


def parse_fasta(text: str) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs from a FASTA-formatted string.

    The header is everything after the leading '>' on the header line,
    stripped of trailing whitespace. The sequence is the concatenation of
    all subsequent non-header lines until the next header (or EOF),
    with whitespace removed and uppercase normalisation applied.

    A FASTA file with no records yields nothing.

    Examples:
        >>> list(parse_fasta(">a\\nACGT\\n>b\\nGGGG"))
        [('a', 'ACGT'), ('b', 'GGGG')]
    """
    # TODO: implement this function.
    #
    # Suggested approach:
    #   1. Split the input on newlines.
    #   2. Walk the lines. When you see a line starting with '>', that
    #      header begins a new record. Save the in-progress record (if
    #      any), then start a new one.
    #   3. Otherwise, append the line's contents (stripped, upper-cased)
    #      to the in-progress record's sequence.
    #   4. After the loop, emit the final in-progress record.
    #
    # Edge cases to handle:
    #   - The file may end without a trailing newline.
    #   - Blank lines between records should be ignored.
    #   - Lines may have trailing whitespace; strip it.
    raise NotImplementedError("Implement parse_fasta")


def fasta_stats(records: list[tuple[str, str]]) -> dict[str, float | int]:
    """Compute summary stats for a list of (header, sequence) records.

    Returns a dict with keys:
        - 'n_records'  : int, number of records
        - 'total_bp'   : int, sum of all sequence lengths
        - 'mean_len'   : float, mean sequence length (0.0 if no records)
        - 'gc_percent' : float, percentage of G+C across all sequences
                          (counted over A, C, G, T only — N and other
                          ambiguity codes are excluded from both numerator
                          and denominator), 0.0 if no A/C/G/T bases

    Examples:
        >>> fasta_stats([('a', 'ACGT'), ('b', 'GGGG')])
        {'n_records': 2, 'total_bp': 8, 'mean_len': 4.0, 'gc_percent': 75.0}
    """
    # TODO: implement this function. Do not call Biopython.
    #
    # gc_percent uses only A/C/G/T as the denominator. N's are common
    # in real reference genomes (the human reference has millions of N's
    # in unsequenced regions); including them would deflate the GC%.
    raise NotImplementedError("Implement fasta_stats")


# ----------------------------------------------------------------------
# Self-test. Run this file directly with `python exercise-02-fasta-by-hand.py`.
# ----------------------------------------------------------------------
if __name__ == "__main__":
    records = list(parse_fasta(TEST_FASTA))

    # Basic structure checks.
    assert len(records) == 3, f"expected 3 records, got {len(records)}"
    headers = [h for h, _ in records]
    assert headers[0] == "seq1 a tiny test record"
    assert headers[1] == "seq2 another record, lowercase mixed in"
    assert headers[2] == "seq3 single-line record"

    # Sequence content checks. Multi-line records must be joined.
    seqs = [s for _, s in records]
    assert seqs[0] == "ATGCGTACGT", f"seq1 wrong: {seqs[0]!r}"
    assert seqs[1] == "ATGCGCGC", f"seq2 wrong (case-normalize): {seqs[1]!r}"
    assert seqs[2] == "AAAAAA", f"seq3 wrong: {seqs[2]!r}"

    # Edge cases.
    assert list(parse_fasta("")) == [], "empty input must yield no records"
    assert list(parse_fasta(">only_header\n")) == [("only_header", "")], (
        "a header with no sequence is still a record (empty sequence)"
    )

    # fasta_stats checks.
    stats = fasta_stats(records)
    assert stats["n_records"] == 3
    assert stats["total_bp"] == 24, f"total_bp wrong: {stats['total_bp']}"
    assert abs(stats["mean_len"] - 8.0) < 1e-9, f"mean_len wrong: {stats['mean_len']}"
    # GC content of "ATGCGTACGTATGCGCGCAAAAAA":
    #   counts: A=9, T=3, G=6, C=6. GC = 12 / 24 = 50.0
    assert abs(stats["gc_percent"] - 50.0) < 1e-9, f"gc% wrong: {stats['gc_percent']}"

    # Empty-input stats.
    empty_stats = fasta_stats([])
    assert empty_stats == {"n_records": 0, "total_bp": 0, "mean_len": 0.0, "gc_percent": 0.0}

    print("All assertions passed. Move on to Exercise 3.")
