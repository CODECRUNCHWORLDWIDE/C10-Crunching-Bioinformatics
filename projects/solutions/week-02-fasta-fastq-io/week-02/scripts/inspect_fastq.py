#!/usr/bin/env python3
"""Inspect a FASTQ file: record count, mean length, mean per-read mean Q, encoding.

This is homework Problem 1 (`p1_inspect_fastq.py`) adapted for the mini-project:
the numbers are also returned as a dict so `summarize.py` can reuse them instead
of re-implementing the same pass.

Usage:
    python scripts/inspect_fastq.py path/to/reads.fastq.gz

C10 Crunching Bioinformatics, Week 2.
Pinned: Python 3.11+, Biopython 1.83.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import re
import sys
from pathlib import Path
from typing import IO

from Bio import SeqIO

# Format name -> the Biopython SeqIO format string that decodes it correctly.
BIOPYTHON_FORMAT = {
    "Phred+33": "fastq",            # Sanger / Illumina 1.8+
    "Phred+64": "fastq-illumina",   # Illumina 1.3+ / 1.5+
}

# ENA/SRA run accessions: (E|S|D)RR followed by digits.
ACCESSION_RE = re.compile(r"([EDS]RR\d+)")


def open_maybe_gz(path: Path) -> IO[str]:
    """Return a text-mode handle for a plain or gzipped FASTQ."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def accession_from_name(path: Path) -> str:
    """Best-effort accession from the file name; falls back to the stem."""
    match = ACCESSION_RE.search(path.name)
    if match:
        return match.group(1)
    return path.name.split(".")[0]


def iter_quality_lines(path: Path, limit: int):
    """Yield up to `limit` raw quality strings, counting lines in fours.

    Deliberately raw: encoding detection must look at the *bytes on disk*, and
    `SeqIO.parse` has already decoded them (and would have picked an offset for
    us, which is the very thing we are trying to determine).
    """
    with open_maybe_gz(path) as handle:
        for index, line in enumerate(handle):
            if index % 4 == 3:            # line 4 of each record
                yield line.rstrip("\n")
                limit -= 1
                if limit <= 0:
                    return


def detect_encoding(path: Path, sample: int = 1000, strict: bool = False) -> str:
    """Return "Phred+33", "Phred+64", or "unknown" from the minimum ASCII byte.

    Boundaries are the ones stated in the Week 2 homework:
      min ASCII < 59  -> Phred+33 (a Phred+64 file cannot go below `@` = 64)
      min ASCII >= 64 -> Phred+64
      59..63          -> ambiguous. Default to Phred+33 with a stderr warning;
                         `strict=True` reports "unknown" instead.
    """
    quals = [q for q in itertools.islice(iter_quality_lines(path, sample), sample) if q]
    if not quals:
        return "unknown"
    min_ascii = min(ord(c) for q in quals for c in q)
    if min_ascii < 59:
        return "Phred+33"
    if min_ascii >= 64:
        return "Phred+64"
    print(
        f"warning: minimum ASCII {min_ascii} falls in the ambiguous 59-63 band; "
        f"sampled only {len(quals)} records",
        file=sys.stderr,
    )
    return "unknown" if strict else "Phred+33"


def inspect(path: Path) -> dict:
    """Stream the whole file once and return the summary numbers."""
    encoding = detect_encoding(path)
    fmt = BIOPYTHON_FORMAT.get(encoding, "fastq")

    n_records = 0
    total_len = 0
    sum_of_read_means = 0.0

    with open_maybe_gz(path) as handle:
        for record in SeqIO.parse(handle, fmt):
            quals = record.letter_annotations["phred_quality"]
            if not quals:
                continue
            n_records += 1
            total_len += len(quals)
            # Accumulate sums; one division per read beats statistics.fmean here
            # because fmean builds an intermediate float list per call.
            sum_of_read_means += sum(quals) / len(quals)

    return {
        "accession": accession_from_name(path),
        "n_records": n_records,
        "mean_len": round(total_len / n_records, 2) if n_records else 0.0,
        "mean_of_read_mean_q": round(sum_of_read_means / n_records, 2) if n_records else 0.0,
        "encoding_detected": encoding,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Summarise a FASTQ file.")
    parser.add_argument("path", type=Path, help="FASTQ file (.gz supported)")
    args = parser.parse_args(argv)

    if not args.path.exists():
        print(f"error: no such file: {args.path}", file=sys.stderr)
        return 1

    summary = inspect(args.path)
    print(f"accession:            {summary['accession']}")
    print(f"records:              {summary['n_records']}")
    print(f"mean read length:     {summary['mean_len']} bp")
    print(f"mean per-read mean Q: {summary['mean_of_read_mean_q']}")
    print(f"encoding detected:    {summary['encoding_detected']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
