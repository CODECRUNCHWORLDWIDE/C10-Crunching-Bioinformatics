#!/usr/bin/env python3
"""Emit the mini-project's machine-readable summary.json.

Produces exactly the shape the mini-project specifies, with tool versions read
from the environment rather than typed in by hand — a version string you had to
retype is a version string that will eventually be wrong.

Usage:
    python scripts/summarize.py --raw data/sub.fastq.gz \
        --trimmed qc/trimmed.fastq.gz --out qc/summary.json \
        --sample NA12878 --accession ERR1019034 --chromosome 22

C10 Crunching Bioinformatics, Week 2 mini-project reference implementation.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import IO

import Bio
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from inspect_fastq import detect_encoding  # noqa: E402

# The read positions the mini-project asks for, 1-based as FastQC reports them.
REPORT_POSITIONS = (1, 35, 75)


def open_maybe_gz(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def scan(path: Path) -> tuple[int, float, dict[int, float]]:
    """Return (n_reads, mean_len, {1-based position: mean Phred}) in one pass."""
    n_reads = 0
    total_len = 0
    sums: list[int] = []
    counts: list[int] = []

    with open_maybe_gz(path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            quals = record.letter_annotations["phred_quality"]
            n_reads += 1
            total_len += len(quals)
            while len(sums) < len(quals):
                sums.append(0)
                counts.append(0)
            for position, q in enumerate(quals):
                sums[position] += q
                counts[position] += 1

    per_position = {}
    for wanted in REPORT_POSITIONS:
        index = wanted - 1
        if index < len(counts) and counts[index]:
            per_position[wanted] = round(sums[index] / counts[index], 2)
        else:
            # A read length shorter than the requested position is a fact about
            # the data, not an error. Say so explicitly rather than emitting 0.
            per_position[wanted] = None

    mean_len = round(total_len / n_reads, 2) if n_reads else 0.0
    return n_reads, mean_len, per_position


def tool_version(executable: str, args: list[str]) -> str | None:
    """Return the first line of `executable args`, or None if it is not installed."""
    if shutil.which(executable) is None:
        return None
    try:
        out = subprocess.run(
            [executable, *args], capture_output=True, text=True, timeout=60, check=False
        )
    except OSError:
        return None
    text = (out.stdout or out.stderr).strip().splitlines()
    return text[0].strip() if text else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write the mini-project summary.json.")
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--trimmed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sample", default="NA12878")
    parser.add_argument("--accession", default="ERR1019034")
    parser.add_argument("--chromosome", default="22")
    args = parser.parse_args(argv)

    for path in (args.raw, args.trimmed):
        if not path.exists():
            print(f"error: no such file: {path}", file=sys.stderr)
            return 1

    n_raw, len_raw, q_raw = scan(args.raw)
    n_trim, len_trim, _ = scan(args.trimmed)

    summary = {
        "sample": args.sample,
        "accession": args.accession,
        "chromosome": args.chromosome,
        "n_reads_raw": n_raw,
        "n_reads_trimmed": n_trim,
        "mean_len_raw": len_raw,
        "mean_len_trimmed": len_trim,
        "mean_per_base_q_raw_pos_1": q_raw[1],
        "mean_per_base_q_raw_pos_35": q_raw[35],
        "mean_per_base_q_raw_pos_75": q_raw[75],
        "encoding_detected": detect_encoding(args.raw),
        "tools": {
            "fastqc": tool_version("fastqc", ["--version"]),
            "biopython": Bio.__version__,
            "seqkit": tool_version("seqkit", ["version"]),
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="\n") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
