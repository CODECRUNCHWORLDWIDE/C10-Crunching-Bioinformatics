#!/usr/bin/env python3
"""Per-base quality profile plot for a FASTQ file, computed in one streaming pass.

This is Exercise 2 (`plot_quality.py`) adapted for the mini-project.

Memory is bounded: the accumulator is a per-position histogram over the 0..93
Phred range, so it grows with the *read length* (a few hundred), never with the
read count. That is what lets the same script run on a 100k-read subset and on
a 400-million-read whole-genome run.

Usage:
    python scripts/plot_quality.py reads.fastq.gz out.png [--label "raw"]

C10 Crunching Bioinformatics, Week 2.
Pinned: Python 3.11+, Biopython 1.83, matplotlib.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import IO

import matplotlib

matplotlib.use("Agg")            # headless: no display needed on CI or a server
import matplotlib.pyplot as plt  # noqa: E402  (must follow the backend choice)

from Bio import SeqIO            # noqa: E402

MAX_Q = 93                       # Phred+33 tops out at 93 ('~')


def open_maybe_gz(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def quality_histograms(path: Path) -> tuple[list[list[int]], int]:
    """Return (per-position Q histograms, n_reads). One pass, bounded memory."""
    histograms: list[list[int]] = []
    n_reads = 0
    with open_maybe_gz(path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            n_reads += 1
            quals = record.letter_annotations["phred_quality"]
            while len(histograms) < len(quals):
                histograms.append([0] * (MAX_Q + 1))
            for position, q in enumerate(quals):
                histograms[position][q] += 1
    return histograms, n_reads


def summarise(histograms: list[list[int]]) -> tuple[list[float], list[int], list[int]]:
    """Mean, lower quartile and upper quartile of Q at each read position."""
    means: list[float] = []
    q25: list[int] = []
    q75: list[int] = []
    for hist in histograms:
        total = sum(hist)
        if total == 0:
            means.append(0.0)
            q25.append(0)
            q75.append(0)
            continue
        means.append(sum(q * c for q, c in enumerate(hist)) / total)
        q25.append(_quantile_from_hist(hist, total, 0.25))
        q75.append(_quantile_from_hist(hist, total, 0.75))
    return means, q25, q75


def _quantile_from_hist(hist: list[int], total: int, fraction: float) -> int:
    target = fraction * total
    running = 0
    for q, count in enumerate(hist):
        running += count
        if running >= target:
            return q
    return len(hist) - 1


def plot(path: Path, out_png: Path, label: str) -> int:
    histograms, n_reads = quality_histograms(path)
    if n_reads == 0:
        print(f"error: no reads found in {path}", file=sys.stderr)
        return 1

    means, q25, q75 = summarise(histograms)
    positions = list(range(1, len(means) + 1))    # 1-based, like FastQC

    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=150)
    ax.fill_between(positions, q25, q75, alpha=0.25, label="interquartile range")
    ax.plot(positions, means, linewidth=2, label="mean Phred Q")
    ax.axhline(20, linestyle="--", linewidth=1, color="crimson", label="Q20 (1% error)")

    ax.set_xlabel("Position in read (bp, 1-based)")
    ax.set_ylabel("Phred quality score")
    ax.set_title(f"Per-base quality — {label} — {path.name} (n = {n_reads:,} reads)")
    ax.set_ylim(0, 42)
    ax.set_xlim(1, len(means))
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left")
    fig.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)

    print(
        f"{out_png}: {n_reads} reads, {len(means)} positions, "
        f"mean Q pos1={means[0]:.2f} posN={means[-1]:.2f}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Plot per-base quality of a FASTQ.")
    parser.add_argument("fastq", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--label", default="raw", help="text shown in the title")
    args = parser.parse_args(argv)

    if not args.fastq.exists():
        print(f"error: no such file: {args.fastq}", file=sys.stderr)
        return 1
    return plot(args.fastq, args.out, args.label)


if __name__ == "__main__":
    raise SystemExit(main())
