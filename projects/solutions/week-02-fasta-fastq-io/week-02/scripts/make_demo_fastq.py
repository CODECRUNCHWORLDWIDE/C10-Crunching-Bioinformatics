#!/usr/bin/env python3
"""Generate a deterministic, Illumina-shaped demo FASTQ for offline smoke tests.

The mini-project is graded on a *real* 1000 Genomes subset (ERR1019034). This
generator exists so the pipeline can be exercised end to end with no network,
no 3 GB download, and byte-identical output on every machine — which is what
lets the reference numbers in ../README.md be checked rather than believed.

The quality model is a caricature of a HiSeq 2000 100 bp run: high Q at the 5'
end, a decaying 3' tail, a small population of bad clusters, and Ns injected
where the model draws a very low Q. It is not real data and must never be
presented as such.

Usage:
    python scripts/make_demo_fastq.py data/demo.fastq.gz --reads 100000 --seed 2026

C10 Crunching Bioinformatics, Week 2 mini-project reference implementation.
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

READ_LEN = 100
BASES = "ACGT"


def phred_to_char(q: int) -> str:
    """Phred+33 encoding: one printable ASCII byte per base."""
    return chr(q + 33)


def simulate_read(rng: random.Random) -> tuple[str, str]:
    """Return (sequence, Phred+33 quality string) for one synthetic read."""
    # Per-read offset: most clusters are fine, ~3% are bad.
    read_offset = rng.gauss(0.0, 1.8)
    if rng.random() < 0.03:
        read_offset -= 14.0

    seq_chars: list[str] = []
    qual_chars: list[str] = []
    for i in range(READ_LEN):
        # Positional decay: flat-ish for the first third, steep in the tail.
        positional = 38.0 - 20.0 * (i / (READ_LEN - 1)) ** 2.5
        q = int(round(positional + read_offset + rng.gauss(0.0, 1.5)))
        q = max(2, min(41, q))
        # Real base callers emit N when they have no confidence at all.
        if q < 8 and rng.random() < 0.30:
            seq_chars.append("N")
            q = 2
        else:
            seq_chars.append(rng.choice(BASES))
        qual_chars.append(phred_to_char(q))
    return "".join(seq_chars), "".join(qual_chars)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path, help="output path (.fastq or .fastq.gz)")
    parser.add_argument("--reads", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(args.out).endswith(".gz") else open

    with opener(args.out, "wt", newline="\n") as handle:
        for i in range(args.reads):
            seq, qual = simulate_read(rng)
            # Illumina 1.8+ style header. The read id is the first token.
            handle.write(
                f"@DEMO:1:FLOWCELL:1:1101:{1000 + i}:{2000 + i} 1:N:0:ATCACG\n"
                f"{seq}\n+\n{qual}\n"
            )

    print(f"wrote {args.reads} reads to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
