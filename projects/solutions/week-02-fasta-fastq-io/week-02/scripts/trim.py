#!/usr/bin/env python3
"""Sliding-window quality trimmer with layered configuration and a receipt.

This is homework Problem 4 (`p4_trim.py`) adapted for the mini-project.

Policy, applied in this order to every read:
  1. Sliding window of `--window` bases from the 5' end. Truncate at the first
     window whose mean Phred is below `--window-threshold` (Trimmomatic's
     SLIDINGWINDOW:4:20 semantics, Bolger et al. 2014).
  2. Drop the read if the surviving length is below `--min-len`.
  3. Drop the read if the surviving mean Phred is below `--min-mean-q`.

Configuration precedence, lowest to highest:
    built-in defaults  <  --config JSON  <  explicit command-line flags

Usage:
    python scripts/trim.py --in reads.fastq.gz --out trimmed.fastq.gz
    python scripts/trim.py --in reads.fastq.gz --out trimmed.fastq.gz \
        --config week-02-config.json --min-len 40

C10 Crunching Bioinformatics, Week 2.
Pinned: Python 3.11+, Biopython 1.83.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import IO

from Bio import SeqIO

DEFAULTS = {
    "min_len": 36,
    "min_mean_q": 20,
    "window": 4,
    "window_threshold": 20,
}

# Inner width of the C10 reproducibility receipt, in characters.
BOX_WIDTH = 57


def open_read(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def open_write(path: Path) -> IO[str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", newline="\n")
    return open(path, "w", newline="\n")


def trim_position(quals: list[int], window: int, threshold: int) -> int:
    """Return the index at which to truncate the read.

    Slides a `window`-wide window from the 5' end and returns the start index of
    the first window whose mean Phred falls below `threshold`. If no window
    fails, the read is kept whole. Reads shorter than one window are kept whole
    and left to the length filter.
    """
    if len(quals) < window:
        return len(quals)
    running = sum(quals[:window])
    limit = threshold * window          # integer compare; avoids a float divide
    if running < limit:
        return 0
    for i in range(1, len(quals) - window + 1):
        running += quals[i + window - 1] - quals[i - 1]
        if running < limit:
            return i
    return len(quals)


def force_utf8_console() -> None:
    """Make the receipt printable everywhere.

    The receipt is drawn with U+250C-style box characters. A Windows console
    still defaults to cp1252, where those characters do not exist, and the
    script would die with `UnicodeEncodeError: 'charmap' codec can't encode
    characters in position 0-58` *after* having already written the output
    FASTQ. Re-encoding the streams as UTF-8 costs nothing on POSIX and removes
    the only platform-specific failure in this pipeline.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8")
        except (AttributeError, OSError):   # already-redirected or exotic stream
            pass


def receipt(lines: list[tuple[str, str]]) -> str:
    """Render the C10 reproducibility receipt at a fixed width."""
    out = ["┌" + "─" * BOX_WIDTH + "┐"]
    out.append("│" + "  REPRODUCIBILITY".ljust(BOX_WIDTH) + "│")
    out.append("│" + " " * BOX_WIDTH + "│")
    for label, value in lines:
        body = f"  {label:<13}{value}"
        if len(body) > BOX_WIDTH:       # never let a long path break the box
            body = body[: BOX_WIDTH - 1] + "…"
        out.append("│" + body.ljust(BOX_WIDTH) + "│")
    out.append("└" + "─" * BOX_WIDTH + "┘")
    return "\n".join(out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Quality-trim and filter a FASTQ.")
    # `--in` cannot become `args.in` (a Python keyword), so name the dest.
    parser.add_argument("--in", dest="in_path", type=Path, required=True)
    parser.add_argument("--out", dest="out_path", type=Path, required=True)
    parser.add_argument("--min-len", type=int)
    parser.add_argument("--min-mean-q", type=int)
    parser.add_argument("--window", type=int)
    parser.add_argument("--window-threshold", type=int)
    parser.add_argument("--config", type=Path, help="JSON file of the same keys")
    return parser


def resolve_settings(argv: list[str] | None) -> argparse.Namespace:
    """Parse twice so JSON can sit between the defaults and the flags."""
    parser = build_parser()
    args = parser.parse_args(argv)

    settings = dict(DEFAULTS)
    if args.config is not None:
        if not args.config.exists():
            print(f"error: no such config: {args.config}", file=sys.stderr)
            raise SystemExit(1)
        with open(args.config) as handle:
            loaded = json.load(handle)
        unknown = set(loaded) - set(DEFAULTS)
        if unknown:
            print(f"error: unknown config keys: {sorted(unknown)}", file=sys.stderr)
            raise SystemExit(1)
        settings.update(loaded)

    # Anything the user typed on the command line wins over the JSON.
    for key in DEFAULTS:
        supplied = getattr(args, key)
        if supplied is not None:
            settings[key] = supplied
        else:
            setattr(args, key, settings[key])
    return args


def main(argv: list[str] | None = None) -> int:
    force_utf8_console()
    args = resolve_settings(argv)

    if not args.in_path.exists():
        print(f"error: no such file: {args.in_path}", file=sys.stderr)
        return 1

    n_in = n_out = 0
    dropped_short = dropped_lowq = 0
    bp_in = bp_out = 0

    with open_read(args.in_path) as reader, open_write(args.out_path) as writer:
        for record in SeqIO.parse(reader, "fastq"):
            n_in += 1
            quals = record.letter_annotations["phred_quality"]
            bp_in += len(quals)

            cut = trim_position(quals, args.window, args.window_threshold)
            # Slicing the SeqRecord keeps seq and phred_quality in lockstep.
            # Slicing record.seq alone would leave the qualities stale and
            # Biopython would refuse to write the record.
            trimmed = record[:cut]
            kept_quals = trimmed.letter_annotations["phred_quality"]

            if len(kept_quals) < args.min_len:
                dropped_short += 1
                continue
            if sum(kept_quals) / len(kept_quals) < args.min_mean_q:
                dropped_lowq += 1
                continue

            SeqIO.write(trimmed, writer, "fastq")
            n_out += 1
            bp_out += len(kept_quals)

    settings_line = (
        f"window={args.window} thr={args.window_threshold} "
        f"min_len={args.min_len} min_q={args.min_mean_q}"
    )
    print(
        receipt(
            [
                ("Input:", str(args.in_path)),
                ("Output:", str(args.out_path)),
                ("Records in:", str(n_in)),
                ("Records out:", str(n_out)),
                ("Settings:", settings_line),
            ]
        )
    )
    print(
        f"dropped: {dropped_short} too short, {dropped_lowq} below mean Q; "
        f"bases {bp_in} -> {bp_out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
