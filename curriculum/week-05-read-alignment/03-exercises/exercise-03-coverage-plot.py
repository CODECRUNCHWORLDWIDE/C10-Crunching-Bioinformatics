"""
Exercise 3 - Coverage plot from a sorted, indexed BAM.

Goal: open a sorted+indexed BAM file with pysam, compute per-position
read depth across the reference contig, summarize the distribution
(mean, median, stddev, coefficient of variation, fraction of zero-
coverage positions), and render a windowed coverage plot with
matplotlib. The output PNG is the visual deliverable for the week.

Estimated time: 50 minutes. Pure Python after the BAM is built.

Acceptance criteria:
- `python exercise-03-coverage-plot.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The output PNG at results/lambda_coverage.png exists and is non-empty.
- You implemented four functions: `per_position_depth`,
  `summarize_depth`, `window_depth`, and `plot_coverage`.

Requirements:
    conda install -c bioconda pysam=0.22 samtools=1.19
    pip install matplotlib numpy
    (and a sorted+indexed BAM from Exercise 1)

What you learn:
- How `pysam.AlignmentFile.pileup` walks a BAM column by column.
- The per-position-depth metric and what it looks like distributionally
  (Poisson-ish for uniform coverage, biased by GC content in practice).
- How to window depth values for a plot that fits a sane number of
  pixels (you cannot plot 4.6 million positions; you bin them).
- The mean / median / cv / zero-coverage summary that goes in every
  variant-calling methods section.

TO COMPLETE: implement the four functions below. Run the file; all
assertions must pass.

Tool versions assumed:
- pysam 0.22
- numpy 1.26.4
- matplotlib 3.8+
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import numpy as np


# Path to the BAM file produced by Exercise 1.
# Falls back to a small bundled example if needed.
BAM_PATH = Path(__file__).parent / "aln" / "lambda.sorted.bam"

# Expected reference contig name and length.
EXPECTED_CONTIG = "NC_001416.1"
EXPECTED_REF_LENGTH = 48502

# Window size in bp for the coverage plot. 200 bp gives ~240 windows
# across the 48.5 kb lambda reference - a reasonable plot resolution.
WINDOW_SIZE = 200

# Output PNG path.
OUTPUT_PNG = Path(__file__).parent / "results" / "lambda_coverage.png"


def per_position_depth(bam_path: Path, contig: str, length: int) -> np.ndarray:
    """Compute per-position read depth across the full contig.

    Open the BAM with pysam, walk every pileup column on `contig`, and
    populate an `np.int32` array of length `length` with the per-column
    `nsegments` count (number of aligned reads spanning that position).
    Positions with no reads remain at 0.

    Use `truncate=True` on the pileup call to restrict to the requested
    region; without it, pysam may include columns slightly outside
    [start, stop).

    Returns:
        An np.ndarray of shape (length,) and dtype int32.

    Hint: `pysam.AlignmentFile(bam_path, "rb").pileup(contig, 0, length,
    truncate=True)` yields PileupColumn objects with `.reference_pos`
    (0-based) and `.nsegments` (read count at that column).
    """
    import pysam

    depth = np.zeros(length, dtype=np.int32)
    # TODO: open the BAM in "rb" mode.
    # TODO: iterate over pileup columns; assign depth[column.reference_pos]
    #       = column.nsegments.
    raise NotImplementedError("Walk pileup columns and fill depth array")


def summarize_depth(depth: np.ndarray) -> dict:
    """Compute summary statistics over a per-position depth array.

    Returns a dict with keys:
        "mean":              float, np.mean(depth)
        "median":            float, np.median(depth)
        "stddev":            float, np.std(depth)
        "min":               int, np.min(depth)
        "max":               int, np.max(depth)
        "coefficient_of_variation": float, stddev / mean (NaN if mean=0)
        "n_zero_positions":  int, (depth == 0).sum()
        "fraction_zero":     float, n_zero / len(depth)
        "n_positions_ge_5":  int, (depth >= 5).sum()
        "n_positions_ge_20": int, (depth >= 20).sum()
        "length":            int, len(depth)

    A healthy alignment has:
        coefficient_of_variation < 0.5 (low variance relative to mean)
        fraction_zero very close to 0.0 (no big coverage gaps)
        n_positions_ge_5 == length (every position has at least 5x)
    """
    # TODO: compute and return the dict.
    raise NotImplementedError("Compute depth summary statistics")


def window_depth(depth: np.ndarray, window: int) -> np.ndarray:
    """Average per-position depth into non-overlapping windows.

    For a 48,502-position depth array and window=200, returns an array
    of shape (242,) with the mean depth per 200 bp window. The final
    partial window (positions 48400-48502) is included as its own
    window mean.

    Args:
        depth: per-position depth, shape (n,).
        window: window size in bp.

    Returns:
        np.ndarray of shape (ceil(n / window),) with the per-window
        mean depth.

    Hint: reshape complete windows with depth[:n_complete*window].reshape(
        n_complete, window).mean(axis=1), then append a final partial
    window if there is a remainder.
    """
    # TODO: compute n_complete = len(depth) // window.
    # TODO: average complete windows; if there is a leftover, average
    #       that too and concatenate.
    raise NotImplementedError("Window-average the depth array")


def plot_coverage(
    windowed: np.ndarray,
    window: int,
    contig: str,
    title: str,
    out_path: Path,
) -> None:
    """Render a per-window coverage plot to `out_path`.

    The plot has:
        - x axis: reference position in kb (integer ticks every ~10 kb)
        - y axis: mean coverage per window
        - a horizontal dashed line at the mean coverage
        - a horizontal dotted line at the median coverage
        - a one-line title
        - a legend

    Use matplotlib's default style. Save at 150 DPI for a sharp PNG.

    Args:
        windowed: array of per-window mean depths, shape (n_windows,).
        window: window size in bp (for the x-axis scaling).
        contig: contig name (for the y-axis label).
        title: plot title.
        out_path: where to save the PNG.

    Hint: import matplotlib.pyplot as plt. The x-axis positions in kb
    are np.arange(len(windowed)) * window / 1000.
    """
    import matplotlib

    # Use the Agg backend so we can run headless (no display).
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # TODO: build the figure, plot windowed depth vs x positions in kb.
    # TODO: add mean and median lines.
    # TODO: set labels, title, legend; save with dpi=150.
    raise NotImplementedError("Render the coverage plot")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-coverage-plot.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Toolchain sanity check.
    for tool in ("samtools",):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"[exercise-03] {tool!r} is not on PATH. Install it:\n"
                f"    conda install -c bioconda {tool}=...\n"
                f"and re-run."
            )

    # Confirm pysam imports.
    try:
        import pysam  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "[exercise-03] pysam is not installed. Install with:\n"
            "    conda install -c bioconda pysam=0.22\n"
            "or pip install pysam==0.22"
        ) from exc

    if not BAM_PATH.exists():
        raise SystemExit(
            f"[exercise-03] BAM not found at {BAM_PATH}.\n"
            f"Run Exercise 1 first to produce the sorted+indexed BAM."
        )

    bai_path = Path(str(BAM_PATH) + ".bai")
    if not bai_path.exists():
        raise SystemExit(
            f"[exercise-03] BAM index not found at {bai_path}.\n"
            f"Run `samtools index {BAM_PATH}` to produce it."
        )

    print("[exercise-03] Step 1: computing per-position depth ...")
    depth = per_position_depth(BAM_PATH, EXPECTED_CONTIG, EXPECTED_REF_LENGTH)
    assert depth.shape == (EXPECTED_REF_LENGTH,), (
        f"depth shape {depth.shape} != ({EXPECTED_REF_LENGTH},)"
    )
    assert depth.dtype == np.int32, (
        f"depth dtype {depth.dtype} != int32"
    )

    print(f"[exercise-03] Step 2: summarizing depth ...")
    stats = summarize_depth(depth)
    for key in (
        "mean", "median", "stddev", "min", "max",
        "coefficient_of_variation", "n_zero_positions", "fraction_zero",
        "n_positions_ge_5", "n_positions_ge_20", "length",
    ):
        assert key in stats, f"summary missing key {key!r}"

    # Sanity. With 1000 paired-end 150 bp reads against 48.5 kb lambda
    # at uniform coverage, expected mean coverage =
    #   (1000 pairs * 2 reads/pair * 150 bp/read) / 48502 bp = 6.2x.
    # Allow a wide range to accommodate wgsim's stochasticity.
    assert stats["mean"] > 1.0, (
        f"mean coverage {stats['mean']:.2f}x is implausibly low; "
        f"did BWA align anything?"
    )
    assert stats["length"] == EXPECTED_REF_LENGTH
    assert stats["min"] >= 0
    assert stats["max"] >= stats["mean"]

    print(f"[exercise-03] Step 3: windowing depth into {WINDOW_SIZE}-bp bins ...")
    windowed = window_depth(depth, WINDOW_SIZE)
    expected_windows = (EXPECTED_REF_LENGTH + WINDOW_SIZE - 1) // WINDOW_SIZE
    assert windowed.shape[0] == expected_windows, (
        f"got {windowed.shape[0]} windows; expected {expected_windows}"
    )

    # Window means should be in the same neighborhood as the raw mean.
    assert abs(windowed.mean() - stats["mean"]) < stats["mean"] * 0.1, (
        f"windowed mean {windowed.mean():.2f} differs from "
        f"per-position mean {stats['mean']:.2f} by > 10%"
    )

    print(f"[exercise-03] Step 4: rendering coverage plot ...")
    plot_coverage(
        windowed,
        WINDOW_SIZE,
        EXPECTED_CONTIG,
        title=(
            f"Bacteriophage lambda ({EXPECTED_CONTIG}) coverage from "
            f"BWA-MEM alignment\n"
            f"mean = {stats['mean']:.1f}x  "
            f"median = {stats['median']:.0f}x  "
            f"CV = {stats['coefficient_of_variation']:.2f}"
        ),
        out_path=OUTPUT_PNG,
    )
    assert OUTPUT_PNG.exists(), f"coverage plot not written to {OUTPUT_PNG}"
    assert OUTPUT_PNG.stat().st_size > 1000, (
        f"coverage plot at {OUTPUT_PNG} is suspiciously small "
        f"({OUTPUT_PNG.stat().st_size} bytes)"
    )

    print()
    print("[exercise-03] Coverage summary:")
    print(f"  Contig:               {EXPECTED_CONTIG} ({EXPECTED_REF_LENGTH} bp)")
    print(f"  Mean depth:           {stats['mean']:.2f}x")
    print(f"  Median depth:         {stats['median']:.1f}x")
    print(f"  Stddev:               {stats['stddev']:.2f}")
    print(f"  Coefficient of var.:  {stats['coefficient_of_variation']:.3f}")
    print(f"  Min / Max:            {stats['min']} / {stats['max']}")
    print(f"  Bases >= 5x coverage: {stats['n_positions_ge_5']} "
          f"({100*stats['n_positions_ge_5']/stats['length']:.1f}%)")
    print(f"  Bases >= 20x coverage: {stats['n_positions_ge_20']} "
          f"({100*stats['n_positions_ge_20']/stats['length']:.1f}%)")
    print(f"  Zero-coverage bases:  {stats['n_zero_positions']} "
          f"({100*stats['fraction_zero']:.2f}%)")
    print()
    print(f"[exercise-03] Plot saved to {OUTPUT_PNG}")
    print(f"[exercise-03] All assertions passed.")
    print(f"[exercise-03] You now have the end-to-end FASTQ -> coverage")
    print(f"[exercise-03] plot pipeline running locally. Continue to")
    print(f"[exercise-03] the challenge or mini-project.")
