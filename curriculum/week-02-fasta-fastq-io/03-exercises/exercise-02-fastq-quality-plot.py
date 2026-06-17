"""
Exercise 2 - FASTQ quality plot.

Goal: write code that reads a FASTQ file, computes the mean Phred quality
at each base position across all reads, and saves a PNG line plot of the
result. This is the single most-looked-at QC plot in short-read
bioinformatics. FastQC computes a fancier version (box plots per base);
the mean is the version you can compute in 20 lines of Python.

Estimated time: 45 minutes.

Acceptance criteria:
- `python exercise-02-fastq-quality-plot.py` runs without crashing.
- The functions below pass the assertions at the bottom.
- The script writes a `mean_quality.png` figure when given a FASTQ file
  with at least 1 record. We exercise that path with an inline synthetic
  FASTQ so the test does not depend on the user having downloaded data.

Required:
    python -m pip install biopython==1.83 matplotlib

What you learn:
- How to extract per-base Phred quality from a SeqRecord
  (`record.letter_annotations["phred_quality"]`).
- How to aggregate across reads of *variable* length without a crash on
  the index-out-of-range you'd get from a naive nested loop.
- How to convert Phred Q to error probability and back.
- Why matplotlib is "fine" for QC plots and overkill is unnecessary.

TO COMPLETE: implement the four functions below. Run the file; all
assertions must pass.
"""

from __future__ import annotations

import io
import math
import statistics
import tempfile
from pathlib import Path

import matplotlib

# Use the non-interactive backend so the script works on a headless box
# (a CI runner, a remote shell). Switch to a GUI backend yourself if you
# want to view interactively.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Bio import SeqIO  # noqa: E402


# A small inline FASTQ to test against without needing a download.
# Four reads of varying length and quality. Phred+33 encoding.
# Quality characters: '!' = Q0, '+' = Q10, '5' = Q20, '?' = Q30, 'I' = Q40.
TEST_FASTQ = """@read1 short, high quality throughout
ACGTACGT
+
IIIIIIII
@read2 medium length, decays at 3'
ACGTACGTACGT
+
IIIIIIIII???
@read3 longest, quality drops sharply
ACGTACGTACGTACGT
+
IIIIIIII?????+++
@read4 same length as read1, low quality
ACGTACGT
+
+++++!!!
"""


def phred_to_perror(q: int) -> float:
    """Convert a Phred quality value Q to an error probability P.

    P = 10 ** (-Q / 10)

    Examples:
        >>> phred_to_perror(10)
        0.1
        >>> round(phred_to_perror(30), 4)
        0.001
    """
    # TODO: implement (it is one line).
    raise NotImplementedError("Implement phred_to_perror")


def perror_to_phred(p: float) -> float:
    """Convert an error probability P to a Phred quality value Q.

    Q = -10 * log10(P)

    Raises ValueError if p <= 0 or p > 1 (the formula is not defined there).

    Examples:
        >>> round(perror_to_phred(0.01), 4)
        20.0
        >>> round(perror_to_phred(0.001), 4)
        30.0
    """
    if p <= 0 or p > 1:
        raise ValueError(f"p must be in (0, 1]; got {p}")
    # TODO: implement (it is one line - math.log10).
    raise NotImplementedError("Implement perror_to_phred")


def per_position_mean_quality(fastq_path: Path) -> list[float]:
    """Return mean Phred quality at each base position across the FASTQ file.

    For each position i:
        mean_q[i] = mean of qualities[i] over all reads that are at least
                    (i+1) bases long.

    The length of the returned list equals the longest read in the file.
    Position 0 in the list corresponds to the first base of the read.

    Empty file -> empty list.

    Implementation hint:
      - Open with SeqIO.parse(fastq_path, "fastq").
      - For each record, grab record.letter_annotations["phred_quality"].
      - Accumulate sums and counts into parallel lists, growing them as
        we see longer reads.
      - Divide at the end.
    """
    # TODO: implement.
    raise NotImplementedError("Implement per_position_mean_quality")


def save_quality_plot(mean_quals: list[float], out_path: Path) -> None:
    """Save a line plot of mean Phred quality vs base position.

    The x axis is 1-based base position (so position 1 is the first base,
    matching how FastQC labels reads).
    The y axis is Phred quality, with a horizontal reference line at Q20
    (the "useful for variant calling" rule of thumb).
    """
    if not mean_quals:
        # Nothing to plot. Still create an (empty) figure so callers don't
        # error on missing file.
        fig, ax = plt.subplots()
        ax.set_title("Per-base mean Phred quality (no data)")
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return

    positions = list(range(1, len(mean_quals) + 1))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(positions, mean_quals, linewidth=1.6)
    ax.axhline(20, linestyle="--", linewidth=1, alpha=0.6, label="Q20")
    ax.set_xlabel("Base position (1-based)")
    ax.set_ylabel("Mean Phred quality")
    ax.set_title("Per-base mean Phred quality across all reads")
    ax.set_ylim(0, max(45, max(mean_quals) + 2))
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-02-fastq-quality-plot.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # 1) Math helpers.
    assert math.isclose(phred_to_perror(10), 0.1, abs_tol=1e-9)
    assert math.isclose(phred_to_perror(20), 0.01, abs_tol=1e-9)
    assert math.isclose(phred_to_perror(30), 0.001, abs_tol=1e-9)
    assert math.isclose(phred_to_perror(0), 1.0, abs_tol=1e-9)

    assert math.isclose(perror_to_phred(0.1), 10.0, abs_tol=1e-9)
    assert math.isclose(perror_to_phred(0.01), 20.0, abs_tol=1e-9)
    assert math.isclose(perror_to_phred(0.001), 30.0, abs_tol=1e-9)
    assert math.isclose(perror_to_phred(1.0), 0.0, abs_tol=1e-9)

    try:
        perror_to_phred(0.0)
    except ValueError:
        pass
    else:
        raise AssertionError("perror_to_phred(0) must raise ValueError")

    # 2) per_position_mean_quality on the inline FASTQ.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        fq = tmp_dir / "test.fastq"
        png = tmp_dir / "mean_quality.png"
        fq.write_text(TEST_FASTQ)

        means = per_position_mean_quality(fq)

        # Longest read is 16 bp -> 16 positions.
        assert len(means) == 16, f"expected 16 positions, got {len(means)}"

        # Position 1: all four reads contribute.
        #   read1 Q40, read2 Q40, read3 Q40, read4 Q10
        #   mean = (40+40+40+10) / 4 = 32.5
        assert math.isclose(means[0], 32.5, abs_tol=1e-9), (
            f"position 1 mean wrong: {means[0]} (expected 32.5)"
        )

        # Position 8 (index 7): read1 Q40, read2 Q40, read3 Q40, read4 Q0
        #   mean = (40+40+40+0) / 4 = 30.0
        assert math.isclose(means[7], 30.0, abs_tol=1e-9), (
            f"position 8 mean wrong: {means[7]} (expected 30.0)"
        )

        # Position 9 (index 8): read1 done, others contribute.
        #   read2 Q40, read3 Q30  ->  mean = 35.0
        assert math.isclose(means[8], 35.0, abs_tol=1e-9), (
            f"position 9 mean wrong: {means[8]} (expected 35.0)"
        )

        # Position 16 (index 15): only read3 reached this far.
        #   read3 last base Q10
        assert math.isclose(means[15], 10.0, abs_tol=1e-9), (
            f"position 16 mean wrong: {means[15]} (expected 10.0)"
        )

        # Sanity: the per-position mean is non-increasing in our toy case
        # (we designed the reads to decay at the 3' end).
        # Allow equality at the start because read1/2/3 share their head.
        for i in range(1, len(means)):
            assert means[i] <= means[i - 1] + 1e-9, (
                f"means should not increase at position {i+1}: "
                f"{means[i-1]} -> {means[i]}"
            )

        # 3) The plot path: produce a PNG file > 1 KB.
        save_quality_plot(means, png)
        assert png.exists() and png.stat().st_size > 1000, (
            f"plot not saved or too small: {png.stat().st_size if png.exists() else 'missing'}"
        )

        # 4) The empty-file path: must not crash.
        empty_fq = tmp_dir / "empty.fastq"
        empty_fq.write_text("")
        empty_means = per_position_mean_quality(empty_fq)
        assert empty_means == [], f"empty FASTQ should yield empty list, got {empty_means!r}"

        empty_png = tmp_dir / "empty.png"
        save_quality_plot(empty_means, empty_png)
        assert empty_png.exists(), "empty-data PNG should still be produced"

    # 5) Bio.SeqIO actually parsed FASTQ for us - confirm by re-reading the
    #    inline FASTQ via a string handle.
    handle = io.StringIO(TEST_FASTQ)
    records = list(SeqIO.parse(handle, "fastq"))
    assert len(records) == 4
    # First read: all I's = Q40 across 8 bases.
    quals = records[0].letter_annotations["phred_quality"]
    assert quals == [40] * 8, quals
    # Per-read mean.
    mean_r0 = statistics.fmean(quals)
    assert math.isclose(mean_r0, 40.0, abs_tol=1e-9)

    print("All assertions passed.")
    print("Inspect mean_quality.png yourself on a real FASTQ to see read decay.")
    print("Move on to Exercise 3 - quality filtering.")
