"""
Exercise 3 - Filter FASTQ reads by quality.

Goal: trim each read with a sliding-window quality heuristic, then keep
only reads whose post-trim length is >= 36 bp and whose post-trim mean
quality is >= Q20. Write the survivors out to a new FASTQ. Count how
many we kept vs dropped. This is, in 80 lines, what Trimmomatic does
in 8,000.

Estimated time: 40 minutes.

Acceptance criteria:
- `python exercise-03-filter-by-quality.py` runs without crashing.
- All assertions at the bottom pass.
- Your filter function takes (in_path, out_path, min_len, min_mean_q,
  window, window_threshold) and returns a small report dict.

Required:
    python -m pip install biopython==1.83

What you learn:
- Slicing a SeqRecord with `record[:n]` keeps the per-base qualities in
  sync. Slicing only `record.seq` does NOT. Always slice the record.
- A sliding-window trimmer in 6 lines.
- Why the "min length 36" rule of thumb exists - it is the lower bound
  for confident short-read alignment to a vertebrate genome.

TO COMPLETE: implement `find_trim_position` and `filter_fastq`.
Run the file; all assertions must pass.
"""

from __future__ import annotations

import io
import statistics
import tempfile
from pathlib import Path
from typing import TypedDict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Defaults align with Trimmomatic SLIDINGWINDOW:4:20 + MINLEN:36.
DEFAULT_WINDOW = 4
DEFAULT_WINDOW_THRESHOLD = 20
DEFAULT_MIN_LEN = 36
DEFAULT_MIN_MEAN_Q = 20


class FilterReport(TypedDict):
    n_in: int
    n_kept: int
    n_dropped_short: int
    n_dropped_lowq: int
    mean_len_in: float
    mean_len_out: float


def find_trim_position(quals: list[int], window: int, threshold: int) -> int:
    """Return the position at which to truncate a read.

    Slide a window of size `window` along the qualities list from left
    to right. The first window whose MEAN Phred quality is strictly
    LESS THAN `threshold` defines the trim point: cut everything from
    that window's start (inclusive) onward.

    If no window falls below the threshold, return len(quals) (no trim).

    Edge cases:
      - quals shorter than window -> return len(quals) (cannot evaluate).
      - empty quals -> return 0.

    Examples:
        # First failing window starts at index 3: [40,5,5,5] mean 13.75 < 20.
        >>> find_trim_position([40, 40, 40, 40, 5, 5, 5, 5], 4, 20)
        3
        >>> find_trim_position([40] * 10, 4, 20)
        10
        >>> find_trim_position([], 4, 20)
        0
    """
    # TODO: implement.
    raise NotImplementedError("Implement find_trim_position")


def filter_fastq(
    in_path: Path,
    out_path: Path,
    min_len: int = DEFAULT_MIN_LEN,
    min_mean_q: float = DEFAULT_MIN_MEAN_Q,
    window: int = DEFAULT_WINDOW,
    window_threshold: int = DEFAULT_WINDOW_THRESHOLD,
) -> FilterReport:
    """Trim and filter reads, write survivors, return a report.

    Pipeline per record:
      1) Compute the trim position via `find_trim_position`.
      2) Slice the SeqRecord (NOT just .seq) up to that position.
         The slice syntax `record[:n]` keeps letter_annotations in sync.
      3) If the trimmed length < min_len: drop, count as dropped_short.
      4) Else if the trimmed mean quality < min_mean_q: drop, count as
         dropped_lowq.
      5) Else: emit to the output FASTQ.

    The report dict (TypedDict above) records:
      n_in              total records read
      n_kept            records emitted
      n_dropped_short   records dropped for being too short post-trim
      n_dropped_lowq    records dropped for low mean quality post-trim
      mean_len_in       mean read length BEFORE trimming
      mean_len_out      mean read length of EMITTED reads (0.0 if none)

    Use SeqIO.parse for input and accumulate the survivors before
    calling SeqIO.write. For very large files you would stream the
    write, but for the exercise the simpler API is fine.
    """
    # TODO: implement the pipeline.
    raise NotImplementedError("Implement filter_fastq")


# A tiny inline FASTQ for self-test. Phred+33 encoding. Quality chars:
#   '!' = Q0,  '&' = Q5,  '+' = Q10,  '6' = Q21,  '5' = Q20,
#   '?' = Q30,  'I' = Q40.
#
# We deliberately exercise the filter with NON-default thresholds in the
# self-test below (window_threshold=15, min_mean_q=25, min_len=36) so
# that every branch is reachable. With the production defaults
# (window_threshold=20, min_mean_q=20) the low-quality branch is hard to
# distinguish from the too-short branch because any sub-Q20 plateau will
# also trim the read down to nothing.
#
# read_keep:          40 bp all Q40                  -> kept unmodified
# read_trim_keep:     50 bp, decays to Q5 at pos 40  -> trims to 40 bp Q40, kept
# read_too_short:     20 bp all Q40                  -> dropped (len < 36)
# read_low_q:         40 bp all Q21                  -> survives trim;
#                                                       mean Q21 < min_mean_q=25,
#                                                       dropped for low quality
# read_trim_to_short: 50 bp, decays to Q5 at pos 10  -> trims to 10 bp, short
TEST_FASTQ = (
    "@read_keep len=40, all Q40\n"
    "A" * 40 + "\n"
    "+\n"
    "I" * 40 + "\n"
    "@read_trim_keep len=50, decays at 40\n"
    "C" * 50 + "\n"
    "+\n"
    "I" * 40 + "&" * 10 + "\n"
    "@read_too_short len=20, all Q40\n"
    "G" * 20 + "\n"
    "+\n"
    "I" * 20 + "\n"
    "@read_low_q len=40, all Q21\n"
    "T" * 40 + "\n"
    "+\n"
    "6" * 40 + "\n"
    "@read_trim_to_short len=50, decays at 10\n"
    "A" * 50 + "\n"
    "+\n"
    "I" * 10 + "&" * 40 + "\n"
)


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-filter-by-quality.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # find_trim_position - the unit tests are in the docstring examples,
    # but we re-check a couple here too.
    # [40,40,40,40,5,5,5,5]: first window with mean<20 starts at i=3
    # ([40,5,5,5] mean=13.75).
    assert find_trim_position([40, 40, 40, 40, 5, 5, 5, 5], 4, 20) == 3
    assert find_trim_position([40] * 10, 4, 20) == 10
    assert find_trim_position([], 4, 20) == 0
    # If the threshold is higher than every base, the first window fails.
    assert find_trim_position([30, 30, 30, 30, 30], 4, 35) == 0
    # Window longer than the read - cannot evaluate, no trim.
    assert find_trim_position([40, 40, 40], 4, 20) == 3

    # Now the full filter pipeline on our synthetic FASTQ.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        in_fq = tmp_dir / "in.fastq"
        out_fq = tmp_dir / "out.fastq"
        in_fq.write_text(TEST_FASTQ)

        # Use non-default thresholds so all four branches are reachable.
        # See the TEST_FASTQ comment above for the design.
        report = filter_fastq(
            in_fq,
            out_fq,
            min_len=36,
            min_mean_q=25,
            window=4,
            window_threshold=15,
        )

        # Five reads in.
        assert report["n_in"] == 5, report

        # Two survivors (read_keep, read_trim_keep).
        assert report["n_kept"] == 2, report

        # Two short rejections (read_too_short, read_trim_to_short).
        assert report["n_dropped_short"] == 2, report

        # One low-quality rejection (read_low_q, all Q21, mean < 25).
        assert report["n_dropped_lowq"] == 1, report

        # mean_len_in = (40 + 50 + 20 + 40 + 50) / 5 = 40.0
        assert abs(report["mean_len_in"] - 40.0) < 1e-9, report
        # mean_len_out:
        #   read_keep         -> length 40 (no trim)
        #   read_trim_keep    -> length 39 (sliding window with threshold 15
        #                        triggers when window covers indices 39..42:
        #                        [40,5,5,5] mean 13.75 < 15 -> trim at 39)
        # mean of survivors = (40 + 39) / 2 = 39.5
        assert abs(report["mean_len_out"] - 39.5) < 1e-9, report

        # The output FASTQ exists and has exactly two records.
        survivors = list(SeqIO.parse(out_fq, "fastq"))
        assert len(survivors) == 2, len(survivors)
        survivor_ids = {r.id for r in survivors}
        assert survivor_ids == {"read_keep", "read_trim_keep"}, survivor_ids

        # The trimmed survivor (read_trim_keep) is 39 bp post-trim and is
        # entirely Q40 (the trim point lies inside the I* run).
        trimmed = next(r for r in survivors if r.id == "read_trim_keep")
        assert len(trimmed.seq) == 39, len(trimmed.seq)
        assert all(q == 40 for q in trimmed.letter_annotations["phred_quality"])

        # Round-trip sanity: parse the output ourselves via a StringIO.
        with out_fq.open() as h:
            text = h.read()
        again = list(SeqIO.parse(io.StringIO(text), "fastq"))
        assert len(again) == 2

        # Per-record mean quality of survivors >= the threshold we used (25).
        for r in survivors:
            mq = statistics.fmean(r.letter_annotations["phred_quality"])
            assert mq >= 25, (r.id, mq)

    # Type-check the public function signature is well-formed.
    sample_record = SeqRecord.__init__  # noqa: F841 - just confirms import OK

    print("All assertions passed.")
    print("Your filter behaves like Trimmomatic SLIDINGWINDOW:4:20 MINLEN:36.")
    print("Move on to the Challenge: streaming a large FASTA without loading it.")
