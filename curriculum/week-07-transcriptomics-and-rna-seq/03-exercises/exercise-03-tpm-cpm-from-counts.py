"""
Exercise 3 - Compute CPM and TPM from a per-transcript counts TSV.

Goal: read a kallisto abundance.tsv (or any equivalent per-transcript
counts table with length information), compute CPM and TPM by hand
without relying on a precomputed TPM column, and verify the result
against kallisto's own TPM column. This exercise pins down the
normalization arithmetic that Lecture 3 introduced.

Estimated time: 45 minutes (15 minutes setup, 20 minutes implementation,
10 minutes verification).

Acceptance criteria:
- `python exercise-03-tpm-cpm-from-counts.py` runs without crashing.
- All `assert` checks at the bottom pass.
- The output TSV at results/tpm_cpm.tsv exists with seven columns:
  target_id, length, eff_length, est_counts, kallisto_tpm,
  computed_cpm, computed_tpm.
- The per-row difference between `computed_tpm` and `kallisto_tpm`
  is < 0.5% for all transcripts with `est_counts > 10`.
- You implemented four functions: `load_abundance`, `compute_cpm`,
  `compute_tpm`, and `verify_tpm_identity`.

Requirements:
    pip install pandas numpy
    (no internet connection needed; this is pure arithmetic on a TSV)

What you learn:
- The CPM formula: CPM_g = 10^6 * count_g / library_size.
- The TPM formula: TPM_g = 10^6 * (count_g / eff_length_g)
                          / sum_g'(count_g' / eff_length_g').
- Why TPM divides by eff_length and CPM does not.
- The TPM identity: sum_g(TPM_g) = 10^6, by construction.
- How to verify a hand-computed TPM against kallisto's own TPM column.

TO COMPLETE: implement the four functions below. Run the file; all
assertions must pass.

Tool versions assumed:
- Python 3.11+
- pandas 2.2+
- numpy 1.26+
- kallisto 0.50.1 (only used to produce the input TSV; this script
  does not call kallisto, it just reads the TSV)

Reference: Wagner, Kin, Lynch 2012, "Measurement of mRNA abundance
using RNA-seq data: RPKM measure is inconsistent among samples,"
Theory in Biosciences 131:281. The cleanest TPM derivation in the
literature.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# Path to the kallisto abundance.tsv produced by Exercise 2. If you
# ran Exercise 2 in week-07/exercises/, this path is correct. Adjust
# if your layout differs.
DEFAULT_ABUNDANCE_TSV = Path(__file__).parent / "quant" / "SRR453568" / "abundance.tsv"

# Where the output TSV goes.
DEFAULT_OUTPUT_TSV = Path(__file__).parent / "results" / "tpm_cpm.tsv"

# Floating-point tolerance for the TPM identity check.
TPM_SUM_TOLERANCE = 1.0  # |sum(TPM) - 1e6| must be below this.

# Per-row tolerance for the kallisto-vs-computed-TPM comparison.
# Expressed as a fraction (0.005 = 0.5%).
TPM_ROW_RELATIVE_TOLERANCE = 0.005


def load_abundance(abundance_tsv: Path) -> pd.DataFrame:
    """Read a kallisto abundance.tsv (or any equivalent) into pandas.

    The expected columns are:
        target_id    str    transcript ID
        length       int    transcript length in bp
        eff_length   float  effective length (length - mean_fragment_len + 1)
        est_counts   float  per-transcript estimated counts (may be fractional)
        tpm          float  kallisto's own TPM (the "ground truth" we verify against)

    Args:
        abundance_tsv: path to the kallisto abundance.tsv.

    Returns:
        DataFrame with the five expected columns.

    Raises:
        FileNotFoundError: if the TSV is not present.
        ValueError: if the TSV does not have the expected columns.
    """
    if not abundance_tsv.exists():
        raise FileNotFoundError(
            f"abundance.tsv not found at {abundance_tsv}. Run Exercise 2 first."
        )

    df = pd.read_csv(abundance_tsv, sep="\t")
    required = {"target_id", "length", "eff_length", "est_counts", "tpm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"abundance.tsv missing required columns: {sorted(missing)}; "
            f"present columns are {sorted(df.columns)}."
        )

    df["length"] = df["length"].astype(int)
    df["eff_length"] = df["eff_length"].astype(float)
    df["est_counts"] = df["est_counts"].astype(float)
    df["tpm"] = df["tpm"].astype(float)
    df = df.rename(columns={"tpm": "kallisto_tpm"})
    return df


def compute_cpm(counts: pd.Series) -> pd.Series:
    """Compute counts per million (CPM) for a sample.

    CPM_g = 10^6 * count_g / library_size

    where library_size = sum over all genes of count_g for this sample.

    CPM corrects for library size only. It does NOT correct for gene
    length, so CPM is not appropriate for cross-gene comparison within
    a single sample. Use TPM for that.

    Args:
        counts: per-transcript (or per-gene) raw counts.

    Returns:
        Series of CPM values aligned to the input index.

    Raises:
        ValueError: if `counts` is empty or sum is zero.
    """
    if len(counts) == 0:
        raise ValueError("counts is empty; cannot compute CPM.")
    library_size: float = float(counts.sum())
    if library_size <= 0.0:
        raise ValueError(
            f"library_size = {library_size} is non-positive; "
            "the sample has no assigned reads."
        )
    cpm: pd.Series = 1e6 * counts / library_size
    return cpm


def compute_tpm(counts: pd.Series, eff_length: pd.Series) -> pd.Series:
    """Compute transcripts per million (TPM) for a sample.

    TPM_g = 10^6 * (count_g / eff_length_g)
                  / sum_g'(count_g' / eff_length_g')

    TPM corrects for both library size and effective transcript length.
    The per-sample sum of TPM is always 10^6, by construction. This
    "sum-to-million" property is what makes TPM the correct cross-sample
    visualization metric.

    Args:
        counts:     per-transcript counts (must align with eff_length).
        eff_length: per-transcript effective length in bp.

    Returns:
        Series of TPM values aligned to the input index.

    Raises:
        ValueError: if the inputs do not align or have zeros where they
                    cannot be tolerated.
    """
    if len(counts) != len(eff_length):
        raise ValueError(
            f"counts and eff_length differ in length "
            f"({len(counts)} vs {len(eff_length)})."
        )
    if (eff_length <= 0.0).any():
        # Floor any non-positive effective lengths at a small positive
        # value to avoid division-by-zero. Real kallisto outputs do not
        # have non-positive eff_length, but defensive code is cheap.
        eff_length = eff_length.where(eff_length > 0.0, other=1.0)

    rate: pd.Series = counts / eff_length
    total_rate: float = float(rate.sum())
    if total_rate <= 0.0:
        raise ValueError(
            "Sum of count/eff_length is zero; the sample has no signal."
        )
    tpm: pd.Series = 1e6 * rate / total_rate
    return tpm


def verify_tpm_identity(tpm: pd.Series, tol: float = TPM_SUM_TOLERANCE) -> None:
    """Verify the identity sum(TPM) = 1e6.

    This is true by construction; if it fails, the implementation of
    compute_tpm has a bug.

    Args:
        tpm: Series of TPM values.
        tol: absolute tolerance for |sum(tpm) - 1e6|.

    Raises:
        AssertionError: if the identity does not hold within `tol`.
    """
    s: float = float(tpm.sum())
    diff: float = abs(s - 1e6)
    if diff > tol:
        raise AssertionError(
            f"TPM identity failed: sum(TPM) = {s:.4f}, "
            f"|sum - 1e6| = {diff:.4f}, tolerance = {tol:.4f}."
        )


def main(abundance_tsv: Path, output_tsv: Path) -> dict[str, Any]:
    """End-to-end driver: load abundance, compute CPM and TPM,
    verify, save, return a small summary dict.
    """
    print(f"[exercise-03] Loading {abundance_tsv} ...")
    df = load_abundance(abundance_tsv)
    print(f"[exercise-03]   loaded {len(df):,} transcripts.")

    print("[exercise-03] Computing CPM ...")
    df["computed_cpm"] = compute_cpm(df["est_counts"])

    print("[exercise-03] Computing TPM ...")
    df["computed_tpm"] = compute_tpm(df["est_counts"], df["eff_length"])

    print("[exercise-03] Verifying sum(TPM) = 1e6 identity ...")
    verify_tpm_identity(df["computed_tpm"])

    print("[exercise-03] Comparing computed TPM to kallisto TPM ...")
    expressed = df[df["est_counts"] > 10].copy()
    expressed["abs_diff"] = (expressed["computed_tpm"] - expressed["kallisto_tpm"]).abs()
    # Avoid division-by-zero in the relative tolerance check.
    expressed["rel_diff"] = expressed["abs_diff"] / expressed["kallisto_tpm"].clip(lower=1e-9)
    max_rel_diff: float = float(expressed["rel_diff"].max())
    median_rel_diff: float = float(expressed["rel_diff"].median())

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    out = df[
        [
            "target_id",
            "length",
            "eff_length",
            "est_counts",
            "kallisto_tpm",
            "computed_cpm",
            "computed_tpm",
        ]
    ]
    out.to_csv(output_tsv, sep="\t", index=False, float_format="%.4f")
    print(f"[exercise-03]   wrote {output_tsv}")

    cpm_sum: float = float(df["computed_cpm"].sum())
    tpm_sum: float = float(df["computed_tpm"].sum())

    return {
        "n_transcripts": int(len(df)),
        "n_expressed": int((df["est_counts"] > 10).sum()),
        "cpm_sum": cpm_sum,
        "tpm_sum": tpm_sum,
        "max_rel_diff": max_rel_diff,
        "median_rel_diff": median_rel_diff,
        "output_tsv": str(output_tsv),
    }


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-03-tpm-cpm-from-counts.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    abundance_tsv = DEFAULT_ABUNDANCE_TSV
    output_tsv = DEFAULT_OUTPUT_TSV

    if not abundance_tsv.exists():
        # Allow the user to override via argv. If neither is available,
        # fall back to a tiny synthetic dataset so the script can still
        # exercise the arithmetic without Exercise 2's output.
        if len(sys.argv) >= 2:
            abundance_tsv = Path(sys.argv[1])
        if not abundance_tsv.exists():
            print(
                f"[exercise-03] WARNING: {abundance_tsv} not found. "
                "Falling back to synthetic test data for arithmetic checks."
            )
            synthetic = pd.DataFrame(
                {
                    "target_id": ["t1", "t2", "t3", "t4", "t5"],
                    "length": [5000, 2000, 1000, 500, 300],
                    "eff_length": [4800.0, 1800.0, 800.0, 300.0, 100.0],
                    "est_counts": [1000.0, 2000.0, 1500.0, 500.0, 100.0],
                    "kallisto_tpm": [0.0, 0.0, 0.0, 0.0, 0.0],
                }
            )
            # Build a synthetic ground-truth TPM that matches our formula
            # so the cross-check passes on the fallback.
            rate = synthetic["est_counts"] / synthetic["eff_length"]
            synthetic["kallisto_tpm"] = 1e6 * rate / rate.sum()
            # Inject the synthetic table into a temp TSV so the rest of
            # the script can run unchanged.
            fallback_tsv = Path(__file__).parent / "_fallback_abundance.tsv"
            synthetic.rename(columns={"kallisto_tpm": "tpm"}).to_csv(
                fallback_tsv, sep="\t", index=False
            )
            abundance_tsv = fallback_tsv

    summary = main(abundance_tsv, output_tsv)

    print()
    print("[exercise-03] Summary:")
    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6f}")
        else:
            print(f"  {k}: {v}")

    # ------------------------------------------------------------------
    # Assertions.
    # ------------------------------------------------------------------
    assert summary["n_transcripts"] > 0, "no transcripts loaded"
    assert summary["n_expressed"] > 0, "no transcripts above the 10-count threshold"

    # CPM sum is 1e6 by construction.
    assert abs(summary["cpm_sum"] - 1e6) < TPM_SUM_TOLERANCE, (
        f"sum(CPM) = {summary['cpm_sum']:.4f}, expected ~1e6 "
        f"(tolerance {TPM_SUM_TOLERANCE})."
    )

    # TPM sum is 1e6 by construction.
    assert abs(summary["tpm_sum"] - 1e6) < TPM_SUM_TOLERANCE, (
        f"sum(TPM) = {summary['tpm_sum']:.4f}, expected ~1e6 "
        f"(tolerance {TPM_SUM_TOLERANCE})."
    )

    # Computed TPM matches kallisto's own TPM to within tolerance for
    # well-expressed transcripts.
    assert summary["max_rel_diff"] < TPM_ROW_RELATIVE_TOLERANCE, (
        f"max relative TPM difference = {summary['max_rel_diff']:.6f}, "
        f"tolerance = {TPM_ROW_RELATIVE_TOLERANCE}. "
        "Your compute_tpm() may have a bug, or eff_length is wrong."
    )

    print()
    print("[exercise-03] All assertions passed.")
    print("[exercise-03] Output:")
    print(f"[exercise-03]   {output_tsv}")
    print(
        "[exercise-03] Continue to the SOLUTIONS.md write-up if you want to compare your "
        "implementation to the reference."
    )
