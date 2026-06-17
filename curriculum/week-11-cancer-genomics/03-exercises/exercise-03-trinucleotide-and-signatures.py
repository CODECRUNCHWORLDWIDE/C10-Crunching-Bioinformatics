"""
Exercise 3 - Compute the 96-class trinucleotide spectrum and decompose
into COSMIC v3.3 mutational signatures.

Educational and research use only. Do not use any output of this pipeline
to guide patient care.

Goal: take the FilterMutectCalls PASS VCF from Exercise 2, extract the
trinucleotide context for each SNV (with pyrimidine normalization),
build the 96-class spectrum as a TSV, and run SigProfilerAssignment to
decompose against the COSMIC v3.3 SBS catalog. Emit a Markdown summary
of the top signatures and the cosine similarity.

The exercise covers:

- Reading a VCF with pysam, filtering for SNV PASS variants.
- Reading the reference FASTA with pysam.FastaFile for trinucleotide context.
- Pyrimidine normalization (reverse-complement A/G reference to C/T).
- Building the 96-class spectrum in the canonical Alexandrov order.
- Calling SigProfilerAssignment.cosmic_fit() and parsing the output.
- Computing the cosine similarity by hand as a sanity check.
- Writing a Markdown summary alongside the run-info JSON.

Estimated time: 105 minutes (35 min reading, 55 min implementing,
15 min running).

Acceptance criteria:
- `python exercise-03-trinucleotide-and-signatures.py
    --filtered-vcf results/ex02/filtered.vcf.gz
    --reference data/chr22_GRCh38.fasta
    --out-dir results/ex03` runs end to end when SigProfilerAssignment
  is installed and the GRCh38 SigProfiler reference is set up; gracefully
  skips otherwise.
- `results/ex03/spectrum_96.tsv` exists with all 96 rows.
- `results/ex03/signature_summary.md` exists.
- `results/ex03/run-info.json` records the COSMIC version, the top
  signatures with their fractional contributions, and the cosine similarity.

Requirements:
    conda install -c bioconda sigprofilerassignment=0.1.4
        sigprofilermatrixgenerator=1.2.26 pysam=0.22.1
    # one-time:
    python -c "from SigProfilerMatrixGenerator import install; install.install('GRCh38')"

What you learn:
- The pyrimidine-normalization convention for mutational signatures.
- The 96-class spectrum data structure and canonical ordering.
- The SigProfilerAssignment.cosmic_fit() Python API.
- The cosine-similarity goodness-of-fit measure.

Tool versions assumed:
- Python 3.11+
- pysam 0.22.1
- SigProfilerAssignment 0.1.4
- SigProfilerMatrixGenerator 1.2.26

References:
- COSMIC signatures: Alexandrov et al. 2020, Nature 578:94
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7054213/
- SigProfilerAssignment:
  https://github.com/AlexandrovLab/SigProfilerAssignment
- deconstructSigs: Rosenthal et al. 2016, Genome Biology 17:31
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4762164/
- COSMIC database: Sondka et al. 2024, NAR 52:D1210
  https://academic.oup.com/nar/article/52/D1/D1210/7416441
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import shutil
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


# ----------------------------------------------------------------------
# Constants.
# ----------------------------------------------------------------------

DEFAULT_GENOME_BUILD: str = "GRCh38"
DEFAULT_COSMIC_VERSION: str = "3.3"
DEFAULT_SAMPLE_NAME: str = "TUMOR"
DEFAULT_TOP_N_SIGNATURES: int = 5

COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# Canonical Alexandrov-lab order: six substitution types in C>A, C>G, C>T,
# T>A, T>C, T>G order; within each, 16 trinucleotides in A < C < G < T order.
SUBSTITUTION_TYPES: list[str] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
FLANK_BASES: list[str] = ["A", "C", "G", "T"]


def canonical_96_classes() -> list[str]:
    """Return the 96 class labels in the canonical Alexandrov order.

    Format: '<left>[<ref>><alt>]<right>', e.g. 'A[C>T]G'.
    """
    classes: list[str] = []
    for sub in SUBSTITUTION_TYPES:
        ref, alt = sub.split(">")
        for left in FLANK_BASES:
            for right in FLANK_BASES:
                classes.append(f"{left}[{ref}>{alt}]{right}")
    return classes


# ----------------------------------------------------------------------
# Provenance.
# ----------------------------------------------------------------------

@dataclass
class ExerciseRunInfo:
    """Provenance metadata for the Exercise 3 run."""
    run_date: str = ""
    filtered_vcf: str = ""
    reference_fasta: str = ""
    out_dir: str = ""
    sample_name: str = DEFAULT_SAMPLE_NAME
    genome_build: str = DEFAULT_GENOME_BUILD
    cosmic_version: str = DEFAULT_COSMIC_VERSION
    pysam_version: str = ""
    sigprofiler_version: str = ""
    n_snvs_used: int = 0
    n_pass_records: int = 0
    n_indels_skipped: int = 0
    top_signatures: list[dict[str, Any]] = field(default_factory=list)
    cosine_similarity: float = 0.0
    skipped: bool = False
    skip_reason: str = ""
    python_version: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------

def pysam_version_string() -> str:
    """Return the installed pysam version, or 'unknown'."""
    try:
        import pysam  # type: ignore[import-not-found]
        return getattr(pysam, "__version__", "unknown")
    except Exception:
        return "unknown"


def sigprofiler_version_string() -> str:
    """Return the installed SigProfilerAssignment version, or 'unknown'."""
    try:
        import SigProfilerAssignment  # type: ignore[import-not-found]
        return getattr(SigProfilerAssignment, "__version__", "unknown")
    except Exception:
        return "unknown"


def python_version_string() -> str:
    """Return the running Python version."""
    return ".".join(str(part) for part in sys.version_info[:3])


# ----------------------------------------------------------------------
# Stage 1 - extract trinucleotide context.
# ----------------------------------------------------------------------

def reverse_complement(s: str) -> str:
    """Return the reverse complement of a DNA string (A/C/G/T/N)."""
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(s))


def normalize_to_pyrimidine(
    ref: str, alt: str, context: str,
) -> tuple[str, str, str]:
    """Normalize a SNV to its pyrimidine-reference representation.

    If ref is C or T: returns (ref, alt, context) unchanged.
    If ref is A or G: complements both alleles and reverse-complements the
    trinucleotide context so the central base is C or T.
    """
    if ref in {"C", "T"}:
        return (ref, alt, context)
    if ref in {"A", "G"}:
        return (
            COMPLEMENT[ref],
            COMPLEMENT[alt],
            reverse_complement(context),
        )
    raise ValueError(f"Unsupported reference base: {ref}")


def trinucleotide_class(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    fasta: Any,
) -> str:
    """Return the 96-class label, e.g. 'A[C>T]G'. 1-based pos.

    fasta must support .fetch(chrom, start_0based, end_0based_exclusive)
    returning a string (pysam.FastaFile.fetch matches this signature).

    Raises ValueError on non-SNV alleles.
    """
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("trinucleotide_class is SNV-only")
    context: str = fasta.fetch(chrom, pos - 2, pos + 1).upper()
    if len(context) != 3:
        raise ValueError(
            f"context fetch returned {len(context)} bases at {chrom}:{pos}"
        )
    if "N" in context:
        raise ValueError(f"context contains N at {chrom}:{pos}: {context}")
    central: str = context[1]
    if central != ref.upper():
        raise ValueError(
            f"reference mismatch at {chrom}:{pos}: BAM says {ref}, FASTA says {central}"
        )
    norm_ref, norm_alt, norm_context = normalize_to_pyrimidine(
        ref.upper(), alt.upper(), context,
    )
    return f"{norm_context[0]}[{norm_ref}>{norm_alt}]{norm_context[2]}"


# ----------------------------------------------------------------------
# Stage 2 - parse PASS SNVs from the filtered VCF.
# ----------------------------------------------------------------------

@dataclass
class SpectrumResult:
    """The 96-class spectrum plus parse counts."""
    counts: dict[str, int]
    n_pass_records: int
    n_snvs_used: int
    n_indels_skipped: int


def build_96_class_spectrum(
    filtered_vcf: Path,
    reference: Path,
) -> SpectrumResult:
    """Walk a filtered Mutect2 VCF and build the 96-class trinucleotide spectrum.

    Skips non-PASS records and indels. Returns SpectrumResult with all 96
    classes present (zeros for unobserved).
    """
    try:
        import pysam  # type: ignore[import-not-found]
    except ImportError as exc:
        raise FileNotFoundError("pysam is not importable") from exc

    counts: dict[str, int] = {cls: 0 for cls in canonical_96_classes()}
    n_pass: int = 0
    n_snvs: int = 0
    n_indels: int = 0

    with pysam.FastaFile(str(reference)) as fasta, pysam.VariantFile(str(filtered_vcf)) as vf:
        for rec in vf:
            filt: list[str] = list(rec.filter.keys())
            if filt and filt != ["PASS"]:
                continue
            n_pass += 1
            for alt in (rec.alts or []):
                if len(rec.ref) != 1 or len(alt) != 1:
                    n_indels += 1
                    continue
                try:
                    cls: str = trinucleotide_class(
                        chrom=rec.chrom,
                        pos=rec.pos,
                        ref=rec.ref,
                        alt=alt,
                        fasta=fasta,
                    )
                except ValueError:
                    continue
                counts[cls] += 1
                n_snvs += 1
    return SpectrumResult(
        counts=counts,
        n_pass_records=n_pass,
        n_snvs_used=n_snvs,
        n_indels_skipped=n_indels,
    )


def write_spectrum_tsv(
    counts: dict[str, int],
    sample_name: str,
    out_path: Path,
) -> None:
    """Write the 96-class spectrum as a SigProfiler-compatible TSV.

    Header: 'Mutation Types\\t<sample_name>'
    Body: one row per of the 96 classes, in canonical order.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write(f"Mutation Types\t{sample_name}\n")
        for cls in canonical_96_classes():
            fh.write(f"{cls}\t{counts.get(cls, 0)}\n")


# ----------------------------------------------------------------------
# Stage 3 - run SigProfilerAssignment.
# ----------------------------------------------------------------------

def run_sigprofiler_assignment(
    spectrum_tsv: Path,
    out_dir: Path,
    genome_build: str = DEFAULT_GENOME_BUILD,
    cosmic_version: str = DEFAULT_COSMIC_VERSION,
) -> Path:
    """Run SigProfilerAssignment.cosmic_fit() against the COSMIC v3.3 catalog.

    Returns the path to the activities TSV (per-signature mutation counts).
    Raises FileNotFoundError if SigProfilerAssignment is not importable.
    """
    try:
        from SigProfilerAssignment import Analyzer  # type: ignore[import-not-found]
    except ImportError as exc:
        raise FileNotFoundError("SigProfilerAssignment is not importable") from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    Analyzer.cosmic_fit(
        samples=str(spectrum_tsv),
        output=str(out_dir),
        input_type="matrix",
        context_type="96",
        genome_build=genome_build,
        cosmic_version=cosmic_version,
        collapse_to_SBS96=True,
        make_plots=False,
    )
    candidates: list[Path] = sorted(out_dir.rglob("Assignment_Solution_Activities.txt"))
    if not candidates:
        raise FileNotFoundError(
            f"SigProfilerAssignment finished but no Activities file found under {out_dir}"
        )
    return candidates[0]


# ----------------------------------------------------------------------
# Stage 4 - parse the activities TSV.
# ----------------------------------------------------------------------

def parse_activities_tsv(
    activities_path: Path,
    top_n: int = DEFAULT_TOP_N_SIGNATURES,
) -> list[dict[str, Any]]:
    """Read the Assignment_Solution_Activities.txt and return top-N signatures.

    The file has samples as rows and signatures as columns:
        Samples  SBS1  SBS2  SBS3  ...
        TUMOR    143   0     38    ...

    Returns a list of dicts with keys: signature, count, fraction.
    """
    if not activities_path.exists():
        return []
    with activities_path.open() as fh:
        header_line: str = fh.readline().rstrip("\n")
        signature_names: list[str] = header_line.split("\t")[1:]
        data_line: str = fh.readline().rstrip("\n")
        parts: list[str] = data_line.split("\t")
        if len(parts) < 2:
            return []
        counts: list[int] = []
        for v in parts[1:]:
            try:
                counts.append(int(round(float(v))))
            except ValueError:
                counts.append(0)
    if not counts:
        return []
    total: int = sum(counts)
    if total <= 0:
        return []
    rows: list[dict[str, Any]] = []
    for name, count in zip(signature_names, counts):
        if count <= 0:
            continue
        rows.append({
            "signature": name,
            "count": count,
            "fraction": count / total,
        })
    rows.sort(key=lambda r: r["count"], reverse=True)
    return rows[:top_n]


# ----------------------------------------------------------------------
# Stage 5 - cosine similarity by hand.
# ----------------------------------------------------------------------

def cosine_similarity(a: list[float], b: list[float]) -> float:
    """Return the cosine similarity between two vectors of equal length.

    Returns 0.0 if either vector is all-zero or has unequal length.
    """
    if len(a) != len(b):
        return 0.0
    dot: float = sum(x * y for x, y in zip(a, b))
    norm_a: float = math.sqrt(sum(x * x for x in a))
    norm_b: float = math.sqrt(sum(y * y for y in b))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return dot / (norm_a * norm_b)


def parse_sample_stats(stats_path: Path) -> float:
    """Read SigProfilerAssignment's Samples_Stats and return the cosine similarity.

    The Solution_Stats file is a TSV with one row per sample; the
    cosine-similarity column varies in name across versions ('Cosine Similarity'
    or 'Cosine_similarity'). Returns 0.0 if the file or column is missing.
    """
    if not stats_path.exists():
        return 0.0
    with stats_path.open() as fh:
        header_line: str = fh.readline().rstrip("\n")
        headers: list[str] = header_line.split("\t")
        cos_idx: int = -1
        for i, h in enumerate(headers):
            if "cosine" in h.lower():
                cos_idx = i
                break
        if cos_idx < 0:
            return 0.0
        for line in fh:
            parts: list[str] = line.rstrip("\n").split("\t")
            if len(parts) > cos_idx:
                try:
                    return float(parts[cos_idx])
                except ValueError:
                    continue
    return 0.0


def find_stats_file(sigprofiler_out_dir: Path) -> Path | None:
    """Locate the Samples_Stats file under the SigProfiler output tree."""
    candidates: list[Path] = sorted(
        sigprofiler_out_dir.rglob("*Samples_Stats*")
    )
    return candidates[0] if candidates else None


# ----------------------------------------------------------------------
# Stage 6 - render Markdown.
# ----------------------------------------------------------------------

def render_signature_markdown(
    out_path: Path,
    n_snvs_used: int,
    top_signatures: list[dict[str, Any]],
    cosine: float,
    cosmic_version: str,
    genome_build: str,
) -> None:
    """Write a Markdown summary of the signature decomposition."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "# Exercise 3 - Mutational Signature Decomposition",
        "",
        "**Educational and research use only.** Not for clinical decisions.",
        "",
        f"- SNVs used (PASS, SNV-only): {n_snvs_used}",
        f"- COSMIC catalog version: v{cosmic_version}",
        f"- Genome build: {genome_build}",
        f"- Reconstructed-spectrum cosine similarity: {cosine:.3f}",
        "",
    ]
    if n_snvs_used < 50:
        lines.extend([
            "> **Warning.** Fewer than 50 SNVs were available. The Alexandrov-lab "
            "recommendation is at least 100 SNVs for a stable decomposition. The "
            "signature attribution below is unstable and should not be quoted "
            "without acknowledging this limitation.",
            "",
        ])
    if cosine < 0.85 and cosine > 0:
        lines.extend([
            "> **Warning.** Cosine similarity below 0.85 indicates a poor fit. The "
            "reported signatures may not adequately explain the observed spectrum.",
            "",
        ])
    lines.append("## Top contributing signatures")
    lines.append("")
    lines.append("| Signature | Mutations | Fraction |")
    lines.append("|-----------|----------:|---------:|")
    for sig in top_signatures:
        lines.append(
            f"| {sig['signature']:<10s} | {sig['count']:>9d} | {sig['fraction']:>7.3f} |"
        )
    lines.extend([
        "",
        "Top-line aetiology hints (consult COSMIC documentation for full details):",
        "",
        "- SBS1: age-related (5-methylcytosine deamination)",
        "- SBS5: clock-like (unknown source)",
        "- SBS3: homologous-recombination deficiency (BRCA-pathway loss)",
        "- SBS4: tobacco smoking",
        "- SBS6 / SBS15 / SBS20 / SBS26: mismatch-repair deficiency",
        "- SBS7a-d: UV exposure",
        "- SBS18: reactive oxygen damage",
        "",
        f"Reference: Alexandrov et al. 2020, Nature 578:94 (PMC: 7054213). "
        f"Reconstructed-spectrum cosine: {cosine:.3f}.",
        "",
    ])
    out_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------
# Stage 7 - run-info JSON.
# ----------------------------------------------------------------------

def write_run_info(run_info: ExerciseRunInfo, out_path: Path) -> None:
    """Write the run-info JSON. Raises ValueError on empty run_date."""
    if not run_info.run_date:
        raise ValueError("run_info.run_date is empty; refusing to write.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(run_info.to_dict(), fh, indent=2, sort_keys=True)
        fh.write("\n")


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------

def run_exercise(
    filtered_vcf: Path,
    reference: Path,
    out_dir: Path,
    sample_name: str = DEFAULT_SAMPLE_NAME,
    genome_build: str = DEFAULT_GENOME_BUILD,
    cosmic_version: str = DEFAULT_COSMIC_VERSION,
    skip_if_missing: bool = True,
) -> Path:
    """Run the full exercise. Returns the path to run-info.json."""
    out_dir.mkdir(parents=True, exist_ok=True)
    spectrum_tsv: Path = out_dir / "spectrum_96.tsv"
    sigprofiler_out: Path = out_dir / "sigprofiler_out"
    summary_md: Path = out_dir / "signature_summary.md"
    run_info_path: Path = out_dir / "run-info.json"

    info = ExerciseRunInfo(
        run_date=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        filtered_vcf=str(filtered_vcf),
        reference_fasta=str(reference),
        out_dir=str(out_dir),
        sample_name=sample_name,
        genome_build=genome_build,
        cosmic_version=cosmic_version,
        pysam_version=pysam_version_string(),
        sigprofiler_version=sigprofiler_version_string(),
        python_version=python_version_string(),
    )

    if not filtered_vcf.exists() or not reference.exists():
        if skip_if_missing:
            info.skipped = True
            info.skip_reason = "filtered VCF or reference not found"
            write_run_info(info, run_info_path)
            print("[ex03] skipped: missing inputs", file=sys.stderr)
            return run_info_path
        raise FileNotFoundError("filtered VCF or reference not found")

    # Stage 1+2 - build the 96-class spectrum.
    try:
        spectrum = build_96_class_spectrum(filtered_vcf, reference)
    except FileNotFoundError as exc:
        if skip_if_missing:
            info.skipped = True
            info.skip_reason = str(exc)
            write_run_info(info, run_info_path)
            print(f"[ex03] skipped: {exc}", file=sys.stderr)
            return run_info_path
        raise

    info.n_pass_records = spectrum.n_pass_records
    info.n_snvs_used = spectrum.n_snvs_used
    info.n_indels_skipped = spectrum.n_indels_skipped

    write_spectrum_tsv(spectrum.counts, sample_name, spectrum_tsv)
    print(
        f"[ex03] built spectrum: {spectrum.n_snvs_used} SNVs across "
        f"{sum(1 for v in spectrum.counts.values() if v > 0)} non-empty classes",
        file=sys.stderr,
    )

    # Stage 3 - run SigProfilerAssignment (optional skip).
    try:
        activities_path: Path = run_sigprofiler_assignment(
            spectrum_tsv=spectrum_tsv,
            out_dir=sigprofiler_out,
            genome_build=genome_build,
            cosmic_version=cosmic_version,
        )
    except FileNotFoundError as exc:
        if skip_if_missing:
            info.skipped = True
            info.skip_reason = str(exc)
            render_signature_markdown(
                out_path=summary_md,
                n_snvs_used=info.n_snvs_used,
                top_signatures=[],
                cosine=0.0,
                cosmic_version=cosmic_version,
                genome_build=genome_build,
            )
            write_run_info(info, run_info_path)
            print(f"[ex03] skipped SigProfiler: {exc}", file=sys.stderr)
            return run_info_path
        raise
    except Exception as exc:
        info.skipped = True
        info.skip_reason = f"SigProfilerAssignment failed: {exc}"
        write_run_info(info, run_info_path)
        print(f"[ex03] SigProfilerAssignment failed: {exc}", file=sys.stderr)
        return run_info_path

    # Stage 4 - parse activities and stats.
    info.top_signatures = parse_activities_tsv(activities_path)
    stats_path = find_stats_file(sigprofiler_out)
    if stats_path is not None:
        info.cosine_similarity = parse_sample_stats(stats_path)

    # Stage 5 - render Markdown.
    render_signature_markdown(
        out_path=summary_md,
        n_snvs_used=info.n_snvs_used,
        top_signatures=info.top_signatures,
        cosine=info.cosine_similarity,
        cosmic_version=cosmic_version,
        genome_build=genome_build,
    )

    write_run_info(info, run_info_path)
    print(
        f"[ex03] wrote {summary_md} with {len(info.top_signatures)} top signatures, "
        f"cosine={info.cosine_similarity:.3f}",
        file=sys.stderr,
    )
    return run_info_path


# ----------------------------------------------------------------------
# CLI.
# ----------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exercise 3 - 96-class spectrum + COSMIC v3 signature decomposition.",
    )
    parser.add_argument("--filtered-vcf", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sample-name", type=str, default=DEFAULT_SAMPLE_NAME)
    parser.add_argument("--genome-build", type=str, default=DEFAULT_GENOME_BUILD)
    parser.add_argument("--cosmic-version", type=str, default=DEFAULT_COSMIC_VERSION)
    parser.add_argument("--no-skip-if-missing", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run_info_path: Path = run_exercise(
        filtered_vcf=args.filtered_vcf,
        reference=args.reference,
        out_dir=args.out_dir,
        sample_name=args.sample_name,
        genome_build=args.genome_build,
        cosmic_version=args.cosmic_version,
        skip_if_missing=not args.no_skip_if_missing,
    )
    print(f"[ex03] wrote {run_info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
