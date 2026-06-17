"""
Exercise 1 - Build a minimal Snakemake rule for the Salmon-quant step
and render the DAG.

Educational and research use only. Outputs are not clinical software.

Goal: take a sample-sheet TSV and a directory of paired-end FASTQ files,
write a Snakefile with three rules (salmon_index, salmon_quant, aggregate),
and produce (a) a per-sample Salmon quant directory and (b) a wide
transcripts-per-million (TPM) matrix as a TSV. The exercise wraps the
Snakefile with a small Python driver that renders the DAG via subprocess
and writes a run-info.json.

The exercise covers:

- Writing a minimal Snakefile with `rule all`, wildcards, and `expand`.
- Calling `snakemake --dag | dot -Tsvg > dag.svg` via subprocess and
  saving the result.
- Running Snakemake in dry-run mode (`-n`) to verify the DAG without
  executing the shell commands.
- Capturing the Snakemake version, the rule count, and the DAG SVG
  hash in a run-info.json.

Estimated time: 90 minutes (30 min reading, 50 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-01-snakemake-rna-seq.py
    --sample-sheet data/samples.tsv
    --fastq-dir data/fastq
    --transcriptome ref/transcripts.fa.gz
    --out-dir results/ex01` runs end to end without errors when
  snakemake and graphviz `dot` are on the PATH; gracefully skips with
  a clear message otherwise.
- `results/ex01/Snakefile` exists with at least three rules.
- `results/ex01/dag.svg` exists and is a valid SVG file (XML root tag
  is `<svg ...>`).
- `results/ex01/run-info.json` exists with the Snakemake version,
  the rule count, the wildcard values resolved from the sample sheet,
  the DAG SVG SHA-256, and the run date.
- Re-running the script writes the same `dag.svg` byte-for-byte (the
  DAG is deterministic given the same inputs).

Requirements:
    conda install -c bioconda snakemake=7.32.4 salmon=1.10.2
    conda install -c conda-forge graphviz=10.0.1

What you learn:
- The Snakemake rule grammar at the minimum-viable level.
- `expand()` for wildcards, `rule all` as the default target.
- DAG rendering via subprocess: `snakemake --dag | dot -Tsvg > dag.svg`.
- The "verify with a dry-run before you execute" reproducibility pattern.

Tool versions assumed:
- Python 3.11+
- snakemake 7.32.4 (CLI tool; subprocess.run it)
- salmon 1.10.2 (CLI; not actually run by this exercise, only declared
  in the rule's shell line so the DAG resolves)
- graphviz 10.0.1 (for the `dot` binary)

References:
- Snakemake: Mölder et al. 2021, F1000Research 10:33
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8114187/
- Salmon: Patro et al. 2017, Nat Methods 14:417
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/

Author: C10 Crunching Bioinformatics, Week 12
License: MIT (code), CC-BY-4.0 (documentation), see repository LICENSE
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


SNAKEFILE_TEMPLATE: str = '''\
# Snakefile - C10 Week 12 Exercise 1 (minimal RNA-seq quantification)
#
# Educational and research use only. Outputs are not clinical software.
#
# Rules: salmon_index, salmon_quant, aggregate_tpm.
# Inputs: paired-end FASTQ files in {fastq_dir}, transcriptome FASTA
#         at {transcriptome}, sample sheet at {sample_sheet}.
# Output: per-sample Salmon quant directories under quants/, wide TPM
#         matrix at results/tpm_matrix.tsv.

import csv

SAMPLES = []
with open("{sample_sheet}") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        SAMPLES.append(row["sample_id"])

rule all:
    input:
        "results/tpm_matrix.tsv"

rule salmon_index:
    input:
        transcriptome = "{transcriptome}"
    output:
        index_dir = directory("ref/salmon_index")
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        "salmon index -t {{input.transcriptome}} -i {{output.index_dir}} -k 31 -p {{threads}}"

rule salmon_quant:
    input:
        index_dir = "ref/salmon_index",
        r1 = "{fastq_dir}/{{sample}}_R1.fastq.gz",
        r2 = "{fastq_dir}/{{sample}}_R2.fastq.gz"
    output:
        quant_dir = directory("quants/{{sample}}")
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        "salmon quant -i {{input.index_dir}} -l A "
        "-1 {{input.r1}} -2 {{input.r2}} "
        "-o {{output.quant_dir}} -p {{threads}} "
        "--seqBias --gcBias --validateMappings"

rule aggregate_tpm:
    input:
        quants = expand("quants/{{sample}}", sample=SAMPLES)
    output:
        tpm_matrix = "results/tpm_matrix.tsv"
    run:
        import csv, os
        tpm = {{}}
        for quant_dir in input.quants:
            sample = os.path.basename(quant_dir)
            path = os.path.join(quant_dir, "quant.sf")
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\\t")
                for row in reader:
                    tpm.setdefault(row["Name"], {{}})[sample] = row["TPM"]
        with open(output.tpm_matrix, "w") as fh:
            writer = csv.writer(fh, delimiter="\\t")
            writer.writerow(["transcript_id"] + SAMPLES)
            for tx in sorted(tpm):
                writer.writerow([tx] + [tpm[tx].get(s, "0.0") for s in SAMPLES])
'''


def sha256_of(path: Path) -> str:
    """Return the SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def check_tool_available(name: str) -> bool:
    """Return True if the named CLI tool is on the PATH."""
    return shutil.which(name) is not None


def read_sample_sheet(path: Path) -> list[str]:
    """Read sample_id values from a TSV with a sample_id header column."""
    import csv

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "sample_id" not in reader.fieldnames:
            raise ValueError(
                f"Sample sheet at {path} is missing a 'sample_id' column. "
                f"Saw columns: {reader.fieldnames!r}."
            )
        samples = [row["sample_id"] for row in reader if row["sample_id"].strip()]
    if not samples:
        raise ValueError(f"Sample sheet at {path} produced zero sample IDs.")
    return samples


def write_snakefile(
    out_dir: Path,
    sample_sheet: Path,
    fastq_dir: Path,
    transcriptome: Path,
) -> Path:
    """Render the Snakefile template into out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    snakefile = out_dir / "Snakefile"
    text = SNAKEFILE_TEMPLATE.format(
        sample_sheet=sample_sheet,
        fastq_dir=fastq_dir,
        transcriptome=transcriptome,
    )
    snakefile.write_text(text)
    return snakefile


def snakemake_version() -> str:
    """Capture the snakemake --version output, or 'unavailable' if missing."""
    if not check_tool_available("snakemake"):
        return "unavailable"
    try:
        proc = subprocess.run(
            ["snakemake", "--version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return proc.stdout.strip() or proc.stderr.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        return f"unavailable (error: {exc})"


def render_dag(snakefile: Path, out_dir: Path) -> Path:
    """Render the DAG via `snakemake --dag | dot -Tsvg > dag.svg`."""
    dag_svg = out_dir / "dag.svg"
    cmd_snakemake = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--dag",
        "--directory",
        str(out_dir),
    ]
    cmd_dot = ["dot", "-Tsvg"]
    snake_proc = subprocess.run(
        cmd_snakemake,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    if not snake_proc.stdout.strip():
        raise RuntimeError(
            "snakemake --dag produced empty output; the Snakefile may not parse."
        )
    dot_proc = subprocess.run(
        cmd_dot,
        check=True,
        input=snake_proc.stdout,
        capture_output=True,
        text=True,
        timeout=60,
    )
    dag_svg.write_text(dot_proc.stdout)
    return dag_svg


def dry_run_snakemake(snakefile: Path, out_dir: Path) -> str:
    """Run `snakemake -n` and return its stdout for inspection."""
    proc = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(out_dir),
            "-n",
            "--quiet",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    return proc.stdout + proc.stderr


def count_rules(snakefile: Path) -> int:
    """Count `rule <name>:` lines in the Snakefile."""
    text = snakefile.read_text()
    return sum(1 for line in text.splitlines() if line.startswith("rule "))


def write_run_info(
    out_dir: Path,
    snakefile: Path,
    dag_svg: Path | None,
    samples: list[str],
    sample_sheet: Path,
    fastq_dir: Path,
    transcriptome: Path,
    snakemake_ver: str,
    rule_count: int,
    dry_run_log: str,
    skip_reason: str | None,
) -> Path:
    """Write the canonical run-info.json next to the Snakefile."""
    info: dict[str, Any] = {
        "project": "C10-Week12-Exercise-01",
        "version": "v1.0",
        "run_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "workflow_manager": "snakemake",
        "workflow_manager_version": snakemake_ver,
        "workflow_file": str(snakefile.name),
        "workflow_file_sha256": sha256_of(snakefile),
        "rule_count": rule_count,
        "samples_resolved": samples,
        "inputs": {
            "sample_sheet": str(sample_sheet),
            "fastq_dir": str(fastq_dir),
            "transcriptome": str(transcriptome),
        },
        "dag_svg": str(dag_svg.name) if dag_svg is not None else None,
        "dag_svg_sha256": sha256_of(dag_svg) if dag_svg is not None else None,
        "dry_run_log_lines": dry_run_log.count("\n"),
        "skip_reason": skip_reason,
        "host": {
            "hostname": os.uname().nodename,
            "system": os.uname().sysname,
            "release": os.uname().release,
        },
        "license": "MIT for code, CC-BY-4.0 for documentation",
        "disclaimer": "Educational and research use only. Not validated for clinical use.",
    }
    path = out_dir / "run-info.json"
    path.write_text(json.dumps(info, indent=2, sort_keys=True))
    return path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Generate a minimal Snakefile for paired-end RNA-seq Salmon "
            "quantification, render the DAG to SVG, and emit a run-info JSON."
        )
    )
    p.add_argument(
        "--sample-sheet",
        type=Path,
        required=True,
        help="TSV with at minimum a 'sample_id' column.",
    )
    p.add_argument(
        "--fastq-dir",
        type=Path,
        required=True,
        help="Directory containing {sample}_R1.fastq.gz / {sample}_R2.fastq.gz.",
    )
    p.add_argument(
        "--transcriptome",
        type=Path,
        required=True,
        help="Transcriptome FASTA (e.g. GENCODE v44 transcripts.fa.gz).",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Directory to receive Snakefile, dag.svg, and run-info.json.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: read the sample sheet.
    if not args.sample_sheet.exists():
        print(
            f"[ex01] Sample sheet does not exist: {args.sample_sheet}",
            file=sys.stderr,
        )
        return 2
    samples = read_sample_sheet(args.sample_sheet)
    print(f"[ex01] Resolved {len(samples)} samples from {args.sample_sheet}.")
    for s in samples:
        print(f"[ex01]   sample_id = {s}")

    # Step 2: write the Snakefile.
    snakefile = write_snakefile(
        out_dir=out_dir,
        sample_sheet=args.sample_sheet,
        fastq_dir=args.fastq_dir,
        transcriptome=args.transcriptome,
    )
    rule_count = count_rules(snakefile)
    print(f"[ex01] Wrote {snakefile} with {rule_count} rules.")

    # Step 3: render the DAG (skip gracefully if tools are missing).
    skip_reason: str | None = None
    dag_svg: Path | None = None
    dry_run_log: str = ""

    if not check_tool_available("snakemake"):
        skip_reason = "snakemake CLI not on PATH"
        print(f"[ex01] Skipping DAG render: {skip_reason}.")
    elif not check_tool_available("dot"):
        skip_reason = "graphviz `dot` not on PATH"
        print(f"[ex01] Skipping DAG render: {skip_reason}.")
    else:
        try:
            dry_run_log = dry_run_snakemake(snakefile, out_dir)
            print("[ex01] Dry-run output (truncated):")
            for line in dry_run_log.splitlines()[:20]:
                print(f"[ex01]   {line}")
            dag_svg = render_dag(snakefile, out_dir)
            svg_text = dag_svg.read_text()
            if "<svg" not in svg_text:
                raise RuntimeError(
                    "Rendered DAG does not appear to be a valid SVG."
                )
            print(f"[ex01] Wrote {dag_svg} (size {len(svg_text)} chars).")
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            skip_reason = f"DAG render failed: {exc}"
            print(f"[ex01] {skip_reason}")
            dag_svg = None

    # Step 4: write run-info.json.
    snakemake_ver = snakemake_version()
    info_path = write_run_info(
        out_dir=out_dir,
        snakefile=snakefile,
        dag_svg=dag_svg,
        samples=samples,
        sample_sheet=args.sample_sheet,
        fastq_dir=args.fastq_dir,
        transcriptome=args.transcriptome,
        snakemake_ver=snakemake_ver,
        rule_count=rule_count,
        dry_run_log=dry_run_log,
        skip_reason=skip_reason,
    )
    print(f"[ex01] Wrote {info_path}.")

    # Step 5: summary.
    print("[ex01] Summary:")
    print(f"[ex01]   snakemake version: {snakemake_ver}")
    print(f"[ex01]   rule count:        {rule_count}")
    print(f"[ex01]   samples:           {len(samples)}")
    print(f"[ex01]   dag.svg:           {'OK' if dag_svg else 'skipped'}")
    print(f"[ex01]   run-info.json:     {info_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
