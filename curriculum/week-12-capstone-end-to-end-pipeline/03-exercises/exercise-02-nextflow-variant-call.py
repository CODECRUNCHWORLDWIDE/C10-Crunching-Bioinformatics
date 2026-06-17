"""
Exercise 2 - Write a Nextflow process for the bcftools-call step
and capture the channel topology.

Educational and research use only. Outputs are not clinical software.

Goal: take a sample sheet listing BAM files and a reference FASTA,
write a Nextflow `main.nf` with two processes (CALL_VARIANTS,
MERGE_VCFS), drive Nextflow's `-with-dag dag.svg` to render the DAG,
and produce a run-info.json that captures the Nextflow version, the
channel-to-process wiring, and the DAG hash.

The exercise covers:

- Writing a minimal Nextflow workflow with channels and processes.
- The Nextflow channel grammar: `Channel.fromPath`, `splitCsv`, `map`.
- DAG rendering via `nextflow run main.nf -with-dag dag.svg`.
- Comparing Snakemake and Nextflow at the DAG level: the two managers
  should produce isomorphic DAGs for the same pipeline.
- The "dry-run with `-preview`" idiom (Nextflow 22.10+).

Estimated time: 90 minutes (40 min reading, 40 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-02-nextflow-variant-call.py
    --sample-sheet data/samples.tsv
    --reference ref/chr22.fasta
    --out-dir results/ex02` runs end to end without errors when
  nextflow is on the PATH; gracefully skips with a clear message
  otherwise.
- `results/ex02/main.nf` exists with at least two `process` blocks.
- `results/ex02/nextflow.config` exists with a `process.container`
  directive and a `conda` profile.
- `results/ex02/dag.svg` exists (when Nextflow is available) and
  passes the `<svg` sniff test.
- `results/ex02/run-info.json` exists with the Nextflow version,
  process count, channel count, and DAG SHA-256.
- Re-running the script produces the same `main.nf` byte-for-byte.

Requirements:
    conda install -c bioconda nextflow=23.10.1 bcftools=1.20 samtools=1.20

What you learn:
- The Nextflow process grammar: input, output, script, tag, cpus.
- Channels as the dataflow primitive: how `fromPath`, `splitCsv`, and
  `map` compose to a sample channel.
- The Snakemake vs Nextflow comparison: same DAG, different syntax.
- The `nextflow.config` patterns: profiles, containers, conda.

Tool versions assumed:
- Python 3.11+
- nextflow 23.10.1 (CLI; subprocess.run it)
- bcftools 1.20 (CLI; declared in process script but not invoked here)

References:
- Nextflow: Di Tommaso et al. 2017, Nat Biotechnol 35:316
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8224876/
- bcftools: Danecek et al. 2021, GigaScience 10:giab008
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7931819/

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


MAIN_NF_TEMPLATE: str = r'''/*
 * main.nf - C10 Week 12 Exercise 2 (minimal variant-calling demo)
 *
 * Educational and research use only. Outputs are not clinical software.
 *
 * Processes:
 *   CALL_VARIANTS: bcftools mpileup + bcftools call per sample.
 *   MERGE_VCFS:    bcftools merge across all samples.
 *
 * Inputs:
 *   --sample_sheet : TSV with columns sample_id, bam_path.
 *   --reference    : reference FASTA.
 *
 * Output:
 *   results/merged.vcf.gz
 */

nextflow.enable.dsl = 2

params.sample_sheet = "samples.tsv"
params.reference    = "reference.fasta"
params.outdir       = "results"

workflow {
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.bam_path)) }

    reference_ch = Channel.value(file(params.reference))

    CALL_VARIANTS(samples_ch, reference_ch)
    MERGE_VCFS(CALL_VARIANTS.out.collect())
}

process CALL_VARIANTS {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'
    time '2.h'

    input:
    tuple val(sample_id), path(bam)
    path reference

    output:
    path "${sample_id}.vcf.gz"

    script:
    """
    bcftools mpileup -f ${reference} -Ou ${bam} \
        | bcftools call --threads ${task.cpus} -mv -Oz -o ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """
}

process MERGE_VCFS {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 4
    memory '8 GB'
    time '1.h'

    input:
    path vcfs

    output:
    path "merged.vcf.gz"

    script:
    """
    bcftools merge --threads ${task.cpus} -Oz -o merged.vcf.gz ${vcfs.join(' ')}
    bcftools index merged.vcf.gz
    """
}
'''


NEXTFLOW_CONFIG: str = '''\
// nextflow.config - C10 Week 12 Exercise 2

manifest {
    name        = 'c10-week12-exercise2'
    description = 'Minimal Nextflow variant-calling DAG for the C10 Week 12 capstone exercises.'
    version     = '1.0.0'
    nextflowVersion = '>=22.10.0'
    author      = 'C10 Crunching Bioinformatics'
}

process {
    container = 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    errorStrategy = 'retry'
    maxRetries = 1
}

profiles {

    standard {
        process.executor = 'local'
    }

    conda {
        conda.enabled = true
        process.conda = "bioconda::bcftools=1.20 bioconda::samtools=1.20"
    }

    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        singularity.cacheDir = "${HOME}/.singularity/cache"
    }

    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    slurm {
        process.executor = 'slurm'
        process.queue = 'short'
    }
}

report {
    enabled = true
    file = 'reports/exec_report.html'
}

dag {
    enabled = true
    file = 'reports/dag.svg'
}
'''


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def check_tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def nextflow_version() -> str:
    if not check_tool_available("nextflow"):
        return "unavailable"
    try:
        proc = subprocess.run(
            ["nextflow", "-version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        # Nextflow prints to stdout or stderr depending on the version.
        return (proc.stdout + proc.stderr).strip().splitlines()[0]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        return f"unavailable (error: {exc})"


def write_workflow_files(out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    main_nf = out_dir / "main.nf"
    main_nf.write_text(MAIN_NF_TEMPLATE)
    config_nf = out_dir / "nextflow.config"
    config_nf.write_text(NEXTFLOW_CONFIG)
    return main_nf, config_nf


def count_processes(main_nf: Path) -> int:
    return sum(1 for line in main_nf.read_text().splitlines() if line.startswith("process "))


def count_channels(main_nf: Path) -> int:
    """Approximate the number of distinct Channel.* references."""
    text = main_nf.read_text()
    return text.count("Channel.")


def render_dag(out_dir: Path, sample_sheet: Path, reference: Path) -> Path | None:
    """Render the Nextflow DAG using `-with-dag` (requires Nextflow on PATH)."""
    dag_svg = out_dir / "dag.svg"
    cmd = [
        "nextflow",
        "run",
        str(out_dir / "main.nf"),
        "-with-dag",
        str(dag_svg),
        "-preview",
        "--sample_sheet",
        str(sample_sheet),
        "--reference",
        str(reference),
        "--outdir",
        str(out_dir / "results"),
    ]
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            cwd=str(out_dir),
        )
        if dag_svg.exists() and dag_svg.stat().st_size > 0:
            return dag_svg
        # If `-preview` is unsupported, retry without it.
        cmd_no_preview = [c for c in cmd if c != "-preview"]
        cmd_no_preview.extend(["-stub"])
        proc = subprocess.run(
            cmd_no_preview,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            cwd=str(out_dir),
        )
        if dag_svg.exists() and dag_svg.stat().st_size > 0:
            return dag_svg
        return None
    except subprocess.TimeoutExpired:
        return None


def write_run_info(
    out_dir: Path,
    main_nf: Path,
    config_nf: Path,
    dag_svg: Path | None,
    nextflow_ver: str,
    process_count: int,
    channel_count: int,
    skip_reason: str | None,
    sample_sheet: Path,
    reference: Path,
) -> Path:
    info: dict[str, Any] = {
        "project": "C10-Week12-Exercise-02",
        "version": "v1.0",
        "run_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "workflow_manager": "nextflow",
        "workflow_manager_version": nextflow_ver,
        "workflow_file": str(main_nf.name),
        "workflow_file_sha256": sha256_of(main_nf),
        "config_file": str(config_nf.name),
        "config_file_sha256": sha256_of(config_nf),
        "process_count": process_count,
        "channel_count_approx": channel_count,
        "inputs": {
            "sample_sheet": str(sample_sheet),
            "reference": str(reference),
        },
        "dag_svg": str(dag_svg.name) if dag_svg is not None else None,
        "dag_svg_sha256": sha256_of(dag_svg) if dag_svg is not None else None,
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
            "Generate a minimal Nextflow main.nf + nextflow.config for the "
            "bcftools variant-calling step, render the DAG, and emit run-info."
        )
    )
    p.add_argument("--sample-sheet", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    main_nf, config_nf = write_workflow_files(out_dir)
    process_count = count_processes(main_nf)
    channel_count = count_channels(main_nf)
    print(f"[ex02] Wrote {main_nf}.")
    print(f"[ex02] Wrote {config_nf}.")
    print(f"[ex02] Process count:  {process_count}")
    print(f"[ex02] Channel count:  {channel_count}")

    skip_reason: str | None = None
    dag_svg: Path | None = None

    if not check_tool_available("nextflow"):
        skip_reason = "nextflow CLI not on PATH"
        print(f"[ex02] Skipping DAG render: {skip_reason}.")
    elif not args.sample_sheet.exists():
        skip_reason = f"sample sheet not present at {args.sample_sheet}"
        print(f"[ex02] Skipping DAG render: {skip_reason}.")
    elif not args.reference.exists():
        skip_reason = f"reference not present at {args.reference}"
        print(f"[ex02] Skipping DAG render: {skip_reason}.")
    else:
        dag_svg = render_dag(out_dir, args.sample_sheet, args.reference)
        if dag_svg is None:
            skip_reason = "Nextflow DAG render did not produce a non-empty SVG"
            print(f"[ex02] {skip_reason}")
        else:
            print(f"[ex02] Wrote {dag_svg} (size {dag_svg.stat().st_size} bytes).")

    nf_ver = nextflow_version()
    info_path = write_run_info(
        out_dir=out_dir,
        main_nf=main_nf,
        config_nf=config_nf,
        dag_svg=dag_svg,
        nextflow_ver=nf_ver,
        process_count=process_count,
        channel_count=channel_count,
        skip_reason=skip_reason,
        sample_sheet=args.sample_sheet,
        reference=args.reference,
    )
    print(f"[ex02] Wrote {info_path}.")

    print("[ex02] Summary:")
    print(f"[ex02]   nextflow version: {nf_ver}")
    print(f"[ex02]   processes:        {process_count}")
    print(f"[ex02]   channels (approx): {channel_count}")
    print(f"[ex02]   dag.svg:          {'OK' if dag_svg else 'skipped'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
