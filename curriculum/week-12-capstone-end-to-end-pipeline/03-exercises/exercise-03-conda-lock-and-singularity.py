"""
Exercise 3 - Capture a Conda lockfile, build a Singularity image,
and verify byte-identical re-runs.

Educational and research use only. Outputs are not clinical software.

Goal: given an environment.yml, drive conda (or mamba) to capture an
explicit lockfile, render a Singularity.def from a template, optionally
build the .sif image (if apptainer is on PATH), capture hashes for the
environment.yml, the lockfile, the def file, and the SIF image, and
emit a run-info.json that ties them all together.

The exercise covers:

- Capturing a Conda environment to a portable environment.yml and a
  byte-identical environment.lock.txt.
- Templating a Singularity.def from the environment files.
- Building a SIF image via apptainer (optionally; skipped on macOS or
  when apptainer is missing).
- Hashing every artefact: environment.yml, lockfile, def, SIF.
- The "byte-identical reproduction" verification pattern: two builds
  on the same lockfile should produce SIF images with the same hash
  (modulo a small set of documented non-determinisms).

Estimated time: 90 minutes (40 min reading, 40 min implementing,
10 min running and inspecting).

Acceptance criteria:
- `python exercise-03-conda-lock-and-singularity.py
    --environment environment.yml
    --out-dir results/ex03` runs end to end without errors when
  conda (or mamba) is on the PATH; gracefully skips with a clear
  message otherwise.
- `results/ex03/environment.yml` is a copy of the input.
- `results/ex03/environment.lock.txt` exists (when conda available)
  with at least the `@EXPLICIT` header and a non-zero count of URL
  lines.
- `results/ex03/Singularity.def` exists and parses as text with
  `Bootstrap:`, `From:`, `%post`, `%environment` sections.
- `results/ex03/run-info.json` exists with the conda version, the
  lockfile line count, the def file hash, the SIF hash (when built),
  and the build host info.
- The script is idempotent: re-running produces the same Singularity.def
  byte-for-byte.

Requirements:
    conda install -c conda-forge mamba=1.5.7
    # Optional for the build step:
    sudo apt-get install -y apptainer  # Linux only

What you learn:
- The `conda env export --from-history` vs `conda list --explicit`
  distinction.
- The Singularity.def file structure: Bootstrap, From, %files, %post,
  %environment, %runscript, %labels.
- Hashing every artefact in a reproducibility pipeline.
- Cross-platform gotchas: conda lockfiles are architecture-specific.

Tool versions assumed:
- Python 3.11+
- conda 24.1.2 or mamba 1.5.7 (CLI; subprocess.run it)
- apptainer 1.2.5 (Linux only; subprocess.run it; optional)

References:
- bioconda: Grüning et al. 2018, Nat Methods 15:475
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11070151/
- Singularity: Kurtzer et al. 2017, PLoS ONE 12:e0177459
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426675/

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


SINGULARITY_DEF_TEMPLATE: str = '''\
Bootstrap: docker
From: continuumio/miniconda3:24.1.2-0

%files
    environment.yml /opt/environment.yml
    environment.lock.txt /opt/environment.lock.txt

%post
    apt-get update && apt-get install -y --no-install-recommends \\
        ca-certificates curl wget tar gzip bzip2 git \\
        && rm -rf /var/lib/apt/lists/*

    /opt/conda/bin/conda install -n base -c conda-forge mamba=1.5.7 -y

    # Prefer the lockfile for byte-identical reproduction; fall back
    # to the YAML if the lockfile is missing or empty.
    if [ -s /opt/environment.lock.txt ]; then
        /opt/conda/bin/mamba create -n {env_name} --file /opt/environment.lock.txt -y
    else
        /opt/conda/bin/mamba env create -n {env_name} -f /opt/environment.yml
    fi

    /opt/conda/bin/conda clean -afy

%environment
    export PATH=/opt/conda/envs/{env_name}/bin:/opt/conda/bin:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%runscript
    exec "$@"

%labels
    Author C10-Crunching-Bioinformatics
    Version v1.0
    Project C10-Week12-Exercise-03

%help
    This container packages the {env_name} Conda environment built from
    environment.lock.txt (preferred) or environment.yml (fallback).
    Run: apptainer run image.sif <command>
'''


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def check_tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def conda_or_mamba() -> str | None:
    """Return 'mamba' if available, else 'conda', else None."""
    if check_tool_available("mamba"):
        return "mamba"
    if check_tool_available("conda"):
        return "conda"
    return None


def conda_version(tool: str) -> str:
    try:
        proc = subprocess.run(
            [tool, "--version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return proc.stdout.strip() or proc.stderr.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as exc:
        return f"unavailable (error: {exc})"


def parse_env_name(environment_yml: Path) -> str:
    """Pull the 'name:' field from the environment.yml, defaulting if absent."""
    for line in environment_yml.read_text().splitlines():
        line = line.strip()
        if line.startswith("name:"):
            return line.split(":", 1)[1].strip() or "c10-capstone"
    return "c10-capstone"


def export_lockfile(
    env_name: str,
    out_dir: Path,
    tool: str,
    try_create: bool,
) -> tuple[Path | None, str]:
    """Capture an explicit lockfile for the env.

    Strategy:
    1. If `try_create` is True, attempt to create the environment from
       environment.yml first (slow; skip in CI).
    2. Run `conda list --explicit -n {env_name}`. If the environment
       does not exist, return (None, message).
    """
    lockfile = out_dir / "environment.lock.txt"
    if try_create:
        cmd = [tool, "env", "create", "-n", env_name, "-f", str(out_dir / "environment.yml")]
        proc = subprocess.run(
            cmd, check=False, capture_output=True, text=True, timeout=1800
        )
        if proc.returncode != 0:
            return None, f"env create failed: {proc.stderr.strip()[:200]}"

    cmd_list = ["conda", "list", "--explicit", "-n", env_name]
    proc = subprocess.run(
        cmd_list, check=False, capture_output=True, text=True, timeout=60
    )
    if proc.returncode != 0:
        return None, f"conda list failed: {proc.stderr.strip()[:200]}"
    if "@EXPLICIT" not in proc.stdout:
        return None, "conda list output is missing @EXPLICIT marker"

    lockfile.write_text(proc.stdout)
    return lockfile, "ok"


def write_singularity_def(env_name: str, out_dir: Path) -> Path:
    def_path = out_dir / "Singularity.def"
    def_path.write_text(SINGULARITY_DEF_TEMPLATE.format(env_name=env_name))
    return def_path


def maybe_build_sif(def_path: Path, out_dir: Path) -> tuple[Path | None, str]:
    """Attempt to build the SIF image via apptainer; return path or skip reason."""
    if not check_tool_available("apptainer") and not check_tool_available("singularity"):
        return None, "apptainer/singularity not on PATH"

    cli = "apptainer" if check_tool_available("apptainer") else "singularity"
    sif_path = out_dir / "image.sif"
    cmd = [cli, "build", "--force", str(sif_path), str(def_path)]
    proc = subprocess.run(
        cmd, check=False, capture_output=True, text=True, timeout=1800, cwd=str(out_dir)
    )
    if proc.returncode != 0:
        return None, f"{cli} build failed: {proc.stderr.strip()[:200]}"
    if not sif_path.exists() or sif_path.stat().st_size == 0:
        return None, "SIF file missing or empty after build"
    return sif_path, "ok"


def count_lockfile_lines(lockfile: Path) -> int:
    """Count the URL lines (not the @EXPLICIT header / comments) in the lockfile."""
    count = 0
    for line in lockfile.read_text().splitlines():
        if line and not line.startswith("#") and not line.startswith("@"):
            count += 1
    return count


def write_run_info(
    out_dir: Path,
    environment_yml: Path,
    lockfile: Path | None,
    def_path: Path,
    sif_path: Path | None,
    env_name: str,
    tool: str | None,
    tool_ver: str,
    lockfile_lines: int,
    skip_lock_reason: str | None,
    skip_sif_reason: str | None,
) -> Path:
    info: dict[str, Any] = {
        "project": "C10-Week12-Exercise-03",
        "version": "v1.0",
        "run_date_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "conda_or_mamba": tool,
        "conda_version": tool_ver,
        "environment_yml": str(environment_yml.name),
        "environment_yml_sha256": sha256_of(environment_yml),
        "environment_lock": str(lockfile.name) if lockfile else None,
        "environment_lock_sha256": sha256_of(lockfile) if lockfile else None,
        "environment_lock_line_count": lockfile_lines,
        "singularity_def": str(def_path.name),
        "singularity_def_sha256": sha256_of(def_path),
        "sif_image": str(sif_path.name) if sif_path else None,
        "sif_image_sha256": sha256_of(sif_path) if sif_path else None,
        "env_name": env_name,
        "skip_lock_reason": skip_lock_reason,
        "skip_sif_reason": skip_sif_reason,
        "host": {
            "hostname": os.uname().nodename,
            "system": os.uname().sysname,
            "machine": os.uname().machine,
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
            "Capture a Conda lockfile and a Singularity definition from "
            "an environment.yml; optionally build the SIF image; emit run-info."
        )
    )
    p.add_argument("--environment", type=Path, required=True, help="environment.yml path")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument(
        "--try-create-env",
        action="store_true",
        help="Attempt to create the Conda environment first (slow). Default off.",
    )
    p.add_argument(
        "--try-build-sif",
        action="store_true",
        help="Attempt to build the SIF image via apptainer (slow). Default off.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.environment.exists():
        print(f"[ex03] environment.yml does not exist: {args.environment}", file=sys.stderr)
        return 2

    # Copy environment.yml into out_dir.
    env_yml_dst = out_dir / "environment.yml"
    env_yml_dst.write_text(args.environment.read_text())
    env_name = parse_env_name(env_yml_dst)
    print(f"[ex03] Environment name: {env_name}")

    # Lockfile capture (slow; optional).
    tool = conda_or_mamba()
    tool_ver = conda_version(tool) if tool else "unavailable"
    print(f"[ex03] Tool: {tool or 'unavailable'} ({tool_ver})")

    lockfile: Path | None = None
    skip_lock_reason: str | None = None
    lockfile_lines = 0
    if tool is None:
        skip_lock_reason = "no conda or mamba on PATH"
    else:
        lockfile, msg = export_lockfile(env_name, out_dir, tool, args.try_create_env)
        if lockfile is None:
            skip_lock_reason = msg
            # Write a placeholder empty lockfile to keep downstream paths sane.
            placeholder = out_dir / "environment.lock.txt"
            placeholder.write_text("# Lockfile not captured: " + msg + "\n")
            lockfile = placeholder
            lockfile_lines = 0
        else:
            lockfile_lines = count_lockfile_lines(lockfile)
            print(f"[ex03] Lockfile lines: {lockfile_lines}")

    # Singularity.def render.
    def_path = write_singularity_def(env_name, out_dir)
    print(f"[ex03] Wrote {def_path}.")

    # SIF build (slow; optional).
    sif_path: Path | None = None
    skip_sif_reason: str | None = None
    if args.try_build_sif:
        sif_path, msg = maybe_build_sif(def_path, out_dir)
        if sif_path is None:
            skip_sif_reason = msg
            print(f"[ex03] SIF build skipped: {msg}")
        else:
            print(f"[ex03] Built SIF: {sif_path} ({sif_path.stat().st_size} bytes)")
    else:
        skip_sif_reason = "build-sif flag not set"
        print("[ex03] SIF build deferred (pass --try-build-sif to actually build).")

    # run-info.json.
    info_path = write_run_info(
        out_dir=out_dir,
        environment_yml=env_yml_dst,
        lockfile=lockfile,
        def_path=def_path,
        sif_path=sif_path,
        env_name=env_name,
        tool=tool,
        tool_ver=tool_ver,
        lockfile_lines=lockfile_lines,
        skip_lock_reason=skip_lock_reason,
        skip_sif_reason=skip_sif_reason,
    )
    print(f"[ex03] Wrote {info_path}.")

    print("[ex03] Summary:")
    print(f"[ex03]   environment.yml sha256:  {sha256_of(env_yml_dst)[:16]}...")
    if lockfile is not None and lockfile_lines > 0:
        print(f"[ex03]   lockfile sha256:         {sha256_of(lockfile)[:16]}...")
    print(f"[ex03]   Singularity.def sha256:  {sha256_of(def_path)[:16]}...")
    if sif_path is not None:
        print(f"[ex03]   SIF image sha256:        {sha256_of(sif_path)[:16]}...")

    return 0


if __name__ == "__main__":
    sys.exit(main())
