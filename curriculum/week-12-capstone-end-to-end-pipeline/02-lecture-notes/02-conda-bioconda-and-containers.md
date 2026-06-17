# Lecture 2 — Conda, bioconda, and Containers

> **Educational and research use only.** Conda, bioconda, Singularity / Apptainer, and Docker are all free open-source tools. The container images you build for the capstone are educational artefacts; outputs are not clinical software.

Lecture 1 covered the workflow-manager layer. Lecture 2 covers the two layers below it: the **package manager** (Conda + bioconda) that resolves and installs tool dependencies, and the **container** (Singularity / Apptainer or Docker) that wraps the whole environment for bootstrapping reproducibility.

The three-layer model:

```
+-----------------------------------------------+
| Workflow manager (Snakemake / Nextflow)        |  ← rules / processes; DAG; --use-conda --use-singularity
+-----------------------------------------------+
| Package manager (Conda + bioconda)             |  ← environment.yml; environment.lock.txt
+-----------------------------------------------+
| Container runtime (Singularity / Apptainer)    |  ← .sif image; Singularity.def
+-----------------------------------------------+
| OS (Linux kernel, glibc, etc.)                 |  ← the laptop or the cluster node
+-----------------------------------------------+
```

Each layer pins the layer above it. The OS pins glibc / kernel; the container pins the OS subset; the package manager pins the tool versions; the workflow manager pins the rule graph and the input file paths. The capstone targets reproducibility at all four layers, but the easiest engineering wins are at the package manager and the container.

## 2.1 — What Conda does

Conda is a cross-platform package and environment manager. It is platform-agnostic (Linux x86_64, Linux aarch64, macOS x86_64, macOS arm64, Windows x86_64). It installs binary packages plus their dependencies. It manages multiple isolated environments per machine.

The four commands you need:

```bash
# Create an environment from a YAML file:
conda env create -f environment.yml -n c10-capstone

# Activate the environment:
conda activate c10-capstone

# List the installed packages:
conda list

# Export the resolved environment:
conda env export --from-history > environment.yml
conda list --explicit > environment.lock.txt
```

The two export commands produce two different files. **Why both?**

- `conda env export --from-history` writes a human-readable YAML listing the packages you explicitly asked for, with their version constraints. This is the *intention*: "I want Salmon 1.10 and DESeq2 1.42." It is portable across architectures (a macOS arm64 user can create an environment from a Linux x86_64 user's `environment.yml`).
- `conda list --explicit` writes a flat file of resolved package URLs with SHA-256 hashes. This is the *byte-identical reproduction*: every URL and every hash is pinned. The lockfile is architecture-specific (a macOS arm64 lockfile cannot install on Linux x86_64).

The capstone commits **both**. The `environment.yml` is the readable specification; the `environment.lock.txt` is the byte-identical reproduction.

## 2.2 — Channels and channel order

A Conda channel is a hosted repository of packages. Conda resolves a package by searching the channels in priority order, top to bottom, picking the first version-matching package.

For bioinformatics, the canonical channel order is:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

(Adding pushes to the top; the resulting priority order is **conda-forge > bioconda > defaults**.)

**`conda-forge`** is the community channel for general-purpose packages (Python, numpy, pandas, scipy, R, etc.). It is the largest Conda channel and has 25k+ packages as of 2026.

**`bioconda`** (Grüning et al. 2018, *Nature Methods* 15:475; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11070151/>) is the bioinformatics channel. 8000+ packages. Maintained by a continuous-integration build farm. Every package's recipe is in the bioconda GitHub repo and every package's build is automated.

**`defaults`** is the Anaconda Inc. commercial channel. Free for academic use; recent licence changes restrict commercial use at large enterprises (over 200 employees). The capstone uses bioconda and conda-forge as the primary sources and treats `defaults` as a fallback.

**`channel_priority: strict`** is the canonical setting. Strict priority means Conda will not mix packages from a lower-priority channel to satisfy a dependency declared in a higher-priority channel. Without strict priority, you get hard-to-debug installs where (say) numpy comes from conda-forge but its scipy dependency comes from defaults, and the two channels' builds are subtly incompatible.

## 2.3 — A minimal `environment.yml`

For an RNA-seq capstone:

```yaml
name: c10-capstone-rnaseq
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - snakemake=7.32.4
  - salmon=1.10.2
  - samtools=1.20
  - fastp=0.23.4
  - bioconductor-tximport=1.30.0
  - bioconductor-deseq2=1.42.0
  - r-base=4.3.2
  - r-ggplot2=3.5.0
  - pandas=2.2.1
  - matplotlib=3.8.3
```

Read it top to bottom:

- **`name:`** — the environment name. `conda env create` uses this if you do not pass `-n`.
- **`channels:`** — the channel order, mirroring what `conda config` would set.
- **`dependencies:`** — the packages. Each line is `package=version`. Without `=version`, Conda picks the latest; **always pin** for the capstone.

Notice that R packages live in bioconda under the `bioconductor-` prefix (`bioconductor-tximport`, `bioconductor-deseq2`) and that R itself comes from conda-forge (`r-base`). The naming convention is `r-<cran-package>` for CRAN packages and `bioconductor-<bioconductor-package>` for Bioconductor packages.

## 2.4 — The lockfile

Run:

```bash
conda env create -f environment.yml -n c10-capstone-rnaseq
conda activate c10-capstone-rnaseq
conda list --explicit > environment.lock.txt
```

The `environment.lock.txt` is a flat list of URLs with SHA-256 hashes:

```
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/python-3.11.8-hab00c5b_0_cpython.conda#abc...
https://conda.anaconda.org/bioconda/linux-64/salmon-1.10.2-h6dccd9a_1.tar.bz2#def...
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-13.2.0-h807b86a_5.conda#ghi...
...
```

Re-creating the environment from the lockfile:

```bash
conda create -n c10-capstone-rnaseq --file environment.lock.txt
```

This is byte-identical to the original (modulo a few transient install logs). It is the closest you get to "the exact environment" on the same architecture.

**Caveat:** the lockfile is **architecture-specific**. A lockfile generated on Linux x86_64 does not install on macOS arm64. The convention is to commit a `environment.lock.linux-64.txt` and (if relevant) a `environment.lock.osx-arm64.txt`. The capstone targets Linux x86_64 because that is what Snakemake / Nextflow clusters and Singularity assume.

## 2.5 — mamba and micromamba

Conda's dependency solver is a Python implementation of pseudo-boolean SAT. For small environments (10-30 packages) it is fast; for large environments (100+ packages with deep dependency trees), it can take 5-30 minutes. This is a known pain point.

`mamba` is a drop-in replacement written in C++. It uses the libsolv library (the same solver as Fedora's DNF and OpenSUSE's zypper). It is 10-100x faster. Install it once into the base environment:

```bash
conda install -n base -c conda-forge mamba
```

After that, replace every `conda` command with `mamba` (or use Snakemake's `--conda-frontend mamba` flag) for the speed-up. `mamba env create -f environment.yml` is the canonical fast install.

`micromamba` is a smaller static binary of mamba (no Python required at the runtime). It is ideal for Dockerfiles and Singularity definitions where you do not want to bake a full Conda installation into the container.

The 2026 convention: use mamba (or micromamba) as the solver. Snakemake's `--conda-frontend mamba` honours this; Nextflow's `conda.useMamba = true` in `nextflow.config` honours this.

## 2.6 — What a container is

A container is a Linux process tree running in a namespace-isolated subset of the host. The container sees its own root filesystem (read-only or copy-on-write), its own process IDs, its own network stack (or shares the host's), and its own user IDs. From the container's perspective it is a self-contained Linux system; from the host's perspective the container's processes appear in `ps` with the host's PIDs.

Containers are not virtual machines. A VM emulates the hardware (CPU, RAM, disk, network) and runs a guest kernel; a container shares the host kernel and isolates only the process tree. This is why containers are lighter than VMs (no boot, no kernel duplication, no hypervisor overhead) and why they can only run Linux processes on a Linux host. The "Docker for Mac" experience involves a hidden Linux VM that hosts the containers.

For bioinformatics reproducibility, the relevant property of a container is: **the same container image, run on different hosts, produces the same output**. The container freezes the OS, the libc, the Python version, the tool versions, every shared library, every random configuration file in `/etc`. Two hosts with different OS versions still see the same container.

## 2.7 — Singularity / Apptainer vs Docker

| Property | Singularity / Apptainer | Docker |
|----------|--------------------------|--------|
| Daemon | None (runs as plain process) | Required (dockerd) |
| Root | Optional; rootless by default | Required (or rootless-Docker with caveats) |
| Image format | SIF (Singularity Image Format), single file | Layered images stored in `/var/lib/docker` |
| Image distribution | Single `.sif` file; `apptainer pull`, `cp`, `rsync` | Image registry (Docker Hub, quay.io); `docker pull` |
| HPC compatibility | High (no daemon, rootless) | Low (daemon, root usually required) |
| Cloud compatibility | Via Docker conversion | High |
| Build host requirement | Can build from Docker | Cannot read SIF natively |

For **HPC / academic-cluster** capstone deposits, **Singularity / Apptainer is the canonical choice**. The single `.sif` file deposits cleanly, the rootless execution model works in shared-user environments, and the HPC sysadmins do not need to run a Docker daemon.

For **cloud-only** deposits (AWS Batch, Google Cloud Run, Kubernetes), **Docker is the canonical choice**. The cloud orchestrators are Docker-native and pulling from a registry is the assumed workflow.

The capstone defaults to Singularity / Apptainer because the `.sif` file is a single deposit-friendly artefact that Zenodo can accept directly.

## 2.8 — A minimal `Singularity.def`

```
Bootstrap: docker
From: continuumio/miniconda3:24.1.2-0

%files
    environment.yml /opt/environment.yml
    environment.lock.txt /opt/environment.lock.txt

%post
    # Update the base image
    apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget tar gzip bzip2 git \
        && rm -rf /var/lib/apt/lists/*

    # Install mamba into base
    /opt/conda/bin/conda install -n base -c conda-forge mamba -y

    # Create the capstone environment from the lockfile
    /opt/conda/bin/mamba create -n c10-capstone --file /opt/environment.lock.txt -y

    # Clean up
    /opt/conda/bin/conda clean -afy

%environment
    export PATH=/opt/conda/envs/c10-capstone/bin:/opt/conda/bin:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%runscript
    exec "$@"

%labels
    Author your_name
    Version v1.0
    Project C10-Capstone

%help
    This container packages the C10 capstone environment.
    Run: apptainer run c10-capstone.sif snakemake --cores 4 --use-conda
```

Read it section by section:

- **`Bootstrap: docker` + `From: continuumio/miniconda3:24.1.2-0`** — start from the official miniconda Docker image (a thin Linux + miniconda installation). The version (`24.1.2-0`) is pinned.
- **`%files`** — copy the `environment.yml` and lockfile into the container at build time.
- **`%post`** — commands to run inside the container at build time. Install system dependencies, install mamba, create the Conda environment from the lockfile.
- **`%environment`** — environment variables set at every `apptainer exec` / `apptainer run`. Put the capstone Conda env's `bin/` on the PATH.
- **`%runscript`** — the default entry point. `apptainer run c10-capstone.sif foo bar` runs `foo bar` inside the container.
- **`%labels`** — metadata. Surfaced by `apptainer inspect c10-capstone.sif`.
- **`%help`** — help text. Surfaced by `apptainer run-help c10-capstone.sif`.

Build the image:

```bash
apptainer build c10-capstone.sif Singularity.def
```

Run it:

```bash
apptainer run c10-capstone.sif snakemake --cores 4
```

Or interactively:

```bash
apptainer shell c10-capstone.sif
```

## 2.9 — A minimal `Dockerfile`

The same recipe as a Dockerfile:

```dockerfile
FROM continuumio/miniconda3:24.1.2-0

COPY environment.yml /opt/environment.yml
COPY environment.lock.txt /opt/environment.lock.txt

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget tar gzip bzip2 git \
    && rm -rf /var/lib/apt/lists/*

RUN /opt/conda/bin/conda install -n base -c conda-forge mamba -y \
    && /opt/conda/bin/mamba create -n c10-capstone --file /opt/environment.lock.txt -y \
    && /opt/conda/bin/conda clean -afy

ENV PATH=/opt/conda/envs/c10-capstone/bin:/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /work
ENTRYPOINT []
CMD ["snakemake", "--cores", "4"]
```

Build:

```bash
docker build -t c10-capstone:v1.0 .
```

Run:

```bash
docker run --rm -v "$PWD:/work" c10-capstone:v1.0 snakemake --cores 4
```

Convert to Apptainer / Singularity:

```bash
apptainer pull docker-daemon://c10-capstone:v1.0
```

(The `apptainer pull` command works on locally-built Docker images via the Docker daemon, or on Docker Hub / quay.io images via `docker://` URLs.)

## 2.10 — Snakemake's container integration

Snakemake natively wires Conda environments per rule:

```python
rule align_reads:
    input: ...
    output: ...
    conda: "envs/bwa.yaml"
    shell: "bwa mem ..."
```

Run with `snakemake --use-conda`. Snakemake creates a Conda environment per rule on first invocation (cached in `.snakemake/conda/`) and activates it before each shell command.

Snakemake also natively wires Singularity containers per rule:

```python
rule align_reads:
    input: ...
    output: ...
    singularity: "docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    shell: "bwa mem ..."
```

Run with `snakemake --use-singularity`. Snakemake pulls the container on first invocation, mounts the working directory, and runs the shell command inside.

You can combine them: `--use-conda --use-singularity` runs each rule inside a per-rule Conda environment inside a per-rule container. The container provides the base OS; the Conda environment inside it pins the tools.

For the capstone, the typical configuration is:

- **One Conda environment** for the whole pipeline (the easiest to manage), specified in `environment.yml` at the workflow root. The `Snakefile`'s rules do not need per-rule `conda:` directives.
- **One Singularity container** for the whole pipeline, built from the same `environment.yml` via `Singularity.def`. The workflow runs with `--use-singularity` and the container's `PATH` provides the tools.

This is the "one environment, one container" model. Simpler to reason about than the "one environment per rule" model. The capstone evaluation accepts either.

## 2.11 — Nextflow's container integration

Nextflow declares the container in `nextflow.config`:

```groovy
process {
    container = 'c10-capstone:v1.0'
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = "${HOME}/.singularity/cache"
}

// Or for Docker:
// docker {
//     enabled = true
// }
```

Each process inherits the container by default. Override per-process if needed:

```groovy
process CUSTOM_RULE {
    container = 'quay.io/biocontainers/somespecialtool:1.2.3--pyh5e36f6f_0'
    script: "somespecialtool --input ..."
}
```

Run with `nextflow run main.nf -profile singularity` (assumes `profiles.singularity` is defined in `nextflow.config`).

## 2.12 — BioContainers: the Conda-Docker bridge

Every bioconda package is automatically built into a Docker image hosted at quay.io under `quay.io/biocontainers/<package>:<version>--<build-hash>`. This is the BioContainers system (<https://biocontainers.pro/>). The build is automated; the image is reproducible.

You can grab a single-tool container per rule:

```python
rule align_reads:
    container: "docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    shell: "bwa mem ..."
```

This is the "per-rule container" model. It is more granular than the "one container per pipeline" model and is the convention in nf-core. For a capstone, the one-container model is simpler; for a production pipeline, the per-rule model is cleaner.

## 2.13 — Hashing and verification

The capstone's reproducibility claim is verifiable. After building the container, capture its hash:

```bash
sha256sum c10-capstone.sif
```

Embed the hash in `run-info.json`:

```json
{
  "singularity_image": "c10-capstone.sif",
  "singularity_image_sha256": "abc123...",
  "build_date": "2026-05-14",
  "git_commit": "def456",
  "snakemake_version": "7.32.4"
}
```

After running the pipeline, capture each output file's hash:

```bash
sha256sum results/*.tsv results/*.vcf.gz > results/hashes.txt
```

Re-running on a clean machine should produce the same hashes for deterministic rules. The verification step:

```bash
sha256sum -c results/hashes.txt
```

This is the byte-identical-reproduction check. The Week 12 challenges include one that walks through this end to end.

## 2.14 — Common pitfalls

**Pinning the lockfile to the wrong architecture.** A Linux x86_64 lockfile does not install on macOS arm64. Commit one lockfile per architecture you support, or pick Linux x86_64 as the canonical architecture and document it.

**Trusting `defaults` for bioconda packages.** The Anaconda commercial channel does not have most bioinformatics tools. Always set `channel_priority: strict` and put bioconda above defaults.

**Forgetting to clean Conda's cache before building the container.** `conda clean -afy` shrinks the image by 200-500 MB. Without it, the SIF file ships every cached package tarball.

**Building the container on macOS for Linux runtime.** Apptainer / Singularity images must be built on Linux (or in a Linux VM). Docker Desktop on macOS does this transparently. Snakemake / Nextflow's `--use-singularity` mode on macOS requires a Linux VM (via Lima, Multipass, or Docker Desktop).

**Pulling a Docker image with `:latest`.** `:latest` is not a version pin; it is "whatever the registry currently calls latest." Always pin to a specific tag (`:v1.0.0`) or, better, a specific digest (`@sha256:abc...`).

**Mounting the wrong working directory.** Singularity's default bind list does not include `/scratch`, `/projects`, or other site-specific paths. Use `-B /scratch:/scratch` (or set `SINGULARITY_BINDPATH`) explicitly. Snakemake's `--singularity-args "-B /scratch"` is the canonical way to add binds.

**Running the container as root.** Singularity / Apptainer's design point is rootless execution. If you run as root (or with `--fakeroot`) inside the container, the output files are owned by root and the host user cannot read them without a `chown`.

## 2.15 — Recap

- Conda + bioconda + conda-forge is the standard package-manager stack for bioinformatics. Pin every package to a specific version in `environment.yml`. Commit both `environment.yml` (human-readable, portable) and `environment.lock.txt` (byte-identical, architecture-specific).
- Mamba is the faster drop-in solver. Use it everywhere.
- Singularity / Apptainer is the HPC-friendly container runtime. The SIF file is a single-file deposit-friendly artefact. Build from a `Singularity.def` that bootstraps off the miniconda Docker base.
- Snakemake and Nextflow natively integrate Conda and Singularity. The recommended capstone configuration: one Conda environment, one Singularity image, the workflow manager wires the environment into every rule.
- Hash everything: the container, the inputs, the outputs. The byte-identical reproduction check is `sha256sum -c results/hashes.txt`.

Lecture 3 covers the four capstone project tracks (variant discovery, RNA-seq DE, MSA + phylogeny, long-read assembly), the `run-info.json` pattern, the Zenodo DOI deposit, and the wrap-up sidebar that closes C10.
