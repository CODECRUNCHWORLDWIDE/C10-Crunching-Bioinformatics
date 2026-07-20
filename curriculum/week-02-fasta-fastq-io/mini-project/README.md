# Mini-Project — QC Report on a 1000 Genomes Subset

> Produce a clean, reproducible quality-control (QC) report for a real public FASTQ subset from the 1000 Genomes Project. The deliverable is a one-page Markdown report, a clean `qc/` directory of reproducible outputs (FastQC HTML, your own quality plots, summary JSON), and a `run.sh` script that re-creates the whole thing in one command.

This is the first C10 mini-project that produces a **runnable analytical artifact**. By the end of it you will have parsed a real Illumina FASTQ, run FastQC, written your own per-base quality plot, applied a documented trimming policy, and produced a methods-section-quality report comparing pre- and post-trim. Every later week of C10 produces a report in roughly this shape; you are learning the format now.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-02/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 2 section)
└── week-02/
    ├── README.md                   one-page QC report (~600 words)
    ├── run.sh                      one-command reproduction script
    ├── env.yml                     conda environment file pinning all tool versions
    ├── data/                       (.gitignored — never commit raw FASTQ)
    │   └── .gitkeep
    ├── qc/
    │   ├── fastqc_raw/             FastQC report on raw reads (HTML + zip)
    │   ├── fastqc_trimmed/         FastQC report on trimmed reads
    │   ├── mean_quality_raw.png    your own per-base mean-Q plot
    │   ├── mean_quality_trimmed.png
    │   └── summary.json            machine-readable summary (your code)
    └── scripts/
        ├── inspect_fastq.py        from homework P1, adapted
        ├── trim.py                 from homework P4, adapted
        └── plot_quality.py         from exercise 2, adapted
```

By the end you will have a clean, reproducible Week 2 directory you can point a recruiter at — and the report inside is the kind of methods-section reading that opens conversations with working bioinformaticians.

---

## The dataset

We use the 1000 Genomes Project (phase 3) and the sample **NA12878** (CEPH/Utah, female, the most-sequenced human genome on Earth, the canonical reference sample for short-read pipeline development).

You will work with **chromosome 22 only** to keep run times reasonable. Chr22 is the smallest human autosome (~50 Mb) and a long-standing tradition for bioinformatics teaching examples.

Pick **exactly one** of the following accession options:

- **Option A — Direct FASTQ.** `ERR1019034` from ENA, paired-end Illumina HiSeq, 100 bp reads. Download just R1 for this project (~3 GB compressed, ~10 GB uncompressed — work with the .gz).
- **Option B — Subset.** If 3 GB is too much, generate a 1-million-read subset with `seqkit head -n 1000000` after download.
- **Option C — Conservative.** A 100k-read subset is plenty for the QC patterns to be visible and runs in seconds.

Pick the option that fits your time and disk budget. Document which one you used.

---

## Rules

- **You may** use Biopython 1.83, FastQC 0.12.1, seqkit 2.8, gzip, matplotlib, and the standard library.
- **You may** consult Biopython's tutorial, FastQC's docs, and the lecture notes.
- **You may NOT** copy a pre-written QC pipeline. The point is to *write* the pipeline. Looking at Snakemake / Nextflow QC workflows for inspiration is fine; copying them is not.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file.

---

## Acceptance criteria

- [ ] `week-02/README.md` is a one-page (≤ 700 word) report containing:
  - One-sentence description of dataset, sample, accession, and chromosome.
  - Methods section in C10 voice: every tool pinned ("FastQC 0.12.1", "Biopython 1.83", "seqkit 2.8"), every parameter explicit ("`--min-len 36 --min-mean-q 20 --window 4 --window-threshold 20`").
  - Quantitative findings: read count, mean read length, mean per-base quality at positions 1, 35, 75 (or whatever your read length is), pre- and post-trim.
  - At least one numbered observation from the FastQC report that you investigated — was it a real problem? A library-type artifact?
  - A reproducibility receipt block (see below).
- [ ] `week-02/run.sh` runs end-to-end on a fresh clone:
  - Activates the conda env from `env.yml`.
  - Downloads the FASTQ (or accepts it as the first argument if already on disk).
  - Runs `fastqc` on raw.
  - Runs your trim script.
  - Runs `fastqc` on trimmed.
  - Runs your `plot_quality.py` for both pre- and post-trim.
  - Writes `summary.json` with the key numbers.
- [ ] `week-02/env.yml` pins every tool to an exact version. Use `conda env export --no-builds` after manual install if you want a starting point.
- [ ] `week-02/qc/summary.json` is valid JSON with the shape:
  ```json
  {
    "sample": "NA12878",
    "accession": "ERR1019034",
    "chromosome": "22",
    "n_reads_raw": 100000,
    "n_reads_trimmed": 97834,
    "mean_len_raw": 100.0,
    "mean_len_trimmed": 91.4,
    "mean_per_base_q_raw_pos_1": 36.2,
    "mean_per_base_q_raw_pos_35": 35.1,
    "mean_per_base_q_raw_pos_75": 29.3,
    "encoding_detected": "Phred+33",
    "tools": {
      "fastqc": "0.12.1",
      "biopython": "1.83",
      "seqkit": "2.8.x"
    }
  }
  ```
- [ ] `week-02/qc/fastqc_raw/` and `week-02/qc/fastqc_trimmed/` contain the FastQC HTML reports.
- [ ] `week-02/qc/mean_quality_*.png` are produced by your own code (not screenshots of FastQC).
- [ ] `data/` is `.gitignored`. **You must not commit raw FASTQ files.** They are huge and they belong on the original archive.
- [ ] The repo passes a fresh-clone test: `git clone`, `cd week-02`, `bash run.sh ./data/your.fastq.gz` reproduces everything in `qc/`.

---

## Suggested order of operations

### Phase 1 — Environment setup (~30 min)

1. Create `week-02/env.yml` with:
   ```yaml
   name: c10-week-02
   channels:
     - conda-forge
     - bioconda
   dependencies:
     - python=3.11
     - biopython=1.83
     - matplotlib
     - fastqc=0.12.1
     - seqkit=2.8
     - pip
   ```
2. `conda env create -f week-02/env.yml`. Activate it. Confirm: `fastqc --version` → `FastQC v0.12.1`.
3. Commit: `Week 2 env.yml`.

### Phase 2 — Download a subset (~30 min)

1. Pick your accession (see above). Download from ENA's direct FASTQ URL. Example for ERR1019034 R1:
   ```bash
   curl -L 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1019034/ERR1019034_1.fastq.gz' \
       -o data/ERR1019034_1.fastq.gz
   ```
2. Make a subset if needed: `seqkit head -n 100000 data/ERR1019034_1.fastq.gz | gzip > data/sub.fastq.gz`.
3. Add `data/` to `.gitignore`. Add the file path to `run.sh` as an environment variable.
4. Commit: `Subset download script for ERR1019034`.

### Phase 3 — Adapt your scripts (~1 hour)

1. Copy your homework P1 (`inspect_fastq.py`), homework P4 (`trim.py`), and exercise 2 (`plot_quality.py`) into `week-02/scripts/`.
2. Make them accept a path as their first argument and write into the `qc/` directory.
3. Run each manually once to confirm.
4. Commit: `Mini-project scripts wired up`.

### Phase 4 — Wire up `run.sh` (~1 hour)

Your `run.sh` should be a straight shell script. Aim for ~30 lines. Pattern:

```bash
#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-data/sub.fastq.gz}"

mkdir -p qc/fastqc_raw qc/fastqc_trimmed

# 1) FastQC on raw
fastqc "${INPUT}" -o qc/fastqc_raw/

# 2) Your own quality plot on raw
python scripts/plot_quality.py "${INPUT}" qc/mean_quality_raw.png

# 3) Trim
python scripts/trim.py \
    --in "${INPUT}" \
    --out qc/trimmed.fastq.gz \
    --min-len 36 --min-mean-q 20 --window 4 --window-threshold 20

# 4) FastQC on trimmed
fastqc qc/trimmed.fastq.gz -o qc/fastqc_trimmed/

# 5) Your own quality plot on trimmed
python scripts/plot_quality.py qc/trimmed.fastq.gz qc/mean_quality_trimmed.png

# 6) Summary JSON
python scripts/summarize.py \
    --raw "${INPUT}" \
    --trimmed qc/trimmed.fastq.gz \
    --out qc/summary.json

echo "Done. See qc/summary.json and week-02/README.md."
```

Commit: `run.sh end-to-end`.

### Phase 5 — Write the report (~2 hours)

Open `week-02/README.md`. Structure:

```
# Week 2 QC report — ERR1019034 chr22 subset

## Sample and accession
- Sample: NA12878 (1000 Genomes phase 3, CEPH/Utah)
- Accession: ERR1019034
- Chromosome: 22 (subset, n=100,000 reads from R1)

## Methods
[150 words. Tools pinned, parameters explicit. Cite Biopython 1.83,
FastQC 0.12.1, seqkit 2.8, your `trim.py` with parameters.]

## Results
- Reads in: 100,000. Reads after trim: 97,834.
- Mean read length: 100.0 -> 91.4 bp.
- Per-base mean Phred at positions 1 / 35 / 75: 36.2 / 35.1 / 29.3 (raw).
- ...

## One observation worth investigating
[100 words on a specific FastQC module — a "warn" or "fail" you looked
at — and your interpretation. Cite exact numbers.]

## Reproducibility
[the receipt block, below]
```

Aim for **clear, quantitative, undramatic**. C10 voice.

Commit: `Week 2 report v1`.

### Phase 6 — Reproducibility receipt and polish (~1 hour)

Add the receipt block at the bottom of `week-02/README.md`:

```
┌───────────────────────────────────────────────────────────────┐
│  REPRODUCIBILITY                                              │
│                                                               │
│  Data source:   1000 Genomes phase 3, ERR1019034, chr22 R1    │
│                  subset (100k reads, seqkit head -n 100000)   │
│  Pipeline:      bash run.sh                                   │
│  Container:     conda env (env.yml pinned to FastQC 0.12.1,   │
│                  Biopython 1.83, seqkit 2.8)                  │
│  Command:       bash run.sh data/sub.fastq.gz                 │
│  Wall time:     ~6 min (M1, 16 GB RAM)                        │
└───────────────────────────────────────────────────────────────┘
```

Sanity-check: every file path in the README is correct, every URL works, `bash run.sh` succeeds on a fresh checkout. Commit: `Mini-project reproducibility receipt`.

### Phase 7 — Update the repo root README (~30 min)

Add a Week 2 section to your portfolio repo's top-level `README.md` linking to `week-02/README.md` with a one-paragraph summary. Commit: `Week 2 entry in portfolio README`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Reproducibility | 25% | `bash run.sh` works on a fresh checkout. Versions pinned. `data/` excluded. |
| Quantitative report | 20% | Every claim has a number. C10 voice throughout. |
| FastQC interpretation | 15% | At least one module discussed thoughtfully — not just "looked good." |
| Code quality | 15% | Scripts have docstrings, accept CLI args, fail loudly on missing input. |
| Voice and precision | 10% | Reads like a methods section. No determinism. No "looked clean." |
| Repo hygiene | 10% | env.yml present, no committed raw FASTQ, commit history meaningful. |
| Plot quality | 5% | mean_quality_*.png are labelled axes, Q20 reference line, no default matplotlib title. |

---

## What this prepares you for

- **Week 3** assumes you can read a FASTA and bring sequence pairs into NumPy.
- **Week 5** assumes a FASTQ has been QC'd and trimmed before alignment — this is exactly the QC step in front of every aligner.
- **Week 12 capstone** has a "Data QC" section. It will look like this report, scaled up.

---

## Submission

When done:

1. Confirm `bash run.sh` works on a fresh clone with only `env.yml` and a small input FASTQ.
2. Confirm `week-02/README.md` renders cleanly on GitHub.
3. Confirm `data/` is `.gitignored` and you have not accidentally committed raw FASTQ.
4. Push.
5. Open Week 3 — only after the report is committed and the receipt block has been filled in honestly.

---

## Resources

- [1000 Genomes data portal](https://www.internationalgenome.org/data-portal/sample) — where the FASTQ lives.
- [ENA FASTQ download](https://www.ebi.ac.uk/ena/browser/view/ERR1019034) — direct gzipped FASTQ URLs.
- [FastQC docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) — read each module before interpreting.
- [seqkit](https://bioinf.shenwei.me/seqkit/) — for `head`, `stats`, sanity checks.
- The C10 [`branding/BRAND.md`](../../../branding/BRAND.md) — voice and receipt format.
- [`resources.md`](../resources.md) — full week resource list.
