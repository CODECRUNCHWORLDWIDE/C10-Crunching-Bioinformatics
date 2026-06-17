# Exercise 1 — Trim a paired-end RNA-seq FASTQ with fastp

> **Estimated time:** 45 minutes (15 minutes setup + download, 15 minutes running fastp, 15 minutes reading the report).
> **Goal:** Download a small yeast RNA-seq FASTQ pair from the SRA, run `fastp` end to end with the canonical paired-end flags, read both the HTML and JSON QC reports, and decide whether the sample is good enough to pseudoalign. Save the QC numbers — the mini-project methods section reuses them.

This exercise is the first hands-on step of the week. You produce two trimmed FASTQs and one QC report. You will reuse these outputs in Exercise 2 (kallisto pseudoalignment) and in the mini-project (3-sample matrix).

---

## Background

Raw Illumina RNA-seq reads need adapter trimming, low-quality 3'-tail trimming, polyG-tail trimming (NovaSeq artifact), and short-read filtering before pseudoalignment or alignment. `fastp` (Chen et al. 2018, *Bioinformatics* 34:i884) does all four in a single pass, in ~30 seconds per million paired reads, and emits an HTML + JSON QC report alongside the trimmed FASTQs. It is the modern replacement for the `Trimmomatic` + `FastQC` two-step.

Lecture 1 Section 4 introduced the canonical command:

```bash
fastp \
    -i raw/SRR12345_1.fq.gz \
    -I raw/SRR12345_2.fq.gz \
    -o trim/SRR12345_1.trim.fq.gz \
    -O trim/SRR12345_2.trim.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --trim_poly_g \
    -h qc/SRR12345.fastp.html \
    -j qc/SRR12345.fastp.json \
    -w 4
```

This exercise has you run it on one yeast sample (`SRR453568`, galactose growth, replicate 1), read the report, and report the numbers.

---

## Prerequisites

You have:

- A conda environment with `fastp 0.23.4` and the `sra-toolkit` (for FASTQ download) installed. The canonical install: `conda install -c bioconda fastp=0.23.4 sra-tools=3.0.10`.
- Python 3.11+ with `pandas` and the standard library, for reading the JSON report in Step 5.
- An internet connection capable of pulling ~250 MB from NCBI SRA.
- ~500 MB of free disk for the raw + trimmed FASTQs.
- Roughly 30 minutes of patience for the SRA download.

If you cannot install the SRA toolkit (e.g. on a locked-down corporate machine), use the bundled subsetted FASTQs that the mini-project README points at; they are ~150 MB per mate rather than ~1 GB.

---

## Step 1 — Layout the directory

In your portfolio repo, create:

```
crunch-bio-portfolio-<yourhandle>/
└── week-07/
    └── exercises/
        ├── raw/         (will contain SRR453568_1.fastq.gz, _2.fastq.gz)
        ├── trim/        (will contain SRR453568_1.trim.fq.gz, _2.trim.fq.gz)
        └── qc/          (will contain SRR453568.fastp.html, .json)
```

Run all commands from `week-07/exercises/`.

---

## Step 2 — Download the FASTQ

The yeast SRA accession we use is `SRR453568`. It is one of the Gierlinski et al. 2015 (*Bioinformatics* 31:3625) replicate samples — 3.2 M paired reads, 75 bp each, paired-end, polyA-selected, from *Saccharomyces cerevisiae* grown on galactose. We will reuse this sample throughout the week.

```bash
# Pull the SRA archive (~250 MB).
prefetch SRR453568 --output-directory raw/

# Extract paired-end FASTQ.
fasterq-dump raw/SRR453568/SRR453568.sra \
    --split-files \
    -O raw/

# Compress (fastp accepts gzipped input).
gzip raw/SRR453568_1.fastq raw/SRR453568_2.fastq
```

Verify:

```bash
ls -lh raw/
# Expected: SRR453568_1.fastq.gz (~150 MB), SRR453568_2.fastq.gz (~150 MB)

zcat raw/SRR453568_1.fastq.gz | head -4
# Expected: a FASTQ block - @SRR453568.1 ..., DNA seq, +, quality string
```

If `prefetch` is unavailable or too slow, the mini-project README also provides a no-cost direct HTTP mirror.

---

## Step 3 — Run fastp

The canonical paired-end fastp call:

```bash
mkdir -p trim qc

fastp \
    -i raw/SRR453568_1.fastq.gz \
    -I raw/SRR453568_2.fastq.gz \
    -o trim/SRR453568_1.trim.fq.gz \
    -O trim/SRR453568_2.trim.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --trim_poly_g \
    -h qc/SRR453568.fastp.html \
    -j qc/SRR453568.fastp.json \
    -w 4
```

Expected runtime: 25-45 seconds on a 4-core laptop, depending on disk speed.

Expected console output (the tail of the fastp log):

```
Read1 before filtering:
total reads: 3200000
total bases: 240000000
Q20 bases: 234567890(97.74%)
Q30 bases: 227890123(94.95%)
...
Read1 after filtering:
total reads: 3072412
...
Filtering result:
reads passed filter: 6144824
reads failed due to low quality: 81230
reads failed due to too many N: 412
reads failed due to too short: 173534
reads with adapter trimmed: 2256312
```

Verify outputs exist:

```bash
ls -lh trim/ qc/
# trim/SRR453568_1.trim.fq.gz  (~130 MB)
# trim/SRR453568_2.trim.fq.gz  (~130 MB)
# qc/SRR453568.fastp.html      (~600 KB)
# qc/SRR453568.fastp.json      (~200 KB)
```

---

## Step 4 — Read the HTML report

Open `qc/SRR453568.fastp.html` in a browser.

Walk through the sections in order:

1. **Summary** at the top. Record:
   - `total reads`: the number of reads in the raw FASTQ (forward + reverse, so for 3.2 M pairs you expect ~6.4 M).
   - `Q30 bases`: the percentage of bases with Q ≥ 30.
   - `duplication`: the estimated PCR/optical duplication rate.
   - `insert size peak`: the most-common fragment length.

2. **Filtering result**. Record:
   - `reads passed filter`: the number that made it through.
   - `reads with adapter trimmed`: how many had adapter.
   - Compute: `pct_retained = reads_passed / total_reads`.

3. **Quality (before/after)**. Look at the per-base quality boxplot. You want:
   - Median Q ≥ 30 across most of the read.
   - A slight drop in the last ~10 bp (this is normal).
   - **Bad signs**: Q crashes below 20 before position 60; the last 30 bp are essentially all low quality. If you see this, your sample is over-cycled and you should retrim with more aggressive `--cut_tail_window_size 5 --cut_tail_mean_quality 25`.

4. **Adapter cutting**. The fraction of reads where adapter was detected. For 75 bp reads with ~200 bp inserts, expect ~30%. If you see ≥ 70%, your insert size is shorter than the read length and many reads ran off the end of the cDNA — sample-prep issue.

5. **Insert size**. The fragment-length distribution. Should be unimodal with a peak at ~200-400 bp. A bimodal distribution suggests a library issue (e.g. two distinct fragmentation populations).

6. **Overrepresented sequences**. If you see rRNA k-mers in the top 10, your rRNA depletion / polyA selection failed and many reads are on rRNA. Document this; kallisto will still quantify the non-rRNA transcripts correctly, but your effective library size is much smaller than the raw read count suggests.

---

## Step 5 — Parse the JSON report

For programmatic QC, the JSON report is the right input. Create `week-07/exercises/fastp_summary.py`:

```python
"""Parse a fastp JSON report and print the QC fields that matter."""

from __future__ import annotations
import json
from pathlib import Path


def summarize_fastp_json(json_path: Path) -> dict[str, float | int | str]:
    """Read a fastp JSON report and return key QC fields."""
    with json_path.open() as f:
        data = json.load(f)

    summary = data["summary"]
    before = summary["before_filtering"]
    after = summary["after_filtering"]

    reads_in: int = before["total_reads"]
    reads_out: int = after["total_reads"]
    pct_retained: float = 100.0 * reads_out / max(reads_in, 1)
    pct_q30: float = 100.0 * after["q30_bases"] / max(after["total_bases"], 1)

    adapter = data.get("adapter_cutting", {})
    adapter_pct: float = 100.0 * adapter.get("adapter_trimmed_reads", 0) / max(reads_in, 1)

    dup = data.get("duplication", {})
    duplication_rate: float = dup.get("rate", 0.0)

    insert = data.get("insert_size", {})
    insert_size_peak: int = insert.get("peak", 0)

    return {
        "reads_in": reads_in,
        "reads_out": reads_out,
        "pct_retained": pct_retained,
        "pct_q30": pct_q30,
        "adapter_pct": adapter_pct,
        "duplication_rate": duplication_rate,
        "insert_size_peak": insert_size_peak,
    }


if __name__ == "__main__":
    summary = summarize_fastp_json(Path("qc/SRR453568.fastp.json"))
    for k, v in summary.items():
        print(f"  {k}: {v}")
```

Run with `python fastp_summary.py`. Expected output for a healthy SRR453568:

```
  reads_in: 6400000
  reads_out: 6144824
  pct_retained: 96.0
  pct_q30: 93.5
  adapter_pct: 35.2
  duplication_rate: 0.42
  insert_size_peak: 210
```

---

## Step 6 — Decide

Apply the Lecture 1 §7 health thresholds:

| Field | Healthy | Marginal | Bad |
|-------|---------|----------|-----|
| `pct_retained` | ≥ 90% | 75-90% | < 75% |
| `pct_q30` | ≥ 90% | 85-90% | < 85% |
| `adapter_pct` | 20-50% | 50-70% | > 70% |
| `duplication_rate` | < 0.60 | 0.60-0.80 | > 0.80 |
| `insert_size_peak` | 200-400 | 100-200 or 400-600 | < 100 or > 600 |

For SRR453568, all five fields should be in "Healthy." Record the numbers in `week-07/exercises/notes/e1-fastp.md`:

```markdown
# Exercise 1 — fastp QC for SRR453568

- Tool: fastp 0.23.4
- Input: SRR453568_1.fastq.gz + SRR453568_2.fastq.gz (3.2 M pairs)
- Output: SRR453568_1.trim.fq.gz + SRR453568_2.trim.fq.gz, plus HTML and JSON reports.

## QC summary

| Field               | Value     | Verdict |
|---------------------|----------:|---------|
| reads_in            | 6,400,000 | —       |
| reads_out           | 6,144,824 | —       |
| pct_retained        | 96.0%     | healthy |
| pct_q30             | 93.5%     | healthy |
| adapter_pct         | 35.2%     | healthy |
| duplication_rate    | 0.42      | healthy |
| insert_size_peak    | 210       | healthy |

Verdict: Sample is healthy. Proceed to Exercise 2 (kallisto pseudoalignment).
```

Commit `notes/e1-fastp.md`, `qc/SRR453568.fastp.html`, and `qc/SRR453568.fastp.json`. Gitignore `raw/` and `trim/` — they are too big and reconstructable from the SRA.

---

## Acceptance criteria

- [ ] `raw/SRR453568_1.fastq.gz` and `raw/SRR453568_2.fastq.gz` exist (or you used the no-cost mirror to fetch equivalent files).
- [ ] `trim/SRR453568_1.trim.fq.gz` and `trim/SRR453568_2.trim.fq.gz` exist.
- [ ] `qc/SRR453568.fastp.html` and `qc/SRR453568.fastp.json` exist.
- [ ] `fastp_summary.py` runs and prints the 7-field summary dict.
- [ ] `notes/e1-fastp.md` is committed with the QC table and verdict.
- [ ] Commit message like `e1: fastp trim SRR453568, 96.0% retained, Q30 93.5%`.

---

## Common pitfalls

**`prefetch` is very slow or hangs.** SRA Toolkit prefetch can be flaky depending on geography and time of day. Workarounds:
- Use `fasterq-dump` directly with the SRA accession (older fasterq-dump versions auto-prefetch): `fasterq-dump SRR453568 --split-files -O raw/`.
- Use the no-cost mirror linked from the mini-project README.
- Use a smaller subset: `head -n 800000 raw/SRR453568_1.fastq | gzip > raw/SRR453568_1.fastq.gz` (and matching for `_2`) — gives you 200k pairs, enough for the exercise but much smaller than the full sample.

**fastp emits "fastp tries to detect the adapter sequence, but failed."** This is fine — `--detect_adapter_for_pe` uses overlap analysis rather than sequence-matching, and the failure message just means the sequence-matching fallback would not be used. Trimming still works.

**The HTML report does not open.** Some browsers block local HTML files with JavaScript. Try opening from a file:// URL or copy the HTML to a temp directory served by `python -m http.server`.

**Q30 is 99%.** Check the input. Some SRA archives are already quality-trimmed; fastp on a pre-trimmed FASTQ trivially passes everything. This is fine but means your trimming is a no-op.

**Output FASTQs are smaller than expected.** Check `pct_retained`. A drop to 60-70% indicates a degraded sample or wrong parameters. Re-check `--qualified_quality_phred` (try 15 if 20 is too aggressive).

---

## What you learned

- `fastp` is the modern one-pass adapter + quality + polyG trimmer for Illumina paired-end RNA-seq.
- The seven QC fields (`reads_in`, `pct_retained`, `pct_q30`, `adapter_pct`, `duplication_rate`, `insert_size_peak`, plus over-represented sequences) decide whether a sample is good enough to quantify.
- The HTML report is for human eyes; the JSON report is for scripting and for the mini-project methods section.
- SRA download with `prefetch` + `fasterq-dump` is slow but free; for big projects, plan a multi-hour download window.
- Always pin tool versions: "fastp 0.23.4" in the methods, not "fastp."

Continue to [Exercise 2 — Pseudoalign with kallisto](./exercise-02-pseudoalign-with-kallisto.md).
