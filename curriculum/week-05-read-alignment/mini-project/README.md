# Mini-Project — Align a Public Small-Genome Dataset and Produce Coverage Plots

> Build a reproducible FASTQ → sorted, duplicate-marked BAM → coverage-plot pipeline on a real public SRA dataset (Illumina paired-end resequencing of *E. coli* K-12 MG1655, SRA accession `SRR1770413` aligned against `NC_000913.3`). End with a duplicate-marked sorted BAM, per-position coverage statistics, a windowed coverage plot, and a short written interpretation of the coverage distribution.

This is the first C10 mini-project that produces a **methods-section-quality alignment with measured QC**, not just a single-read demonstration. By the end of it you will have a `align_and_qc.py` script and a `run.sh` wrapper you can point a recruiter at, a results directory with a coverage plot and `flagstat`/`markdup` summaries, and a write-up that defends each pipeline parameter and names the failure modes you observed.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-05/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 5 section)
└── week-05/
    ├── README.md                   one-page report (~800 words)
    ├── run.sh                      one-command reproduction script
    ├── env.yml                     conda environment file pinning all tool versions
    ├── data/
    │   ├── ecoli.fa                NC_000913.3 reference FASTA
    │   ├── ecoli.fa.{amb,ann,bwt,pac,sa,fai}   BWA + samtools indexes
    │   ├── SRR1770413_1.fq.gz      paired-end R1 reads (gitignored, large)
    │   └── SRR1770413_2.fq.gz      paired-end R2 reads (gitignored, large)
    ├── align_and_qc.py             the alignment + QC pipeline
    ├── plot_coverage.py            the coverage-plot renderer
    ├── aln/
    │   ├── SRR1770413.sorted.bam       coordinate-sorted (gitignored)
    │   ├── SRR1770413.markdup.bam      duplicate-marked (gitignored)
    │   └── SRR1770413.markdup.bam.bai  index for the markdup BAM
    └── results/
        ├── flagstat.txt            samtools flagstat output
        ├── markdup_stats.txt       samtools markdup -s summary
        ├── coverage_summary.tsv    samtools coverage output
        ├── depth.tsv.gz            samtools depth -a output (compressed)
        ├── coverage_plot.png       windowed coverage plot
        ├── mapq_histogram.png      MAPQ distribution histogram
        └── insert_size_histogram.png  paired-end insert size histogram
```

By the end you will have a clean, reproducible Week 5 directory you can point a recruiter at — and `align_and_qc.py` is the kind of pipeline that opens conversations with working bioinformaticians and clinical-sequencing shops.

---

## The dataset

You will work with the **`SRR1770413`** SRA run: Illumina HiSeq 2500 paired-end resequencing of *E. coli* K-12 MG1655, 2 × 150 bp reads, ~5 GB compressed. The reference is **`NC_000913.3`** (the *E. coli* K-12 MG1655 complete genome from RefSeq, 4,641,652 bp). This is the same reference you used in Week 4 Exercise 2, so the FASTA fetch is already familiar.

The run was uploaded to SRA in 2014 as part of a comparative-genomics study; the exact biology is not important for our purposes — we treat it as "generic well-behaved Illumina paired-end data from a known bacterial reference." Mean expected coverage when aligning all reads is ~200x (far more than needed); for the mini-project we subsample to ~50x to keep the BAM under 1 GB.

Fetch:

```bash
fasterq-dump --split-files -p --threads 4 -X 750000 SRR1770413 -O data/
gzip data/SRR1770413_1.fastq data/SRR1770413_2.fastq
# Rename to match the pipeline:
mv data/SRR1770413_1.fastq.gz data/SRR1770413_1.fq.gz
mv data/SRR1770413_2.fastq.gz data/SRR1770413_2.fq.gz
```

750,000 read pairs × 300 bp = 225 Mb of read data against a 4.6 Mb reference = ~49x mean coverage. That is the working dataset.

### Why E. coli and not human

*E. coli* MG1655 is the standard small-bacterial-reference benchmark for three reasons: (1) the genome is small (4.6 Mb) so the full pipeline fits on a laptop and runs in under 10 minutes, (2) the reference is well-curated and the biology is well-studied (every gene is named), and (3) the rRNA-operon repetitive regions produce a characteristic coverage signature that lets you check whether your pipeline is behaving correctly. Week 6 will move to a *human* small-region dataset (a 1 Mb slice of chr22 from 1000 Genomes); Week 5's *E. coli* is the warm-up.

---

## Rules

- **You may** use BWA 0.7.17, minimap2 2.26, samtools 1.19, pysam 0.22, Biopython 1.83, pandas, matplotlib, numpy, and the standard library.
- **You may** consult Lectures 1 and 2, the BWA-MEM paper (Li 2013), the SAM/BAM spec, the samtools documentation, and your Week 5 exercises and challenge.
- **You may NOT** copy a pre-written alignment pipeline from the internet. The point is to *build* the pipeline. Reading the GATK Best Practices for inspiration is fine; copy-pasting a `nf-core/sarek` config is not.
- **You must** cache the reference FASTA, BWA index, and BAM files to disk. A second run of `bash run.sh` should not re-fetch the reference or re-build the index.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file. The only network access on a fresh run should be the initial reference fetch and SRA download.

---

## Acceptance criteria

- [ ] `week-05/align_and_qc.py` exports a function `align(ref_fa, r1, r2, out_dir, *, threads=4, read_group=None) -> Path` that runs the full FASTQ → markdup BAM pipeline and returns the path to the duplicate-marked sorted+indexed BAM.
- [ ] The pipeline implements **five stages**:
  1. `bwa index` (skip if index files exist).
  2. `bwa mem -t N -R '@RG\tID:<id>\tSM:<sample>\tLB:<library>\tPL:ILLUMINA' ref r1 r2 | samtools sort -@ N -o sorted.bam -`.
  3. `samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup -s` (the four-step duplicate-marking idiom).
  4. `samtools index` on the markdup BAM.
  5. `samtools flagstat` + `samtools coverage` + `samtools depth -a` for QC summaries.
- [ ] `week-05/plot_coverage.py` produces three PNGs:
  - `results/coverage_plot.png` — windowed coverage across the reference, with mean and median lines.
  - `results/mapq_histogram.png` — distribution of MAPQ values (expect bimodal at 0 and 60).
  - `results/insert_size_histogram.png` — distribution of `TLEN` for properly-paired reads (expect a clean peak around the library's median insert size, typically 300–500 bp for Illumina TruSeq).
- [ ] `week-05/README.md` is a one-page (≤ 1,000 word) report containing:
  - One-sentence description of the dataset, the reference, and the pipeline's five stages.
  - Methods section in C10 voice: every tool pinned ("BWA 0.7.17", "samtools 1.19", "minimap2 2.26", "pysam 0.22"), every parameter explicit ("`bwa mem -t 4 -R '@RG\\tID:SRR1770413\\tSM:ecoli\\tLB:lib1\\tPL:ILLUMINA'`", "`samtools markdup -s`").
  - Quantitative findings: mean coverage, median coverage, CV, mapping rate, properly-paired rate, duplicate rate, MAPQ-60 fraction. Each value to one decimal place.
  - A coverage-plot interpretation: describe the overall shape, identify any sharp dips or spikes, name the biological cause if known (rRNA operons, prophage insertions).
  - A failure-modes section listing every QC signal that surprised you, with one-sentence diagnoses.
  - A reproducibility receipt block.
- [ ] `week-05/run.sh` runs end-to-end on a fresh clone:
  - Activates the conda env from `env.yml`.
  - Downloads the *E. coli* reference via `Bio.Entrez.efetch` if not already cached on disk.
  - Builds the BWA index if not already built.
  - Downloads the SRA reads via `fasterq-dump` if not already on disk.
  - Runs `align_and_qc.py` to align, sort, mark duplicates, and produce QC summaries.
  - Runs `plot_coverage.py` to render the three PNGs.
- [ ] `week-05/env.yml` pins every tool to an exact version:
  ```yaml
  name: c10-week-05
  channels:
    - conda-forge
    - bioconda
  dependencies:
    - python=3.11
    - numpy=1.26.4
    - matplotlib
    - pandas
    - biopython=1.83
    - bwa=0.7.17
    - minimap2=2.26
    - samtools=1.19
    - pysam=0.22
    - sra-tools=3.0.10
    - pip
  ```
- [ ] `week-05/results/flagstat.txt` contains a `samtools flagstat` summary with > 95% mapped.
- [ ] `week-05/results/markdup_stats.txt` contains a `samtools markdup -s` summary with a duplicate rate < 25%.
- [ ] `week-05/results/coverage_summary.tsv` contains a `samtools coverage` per-contig summary with `coverage` (the column) > 99%.
- [ ] `week-05/results/coverage_plot.png` is a matplotlib plot with axis labels, a title naming the SRA run and reference, and visible mean/median reference lines.
- [ ] The repo passes a fresh-clone test: `git clone`, `cd week-05`, `bash run.sh` reproduces everything in `results/` (modulo a ±0.1% difference in counts if NCBI's reference has been updated).

---

## Suggested order of operations

### Phase 1 — Environment setup (~30 min)

1. Create `week-05/env.yml` (see acceptance criteria).
2. `conda env create -f week-05/env.yml`. Activate it. Confirm: `python -c "import pysam; print(pysam.__version__)"` → `0.22`, `bwa 2>&1 | head -3` → `Version: 0.7.17`, `samtools --version | head -1` → `samtools 1.19`, `minimap2 --version` → `2.26-r1175`.
3. Commit: `Week 5 env.yml`.

### Phase 2 — Fetch the reference and reads (~30–60 min)

1. Write `scripts/fetch_reference.py` that calls `Bio.Entrez.efetch(db="nuccore", id="NC_000913.3", rettype="fasta")` and writes `data/ecoli.fa`. Set `Bio.Entrez.email`. If `data/ecoli.fa` already exists, skip.
2. Write `scripts/fetch_reads.sh` that runs `fasterq-dump --split-files -p --threads 4 -X 750000 SRR1770413 -O data/` and gzips the output. If `data/SRR1770413_1.fq.gz` already exists, skip.
3. Verify file sizes: reference ~5 MB, R1 + R2 together ~280 MB compressed.
4. Commit: `Reference and reads fetched`.

### Phase 3 — Build the alignment pipeline (~2 hours)

1. Write `align_and_qc.py` with the `align(...)` function. Outline:
   ```python
   def align(ref_fa: Path, r1: Path, r2: Path, out_dir: Path, *,
             threads: int = 4, sample: str = "ecoli_K12",
             library: str = "lib1", run_id: str = "SRR1770413") -> Path:
       out_dir.mkdir(parents=True, exist_ok=True)
       # Stage 1: index reference (idempotent).
       if not Path(f"{ref_fa}.bwt").exists():
           subprocess.run(["bwa", "index", str(ref_fa)], check=True)
       if not Path(f"{ref_fa}.fai").exists():
           subprocess.run(["samtools", "faidx", str(ref_fa)], check=True)
       # Stage 2: bwa mem | samtools sort -.
       sorted_bam = out_dir / f"{run_id}.sorted.bam"
       rg = (f"@RG\\tID:{run_id}\\tSM:{sample}\\tLB:{library}"
             "\\tPL:ILLUMINA")
       cmd = (f"bwa mem -t {threads} -R '{rg}' {ref_fa} {r1} {r2} | "
              f"samtools sort -@ {threads} -o {sorted_bam} -")
       subprocess.run(cmd, shell=True, check=True)
       # Stages 3-4: name-sort | fixmate -m | sort | markdup.
       markdup_bam = out_dir / f"{run_id}.markdup.bam"
       cmd = (f"samtools sort -n -@ {threads} {sorted_bam} | "
              f"samtools fixmate -m - - | "
              f"samtools sort -@ {threads} - | "
              f"samtools markdup -s - {markdup_bam}")
       subprocess.run(cmd, shell=True, check=True)
       subprocess.run(["samtools", "index", str(markdup_bam)], check=True)
       return markdup_bam
   ```
2. Add input validation: every path must exist; threads must be ≥ 1; the read-group fields must be non-empty.
3. Test on a tiny subset (1000 reads from `head -4000 SRR1770413_1.fq | gzip > test_R1.fq.gz`) to verify the pipeline works before running on the full 750k pairs.
4. Commit: `align_and_qc.py end to end on test subset`.

### Phase 4 — Run on the full dataset (~10–30 min wall-clock)

1. `python align_and_qc.py data/ecoli.fa data/SRR1770413_1.fq.gz data/SRR1770413_2.fq.gz aln/`.
2. Confirm `aln/SRR1770413.markdup.bam` and `aln/SRR1770413.markdup.bam.bai` exist.
3. Confirm via `samtools flagstat aln/SRR1770413.markdup.bam`:
   - Total reads ≈ 1.5 M (2 × 750k).
   - Mapped > 95%.
   - Properly paired > 90%.
   - Duplicates 5–25% (varies by run).
4. Commit: `Aligned SRR1770413 to NC_000913.3`.

### Phase 5 — QC summaries (~30 min)

1. Run `samtools flagstat aln/SRR1770413.markdup.bam > results/flagstat.txt`.
2. Run `samtools coverage aln/SRR1770413.markdup.bam > results/coverage_summary.tsv`.
3. Run `samtools depth -a aln/SRR1770413.markdup.bam | gzip > results/depth.tsv.gz`.
4. Run `samtools markdup -s` once more on a known input to capture the stats text and save to `results/markdup_stats.txt` (or capture from the `align_and_qc.py` run via `subprocess.run(..., capture_output=True)` — cleaner).
5. Commit: `QC summaries`.

### Phase 6 — Plots (~1 hour)

1. Write `plot_coverage.py` that takes `results/depth.tsv.gz` as input and produces three PNGs:
   - **Coverage plot**: read the depth TSV, bin into 5 kb windows, plot the per-window mean depth. Add horizontal lines at mean and median. The plot should be 12 inches wide × 4 inches tall at 150 DPI.
   - **MAPQ histogram**: open the BAM with pysam, iterate over reads, collect `mapping_quality`, plot as `plt.hist(mapqs, bins=range(0, 62))`. Expect bimodal: most reads at 60, a few at 0, and a thin tail in between.
   - **Insert-size histogram**: filter to properly-paired primary reads, collect `template_length` (taking absolute value), plot as `plt.hist(insert_sizes, bins=200, range=(0, 1000))`. Expect a tight peak around 300–500 bp (the library's median fragment size).
2. Commit: `Three QC plots`.

### Phase 7 — Write the report (~2 hours)

Open `week-05/README.md`. Structure:

```
# Week 5 — Read alignment of SRR1770413 against NC_000913.3

## Dataset and biological question
- Read set: SRR1770413, 2 × 150 bp Illumina paired-end resequencing of
  E. coli K-12 MG1655, 750,000 read pairs (~225 Mb of read data,
  expected ~49x mean coverage).
- Reference: NC_000913.3, E. coli K-12 MG1655 complete genome,
  4,641,652 bp.
- Question: can a standard bwa mem + markdup pipeline produce a uniform-
  coverage BAM suitable for downstream variant calling on a small
  well-behaved bacterial dataset?

## Methods
[200 words. BWA 0.7.17 with default parameters, -R read group as
specified. samtools 1.19 sort and fixmate -m -> markdup four-step
idiom. samtools depth -a with positions reported including zero-
coverage bases. Three plots rendered in matplotlib at 150 DPI.
Implementation in pure Python via subprocess; no nf-core or
snakemake.]

## Findings
- Mapping rate: 99.X% (from flagstat)
- Proper-pair rate: 9X.X%
- Duplicate rate: X.X% (from markdup -s)
- Mean coverage: 4X.Xx (samtools coverage `meandepth` column)
- Median coverage: 4X.Xx (from depth.tsv.gz)
- Coefficient of variation: 0.XX
- MAPQ-60 fraction: 9X.X% of all primary alignments
- Fraction of reference bases at >= 20x: 9X.X%

## Coverage interpretation
[150 words. The overall shape is roughly uniform around the mean.
Six characteristic dips at the rRNA operons rrnA-rrnH (E. coli has
7 rRNA operons at canonical positions on the chromosome) drop
coverage to ~5x. No other large-scale features. The MAPQ
distribution is bimodal at 0 (multi-mapped reads, ~3% of reads,
primarily in the rRNA operons and the IS1 element) and 60
(unique placement, ~97%).]

## Failure modes observed
- rRNA operons (positions ~4.0 Mb and ~4.2 Mb): MAPQ-0 multimappers
  cause the characteristic coverage dip; this is expected.
- IS1 transposon copies: similar multimapper signature at smaller
  positions; documented in the K-12 reference.
- [Any other surprises]

## Reproducibility
[the receipt block, below]
```

Aim for **clear, quantitative, undramatic**. C10 voice. No "the coverage looks good"; give numbers.

Commit: `Week 5 report v1`.

### Phase 8 — Reproducibility receipt and polish (~30 min)

Add the receipt block at the bottom of `week-05/README.md`:

```
+--------------------------------------------------------------------+
|  REPRODUCIBILITY                                                   |
|                                                                    |
|  Data source:   SRA SRR1770413 (downloaded YYYY-MM-DD)              |
|                  + RefSeq NC_000913.3 (E. coli K-12 MG1655)        |
|  Pipeline:      bash run.sh                                        |
|  Container:     conda env (env.yml pinned to BWA 0.7.17,           |
|                  samtools 1.19, minimap2 2.26, pysam 0.22,         |
|                  Python 3.11)                                      |
|  Command:       bash run.sh                                        |
|  Wall time:     ~12 min (M1, 16 GB RAM, 4 threads, network for     |
|                  first run only)                                   |
+--------------------------------------------------------------------+
```

Sanity-check: every file path in the README is correct, every URL works, `bash run.sh` succeeds on a fresh checkout. Commit: `Mini-project reproducibility receipt`.

### Phase 9 — Update the repo root README (~15 min)

Add a Week 5 section to your portfolio repo's top-level `README.md` linking to `week-05/README.md` with a one-paragraph summary. Commit: `Week 5 entry in portfolio README`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Correctness | 25% | Mapping rate > 99%, properly-paired > 95%, duplicate rate < 20%, mean coverage within 10% of expected. Pipeline runs end to end without manual intervention. |
| Reproducibility | 20% | `bash run.sh` works on a fresh clone. Versions pinned. Reference and reads cached. |
| Code quality | 15% | `align_and_qc.py` has module docstrings, structured return types, input validation, subprocess errors propagated. Idempotent on re-run. |
| Quantitative report | 15% | Every claim has a number. C10 voice throughout. Coverage interpretation names specific biological features. |
| Plot quality | 15% | All three plots are labelled, captioned, biologically interpretable. Mean/median reference lines on the coverage plot. Bimodal MAPQ shows clearly. |
| Voice and precision | 10% | Reads like a methods section. No "fast" without seconds, no "good" without numbers. |

---

## What this prepares you for

- **Week 6** runs `bcftools mpileup + bcftools call` against your Week 5 sorted+markdup BAM to produce a VCF of variants vs the reference. The coverage and duplicate-rate signals from this week directly drive the hard-filter thresholds in Week 6.
- **Week 7** counts reads against gene-feature intervals (RNA-seq), reusing the same pysam patterns from this week.
- **Week 11** asks you to convert this pipeline to Snakemake. The current `run.sh` is the rehearsal: every step is a `subprocess.run(...)` invocation that can become a Snakemake rule.
- **Week 12 capstone** will likely include some form of alignment. The honest QC methodology you build here is the template.

---

## Common pitfalls

A short list, from instructor experience:

- **Forgetting `-R` read group.** Variant callers in Week 6 will refuse to run, or will silently assume one sample, on a BAM without read groups. Set it now and you will not have to come back.
- **Using `samtools markdup` without `samtools fixmate -m` first.** It will run but silently produce wrong duplicate calls (missing the `ms` tag means markdup falls back to a less-accurate heuristic). Always sort-by-name → fixmate -m → sort-by-coord → markdup.
- **Forgetting `samtools depth -a`.** Without `-a`, zero-coverage positions are silently dropped from the output, biasing your mean coverage upward by `~fraction_of_zero_bases · current_mean`. For a healthy E. coli alignment this is < 1% but on a poorly-covered region it can be 20%+.
- **Hardcoding paths.** `run.sh` should compute paths relative to its own location (`SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)`) so a fresh clone in any directory works.
- **Committing the 2 GB BAM.** Add `*.bam`, `*.bai`, `*.fq.gz` to `.gitignore`. Commit the smaller text outputs (`flagstat.txt`, `coverage_summary.tsv`, `markdup_stats.txt`) and the PNGs.

---

## Submission

When done:

1. Confirm `bash run.sh` works on a fresh clone (delete `aln/` and re-run; it should rebuild from cached reference + reads).
2. Confirm `week-05/README.md` renders cleanly on GitHub.
3. Confirm `results/coverage_plot.png`, `results/mapq_histogram.png`, `results/insert_size_histogram.png` are committed and meaningful.
4. Confirm `results/flagstat.txt`, `results/coverage_summary.tsv`, `results/markdup_stats.txt` are committed.
5. Push.
6. Open Week 6 — only after the report is committed and the failure-modes section explains every QC surprise.

---

## Resources

- [BWA manual](https://bio-bwa.sourceforge.net/bwa.shtml) — the command-line reference.
- [BWA-MEM paper](https://arxiv.org/abs/1303.3997) — Li 2013.
- [minimap2 paper](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778) — Li 2018.
- [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) — the canonical 25-page reference.
- [samtools documentation](http://www.htslib.org/doc/samtools.html) — every command with examples.
- [pysam documentation](https://pysam.readthedocs.io/) — the Python BAM-reader API.
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912) — the industry-standard pipeline for short-variant discovery (the duplicate-marking step in particular is worth reading).
- [Lecture 1 — From FASTQ to BAM](../lecture-notes/01-from-fastq-to-bam.md)
- [Lecture 2 — SAM/BAM format and samtools](../lecture-notes/02-sam-bam-format-and-samtools.md)
- [Challenge 1 — Detect duplicates by hand](../challenges/challenge-01-detect-duplicates.md)
- [resources.md](../resources.md) — full week resource list.
