# Mini-Project — Pure-NumPy Smith-Waterman Benchmarked vs Biopython

> Build a production-quality, fully-traced, affine-gap-capable Smith-Waterman local aligner in pure NumPy. Benchmark it against `Bio.Align.PairwiseAligner` on a real DNA case study — a pair of SARS-CoV-2 spike-gene fragments from two patient isolates separated by a few thousand viral generations — and produce a one-page report with timings, peak memory, and a brief biological discussion.

This is the first C10 mini-project that produces a **runnable Python library**, not just an analytical artifact. By the end of it you will have a `smith_waterman.py` module you can import in the rest of the course, a benchmark harness that produces reproducible numbers, and a methods-section-quality report that compares your implementation to Biopython on a real DNA pair. Every later week of C10 that needs alignment (BLAST in Week 4, MSA in Week 9) will reach into this code.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-03/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 3 section)
└── week-03/
    ├── README.md                   one-page report (~700 words)
    ├── run.sh                      one-command reproduction script
    ├── env.yml                     conda environment file pinning all tool versions
    ├── data/
    │   ├── spike_isolate_a.fasta   (small enough to commit: ~3.8 kb)
    │   └── spike_isolate_b.fasta
    ├── smith_waterman.py           the pure-NumPy aligner module
    ├── benchmark.py                the benchmark harness
    ├── case_study.py               the SARS-CoV-2 spike pair script
    └── results/
        ├── benchmark.json          machine-readable timings
        ├── case_study.txt          the alignment with score and percent identity
        └── case_study.png          a dot-plot or score-trace figure
```

By the end you will have a clean, reproducible Week 3 directory you can point a recruiter at — and the `smith_waterman.py` inside is the kind of code that opens conversations with working bioinformaticians.

---

## The dataset

We use the SARS-CoV-2 spike gene (`S`) — the part of the viral genome that encodes the surface glycoprotein. Two isolates:

- **Isolate A — Wuhan-Hu-1 reference.** GenBank accession `NC_045512.2`, position `21563..25384` (3822 bp). The original 2020-01-05 reference.
- **Isolate B — a later Omicron BA.1 isolate.** GenBank accession `OL672836.1`, spike CDS. Roughly 50 substitutions and a couple of in-frame deletions across the ~3.8 kb spike — substantial divergence, but still alignable end to end.

You will work with the **spike gene only** to keep run times reasonable. Full-genome alignment of two 30 kb viral genomes is a Week-5 problem (it needs an indexed seed-and-extend aligner). 3.8 kb x 3.8 kb is the ceiling of "things you can do with pure NumPy Smith-Waterman in well under a minute on a laptop."

Both FASTA files are small enough (< 4 kb each) that they live inside the repo — no large-file workarounds, no `.gitignore` rules. Commit them.

### How to fetch the sequences

```bash
# Isolate A: spike gene (CDS) from the SARS-CoV-2 reference.
curl -L 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&seq_start=21563&seq_stop=25384' \
    -o data/spike_isolate_a.fasta

# Isolate B: spike gene from an Omicron BA.1 isolate.
curl -L 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=OL672836.1&rettype=fasta' \
    -o data/spike_isolate_b.fasta
```

Confirm with `seqkit stats data/spike_isolate_*.fasta` — each file should report 1 record, ~3.8 kb.

---

## Rules

- **You may** use NumPy 1.26, Biopython 1.83, matplotlib, and the standard library.
- **You may** consult Biopython's tutorial, the Smith-Waterman 1981 paper, and the lecture notes.
- **You may NOT** copy a pre-written aligner from the internet. The point is to *write* the algorithm. Reading scikit-bio's source for inspiration is fine; copy-pasting it is not.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file.

---

## Acceptance criteria

- [ ] `week-03/smith_waterman.py` exports a function `align(seq_a, seq_b, *, match=1, mismatch=-1, gap_open=-2, gap_extend=-1) -> AlignResult` where `AlignResult` is a `NamedTuple` or `dataclass` with fields `score: int`, `aligned_a: str`, `aligned_b: str`, `a_start: int`, `a_end: int`, `b_start: int`, `b_end: int`. Coordinates are 1-based, inclusive (BLAST convention).
- [ ] The implementation uses the Gotoh three-matrix recurrence (challenge 1) so it handles affine gaps correctly. If you only completed the linear-gap version of Smith-Waterman from Exercise 2, that is acceptable for the score-only path **but** your report must call out the limitation explicitly.
- [ ] `week-03/README.md` is a one-page (≤ 800 word) report containing:
  - One-sentence description of the dataset, accessions, and the biological question ("how diverged is BA.1 spike from the original reference?").
  - Methods section in C10 voice: every tool pinned ("Biopython 1.83", "NumPy 1.26"), every parameter explicit ("`match=+1, mismatch=-1, gap_open=-5, gap_extend=-1` under `Bio.Align.PairwiseAligner` mode='local'").
  - Quantitative findings: local-alignment score, length of aligned region, percent identity over the aligned region, number of gaps.
  - A benchmark table: wall-clock for NumPy implementation vs Biopython on three input sizes (500 bp, 1000 bp, 3800 bp), median of three runs each.
  - A reproducibility receipt block (see below).
- [ ] `week-03/run.sh` runs end-to-end on a fresh clone:
  - Activates the conda env from `env.yml`.
  - Fetches both FASTA files (or accepts them as already on disk).
  - Runs `case_study.py` to align the two spike sequences, writing `results/case_study.txt` and `results/case_study.png`.
  - Runs `benchmark.py` to time the NumPy and Biopython implementations on three sizes, writing `results/benchmark.json`.
- [ ] `week-03/env.yml` pins every tool to an exact version:
  ```yaml
  name: c10-week-03
  channels:
    - conda-forge
    - bioconda
  dependencies:
    - python=3.11
    - numpy=1.26.4
    - biopython=1.83
    - matplotlib
    - pip
  ```
- [ ] `week-03/results/benchmark.json` is valid JSON with the shape:
  ```json
  {
    "machine": "M1 MacBook Pro, 16 GB",
    "python": "3.11.x",
    "biopython": "1.83",
    "numpy": "1.26.4",
    "sizes": [500, 1000, 3800],
    "numpy_seconds_median": [0.04, 0.18, 2.4],
    "biopython_seconds_median": [0.001, 0.003, 0.04],
    "ratio_numpy_over_biopython": [40, 60, 60]
  }
  ```
- [ ] `week-03/results/case_study.txt` contains the EMBOSS-style three-line alignment plus a header block with score, percent identity, and aligned coordinates.
- [ ] `week-03/results/case_study.png` is a figure produced by your own matplotlib code — either a dot-plot, a score-trace, or a per-position match indicator. Labelled axes; no default matplotlib title.
- [ ] The repo passes a fresh-clone test: `git clone`, `cd week-03`, `bash run.sh` reproduces everything in `results/`.

---

## Suggested order of operations

### Phase 1 — Environment setup (~30 min)

1. Create `week-03/env.yml` (see acceptance criteria above for the exact contents).
2. `conda env create -f week-03/env.yml`. Activate it. Confirm: `python -c "import numpy; print(numpy.__version__)"` → `1.26.4`.
3. Commit: `Week 3 env.yml`.

### Phase 2 — Port and harden `smith_waterman.py` (~2 hours)

1. Start from your Exercise 2 implementation.
2. Add affine-gap support (the Gotoh three-matrix recurrence from Challenge 1) if you have time; otherwise document the linear-gap restriction in the report.
3. Add a proper module docstring citing Smith & Waterman 1981 and the BLAST defaults.
4. Add a `AlignResult` `NamedTuple` so callers get a structured return rather than a tuple.
5. Add input validation: reject empty sequences (or handle them gracefully — return a zero-score result with empty alignment), uppercase the inputs, raise on non-DNA characters (or accept `N` and treat it as a wildcard with score 0).
6. Write a few inline unit tests at the bottom of the file (`if __name__ == "__main__":`) that confirm the score against the Lecture 1 §4.7 by-hand example and a couple more.
7. Commit: `smith_waterman.py with affine gap support`.

### Phase 3 — Case study script (~1 hour)

1. Write `case_study.py` that:
   - Loads `data/spike_isolate_a.fasta` and `data/spike_isolate_b.fasta` with `Bio.SeqIO.read`.
   - Aligns them with `smith_waterman.align(...)` using `match=+1, mismatch=-1, gap_open=-5, gap_extend=-1`.
   - Computes percent identity over the aligned region (count of match columns / total aligned columns excluding gaps).
   - Writes the alignment to `results/case_study.txt` in EMBOSS-style three-line format (top sequence, match indicator with `|` for match and `.` for mismatch, bottom sequence — wrap at 60 cols).
   - Produces a dot-plot of the alignment with matplotlib (one dot per match column at coordinates `(a_pos, b_pos)`; mismatches and gaps absent). Save as `results/case_study.png`.
2. Run it manually once. Confirm the percent identity is in the 95-99% range (BA.1 spike is mostly conserved against the reference at the nucleotide level; only ~50 of ~3800 positions differ).
3. Commit: `case_study.py for spike pair`.

### Phase 4 — Benchmark harness (~1.5 hours)

1. Write `benchmark.py` that:
   - Generates three pairs of random DNA sequences at lengths 500, 1000, and 3800 (use a seeded `random.Random(...)` for reproducibility — your numbers should be the same on a re-run).
   - For each pair, times `smith_waterman.align(...)` and `Bio.Align.PairwiseAligner(...).score(...)` three times each, takes the median.
   - Writes `results/benchmark.json` in the shape shown above.
   - Prints a Markdown table to stdout.
2. Run it. Expect your NumPy implementation to be ~40-100x slower than Biopython's C extension. **This is fine.** The point of the benchmark is not to win against Biopython; it is to *know* the ratio honestly.
3. Commit: `benchmark.py with three-size ratio`.

### Phase 5 — Wire up `run.sh` (~30 min)

Aim for ~25 lines of straight bash:

```bash
#!/usr/bin/env bash
set -euo pipefail

# 1) Confirm data is on disk; fetch if not.
if [ ! -f data/spike_isolate_a.fasta ]; then
    bash scripts/fetch_data.sh   # or curl one-liners here
fi

mkdir -p results

# 2) Run the case study.
python case_study.py

# 3) Run the benchmark.
python benchmark.py

echo "Done. See results/ for outputs."
```

Commit: `run.sh end-to-end`.

### Phase 6 — Write the report (~2 hours)

Open `week-03/README.md`. Structure:

```
# Week 3 — Pairwise alignment: SARS-CoV-2 spike-gene case study

## Sequences and accessions
- Isolate A: NC_045512.2, spike CDS, positions 21563..25384 (3822 bp)
- Isolate B: OL672836.1, spike CDS (3819 bp)

## Methods
[150 words. Smith-Waterman, mode='local', match +1, mismatch -1,
gap_open -5, gap_extend -1. Implementation: pure NumPy with the
Gotoh three-matrix affine-gap recurrence. Cite Biopython 1.83 for the
verification path.]

## Findings
- Local-alignment score: 3712
- Aligned length: 3819 columns (essentially end-to-end on isolate B)
- Percent identity over aligned region: 96.8%
- Number of gaps: 2 (a 9 bp deletion at position ~21000 and a 3 bp
  deletion at position ~22300, both in-frame; consistent with the
  spike's NTD-loop variability in Omicron BA.1)

## Benchmark
[the three-size table, median wall-clock]

## Reproducibility
[the receipt block, below]
```

Aim for **clear, quantitative, undramatic**. C10 voice. No "the alignment looks good"; give numbers.

Commit: `Week 3 report v1`.

### Phase 7 — Reproducibility receipt and polish (~30 min)

Add the receipt block at the bottom of `week-03/README.md`:

```
┌───────────────────────────────────────────────────────────────┐
│  REPRODUCIBILITY                                              │
│                                                               │
│  Data source:   NC_045512.2 spike CDS, OL672836.1 spike CDS   │
│  Pipeline:      bash run.sh                                   │
│  Container:     conda env (env.yml pinned to NumPy 1.26.4,    │
│                  Biopython 1.83, Python 3.11)                 │
│  Command:       bash run.sh                                   │
│  Wall time:     ~12 s (M1, 16 GB RAM)                         │
└───────────────────────────────────────────────────────────────┘
```

Sanity-check: every file path in the README is correct, every URL works, `bash run.sh` succeeds on a fresh checkout. Commit: `Mini-project reproducibility receipt`.

### Phase 8 — Update the repo root README (~30 min)

Add a Week 3 section to your portfolio repo's top-level `README.md` linking to `week-03/README.md` with a one-paragraph summary. Commit: `Week 3 entry in portfolio README`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Correctness | 25% | NumPy scores match Biopython to the integer on all benchmark sizes |
| Reproducibility | 20% | `bash run.sh` works on a fresh checkout. Versions pinned. Data fetched or committed. |
| Code quality | 15% | `smith_waterman.py` has a module docstring, structured return type, input validation, inline tests |
| Quantitative report | 15% | Every claim has a number. C10 voice throughout. |
| Benchmark honesty | 10% | Real timings from your machine; ratio explained, not glossed |
| Voice and precision | 10% | Reads like a methods section. No "fast" without seconds, no "good" without numbers |
| Plot quality | 5% | `case_study.png` is labelled, captioned, and biologically interpretable |

---

## What this prepares you for

- **Week 4** assumes you can score an alignment with `Bio.Align.PairwiseAligner` and read a BLAST tabular hit — your `case_study.py` is the rehearsal for that.
- **Week 9** revisits alignment in the form of *multiple* sequence alignment (MAFFT, MUSCLE). The pairwise core is the same recurrence.
- **Week 12 capstone** will benchmark whatever pipeline you build against an established baseline. The honest-ratio table you produce here is the template.

---

## Submission

When done:

1. Confirm `bash run.sh` works on a fresh clone.
2. Confirm `week-03/README.md` renders cleanly on GitHub.
3. Confirm `results/benchmark.json` is valid JSON.
4. Confirm `results/case_study.png` is committed and is a meaningful figure.
5. Push.
6. Open Week 4 — only after the report is committed and the receipt block has been filled in honestly.

---

## Resources

- [NCBI eutils efetch](https://www.ncbi.nlm.nih.gov/books/NBK25497/) — for the data fetch.
- [Bio.Align.PairwiseAligner docs](https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.PairwiseAligner) — the verification path.
- [Smith & Waterman 1981 paper](https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/SmithWaterman1981.pdf) — four pages, free PDF.
- [Gotoh 1982 affine-gap algorithm](https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf) — lecture notes citing the original.
- The C10 [`branding/BRAND.md`](../../../branding/BRAND.md) — voice and receipt format.
- [`resources.md`](../resources.md) — full week resource list.
