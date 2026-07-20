# Week 2 — Exercises

Three focused Python drills. Each one runs from the command line, finishes inside 45 minutes, and produces output you can keep.

## Index

1. **[Exercise 1 — Parse FASTA with Biopython](exercise-01-parse-fasta-biopython.py)** — re-do Week 1's FASTA parser using `Bio.SeqIO` and compare. (~35 min)
2. **[Exercise 2 — FASTQ quality plot](exercise-02-fastq-quality-plot.py)** — compute and plot per-base mean quality across a real FASTQ. (~45 min)
3. **[Exercise 3 — Filter by quality](exercise-03-filter-by-quality.py)** — trim and filter a FASTQ on length and mean quality. (~40 min)

## How to work the exercises

- Install Biopython 1.83 before starting: `python -m pip install biopython==1.83 matplotlib`.
- **Type the code yourself.** Copy-paste defeats the muscle-memory point.
- Each file is runnable: `python exercise-XX.py`. All assertions must pass.
- If you get stuck for more than 10 minutes, peek at the inline hints.
- No solutions are checked in. The exercises are self-checking.
