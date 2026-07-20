# Week 5 — Challenges

One challenge this week. It is the data-archeology problem that separates "I can run `samtools markdup`" from "I understand what `samtools markdup` is doing and can defend my duplicate rate to a reviewer." Budget ~90 minutes.

## Index

1. **[Challenge 1 — Detect duplicates by hand](challenge-01-detect-duplicates.md)** — implement a positional duplicate detector from scratch in Python: for each read pair in a sorted BAM, group by (5' position, strand, mate 5' position, mate strand), and call any group of size ≥ 2 duplicates. Compare your duplicate rate to `samtools markdup`'s answer on the same BAM. Explain any discrepancies. (~90 min)

## How to work the challenge

- Read the prompt in full before writing any code. Sketch the data flow on paper: BAM → reads → (5' position keys) → grouping → duplicate calls.
- **Use the work from Exercises 1 and 3.** The challenge is the next layer on top, not a redo of the exercises. You should have a sorted+indexed BAM from Exercise 1 and the pysam patterns from Exercise 3.
- The challenge intentionally avoids re-implementing optical-duplicate detection (which requires the read names to be parsed for flow-cell coordinates). Stick to PCR-duplicate detection from alignment positions and orientations.
- Be honest about discrepancies. `samtools markdup` handles a few edge cases your hand-rolled detector probably will not — supplementary alignments, hard-clipping adjustments to the 5' position, the unmapped mate case. List them in your writeup; do not pretend your answer should match exactly.
- The point of the challenge is to internalize that **duplicate detection is positional, not sequence-based**. Once you have written 30 lines of Python to detect duplicates yourself, you will never again confuse "duplicate" with "identical read sequence."
