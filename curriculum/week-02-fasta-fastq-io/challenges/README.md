# Week 2 — Challenges

One challenge this week. It is the IO problem that separates "I parsed a small file" from "I can process production data." Budget 60–90 minutes.

## Index

1. **[Challenge 1 — Streaming a large FASTA](challenge-01-streaming-large-fasta.md)** — parse a multi-gigabyte gzipped FASTA without loading it into RAM, then compute three summary statistics in a single pass. (~75 min)

## How to work the challenge

- Read the prompt fully before writing any code. Pseudocode the data flow on paper.
- **No loading the whole file.** That is the entire point of the exercise.
- Measure your peak RSS (`/usr/bin/time -v` on Linux, Activity Monitor on macOS, or `tracemalloc` in pure Python). A correct solution stays under ~150 MB regardless of the input file size.
- When you finish, compare your solution against `seqkit stats` on the same input.
