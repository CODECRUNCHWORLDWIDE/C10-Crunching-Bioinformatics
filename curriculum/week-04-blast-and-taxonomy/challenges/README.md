# Week 4 — Challenges

One challenge this week. It is the design problem that separates "I can call `blastn` and parse the output" from "I can build a classifier and defend its choices on a labelled test set." Budget ~120 minutes.

## Index

1. **[Challenge 1 — Taxonomy classifier](challenge-01-taxonomy-classifier.md)** — implement two competing classifier strategies on the parsed BLAST output (top-hit and majority-among-top-N), evaluate them on a small labelled dataset, and write up which one wins and why. (~120 min)

## How to work the challenge

- Read the prompt in full before writing any code. Sketch the data flow on paper: query → BLAST hits → taxonomy lookups → classification → metric.
- **Use the work from Exercises 2 and 3.** The challenge is the next layer on top, not a redo of the exercises.
- Pick your evaluation metric *before* you compute results. Precision, recall, and accuracy can all be computed on a multi-class classifier; the right one to report depends on whether you care about coverage or correctness.
- Be honest about failure modes. If your classifier confidently mislabels a query, the writeup should name the failure mode (low complexity? paralog? chimeric reference?) and propose a mitigation. The mini-project will revisit each failure mode at scale.
