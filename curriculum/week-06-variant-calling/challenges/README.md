# Week 6 — Challenges

One challenge this week. It is the cross-tool comparison problem that separates "I can run `bcftools call`" from "I can defend a variant-set difference between two production callers to a reviewer." Budget ~90 minutes.

## Index

1. **[Challenge 1 — Compare bcftools and GATK on the same BAM](challenge-01-compare-callers.md)** — call variants on the Week 5 mini-project BAM with both `bcftools call -m` and GATK `HaplotypeCaller`, normalize both VCFs with `bcftools norm`, compute the intersection and private sets with `bcftools isec`, and explain the disagreements at the per-variant level. Identify at least two specific variants where the callers disagree and explain why. (~90 min)

## How to work the challenge

- Read the prompt in full before writing any code. Sketch the data flow on paper: BAM → bcftools VCF, BAM → GATK VCF, both VCFs → `bcftools norm` → `bcftools isec` → comparison.
- **Use the work from Exercise 1.** The bcftools side of the comparison is the same pipeline you ran in Exercise 1, applied to the mini-project's larger BAM instead of the lambda toy BAM.
- The challenge intentionally avoids tuning the callers' parameters. Run both with their defaults. The point is to understand where defaults disagree, not to engineer agreement.
- Be honest about discrepancies. The two callers will disagree on ~5% of SNPs and ~15% of indels on a well-behaved BAM. List the cases in your writeup; do not pretend the disagreements are a bug.
- The point of the challenge is to internalize that **variant calling is model-dependent**. The Bayesian column-by-column model in `bcftools` and the local-reassembly model in GATK produce different calls in different corner cases. Once you have seen 10 specific examples of how, you will never again ask "which caller is right" without first asking "which caller has the right model for this region."
