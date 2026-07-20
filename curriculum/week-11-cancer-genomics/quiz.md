# Week 11 — Quiz

> **Educational and research use only.** This quiz tests your knowledge of the matched-pair model, Mutect2, FilterMutectCalls, mutational signatures, and the COSMIC / OncoKB / CIViC interpretation layer. Nothing here should be applied to clinical care without an accredited laboratory pipeline.

Ten multiple-choice questions on tumor-normal somatic calling, Mutect2 mechanics, filter sets, contamination, mutational signatures, and clinical-interpretation databases. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom; do not peek.

---

**Q1.** The matched tumor-normal model of somatic variant calling assumes:

- A) The tumor and normal samples were sequenced on different platforms to maximize signal diversity.
- B) The tumor and normal samples come from the same patient; the normal is genuinely non-tumor tissue; the tumor is enriched for tumor cells; both samples were processed equivalently. Somatic variants are present in tumor but absent (or at very low AF) in normal.
- C) The normal sample is sequenced first and the tumor sample is sequenced only if specific candidate variants are identified.
- D) The model produces a list of variants present in the normal but not in the tumor.

---

**Q2.** A tumor sample with 40% purity carries a heterozygous truncal somatic variant. The expected observed allele frequency in the BAM at 50x coverage is approximately:

- A) 50% (the variant is heterozygous in the tumor cells).
- B) 40% (the purity is 40%).
- C) 20% (0.5 × purity = 0.5 × 0.4 = 0.2).
- D) 100% (somatic variants are always at AF 1.0).

---

**Q3.** Mutect2 in tumor-normal mode is invoked with `-tumor TUMOR_X -normal NORMAL_Y`. The names `TUMOR_X` and `NORMAL_Y` must match:

- A) The filenames of the BAM files.
- B) The values of the BAM `@RG SM:` headers of the tumor and normal BAMs respectively.
- C) The patient's name in the EHR.
- D) The reference build name.

---

**Q4.** The Panel of Normals (PON) for Mutect2 is:

- A) A list of patient names whose samples have been sequenced.
- B) A multi-sample VCF distilled from many technically-similar normal samples; variants present in the PON at non-trivial frequency are flagged as recurrent technical artifacts and Mutect2 down-weights or filters them.
- C) A copy of gnomAD restricted to common variants.
- D) An optional output file that records normal-sample calls.

---

**Q5.** `CalculateContamination` followed by `FilterMutectCalls --contamination-table` does what:

- A) Computes the tumor purity and adjusts the allele frequencies to the pure-tumor value.
- B) Computes the cross-sample contamination fraction (the fraction of DNA from a different individual) and flags variants whose AF is consistent with the contamination as the `contamination` filter reason.
- C) Removes all low-AF variants below a fixed threshold.
- D) Identifies copy-number variants in the tumor.

---

**Q6.** For a somatic SNV at chromosome chr22, position 23,456,789, with REF=G and ALT=A in a reference where the 3-mer at positions 23,456,788-23,456,790 is `CGT`, the pyrimidine-normalized 96-class label is:

- A) `C[G>A]T`.
- B) `A[C>T]G` (after reverse-complementing the context and complementing the alleles, since REF is a purine).
- C) `T[G>A]C`.
- D) `G[A>T]C`.

---

**Q7.** A SigProfilerAssignment decomposition reports SBS1 (35%), SBS5 (40%), SBS3 (15%), SBS39 (8%), residual (2%) with a reconstructed-spectrum cosine similarity of 0.79. The correct interpretation is:

- A) The decomposition is clean; report SBS1 / SBS5 / SBS3 as the top three signatures with high confidence.
- B) The decomposition is poor (cosine < 0.85) and SBS3 / SBS39 are degenerate; report the result with a flag for the cosine, name the SBS3 / SBS39 ambiguity, and consider that the spectrum may be too noisy for a stable fit.
- C) The decomposition is exclusively explained by SBS5 because it has the highest contribution.
- D) The decomposition has failed completely; discard the result.

---

**Q8.** OncoKB FDA Evidence Level 1 means:

- A) The variant has any biological evidence published in any tumor type.
- B) FDA-recognized biomarker predictive of response to an FDA-approved drug in this specific tumor type.
- C) The variant is found at frequency 1% in COSMIC.
- D) The variant has been mentioned in one CIViC evidence item.

---

**Q9.** The standard Mutect2 → FilterMutectCalls pipeline emits a filtered VCF where most variants are not PASS. The largest non-PASS filter category on a typical real-data run is usually:

- A) `contamination`.
- B) `slippage`.
- C) `germline` (the matched normal is doing its job: most candidate somatic variants turn out to be inherited germline).
- D) `strand_bias`.

---

**Q10.** Two Mutect2 runs on the same tumor-normal BAM pair, same reference, same PON, and same germline-resource produce VCFs that differ at 30% of low-AF variants. The most likely cause is:

- A) Mutect2 is non-deterministic by design.
- B) One of the runs used a different reference build (GRCh37 vs GRCh38), a different gnomAD release, a different PON, or a different GATK version. Mutect2 itself is deterministic given identical inputs.
- C) The variants are randomly sampled.
- D) Mutect2 emits a different sub-set on every run.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — Matched-pair model. Lecture 1 §1-§3. The model assumes equal processing of both samples; the somatic call is the set difference.

2. **C** — `observed_AF = purity * tumor_cell_AF + (1 - purity) * normal_AF`. For purity 0.4, tumor_cell_AF 0.5, normal_AF 0: observed_AF = 0.5 * 0.4 = 0.2. Lecture 1 §4.

3. **B** — The values of the BAM `@RG SM:` headers. Lecture 2 §3. Mis-tagging the names produces silent inverted calls; the exercise wrapper verifies them via `samtools view -H` before calling Mutect2.

4. **B** — Multi-sample VCF of recurrent normal-sample technical artifacts. Lecture 1 §9; Lecture 2 §2. Platform-specific; the Broad publishes a free GRCh38 PON.

5. **B** — Cross-sample contamination estimation and filter. Lecture 1 §5; Lecture 2 §6. Note: purity and contamination are different concepts; contamination is corrected by `--contamination-table`, purity is not.

6. **B** — `A[C>T]G`. The original REF `G` is a purine, so we reverse-complement the context `CGT` to `ACG`; the alleles `G>A` complement to `C>T`. The result is `A[C>T]G`. Lecture 3 §2; Exercise 3 helper `normalize_to_pyrimidine`.

7. **B** — Poor cosine and known degeneracy. Lecture 3 §7-§8. The cosine of 0.79 is below the 0.85 threshold; SBS3 and SBS39 are known to be co-linear. Report the ambiguity, do not collapse it.

8. **B** — FDA-recognized biomarker predictive of response in this tumor type. Lecture 3 §9 (the OncoKB section). The OncoKB levels are public at <https://www.oncokb.org/levels>.

9. **C** — `germline`. Lecture 2 §5 and §10. On real-data runs, hundreds to thousands of low-AF candidates turn out to be inherited variants that the matched normal also carries.

10. **B** — One of the upstream inputs differed. Lecture 2 §13-§14. Mutect2 itself is deterministic. The most common silent input differences are reference-build mismatch, gnomAD-release version, PON release, and GATK version.

</details>

---

If you scored under 7, re-read Lecture 1 §3-§5 (matched normal, purity, contamination), Lecture 2 §2-§3 (Mutect2 internals, sample-name verification), Lecture 2 §5-§6 (filter set and contamination flow), and Lecture 3 §2 and §7-§9 (pyrimidine normalization, the signature degeneracy, and the OncoKB level definitions). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
