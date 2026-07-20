# Week 8 — Quiz

> **Educational and research use only.** This quiz tests your knowledge of the mechanics of variant annotation and the ACMG/AMP framework. Knowing the mechanics is not the same as making a clinical interpretation. The disclaimer that opens every Week 8 file applies here.

Ten multiple-choice questions on VEP, SnpEff, gnomAD, ClinVar, the ACMG/AMP criteria, SIFT, PolyPhen-2, and PharmGKB. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** The three orthogonal axes of variant interpretation, in the Week 8 framework, are:

- A) Reference, alternate, quality.
- B) Functional consequence (VEP/SnpEff), population frequency (gnomAD), and clinical knowledge (ClinVar/dbSNP).
- C) Single-nucleotide variant, insertion, deletion.
- D) Pathogenic, benign, uncertain.

---

**Q2.** The VEP `CSQ` INFO field is structured as:

- A) One key-value pair per variant.
- B) A single comma-separated list of consequence strings.
- C) Pipe-separated fields per transcript record; multiple transcript records per variant are separated by commas; the field order is documented in the `##INFO=<ID=CSQ` header line.
- D) JSON embedded in the INFO column.

---

**Q3.** The ACMG **BA1** criterion, in the 2015 framework, fires when:

- A) The variant is absent from all databases.
- B) gnomAD popmax allele frequency exceeds 5% in any general population.
- C) SIFT and PolyPhen agree that the variant is deleterious.
- D) The variant is a frameshift.

---

**Q4.** The ACMG **PVS1** criterion is the only "Pathogenic Very Strong" criterion and requires:

- A) A null variant (nonsense, frameshift, canonical +/-1,2 splice site) AND the gene must have loss-of-function as a known disease mechanism.
- B) A missense variant in any disease gene.
- C) A variant absent from gnomAD.
- D) A ClinVar Pathogenic assertion with 4-star review status.

---

**Q5.** A SIFT score of 0.02 with a PolyPhen-2 score of 0.94 most likely indicates:

- A) The variant is benign.
- B) SIFT and PolyPhen disagree; manual review needed.
- C) Both predictors agree the variant is deleterious; this supports the ACMG PP3 criterion.
- D) The annotation failed.

---

**Q6.** A variant in ClinVar with `CLNREVSTAT=criteria_provided,_conflicting_classifications` should be:

- A) Treated as Pathogenic (the most extreme classification submitted).
- B) Treated as Benign (the least extreme classification submitted).
- C) Flagged as conflicting; the per-submission classifications should be reviewed by a human; no automated pipeline should resolve the conflict.
- D) Dropped from the report entirely.

---

**Q7.** The gnomAD "popmax" allele frequency is:

- A) The global allele frequency across all sub-populations.
- B) The maximum allele frequency across non-bottlenecked sub-populations (Finnish, Ashkenazi Jewish, Other, and Middle Eastern are excluded from popmax by gnomAD convention).
- C) The allele frequency in the largest sub-population only.
- D) The same as the dbSNP frequency.

---

**Q8.** Of the 28 ACMG/AMP 2015 criteria, the ones that are cleanly mechanically computable from a VCF plus VEP/gnomAD/ClinVar/SIFT/PolyPhen are approximately:

- A) All 28.
- B) About 8: PVS1 (with a gene list), PM2, PM4, PP3, BA1, BS2, BP4, BP7. The other ~20 require functional studies, pedigree segregation, or patient phenotype.
- C) Only PVS1 and PM2.
- D) None; clinical judgment is always required.

---

**Q9.** Pharmacogenomics differs from disease-variant interpretation in that:

- A) The unit of analysis is the star allele (a defined haplotype), not the individual variant; the output is a metabolizer phenotype plus a per-drug recommendation; the canonical evidence framework is CPIC, not ACMG.
- B) PharmGKB is a paid service.
- C) Pharmacogenomic variants are always pathogenic.
- D) The two are identical except for the disease.

---

**Q10.** When generating a variant annotation report, the most important reproducibility metadata to record is:

- A) The Python version.
- B) The Linux kernel version.
- C) The exact versions of every annotation database queried (VEP cache version, SnpEff database, ClinVar release date, gnomAD version) plus the run date and the assembly.
- D) The disk size.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — The three axes are functional consequence (VEP, SnpEff), population frequency (gnomAD), and clinical knowledge (ClinVar, dbSNP). Each axis has a canonical free database. Reporting only one or two of the three produces incomplete interpretations and is the most common Week 8 error. Lecture 1 §2.

2. **C** — Pipe-separated fields per transcript record; comma-separated across records; field order from the header. The field order is `Allele|Consequence|IMPACT|SYMBOL|...|SIFT|PolyPhen|...`. To parse: split on `,` then on `|`, then zip with the header field names. Lecture 1 §3, §8.

3. **B** — Stand-alone benign at popmax > 5%. The 5% threshold reflects the empirical observation that a variant common in healthy populations is essentially incompatible with a rare Mendelian-disease causal interpretation. BA1 overrides any other evidence (a variant cannot be both Pathogenic and BA1). The 2019 SVI refinements suggested context-specific thresholds (e.g. 0.02 for some specific disorders), but 0.05 is still the default. Lecture 3 §2, §3.

4. **A** — Null variant in a gene where LOF is a known disease mechanism. The "gene where LOF is a known mechanism" qualifier is essential — many genes are haplosufficient and losing one copy is tolerated. The ClinGen Haploinsufficiency Tier 3 list is the canonical reference. The 2023 SVI update (Walsh et al. 2023) refined PVS1 with NMD-rule considerations and exon location, but the core criterion is unchanged. Lecture 3 §2, §3.

5. **C** — Both predictors agree deleterious; this supports PP3. SIFT < 0.05 is "deleterious"; PolyPhen > 0.85 is "probably damaging." Together, two lines of computational evidence pointing the same way is the ACMG PP3 criterion. Note that PP3 is one of the *supporting* criteria, not a strong one; you need multiple PP criteria or one or more PS/PM criteria to reach Likely Pathogenic or Pathogenic. Lecture 3 §3, §5.

6. **C** — Flag and review. ~3% of ClinVar records have conflicting interpretations. The ClinVar review-status field for conflicts (`criteria_provided,_conflicting_classifications`) explicitly tells you "do not resolve this automatically." A pipeline that picks one submission and reports it as the consensus is overclaiming. The Week 8 reports always flag conflicts as a separate category. Lecture 2 §2.

7. **B** — Max AF across non-bottlenecked sub-populations. Finnish, Ashkenazi Jewish, Other, and Middle Eastern are excluded from popmax because they are bottlenecked or small. The popmax is what most ACMG-flavored filters pivot on: PM2 (popmax < 0.0001) and BA1 (popmax > 0.05). For an individual interpretation you should also check the patient's ancestry and the matched sub-population AF. Lecture 2 §1.

8. **B** — About 8 of 28 are mechanically computable. The full list: PVS1 (with a curated gene list), PM2 (gnomAD popmax), PM4 (VEP consequence), PP3 (SIFT + PolyPhen agreement), BA1 (gnomAD popmax > 0.05), BS2 (gnomAD homozygous count), BP4 (SIFT + PolyPhen agreement of "tolerated"), BP7 (synonymous without splice impact). A few more (PS1, PM5, BS1, BP3) are computable with additional curated tables. The remaining ~14-17 criteria require functional studies, pedigree segregation, or patient phenotype and cannot be mechanized. Lecture 3 §2, §3.

9. **A** — Star alleles, metabolizer phenotype, CPIC framework. Pharmacogenomic variants act on drug-metabolism enzymes; the relevant unit is the haplotype, not the variant; the output is a phenotype (ultrarapid / normal / intermediate / poor metabolizer) plus the CPIC-published per-drug recommendation. PharmGKB is free. Lecture 3 §6.

10. **C** — Database versions, run date, assembly. Without these, the same VCF run a year later will produce different annotations and the difference cannot be debugged. The run-info JSON should include the VEP cache version, the SnpEff database version, the ClinVar release date, the gnomAD version, and the date of the run. Python and kernel versions are less important; the database versions are critical. Lecture 1 §6, Lecture 2 §5, Resources style guide.

</details>

---

If you scored under 7, re-read Lecture 1 §2 (the three axes) and Lecture 3 §2-3 (the ACMG criteria and the mechanically computable subset). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
