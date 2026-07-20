# Week 11 — Homework

> **Educational and research use only.** All problems below use simulated or anonymized public data. None of the outputs are clinical results. Do not apply the methods to identifiable patient data.

Six practice problems for the week. Each problem builds on a single Week 11 lecture or exercise. Aim to complete at least four; the mini-project assumes Problems 1, 2, and 3 are done.

---

## Problem 1 — Sample-name verification on simulated mis-tagged BAMs

**Source:** Lecture 2 §3; Exercise 1.

**Goal:** demonstrate that the sample-name verification step catches a mis-tagged BAM pair before Mutect2 runs.

**Steps:**

1. Copy `data/tumor_chr22.bam` and `data/normal_chr22.bam` to `homework/p01/`.
2. Use `samtools addreplacerg -r 'ID:SWAP\tSM:TUMOR_SAMPLE'` to *re-tag* the normal BAM with the tumor's sample name — simulating a mis-tag at library prep.
3. Run the Exercise 1 script on the mis-tagged pair. It should error out at the sample-name verification step with `ValueError: Tumor and normal samples share the same name 'TUMOR_SAMPLE'`.
4. Untag and re-run to confirm the script PASSes the verification when the BAMs are correct.

**Deliverable:** `homework/p01/notes.md` with the command sequences, the error message you observed, and a one-sentence explanation of why this catch matters (silent inverted calls).

**Time estimate:** 30 minutes.

---

## Problem 2 — Inspect the FILTER tally on a contaminated-vs-clean pair

**Source:** Lecture 2 §5-§6; Exercise 2.

**Goal:** observe how the FilterMutectCalls FILTER distribution changes between a clean and a contaminated tumor sample.

**Steps:**

1. Run Exercise 2 with the standard didactic pair. Record the FILTER tally and the contamination estimate in `homework/p02/clean.json`.
2. Use `samtools merge` to combine 1% of an unrelated patient's BAM into the tumor BAM (simulated contamination). The script `data/scripts/make_contaminated.sh` performs this if you do not want to do it by hand.
3. Re-run Exercise 2 on the contaminated pair. Record the FILTER tally and the contamination estimate in `homework/p02/contaminated.json`.
4. Compare: how did the `contamination` filter count change? How did the contamination fraction estimate change? Did the PASS variant count drop?

**Deliverable:** `homework/p02/notes.md` with the two filter tables side by side, the two contamination estimates, and a paragraph explaining what changed and why.

**Time estimate:** 60 minutes.

---

## Problem 3 — Verify the canonical 96-class order

**Source:** Lecture 3 §2; Exercise 3.

**Goal:** verify your understanding of the 96-class Alexandrov-lab ordering by hand.

**Steps:**

1. Without looking at the Exercise 3 source, write a Python function `canonical_96_classes() -> list[str]` that emits the 96 class labels in the canonical order. The first six entries should be `A[C>A]A, A[C>A]C, A[C>A]G, A[C>A]T, C[C>A]A, C[C>A]C`. The last entry should be `T[T>G]T`.
2. Verify by importing the Exercise 3 helper and comparing the two lists with `assert mine == reference`.
3. Write a function `class_index(cls: str) -> int` that returns the 0-based position of a class label in the canonical order. Verify that `class_index('A[C>A]A') == 0` and `class_index('T[T>G]T') == 95`.

**Deliverable:** `homework/p03/p03.py` with the two functions and the assertions.

**Time estimate:** 45 minutes.

---

## Problem 4 — Build a single PON from three control BAMs

**Source:** Lecture 1 §9; Lecture 2 §2.

**Goal:** build a small panel-of-normals from three control BAMs and observe how it affects FilterMutectCalls.

**Steps:**

1. The didactic data folder ships three "control normal" BAMs at `data/control_normal_1.bam`, `data/control_normal_2.bam`, `data/control_normal_3.bam`. Each is a small chr22-subset BAM from a non-cancer donor.
2. Run `gatk Mutect2 --tumor` on each control BAM in turn (no matched normal); this produces three single-sample VCFs of "candidate somatic" calls. On normal-tissue BAMs, these are actually a list of recurrent technical artifacts and rare germline variants.
3. Combine the three VCFs into a PON with `gatk CreateSomaticPanelOfNormals --vcfs vcf1 --vcfs vcf2 --vcfs vcf3 -O custom_pon.vcf.gz`.
4. Re-run Exercise 1 with `--pon custom_pon.vcf.gz` instead of the Broad-public PON; record the PASS count.
5. Compare with the Broad-public-PON run from Exercise 1.

**Deliverable:** `homework/p04/notes.md` with the PASS counts under each PON, the size of the custom PON, and a paragraph on what the comparison tells you about PON specificity.

**Time estimate:** 90 minutes.

---

## Problem 5 — Reproduce a published smoking signature

**Source:** Lecture 3 §4 and §8.

**Goal:** verify that a tobacco-smoking signature (SBS4) is recovered cleanly on a known smoking-associated lung cancer sample.

**Steps:**

1. Download one of the PCAWG sample profiles publicly distributed by COSMIC for a tobacco-associated lung cancer (the COSMIC signature website at <https://cancer.sanger.ac.uk/signatures/sbs/> has profile data; or use the synapse-distributed PCAWG sample `SP116816` if you have Synapse access).
2. Build the 96-class spectrum from the sample's PASS VCF using the Exercise 3 helper.
3. Run SigProfilerAssignment against COSMIC v3.3 SBS catalog.
4. Verify that SBS4 (tobacco smoking) is in the top three signatures and has fractional contribution >= 0.3. If not, investigate why (low mutation count? Reference-build mismatch?).

**Deliverable:** `homework/p05/notes.md` with the input sample identifier, the top three signatures with fractional contributions, the cosine similarity, and a one-paragraph interpretation. Cite COSMIC v3.3 and Alexandrov et al. 2020 (PMC 7054213).

**Time estimate:** 60 minutes.

---

## Problem 6 — Look up TP53 R175H in three databases

**Source:** Lecture 3 §9-§11; Challenge 2.

**Goal:** practice the variant-level interpretation workflow on a well-characterized hotspot mutation.

**Steps:**

1. TP53 R175H is the most-frequently observed TP53 mutation in pan-cancer. Look it up in COSMIC at <https://cancer.sanger.ac.uk/cosmic/search?q=TP53+R175H> and record: the number of tumor samples it has been observed in, the top three tumor types, the COSMIC Cancer Gene Census classification of TP53 (oncogene / tumor suppressor / both).
2. Look it up in OncoKB at <https://www.oncokb.org/gene/TP53/R175H>. Record: the Mutation Effect, the Oncogenic classification, the highest FDA evidence level, and any linked drugs (with their development status).
3. Look it up in CIViC at <https://civicdb.org/genes/2/summary>. Record: the number of evidence items, the highest star-rating, the most-cited paper.
4. Combine into a one-page Markdown summary.

**Deliverable:** `homework/p06/tp53_r175h.md` with the three database entries and a final paragraph that summarizes the variant's clinical-actionability profile in 200-300 words. Note that R175H is a *loss-of-function* mutation in a tumor suppressor, so the "actionable therapy" status is more nuanced than for an oncogene gain-of-function hotspot like KRAS G12C.

**Time estimate:** 60 minutes.

---

## Optional problem — Compare per-AF-bin agreement between Mutect2 and Strelka2

**Source:** Challenge 1.

**Goal:** quantify how the caller-agreement depends on tumor allele frequency.

**Steps:**

1. Complete Challenge 1 first.
2. Bin the union of Mutect2 PASS and Strelka2 PASS variants by tumor AF: `[0, 0.05), [0.05, 0.1), [0.1, 0.2), [0.2, 0.3), [0.3, 1.0]`.
3. For each bin, compute the Jaccard index restricted to that bin.
4. Plot or table the Jaccard-vs-AF curve.

**Deliverable:** `homework/optional/jaccard_by_af.md` with the table and a one-paragraph interpretation. Expected: Jaccard rises with AF; at AF >= 30% it should be > 0.95.

**Time estimate:** 90 minutes.

---

## Submission checklist

- [ ] Each problem's `homework/pNN/` directory contains the source files and a `notes.md` write-up.
- [ ] The write-ups are in lab-notebook voice: pin tool versions, name parameters, report specific numbers.
- [ ] No problem references identifiable patient data; all inputs are simulated or public.
- [ ] The `homework/` directory is committed to your portfolio repo before the mini-project deadline.

## A note on the source of the simulated control BAMs

The "control normal" BAMs in `data/control_normal_*.bam` are synthetic; each was simulated from the GRCh38 chr22 reference with a small number of random heterozygous variants at gnomAD-frequency-distributed sites. They are not derived from any real patient sample. The PON built from them (Problem 4) is a pedagogical toy; a production PON should be built from real normal samples sequenced on the target platform.
