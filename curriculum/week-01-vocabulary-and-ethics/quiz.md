# Week 1 — Quiz

Ten multiple-choice questions, mixed biology and ethics. Take it with your lecture notes closed. Aim for 9/10 before moving to Week 2. Answer key at the bottom — do not peek.

---

**Q1.** How many distinct codons are there in the standard genetic code, and how many standard amino acids do they encode?

- A) 20 codons, 64 amino acids
- B) 64 codons, 20 amino acids
- C) 32 codons, 20 amino acids
- D) 64 codons, 64 amino acids

---

**Q2.** Which of the following best describes the difference between a **gene** and a **transcript**?

- A) They are synonyms.
- B) A gene is a region of DNA; a transcript is a specific RNA produced from a gene, and one gene can produce multiple transcripts.
- C) A transcript is the DNA copy of a gene; a gene is the RNA copy of a transcript.
- D) A gene is a single chromosome; a transcript is a region within a chromosome.

---

**Q3.** Which file format is used to store **aligned sequencing reads** in a binary, compressed form?

- A) FASTA
- B) FASTQ
- C) BAM
- D) VCF

---

**Q4.** You read in a press release: *"Scientists have discovered the gene for musical talent."* What is the most accurate critique of this sentence?

- A) The press release should have used "transcript" instead of "gene."
- B) Complex behavioural traits are polygenic and gene-environment interactive; framing it as a single "gene for X" is essentially always wrong.
- C) The press release should have specified a reference genome version.
- D) The sentence is correct; modern genomics has identified many "the gene for X" associations.

---

**Q5.** Gymrek et al. (2013) showed that public 1000 Genomes data could be partially re-identified by:

- A) Comparing variants against publicly listed missing-persons databases.
- B) Inferring surnames from Y-chromosome short tandem repeats cross-referenced with public genealogy databases.
- C) Recovering names embedded in the FASTQ read headers.
- D) Matching variants against social-media profile photos.

---

**Q6.** Which is the **correct** way to cite a tool in a methods section?

- A) "We used Biopython."
- B) "We used the latest Biopython."
- C) "We used Biopython 1.83 with Python 3.11.5."
- D) "We used a bioinformatics library."

---

**Q7.** What is the relationship between a **chromosome** and a **genome**?

- A) A chromosome is larger than a genome.
- B) A genome is the entire DNA content of an organism; a chromosome is one of the large DNA molecules the genome is partitioned into.
- C) They are different names for the same thing.
- D) A chromosome is a subset of a gene; a genome is a chromosome of chromosomes.

---

**Q8.** Your friend, a Ph.D. learner, offers you a small de-identified VCF file from their lab for "just practice." Their lab has not given public consent for learner use. According to the C10 data-ethics rules, what should you do?

- A) Accept the file — it is de-identified, so the donors are protected.
- B) Accept the file but commit it to a private repo, not a public one.
- C) Decline and use a public consent-cleared dataset like 1000 Genomes instead.
- D) Accept the file but only analyze two samples, not all of them.

---

**Q9.** Which sequence is the correct **reverse complement** of `ATGCGT`?

- A) `TACGCA`
- B) `ACGCAT`
- C) `TGCGTA`
- D) `GCATGC`

---

**Q10.** Which of the following statements about reference genomes is **most accurate**?

- A) Every member of a species has a genome identical to the reference.
- B) The reference is the "best" genome of a species; everyone else has errors.
- C) The reference is a community-curated coordinate system; individuals differ from it at millions of positions, which is normal.
- D) There is only one reference genome per species, and it never changes.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — 4 nucleotides taken 3 at a time gives 4³ = 64 codons; these encode 20 amino acids plus a stop signal. The redundancy is the source of the "wobble" at the third position.
2. **B** — Genes are loci on DNA; transcripts are RNA copies of genes. Alternative splicing means one gene routinely produces multiple isoform transcripts.
3. **C** — BAM is the binary, compressed form of SAM. FASTA is unaligned sequences; FASTQ adds quality scores; VCF stores variants.
4. **B** — "The gene for X" framing implies single-gene determinism that almost never applies to complex traits. This is the language pattern C10 trains you to recognize and avoid.
5. **B** — The paper inferred surnames from Y-STRs cross-referenced with public genealogy databases. It became a foundational re-identification result.
6. **C** — Versions are cited; "latest" is meaningless six months later. Tool name, version, and language version are the minimum.
7. **B** — The human genome is the totality of DNA; it is partitioned into 23 chromosome pairs (24 if you count X and Y separately). Conflating the two is one of the most common vocabulary mistakes.
8. **C** — C10 Rule 1: use public consent-cleared datasets for coursework. The friend's offer is well-intended but does not have the right consent provenance for course use. A public dataset has clean provenance.
9. **B** — Complement of ATGCGT is TACGCA; reverse it for ACGCAT. (Or: reverse first to TGCGTA, then complement to ACGCAT — the operations commute.)
10. **C** — The reference is a coordinate system, not a "best" genome. Individuals carry millions of variants relative to GRCh38; that is normal human variation, not error.

</details>

---

If you scored under 7, re-read the lectures for the questions you missed — especially any ethics or vocabulary questions. If you scored 9 or 10, you are ready to start the [homework](./homework.md).
