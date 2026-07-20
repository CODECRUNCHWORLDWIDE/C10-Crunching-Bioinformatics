# Week 6 — Quiz

Ten multiple-choice questions on the genotype-likelihood model, the `bcftools mpileup` + `call` pipeline, the VCF format, the GATK Best Practices hard filters, VEP annotation, and the `bcftools` toolchain. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** The `bcftools call -m` model computes the per-genotype posterior probability `P(genotype | reads)` using:

- A) Counting the number of reads carrying each allele and dividing by the total.
- B) The binomial likelihood `P(reads | genotype)` times a prior `P(genotype)` from Hardy-Weinberg or a coalescent model, then normalizing.
- C) A neural network trained on the 1000 Genomes Project data.
- D) A non-parametric bootstrap over the reads.

---

**Q2.** Which of the following is **not** an `INFO` field emitted by `bcftools mpileup` with `-a 'AD,DP,SP'`?

- A) `DP` — total read depth at the position.
- B) `AD` — per-allele read depths.
- C) `SP` — Phred-scaled strand-bias p-value.
- D) `SOR` — symmetric strand-odds-ratio.

---

**Q3.** For variant calling on a single *E. coli* sample with `bcftools call`, the correct ploidy flag is:

- A) `--ploidy 1` (haploid).
- B) `--ploidy 2` (the default, diploid).
- C) `--ploidy 4` (some bacteria are polyploid).
- D) ploidy is auto-detected and the flag does not need to be set.

---

**Q4.** A VCF record reads `NC_000913.3 150123 . A G 227 . DP=42;MQ=60 GT:PL:DP:AD 1:255,0:42:0,42`. This is:

- A) A heterozygous SNP in a diploid sample.
- B) A homozygous-alt SNP in a haploid sample, depth 42, all 42 reads carry the alt allele.
- C) A 1-bp deletion.
- D) A multiallelic site with three possible genotypes.

---

**Q5.** The CIGAR-equivalent for VCF — the operation that left-aligns indels and splits multiallelic records — is performed by:

- A) `bcftools sort`.
- B) `bcftools norm -f ref.fa -m -any`.
- C) `bcftools filter -e 'INDEL'`.
- D) `bcftools index`.

---

**Q6.** A variant has `INFO/QD = 1.2`. According to the GATK Best Practices hard filters for SNPs:

- A) PASS (QD < 2 is fine).
- B) Filtered out (QD < 2 is below the QD2 threshold).
- C) Cannot tell without knowing the indel-or-SNP status.
- D) QD is not in the GATK hard-filter recipe.

---

**Q7.** The GATK Best Practices hard filter for indel strand bias (`FS`) uses a threshold of **200**, while for SNPs it is **60**. Why?

- A) Indels are intrinsically more strand-biased; soft-clipping at indel breakpoints preferentially affects one strand.
- B) GATK has a bug; the indel threshold should be 60.
- C) Indels are usually called from a single strand, so any strand-bias signal is meaningless.
- D) The thresholds are arbitrary; both could be 100.

---

**Q8.** Ensembl VEP annotates a variant `chr17:7674220 G>A` with consequence `missense_variant` on transcript `ENST00000269305` (TP53). The HGVSp string is `p.Arg175His`. This means:

- A) The amino acid at position 175 of the TP53 protein changes from arginine to histidine.
- B) The DNA at position 175 changes from arginine codon to histidine codon (this is nonsensical — DNA is not amino acids).
- C) The variant is at chromosome 17, position 175.
- D) The IMPACT category is HIGH.

---

**Q9.** Which `pysam` call correctly opens a bgzipped, tabix-indexed VCF and iterates over variants in a region?

- A) `pysam.AlignmentFile("calls.vcf.gz", "rb").fetch("chr1", 1000, 2000)`.
- B) `pysam.VariantFile("calls.vcf.gz").fetch("chr1", 1000, 2000)`.
- C) `Bio.SeqIO.parse("calls.vcf.gz", "vcf")`.
- D) `pysam.read_vcf("calls.vcf.gz", contig="chr1", start=1000, stop=2000)`.

---

**Q10.** Which of the following is **not** a well-known failure mode of single-sample variant calling?

- A) Low-coverage positions (< 5x) producing false negatives because the likelihood model has too few observations.
- B) Unmarked PCR-duplicate stacks producing false positives because the same fragment is counted as independent observations.
- C) MAPQ-0 multimappers in repetitive regions producing false positives at every paralogous copy.
- D) The `bcftools call -m` model rejecting variants in genes longer than 10 kb.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — Bayesian inference: per-genotype binomial likelihood times prior, normalize to get posterior. The "count and divide" approach (A) is what happens at infinite depth; with finite depth, you need the likelihood model to weight observations by base quality and account for sampling noise. Lecture 1 §2.

2. **D** — `SOR` is a GATK metric, not a `bcftools` one. `bcftools mpileup -a 'AD,DP,SP'` emits `AD`, `DP`, and `SP`. GATK emits `SOR` (and `FS`, `MQRankSum`, `ReadPosRankSum`, `QD`) but `bcftools` does not. This is why the bacterial mini-project hard-filter recipe uses `INFO/SP` rather than `SOR`. Lecture 1 §3, Lecture 2 §4.3.

3. **A** — `--ploidy 1`. *E. coli* is haploid; the default `--ploidy 2` would call every variant as heterozygous (`GT=0/1`), which is biologically nonsensical. Lecture 1 §4, §8.

4. **B** — Haploid homozygous-alt: `GT=1` is the haploid genotype "the alt allele," `AD=0,42` is "0 reads reference, 42 reads alt," `DP=42` is total depth. A diploid heterozygous would be `GT=0/1` and `AD≈21,21`. Lecture 1 §5.2.

5. **B** — `bcftools norm -f ref.fa -m -any`. `-f ref.fa` enables left-alignment of indels; `-m -any` splits multiallelic records. Without normalization, two VCFs of the same biological event will appear to disagree because of representation differences. Lecture 2 §3.

6. **B** — Filtered out. `QD < 2` is the `QD2` filter in the Best Practices recipe; any SNP (or indel) with `QD < 2` is removed. The threshold is the same for both, even though the other thresholds (`FS`, `SOR`) differ. Lecture 2 §4.

7. **A** — Indels are intrinsically more strand-biased. The aligner soft-clips reads at indel breakpoints, and the soft-clipping is asymmetric across strands; the result is that real indels look strand-biased even when they are not artifactual. The looser threshold (200 vs 60) accounts for this. Lecture 2 §4.2.

8. **A** — `p.Arg175His` reads "at protein position 175, arginine changes to histidine." This is the HGVS protein notation. Codons are 3 bases; the underlying coding-sequence (`HGVSc`) for this is `c.524G>A`. Lecture 2 §5.3.

9. **B** — `pysam.VariantFile("calls.vcf.gz").fetch("chr1", 1000, 2000)`. The class is `VariantFile`, not `AlignmentFile` (that is for BAM). `Bio.SeqIO` does not handle VCF. Lecture 1 §7, Resources.

10. **D** — There is no such failure mode. The likelihood model in `bcftools call -m` operates per position; the length of the surrounding gene is irrelevant to the model. A, B, and C are real failure modes covered in Lecture 1 §1 and Lecture 2 §8.

</details>

---

If you scored under 7, re-read Lecture 1 §2-4 (genotype model + bcftools pipeline) and Lecture 2 §1-5 (GATK, hard filters, VEP). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
