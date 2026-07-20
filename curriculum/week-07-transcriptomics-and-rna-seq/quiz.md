# Week 7 — Quiz

Ten multiple-choice questions on RNA-seq biology, fastp trimming, the k-mer pseudoalignment idea, the EM step, TPM/CPM normalization, and the alignment-vs-pseudoalignment trade-off. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** A short-read aligner like `bwa mem` does not work for RNA-seq on the human genome because:

- A) The reads are too short.
- B) Many RNA-seq reads span exon-exon junctions, and `bwa mem` cannot align a read whose two halves are separated by a kilobase of intron on the genome.
- C) The reference is too large.
- D) `bwa mem` does not support paired-end reads.

---

**Q2.** The default k-mer size in `kallisto index` is:

- A) 21.
- B) 25.
- C) 31.
- D) 51.

---

**Q3.** kallisto's "compatibility class" for a read is:

- A) The single transcript the read aligns to.
- B) The set of transcripts that share at least one k-mer with the read, intersected across all k-mers of the read.
- C) The GTF annotation features that overlap the read.
- D) The list of FASTQ files the read came from.

---

**Q4.** TPM and CPM differ in that:

- A) TPM is computed per kilobase, CPM is per million.
- B) CPM corrects for library size only; TPM corrects for both library size and effective transcript length.
- C) They are the same thing under different names.
- D) TPM is for paired-end data, CPM is for single-end.

---

**Q5.** The identity `sum_g(TPM_g) = 10^6` holds:

- A) Only on samples with library size exactly 10^6.
- B) Only after dropping low-count transcripts.
- C) Always, by construction.
- D) Only for kallisto, not for Salmon.

---

**Q6.** A 1,000 bp gene with 1,000 reads in a sample of library size 10^6 reads has:

- A) CPM = 1,000, TPM ≈ 1,250 (with eff_length 800).
- B) CPM = 1, TPM ≈ 1.25.
- C) CPM = 1,000, RPKM = 1.
- D) Cannot be computed without the fragment length.

---

**Q7.** A `fastp` HTML report shows `pct_retained = 65%`, `pct_q30 = 78%`, `duplication_rate = 0.89`, `adapter_pct = 82%`. The most likely diagnosis is:

- A) The sample is fine; modern QC tools are conservative.
- B) The sample is degraded or PCR over-amplified; short inserts caused adapter readthrough, duplication is from PCR, and the run should be repeated.
- C) The wrong reference was used.
- D) `fastp` is misconfigured.

---

**Q8.** Salmon's `--validateMappings` flag enables:

- A) Re-validating the FASTQ headers.
- B) A selective-alignment step that score-checks each pseudoaligned read against its candidate transcripts and drops poor matches.
- C) Checking that the index is fresh.
- D) Re-running kallisto on the same data to confirm Salmon's output.

---

**Q9.** For differential expression with DESeq2 or edgeR, the right input is:

- A) TPM-normalized counts.
- B) CPM-normalized counts.
- C) RPKM-normalized counts.
- D) Raw integer counts (rounded if the source produced fractional values).

---

**Q10.** kallisto/Salmon are typically much faster than STAR/HISAT2 + featureCounts because:

- A) They use more cores.
- B) They skip the per-base alignment step entirely and use k-mer compatibility classes instead.
- C) They run on the GPU.
- D) They use a smaller reference.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — Splicing. A read that spans an exon-exon junction has its two halves separated on the genome by an intron that is not in the read. `bwa mem` cannot infer this jump; STAR and HISAT2 can. The alternative, pseudoalignment to a transcriptome FASTA, sidesteps the problem because the transcriptome already has the splicing resolved. Lecture 1 §2, Lecture 2 §1.

2. **C** — k=31. Bray et al. 2016 §Online Methods establishes this empirically: smaller k (21) gives too many chance matches, larger k (51) is too sensitive to single-nucleotide errors. k=31 is the default in kallisto, Salmon, and most other k-mer-based bioinformatics tools. Lecture 2 §2.

3. **B** — Intersection of color sets. For each k-mer in a read, kallisto looks up the set of transcripts that k-mer is known to appear in (the "color set"), and intersects across all k-mers of the read. The result is the compatibility class — the set of transcripts the read could have come from. The read is never aligned base-by-base. Lecture 2 §2.

4. **B** — CPM corrects for library size only; TPM also corrects for effective transcript length. The TPM definition is `TPM_g = 10^6 × (c_g / l_g) / sum_g'(c_g' / l_g')`. The TPM denominator is the per-sample sum of per-gene rates, which is what makes TPM sum to 10^6 per sample. Lecture 3 §2, §4.

5. **C** — Always, by construction. The proof: `sum_g(TPM_g) = sum_g(10^6 × r_g / S) = (10^6/S) × sum_g(r_g) = (10^6/S) × S = 10^6`. The identity holds for both kallisto and Salmon, on any sample, regardless of filtering. Lecture 3 §4.

6. **A** — CPM = 10^6 × 1000 / 10^6 = 1,000. TPM ≈ 10^6 × (1000/800) / sum_g(c_g/eff_length_g) where the sum is dominated by the rest of the library; if we assume the rest of the library averages roughly the same rate (1.25 per gene), the per-sample sum is ~800 and TPM(this gene) ≈ 1.25 × 10^6 / 800 ≈ 1,562 — within an order of magnitude of the "1,250" estimate in option A. The exact number depends on the rest of the sample, but the algebra gives a value in the 10^3 range. C is wrong because RPKM = 10^9 × 1000 / (10^6 × 1000) = 1, not 1.0 by way of CPM = 1,000. Lecture 3 §2 and §4.

7. **B** — Degraded or over-PCR'd sample. All four numbers are out of range: pct_retained < 75% (Lecture 1 §7), pct_q30 < 85%, duplication > 0.80, adapter_pct > 70% (short inserts mean lots of adapter readthrough). Together these indicate a sample that has lost complexity to PCR over-amplification on short fragments. Re-running with fewer PCR cycles and longer inserts would help; with this sample as-is, downstream quantification is partially compromised. Lecture 1 §4, §5.

8. **B** — Selective alignment score check. After pseudoalignment identifies a candidate set of transcripts, Salmon does a fast Smith-Waterman-like alignment score check between the read and each candidate. Transcripts where the actual alignment score is below threshold are dropped from the read's compatibility class. This recovers ~1-2% accuracy at the cost of ~2x runtime. Default in Salmon ≥ 1.0. Lecture 2 §4.1.

9. **D** — Raw integer counts. DESeq2 and edgeR apply their own size-factor / TMM normalization internally and want integer counts on input. Feeding them TPM (already normalized) double-normalizes and produces wrong p-values. Round any fractional counts (from kallisto/Salmon's EM) before passing in. Lecture 3 §1, §11.

10. **B** — They skip per-base alignment. Pseudoalignment is O(read_length) hash lookups; alignment is O(read_length × log(genome)) Smith-Waterman per read. For 100 bp reads, the per-read cost is ~7 μs vs ~100-1,000 μs. For a 3 M-read sample, ~30 seconds vs ~30 minutes. Lecture 2 §2.

</details>

---

If you scored under 7, re-read Lecture 1 §2-4 (RNA-seq biology + fastp) and Lecture 2 §2-3 (k-mer compatibility + EM). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
