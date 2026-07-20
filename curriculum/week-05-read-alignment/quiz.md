# Week 5 — Quiz

Ten multiple-choice questions on the FM-index, the BWA-MEM and minimap2 alignment pipelines, SAM/BAM format, the FLAG bit field, the CIGAR alphabet, MAPQ, and the `samtools` toolchain. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** The FM-index over a reference of length `n` answers "how many times does a length-`m` pattern occur in the reference?" in time:

- A) `O(n)`.
- B) `O(m · n)`.
- C) `O(m)`, independent of `n`.
- D) `O(log n)`.

---

**Q2.** Which of the following is **not** true about the Burrows-Wheeler Transform?

- A) It is a reversible permutation of the input string.
- B) It produces a string with the same multiset of characters as the input.
- C) It compresses well because similar characters tend to cluster together.
- D) It computes Smith-Waterman in `O(n log n)` time.

---

**Q3.** For aligning 150 bp paired-end Illumina reads against the *E. coli* K-12 MG1655 reference (`NC_000913.3`), the most appropriate aligner choice is:

- A) `bwa aln` + `bwa sampe` (legacy backtrack mode).
- B) `bwa mem` or `minimap2 -ax sr`.
- C) `minimap2 -ax map-ont`.
- D) Smith-Waterman against every position.

---

**Q4.** The SAM `FLAG = 99` decodes to which combination of bits?

- A) `0x80 + 0x10 + 0x2 + 0x1` = READ2 + REVERSE + PROPER_PAIR + PAIRED.
- B) `0x40 + 0x20 + 0x2 + 0x1` = READ1 + MREVERSE + PROPER_PAIR + PAIRED.
- C) `0x400 + 0x40 + 0x2 + 0x1` = DUP + READ1 + PROPER_PAIR + PAIRED.
- D) `0x4 + 0x10 + 0x40 + 0x80` = UNMAP + REVERSE + READ1 + READ2.

---

**Q5.** The CIGAR string `120M2I20M3S` decodes to a read with implied query length and reference span of:

- A) Query 145, ref 142.
- B) Query 142, ref 145.
- C) Query 145, ref 140.
- D) Query 142, ref 140.

---

**Q6.** A read has `MAPQ = 0` in a `bwa mem` BAM. The correct interpretation is:

- A) The read did not align — it is unmapped.
- B) The alignment has a 100% probability of being wrong.
- C) The read aligned equally well to ≥ 2 reference positions; the aligner cannot pick one.
- D) The read had Phred quality 0 at every position.

---

**Q7.** The canonical samtools pipeline for marking PCR duplicates from a freshly-aligned BAM is:

- A) `samtools markdup` directly on the coordinate-sorted BAM.
- B) `samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup`.
- C) `samtools dedup | samtools sort | samtools index`.
- D) `samtools view -F 0x400 | samtools sort | samtools index`.

---

**Q8.** What does `samtools depth -a aln.bam` do that `samtools depth aln.bam` does not?

- A) Outputs alignment-quality information.
- B) Includes reference positions with zero read coverage.
- C) Includes secondary alignments in the depth calculation.
- D) Auto-indexes the BAM if no `.bai` is present.

---

**Q9.** Which Biopython / pysam call correctly reads a sorted, indexed BAM and iterates over reads overlapping a region?

- A) `pysam.AlignmentFile("aln.bam", "rb").fetch("chr1", 1000, 2000)`.
- B) `Bio.SeqIO.parse("aln.bam", "bam")`.
- C) `pysam.AlignmentFile("aln.bam").read_region(1000, 2000)`.
- D) `pysam.read("aln.bam", contig="chr1", start=1000, stop=2000)`.

---

**Q10.** Which of the following is **not** a well-known failure mode of short-read alignment?

- A) Low-complexity reference regions (microsatellites, homopolymer runs) producing high multimapper rates and MAPQ = 0 reads.
- B) Large structural variants (insertions / deletions > 50 bp) exceeding the aligner's affine-gap penalty and producing soft-clipped reads.
- C) Adapter contamination producing reads that fail to align or align to unexpected references.
- D) BWA-MEM's MAPQ cap of 60 causing false-positive variant calls in homozygous regions.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **C** — `O(m)`, independent of reference length. This is the FM-index's defining property and the reason short-read aligners scale to 3 Gb references. The seed step is `O(m)` per read; the extension step adds linear-in-read-length cost. Total per-read cost is `O(m)`, total run cost is `O(read_count · m)`.

2. **D** — The BWT does not compute Smith-Waterman. It is a *transform* of the reference that enables efficient indexing (when paired with rank/select structures into an FM-index); the alignment itself is still Smith-Waterman on the windowed extension. A, B, and C are all true.

3. **B** — Either `bwa mem` or `minimap2 -ax sr` works. Both are tuned for 150 bp paired-end Illumina against bacterial-to-mammalian references. `bwa aln` is for legacy reads < 70 bp. `minimap2 -ax map-ont` is for Oxford Nanopore long reads. Smith-Waterman against every position is the brute-force baseline that BWA replaces.

4. **B** — FLAG 99 = `0x40 + 0x20 + 0x2 + 0x1` = READ1 + MREVERSE + PROPER_PAIR + PAIRED. This is the standard FLAG for the forward read of a healthy paired-end pair. The reverse mate carries FLAG 147 = `0x80 + 0x10 + 0x2 + 0x1` = READ2 + REVERSE + PROPER_PAIR + PAIRED.

5. **C** — Query 145, ref 140. Breaking down `120M2I20M3S`: 120 M's consume both query and reference (+120 to each), 2 I's consume query only (+2 query), 20 M's add +20 to each, 3 S's consume query only (+3 query). Query total = 120 + 2 + 20 + 3 = 145. Reference total = 120 + 0 + 20 + 0 = 140. M, I, S, =, X consume query; M, D, N, =, X consume reference.

6. **C** — MAPQ 0 means the read aligned equally well to ≥ 2 positions. The aligner declines to pick one. The read *is* aligned (check the `0x4` flag, not MAPQ, for unmapped status); the *placement* is ambiguous. Variant callers typically discard MAPQ-0 reads because the position is uncertain.

7. **B** — `sort -n | fixmate -m | sort | markdup`. `samtools markdup` requires that mate-coordinate fields (`PNEXT`, `MS`, `ms`) are populated, which `samtools fixmate -m` does. `fixmate` requires name-sorted input; `markdup` requires coordinate-sorted input. Hence the four-step idiom.

8. **B** — `-a` includes zero-coverage positions in the output. Without `-a`, positions with no reads are silently omitted, biasing any depth distribution computed from the output upward. Always use `-a` for QC.

9. **A** — `pysam.AlignmentFile("aln.bam", "rb").fetch("chr1", 1000, 2000)` is the canonical pattern. The `"rb"` mode opens BAM (binary). `Bio.SeqIO` does not handle BAM (it handles FASTA, FASTQ, GenBank, etc.); the other API surfaces are imaginary.

10. **D** — BWA-MEM's MAPQ cap of 60 is by design, not a failure mode. It means "P_misalignment ≤ 10^-6"; the cap exists because going higher would imply unrealistic precision given finite read information. It does not cause false-positive variant calls — the cap is at the top of the scale, not the bottom. A, B, and C are real failure modes covered in Lecture 1 §8.

</details>

---

If you scored under 7, re-read Lecture 1 §2–3 (FM-index and BWA-MEM) and Lecture 2 §1–4 (SAM, FLAG, CIGAR, MAPQ). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
