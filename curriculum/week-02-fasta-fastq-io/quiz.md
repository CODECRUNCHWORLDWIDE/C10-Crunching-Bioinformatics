# Week 2 — Quiz

Ten multiple-choice questions on FASTA, FASTQ, Phred quality, and Biopython IO. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** Which line of a FASTQ record carries the per-base quality string?

- A) Line 1 (the `@` header).
- B) Line 2 (the sequence).
- C) Line 3 (the `+` separator).
- D) Line 4.

---

**Q2.** A FASTQ quality character is `?`. ASCII for `?` is 63. Assuming Phred+33 encoding, what Phred quality does this character represent?

- A) Q63.
- B) Q30.
- C) Q33.
- D) Q96.

---

**Q3.** Which Biopython function should you reach for when you have a FASTA file that contains **exactly one record** and you want a hard error if it contains more or fewer?

- A) `Bio.SeqIO.parse`
- B) `Bio.SeqIO.read`
- C) `Bio.SeqIO.to_dict`
- D) `Bio.SeqIO.index`

---

**Q4.** You inspect a FASTQ file and the lowest ASCII character anywhere in the quality strings is `B` (ASCII 66). The most likely encoding is:

- A) Sanger / Phred+33.
- B) Illumina 1.8+ / Phred+33.
- C) Solexa+64 or Illumina 1.3+/1.5+ / Phred+64.
- D) Unencoded — the file is malformed.

---

**Q5.** A Phred quality score of Q40 corresponds to an error probability of:

- A) 1 in 10.
- B) 1 in 100.
- C) 1 in 1,000.
- D) 1 in 10,000.

---

**Q6.** Why does the implicit FASTA rule "all lines after the header until the next `>` line are sequence" matter operationally?

- A) It lets you `cat` two FASTA files together and get a valid FASTA.
- B) It is the reason FASTA records cannot be parsed in parallel.
- C) It is what enables compression of FASTA below 50% of plain-text size.
- D) It is the reason FASTA cannot represent protein sequences.

---

**Q7.** You read a 30x-coverage human-genome FASTQ file (35 GB gzipped) and want the mean per-base quality at each read position. The correct approach is:

- A) `records = list(SeqIO.parse(file, "fastq"))` then iterate.
- B) Decompress to disk, load into a pandas DataFrame, groupby position.
- C) Iterate with `SeqIO.parse` directly over `gzip.open(file, "rt")`, accumulating per-position sums and counts, never holding more than one record.
- D) Load into Biopython's `SeqIO.to_dict` and index by read id.

---

**Q8.** A `SeqRecord` parsed from a FASTQ file stores the per-base Phred qualities at which attribute?

- A) `record.qualities`
- B) `record.seq.qualities`
- C) `record.letter_annotations["phred_quality"]`
- D) `record.annotations["quality"]`

---

**Q9.** Which of the following is the **correct** way to trim a `SeqRecord` to its first 50 bases while keeping the per-base quality string in sync?

- A) `record.seq = record.seq[:50]`
- B) `trimmed = record[:50]`
- C) `record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][:50]`
- D) `record.seq = record.seq[:50]; record.qual = record.qual[:50]`

---

**Q10.** Which statement about FastQC's pass / warn / fail icons is **most accurate**?

- A) A "fail" icon means the data is unusable and must be discarded.
- B) The icons are calibrated for standard Illumina whole-genome shotgun and can produce expected "fails" on RNA-seq, amplicon, or small-RNA libraries.
- C) The icons are derived from Phred quality scores alone.
- D) A "pass" icon means the library is publication-ready.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **D** — Line 4 is the quality string. Lines 1, 2, 3, 4 are header, sequence, `+` separator, quality. A robust parser counts lines in groups of four; it never splits on a delimiter, because `@` and `+` are both legitimate quality characters.

2. **B** — Phred+33: `Q = ASCII - 33 = 63 - 33 = 30`. Q30 corresponds to a 1-in-1,000 error rate and is the canonical Illumina quality target.

3. **B** — `SeqIO.read` asserts exactly one record. `SeqIO.parse` is a generator (will not fail on multi-record files), `to_dict` and `index` are for keyed access. The fail-loud helper is the right tool when "single record" is a precondition.

4. **C** — Phred+33 quality characters range from `!` (33) upward. Anything in the `@`–`B` region (64–66) cannot be a Phred+33 Q0–Q2 (which would be `!`–`#`), so the file is almost certainly Phred+64. The Illumina 1.5+ quirk explicitly marks low-confidence bases with `B`.

5. **D** — `Q = -10 log10(P)`, so Q40 → P = 10^(-4) = 0.0001 = 1 in 10,000. Memorize Q20=1%, Q30=0.1%, Q40=0.01%.

6. **A** — The "all lines until next `>` are sequence" rule is what lets you concatenate FASTA files without any glue logic. It is also why grep-style processing of FASTA is so easy. (Parallel parsing is possible by scanning for `^>` markers; compression isn't structurally enabled by the rule; protein FASTA is unaffected.)

7. **C** — Stream. `SeqIO.parse` over a gzip handle is bounded-memory regardless of file size. The other options either load the whole file into RAM (A, D) or waste disk and add an unnecessary pandas dependency (B).

8. **C** — `record.letter_annotations["phred_quality"]`. Biopython stores per-letter info (quality, secondary structure, etc.) in this dict for any record where the annotation applies per-base. The list length is always equal to `len(record.seq)`.

9. **B** — Slice the whole `SeqRecord` with `record[:50]`. This is the only option that keeps `seq` and `letter_annotations` in lockstep. Option A leaves the quality unchanged (and Biopython will refuse to write the record). Option C is dangerous — you must update both halves. Option D uses a `.qual` attribute that doesn't exist on a Biopython `SeqRecord`.

10. **B** — The icons are heuristics calibrated for standard Illumina whole-genome shotgun. Common "fails" on RNA-seq (random-hexamer priming bias at positions 1–15), amplicon (intentional duplication), and small-RNA (intentionally short reads) libraries are *expected*, not red flags. Read the modules; don't read the icons in isolation.

</details>

---

If you scored under 7, re-read the lectures for the questions you missed — especially anything involving Phred encoding or `SeqRecord` slicing. If you scored 9 or 10, you are ready to move to the [homework](./homework.md).
