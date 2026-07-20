# Week 4 — Quiz

Ten multiple-choice questions on BLAST mechanics, E-values, the BLAST family, output formats, and the Biopython API. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** The BLAST seed-and-extend heuristic has four steps. Which of the following lists them in the correct order?

- A) Score → seed → extend → report.
- B) Extend → seed → score → report.
- C) Seed → extend → score → report.
- D) Seed → score → extend → report.

---

**Q2.** For default `blastp` against a protein database with BLOSUM62, the word size `W` and neighborhood threshold `T` are:

- A) `W = 11`, `T = 0` (exact 11-mer required).
- B) `W = 3`, `T = 11` (any 3-mer with similarity score ≥ 11 against the query 3-mer is a seed).
- C) `W = 28`, `T = 22` (megablast-style, exact 28-mer required).
- D) `W = 1`, `T = 4` (single-residue seeding).

---

**Q3.** The Karlin-Altschul E-value formula is `E = K · m · n · exp(-λ · S)`. Which symbol does *not* depend on the database being searched?

- A) `n`.
- B) `m`.
- C) `K`.
- D) `λ`.

---

**Q4.** A BLAST hit has E-value `4e-87`. Which interpretation is correct?

- A) There is a 4×10⁻⁸⁷ probability the hit is biologically meaningful.
- B) Under the K-A null (random i.i.d. residues), the expected number of hits at this score or better in a database this size is 4×10⁻⁸⁷ — effectively zero, so the observed hit is not chance.
- C) The hit covers 4×10⁻⁸⁷ percent of the query.
- D) The hit took 4×10⁻⁸⁷ seconds to find.

---

**Q5.** You have an unannotated 250 kb DNA contig from a metagenomics assembly. You want to identify which proteins it encodes. The right BLAST family member is:

- A) `blastn` against `nt`.
- B) `blastp` against `nr`.
- C) `blastx` against `nr` (translates the query in all 6 frames, searches a protein database).
- D) `tblastn` against `nt`.

---

**Q6.** The two-hit rule (BLAST 2.0) requires:

- A) The query to have at least two amino acid types represented.
- B) Two non-overlapping seed words within distance `A` on the same diagonal before triggering an extension.
- C) Two complete passes over the database, one forward and one reverse-complement.
- D) Two separate E-value cutoffs, one for the seed phase and one for the extend phase.

---

**Q7.** Which BLAST output format is best suited for parsing into a pandas DataFrame?

- A) `-outfmt 0` (human-readable pairwise alignments).
- B) `-outfmt 5` (XML — readable, but verbose; better with `Bio.Blast.NCBIXML`).
- C) `-outfmt 6` (tab-separated, twelve standard columns, optionally with custom fields).
- D) `-outfmt 11` (SAM/BAM).

---

**Q8.** The 16S rRNA percent-identity threshold conventionally used to discriminate **same species** is approximately:

- A) ~ 80%.
- B) ~ 90%.
- C) ~ 95%.
- D) ~ 98.7% (Stackebrandt & Ebers 2006).

---

**Q9.** Which of the following Biopython calls correctly submits a `blastn` query to NCBI against the `16S_ribosomal_RNA` database?

- A) `Bio.Blast.NCBIWWW.qblast("blastn", "16S_ribosomal_RNA", str(query.seq), expect=1e-50)`.
- B) `Bio.Blast.run("blastn", db="16S_ribosomal_RNA", query=query, evalue=1e-50)`.
- C) `Bio.Blast.NCBIWWW.run(program="blastn", database="16S_ribosomal_RNA", sequence=query)`.
- D) `Bio.Entrez.efetch(db="blast", program="blastn", query=str(query.seq))`.

---

**Q10.** Which of the following is **not** a known failure mode of a BLAST-based top-hit taxonomy classifier?

- A) Low-complexity regions (e.g., `AAAAAAAAA...` runs) matching unrelated database entries.
- B) A chimeric reference database entry, splicing together sequences from two organisms.
- C) Paralog/ortholog confusion — a paralog from a different species with higher similarity than the true ortholog.
- D) E-value cutoff being applied at the reporting phase rather than the search phase, biasing toward shorter alignments.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **C** — Seed, extend, score, report. Seed finds candidate matches via short k-mer hits; extend grows them outward via the X-drop rule; score runs a gapped alignment on the surviving HSPs; report ranks by E-value and prints the survivors above the cutoff.

2. **B** — `W = 3` and `T = 11`. The neighborhood threshold means BLAST does *not* require the seed 3-mer to match the query 3-mer exactly; any 3-mer with BLOSUM62 similarity ≥ 11 to the query 3-mer is a seed. This is why protein BLAST is more sensitive than nucleotide BLAST at the seed step. (A) and (C) describe `blastn`-style word sizes — wrong family. (D) is too low to be useful.

3. **B** — `m` (the query length) does not depend on the database. `n` is the (effective) database length and obviously depends on the database. `K` and `λ` depend on the substitution matrix *and* on the residue composition of the database, so they shift slightly when you change databases. The cleanest answer is that `m` is the only one that is purely a property of the query.

4. **B** — The E-value is the expected number of chance hits at this score or better under the K-A null, *not* a posterior probability of the hit being real. Choice (A) is the most common misinterpretation in introductory courses and gets people fired in industry — calibrate now.

5. **C** — `blastx`. The query is DNA (so you need to translate something), the target is a protein database, so translate the query in all six frames and search against `nr`. `tblastn` is the mirror image (protein query, translated nucleotide DB) and is appropriate when you have a protein and want to find where it is encoded in an unannotated genome.

6. **B** — Two non-overlapping seeds within distance `A` on the same diagonal. The 1997 BLAST 2.0 paper introduced this to cut the false-positive extension rate by ~50x at minimal sensitivity cost.

7. **C** — `-outfmt 6` is the tabular format with twelve fixed columns (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore). `pandas.read_csv(..., sep="\t", names=[...])` handles it directly. `-outfmt 5` (XML) is also fine but takes more code. `-outfmt 0` is for humans; `-outfmt 11` is not a BLAST output format (it does not exist).

8. **D** — ~98.7%. Stackebrandt & Ebers 2006 proposed this threshold for the 16S gene; recent literature debates 97–99%. Below 98.7% you are typically below species; above 98.7% you are within species. Note this is a 16S-specific threshold — for whole-genome ANI the species threshold is ~95%.

9. **A** — `Bio.Blast.NCBIWWW.qblast(program, database, sequence, ...)` is the canonical call. The positional argument order is `(program, database, sequence)`. (B), (C), and (D) reference nonexistent APIs or wrong call signatures.

10. **D** — The E-value cutoff in BLAST is applied at the *reporting* phase: BLAST does the full search and then filters to hits below the user's cutoff. This is by design and not a failure mode. (A), (B), and (C) are all real failure modes covered in Lecture 1 §6 and §4.1, and mitigated by techniques covered in the mini-project: `dustmasker`/`seg` for low complexity, database curation for chimeras, and orthology-aware methods like OrthoFinder for paralog/ortholog distinction.

</details>

---

If you scored under 7, re-read Lecture 1 for the algorithm and statistics questions and Lecture 2 for the BLAST family and Biopython API questions. If you scored 9 or 10, you are ready to start the [homework](./homework.md).
