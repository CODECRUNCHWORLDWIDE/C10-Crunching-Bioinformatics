# Week 10 — Quiz

> **Reproducibility note.** This quiz tests your knowledge of the long-read sequencing chemistries, the OLC assembly paradigm, the polishing tools, and the assembly QC metrics. Knowing the mechanics is the difference between a defensible assembly and a FASTA-shaped picture; pin every parameter and the result is reproducible.

Ten multiple-choice questions on ONT R10.4 chemistry, PacBio HiFi, Dorado, Flye, Canu, Hifiasm, Medaka, BUSCO, N50 / L50, and Bandage. Take it with the lecture notes closed. Aim for 9/10 before the mini-project. Answer key at the bottom — do not peek.

---

**Q1.** Oxford Nanopore R10.4.1 sequencing measures:

- A) The fluorescence emission of a single nucleotide incorporated by a polymerase in a zero-mode waveguide.
- B) The ionic current through a protein pore as a single strand of DNA translocates through it; the signal is called the squiggle and is recorded at ~4 kHz.
- C) The mass-spectrometry signature of each base after enzymatic digestion.
- D) The light scattering off fluorescent dye-labelled bases attached to a glass slide.

---

**Q2.** PacBio HiFi reads achieve ~99.8% per-base accuracy because:

- A) The PacBio polymerase has a higher fidelity than other sequencing polymerases.
- B) The PacBio instrument uses a four-colour fluorescent imaging system at low temperature.
- C) The HiFi protocol records multiple polymerase passes over the same circular template and computes a per-base consensus across passes (the CCS algorithm); errors that are not shared across passes are corrected.
- D) The HiFi protocol uses electrochemical detection with an Au surface.

---

**Q3.** The canonical Flye call for modern Dorado-SUP-basecalled R10.4.1 reads is:

- A) `flye --pacbio-hifi reads.fastq --genome-size 1m --out-dir out --threads 4`
- B) `flye --nano-raw reads.fastq --genome-size 1m --out-dir out --threads 4`
- C) `flye --nano-hq reads.fastq --genome-size 1m --out-dir out --threads 4`
- D) `flye --auto reads.fastq --out-dir out`

---

**Q4.** Flye's repeat graph differs from Canu's best-overlap graph in that:

- A) Flye uses a de Bruijn graph; Canu uses an overlap graph.
- B) Flye represents each repeated region as a single node tagged as a repeat, with the non-repeat flanks as edges entering and exiting; Canu keeps each read as a node and lets the best-pair-of-overlaps speak for itself.
- C) Flye is unable to handle repetitive genomes; Canu can.
- D) Flye and Canu produce identical graphs on identical inputs.

---

**Q5.** Hifiasm is the preferred assembler for PacBio HiFi input because:

- A) It has a faster I/O implementation than Flye.
- B) HiFi reads are already ~99.8% accurate, so Hifiasm skips the error-correction stage that Flye and Canu spend most of their compute on; it builds the assembly graph directly from raw HiFi reads and can emit haplotype-resolved primary and alternate contig sets in a single pass.
- C) Hifiasm includes a built-in BUSCO step.
- D) Hifiasm is the only free assembler that handles HiFi input.

---

**Q6.** Medaka polishes an ONT-derived assembly. The most common failure mode of Medaka is:

- A) The neural network runs out of memory on contigs longer than 5 Mb.
- B) Medaka and the basecaller use different output FASTA conventions, so the polished file has the wrong sequence IDs.
- C) The user passes a Medaka model that does not match the basecaller model that produced the reads (e.g. `r941_min_sup_g507` applied to R10.4.1 reads); the resulting polish is *worse* than no polish because the model is trained on a different signal distribution.
- D) Medaka does not work on circular contigs.

---

**Q7.** For a contig set with lengths (in bp, sorted descending) 600,000 / 400,000 / 200,000 / 100,000, the N50 and L50 are:

- A) N50 = 600,000; L50 = 1.
- B) N50 = 400,000; L50 = 2.
- C) N50 = 300,000; L50 = 2.5.
- D) N50 = 200,000; L50 = 3.

---

**Q8.** A BUSCO short_summary line `C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124` on a *bacteria_odb10* run means:

- A) The assembly is 99.2% identical to a reference at the nucleotide level.
- B) Out of 124 single-copy ortholog gene models in `bacteria_odb10`, 99.2% are present as complete-length genes (99.2% as single-copy, 0.0% as duplicated), 0.0% are fragmented, and 0.8% (one gene) is missing.
- C) BUSCO is 99.2% confident in the lineage assignment.
- D) The assembly has 124 contigs of which 99.2% are at least 1 kb long.

---

**Q9.** A clean bacterial assembly graph viewed in Bandage looks like:

- A) A tree with many leaves.
- B) A grid pattern of segments.
- C) A single self-loop (one segment with one link that returns to itself); for a closed circular chromosome.
- D) A scatter of disconnected dots.

---

**Q10.** When you ship a long-read assembly to a collaborator, the most important reproducibility metadata to include in the run-info JSON is:

- A) The size of the input FASTQ in bytes.
- B) The Python version.
- C) The basecaller name and version (Dorado 0.7.2), the basecaller model (`dna_r10.4.1_e8.2_400bps_sup@v4.3.0`), the read filter parameters (min length, min QV), the assembler and version (Flye 2.9.5), the assembler input mode (`--nano-hq`), the genome-size hint, the thread count, the polish tool and model (Medaka 1.12.0 / `r1041_e82_400bps_sup_v4.3.0`), the BUSCO version and lineage dataset name + creation date, and the run date.
- D) A copy of the assembly FASTA pasted into the email body.

---

## Answer key

<details>
<summary>Click to reveal answers</summary>

1. **B** — Ionic current through a protein pore; the squiggle. Lecture 1 §3. The R10.4.1 pore is a double-pore design that improves the homopolymer signal; the basecaller is a neural network that maps the time series of currents to a base sequence.

2. **C** — Multiple passes over the same circular template; CCS consensus. Lecture 1 §4. The HiFi protocol's accuracy comes from the consensus, not from a higher-fidelity polymerase or a better detector; with 5+ passes the error rate drops below 1%; with 10+ passes below 0.1%.

3. **C** — `--nano-hq` for Dorado-SUP-basecalled R10.4.1. Lecture 2 §4. `--nano-raw` is for legacy Guppy-basecalled R9.4 reads; using the wrong input mode flag loosens or tightens the overlap thresholds in ways that produce defensible-looking but worse assemblies.

4. **B** — Flye's repeat graph has explicit repeat nodes; Canu's best-overlap graph keeps every read as a node and edges only the best overlap on each end. Lecture 2 §3. Both are OLC; the difference is how the graph data structure encodes repeats.

5. **B** — HiFi reads are already accurate; Hifiasm skips correction and emits haplotype-resolved contigs. Lecture 2 §6. The wall-clock advantage on small genomes is ~5x; on diploid large genomes the haplotype-resolution advantage is decisive.

6. **C** — Model mismatch produces a polish worse than no polish. Lecture 3 §3. The Medaka model is trained on a specific basecaller's error distribution; applying it to reads from a different basecaller asks the model to correct errors that are not there.

7. **B** — Total = 1,300,000; half = 650,000. Cumulative pass: 600,000 (under), 600,000 + 400,000 = 1,000,000 (over). So N50 = 400,000, L50 = 2. Lecture 3 §5.

8. **B** — 99.2% complete out of 124 orthologs; 0.0% duplicated; 0.0% fragmented; 0.8% missing. Lecture 3 §6. The C number is single-copy + duplicated; the F number is fragments; the M number is missing.

9. **C** — Single self-loop for a clean circular bacterial chromosome. Lecture 3 §7. A linear path is also clean (linear chromosome); a tangle or many bubbles is a flag.

10. **C** — Everything: basecaller, model, filter, assembler, mode, polish model, BUSCO lineage + date, seed, run date. Resources style guide; Lecture 3 §10. Without these the same FASTQ can produce different assemblies on different days and the difference is impossible to debug. The run-info JSON is the most important artefact in any assembly output directory.

</details>

---

If you scored under 7, re-read Lecture 1 §2-§4 (the two platforms), Lecture 2 §3-§6 (the three assemblers), and Lecture 3 §3 and §10 (Medaka, failure modes). If you scored 9 or 10, you are ready to start the [homework](./homework.md).
