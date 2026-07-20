# Week 10 Homework

> **Reproducibility note.** Every homework problem produces a file the grader can read alongside a `run-info.json`. The numbers below are illustrative; your numbers will differ slightly with different Flye, Medaka, or BUSCO versions. Pin the versions, pin the seed, and commit both the result file and the run-info.

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-10/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Coverage sweep on the 1 Mb reference

**Problem statement.** Re-run the Exercise 1 pipeline (badread + Flye) at three coverage levels: 20x, 50x, and 100x. Compare the resulting assemblies on (a) contig count, (b) N50, (c) total assembled length, (d) QV against the reference. Plot the four metrics as a function of coverage.

Answer in `homework/notes/p1-coverage-sweep.md`:

1. At which coverage level does the assembly first become a single contig?
2. Does the QV continue to improve with coverage beyond 50x?
3. At what coverage level does the assembly's total length first reach >= 99% of the reference length?
4. If you had to pick one coverage level for a real bacterial-genome project, what would you pick and why?

**Acceptance criteria.**

- `homework/p1_coverage_sweep.py` runs end to end (orchestrates three Flye runs at 20x, 50x, 100x; one minimap2 + QV calculation per run).
- `homework/results/p1_metrics.tsv` exists with columns: `coverage`, `n_contigs`, `total_bp`, `n50_bp`, `qv`.
- `homework/results/p1_metrics_plot.png` exists with four panels (one per metric) vs coverage.
- `homework/notes/p1-coverage-sweep.md` contains four numbered answers.
- Commit message like `p1: coverage sweep on 1Mb, single contig at 30x, QV plateau at 50x, 99% length at 25x`.

**Hint.** Use `subprocess.run` to call Flye three times in a loop; reuse the Exercise 1 orchestrator's `run_exercise` function with different `--coverage` arguments. The matplotlib plot is a 2x2 grid of `axes.plot(coverages, metric_values)` calls.

**Estimated time.** 75 minutes.

---

## Problem 2 — Medaka model sensitivity

**Problem statement.** Take the Exercise 1 Flye assembly. Run Medaka three times: once with the correct model (`r1041_e82_400bps_sup_v4.3.0`), once with a deliberate mismatch (`r941_min_sup_g507`, which is for R9.4.1 Guppy SUP), and once with a different correct-chemistry variant (`r1041_e82_400bps_hac_v4.3.0`). Compute QV for each polish.

Answer in `homework/notes/p2-medaka-models.md`:

1. What was the QV for each of the three polishes?
2. Did the mismatched-chemistry polish produce a worse QV than no polish?
3. Did the SUP vs HAC variants produce different QVs? By how much?
4. Why should the run-info JSON record the exact Medaka model string, not just "Medaka"?

**Acceptance criteria.**

- `homework/p2_medaka_models.py` runs three Medaka polishes (the script must call `assert_medaka_chemistry_matches` with `expected_chemistry_prefix="r1041"` for the first and third runs, and bypass the check for the second deliberate-mismatch run by passing `expected_chemistry_prefix="r941"` and noting it as a "did this on purpose" comment).
- `homework/results/p2_qv_comparison.tsv` exists with columns: `model`, `qv`, `mismatches`, `indels`.
- `homework/notes/p2-medaka-models.md` contains four numbered answers.
- Commit message like `p2: medaka model sensitivity, correct R10.4 SUP gave QV 38.2, R9.4 SUP gave QV 22.1 (worse than no polish), R10.4 HAC gave QV 35.7`.

**Hint.** The Medaka models are at <https://github.com/nanoporetech/medaka/blob/master/medaka/options.py>. Download them ahead of time with `medaka tools download_models r941_min_sup_g507 r1041_e82_400bps_sup_v4.3.0 r1041_e82_400bps_hac_v4.3.0` to avoid network hits during the script.

**Estimated time.** 60 minutes.

---

## Problem 3 — Flye vs Canu on the 5 Mb repeat-rich reference

**Problem statement.** Take the Challenge 2 5 Mb repeat-rich reference. Simulate 50x of badread reads. Run Flye and Canu on the same FASTQ. Compare on: (a) wall-clock time, (b) peak memory (use `/usr/bin/time -v` on Linux or `gtime -v` from `brew install gnu-time` on macOS), (c) N50, (d) contig count, (e) how each handled the three deliberate repeats (tandem 10x, dispersed 5x, inverted 3 kb).

Answer in `homework/notes/p3-flye-vs-canu.md`:

1. Which assembler was faster on wall-clock? By how much?
2. Which assembler used more peak memory? By how much?
3. Did Canu resolve any of the three repeats that Flye collapsed? Did Flye resolve any that Canu collapsed?
4. Which assembler would you pick if your job was "assemble a heavily repetitive Archaeon genome with limited memory"?

**Acceptance criteria.**

- `homework/p3_flye_vs_canu.py` runs both assemblers in subprocess and captures the wall-clock and memory.
- `homework/results/p3_comparison.tsv` exists with one row per assembler (columns: `assembler`, `wall_clock_sec`, `peak_memory_mb`, `n_contigs`, `total_bp`, `n50_bp`).
- `homework/notes/p3-flye-vs-canu.md` contains four numbered answers.
- Commit message like `p3: Flye vs Canu on 5Mb repeat-rich, Flye 3min/2GB, Canu 18min/8GB, Flye lost tandem 10x, Canu resolved it`.

**Hint.** Use `time.perf_counter()` for wall-clock; for memory, parse the `Maximum resident set size` line from `/usr/bin/time -v` stderr. On macOS install GNU time: `brew install gnu-time` and call `gtime` instead of `/usr/bin/time`.

**Estimated time.** 90 minutes (Canu alone is ~15-20 minutes of wall-clock).

---

## Problem 4 — Hifiasm on a simulated diploid input

**Problem statement.** Simulate two "haplotypes" by taking the 1 Mb reference and introducing 50 random SNPs (one every 20 kb) plus 5 random 1 kb structural variants into a copy. Concatenate the two FASTAs and call badread on the result at 50x coverage with HiFi-style parameters. Run Hifiasm with `--n-hap 2`. Compare the primary and alternate contig sets to the two source haplotypes.

Answer in `homework/notes/p4-hifiasm-diploid.md`:

1. How many contigs in the primary set? In the alternate set?
2. What is the N50 of each set?
3. Did Hifiasm correctly assign reads to the two haplotypes? (Compute the fraction of SNPs that appear correctly resolved between the two contig sets.)
4. What is the most important reproducibility caveat for diploid HiFi assembly?

**Acceptance criteria.**

- `homework/p4_diploid_hifiasm.py` runs end to end, including the haplotype-construction step.
- `homework/results/p4_hap_a.fasta`, `p4_hap_b.fasta` (the constructed haplotypes), `p4_reads.fastq`, and the Hifiasm output FASTAs all exist.
- `homework/notes/p4-hifiasm-diploid.md` contains four numbered answers.
- Commit message like `p4: Hifiasm diploid sim, 1 primary + 1 alt contig, 47/50 SNPs correctly phased`.

**Hint.** Random-mutation step: read the reference, walk every 20 kb, pick a random non-original base, write the mutated FASTA. Use `random.seed(42)` so the mutations are reproducible. For SV insertion, pick 5 random positions and insert a 1 kb run of random bases.

**Estimated time.** 75 minutes.

---

## Problem 5 — BUSCO on a real E. coli assembly

**Problem statement.** Take a real *E. coli* K-12 MG1655 ONT R10.4 dataset from the SRA / ENA (run accession `SRR23984905` or similar, or use a pre-bundled `data/ecoli_subset.fastq` of ~10x coverage). Run Flye + Medaka and then BUSCO with `bacteria_odb10` and again with `enterobacterales_odb10`. Compare the two BUSCO scores.

Answer in `homework/notes/p5-busco-ecoli.md`:

1. What was the BUSCO C, S, D, F, M for `bacteria_odb10`?
2. What was the BUSCO C, S, D, F, M for `enterobacterales_odb10`?
3. Did the more-specific lineage produce more or fewer total orthologs? More or fewer Missing?
4. Which lineage would you use in the published methods section of an *E. coli* genome paper, and why?

**Acceptance criteria.**

- `homework/p5_busco_ecoli.py` runs the full pipeline.
- `homework/results/p5_busco_bacteria.txt` and `p5_busco_enterobacterales.txt` (the BUSCO short_summary files) exist.
- `homework/notes/p5-busco-ecoli.md` contains four numbered answers.
- Commit message like `p5: BUSCO on E coli, bacteria_odb10 C=99.2% F=0% M=0.8%, enterobacterales_odb10 C=98.4% F=0.5% M=1.1%`.

**Hint.** Pre-cache the BUSCO downloads: `busco --download bacteria_odb10 enterobacterales_odb10`. The `enterobacterales_odb10` dataset has ~440 orthologs (vs 124 for `bacteria_odb10`); the percentages will be similar but on much larger denominators, which makes small differences statistically meaningful.

**Estimated time.** 75 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 400-500 word reflection at `homework/notes/week-10-reflection.md` answering:

1. Before Week 10, what did you think a "long read" was? After Week 10, what is it actually? Pick one stage of the pipeline (basecalling, assembly, polishing, QC) and say what surprised you about how the canonical defaults are set.
2. The Flye repeat graph (Kolmogorov et al. 2019) is the technical heart of modern OLC assembly. After reading the paper (or the lecture's summary), what is the single most important conceptual move it made over the older string-graph approach?
3. BUSCO is widely reported in published assemblies. After Week 10, what does a BUSCO score actually tell you, and what does it *not* tell you? Give two concrete scenarios where a high BUSCO would mislead you.
4. The mini-project produces a polished assembly FASTA plus a `run-info.json`. Imagine you hand the FASTA to a non-bioinformatician colleague who asks "what does this assembly mean?". What is the *most important* sentence you would say first? Why?

**Acceptance criteria.**

- File exists, 400-500 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. The boundaries you note here are what will keep you out of trouble in any future job that touches assembly.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 75 min |
| 2 | 60 min |
| 3 | 90 min |
| 4 | 75 min |
| 5 | 75 min |
| 6 | 30 min |
| **Total** | **~6.7 h** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
