# Mini-Project — End-to-end long-read assembly pipeline

> **Reproducibility note.** This mini-project produces a polished long-read assembly from a simulated FASTQ. The assembly is reproducible only if you ship the inputs, the parameters, the seed, and the tool versions alongside the FASTA. Without the `run-info.json` the assembly is an opinion; with it, the assembly is a reproducible result. **The output of this pipeline must travel with its run-info, every time.**

Build a reproducible long-read assembly pipeline that takes a reference FASTA (the bundled `data/reference_1mb.fasta`, a 1 Mb scaffold), simulates a 50x nanopore-style read set with `badread`, filters the reads by length and quality, runs **Flye** to assemble them into contigs, polishes the contigs with **Medaka**, computes assembly statistics (N50, L50, total length, GC), runs **BUSCO** for gene-content completeness, computes QV against the reference, renders the assembly graph with **Bandage CLI**, and emits a `run-info.json` recording every parameter. End with: the simulated FASTQ, the Flye assembly FASTA + GFA, the Medaka-polished FASTA, the BUSCO summary, the QV report, the rendered graph PNG, and a 600-800 word write-up that defends every parameter and names every limit.

This is the C10 mini-project that produces a **methods-section-quality assembly artefact with measured provenance**, not just a one-shot demonstration. By the end of it you will have a `long_read_pipeline.py` script and a `run.sh` wrapper, a results directory with the assembly and the run-info, and a write-up that defends every parameter and explicitly names what the assembly cannot reliably claim.

**Estimated time:** 7 hours (split across Wednesday, Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-10/mini-project/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 10 section)
└── week-10/
    └── mini-project/
        ├── README.md               one-page report (~600-800 words)
        ├── run.sh                  one-command reproduction script
        ├── env.yml                 conda environment file pinning all tool versions
        ├── data/
        │   └── reference_1mb.fasta input reference (1 Mb)
        ├── long_read_pipeline.py   the orchestration script
        ├── starter.py              skeleton implementation with TODOs
        └── results/
            ├── reads.fastq                       badread simulated FASTQ
            ├── reads_filtered.fastq              after length / QV filter
            ├── flye_out/
            │   ├── assembly.fasta                Flye draft contigs
            │   ├── assembly_info.txt             per-contig length / coverage
            │   ├── assembly_graph.gfa            assembly graph (Bandage-compatible)
            │   └── flye.log
            ├── medaka_out/consensus.fasta        Medaka-polished assembly
            ├── polished.fasta                    final assembly (copy of consensus.fasta)
            ├── asmstats.tsv                      N50, L50, contig table
            ├── qc_report.md                      asmstats + BUSCO + QV Markdown report
            ├── busco_out/                        BUSCO output directory
            ├── bandage_graph.png                 graph rendering
            ├── polished_vs_ref.sam               minimap2 alignment to reference
            └── run-info.json                     run provenance
```

By the end you will have a clean, reproducible Week 10 directory you can point a recruiter at — and `long_read_pipeline.py` is the kind of pipeline that opens conversations with working bioinformaticians and biotech / academic shops, *as long as* you can speak to its limits.

---

## The dataset

You will work with the bundled `data/reference_1mb.fasta` — a 1.04 Mb synthetic reference plus three deliberate repeat features (a tandem 5x repeat of a 100 bp unit; a dispersed 3x repeat of a 1 kb unit; one inverted 2 kb repeat). The repeats are smaller than the Challenge 2 case and the assembler should resolve all three at 50x coverage with 15 kb mean read length. If your assembly collapses one of them, the mini-project write-up flags it as the expected failure mode and explains the cause.

The choice of "synthetic with deliberate repeats" rather than "real E. coli" is for two reasons:

- **Reproducibility.** Real sequencer data has a non-reproducible random-error component that the simulator does not; pinning the simulator seed gives byte-identical reads across runs and lets you focus on the assembly side.
- **Pedagogical clarity.** The deliberate repeats let the write-up name what the assembler did and did not handle, with the truth visible (you constructed the repeats yourself).

Expected coverage by the pipeline:

- 50x simulated reads should produce ~3,500 records with mean length ~15 kb.
- Flye should produce 1-2 contigs totaling ~1.04 Mb; the longest contig should be > 500 kb.
- Medaka should improve the QV from ~25 (raw Flye) to ~38 (after polish).
- BUSCO on the synthetic reference is not informative (no real genes); skip BUSCO with `--skip-busco` or use the *E. coli* extension below.

### Optional E. coli extension

If you want a BUSCO-meaningful run, swap the reference for a real *E. coli* K-12 MG1655 sequence (NCBI accession `U00096.3`). The pipeline should still work without modification (Flye, Medaka, BUSCO all run on real data identically); the write-up will then have a meaningful BUSCO score to report (expected: C >= 99% on `bacteria_odb10`).

---

## Rules

- **You may** use Flye 2.9.5, Medaka 1.12.0, BUSCO 5.7.1, minimap2 2.28, badread 0.4.1, Bandage 0.9.0, Biopython 1.84, and the standard library.
- **You may** consult Lectures 1, 2, 3, the Flye paper (Kolmogorov et al. 2019), the Canu paper (Koren et al. 2017), the Hifiasm paper (Cheng et al. 2021), the HiFi paper (Wenger et al. 2019; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6776680/>), the BUSCO paper (Manni et al. 2021), the Bandage paper (Wick et al. 2015), and the Week 10 exercises and challenges.
- **You may NOT** copy a pre-written assembly pipeline from the internet. The point is to *build* the pipeline. Reading the Flye GitHub README for inspiration is fine; copy-pasting a complete pipeline is not.
- **You must** pin every tool version, every input mode flag, the read filter parameters, the Medaka model, the BUSCO lineage, and the seed in the `run-info.json`.
- **You must** verify that the Medaka model's chemistry prefix matches the basecaller's chemistry (here: simulated `nanopore2023`, which maps to `r1041`).
- **You should** ship a Markdown QC report alongside the run-info JSON; the JSON is for machines and the Markdown is for humans.
- **You must** commit the input reference, the simulated FASTQ (gzipped), the Flye assembly FASTA + GFA, the polished FASTA, the QC report, the Bandage PNG, and the `run-info.json`. Gitignore the Flye intermediate stages, the Medaka intermediate BAMs, and the BUSCO download cache.

---

## Acceptance criteria

- [ ] `mini-project/long_read_pipeline.py` exports a function `assemble_genome(reference_fasta: Path, out_dir: Path, seed: int = 42) -> Path` that runs the full pipeline and returns the path to the final polished FASTA.
- [ ] The pipeline implements **eight stages**:
  1. **Validate input.** Confirm the reference FASTA exists, has at least one record, total length >= 1,000 bp.
  2. **Simulate reads.** badread at 50x nanopore-style coverage with the pinned seed.
  3. **Filter reads.** Drop reads shorter than 1,000 bp or with mean Phred < 10.
  4. **Run Flye.** `--nano-hq --genome-size <auto> --threads 4 --iterations 1`.
  5. **Polish with Medaka.** `medaka_consensus -m r1041_e82_400bps_sup_v4.3.0 -t 4`.
  6. **Compute asmstats.** N50, L50, total length, GC fraction.
  7. **QV against reference.** minimap2 `-x asm5` then SAM tally.
  8. **Render the graph.** `Bandage image ... --height 800`.
- [ ] `long_read_pipeline.py` produces:
  - `results/reads.fastq` — badread output.
  - `results/reads_filtered.fastq` — after the read filter.
  - `results/flye_out/assembly.fasta`, `assembly_info.txt`, `assembly_graph.gfa`.
  - `results/medaka_out/consensus.fasta` and `results/polished.fasta`.
  - `results/asmstats.tsv` — per-contig + summary.
  - `results/qc_report.md` — Markdown report.
  - `results/polished_vs_ref.sam` — minimap2 alignment.
  - `results/bandage_graph.png` — rendered graph (or a `.txt` describing the topology if Bandage is unavailable).
  - `results/run-info.json` — versions, dates, parameters.
- [ ] `mini-project/README.md` is a one-page (~600-800 word) report containing:
  - One-sentence description of the input, the simulator, the assembler, and the polisher.
  - Methods section in C10 voice: every tool pinned ("Flye 2.9.5", "Medaka 1.12.0", "minimap2 2.28"), every parameter explicit, the seed stated.
  - Results section in C10 voice: read count after filter, contig count, N50, L50, longest contig, raw QV, polished QV, the deepest-supported assembly feature (the longest contig).
  - Discussion section: 150-250 words on the limits of the pipeline. What is *not* in the assembly? What does the N50 *not* tell you? Which standard failure modes (collapsed repeats, basecaller-model drift, polish-model mismatch) is the pipeline plausibly susceptible to?
- [ ] `run.sh` is a single bash script that, given a fresh checkout + `conda env create -f env.yml`, reproduces the entire pipeline from scratch in under ten minutes.
- [ ] The repo is **public** and at least one classmate or instructor has been added as a collaborator.
- [ ] **Most importantly**: the run-info JSON is complete and the discussion section is honest about the limits.

---

## Suggested approach (rough timeline)

### Wednesday (1 hour)

1. (15 min) `git clone`, set up `mini-project/` directory.
2. (30 min) Write `env.yml` with pinned tool versions; create conda env; verify each tool is on the PATH (badread, Flye, Medaka, minimap2, Biopython, and optionally BUSCO and Bandage).
3. (15 min) Read this README end to end and the Challenge 1 / Challenge 2 READMEs. Sketch the eight-stage flow on paper.

### Thursday (2 hours)

1. (45 min) Implement Stages 1-3 (validate, simulate, filter) in `long_read_pipeline.py`. Reuse the Exercise 1 helpers; they cover most of this.
2. (45 min) Implement Stages 4 (Flye) and 6 (asmstats). Reuse the Exercise 1 and Exercise 2 helpers.
3. (30 min) Implement Stage 5 (Medaka polish). Reuse the Exercise 3 helpers; include the chemistry-prefix assertion.

### Friday (2 hours)

1. (45 min) Implement Stage 7 (QV via minimap2) and Stage 8 (Bandage). Both have graceful-skip paths.
2. (30 min) Write `run.sh` and `env.yml`. Do a full `bash run.sh` on a fresh `conda env create -f env.yml` to verify end-to-end reproducibility.
3. (45 min) Draft the README.md write-up (methods + results + 200-word discussion).

### Saturday (2 hours)

1. (30 min) Optional: swap the synthetic reference for *E. coli* and re-run; the BUSCO score becomes meaningful.
2. (30 min) Polish the QC report (add the Bandage PNG inline; format the asmstats table; add the limits section).
3. (30 min) Final edit of the write-up.
4. (30 min) Push, add a classmate as a collaborator, write the commit messages with specific numbers.

---

## Tips for the write-up

- **Lead with numbers.** "badread simulated 3,457 nanopore-style reads at 50x coverage on the 1.04 Mb reference; after the length / QV filter (min 1,000 bp, min Q10), 3,289 reads remained. Flye 2.9.5 in --nano-hq mode assembled them into 1 circular contig of length 1,038,203 bp in 47 seconds; the N50 is 1,038,203 bp and the L50 is 1. Medaka 1.12.0 with the `r1041_e82_400bps_sup_v4.3.0` model polished the draft to a QV of 38.2 (mismatches 67, indels 39, error rate 1.0e-4); the raw Flye QV was 25.4 (mismatches 1,103, indels 287). The Bandage rendering shows one self-loop on the primary contig: the reference closed cleanly."
- **Defend every parameter.** "We used Flye `--nano-hq` rather than `--nano-raw` because the badread `nanopore2023` error model corresponds to Dorado SUP / R10.4.1 reads, not R9.4 Guppy reads; we used Medaka with the `r1041_e82_400bps_sup_v4.3.0` model to match the basecaller chemistry (a mismatched model would degrade the polish); we used 50x coverage because the empirical floor for nanopore Flye is ~25x and 50x leaves comfortable headroom; the random seed is 42 by curriculum convention."
- **Name two failure modes the pipeline is susceptible to.** "The synthetic reference contains three deliberate repeats (tandem 5x, dispersed 3x, inverted 2x); on this run all three resolved correctly, but at lower coverage (20x) the tandem 5x repeat collapses into a single high-coverage segment. The Medaka model mismatch is the second failure mode the pipeline is structurally vulnerable to; we mitigate with the `assert_medaka_chemistry_matches` check at the start of Stage 5."
- **Acknowledge what is missing.** "BUSCO is not informative on the synthetic reference (no real genes); for a real-data validation we would re-run on an *E. coli* reference and report the `bacteria_odb10` score. Real ONT data also has correlated errors that the simulator's i.i.d. model does not capture; a small fraction of true positives may be artificially easy in this pipeline."

---

## Stretch goals (optional)

- Add the Canu branch (`canu -nanopore reads.fastq genomeSize=1m useGrid=false maxThreads=4`) and compare to Flye on the same FASTQ. Use the Problem 3 comparison structure.
- Add the Hifiasm-on-HiFi branch (simulate HiFi-style reads from the same reference; run Hifiasm; compare). This is essentially Challenge 1 inlined into the mini-project.
- Compute the per-contig coverage profile by aligning the reads back to the assembly (`minimap2 -x map-ont`) and reporting any contig with coverage > 1.5x the median. This is the "did the assembler collapse a repeat?" check.
- Render the Bandage graph in both a basic CLI image and the Bandage GUI (taking a screenshot manually). The GUI view shows the topology more clearly for small graphs; the CLI image is what your run.sh produces reproducibly.
- Reproduce a small published long-read assembly end-to-end. The *E. coli* K-12 MG1655 ONT R10.4 dataset on the EBI ENA (e.g. `SRR23984905`) is a common target.

---

## What to commit

By the end of the mini-project your `week-10/mini-project/` should contain:

```
mini-project/
    README.md
    env.yml
    run.sh
    long_read_pipeline.py
    starter.py
    data/reference_1mb.fasta
    results/
        reads.fastq.gz                            (gzipped)
        reads_filtered.fastq.gz                   (gzipped)
        flye_out/{assembly.fasta, assembly_info.txt, assembly_graph.gfa, flye.log}
        medaka_out/consensus.fasta
        polished.fasta
        asmstats.tsv
        qc_report.md
        polished_vs_ref.sam
        bandage_graph.png                         (or topology.txt if Bandage unavailable)
        run-info.json
```

Gitignore the Flye intermediate stages (`00-assembly/`, `10-consensus/`, `20-repeat/`, `30-contigger/`, `40-polishing/`), the Medaka intermediate BAMs, the BUSCO download cache, and any uncompressed FAST5 / POD5 archives (not used in the simulated case; relevant for the *E. coli* extension). The commit message for the final pipeline run should be specific, e.g. `mini-project: long-read pipeline on 1Mb sim, 1 circular contig, raw QV 25.4 -> polished QV 38.2`.
