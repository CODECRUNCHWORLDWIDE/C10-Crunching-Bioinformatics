# Week 10 — Exercise Solutions

> **Reproducibility note.** The expected outputs below assume Flye 2.9.5, Canu 2.2, Hifiasm 0.19.9, Medaka 1.12.0, BUSCO 5.7.1, badread 0.4.1, minimap2 2.28, and Biopython 1.84. Numbers will differ slightly on older or newer tool versions. The shape of the answers will not.

Each solution names the file you should write, the function bodies the reference implementation expects, and the expected numbers on the demo 1 Mb dataset.

---

## Solution 1 — badread + Flye end to end

**File:** `exercise-01-flye-via-subprocess.py`. The starter is already runnable; the work is internalizing what each helper function does and why each parameter is pinned.

### Expected output

Running:

```bash
python exercise-01-flye-via-subprocess.py \
    --reference data/reference_1mb.fasta \
    --out-dir results/ex01 \
    --seed 42
```

produces:

```
results/ex01/reads.fastq              ~50 MB (3,500 records, mean ~15 kb)
results/ex01/flye_out/assembly.fasta  1 contig, ~1.04 Mb
results/ex01/flye_out/assembly_info.txt
                                       contig_1   length=1038203   coverage=49   circular=Y
results/ex01/flye_out/assembly_graph.gfa
                                       1 S line, 1 L line (self-loop), 1 P line
results/ex01/run-info.json             { run_date, flye_version: 2.9.5, ... }
```

The badread simulation produces ~3,500 reads on a 1 Mb reference at 50x coverage (mean read length 15 kb means 1,000,000 / 15,000 reads per coverage unit, times 50 coverage = ~3,300 reads, plus the length-distribution tail). Flye in `--nano-hq` mode assembles them into a single ~1.04 Mb circular contig in ~30-60 seconds.

### Discussion points

1. **Why `--seed 42` is on the badread call but not (directly) on the Flye call?** badread is the explicit randomness source: the same seed always produces the same FASTQ. Flye is mostly deterministic given the same input and the same thread count; the only non-determinism is the thread-merge wobble we discussed in Lecture 2 §10. Pinning `--threads 4` is the Flye-side reproducibility lever.
2. **Why `--genome-size 1m` and not auto-detection?** Flye uses the genome size hint for coverage normalization in the overlap stage. The hint can be off by 2-3x without much impact; off by 10x degrades the assembly. For a 1 Mb reference passing `1m` is the obvious choice.
3. **Why `--nano-hq` rather than `--nano-raw`?** The badread `nanopore2023` error model simulates the ~5% per-base error rate characteristic of Dorado-SUP-basecalled R10.4 reads. That is what `--nano-hq` is tuned for. Using `--nano-raw` would tell Flye to expect ~10-15% per-base error and the overlap thresholds would loosen, producing more false overlaps.

### Sanity check: same seed + same input + same thread count = byte-identical output

```bash
python exercise-01-flye-via-subprocess.py --reference data/reference_1mb.fasta --out-dir results/ex01a --seed 42
python exercise-01-flye-via-subprocess.py --reference data/reference_1mb.fasta --out-dir results/ex01b --seed 42
diff results/ex01a/reads.fastq results/ex01b/reads.fastq && echo "byte-identical fastq"
diff <(grep -v "^#" results/ex01a/flye_out/assembly_info.txt) \
     <(grep -v "^#" results/ex01b/flye_out/assembly_info.txt) \
     && echo "byte-identical assembly_info"
```

The FASTQ diff should be empty. The Flye assembly may differ on a handful of bases in the contig FASTA due to multi-threaded merge order; the `assembly_info.txt` should report the same contig count, length, and coverage.

---

## Solution 2 — asmstats + BUSCO

**File:** `exercise-02-asmstats-n50-busco.py`. The reference implementation contains the N50 / L50 walk, the BUSCO parser, and the Markdown report writer.

### Expected output

```bash
python exercise-02-asmstats-n50-busco.py \
    --assembly results/ex01/flye_out/assembly.fasta \
    --out-dir results/ex02 \
    --lineage bacteria_odb10 \
    --skip-busco
```

produces:

```
results/ex02/asmstats.tsv             2 lines (one per contig + a summary comment)
results/ex02/qc_report.md             Markdown with asmstats table, no BUSCO block
results/ex02/run-info.json            { run_date, asmstats: {...}, busco_skipped: true }
```

Expected asmstats for the demo 1 Mb assembly:

| Metric              | Value         |
|---------------------|---------------|
| n_contigs           | 1             |
| total_length_bp     | 1,038,203     |
| longest_contig_bp   | 1,038,203     |
| shortest_contig_bp  | 1,038,203     |
| N50                 | 1,038,203     |
| L50                 | 1             |
| GC fraction         | 0.508         |

For a real bacterial genome with BUSCO enabled (`--lineage bacteria_odb10`, BUSCO 5.7.1, no `--skip-busco`):

| BUSCO metric        | Expected      |
|---------------------|---------------|
| Complete (C)        | 98-99%        |
| Single-copy (S)     | 98-99%        |
| Duplicated (D)      | 0-1%          |
| Fragmented (F)      | 0-2%          |
| Missing (M)         | 0-2%          |

### Discussion points

1. **Why does the synthetic 1 Mb reference score BUSCO C ~= 0%?** BUSCO looks for orthologs from a curated gene set; the synthetic FASTA is random nucleotides with no real genes. To exercise BUSCO realistically, point the script at an *E. coli* assembly (use `data/ecoli_reference.fasta` if bundled, or download `NC_000913.3` from GenBank).
2. **Why pin the BUSCO lineage creation date in the run-info?** Two `bacteria_odb10` releases (say 2023-02-01 and 2024-01-08) have slightly different ortholog sets; the same assembly can score differently. The parser extracts the creation date from the summary header and the run-info JSON records it.
3. **Why the walk-cumulative-then-break loop for N50?** Stopping early avoids walking the full sorted list when L50 is small. For a single-contig assembly the loop runs once.

### Sanity check: re-running on the same assembly is idempotent

```bash
python exercise-02-asmstats-n50-busco.py --assembly results/ex01/flye_out/assembly.fasta --out-dir results/ex02a --skip-busco
python exercise-02-asmstats-n50-busco.py --assembly results/ex01/flye_out/assembly.fasta --out-dir results/ex02b --skip-busco
diff results/ex02a/asmstats.tsv results/ex02b/asmstats.tsv && echo "byte-identical"
```

The TSV diff should be empty.

---

## Solution 3 — Medaka polish + QV comparison

**File:** `exercise-03-medaka-polish.py`. The reference implementation runs medaka_consensus, aligns the draft and the polished assemblies to the reference with minimap2 (`-x asm5`), parses the resulting SAM to tally edit distance, and emits a Markdown comparison report.

### Expected output

```bash
python exercise-03-medaka-polish.py \
    --reads results/ex01/reads.fastq \
    --draft results/ex01/flye_out/assembly.fasta \
    --reference data/reference_1mb.fasta \
    --out-dir results/ex03 \
    --medaka-model r1041_e82_400bps_sup_v4.3.0
```

produces:

```
results/ex03/medaka_out/consensus.fasta   the polished assembly
results/ex03/polished.fasta               copy of consensus.fasta at the top level
results/ex03/raw_vs_ref.sam               draft -> reference alignment
results/ex03/polished_vs_ref.sam          polished -> reference alignment
results/ex03/qv_report.md                 Markdown report with QV comparison
results/ex03/run-info.json                { medaka_model, raw_qv, polished_qv, ... }
```

Expected QV on the demo 1 Mb assembly:

| Metric          | Raw (Flye)   | Polished (Medaka) |
|-----------------|--------------|-------------------|
| Aligned bp      | 1,037,901    | 1,038,043         |
| Mismatches      | ~700-1,200   | ~50-150           |
| Insertions      | ~80-200      | ~30-80            |
| Deletions       | ~100-300     | ~40-120           |
| Error rate      | ~0.001-0.002 | ~0.0001-0.0003    |
| QV (Phred)      | ~25-30       | ~35-40            |

The polish improvement on simulated nanopore reads is typically 5-10x in error rate (a +5-10 in QV). On real R10.4 data the improvement is usually closer to 3-5x; the simulator is a little easier than reality.

### Discussion points

1. **Why does Medaka require the basecaller-matched model?** Medaka's neural network is trained on a specific basecaller's error distribution. Applying an R9.4-trained model to R10.4 reads asks the model to correct errors that are not there and skip errors that are; the polish degrades the assembly. The `assert_medaka_chemistry_matches` helper enforces this at the chemistry prefix.
2. **Why minimap2 `-x asm5` and not `-x map-ont`?** `asm5` is the assembly-to-reference preset, expecting < 5% divergence and end-to-end mapping; `map-ont` is the long-read-to-reference preset, expecting many overlapping reads. For a single-contig assembly aligned to its source reference, `asm5` is correct.
3. **Why prefer the NM tag over CIGAR's `=` / `X` counts when both are present?** Minimap2 with default settings writes M (match-or-mismatch); the NM tag is the explicit edit distance and is set even without `=` / `X` mode. Using NM gives us a single source of truth across minimap2 invocations.
4. **Why is QV capped at 60?** A zero-error alignment has `error_rate = 0` and the log goes to negative infinity. Q60 is the conventional ceiling for "effectively zero error rate"; for QV > 60 the error rate is below the resolution of the alignment-based estimator.

### Sanity check: skip-medaka path

```bash
python exercise-03-medaka-polish.py --reads results/ex01/reads.fastq \
    --draft results/ex01/flye_out/assembly.fasta \
    --reference data/reference_1mb.fasta \
    --out-dir results/ex03_skip \
    --skip-medaka
```

This should produce the same raw QV in both `raw_qv` and `polished_qv` (the script copies the draft as the polished output when Medaka is skipped). The QV report's "Polish improved QV" line should be absent.

---

## A note on running these without the external tools installed

The Python files compile under `python3 -m py_compile` regardless of which tools are on the PATH (Biopython, badread, Flye, Medaka, BUSCO, minimap2 are all imported lazily inside the functions that need them; the external tools are called via `subprocess.run`). To *run* the exercises you need:

- **Exercise 1:** badread and Flye on the PATH. Biopython for parsing.
- **Exercise 2:** Biopython. BUSCO is optional (use `--skip-busco`).
- **Exercise 3:** minimap2 for QV. Medaka is optional (use `--skip-medaka`).

Conda one-liner:

```bash
conda install -c bioconda -c conda-forge \
    python=3.11 flye=2.9.5 badread=0.4.1 medaka=1.12.0 \
    busco=5.7.1 minimap2=2.28 biopython=1.84
```

The Medaka install on conda has occasionally been finicky on Apple silicon; if it fails, install via pip in a fresh venv (`pip install medaka==1.12.0`). The BUSCO download for `bacteria_odb10` is ~50 MB; pre-cache it with `busco --download bacteria_odb10` to speed up offline runs.

---

## What to commit

After running all three exercises, your `week-10/exercises/` directory should contain:

```
exercises/
    exercise-01-flye-via-subprocess.py
    exercise-02-asmstats-n50-busco.py
    exercise-03-medaka-polish.py
    SOLUTIONS.md
    results/
        ex01/{reads.fastq, flye_out/{assembly.fasta, assembly_graph.gfa,
              assembly_info.txt, flye.log}, run-info.json}
        ex02/{asmstats.tsv, qc_report.md, run-info.json,
              [busco_out/busco_run/short_summary*.txt]}
        ex03/{polished.fasta, raw_vs_ref.sam, polished_vs_ref.sam,
              qv_report.md, run-info.json}
```

Gitignore the Flye intermediate subdirectories (`00-assembly/`, `10-consensus/`, `20-repeat/`, `30-contigger/`, `40-polishing/`), the Medaka intermediate BAMs, and the BUSCO download cache. The deliverables to keep are: the run-info JSONs, the QC reports, the asmstats TSV, and the polished assembly FASTA.
