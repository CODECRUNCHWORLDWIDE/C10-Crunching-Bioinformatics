# Challenge 1 — Hifiasm on simulated HiFi reads, no post-assembly polish

> **Reproducibility note.** Hifiasm is mostly deterministic given the same input and thread count, but the haplotype-phasing stage uses a stable-sort tie-breaker that can produce slightly different primary / alternate splits on different machines. Pin the version (`Hifiasm 0.19.9`), pin the thread count (`-t 4`), and pin the badread seed (`--seed 42`) and the result is reproducible across runs on the same machine. Compare across machines only if you re-run the simulation locally.

**Estimated time:** 2.5 hours.
**Goal:** Swap Flye for Hifiasm on simulated PacBio-HiFi-style reads. Hifiasm requires no post-assembly polish (HiFi reads are already QV30+), so the pipeline is shorter than the ONT Flye + Medaka path. Compare the Hifiasm output to a Flye-on-HiFi run on the same FASTQ; report the contiguity, the QV, and the wall-clock time for each.

This challenge is the bridge between "I can assemble ONT reads" and "I can pick the right assembler for the input chemistry." The conceptual lift is: HiFi reads are different from ONT reads at every layer (length distribution, error mode, per-base accuracy), and the assembler choice should follow.

---

## Background — Why a HiFi-specialized assembler beats a general-purpose OLC

A Flye assembly of HiFi reads (`--pacbio-hifi`) works but is suboptimal:

1. **Flye's error-correction stage is overkill on HiFi.** HiFi reads are ~99.8% accurate at the read level; Flye spends much of its compute correcting errors that are not there.
2. **Flye collapses haplotype-distinct regions.** A 1-SNP-different region between two haplotypes is folded into a single consensus contig. For haploid organisms this is fine; for diploid this loses information.
3. **Flye's output is unphased.** Even on haploid input, Hifiasm's primary / alternate contig structure is cleaner than Flye's single FASTA when minor heterozygous events are present (e.g. a contaminating second strain).

Hifiasm (Cheng et al. 2021, *Nat Methods* 18:170; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8005227/>) skips the correction stage, builds the string graph directly from raw HiFi reads, and emits haplotype-resolved contigs in a single pass. On HiFi input it is ~5x faster than Flye and produces a more interpretable graph.

---

## Task

Build a Python wrapper `hifiasm_pipeline.py` that:

1. Simulates HiFi-style reads from the same 1 Mb reference used in Exercise 1, but with PacBio-tuned badread parameters.
2. Runs Hifiasm with the pinned thread count and haploid expectation.
3. Converts the Hifiasm GFA output to FASTA.
4. Computes asmstats (reusing the Exercise 2 helpers).
5. Aligns the assembly to the reference with minimap2 `-x asm5` and computes QV (reusing the Exercise 3 helpers).
6. Re-runs the same FASTQ through Flye `--pacbio-hifi` for comparison.
7. Emits a Markdown comparison report and a `run-info.json` recording both assembler runs.

### Layout

```
crunch-bio-portfolio-<yourhandle>/
└── week-10/
    └── challenge-01/
        ├── README.md                how-to-run + write-up
        ├── env.yml                  conda env file
        ├── hifiasm_pipeline.py      the orchestration script
        ├── data/
        │   └── reference_1mb.fasta  shared with Exercise 1
        └── results/
            ├── reads_hifi.fastq                       badread HiFi-style simulation
            ├── hifiasm_out/ecoli.bp.p_ctg.gfa         Hifiasm primary contig graph
            ├── hifiasm_out/ecoli.fasta                converted FASTA
            ├── flye_hifi_out/assembly.fasta           Flye --pacbio-hifi assembly
            ├── hifiasm_vs_ref.sam                     minimap2 alignment
            ├── flye_hifi_vs_ref.sam
            ├── comparison_report.md                   side-by-side QC
            └── run-info.json
```

### badread parameters for HiFi-style reads

```bash
badread simulate \
    --reference data/reference_1mb.fasta \
    --quantity 50x \
    --length 18000,5000 \
    --identity 99.5,0.3,99.95 \
    --error_model pacbio2021 \
    --qscore_model pacbio2021 \
    --seed 42 \
    > results/reads_hifi.fastq
```

The narrower identity distribution (`99.5,0.3,99.95` vs `95,3,99` for nanopore) reflects HiFi's tighter per-read accuracy distribution. The `pacbio2021` error model places errors in the random-distributed pattern characteristic of HiFi reads (no homopolymer concentration; no methylation-related errors).

### Hifiasm call

```bash
hifiasm \
    -o results/hifiasm_out/ecoli \
    -t 4 \
    --n-hap 1 \
    results/reads_hifi.fastq
```

`--n-hap 1` is the haploid expectation; for diploid organisms you would set `--n-hap 2` and Hifiasm would emit a distinct alternate-contig FASTA.

### GFA-to-FASTA conversion

Hifiasm emits a GFA, not a FASTA. The conversion is a one-line awk:

```bash
awk '/^S/ {print ">"$2"\n"$3}' \
    results/hifiasm_out/ecoli.bp.p_ctg.gfa \
    > results/hifiasm_out/ecoli.fasta
```

Or in Python (use this for the orchestration script):

```python
from __future__ import annotations

from pathlib import Path


def gfa_to_fasta(gfa_path: Path, fasta_path: Path) -> int:
    """Extract 'S' lines from a GFA and write as FASTA. Returns the segment count."""
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    n_segments: int = 0
    with gfa_path.open() as gfa_fh, fasta_path.open("w") as fasta_fh:
        for line in gfa_fh:
            if not line.startswith("S"):
                continue
            parts: list[str] = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seg_id: str = parts[1]
            seg_seq: str = parts[2]
            fasta_fh.write(f">{seg_id}\n{seg_seq}\n")
            n_segments += 1
    return n_segments
```

### Flye-on-HiFi call (for comparison)

```bash
flye \
    --pacbio-hifi results/reads_hifi.fastq \
    --genome-size 1m \
    --out-dir results/flye_hifi_out \
    --threads 4 \
    --iterations 0
```

`--iterations 0` skips Flye's internal polish (HiFi is already accurate enough that an extra polish round has no effect; we save the wall-clock time).

---

## Acceptance criteria

- [ ] `hifiasm_pipeline.py` is a single CLI script that runs end to end with `python hifiasm_pipeline.py --reference data/reference_1mb.fasta --out-dir results --seed 42`.
- [ ] The script implements **six stages**: badread HiFi simulation; Hifiasm; GFA -> FASTA; Flye --pacbio-hifi (comparison); minimap2 alignment of both assemblies; QV tally and comparison report.
- [ ] `results/hifiasm_out/ecoli.fasta` exists with at least one contig totaling ~1.04 Mb.
- [ ] `results/flye_hifi_out/assembly.fasta` exists for the comparison.
- [ ] `results/comparison_report.md` is a one-page report with:
  - Side-by-side asmstats table for Hifiasm and Flye (n_contigs, total_length, N50, L50).
  - Side-by-side QV table for both assemblies.
  - Wall-clock time for each assembler (use `time.perf_counter()` in the script).
  - One paragraph naming which assembler "won" on each axis (contiguity, QV, wall-clock).
- [ ] `results/run-info.json` records: badread version + parameters, Hifiasm version + thread count, Flye version + input mode, minimap2 version, seed, run date, both QV results.
- [ ] The script gracefully skips Hifiasm or Flye if the binary is missing (print a warning, continue with the other assembler).

---

## Expected numbers

On the demo 1 Mb reference at 50x HiFi-style coverage (seed 42):

| Metric              | Hifiasm   | Flye --pacbio-hifi |
|---------------------|-----------|--------------------|
| n_contigs           | 1         | 1                  |
| total_length_bp     | ~1,040,000| ~1,040,000         |
| N50                 | ~1,040,000| ~1,040,000         |
| L50                 | 1         | 1                  |
| Mismatches vs ref   | < 50      | < 100              |
| Indels vs ref       | < 30      | < 50               |
| QV                  | 45-50     | 40-45              |
| Wall clock          | ~30 sec   | ~60-90 sec         |

Hifiasm should match or beat Flye on every axis on HiFi input. The QV advantage comes from Hifiasm's HiFi-tuned overlap thresholds; the wall-clock advantage comes from skipping the correction stage. On a diploid input the haplotype-resolution advantage would be much larger; for this haploid simulated dataset the two assemblers are very close.

---

## Write-up

In `challenge-01/README.md`, answer in ~400 words:

1. **What did Hifiasm do that Flye did not?** Name the algorithmic step (skipped the error-correction stage; built the string graph directly from raw HiFi reads).
2. **On the demo dataset, which assembler "won"?** Report the numbers and call out the win.
3. **For a real human genome (3 Gb, diploid, ~30x HiFi coverage on a Revio cell), would your answer change?** Hint: the haplotype-resolution gap widens dramatically on diploid input.
4. **What is the most important reproducibility caveat in this pipeline?** Hint: it is not the seed; it is the input mode flag.

---

## Stretch goals (optional)

- Try `--n-hap 2` on the same HiFi FASTQ and see how the primary and alternate contig sets differ. On a haploid simulated dataset both will be ~1 Mb; the difference is illustrative of the diploid pipeline.
- Add the `--ul` flag to Hifiasm with a simulated ONT ultra-long read set and assemble the same reference. The combined output should be more contiguous than HiFi alone on a repeat-rich region.
- Compute a per-contig coverage profile (using `samtools depth` on the read-to-assembly BAM) and check whether any contig has coverage > 1.5x the median (the collapsed-repeat warning sign).
- Reproduce a small published HiFi assembly. The *Saccharomyces cerevisiae* S288C HiFi dataset on the SRA (run accession `SRR15275213` or similar) is a common target; your Hifiasm assembly should agree with the GenBank reference (`U00091` to `U00103` for the 16 chromosomes) on >= 99.9% of bases.

---

## What to commit

```
challenge-01/
    README.md
    env.yml
    hifiasm_pipeline.py
    data/reference_1mb.fasta
    results/
        reads_hifi.fastq
        hifiasm_out/ecoli.bp.p_ctg.gfa
        hifiasm_out/ecoli.fasta
        flye_hifi_out/assembly.fasta
        flye_hifi_out/assembly_info.txt
        flye_hifi_out/assembly_graph.gfa
        hifiasm_vs_ref.sam
        flye_hifi_vs_ref.sam
        comparison_report.md
        run-info.json
```

Gitignore: Flye intermediate stages (`00-assembly/`, `10-consensus/`, etc.) and Hifiasm intermediate files (`*.bin`, `*.ec.fa`). The deliverables to keep are the GFA, the FASTA, the SAM, and the run-info JSON.

Commit message like `challenge-01: Hifiasm vs Flye on 1Mb HiFi sim, Hifiasm QV 47.3 / 30s, Flye 43.1 / 75s`.
