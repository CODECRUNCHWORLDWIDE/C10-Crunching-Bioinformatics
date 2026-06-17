# Challenge 2 — Bandage graph inspection on a repeat-rich assembly

> **Reproducibility note.** Bandage CLI image output is deterministic given the same input GFA and the same Bandage version; the GUI layout is stochastic and a manual screenshot will not be byte-identical across runs. For this challenge we use `Bandage image` (the headless CLI) so the rendered PNG is reproducible, and we describe the graph topology in words so the description survives layout changes.

**Estimated time:** 2 hours.
**Goal:** Run Flye on a deliberately repeat-rich 5 Mb reference, open the resulting `assembly_graph.gfa` in Bandage, classify the graph topology (linear / circular / bubble / tangle), identify which repeats the assembler resolved and which it collapsed, and write a one-page report defending your reading of the graph.

This challenge is the bridge between "I have an N50 number" and "I understand what my assembly looks like at the graph level." The conceptual lift is: the FASTA is the layout's output, but the GFA is the layout's reasoning, and reading the GFA is the only way to diagnose whether a clean-looking FASTA hides a collapsed repeat.

---

## Background — Why graph inspection matters

An assembly FASTA can lie about the underlying genome:

- **Overcollapsed repeats.** N tandem copies of a 1 kb unit are emitted as a single contig of length 1 kb (not N kb) with coverage N times the genome median. The N50 looks fine; the genome length looks short; the per-contig coverage gives the game away.
- **Undercollapsed repeats.** N copies are emitted as N separate contigs; the FASTA looks fragmented but the underlying biology is one cluster.
- **Chimeric contigs.** Two distant regions joined by a low-confidence overlap. The FASTA shows one contig; the GFA shows a bubble where the assembler chose one path.

Bandage's job is to make these visible. A clean bacterial assembly looks like one self-loop (circular chromosome) or one linear path; a problematic assembly has bubbles, tangles, or tips that the FASTA alone does not reveal.

---

## Task

Build a Python script `bandage_inspect.py` that:

1. Simulates reads from `data/reference_5mb_with_repeats.fasta` (provided; see "Reference construction" below).
2. Runs Flye in `--nano-hq` mode at 50x coverage.
3. Calls `Bandage image` to render the assembly graph as a PNG.
4. Parses the GFA file and reports:
   - The number of `S` (segment) lines.
   - The number of `L` (link) lines.
   - The number of `P` (path) lines.
   - For each segment: length, coverage (from the `dp:f:` tag), in-degree, out-degree.
   - A topology classification: `linear` (no links), `circular` (self-loop), `bubble` (multiple parallel paths), `tangle` (high-degree subgraph).
5. Writes a Markdown report describing the topology in words.

### Reference construction

The repeat-rich 5 Mb reference is built by appending three deliberate repeats to a random 4.7 Mb scaffold:

- **Tandem 10x repeat.** 10 head-to-tail copies of a 200 bp unit, inserted at position 1,500,000. Total 2 kb.
- **Dispersed 5x repeat.** 5 copies of a 1,500 bp unit, scattered at positions 2,000,000 / 2,500,000 / 3,000,000 / 3,500,000 / 4,000,000. Total 7.5 kb.
- **Inverted 3 kb repeat.** One forward copy and one reverse-complement copy of a 3,000 bp unit, at positions 4,500,000 and 4,700,000. Total 6 kb.

The construction script `data/build_reference_with_repeats.py` is provided. A 50x ONT-style simulation produces ~25,000 reads (mean 15 kb) and runs through Flye in ~3-5 minutes.

### Bandage CLI

```bash
Bandage image \
    results/flye_5mb/assembly_graph.gfa \
    results/bandage_graph.png \
    --height 1200
```

The output PNG can be rendered on a headless server. For the GUI inspection (optional; do this in the report write-up), open the GFA in the Bandage GUI and look at the topology around each repeat region.

---

## Acceptance criteria

- [ ] `bandage_inspect.py` runs end to end with `python bandage_inspect.py --reference data/reference_5mb_with_repeats.fasta --out-dir results --seed 42`.
- [ ] `results/flye_5mb/assembly.fasta` and `results/flye_5mb/assembly_graph.gfa` exist.
- [ ] `results/bandage_graph.png` exists (or a warning if `Bandage` is not on the PATH).
- [ ] `results/gfa_summary.tsv` lists every segment with its length, coverage, in-degree, out-degree.
- [ ] `results/topology_report.md` is a one-page report that:
  - States the overall topology classification.
  - For each of the three deliberate repeats, names whether Flye resolved it (separate segments in the GFA) or collapsed it (single segment with anomalously high coverage).
  - Includes the Bandage PNG (embedded via `![graph](bandage_graph.png)`).
- [ ] `results/run-info.json` records: Flye version, badread version, Bandage version (or `unknown` if not present), seed, run date.

---

## Expected output

On the demo 5 Mb-with-repeats reference at 50x:

- **Tandem 10x repeat (200 bp x 10):** Almost always collapsed by Flye. Visible as a single segment of length ~200 bp with coverage 10x the genome median. The Bandage view shows a small node with a self-loop or with two flanking edges to the upstream and downstream context.
- **Dispersed 5x repeat (1,500 bp x 5):** Sometimes resolved (if the flanking reads are long enough), sometimes collapsed. Visible as either five separate segments (each ~1.5 kb, normal coverage) or a single repeat segment with 5x coverage.
- **Inverted 3 kb repeat (3 kb x 2):** Usually resolved at 15 kb mean read length (the read spans the repeat plus the flank). Visible as two separate segments with normal coverage and an inversion mark in the GFA's `L` lines (one of the two endpoints is in the opposite orientation).

The overall N50 is typically 1-3 Mb (Flye breaks the assembly at the collapsed-tandem-repeat boundary; the assembly fragments into 2-4 contigs covering ~95-99% of the reference).

The topology classification is `tangle` (because of the collapsed repeats) or `linear+tangle` (the bulk is linear; the repeats form a tangle).

---

## The GFA parser

```python
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class GfaSegment:
    seg_id: str
    length: int
    coverage: float
    in_degree: int = 0
    out_degree: int = 0


def parse_gfa(gfa_path: Path) -> tuple[list[GfaSegment], int, int]:
    """Parse a GFA. Returns (segments, n_links, n_paths)."""
    segments: dict[str, GfaSegment] = {}
    n_links: int = 0
    n_paths: int = 0
    coverage_tag = re.compile(r"dp:f:([\d.]+)")
    with gfa_path.open() as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts: list[str] = line.rstrip("\n").split("\t")
                seg_id: str = parts[1]
                seq_len: int = len(parts[2]) if parts[2] != "*" else 0
                cov_match = coverage_tag.search(line)
                cov: float = float(cov_match.group(1)) if cov_match else 0.0
                segments[seg_id] = GfaSegment(
                    seg_id=seg_id,
                    length=seq_len,
                    coverage=cov,
                )
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    a: str = parts[1]
                    b: str = parts[3]
                    if a in segments:
                        segments[a].out_degree += 1
                    if b in segments:
                        segments[b].in_degree += 1
                n_links += 1
            elif line.startswith("P\t"):
                n_paths += 1
    return list(segments.values()), n_links, n_paths


def classify_topology(
    segments: list[GfaSegment],
    n_links: int,
) -> str:
    """Return one of: 'linear', 'circular', 'bubble', 'tangle', 'mixed'.

    Heuristics:
    - linear: 0 links and 1+ segments.
    - circular: 1 link per segment, all self-loops.
    - bubble: at most one segment has in_degree > 1 OR out_degree > 1.
    - tangle: multiple segments with in_degree > 1 OR out_degree > 1.
    """
    if not segments:
        return "empty"
    if n_links == 0:
        return "linear"
    branching: list[GfaSegment] = [
        s for s in segments if s.in_degree > 1 or s.out_degree > 1
    ]
    self_loops: int = sum(
        1 for s in segments if s.in_degree >= 1 and s.out_degree >= 1
    )
    if not branching and self_loops == len(segments):
        return "circular"
    if len(branching) <= 2:
        return "bubble"
    return "tangle"
```

---

## Write-up

In `challenge-02/results/topology_report.md`, answer in ~400 words:

1. **What is the overall topology of the assembly graph?** Name the classification and the segment counts.
2. **For each of the three deliberate repeats (tandem 10x, dispersed 5x, inverted 3 kb), state whether Flye resolved or collapsed it.** Use the coverage on the candidate segment as the discriminator: coverage > 1.5x the median is a collapse.
3. **If you had to assemble this genome again with the *same* reads, what would you change?** Hint: longer reads (simulated ONT ultra-long with `--length 100000,50000`) or a different assembler (Canu with `corOutCoverage=80`).
4. **Reading the Bandage PNG, identify two graph features that the FASTA alone does not reveal.** Be specific: position in the graph, segment IDs, what it tells you about the assembly.

---

## Stretch goals (optional)

- Re-run Flye with `--nano-raw` instead of `--nano-hq` and see how the graph topology shifts. The looser overlap thresholds usually produce more tangles.
- Run Canu on the same reads and compare the Canu GFA topology to the Flye topology. Canu is typically more conservative at repeats; expect Canu to resolve more repeats but produce smaller individual contigs in the resolved regions.
- Use `bandage info` (the CLI subcommand for graph statistics) to extract per-segment and per-edge data, then plot the segment-length-vs-coverage scatter as a PNG. The collapsed repeats appear as outliers at the top of the scatter (short segments with anomalously high coverage).

---

## What to commit

```
challenge-02/
    README.md
    bandage_inspect.py
    data/{reference_5mb_with_repeats.fasta, build_reference_with_repeats.py}
    results/
        reads.fastq
        flye_5mb/{assembly.fasta, assembly_graph.gfa, assembly_info.txt, flye.log}
        bandage_graph.png
        gfa_summary.tsv
        topology_report.md
        run-info.json
```

Gitignore: Flye intermediate stages, the badread quality-distribution PNGs, and any Bandage GUI screenshots (the CLI-rendered PNG is the reproducible artifact).

Commit message like `challenge-02: Bandage inspect on 5Mb repeat-rich sim, Flye resolved inverted 3kb, collapsed tandem 10x and 1 dispersed`.
