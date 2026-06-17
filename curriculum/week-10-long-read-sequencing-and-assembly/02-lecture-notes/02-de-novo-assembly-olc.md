# Lecture 2 — De Novo Assembly and the OLC Graph

> **Reproducibility note.** Long-read assemblers are heuristics for an NP-hard problem (the simultaneous overlap-and-layout of N reads). They are mostly deterministic given the same input and version, but multi-threading can introduce order-of-merge wobble on equivalent overlaps, and a small change in the input (a single read added or removed) can change the contig boundaries by several kilobases. Pin the tool, pin the version, pin the input mode flag (`--nano-hq` vs `--pacbio-hifi`), pin the thread count, and the assembly becomes reproducible. Skip any of those and you have a moving target.

> **Duration:** ~2.5 hours of reading + a Flye run on the simulated dataset.
> **Outcome:** You can describe the overlap-layout-consensus (OLC) paradigm, name Flye, Canu, and Hifiasm as the three canonical free long-read assemblers, state when each is preferred, describe what a repeat graph is and why it differs from the de Bruijn graph used for short reads, and call Flye from Python with `subprocess.run(..., check=True)`.

If you only remember one thing from this lecture, remember this:

> **Long-read assembly is an Overlap-Layout-Consensus problem solved on a graph data structure. Flye builds a repeat graph; Canu builds a best-overlap graph; Hifiasm builds a phased assembly graph. The choice of tool depends on the input chemistry (ONT vs HiFi) and the assembly goal (haploid contigs vs haplotype-resolved diploid). All three are free, open source, and produce a FASTA of contigs plus a GFA of the graph.**

Lecture 1 left us with a cleaned FASTQ of long reads. Lecture 2 takes the FASTQ to a FASTA of contigs.

---

## 1. Where we are in the pipeline

```
FASTQ (cleaned long reads) ->
        OLC assembler (Flye for ONT --nano-hq; Hifiasm for HiFi; Canu as the slower alternative) ->
        FASTA (M contigs)
        + GFA (the assembly graph)
        + assembly_info.txt (per-contig length, coverage, circularity).
```

Lecture 2 covers the middle box — the OLC assembler. Lecture 3 takes the FASTA and the GFA through polishing and QC.

The unit of work for an assembler is the **overlap**: a pair of reads that share a long enough end-to-end match for the assembler to believe they originated from the same region of the source DNA. With N reads the naive number of pairwise overlap candidates is N(N-1)/2, which is impractical for N > ~10,000. Every modern long-read assembler reduces the search space with one of two approaches:

- **Hash-based seed-and-extend** (Flye, Canu, minimap2): pick a small set of "seed" k-mers in each read, store them in a hash table, find read pairs that share enough seeds, and extend the seed matches to full overlaps.
- **Sparse-graph k-mer matching** (Hifiasm): use a much sparser k-mer fingerprint and rely on the higher per-base accuracy of HiFi reads to make the matches reliable.

Both approaches reduce the practical compute from O(N^2) to roughly O(N log N) per read, which is the difference between "runs on a laptop" and "needs a cluster."

---

## 2. The OLC paradigm in one diagram

The classical overlap-layout-consensus pipeline has three stages:

```
            +------------+         +------------+         +-----------+
reads --->  | (1) OVERLAP|  --->   | (2) LAYOUT |  --->   |(3)CONSEN- |  ---> contigs
            +------------+         +------------+         |   SUS     |
                                                          +-----------+
```

- **(1) Overlap.** For each pair of reads, find the alignment that maximizes a scoring function under a minimum-length-and-identity threshold. Output: an "overlap set" — a list of edges between reads with offsets and identity scores.
- **(2) Layout.** Treat the overlap set as a graph (nodes = reads, edges = overlaps) and clean it: remove transitive edges, remove low-confidence edges, identify bubbles (alternative paths through homologous regions, which usually mean sequencing errors or heterozygosity), and identify tangles (regions where the graph branches because of repeats). Output: a smaller graph whose paths are candidate contigs.
- **(3) Consensus.** For each path in the cleaned graph, walk the path read-by-read and call a per-position consensus from the overlapping reads. Output: a FASTA of contigs.

The three stages are mostly independent in modern assemblers (the consensus does not change which reads overlap), so an assembler is more or less a pipeline of three heuristics chained together. The choice of heuristic at each stage is where the assemblers differ.

---

## 3. The repeat graph (Flye) vs the de Bruijn graph (short-read) vs the best-overlap graph (Canu)

The "graph" in OLC is a small zoo of data structures. The three that matter for long-read assembly:

### De Bruijn graph (short reads only)

Used by short-read assemblers (SPAdes, Velvet) and not directly by long-read assemblers. The de Bruijn graph has one node per k-mer (k typically 21-65) and one edge per consecutive k-mer pair. With short reads it scales beautifully because the k-mer count is bounded by the genome size. With long reads it fails: the k-mer error rate is too high (a 5% per-base error rate with k=21 means most k-mers have an error somewhere, and the graph becomes mostly noise nodes).

We mention de Bruijn graphs here only to contrast with the long-read graphs below; you do not build them in Week 10.

### Best-overlap graph (Canu)

Canu's representation. Each read is a node. For each read, keep only the *single* best overlap on each end (in the orientation that extends the read). The result is a graph where each node has degree at most 2 (one best overlap on each end). The graph is a set of paths and cycles; each path is a candidate contig.

The best-overlap-graph is conservative: a read with two roughly-equal candidate overlaps on one end (a "repeat boundary") is conservatively represented as a single edge to whichever overlap wins by a hair, with the loser available as an alternate path that the post-processing can choose. This works well for highly repetitive genomes (the original Canu motivation was the *Lemur catta* and *Tarsius syrichta* primate genomes) but is slow because the per-read overlap computation is expensive.

### Repeat graph (Flye)

Flye's representation. Each repeated region in the reads is collapsed into a single "repeat node" and the non-repeated flanks are represented as edges that enter and exit the repeat node. The output looks like a string graph but with the repeat regions explicitly tagged. The repeat-aware structure means Flye can recognize when an overlap is repeat-internal (and therefore unreliable) versus boundary-spanning (and therefore high-confidence), and the layout stage can prefer the boundary-spanning overlaps.

The repeat graph is the technical contribution of the Flye paper (Kolmogorov et al. 2019; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699608/>); §2 of the paper is the canonical reference. The intuition is: repeats are the part of the genome where OLC fails most often, so represent them explicitly and the failure becomes easier to diagnose.

For a clean bacterial genome with few repeats the repeat graph and the best-overlap graph produce essentially the same contigs. The differences appear on repeat-rich inputs (the rRNA operon copies in *E. coli*; the CRISPR arrays in some Archaea; large insertion-sequence families). Flye usually handles these better than Canu, in less compute time.

### Phased assembly graph (Hifiasm)

Hifiasm's representation. Specialized for diploid organisms: each pair of haplotype-distinct contigs (the maternal and paternal copies of a region) is kept separate in the graph, rather than collapsed into a single consensus contig. The output is two contig sets: a **primary** set (the longer haplotype per region) and an **alternate** set (the other haplotype).

For haploid organisms (most bacteria) the primary set is "the" assembly and the alternate set is empty. For diploid organisms (humans, plants) both sets are useful: the primary is the conventional assembly, the alternate is the haplotype-resolved counterpart. This is the main reason Hifiasm has displaced Canu for human and plant HiFi assemblies — the haplotype resolution comes for free.

---

## 4. Flye in detail

Citation: Kolmogorov M, Yuan J, Lin Y, Pevzner PA. "Assembly of long, error-prone reads using repeat graphs." *Nat Biotechnol* 37:540 (2019). Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699608/>. Tool docs at <https://github.com/fenderglass/Flye>.

Flye is the Week 10 default for ONT input. It is a C++ program with a Python front-end and a permissive license (BSD-3-Clause).

### The Flye pipeline

Internally Flye runs a sequence of stages, each writing into a numbered subdirectory of the output:

1. `00-assembly/` — Builds a repeat graph from the raw reads. Uses a custom k-mer-frequency-aware seed-and-extend algorithm to find overlaps without an explicit pairwise comparison stage.
2. `10-consensus/` — Computes a first-pass consensus per contig using the read overlaps.
3. `20-repeat/` — Resolves repeat regions: for each repeat node in the graph, count the read coverage and the read-flank patterns to estimate the repeat copy number; if it can be resolved (reads spanning two copies are present), unfold the node; otherwise leave it as a single "compressed" representation.
4. `30-contigger/` — Walks the resolved graph and emits contig sequences.
5. `40-polishing/` — Polishes the contig sequences using the read overlaps (this is Flye's internal polish; Medaka in Lecture 3 is a separate, more thorough polish).

The output files of interest are:

- `assembly.fasta` — the contigs.
- `assembly_graph.gfa` — the final repeat graph (after resolution) in GFA format; viewable in Bandage.
- `assembly_info.txt` — per-contig summary table: `seq_name`, `length`, `coverage`, `circular`, `repeat`, `mult`, `alt_group`, `graph_path`.
- `flye.log` — the full run log with timestamps, k-mer statistics, repeat resolution decisions, and a final summary.

### The canonical Flye call

```bash
flye \
    --nano-hq reads.fastq \
    --genome-size 1m \
    --out-dir flye_out \
    --threads 4 \
    --iterations 1
```

- `--nano-hq` — Modern Dorado-SUP-basecalled R10.4 reads. Use `--nano-raw` only for legacy Guppy-basecalled R9.4 reads; the parameter affects overlap thresholds.
- `--genome-size 1m` — Estimated genome size; used for coverage normalization. Off by 2-3x is fine; off by 10x degrades the assembly.
- `--threads 4` — Pin the thread count for reproducibility (see the caveat in §10).
- `--iterations 1` — Number of internal polishing rounds. Default 1 is usually sufficient; the external Medaka polish in Lecture 3 picks up the rest.

For the 50x simulated nanopore dataset on the 1 Mb reference, Flye produces a typical output of:

```
contig_1    length=1038203    coverage=49    circular=Y    repeat=N    multiplicity=1
```

A single contig, ~1.04 Mb, circular (Flye detected the circularity from the read overlap pattern), and the coverage matches the input 50x. The N50 is 1,038,203; the L50 is 1.

### Calling Flye from Python

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def run_flye(
    reads_fastq: Path,
    out_dir: Path,
    genome_size: str = "1m",
    input_mode: str = "--nano-hq",
    threads: int = 4,
    iterations: int = 1,
) -> Path:
    """Run Flye. Returns the path to the assembly FASTA.

    Raises subprocess.CalledProcessError if Flye fails or is not on the PATH.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "flye",
        input_mode, str(reads_fastq),
        "--genome-size", genome_size,
        "--out-dir", str(out_dir),
        "--threads", str(threads),
        "--iterations", str(iterations),
    ]
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return out_dir / "assembly.fasta"
```

The function is straightforward: build the command, run it, return the expected output path. Flye writes its log to `out_dir/flye.log` automatically; we do not need to capture it on stdout.

### Reading the assembly_info.txt

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class ContigInfo:
    seq_name: str
    length: int
    coverage: int
    circular: bool
    repeat: bool
    multiplicity: int


def parse_flye_assembly_info(info_path: Path) -> list[ContigInfo]:
    """Parse Flye's assembly_info.txt into a list of ContigInfo records.

    Skips the header line (which starts with '#') and any blank lines.
    """
    rows: list[ContigInfo] = []
    with info_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts: list[str] = line.split("\t")
            rows.append(ContigInfo(
                seq_name=parts[0],
                length=int(parts[1]),
                coverage=int(parts[2]),
                circular=(parts[3] == "Y"),
                repeat=(parts[4] == "Y"),
                multiplicity=int(parts[5]),
            ))
    return rows
```

The `assembly_info.txt` is a TSV with one row per contig and one of the more useful summaries the assembler emits. The `circular=Y` flag is the load-bearing field for bacterial genomes: a closed circular chromosome should always have `circular=Y` on its primary contig.

---

## 5. Canu in detail

Citation: Koren S, Walenz BP, Berlin K, Miller JR, Bergman NH, Phillippy AM. "Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation." *Genome Research* 27:722 (2017). Free at <https://genome.cshlp.org/content/27/5/722> or <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411767/>. Tool docs at <https://canu.readthedocs.io/>.

Canu is the older OLC assembler and the direct descendant of Celera Assembler (the assembler used for the original human genome). Its design priorities are accuracy on repetitive genomes and clean handling of the assembly graph; the cost is wall-clock time.

### The Canu pipeline

Canu runs three sub-pipelines, each a separate executable:

1. **Correction** (`canu -correct`): for each read, identify the high-quality overlapping reads and use them to error-correct the original read. The output is a "corrected" FASTQ in which the per-base error rate is reduced to <1%. This stage is the most compute-intensive in Canu and the one Flye and Hifiasm explicitly try to skip.
2. **Trimming** (`canu -trim`): for each corrected read, identify the high-confidence span (the inner region of the read where flanking adapters or low-quality ends have been removed) and emit only that.
3. **Assembly** (`canu -assemble`): build the best-overlap graph from the trimmed corrected reads, walk the graph, and emit contigs.

The wrapper `canu` (without an explicit subcommand) runs all three in sequence.

### The canonical Canu call

```bash
canu \
    -p ecoli \
    -d canu_out \
    genomeSize=1m \
    -nanopore reads.fastq \
    useGrid=false \
    maxThreads=4 \
    maxMemory=16g
```

For the same 50x / 1 Mb simulated dataset, Canu takes ~5-10 minutes on a 4-core laptop (compared to ~30-60 seconds for Flye). The output is `canu_out/ecoli.contigs.fasta` (typically 1-2 contigs on a clean dataset; very similar to Flye on the same input). The `canu_out/ecoli.report` file has the run statistics; the `canu_out/ecoli.contigs.gfa` is the graph.

### When to prefer Canu over Flye

The conventional wisdom (as of 2024-2025) is: **default to Flye; use Canu only when you specifically need its features**. The features:

- **More conservative repeat handling.** Canu's best-overlap-graph is more conservative on repeats than Flye's repeat graph; on heavily repetitive genomes (Archaea with many CRISPR arrays; some plant chloroplast genomes) Canu sometimes resolves repeats Flye collapses.
- **Better-documented graph.** The Canu GFA is more granular than the Flye GFA; if you intend to manually inspect the graph in Bandage and rewire it, Canu's representation is easier to work with.
- **Production-tested workflow.** Canu has been the backbone of many published assemblies since 2017; its behavior on edge-case inputs is well-characterized.

The Canu cost is wall-clock time (~10x Flye) and memory (Canu wants 16-32 GB even on small genomes because of the correction stage). For most Week 10 / mini-project inputs Flye is the right default and Canu is the comparison; we run Canu in the exercises as a "second opinion" rather than the primary tool.

---

## 6. Hifiasm in detail

Citation: Cheng H, Concepcion GT, Feng X, Zhang H, Li H. "Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm." *Nat Methods* 18:170 (2021). Free at <https://www.nature.com/articles/s41592-020-01056-5> or <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8005227/>. Tool docs at <https://github.com/chhylp123/hifiasm>.

Hifiasm is specialized for PacBio HiFi input. The conceptual reframe is:

- **HiFi reads are already accurate.** A 15-20 kb read at QV30-50 does not need a correction stage; the assembler can use the raw reads directly.
- **HiFi reads are length-limited.** They are shorter than the longest ONT ultra-long reads, so the assembler cannot rely on ultra-long-span-the-repeat tricks.
- **Diploid organisms benefit from haplotype resolution.** Two HiFi reads from the same region of two different haplotypes are nearly identical (~99.9% identity) but differ at the SNPs between them; an assembler that resolves the two as separate contigs preserves the heterozygosity.

Hifiasm's algorithm is a string graph (similar to a best-overlap graph) built on k-mers, with an additional **phasing** stage that uses the SNP differences between reads to assign each read to a haplotype.

### The Hifiasm pipeline

1. **All-vs-all overlap.** Use a sparse minimizer-based scheme to find candidate overlaps among the HiFi reads.
2. **String graph construction.** Build a string graph from the high-quality overlaps; clean it by removing transitive edges and low-confidence tips.
3. **Haplotype phasing.** For each "bubble" in the string graph (a divergence corresponding to a haplotype-distinct region), assign the two branches to haplotypes using the SNP markers shared with adjacent unitigs.
4. **Output two contig sets.** The primary contigs (`*.bp.p_ctg.gfa`) and the alternate contigs (`*.bp.a_ctg.gfa`).

For haploid input (the Week 10 default) Hifiasm runs the same pipeline; the alternate contig set is empty.

### The canonical Hifiasm call

```bash
hifiasm \
    -o hifiasm_out/ecoli \
    -t 4 \
    --n-hap 1 \
    reads.hifi.fastq
```

- `-o hifiasm_out/ecoli` — output prefix; the primary contigs are at `hifiasm_out/ecoli.bp.p_ctg.gfa`.
- `-t 4` — thread count.
- `--n-hap 1` — haploid expectation (use `--n-hap 2` for diploid).

The primary output is a GFA, not a FASTA. To convert:

```bash
awk '/^S/ {print ">"$2"\n"$3}' hifiasm_out/ecoli.bp.p_ctg.gfa > hifiasm_out/ecoli.fasta
```

(The GFA's `S` lines are segments; `$2` is the segment ID and `$3` is the segment sequence.)

For a 50x simulated HiFi dataset on a 1 Mb reference, Hifiasm runs in ~30 seconds and produces a single contig of ~1.04 Mb. The N50 is the contig length; the L50 is 1.

### When to prefer Hifiasm

- **HiFi input.** Always Hifiasm. Flye works on HiFi (with `--pacbio-hifi`) but is slower and emits unphased contigs.
- **Diploid organisms.** Always Hifiasm. The haplotype-resolved output is the entire reason Hifiasm exists.
- **Mixed HiFi + ONT input.** Hifiasm has experimental support for ONT ultra-long as a "scaffolding" track (`--ul reads_ont.fastq`); the combined output is more contiguous than HiFi alone. This is the approach Verkko uses for the human reference assembly.

For ONT-only input Hifiasm is not appropriate; the algorithm assumes per-read accuracy that ONT does not provide.

---

## 7. Choosing the right assembler for the input

A flowchart:

```
What is the input platform?
    ONT (R10.4.1, Dorado SUP)         -> Flye --nano-hq
    ONT (R9.4.1, Guppy SUP)           -> Flye --nano-raw
    PacBio HiFi                       -> Hifiasm
    PacBio CLR (pre-HiFi)             -> Canu -pacbio (legacy)
    Mixed ONT + HiFi                  -> Verkko (out of scope for Week 10)

What is the genome property?
    Bacterial (1-10 Mb, low repeat)   -> any; Flye fastest
    Archaeal (1-5 Mb, CRISPR repeats) -> Canu more conservative
    Small eukaryote (10-100 Mb, diploid) -> Hifiasm on HiFi; Verkko on hybrid
    Large eukaryote (> 100 Mb)        -> Verkko on hybrid; not Week 10

What is the analysis goal?
    Quick draft                       -> Flye
    Publication-quality haploid       -> Flye + Medaka, or Hifiasm
    Haplotype-resolved diploid        -> Hifiasm
    Maximum interpretability of graph -> Canu
```

The Week 10 defaults are:

- **Exercises and mini-project:** Flye --nano-hq on simulated nanopore reads.
- **Challenge 1:** Hifiasm on simulated HiFi reads.
- **Comparison in exercise 1:** optional Canu run.

---

## 8. Reading the GFA: bubbles, tangles, and clean linear paths

The Graphical Fragment Assembly (GFA) format is the standard for assembly graphs. The spec is at <https://gfa-spec.github.io/GFA-spec/GFA1.html>. Each line is one record:

- `H` — header (the version line: `H VN:Z:1.0`).
- `S <id> <sequence> [<tags>]` — a segment (a sequence node).
- `L <from> <from_orient> <to> <to_orient> <overlap_cigar> [<tags>]` — a link (an edge).
- `P <name> <segments> <overlaps>` — a path (a named walk through segments).

A clean linear assembly looks like:

```
H VN:Z:1.0
S contig_1 ACGTACGT...
P contig_1 contig_1+ *
```

One `S` line, one trivial `P` line, no `L` lines. The contig is the only path; there are no branching points.

A clean circular assembly looks like:

```
H VN:Z:1.0
S contig_1 ACGTACGT...
L contig_1 + contig_1 + 0M
P contig_1 contig_1+ *
```

One `S` line, one `L` line linking the segment to itself (the circular closure), one trivial `P` line.

A heterozygous-bubble assembly (typical of a diploid input):

```
H VN:Z:1.0
S contig_1a ACGTACGT...
S contig_1b ACGTAAGT...    (one SNP different)
L upstream + contig_1a + 0M
L upstream + contig_1b + 0M
L contig_1a + downstream + 0M
L contig_1b + downstream + 0M
```

Two segments at the bubble, two incoming edges, two outgoing edges. Either branch is a valid path; the assembler emits one as primary and the other as alternate (Hifiasm) or collapses them (Flye unless `--keep-haplotypes`).

A tangle (typical of an unresolved repeat):

```
S repeat_node ACGTAC...     (the repeat unit, multiplicity > 1)
L flank_1a + repeat_node + 0M
L flank_2a + repeat_node + 0M
L repeat_node + flank_1b + 0M
L repeat_node + flank_2b + 0M
```

One segment with multiple incoming and multiple outgoing edges; the assembler could not determine which flank-pair was the correct walk. The output FASTA breaks at the tangle (the assembler emits the flanks as separate contigs and leaves the repeat as a separate fragment).

The art of looking at a GFA in Bandage is recognizing these topologies. Lecture 3 ties this to the QC stage; Challenge 2 puts it into practice.

---

## 9. Comparing assemblies: when to declare them "the same"

Two assemblies of the same input that differ in any of:

- The contig count.
- The total assembled length (by more than ~1% — under 1% is consensus-stage rounding).
- The N50 (by more than ~5%).
- The presence or absence of a circular contig flag.
- The number of `L` lines in the GFA.

…are different assemblies that warrant investigation. Differences below those thresholds are typically thread-merge wobble or consensus-stage rounding and not worth chasing.

For the mini-project, the conventional "is my assembly stable?" test is to run Flye twice on the same FASTQ with the same thread count and `diff` the output FASTA. They should be byte-identical or differ on at most ~50 bp out of ~1 Mb (rounding in the per-position consensus on a few base-equivalent positions). If they differ by more, the assembler has a non-determinism source you have not pinned.

---

## 10. Standard failure modes at the assembly stage

- **Insufficient coverage.** Below ~25-30x for ONT or ~15-20x for HiFi, the assembler cannot find enough overlaps to bridge most regions; the assembly fragments into many small contigs. The fix is more reads.
- **Unresolved repeats.** A repeat region wider than the read N50 cannot be threaded; the assembler breaks the assembly at the repeat boundaries. The fix is longer reads (ONT ultra-long, or Hi-C scaffolding after assembly).
- **Overcollapsed repeats.** N tandem copies of a repeat unit are collapsed into one contig. Visible as a coverage anomaly: the collapsed contig has K-times-higher read coverage than the genome average. The mini-project pipeline checks the per-contig coverage and flags any contig with coverage > 1.5x the genome median.
- **Adapter contamination misassembly.** Adapter sequences inside reads cause the assembler to join two distant regions on a phantom overlap. The fix is an adapter-trimming pre-pass with `Trimmomatic` or `Porechop`.
- **Chimeric reads creating chimeric contigs.** A 1-3% chimera rate in ONT reads can produce false joins in the contigger; the result is a chimeric contig with two regions of obviously-mismatched coverage. The mini-project pipeline flags coverage discontinuities within a single contig.
- **Wrong input mode.** Running Flye `--nano-raw` on R10.4.1 reads (or `--nano-hq` on R9.4.1 reads) produces a measurably worse assembly because the overlap thresholds are tuned to the input chemistry's error rate. Pin the input mode in the run-info JSON.
- **Thread-merge wobble.** Running Flye `--threads 4` twice produces byte-identical output; running `--threads 4` and `--threads 8` on the same input produces assemblies that differ on ~50-200 bp out of ~1 Mb. The differences are small but reproducibility requires pinning the thread count.

The mini-project write-up names which of these the assembly can plausibly suffer from. For the simulated 1 Mb / 50x case the answer is short: overcollapsed repeats (the mini-project's synthetic 5 Mb reference deliberately includes a tandem repeat to test for this) and thread-merge wobble are the realistic risks; coverage and adapter-contamination are non-issues by construction.

---

## 11. Lecture 2 wrap-up

We named three long-read assemblers (Flye, Canu, Hifiasm), described what makes each different (the graph data structure each uses), stated when each is preferred (ONT vs HiFi; haploid vs diploid; speed vs interpretability), walked through the GFA file format and the topologies (linear, circular, bubble, tangle) it can express, and listed the standard failure modes at the assembly stage.

Lecture 3 takes the FASTA + GFA the assembler emits and runs the polishing (Medaka for ONT; usually nothing for HiFi) and QC (N50 / L50 by hand; BUSCO for gene content; Bandage for graph topology). Exercise 1 ties Lectures 1 and 2 together: simulate reads with badread, run Flye, parse the output, and verify the assembly_info.txt.

---

## Lecture 2 sanity check (15 minutes)

If you have `badread` and `flye` installed on the local PATH:

```bash
# 1. Generate a 100 kb random reference (Flye runs are faster on a smaller reference).
python -c "
import random
random.seed(1)
seq = ''.join(random.choice('ACGT') for _ in range(100000))
print('>chr1')
print(seq)
" > /tmp/ref.fasta

# 2. Simulate 50x of nanopore reads.
badread simulate --reference /tmp/ref.fasta --quantity 50x \
    --length 15000,13000 --identity 95,3,99 \
    --error_model nanopore2023 --qscore_model nanopore2023 \
    --seed 42 > /tmp/reads.fastq

# 3. Run Flye.
flye --nano-hq /tmp/reads.fastq --genome-size 100k \
    --out-dir /tmp/flye_out --threads 2

# 4. Look at the assembly summary.
cat /tmp/flye_out/assembly_info.txt
```

Expected output: a single contig in the `assembly_info.txt` with length ~100,000 (within ~5%), coverage ~50, circular=Y (Flye usually detects circularity even on synthetic input because the read overlaps wrap around). The `flye.log` should end with a "Total length: 100,xxx" summary near the bottom.

If the sanity check passes, you are ready for Lecture 3.

---

## Further reading

- Kolmogorov M et al. "Assembly of long, error-prone reads using repeat graphs." *Nat Biotechnol* 37:540 (2019). The Flye paper.
- Koren S et al. "Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation." *Genome Res* 27:722 (2017). The Canu paper.
- Cheng H et al. "Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm." *Nat Methods* 18:170 (2021). The Hifiasm paper.
- Myers EW. "The fragment assembly string graph." *Bioinformatics* 21 Suppl 2:ii79 (2005). The original string-graph paper that the modern OLC graphs descend from.
- Li H. "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics* 34:3094 (2018). The long-read mapper underneath Hifiasm's overlap stage.

Continue to [Lecture 3 — Polishing and Assembly QC](./03-polishing-and-qc.md).
