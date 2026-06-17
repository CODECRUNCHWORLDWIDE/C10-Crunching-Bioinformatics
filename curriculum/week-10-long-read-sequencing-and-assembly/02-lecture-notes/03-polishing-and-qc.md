# Lecture 3 — Polishing and Assembly QC

> **Reproducibility note.** Polishing is a neural-network-driven correction step. Medaka and DeepConsensus are deterministic given a fixed model and a fixed input, but the model number is load-bearing: a Medaka model trained for R9.4.1 Guppy SUP, applied to R10.4.1 Dorado SUP reads, produces a *worse* polish than no polish at all because the model expects a different signal distribution. Always pin the polish model in the run-info JSON. The same logic applies to BUSCO: pin the BUSCO version, the lineage dataset name, and the dataset creation date; an updated `bacteria_odb10` can score the same assembly differently.

> **Duration:** ~2 hours of reading + a Medaka and a BUSCO run on the simulated dataset.
> **Outcome:** You can describe what Medaka and DeepConsensus do (neural-network polishers operating at different stages of the pipeline), compute N50 and L50 by hand from a contig FASTA, run BUSCO on the polished assembly and read its summary, classify a Bandage graph view (linear / circular / bubble / tangle), and write a QC report in Markdown.

If you only remember one thing from this lecture, remember this:

> **The polishing stage corrects systematic per-base errors that the assembler's internal consensus cannot. The QC stage measures the assembly along two orthogonal axes: contiguity (N50, L50, total length) and content correctness (BUSCO gene completeness). Both axes can lie: a high N50 can come from collapsed repeats; a high BUSCO can co-exist with a wrong-sized assembly. The honest QC report names both numbers, the limit of each, and the failure modes the assembly could plausibly suffer from.**

Lecture 2 left us with a FASTA of contigs and a GFA of the assembly graph. Lecture 3 polishes the contigs and characterizes the assembly's quality.

---

## 1. Where we are in the pipeline

```
FASTA (Flye contigs, ~99% accurate per base) ->
        Medaka (or DeepConsensus upstream of the assembly for HiFi) ->
        FASTA (polished, ~99.9% accurate per base) ->
        QC stage:
            asmstats (N50, L50, total length, contig count, GC fraction) ->
            BUSCO (gene-content completeness against a lineage dataset) ->
            Bandage (graph topology inspection from the GFA) ->
            QV computation (if a reference is available) ->
        QC report (Markdown) +
        run-info.json (versions, models, parameters, run date).
```

Lecture 3 covers the polish-and-QC tail. The mini-project ties Lectures 1, 2, and 3 together into a one-command pipeline.

---

## 2. Why polishing matters

A Flye assembly of 50x nanopore reads on a 1 Mb reference typically has:

- ~99% per-base accuracy (~10 errors per 1,000 bp, QV20).
- Errors concentrated at homopolymer runs (a stretch of N A's may be called as N-1 or N+1 A's) and at methylated bases (the basecaller's model is less confident on the modified-base sites).
- A few hundred ambiguous-base calls (`N` in the contig) at low-coverage positions.

A QV20 assembly is fine for many uses (large-scale structural comparison, gene prediction down to 80-90% accuracy), but it is **not publishable** as a reference genome. For publication you want QV35+ (about one error per 3,000 bp; or one error every 10,000 bp on average for QV40), which is what a single round of Medaka polishing on Flye-on-nanopore typically reaches.

Polishing has two flavors:

- **Post-assembly polishing.** Take the draft contigs and the original reads; re-align the reads to the draft; for each position, call a per-position consensus from the read pileup, optionally using a neural network. **Medaka** is the ONT-specialized tool here.
- **Upstream consensus polishing.** Before assembly, take the raw subreads (HiFi only) and compute a higher-quality consensus per ZMW; then assemble. **DeepConsensus** is the HiFi-specialized tool here.

The two flavors are not interchangeable: Medaka is trained on the ONT basecaller's error profile and applied to assembled contigs; DeepConsensus is trained on the HiFi subread error profile and applied to per-read consensus.

---

## 3. Medaka in detail

Tool: <https://github.com/nanoporetech/medaka>. License: Mozilla Public License 2.0 (open source).

Medaka is the Oxford Nanopore Technologies polisher. It is a neural-network-based per-position consensus caller; the input is (a) a draft assembly FASTA, (b) the original nanopore reads, and (c) a Medaka model file matched to the basecaller model that produced the reads. The output is a polished FASTA where each base has been re-called from the read pileup using the model's per-position probabilities.

### The Medaka pipeline

Internally Medaka runs:

1. **Read alignment.** `minimap2 -a -x map-ont reads.fastq draft.fasta` to produce a BAM of reads aligned to the draft.
2. **Pileup featurization.** For each draft position, extract a fixed-size feature vector from the read pileup (the bases under the pileup, their qualities, the local alignment context).
3. **Neural-network consensus.** Pass the feature vectors through a recurrent neural network (older models) or a transformer (newer models) trained on a large nanopore consensus training set. The output is a per-position probability over A / C / G / T / `-` (gap, meaning "delete this base") / new bases (insertions).
4. **Polished sequence emission.** Apply the per-position predictions to the draft to produce the polished FASTA.

### The canonical Medaka call

```bash
medaka_consensus \
    -i reads.fastq \
    -d flye_out/assembly.fasta \
    -o medaka_out \
    -m r1041_e82_400bps_sup_v4.3.0 \
    -t 4
```

For the 1 Mb / 50x simulated dataset, Medaka takes ~2-5 minutes on a 4-core CPU. The output is `medaka_out/consensus.fasta` (the polished FASTA). On simulated data, the QV typically improves from ~20 (raw Flye) to ~35-40 (after Medaka).

### Model selection (the load-bearing parameter)

Medaka's model list is maintained at <https://github.com/nanoporetech/medaka/blob/master/medaka/options.py>. The recipe is:

- `r1041_e82_400bps_sup_v4.3.0` — current default for R10.4.1 + Dorado SUP + Dorado v0.7.x.
- `r1041_e82_400bps_hac_v4.3.0` — for HAC-basecalled R10.4.1 reads.
- `r941_min_sup_g507` — legacy: R9.4.1 + Guppy SUP + Guppy 5.0.7 (the most-cited deprecated model).
- `r941_min_hac_g507` — legacy: R9.4.1 + Guppy HAC.

The model name encodes the basecaller chemistry (`r1041` = R10.4.1), the model variant (`sup` = super-accuracy), and the version. Mismatching is a frequent foot-gun: applying `r941_min_sup_g507` to R10.4.1 reads produces a polish that is *worse than no polish*, because the model is using its training-data prior about which errors to correct on a dataset whose errors are differently distributed.

The defensive pattern: pin the basecaller model in the `run-info.json` for the read stage; pin the Medaka model in the `run-info.json` for the polish stage; in the polish-stage code, assert that the Medaka model matches the basecaller. Example:

```python
def assert_medaka_model_matches_basecaller(
    medaka_model: str,
    basecaller_model: str,
) -> None:
    """Soft check: warn if the chemistry prefix of the two models disagrees.

    Hard chemistry mismatches (r1041 vs r941) almost always produce a
    worse polish than skipping the polish.
    """
    medaka_chem: str = medaka_model.split("_")[0]
    basecaller_chem: str = (
        "r1041" if "r10.4" in basecaller_model
        else "r941" if "r9.4" in basecaller_model
        else "unknown"
    )
    if medaka_chem != basecaller_chem:
        raise ValueError(
            f"Medaka model {medaka_model} and basecaller {basecaller_model} "
            f"do not share a chemistry prefix; mismatch will degrade the polish."
        )
```

### Calling Medaka from Python

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def run_medaka(
    reads_fastq: Path,
    draft_fasta: Path,
    out_dir: Path,
    medaka_model: str = "r1041_e82_400bps_sup_v4.3.0",
    threads: int = 4,
) -> Path:
    """Run medaka_consensus. Returns the polished FASTA path.

    Raises subprocess.CalledProcessError if Medaka fails.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "medaka_consensus",
        "-i", str(reads_fastq),
        "-d", str(draft_fasta),
        "-o", str(out_dir),
        "-m", medaka_model,
        "-t", str(threads),
    ]
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return out_dir / "consensus.fasta"
```

The function returns the path to the polished FASTA. The Medaka log is at `out_dir/medaka.log` and includes the model name, the number of windows polished, and the per-window confidence statistics.

### When to skip Medaka

- **HiFi input.** HiFi assemblies are typically QV40+ from the assembler; Medaka was not trained on HiFi error profiles and will not help (and may hurt). Skip Medaka for HiFi.
- **Low coverage (< 20x).** Medaka relies on a confident read pileup at each position; below ~20x coverage the pileup is too thin and the polish can introduce errors. For very low coverage, accept the QV20 assembly and note the limit.
- **Wrong basecaller model.** As above; mismatched models produce worse polishes than no polish. If you cannot match the model, skip.

---

## 4. DeepConsensus in detail

Citation: Baid G, Cook DE, Shafin K, Yun T, Llinares-Lopez F, Berthet Q, Belyaeva A, Topfer A, Wenger AM, Rowell WJ, Yang H, et al. "DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer." *Nature Biotechnology* 41:232 (2023). Free at <https://www.nature.com/articles/s41587-022-01435-7>. Tool docs at <https://github.com/google/deepconsensus>.

DeepConsensus is the HiFi-specialized polisher, but unlike Medaka it operates **upstream of the assembly**: it polishes the per-read consensus that the standard PacBio CCS algorithm produces.

### What DeepConsensus does

The standard PacBio CCS algorithm builds a HiFi read from multiple subreads of the same circular template using a hidden-Markov-model-based consensus. The HMM is fast and produces QV30-35 reads on a typical Sequel II run. DeepConsensus replaces the HMM with a transformer-based model trained on a large HiFi training set: same input (subreads), more accurate output (HiFi reads at QV35-45).

The improvement at the read level translates to the assembly level: a Hifiasm assembly of DeepConsensus reads is typically QV45-50, versus QV40-45 for standard-CCS HiFi reads. The trade-off is compute: DeepConsensus is ~10-100x slower than the standard CCS algorithm and is typically run on a GPU.

### The canonical DeepConsensus call

DeepConsensus is invoked via the `deepconsensus run` subcommand. Inputs: the subreads BAM, the standard-CCS HiFi BAM (used as a draft). Output: a DeepConsensus-polished HiFi BAM.

```bash
deepconsensus run \
    --subreads_to_ccs subreads_to_ccs.bam \
    --ccs_bam hifi_reads.bam \
    --checkpoint deepconsensus_model/checkpoint \
    --output deepconsensus_hifi.fastq
```

The Week 10 exercises and mini-project do not invoke DeepConsensus (it requires GPU access most students do not have); the Hifiasm challenge runs on `badread`-simulated HiFi-style reads whose accuracy already matches DeepConsensus output. The call is documented here for production reference.

### When to use DeepConsensus

- **Publication-quality HiFi assemblies of large genomes.** Worth the ~1-2 GPU hours per sample.
- **Variant calling on HiFi reads.** A HiFi consensus that is QV45 instead of QV35 reduces the false-positive variant rate by ~10x.
- **Comparative analyses where every base matters.** Pangenome graphs, structural-variant calling, allele-specific expression.

For exploratory analyses, the standard CCS HiFi reads are usually fine.

---

## 5. Assembly contiguity statistics: N50, L50, total length

The conventional contiguity statistics for a contig set:

```
n_contigs       = number of contigs
total_length_bp = sum of contig lengths
longest_bp      = max(contig lengths)
N50             = length of the contig at which the cumulative length (descending) first reaches >= total_length_bp / 2
L50             = the rank (1-indexed) of that contig
```

Implementing N50 / L50 by hand is a five-line exercise. Every Week 10 deliverable computes them in Python rather than calling out to a tool like `seqkit stats` or `assembly-stats`:

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass
class AssemblyStats:
    n_contigs: int
    total_length_bp: int
    longest_contig_bp: int
    shortest_contig_bp: int
    n50_bp: int
    l50: int
    gc_fraction: float


def compute_assembly_stats(fasta_path: Path) -> AssemblyStats:
    """Compute N50, L50, total length, contig count, GC fraction.

    For empty assemblies, returns AssemblyStats with all-zero fields.
    """
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        return AssemblyStats(
            n_contigs=0,
            total_length_bp=0,
            longest_contig_bp=0,
            shortest_contig_bp=0,
            n50_bp=0,
            l50=0,
            gc_fraction=0.0,
        )

    lengths: list[int] = sorted((len(r.seq) for r in records), reverse=True)
    total: int = sum(lengths)

    cumulative: int = 0
    n50_bp: int = 0
    l50: int = 0
    for rank, length in enumerate(lengths, start=1):
        cumulative += length
        if cumulative >= total / 2:
            n50_bp = length
            l50 = rank
            break

    gc_bases: int = 0
    for r in records:
        seq_upper: str = str(r.seq).upper()
        gc_bases += seq_upper.count("G") + seq_upper.count("C")
    gc_fraction: float = gc_bases / total if total else 0.0

    return AssemblyStats(
        n_contigs=len(records),
        total_length_bp=total,
        longest_contig_bp=lengths[0],
        shortest_contig_bp=lengths[-1],
        n50_bp=n50_bp,
        l50=l50,
        gc_fraction=gc_fraction,
    )
```

### Interpreting N50 and L50

- **N50 high, L50 = 1.** Single-contig assembly. Best case for a small bacterial genome.
- **N50 high, L50 small (2-3).** A few large contigs. Typical for a clean bacterial-sized assembly with a small number of unresolved repeats.
- **N50 small, L50 large.** Many small contigs. Fragmented assembly; under-covered or wide-repeat-blocked.
- **N50 high, L50 = 1, but total_length_bp << expected_genome_size.** Suspicious: the longest contig may be an overcollapsed repeat region. Check the per-contig coverage; if it is N times the genome median, you have a collapse.

The mini-project pipeline reports all four numbers and flags any combination that warrants investigation.

### A worked example

Assume an assembly with contigs of lengths (sorted descending): 500k, 300k, 200k, 100k. Total = 1,100,000.

Cumulative pass:

- Rank 1 (500k): cumulative = 500,000. Less than 1,100,000 / 2 = 550,000.
- Rank 2 (300k): cumulative = 800,000. Greater than 550,000. Stop.

So N50 = 300,000 and L50 = 2.

Now assume an assembly with a single contig of length 1,000,000. Total = 1,000,000.

- Rank 1 (1M): cumulative = 1,000,000. Greater than 500,000. Stop.

N50 = 1,000,000 and L50 = 1. The maximum-contiguity case.

---

## 6. BUSCO: gene-content completeness

Citation: Manni M, Berkeley MR, Seppey M, Simao FA, Zdobnov EM. "BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes." *Molecular Biology and Evolution* 38:4647 (2021). Free at <https://academic.oup.com/mbe/article/38/10/4647/6329644>. Tool docs at <https://busco.ezlabgyk.unige.ch/>.

BUSCO (Benchmarking Universal Single-Copy Orthologs) scores an assembly against a curated set of genes that are expected to be present in single copy in every member of a clade. The output is a fraction:

- **C** (complete) — fraction of orthologs present as a complete-length gene.
- **S** (complete and single-copy) — present in exactly one copy.
- **D** (complete and duplicated) — present in more than one copy.
- **F** (fragmented) — present but split across multiple contigs or interrupted.
- **M** (missing) — not detected.

For a well-assembled genome, S should be high (~98-99%), D should be low (~0-1%; high D suggests the assembler did not collapse haplotypes that should have been collapsed), F should be low (~0-2%), and M should be low (~0-3%; high M suggests fragmentation).

### Picking a lineage

The lineage list is at <https://busco-data.ezlab.org/v5/data/lineages/>. The general rule: pick the most specific lineage that covers your organism, because more specific lineages have more single-copy orthologs and produce a higher-resolution score. Examples:

- **Bacteria** (most uses): `bacteria_odb10`, 124 orthologs.
- **E. coli** specifically: `enterobacterales_odb10`, ~440 orthologs.
- **Vertebrates**: `vertebrata_odb10`, 3,354 orthologs.
- **Mammals**: `mammalia_odb10`, 9,226 orthologs.
- **Fungi**: `fungi_odb10`, 758 orthologs.
- **Plants**: `viridiplantae_odb10`, 425 orthologs.

For the Week 10 simulated dataset (a synthetic 1 Mb FASTA with no real genes) BUSCO will report C ~= 0, because the synthetic reference has no orthologs to detect. For the *E. coli* mini-project case, `bacteria_odb10` is the right choice and a well-polished assembly should score C >= 98%.

### The canonical BUSCO call

```bash
busco \
    -i medaka_out/consensus.fasta \
    -l bacteria_odb10 \
    -o busco_run \
    --out_path busco_out \
    -m genome \
    -c 4 \
    --offline \
    --download_path busco_downloads
```

For a 1 Mb bacterial assembly this runs in ~30-60 seconds. The result of interest is `busco_out/busco_run/short_summary.specific.bacteria_odb10.busco_run.txt`:

```
# BUSCO version is: 5.7.1
# The lineage dataset is: bacteria_odb10 (Creation date: 2024-01-08, number of genomes: 4085, number of BUSCOs: 124)
        C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124
        123     Complete BUSCOs (C)
        123     Complete and single-copy BUSCOs (S)
        0       Complete and duplicated BUSCOs (D)
        0       Fragmented BUSCOs (F)
        1       Missing BUSCOs (M)
        124     Total BUSCO groups searched
```

Read the `C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124` line as: 99.2% complete (99.2% single-copy plus 0.0% duplicated), 0.0% fragmented, 0.8% missing, out of 124 orthologs total.

### Parsing the BUSCO summary in Python

```python
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class BuscoSummary:
    busco_version: str
    lineage_dataset: str
    n_total: int
    n_complete: int
    n_single_copy: int
    n_duplicated: int
    n_fragmented: int
    n_missing: int
    pct_complete: float
    pct_single_copy: float
    pct_duplicated: float
    pct_fragmented: float
    pct_missing: float


def parse_busco_summary(summary_path: Path) -> BuscoSummary:
    """Parse a BUSCO short_summary.specific.*.txt file."""
    text: str = summary_path.read_text()
    version_match = re.search(r"BUSCO version is:\s+(\S+)", text)
    lineage_match = re.search(r"lineage dataset is:\s+(\S+)", text)
    pct_match = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    n_complete_match = re.search(r"(\d+)\s+Complete BUSCOs \(C\)", text)
    n_single_match = re.search(r"(\d+)\s+Complete and single-copy BUSCOs \(S\)", text)
    n_dup_match = re.search(r"(\d+)\s+Complete and duplicated BUSCOs \(D\)", text)
    n_frag_match = re.search(r"(\d+)\s+Fragmented BUSCOs \(F\)", text)
    n_missing_match = re.search(r"(\d+)\s+Missing BUSCOs \(M\)", text)

    if not (
        version_match
        and lineage_match
        and pct_match
        and n_complete_match
        and n_single_match
        and n_dup_match
        and n_frag_match
        and n_missing_match
    ):
        raise ValueError(f"Could not parse BUSCO summary at {summary_path}")

    return BuscoSummary(
        busco_version=version_match.group(1),
        lineage_dataset=lineage_match.group(1),
        n_total=int(pct_match.group(6)),
        n_complete=int(n_complete_match.group(1)),
        n_single_copy=int(n_single_match.group(1)),
        n_duplicated=int(n_dup_match.group(1)),
        n_fragmented=int(n_frag_match.group(1)),
        n_missing=int(n_missing_match.group(1)),
        pct_complete=float(pct_match.group(1)),
        pct_single_copy=float(pct_match.group(2)),
        pct_duplicated=float(pct_match.group(3)),
        pct_fragmented=float(pct_match.group(4)),
        pct_missing=float(pct_match.group(5)),
    )
```

The parser is robust to BUSCO format minor changes (the percentages are extracted as a single regex; the counts as five separate regexes). If BUSCO's output format changes more substantially, the parser raises a clear error rather than silently producing wrong numbers.

### Interpreting BUSCO scores

- **C >= 98%, F < 2%, M < 2%.** Well-assembled genome. Publication-quality on the gene-content axis.
- **C = 80-98%, F >= 5%.** Fragmented assembly; many genes are split across contigs. The N50 is also likely low.
- **C >= 95%, D >= 5%.** Duplicated orthologs. Either the input was diploid and the assembler did not collapse haplotypes (Hifiasm in `--n-hap 2` mode does this by design; Flye should collapse), or the assembler emitted spurious near-duplicates.
- **C < 80%.** Either the assembly is severely incomplete, or the wrong lineage was picked. Check the lineage choice and re-run.

BUSCO is a complementary measure to N50: a high N50 with a low BUSCO is a structurally-contiguous assembly that is missing genes (suggests a wrong-species reference or a high error rate that prevents gene prediction); a low N50 with a high BUSCO is a fragmented but gene-complete assembly (the genes are intact but distributed across many small contigs). Both numbers are needed for a complete QC story.

---

## 7. Bandage: visual inspection of the assembly graph

Citation: Wick RR, Schultz MB, Zobel J, Holt KE. "Bandage: interactive visualization of de novo genome assemblies." *Bioinformatics* 31:3350 (2015). Free at <https://academic.oup.com/bioinformatics/article/31/20/3350/196114>. Tool docs at <https://rrwick.github.io/Bandage/>.

Bandage is a GUI for viewing GFA files. The CLI (`Bandage image input.gfa output.png`) renders the graph to a PNG without opening the GUI.

### Topologies to recognize

- **Linear path.** A single segment with no incoming or outgoing edges. The most common topology for a fragmented assembly contig.
- **Self-loop (circular).** A single segment with one self-loop edge. The healthy topology for a bacterial chromosome or a plasmid.
- **Bubble.** Two parallel segments between the same endpoints. A heterozygous-diploid region (kept in Hifiasm) or a sequencing-error artefact (kept in Flye `--keep-haplotypes`; collapsed otherwise).
- **Tangle.** A subgraph with many segments and many edges, no clear path through. An unresolved repeat region.
- **Tip / dead-end.** A segment with one connection at one end and no connection at the other. Often a low-confidence read fragment.

The mini-project pipeline (and Challenge 2) renders the GFA via:

```bash
Bandage image flye_out/assembly_graph.gfa graph.png --height 800
```

…and the write-up describes which topologies are present. A clean bacterial assembly should be one self-loop; a tangle is a flag.

### Bandage CLI from Python

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def render_bandage_image(
    gfa_path: Path,
    output_png: Path,
    height: int = 800,
) -> None:
    """Render a GFA graph to a PNG via the Bandage CLI.

    Falls back silently if Bandage is not on the PATH (the rendering is
    not strictly required; the assembly is judged on the contig FASTA).
    """
    output_png.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "Bandage", "image",
        str(gfa_path),
        str(output_png),
        "--height", str(height),
    ]
    try:
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        print("[bandage] Bandage CLI not on PATH; skipping graph render.")
    except subprocess.CalledProcessError as exc:
        print(f"[bandage] Bandage failed: {exc.stderr}")
```

The graceful-skip pattern is important: Bandage is a binary download from GitHub, not a pip-installable Python package, so we should not hard-fail the pipeline when it is absent. The same pattern applies to Medaka and BUSCO in the mini-project's optional-tool branches.

---

## 8. QV: Phred-scaled per-base error rate

When a reference genome is available (the simulated case in Week 10; the validated-reference case for many published bacterial genomes), the assembly QV can be computed exactly by aligning the assembly to the reference and tallying matches / mismatches / indels:

```
QV = -10 * log10(error_rate)
where error_rate = (mismatches + indels) / aligned_bases
```

The alignment is done with `minimap2`:

```bash
minimap2 -a -x asm5 reference.fasta assembly.fasta > assembly_to_ref.sam
```

(`asm5` is the minimap2 preset for assembly-to-reference alignment at 5% divergence; for higher-divergence comparisons use `asm10` or `asm20`.)

A Python helper to tally the SAM and compute QV:

```python
from __future__ import annotations

import math
import re
from pathlib import Path


def compute_assembly_qv(sam_path: Path) -> tuple[int, int, int, float]:
    """Compute (aligned_bp, mismatches, indels, QV) from a SAM file.

    Assumes the SAM contains one alignment per contig (the typical
    asm5 / asm10 output from minimap2). Skips supplementary and
    unmapped records.
    """
    aligned_bp: int = 0
    mismatches: int = 0
    indels: int = 0
    cigar_token = re.compile(r"(\d+)([MIDNSHP=X])")
    with sam_path.open() as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            parts: list[str] = line.split("\t")
            flag: int = int(parts[1])
            if flag & 0x4 or flag & 0x800:  # unmapped or supplementary
                continue
            cigar: str = parts[5]
            for length_str, op in cigar_token.findall(cigar):
                length: int = int(length_str)
                if op in ("M", "=", "X"):
                    aligned_bp += length
                if op == "X":
                    mismatches += length
                if op in ("I", "D"):
                    indels += length
            # The NM tag gives total edit distance; prefer it if present.
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    nm: int = int(tag.split(":")[2])
                    if "=" not in cigar and "X" not in cigar:
                        # CIGAR uses M (match-or-mismatch); NM gives mismatches + indels.
                        # Approximate split: assume mismatch:indel ratio 4:1 (typical).
                        mismatches = int(nm * 0.8)
                        indels = int(nm * 0.2)
                    break
    if aligned_bp == 0:
        return 0, 0, 0, 0.0
    error_rate: float = (mismatches + indels) / aligned_bp
    if error_rate <= 0:
        qv: float = 60.0  # cap; effectively zero error
    else:
        qv = -10.0 * math.log10(error_rate)
    return aligned_bp, mismatches, indels, qv
```

For the 1 Mb simulated / 50x case after one round of Medaka polish, the QV is typically 35-40 (about 1 error per 3,000-10,000 bp). Raw Flye output (before Medaka) is typically QV20-25 (about 1 error per 100-300 bp). The polish gain is a factor of ~5-10x in error rate, which is the headline number worth reporting.

For real-data cases without a reference, QV is estimated from k-mer methods (Merqury; Rhie et al. 2020, *Genome Biology* 21:245, free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7488247/>) using high-accuracy short reads as the ground truth. Merqury is out of scope for Week 10 but is the right tool for production assembly QC.

---

## 9. Writing the QC report

The Week 10 deliverable is a Markdown QC report rendered from the various statistics:

```markdown
# Assembly QC Report

Generated: 2025-05-14T13:42:11Z

## Input
- Platform: nanopore_simulated
- Reads: 3,289 reads after filter (length >= 1000, mean QV >= 10)
- Coverage: 49.4x

## Assembly
- Tool: Flye 2.9.5 (--nano-hq --genome-size 1m --threads 4 --iterations 1)
- Contigs: 1
- Total length: 1,038,203 bp
- Longest contig: 1,038,203 bp
- N50: 1,038,203 bp
- L50: 1
- GC fraction: 0.508
- Circular contigs: 1 (contig_1)

## Polish
- Tool: Medaka 1.12.0
- Model: r1041_e82_400bps_sup_v4.3.0
- Polished assembly QV (vs reference): 38.2

## BUSCO
- Version: 5.7.1
- Lineage: bacteria_odb10 (created 2024-01-08, 124 orthologs)
- Complete (C): 99.2% (single-copy 99.2%, duplicated 0.0%)
- Fragmented (F): 0.0%
- Missing (M): 0.8%

## Graph topology (Bandage)
- One self-loop on contig_1: the bacterial chromosome closed cleanly.
- No tangles, no bubbles, no tips.

## Limits and caveats
- Validated on simulated reads only; real-data validation requires a public ENA flow cell.
- Synthetic reference does not contain real genes; BUSCO score on the demo reference is artificially zero (the bacteria_odb10 result above is from the *E. coli* mini-project case).
- Thread count pinned at 4; varying the thread count can produce contigs differing on ~50 bp.

## Provenance
See `run-info.json` for full version pinning and the input MD5.
```

The report is a literal Markdown file; the pipeline writes it as part of the post-run summary. The `run-info.json` is the machine-readable counterpart.

---

## 10. Standard failure modes at the polish-and-QC stage

- **Wrong Medaka model.** As above; produces worse than no polish. Pin the model in the run-info JSON; assert the chemistry prefix matches the basecaller.
- **BUSCO with the wrong lineage.** A bacterial assembly scored against `vertebrata_odb10` will report C ~= 0% and look catastrophic. Pin the lineage; verify it makes sense for the organism.
- **BUSCO dataset version drift.** `bacteria_odb10` released in 2024-01-08 has 124 orthologs; an updated release may have 125 or 123. Always record the lineage creation date alongside the lineage name.
- **N50 inflation from collapsed repeats.** A 1 Mb assembly that should be 1.2 Mb but reports N50 = 1 Mb is overcollapsed; the apparent contiguity is misleading. Check the per-contig coverage against the genome median; collapsed contigs have coverage > 1.5x the median.
- **Wrong reference for QV computation.** Aligning the assembly to a closely-related but not-the-same reference produces a low QV (the differences are real biology, not assembly errors). Verify the reference identity (NCBI accession; MD5 of the FASTA) before reporting QV.
- **Polish-introduced indels at homopolymers.** Even with the correct Medaka model, a small fraction of homopolymer positions (~0.01%) are polished incorrectly (e.g. an A8 reference becomes an A7 in the polished contig). The error is platform-systematic and reduces with newer chemistries; pinning the chemistry version helps document the residual rate.

The mini-project write-up names which of these the assembly can plausibly suffer from. The N50-inflation-from-collapse and the homopolymer indels are the two ONT-specific failure modes that the write-up should always address.

---

## 11. Lecture 3 wrap-up

We described post-assembly polishing with Medaka (ONT) and upstream consensus polishing with DeepConsensus (HiFi); computed N50 and L50 by hand; ran BUSCO on a polished assembly and interpreted the summary line; named the Bandage topologies (linear, circular, bubble, tangle) and matched them to assembly quality outcomes; computed QV from a minimap2-aligned SAM; and wrote a QC report in Markdown.

The exercises pick up here. Exercise 1 simulates reads with badread and runs Flye. Exercise 2 computes N50 / L50 and parses a BUSCO summary. Exercise 3 polishes with Medaka and quantifies the QV improvement. The mini-project ties everything from Lectures 1-3 into a one-command end-to-end pipeline with a `run.sh` reproducer.

---

## Lecture 3 sanity check (10 minutes)

If you have an assembly FASTA from the Lecture 2 sanity check, compute its N50 by hand:

```bash
python -c "
from Bio import SeqIO
records = list(SeqIO.parse('/tmp/flye_out/assembly.fasta', 'fasta'))
lengths = sorted((len(r.seq) for r in records), reverse=True)
total = sum(lengths)
cumulative = 0
for rank, length in enumerate(lengths, start=1):
    cumulative += length
    if cumulative >= total / 2:
        print(f'n_contigs = {len(records)}')
        print(f'total = {total}')
        print(f'N50 = {length}')
        print(f'L50 = {rank}')
        break
"
```

Expected output on the Lecture 2 100 kb sanity check: `n_contigs = 1`, `total = ~100,000`, `N50 = ~100,000`, `L50 = 1`. If the assembly fragmented, you may see 2-3 contigs with a smaller N50 — that is also a valid sanity-check outcome.

---

## Further reading

- Wenger AM et al. "Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome." *Nat Biotechnol* 37:1155 (2019). The HiFi paper.
- Baid G et al. "DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer." *Nat Biotechnol* 41:232 (2023). The DeepConsensus paper.
- Manni M et al. "BUSCO update: novel and streamlined workflows ..." *Mol Biol Evol* 38:4647 (2021). The BUSCO paper.
- Wick RR et al. "Bandage: interactive visualization of de novo genome assemblies." *Bioinformatics* 31:3350 (2015). The Bandage paper.
- Rhie A, Walenz BP, Koren S, Phillippy AM. "Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies." *Genome Biology* 21:245 (2020). The reference-free QV estimator; production tool for QC when no reference is available.

Continue to the [exercises](../03-exercises/) for the hands-on work.
