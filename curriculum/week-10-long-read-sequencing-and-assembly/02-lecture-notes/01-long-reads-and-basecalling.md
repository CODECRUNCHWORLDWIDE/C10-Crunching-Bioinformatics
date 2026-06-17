# Lecture 1 — Long-Read Sequencing and Basecalling

> **Reproducibility note.** Long-read data is a moving target. The Oxford Nanopore Technologies (ONT) basecaller models change every few months; the PacBio HiFi consensus algorithm has a documented version history; and the same raw signal archive can produce different basecalls depending on which Dorado release you run. Pin the chemistry (`R10.4.1`), pin the basecaller model (`dna_r10.4.1_e8.2_400bps_sup@v4.3.0`), and pin the Dorado version (`0.7.2`) in any `run-info.json`. Skip any of those and the FASTQ you analyse becomes a moving target.

> **Duration:** ~3 hours of reading + a short hands-on Dorado / badread sanity check.
> **Outcome:** You can describe the physical principle of nanopore sequencing (current trace through a pore, neural network basecall) and PacBio HiFi (circular consensus over multiple polymerase passes), state the read-length and per-base error-rate ranges of each, pick the correct Flye input mode (`--nano-hq` vs `--pacbio-hifi`), and call `badread simulate` from Python with `subprocess.run(..., check=True)` to produce a synthetic FASTQ.

If you only remember one thing from this lecture, remember this:

> **A long read is a long string of bases with a per-base error rate that depends on the platform (ONT or HiFi) and the chemistry / model version. The error mode is platform-specific: ONT errors are systematic at homopolymers; HiFi errors are randomly distributed at the per-read level but the consensus over multiple passes removes most of them. Pin the chemistry, pin the basecaller model, and the downstream assembly becomes reproducible. Skip either and you are doing exploratory work, not engineering work.**

Week 9's MSA-and-tree pipeline is parked. Week 10's input is a FASTQ of long reads from a sequencer. Week 10's intermediate artefact is a set of assembled contigs. The QC and the polish are downstream.

---

## 1. Where we are in the pipeline

The end-to-end Week 10 pipeline runs:

```
Raw signal (POD5 from ONT; subreads BAM from PacBio) ->
        basecaller (Dorado for ONT; CCS / DeepConsensus for HiFi) ->
        FASTQ (N long reads, mean ~15 kb, per-base QV variable) ->
        read filter (drop reads < 1,000 bp; optionally drop reads < QV10) ->
        FASTQ (clean) ->
        assembler (Flye for ONT; Hifiasm for HiFi; Canu as the alternative) ->
        FASTA (M contigs, total length ~ genome size) ->
        polisher (Medaka for ONT; usually nothing for HiFi) ->
        FASTA (polished) ->
        QC (asmstats / N50, BUSCO, Bandage on the GFA).
```

Lecture 1 covers the first three boxes — raw signal, basecaller, FASTQ — and stops at the cleaned reads. Lectures 2 and 3 take the cleaned reads through assembly, polishing, and QC.

The single most consequential step on the read side is the basecaller. A poor basecall produces poor reads produces a poor assembly. The assembler cannot recover from miscalled bases, particularly miscalls at homopolymers (a run of A's in the reference may be miscalled as a slightly shorter or longer A-run; if every read miscalls the homopolymer, the consensus is wrong by one base regardless of read depth). The polishing stage helps, but the basecaller sets the floor on per-base accuracy.

---

## 2. The two long-read platforms in one paragraph each

### Oxford Nanopore Technologies (ONT)

ONT drives a single strand of DNA through a protein pore embedded in a membrane and measures the ionic current as each base passes. The current trace is called the **squiggle** and a neural network (the basecaller) translates the squiggle into a base sequence. The pores live on a flow cell (a hardware unit; the current generation is the `R10.4.1` chemistry on a `MinION`, `GridION`, or `PromethION` instrument). A single flow cell produces ~10-100 Gbp of basecalled sequence over a 24-72 hour run. The read length distribution is wide: the median is ~15-25 kb, the N50 is ~25-40 kb, and the longest reads can exceed 1 Mb on libraries prepared for "ultra-long" mode. The per-base error rate after Dorado SUP basecalling is ~1-2% on most regions and slightly higher in homopolymer runs and methylated regions.

### PacBio Single Molecule Real Time (SMRT) sequencing

PacBio attaches a polymerase to the bottom of a zeptolitre well (a "Zero Mode Waveguide", ZMW) and watches the polymerase incorporate fluorescently-labelled nucleotides. The raw reads (continuous long reads, CLR) are noisy (~10-15% per-base error rate) but the same circular template is read multiple times by the same polymerase. The **HiFi** protocol (the dominant mode since 2019; Wenger et al. 2019, *Nat Biotechnol* 37:1155) takes the multiple passes and computes a per-base consensus across them; the resulting **HiFi read** is ~15-20 kb long and ~99.8-99.9% accurate. The instrument is the `Sequel II` (older) or the `Revio` (current; 2022-onwards). A single Revio SMRT Cell produces ~30 Gbp of HiFi reads in ~24 hours.

### The trade-off

| Property | ONT R10.4.1 (Dorado SUP) | PacBio HiFi (Revio) |
|----------|--------------------------|---------------------|
| Mean read length | ~20 kb | ~15-20 kb |
| Read length N50 | ~30 kb | ~18 kb |
| Longest reads | > 1 Mb (ultra-long) | ~30 kb |
| Per-base error rate | ~1-2% (SUP); ~5% (HAC) | ~0.1-0.2% |
| Error mode | systematic at homopolymers and methylated bases | mostly random; rare CCS-pass-dependent errors |
| Cost per Gbp (2024) | ~$10-30 | ~$30-50 |
| Run time | 24-72 hours | ~24 hours |
| Instrument footprint | desktop (MinION) to chassis (PromethION) | chassis only (Revio) |
| Open-source toolchain | Dorado (open) | CCS, DeepConsensus, Hifiasm (open) |
| Best assembler | Flye (`--nano-hq`); Canu | Hifiasm; Flye (`--pacbio-hifi`) |
| Best polisher (post-assembly) | Medaka | usually none (polish is upstream at the CCS step) |

The summary: ONT optimizes for read length and run-time flexibility at the cost of per-base error; HiFi optimizes for per-base accuracy at the cost of read length and per-read cost. Both produce reads long enough to span most repeats in a typical bacterial or small-eukaryotic genome.

---

## 3. The nanopore squiggle in detail

Citation context: Jain M, Olsen HE, Paten B, Akeson M. "The Oxford Nanopore MinION: delivery of nanopore sequencing to the genomics community." *Genome Biology* 17:239 (2016). Free at <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1103-0>.

When a DNA strand translocates through the R10.4.1 nanopore, the ionic current through the pore changes depending on the identity of the bases inside the pore at that instant. The pore has a "reader head" that is approximately 5 bases wide; the current at any instant is determined by all 5 bases inside the reader head, not just the central base. So the basecaller is not a per-base classifier — it is a sequence-to-sequence model that maps the time series of currents (the squiggle) to a base string.

The squiggle is recorded as a time series at ~4 kHz sampling rate: ~4,000 current measurements per second of translocation. A 15 kb read at ~400 bases per second of translocation produces ~150,000 current measurements. The raw file format is **POD5** (an Arrow-based columnar format introduced in 2022) or the older **FAST5** (an HDF5-based format from 2014-2021). POD5 is ~3-5x smaller than FAST5 for the same data and reads ~10x faster on modern hardware; new ONT pipelines should default to POD5.

The basecaller (Dorado, the current ONT default; replacing the older Guppy from 2018-2022) is a recurrent neural network in the older models (LSTM-based) and a transformer in the newer models. The model is selected at the command line. The current default for R10.4.1 with super-accuracy mode is `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`, parsed as:

- `dna` — DNA mode (RNA models also exist).
- `r10.4.1` — the chemistry version.
- `e8.2` — the protein-engineering revision of the pore.
- `400bps` — the translocation speed (400 bases per second; older chemistries used 70 bps and 260 bps).
- `sup` — super-accuracy variant; the slowest and most accurate. Other variants are `hac` (high-accuracy; ~3x faster than `sup`) and `fast` (~10x faster than `sup`).
- `v4.3.0` — the model version.

Every run-info JSON for a nanopore experiment must record the full model string. Two runs with `dna_r10.4.1_e8.2_400bps_sup@v4.2.0` and `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` on the same POD5 archive will produce slightly different FASTQ — a few percent of bases will differ — and the downstream assembly will differ correspondingly. Pinning the model version is non-negotiable in production.

### Dorado from the command line

```bash
dorado basecaller \
    dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    pod5/ \
    --emit-fastq \
    > reads.fastq
```

Dorado requires a GPU for any non-toy run. On a consumer GPU (NVIDIA RTX 3060 or similar) it basecalls ~5-10x faster than the older Guppy. CPU-only basecalling is supported but ~100x slower; do not basecall a real flow cell on CPU.

Week 10's exercises and mini-project bypass Dorado entirely: we use `badread simulate` to produce a synthetic FASTQ from a known reference, with a tunable error profile that mimics R10.4.1 SUP. The Dorado call is documented here for completeness and for anyone with access to a real sequencer.

---

## 4. The PacBio HiFi protocol in detail

Citation: Wenger AM, Peluso P, Rowell WJ, Chang P-C, Hall RJ, Concepcion GT, Ebler J, Fungtammasan A, Kolesnikov A, Olson ND, et al. "Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome." *Nature Biotechnology* 37:1155 (2019). Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6776680/>.

The PacBio HiFi protocol works like this:

1. **Library prep.** The DNA is sheared to a target insert size of ~15-20 kb. Each insert is ligated to **SMRTbell adapters**, which form a hairpin at each end. The resulting molecule is a **circular template**: a single-stranded loop that contains the insert in both directions plus the two adapters.
2. **Polymerase loading.** A polymerase is loaded into each ZMW well and binds to one of the adapters.
3. **SMRT sequencing.** The polymerase circles the template, incorporating one fluorescently-labelled nucleotide at a time. Each incorporation event produces a flash of fluorescence that the instrument records. The polymerase makes ~5-15 passes around a 15 kb template before falling off or the run ending.
4. **Subread output.** Each pass produces a **subread**: a noisy long read of the insert in one direction, separated from the next subread by a stretch of adapter sequence. The subreads are emitted as an unaligned BAM.
5. **CCS / HiFi consensus.** The PacBio software (the `ccs` tool, open source; or the more recent `DeepConsensus`, Baid et al. 2023) takes the subreads from one ZMW and computes a consensus across them. With 5+ passes the consensus error rate drops below 1%; with 10+ passes it drops below 0.1%. The output is the **HiFi read**: one record per ZMW, ~15-20 kb long, ~99.8-99.9% accurate, with per-base quality scores that are calibrated against the consensus depth.

The HiFi read inherits the random-error advantage of consensus: any error that is not shared across multiple passes is corrected. The remaining errors are concentrated at the positions where the CCS algorithm itself is uncertain (very high or very low GC content; certain modified bases).

The output FASTQ is conventionally named `*.hifi_reads.fastq.gz`; the quality scores are Phred+33 with per-base values typically in Q30-Q50.

### CCS from the command line

```bash
ccs subreads.bam hifi_reads.bam \
    --min-passes 3 \
    --min-rq 0.99 \
    --threads 16
```

This produces an unaligned BAM of HiFi reads. To convert to FASTQ:

```bash
samtools fastq hifi_reads.bam > hifi_reads.fastq
```

DeepConsensus (Baid et al. 2023) is a drop-in replacement for the standard `ccs` consensus that uses a transformer-based model trained on a large HiFi training set. On the same subreads it typically produces a HiFi consensus that is ~0.05-0.1 log-error better (a QV40 read becomes QV43-45). The added compute is significant — ~1-2 GPU-hours per Gbp of subreads — so DeepConsensus is run for the final, high-quality analyses rather than the exploratory runs.

### Why HiFi assemblies often need no post-assembly polish

Because the consensus is computed *upstream* at the subread step, the HiFi reads themselves are already accurate enough that the assembly contigs inherit the accuracy. A Hifiasm assembly from HiFi input is typically QV40+; running Medaka or a similar post-assembly polisher on it would either not improve the contigs or actively degrade them (the polishers are trained on different read types). The conventional HiFi pipeline is `subreads -> CCS (or DeepConsensus) -> Hifiasm -> done`; the polishing step is absent.

---

## 5. Read quality, read length, and the assembly graph

The two read properties that dominate assembly outcomes are **per-base error rate** and **read length**. They affect different stages of the OLC pipeline (Lecture 2 unpacks the stages):

- **Per-base error rate** affects the overlap stage. Two reads from the same region must agree on a long enough stretch of bases to be detected as overlapping. With 1% per-base error (HiFi), a 5,000-bp overlap is easy. With 5% per-base error (R10.4.1 HAC), the same overlap is harder; the assembler typically requires longer overlaps (10-20 kb minimum) to achieve the same false-overlap rate.
- **Read length** affects the layout stage. To span a repeat of length L, you need reads of length > L. A 1 kb repeat is spanned by typical 15 kb reads; a 30 kb repeat (e.g. the human rDNA cluster) is rarely spanned. Repeats wider than the read length cannot be resolved; the assembler either collapses them or breaks the assembly there.

The two properties have a partial trade-off in real instruments: ONT pushes for length at the cost of accuracy; HiFi pushes for accuracy at the cost of length. For most bacterial genomes (4-7 Mb, repeats up to ~3 kb), both platforms produce reads long enough to span all repeats; the choice between them is set by per-base accuracy and per-Gbp cost.

For larger genomes with longer repeats (large eukaryotes, mammals with their long retrotransposon families), neither platform alone is sufficient and the conventional workflow is **ONT ultra-long + HiFi**: HiFi for per-base accuracy across most of the genome, ONT ultra-long for spanning the largest repeats. The combined assembler is Verkko (Rautiainen et al. 2023, *Nat Biotechnol* 41:1474, free at <https://www.nature.com/articles/s41587-023-01662-6>); it is out of scope for Week 10 but worth knowing exists.

---

## 6. Read simulation with badread

Week 10 cannot run on real sequencer output: a real run is 10-100 GB of FASTQ and the basecalling step requires GPU access most students do not have. The workaround is to **simulate** reads from a known reference using `badread` (Wick 2019, *Journal of Open Source Software* 4:1316; <https://github.com/rrwick/Badread>).

`badread simulate` takes a reference FASTA and emits a FASTQ that mimics a real long-read run in three properties:

1. **Length distribution.** Long-tailed; the user provides a mean and standard deviation; the empirical fit is approximately a gamma distribution.
2. **Per-read identity distribution.** Each read has its own per-base error rate sampled from a Beta-like distribution; some reads are near-perfect and others are very noisy.
3. **Error model.** Within a read, errors are placed with a position-dependent probability that mimics real ONT or PacBio data (homopolymer-concentrated for ONT; mostly random for PacBio).

The canonical Week 10 invocation for an ONT-like read set:

```bash
badread simulate \
    --reference data/reference_1mb.fasta \
    --quantity 50x \
    --length 15000,13000 \
    --identity 95,3,99 \
    --error_model nanopore2023 \
    --qscore_model nanopore2023 \
    --seed 42 \
    > data/reads_nanopore.fastq
```

The `--seed 42` makes the run deterministic: two invocations with the same seed and the same reference produce byte-identical FASTQ.

The trade-offs of working with simulated data:

- **Pro.** The "true" sequence is known (it is the reference). Assembly QC against the reference is trivial; the Phred-scaled error rate (QV) can be computed exactly.
- **Pro.** The dataset is small (a 1 Mb reference at 50x is ~50 MB FASTQ). It runs through Flye + Medaka in under 10 minutes on a laptop.
- **Con.** The error model is a model. Real R10.4.1 SUP reads have correlated errors (a tricky region miscalls in many reads at once) that the simulator's i.i.d. error model does not capture. Assemblies that look easy on simulated data can be harder on real data.

We use simulated data for the exercises and the mini-project because the alternative (downloading and basecalling a real flow cell) is out of scope for a week-long unit. The mini-project write-up names the limit explicitly: "this pipeline was validated on simulated reads only; real-data validation requires running on a publicly available ENA flow cell."

### Calling badread from Python

The pattern is the same as the Week 5 / Week 8 / Week 9 `subprocess.run` template:

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def simulate_nanopore_reads(
    reference: Path,
    output_fastq: Path,
    coverage: str = "50x",
    length: str = "15000,13000",
    identity: str = "95,3,99",
    seed: int = 42,
) -> None:
    """Run badread simulate. Writes the FASTQ to output_fastq."""
    output_fastq.parent.mkdir(parents=True, exist_ok=True)
    cmd: list[str] = [
        "badread", "simulate",
        "--reference", str(reference),
        "--quantity", coverage,
        "--length", length,
        "--identity", identity,
        "--error_model", "nanopore2023",
        "--qscore_model", "nanopore2023",
        "--seed", str(seed),
    ]
    with output_fastq.open("w") as fh:
        result = subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
    # badread writes progress to stderr; surface it on demand.
    if result.stderr:
        for line in result.stderr.splitlines():
            if line.strip().startswith("Generating"):
                print(f"[badread] {line.strip()}")
```

The seed argument is the load-bearing parameter for reproducibility. Without it the second run on the same reference produces a different FASTQ and the downstream assembly will differ.

---

## 7. Read filtering before assembly

Modern Flye and Hifiasm both apply internal read filters, but a pre-assembly filter pass is good hygiene because it makes the run-time predictable and the run-info JSON honest. The Week 10 default filter:

- Drop reads shorter than **1,000 bp**. Very short reads waste compute in the overlap stage and rarely contribute useful overlaps; the typical long-read assembler benefits from a length floor.
- Drop reads with mean Phred quality below **Q10** for ONT data; not used for HiFi (where every read is Q30+).
- Drop reads with adapter contamination on either end (a separate step in production; out of scope for Week 10).

The filter can be implemented in 30 lines of Biopython:

```python
from __future__ import annotations

from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def filter_long_reads(
    input_fastq: Path,
    output_fastq: Path,
    min_length: int = 1000,
    min_mean_qv: float = 10.0,
) -> tuple[int, int]:
    """Filter reads by length and mean Phred quality.

    Returns (n_in, n_out).
    """
    output_fastq.parent.mkdir(parents=True, exist_ok=True)
    n_in: int = 0
    n_out: int = 0
    with output_fastq.open("w") as fh:
        for record in SeqIO.parse(str(input_fastq), "fastq"):
            n_in += 1
            if len(record.seq) < min_length:
                continue
            quals: list[int] = record.letter_annotations["phred_quality"]
            if quals and (sum(quals) / len(quals)) < min_mean_qv:
                continue
            SeqIO.write([record], fh, "fastq")
            n_out += 1
    return n_in, n_out
```

For the 1 Mb / 50x simulated dataset, the filter typically drops ~5% of reads (the short tail of the length distribution) and reduces the FASTQ size from ~50 MB to ~45 MB. The assembly quality is unchanged within a few hundred bp of N50; the run time is reduced by a similar percentage.

**Stronger filters are not always better.** Dropping aggressively (e.g. requiring length > 10 kb and mean QV > 15) can reduce coverage below the assembly-viable threshold (~25-30x for nanopore Flye; ~15-20x for HiFi Hifiasm) and produce a worse assembly than the unfiltered run. The mini-project README pins the filter to "min length 1,000, min QV 10" as a compromise.

---

## 8. The run-info JSON for the read stage

Every Week 10 deliverable emits a `run-info.json` recording the parameters of every stage. The read-stage fields are:

```json
{
    "run_date": "2025-05-14T13:42:11Z",
    "platform": "nanopore_simulated",
    "reference_fasta": "data/reference_1mb.fasta",
    "reference_md5": "a4b5c6...",
    "reference_length_bp": 1040000,
    "basecaller": "badread 0.4.1",
    "basecaller_model": "nanopore2023",
    "read_simulation_seed": 42,
    "read_simulation_coverage": "50x",
    "read_length_mean_sd": "15000,13000",
    "read_identity_mean_sd_max": "95,3,99",
    "read_filter_min_length_bp": 1000,
    "read_filter_min_mean_qv": 10.0,
    "n_reads_simulated": 3457,
    "n_reads_after_filter": 3289,
    "fastq_bytes": 47892013
}
```

For a real-data pipeline the same fields exist with different values:

```json
{
    "platform": "nanopore_R10.4.1",
    "basecaller": "Dorado 0.7.2",
    "basecaller_model": "dna_r10.4.1_e8.2_400bps_sup@v4.3.0",
    "instrument": "PromethION P24",
    "flow_cell_id": "PAU12345",
    ...
}
```

Pin the fields. Two assemblies built from FASTQs without matching basecaller models or matching simulation seeds are not directly comparable; making them comparable requires re-basecalling or re-simulating from the same parameters, and the run-info is the only way to know what those parameters were.

---

## 9. Standard failure modes at the read stage

- **Basecaller drift.** A POD5 archive basecalled with Dorado v0.5.0 and the same archive basecalled with Dorado v0.7.0 produce FASTQs that differ on ~1-3% of bases; the downstream assembly differs correspondingly. Always pin the basecaller version in the run-info.
- **Wrong basecaller model.** Using the R9.4-trained `dna_r9.4.1_450bps_sup` model on R10.4.1 data is the canonical foot-gun: the basecaller produces a plausible-looking FASTQ that is systematically wrong (the model was not trained on the new chemistry's signal characteristics). The FASTQ headers do not encode the chemistry; you have to track it externally.
- **Adapter contamination.** A small fraction of reads contains the SMRTbell or ONT adapter sequence at one end or in the middle; the assembler can mistake the adapter for a real overlap and produce chimeric contigs. The Week 10 didactic dataset is adapter-free (`badread` does not insert adapters by default); real data needs an adapter-trimming pass with `Trimmomatic` or `Porechop`.
- **Chimeric reads.** Both ONT (multiple molecules through the same pore in succession) and HiFi (multiple inserts on the same SMRTbell adapter) can produce chimeric reads: a single read that spans two distinct regions of the source DNA. Modern HiFi has < 0.1% chimera rate; ONT has 1-3% depending on the library prep. The assembler's overlap stage usually rejects chimeras, but not always.
- **Quality-score miscalibration.** Old ONT (Guppy-basecalled, pre-2022) had quality scores that were not honestly Phred-scaled: a "Q10" base was empirically much worse than 10% error rate. Filtering by quality on miscalibrated scores can drop the wrong reads or keep too many. Dorado SUP from 2023 onwards has well-calibrated quality scores.

The mini-project write-up names which of these the pipeline can plausibly suffer from on its input. For the simulated-data case the answer is short: only basecaller-drift is irrelevant (there is no real basecaller); the others have analogues in the simulator's parameters.

---

## 10. Lecture 1 wrap-up

We named two long-read platforms (ONT R10.4.1; PacBio HiFi on Sequel II / Revio), one basecaller per platform (Dorado for ONT; CCS / DeepConsensus for HiFi), the trade-offs between them (length vs accuracy vs cost), the simulator we will actually use in the exercises (badread), the read filter we apply before assembly (min length 1,000; min Phred 10 for ONT), and the run-info JSON fields the FASTQ stage must record.

Lecture 2 picks up the cleaned reads and runs them through the OLC assembly graph in Flye, Canu, and Hifiasm. The lecture is mechanism-focused: how each tool builds the graph, how each tool handles repeats, and why the choice of tool depends on the input chemistry. Lecture 3 closes the loop with polishing (Medaka, DeepConsensus) and QC (BUSCO, Bandage, N50 / L50).

Before moving on, do the Lecture 1 sanity check below.

---

## Lecture 1 sanity check (15 minutes)

If you have `badread` and `biopython` installed on the local PATH:

```bash
# 1. Download or generate a small reference (10 kb of random nucleotides is fine for a sanity check).
python -c "
import random
random.seed(1)
seq = ''.join(random.choice('ACGT') for _ in range(10000))
print('>chr1')
print(seq)
" > /tmp/ref.fasta

# 2. Simulate 30x of badread reads.
badread simulate --reference /tmp/ref.fasta --quantity 30x \
    --length 5000,3000 --identity 95,3,99 \
    --error_model nanopore2023 --qscore_model nanopore2023 \
    --seed 42 > /tmp/reads.fastq

# 3. Count the reads and look at the first record.
python -c "
from Bio import SeqIO
recs = list(SeqIO.parse('/tmp/reads.fastq', 'fastq'))
print(f'n_reads = {len(recs)}')
print(f'first_id = {recs[0].id}')
print(f'first_len_bp = {len(recs[0].seq)}')
print(f'mean_len_bp = {sum(len(r.seq) for r in recs) / len(recs):.0f}')
"
```

Expected output: ~60 reads (30x coverage on 10 kb with mean read length 5 kb), the first read's ID starting with `read_`, the first read between 1,000 and 15,000 bp, and the mean read length close to 5,000 bp. The exact numbers depend on the badread version; the order of magnitude is robust.

If the sanity check passes, you are ready for Lecture 2.

---

## Further reading

- Jain M, Olsen HE, Paten B, Akeson M. "The Oxford Nanopore MinION: delivery of nanopore sequencing to the genomics community." *Genome Biology* 17:239 (2016). The classic introduction to nanopore sequencing.
- Wenger AM et al. "Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome." *Nature Biotechnology* 37:1155 (2019). The HiFi paper.
- Rautiainen M, Nurk S, Walenz BP, Logsdon GA, Porubsky D, Rhie A, Eichler EE, Phillippy AM, Koren S. "Telomere-to-telomere assembly of diploid chromosomes with Verkko." *Nat Biotechnol* 41:1474 (2023). The combined ONT-ultra-long + HiFi assembler; out of scope for Week 10 but the state of the art for large eukaryotic genomes.
- Nurk S et al. "The complete sequence of a human genome." *Science* 376:44 (2022). The T2T-CHM13 paper; the first end-to-end gapless human assembly, built from ONT + HiFi.
- Wick RR. Badread documentation. <https://github.com/rrwick/Badread>. The simulator we use in the exercises.

Continue to [Lecture 2 — De Novo Assembly and the OLC Graph](./02-de-novo-assembly-olc.md).
