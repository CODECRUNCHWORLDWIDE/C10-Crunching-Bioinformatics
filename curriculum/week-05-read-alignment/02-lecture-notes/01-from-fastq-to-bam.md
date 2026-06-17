# Lecture 1 — From FASTQ to BAM

> **Duration:** ~3 hours of reading + paper-and-pencil + a brief Python sanity check.
> **Outcome:** You can describe the FM-index in two paragraphs, name the four parameters that control BWA-MEM's seeding (`-k`, `-r`, `-c`, `-T`), choose between `bwa mem` and `minimap2` for a given read profile, and run the FASTQ → SAM → sorted BAM pipeline end to end on a small reference without copying commands from a tutorial.

If you only remember one thing from this lecture, remember this:

> **Modern short-read alignment is a two-phase algorithm: an FM-index over the reference lets you locate every exact occurrence of a short query substring in `O(m)` time independent of the reference's length `n`; a banded dynamic-programming extension then fills in the gaps around each seed. BWA-MEM (Li 2013) and minimap2 (Li 2018) are the production implementations. They are both fast because the seed lookup is sub-linear, not because the alignment is heuristic — the final alignment is still Smith-Waterman-style on the windowed extension, just on a tiny fraction of the reference.**

Week 4's BLAST was seed-and-extend at a 30,000-entry curated database scale. This week's BWA-MEM is seed-and-extend at a 1-entry-but-3-gigabase-long database scale. The seed step is the same idea in both — find short exact (or near-exact) matches first, then run optimal alignment on what survives. The data structure underneath is what differs: BLAST uses a hash table of word positions in the database; BWA uses an FM-index over the BWT of the reference. The FM-index buys you O(m) per-read search time instead of O(m · n_total), and that is the difference between "weeks per genome" and "minutes per genome."

---

## 1. The problem `bwa mem` and `minimap2` solve

A typical Illumina NovaSeq run produces 1–10 billion paired-end reads of 150 bp each. The reference is, in the human case, a 3 Gb (`3 × 10^9` bp) genome assembly. The naive approach — Smith-Waterman of every read against every position in the reference — costs:

```
T_SW ≈ O(read_count × m × n)
     ≈ 10^9 × 150 × 3·10^9
     ≈ 4.5 × 10^20 cell updates.
```

At ~10^9 cell updates per second on a single optimistic SIMD-vectorized core (Farrar's 2007 SW), one full alignment run takes about 14,000 years. The order-of-magnitude conclusion is the same as Week 4's: brute-force Smith-Waterman is not just slow, it is *centuries* slow for any genome larger than a phage.

The trick that makes alignment feasible is the same trick BLAST used: **filter first, align second**. For each read, locate the small number of reference positions where a meaningful alignment is *possible* (the seed step), then run Smith-Waterman-style extension only on those positions (the extension step). The seed step is sub-linear in the reference size when the reference is indexed correctly. The extension step is linear in the read length. The overall cost is dominated by the seed step, which is `O(m)` per read with an FM-index, giving a per-read alignment cost of `O(m)` *independent of reference length*. For 10^9 reads at 150 bp each, total cost is `10^9 × 150 = 1.5 × 10^11` operations — about 2.5 minutes at 10^9 ops/sec on a single core, or ~10 seconds on 16 cores. This is the speed that made whole-genome sequencing a routine clinical workflow.

---

## 2. The FM-index in three concepts

The FM-index (Ferragina & Manzini, *FOCS 2000*) is a compressed self-index built on top of the Burrows-Wheeler Transform of the reference. You do not need to implement it to use BWA — but a working mental model of its three core operations is the difference between "BWA is magic" and "I can predict what BWA will and will not do well."

### 2.1 The Burrows-Wheeler Transform

The BWT of a string `T` is a reversible permutation of `T` produced by sorting all rotations of `T$` (where `$` is a sentinel less than every other character) lexicographically and reading off the last column of the sorted matrix.

Example. Let `T = ACGTACGT$`. The nine rotations of `T` (treating `$` as the wraparound) are:

```
ACGTACGT$
CGTACGT$A
GTACGT$AC
TACGT$ACG
ACGT$ACGT
CGT$ACGTA
GT$ACGTAC
T$ACGTACG
$ACGTACGT
```

Sorted lexicographically (with `$ < A < C < G < T`):

```
$ACGTACGT     last column: T
ACGT$ACGT     last column: T
ACGTACGT$     last column: $
CGT$ACGTA     last column: A
CGTACGT$A     last column: A
GT$ACGTAC     last column: C
GTACGT$AC     last column: C
T$ACGTACG     last column: G
TACGT$ACG     last column: G
```

The BWT of `T` is the last column, read top to bottom: `TT$AACCGG`. This is a permutation of the original `T` — same characters, different order. The key empirical fact is that the BWT tends to group similar characters together (every `T` is followed by either a `$` or an `A`, in this toy case both `T`s landed adjacent to each other), which makes the BWT highly compressible by run-length encoding. For real genomes (3 Gb of `A/C/G/T/N`), the BWT compresses to ~1 GB; the FM-index of the human genome is roughly 3 GB on disk including auxiliary structures.

### 2.2 Rank and select

The FM-index supports two queries in O(1) time after `O(n)` preprocessing:

- **`rank(c, i)`** — the number of occurrences of character `c` in the BWT up to position `i`. For our example BWT `TT$AACCGG`, `rank('A', 5)` (count of `A` in positions 0..4 inclusive) = 2.
- **`select(c, k)`** — the position in the BWT of the `k`-th occurrence of `c`. `select('A', 1)` (position of the 1st `A`) = 3.

Both operations are O(1) with the right precomputed `Occ[c, i]` arrays (and a few bit-vector tricks for compression). The details are in the Ferragina-Manzini 2000 paper, sections 3–4.

### 2.3 Backward search

Given a query pattern `P` of length `m`, the FM-index answers "how many times does `P` occur in `T`?" and "where exactly?" in `O(m)` time. The algorithm is **backward search**: process `P` from right to left, maintaining an interval `[lo, hi)` of the suffix array that contains all suffixes prefixed by the suffix of `P` processed so far.

Pseudocode:

```
function backward_search(P, C, Occ):
    m = len(P)
    c = P[m-1]                       # rightmost char
    lo = C[c]                        # first row in sorted matrix starting with c
    hi = C[c+1]                      # one past last row starting with c
    for i = m-2 down to 0:
        c = P[i]
        lo = C[c] + Occ[c, lo]
        hi = C[c] + Occ[c, hi]
        if lo >= hi: return []       # P does not occur
    return suffix_array[lo .. hi]    # all occurrence positions
```

Two arrays at play: `C[c]` is the number of characters in `T` lexicographically smaller than `c` (a length-σ array, where σ is the alphabet size — 5 for DNA with N). `Occ[c, i]` is the rank-1 above. Each iteration of the loop is two `Occ` lookups and two additions. The total cost is `O(m)` — *independent of the reference length n*. That is the magic.

A read of length 150 bp, backward-searched against a 3 Gb human reference, takes ~150 array accesses. On modern hardware that is ~1 microsecond per exact-match query. A billion such queries — once per read — is ~1000 seconds, or 17 minutes single-core. The whole alignment pipeline (which adds chain-and-extend on top) is ~4x slower than this, so a billion 150 bp reads cost ~1 hour single-core or ~5 minutes on 16 cores.

---

## 3. BWA-MEM in five steps

BWA-MEM (Li 2013, `arXiv:1303.3997`) is the BWA mode used for reads ≥ 70 bp — which is essentially all modern Illumina output. The algorithm:

1. **Seed.** For each position in the read, find the **Maximal Exact Match** against the reference using the FM-index — the longest substring starting at this position that occurs at least once in the reference, such that extending by one more base would mismatch. MEMs are found via backward search variants in `O(m + occurrences)`.
2. **Reseed.** If a MEM is longer than a threshold (`-r` flag, default 1.5×k), split it into shorter overlapping seeds. Long MEMs that come from highly repetitive regions can produce thousands of occurrences and dominate the chain step; reseeding bounds this.
3. **Chain.** Group MEMs that are colinear on the same reference contig and within an expected insert-size window into "chains." Chains with the same set of MEMs but different colinear ordering are different chains. Drop chains whose total seed coverage is below a threshold (`-c` flag, default 10000).
4. **Extend.** For each surviving chain, run a banded Smith-Waterman extension across the gaps between MEMs. The banded width is dynamic, set by the chain's gap-and-mismatch profile. This is where the affine-gap penalty matters: default `-O 6 -E 1` is gap-open 6 and gap-extend 1, tuned for Illumina error profiles.
5. **Score.** The final score is the banded SW score across the chain. The top alignment is the primary; alternatives within `-T` of the primary score are reported as secondary. Compute MAPQ from the gap between primary and second-best scores.

Every step is a filter. The seed step drops ~99.9% of the reference (anything that does not match any MEM in the read). The chain step drops most surviving seeds (anything not colinear with at least one other). The extension step runs Smith-Waterman on a tiny fraction of the reference per read.

### 3.1 The four key BWA-MEM parameters

- **`-k INT`** — minimum seed length. Default 19. Smaller `k` means more seeds per read, slower, more sensitive. Larger `k` means fewer seeds, faster, less sensitive. The default is well-tuned for 150 bp Illumina; for shorter reads (50 bp) you may want `-k 15`. Do not lower below 12 unless you have a specific reason.
- **`-r FLOAT`** — reseed trigger. Default 1.5. If a MEM is at least `1.5 × k` long, split into overlapping shorter seeds at that position. Lowering to `1.0` is more aggressive (more reseeding, slower, slightly more sensitive in repeats); raising to `2.0` is faster but worse in repetitive regions.
- **`-c INT`** — skip a MEM if it has more than `INT` occurrences in the reference. Default 500. Lower `-c` (e.g. 100) is faster but loses sensitivity in repeats. Higher `-c` is slower but better in repeats. Set to a much larger value (~10000) when you specifically care about repetitive regions.
- **`-T INT`** — output threshold. Default 30. Alignments with final score below `INT` are dropped. Lowering produces more alignments (including poor ones); raising produces stricter output. This is the closest analog to BLAST's `-evalue` and is a *reporting* parameter, not a search parameter — the search is the same regardless of `-T`.

There are ~30 more flags in `bwa mem --help`. The four above are the ones you will tune. Read the manual.

### 3.2 The `-R` read group is non-negotiable

```bash
bwa mem -t 4 -R "@RG\tID:run1\tSM:ecoli_K12\tLB:lib1\tPL:ILLUMINA" \
  ref/ecoli reads/R1.fq.gz reads/R2.fq.gz > aln/sample.sam
```

The `-R` flag emits a `@RG` (read group) header into the SAM/BAM. **Always set it.** Variant callers require read groups to distinguish samples (different `SM` values) and libraries (different `LB` values, used by duplicate markers to detect PCR duplicates within a library). A BAM without a read group is treated as "unknown sample" by `bcftools mpileup` and `gatk HaplotypeCaller` and will produce confusing output or hard errors downstream.

The fields:

- `ID` — unique within a BAM. The lane or run identifier.
- `SM` — sample name. The biological sample; multiple libraries from the same sample share an `SM`.
- `LB` — library. The DNA library prep; duplicates are detected within a library.
- `PL` — platform. `ILLUMINA`, `PACBIO`, `ONT` (Oxford Nanopore), `MGI`. Drives the duplicate-detection heuristics.
- `PU` (optional) — platform unit. Run + lane + barcode. Required by some pipelines (GATK).

---

## 4. minimap2 in four steps

minimap2 (Li 2018, *Bioinformatics* 34:3094) is the modern replacement for BWA-MEM that handles short reads, long noisy reads (Nanopore, PacBio CLR), long accurate reads (PacBio HiFi), and splice-aware RNA alignment from a single binary. The algorithm differs from BWA-MEM in three places: minimizer seeding instead of MEM seeding, a different chaining strategy, and per-preset parameter tuning.

1. **Minimizer index.** A minimizer is the lexicographically smallest k-mer in a sliding window of `w` consecutive k-mers (Roberts et al. 2004). For default short-read settings (`-k 21 -w 11`), this samples roughly `1/(w+1)` ≈ 8% of the k-mers in the reference, producing a much smaller index than BWA's full FM-index — ~6 GB for a 3 Gb reference vs ~3 GB for BWA, but with simpler access patterns. The minimap2 index format (`.mmi`) is the on-disk version.
2. **Seed.** For each read, extract its minimizers and look up each in the index. This gives a list of (read minimizer position, reference minimizer positions) hits. The number of hits per read is typically 10–50 for short reads, hundreds-to-thousands for long reads.
3. **Chain.** Find the highest-scoring colinear chain of minimizer hits using a 2D dynamic-programming chain algorithm. Long-read chains can span 10s of kilobases; short-read chains span 100–200 bp. The chain score function is tunable per preset.
4. **Extend.** Run a banded SSE2/AVX2-vectorized Smith-Waterman extension across the gaps in the top chain. Output as SAM (`-a` flag) or PAF (default). For long reads, the extension band can be wide (10s of bp); for short reads it is narrow.

### 4.1 minimap2 vs BWA-MEM: when to use each

| Question | Use |
|----------|-----|
| Illumina paired-end 150 bp DNA-seq | Either. `bwa mem` is the historical default and produces near-identical output. `minimap2 -ax sr` is slightly faster (~10%) on modern CPUs. |
| Oxford Nanopore long reads (10 kb median, 5–15% error) | `minimap2 -ax map-ont` |
| PacBio CLR long reads (10 kb median, 10–15% error) | `minimap2 -ax map-pb` |
| PacBio HiFi reads (15 kb median, <1% error) | `minimap2 -ax map-hifi` |
| Spliced cDNA / RNA reads (Iso-Seq, ONT direct cDNA) | `minimap2 -ax splice` |
| Whole-genome assembly-vs-assembly | `minimap2 -ax asm5` / `asm10` / `asm20` |
| Reads < 50 bp (legacy 36–50 bp Illumina) | Neither is ideal. `bwa aln` (legacy) or `bowtie2` are better tuned. |

**Default to minimap2.** It is faster, handles more read types, and is more actively developed (Heng Li still commits to it weekly; BWA gets bug fixes but no new features since 2017). The main reason to use BWA-MEM is reproducibility — large existing pipelines (GATK Best Practices, the 1000 Genomes pipeline) standardize on BWA-MEM, and switching aligners changes the alignments enough to perturb downstream variant calls by a fraction of a percent.

---

## 5. The FASTQ → SAM → sorted BAM pipeline

Here is the canonical short-read alignment pipeline, end to end. You will run it dozens of times this semester; memorize the steps.

### 5.1 Step 1 — Index the reference (one-time)

```bash
# E. coli K-12 MG1655 reference (4.6 Mb, ~5 seconds to index).
bwa index ref/ecoli.fa
# Produces: ref/ecoli.fa.amb, ref/ecoli.fa.ann, ref/ecoli.fa.bwt,
#           ref/ecoli.fa.pac, ref/ecoli.fa.sa.

# Or with minimap2 (faster indexing; smaller .mmi for short-read settings).
minimap2 -d ref/ecoli.mmi -x sr ref/ecoli.fa
```

The BWA index is ~3.5x the size of the input FASTA. For a 4.6 Mb FASTA, the index is ~12 MB total across the 5 files. For a 3 Gb human reference, the index is ~10 GB and takes ~90 minutes to build (`bwa index -a bwtsw` is required for genomes ≥ 2 GB).

### 5.2 Step 2 — Align and pipe to sort

```bash
bwa mem -t 4 \
  -R "@RG\tID:SRR1770413\tSM:ecoli_K12\tLB:lib1\tPL:ILLUMINA" \
  ref/ecoli.fa \
  reads/SRR1770413_1.fq.gz reads/SRR1770413_2.fq.gz \
| samtools sort -@ 4 -o aln/SRR1770413.sorted.bam -
```

The pipe (`|`) is intentional. `bwa mem` writes SAM to stdout, `samtools sort` reads it and writes coordinate-sorted BAM. Avoiding an intermediate SAM file on disk saves I/O and ~15 GB of temporary space.

### 5.3 Step 3 — Index the BAM

```bash
samtools index aln/SRR1770413.sorted.bam
# Produces aln/SRR1770413.sorted.bam.bai
```

The `.bai` is a binary index over the BAM file letting you do range queries (`samtools view aln.bam chr1:1000-2000`) in milliseconds instead of streaming the whole file. **Every downstream tool requires this index.** If you forget it, the next call will yell at you with a moderately helpful error message.

### 5.4 Step 4 — Sanity check with flagstat and depth

```bash
samtools flagstat aln/SRR1770413.sorted.bam
# 2,500,000 + 0 in total (QC-passed reads + QC-failed)
# 0       + 0 secondary
# 12,345  + 0 supplementary
# 0       + 0 duplicates (we have not marked them yet)
# 2,488,765 + 0 mapped (99.55%)
# 2,500,000 + 0 paired in sequencing
# 1,250,000 + 0 read1
# 1,250,000 + 0 read2
# 2,460,123 + 0 properly paired (98.40%)
# ...

samtools depth -a aln/SRR1770413.sorted.bam | awk '{sum+=$3} END {print sum/NR}'
# 47.3   (mean coverage across all reference positions)
```

If `flagstat` shows `> 90%` mapped and `> 95%` properly paired, the alignment is healthy. If mapped drops below 80%, something is wrong — usually wrong reference, wrong organism, or massive adapter contamination. If properly-paired drops below 90% but mapped is fine, the insert-size distribution is unusual (very short or very long fragments, or chimeric library prep).

### 5.5 Step 5 — Mark duplicates (covered in Lecture 2)

```bash
# Required preparation: ensure mate-coordinate fields are filled.
samtools sort -n -@ 4 aln/SRR1770413.sorted.bam | \
  samtools fixmate -m - - | \
  samtools sort -@ 4 -o aln/SRR1770413.fixmate.bam - && \
  samtools markdup aln/SRR1770413.fixmate.bam aln/SRR1770413.markdup.bam && \
  samtools index aln/SRR1770413.markdup.bam
```

The pipeline `sort -n | fixmate -m | sort | markdup` is the canonical samtools idiom for duplicate marking. `samtools markdup` itself requires that mates have been paired up by `fixmate` first. The result is a BAM where duplicates carry the `0x400` flag but are *not* removed — variant callers will ignore them, coverage plots will count or skip them per your choice. Lecture 2 covers the details.

---

## 6. A worked sanity check in Python

You will not implement BWA from scratch — that is a semester project. But you should *run* the pipeline end to end on a tiny example so the abstract algorithm becomes concrete. Here is the smallest possible Python+shell sanity check on the bacteriophage lambda reference (48.5 kb) and 1000 simulated 150 bp reads:

```python
from __future__ import annotations
import subprocess
from pathlib import Path

import pysam

# Assume ref/lambda.fa exists; if not, fetch it once:
#   wget -O ref/lambda.fa.gz \
#     https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=fasta
#   gunzip ref/lambda.fa.gz

ref = Path("ref/lambda.fa")
r1 = Path("reads/lambda_R1.fq.gz")
r2 = Path("reads/lambda_R2.fq.gz")
bam = Path("aln/lambda.sorted.bam")
bam.parent.mkdir(parents=True, exist_ok=True)

# Step 1: index the reference (if not already done).
if not (ref.parent / f"{ref.name}.bwt").exists():
    subprocess.run(["bwa", "index", str(ref)], check=True)

# Step 2: align + pipe to sort. Use shell=True for the pipe.
cmd = (
    f"bwa mem -t 2 "
    f"-R '@RG\\tID:lambda_test\\tSM:lambda\\tLB:lib1\\tPL:ILLUMINA' "
    f"{ref} {r1} {r2} "
    f"| samtools sort -@ 2 -o {bam} -"
)
subprocess.run(cmd, shell=True, check=True)

# Step 3: index the BAM.
subprocess.run(["samtools", "index", str(bam)], check=True)

# Step 4: sanity check with pysam.
af = pysam.AlignmentFile(str(bam), "rb")
n_total = af.count()
af.close()
af = pysam.AlignmentFile(str(bam), "rb")
n_mapped = sum(1 for read in af if not read.is_unmapped)
af.close()
print(f"Total reads: {n_total}")
print(f"Mapped: {n_mapped} ({100*n_mapped/n_total:.1f}%)")
```

Expected output for a clean 1000-read paired-end simulation against lambda:

```
Total reads: 2000
Mapped: 1998 (99.9%)
```

If your numbers are wildly different (mapped < 90%), the typical causes are: (a) the reference is wrong (you indexed `lambda.fa` but the reads are from `E. coli`), (b) the reads are corrupted (truncated FASTQ), or (c) the FASTQ encoding is non-standard (some old datasets are Phred+64 instead of Phred+33).

---

## 7. Choice criteria: which aligner for which read profile

Once you can run both `bwa mem` and `minimap2`, the next question is *when* to use each. The decision tree:

```
Read length?
├── < 50 bp (legacy)               → bowtie2 (out of scope; not in this course)
├── 50–500 bp paired Illumina      → bwa mem  OR  minimap2 -ax sr
│                                    (pick by team convention; both work)
├── 1–100 kb Oxford Nanopore       → minimap2 -ax map-ont
├── 1–100 kb PacBio CLR (legacy)   → minimap2 -ax map-pb
├── 1–100 kb PacBio HiFi (modern)  → minimap2 -ax map-hifi
└── Spliced cDNA / Iso-Seq / RNA   → minimap2 -ax splice  OR  STAR (for short-read RNA)
```

If you have *short-read DNA* and no strong reason to choose otherwise, the field's default is still `bwa mem` for legacy/pipeline reasons. New code should prefer `minimap2 -ax sr` — the output is essentially equivalent (97%+ of reads receive identical alignments), it is ~10% faster, and it is one binary for everything you might align next.

The rest of this lecture series will use `bwa mem` in the mini-project (to match the field's variant-calling convention, since Week 6 will run `bcftools` against a BWA-MEM BAM exactly as the GATK Best Practices do) and `minimap2 -ax sr` in the homework (so you have hands-on experience with both).

---

## 8. Common misconceptions

A short list of "things that seem right but are not":

- **"`bwa mem` finds the optimal alignment."** No. `bwa mem` finds optimal-or-near-optimal alignments via a heuristic seed-and-extend filter. The extension itself is Smith-Waterman on the windowed band, but the *seed* step can miss a true alignment if the read's MEMs all fail (every k-mer in the read happens to have many false hits in the reference). In practice, for reads of ≥ 70 bp at < 5% error, BWA-MEM misses essentially nothing.
- **"A MAPQ of 60 means the alignment is perfect."** No. MAPQ = 60 means *the placement* is confidently unique (P[misalignment] ≤ 10^-6), but the *alignment* may still contain mismatches and indels. The CIGAR string tells you about the alignment; MAPQ tells you about the position. These are different concepts.
- **"Soft clips mean the read is bad."** No. Soft clips at read ends are normal and usually mean (a) the read partially extends past an indel that bwa-mem chose not to align through, or (b) the read partially extends into adapter sequence not trimmed upstream. Pervasive soft-clipping at internal positions is a red flag, but soft clips at ends are routine.
- **"More coverage is always better."** Up to a point. For variant calling, 30x mean coverage is the standard for germline calling; somatic calling wants 100x+. Above 100x, additional coverage produces diminishing returns and can mask PCR-duplicate-driven artifacts. For 1000x+ panel sequencing, you genuinely need that depth (low-allele-fraction variant detection); for whole-genome germline, you do not.
- **"Duplicate marking removes reads."** No, by default `samtools markdup` and `picard MarkDuplicates` *mark* duplicates with flag `0x400` and leave them in the BAM. Variant callers then choose to ignore them at the SNP-calling step. Use `samtools markdup -r` if you actually want to remove them (and lose information you might want later — e.g. for coverage QC, which should count duplicates and non-duplicates separately).

If any of these surprised you, re-read sections 3, 4, and 5.

---

## 9. Where this lecture lands you for Lecture 2

After this lecture you should be able to:

- Describe the FM-index in two paragraphs (BWT + rank/select + backward search).
- Run `bwa index ref.fa`, `bwa mem ref.fa r1.fq.gz r2.fq.gz`, and `samtools sort` end to end.
- Choose between `bwa mem` and `minimap2` for a given read profile.
- Name the four BWA-MEM seed-and-extend parameters (`-k`, `-r`, `-c`, `-T`) and what each controls.
- Distinguish the *seed step* (FM-index lookup, sub-linear) from the *extension step* (windowed Smith-Waterman, linear in read length).

Lecture 2 takes the SAM output and zooms into the format: the eleven mandatory columns, the FLAG bit field, the CIGAR string, the MAPQ scale, and the `samtools` toolchain for sorting, indexing, filtering, and computing coverage and duplicates. By the end of Lecture 2 you will be able to read a SAM line aloud and translate it into "this 150 bp read from sample S1 mapped uniquely (MAPQ 60) to position 1234567 on the forward strand, with a 2 bp insertion at position 67 in the read."

---

## Self-check questions

Before you move on, answer these without looking. If you cannot answer one, re-read the relevant section.

1. State the four steps of the seed-and-extend alignment pipeline as implemented by BWA-MEM. (§3)
2. The FM-index lets you locate every occurrence of a length-`m` query in a length-`n` reference in `O(?)` time. Fill in the question mark and explain why. (§2.3)
3. What is a Maximal Exact Match (MEM)? How is it found in the FM-index? (§3, §2.3)
4. For a 150 bp paired-end Illumina dataset against the *E. coli* MG1655 reference, which aligner would you reach for first? Why? (§4.1)
5. Decode the BWA-MEM command `bwa mem -t 4 -R '@RG\tID:run1\tSM:S1\tLB:lib1\tPL:ILLUMINA' ref/ecoli.fa r1.fq.gz r2.fq.gz`. What does each flag do? What goes wrong if `-R` is omitted? (§3.2)
6. What does `MAPQ = 60` mean? `MAPQ = 0`? Are they on a linear or a log scale? (§7, see Lecture 2 §3)
7. Why does the canonical alignment pipeline use a pipe (`bwa mem ... | samtools sort -`) rather than writing an intermediate SAM file? (§5.2)
8. What is the difference between `bwa index` and `samtools faidx`? When do you use each? (§5.1, §5.4)
9. Name two read-type scenarios where you must use `minimap2` instead of `bwa mem`. (§4.1)
10. Under what biological scenario would the seed step of BWA-MEM miss a real homologous alignment? (§3, §8)

Answers are not provided. If you struggle, the answers are in the section references; do the work.

---

## Further reading

- Li, H. & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* 25(14):1754–1760.
- Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv:1303.3997*.
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094–3100.
- Ferragina, P. & Manzini, G. (2000). Opportunistic data structures with applications. *FOCS 2000*:390–398.
- Burrows, M. & Wheeler, D. J. (1994). A block sorting lossless data compression algorithm. *Digital SRC Research Report* 124.
- The SAM/BAM specification: <https://samtools.github.io/hts-specs/SAMv1.pdf>.
- The BWA source code: <https://github.com/lh3/bwa>.
- The minimap2 source code: <https://github.com/lh3/minimap2>.

---

*Continue to [Lecture 2 — SAM/BAM format and samtools](./02-sam-bam-format-and-samtools.md) once you have answered the self-check questions.*
