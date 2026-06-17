# Week 5 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (BWA 0.7.17, minimap2 2.26, samtools 1.19, pysam 0.22), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Li & Durbin (2009)** — the original BWA paper. *Bioinformatics* 25:1754. Free full text:
  <https://academic.oup.com/bioinformatics/article/25/14/1754/225615>
- **Li (2013)** — BWA-MEM, the "Maximal Exact Match" seed-and-extend mode that replaced BWA-backtrack for reads ≥ 70 bp. *arXiv:1303.3997*. Free preprint:
  <https://arxiv.org/abs/1303.3997>
- **Li (2018)** — minimap2, the long-read and split-read replacement that subsumes BWA-MEM for reads of any length. *Bioinformatics* 34:3094. Free full text:
  <https://academic.oup.com/bioinformatics/article/34/18/3094/4994778>
- **Li et al. (2009)** — the SAM/BAM format paper. *Bioinformatics* 25:2078. Free full text:
  <https://academic.oup.com/bioinformatics/article/25/16/2078/204688>
- **Ferragina & Manzini (2000)** — the FM-index data structure. *FOCS 2000*. Free preprint widely available:
  <https://people.unipmn.it/manzini/papers/focs00draft.pdf>
- **Burrows & Wheeler (1994)** — the BWT, the transform that makes the FM-index possible. *Digital SRC Research Report* 124. Free PDF:
  <https://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.pdf>
- **The SAM/BAM/CRAM specification** — the canonical reference. ~25 pages, surprisingly readable:
  <https://samtools.github.io/hts-specs/SAMv1.pdf>

## Tool reference (the command-line surface)

### BWA 0.7.17

Cheat-sheet of the commands you will use this week. All ship with BWA 0.7.17.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `bwa index` | Build an FM-index from a reference FASTA | `-p` (output prefix), `-a` (`bwtsw` for genomes ≥ 2 GB, `is` for smaller) |
| `bwa mem` | Align reads ≥ 70 bp against an indexed reference | `-t` (threads), `-R` (read-group tag), `-M` (flag secondary as supplementary), `-K` (chunk size for reproducibility) |
| `bwa aln` + `bwa samse` / `bwa sampe` | Legacy "backtrack" mode for reads < 70 bp. Not used in this course. | — |
| `bwa fastmap` | Find super-maximal exact matches without alignment | — |

#### Building an index

```bash
bwa index -p ref/ecoli ref/ecoli.fa
# Produces ref/ecoli.amb, ref/ecoli.ann, ref/ecoli.bwt, ref/ecoli.pac, ref/ecoli.sa.
```

For a 4.6 Mb E. coli reference, this takes ~5 seconds. For a 3 Gb human reference, it takes ~90 minutes.

#### Aligning paired-end reads

```bash
bwa mem -t 4 \
  -R "@RG\tID:run1\tSM:ecoli_K12\tLB:lib1\tPL:ILLUMINA" \
  ref/ecoli \
  reads/R1.fq.gz reads/R2.fq.gz \
  > aln/sample.sam
```

The `-R` flag emits a `@RG` (read group) header. **Always set it.** Downstream tools (variant callers, duplicate markers) require it to distinguish samples and libraries.

### minimap2 2.26

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `minimap2 -d` | Build a minimizer index from a reference (optional — minimap2 can index on the fly) | `-k` (k-mer size), `-w` (minimizer window) |
| `minimap2 -ax <preset>` | Align reads with output as SAM | `-ax sr` (short reads), `-ax map-ont` (Oxford Nanopore), `-ax map-pb` (PacBio CLR), `-ax map-hifi` (PacBio HiFi), `-ax splice` (RNA splice-aware) |

#### Presets — read this table once, refer to it always

| Preset | Read type | Expected error rate | Use case |
|--------|-----------|--------------------:|----------|
| `-ax sr` | Illumina paired-end short reads | < 1% | DNA-seq, equivalent to bwa mem |
| `-ax map-ont` | Oxford Nanopore | 5–15% | Long-read DNA |
| `-ax map-pb` | PacBio CLR | 10–15% | Long-read DNA, older PacBio |
| `-ax map-hifi` | PacBio HiFi | < 1% | Long-accurate-read DNA, recent |
| `-ax splice` | RNA, splice-aware | varies | Spliced cDNA / Iso-Seq |
| `-ax asm5` / `-ax asm10` / `-ax asm20` | Assembly-vs-assembly | 5/10/20% divergence | Whole-genome comparisons |

#### Aligning long reads

```bash
minimap2 -ax map-ont -t 4 ref.fa nanopore_reads.fq.gz > aln.sam
```

minimap2 is roughly 4–5x faster than bwa mem on short reads and ~30x faster than BWA-SW on long reads at equal or better accuracy. It is the modern default for any read length.

### samtools 1.19

The `samtools` suite is ~30 subcommands that together let you read, write, sort, index, filter, and summarize SAM/BAM/CRAM files.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `samtools view` | Read/filter SAM/BAM/CRAM | `-b` (output BAM), `-S` (input SAM, optional in recent versions), `-f`/`-F` (require/exclude flag bits), `-q` (MAPQ threshold), `-h` (include header in text output) |
| `samtools sort` | Sort BAM by coordinate (or `-n` by name) | `-o` (output), `-@` (threads), `-n` (name-sort), `-T` (temp file prefix) |
| `samtools index` | Build a `.bai` index for a coordinate-sorted BAM | — |
| `samtools fixmate` | Fill in mate-coordinate fields (needed before markdup) | `-m` (add MS/ms tags for markdup) |
| `samtools markdup` | Mark PCR/optical duplicates with flag `0x400` | `-r` (remove instead of mark), `-s` (print stats) |
| `samtools flagstat` | Print SAM flag summary | — |
| `samtools depth` | Per-position read depth | `-a` (include zero-coverage positions), `-r` (region) |
| `samtools coverage` | Per-contig coverage summary | `-r` (region) |
| `samtools idxstats` | Per-contig read count from BAM index | — |
| `samtools faidx` | Index a FASTA for random access | — |
| `samtools tview` | Terminal text viewer for a BAM file | `-d T` (text mode) |

#### The canonical FASTQ-to-BAM one-liner

```bash
bwa mem -t 4 ref/ecoli reads/R1.fq.gz reads/R2.fq.gz \
| samtools sort -@ 4 -o aln/sample.sorted.bam -
samtools index aln/sample.sorted.bam
```

That pipe is the most-typed sequence of commands in short-read bioinformatics. Memorize it.

### pysam 0.22

`pysam` is the Python wrapper around the C `htslib` library. It is the way to read BAM files programmatically without parsing text. The API is small:

| Class / function | Purpose |
|------------------|---------|
| `pysam.AlignmentFile(path, mode)` | Open a SAM/BAM/CRAM file. Mode `"rb"` for BAM, `"r"` for SAM. |
| `AlignmentFile.fetch(contig, start, stop)` | Iterate over reads overlapping a region (requires `.bai`) |
| `AlignmentFile.pileup(contig, start, stop)` | Iterate over pileup columns (one per reference position) |
| `AlignedSegment.query_name`, `.flag`, `.reference_name`, `.reference_start`, `.mapping_quality`, `.cigarstring`, `.cigartuples`, `.query_sequence`, `.query_qualities` | The eleven SAM fields plus parsed CIGAR. |
| `AlignedSegment.is_paired`, `.is_unmapped`, `.is_reverse`, `.is_duplicate`, `.is_secondary`, `.is_supplementary` | Boolean accessors for common FLAG bits. |

#### A 10-line coverage computation

```python
import pysam
import numpy as np

bam = pysam.AlignmentFile("aln/sample.sorted.bam", "rb")
ref_length = bam.lengths[0]
coverage = np.zeros(ref_length, dtype=np.int32)
for column in bam.pileup("NC_000913.3", 0, ref_length, truncate=True):
    coverage[column.reference_pos] = column.nsegments
print(f"Mean coverage: {coverage.mean():.1f}x, median: {np.median(coverage):.0f}x")
```

On a 4.6 Mb E. coli BAM at 40x coverage, this runs in ~30 seconds. For a human-genome BAM it is the wrong approach (too slow); use `samtools depth` and read the TSV.

## Reference dataset accessions

Cited by NCBI / SRA accession so you can verify your data is the same as the curriculum's:

- **`NC_000913.3`** — *Escherichia coli* str. K-12 substr. MG1655 complete genome (4,641,652 bp). The canonical bacterial reference. From Week 4; we re-use the FASTA.
- **`SRR1770413`** — Illumina HiSeq 2500 paired-end resequencing of *E. coli* K-12 MG1655, ~5 Gb compressed. The mini-project read set.
- **`SRR2014925`** — Oxford Nanopore MinION reads of *E. coli* K-12 MG1655. Used in the homework Problem 4 minimap2 comparison.
- **`SRR8482586`** — Illumina paired-end reads of *Mycobacterium tuberculosis* H37Rv, ~1 Gb. Used in homework Problem 5 (low-complexity / GC-skewed reference).
- **`NC_001416.1`** — Bacteriophage lambda complete genome (48,502 bp). A toy reference used in exercises so a full pipeline runs in seconds.

If any of these accessions have been retired or updated by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt.

## SAM flag bits — the ones you must memorize

| Hex | Decimal | Name | Meaning |
|----:|--------:|------|---------|
| `0x1` | 1 | `PAIRED` | Read is paired (the run is paired-end) |
| `0x2` | 2 | `PROPER_PAIR` | Read is mapped in a proper pair (orientation + insert size sensible) |
| `0x4` | 4 | `UNMAP` | Read is unmapped |
| `0x8` | 8 | `MUNMAP` | Mate is unmapped |
| `0x10` | 16 | `REVERSE` | Read is reverse-complemented (mapped to the minus strand) |
| `0x20` | 32 | `MREVERSE` | Mate is reverse-complemented |
| `0x40` | 64 | `READ1` | This is the first read of the pair |
| `0x80` | 128 | `READ2` | This is the second read of the pair |
| `0x100` | 256 | `SECONDARY` | Secondary alignment (read multi-mapped, this is not the primary) |
| `0x200` | 512 | `QCFAIL` | Read failed vendor QC |
| `0x400` | 1024 | `DUP` | PCR or optical duplicate |
| `0x800` | 2048 | `SUPPLEMENTARY` | Supplementary alignment (chimeric — read spans two locations) |

Examples to internalize:

- `FLAG = 99` = `0x1 + 0x2 + 0x20 + 0x40` = paired, proper-paired, mate-reverse, first-in-pair. The most common flag for the forward read of a healthy paired-end pair.
- `FLAG = 147` = `0x1 + 0x2 + 0x10 + 0x80` = paired, proper-paired, read-reverse, second-in-pair. The most common flag for the reverse read of a healthy paired-end pair.
- `FLAG = 4` = `0x4` only = unpaired, unmapped read. Common in a small fraction of any real dataset.
- `FLAG = 1024` = `0x400` only = unpaired duplicate (rare). Combined with other bits is common in deduped paired-end data (`FLAG = 1123` = `0x1 + 0x2 + 0x40 + 0x400 + 0x80` for example).

A useful interactive decoder lives at <https://broadinstitute.github.io/picard/explain-flags.html>. Use it to check yourself.

## CIGAR operations — the table

| Op | Long name | Consumes query? | Consumes reference? | Notes |
|----|-----------|:----------------:|:-----------------:|-------|
| `M` | Alignment match | yes | yes | Could be sequence match *or* mismatch. Common in older SAM. |
| `I` | Insertion to reference | yes | no | Bases in read not in reference. |
| `D` | Deletion from reference | no | yes | Bases in reference not in read. |
| `N` | Skipped region | no | yes | E.g. an intron in RNA-seq. |
| `S` | Soft clip | yes | no | Bases at read ends not aligned; still in `SEQ`. |
| `H` | Hard clip | no | no | Bases not in `SEQ` at all (only legal at read ends). |
| `P` | Padding | no | no | Silent deletion against padded reference (rare). |
| `=` | Sequence match | yes | yes | Modern aligners emit this when asked. |
| `X` | Sequence mismatch | yes | yes | Modern aligners emit this when asked. |

A CIGAR of `36M2I12M1D80M` decoded:

- `36M` — 36 aligned bases (matches or mismatches).
- `2I` — 2 inserted bases (in read, not in reference).
- `12M` — 12 aligned bases.
- `1D` — 1 deleted base (in reference, not in read).
- `80M` — 80 aligned bases.

Implied query length = `36 + 2 + 12 + 0 + 80 = 130` bases. Implied reference span = `36 + 0 + 12 + 1 + 80 = 129` bases. (Soft clips would add to query but not reference; hard clips add to neither.)

## MAPQ — the Phred-scaled mapping quality

MAPQ is `-10·log10(P[read is misaligned to this position])`. The integer MAPQ corresponds to a misalignment probability:

| MAPQ | P_wrong | Meaning |
|-----:|--------:|---------|
| 60 | 10^-6 | BWA-MEM ceiling — unique, high-confidence placement |
| 40 | 10^-4 | Good placement; second-best position significantly worse |
| 30 | 10^-3 | Reasonable placement; some ambiguity |
| 20 | 10^-2 | Marginal placement; second-best within ~10 of best |
| 10 | 10^-1 | Poor placement; many candidates |
| 0 | 1 | Multimapper or unaligned; no unique placement |

**Critical caveats**:

- Different aligners use different MAPQ scales. `bwa mem` caps at 60, `minimap2` at 60, `bowtie2` at 42, `STAR` at 255 (with 255 = "this many alignments, unique" in STAR's scheme, *not* infinity). When reporting MAPQ thresholds, name the aligner.
- `MAPQ = 0` does not mean "no alignment." It usually means "this read aligned equally well to ≥ 2 positions, so I cannot pick one." Multimapper status, not failure.
- `MAPQ = 255` historically meant "unavailable" in the SAM spec but some aligners (STAR) overloaded it. Read the aligner manual.

## Free books, chapters, and tutorials

- **Heng Li's blog** — the author of BWA, minimap2, samtools. His posts often explain corner cases that the papers do not:
  <https://lh3.github.io/>
- **The samtools-bcftools-htslib tutorials** at <http://www.htslib.org/doc/> — official, version-pinned, kept up to date with each release.
- **GATK Best Practices for short variant discovery** — the Broad's canonical pipeline reference. The "Mark Duplicates" and "Base Quality Score Recalibration" pages are essential reading even if you do not plan to use GATK:
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035535912>
- **Galaxy Training Network — "Mapping" tutorial** — a hands-on web-based introduction with diagrams:
  <https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html>
- **Bioinformatics Algorithms: An Active Learning Approach (Compeau & Pevzner)** — Chapter 9 covers the BWT and FM-index with diagrams and worked examples; free chapters at:
  <https://www.bioinformaticsalgorithms.org/>

## Public BAM datasets to read for practice

- **The 1000 Genomes Project phase 3 BAMs** — high-coverage and low-coverage human BAMs from publicly consented samples. The HG00096 individual is the textbook "first sample" — try a 1 Mb region from chromosome 22:
  <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00096/alignment/>
- **GIAB (Genome in a Bottle)** — high-confidence BAMs with truth-set variants. Sample NA12878 (HG001) is the gold standard:
  <https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/>
- **The Heng Li example BAM** — a 100-read SAM file Heng distributes for tutorials, useful for SAM-by-hand parsing:
  <https://github.com/samtools/samtools/blob/develop/examples/toy.sam>

## Open-source code to read this week

You will learn more from one hour reading the BWA-MEM source than from three hours of tutorials. The BWA codebase is C; the parts worth reading are small:

- **`bwa/bwamem.c`** — the BWA-MEM main loop: seed (MEM finding via the FM-index), chain, extend. ~3,000 lines, well-commented:
  <https://github.com/lh3/bwa/blob/master/bwamem.c>
- **`bwa/bwt.c`** — the BWT-based suffix array index. ~500 lines:
  <https://github.com/lh3/bwa/blob/master/bwt.c>
- **`minimap2/minimap.h` + `minimap2/map.c`** — the minimap2 chaining and extension. Read `minimap.h` first for the data structures, then `map.c` for the algorithm:
  <https://github.com/lh3/minimap2>
- **`samtools/bam_markdup.c`** — the markdup implementation. ~700 lines. The Picard MarkDuplicates equivalent is much longer but the logic is the same:
  <https://github.com/samtools/samtools/blob/develop/bam_markdup.c>
- **`pysam`** — the Python bindings to htslib. The `pysam/calignmentfile.pyx` Cython file is the surface you call from Python:
  <https://github.com/pysam-developers/pysam>

## Coverage thresholds — a rough rule-of-thumb table

What "enough" coverage looks like depends on the question. A rough guide:

| Mean coverage | Use cases | Notes |
|--------------:|-----------|-------|
| 1–5x | Population SNP detection (low-coverage 1000 Genomes) | Per-position genotypes unreliable; rely on linkage and statistical phasing |
| 10–15x | Germline SNP calling, well-behaved reference | Below the practical floor for indels |
| 30x | Germline whole-genome variant calling | The "industry standard" for clinical sequencing |
| 100x+ | Somatic variant calling (tumor) | Tumor heterogeneity demands depth to detect low-allele-fraction variants |
| 500x+ | Targeted panels, ctDNA | Detection of variants present at < 1% allele fraction |

For an E. coli reference at 4.6 Mb with 1 Gb of 150 bp paired-end reads, mean coverage = `(1·10^9 bases) / (4.6·10^6 bp)` ≈ 217x. That is excessive for the mini-project; we will subsample.

## Tools you will install this week

- **BWA 0.7.17** — `conda install -c bioconda bwa=0.7.17` (or `brew install bwa` on macOS). Adds `bwa` to your PATH.
- **minimap2 2.26** — `conda install -c bioconda minimap2=2.26`. Adds `minimap2` to your PATH.
- **samtools 1.19** — `conda install -c bioconda samtools=1.19`. Adds `samtools` and `htsfile` to your PATH.
- **pysam 0.22** — `pip install pysam==0.22` or `conda install -c bioconda pysam=0.22`. Python bindings to htslib.
- **(SRA toolkit)** — `conda install -c bioconda sra-tools=3.0.10`. Adds `prefetch`, `fasterq-dump`. Needed only for downloading SRA datasets; if you already have the FASTQ.gz files, skip this.

### A complete environment file for the week

```yaml
name: c10-week-05
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - numpy=1.26.4
  - matplotlib
  - pandas
  - biopython=1.83
  - bwa=0.7.17
  - minimap2=2.26
  - samtools=1.19
  - bcftools=1.19
  - pysam=0.22
  - sra-tools=3.0.10
  - pip
```

`conda env create -f env.yml && conda activate c10-week-05`. Total size on disk after install: ~1.2 GB.

---

*If a link 404s, please open an issue so we can replace it.*
