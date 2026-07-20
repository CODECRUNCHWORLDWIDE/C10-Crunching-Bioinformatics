# Week 10 — Resources

> **Reproducibility note.** Every tool, file format, and reference paper on this page is free and publicly accessible. Where we name a version (Flye 2.9.5, Canu 2.2, Hifiasm 0.19.9, Medaka 1.12.0, BUSCO 5.7.1, badread 0.4.1, Dorado 0.7.2, Bandage 0.9.0), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Kolmogorov, Mikhail; Yuan, Jeffrey; Lin, Yu; Pevzner, Pavel A. (2019)** — "Assembly of long, error-prone reads using repeat graphs." The Flye paper. *Nature Biotechnology* 37:540. Free full text via PubMed Central:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699608/>
  Tool documentation:
  <https://github.com/fenderglass/Flye>
  The repeat-graph construction in §2 is the technical core. Read it once end-to-end; it is ~10 pages and is the most readable modern long-read assembly paper.
- **Koren, Sergey; Walenz, Brian P.; Berlin, Konstantin; Miller, Jason R.; Bergman, Nicholas H.; Phillippy, Adam M. (2017)** — "Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation." The Canu paper. *Genome Research* 27:722. Free full text:
  <https://genome.cshlp.org/content/27/5/722.full>
  Free PMC mirror:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411767/>
  Tool documentation:
  <https://canu.readthedocs.io/en/latest/>
- **Cheng, Haoyu; Concepcion, Gregory T.; Feng, Xiaowen; Zhang, Haowen; Li, Heng (2021)** — "Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm." The Hifiasm paper. *Nature Methods* 18:170. Free full text:
  <https://www.nature.com/articles/s41592-020-01056-5>
  Free PMC mirror:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8005227/>
  Tool documentation:
  <https://github.com/chhylp123/hifiasm>
- **Wenger, Aaron M.; Peluso, Paul; Rowell, William J.; Chang, Pi-Chuan; Hall, Richard J.; Concepcion, Gregory T.; Ebler, Jana; Fungtammasan, Arkarachai; Kolesnikov, Alexey; Olson, Nathan D.; et al. (2019)** — "Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome." The HiFi paper. *Nature Biotechnology* 37:1155. Free full text via PubMed Central:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6776680/>
  The CCS protocol section explains the chemistry behind PacBio HiFi reads and is the canonical reference for "what is a HiFi read."
- **Manni, Mose; Berkeley, Matthew R.; Seppey, Mathieu; Simao, Felipe A.; Zdobnov, Evgeny M. (2021)** — "BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes." The BUSCO paper. *Molecular Biology and Evolution* 38:4647. Free full text:
  <https://academic.oup.com/mbe/article/38/10/4647/6329644>
  Tool documentation:
  <https://busco.ezlabgyk.unige.ch/> and <https://busco.ezlabgyk.unige.ch/busco_userguide.html>
- **Wick, Ryan R.; Schultz, Mark B.; Zobel, Justin; Holt, Kathryn E. (2015)** — "Bandage: interactive visualization of de novo genome assemblies." The Bandage paper. *Bioinformatics* 31:3350. Free full text:
  <https://academic.oup.com/bioinformatics/article/31/20/3350/196114>
  Tool documentation:
  <https://rrwick.github.io/Bandage/>
- **Baid, Gunjan; Cook, Daniel E.; Shafin, Kishwar; Yun, Taedong; Llinares-Lopez, Felipe; Berthet, Quentin; Belyaeva, Anastasiya; Topfer, Armin; Wenger, Aaron M.; Rowell, William J.; Yang, Howard; et al. (2023)** — "DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer." The DeepConsensus paper. *Nature Biotechnology* 41:232. Free full text:
  <https://www.nature.com/articles/s41587-022-01435-7>
  Tool documentation:
  <https://github.com/google/deepconsensus>
- **Wick, Ryan R. (2019)** — "Badread: simulation of error-prone long reads." The badread paper. *Journal of Open Source Software* 4:1316. Free full text:
  <https://joss.theoj.org/papers/10.21105/joss.01316>
  Tool documentation:
  <https://github.com/rrwick/Badread>
- **Cock, Peter J. A.; Antao, Tiago; Chang, Jeffrey T.; Chapman, Brad A.; Cox, Cymon J.; Dalke, Andrew; Friedberg, Iddo; Hamelryck, Thomas; Kauff, Frank; Wilczynski, Bartek; de Hoon, Michiel J. L. (2009)** — "Biopython: freely available Python tools for computational molecular biology and bioinformatics." The Biopython paper. *Bioinformatics* 25:1422. Free full text:
  <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>
  Tutorial (the canonical reference for `Bio.SeqIO` and `Bio.SeqRecord` on long-read FASTQ):
  <https://biopython.org/docs/latest/Tutorial/index.html>
- **Li, Heng (2018)** — "Minimap2: pairwise alignment for nucleotide sequences." The minimap2 paper. *Bioinformatics* 34:3094. Free full text:
  <https://academic.oup.com/bioinformatics/article/34/18/3094/4994778>
  Tool documentation:
  <https://github.com/lh3/minimap2>
  Minimap2 is the long-read mapper underneath Medaka, Hifiasm's overlap stage, and the assembly QC pipeline that aligns the assembly back to the reference.
- **Walker, Bruce J.; Abeel, Thomas; Shea, Terrance; Priest, Margaret; Abouelliel, Amr; Sakthikumar, Sharadha; Cuomo, Christina A.; Zeng, Qiandong; Wortman, Jennifer; Young, Sarah K.; Earl, Ashlee M. (2014)** — "Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement." The Pilon paper. *PLoS ONE* 9:e112963. Free full text via PubMed Central:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4237348/>
  Tool documentation:
  <https://github.com/broadinstitute/pilon>

## Tool reference (the command-line surface)

### Flye 2.9.5

Flye is the default modern long-read assembler. It is a C++ program with a Python front-end; the conda package installs the `flye` command on the PATH.

| Flag | Purpose |
|------|---------|
| `--nano-hq <reads.fastq>` | High-quality nanopore input (Dorado SUP / HAC basecalled R10.4.1; QV ~17+). The current default. |
| `--nano-raw <reads.fastq>` | Legacy raw nanopore input (Guppy-basecalled R9.4.1; QV ~10-15). |
| `--pacbio-hifi <reads.fastq>` | PacBio HiFi input (QV ~30+). Flye-on-HiFi works but Hifiasm is preferred. |
| `--pacbio-raw <reads.fastq>` | PacBio CLR (continuous long reads; pre-HiFi era). |
| `--genome-size <size>` | Estimated genome size, e.g. `1m`, `5m`, `100m`. Used for read coverage estimation. |
| `--out-dir <path>` | Output directory (created if missing). |
| `--threads N` | Worker thread count. |
| `--iterations N` | Polishing iterations Flye runs internally (default 1; 0 disables Flye's polish). |
| `--min-overlap N` | Minimum read-pair overlap length (default depends on input mode). |
| `--meta` | Metagenome assembly mode (handles uneven coverage across multiple species). |
| `--keep-haplotypes` | Keep heterozygous bubbles instead of collapsing (diploid mode). |

#### The canonical Flye call (Week 10 default)

```bash
flye \
    --nano-hq reads.fastq \
    --genome-size 1m \
    --out-dir flye_out \
    --threads 4 \
    --iterations 1
```

For a 50x simulated nanopore read set on a 1 Mb reference this runs in ~30-60 seconds on a 4-core laptop. The output files are: `flye_out/assembly.fasta` (the contigs), `flye_out/assembly_graph.gfa` (the assembly graph in GFA format, viewable in Bandage), `flye_out/assembly_info.txt` (the per-contig summary with length, coverage, circularity), and `flye_out/flye.log` (the run log).

**Reproducibility caveat.** Flye is mostly deterministic given the same input, the same version, and the same thread count, but the multi-threaded overlap stage can produce slightly different contig sets on different thread counts because of order-of-merge effects on equivalent overlaps. Pin the thread count (`--threads 4`) in any script that needs byte-identical output.

### Canu 2.2

Canu is the older OLC assembler; the CLI is split into three stages but the wrapper `canu` runs them all.

| Flag | Purpose |
|------|---------|
| `-p <prefix>` | Output file prefix. |
| `-d <dir>` | Output directory. |
| `genomeSize=<size>` | Required estimated genome size, e.g. `1m`. |
| `-nanopore <reads.fastq>` | Nanopore input. |
| `-pacbio <reads.fastq>` | PacBio CLR input. |
| `-pacbio-hifi <reads.fastq>` | PacBio HiFi input. |
| `useGrid=false` | Run on the local machine (do not submit to SGE / Slurm). |
| `maxThreads=N` | Worker thread count. |
| `maxMemory=Ng` | Memory cap. |
| `corOutCoverage=N` | Target read coverage after correction (default 40). |

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

For a 50x simulated nanopore read set on a 1 Mb reference, Canu takes ~5-10 minutes on a 4-core laptop (significantly slower than Flye; this is a deliberate trade-off for the more-conservative graph representation). The output files are: `canu_out/ecoli.contigs.fasta` (the contigs), `canu_out/ecoli.contigs.gfa` (the assembly graph), `canu_out/ecoli.report` (the run report).

### Hifiasm 0.19.9

Hifiasm is specialized for PacBio HiFi input. The CLI is minimal.

| Flag | Purpose |
|------|---------|
| `-o <prefix>` | Output file prefix. |
| `-t N` | Worker thread count. |
| `--n-hap N` | Expected ploidy (1 for haploid, 2 for diploid). |
| `<reads.fastq>` | Input HiFi FASTQ (positional). |

```bash
hifiasm \
    -o hifiasm_out/ecoli \
    -t 4 \
    --n-hap 1 \
    reads.hifi.fastq
```

Hifiasm emits a GFA (`hifiasm_out/ecoli.bp.p_ctg.gfa` for the primary contigs); to convert to FASTA: `awk '/^S/ {print ">"$2"\n"$3}' hifiasm_out/ecoli.bp.p_ctg.gfa > hifiasm_out/ecoli.fasta`. For a 50x simulated HiFi read set on a 1 Mb reference this runs in ~30 seconds.

### Medaka 1.12.0

Medaka polishes ONT assemblies. It is a neural-network-based polisher whose models are trained per-basecaller-model.

| Flag | Purpose |
|------|---------|
| `-i <reads.fastq>` | Input nanopore reads (the same reads that produced the draft). |
| `-d <draft.fasta>` | Draft assembly to polish. |
| `-o <out_dir>` | Output directory. |
| `-m <model>` | Medaka model name. Must match the basecaller model. |
| `-t N` | Worker thread count. |
| `-b N` | Batch size (GPU only). |

#### The canonical Medaka call (Week 10 default)

```bash
medaka_consensus \
    -i reads.fastq \
    -d flye_out/assembly.fasta \
    -o medaka_out \
    -m r1041_e82_400bps_sup_v4.3.0 \
    -t 4
```

For a 50x nanopore read set polishing a 1 Mb Flye assembly this runs in ~2-5 minutes on a 4-core CPU; ~30 seconds on a consumer GPU. The output is `medaka_out/consensus.fasta` (the polished assembly). The Medaka model name encodes the basecaller chemistry (`r1041` = R10.4.1), the model variant (`sup` = super-accuracy), and the version (`v4.3.0`); always pin the exact string in the run-info JSON.

**Model selection.** The full list is at <https://github.com/nanoporetech/medaka/blob/master/medaka/options.py>. The default-for-R10.4.1-SUP recipe (as of 2024-2025) is `r1041_e82_400bps_sup_v4.3.0`. For R9.4.1 reads basecalled with Guppy SUP, the recipe was `r941_min_sup_g507`. Mismatching the model and the basecaller produces *worse* output than skipping the polish.

### BUSCO 5.7.1

BUSCO scores gene-content completeness against a curated lineage dataset.

| Flag | Purpose |
|------|---------|
| `-i <assembly.fasta>` | Input assembly. |
| `-l <lineage>` | Lineage dataset name, e.g. `bacteria_odb10`, `vertebrata_odb10`, `eukaryota_odb10`. |
| `-o <prefix>` | Output prefix. |
| `--out_path <dir>` | Output directory. |
| `-m <mode>` | Mode: `genome` (default), `transcriptome`, `proteins`. |
| `-c N` | Worker thread count. |
| `--offline` | Use locally-downloaded lineage data (skip the BUSCO download). |
| `--download_path <dir>` | Path to the downloaded BUSCO datasets. |

#### The canonical BUSCO call (Week 10 default)

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

For a 1 Mb bacterial assembly against `bacteria_odb10` (124 single-copy orthologs) this runs in ~30-60 seconds. The output of interest is `busco_out/busco_run/short_summary.specific.bacteria_odb10.busco_run.txt`, which reports the percentage of complete single-copy, complete duplicated, fragmented, and missing orthologs.

**Lineage selection.** The full lineage list is at <https://busco-data.ezlab.org/v5/data/lineages/>. Pick the smallest lineage that covers your organism; a smaller lineage runs faster and the percentages are easier to compare across runs.

### badread 0.4.1

badread simulates error-prone long reads from a reference FASTA.

| Flag | Purpose |
|------|---------|
| `simulate` | The subcommand. |
| `--reference <ref.fasta>` | Input reference. |
| `--quantity <coverage>` | Target read coverage, e.g. `50x`. |
| `--length <mean,sd>` | Read length distribution: comma-separated mean and standard deviation, e.g. `15000,13000`. |
| `--identity <mean,sd,max>` | Read identity distribution. For Dorado-style R10.4.1 SUP reads: `95,3,99`. |
| `--error_model <name>` | Error model: `nanopore2023`, `pacbio2021`. |
| `--qscore_model <name>` | Quality score model: `nanopore2023`, `pacbio2021`. |
| `--chimera <pct>` | Chimera fraction (default 1.0). |
| `--seed N` | Random seed. |

#### The canonical badread call (Week 10 default)

```bash
badread simulate \
    --reference reference.fasta \
    --quantity 50x \
    --length 15000,13000 \
    --identity 95,3,99 \
    --error_model nanopore2023 \
    --qscore_model nanopore2023 \
    --seed 42 \
    > reads.fastq
```

For a 1 Mb reference at 50x coverage this produces ~3,500 reads (mean 15 kb) in ~30 seconds. The reads have a realistic length distribution (some short, some up to 80+ kb) and a realistic per-read error profile (mean ~5% error rate, range ~1-15%).

### Dorado 0.7.2

Dorado is the Oxford Nanopore basecaller. It is not on Bioconda; download the binary from <https://github.com/nanoporetech/dorado/releases>. The CLI is short.

```bash
dorado basecaller \
    dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    pod5/ \
    --emit-fastq \
    > reads.fastq
```

Dorado requires a GPU for any non-toy run; on CPU it is ~100x slower. The Week 10 exercises and mini-project use `badread`-simulated reads instead of real POD5 + Dorado output; the Dorado call is documented here for completeness and for anyone working with real sequencer data.

### Bandage 0.9.0

Bandage is a GUI for viewing assembly graphs in GFA format. Download from <https://github.com/rrwick/Bandage/releases>. The CLI is:

```bash
Bandage load flye_out/assembly_graph.gfa
```

This opens the graph in a window. Right-click in the window to lay out the graph (`Graph drawing` -> `Random layout`); the graph topology becomes visible. For a clean bacterial assembly the graph is a single circle. For a problematic assembly the graph has bubbles (heterozygosity or sequencing errors) or tangles (unresolved repeats).

To export a screenshot from the CLI:

```bash
Bandage image flye_out/assembly_graph.gfa graph.png --height 800
```

This writes a PNG without opening the GUI; the CLI export is useful for headless servers and CI pipelines.

### Biopython 1.84

Biopython is the Python lab assistant. The relevant Week 10 modules:

- `Bio.SeqIO` — FASTA / FASTQ read and write. (Carried over from Weeks 2 and 5.)
- `Bio.SeqRecord` — the record object with `.id`, `.seq`, `.letter_annotations["phred_quality"]`.

Week 10 uses Biopython for FASTQ parsing (the simulated read set), FASTA parsing (the assembly contigs), and per-record length / quality summaries. The assembly itself is produced by external tools; Biopython is not a Python-native assembler.

```python
from Bio import SeqIO

records = list(SeqIO.parse("flye_out/assembly.fasta", "fasta"))
lengths = sorted((len(r.seq) for r in records), reverse=True)
total = sum(lengths)
# N50 calculation by walking the sorted lengths until cumulative >= total/2.
```

## File formats

### FASTQ (the input)

Same FASTQ format as Week 5, but with much longer reads (mean ~15 kb for ONT R10.4 vs ~150 bp for Illumina) and a different quality-score interpretation (the ONT base-quality scores from Dorado SUP are reasonably calibrated; the older Guppy SUP scores were not). Each record is four lines:

```
@read_001
ACGT...
+
qual...
```

The quality string is Phred+33; for long reads the per-base scores typically range from Q5 (mediocre) to Q30 (excellent), with the bulk in Q15-Q20 for Dorado SUP-basecalled R10.4 reads.

### FASTA (the assembly)

Standard FASTA. Flye and Hifiasm emit contigs with names like `contig_1`, `contig_2`, etc.; Canu uses `tig00000001`, `tig00000002`, etc. The assembly is conventionally one record per contig.

### GFA (the assembly graph)

**Graphical Fragment Assembly** format. The GFA 1.0 spec is at <https://gfa-spec.github.io/GFA-spec/GFA1.html>. Each line is one record:

- `S` (segment) — a sequence node with an ID and a base string.
- `L` (link) — an edge between two segments with orientation and an overlap CIGAR.
- `P` (path) — a named path through segments (the assembled contig).

For a clean bacterial assembly the GFA has one `S` line per contig, no `L` lines (no unresolved overlaps), and one `P` line that is just the segment. For a problematic assembly the GFA has multiple `S` and many `L` lines; Bandage renders the graph and the topology becomes visible.

### BUSCO summary

The `short_summary.specific.<lineage>.<run_name>.txt` file is the human-readable BUSCO output:

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

The `C` is the headline number: complete single-copy + duplicated as a percentage of the lineage gene set. A well-assembled bacterial genome typically hits C >= 99%; below 95% is a warning sign (either the assembly is fragmented or the wrong lineage was picked).

## Read-quality metrics

### N50, L50

The conventional contiguity statistics for an assembly:

```
Sort the contigs by length, descending.
Walk the sorted list and accumulate lengths.
N50 = the length of the contig at which the cumulative length first reaches >= total_length / 2.
L50 = the rank (1-indexed) of that contig.
```

For a 1 Mb single-circular-chromosome bacterial assembly with one contig of length 1,038,420 bp, N50 = 1,038,420 and L50 = 1. For a fragmented assembly with contigs of lengths 500k, 300k, 200k, 100k (total 1,100,000), the cumulative passes 550k at the 500k contig (still under) and at the 300k contig (cumulative 800k, over 550k), so N50 = 300,000 and L50 = 2.

A higher N50 means a more contiguous assembly; a lower L50 means the same.

### QV (Phred-scaled error rate)

```
QV = -10 * log10(error_rate)
```

A QV30 assembly has one error per 1,000 bp (10^-3); a QV40 has one per 10,000; a QV50 has one per 100,000. HiFi assemblies polished with DeepConsensus typically reach QV45-50; ONT assemblies polished with one round of Medaka typically reach QV35-40. Computing QV requires a ground-truth reference (the input to `badread`, for the simulated case) or a high-quality short-read alignment (for the real-data case).

## Compute requirements

| Step | Demo (1 Mb, 50x) | Real bacterial (5 Mb, 100x) |
|------|------------------|---------------------------------|
| badread simulate | ~30 sec | ~5 min |
| Flye `--nano-hq` | ~30-60 sec | ~5-10 min |
| Canu | ~5-10 min | ~30-60 min |
| Hifiasm (HiFi input) | ~30 sec | ~3 min |
| Medaka polish | ~2-5 min | ~10-20 min |
| BUSCO `bacteria_odb10` | ~30-60 sec | ~2-3 min |
| Bandage CLI export | ~5 sec | ~10 sec |

A laptop with 16 GB of RAM is enough for everything in Week 10 on the 1 Mb demo. The mini-project ramps to ~5 Mb, which is comfortable on 32 GB. Real human-genome assembly requires a workstation or cluster; that is out of scope.

## Datasets

All datasets used in Week 10 are free and publicly distributed, and most are *simulated* from a reference to keep the file sizes manageable.

### The 1 Mb demo reference (used by Exercise 1, 2, 3 and the mini-project)

A 1.04 Mb fragment of the *E. coli* K-12 MG1655 reference (NCBI accession `U00096.3`), or a synthetic 1 Mb FASTA built from random nucleotides plus inserted tandem repeats for didactic purposes. The file is bundled at `data/reference_1mb.fasta` (~1 MB). The reads are simulated by `badread` (the seed pins the output exactly):

```bash
badread simulate --reference data/reference_1mb.fasta --quantity 50x \
    --length 15000,13000 --identity 95,3,99 \
    --error_model nanopore2023 --qscore_model nanopore2023 \
    --seed 42 > data/reads_nanopore.fastq
```

The result is ~3,500 reads, ~50 MB FASTQ.

### The 5 Mb mini-project reference

A 5 Mb bacterial-sized synthetic genome with three repeat regions of varying difficulty (a tandem 10x repeat of a 200 bp unit; a dispersed 5x repeat of a 1,500 bp unit; a single 3 kb inverted repeat). File: `data/reference_5mb.fasta` (~5 MB). The reads are simulated similarly at 50x coverage. The presence of repeats is intentional — the mini-project write-up has to identify which repeats the assembler handled cleanly and which it collapsed.

### The HiFi panel (Challenge 1)

The same 1 Mb reference, simulated with badread's PacBio HiFi-style parameters: `--length 18000,5000 --identity 99.5,0.3,99.95 --error_model pacbio2021`. The HiFi reads have a much narrower error distribution; the assembly with Hifiasm is correspondingly cleaner.

## Free-tier compute and storage

- **GitHub Codespaces** (free 60 hours/month): a 2-core / 8 GB instance is enough for the 1 Mb demo. The mini-project's 5 Mb run is borderline; if the codespace runs out of memory during Flye, drop the `--threads` from 4 to 2.
- **Google Colab** (free tier): also enough for the 1 Mb demo. Flye and Medaka install cleanly via `apt install flye medaka` after enabling the bioconda channel. Colab's CPU is comparable to a 2018-era laptop; a 1 Mb run is comfortable.
- **Local laptop**: any 16 GB RAM laptop is enough for the demo and most of the mini-project. Real long-read assembly on a real bacterial or eukaryotic genome will want 32 GB and at least 8 cores.

## Style guide for this week

- **Pin tool versions.** Always write "Flye 2.9.5", "Canu 2.2", "Hifiasm 0.19.9", "Medaka 1.12.0", "BUSCO 5.7.1", "badread 0.4.1", "Dorado 0.7.2", "Bandage 0.9.0" — never just "Flye" or "the latest Medaka." Assemblies drift across versions; reproducibility requires pinning.
- **Pin the Flye input mode.** `--nano-hq` for modern Dorado-SUP-basecalled R10.4 reads; `--nano-raw` only for legacy Guppy-basecalled R9.4 reads; `--pacbio-hifi` for HiFi. Picking the wrong mode produces defensible-looking garbage.
- **Pin the Medaka model.** `r1041_e82_400bps_sup_v4.3.0` for the current Dorado SUP / R10.4.1 default. Mismatching the model and the basecaller is the most common "polishing makes things worse" failure mode in this domain.
- **Pin the BUSCO lineage dataset.** `bacteria_odb10` for bacteria; `vertebrata_odb10` for vertebrates; etc. Always record the BUSCO version *and* the dataset creation date (which BUSCO prints in the summary header).
- **Pin the seed.** badread takes `--seed`; Flye is mostly deterministic but the multi-thread merge can wobble — pin the thread count too; Medaka is deterministic given a fixed model. Bootstrap-style randomness is rare in this pipeline but the rule is the same: pin everything you can.
- **Cite tools by paper.** "Flye (Kolmogorov et al. 2019)" or "Hifiasm (Cheng et al. 2021)" — not just "Flye" or "Hifiasm." This is standard practice in publications and what reviewers expect.
- **Report numbers, not adjectives.** "Flye assembled the 50x simulated nanopore reads on the 1.04 Mb reference into 3 contigs totaling 1.038 Mb in 47 seconds; the N50 is 521,883 bp; after one round of Medaka polishing with the `r1041_e82_400bps_sup_v4.3.0` model, BUSCO reports 99.2% complete single-copy from `bacteria_odb10`" is a sentence. "The assembly was good" is not.
- **Always print the basecaller model, the assembler version, the polish model, the BUSCO lineage and dataset date, the seed, and the run date in the report header.** Without these you have an assembly that nobody can reproduce.
- **Always explain how repeats were handled.** A clean N50 number on a repeat-rich genome is a warning sign, not a brag; the assembler may have collapsed N copies into one. Bandage the graph and check.

## Common questions

**Q. Flye and Canu produce different assemblies on the same FASTQ. Which is right?**

Neither is uniquely "right" — they are two implementations of the same NP-hard problem (the simultaneous overlap-and-layout of N reads) with different heuristics. Disagreements almost always trace to repetitive regions where the two tools resolve the graph differently. For most downstream uses (gene prediction, comparative analysis), the two assemblies produce nearly identical results on the unambiguous regions. Treat persistent disagreement on a specific region (e.g. the rRNA operon in a bacterial genome) as a flag for manual inspection via Bandage.

**Q. My Flye assembly has 3 contigs but the reference is one circular chromosome. Why?**

Three possibilities. (1) Coverage gaps: a stretch of the reference is not covered by enough reads to span the overlap threshold; Flye breaks the assembly at that gap. (2) Unresolved repeats: an internal repeat region wider than the read N50 cannot be threaded by the OLC graph; Flye breaks the assembly at the repeat boundaries. (3) Adapter or chimeric reads inside the panel: a small fraction of reads contains junctions that fool the overlap stage. In every case the `flye.log` plus the GFA (viewed in Bandage) makes the cause visible.

**Q. Should I polish my Flye assembly with Medaka before or after BUSCO?**

After Flye, before BUSCO. The order matters: Medaka improves the per-base error rate, which improves gene-content scoring by BUSCO (a fragmented gene model from a few stray errors can score as `M` (missing) when it should be `F` (fragmented) or `S` (complete)). The mini-project pipeline runs `flye -> medaka -> busco` in that order.

**Q. Hifiasm on HiFi vs Flye on HiFi — which is better?**

Hifiasm on HiFi. Flye works on HiFi input (`--pacbio-hifi`) but Hifiasm is purpose-built for HiFi: it skips the read-correction stage, it emits haplotype-resolved contigs, and it runs ~5x faster. For diploid organisms (humans, plants) the haplotype resolution is decisive; for haploid bacteria the two assemblers produce very similar contigs but Hifiasm finishes faster. Default to Hifiasm on HiFi input; default to Flye on ONT input.

**Q. The Medaka model I want is not in the official list. What do I do?**

Medaka's model list is updated periodically. If your basecaller model is newer than your Medaka install, either (a) upgrade Medaka (`pip install -U medaka`), (b) use the closest model in the official list and accept a slightly worse polish, or (c) skip the polish and accept the raw Flye QV. Mismatching a model (e.g. using `r941_min_sup_g507` on R10.4 reads) is *worse than no polish*; do not do that.

**Q. Do I need a GPU?**

For Dorado basecalling, yes (CPU is ~100x slower than GPU). For everything else in Week 10, no — Flye, Canu, Hifiasm, BUSCO, and Bandage are all CPU tools. Medaka has a GPU code path that is ~10x faster than the CPU path but the CPU path works fine on the 1 Mb demo. The mini-project deliverables run end-to-end on a CPU-only laptop.

**Q. Do I need to commit the Flye intermediate files?**

No. The `flye_out/00-assembly/`, `flye_out/10-consensus/`, etc. directories contain intermediates that are large (often >> the final assembly) and reproducible from the inputs and the version. Commit `flye_out/assembly.fasta`, `flye_out/assembly_graph.gfa`, `flye_out/assembly_info.txt`, and `flye_out/flye.log`. Gitignore the rest.
