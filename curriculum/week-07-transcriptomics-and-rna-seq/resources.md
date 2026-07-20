# Week 7 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (fastp 0.23.4, kallisto 0.50.1, salmon 1.10.2, hisat2 2.2.1, subread 2.0.6, STAR 2.7.11a, pysam 0.22), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Conesa et al. (2016)** — "A survey of best practices for RNA-seq data analysis," the canonical RNA-seq tutorial review. The "Hitchhiker's Guide" that every modern RNA-seq paper cites. *Genome Biology* 17:13. Free full text:
  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8>
- **Bray, Pimentel, Melsted, Pachter (2016)** — kallisto's "near-optimal probabilistic RNA-seq quantification." The pseudoalignment paper. *Nature Biotechnology* 34:525. Free preprint:
  <https://arxiv.org/abs/1505.02710>
  Paper landing page (paywalled, but the preprint is the same content):
  <https://www.nature.com/articles/nbt.3519>
- **Patro, Duggal, Love, Irizarry, Kingsford (2017)** — Salmon: "fast and bias-aware quantification of transcript expression." The selective-alignment refinement of kallisto's idea. *Nature Methods* 14:417. Free full text:
  <https://www.nature.com/articles/nmeth.4197>
  Open preprint:
  <https://www.biorxiv.org/content/10.1101/021592v3>
- **Chen, Zhou, Chen, Gu (2018)** — `fastp`: a paper-length tool description for the adapter/quality trimmer used in this week's exercises. *Bioinformatics* 34:i884. Free full text:
  <https://academic.oup.com/bioinformatics/article/34/17/i884/5093234>
- **Kim, Paggi, Park, Bennett, Salzberg (2019)** — HISAT2: "graph-based alignment of next-generation sequencing reads to a population of genomes." *Nature Biotechnology* 37:907. Free full text:
  <https://www.nature.com/articles/s41587-019-0201-4>
- **Dobin et al. (2013)** — STAR: "ultrafast universal RNA-seq aligner." *Bioinformatics* 29:15. Free full text:
  <https://academic.oup.com/bioinformatics/article/29/1/15/272537>
- **Liao, Smyth, Shi (2014)** — `featureCounts`: a paper-length tool description. *Bioinformatics* 30:923. Free full text:
  <https://academic.oup.com/bioinformatics/article/30/7/923/232889>
- **Wagner, Kin, Lynch (2012)** — the original TPM paper. *Theory in Biosciences* 131:281. Free preprint:
  <https://link.springer.com/article/10.1007/s12064-012-0162-3>
- **Pachter (2011)** — "Models for transcript quantification from RNA-Seq." The cleanest mathematical derivation of TPM, FPKM, and the EM step. *arXiv*:1104.3889. Free preprint:
  <https://arxiv.org/abs/1104.3889>
- **Soneson, Love, Robinson (2015)** — `tximport` and the transcript-to-gene aggregation question. *F1000Research* 4:1521. Free full text:
  <https://f1000research.com/articles/4-1521/v2>
- **The GFF3 specification** — the canonical annotation format. ~15 pages, accessible:
  <http://gmod.org/wiki/GFF3>
- **The GENCODE project** — the human/mouse reference annotation. Browse and download cDNA + GTF:
  <https://www.gencodegenes.org/>
- **Ensembl FTP** — cross-species transcriptome FASTAs and GTF annotations. The mini-project's yeast files come from here:
  <http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/>
  <http://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/>

## Tool reference (the command-line surface)

### fastp 0.23.4

`fastp` is a single-binary adapter/quality trimmer plus QC reporter. It auto-detects Illumina adapters, trims low-quality 3' tails, removes polyG NovaSeq artifacts, and emits an HTML + JSON QC report in one pass. It is the modern replacement for the `Trimmomatic` + `FastQC` two-step.

| Flag | Purpose |
|------|---------|
| `-i R1.fq.gz -I R2.fq.gz` | Input paired-end FASTQs |
| `-o R1.trim.fq.gz -O R2.trim.fq.gz` | Output trimmed paired-end FASTQs |
| `--detect_adapter_for_pe` | Use overlap-based adapter detection (vs declared sequences) |
| `--qualified_quality_phred 20` | Per-base quality threshold (default 15; bump to 20 for stringent) |
| `--length_required 36` | Drop reads shorter than this after trimming |
| `--trim_poly_g` | Strip polyG tails (NovaSeq dark-base artifact) |
| `--cut_tail` | Sliding-window 3'-end quality trimming |
| `-h fastp.html -j fastp.json` | HTML report + JSON for downstream parsing |
| `-w 4` | Threads |

#### The canonical paired-end fastp call

```bash
fastp \
    -i raw/SRR12345_1.fq.gz \
    -I raw/SRR12345_2.fq.gz \
    -o trim/SRR12345_1.trim.fq.gz \
    -O trim/SRR12345_2.trim.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    -h qc/SRR12345.fastp.html \
    -j qc/SRR12345.fastp.json \
    -w 4
```

For a 3 M-paired-read yeast RNA-seq sample, this runs in ~30 seconds on a laptop and produces a ~600 KB HTML report.

### kallisto 0.50.1

`kallisto` is the original pseudoalignment quantifier. Its index is a colored de Bruijn graph of k-mers (k=31 default), one path per transcript. Its quantification step is EM over the per-read k-mer compatibility classes.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `kallisto index` | Build an index from a transcriptome FASTA | `-i index.idx` (output), `-k 31` (k-mer size, default 31) |
| `kallisto quant` | Quantify a paired-end sample against the index | `-i index.idx`, `-o out_dir/`, `-t 4` (threads), `-b 100` (bootstrap samples for uncertainty), `--single` (single-end mode), `-l 200 -s 30` (fragment-length mean+stddev for single-end) |
| `kallisto inspect` | Print index metadata (k-mer count, n_targets) | `index.idx` |
| `kallisto bus` | Single-cell mode (out of scope for Week 7) | — |

#### The canonical kallisto pipeline

```bash
# 1. Download transcriptome (yeast example).
curl -sLO http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# 2. Build index. ~10 sec for yeast, ~3 minutes for human.
kallisto index -i index/sce.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# 3. Quantify. ~30 seconds for 3M yeast pairs.
kallisto quant \
    -i index/sce.idx \
    -o quant/SRR12345/ \
    -t 4 \
    -b 100 \
    trim/SRR12345_1.trim.fq.gz \
    trim/SRR12345_2.trim.fq.gz
```

Output of `kallisto quant`:

- `abundance.tsv` — per-transcript counts and TPM (tab-separated, columns: `target_id`, `length`, `eff_length`, `est_counts`, `tpm`).
- `abundance.h5` — the same data plus the bootstrap matrix, in HDF5 format. Used by downstream tools like `sleuth` for uncertainty-aware DE.
- `run_info.json` — metadata: total reads, n_pseudoaligned, kallisto version, command line.

### salmon 1.10.2

`salmon` is the Patro et al. 2017 refinement of kallisto's pseudoalignment with selective alignment, richer bias models, and variational Bayesian EM. Its CLI is similar but not identical to kallisto's.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `salmon index` | Build an index from a transcriptome FASTA. With `--decoys`, also incorporate genomic decoy sequences for better specificity | `-t transcripts.fa`, `-i index/`, `-k 31`, `--decoys decoys.txt`, `--gencode` (for GENCODE-style FASTA headers) |
| `salmon quant` | Quantify a paired-end sample | `-i index/`, `-l A` (auto-detect library type), `-1 R1.fq.gz`, `-2 R2.fq.gz`, `-o out_dir/`, `-p 4` (threads), `--validateMappings` (enable selective alignment), `--gcBias`, `--seqBias`, `--posBias` (bias-model flags) |

#### The canonical Salmon pipeline

```bash
# 1. Build index (no decoys, simplest case).
salmon index \
    -t Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
    -i index/sce_salmon/ \
    -k 31

# 2. Quantify.
salmon quant \
    -i index/sce_salmon/ \
    -l A \
    -1 trim/SRR12345_1.trim.fq.gz \
    -2 trim/SRR12345_2.trim.fq.gz \
    -o quant_salmon/SRR12345/ \
    -p 4 \
    --validateMappings \
    --gcBias --seqBias
```

Output: `quant.sf` (per-transcript counts and TPM), `cmd_info.json`, `lib_format_counts.json`, and a `logs/` directory. The `quant.sf` columns are `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` — same content as kallisto's `abundance.tsv`, different column names.

### HISAT2 2.2.1

`HISAT2` is a splice-aware short-read aligner that uses a graph FM-index of the genome plus known splice sites. It is the lightweight alternative to STAR.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `hisat2-build` | Build a HISAT2 index | `<ref.fa> <index_prefix>` |
| `hisat2` | Align paired-end reads | `-x <index_prefix>`, `-1 R1.fq.gz -2 R2.fq.gz`, `-p 4`, `--dta` (downstream-transcript-assembly compatible), `-S out.sam` |

#### The canonical HISAT2 alignment

```bash
# 1. Build index. ~5 min for yeast, ~30 min for human.
hisat2-build ref/sce.fa index/sce_hisat2

# 2. Align. ~3 min for 3M yeast pairs.
hisat2 \
    -x index/sce_hisat2 \
    -1 trim/SRR12345_1.trim.fq.gz \
    -2 trim/SRR12345_2.trim.fq.gz \
    -p 4 \
    --dta \
| samtools sort -o aln/SRR12345.bam -

samtools index aln/SRR12345.bam
```

### STAR 2.7.11a

`STAR` is the Broad/cancer-genomics-standard splice-aware aligner. It is faster than HISAT2 but its index is large (~28 GB for human; ~1 GB for yeast).

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `STAR --runMode genomeGenerate` | Build a STAR index | `--genomeDir`, `--genomeFastaFiles`, `--sjdbGTFfile`, `--sjdbOverhang 99` (typically read_length - 1) |
| `STAR --runMode alignReads` | Align reads | `--genomeDir`, `--readFilesIn R1 R2`, `--readFilesCommand zcat` (for gzipped FASTQ), `--outSAMtype BAM SortedByCoordinate`, `--quantMode GeneCounts` (also produce a featureCounts-style table) |

#### The canonical STAR alignment

```bash
STAR \
    --runMode alignReads \
    --genomeDir index/sce_star/ \
    --readFilesIn trim/SRR12345_1.trim.fq.gz trim/SRR12345_2.trim.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix aln/SRR12345. \
    --runThreadN 4
```

### subread / featureCounts 2.0.6

`featureCounts` (part of the `subread` suite, Liao et al. 2014) is the standard alignment-based per-gene counter. It walks a sorted BAM, intersects each read with each annotated feature in a GTF/GFF3, and tallies reads.

| Flag | Purpose |
|------|---------|
| `-a annotation.gtf` | The annotation file |
| `-t exon` | Feature type to count (usually `exon` for gene-level counts) |
| `-g gene_id` | Group features by this attribute (sums all exons of a gene) |
| `-p --countReadPairs` | Count read pairs, not single reads (for paired-end RNA-seq) |
| `-o counts.tsv` | Output table |
| `-T 4` | Threads |
| `-s 2` | Strandedness (0=unstranded, 1=stranded, 2=reverse-stranded; standard for Illumina dUTP libraries) |

#### The canonical featureCounts call

```bash
featureCounts \
    -a ref/sce.gtf \
    -t exon \
    -g gene_id \
    -p --countReadPairs \
    -s 2 \
    -T 4 \
    -o counts/SRR12345.counts.tsv \
    aln/SRR12345.bam
```

The output is a 7-column TSV: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, `<sample_bam_name>`. The last column is the integer count for the sample. For a multi-sample run, pass multiple BAMs and you get one column per BAM.

### pysam 0.22 (BAM and FASTA access)

For Week 7 you mostly use pysam to read the HISAT2/STAR BAMs and the transcriptome FASTA. Most parsing of RNA-seq quantifier outputs is done with pandas, not pysam, because the quantifier outputs are plain TSV.

### pandas 2.2

Read kallisto's `abundance.tsv` and Salmon's `quant.sf` with `pd.read_csv(..., sep='\t')`. Aggregate per-transcript counts to per-gene counts with `df.groupby('gene_id')[['est_counts','tpm']].sum()`. Build the 3-sample mini-project counts matrix by joining three per-sample frames on `gene_id`.

## Compute requirements

| Step | Yeast (3 M pairs) | Human (30 M pairs) |
|------|-------------------|---------------------|
| `fastp` | ~30 sec | ~5 min |
| `kallisto index` (one-time) | ~10 sec | ~3 min |
| `kallisto quant` | ~30 sec | ~3 min |
| `salmon index` (one-time) | ~30 sec | ~10 min |
| `salmon quant --validateMappings` | ~60 sec | ~10 min |
| `STAR --runMode genomeGenerate` (one-time) | ~3 min, ~1 GB index | ~45 min, ~28 GB index |
| `STAR alignReads` | ~3 min | ~30 min |
| `hisat2-build` (one-time) | ~5 min | ~30 min |
| `hisat2` alignment | ~3 min | ~20 min |
| `featureCounts` | ~10 sec | ~2 min |

A laptop with 16 GB RAM is enough for everything except STAR on human (28 GB index, needs a workstation). The yeast mini-project fits comfortably on any laptop.

## Datasets

All datasets used in Week 7 are no-cost and publicly distributed.

### Yeast RNA-seq toy set (used by exercises 1-3 and the mini-project)

Three *Saccharomyces cerevisiae* paired-end RNA-seq samples from SRA, ~3 M paired reads each. These approximate a glucose-vs-galactose growth-condition comparison (the canonical yeast carbon-source experiment, expected to show strong differential expression at the `GAL1/GAL7/GAL10` regulon).

| SRA accession | Condition | Source |
|---------------|-----------|--------|
| `SRR453566` | Glucose, replicate 1 | Gierlinski et al. 2015, *Bioinformatics* 31:3625 |
| `SRR453567` | Glucose, replicate 2 | Gierlinski et al. 2015 |
| `SRR453568` | Galactose, replicate 1 | Gierlinski et al. 2015 |

Download with `prefetch SRR453566 && fasterq-dump SRR453566 --split-files` from the SRA Toolkit (free, NCBI). Or use the smaller subsetted FASTQs that the mini-project README points at — they are ~150 MB each rather than the ~1 GB full sample.

### Yeast transcriptome and annotation

```bash
curl -sLO http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
curl -sLO http://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
```

The cDNA FASTA is ~3.5 MB compressed; the GTF is ~600 KB compressed. Both contain ~7,000 transcripts of ~6,000 genes.

### Human transcriptome and annotation (for stretch goals)

```bash
# GENCODE 44 basic.
curl -sLO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
curl -sLO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
```

The transcripts FASTA is ~110 MB compressed; the GTF is ~50 MB compressed. ~250,000 transcripts of ~62,000 genes.

## Free-tier compute and storage

- **GitHub Codespaces** (free 60 hours/month for personal accounts): a 4-core / 16 GB instance is enough for the yeast mini-project end to end. Use the Codespaces-provided VS Code for the homework.
- **Google Colab** (free tier): also enough for the yeast mini-project. Install the conda tools in the first cell with `!conda install -c bioconda fastp=0.23 kallisto=0.50 -y`.
- **Local laptop**: any 16 GB RAM laptop is enough for the yeast pipeline. The exercises run end to end in ~5 minutes.

## Style guide for this week

- **Pin tool versions.** Always write "fastp 0.23.4", "kallisto 0.50.1", "salmon 1.10.2", "hisat2 2.2.1", "subread 2.0.6", "STAR 2.7.11a" — never just "fastp" or "kallisto". Other versions exist, behave subtly differently, and matter for reproducibility.
- **Pin annotations.** Always write "Ensembl release 110" or "GENCODE v44", not "the Ensembl annotation" or "the GENCODE GTF". The annotation drifts ~5% per release; reproducibility requires pinning.
- **Cite tools by paper.** "kallisto (Bray et al. 2016)" or "Salmon (Patro et al. 2017)" — not just "kallisto" or "Salmon". This is standard practice in publications and what reviewers expect.
- **Report counts and TPMs as numbers, not adjectives.** "GAL1 has 8,432 reads and TPM 1,247.3 in SRR453568" is a sentence. "GAL1 is highly expressed in galactose" is not — at minimum, give the number.
- **Use TPM for visualization, raw counts for differential expression.** The Week 8 DESeq2/edgeR pipeline takes raw integer counts; the heatmaps and PCAs you make this week use log-TPM. Mixing the two will produce wrong answers and is a common beginner mistake.

## Common questions

**Q. Why kallisto/Salmon and not just STAR + featureCounts?**

For most modern RNA-seq projects, pseudoalignment is ~30-100x faster than alignment + counting, recovers ~99% of the gene-level counts, and does not require disk for a BAM. Use STAR + featureCounts only when you specifically need a BAM (e.g. for splice-junction discovery, fusion calling, or variant calling on transcripts).

**Q. Should I use kallisto or Salmon?**

Both work. Salmon's `--validateMappings` mode (selective alignment) is slightly more accurate, especially when the sample has many sequence variants relative to the reference. Kallisto is simpler and slightly faster. For Week 7 we use kallisto in the main exercises and Salmon in Challenge 2 so you see both.

**Q. Why is kallisto's `est_counts` not always an integer?**

Multi-mapping reads are assigned fractionally by EM. A read that could come from `MYC`, `MYCL1`, or `MYCN` gets ~0.33 fractional count assigned to each, weighted by the current EM estimate of each transcript's abundance. The resulting `est_counts` are real numbers, not integers. Downstream tools like `tximport` round or scale appropriately.

**Q. Why is the human transcriptome only 300 MB when the genome is 3 GB?**

Because the transcriptome skips introns (which are ~95% of the genome by base count), intergenic regions, and the heterochromatic regions. The processed mRNA is ~300 MB; the unprocessed genomic DNA is 10x bigger.

**Q. What does `--validateMappings` actually validate?**

After the pseudoalignment step identifies a candidate set of transcripts each read could have come from, Salmon does a quick Smith-Waterman-like alignment score check between the read and each candidate. Transcripts where the actual alignment score is below threshold are removed from the read's compatibility class. This recovers ~1-2% accuracy on samples with high sequence-error or high-divergence-from-reference rates.

**Q. Do I need to commit the BAMs to git?**

No. BAMs are tens to hundreds of MB per sample. Gitignore them. Commit the `abundance.tsv`/`quant.sf` (~10 KB each) and the final counts matrix (~50 KB). These are all you need to reproduce the analysis from the trimmed FASTQs.
