# Lecture 1 — From Reads to Transcripts

> **Duration:** ~3 hours of reading + a brief Python sanity check.
> **Outcome:** You can describe the RNA-seq experimental workflow from cell lysate to FASTQ in one paragraph, explain why a multi-exon mRNA cannot be aligned to the genome as a single contiguous block, run `fastp` end to end on a paired-end RNA-seq FASTQ pair, read the resulting HTML report, and decide whether your sample is good enough to quantify.

If you only remember one thing from this lecture, remember this:

> **RNA-seq is sampling, not coverage. A DNA-seq read at coverage 30x is one of 30 observations of a single static template (the diploid genome); an RNA-seq read is one of millions of observations from a distribution of transcripts whose abundances vary by gene, by isoform, and by condition over six orders of magnitude. Every downstream choice in Week 7 — splice-aware alignment, multi-mapping resolution by EM, library-size and transcript-length normalization — falls out of that one fact. The reference is no longer the chromosome the reads came from; the reference is the transcriptome the reads were drawn from.**

Week 6's BAM and VCF are parked. Week 7's input is a FASTQ pair, just like Week 5's. The output is no longer a per-position genotype call; it is a per-gene count.

---

## 1. The RNA-seq workflow in one page

A typical RNA-seq experiment, from cells to FASTQ, has the following steps:

1. **Cell or tissue lysis.** Disrupt the membrane in a chaotropic buffer (e.g. guanidinium thiocyanate in Trizol). RNA is released along with DNA and protein.
2. **Total RNA purification.** Phenol/chloroform extraction or silica-column extraction (Qiagen RNeasy, etc.). Output: ~5-20 μg of total RNA per million cells. Of that total, **~85% is ribosomal RNA** (rRNA, the structural RNA of the ribosome); ~10% is transfer RNA, small nuclear RNA, etc.; **~1-5% is mRNA**, the protein-coding messenger RNA that RNA-seq actually wants to measure.
3. **mRNA enrichment.** Two competing strategies:
   - **PolyA selection** (oligo-dT bead pulldown): mature mRNAs have a polyadenylated 3' tail; oligo-dT beads bind that tail and pull mRNA out of the total RNA pool. ~80-90% of reads end up on mRNA. **Caveat**: introduces a 3' coverage bias because the pulldown is anchored at the 3' end and fragments closer to that end survive better. Standard for differential expression on protein-coding genes; the default for the GTEx project and most commercial Illumina kits.
   - **Ribosomal-RNA depletion** (RiboZero, RiboMinus, Ribo-off): hybridization-based pull-down of rRNA, leaving everything else (mRNA + lncRNA + pre-mRNA + unprocessed mRNA). Better for non-polyadenylated transcripts (long noncoding RNAs, histone mRNAs), better for degraded/FFPE samples where polyA tails are gone, but slightly more rRNA carryover.
4. **Reverse transcription.** Convert the RNA to cDNA using a reverse transcriptase (e.g. SuperScript IV). Random hexamer or polyT priming. Output: double-stranded cDNA.
5. **Fragmentation.** Acoustic shearing (Covaris) or enzymatic fragmentation. Target insert size ~200-400 bp. Output: short cDNA fragments.
6. **Adapter ligation, PCR amplification, sequencing.** Same chemistry as DNA-seq. ~10-15 PCR cycles. Sequence ~75-150 bp from each end on an Illumina NovaSeq, NextSeq, or MiSeq.
7. **Demultiplexing.** Split the run by sample-specific barcode. Output: paired-end FASTQ files, one pair per sample.

Reference for the full workflow: **Conesa et al. 2016**, "A survey of best practices for RNA-seq data analysis," *Genome Biology* 17:13 — the open-access "Hitchhiker's Guide" that every modern RNA-seq paper cites. Read sections 1-2 this week; we will revisit sections 3-5 in Week 8.

Every step has a characteristic QC signature in the resulting FASTQ. PCR over-amplification shows up as high duplication (`fastp` reports it). Adapter readthrough at short inserts shows up as adapter content at the 3' end (`fastp` trims it). Failed rRNA depletion shows up as ~80% of reads pseudoaligning to the ~6 rRNA loci (`kallisto`'s `n_pseudoaligned` is high but per-non-rRNA-gene counts are low). 3' polyA-selection bias shows up as a steep drop in per-base coverage at the 5' end of long transcripts (a `qualimap rnaseq` plot makes it obvious; `fastp` does not). You will see all four signatures in the mini-project, with varying severity.

---

## 2. The mRNA and why splice-aware alignment is mandatory on vertebrate RNA

A typical human protein-coding gene is **~30 kb** long on the genome. Of that, **~95% is intron** — non-coding spacer sequence that is transcribed and then spliced out before the mRNA leaves the nucleus. The mature mRNA is **~2 kb** long and consists of ~8 exons joined end to end at splice junctions.

```
genome:    [exon1]──intron──[exon2]──intron──[exon3]──...──[exon8]      ~30 kb total
                                ↓
                  transcription + splicing
                                ↓
mRNA:      [exon1][exon2][exon3]...[exon8]                              ~2 kb total
```

An RNA-seq read of length 100 bp drawn from this mRNA might:

- Fall entirely within one exon (most common — exons average ~150 bp). Alignment to the genome is a contiguous match; ordinary `bwa mem` would work.
- Span an exon-exon junction (~10% of reads in a multi-exon gene). The read has, say, 60 bp from the end of exon 3 followed by 40 bp from the start of exon 4. On the genome, those 60 bp and 40 bp are separated by **a kilobase of intron** that is not in the read. An aligner that does not know about splicing will either fail to align this read, soft-clip it, or place it at a wrong intron-less paralog.

This is why `bwa mem` (the Week 5 short-read aligner) does not work for RNA-seq on a vertebrate genome. You need a **splice-aware** aligner that can recognize "the next 1,000 bp of the read are 1,000 bp away on the genome, separated by an intron." The two production-grade options for vertebrate RNA-seq are:

- **STAR** (Dobin et al. 2013, *Bioinformatics* 29:15) — uses a suffix-array index of the genome and a two-pass strategy: first pass finds anchor seeds, second pass uses the splice-junction database (from the GTF) to extend across junctions. Index is ~28 GB for human. Fast (~30 min per 30 M-read sample on a workstation), highly accurate, the cancer-genomics standard.
- **HISAT2** (Kim et al. 2019, *Nature Biotechnology* 37:907) — uses a graph FM-index of the genome plus known splice sites. Index is ~4 GB for human (much smaller than STAR), runs in ~20 min per sample, slightly less sensitive than STAR but light enough to run on a laptop. The "kallisto-of-classical-alignment" in terms of accessibility.

For yeast and most bacteria, splicing is rare (yeast has ~280 introns total across 6,000 genes; *E. coli* has zero). On those references you could in principle use `bwa mem`, but you should still use HISAT2 or STAR for consistency with the rest of the toolchain.

For the **transcriptome FASTA** (one record per transcript, splicing already resolved), splice-awareness is unnecessary by construction: each read maps to a transcript as a contiguous match. This is why pseudoalignment works on the transcriptome: the exon-exon junction problem disappears when you align to the transcriptome rather than the genome. We come back to this in Lecture 2.

---

## 3. The transcriptome as a reference

A **transcriptome FASTA** has one record per transcript. The header is typically a transcript ID (`ENST00000456328.2` for Ensembl, `YAL003W` for yeast) and the sequence is the mature mRNA (5' UTR + CDS + 3' UTR, exons concatenated, no introns).

For human GENCODE v44 basic:

- ~250,000 transcripts of ~62,000 genes (most genes have 2-5 isoforms; some pseudogenes have only one).
- Total FASTA size: ~300 MB uncompressed, ~110 MB gzipped.
- Median transcript length: ~1,400 bp. Longest transcript: ~110 kb (`TTN` — titin, the giant muscle protein gene).

For yeast Ensembl release 110:

- ~7,000 transcripts of ~6,000 genes.
- Total FASTA size: ~12 MB uncompressed, ~3.5 MB gzipped.
- Median transcript length: ~1,200 bp.

For each transcript, the natural reference for an RNA-seq read is the transcript itself, not the genomic locus the transcript came from. This is the philosophical inversion at the heart of pseudoalignment: **a read is an observation of a transcript, not of a genomic position.**

The annotation that ties transcripts back to genomic positions is the **GTF/GFF3** file. It has 9 tab-separated columns per row:

```
seqid   source   type        start    end       score  strand  phase  attributes
chr1    HAVANA   gene        11869    14409     .      +       .      gene_id=ENSG00000223972.5;...
chr1    HAVANA   transcript  11869    14409     .      +       .      gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;...
chr1    HAVANA   exon        11869    12227     .      +       .      gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;exon_number=1;...
chr1    HAVANA   exon        12613    12721     .      +       .      ...exon_number=2;...
chr1    HAVANA   exon        13221    14409     .      +       .      ...exon_number=3;...
```

Each `transcript` row is followed by one or more `exon` rows that share the same `transcript_id`. The `gene_id` attribute groups transcripts into genes. Most counting tools (`featureCounts`, `htseq-count`) consume the GTF and produce per-gene counts by summing reads over all exons of each gene.

GFF3 specification: <http://gmod.org/wiki/GFF3>. The format is plain text, tab-separated, and the schema fits on one page. You should be comfortable reading it column by column by the end of the week.

---

## 4. Read trimming with `fastp`

Before any alignment or pseudoalignment, raw FASTQ reads need:

- **Adapter trimming.** If the cDNA fragment was shorter than the sequencing read length (e.g. a 150 bp read on a 130 bp insert), the read will "run off the end" of the cDNA and into the adapter. The adapter sequence is at the 3' end of the read and must be removed before mapping.
- **Quality trimming.** Illumina quality drops at the 3' end of long reads (the polymerase is exhausted, the cluster is dephased). Bases below Q20 in the last 20 bp are usually trimmed.
- **PolyG removal.** On NovaSeq and NextSeq instruments (two-color chemistry), a missing signal is read as G. Reads that ran off the end of a short fragment have a polyG tail (often 20-50 G's). These G's are spurious and must be stripped.
- **Short-read filtering.** After trimming, very short reads (< 30-36 bp) align ambiguously to many transcripts. They are usually discarded.

The modern one-stop tool for all four is **`fastp`** (Chen et al. 2018, *Bioinformatics* 34:i884). It is a single binary that takes a FASTQ pair, auto-detects adapters by overlap analysis (no need to declare them), trims by quality and by polyG, drops short reads, and emits an HTML + JSON QC report alongside the trimmed FASTQs.

The canonical paired-end fastp call:

```bash
fastp \
    -i raw/SRR12345_1.fq.gz \
    -I raw/SRR12345_2.fq.gz \
    -o trim/SRR12345_1.trim.fq.gz \
    -O trim/SRR12345_2.trim.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --trim_poly_g \
    -h qc/SRR12345.fastp.html \
    -j qc/SRR12345.fastp.json \
    -w 4
```

What each flag does:

- `-i / -I`: paired-end input FASTQs.
- `-o / -O`: paired-end output FASTQs (trimmed and filtered).
- `--detect_adapter_for_pe`: use the paired-end-overlap method to detect adapters automatically. Robust against unusual library prep.
- `--qualified_quality_phred 20`: a base is "qualified" if its Q score is ≥ 20. Used by the per-window trimming algorithm.
- `--length_required 36`: discard reads shorter than 36 bp after trimming.
- `--trim_poly_g`: strip 3'-end polyG tails. Important for NovaSeq/NextSeq; harmless for HiSeq/MiSeq.
- `-h / -j`: emit HTML and JSON reports.
- `-w 4`: 4 worker threads.

For a 3 M-paired-read yeast sample, this runs in ~30 seconds on a laptop. For a 30 M-paired-read human sample, ~5 minutes.

### Reading the `fastp` HTML report

Open `qc/SRR12345.fastp.html` in a browser. The sections you care about, in order of diagnostic value:

1. **Summary** — total reads in, total reads out, % retained, % adapter trimmed, % low-quality trimmed. A healthy run keeps ≥ 90% of reads.
2. **Quality** — per-base quality boxplot. Should be Q30 across most of the read, with a slight drop at the 3' end. A drop to Q20 before position 50 indicates a degraded run.
3. **Adapter cutting** — how many reads had detectable adapter. ≥ 30% is normal for short-insert libraries; 0% means the auto-detection failed and you may need to declare adapters manually.
4. **Duplication** — estimated duplication rate. For RNA-seq, 30-60% is normal (highly expressed genes generate many copies of the same fragment). > 80% indicates PCR over-amplification or extreme library complexity loss.
5. **Insert size** — fragment-length distribution. Should be unimodal at ~200-400 bp. Bimodal indicates a library prep issue.
6. **K-mer / overrepresented sequences** — top overrepresented k-mers. If you see rRNA k-mers here, your rRNA depletion failed.

If the summary numbers look healthy (≥ 90% retained, Q30 median, < 60% duplication), proceed to pseudoalignment in Lecture 2. If any number is anomalous, debug *before* you spend an hour quantifying garbage.

---

## 5. PolyA selection, rRNA depletion, and the duplication question

A subtle but consequential decision in RNA-seq library prep is the mRNA enrichment strategy:

- **PolyA selection** (oligo-dT pulldown): pulls mature, polyadenylated mRNAs. Typical kit: Illumina TruSeq Stranded mRNA. Output: ~85-95% reads on mRNA, ~5-15% on rRNA carryover.
  - Pros: standard, well-validated, the default for most differential-expression projects.
  - Cons: misses non-polyadenylated transcripts (histone mRNAs, many long noncoding RNAs, bacterial mRNAs which have no polyA tail). Introduces a 3' coverage bias. Does not work on degraded RNA (RNA integrity number, RIN, < 7) where polyA tails are partially digested.
- **rRNA depletion** (RiboZero / RiboCop / Ribo-off): pulls rRNA out, leaves everything else. Typical kit: Illumina TruSeq Stranded Total RNA.
  - Pros: captures non-polyadenylated transcripts, works on degraded RNA, gives a more uniform per-position coverage profile.
  - Cons: more rRNA carryover (5-20%, vs 1-5% for polyA), more expensive kit, more reads "wasted" on non-mRNA non-coding RNAs.

If you are doing a vanilla "compare condition A to condition B on protein-coding genes" experiment, polyA selection is the conservative default. If you care about lncRNA, histone mRNA, FFPE samples, or anything bacterial, use rRNA depletion.

The **duplication question** comes up here because RNA-seq duplication can be either:

- **Optical / PCR duplication** (technical, bad): the same fragment sequenced twice because of optical bleed-through on the flowcell or because PCR over-amplified one molecule. These are *true duplicates* and should ideally be removed.
- **Biological duplication** (real signal, good): a highly expressed gene like `GAPDH` produces millions of copies of the same mRNA, which fragment into many copies of the same short fragment. The "duplicates" are independent observations of the same biological abundance.

DNA-seq distinguishes these two by removing duplicates (Week 5's `markdup` step). RNA-seq generally **does not** mark duplicates, because biological duplication is the signal and you cannot tell it from technical duplication based on coordinates alone (the fragmentation is not random enough). Modern best practice is to *not* mark RNA-seq duplicates and instead rely on the quantification step's per-transcript abundance estimate, which is robust to having many copies of the same fragment.

The exception: if `fastp`'s duplication estimate is > 80%, you have a complexity problem (PCR over-amplification or library exhaustion) and the sample is partially compromised regardless of which interpretation you favor.

---

## 6. Strandedness: a small but critical detail

Most modern Illumina RNA-seq libraries are **stranded**, meaning the orientation of the read tells you which strand of the genome the mRNA came from. This matters for genes on opposite strands that overlap (~5% of human genes have antisense partners). The standard Illumina dUTP protocol produces a **reverse-stranded** library: read 1 is reverse-complementary to the mRNA, read 2 matches the mRNA.

When you run `featureCounts`, pass `-s 2` for reverse-stranded libraries:

```bash
featureCounts -a annotation.gtf -t exon -g gene_id -p --countReadPairs -s 2 -o counts.tsv aln.bam
```

`-s 0` for unstranded, `-s 1` for stranded ("forward"), `-s 2` for reverse-stranded. The wrong setting silently drops or doubles your counts. Always check the library kit documentation. Salmon's `-l A` flag auto-detects strandedness; `kallisto` accepts `--fr-stranded` or `--rf-stranded` flags but does not auto-detect.

For yeast and bacteria, strandedness is less critical because overlapping genes are rarer. But you should still use the correct setting; it costs nothing and prevents subtle bugs.

---

## 7. A Python sanity check: parsing the `fastp` JSON

`fastp` emits a JSON report alongside the HTML. The JSON is machine-readable and the fields you care about for the mini-project are easy to extract. Example:

```python
from __future__ import annotations
import json
from pathlib import Path


def summarize_fastp_json(json_path: Path) -> dict:
    """Read a fastp JSON report and return key QC fields.

    Returns a dict with:
        reads_in, reads_out, pct_retained, pct_q30,
        adapter_pct, duplication_rate, insert_size_peak
    """
    with json_path.open() as f:
        data = json.load(f)

    summary = data["summary"]
    before = summary["before_filtering"]
    after = summary["after_filtering"]
    filt = data["filtering_result"]

    reads_in = before["total_reads"]
    reads_out = after["total_reads"]
    pct_retained = 100.0 * reads_out / max(reads_in, 1)
    pct_q30 = 100.0 * after["q30_bases"] / max(after["total_bases"], 1)

    adapter_cut = data.get("adapter_cutting", {})
    adapter_pct = 100.0 * adapter_cut.get("adapter_trimmed_reads", 0) / max(reads_in, 1)

    dup = data.get("duplication", {})
    duplication_rate = dup.get("rate", 0.0)

    insert = data.get("insert_size", {})
    insert_size_peak = insert.get("peak", 0)

    return {
        "reads_in": reads_in,
        "reads_out": reads_out,
        "pct_retained": pct_retained,
        "pct_q30": pct_q30,
        "adapter_pct": adapter_pct,
        "duplication_rate": duplication_rate,
        "insert_size_peak": insert_size_peak,
    }
```

A healthy yeast RNA-seq sample after `fastp` returns something like:

```
reads_in:           3,200,000
reads_out:          3,072,000
pct_retained:       96.0%
pct_q30:            93.5%
adapter_pct:        35.2%
duplication_rate:   0.42
insert_size_peak:   210
```

A degraded or compromised sample looks like:

```
reads_in:           3,200,000
reads_out:          2,100,000
pct_retained:       65.6%      <-- 30% drop, suspicious
pct_q30:            81.4%      <-- quality crashed
adapter_pct:        78.0%      <-- short inserts, library exhausted
duplication_rate:   0.89       <-- PCR over-amplification
insert_size_peak:   75         <-- much shorter than target ~300
```

Use these thresholds as a sanity check before pseudoalignment. If `pct_retained < 75%` or `pct_q30 < 85%` or `duplication_rate > 0.80`, flag the sample and revisit the library prep before running kallisto.

---

## 8. Worked example: trimming the yeast mini-project sample

We take SRR453568 (yeast, galactose growth, replicate 1) from the SRA. The full FASTQ is ~1 GB per mate; the exercise uses a subset of 800,000 paired reads (~150 MB per mate).

```bash
# 1. Download (or use the bundled subset).
prefetch SRR453568
fasterq-dump SRR453568 --split-files -O raw/

# 2. Trim with fastp.
fastp \
    -i raw/SRR453568_1.fastq \
    -I raw/SRR453568_2.fastq \
    -o trim/SRR453568_1.trim.fq.gz \
    -O trim/SRR453568_2.trim.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --trim_poly_g \
    -h qc/SRR453568.fastp.html \
    -j qc/SRR453568.fastp.json \
    -w 4
```

Expected runtime on a 4-core laptop: ~25 seconds. Expected output:

- `trim/SRR453568_1.trim.fq.gz`, `trim/SRR453568_2.trim.fq.gz`: trimmed FASTQs, ~130 MB each.
- `qc/SRR453568.fastp.html`: 600 KB HTML report.
- `qc/SRR453568.fastp.json`: 200 KB JSON report.

A healthy SRR453568 run produces: ~3.20 M reads in, ~3.07 M reads out (96% retained), median Q ~33, ~33% with adapter trimmed, duplication ~0.40, insert peak ~205 bp.

Save these numbers in your lab notebook. They are the numbers you will write into the mini-project methods section.

---

## 9. What to remember

- **RNA-seq is sampling, not coverage.** The depth at a gene is proportional to the gene's expression level, not to the gene's copy number. Six orders of magnitude in expression means six orders of magnitude in depth.
- **The reference is the transcriptome, not the genome.** A read came from a transcript. Map it to a transcript. This is what makes pseudoalignment work.
- **Splice-aware alignment is mandatory on vertebrate genomes.** STAR or HISAT2; never `bwa mem`. Or skip alignment entirely and use pseudoalignment against the transcriptome.
- **`fastp` is the modern Trimmomatic + FastQC.** One pass, one binary, HTML report, ~30 seconds on a yeast sample. Run it on every sample before quantification.
- **Strandedness matters.** Most modern Illumina libraries are reverse-stranded. Pass `-s 2` to featureCounts.
- **Most RNA-seq projects do not mark duplicates.** Biological duplication is the signal; you cannot tell it from technical duplication by coordinates. Trust the quantification step's abundance estimate instead.
- **PolyA selection is the default; rRNA depletion is the alternative.** Both are valid; the choice depends on whether you need non-polyadenylated transcripts.

Continue to **Lecture 2 — Alignment vs Pseudoalignment** to see how the transcriptome reference plus k-mer compatibility classes plus EM produces a counts matrix in 30 seconds instead of 30 minutes.
