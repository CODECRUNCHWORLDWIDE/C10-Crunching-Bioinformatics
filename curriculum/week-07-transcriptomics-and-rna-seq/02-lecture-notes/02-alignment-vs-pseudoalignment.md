# Lecture 2 — Alignment vs Pseudoalignment

> **Duration:** ~3 hours of reading + paper-and-pencil k-mer practice + a kallisto sanity-check run.
> **Outcome:** You can describe the k-mer compatibility class abstraction from Bray et al. 2016 in two paragraphs, run `kallisto index` + `kallisto quant` end to end on a yeast RNA-seq sample, run `salmon quant --validateMappings` on the same sample, and explain why both produce essentially the same gene-level counts as a STAR/HISAT2 + featureCounts pipeline at 1/30 the wall-clock time.

If you only remember one thing from this lecture, remember this:

> **Pseudoalignment never aligns. It looks up each read's k-mers in a precomputed index of the transcriptome, intersects the set of transcripts each k-mer is known to come from, and produces the read's "compatibility class" — the set of transcripts that could have generated the read. An EM algorithm then estimates per-transcript abundance from the empirical distribution of compatibility classes. The output is a counts table indistinguishable in quality from the alignment-based pipeline, produced in ~3 minutes per human sample instead of ~30 minutes. The trick is that you only need an alignment when you care about the exact base positions; for counting reads per transcript, you only need to know which transcripts the read came from. Bray et al. 2016 is the canonical reference.**

Lecture 1's trimmed FASTQs are the input. This lecture's `abundance.tsv` (kallisto) and `quant.sf` (Salmon) are the output. The transformation is the k-mer index plus the EM step in between.

---

## 1. Two architectures, one counts matrix

There are two production-quality ways to go from a trimmed RNA-seq FASTQ pair to a per-gene counts matrix:

**Architecture A — Classical alignment + counting.**

```
trimmed FASTQ ──> splice-aware aligner (STAR / HISAT2) ──> sorted BAM ──> featureCounts ──> counts.tsv
                  (~3-30 min per sample)                                  (~10 sec)
```

The aligner walks each read's bases against the **genome** index, finds the splice junctions if needed, and outputs a sorted BAM with the position-and-CIGAR of every read. `featureCounts` then intersects each read with each gene annotation in the GTF and tallies. The BAM is preserved and can be used for downstream applications: variant calling on RNA, splice-junction discovery, fusion detection, etc.

**Architecture B — Pseudoalignment + EM.**

```
trimmed FASTQ ──> kallisto / Salmon ──> abundance.tsv / quant.sf
                  (~30 sec - 3 min per sample)
```

No BAM, no alignment. The quantifier walks each read's k-mers against the **transcriptome** index, builds the compatibility class, and runs EM. The output is per-transcript estimated counts and TPM. Sub-second per million reads.

Both architectures produce a per-gene counts matrix. On a typical RNA-seq sample, the two matrices agree on **~99% of the top 5,000 expressed genes** (Patro et al. 2017 Figure 2, Soneson et al. 2015 Figure 3). The disagreements are concentrated in:

- **Paralog families** where multi-mapping is rampant (e.g. ribosomal protein genes, histone genes). The two methods resolve multi-mapping differently and disagree at the per-paralog level (though the family-level total is similar).
- **Lowly expressed genes** (< 10 reads). Stochastic differences between the two methods' algorithms become visible at low counts.
- **Genes with extensive intronic retention** (unspliced reads). Alignment captures them; pseudoalignment to a mature-mRNA transcriptome does not.

For 95% of RNA-seq projects (differential expression on well-expressed protein-coding genes), the two architectures are interchangeable. Pseudoalignment wins on wall-clock time and disk usage. The classical path wins when you need a BAM for downstream non-quantification work.

The rest of this lecture is **how** Architecture B works, why it produces the right answer, and how it compares mathematically to Architecture A. Section 7 then walks through Architecture A's end-to-end commands so you have both in muscle memory.

---

## 2. The k-mer compatibility class (Bray et al. 2016)

The central abstraction of kallisto is the **k-mer compatibility class** of a read.

Take a transcriptome of `T` transcripts, e.g. yeast Ensembl release 110 with `T = 6,975`. Pick a k-mer size, e.g. `k = 31` (kallisto's default). For each transcript, enumerate all `k`-mers (a transcript of length `L` has `L - k + 1` k-mers). For each unique k-mer `m`, record the **color set** `C(m) = { t : t is a transcript containing m }` — the set of transcripts that k-mer appears in.

```
transcript A: ATGCATGCATGCATGC...   k-mers: ATGCATGCATGCATGCAT...
transcript B: ATGCATGCAAGGAATT...   k-mers: ATGCATGCAAGGAATTAA...
transcript C: GTTAATGCATGCATGC...   k-mers: GTTAATGCATGCATGCAT...

k-mer 'ATGCATGCATGCATGCAT': C = {A, C}      (in both A and C)
k-mer 'TGCATGCAAGGAATT':   C = {B}           (only in B)
k-mer 'GTTAATGCATGCATG':   C = {C}           (only in C)
```

The full index is a **colored de Bruijn graph**: nodes are unique k-mers, edges connect k-mers that differ by one nucleotide and appear adjacent in some transcript, and each node carries the color set of which transcripts contain it. For yeast at k=31, the graph has ~6 M nodes and the entire index fits in ~50 MB on disk; for human GENCODE v44 basic, ~200 M nodes and ~3 GB on disk.

Now take a read of length `L_r = 100`. It has `L_r - k + 1 = 70` k-mers. Look up each in the index and intersect the color sets:

```
read k-mer 1: ATGCATGCATGCATGCAT   color set = {A, C}
read k-mer 2: TGCATGCATGCATGCATG   color set = {A, C}
read k-mer 3: GCATGCATGCATGCATGC   color set = {A, C}
...
read k-mer 70: ...                 color set = {A, C}

Compatibility class of read = intersection = {A, C}
```

The compatibility class is the set of transcripts the read could have come from. **The read has not been aligned.** The k-mer lookups are O(1) per k-mer (hash table) and the intersection is fast because color sets are small (most k-mers appear in 1-3 transcripts; even paralog-family k-mers appear in ~10 transcripts).

Total pseudoalignment cost per read: **`O(L_r)`** k-mer hash lookups. For a 100 bp read, ~70 hash lookups, each ~100 ns on a modern CPU. **~7 μs per read.** A sample of 3 M reads pseudoaligns in 3M × 7μs = **21 seconds**.

Compare to STAR/HISAT2, which does base-by-base alignment and splice-junction inference per read: ~100-1,000 μs per read, or ~5-50 minutes for the same sample. Pseudoalignment is **~100x faster** because it skips the alignment, and ~100x faster is the difference between "RNA-seq is a project" and "RNA-seq is a routine query."

### Why k=31?

Statistically, a k-mer of length 31 occurs by chance in a random 4 GB genome once every `4^31 ≈ 4.6 × 10^18` bases — i.e. essentially never. So any k-mer in a read that is in the index almost certainly came from the transcript(s) it points to. Smaller k (~21) gives more false matches and more ambiguous compatibility classes; larger k (~51) gives fewer false matches but is more sensitive to single-nucleotide errors (one base error in a 51-mer destroys 51 k-mer matches centered on that base).

k=31 is the sweet spot empirically (Bray et al. 2016 §Online Methods) and is the default in kallisto, Salmon, and most other k-mer-based bioinformatics tools (`mash`, `Mash Screen`, `KMC`).

---

## 3. The EM step (kallisto's quantification)

The compatibility class tells you the set of transcripts each read could have come from. It does **not** tell you which transcript the read actually came from when the compatibility class has more than one transcript. For that, you need to know the relative abundance of the transcripts — but the abundance is exactly what you are trying to estimate. The classical resolution is the **Expectation-Maximization (EM)** algorithm.

Let `θ_t` be the unknown abundance of transcript `t` (the fraction of all reads in the sample that came from `t`). Let `r` index reads. Let `C_r` be the compatibility class of read `r`. The likelihood of the observed reads given `θ` is:

```
L(θ) = ∏_r ∑_{t ∈ C_r} θ_t / l_t
```

where `l_t` is the effective length of transcript `t` (transcript length minus the mean fragment length — long transcripts produce more reads, so the abundance must be scaled by length to compare across transcripts).

The EM iteration:

- **E-step.** For each read `r` with compatibility class `C_r`, compute the expected fractional assignment of `r` to each transcript `t ∈ C_r` given the current `θ`:
  ```
  p(t | r, θ) = (θ_t / l_t) / ∑_{t' ∈ C_r} (θ_{t'} / l_{t'})
  ```
- **M-step.** Update `θ` as the normalized sum of expected fractional assignments:
  ```
  θ_t ← (∑_r p(t | r, θ)) / N
  ```
  where `N` is the total number of reads.

Iterate until `θ` changes by less than a small threshold. Typically converges in 50-200 iterations on a real sample.

The resulting `θ` is the maximum-likelihood per-transcript abundance vector. Multiplying by `N · l_t` gives the **estimated count** for transcript `t` — the number of reads kallisto thinks came from that transcript. The estimated counts are not integers because the fractional EM assignment is continuous; this is correct, and `tximport` or `DESeq2` handles the rounding/scaling correctly downstream.

Kallisto reports estimated counts in the `est_counts` column of `abundance.tsv` and **TPM** in the `tpm` column. We define TPM in Lecture 3.

### Bootstrap for uncertainty

`kallisto quant -b 100` runs 100 bootstrap replicates over the read set. Each bootstrap resamples the reads with replacement and reruns EM, producing 100 independent abundance estimates per transcript. The variance across bootstraps is a measure of quantification uncertainty: low for confidently quantified transcripts, high for transcripts with ambiguous multi-mapping. Downstream tools like `sleuth` (Pimentel et al. 2017, *Nature Methods* 14:687) use the bootstrap variance to improve differential-expression calls. For Week 7 we run `-b 100` for the mini-project but mostly read the point estimates; `sleuth` is a Week 8/9 topic.

---

## 4. Salmon: selective alignment (Patro et al. 2017)

`Salmon` (Patro et al. 2017, *Nature Methods* 14:417) is the second-generation pseudoaligner. It builds on kallisto's k-mer compatibility class idea and adds three refinements:

### 4.1 Selective alignment

Pure k-mer pseudoalignment can be fooled by reads with sequencing errors. A read with a single base error has its k-mers shifted: instead of 70 perfect k-mer matches to the true transcript, you might see 30 matches (the k-mers that span the error) and 40 spurious matches (the k-mers that happen to also exist in some unrelated transcript). The compatibility class becomes too inclusive, and the EM step underweights the true transcript.

Salmon's fix: after pseudoalignment identifies a candidate set of transcripts, run a fast **selective alignment** (a Smith-Waterman-like score check) between the read and each candidate. Drop candidates where the actual alignment score is below threshold. This recovers ~1-2% accuracy at the cost of ~2x runtime. Enable with `--validateMappings` (Salmon 1.x default since 1.0).

### 4.2 Fragment-bias models

RNA-seq libraries have well-known position biases that pseudoalignment by default ignores:

- **GC bias**: GC-rich fragments are over-represented (PCR favors them). Salmon's `--gcBias` flag learns a per-fragment GC correction factor.
- **Sequence bias**: hexamer priming during library prep is non-uniform. Salmon's `--seqBias` flag learns a per-position sequence-context correction.
- **Positional bias**: polyA selection biases reads toward the 3' end of transcripts. Salmon's `--posBias` flag learns a per-transcript-position correction.

Each bias model adds ~30 sec to the quantification step and improves accuracy by 0.5-2% depending on the sample.

### 4.3 Variational Bayesian EM

Salmon optionally replaces the standard EM with a **variational Bayesian** inference (`--useVBOpt`, default in Salmon 1.x) that is more robust on low-count transcripts and produces marginally better calibrated abundance estimates. The math is heavier; the user-visible effect is a slightly different `est_counts` value in low-count transcripts.

### 4.4 The canonical Salmon call

```bash
# Build index (no decoys, simple version).
salmon index \
    -t Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
    -i index/sce_salmon/ \
    -k 31

# Quantify with all bias corrections on.
salmon quant \
    -i index/sce_salmon/ \
    -l A \
    -1 trim/SRR453568_1.trim.fq.gz \
    -2 trim/SRR453568_2.trim.fq.gz \
    -o quant_salmon/SRR453568/ \
    -p 4 \
    --validateMappings \
    --gcBias --seqBias
```

Output: `quant.sf` (per-transcript counts and TPM), `cmd_info.json`, `lib_format_counts.json`, and a `logs/` directory. The `quant.sf` columns are `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` — same content as kallisto's `abundance.tsv`, different column names.

### Should I use kallisto or Salmon?

- **kallisto** is simpler, slightly faster, less RAM, and has a cleaner CLI. The pseudoalignment-only architecture is conceptually pure (Bray et al. 2016 is the paper that started the field).
- **Salmon** is more accurate on edge cases (samples with high variant load, samples with strong GC bias, low-count transcripts), and its `--decoy` option (using genomic sequence as decoy) improves specificity at low-expression genes by another ~1-2%.

In practice, both are widely used in production. Most published RNA-seq papers cite one or the other; some pipelines (like `nf-core/rnaseq`) run both and report the agreement. For Week 7 we use kallisto in Exercise 2 (because the CLI is simpler) and Salmon in Challenge 2 (so you see both).

---

## 5. Decoy-aware indexing

A subtle accuracy issue with both kallisto and Salmon is that the transcriptome FASTA only contains *known* transcripts. A read from an unannotated genomic region (an intergenic region, a novel transcript, a piece of contaminating DNA) may still find a k-mer match in the transcriptome by chance — and be assigned to some random transcript with non-zero probability.

The **decoy-aware index** (Patro et al. 2017 §Methods, Srivastava et al. 2020, *Genome Biology* 21:239) fixes this by adding the **genome** to the index alongside the transcriptome, marked as "decoy." Any read that prefers a decoy over a transcript is dropped from the EM. The result: transcripts with low abundance are no longer inflated by intergenic noise.

To build a decoy-aware Salmon index:

```bash
# 1. Build decoys.txt from chromosome names.
grep "^>" GRCh38.primary_assembly.genome.fa | sed 's/>//;s/ .*//' > decoys.txt

# 2. Concatenate transcriptome + genome.
cat gencode.v44.transcripts.fa GRCh38.primary_assembly.genome.fa > gentrome.fa

# 3. Build index with decoys.
salmon index \
    -t gentrome.fa \
    -d decoys.txt \
    -i index/human_salmon_decoy/ \
    -k 31 \
    --gencode
```

For the yeast mini-project, the small genome makes decoy-aware indexing not critical. For human/mouse, it is the standard since ~2020.

Kallisto does not natively support decoys until v0.50; if you need decoy-aware pseudoalignment, Salmon is the more mature choice.

---

## 6. Classical alignment: STAR and HISAT2

Both `STAR` (Dobin et al. 2013, *Bioinformatics* 29:15) and `HISAT2` (Kim et al. 2019, *Nature Biotechnology* 37:907) are splice-aware short-read aligners. The Week 7 mini-project uses HISAT2 because its index is smaller (~4 GB for human, ~50 MB for yeast) and its memory footprint fits on any laptop. STAR is faster on a server with adequate RAM and is the cancer-genomics standard.

### 6.1 HISAT2 index and alignment

```bash
# 1. Build index. ~5 sec for yeast.
hisat2-build ref/sce.fa index/sce_hisat2

# 2. Align. ~3 min for 3M yeast pairs.
hisat2 \
    -x index/sce_hisat2 \
    -1 trim/SRR453568_1.trim.fq.gz \
    -2 trim/SRR453568_2.trim.fq.gz \
    -p 4 \
    --dta \
| samtools sort -o aln/SRR453568.bam -

samtools index aln/SRR453568.bam
```

The `--dta` flag emits XS tags for downstream-transcript-assembly tools like StringTie. It is harmless if you do not use those tools.

### 6.2 STAR index and alignment

```bash
# 1. Build index. ~3 min for yeast with default sjdbOverhang.
STAR \
    --runMode genomeGenerate \
    --genomeDir index/sce_star/ \
    --genomeFastaFiles ref/sce.fa \
    --sjdbGTFfile ref/sce.gtf \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 10 \
    --runThreadN 4

# 2. Align. ~2 min for 3M yeast pairs.
STAR \
    --runMode alignReads \
    --genomeDir index/sce_star/ \
    --readFilesIn trim/SRR453568_1.trim.fq.gz trim/SRR453568_2.trim.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix aln/SRR453568. \
    --runThreadN 4
```

STAR's `--quantMode GeneCounts` mode produces a `ReadsPerGene.out.tab` file alongside the BAM — a per-gene counts table without needing a separate featureCounts step. The counts match `featureCounts -s 2` to within 1%.

### 6.3 featureCounts on the HISAT2/STAR BAM

If you used HISAT2 (no built-in counting) or want the featureCounts table for consistency:

```bash
featureCounts \
    -a ref/sce.gtf \
    -t exon \
    -g gene_id \
    -p --countReadPairs \
    -s 2 \
    -T 4 \
    -o counts/SRR453568.counts.tsv \
    aln/SRR453568.bam
```

Output: a 7-column TSV: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, `<sample_bam_name>`. The last column is the integer count for the sample.

### 6.4 Comparing alignment and pseudoalignment

On the yeast mini-project, you should see:

- **kallisto**: ~2.91 M of 3.07 M trimmed reads pseudoalign (~95%). ~3,800 transcripts with `est_counts > 10`.
- **Salmon**: ~2.93 M of 3.07 M reads assigned (~95%). Similar transcript-level abundance vector to kallisto, correlation > 0.99 on the top 3,000 transcripts.
- **HISAT2 + featureCounts -s 2**: ~2.85 M of 3.07 M reads assigned to a gene (~93%). ~3,750 genes with count > 10. Slightly lower assignment rate because alignment imposes a stricter base-by-base match requirement.

Gene-level correlation between kallisto-summed-to-gene and featureCounts: typically Pearson r > 0.99 on the top 3,000 genes, with the residual variance concentrated in paralog families and lowly expressed genes. The three pipelines agree.

---

## 7. Multi-mapping reads and how the three pipelines resolve them

The most consequential difference between the three pipelines is how they handle reads that map to multiple transcripts (or multiple genes).

**Multi-mapping in yeast**: ~3% of reads. The mappable transcriptome has few paralogs; multi-mapping is rare.

**Multi-mapping in human**: ~15-25% of reads. Many paralog families (ribosomal proteins, histones, olfactory receptors, immunoglobulins) generate reads that align to multiple paralogous transcripts.

### 7.1 STAR / HISAT2 + featureCounts (default)

`featureCounts` by default **discards multi-mappers** (reads with NH > 1 in the BAM). This is the conservative choice: every reported count is from a uniquely-mapping read. The downside: paralogous genes are systematically undercounted, and small gene families with high sequence similarity (e.g. some histone clusters) may have near-zero counts even when highly expressed.

Override with `-M --fraction`: assign multi-mappers fractionally (one read mapping to 4 transcripts contributes 0.25 to each). This is closer to what kallisto/Salmon do.

### 7.2 kallisto / Salmon EM

Both kallisto and Salmon resolve multi-mapping via the EM step (§3). A read with compatibility class `{A, B, C}` is fractionally assigned to A, B, C in proportion to their current abundance estimates. After convergence, the fractional counts reflect the best estimate of which transcript each read came from given the joint abundance vector.

This is mathematically the right thing to do — multi-mappers carry real information about the relative abundance of paralogs — but it produces non-integer counts that confuse some downstream tools. The `tximport` package (Soneson et al. 2015) provides the standard transcript-to-gene aggregation that handles fractional counts correctly.

### 7.3 The practical result

For a project that compares two conditions and asks "which genes change?", the EM approach (kallisto/Salmon) is more sensitive at paralog families. For a project that needs integer counts (e.g. a counting-based variant caller, or a downstream tool that errors on non-integers), the unique-mapper approach (`featureCounts` default) is the safer choice.

For Week 7's mini-project, we use kallisto's EM-resolved counts as the primary table.

---

## 8. A toy worked example

To make the pseudoalignment idea concrete, here is a fully worked toy example with three transcripts, two reads, k=4.

```
Transcriptome:
  t1 = ATGCATGCATGC
  t2 = ATGCATGCAAGG
  t3 = GTTAATGCATGC

k-mer index (k=4):
  ATGC -> {t1, t2, t3}  (in all three)
  TGCA -> {t1, t2, t3}
  GCAT -> {t1, t3}      (in t1, t3 only)
  CATG -> {t1, t3}
  TGCA (second occurrence) -> already covered
  CAAG -> {t2}          (only in t2)
  AAGG -> {t2}
  GTTA -> {t3}
  TTAA -> {t3}
  TAAT -> {t3}
  AATG -> {t3}

Read 1 = ATGCATGCAA (10 bp, k-mers: ATGC, TGCA, GCAT, CATG, ATGC, TGCA, GCAA)
  Compatibility class = intersection of color sets
                      = {t1,t2,t3} ∩ {t1,t2,t3} ∩ {t1,t3} ∩ {t1,t3} ∩ {t1,t2,t3} ∩ {t1,t2,t3} ∩ {}  (GCAA not in index)
                      = {t1, t3}  (treating GCAA as a sequencing error; in the real index this would push us to {t2} via partial matching, but in the toy we just intersect non-empty sets)
                      ≈ {t1, t3}

Read 2 = ATGCAAGG (k-mers: ATGC, TGCA, GCAA, CAAG, AAGG)
  Compatibility class = {t1,t2,t3} ∩ {t1,t2,t3} ∩ {} ∩ {t2} ∩ {t2}
                      = {t2}  (uniquely assigned)
```

After all reads have a compatibility class, the EM step estimates the abundances `θ_1, θ_2, θ_3`. If the dataset has 100 reads, of which:

- 50 have compatibility class `{t1}` (uniquely from t1).
- 30 have compatibility class `{t2}` (uniquely from t2).
- 10 have compatibility class `{t3}` (uniquely from t3).
- 10 have compatibility class `{t1, t3}` (could be from either).

Then after EM, the rough abundances are: `θ_1 ≈ (50 + 10 · α) / 100`, `θ_2 ≈ 30 / 100`, `θ_3 ≈ (10 + 10 · (1-α)) / 100`, where `α` is the fraction of the ambiguous reads EM assigns to t1. The EM solves for `α` self-consistently: given current `θ_1 / θ_3`, assign each ambiguous read in that ratio, then update `θ` to match. Converges in ~10 iterations on the toy.

Real samples have 3M reads and 7,000 transcripts; same algorithm, just scaled.

---

## 9. A Python sanity check: reading kallisto's `abundance.tsv`

After `kallisto quant -o quant/SRR453568/`, the output table is `quant/SRR453568/abundance.tsv`:

```
target_id            length   eff_length   est_counts   tpm
YAL001C              3483     3265.7       243.1        56.43
YAL002W              3825     3607.7       189.7        39.82
YAL003W              621      403.7        4521.8       8472.13
YAL005C              1929     1711.7       1212.4       536.41
...
```

Five columns:

- `target_id`: the transcript ID (matches the FASTA header).
- `length`: the actual transcript length in bp.
- `eff_length`: the effective length, computed as `length - mean_fragment_length + 1`. This is what abundance is normalized to.
- `est_counts`: the estimated count for this transcript (may be non-integer).
- `tpm`: transcripts per million.

A quick pandas summary:

```python
from __future__ import annotations
import pandas as pd
from pathlib import Path


def summarize_kallisto_quant(abundance_tsv: Path) -> dict:
    """Read kallisto abundance.tsv and return summary stats.

    Returns a dict with:
        n_transcripts, n_quantified (est_counts > 0),
        total_est_counts, top_transcript, top_tpm,
        mean_eff_length
    """
    df = pd.read_csv(abundance_tsv, sep="\t")

    n_transcripts: int = len(df)
    n_quantified: int = int((df["est_counts"] > 0).sum())
    total_est_counts: float = float(df["est_counts"].sum())
    top_row = df.sort_values("tpm", ascending=False).iloc[0]
    top_transcript: str = str(top_row["target_id"])
    top_tpm: float = float(top_row["tpm"])
    mean_eff_length: float = float(df["eff_length"].mean())

    return {
        "n_transcripts": n_transcripts,
        "n_quantified": n_quantified,
        "total_est_counts": total_est_counts,
        "top_transcript": top_transcript,
        "top_tpm": top_tpm,
        "mean_eff_length": mean_eff_length,
    }
```

A healthy yeast SRR453568 run returns something like:

```
n_transcripts:    6,975
n_quantified:     4,832          (~69%; lowly expressed transcripts have est_counts ≈ 0)
total_est_counts: 2,914,012      (close to but less than the trimmed read count; the difference is unmapped reads)
top_transcript:   YGL103W        (RPL28; ribosomal protein, always one of the top expressed)
top_tpm:          18,427.6
mean_eff_length:  1,256.4
```

Save these numbers. They are the numbers you will write into the mini-project methods section.

---

## 10. What to remember

- **Pseudoalignment is O(read_length) hash lookups + EM. Alignment is O(read_length × log(genome)) Smith-Waterman per read.** Pseudoalignment is ~100x faster because it skips the expensive part.
- **The compatibility class is the set of transcripts a read could have come from.** It is the only output of the "alignment-like" step. The EM step uses the empirical distribution of compatibility classes to estimate per-transcript abundance.
- **kallisto and Salmon agree to ~99% on well-expressed genes.** Differences appear at paralogs, low-count transcripts, and samples with high error rates. Use `--validateMappings` in Salmon if you have a noisy sample.
- **STAR/HISAT2 + featureCounts produce nearly the same gene-level counts.** Pseudoalignment wins on speed; alignment wins when you need a BAM.
- **Multi-mapping is handled differently.** `featureCounts` default drops multi-mappers; kallisto/Salmon assign them fractionally via EM. The fractional approach is more accurate at paralog families.
- **k=31 is the default for a reason.** Smaller k is too sensitive to chance matches; larger k is too sensitive to sequencing errors. Bray et al. 2016 §Methods empirically justifies the choice.
- **Decoy-aware indexing matters on human/mouse.** Add the genome as decoys in Salmon's index to suppress noise from unannotated genomic regions. Less critical on yeast/bacteria.
- **No BAM is needed for differential expression.** If you only need a counts matrix, pseudoalignment is the right tool. Reserve alignment for projects that need the per-read coordinate (variant calling on RNA, fusion detection, novel-isoform discovery).

Continue to **Lecture 3 — Counting and Normalization (TPM, CPM)** to turn the per-transcript counts into a per-gene counts matrix and to choose the right normalization for each downstream use.
