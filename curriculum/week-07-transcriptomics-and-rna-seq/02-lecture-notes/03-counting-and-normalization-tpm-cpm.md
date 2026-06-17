# Lecture 3 — Counting and Normalization (TPM, CPM)

> **Duration:** ~2 hours of reading + paper-and-pencil arithmetic + a Python verification of the TPM identity.
> **Outcome:** You can write the formulas for CPM, RPKM, FPKM, and TPM from memory, explain when each one is the right normalization, aggregate per-transcript kallisto counts to per-gene counts using a transcript-to-gene map, build a 3-sample counts matrix, and compute log-TPM for use in a PCA or heatmap.

If you only remember one thing from this lecture, remember this:

> **Use raw integer counts going into a differential-expression tool (DESeq2, edgeR, limma-voom). Use TPM going into a heatmap, a PCA, or any human-eye comparison. CPM is what differential-expression tools compute internally as part of their normalization, but you generally do not work with CPM directly. RPKM/FPKM are deprecated relics whose only remaining value is that you may encounter them in older papers. The TPM/CPM distinction is the single most common beginner mistake in RNA-seq, and getting it right is the difference between a methods section that survives review and one that does not.**

Lecture 2's per-transcript counts (`abundance.tsv`, `quant.sf`, or featureCounts output) are the input. This lecture's per-gene, normalized counts matrix is the output. Week 8 takes that matrix and runs differential expression.

---

## 1. Why normalize at all?

Raw RNA-seq counts have two confounds that make cross-sample and cross-gene comparison nonsensical:

**Confound 1: library size.** Sample A was sequenced to 30 M reads; Sample B was sequenced to 50 M reads. A gene that has 100 reads in A and 100 reads in B is *not* equally expressed — A has many fewer total reads, so 100 of them is a much larger fraction. Without normalizing, "is this gene up or down?" cannot be answered.

**Confound 2: gene length.** A 10 kb gene expressed at a given level (say, 100 mRNA copies per cell) generates many more reads than a 1 kb gene expressed at the same level (also 100 copies per cell), because RNA-seq fragments the cDNA and each kilobase of cDNA generates one fragment. Without length-normalizing, "is gene X more expressed than gene Y?" cannot be answered.

The four canonical normalizations address these confounds with different conventions:

| Normalization | Corrects for library size? | Corrects for gene length? | Sum-to-million across samples? | Use case |
|---------------|:--------------------------:|:-------------------------:|:------------------------------:|----------|
| **Raw counts** | No  | No  | No  | Differential expression (DESeq2, edgeR, voom) |
| **CPM**        | Yes | No  | No (sums to library_size × 10^6 / library_size = 10^6) | Internal to DE tools; sometimes for visualization |
| **RPKM/FPKM**  | Yes | Yes | No  | Legacy (Mortazavi 2008); avoid in new work |
| **TPM**        | Yes | Yes | Yes (sums to 10^6 per sample, always) | Heatmaps, PCAs, cross-sample visualization, ratio comparisons |

The "sum-to-million-per-sample" property of TPM is what makes it the right normalization for visualization: when you take `log2(TPM + 1)` and plot a heatmap, the rows are directly comparable across samples because each column sums to the same constant. RPKM does *not* sum to a constant across samples; if Sample A has a few very highly expressed long genes and Sample B has the same expression distributed over more genes, the per-sample RPKM totals can differ by 30%, and a "Sample A has higher RPKM at gene X" statement is ambiguous.

References:
- **CPM**: introduced informally by edgeR; the cleanest derivation is in Robinson, McCarthy, Smyth 2010 (*Bioinformatics* 26:139), <https://doi.org/10.1093/bioinformatics/btp616>.
- **RPKM / FPKM**: Mortazavi et al. 2008, *Nature Methods* 5:621, "Mapping and quantifying mammalian transcriptomes by RNA-Seq."
- **TPM**: Wagner, Kin, Lynch 2012, *Theory in Biosciences* 131:281, "Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples."
- **The cleanest mathematical derivation of all four**: Pachter 2011, *arXiv*:1104.3889, "Models for transcript quantification from RNA-Seq."

---

## 2. CPM: counts per million

CPM corrects for library size only:

```
CPM_g = 10^6 × count_g / library_size

   where library_size = sum over all genes of count_g for this sample
```

Worked example. Sample SRR453568 has:

- 2,914,012 total assigned reads (the kallisto `n_pseudoaligned`).
- 8,432 reads on gene `GAL1`.

```
CPM(GAL1) = 10^6 × 8432 / 2914012 = 2,894.0
```

In words: GAL1 captured about 0.29% of all reads in this sample (2,894 of every million).

Multi-sample example. Three samples with library sizes 2.9 M, 3.1 M, 2.7 M, and per-sample gene counts:

| Gene  | counts in S1 | counts in S2 | counts in S3 | CPM in S1 | CPM in S2 | CPM in S3 |
|-------|-------------:|-------------:|-------------:|----------:|----------:|----------:|
| GAL1  | 200          | 220          | 8,432        | 69.0      | 71.0      | 3,123.0   |
| ACT1  | 5,000        | 5,100        | 4,800        | 1,724.1   | 1,645.2   | 1,777.8   |
| RPL28 | 30,000       | 31,000       | 28,000       | 10,344.8  | 10,000.0  | 10,370.4  |

CPM normalizes for library size: ACT1 is similar across samples in CPM (~1,700), even though the raw counts differ. RPL28 is similar across samples in CPM (~10,000), as expected for a housekeeping gene. GAL1 is much higher in S3 (galactose-induced) — *not* a library-size artifact but real biology.

CPM does **not** correct for gene length, so comparing CPM across genes is misleading. CPM(GAL1) = 3,123 in S3 does not mean GAL1 has 3 mRNA copies per 1,000 mRNAs in the cell — it means GAL1 fragments make up 3,123 per million of the read library, which depends on both the gene's mRNA abundance and the gene's length.

When you would use CPM in a script:

- Inside `edgeR::cpm()` for filtering low-expression genes before a DE test: `keep <- rowSums(cpm(counts) > 1) >= 3` keeps genes with CPM > 1 in at least 3 samples.
- For a quick "what fraction of reads is on chrM?" sanity check.

For visualization or cross-gene comparison, prefer TPM.

---

## 3. RPKM and FPKM: deprecated, but you will encounter them

RPKM (Reads Per Kilobase per Million reads) and FPKM (Fragments Per Kilobase per Million fragments, the paired-end variant) were introduced in Mortazavi et al. 2008. The formula:

```
RPKM_g = 10^9 × count_g / (library_size × gene_length_bp)
       = (CPM_g × 1000) / gene_length_bp
```

The "kilobase" comes from dividing by `gene_length_bp / 1000`, hence the `10^9` numerator instead of `10^6`. FPKM is identical mathematically; the rename to "Fragments" reflects that paired-end reads are counted as one fragment (one pair) rather than two reads.

**Why RPKM/FPKM is deprecated.** The Wagner et al. 2012 paper showed that RPKM is **not invariant** across samples: the same gene at the same biological abundance in two samples can have very different RPKM values if the two samples differ in their highly expressed long-gene contribution. The fix is to normalize the per-gene "RPK" (reads per kilobase) by the **per-sample sum of RPKs** rather than by library size. That fix produces TPM.

You will encounter RPKM/FPKM in:

- Older papers (pre-2015) that report RPKM tables.
- The output of `cufflinks` (the 2012-2018 standard isoform quantifier, now superseded by Salmon and StringTie).
- The legacy `--quantMode GeneCounts` output of older STAR versions.

In new work, do not report RPKM/FPKM. Report TPM and raw counts. If you need to compare to a published RPKM table, the conversion is:

```
TPM_g = (RPK_g / sum_g'(RPK_g')) × 10^6
       where RPK_g = count_g / (gene_length_bp / 1000) = RPKM_g × library_size / 10^6
```

---

## 4. TPM: transcripts per million

TPM is the modern standard for cross-sample, cross-gene RNA-seq visualization.

The Wagner et al. 2012 formula:

```
TPM_g = 10^6 × (count_g / eff_length_g) / sum_g'(count_g' / eff_length_g')
```

The derivation in two steps:

- **Step 1**: per-gene rate `r_g = count_g / eff_length_g`. This is the per-base read rate at gene `g`, accounting for the fact that longer genes generate more reads per molecule.
- **Step 2**: normalize so the per-sample sum is 10^6:
  ```
  TPM_g = 10^6 × r_g / sum_g'(r_g')
  ```

The **effective length** `eff_length_g` is the gene length minus the mean fragment length (plus 1). This is because a fragment of length `f` can start anywhere in a transcript of length `L` from position 1 to position `L - f + 1`, so the number of fragment-start positions is `L - f + 1`, and that is what scales the count. For short transcripts where `L < f`, the effective length is conventionally floored at a small positive number to avoid division by zero. kallisto and Salmon report `eff_length` directly in their output tables.

For long genes where `L >> f`, `eff_length ≈ L - f + 1 ≈ L`. For typical Illumina libraries with `f ≈ 200`, the correction matters only for genes < ~500 bp.

### Worked example: TPM by hand

Sample with three genes, library size of 1,000,000 reads:

| Gene  | count   | length_bp | eff_length (~length - 200) | count / eff_length |
|-------|--------:|----------:|---------------------------:|-------------------:|
| A     | 1,000   | 5,000     | 4,800                      | 0.2083             |
| B     | 1,000   | 1,000     | 800                        | 1.2500             |
| C     | 1,000   | 500       | 300                        | 3.3333             |

Sum of rates: `0.2083 + 1.2500 + 3.3333 = 4.7917`.

```
TPM(A) = 10^6 × 0.2083 / 4.7917 = 43,475
TPM(B) = 10^6 × 1.2500 / 4.7917 = 260,870
TPM(C) = 10^6 × 3.3333 / 4.7917 = 695,655
```

Check: sum = 43,475 + 260,870 + 695,655 = 1,000,000. **TPM always sums to 10^6 per sample, by construction.**

This is the property that makes TPM the right normalization for visualization. If you have three samples and you plot `log2(TPM + 1)` for each gene across samples, the rows are directly comparable: the row mean is the "average expression" of the gene, the row variance is the "variability" of the gene, and the column sums are all 10^6 so the columns are calibrated to each other.

### The TPM identity in code

The TPM identity is `sum_g(TPM_g) = 10^6`. Always. The implementation:

```python
from __future__ import annotations
import pandas as pd
import numpy as np


def counts_to_tpm(counts: pd.Series, eff_length: pd.Series) -> pd.Series:
    """Convert per-gene counts to TPM given per-gene effective lengths.

    counts:      raw counts per gene
    eff_length:  effective length per gene (kallisto eff_length column or
                 gene_length_bp - mean_fragment_length + 1)

    Returns: TPM Series; sum will be 1e6 (modulo float precision).
    """
    rate: pd.Series = counts / eff_length
    tpm: pd.Series = 1e6 * rate / rate.sum()
    return tpm
```

Test:

```python
>>> counts = pd.Series([1000, 1000, 1000], index=["A", "B", "C"])
>>> eff_length = pd.Series([4800, 800, 300], index=["A", "B", "C"])
>>> tpm = counts_to_tpm(counts, eff_length)
>>> tpm
A     43475.131
B    260870.787
C    695654.082
dtype: float64
>>> tpm.sum()
999999.99999...   # 1e6 modulo float precision
```

Verify against kallisto's own `tpm` column from `abundance.tsv`. The numbers should match to within 0.1% (floating-point rounding); if they do not, your `eff_length` source is wrong.

---

## 5. Transcript-to-gene aggregation

Both kallisto and Salmon produce **per-transcript** counts. Most downstream analyses (DE, GO enrichment, pathway analysis) work at the **gene** level. The aggregation step is straightforward: for each gene, sum the counts across its transcripts.

The mapping from transcript ID to gene ID lives in the GTF/GFF3 file. A typical row:

```
transcript YAL001C; gene YAL001C  (yeast: most genes have one transcript)
transcript ENST00000456328.2; gene ENSG00000223972.5  (human: gene has 8 transcripts)
```

For yeast, the transcript-to-gene mapping is almost trivial (one transcript per gene; the transcript ID and gene ID are often identical or differ only by suffix). For human, the mapping is many-to-one and the aggregation matters.

### Naive sum-by-gene

```python
from __future__ import annotations
import pandas as pd
from pathlib import Path


def parse_t2g_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Read a GTF and extract the transcript_id -> gene_id mapping.

    Returns a DataFrame with columns ['transcript_id', 'gene_id'].
    """
    rows: list[dict[str, str]] = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            t_id: str = ""
            g_id: str = ""
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id"):
                    t_id = kv.split('"')[1] if '"' in kv else kv.split()[1]
                elif kv.startswith("gene_id"):
                    g_id = kv.split('"')[1] if '"' in kv else kv.split()[1]
            if t_id and g_id:
                rows.append({"transcript_id": t_id, "gene_id": g_id})
    return pd.DataFrame(rows)


def aggregate_transcript_to_gene(
    abundance: pd.DataFrame, t2g: pd.DataFrame
) -> pd.DataFrame:
    """Sum per-transcript counts and TPM to per-gene.

    abundance:  pandas frame from kallisto abundance.tsv with columns
                target_id, eff_length, est_counts, tpm
    t2g:        DataFrame with columns transcript_id, gene_id

    Returns: per-gene DataFrame with gene_id, est_counts, tpm.
    """
    df = abundance.merge(t2g, left_on="target_id", right_on="transcript_id")
    by_gene = df.groupby("gene_id", as_index=False)[["est_counts", "tpm"]].sum()
    return by_gene
```

The naive sum-by-gene is the right thing to do for `tpm` (TPM is additive across transcripts of the same gene because the per-million normalization is at the sample level) and **almost** right for `est_counts`. The Soneson et al. 2015 paper notes that the better aggregation for `est_counts` is "lengthScaledTPM": convert TPM to a count-scale by multiplying by `(library_size / 1e6)`, which correctly handles the fact that different transcripts of the same gene have different effective lengths. The `tximport` R package implements this.

For yeast (one transcript per gene), the distinction does not matter; naive sum and lengthScaledTPM give identical results. For human (multi-isoform genes), lengthScaledTPM is the recommended aggregation.

---

## 6. featureCounts: alignment-based counting

The classical alignment path (HISAT2 / STAR + featureCounts) bypasses transcript-level quantification and counts directly at the gene level. The `featureCounts` tool walks the sorted BAM, intersects each read (or read pair) with the annotated features in the GTF, and tallies.

The canonical featureCounts call:

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

Flag-by-flag:

- `-a annotation.gtf`: the annotation file.
- `-t exon`: count reads that overlap features of this type. For gene-level counts, use `exon` (each exon belongs to a gene; reads on any exon of a gene are tallied to that gene).
- `-g gene_id`: group features by this attribute. With `-t exon -g gene_id`, all exons of the same gene contribute to a single per-gene count.
- `-p --countReadPairs`: paired-end mode; count read pairs as single observations (not double).
- `-s 2`: reverse-stranded library (Illumina dUTP). `-s 0` for unstranded, `-s 1` for forward-stranded.
- `-T 4`: threads.
- `-o counts/...`: output file.

Output format. The TSV has a header comment line, then a column header, then one row per gene:

```
# Program:featureCounts v2.0.6; Command:featureCounts ...
Geneid    Chr     Start    End      Strand   Length   aln/SRR453568.bam
YAL001C   chrI    151169   147594   -        3477     243
YAL002W   chrI    143707   147531   +        3825     189
YAL003W   chrI    142255   142875   +        621      4521
YAL005C   chrI    139503   141431   +        1929     1212
...
```

The first six columns are gene annotation; the seventh column is the count for this sample (the BAM file path becomes the column name). For multi-sample runs, pass multiple BAMs and the output has one count column per BAM.

featureCounts also produces a `<output>.summary` file listing the read-assignment breakdown:

```
Status                                  aln/SRR453568.bam
Assigned                                2854312
Unassigned_Unmapped                     46221
Unassigned_Read_Type                    0
Unassigned_Singleton                    12041
Unassigned_MappingQuality               0
Unassigned_Chimera                      8421
Unassigned_FragmentLength               0
Unassigned_Duplicate                    0
Unassigned_MultiMapping                 102434
Unassigned_Secondary                    0
Unassigned_NonSplit                     0
Unassigned_NoFeatures                   58219
Unassigned_Overlapping_Length           0
Unassigned_Ambiguity                    14982
```

The `Assigned` count is your effective library size for downstream normalization. The `Unassigned_*` rows are useful diagnostics: high `Unassigned_NoFeatures` indicates many reads on intergenic regions (failed rRNA depletion?); high `Unassigned_MultiMapping` indicates a paralog-heavy sample (or a wrong `-s` strandedness flag).

### Multi-sample featureCounts

For a 3-sample mini-project, run featureCounts once with all three BAMs:

```bash
featureCounts \
    -a ref/sce.gtf \
    -t exon \
    -g gene_id \
    -p --countReadPairs \
    -s 2 \
    -T 4 \
    -o counts/all_samples.counts.tsv \
    aln/SRR453566.bam aln/SRR453567.bam aln/SRR453568.bam
```

Output: a TSV with columns `Geneid, Chr, Start, End, Strand, Length, aln/SRR453566.bam, aln/SRR453567.bam, aln/SRR453568.bam`. Strip the first six columns and rename the BAM-path columns to sample names; you have your 3-sample counts matrix:

```python
from __future__ import annotations
import pandas as pd
from pathlib import Path


def featurecounts_to_matrix(fc_tsv: Path, sample_map: dict[str, str]) -> pd.DataFrame:
    """Read featureCounts output and return a per-gene counts matrix.

    sample_map: dict from BAM-path column name to sample name, e.g.
                {"aln/SRR453566.bam": "S1_glucose_rep1", ...}

    Returns: DataFrame indexed by gene_id, columns = sample names.
    """
    df = pd.read_csv(fc_tsv, sep="\t", comment="#")
    df = df.rename(columns=sample_map)
    return df.set_index("Geneid")[list(sample_map.values())]
```

This is the counts matrix you feed into Week 8's `DESeq2` / `edgeR`.

---

## 7. The kallisto-to-counts-matrix pipeline

For a 3-sample mini-project using kallisto, the equivalent pipeline:

```python
from __future__ import annotations
import pandas as pd
from pathlib import Path


def build_counts_matrix_from_kallisto(
    quant_dirs: dict[str, Path],
    t2g: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate per-sample kallisto abundance.tsv to a per-gene counts matrix.

    quant_dirs: dict from sample name to kallisto quant output directory.
    t2g:        transcript-to-gene mapping (output of parse_t2g_from_gtf).

    Returns: DataFrame indexed by gene_id, columns = sample names.
    """
    per_sample: list[pd.Series] = []
    for sample, quant_dir in quant_dirs.items():
        abundance = pd.read_csv(quant_dir / "abundance.tsv", sep="\t")
        df = abundance.merge(t2g, left_on="target_id", right_on="transcript_id")
        gene_counts = df.groupby("gene_id")["est_counts"].sum()
        gene_counts.name = sample
        per_sample.append(gene_counts)

    matrix = pd.concat(per_sample, axis=1).fillna(0.0)
    return matrix
```

A 3-sample yeast matrix has ~6,000 rows and 3 columns. Typical size: ~150 KB as TSV. Commit it.

---

## 8. log-TPM for visualization

For heatmaps, PCAs, and any human-eye comparison, the standard transformation is `log2(TPM + pseudocount)` where the pseudocount is conventionally 1.

```python
from __future__ import annotations
import numpy as np
import pandas as pd


def log_tpm_matrix(counts_matrix: pd.DataFrame, eff_lengths: pd.Series) -> pd.DataFrame:
    """Convert per-gene counts matrix to log2(TPM + 1).

    counts_matrix: rows = genes, columns = samples, integer or float counts.
    eff_lengths:   per-gene effective length, indexed by gene_id.

    Returns: log2(TPM + 1) matrix, same shape as counts_matrix.
    """
    rate_matrix = counts_matrix.div(eff_lengths, axis=0)
    tpm_matrix = rate_matrix * (1e6 / rate_matrix.sum(axis=0))
    return np.log2(tpm_matrix + 1.0)
```

The `+1` pseudocount avoids `-inf` at zero-count genes. The `log2` compresses the dynamic range (TPM spans 0 to ~50,000, log2(TPM) spans 0 to ~16) so that heatmaps and PCAs are not dominated by a handful of very highly expressed genes.

For PCA on a 3-sample yeast log-TPM matrix:

```python
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def pca_on_log_tpm(log_tpm: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """Run PCA on a log-TPM matrix, return per-sample PC scores.

    log_tpm: rows = genes, columns = samples.
    Returns: DataFrame with rows = samples, columns = PC1, PC2, ...
    """
    pca = PCA(n_components=n_components)
    # PCA expects rows = observations; here observations are samples.
    scores = pca.fit_transform(log_tpm.T.values)
    return pd.DataFrame(
        scores,
        index=log_tpm.columns,
        columns=[f"PC{i+1}" for i in range(n_components)],
    )
```

A healthy 3-sample yeast PCA shows the galactose sample separated from the two glucose samples on PC1 (which captures the carbon-source-induced expression program), with the two glucose replicates clustering tightly. The mini-project asks you to produce this plot and write up its interpretation.

---

## 9. Common normalization pitfalls

**Pitfall 1: TPM on the wrong axis.**
TPM normalizes within a sample. You **cannot** compute "TPM across samples" — that does not exist. Each sample has its own TPM column. For cross-sample comparison, compare a gene's TPM across columns; do not normalize the matrix a second time.

**Pitfall 2: Using TPM for differential expression.**
DESeq2 and edgeR want raw integer counts and apply their own size-factor / TMM normalization internally. If you feed them TPM, you double-normalize and the p-values are wrong. Always pass raw counts to a DE tool.

**Pitfall 3: Not rounding non-integer counts.**
kallisto/Salmon emit non-integer `est_counts` for multi-mapping reads. DESeq2 wants integers. Round (`np.round`) or, more correctly, use `tximport(..., type="kallisto")` which applies lengthScaledTPM and produces integer-suitable counts.

**Pitfall 4: Computing TPM with raw transcript length instead of effective length.**
The denominator in TPM is the effective length (transcript_length - mean_fragment_length + 1), not the transcript length itself. For 200 bp fragments and a 500 bp gene, raw_length = 500 and eff_length = 301; using 500 introduces a 40% error on short transcripts. kallisto and Salmon report `eff_length` directly; use it.

**Pitfall 5: Forgetting to filter low-count genes before PCA.**
A 3-sample yeast counts matrix has ~6,000 rows. Of those, ~2,000 have zero or near-zero expression in all samples (lowly expressed transcripts, pseudogenes, etc.). PCA on the full matrix is dominated by noise from these low-count rows. Standard practice: filter to genes with `CPM > 1 in ≥ 2 samples` before PCA. This typically retains ~3,000-4,000 yeast genes.

---

## 10. Choosing the right normalization for the question

| Question | Right normalization |
|----------|---------------------|
| Is gene X up in condition A vs condition B? | Raw counts → DESeq2 → log2 fold change + adjusted p-value |
| What does the per-sample expression profile look like? | log2(TPM + 1) heatmap or PCA |
| What fraction of reads is on chrM? | CPM, summed over chrM genes |
| Should I include gene X in the DE analysis (is it expressed at all)? | CPM > 1 in ≥ N samples filter |
| What is the expression of gene X in absolute terms? | TPM, with a note that "TPM 100 means X captures 0.01% of the transcript mass" |
| What is the ratio of gene X to gene Y in this sample? | TPM(X) / TPM(Y), because TPM is normalized for length |
| What is the same ratio across two samples? | TPM(X)_S1 / TPM(Y)_S1 vs TPM(X)_S2 / TPM(Y)_S2 |
| Compare to a 2010 paper that reports RPKM | Convert via TPM = RPKM × 10^6 / sum(RPKMs) |

---

## 11. What to remember

- **Raw counts in, normalized values out.** DE tools want raw integers. Visualization wants log-TPM. Filtering wants CPM.
- **TPM sums to 10^6 per sample, by construction.** That is what makes it the right cross-sample normalization. RPKM does not have this property.
- **Effective length, not transcript length.** The TPM denominator is `length - mean_fragment_length + 1`, not `length`. kallisto and Salmon report it; use it.
- **Aggregate transcripts to genes before downstream analysis.** Most DE tools, GO tools, and pathway tools operate at the gene level. Sum across transcripts of each gene; or use `tximport` for the lengthScaledTPM-aware aggregation on human/mouse.
- **Filter low-count genes before PCA.** Standard threshold: CPM > 1 in ≥ N samples. Reduces noise from pseudogenes and lowly expressed transcripts.
- **Never feed TPM into a DE tool.** DESeq2 / edgeR want raw integer counts and apply their own normalization. Feeding TPM produces wrong p-values.
- **RPKM/FPKM is deprecated.** Report TPM and raw counts. If a reviewer asks for FPKM (some still do), you can compute it on the fly, but flag the deprecation in the methods.

This closes the lecture trio. **Continue to the exercises** to put the toolchain to work on a real yeast RNA-seq sample.
