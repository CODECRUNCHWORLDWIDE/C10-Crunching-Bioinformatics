# Challenge 2 — Compare quantifiers on the same sample

> **Estimated time:** 2 hours.
> **Goal:** Run kallisto, Salmon, and HISAT2 + featureCounts on the same yeast RNA-seq sample (`SRR453568`), normalize their per-gene outputs to a comparable scale, compute pairwise concordance metrics, and produce a written interpretation of where the three tools agree and where they disagree. The point is not to declare a winner; the point is to map the corners where alignment, pseudoalignment, and selective-alignment-pseudoalignment diverge — and to defend that map in writing.

Lecture 2 introduced the three architectures and asserted "they agree to ~99% on the top 5,000 expressed genes." This challenge has you measure that agreement on a real sample and produce a per-gene comparison you can hand to a reviewer.

---

## Background — the three pipelines in one paragraph each

**`kallisto quant`** (Bray et al. 2016) builds a colored de Bruijn graph of the transcriptome, pseudoaligns each read to a compatibility class of transcripts via k-mer lookup, and runs EM to estimate per-transcript abundance. ~30 sec per sample, no BAM, no alignment.

**`salmon quant --validateMappings`** (Patro et al. 2017) does the same pseudoalignment + EM as kallisto, then adds a selective-alignment score check between each read and each candidate transcript to drop spurious k-mer matches. Adds bias models (`--gcBias`, `--seqBias`). ~60 sec per sample, marginally more accurate on noisy samples.

**`HISAT2` + `samtools sort` + `featureCounts`** (Kim et al. 2019 + Liao et al. 2014) is the classical alignment path: a splice-aware FM-index-based aligner produces a sorted BAM, and `featureCounts` walks the BAM to tally reads per gene against a GTF annotation. ~3-5 min per sample (most in alignment), produces a BAM that downstream tools can reuse.

---

## Prerequisites

- Exercises 1 and 2 are done. You have `trim/SRR453568_{1,2}.trim.fq.gz` and `quant/SRR453568/abundance.tsv` (kallisto).
- A conda environment with `kallisto 0.50.1`, `salmon 1.10.2`, `hisat2 2.2.1`, `subread 2.0.6`, `samtools 1.19`, `pandas 2.2`. Canonical install:
  ```bash
  conda install -c bioconda kallisto=0.50.1 salmon=1.10.2 hisat2=2.2.1 \
      subread=2.0.6 samtools=1.19
  ```
- The yeast Ensembl release 110 cDNA FASTA from Exercise 2.
- The yeast Ensembl release 110 **genome FASTA** and **GTF** for HISAT2 + featureCounts. Download:
  ```bash
  curl -sLO http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
  gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
  curl -sLO http://ftp.ensembl.org/pub/release-110/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
  gunzip Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz
  ```
- ~500 MB free disk for the BAM and indices.

---

## Step 1 — Run kallisto (already done in Exercise 2)

You already have `quant/kallisto/SRR453568/abundance.tsv`. Symlink or copy to a standard location:

```bash
mkdir -p compare/kallisto/SRR453568
cp quant/SRR453568/abundance.tsv compare/kallisto/SRR453568/
cp quant/SRR453568/run_info.json compare/kallisto/SRR453568/
```

---

## Step 2 — Run Salmon

Build the Salmon index:

```bash
mkdir -p index
salmon index \
    -t ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
    -i index/sce_salmon/ \
    -k 31
```

Expected runtime: ~30 sec, index size ~10 MB.

Quantify SRR453568:

```bash
mkdir -p compare/salmon
salmon quant \
    -i index/sce_salmon/ \
    -l A \
    -1 trim/SRR453568_1.trim.fq.gz \
    -2 trim/SRR453568_2.trim.fq.gz \
    -o compare/salmon/SRR453568/ \
    -p 4 \
    --validateMappings \
    --gcBias --seqBias
```

Expected runtime: ~60 sec.

Output: `compare/salmon/SRR453568/quant.sf` with columns `Name, Length, EffectiveLength, TPM, NumReads`.

Spot-check:

```bash
head -5 compare/salmon/SRR453568/quant.sf
# Expected (column names are different from kallisto, but content is parallel):
# Name      Length  EffectiveLength  TPM       NumReads
# YHR055C   465     257.0            234.5     8.0
# ...
```

---

## Step 3 — Run HISAT2 + featureCounts

### Build the HISAT2 index

```bash
hisat2-build ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa index/sce_hisat2
```

Expected runtime: ~30 sec, index size ~30 MB (eight `.ht2` files).

### Align

```bash
mkdir -p compare/hisat2

hisat2 \
    -x index/sce_hisat2 \
    -1 trim/SRR453568_1.trim.fq.gz \
    -2 trim/SRR453568_2.trim.fq.gz \
    -p 4 \
    --dta \
    --summary-file compare/hisat2/SRR453568.hisat2.log \
| samtools sort -o compare/hisat2/SRR453568.bam -

samtools index compare/hisat2/SRR453568.bam
```

Expected runtime: ~3 min. Expected output: a sorted, indexed BAM (~250 MB) and a HISAT2 summary log.

The summary log reports the alignment rate. For a clean yeast sample, expect ~95% overall alignment rate, with ~85% of pairs concordantly aligned exactly once.

### Count with featureCounts

```bash
featureCounts \
    -a ref/Saccharomyces_cerevisiae.R64-1-1.110.gtf \
    -t exon \
    -g gene_id \
    -p --countReadPairs \
    -s 2 \
    -T 4 \
    -o compare/hisat2/SRR453568.counts.tsv \
    compare/hisat2/SRR453568.bam
```

Expected runtime: ~10 sec. Expected output: a 7-column TSV with one row per yeast gene (~6,000 rows).

Inspect:

```bash
head -3 compare/hisat2/SRR453568.counts.tsv
# Geneid    Chr   Start   End   Strand   Length   compare/hisat2/SRR453568.bam
# YAL001C   I     147594  151169   -   3477     198
# YAL002W   I     143707  147531   +   3825     156

cat compare/hisat2/SRR453568.counts.tsv.summary
# Status                          compare/hisat2/SRR453568.bam
# Assigned                        2741234
# Unassigned_Unmapped             ...
# Unassigned_NoFeatures           ...
# ...
```

The `Assigned` count is the per-sample effective library size for featureCounts. It is typically 5-10% lower than kallisto's `n_pseudoaligned` because featureCounts (with default settings) drops multi-mappers, whereas kallisto's EM assigns them fractionally.

---

## Step 4 — Build a per-gene, per-tool comparison table

Write `compare/compare_tools.py`:

```python
"""Compare per-gene counts and TPM across kallisto, Salmon, and HISAT2+featureCounts.

For yeast, the transcript-to-gene mapping is essentially 1:1
(target_id == gene_id for most rows), so we treat per-transcript and
per-gene as the same.
"""

from __future__ import annotations
import pandas as pd
import numpy as np
from pathlib import Path


def load_kallisto(quant_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(quant_dir / "abundance.tsv", sep="\t")
    df = df.rename(columns={"target_id": "gene_id"})
    return df[["gene_id", "est_counts", "tpm"]].rename(
        columns={"est_counts": "counts_kallisto", "tpm": "tpm_kallisto"}
    )


def load_salmon(quant_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(quant_dir / "quant.sf", sep="\t")
    df = df.rename(columns={"Name": "gene_id"})
    return df[["gene_id", "NumReads", "TPM"]].rename(
        columns={"NumReads": "counts_salmon", "TPM": "tpm_salmon"}
    )


def load_featurecounts(counts_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(counts_tsv, sep="\t", comment="#")
    bam_col = [c for c in df.columns if c.endswith(".bam")][0]
    df = df.rename(columns={"Geneid": "gene_id", bam_col: "counts_fc", "Length": "length_fc"})
    df["counts_fc"] = df["counts_fc"].astype(float)
    # featureCounts does not emit TPM; compute it from counts and length.
    eff_length = (df["length_fc"] - 200.0).clip(lower=1.0)
    rate = df["counts_fc"] / eff_length
    df["tpm_fc"] = 1e6 * rate / rate.sum()
    return df[["gene_id", "counts_fc", "tpm_fc"]]


def merge_three(kal: pd.DataFrame, sal: pd.DataFrame, fc: pd.DataFrame) -> pd.DataFrame:
    df = kal.merge(sal, on="gene_id", how="outer").merge(fc, on="gene_id", how="outer")
    df = df.fillna(0.0)
    return df


def pairwise_concordance(merged: pd.DataFrame, min_tpm: float = 1.0) -> dict[str, float]:
    """Compute Pearson r on log2(TPM+1) for each pair, on the well-expressed subset."""
    expressed = merged[
        (merged["tpm_kallisto"] >= min_tpm)
        | (merged["tpm_salmon"] >= min_tpm)
        | (merged["tpm_fc"] >= min_tpm)
    ].copy()
    log = np.log2(expressed[["tpm_kallisto", "tpm_salmon", "tpm_fc"]] + 1.0)
    pairs = {
        "kallisto_vs_salmon": float(log["tpm_kallisto"].corr(log["tpm_salmon"])),
        "kallisto_vs_fc": float(log["tpm_kallisto"].corr(log["tpm_fc"])),
        "salmon_vs_fc": float(log["tpm_salmon"].corr(log["tpm_fc"])),
    }
    return pairs


if __name__ == "__main__":
    kal = load_kallisto(Path("compare/kallisto/SRR453568"))
    sal = load_salmon(Path("compare/salmon/SRR453568"))
    fc = load_featurecounts(Path("compare/hisat2/SRR453568.counts.tsv"))

    merged = merge_three(kal, sal, fc)
    merged.to_csv("compare/three_tools.tsv", sep="\t", index=False, float_format="%.4f")

    pairs = pairwise_concordance(merged, min_tpm=1.0)
    for k, v in pairs.items():
        print(f"  {k}: Pearson r on log2(TPM+1) = {v:.4f}")

    # Top 10 most discordant genes (by TPM ratio kallisto / featureCounts).
    expressed = merged[
        (merged["tpm_kallisto"] > 10) & (merged["tpm_fc"] > 10)
    ].copy()
    expressed["log2_ratio_kal_fc"] = np.log2(
        (expressed["tpm_kallisto"] + 1.0) / (expressed["tpm_fc"] + 1.0)
    )
    print()
    print("Top 10 most discordant (kallisto vs featureCounts):")
    print(
        expressed.reindex(
            expressed["log2_ratio_kal_fc"].abs().sort_values(ascending=False).index
        )
        .head(10)[["gene_id", "tpm_kallisto", "tpm_salmon", "tpm_fc", "log2_ratio_kal_fc"]]
        .to_string(index=False)
    )
```

Run:

```bash
python compare/compare_tools.py
```

Expected output for SRR453568:

```
  kallisto_vs_salmon: Pearson r on log2(TPM+1) = 0.9912
  kallisto_vs_fc: Pearson r on log2(TPM+1) = 0.9871
  salmon_vs_fc: Pearson r on log2(TPM+1) = 0.9883

Top 10 most discordant (kallisto vs featureCounts):
 gene_id  tpm_kallisto  tpm_salmon  tpm_fc  log2_ratio_kal_fc
 RDN37-1     2341.21      2401.55    14.32       7.34
 RDN37-2     2102.45      2154.21    11.87       7.46
 RDN18-1     1832.11      1875.23     8.92       7.68
 ...
```

The top discordant genes are essentially always **rRNA loci** (RDN37, RDN18, RDN58, RDN25, RDN5) or **paralogs** (some tRNA clusters, some histone genes). The discordance is from multi-mapping treatment: kallisto/Salmon assign reads fractionally across the paralogous copies; featureCounts (default) drops multi-mappers entirely. Both behaviors are defensible; they answer slightly different questions.

---

## Step 5 — Write up

Create `compare/README.md`:

```markdown
# Tool comparison — SRR453568

## Methods

- Sample: SRR453568, yeast galactose-growth replicate from Gierlinski et al. 2015.
- Trimmed with fastp 0.23.4 (Exercise 1).
- Quantified with:
  - kallisto 0.50.1 against Ensembl release 110 yeast cDNA (Exercise 2 output).
  - Salmon 1.10.2 against the same cDNA, with --validateMappings --gcBias --seqBias.
  - HISAT2 2.2.1 against the Ensembl release 110 yeast genome, then featureCounts (subread 2.0.6) -t exon -g gene_id -p --countReadPairs -s 2 against the matching GTF.

## Results

| Metric                                | kallisto  | Salmon   | HISAT2+featureCounts |
|---------------------------------------|----------:|---------:|---------------------:|
| Effective library size (assigned reads)| 2,914,012 | 2,927,431 | 2,741,234            |
| Assignment rate                       | 94.8%     | 95.3%    | 89.2%                |
| Wall-clock (excluding index build)    | 38 sec    | 62 sec   | 195 sec              |
| Genes with TPM > 10                   | 3,812     | 3,847    | 3,724                |

## Pairwise log2(TPM+1) Pearson correlations

| Pair                  | Pearson r |
|-----------------------|----------:|
| kallisto vs Salmon    | 0.9912    |
| kallisto vs featureCounts | 0.9871|
| Salmon vs featureCounts | 0.9883  |

## Where they disagree

The top 20 most discordant genes between kallisto and featureCounts (by log2 TPM ratio) are dominated by rRNA loci (RDN37-1, RDN37-2, RDN18-1, RDN5-1, RDN58-1, RDN25-1) and a handful of paralogous tRNAs. featureCounts (default settings) drops multi-mapping reads; kallisto/Salmon assign them fractionally via EM. On rRNA, where the multi-mapping is rampant (multiple identical copies of the rRNA operon in the yeast genome), the discordance is large: kallisto reports TPM ~2,000 for RDN37-1; featureCounts reports TPM ~14 because most reads on RDN37 are multi-mappers and discarded.

Outside the multi-mapping subset (>95% of genes), the three tools agree to within ~5% on log2(TPM+1). For a vanilla differential expression analysis on non-rRNA protein-coding genes, the choice between the three is essentially a matter of wall-clock time and disk usage; the biology is the same.

## Verdict

For this project's downstream use case (differential expression on protein-coding genes), all three tools produce equivalent results. kallisto is the fastest and the simplest; Salmon is marginally more accurate but slower; HISAT2+featureCounts is the slowest, the largest on disk, and the only one that produces a BAM. For the mini-project, we proceed with kallisto.
```

---

## Acceptance criteria

- [ ] `compare/kallisto/SRR453568/abundance.tsv`, `compare/salmon/SRR453568/quant.sf`, and `compare/hisat2/SRR453568.counts.tsv` all exist.
- [ ] `compare/three_tools.tsv` exists with one row per gene and columns `gene_id, counts_kallisto, tpm_kallisto, counts_salmon, tpm_salmon, counts_fc, tpm_fc`.
- [ ] `compare/compare_tools.py` runs and prints the three pairwise Pearson correlations.
- [ ] All three correlations are ≥ 0.97 on `log2(TPM+1)` for the well-expressed subset.
- [ ] `compare/README.md` is committed with a 200-400 word write-up.
- [ ] Commit message like `c2: three-tool comparison on SRR453568, all pairwise r > 0.99`.

---

## Stretch goals

- **Scatter plots.** Make matplotlib scatters of `log2(TPM+1)` for each pair of tools, color-coded by whether the gene is in a paralog family. Save to `compare/scatter_kal_vs_fc.png` etc.
- **Decoy-aware Salmon.** Build a Salmon index with the yeast genome as decoys (Lecture 2 §5) and re-run. The kallisto-vs-Salmon-decoy correlation should be marginally higher than kallisto-vs-Salmon-no-decoy on this sample (yeast has too small a genome for decoys to matter much; on human the effect is ~1-2%).
- **STAR.** Re-do the alignment-based path with STAR 2.7.11a instead of HISAT2. STAR's `--quantMode GeneCounts` produces a featureCounts-equivalent table directly. Compare to the HISAT2+featureCounts result; expect agreement to within ~1%.
- **Run all three on a glucose sample (SRR453566) and compare to the galactose result.** Compute differential expression by hand: `log2(TPM_galactose + 1) - log2(TPM_glucose + 1)` per gene. The top up-regulated genes in galactose should include GAL1, GAL7, GAL10, GAL2, GAL3 with log2 fold change > 5. This is a Week 8 preview.

---

## What you learned

- All three quantification architectures agree on ~99% of well-expressed protein-coding genes.
- The disagreements are concentrated in multi-mapping subsets: rRNA, paralog families, repetitive regions. Choice of tool matters for these subsets.
- kallisto is the fastest; Salmon is marginally more accurate; HISAT2+featureCounts is the only one that produces a BAM. Pick based on what the downstream analysis needs.
- featureCounts default drops multi-mappers; kallisto/Salmon assign them fractionally. The two behaviors answer slightly different biological questions and produce systematically different counts on paralogs.
- Pearson r on log2(TPM+1) is the standard concordance metric for quantifier comparison. ≥ 0.99 on the well-expressed subset is the bar.
- For methods sections, report the tool, the version, the index source, the alignment/assignment rate, and the wall-clock time. Reviewers want all four.
