# Exercise 2 — Pseudoalign a yeast RNA-seq sample with kallisto

> **Estimated time:** 45 minutes (10 minutes download, 5 minutes index, 5 minutes quant, 25 minutes interpretation).
> **Goal:** Download a yeast transcriptome FASTA from Ensembl, build a kallisto index, run `kallisto quant` on the trimmed FASTQs from Exercise 1, read the resulting `abundance.tsv`, and verify the top expressed transcripts make biological sense.

This is the central RNA-seq exercise. By the end you have a per-transcript counts and TPM table for one yeast sample. You reuse this output in Exercise 3 (TPM/CPM computation) and the mini-project (three-sample matrix).

---

## Background

Lecture 2 introduced the kallisto pipeline:

```bash
# 1. Build index from a transcriptome FASTA.
kallisto index -i index.idx Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# 2. Quantify a paired-end sample.
kallisto quant \
    -i index.idx \
    -o quant/SRR453568/ \
    -t 4 \
    -b 100 \
    trim/SRR453568_1.trim.fq.gz \
    trim/SRR453568_2.trim.fq.gz
```

The index is a colored de Bruijn graph of k-mers (k=31 by default) over the transcriptome. The quant step pseudoaligns each read to a compatibility class of transcripts and runs EM to estimate per-transcript abundance.

This exercise walks through both steps on a real yeast sample.

---

## Prerequisites

- Exercise 1 is done. You have `trim/SRR453568_1.trim.fq.gz` and `trim/SRR453568_2.trim.fq.gz`.
- A conda environment with `kallisto 0.50.1` installed. The canonical install: `conda install -c bioconda kallisto=0.50.1`.
- ~50 MB of free disk for the index + ~30 MB for the quant output.
- An internet connection to pull the Ensembl yeast cDNA FASTA (~3.5 MB).

---

## Step 1 — Download the yeast transcriptome

The yeast (Saccharomyces cerevisiae) cDNA FASTA from Ensembl release 110:

```bash
mkdir -p ref
cd ref

curl -sLO http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# Sanity check.
zcat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz | head -4
# Expected:
# >YHR055C cdna chromosome:R64-1-1:VIII:215248:215820:-1 gene:YHR055C ...
# ATGCAACAATTAAATGCATTATATATAGGGGCAAGCTTGCCCATCCAAAGTAACAATAAA
# ...

zcat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz | grep -c "^>"
# Expected: ~6,975 transcripts

cd ..
```

Each FASTA record is one yeast transcript. The header contains the transcript ID (e.g. `YHR055C`), the genomic coordinates, and the gene metadata. For yeast, the transcript ID and the gene ID are usually identical (one transcript per gene) — this simplifies the transcript-to-gene aggregation in Exercise 3.

---

## Step 2 — Build the kallisto index

```bash
mkdir -p index

kallisto index \
    -i index/sce.idx \
    -k 31 \
    ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
```

Expected runtime: ~10 seconds on a laptop. Expected console output:

```
[build] loading fasta file ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
[build] k-mer length: 31
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 6975 contigs and contains 8456321 k-mers
```

The index is ~50 MB on disk. Verify:

```bash
ls -lh index/sce.idx
# Expected: ~50 MB

kallisto inspect index/sce.idx
# Expected: k-mer length 31, number of targets 6975, ...
```

---

## Step 3 — Quantify the trimmed FASTQs

```bash
mkdir -p quant

kallisto quant \
    -i index/sce.idx \
    -o quant/SRR453568/ \
    -t 4 \
    -b 100 \
    trim/SRR453568_1.trim.fq.gz \
    trim/SRR453568_2.trim.fq.gz
```

Expected runtime: ~30-60 seconds (most of which is the 100 bootstrap replicates). Without `-b 100`, ~15 seconds.

Expected console output:

```
[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 6,975
[index] number of k-mers: 8,456,321
[index] number of equivalence classes: 9,873
[quant] running in paired-end mode
[quant] will process pair 1: trim/SRR453568_1.trim.fq.gz
                              trim/SRR453568_2.trim.fq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 3,072,412 reads, 2,914,012 reads pseudoaligned
[quant] estimated average fragment length: 208.45
[em] quantifying the abundances ... done
[em] the Expectation-Maximization algorithm ran for 178 rounds
[bstrp] running 100 bootstrap replicates ... done
```

Record:

- `processed`: total trimmed reads.
- `pseudoaligned`: reads that matched the transcriptome.
- `pseudoalignment rate`: pseudoaligned / processed, in percent.
- `mean fragment length`: kallisto's empirical estimate.
- `EM rounds`: convergence iterations.

For a healthy yeast sample, expect pseudoalignment rate **≥ 90%**. A rate of 75-90% is marginal (often rRNA contamination or genomic DNA contamination). A rate < 75% suggests a wrong reference or a severely contaminated sample.

---

## Step 4 — Read the abundance.tsv

The quant output directory contains:

```
quant/SRR453568/
├── abundance.h5       (HDF5: counts + bootstrap matrix)
├── abundance.tsv      (TSV: counts + TPM, the human-readable form)
└── run_info.json      (metadata: command, version, read counts)
```

The TSV has five columns: `target_id, length, eff_length, est_counts, tpm`. Open it with pandas:

```python
"""Inspect kallisto's abundance.tsv for SRR453568."""

from __future__ import annotations
import pandas as pd
from pathlib import Path


def load_abundance(quant_dir: Path) -> pd.DataFrame:
    """Read kallisto abundance.tsv into a pandas DataFrame."""
    return pd.read_csv(quant_dir / "abundance.tsv", sep="\t")


if __name__ == "__main__":
    df = load_abundance(Path("quant/SRR453568"))

    print(f"  n_transcripts:    {len(df):,}")
    print(f"  n_quantified:     {(df['est_counts'] > 0).sum():,}")
    print(f"  total_est_counts: {df['est_counts'].sum():,.0f}")
    print(f"  mean_eff_length:  {df['eff_length'].mean():.1f}")

    print()
    print("Top 10 transcripts by TPM:")
    print(df.sort_values("tpm", ascending=False).head(10)[
        ["target_id", "length", "est_counts", "tpm"]
    ].to_string(index=False))
```

Expected output:

```
  n_transcripts:    6,975
  n_quantified:     4,832
  total_est_counts: 2,914,012
  mean_eff_length:  1,256.4

Top 10 transcripts by TPM:
 target_id  length  est_counts        tpm
   YGL103W     465     5821.34   18427.61
   YJL177W     477     5443.87   16842.55
   YBR011C     303     3201.21   15673.18
   YML063W     288     2987.45   15203.97
   YGL031C     447     4123.78   14012.83
   YDR450W     297     3001.12   13987.62
   YOL127W     387     3567.21   13452.31
   YGR234W     507     3987.66   13201.45
   YDR133C     321     2876.21   12987.32
   YHR141C     369     2901.45   12873.85
```

Most of these are **ribosomal protein genes** (`YGL103W` is `RPL28`, `YJL177W` is `RPL17B`, etc.). This is the expected signature of any RNA-seq sample: ribosomal proteins dominate the top of the expression distribution because they are highly expressed in every growth condition.

---

## Step 5 — Spot-check the galactose signature

SRR453568 is the **galactose-growth** replicate. The galactose-induced genes are well known: `GAL1, GAL2, GAL3, GAL7, GAL10, GCY1`. In a galactose sample, these should be in the top ~100 expressed; in a glucose sample, they should be near the bottom (Mig1-repressed).

```python
"""Check the GAL regulon's expression in SRR453568."""

from __future__ import annotations
import pandas as pd
from pathlib import Path


GAL_GENES = {
    "YBR020W": "GAL1",
    "YLR081W": "GAL2",
    "YDR009W": "GAL3",
    "YBR018C": "GAL7",
    "YBR019C": "GAL10",
    "YOR120W": "GCY1",
}


if __name__ == "__main__":
    df = pd.read_csv("quant/SRR453568/abundance.tsv", sep="\t")
    gal = df[df["target_id"].isin(GAL_GENES.keys())].copy()
    gal["gene_name"] = gal["target_id"].map(GAL_GENES)
    gal = gal.sort_values("tpm", ascending=False)

    print("Galactose regulon expression in SRR453568:")
    print(gal[["target_id", "gene_name", "est_counts", "tpm"]].to_string(index=False))
```

Expected output for SRR453568 (galactose):

```
Galactose regulon expression in SRR453568:
 target_id gene_name  est_counts      tpm
   YBR020W     GAL1    8432.13   2894.42
   YBR018C     GAL7    7102.45   2547.81
   YBR019C    GAL10    8765.12   2487.93
   YLR081W     GAL2    4321.87   1542.27
   YDR009W     GAL3    1873.45    872.16
   YOR120W     GCY1    1102.34    523.61
```

All six are TPM > 500, with `GAL1/7/10` (the canonical operon) in the top ~50 expressed transcripts of the sample. This is the **textbook galactose induction signature**.

For comparison, in a glucose sample (SRR453566 or SRR453567), `GAL1` would have TPM ~5 and `GAL7/10` would have TPM ~10-30 (the residual basal expression). The fold change between glucose and galactose for these genes is ~500x, which is what Week 8's differential expression analysis will recover.

---

## Step 6 — Read the run_info.json

```python
"""Parse kallisto's run_info.json for QC numbers."""

from __future__ import annotations
import json
from pathlib import Path


def kallisto_summary(quant_dir: Path) -> dict[str, float | int | str]:
    """Read kallisto run_info.json and return a QC summary."""
    with (quant_dir / "run_info.json").open() as f:
        data = json.load(f)
    return {
        "n_targets": int(data["n_targets"]),
        "n_processed": int(data["n_processed"]),
        "n_pseudoaligned": int(data["n_pseudoaligned"]),
        "p_pseudoaligned": float(data["p_pseudoaligned"]),
        "n_unique": int(data["n_unique"]),
        "p_unique": float(data["p_unique"]),
        "kallisto_version": str(data["kallisto_version"]),
        "n_bootstraps": int(data["n_bootstraps"]),
        "frag_length_mean": float(data.get("frag_length_mean", 0)),
        "frag_length_sd": float(data.get("frag_length_sd", 0)),
    }


if __name__ == "__main__":
    s = kallisto_summary(Path("quant/SRR453568"))
    for k, v in s.items():
        print(f"  {k}: {v}")
```

Expected output:

```
  n_targets: 6975
  n_processed: 3072412
  n_pseudoaligned: 2914012
  p_pseudoaligned: 94.84
  n_unique: 2682431
  p_unique: 87.31
  kallisto_version: 0.50.1
  n_bootstraps: 100
  frag_length_mean: 208.45
  frag_length_sd: 42.13
```

The `p_pseudoaligned` of ~95% is healthy. The `p_unique` of ~87% (reads with a single-transcript compatibility class) means ~13% of pseudoaligned reads were multi-mapping and resolved fractionally by EM.

---

## Step 7 — Save the QC notes

Create `week-07/exercises/notes/e2-kallisto.md`:

```markdown
# Exercise 2 — kallisto quant for SRR453568

- Tool: kallisto 0.50.1
- Index: Ensembl release 110 yeast cDNA (6,975 transcripts, k=31)
- Input: SRR453568_1.trim.fq.gz + SRR453568_2.trim.fq.gz (3.07 M trimmed pairs)
- Bootstrap replicates: 100

## Summary

| Field             | Value      |
|-------------------|-----------:|
| n_processed       | 3,072,412  |
| n_pseudoaligned   | 2,914,012  |
| p_pseudoaligned   | 94.84%     |
| p_unique          | 87.31%     |
| frag_length_mean  | 208.45     |
| frag_length_sd    | 42.13      |
| EM rounds         | 178        |

## Top expressed transcripts

| ID        | gene  | est_counts | TPM      |
|-----------|-------|-----------:|---------:|
| YGL103W   | RPL28 | 5,821      | 18,427.6 |
| YJL177W   | RPL17B| 5,443      | 16,842.5 |
| YBR020W   | GAL1  | 8,432      |  2,894.4 |
| YBR018C   | GAL7  | 7,102      |  2,547.8 |
| YBR019C   | GAL10 | 8,765      |  2,487.9 |

## Verdict

Sample is healthy. Pseudoalignment rate 94.8% is well above the 90% threshold. The GAL1/7/10 regulon is induced (TPM > 2,000 vs the expected ~10 in glucose), consistent with the galactose growth condition declared in the metadata. Proceed to Exercise 3 (TPM/CPM computation).
```

Commit:

- `quant/SRR453568/abundance.tsv`
- `quant/SRR453568/run_info.json`
- `notes/e2-kallisto.md`

Gitignore:

- `ref/` (the transcriptome FASTA is easy to re-download)
- `index/` (the index is fast to rebuild)
- `quant/SRR453568/abundance.h5` (it is ~10 MB and the bootstrap matrix is reconstructable)

---

## Acceptance criteria

- [ ] `ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz` exists.
- [ ] `index/sce.idx` exists, ~50 MB.
- [ ] `quant/SRR453568/abundance.tsv` exists with 6,975 rows.
- [ ] `quant/SRR453568/run_info.json` exists.
- [ ] You ran the three Python snippets above and got output matching expectations (within natural sample-to-sample variation).
- [ ] `notes/e2-kallisto.md` is committed with the QC table and verdict.
- [ ] Commit message like `e2: kallisto quant SRR453568, 94.8% pseudoaligned, GAL1/7/10 induced`.

---

## Common pitfalls

**`kallisto index` complains about `unknown nucleotide`.** The transcriptome FASTA contains a non-ACGT character (usually N or X). Most Ensembl FASTAs are clean; if you see this, you may be using a custom-annotated reference. Filter out Ns with `awk '/^>/{p=1} !/^>/{if($0 !~ /N/){p=p&&1}else{p=0}} p' ref.fa > ref.clean.fa` — but check biologically first; some Ns are real (gaps in the assembly).

**`p_pseudoaligned` is 30%.** Three common causes:
1. You used the genome FASTA instead of the cDNA FASTA. Kallisto needs transcripts, not chromosomes.
2. The FASTQs are not actually paired (or are paired in the wrong order). Verify with `zcat trim/SRR453568_1.trim.fq.gz | head -4` and same for `_2`.
3. The sample is from a different organism than the index. Double-check the SRA metadata.

**Fragment-length estimate is wrong on single-end.** For single-end data, kallisto cannot estimate the fragment length from the data; you must pass `--single -l 200 -s 30` (mean=200, sd=30, or whatever the library was). The mini-project is paired-end so this does not arise.

**`est_counts` are not integers.** This is correct (Lecture 2 §3). Multi-mapping reads are fractionally assigned by EM. DESeq2 wants integers; round downstream.

**Bootstrap takes forever.** Drop `-b 100` if you only need point estimates; `-b 0` is the default. The Week 8 `sleuth` differential expression tool uses bootstrap variance, but for this exercise the point estimate is enough.

---

## What you learned

- The kallisto index is a colored de Bruijn graph of k-mers (k=31) over the transcriptome. Building it takes ~10 seconds for yeast.
- `kallisto quant` pseudoaligns reads and runs EM to estimate per-transcript abundance, in ~30 seconds per million paired reads on a laptop.
- `abundance.tsv` has five columns: `target_id, length, eff_length, est_counts, tpm`. The `tpm` column is already normalized; you can use it directly for cross-sample visualization.
- `run_info.json` is the QC source-of-truth for kallisto. `p_pseudoaligned ≥ 90%` is the health threshold.
- The GAL1/7/10 regulon shows up at TPM > 2,000 in galactose samples and TPM < 20 in glucose samples. This is a fold change of ~100-500x that Week 8's DE analysis will recover.
- Multi-mapping reads (~13% in yeast, ~25% in human) are resolved fractionally by EM. The `p_unique` field tells you what fraction was unique.

Continue to [Exercise 3 — TPM and CPM from counts (Python)](./exercise-03-tpm-cpm-from-counts.py).
