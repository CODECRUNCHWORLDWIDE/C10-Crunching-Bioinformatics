# Exercise Solutions — Week 7

> Reference implementations and worked numerical expectations for Exercises 1-3. Read these only after you have written your own answers. The point of an exercise is the friction of doing it the first time; the solution exists to clarify after the fact.

---

## Exercise 1 — Trim a paired-end RNA-seq FASTQ with fastp

### Reference command

```bash
fastp \
    -i raw/SRR453568_1.fastq.gz \
    -I raw/SRR453568_2.fastq.gz \
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

### Reference fastp_summary.py

```python
from __future__ import annotations
import json
from pathlib import Path


def summarize_fastp_json(json_path: Path) -> dict[str, float | int | str]:
    """Read a fastp JSON report and return key QC fields."""
    with json_path.open() as f:
        data = json.load(f)

    summary = data["summary"]
    before = summary["before_filtering"]
    after = summary["after_filtering"]

    reads_in: int = before["total_reads"]
    reads_out: int = after["total_reads"]
    pct_retained: float = 100.0 * reads_out / max(reads_in, 1)
    pct_q30: float = 100.0 * after["q30_bases"] / max(after["total_bases"], 1)

    adapter = data.get("adapter_cutting", {})
    adapter_pct: float = 100.0 * adapter.get("adapter_trimmed_reads", 0) / max(reads_in, 1)

    dup = data.get("duplication", {})
    duplication_rate: float = dup.get("rate", 0.0)

    insert = data.get("insert_size", {})
    insert_size_peak: int = insert.get("peak", 0)

    return {
        "reads_in": reads_in,
        "reads_out": reads_out,
        "pct_retained": pct_retained,
        "pct_q30": pct_q30,
        "adapter_pct": adapter_pct,
        "duplication_rate": duplication_rate,
        "insert_size_peak": insert_size_peak,
    }


if __name__ == "__main__":
    summary = summarize_fastp_json(Path("qc/SRR453568.fastp.json"))
    for k, v in summary.items():
        print(f"  {k}: {v}")
```

### Expected numbers for SRR453568

| Field                | Value      |
|----------------------|-----------:|
| reads_in             | 6,400,000  |
| reads_out            | 6,144,824  |
| pct_retained         | 96.0%      |
| pct_q30              | 93.5%      |
| adapter_pct          | 35.2%      |
| duplication_rate     | 0.42       |
| insert_size_peak     | 210        |

All five gating fields are in the "Healthy" band, so the sample is OK to pseudoalign. The methods-section summary line: *"Reads were adapter- and quality-trimmed with fastp 0.23.4 (--qualified_quality_phred 20 --length_required 36 --trim_poly_g), retaining 6,144,824 of 6,400,000 reads (96.0%). Median Q30 was 93.5%; duplication rate 0.42; insert-size peak 210 bp."*

### Common mistakes

- **Forgetting `--detect_adapter_for_pe`** — fastp will fall back to declared adapter sequences and miss non-standard kits. For paired-end Illumina, always use the overlap method.
- **Setting `--qualified_quality_phred 15`** — the fastp default. For Q30-targeted libraries (most modern Illumina), 20 is a better threshold; 15 lets through more low-quality data.
- **Not committing the HTML report.** It is ~600 KB, version-controllable, and the reviewer-readable proof that your trimming was reasonable. Commit it.

---

## Exercise 2 — Pseudoalign with kallisto

### Reference commands

```bash
# Build index.
kallisto index -i index/sce.idx -k 31 ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# Quantify.
kallisto quant \
    -i index/sce.idx \
    -o quant/SRR453568/ \
    -t 4 \
    -b 100 \
    trim/SRR453568_1.trim.fq.gz \
    trim/SRR453568_2.trim.fq.gz
```

### Expected numbers for SRR453568

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

### Top expressed transcripts (galactose-induced sample)

| target_id | gene (Ensembl) | est_counts | TPM     |
|-----------|----------------|-----------:|--------:|
| YGL103W   | RPL28          |  5,821     | 18,427.6 |
| YJL177W   | RPL17B         |  5,443     | 16,842.5 |
| YBR011C   | IPP1           |  3,201     | 15,673.2 |
| YBR020W   | GAL1           |  8,432     |  2,894.4 |
| YBR018C   | GAL7           |  7,102     |  2,547.8 |
| YBR019C   | GAL10          |  8,765     |  2,487.9 |

The dominant species (`RPL*`, `IPP1`) are housekeeping; the diagnostic biological signature is the `GAL1/7/10` regulon at TPM > 2,000, which is what tells you the sample was grown on galactose.

### Common mistakes

- **Confusing the cDNA FASTA with the genome FASTA.** kallisto wants transcripts, not chromosomes. The Ensembl URL is `/cdna/`, not `/dna/`.
- **Single-end-mode FASTQs treated as paired.** Verify `_1` and `_2` are the same read count and are properly paired. `seqkit stats` is a quick sanity check.
- **Skipping `-b 100`.** Without bootstrap, downstream tools like `sleuth` cannot quantify uncertainty. For Week 7 we mostly use point estimates, but the bootstrap is cheap (~30 sec) and worth saving.

### Why `est_counts` is not always an integer

Multi-mapping reads are fractionally assigned by EM. A read whose compatibility class is `{YGL103W, YGL103C-A}` (paralogs with high sequence identity) gets a fractional count assigned to each, weighted by the EM-inferred relative abundance. The downstream `tximport` package handles the fractional counts correctly; DESeq2 wants integers, so round in the aggregation step.

---

## Exercise 3 — TPM and CPM from counts

### Reference compute_cpm

```python
def compute_cpm(counts: pd.Series) -> pd.Series:
    library_size: float = float(counts.sum())
    if library_size <= 0.0:
        raise ValueError("library_size <= 0")
    return 1e6 * counts / library_size
```

### Reference compute_tpm

```python
def compute_tpm(counts: pd.Series, eff_length: pd.Series) -> pd.Series:
    rate = counts / eff_length
    total_rate: float = float(rate.sum())
    if total_rate <= 0.0:
        raise ValueError("total_rate <= 0")
    return 1e6 * rate / total_rate
```

### The TPM identity

`sum(TPM) = 1e6`, always, by construction. The proof:

```
sum_g(TPM_g) = sum_g(1e6 * (c_g / l_g) / S)
             = (1e6 / S) * sum_g(c_g / l_g)
             = (1e6 / S) * S        (where S = sum_g(c_g / l_g))
             = 1e6.
```

The CPM identity is the same: `sum(CPM) = sum(1e6 * c_g / N) = 1e6 * N / N = 1e6` (where N is library size).

### Verifying against kallisto's TPM

On the SRR453568 sample, the computed-vs-kallisto TPM should match to < 0.5% relative error for any transcript with `est_counts > 10`. If it does not:

- Check that you used `eff_length` (not `length`) in the denominator. The 200 bp fragment-length adjustment is what kallisto uses internally.
- Check that the `est_counts` column was loaded as float, not int. The fractional values from multi-mapping resolution matter.
- Check that you did not pre-filter the input. The TPM denominator is the sum across **all** transcripts in the sample; dropping low-count transcripts before computing TPM produces inflated TPMs.

### Why CPM ≠ TPM

CPM normalizes by library size only; TPM normalizes by library size **and** effective transcript length. They agree only when all transcripts have the same effective length (which is essentially never true in real data).

A concrete contrast on a toy 3-transcript example:

| transcript | counts | eff_length | CPM      | TPM      |
|------------|-------:|-----------:|---------:|---------:|
| A (long)   | 1000   | 5000       | 333,333  |  43,475  |
| B (medium) | 1000   | 1000       | 333,333  | 260,870  |
| C (short)  | 1000   | 200        | 333,333  | 696,655  |

CPM says all three are equally expressed (same fraction of the library), which is **wrong** for cross-transcript comparison because the long transcript generates more reads per molecule.

TPM says C is the most abundant in molecules-per-million, which is **right**: the short transcript captures the same number of reads as the long transcript despite being 25x shorter, so its per-molecule abundance is 25x higher.

This is why every modern RNA-seq tutorial uses TPM for visualization.

### Common mistakes

- **Using `length` instead of `eff_length` in TPM denominator.** Introduces a ~10-40% error on short transcripts (eff_length ≈ length - 200; for length = 500 the error is 40%).
- **Computing TPM after filtering low-count transcripts.** TPM is a sum-to-million normalization. If you drop 30% of transcripts (even if they have zero counts), the remaining TPMs sum to 1e6 over a smaller denominator, so the values are 30% inflated. Compute TPM first, then filter.
- **Computing TPM across samples.** TPM is a within-sample normalization. There is no such thing as "TPM normalized across the matrix"; each column has its own sum-to-million constraint.

### Connecting back to Lecture 3

The formula and the identity verification in this exercise are the algebraic core of Lecture 3. If you can write both functions from memory and explain why `sum(TPM) = 1e6` in one breath, you have internalized the most important normalization in modern RNA-seq.

---

## Where to go next

- Three samples? See the **mini-project README** — it builds a 3-sample matrix from these three per-sample `abundance.tsv` outputs, plus the transcript-to-gene aggregation from Lecture 3 §5.
- Want to compare quantifiers? **Challenge 2** runs kallisto, Salmon, and HISAT2 + featureCounts on the same sample and quantifies the agreement.
- Want to scale to a 10-sample pipeline with reproducible orchestration? **Challenge 1** wires fastp + kallisto into a Snakemake workflow.
