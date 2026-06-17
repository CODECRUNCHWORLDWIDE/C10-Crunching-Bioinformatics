# Lecture 1 — Multiple Sequence Alignment and the Progressive Heuristic

> **Reproducibility note.** This lecture describes the *mechanics* of multiple sequence alignment. Every MSA algorithm in widespread use is a heuristic for an NP-hard problem; two algorithms can produce two alignments of the same FASTA and both can be defensible. Pin the tool, pin the algorithm flag, pin the version, and the alignment becomes reproducible. Skip any of those and the alignment becomes a moving target.

> **Duration:** ~3 hours of reading + a brief MAFFT sanity check.
> **Outcome:** You can describe the progressive alignment heuristic, name MAFFT, Clustal Omega, and MUSCLE 5 as the three canonical free aligners, state when each is preferred, run MAFFT from Python with `subprocess.run(..., check=True)`, and parse the result with `Bio.AlignIO`.

If you only remember one thing from this lecture, remember this:

> **A multiple sequence alignment hypothesizes that the residues in each column descend from a single ancestral residue. The alignment is not a "true answer" — it is an estimate, made by a heuristic algorithm, of which residues are homologous. Two algorithms can disagree on the placement of indels, and the disagreements travel downstream into the tree. The fix is to (a) pick one algorithm and version, (b) trim columns with low confidence, (c) record both choices in the run-info, and (d) inspect the alignment by eye for any region your downstream conclusion depends on.**

Week 8's VCF is parked. Week 9's input is a FASTA of homologous sequences. Week 9's intermediate artefact is the alignment. The tree is downstream.

---

## 1. Where we are in the pipeline

The end-to-end Week 9 pipeline runs:

```
FASTA (N homologous sequences) ->
        progressive MSA (MAFFT, Clustal Omega, or MUSCLE 5) ->
        ALIGNED FASTA (N rows x M aligned columns, gaps inserted) ->
        column-trim (drop > 50% gap columns; optionally trimAl) ->
        ALIGNED FASTA (cleaner; M' <= M columns) ->
        distance matrix (Jukes-Cantor, Kimura 2-parameter) ->
        tree (NJ via Biopython; ML via IQ-TREE 2) ->
        Newick / Nexus tree file ->
        rendered PNG / SVG.
```

Lecture 1 covers the first three boxes — input FASTA, MSA, trim — and stops at the cleaned alignment. Lectures 2 and 3 take the cleaned alignment and turn it into a tree.

The single most consequential step is the MSA. A poor alignment produces a poor distance matrix produces a poor tree. The tree builder cannot recover from misaligned columns; if a column is wrongly hypothesized to be homologous, the substitutions in that column are counted as real evolutionary events when they are an alignment artefact. Felsenstein's 1981 ML framework presupposes that the columns are correctly aligned; the entire downstream apparatus inherits that assumption.

---

## 2. What is a multiple sequence alignment

A pairwise alignment between two sequences is the Week 3 problem: Needleman-Wunsch (global) or Smith-Waterman (local) dynamic programming, polynomial in the product of the two lengths. The optimal alignment is well-defined under a scoring matrix and gap penalty; Week 3 implementations recover it exactly.

A multiple sequence alignment between N sequences is the same problem with N dimensions of dynamic programming. The exact algorithm (Carrillo and Lipman 1988) is `O(L^N)` in time and `O(L^N)` in space, where `L` is the sequence length. For three sequences of length 1000 that is `10^9` cells; for ten sequences it is `10^30`. The exact MSA is computable only for the smallest of inputs (N <= 3-5, L <= 100-200) and is essentially never used in practice.

Every practical MSA tool is a heuristic. The dominant heuristic family is **progressive alignment** (Feng and Doolittle 1987; Hogeweg and Hesper 1984): order the sequences by pairwise distance, build a guide tree, then align profiles (alignment columns aggregated across sub-trees) along the guide tree from the leaves to the root. The output is an alignment that approximates the exact MSA without paying the exponential cost.

The progressive heuristic has two important failure modes:

- **Early-fix bias.** Once a gap is inserted in an early profile, the gap stays. Subsequent rounds cannot revisit the placement. This is why "once a gap, always a gap" is a thing people say in phylogenetics circles.
- **Guide-tree dependence.** A wrong guide tree (which is the entire point — we have not built the real tree yet) places gaps in the wrong order. MAFFT's iterative-refinement modes and Clustal Omega's HMM profile mode attempt to mitigate this by iterating; MUSCLE 5's super5 algorithm tackles it with a different ensemble approach.

The Week 9 workhorse, MAFFT, is a progressive aligner with optional iterative refinement.

---

## 3. MAFFT (Katoh and Standley 2013)

Citation: Katoh K, Standley DM. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." *Molecular Biology and Evolution* 30:772 (2013). Free full text at <https://academic.oup.com/mbe/article/30/4/772/1073398>. Tool docs at <https://mafft.cbrc.jp/alignment/software/>.

MAFFT is a C program with a Perl front-end. The conda recipe installs the `mafft` binary on the PATH plus support scripts. The CLI takes a FASTA on stdin or a path, writes the alignment to stdout, and accepts a small number of algorithm-choice flags.

The Week 9 default invocation is:

```bash
mafft --retree 2 --maxiterate 0 --nuc --anysymbol --quiet --thread 4 input.fasta > aligned.fasta
```

- `--retree 2` selects the FFT-NS-2 algorithm: two rounds of progressive alignment, no iterative refinement. Fastest of the MAFFT modes.
- `--maxiterate 0` disables iterative refinement (the default is 0; including it makes the choice explicit and reproducible).
- `--nuc` pins nucleotide input. MAFFT's auto-detection is reliable but pinning removes ambiguity.
- `--anysymbol` permits IUPAC ambiguity codes (N, R, Y, etc.) without warning.
- `--quiet` suppresses progress messages on stderr.
- `--thread 4` enables four worker threads.

For the demo cytochrome b panel (ten sequences, ~1,140 bp each), this runs in ~0.4 seconds and produces a 1,140-column alignment. For an L-INS-i refinement (highest-accuracy mode on small inputs), swap to:

```bash
mafft --localpair --maxiterate 1000 --nuc --anysymbol --quiet --thread 4 input.fasta > aligned.linsi.fasta
```

L-INS-i takes ~3 seconds on the same input and produces a marginally cleaner alignment (typically 2-5 fewer indel columns in a typical vertebrate cytochrome b panel). For Week 9's didactic purposes we use `--retree 2` in the exercises and L-INS-i in the challenge.

### Calling MAFFT from Python

The pattern is the same as the Week 5 / Week 8 `subprocess.run` template:

```python
from __future__ import annotations

import subprocess
from pathlib import Path


def run_mafft(input_fasta: Path, output_fasta: Path, threads: int = 4) -> None:
    """Run MAFFT in FFT-NS-2 mode. Writes the alignment to output_fasta."""
    cmd = [
        "mafft",
        "--retree", "2",
        "--maxiterate", "0",
        "--nuc",
        "--anysymbol",
        "--quiet",
        "--thread", str(threads),
        str(input_fasta),
    ]
    with output_fasta.open("w") as fh:
        result = subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
    if result.stderr:
        print(f"[mafft] {result.stderr.strip()}")
```

Note the deliberate choices:

- `check=True` raises a `subprocess.CalledProcessError` if MAFFT exits non-zero. Always prefer this to checking the return code manually.
- `stdout=fh` writes the alignment directly to disk without round-tripping through a Python string. For a 1,140-column alignment this is irrelevant; for a 50,000-column ribosomal-RNA alignment it matters.
- `text=True` decodes stderr as UTF-8.

### MAFFT versioning and the `--auto` trap

MAFFT's `--auto` flag picks one of `FFT-NS-2`, `FFT-NS-i`, `L-INS-i`, or `G-INS-i` based on the input size and sequence count. The decision boundary has changed between MAFFT 7.310, 7.490, and 7.526 — meaning that the same FASTA can pick different algorithms in different MAFFT versions and produce different alignments.

For Week 9 we **never** use `--auto`. We pin one algorithm per script and record it in the run-info JSON. The same is true of any production phylogenetics pipeline.

---

## 4. Clustal Omega (Sievers et al. 2011)

Citation: Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG. "Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega." *Molecular Systems Biology* 7:539 (2011). Free full text at <https://www.embopress.org/doi/full/10.1038/msb.2011.75>. Tool docs at <http://www.clustal.org/omega/>.

Clustal Omega is the modern successor to ClustalW (Thompson, Higgins, Gibson 1994). It uses HMM profiles and a modified guide tree to scale to very large inputs. The CLI:

```bash
clustalo -i input.fasta -o aligned.clustalo.fasta --seqtype=DNA --iter=2 --threads=4 --force
```

Clustal Omega's strengths:

- **Scales to ~190,000 sequences** in published benchmarks. MAFFT begins to slow noticeably beyond ~10,000.
- **High-quality on highly divergent inputs** because the HMM profile captures long-range conservation patterns better than the simple position-specific scoring matrix that MAFFT uses.
- **Free, open source, conda-installable**.

Clustal Omega's weaknesses for Week 9:

- **Slower than MAFFT on small inputs.** For the ten-sequence cytochrome b panel, MAFFT runs in 0.4 seconds and Clustal Omega runs in 1.2 seconds. The difference is invisible at this scale; it matters at the 1,000-sequence scale, where MAFFT is twice as fast.

We use MAFFT as the default in Week 9 exercises and the mini-project; we mention Clustal Omega as the right tool when the input size exceeds a few thousand sequences.

---

## 5. MUSCLE 5 (Edgar 2022)

Citation: Edgar RC. "Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny." *Nature Communications* 13:6968 (2022). Free full text at <https://www.nature.com/articles/s41467-022-34630-w>. Tool docs at <https://drive5.com/muscle/> and the source at <https://github.com/rcedgar/muscle>.

MUSCLE 5 is the 2022 rewrite of MUSCLE (Edgar 2004). Two algorithms in the same binary:

- `muscle -align input.fasta -output aligned.fasta`: the PPP algorithm, highest-accuracy on small-to-medium inputs.
- `muscle -super5 input.fasta -output aligned.fasta`: the super5 algorithm, faster on large inputs.

In the published *Nature Communications* benchmark (Edgar 2022), MUSCLE 5 reports the highest per-second accuracy on the BAliBASE and BENCH MSA test sets. For Week 9's small inputs, MUSCLE 5 and MAFFT produce essentially identical alignments; we mention MUSCLE 5 because it is a strong free alternative and you should know it exists.

The Week 9 exercises stick to MAFFT; one homework problem repeats Exercise 1 with MUSCLE 5 and compares the alignment column counts.

---

## 6. When to prefer which aligner

| Input size | Sequence divergence | Recommended | Why |
|------------|---------------------|-------------|-----|
| 10-200 sequences | Low-to-medium | **MAFFT L-INS-i** | Highest-accuracy MAFFT mode; small enough to afford the iterative refinement. |
| 10-200 sequences | High (< 30% identity) | **MUSCLE 5 PPP** | MUSCLE 5's PPP algorithm is the most accurate on highly divergent small inputs. |
| 200-10,000 sequences | Any | **MAFFT FFT-NS-i** | Fastest among the iterative-refinement modes. |
| > 10,000 sequences | Any | **Clustal Omega** | Scales gracefully to ~190,000 sequences; MAFFT slows. |
| Protein alignments where you suspect remote homology | Any size | **MAFFT L-INS-i + HHpred** | For the hardest cases. Out of Week 9 scope. |

The takeaway: MAFFT covers most cases. Clustal Omega scales further. MUSCLE 5 is highest-accuracy on the hardest small cases. All three are free.

---

## 7. Parsing the alignment with Biopython

After MAFFT writes the aligned FASTA, the Python-side work is parsing and inspecting. The relevant Biopython modules:

- `Bio.SeqIO` — generic FASTA read/write (Week 2).
- `Bio.AlignIO` — multiple sequence alignment read/write. Formats: `fasta`, `phylip`, `clustal`, `stockholm`, `nexus`.
- `Bio.Align.MultipleSeqAlignment` — the in-memory alignment object. Indexable rows (records) and columns.

The canonical Week 9 read:

```python
from __future__ import annotations

from pathlib import Path
from typing import Iterable


def load_alignment(aligned_fasta: Path) -> "MultipleSeqAlignment":
    """Read the MAFFT output and return a Biopython MultipleSeqAlignment.

    Lazy-imports AlignIO inside the function so the file compiles
    even on machines where Biopython is not yet installed.
    """
    from Bio import AlignIO
    return AlignIO.read(aligned_fasta, "fasta")


def column_gap_fractions(alignment: "MultipleSeqAlignment") -> list[float]:
    """Return the fraction of '-' characters per alignment column."""
    n_rows: int = len(alignment)
    n_cols: int = alignment.get_alignment_length()
    fractions: list[float] = []
    for j in range(n_cols):
        col: str = alignment[:, j]
        gap_count: int = col.count("-")
        fractions.append(gap_count / n_rows)
    return fractions


def trim_gappy_columns(
    alignment: "MultipleSeqAlignment",
    max_gap_fraction: float = 0.5,
) -> "MultipleSeqAlignment":
    """Drop alignment columns where the gap fraction exceeds max_gap_fraction.

    Returns a new MultipleSeqAlignment with the kept columns. The default
    threshold (0.5) is the conventional starting point; tighter thresholds
    (e.g. 0.2) produce smaller, cleaner alignments at the cost of
    information loss in moderately-conserved regions.
    """
    from Bio.Align import MultipleSeqAlignment

    n_cols: int = alignment.get_alignment_length()
    keep_cols: list[int] = [
        j for j in range(n_cols)
        if alignment[:, j].count("-") / len(alignment) <= max_gap_fraction
    ]
    if not keep_cols:
        raise ValueError("No columns survive the gap threshold; alignment is uniformly gappy.")
    # Build new records by slicing on the kept columns.
    new_records = []
    for record in alignment:
        kept_seq: str = "".join(str(record.seq)[j] for j in keep_cols)
        new_record = record[:0]  # empty SeqRecord; reuse the metadata
        new_record.seq = type(record.seq)(kept_seq)
        new_records.append(new_record)
    return MultipleSeqAlignment(new_records)
```

The trim step is conventional but not free. Dropping columns where more than half the rows are gaps removes mostly-uninformative columns (where the homology hypothesis is weak), but it also discards a small amount of signal in moderately-divergent regions. For a 1,140-column cytochrome b alignment, the 50% trim typically drops 30-50 columns and leaves ~1,090. For protein alignments with longer indel tracts, the trim can drop 5-15% of the columns.

A more principled trimmer is `trimAl` (Capella-Gutierrez et al. 2009, *Bioinformatics* 25:1972, free at <https://academic.oup.com/bioinformatics/article/25/15/1972/213148>). Its `-automated1` mode picks a per-alignment threshold from a small grid of options based on the alignment's overall gap distribution. We use the simple "drop > 50% gap columns" rule in the exercises for transparency; the mini-project optionally switches to `trimAl -automated1`.

---

## 8. Inspecting the alignment by eye

Every alignment that drives a downstream tree should be **looked at**. The lab notebook convention is:

1. Open the aligned FASTA in **Jalview** (free, Java GUI) or **AliView** (free, Java GUI). The web-based **MView** is the no-install alternative.
2. Scan for obvious misalignments: a row that consists mostly of gaps; a region where one row's residues do not match the column consensus; an indel that splits a conserved motif.
3. If the alignment looks wrong in a region your downstream conclusion depends on, do not "fix" it by hand. Re-run with a different algorithm (L-INS-i instead of FFT-NS-2; MUSCLE 5 instead of MAFFT) and compare. If the conclusion changes between aligners, your conclusion is not robust and you should say so in the write-up.

For headless inspection we provide a `summarize_alignment` helper that reports:

- Row count and column count.
- Per-column conservation (fraction of rows matching the column consensus).
- Per-row sequence-level metadata (length, gap fraction, GC content for nucleotide alignments).
- The five most-gappy columns.

The Week 9 exercises run this summarizer on every alignment and commit the output as a sanity-check artefact.

---

## 9. The seed and reproducibility

Pure progressive MAFFT (`--retree 2 --maxiterate 0`) is deterministic: same input, same version, same algorithm flag = byte-identical output. There is no random component. Adding iterative refinement (`--maxiterate 1000`) introduces a small amount of non-determinism via the iteration order (tie-breaking in the score comparison). MAFFT does not expose a seed flag; for strict reproducibility, pin the algorithm to a fully-deterministic mode (`--retree 2`) or accept that two L-INS-i runs may disagree on ~0.1% of columns.

Clustal Omega is deterministic by default. MUSCLE 5's PPP is deterministic; super5 is randomized and accepts a `-randseed N` flag — pin it.

Bootstrap resampling, NJ tie-breaking, and ML starting-tree generation are all randomized; we cover seed management for those in Lecture 2 and Lecture 3.

---

## 10. Common errors and how to spot them

- **MAFFT silently swaps `T` for `U` in RNA input.** If you feed RNA sequences without `--nuc`, MAFFT may convert. Pin `--nuc` and pre-convert your RNA to DNA (`U -> T`) explicitly.
- **Misordered FASTA records.** MAFFT preserves input order in the output. If your downstream code assumes alphabetical order (e.g. for a distance matrix), explicitly sort the records after reading.
- **Sequences of very different lengths.** A 100-bp read aligned with 1,500-bp 16S sequences will produce a column-mostly-gaps row. Drop too-short sequences before alignment.
- **Identical sequences.** MAFFT does not deduplicate; identical sequences produce zero-distance rows in the distance matrix and a zero-length branch in the NJ tree. Deduplicate before alignment.
- **Non-IUPAC characters.** A `*` (sometimes used for stop codons) or `?` (unknown) without `--anysymbol` may cause MAFFT to raise a warning or substitute silently. Either pre-clean the input or pass `--anysymbol`.

In Exercise 1 you run a `validate_input_fasta` function that catches the above before invoking MAFFT. The pattern matches the Week 8 `validate_input` function for VCFs.

---

## 11. What to take to the rest of the week

By the end of Lecture 1 you should be able to:

- Run MAFFT from Python on a FASTA, with the algorithm flag pinned, and capture the alignment in a `MultipleSeqAlignment` object.
- Compute per-column gap fractions, trim columns with > 50% gaps, and report the row count, column count, and conservation summary.
- State the differences between MAFFT, Clustal Omega, and MUSCLE 5 and pick the right one for an input size.
- Pin the MAFFT version (`mafft --version`), the algorithm flag (`--retree 2`), and the column-trim threshold (`--max-gap-fraction 0.5`) in any run-info JSON.

Lecture 2 takes the cleaned alignment and turns it into a distance matrix and a tree.

## References

- Carrillo H, Lipman D. "The multiple sequence alignment problem in biology." *SIAM Journal on Applied Mathematics* 48:1073 (1988). The original NP-hardness result for exact MSA. Behind a paywall in the *SIAM* archives; the result is standard and you do not need to read the paper to know the consequence.
- Feng DF, Doolittle RF. "Progressive sequence alignment as a prerequisite to correct phylogenetic trees." *Journal of Molecular Evolution* 25:351 (1987). The progressive alignment paper.
- Hogeweg P, Hesper B. "The alignment of sets of sequences and the construction of phyletic trees: an integrated method." *Journal of Molecular Evolution* 20:175 (1984). The independent progressive alignment paper.
- Thompson JD, Higgins DG, Gibson TJ. "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice." *Nucleic Acids Research* 22:4673 (1994). The ClustalW paper.
- Edgar RC. "MUSCLE: multiple sequence alignment with high accuracy and high throughput." *Nucleic Acids Research* 32:1792 (2004). The original MUSCLE paper, superseded by the MUSCLE 5 paper above.
- Katoh K, Standley DM. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." *Molecular Biology and Evolution* 30:772 (2013). Free at <https://academic.oup.com/mbe/article/30/4/772/1073398>.
- Sievers F et al. "Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega." *Molecular Systems Biology* 7:539 (2011). Free at <https://www.embopress.org/doi/full/10.1038/msb.2011.75>.
- Edgar RC. "Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny." *Nature Communications* 13:6968 (2022). Free at <https://www.nature.com/articles/s41467-022-34630-w>.
- Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T. "trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses." *Bioinformatics* 25:1972 (2009). Free at <https://academic.oup.com/bioinformatics/article/25/15/1972/213148>.
- Cock PJA et al. "Biopython: freely available Python tools for computational molecular biology and bioinformatics." *Bioinformatics* 25:1422 (2009). Free at <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>. Tutorial at <https://biopython.org/docs/latest/Tutorial/index.html>.
