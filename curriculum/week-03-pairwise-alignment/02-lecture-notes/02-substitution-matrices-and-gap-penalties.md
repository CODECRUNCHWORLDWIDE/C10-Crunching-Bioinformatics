# Lecture 2 — Substitution Matrices and Gap Penalties

> **Duration:** ~2 hours of reading + hands-on.
> **Outcome:** You can choose between BLOSUM62, BLOSUM45, BLOSUM80, and PAM250 for a given pair of proteins and justify the choice. You can choose between a linear and an affine gap penalty and justify that. You can read the standard NCBI matrix-file format off disk and load a matrix into a NumPy lookup. You can articulate, in one sentence each, why the PAM and BLOSUM scales run in opposite directions.

The Needleman-Wunsch and Smith-Waterman recurrences from Lecture 1 are content-free without a **scoring function**. The scoring function is two things: a **substitution matrix** `s(a, b)` and a **gap penalty** `g` (or, for affine gaps, a pair `(open, extend)`). This lecture is about choosing both.

If you only remember one thing:

> **The substitution matrix encodes evolutionary distance.** A protein pair that diverged 100 million years ago needs a different matrix than a pair that diverged 10 million years ago. BLOSUM62 is the safe default for "moderate" divergence (the BLAST default). BLOSUM45 is for distant pairs. BLOSUM80 is for close pairs. PAM250 is the older alternative, still seen in legacy tools. Pick deliberately, document the choice, and never silently accept a tool default.

---

## 1. Why we need a substitution matrix

In DNA alignment, the simplest scoring scheme is `match = +1, mismatch = -1`. There is no biological reason to score `A vs G` differently from `A vs T` — both are point substitutions, equally "wrong." A bit more refinement: distinguish **transitions** (purine ↔ purine `A ↔ G`, or pyrimidine ↔ pyrimidine `C ↔ T`) from **transversions** (purine ↔ pyrimidine). Transitions occur about twice as often as transversions in most organisms, so a slightly refined DNA matrix might score `match = +1, transition = -1, transversion = -2`. This is the kind of thing a phylogenetics package (Week 9) will hand you, but for the work this week we will stick with `+1 / -1`.

Protein alignment is fundamentally different. The 20 amino acids vary on multiple axes — hydrophobicity, charge, size, aromaticity, structural role. A substitution from leucine (`L`) to isoleucine (`I`) is biologically nearly silent: both are aliphatic, both pack into the hydrophobic core of a protein, the protein's fold is preserved. A substitution from leucine to glutamate (`E`, negatively charged) is biologically catastrophic: a charged residue cannot occupy a hydrophobic core. The scoring matrix should reflect this — `L ↔ I` should score *almost as well as* a match; `L ↔ E` should score worse than a random mismatch.

The two great families of protein substitution matrices that encode this biology are **PAM** (Dayhoff et al., 1978) and **BLOSUM** (Henikoff & Henikoff, 1992).

---

## 2. PAM matrices (Dayhoff, 1978)

Margaret Dayhoff and colleagues at the National Biomedical Research Foundation built the **Point Accepted Mutation** (PAM) matrices from the manually curated *Atlas of Protein Sequence and Structure*. The construction:

1. Manually align closely related protein families (cytochrome c, hemoglobin, etc.).
2. Build a phylogenetic tree for each family.
3. Count, on each branch of each tree, the substitutions that were **accepted** by natural selection (i.e., observed and not eliminated).
4. Normalize: PAM1 is the matrix you get for an evolutionary distance such that 1% of residues have been substituted.
5. Extrapolate: PAM*N* = PAM1 raised to the power *N*. PAM250 represents the distance at which 250 substitutions per 100 residues have accumulated (with back-mutations, many positions are substituted more than once, so PAM250 is *not* simply "75% identity").

**PAM250** is the historically dominant choice for distantly related protein pairs. PAM30 and PAM70 are tighter matrices for short, similar peptides — they were the default for BLAST short queries for many years.

> **Critique.** The PAM matrices' extrapolation step assumes the substitution process is Markovian and that the closely-related-family substitution rates extrapolate to distant relationships. Neither is fully true. By the early 1990s, the field had data on a much larger sample of distantly related proteins, and PAM250 was visibly suboptimal for the cases it was designed for. Enter BLOSUM.

---

## 3. BLOSUM matrices (Henikoff & Henikoff, 1992)

Steven and Jorja Henikoff, *PNAS* 89:10915, 1992. They built the **BLOck SUbstitution Matrices** from the BLOCKS database — a curated collection of ungapped, highly conserved local protein-alignment blocks. The construction:

1. Take the BLOCKS database (thousands of ungapped multiple alignments of conserved protein regions).
2. Cluster sequences within each block by percent identity at threshold *N* — sequences more than *N%* identical to each other are merged into a single "cluster" so they don't dominate the count.
3. Across all blocks, count substitution frequencies in the clustered alignments.
4. The substitution matrix entry `s(a, b)` is the log-odds ratio of `observed frequency of a↔b` to `expected frequency assuming independence`, multiplied by `2/log(2)` (a scaling factor) and rounded to an integer.

The cluster threshold gives the matrix its name:

- **BLOSUM80** clusters at 80% identity — sequences are weighted toward close pairs. Best for **closely related** protein alignment.
- **BLOSUM62** clusters at 62% identity — the BLAST default. Best for **moderate divergence**.
- **BLOSUM45** clusters at 45% identity — sequences in clusters are quite different. Best for **distant homology**.

> **The scale runs opposite to PAM.** A high BLOSUM number means close, a low BLOSUM number means distant. A high PAM number means distant, a low PAM number means close. Memorize this — it will trip you up in code-review.

The rough mapping (not exact, the construction methodologies differ): BLOSUM45 ≈ PAM250, BLOSUM62 ≈ PAM160, BLOSUM80 ≈ PAM120, BLOSUM90 ≈ PAM100. NCBI's documentation will sometimes use these mappings when discussing BLAST defaults.

### BLOSUM62 in detail

This is the matrix you will see most often. Some entries to internalize:

| Pair | Score | Biological reading |
|------|------:|---------------------|
| `W vs W` | +11 | Tryptophan is rare and structurally distinctive — matches are very informative |
| `C vs C` | +9 | Cysteines often form disulfide bridges — matches indicate strong structural conservation |
| `A vs A` | +4 | Alanine is small, common, and freely substitutable — matches are unremarkable |
| `L vs I` | +2 | Both aliphatic; substitution is nearly silent |
| `K vs R` | +2 | Both positively charged; substitution is conservative |
| `D vs E` | +2 | Both negatively charged; substitution is conservative |
| `L vs E` | -3 | Hydrophobic ↔ charged; substitution is disruptive |
| `W vs A` | -3 | Bulky aromatic ↔ small aliphatic; substitution is disruptive |
| `*` vs anything | varies | `*` is the stop codon; some matrices include it, some do not |

A few quick rules of thumb:

- **Aromatics (W, Y, F)** score very well against each other (+2 to +3) and poorly against everything else.
- **Aliphatics (L, I, V, M)** score well against each other (+1 to +2).
- **Cysteine** is its own residue family — a `C ↔ C` match is the highest off-diagonal-of-the-class score, but `C` against anything else is around -1 to -3.
- **Charge** matters: `K ↔ R` and `D ↔ E` are mild, `K ↔ D` is harsh.

This is the structure that PAM and BLOSUM share with the actual biochemistry. Use it to sanity-check your matrix loader — if `s('W', 'W')` is not greater than 8, you have loaded the matrix wrong.

---

## 4. Choosing a matrix in practice

This is the workflow that production pipelines follow, and the one we will use this week:

1. **Is it DNA or protein?** DNA: use `+1 / -1` or `+1 / -2 / -2` (transition-transversion aware) for teaching; use minimap2's defaults (`+1 / -19 / -39 ...`) for production read alignment in Week 5. Protein: continue to step 2.
2. **Do you know the divergence?** If the pair are obvious orthologs from closely related species (human vs chimp), use **BLOSUM80** or PAM30. If they are distantly homologous (human vs *E. coli*), use **BLOSUM45** or PAM250.
3. **Don't know the divergence?** Use **BLOSUM62**. It is the BLAST default for a reason — it is the best single-matrix compromise across the range of biologically interesting pairs.
4. **Document the choice.** "Aligned with BLOSUM62, affine gap open `-11` / extend `-1`" is the methods-section line. "Aligned with default settings" is not.

Two operational tips:

- **Always check the matrix you actually got loaded.** Biopython's `substitution_matrices.load("BLOSUM62")` returns a `SubstitutionMatrix` object you can index as `m['W', 'W']`. Print a few entries and confirm they match the published values. The number of times a colleague has run BLAST with `BLOSUM62` but loaded `BLOSUM45` from a stale file is non-zero.
- **NCBI's matrix files** (`https://ftp.ncbi.nih.gov/blast/matrices/`) are the canonical ASCII source. Each is a small text file with a 20x20 grid (plus rows/columns for `B`, `Z`, `X`, `*` which encode ambiguity codes — average, stop codon). Read the file once with your eyes before trusting any parser.

---

## 5. Gap penalties — linear vs affine

A gap of length `k` corresponds to a single biological event: an **insertion** or a **deletion** of a contiguous stretch of `k` residues. The biology of indels is highly skewed — small indels are common, large indels are rare, but a large indel is *much* more likely than `k` independent single-residue indels.

The **linear** gap penalty `-g * k` does not capture this. Under a linear penalty, the cost of one 10-residue gap is identical to the cost of ten 1-residue gaps in different places. The alignment that the algorithm picks is therefore biased toward many small gaps when one large gap is more likely.

The **affine** gap penalty fixes this:

```
gap_cost(k) = open + extend * (k - 1)
```

where `open` is the cost of *starting* a gap and `extend` is the cost of *each additional residue* in the gap. With `open = -11` and `extend = -1` (the BLAST protein defaults with BLOSUM62), a gap of length 1 costs `-11`, a gap of length 2 costs `-12`, ..., a gap of length 10 costs `-20`. The "first residue is expensive, subsequent residues are cheap" structure matches the biology of a single indel event.

Production aligners — BLAST, BWA, minimap2, MAFFT, EMBOSS `needle`/`water` — all use affine gaps. The linear gap penalty survives mostly in teaching examples (including this week's exercises, until Friday).

### The Gotoh 1982 trick

Implementing affine gaps naively requires examining all `k` for each cell, blowing up the time complexity to O(mn(m + n)). Osamu Gotoh's 1982 paper showed how to keep the algorithm O(mn) by introducing **three matrices** instead of one:

- `M[i, j]` — best score ending in a match/mismatch at `(i, j)`.
- `X[i, j]` — best score ending in a gap in `A` (i.e., gap aligned to `B[j]`).
- `Y[i, j]` — best score ending in a gap in `B`.

Each matrix has its own recurrence; transitions between matrices encode "open a gap" vs "extend a gap." The end-of-alignment score is `max(M[m, n], X[m, n], Y[m, n])`. The challenge problem for this week (challenge 1) walks you through the Gotoh recurrence in detail.

---

## 6. A worked example: BLOSUM62 in action

Consider this pair of protein fragments:

```
A:  W F Q L R G F E
B:  W F K L R - F E
```

Score under BLOSUM62, with affine gap open `-11`, extend `-1`:

| Column | A | B | s(A, B) | Running |
|-------:|---|---|--------:|--------:|
| 1 | W | W | +11 | +11 |
| 2 | F | F | +6 | +17 |
| 3 | Q | K | +1 | +18 |
| 4 | L | L | +4 | +22 |
| 5 | R | R | +5 | +27 |
| 6 | G | - | -11 (gap open) | +16 |
| 7 | F | F | +6 | +22 |
| 8 | E | E | +5 | +27 |

The alignment scores `+27` under BLOSUM62 with a single 1-residue gap. Note how the conservative substitution `Q vs K` (both polar, K positively charged) still scores `+1` — better than the `-1` we would have given it with the toy DNA scheme.

If we had used a linear gap penalty of `-11`, the score would be identical (the gap is only 1 residue long, so `open == open + extend * 0`). The difference between linear and affine only appears for gaps of length ≥ 2.

---

## 7. Loading a substitution matrix in Python

Three ways, from most to least dependency:

### Biopython (the easy way)

```python
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")
print(blosum62['W', 'W'])    # 11.0
print(blosum62['L', 'I'])    # 2.0
print(blosum62['L', 'E'])    # -3.0
```

The matrix is a `Bio.Align.substitution_matrices.Array` — a NumPy ndarray subclass with letter indexing. You can call `.alphabet` to see the residues it covers (typically `'ARNDCQEGHILKMFPSTWYVBZX*'` — 20 amino acids + 4 ambiguity codes).

### Loading the NCBI matrix file directly

The file format is a column-aligned ASCII grid. The header row lists the alphabet, each subsequent row gives the scores for one residue against all others:

```
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
...
```

Parsing it yourself is a 15-line script (and it makes a good exercise — but Biopython's loader is faster and equally available).

### Building it from scratch

For the toy DNA case in the exercises, we build it in NumPy:

```python
import numpy as np

ALPHABET = "ACGT"
MATCH = 1
MISMATCH = -1

def build_dna_matrix(match=MATCH, mismatch=MISMATCH):
    n = len(ALPHABET)
    M = np.full((n, n), mismatch, dtype=np.int32)
    np.fill_diagonal(M, match)
    return M, {c: i for i, c in enumerate(ALPHABET)}

S, IDX = build_dna_matrix()
score = S[IDX['A'], IDX['A']]   # +1
score = S[IDX['A'], IDX['G']]   # -1
```

This is the form your NumPy implementations will use this week. The lookup `S[IDX[a], IDX[b]]` is `O(1)` and vectorizes cleanly if you precompute integer indices for whole sequences.

---

## 8. Common bugs in scoring

**Mismatching alphabets.** A protein sequence with an `X` (unknown residue) or a `*` (stop codon) handed to a matrix without those columns will crash, or worse, silently use a wrong entry. Biopython's BLOSUM62 includes `BZX*`. Your hand-built lookup likely does not. Either expand the alphabet or filter the input.

**Case sensitivity.** Most matrices are upper-case. If you forget to upper-case your input, `s('a', 'a')` will key-error. Cast at the boundary: `seq = str(record.seq).upper()`.

**Wrong matrix loaded.** BLOSUM50 and BLOSUM62 look similar in print; a typo in a config file is hard to spot. Always log the matrix name and a sanity check (`s('W', 'W')`).

**Linear vs affine confusion.** A tool that advertises "gap penalty -10" may mean linear `-10` per residue, or affine `open = -10, extend = -1`, or affine `open = -10, extend = -10` (which is degenerate to linear). Read the docs. EMBOSS `needle`: `open = 10, extend = 0.5` by default. BLAST protein: `open = 11, extend = 1`. Biopython 1.83: **no defaults**, you must set them.

**Sign convention.** Some literature writes gap penalties as positive numbers (and subtracts), others as negative numbers (and adds). Biopython's `PairwiseAligner` uses signed values added to the score — gaps must be **negative**. Get this wrong and your scores are nonsensically high.

---

## 9. Recap and next lecture

You should now be able to:

- Choose between PAM250, BLOSUM45, BLOSUM62, BLOSUM80 for a given protein pair, justified by expected divergence.
- Explain why the PAM and BLOSUM scales run in opposite directions.
- Sketch why an affine gap penalty matches the biology of single-event indels better than a linear penalty.
- State the BLAST-default gap-open and gap-extend values for protein alignment with BLOSUM62 (`-11`, `-1`).
- Load BLOSUM62 from Biopython, index it, and confirm `s('W', 'W') = 11`.
- Build a small DNA substitution matrix in NumPy and use it for `O(1)` lookups in a tight loop.
- Diagnose the four common scoring bugs (alphabet mismatch, case, wrong matrix, sign convention).

The exercises this week use the toy `+1 / -1` DNA scheme so you can hand-verify your NumPy code against the by-hand 5x5 example from Lecture 1. The challenge problem moves to affine gaps. The mini-project uses BLOSUM62 (via Biopython) on a real DNA pair for the comparison run.

In Week 4 we will use these same alignment ideas — but turned into a database search via BLAST. BLAST is, at its core, a heuristic Smith-Waterman: it indexes the database with short word seeds, extends each seed into a full SW alignment, and returns the highest-scoring hits with E-values. Everything you learned this week is what BLAST is doing under the hood.

---

*Return to the [Week 3 README](../00-overview.md), or jump to the [Exercises](../03-exercises/00-overview.md).*
