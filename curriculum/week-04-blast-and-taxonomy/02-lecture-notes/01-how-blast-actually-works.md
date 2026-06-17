# Lecture 1 — How BLAST Actually Works

> **Duration:** ~3 hours of reading + paper-and-pencil + a brief Python sanity check.
> **Outcome:** You can describe the BLAST seed-and-extend heuristic in two paragraphs, name the four parameters that control it (`W`, `T`, `A`, `X`), state the Karlin-Altschul E-value formula and what each symbol means, and explain in one sentence each how BLAST 2.0's two-hit rule and X-drop extension trade optimality for speed.

If you only remember one thing from this lecture, remember this:

> **BLAST is a heuristic seed-and-extend algorithm with a statistical scoring framework on top.** It does *not* compute the optimal local alignment between a query and every entry in the database. It computes the optimal local alignment between the query and only those database entries that share at least one (for BLAST 1.0) or two (for BLAST 2.0) short, high-scoring **seed words** with the query. The seed step is sub-linear in the database size when the database is indexed correctly. The extension step is linear in the matched region. The E-value framework on top tells you, for any given final score, how often a hit that good would arise by pure chance.

This is the lecture where you build intuition for *why* BLAST is fast. Lecture 2 is where you learn to *use* it. Without Lecture 1 the parameters are magic numbers; with Lecture 1 they are knobs whose effect you can predict.

---

## 1. The problem Smith-Waterman cannot solve at scale

In Week 3 you implemented Smith-Waterman. The recurrence is:

```
H[i, j] = max(  0,
                H[i-1, j-1] + s(A[i], B[j]),
                H[i-1, j]   + g,
                H[i,   j-1] + g  )
```

The time complexity is `O(mn)` for a single pair of sequences of length `m` and `n`. The optimal local alignment is *guaranteed*: there is no input under which Smith-Waterman misses the highest-scoring local match.

Now consider the database-search version of the problem. The query is a single sequence of length `m` (typically 100–10,000 residues). The database is a collection of `N` sequences with total length `n_total` (the NCBI `nr` protein database has `N ≈ 6 × 10^8` and `n_total ≈ 2 × 10^11`). Smith-Waterman of the query against every database entry costs:

```
T_SW ≈ O(m × n_total)
     ≈ 10^4 × 2 × 10^11
     ≈ 2 × 10^15 cell updates per query.
```

At ~10^9 cell updates per second on a single CPU core (an optimistic estimate for a well-vectorized SW implementation; Farrar's SIMD SW from 2007 hits this), one query takes ~2 × 10^6 seconds, or about 23 days. The exact figures shift with the year and the database, but the order-of-magnitude conclusion does not: **doing Smith-Waterman against `nr` for a single query takes weeks on a single core**. That is the wall. BLAST exists to climb it.

The trick BLAST uses is not algorithmic in the "find a better recurrence" sense. It is *filtering*. Of those `N ≈ 6 × 10^8` database entries, the *vast majority share no significant similarity with the query at all*. If we can cheaply identify the ones that *might* be similar — and only run Smith-Waterman on those — we get a `10^4`-to-`10^6`-fold speedup at near-zero loss of sensitivity. That filter is the seed-and-extend heuristic.

---

## 2. The seed-and-extend pattern

The BLAST algorithm, distilled to four steps:

1. **Seed.** Find short exact (for `blastn`) or high-scoring (for `blastp`) word matches between the query and the database.
2. **Extend.** From each seed, extend the alignment outward in both directions without gaps, stopping when the running score drops by more than a threshold from its maximum seen so far (the X-drop rule).
3. **Score.** On the seed-extended HSPs (high-scoring segment pairs) that survive the seed-and-extend filter, run a gapped Smith-Waterman-style extension to produce the final alignment.
4. **Report.** Compute an E-value for each surviving alignment via the Karlin-Altschul formula and report hits below the user's E-value cutoff, sorted by E-value ascending.

Each step is a filter: the seeds drop most database entries, the X-drop extension drops most seeds, and the gapped extension drops most extensions. By the time you get to the final scoring, you are working on a few thousand candidate alignments, not a few hundred million.

### 2.1 Seeds for `blastp` (protein)

For protein search, BLAST splits the query into all overlapping **k-mers** of length `W = 3` (default). For each query k-mer, it enumerates the set of all 3-letter strings whose `BLOSUM62` similarity score against the query k-mer is at least the **neighborhood threshold** `T = 11` (default). This **neighborhood** typically contains 50–200 words per query k-mer. The union of all neighborhood words across the query is the BLAST **seed lookup table** for this query.

Example. Suppose the query contains the 3-mer `WFE`. The BLOSUM62 self-score of `WFE` is `s(W,W) + s(F,F) + s(E,E) = 11 + 6 + 5 = 22`. We enumerate all 3-mers `xyz` with `s(W,x) + s(F,y) + s(E,z) ≥ 11`. The set includes `WFE` itself (score 22), close variants like `WYE` (score `11 + 3 + 5 = 19`) and `WFD` (`11 + 6 + 2 = 19`), and a sprinkling of more divergent 3-mers that still happen to clear 11. Words like `AAA` (`s(W,A) + s(F,A) + s(E,A) = -3 + -2 + -1 = -6`) do not clear; they are not seeds.

BLAST then scans the database for any occurrence of any neighborhood word. The scan is `O(database length)` using a precomputed index (a hash table from words to positions), and is the dominant cost.

### 2.2 Seeds for `blastn` (nucleotide)

For nucleotide search, the alphabet is only 4 letters, so `W = 3` would give `4^3 = 64` possible words — too few to discriminate, every database position would be a seed. BLAST nucleotide search uses much larger word sizes:

- `blastn -task megablast`: `W = 28` (default; "exact 28-mer match required").
- `blastn -task blastn`: `W = 11` (more sensitive, slower).
- `blastn -task blastn-short`: `W = 7` (for queries < 30 bp like primers).
- `blastn -task dc-megablast`: discontiguous megablast, `W = 11` with a "spaced seed" template that allows mismatches at specific positions.

For nucleotide seeds, BLAST does *not* use a neighborhood. The seed must be an exact match. The reasoning: at `W = 28`, the chance of an exact 28-letter match by chance in a 10^11-base database is `10^11 / 4^28 ≈ 1.3 × 10^-6`, low enough that exact matching is sensitive *and* specific.

### 2.3 The two-hit rule (BLAST 2.0)

The 1990 BLAST algorithm extended every seed it found. In practice this still produced too many false-positive extensions to be efficient on large databases. Altschul et al. 1997 introduced the **two-hit rule**: only extend a seed if there is a **second** seed within `A = 40` residues of it, on the same diagonal of the alignment matrix, and the two seeds do not overlap. This costs almost no sensitivity (homologous regions typically contain *several* seeds, not just one) and cuts the number of triggered extensions by roughly 50x.

Geometric intuition: a "diagonal" of the dot-plot matrix is the set of cells `(i, j)` with `i - j = d` for some fixed `d`. Two seeds on the same diagonal at offsets `(i1, j1)` and `(i2, j2)` represent two short matching segments at the same relative shift between query and subject. If both segments matched by chance independently, the joint probability is squared — that is the source of the false-positive suppression.

### 2.4 The X-drop extension

Given a seed pair that has passed the two-hit filter, BLAST extends the alignment outward to the left and to the right *one residue at a time*, accumulating a running score. It does **not** allow gaps in this phase. The extension stops as soon as the running score drops by more than `X` below the maximum score seen so far on this extension. Default `X` values depend on the task; for `blastp` with BLOSUM62, `X ≈ 22` bits (so the running score is allowed to drop by ~22 before the extension is killed).

The X-drop rule is the second filter. Most extensions die fast: a seed that landed by chance has nothing on either side and the score crashes within 5–10 residues. The extensions that survive are the **HSPs** — high-scoring segment pairs — and they are passed to the final gapped extension.

### 2.5 Gapped extension

Once an HSP has survived the X-drop ungapped extension, BLAST runs a *gapped* extension on it — essentially a windowed Smith-Waterman around the HSP, allowing affine gaps. This is the only Smith-Waterman computation in the pipeline, and it runs on a tiny fraction of the database. Defaults: `gap_open = 11`, `gap_extend = 1` for protein with BLOSUM62; `gap_open = 5`, `gap_extend = 2` for `blastn`.

The final alignment that BLAST reports is the result of the gapped extension. Its score is what gets fed into the K-A formula to compute the E-value.

---

## 3. Karlin-Altschul statistics: what an E-value is

The 1990 BLAST paper introduced the seed-and-extend heuristic. The 1990 Karlin-Altschul paper (same year, same lab) introduced the statistical framework for interpreting the resulting scores. The two papers are paired and should be read together.

### 3.1 The null hypothesis

The K-A null model: the query and the database are independent sequences of i.i.d. residues drawn from the background distribution of amino acid frequencies (or nucleotide frequencies). Under this model, what is the distribution of *the maximum local alignment score* between the query and the database?

The answer (Karlin & Altschul 1990, theorem 1): under the null, the maximum local alignment score follows a **Gumbel extreme-value distribution** with parameters `K` and `λ`:

```
P(S >= x) ≈ 1 - exp( -K · m · n · exp(-λx) )
```

where `m` is the query length, `n` is the database length (in residues — for protein) or some effective length adjustment for nucleotide, and `K`, `λ` are constants that depend on the scoring system (substitution matrix and gap penalty).

For large `x` this simplifies to:

```
P(S >= x) ≈ K · m · n · exp(-λx)
```

This is the probability that a hit with score `>= x` arises *somewhere* in a search of size `m × n` under the null. Multiply through by the database size and you get the **expected number** of such hits — the **E-value**:

```
E(x) = K · m · n · exp(-λ · x)
```

That is the canonical Karlin-Altschul E-value formula. Memorize it.

### 3.2 Reading an E-value correctly

An E-value of `1e-50` for a hit at score `S = 312` in a search against `nr` means: *under the null hypothesis that the query and database are random sequences of i.i.d. residues, the expected number of hits at score 312 or better in a database this size is `1e-50` — effectively zero, so the observed hit is not chance.*

An E-value of `5` for a hit at score `S = 27` means: *we would expect about 5 hits at score 27 or better even if there were no real homology in the database*. So a single such hit is not informative on its own.

Note the asymmetry. **A small E-value is evidence against the null** (the hit is unlikely under chance). A small E-value is *not* "the probability the hit is wrong" — that is a Bayesian claim that requires a prior, which the E-value does not give you. This nuance routinely trips up biology students learning BLAST for the first time. Re-read the previous paragraph until it sticks.

### 3.3 Bit scores

The raw alignment score `S` is not directly comparable across searches with different scoring schemes — BLOSUM62 and PAM250 are on different scales — or across databases of different sizes. The **bit score** normalizes:

```
bit_score = (λ · S − ln K) / ln 2
```

The bit score is the negative log-base-2 of the chance probability. It is on a scale where every additional bit halves the chance probability. A bit score of 30 means the hit is ~10^-9 likely by chance (`2^-30 ≈ 9.3 × 10^-10`). A bit score of 80 means ~10^-24 likely. These numbers do not depend on the database size, only on the score and the scoring system. So if you want to ask "how good is this alignment biologically?" use the bit score. If you want to ask "how many hits this good would I see in this database under the null?" use the E-value (which combines the bit score with the database size).

### 3.4 The effective database length

You will sometimes see BLAST report `Effective database length` and `Effective query length`. These are corrections for edge effects: a hit cannot start in the last `m - 1` positions of a database sequence (it would not have room to extend), so the *effective* search space is slightly smaller than `m × n_total`. K-A's `m` and `n` in the formula above are typically the effective lengths, not the raw lengths. The correction is small for long sequences and you can ignore it for the lecture; it matters when you are computing E-values by hand to debug a BLAST run.

---

## 4. Seed-and-extend vs Smith-Waterman: the tradeoff

You should now be able to fill in the following table without looking. If you cannot, re-read sections 1–3.

| Property | Smith-Waterman | BLAST |
|----------|----------------|-------|
| Algorithm class | Optimal dynamic programming | Heuristic seed-and-extend |
| Time complexity per pair | `O(mn)` | Effectively `O(m + hits)` |
| Time complexity per query against a database of size `N` | `O(N · m · avg_seq_len)` | Roughly `O(N · m / k)` for some `k` ≫ 1 dependent on seed parameters |
| Sensitivity | Guaranteed optimal | High but not guaranteed; can miss low-similarity matches |
| Specificity | Same | Same — final scoring is Smith-Waterman-equivalent on the surviving pairs |
| Tunable parameters | Substitution matrix, gap penalties | Above, plus `W`, `T`, `A`, `X` |
| Statistical scoring | None built in | Karlin-Altschul E-value and bit score |
| When to use | Two specific sequences; reference-grade alignment | Database search; you need an answer in seconds, not weeks |

The fundamental tradeoff: BLAST gives up the *guarantee* that it will find every homologous match. In exchange, it runs ~10^4 to ~10^6 times faster. For nearly every database-search use case, that is the right trade. For two specific reference sequences you care deeply about (the SARS-CoV-2 spike pair in the Week 3 mini-project, or a curated pair of orthologs in a phylogenetic study), use Smith-Waterman.

### 4.1 What BLAST can miss

The honest answer: BLAST's seed step requires *at least one short, high-scoring word match* between query and subject (BLAST 1) or *at least two such matches close together on the same diagonal* (BLAST 2). For very distantly related proteins (< 20% identity), there may be no such word. Such homologies — when they exist — are typically found by structural alignment (DALI, FoldSeek) or PSI-BLAST iterative profile search, not by single-pass BLAST. The 20% identity floor is sometimes called the "twilight zone" of sequence homology.

For closer homologs (> 30% identity, and certainly > 50% identity), BLAST misses essentially nothing of practical interest. The seed-and-extend heuristic was tuned with care, and its sensitivity is excellent in the regime where most working bioinformaticians actually operate.

### 4.2 The parameter knobs

Here is what to know about each parameter, in the order they affect the algorithm:

- **Word size `W`** — Lower `W` means more seeds, slower, more sensitive. Higher `W` means fewer seeds, faster, less sensitive. For `blastn`, going from `W = 28` (megablast, optimized for highly similar sequences) to `W = 11` (default `blastn`) increases run time by ~5–10x but extends sensitivity into the 70–90% identity range. For `blastp`, `W = 3` is essentially fixed — going lower makes the neighborhood explode and going higher loses too much sensitivity.
- **Neighborhood threshold `T`** (protein only) — Lower `T` means a larger neighborhood per k-mer, more seeds, slower, more sensitive. Default `T = 11` for BLOSUM62 is a good compromise. Reducing to `T = 9` is the "high sensitivity" setting.
- **Two-hit window `A`** — Lower `A` means seeds must be very close to trigger extension; higher `A` (default 40) is more permissive. Reducing `A` is a way to speed up the algorithm at modest sensitivity cost.
- **X-drop threshold `X`** — Lower `X` (more aggressive cutoff) means extensions die faster, fewer HSPs, faster; higher `X` means more permissive extension, more HSPs, slower. Default `X = 22` bits for `blastp`.
- **E-value cutoff `-evalue`** — The reporting threshold. *Does not affect the search itself*, only which results are shown. Default `1e-3` for command-line BLAST+ and `10` for the web interface.

The week's homework problem 4 asks you to vary `W` for `blastn` on a real query and observe the sensitivity-vs-time tradeoff. Do it. The intuition is worth more than the numbers.

---

## 5. A worked sanity check in Python

You will not implement BLAST from scratch this week — that is a full semester project. But you should *run* it end to end on a tiny example so the abstract algorithm becomes concrete. Lecture 2 has the full how-to. For now, the smallest possible Python sanity check using only `Bio.Blast.NCBIWWW`:

```python
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "you@example.com"

# Fetch the human BRCA1 mRNA reference sequence (~7.2 kb).
handle = Entrez.efetch(db="nuccore", id="NM_007294.4", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
print(f"Query: {record.id}, length: {len(record.seq)} bp")

# Submit to NCBI BLAST against refseq_rna with megablast defaults.
result_handle = NCBIWWW.qblast(
    program="blastn",
    database="refseq_rna",
    sequence=str(record.seq),
    expect=1e-50,
    hitlist_size=5,
    megablast=True,
)

# Parse the XML.
blast_record = NCBIXML.read(result_handle)
print(f"Top hits for {record.id}:")
for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    print(
        f"  {alignment.accession}  {alignment.title[:60]:60s}"
        f"  E={hsp.expect:.2e}  bit={hsp.bits:.0f}"
        f"  id={hsp.identities}/{hsp.align_length}"
    )
```

Expected output (your numbers will vary slightly if the NCBI refseq_rna database has been updated):

```
Query: NM_007294.4, length: 7224 bp
Top hits for NM_007294.4:
  NM_007294    Homo sapiens BRCA1 DNA repair associated  E=0.00e+00  bit=13340  id=7224/7224
  NM_007297    Homo sapiens BRCA1 DNA repair associated  E=0.00e+00  bit=11900  id=7106/7224
  ...
```

The top hit is the query against itself — accession `NM_007294`, identity `7224/7224`, E-value zero (literally; the score is so high the K-A formula underflows). The next hits are alternative splice variants of the same gene. A perfectly-behaved sanity check.

Once you have run this once, two things should be clear:

1. **The Biopython BLAST API is small.** Three functions (`qblast` to submit, `read`/`parse` to parse, `NCBIXML.BlastRecord` to navigate). Get comfortable with them.
2. **NCBI BLAST takes 10–60 seconds even for a small query.** Plan for the wait. Cache results to disk.

If your code fails with a network or rate-limit error, that is the wrong kind of teaching moment but a normal kind of operational moment. Set `Entrez.email`, retry with backoff, and try again at 3 a.m. Eastern time when the queue is empty.

---

## 6. Common misconceptions

A short list of "things that seem right but are not":

- **"BLAST finds the optimal alignment."** No. BLAST finds *highly probable* optimal-or-near-optimal alignments via a heuristic seed-and-extend filter. The optimal alignment under Smith-Waterman is guaranteed; BLAST's is not. In practice, for sequences with > 30% identity, BLAST misses essentially nothing. For < 20% identity, BLAST can miss real homology that PSI-BLAST or structural alignment would find.
- **"A low E-value means the hit is real."** A low E-value is evidence against the null hypothesis of chance similarity. It is *not* a posterior probability and is *not* a guarantee of biological homology. A low E-value to a *contaminated* or *misannotated* database entry is still a low E-value — it just means your hit is to a contaminated entry. Always check what the subject actually is.
- **"Bit score and E-value mean the same thing."** Bit score is a database-independent normalization of the raw score. E-value is the bit score combined with the database size to give an expected number of chance hits. Bit score = "how good is this alignment intrinsically?"; E-value = "how surprised should I be to find a hit this good in *this* database?"
- **"`blastn` with default settings is always the right starting point."** Default `blastn -task megablast` is optimized for highly similar sequences (> 95% identity). For more diverged sequences, switch to `-task blastn` (`W = 11`) or `-task dc-megablast`. Read the BLAST+ user manual before trusting any default.
- **"BLAST against `nr` is always more sensitive than BLAST against a curated database."** It is more *exhaustive*, but `nr` contains many low-quality and redundant entries that can dilute the signal. For organism identification, querying the much smaller curated `16S_ribosomal_RNA` database is both faster and produces more interpretable hits.

If any of these surprised you, re-read sections 2 and 3.

---

## 7. Where this lecture lands you for Lecture 2

After this lecture you should be able to:

- Describe the BLAST seed-and-extend pipeline in four steps (seed, extend, score, report).
- Name the four key parameters (`W`, `T`, `A`, `X`) and what each controls.
- Write down the K-A E-value formula and state in one sentence what each symbol means.
- Distinguish "raw score", "bit score", and "E-value" — including which depends on database size and which does not.
- Explain in one sentence what BLAST gives up relative to Smith-Waterman, and what it gains.

Lecture 2 takes this conceptual material and turns it into the day-to-day operating manual: how to run `blastn`/`blastp`/`tblastn` from the command line, how to build a local database with `makeblastdb`, how to submit a query to NCBI via `Bio.Blast.NCBIWWW`, and how to parse the resulting XML or tabular output into a pandas DataFrame ready to feed into your taxonomy classifier on Thursday.

---

## Self-check questions

Before you move on, answer these without looking. If you cannot answer one, re-read the relevant section.

1. What are the four steps of the BLAST seed-and-extend pipeline? (§2)
2. For `blastp`, what is `W`? What is `T`? Why do they have different values than the `blastn` defaults? (§2.1, §2.2)
3. State the K-A E-value formula. (§3.1)
4. A hit has E-value `4e-87`. Is the hit "real"? What is the precise statistical claim? (§3.2)
5. A query of 1,000 bp is searched against two databases — `refseq_rna` (~10^9 bp) and `nt` (~10^11 bp). The same alignment in both cases gives the same raw score, same bit score, but the E-value in `nt` is ~100x larger. Why? (§3.1, §3.3)
6. What is the two-hit rule, and what does it gain? (§2.3)
7. Name two reasons BLAST can miss a real homologous match. (§4.1)
8. Why does `blastn` use `W = 28` (megablast) by default rather than `W = 3`? (§2.2)
9. Is the E-value cutoff a search parameter or a reporting parameter? Explain the difference. (§4.2)
10. Under what biological scenario would you reach for Smith-Waterman over BLAST despite BLAST being faster? (§4)

Answers are not provided. If you struggle, the answers are in the section references; do the work.

---

## Further reading

- Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. *Journal of Molecular Biology*, 215(3), 403–410.
- Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., & Lipman, D. J. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research*, 25(17), 3389–3402.
- Karlin, S., & Altschul, S. F. (1990). Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. *PNAS*, 87(6), 2264–2268.
- NCBI BLAST+ user manual: <https://www.ncbi.nlm.nih.gov/books/NBK279690/>.
- Biopython tutorial, Chapter 7: <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>.

---

*Continue to [Lecture 2 — Running BLAST locally and via NCBI](./02-running-blast-locally-and-via-ncbi.md) once you have answered the self-check questions.*
