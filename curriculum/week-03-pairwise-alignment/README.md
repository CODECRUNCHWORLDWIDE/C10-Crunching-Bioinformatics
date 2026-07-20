# Week 3 — Pairwise Alignment

In Week 2 you learned to parse, inspect, and quality-control sequence files. In Week 3 we **do something with the sequences**. The single most-used operation in all of bioinformatics is **pairwise alignment**: given two sequences `A` and `B`, find the arrangement of matches, mismatches, and gaps that maximizes a biologically meaningful score. By Friday of Week 3 you will be able to derive Needleman-Wunsch on paper for a 5x5 matrix, implement both Needleman-Wunsch and Smith-Waterman in pure NumPy, choose between BLOSUM62 and PAM250 for a protein alignment, and explain why an affine gap penalty exists and when to reach for it.

The other half of the week is **scoring**. Alignment is not a string-edit problem; it is an *optimization* problem over a scoring function, and the scoring function encodes biology. We will look at how Margaret Dayhoff built the original PAM matrices from manually curated protein families in the late 1970s, why Henikoff & Henikoff's BLOSUM matrices (1992) replaced them for distant homology, and why the gap-open / gap-extend split (the **affine gap penalty**) more accurately models the biology of insertions and deletions than a single linear penalty. The mini-project is a from-scratch NumPy Smith-Waterman benchmarked against Biopython's `PairwiseAligner` on a small DNA case study.

## Learning objectives

By the end of this week, you will be able to:

- **Fill in** the Needleman-Wunsch dynamic-programming matrix by hand for a 5-residue x 5-residue example, including the initial gap row and column, and read off the optimal global alignment by traceback.
- **Fill in** the Smith-Waterman matrix by hand for the same example, including the zero-floor rule, and read off the optimal local alignment.
- **Implement** Needleman-Wunsch in pure NumPy with linear gap penalty, returning both the score and the aligned strings.
- **Implement** Smith-Waterman in pure NumPy with linear gap penalty, returning the local alignment plus its start and end coordinates in both sequences.
- **Choose** between BLOSUM62 and PAM250 for a protein alignment based on the expected evolutionary distance of the pair, and justify the choice in one sentence.
- **Choose** a linear vs affine gap penalty for a given biological context (DNA vs protein, close vs distant homology).
- **Use** Biopython's `Bio.Align.PairwiseAligner` to verify your NumPy implementation, and explain at least one case where the two implementations should differ at the byte level but agree on score.
- **Benchmark** your NumPy implementation against Biopython on a pair of ~500 bp sequences and report wall-clock and peak memory in the C10 reproducibility-receipt format.

## Prerequisites

This week assumes Weeks 1 and 2 are **done and committed**. Specifically:

- You can parse a FASTA file with `Bio.SeqIO.parse` and pull two sequences into Python as strings (Week 2 Exercise 1).
- You are comfortable enough with NumPy to build a 2D array, index it as `M[i, j]`, and walk a traceback by integer offsets. If you have never written `np.zeros((m+1, n+1), dtype=np.int32)` before, do the NumPy quickstart now.
- You have Python 3.11+, Biopython 1.83, and NumPy 1.26 installed. We will pin those exact versions in `env.yml` and in the mini-project's reproducibility receipt.

You do not need biology beyond the central dogma. You do need a quiet hour to fill in a 5x5 matrix by hand; this is the lecture that everyone skips and then re-learns when their NumPy traceback comes out wrong.

## Topics covered

- The pairwise-alignment problem statement: given `A`, `B`, a substitution matrix `s(a, b)`, and a gap penalty `g`, find the alignment with maximum score
- The Needleman-Wunsch recurrence (1970) for global alignment — initialization of row 0 and column 0 with gap-penalty multiples, the `max` of three predecessors, the traceback from `M[m, n]` back to `M[0, 0]`
- The Smith-Waterman recurrence (1981) for local alignment — the zero-floor rule, the traceback from the global maximum cell back to the first zero, why local alignment is the right default for searching a protein against a database
- Substitution matrices: **PAM250** (Dayhoff, 1978) for closely related proteins; **BLOSUM62** (Henikoff & Henikoff, 1992) for moderate divergence; **BLOSUM45** for distant homology; **BLOSUM80** for close homology; the simple `+1 / -1` match-mismatch scheme for DNA
- Gap penalties: linear (`g * k` for a gap of length `k`), affine (`open + extend * (k - 1)`), and why affine matches the biology of single-event indels of arbitrary length
- Time and space complexity: O(mn) time and O(mn) space for the naive implementation; the Hirschberg trick for O(min(m, n)) space; why we do not implement Hirschberg this week
- Biopython 1.83 `Bio.Align.PairwiseAligner` — the modern replacement for the deprecated `Bio.pairwise2`, what its `mode`, `substitution_matrix`, `open_gap_score`, and `extend_gap_score` attributes do
- Edge cases: zero-length input, sequences of very different length, ambiguity codes (`N`, `X`), case sensitivity

## Weekly schedule

The schedule below adds up to approximately **36 hours**. Treat it as a target. The by-hand matrix on Monday is the single highest-leverage hour of the week — do not skip it.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | Needleman-Wunsch by hand, the recurrence           |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | Smith-Waterman by hand, local vs global            |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Wednesday | Substitution matrices: PAM, BLOSUM                 |    1h    |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     6.5h    |
| Thursday  | NumPy implementations, traceback                   |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Affine gaps, Biopython PairwiseAligner             |    0h    |    1h     |     1h     |    0.5h   |   1h     |     2h       |    0h      |     5.5h    |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **7.5h**  | **3h**     | **3h**    | **6h**   | **8h**       | **2h**     | **35.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | Needleman-Wunsch and Smith-Waterman papers, BLOSUM and PAM matrices, Biopython PairwiseAligner docs |
| [lecture-notes/01-needleman-wunsch-and-smith-waterman-by-hand.md](./lecture-notes/01-needleman-wunsch-and-smith-waterman-by-hand.md) | The two recurrences, derived and walked through on a 5x5 worked example |
| [lecture-notes/02-substitution-matrices-and-gap-penalties.md](./lecture-notes/02-substitution-matrices-and-gap-penalties.md) | PAM250, BLOSUM62, BLOSUM45/80, linear vs affine gaps, how to choose |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-needleman-wunsch-numpy.py](./exercises/exercise-01-needleman-wunsch-numpy.py) | Implement Needleman-Wunsch in pure NumPy |
| [exercises/exercise-02-smith-waterman-numpy.py](./exercises/exercise-02-smith-waterman-numpy.py) | Implement Smith-Waterman in pure NumPy |
| [exercises/exercise-03-biopython-compare.py](./exercises/exercise-03-biopython-compare.py) | Verify your implementations against `Bio.Align.PairwiseAligner` |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-affine-gap-penalty.md](./challenges/challenge-01-affine-gap-penalty.md) | Extend Needleman-Wunsch to affine gaps (Gotoh's three-matrix recurrence) |
| [quiz.md](./quiz.md) | 10 multiple-choice questions on alignment, scoring, and complexity |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | Pure-NumPy Smith-Waterman benchmarked vs Biopython, on a small DNA case study |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("Biopython 1.83," "NumPy 1.26"). We cite matrices by name and citation ("BLOSUM62, Henikoff & Henikoff 1992"). We say "the local alignment score is 47 under BLOSUM62 with gap-open -11 / gap-extend -1," not "the alignment looks reasonable." An alignment score is a number on a known scale, not a vibe. If your benchmark report uses the words "fast" or "slow" without seconds attached, you have not written a benchmark report yet.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Open `Bio/Align/_pairwisealigner.c` in the Biopython source. Read the inner loop. Notice how it folds the substitution matrix lookup into a single index operation and how it batches the traceback array as a packed array of small integers. Compare to your NumPy version.
- Implement Hirschberg's algorithm — the divide-and-conquer trick that reduces Needleman-Wunsch space from O(mn) to O(min(m, n)) while preserving the alignment, not just the score. Hirschberg, *CACM* 1975. You will need it the day you try to align two whole chromosomes.
- Read the Smith-Waterman 1981 paper end to end. It is four pages. The "obviousness in retrospect" of the zero-floor rule is the kind of insight that wins a long-cited paper.
- Build the BLOSUM62 matrix from scratch. Download the BLOCKS database (still hosted on Henikoff's lab page), tabulate substitution frequencies in conserved blocks, take the log-odds ratio, round, and compare to the canonical BLOSUM62. The off-by-one errors in your reconstruction are where the biology lives.

## Up next

Continue to [Week 4 — BLAST and taxonomy](../week-04/) once you have pushed your mini-project to GitHub.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
