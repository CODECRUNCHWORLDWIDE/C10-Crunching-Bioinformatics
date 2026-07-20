# Week 3 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (Biopython 1.83, NumPy 1.26), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Needleman & Wunsch (1970)** — the original global-alignment paper. Free PDF via *Journal of Molecular Biology*:
  <https://www.sciencedirect.com/science/article/pii/0022283670900574>
  Mirror (free PDF, often easier): <https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/needlemanandwunsch1970.pdf>
- **Smith & Waterman (1981)** — the local-alignment paper. Four pages, free PDF:
  <https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/SmithWaterman1981.pdf>
- **Henikoff & Henikoff (1992)** — the BLOSUM paper, *PNAS*. Free full text:
  <https://www.pnas.org/doi/10.1073/pnas.89.22.10915>
- **Biopython 1.83 `Bio.Align.PairwiseAligner` documentation** — the modern replacement for the deprecated `Bio.pairwise2`:
  <https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.PairwiseAligner>
- **Biopython Tutorial and Cookbook — Chapter 6 (Pairwise alignment)** — free PDF:
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>

## Substitution matrices (download the files)

- **NCBI BLOSUM and PAM matrices** — the canonical ASCII matrix files used by every aligner. Free:
  <https://ftp.ncbi.nih.gov/blast/matrices/>
  Files of interest: `BLOSUM62`, `BLOSUM45`, `BLOSUM80`, `PAM30`, `PAM70`, `PAM250`.
- **EBI mirror of the BLOSUM/PAM matrices** — same files, sometimes faster:
  <https://www.ebi.ac.uk/Tools/sss/fasta/matrices/>
- **Biopython's bundled matrices** — `Bio.Align.substitution_matrices`. List with:
  ```python
  from Bio.Align import substitution_matrices
  print(substitution_matrices.load())   # prints all available names
  m = substitution_matrices.load("BLOSUM62")
  ```

## Format references (have these open in tabs)

- **NCBI matrix file format** — the column-aligned ASCII layout used by `BLOSUM62` and friends:
  <https://en.wikipedia.org/wiki/BLOSUM>
- **PAM matrix family — Dayhoff et al. 1978** — the foundational paper, available via the NCBI Bookshelf and several university archives:
  <https://en.wikipedia.org/wiki/Point_accepted_mutation>
- **Affine gap penalty (Gotoh 1982)** — *Journal of Molecular Biology*, the O(mn) three-matrix algorithm:
  <https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf> (lecture notes citing Gotoh)

## Tools you will install this week

- **Biopython 1.83** — `pip install biopython==1.83`. The `Bio.Align.PairwiseAligner` API lives here.
- **NumPy 1.26** — `pip install numpy==1.26.4`. The whole week is NumPy arrays plus integer indexing.
- **matplotlib** — for the dot-plot at the end of the by-hand lecture. Any recent version is fine.

## Free books (chapter-level)

- **Biopython Tutorial and Cookbook** — chapter 6 covers pairwise alignment end to end:
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **Bioinformatics Algorithms: An Active Learning Approach (Compeau & Pevzner)** — chapter 5 is the canonical undergrad treatment of dynamic-programming alignment, with diagrams. Free chapters:
  <https://www.bioinformaticsalgorithms.org/>
- **Durbin, Eddy, Krogh, Mitchison — Biological Sequence Analysis (1998)** — chapter 2 is the rigorous version. The classic graduate text. Not free in print, but most university libraries have it; chapter 2 has been widely uploaded as scanned PDF.

## Worked examples online (free)

- **Rosalind — Edit Distance and Alignment problems** — free interactive problems that drill exactly this week's material:
  <https://rosalind.info/problems/locations/>
  Start with `EDIT`, `GLOB` (global alignment), and `LOCA` (local alignment).
- **EBI EMBOSS `needle` and `water`** — the canonical Unix command-line implementations. Free web tool to sanity-check your by-hand alignments:
  <https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle>
  <https://www.ebi.ac.uk/jdispatcher/psa/emboss_water>

## Open-source code to read this week

You can learn more from one hour reading other people's code than from three hours of tutorials. Pick one:

- **Biopython `Bio.Align._pairwisealigner.c`** — the C extension that backs `PairwiseAligner`. The inner loop is ~200 lines and very readable:
  <https://github.com/biopython/biopython/tree/master/Bio/Align>
- **scikit-bio `skbio.alignment`** — pure-Python implementations of NW and SW, used for teaching:
  <https://github.com/scikit-bio/scikit-bio/tree/main/skbio/alignment>
- **parasail** — vectorized SIMD implementations of Needleman-Wunsch, Smith-Waterman, and semi-global alignment. The SIMD tricks are out of scope this week, but the API surface is the same:
  <https://github.com/jeffdaily/parasail>

## Substitution-matrix cheat sheet

Keep this open while you work the exercises.

| Matrix | Best for | Built from | Citation |
|--------|----------|-----------|----------|
| **BLOSUM62** | Moderately diverged proteins (default for BLAST protein search) | BLOCKS database, sequences <62% identical clustered | Henikoff & Henikoff 1992 |
| **BLOSUM45** | Distantly related proteins | BLOCKS, <45% identical clusters | Henikoff & Henikoff 1992 |
| **BLOSUM80** | Closely related proteins | BLOCKS, <80% identical clusters | Henikoff & Henikoff 1992 |
| **PAM250** | Distantly related proteins (now superseded by BLOSUM45 in practice) | Manually curated families, extrapolated from PAM1 | Dayhoff et al. 1978 |
| **PAM30 / PAM70** | Short, very similar peptides (often used for short-query BLAST) | Same family, less extrapolation | Dayhoff et al. 1978 |
| **+1 / -1 (or +2 / -3)** | DNA, no evolutionary modeling | Toy/teaching | Various |

The intuition: PAM*N* models sequences that have undergone *N* point accepted mutations per 100 residues — so a high PAM number means distantly related. BLOSUM*N* clusters sequences *more than N% identical* before counting substitutions — so a *low* BLOSUM number means distantly related. The two scales run in opposite directions.

## Gap-penalty cheat sheet

| Scheme | Penalty for a gap of length `k` | When to use |
|--------|---------------------------------|-------------|
| Linear | `g * k` | Teaching, very rough alignments, when speed dominates correctness |
| Affine | `open + extend * (k - 1)` | Production protein and DNA alignment — the biological default |
| BLAST default (protein) | `open = 11`, `extend = 1` | The setting BLAST uses with BLOSUM62 |
| EMBOSS `needle` default | `open = 10`, `extend = 0.5` | EMBOSS's global aligner |
| Biopython default | `open = 0`, `extend = 0` | Biopython does **not** set a default — you must pick |

If a tool's documentation is silent on gap penalties, find out before trusting the output.

---

*If a link 404s, please open an issue so we can replace it.*
