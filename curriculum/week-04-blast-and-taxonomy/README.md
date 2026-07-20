# Week 4 — BLAST and Taxonomy

In Week 3 you built a Smith-Waterman aligner from scratch in NumPy and benchmarked it against Biopython on a pair of ~3.8 kb SARS-CoV-2 spike sequences. That gave you a tool that finds the optimal local alignment between *two* sequences in `O(mn)` time. In Week 4 we ask the next question: **what do you do when one of the sequences is a 200-residue protein and the other is a 600-million-protein database?** `O(mn)` over the full database is `O(120 billion)` per query. That is the wall every sequence search hits, and that is the problem BLAST (Altschul, Gish, Miller, Myers, Lipman, 1990) solved by sacrificing optimality for the right kind of speed.

By Friday of Week 4 you will be able to run `blastn`, `blastp`, and `tblastn` from the command line on a local database you built yourself, run the same query via the NCBI web API with `Bio.Blast.NCBIWWW`, parse the resulting XML or tabular hit table with `Bio.Blast.NCBIXML`, read an E-value off a hit and explain in one sentence what null hypothesis it falsifies, and identify the genus of an unknown DNA sequence by taking the lowest-E-value hit's taxonomic lineage from the NCBI Taxonomy database. The mini-project is a BLAST-driven taxonomy classifier for ~20 unknown sequences, with a precision-recall accounting that you will defend in writing.

The other half of the week is **what BLAST gives up to be fast**. Smith-Waterman is *optimal*: it is guaranteed to find the highest-scoring local alignment under the chosen substitution matrix and gap penalty. BLAST is *heuristic*: it scans the database for short, high-scoring **seed words** (the "two-hit" rule for protein, the contiguous-word rule for nucleotide), extends those seeds outward without gaps until the score drops by a threshold, joins compatible seed hits, and only *then* runs a gapped Smith-Waterman-style extension on the surviving pairs. This is the **seed-and-extend** pattern, and once you have internalized it you will see it everywhere — in BWA for short-read alignment (Week 5), in minimap2 for long-read alignment (Week 5), in DIAMOND for protein-vs-protein at scale, in any aligner that ever has to work on a database of more than a few thousand sequences. Karlin & Altschul's 1990 statistical framework (the K-A formula, `E = K·m·n·exp(-λS)`) is the other half — it tells you, given a database size and a raw alignment score, how often you would expect a hit that good to arise *by chance* under a null model of i.i.d. random residues. The E-value is not "how similar these sequences are." It is the *expected number of false positives at this score or better* given the database you searched. Get that distinction right and the rest of the week falls into place.

## Learning objectives

By the end of this week, you will be able to:

- **Describe** the BLAST seed-and-extend heuristic in two paragraphs, naming the word size `W`, the neighborhood threshold `T`, the two-hit window `A`, the X-drop threshold `X`, and the role each plays in trading optimality for speed.
- **State** what an E-value is in one sentence (expected number of chance hits at this score or better, in a database of this size, under a random null model) and what a bit score is in one sentence (a normalized score independent of database size and substitution-matrix scale).
- **Run** `makeblastdb` to build a local nucleotide database from a FASTA file and `blastn` to query it, choosing word size and E-value cutoff explicitly.
- **Run** `blastp` against a local protein database and `tblastn` against a 6-frame translation of a nucleotide database, and explain which of the BLAST family members is appropriate for a given biological question.
- **Run** the same `blastn` query against NCBI's `nt` database via `Bio.Blast.NCBIWWW.qblast` and parse the returned XML with `Bio.Blast.NCBIXML.parse`, extracting hit titles, accessions, scores, E-values, and aligned coordinates.
- **Read** a BLAST tabular output (format `6` or `7`) into a pandas DataFrame, filter by E-value and percent identity, and produce a one-hit-per-query summary by minimum E-value.
- **Fetch** the NCBI Taxonomy lineage for a hit accession using `Bio.Entrez.efetch`, and identify the lowest common ancestor across the top-N hits for a query (the LCA classifier baseline).
- **Classify** ~20 unknown DNA sequences to genus or species using a top-hit-with-E-value-cutoff classifier, report classification confidence per query, and produce a confusion matrix against a ground-truth label set.
- **Identify** at least three failure modes of BLAST-based classification (low-complexity regions giving spurious hits, contaminated reference databases producing high-confidence wrong calls, paralog vs ortholog ambiguity) and the standard mitigation for each.

## Prerequisites

This week assumes Weeks 1, 2, and 3 are **done and committed**. Specifically:

- You can parse a FASTA file with `Bio.SeqIO.parse` and pull a list of sequences into Python (Week 2 Exercise 1).
- You can score a pairwise local alignment with `Bio.Align.PairwiseAligner` and read off the score, the aligned region, and the percent identity (Week 3 mini-project).
- You have a working `crunch-bio-portfolio-<yourhandle>/` repo with a `week-03/` directory committed. The Week 4 mini-project lives in `week-04/` alongside it.
- You have Python 3.11+, Biopython 1.83, and pandas installed. You will need to install BLAST+ (the NCBI command-line BLAST package) this week — `conda install -c bioconda blast=2.15` is the canonical path.

You do not need biology beyond "DNA is made of A/C/G/T, proteins are made of 20 amino acids, the central dogma runs DNA → RNA → protein." You do need patience for the first time you run a `blastn` against `nt` over the network — NCBI's web service rate-limits and queues, and a single query can take 30 seconds to several minutes depending on time of day. Plan around it.

## Topics covered

- The pairwise-search problem at scale: a 30,000-protein query bank against a 600M-protein database is `O(10^15)` Smith-Waterman cell updates; we need a heuristic.
- The seed-and-extend pattern: **seed** (find short exact or near-exact word matches between query and database), **extend** (grow the seed in both directions while score stays above a drop-off threshold), **score** (compute a final gapped local alignment on the surviving pairs), **report** (rank by E-value).
- Word size `W` and the BLAST neighborhood threshold `T`: for `blastp`, `W = 3` and a residue triple is a *seed* if its `BLOSUM62` similarity score against the query triple is at least `T = 11`; for `blastn` (megablast), `W = 28` and the seed must be exact.
- The two-hit rule (BLAST 2.0, Altschul et al. 1997): require *two non-overlapping seeds within distance A of each other on the same diagonal* before triggering extension, which cuts false-positive extensions by ~50x at minimal sensitivity cost.
- Karlin-Altschul statistics: the E-value formula `E = K·m·n·exp(-λS)` where `m` is the query length, `n` is the effective database length, `S` is the raw alignment score, and `K, λ` are parameters of the score matrix and the residue composition.
- Bit scores: `bit = (λS - ln K) / ln 2`, a normalized score on a scale where adding 1 bit halves the expected number of chance hits. Bit scores are comparable across databases of different sizes; raw scores are not.
- BLAST family members: `blastn` (nucleotide query, nucleotide DB), `blastp` (protein, protein), `blastx` (translated nucleotide query, protein DB), `tblastn` (protein query, translated nucleotide DB), `tblastx` (translated nucleotide query, translated nucleotide DB). Most of the week is `blastn` and `blastp`; `tblastn` appears in the mini-project for the gene-finding use case.
- The NCBI Entrez E-utilities: `efetch` (fetch a record by accession), `esearch` (run a database query), `elink` (follow cross-database links), `summary` (compact metadata). All accessible from `Bio.Entrez` with rate-limit handling baked in.
- The NCBI Taxonomy database and the lineage lookup: every GenBank record has a `taxon_id` field that resolves to a hierarchical lineage (species → genus → family → order → class → phylum → kingdom → domain). Cross-walks from accession to lineage are the foundation of every BLAST-based classifier.
- Output formats: BLAST `outfmt=6` (tabular, 12 columns by default), `outfmt=7` (tabular with header comments), `outfmt=5` (XML), `outfmt=0` (the human-readable default). For pipeline use, format 6 or 7. For programmatic parsing in Python, format 5 with `Bio.Blast.NCBIXML`.
- Classifier basics: top-hit classifier (genus of the lowest-E-value hit), top-N consensus (majority genus among the top-5 hits below an E-value cutoff), lowest common ancestor (LCA: the deepest taxonomic node shared by all hits above a confidence threshold). LCA is the de facto standard in metagenomics tools like Kraken and DIAMOND-MEGAN.
- Failure modes: **low-complexity sequence** (`AAAAAAA...` runs match everything; use `dustmasker` for DNA, `seg` for protein), **chimeric reference entries** (curate your database), **paralog vs ortholog confusion** (a top hit to a paralog of a different species is not informative about the species; use orthology-aware methods like OrthoFinder when this matters).

## Weekly schedule

The schedule below adds up to approximately **36 hours**. Treat it as a target. Monday's lecture on how BLAST actually works is the single hour that decides whether the rest of the week makes sense — read it twice if needed.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | How BLAST works: seed-and-extend, K-A statistics   |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | Running blastn / blastp locally and via NCBI       |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Wednesday | E-values, bit scores, output formats               |    1h    |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     6.5h    |
| Thursday  | Parsing BLAST output, taxonomy via Entrez          |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Classifier design, LCA, failure modes              |    0h    |    1h     |     1h     |    0.5h   |   1h     |     2h       |    0h      |     5.5h    |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **7.5h**  | **3h**     | **3h**    | **6h**   | **8h**       | **2h**     | **35.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | BLAST+ docs, NCBI Entrez E-utils, Karlin-Altschul paper, NCBI Taxonomy |
| [lecture-notes/01-how-blast-actually-works.md](./lecture-notes/01-how-blast-actually-works.md) | The seed-and-extend heuristic, the two-hit rule, the K-A E-value formula, and exactly what BLAST trades for its speed |
| [lecture-notes/02-running-blast-locally-and-via-ncbi.md](./lecture-notes/02-running-blast-locally-and-via-ncbi.md) | `blastn` / `blastp` / `tblastn` from the command line, `makeblastdb` to build a local database, `Bio.Blast.NCBIWWW` for remote queries, reading and parsing output |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-blast-an-unknown-sequence.py](./exercises/exercise-01-blast-an-unknown-sequence.py) | Submit one unknown DNA sequence to NCBI `blastn` and parse the top hit |
| [exercises/exercise-02-build-local-db.py](./exercises/exercise-02-build-local-db.py) | Build a small local BLAST database from a FASTA file and query it |
| [exercises/exercise-03-parse-blast-output.py](./exercises/exercise-03-parse-blast-output.py) | Parse BLAST tabular and XML output into a pandas DataFrame and filter hits |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-taxonomy-classifier.md](./challenges/challenge-01-taxonomy-classifier.md) | Build a top-hit and LCA classifier, compare on a labelled test set |
| [quiz.md](./quiz.md) | 10 multiple-choice questions on BLAST mechanics, E-values, and the API |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | BLAST-driven taxonomy classifier for ~20 unknown sequences with precision-recall |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("Biopython 1.83," "BLAST+ 2.15.0"). We cite tools by their paper ("BLAST 2.0, Altschul et al. *Nucleic Acids Res.* 1997"). We say "the top hit has E-value 4e-87 and 99.2% identity over 1,247 aligned positions" not "the match looks good." An E-value is a number on a known scale. If your classifier report uses the words "good hit" or "strong match" without a number, you have not written a methods section yet.

## A note on the network

BLAST against NCBI is a *shared service*. Roughly:

- Each search takes 10–120 seconds depending on database size, query length, and how busy NCBI is.
- NCBI asks that automated callers (1) provide an email address (`Bio.Entrez.email = "you@example.com"`), (2) cap requests at 3/second without an API key or 10/second with one, (3) put bulk jobs on the off-hours queue (`Bio.Blast.NCBIWWW.qblast(..., service="psi")` or run during night/weekend US-east-coast hours), and (4) honor `Retry-After` headers on 429 responses.
- For the mini-project's ~20 sequences, expect the BLAST submission phase to take ~10–30 minutes. Plan around it. The exercises use a local database explicitly so you can iterate on parsing without waiting on NCBI.

If NCBI is down or extremely slow on the day you do the exercises, the alternative is to use the EBI BLAST web service (`https://www.ebi.ac.uk/Tools/sss/ncbiblast/`) which has its own queue and is often less loaded. Your homework will note when this fallback applies.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Read the original BLAST 1990 paper (Altschul, Gish, Miller, Myers, Lipman, *J. Mol. Biol.* 215:403) and the BLAST 2.0 paper (Altschul et al., *Nucleic Acids Res.* 25:3389, 1997) end to end. The 1997 paper introduces the two-hit rule and gapped alignment. Both are ~10 pages.
- Install [`DIAMOND`](https://github.com/bbuchfink/diamond) and benchmark a `blastp` query against the SwissProt database in BLAST+ vs DIAMOND. DIAMOND is ~100x faster than BLAST+ for protein-vs-protein at scale and uses the same K-A statistics on the back end. Read the DIAMOND paper (Buchfink et al., *Nat. Methods* 12:59, 2015).
- Reproduce Karlin & Altschul 1990's E-value formula on a small synthetic dataset: generate 1,000 random 100-residue protein sequences, all-vs-all `blastp`, plot the observed score distribution against the predicted Gumbel extreme-value distribution. The fit should be tight for high scores and noisy in the tail. Cite Karlin & Altschul, *PNAS* 87:2264, 1990.
- Build a **Kraken-style k-mer classifier** on a small dataset and benchmark against your BLAST classifier. Kraken (Wood & Salzberg 2014) trades exact alignment for an exact-k-mer lookup against a precomputed LCA index; it is ~1000x faster than BLAST at similar accuracy for metagenomic classification.

## Up next

Continue to [Week 5 — Read alignment](../week-05/) once you have pushed your mini-project to GitHub. Week 5 takes seed-and-extend from the BLAST setting (one query at a time against a static database) to the read-alignment setting (millions of short reads against a single reference genome), and you will meet the FM-index + Burrows-Wheeler Transform combination that makes BWA possible.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
