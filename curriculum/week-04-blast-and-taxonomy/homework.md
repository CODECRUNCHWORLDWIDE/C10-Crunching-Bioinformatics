# Week 4 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-04/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — BLAST a single sequence against NCBI nt and answer five questions

**Problem statement.** Fetch the *Drosophila melanogaster* alcohol dehydrogenase (ADH) mRNA reference, NCBI accession `NM_001275680.1` (~1,400 bp), via `Bio.Entrez.efetch`. Submit it to NCBI BLAST as a `blastn` query against the `refseq_rna` database with `expect=1e-10` and `hitlist_size=20`. Cache the XML response to `homework/cache/p1_adh.xml`. Then in `homework/p1_adh_blast.py`, parse the cached XML and answer the following questions in `homework/notes/p1-adh-answers.md`:

1. What is the accession and species of the top hit (lowest E-value)?
2. What are the bit score, E-value, and percent identity of the top hit?
3. Among the top 20 hits, how many are from genus *Drosophila*? List the species.
4. Is there a non-*Drosophila* hit in the top 20? If so, what genus, and at what bit score does it appear?
5. Based on the top hit's E-value and bit score, how confident are you that the query is correctly identified? Use the cheat-sheet thresholds in `resources.md`.

**Acceptance criteria.**

- `python homework/p1_adh_blast.py` runs, reads the cache (or fetches+caches), and prints a short report to stdout.
- The XML response is cached to `homework/cache/p1_adh.xml` (committed to the repo).
- `notes/p1-adh-answers.md` contains five numbered, complete-sentence answers.
- Commit message like `p1: NM_001275680 ADH against refseq_rna, top hit at E=...`.

**Hint.** The top hit is almost certainly the query itself (or its RefSeq replacement). The interesting biology is in hits #2 through #20, where you start to see the other *Drosophila* species' ADH genes — *D. simulans*, *D. yakuba*, *D. erecta*, etc., in roughly their phylogenetic distance from *D. melanogaster*.

**Estimated time.** 60 minutes.

---

## Problem 2 — Build and query a local protein database with `blastp`

**Problem statement.** Download the human BRCA1 protein sequence (UniProt accession `P38398`, ~1,863 amino acids) and the mouse BRCA1 protein sequence (UniProt accession `P48754`, ~1,812 aa). Build a local protein BLAST database containing both. Then use the human BRCA1 as the query against the database with `blastp`. Save tabular output as `homework/p2_results.tsv`.

In `homework/notes/p2-brca1-pblast.md`, answer:

1. The top hit is the query against itself. What is its bit score?
2. The second hit is human vs mouse. What is the bit score, E-value, and percent identity?
3. Compute the bit-score *ratio* of the mouse hit to the human-self hit. What does this tell you about how well-conserved BRCA1 is between human and mouse at the protein level?
4. Repeat the search with `-matrix BLOSUM45` instead of the default BLOSUM62. Does the bit score of the human-mouse hit go up or down? Why? (Hint: BLOSUM45 is tuned for *more distant* homology — its scoring matrix has different absolute values, so the bit score will move in a specific direction.)

**Acceptance criteria.**

- `python homework/p2_brca1_pblast.py` runs (it calls subprocess BLAST under the hood).
- `p2_results.tsv` contains the BLAST output for both BLOSUM62 and BLOSUM45 runs.
- `notes/p2-brca1-pblast.md` contains four numbered answers with the specific numbers.
- Commit message like `p2: BRCA1 protein blastp BLOSUM62 vs BLOSUM45`.

**Hint.** To fetch UniProt sequences directly: `wget https://rest.uniprot.org/uniprotkb/P38398.fasta` and `wget https://rest.uniprot.org/uniprotkb/P48754.fasta`. Concatenate the two files into `homework/data/brca1_db.fasta`, then `makeblastdb -in brca1_db.fasta -dbtype prot -out brca1_db`.

**Estimated time.** 60 minutes.

---

## Problem 3 — Word size and the speed-vs-sensitivity tradeoff

**Problem statement.** Lecture 1 §4.2 introduced word size `W` as the primary speed-vs-sensitivity knob. This problem makes the tradeoff concrete. Using the *Drosophila* ADH query from Problem 1, run `blastn` against `refseq_rna` four times with different settings:

1. `-task megablast` (default `W = 28`).
2. `-task blastn` (`W = 11`).
3. `-task blastn -word_size 7`.
4. `-task blastn-short` (`W = 7`, tuned for short queries).

For each run, record (a) wall-clock time using `time.perf_counter`, (b) total number of hits returned (`len(blast_record.alignments)`), and (c) the lowest bit score among the returned hits. Write `homework/p3_word_size.py` and produce a Markdown table in `homework/notes/p3-word-size.md`:

```
| Setting           | Time (s) | # hits | Min bit score |
|-------------------|---------:|-------:|--------------:|
| megablast (W=28)  |          |        |               |
| blastn (W=11)     |          |        |               |
| blastn W=7        |          |        |               |
| blastn-short      |          |        |               |
```

Add a 200-word commentary on the tradeoff. In particular: did lower `W` change the *top* hit, or only add more hits in the long tail? Why?

**Acceptance criteria.**

- Script runs without crashing.
- Four BLAST results cached separately.
- Markdown table is filled with **real numbers from your machine**.
- Commentary addresses both "did the top hit change?" and "how much extra time per extra hit?"

**Hint.** This problem hits NCBI four times. Run it once during off-hours and cache aggressively. If NCBI rate-limits you (unlikely on 4 requests but possible if you debug-iterate), the EBI BLAST web service is a polite fallback.

**Estimated time.** 60 minutes.

---

## Problem 4 — `tblastn` for gene finding

**Problem statement.** You have a small piece of unannotated bacterial DNA: positions 100,000–105,000 of *Escherichia coli* K-12 MG1655 (`NC_000913.3`). You suspect it encodes the gene `lacZ` (β-galactosidase, ~3,075 bp coding sequence in *E. coli*). Use `tblastn` to test this:

1. Fetch the `LacZ` protein sequence (UniProt `P00722`, ~1,024 aa).
2. Fetch the 5 kb DNA fragment (`Entrez.efetch(db="nuccore", id="NC_000913.3", seq_start=100000, seq_stop=105000, rettype="fasta")`).
3. Build a local nucleotide database from the 5 kb fragment.
4. Query the LacZ protein against it with `tblastn`. Save tabular output.

Answer in `homework/notes/p4-tblastn.md`:

1. Did `tblastn` find a hit? If so, what was the bit score, E-value, and percent identity?
2. Which reading frame (`sframe`)? Is the hit on the forward or reverse strand?
3. What are the `sstart` and `send` coordinates of the hit on the 5 kb fragment? Compute the implied location in the full genome (5 kb fragment starts at 100,000 of the full chromosome).
4. Look up the actual `lacZ` coordinates in *E. coli* K-12 MG1655 (the gene is at 365,529–368,603 on the chromosome, or check NCBI for the current annotation). Did your `tblastn` hit fall in the 100,000–105,000 fragment? If not, what gene did you actually identify? (Hint: the 100k–105k fragment is *not* where `lacZ` lives. You will get a "no hits" result, which is itself the answer.)

**Acceptance criteria.**

- Script runs without crashing.
- Tabular output saved.
- `notes/p4-tblastn.md` answers four questions, including the "no hits" finding if applicable.

**Hint.** The point of this problem is to teach you that "no hits" is a *valid and informative* result. Many homework submissions force a hit by relaxing parameters; do not do that. Report what BLAST found, even if it is nothing.

**Estimated time.** 60 minutes.

---

## Problem 5 — Compute an E-value by hand from a raw score

**Problem statement.** A BLAST search returns a hit with raw alignment score `S = 78` against a query of length `m = 280` residues in a database of effective size `n = 1.5 × 10^11` residues. The substitution matrix in use is BLOSUM62, for which the published K-A parameters under the standard BLAST gap-open / gap-extend are approximately `K = 0.041` and `λ = 0.267` (units: per-score-unit).

In `homework/notes/p5-evalue.md`:

1. Compute the E-value by hand, showing every step:
   - `E = K · m · n · exp(-λ · S)`
   - Substitute the numbers; report the intermediate `λ · S` and `exp(-λ · S)` values.
   - Report `E` in scientific notation.
2. Compute the bit score: `bit = (λ · S − ln K) / ln 2`. Report to one decimal place.
3. Verify both in Python: write `homework/p5_evalue.py` that computes the same two numbers and prints them. The numbers should match your by-hand computation to within rounding.
4. If the database size doubled (to `n = 3.0 × 10^11`) but the same hit was found, what would the new E-value be? What would the new bit score be? Explain in one sentence each.

**Acceptance criteria.**

- `notes/p5-evalue.md` shows the by-hand computation step by step.
- `p5_evalue.py` runs and prints the two numbers.
- Numbers agree to 2 significant figures.
- Final question is answered: doubled database → E-value doubles, bit score unchanged.

**Hint.** Use `math.exp` and `math.log` in Python. `λ · S = 0.267 × 78 = 20.826`. `exp(-20.826) ≈ 9.0 × 10⁻¹⁰`. `E ≈ 0.041 × 280 × 1.5×10¹¹ × 9.0×10⁻¹⁰`. You should get `E ≈ 1.55`. That is a borderline hit — not significant at `1e-3` cutoff but suggestive enough to investigate.

**Estimated time.** 45 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300–400 word reflection at `homework/notes/week-04-reflection.md` answering:

1. Before Week 4, what did you think an E-value was? What is it actually? In what way is the distinction operationally important?
2. The first time you waited 60+ seconds for an NCBI BLAST result to come back, what did you do during the wait? Did anything about the experience change how you would *design* a pipeline that does many BLAST queries?
3. If you ran the Problem 3 word-size comparison: in your own words, when should you reach for `megablast` vs `blastn` vs `blastn-short`? What rule of thumb would you give a junior bioinformatician?
4. The mini-project asks you to classify ~20 unknown sequences with a *measurable* precision-recall accounting. What would make that accounting *honest* vs *cherry-picked*? Name one specific thing you will avoid.

**Acceptance criteria.**

- File exists, 300–400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. Be honest. The mistakes you note here are what you will re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 1 h 0 min |
| 2 | 1 h 0 min |
| 3 | 1 h 0 min |
| 4 | 1 h 0 min |
| 5 | 45 min |
| 6 | 30 min |
| **Total** | **~5 h 15 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
