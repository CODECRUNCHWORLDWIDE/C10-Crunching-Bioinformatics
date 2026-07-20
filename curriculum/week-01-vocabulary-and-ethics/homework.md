# Week 1 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours** in total. Work in your Week-1 Git repository so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Read a real methods section

**Problem statement.** Pick one recent (last two years) open-access paper from any of: *Bioinformatics*, *Genome Biology*, *Nature Methods*, *PLOS Computational Biology*. Read **only the Methods section** end-to-end. In `notes/methods-read.md`, write:

- The paper title, authors, journal, year, DOI.
- A bulleted list of every tool the authors used, with the version they pinned.
- A bulleted list of every dataset they used, with the accession or version.
- One paragraph (~150 words) on what the Methods section did *well* and what (if anything) it left ambiguous.

**Acceptance criteria.**

- File exists at `notes/methods-read.md`.
- Tools list has at least 4 entries, all with versions.
- Datasets list has at least 1 accession.
- The paragraph distinguishes biological claim from statistical claim, in the voice of C10.
- Committed.

**Hint.** [PLOS Computational Biology](https://journals.plos.org/ploscompbiol/) is fully open access and tends to have detailed methods sections. *Nature Methods* sometimes paywalls the full paper but the Methods section is usually accessible.

**Estimated time.** 1 hour.

---

## Problem 2 — Hand-translate a sequence

**Problem statement.** Take the SARS-CoV-2 spike protein coding sequence (the first 60 nucleotides are given below). Translate it **by hand**, using the codon table from Lecture 1. Write the protein sequence in `notes/hand-translation.md` along with the codon-by-codon breakdown.

```
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
```

**Acceptance criteria.**

- File `notes/hand-translation.md` exists.
- The 60-nucleotide sequence is split into 20 codons.
- Each codon's amino acid is identified by hand (one-letter code).
- The protein sequence reads `MFVFLVLLPLVSSQCVNLTT`.
- At the bottom, verify your work by running `python challenges/challenge-01-reverse-complement.py` on the same sequence (so your `translate()` agrees with your hand work).
- Committed.

**Hint.** Use the codon table in [Lecture 1 §4](./lecture-notes/01-the-central-dogma-in-90-minutes.md). Go three letters at a time. If your hand answer and your code answer disagree, the bug is one of: an off-by-one in your codon boundaries, a flipped U/T, or a typo in your `CODON_TABLE`.

**Estimated time.** 45 minutes.

---

## Problem 3 — Re-identification thought experiment

**Problem statement.** Read Gymrek et al. 2013, *"Identifying personal genomes by surname inference"* (linked in [resources.md](./resources.md)). It is 6 pages. In `notes/re-id-essay.md`, write a 400–500 word essay answering:

1. What specific data did the authors use as input?
2. What public databases did they cross-reference to recover identities?
3. Which sex was vulnerable to this particular attack, and why?
4. If you were designing a 2026 consent form for a new genomic study, what would you tell donors about re-identification risk?

**Acceptance criteria.**

- File exists at `notes/re-id-essay.md`.
- 400–500 words, four numbered paragraphs.
- Cites Gymrek et al. 2013 by author and year at least once in the body of the essay.
- Uses C10 voice (no determinism language; cites uncertainty honestly).
- Committed.

**Hint.** The paper itself is short and very readable for a *Science* paper. Read the abstract, introduction, and discussion at full attention; the methods can be skimmed.

**Estimated time.** 1 hour 15 minutes.

---

## Problem 4 — Build a tiny FASTA validator

**Problem statement.** Write a Python script `homework/p4_fasta_validate.py` that takes a path to a FASTA file on the command line and prints a one-line summary, plus warnings for common problems. The summary should include the number of records, total base pairs, mean record length, and overall GC%. The warnings should fire for:

- A record with an empty sequence.
- A record with characters outside `{A, C, G, T, N}` (case-insensitive).
- Two records sharing the same header (collision).

**Acceptance criteria.**

- `python homework/p4_fasta_validate.py path/to/file.fasta` runs.
- The script imports only stdlib (or it may import the `parse_fasta` from your `exercise-02-fasta-by-hand.py`).
- On a valid file, prints exactly one summary line and zero warnings.
- On a file with any of the three problems above, prints the summary plus one warning per problem.
- Exits 0 on valid files, exits 1 on any warning.
- Committed.

**Hint.** Reuse your `parse_fasta` from Exercise 2. The "duplicate header" check needs a `set` of headers seen so far.

**Estimated time.** 1 hour 15 minutes.

---

## Problem 5 — Refine your glossary

**Problem statement.** Take the draft glossary you wrote in [Exercise 1](./exercises/exercise-01-glossary.md) and edit it. Specifically:

- Read each entry out loud. If it does not sound like something you would say in a code review, rewrite it.
- For at least 5 of the 20 entries, add a *second* "easy to confuse with" item, this time chosen because someone in the C10 community Slack actually got it wrong.
- Add a final entry **#21: "the vocabulary problem"** — your own definition of what we mean by that phrase in [Lecture 1 §7](./lecture-notes/01-the-central-dogma-in-90-minutes.md).

**Acceptance criteria.**

- The glossary file is edited and committed with a message like `Glossary v2: tightened wording, added vocabulary problem entry`.
- The 21st entry exists.
- 5 entries have a second "easy to confuse with" note.
- Diff is visible in the commit history (i.e. you did not delete the v1 and start over — you *edited*).

**Hint.** The act of editing is the learning. If your v2 looks identical to your v1, you did not actually re-read it.

**Estimated time.** 45 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300–400 word reflection at `notes/week-01-reflection.md` answering:

1. Of vocabulary, central dogma, and data ethics — which felt easiest? Which felt hardest? Why?
2. Did anything you previously believed about genetics turn out to be off this week? If so, what?
3. Which public dataset are you most curious to work with later in C10, and why?
4. What is one thing you would want to learn next that this week did not cover?

**Acceptance criteria.**

- File exists at `notes/week-01-reflection.md`, 300–400 words.
- Each numbered question is addressed in its own paragraph.
- File is committed.

**Hint.** This is for you, not for a grade. Be honest. Future-you reading it after Week 12 will be grateful.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 1 h 0 min |
| 2 | 45 min |
| 3 | 1 h 15 min |
| 4 | 1 h 15 min |
| 5 | 45 min |
| 6 | 30 min |
| **Total** | **~5 h 30 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
