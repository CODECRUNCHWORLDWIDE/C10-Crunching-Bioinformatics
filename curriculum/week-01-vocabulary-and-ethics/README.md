# Week 1 — Vocabulary and Ethics

Welcome to **C10 · Crunching Bioinformatics**. Week 1 is unusual: we will not call a single bioinformatics tool, parse a single FASTQ file, or run a single alignment. Instead, we install the two things every later week depends on — **the vocabulary** and **the ethics**.

Most short bioinformatics tutorials skip both. They open with `pip install biopython`, hand you a sample FASTA, and you walk away two weeks later confidently misusing the word "gene." We are not doing that. By Friday of Week 1 you will be able to **explain the central dogma in five minutes to someone with no biology background, parse a FASTA file in pure Python without Biopython, and articulate the data-ethics rules** that govern every later week of this course.

The other half of the week is a tour of the public data infrastructure of bioinformatics — NCBI, Ensembl, EBI, GISAID, GTEx, 1000 Genomes. You will leave the week with a one-page personal data inventory listing, for each public source you will use in C10, what it contains, what its access policy is, and what its citation requirement is.

## Learning objectives

By the end of this week, you will be able to:

- **Define** every term in the Week 1 glossary in your own words, without copying a Wikipedia sentence.
- **Explain** the central dogma — DNA to RNA to protein — in five minutes, on a whiteboard, to a non-biologist.
- **Distinguish** the words *genome*, *chromosome*, *gene*, *transcript*, *exon*, *codon*, and *nucleotide* — without conflating them.
- **Name** the canonical bioinformatics file formats (FASTA, FASTQ, SAM/BAM, VCF, GFF/GTF) and what each is for, even though we go deep on them only in later weeks.
- **Parse** a FASTA file in pure Python — no Biopython, no regex tricks, just `open()` and a `for` loop.
- **Identify** at least four public consent-cleared datasets and what each is for (1000 Genomes, GTEx-public, GISAID, RefSeq, Ensembl reference genomes).
- **Articulate** the two data-ethics rules of C10 and explain *why* each one exists, not just *that* it exists.
- **Recognize** when a sentence about genetics crosses from "associated with" into genetic determinism — and avoid that language in your own writing.

## Prerequisites

This week assumes you have completed **C1 weeks 1–11**, or have equivalent skill. Specifically:

- Comfortable in a terminal — `cd`, `ls`, `python`, `pip`.
- You can read and write functions, dictionaries, list comprehensions, file IO.
- You know high-school biology — DNA is a double helix, proteins are made of amino acids, cells have a nucleus. We will re-teach what's needed but we will not stop to explain mitosis.
- You have a public GitHub account and have committed code from a terminal at least once.

If any of those are shaky, **stop** and review the relevant C1 week before continuing. C10 is paced for someone who already has Python comfort and is here for the biology.

## Topics covered

- What bioinformatics actually is, and what it is not
- The 90-minute central dogma — DNA, RNA, protein, transcription, translation
- The genetic code and the codon table — 64 codons, 20 amino acids, redundancy and wobble
- Genome vs chromosome vs gene vs transcript vs nucleotide — the vocabulary problem
- Reference genomes — what GRCh38, GRCm39, T2T-CHM13 are and how they relate
- File-format tour: FASTA, FASTQ, SAM, BAM, VCF, GFF, GTF, BED (we go deep on these in Weeks 2–6)
- The "vocabulary problem" — biologists and engineers use the same words differently
- Why genetic data is uniquely sensitive — the re-identification literature
- The two C10 data-ethics rules: public consent-cleared datasets only; never analyze a friend or family member's DNA
- The public-data ecosystem: NCBI, Ensembl, EBI, GISAID, GTEx, TCGA, dbSNP, ClinVar
- IRB, GDPR, HIPAA at a high level — we are not lawyers; we point at the right experts
- The open-science norm and why bioinformatics is structurally open-source-friendly
- Citing versions and accessions like a methods section, not a tweet

## Weekly schedule

The schedule below adds up to approximately **36 hours**. Treat it as a target. Some sections will click in 20 minutes, others will need 3 hours. That is fine.

| Day       | Focus                                              | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|----------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | Central dogma, vocabulary                          |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | File-format tour, glossary work                    |    1h    |    2h     |     1h     |    0.5h   |   1h     |     0h       |    0.5h    |     6h      |
| Wednesday | Data ethics, re-identification, IRB                |    2h    |    1.5h   |     1h     |    0.5h   |   1h     |     0h       |    0.5h    |     6.5h    |
| Thursday  | Public-data tour, data inventory                   |    1h    |    1.5h   |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     7.5h    |
| Friday    | FASTA-by-hand, reverse complement                  |    0h    |    1.5h   |     1h     |    0.5h   |   1h     |     2h       |    0h      |     6h      |
| Saturday  | Mini-project deep work                             |    0h    |    0h     |     0h     |    0h     |   1h     |     3h       |    0h      |     4h      |
| Sunday    | Quiz, review, polish                               |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                    | **6h**   | **8h**    | **4h**     | **3h**    | **6h**   | **7h**       | **2h**     | **36h**     |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | Curated readings, NCBI/Ensembl handbooks, free books, official docs |
| [lecture-notes/01-the-central-dogma-in-90-minutes.md](./lecture-notes/01-the-central-dogma-in-90-minutes.md) | DNA to RNA to protein, vocabulary, file formats |
| [lecture-notes/02-data-ethics-and-public-data-sources.md](./lecture-notes/02-data-ethics-and-public-data-sources.md) | Consent, re-identification, IRB, public datasets |
| [exercises/README.md](./exercises/README.md) | Index of short drills |
| [exercises/exercise-01-glossary.md](./exercises/exercise-01-glossary.md) | Define every term in your own words |
| [exercises/exercise-02-fasta-by-hand.py](./exercises/exercise-02-fasta-by-hand.py) | Parse a FASTA file in pure Python |
| [exercises/exercise-03-data-inventory.md](./exercises/exercise-03-data-inventory.md) | Catalog the public datasets you will use |
| [challenges/README.md](./challenges/README.md) | Index of weekly challenges |
| [challenges/challenge-01-reverse-complement.py](./challenges/challenge-01-reverse-complement.py) | `reverse_complement`, `GC_content`, `translate` from scratch |
| [quiz.md](./quiz.md) | 10 multiple-choice questions, biology + ethics |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | A 1-page glossary + a public-data inventory |

## A note on tone

C10 is written in **lab-notebook voice**. We cite versions ("Biopython 1.83," "1000 Genomes phase 3," "GRCh38.p14"). We distinguish biological claim from statistical claim ("gene X is associated with" — not "gene X causes"). We acknowledge biological uncertainty. We do not use determinism language. This is not stylistic preference; this is how the field actually writes, and learning the voice now will make your Week 12 capstone read like a real methods section instead of a blog post.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Read the *introductions only* of three Nature Methods papers from the last six months. Note the citation style, the version-pinning, and how each paper signals data availability.
- Skim the 1000 Genomes Project consent documents (linked in `resources.md`). Note what donors were and were not consented for.
- Read the ELSI (Ethical, Legal, Social Implications) page on the NHGRI site end-to-end. Bookmark two sources you trust for genetic-ethics questions.
- Browse the [Bioinformatics Data Skills](https://vincebuffalo.com/bds/) free chapters and pick one to read this week.

## Up next

Continue to [Week 2 — FASTA, FASTQ, IO](../week-02-fasta-fastq-io/) once you have pushed your mini-project to GitHub.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
