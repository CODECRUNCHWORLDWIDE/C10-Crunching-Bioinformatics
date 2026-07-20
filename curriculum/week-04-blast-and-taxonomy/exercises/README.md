# Week 4 — Exercises

Three focused Python drills. Each one runs from the command line, finishes inside 60 minutes, and produces output you can keep. Together they take you from "I read about BLAST" to "I can submit a BLAST query, build a local database, and parse the results into pandas." The mini-project then composes these three into a working classifier.

## Index

1. **[Exercise 1 — BLAST an unknown sequence](exercise-01-blast-an-unknown-sequence.py)** — submit one unknown 16S rRNA sequence to NCBI `blastn` via `Bio.Blast.NCBIWWW`, parse the top hit, and identify the genus. (~50 min, network-bound — most of the time is waiting on NCBI.)
2. **[Exercise 2 — Build a local BLAST database](exercise-02-build-local-db.py)** — download a small reference FASTA (the *E. coli* K-12 MG1655 reference genome from NCBI accession `NC_000913.3`), build a local `blastn` database with `makeblastdb`, and query an unknown sequence against it. (~45 min, no network after the initial fetch.)
3. **[Exercise 3 — Parse BLAST output](exercise-03-parse-blast-output.py)** — read both tabular (`outfmt 6`) and XML (`outfmt 5`) BLAST output, filter by E-value and percent identity, and produce a one-hit-per-query summary DataFrame. (~50 min, pure Python — no BLAST runs.)

## How to work the exercises

- Install the toolchain *first*: `conda install -c bioconda blast=2.15 biopython=1.83 pandas`. The exercises assume `blastn`, `makeblastdb`, and `blastdbcmd` are on your `PATH`.
- **Set `Bio.Entrez.email`** to a real address (yours) at the top of every script. NCBI may block your IP if you do not.
- **Type the code yourself.** Copy-paste defeats the muscle-memory point.
- Each file is runnable: `python exercise-XX.py`. All assertions must pass.
- Exercise 1 makes a real network call to NCBI. Expect 30–120 seconds of wait time. Cache the XML response to disk on first run; subsequent runs read from cache.
- If you get stuck for more than 10 minutes, peek at the inline hints, then re-read Lecture 2.
- No solutions are checked in. The exercises are self-checking — every script ends in an assertion block.

## Accession IDs used this week

Cited by NCBI accession so you can verify your data is the same as the curriculum's:

- **`NC_000913.3`** — *Escherichia coli* str. K-12 substr. MG1655 complete genome (4,641,652 bp). The canonical bacterial reference genome. Used in Exercise 2.
- **`NR_117741.1`** — *Staphylococcus aureus* 16S rRNA gene reference sequence (~1,547 bp). Used as the worked example in Exercise 1.
- **`NR_074549.1`** — *Escherichia coli* K-12 MG1655 16S ribosomal RNA reference (~1,541 bp). Used in Exercise 2.
- **`NM_007294.4`** — *Homo sapiens* BRCA1 mRNA (7,224 bp). Used in the Lecture 1 §5 sanity check.

If any of these accessions has been retired or updated by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt.
