# Week 2 Homework

Six practice problems that revisit the week's topics. The full set should take about **6 hours**. Work in your `crunch-bio-portfolio-<yourhandle>/week-02/` directory so each problem produces at least one commit you can point to later.

Each problem includes:

- A short **problem statement**.
- **Acceptance criteria** so you know when you are done.
- A **hint** if you get stuck.
- An **estimated time**.

---

## Problem 1 — Read a real FASTQ from SRA

**Problem statement.** Pick one public Illumina FASTQ from the Sequence Read Archive (SRA) or the European Nucleotide Archive (ENA). Good starter accessions: `ERR1019034` (1000 Genomes phase 3, NA12878 chr22 subset), `SRR12345678` (any random recent run). Download just the first 100,000 reads (use `seqkit head -n 100000` after `prefetch` + `fasterq-dump`, or use ENA's direct FASTQ URLs and `gzip -dc | head`).

Write `homework/p1_inspect_fastq.py` that opens the file and prints:

- The accession (from the file name).
- Number of records.
- Mean read length.
- Mean per-read mean Phred quality (yes, the mean of means — it's the right starting summary).
- The detected encoding (Phred+33 vs Phred+64), using the minimum-ASCII heuristic from [Lecture 2 §4](./lecture-notes/02-fastq-and-quality-scores.md).

**Acceptance criteria.**

- File runs as `python homework/p1_inspect_fastq.py path/to/reads.fastq.gz`.
- Uses `gzip.open` if the path ends in `.gz`.
- Uses `Bio.SeqIO.parse` in streaming mode (no `list()`).
- Detected encoding matches what `fastqc --version` reports for the same file (run FastQC and compare).
- Committed with a message like `p1: FASTQ inspector working on ERR1019034 subset`.

**Hint.** `record.letter_annotations["phred_quality"]` gives you an `int` list per record. `statistics.fmean` is your friend; for many records, accumulate sums and divide at the end rather than calling `fmean` per record (it's faster).

**Estimated time.** 60 minutes.

---

## Problem 2 — Round-trip a multi-FASTA through SeqRecord

**Problem statement.** Take the SARS-CoV-2 reference genome (`NC_045512.2`) as your input. Write `homework/p2_roundtrip.py` that:

1. Parses the file with `Bio.SeqIO.parse`.
2. For each record, builds a *new* `SeqRecord` with the same id, the same sequence, and a description that starts with `"[REPROCESSED] "` followed by the original description.
3. Writes the new records back to a different file.
4. Confirms by re-parsing the output that the round-trip preserved the sequence exactly (string equality) and added the prefix in the description.

**Acceptance criteria.**

- File runs as `python homework/p2_roundtrip.py input.fasta output.fasta`.
- The round-trip is loss-less for the sequence.
- The description prefix is present in every output record.
- All `assert`s in the script pass.
- Committed.

**Hint.** Build the new `SeqRecord` with `SeqRecord(seq=record.seq, id=record.id, description=f"[REPROCESSED] {record.description}")`. Watch the description — Biopython's `record.description` includes the id as its first token; if you concatenate naively you will duplicate the id in the output. Strip the id off first.

**Estimated time.** 45 minutes.

---

## Problem 3 — Detect quality encoding from a sample

**Problem statement.** Implement, in `homework/p3_detect_encoding.py`, a function that takes a path to a FASTQ file and returns one of `"Phred+33"`, `"Phred+64"`, or `"unknown"`. Reason from the minimum ASCII character across **the first 1,000 quality strings** (sample, don't scan the whole file).

Test on:

- One file you know is modern Illumina (Phred+33). Any 1000 Genomes FASTQ works.
- One file you generate yourself by re-encoding a Phred+33 file as Phred+64. Hint:
  ```python
  import gzip
  from Bio import SeqIO
  with gzip.open("input.fastq.gz", "rt") as h_in, open("output.fastq", "w") as h_out:
      records = SeqIO.parse(h_in, "fastq")
      SeqIO.write(records, h_out, "fastq-illumina")
  ```
  The `"fastq-illumina"` format is Phred+64. (This is a synthetic test — you would never re-encode real data this way in production.)

**Acceptance criteria.**

- `python homework/p3_detect_encoding.py path/to/reads.fastq` prints exactly one of the three strings.
- Detector correctly classifies both files above.
- Sampling logic stops after 1,000 records (use `itertools.islice`).
- Committed.

**Hint.** Boundary in the lecture: ASCII < 59 → Phred+33; ≥ 64 → Phred+64; in the gap (59–63), prefer Phred+33 but flag as low-confidence detection. With only 1,000 records sampled you may legitimately get "unknown" on a file whose qualities never hit Q0.

**Estimated time.** 45 minutes.

---

## Problem 4 — Build a FASTQ trimmer with a config

**Problem statement.** Take your Exercise 3 `filter_fastq` and harden it. Move it to `homework/p4_trim.py` and add:

- Command-line flags: `--in`, `--out`, `--min-len`, `--min-mean-q`, `--window`, `--window-threshold`. Use `argparse`.
- A `--config` flag that points at a JSON file with the same keys (overrides defaults; command-line flags override JSON).
- A printed **reproducibility receipt** at the end of the run, matching the C10 brand format:
  ```
  ┌─────────────────────────────────────────────────────────┐
  │  REPRODUCIBILITY                                        │
  │                                                         │
  │  Input:        <path>                                   │
  │  Output:       <path>                                   │
  │  Records in:   <n>                                      │
  │  Records out:  <n>                                      │
  │  Settings:     window=4 thr=20 min_len=36 min_q=20      │
  └─────────────────────────────────────────────────────────┘
  ```

**Acceptance criteria.**

- `python homework/p4_trim.py --in reads.fastq --out trimmed.fastq` works with defaults.
- Same invocation with `--config week-02-config.json` and the right keys produces an equivalent run.
- Reproducibility receipt is printed; exit code is 0 on success, 1 on missing input.
- Committed.

**Hint.** Layer the configuration: read JSON if `--config` is set, then layer command-line args on top. `argparse.ArgumentParser.parse_known_args` is useful for letting JSON keys appear without a flag.

**Estimated time.** 75 minutes.

---

## Problem 5 — FastQC sanity run

**Problem statement.** Install FastQC 0.12.1. Run it on **two** FASTQ files: the raw FASTQ from Problem 1, and the trimmed output from Problem 4. Compare the reports. In `homework/notes/fastqc-comparison.md`, write:

- One paragraph summarizing what changed between the raw and trimmed reports.
- A bulleted list of every module that changed status (pass → warn, fail → pass, etc.).
- One paragraph on whether you believe the trimming was *appropriate* for this dataset — was your `min-mean-q=20` too aggressive, too lenient, about right? Cite specific numbers from the report.

**Acceptance criteria.**

- Both FastQC HTML reports are saved into `homework/fastqc/` (commit the `.zip` archives too, but `.gitignore` the unpacked HTML asset folders).
- `notes/fastqc-comparison.md` is committed with a substantive comparison (≥ 250 words).
- C10 voice: no "the data looks good" — give numbers.

**Hint.** FastQC's "Basic Statistics" and "Per sequence quality scores" modules are where the trim impact is most visible. The "Adapter content" module often *worsens* slightly after trimming because shorter reads have a higher fraction of their length occupied by remaining adapter, even if the absolute adapter count fell.

**Estimated time.** 60 minutes.

---

## Problem 6 — Mini reflection essay

**Problem statement.** Write a 300–400 word reflection at `homework/notes/week-02-reflection.md` answering:

1. Of FASTA, FASTQ, Phred encoding, and streaming IO — which felt easiest? Which felt hardest? Why?
2. Did anything you previously believed about sequencing quality turn out to be off this week? If so, what?
3. After running FastQC for the first time, what one module surprised you the most?
4. What is one thing you would want to learn next that this week did not cover?

**Acceptance criteria.**

- File exists, 300–400 words, four numbered paragraphs.
- Committed.

**Hint.** This is for you, not for a grade. Be honest. The mistakes you note here are what you will go back and re-read after the mini-project.

**Estimated time.** 30 minutes.

---

## Time budget recap

| Problem | Estimated time |
|--------:|--------------:|
| 1 | 1 h 0 min |
| 2 | 45 min |
| 3 | 45 min |
| 4 | 1 h 15 min |
| 5 | 1 h 0 min |
| 6 | 30 min |
| **Total** | **~5 h 15 min** |

When you have finished all six, push your repo and open the [mini-project](./mini-project/README.md).
