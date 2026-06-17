# Lecture 1 — FASTA: The Line Format That Won

> **Duration:** ~2 hours of reading + hands-on.
> **Outcome:** You can read, write, index, and round-trip FASTA files with Biopython 1.83; you can decide between `SeqIO.parse`, `SeqIO.read`, `SeqIO.to_dict`, and the lower-level `SimpleFastaParser`; and you know enough about description-line conventions across NCBI, GenBank, and UniProt to not be tripped up by them at 11 PM the night before a deadline.

If you only remember one thing from this lecture, remember this:

> **FASTA is the simplest possible serialization of a labelled sequence.** A `>` line introduces a record, all following lines are the sequence until the next `>` (or end of file). Forty years of bioinformatics have not replaced it because there is nothing left to optimize at the format level. The complexity lives in the *description line* — and that is where every bug you will hit this week comes from.

---

## 1. Why FASTA won

The FASTA format was introduced in 1985 by the FASTP program (Lipman & Pearson, *PNAS* 82:1435). It is older than the World Wide Web. It predates UTF-8. It is plain ASCII, line-oriented, and trivially parseable with `grep` and `awk`. Every tool in bioinformatics speaks it. Every database in bioinformatics exports it.

There have been many proposals to replace FASTA — XML-based formats in the 2000s, JSON-based formats in the 2010s, schema-typed Protobuf formats more recently. None won. The reason is **operational**: a 3 GB human-genome FASTA is `grep`-able. A 3 GB JSON is not. The cost of "richer" formats was always paid by the human at the terminal at 2 AM, and that human voted with their hands.

Two practical consequences:

1. **You will be reading FASTA in 2026, 2030, and (assuming the field still exists) 2050.** Learn it once, properly.
2. **The format itself has zero ambiguity. The ambiguity is in the description line.** That is where Week 2's bugs come from.

---

## 2. The format, formally

A FASTA file is a sequence of **records**. Each record is:

```
>description line
ACGT
ACGT
ACGT
```

- Exactly one **header line** starting with `>`. Everything after the `>` on that line is the **description**. By convention the first whitespace-delimited token is the **identifier** and everything after it is the **free-text description**.
- One or more **sequence lines**. By convention each line is 60 or 80 characters wide, but parsers must accept any width including a single line of arbitrary length.
- The record ends at the next `>` line or end of file.

That is the entire specification. There is no `END` marker. There is no record count at the top. There is no escape sequence for a literal `>` in a sequence line (because sequence lines do not contain `>`).

> **The implicit rule.** "All lines after the header until the next `>` line are sequence." That is the rule that lets you `cat` two FASTA files together and have a valid FASTA file. It is also the rule that lets a malformed file (one with a stray `>` in a description) silently corrupt your parse.

Multi-FASTA files are just FASTA files with more than one record. The human reference genome (GRCh38) is a multi-FASTA with one record per chromosome — 25 records (chromosomes 1-22, X, Y, MT) plus a few hundred "alternate haplotype" records.

---

## 3. The description-line zoo

The description line is where the format earns its reputation. Three conventions you must recognize:

### NCBI

NCBI uses a space-delimited convention. The identifier is the first whitespace-delimited token, conventionally an **accession** (like `NC_045512.2` — the SARS-CoV-2 reference genome) followed by free-text description.

```
>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCAC
```

The accession `NC_045512.2` tells you exactly which record this is, in any NCBI database, forever. The `.2` is the *version* — `.1` was the original 2020-01-05 deposit, `.2` corrects a minor issue. Pin the version.

### GenBank FASTA

GenBank's downloaded FASTA uses a pipe-delimited "FASTA defline" with multiple ID tokens:

```
>gi|1798174254|ref|NC_045512.2| Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
```

The `gi|...|ref|...` style was officially deprecated in 2016 (GI numbers were phased out) but you will see it in files generated before the change. Biopython parses the first token (`gi|1798174254|ref|NC_045512.2|`) as the record id; you frequently want only the `NC_045512.2` part, which means splitting on `|` and picking the appropriate field.

### UniProt

UniProt's FASTA defline is the most structured of the three:

```
>sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S PE=1 SV=1
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS
NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFAPLLRSYSFNTQRNFLLDIPCYFTNGLT
```

- `sp` means **Swiss-Prot** (manually curated). `tr` means **TrEMBL** (auto-curated).
- `P0DTC2` is the UniProt accession.
- `SPIKE_SARS2` is the entry name.
- The key-value tags `OS=`, `OX=`, `GN=`, `PE=`, `SV=` give Organism (Species), Organism (taXonomy id), Gene Name, Protein Evidence level, and Sequence Version. They are pleasantly structured, but they are not part of the FASTA spec — they are UniProt's local convention.

> **The takeaway.** The format is trivial; the **identifier conventions are not portable** between databases. Always confirm which database your FASTA came from before parsing the description line for metadata.

---

## 4. Biopython 1.83 — the four functions that do 95% of the work

Install:

```bash
python -m pip install biopython==1.83
```

Then in Python:

```python
from Bio import SeqIO
```

Four functions cover almost everything you will do with FASTA this week.

### `SeqIO.parse(handle, "fasta")` — the iterator

```python
from Bio import SeqIO

for record in SeqIO.parse("genome.fasta", "fasta"):
    print(record.id, len(record.seq))
```

`parse` returns a **generator** of `SeqRecord` objects. It does **not** load the file into memory. Use it for any file larger than a few hundred MB, or any time you only need one record at a time. This is the right default.

### `SeqIO.read(handle, "fasta")` — the single-record shortcut

```python
record = SeqIO.read("spike_protein.fasta", "fasta")
print(record.id, len(record.seq))
```

`read` asserts that the file contains **exactly one record** and returns it. If there are zero or more than one, it raises `ValueError`. Use it for single-sequence files (one chromosome, one protein, one reference). If you accidentally point `read` at a multi-FASTA, you want the explicit error rather than the silent "I'll just take the first record" behaviour that hand-rolled parsers tend to default to.

### `SeqIO.to_dict(SeqIO.parse(...))` — the in-memory index

```python
records = SeqIO.to_dict(SeqIO.parse("transcripts.fasta", "fasta"))
print(records["ENST00000288602.11"].seq[:60])
```

`to_dict` materializes the iterator into a dict keyed by `record.id`. **It loads the entire file into RAM.** A human reference genome will fit; a 50 GB long-read FASTQ will not. Use `to_dict` for files of a few hundred MB or smaller.

For larger files, use `SeqIO.index()` (memory-mapped, on-demand) or `SeqIO.index_db()` (SQLite-backed). The interface is the same — dictionary access by `record.id` — but the cost is paid on read, not on load. We use these in Week 5 when we routinely deal with multi-GB files.

### `SeqIO.write(records, handle, "fasta")` — the writer

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

new_records = [
    SeqRecord(Seq("ATGCGTACGT"), id="seq1", description="my first record"),
    SeqRecord(Seq("GGGGGGGGGG"), id="seq2", description="all Gs"),
]
SeqIO.write(new_records, "output.fasta", "fasta")
```

`write` takes an **iterable** of records — a list, a generator, or another `SeqIO.parse` call — and writes them out. By default, sequence lines are wrapped at 60 characters wide (configurable via the lower-level `FastaIO.FastaWriter`).

The pattern `SeqIO.write(SeqIO.parse(input, "fasta"), output, "fasta")` is a streaming copy — useful for renaming records, filtering on length, or converting line widths, without loading the file.

---

## 5. `SeqRecord` vs `Seq` vs raw `str` — when to use which

A `SeqRecord` is a wrapper around a sequence with metadata. Its key attributes:

```python
record.id           # str — the identifier (first token of the description)
record.description  # str — the full description line, including the id
record.seq          # Bio.Seq.Seq — the sequence itself
record.annotations  # dict — format-specific metadata
record.features     # list[SeqFeature] — GenBank-style annotations
record.letter_annotations  # dict — per-letter info (e.g. quality scores)
```

The `Seq` object inside is a string-like wrapper that knows it is biological:

```python
record.seq.reverse_complement()  # returns a new Seq
record.seq.translate(table=1)    # returns a new Seq (protein)
record.seq.transcribe()          # DNA -> RNA (T -> U)
str(record.seq)                  # convert to a plain str
```

The plain `str` is what you reach for in **tight loops** — slicing, counting characters, hashing. Every method call on a `Seq` object goes through Biopython's alphabet-checking machinery, which is the correct trade-off for application code and the wrong one for inner loops.

**Decision rule.**

- Use **`SeqRecord`** whenever you need to round-trip metadata (read in, transform, write out).
- Use **`Seq`** when you want the biological methods (reverse_complement, translate) without dropping back to a plain string.
- Use **`str`** when you are inside a loop that touches every character, or whenever the metadata doesn't matter.

Cast aggressively. `str(record.seq)` is cheap. The mistake students make is wrapping every transient computation in `SeqRecord`, which makes their code 5x slower than necessary and obscures what's going on.

---

## 6. A worked example — extract chr22 from GRCh38

Suppose you have the GRCh38 primary-assembly FASTA (one file, 3.0 GB, all chromosomes). You want a single-chromosome FASTA for chr22 only (the smallest autosome, our standing example for "small enough to play with").

```python
from Bio import SeqIO

INPUT  = "GRCh38.primary_assembly.genome.fa"
OUTPUT = "chr22.fasta"
TARGET_ID = "chr22"

with open(OUTPUT, "w") as out_handle:
    for record in SeqIO.parse(INPUT, "fasta"):
        if record.id == TARGET_ID:
            SeqIO.write([record], out_handle, "fasta")
            break  # chr22 appears once; no need to keep scanning
```

Notes:

- `SeqIO.parse` is a generator — `chr1` is processed before `chr2` is loaded.
- We open `OUTPUT` once and write a single record into it. Writing one record at a time via `SeqIO.write` is fine; the function accepts any iterable.
- The `break` is the optimisation. Without it, the loop would continue scanning the file to end-of-file even though it already had its answer. Reading the rest of GRCh38 from disk is not free.

This pattern — stream in, filter, stream out — is the canonical FASTA-IO idiom in bioinformatics. Memorize it.

---

## 7. Common gotchas

**Lowercase, ambiguity codes, and N's.** Reference genomes routinely contain `N` (any base, used in unsequenced regions) and lowercase letters (which conventionally mark repetitive regions, "soft-masked"). Your parser must accept both. If you `assert seq.upper() in "ACGT"`, you have a bug that will fire on the first real file you hand it.

**Stray `>` in descriptions.** Some pipelines emit descriptions that include `>` characters mid-line. The FASTA spec does not forbid this, but parsers split on `>` at the start of a line; embedded `>` in a description is fine. The bug pattern is: someone wrote a parser that splits the *whole file* on `>` (rather than recognising the `>` only at line start). Don't be that person.

**Mixed line endings.** Files generated on Windows (`\r\n`), Mac classic (`\r`), or Unix (`\n`) all show up in production. Biopython handles all three transparently. Your own pure-Python parser from Week 1 should also — `line.rstrip()` (with no argument) strips all standard whitespace including `\r`.

**Sequence on the header line.** Some tools (notably very old EMBOSS scripts) emit the sequence immediately after `>id`, on the same line, with no newline. This is non-conformant, but you will see it. Biopython will fail; if your data fails to parse, check this.

**Empty records.** A header with no following sequence is a valid record with an empty sequence. Whether you should *treat* it as valid depends on context — in QC, an empty record usually means a bug upstream. Always count.

---

## 8. When to drop down from `SeqIO` to `FastaIO`

`Bio.SeqIO.parse(..., "fasta")` builds a full `SeqRecord` for every record, which means it allocates `Seq` and metadata objects on every iteration. For multi-million-record files, this is the bottleneck.

The lower-level `Bio.SeqIO.FastaIO.SimpleFastaParser`:

```python
from Bio.SeqIO.FastaIO import SimpleFastaParser

with open("transcripts.fasta") as handle:
    for title, sequence in SimpleFastaParser(handle):
        # title: str, the description line without the leading '>'
        # sequence: str, the joined sequence with newlines removed
        pass
```

returns plain `(title, sequence)` tuples of strings. It is roughly **3-5x faster** than the full `SeqIO.parse` route on large transcript files. Use it when:

- You don't need the `Seq` methods.
- You don't need the GenBank-style features.
- You are processing the file with tight string operations of your own.

If you don't know which to use, stay with `SeqIO.parse`. The performance gap matters only when you are processing millions of records and you have already profiled.

---

## 9. Writing FASTA — the parts people get wrong

Three rules for FASTA output:

1. **Line-wrap your sequence.** Default to 60 or 80 columns. Some downstream tools fail (or simply hang) on multi-megabase single-line FASTA, especially text editors and older parsers. Biopython's default is 60.
2. **Do not put a space between `>` and the id.** `>NC_045512.2 description` is correct. `> NC_045512.2 description` will be parsed by some tools as a record with id `""` and description `"NC_045512.2 description"`, which silently breaks everything.
3. **End the file with a newline.** POSIX text files end in `\n`. Some tools strictly require this; some are forgiving; you will not always know which is which.

Biopython's `SeqIO.write` gets all three right by default. The bugs happen when people hand-write FASTA from f-strings without using `SeqIO.write`. Don't.

---

## 10. A note on the GenBank format (cousin to FASTA)

The **GenBank** flat-file format is a richer text serialization used by NCBI to ship annotated reference sequences. A GenBank file contains the same sequence as the equivalent FASTA, plus a header with metadata (organism, references, features) and a feature table (genes, CDS, repeats, etc.). Biopython parses it with the same `SeqIO.parse(handle, "genbank")` interface.

We won't go deep on GenBank this week — Week 6 (annotation) uses it. But know that the relationship is: FASTA is the sequence alone, GenBank is the sequence plus the annotations. If someone hands you a `.gb` file and says "the FASTA," they mean "the sequence from the GenBank file"; you can extract it with `SeqIO.convert("input.gb", "genbank", "output.fasta", "fasta")`.

---

## 11. Recap and next lecture

You should now be able to:

- Parse a FASTA file with `Bio.SeqIO.parse` and explain why it's a generator.
- Read a single-record FASTA with `Bio.SeqIO.read` and articulate when to prefer it over `parse`.
- Write a FASTA file with `Bio.SeqIO.write`, line-wrapped, with correct id and description handling.
- Recognize NCBI, GenBank, and UniProt description-line conventions.
- Choose between `SeqRecord`, `Seq`, and `str` for a given task.
- Identify the lower-level `SimpleFastaParser` and explain when to drop to it.

In [Lecture 2](./02-fastq-and-quality-scores.md) we add the format that every short-read pipeline begins with: **FASTQ**. We will see why the format adds a quality string, what a Phred score actually means, and why the "encoding zoo" (Sanger / Solexa / Illumina 1.3 / 1.5 / 1.8) is a real problem you must be prepared to detect and route around in production.

---

*Continue to [Lecture 2 — FASTQ and Quality Scores](./02-fastq-and-quality-scores.md).*
