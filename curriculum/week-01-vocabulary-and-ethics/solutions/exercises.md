# Week 1 — Exercise Solutions

Answers to all three Week-1 exercises, in assignment order. Read [00-overview.md](./00-overview.md) first if you have not.

---

## Exercise 1 — Glossary in Your Own Words

**Assignment:** [../exercises/exercise-01-glossary.md](../exercises/exercise-01-glossary.md)

### The complete answer

This is a model `notes/week-01-glossary.md`. It is a *first draft* — Homework Problem 5 and the mini-project are where it gets polished. All 20 required terms, each with the three mandatory parts, and no example reused between entries.

```markdown
# Week 1 Glossary — first pass

### 1. Nucleotide

One letter of a DNA or RNA sequence. DNA uses four of them (A, C, G, T); RNA
swaps T for U. It is the smallest unit you will ever index into — everything
above it (genes, chromosomes, genomes) is just a longer string of these.

*Example:* the 60-mer `ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC`
is 60 nucleotides long.
*Easy to confuse with:* **base pair**. A nucleotide is one letter on one
strand; a base pair is the A–T or G–C pairing *across* the two strands.
Genome sizes are quoted in base pairs, read lengths in nucleotides, and they
are not the same unit even though both are "3 billion" for a human genome.

### 2. Genome

The entire DNA content of one organism — every chromosome, plus the
mitochondrial DNA, plus any plasmids if it is a bacterium. For a human that is
about 3.2 billion base pairs across 24 distinct chromosome sequences.

*Example:* the SARS-CoV-2 reference genome `NC_045512.2` is 29,903 nucleotides
of single-stranded RNA — five orders of magnitude smaller than ours.
*Easy to confuse with:* **exome**, the protein-coding subset only, which is
roughly 1–2% of the human genome. "Whole-genome sequencing" and
"whole-exome sequencing" are different experiments with different costs.

### 3. Chromosome

One of the large DNA molecules a genome is partitioned into. Humans have 23
pairs, of which we usually sequence 24 distinct sequences: chr1–chr22, chrX,
chrY. In a FASTA reference file, each chromosome is one record.

*Example:* `chr22` in GRCh38 is about 50.8 Mb — the reference file's chr22
record is one header line followed by ~847,000 lines of 60 bases each.
*Easy to confuse with:* **contig** / **scaffold**. Those are assembly
artefacts — stretches of sequence the assembler is confident about, not yet
placed on a chromosome. A draft genome has thousands of contigs and zero
chromosomes.

### 4. Gene

A stretch of a chromosome that gets transcribed into a functional product,
plus (informally) the regulatory sequence around it. It is a *locus* — a named
region with coordinates — not a physical molecule you can hold.

*Example:* `TP53` sits on human chromosome 17 and encodes a 393-amino-acid
protein in its canonical isoform.
*Easy to confuse with:* **ORF (open reading frame)** — a stretch between a
start codon and an in-frame stop codon. Every protein-coding gene contains an
ORF, but plenty of ORFs occur by chance and code for nothing.

### 5. Transcript

One specific RNA molecule produced from a gene. A gene is a region; a
transcript is a product. Alternative splicing means one gene routinely produces
several different transcripts (isoforms), which is why transcript-level and
gene-level counts differ in RNA-seq.

*Example:* the MANE Select transcript for `TP53` is `ENST00000269305`; the gene
has more than a dozen other annotated isoforms.
*Easy to confuse with:* **gene**. If you say "we quantified 20,000 genes" but
your tool emitted `ENST…` identifiers, you quantified transcripts and someone
downstream is going to double-count.

### 6. Exon

A stretch of a gene that survives splicing and ends up in the mature
transcript. Exons are what gets joined together; introns are what gets cut out.

*Example:* the human `HBB` gene (beta-globin) has three exons separated by two
introns; the mature mRNA is the three exons concatenated.
*Easy to confuse with:* **CDS (coding sequence)**. Exons include the
untranslated regions at both ends of the transcript. The CDS is only the part
that gets translated. Exon 1 of many genes is largely 5' UTR.

### 7. Intron

A stretch inside a gene that is transcribed and then cut out of the RNA before
translation. Introns can be enormous relative to the coding sequence, and they
are where most splice-site variants do their damage.

*Example:* the human `DMD` gene spans about 2.2 Mb of chromosome X but its
mature transcript is only ~14 kb — nearly all of that gene is intron.
*Easy to confuse with:* **intergenic region**. Introns are *inside* a gene and
are transcribed. Intergenic sequence lies between genes and is not part of any
transcript.

### 8. Codon

A three-nucleotide window on an mRNA that the ribosome reads as one unit. Four
letters taken three at a time gives 64 codons, which map onto 20 amino acids
plus a stop signal — so the code is redundant, mostly at the third position.

*Example:* `AUG` codes for methionine and also signals the start of
translation; `UAA`, `UAG`, `UGA` are the three stops.
*Easy to confuse with:* **k-mer**. A k-mer is any length-*k* substring at any
offset. A codon is a 3-mer that is also *in frame* — offset by a multiple of
three from the start codon. Frame is the whole difference.

### 9. Amino acid

The monomer a protein is built from. Twenty standard ones, each with a
one-letter code, chained together in the order the codons specify.

*Example:* leucine is written `L` and is encoded by six different codons —
`UUA`, `UUG`, `CUU`, `CUC`, `CUA`, `CUG` — which is the redundancy of the
genetic code made concrete.
*Easy to confuse with:* **residue**. Once an amino acid is joined into a chain
it loses a water molecule and is called a residue. Protein papers count
residues; the numbers are the same, the word is not.

### 10. Protein

A folded chain of amino acids that does chemical or structural work in a cell.
Computationally it is a string over a 20-letter alphabet — until you care about
the fold, at which point it becomes a 3D-structure problem.

*Example:* the SARS-CoV-2 spike protein is 1,273 amino acids long and begins
`MFVFLVLLPLVSSQCVNLTT`.
*Easy to confuse with:* **peptide**. Same chemistry, shorter — conventionally
under ~50 residues, and usually not independently folded.

### 11. Reference genome

A community-curated example sequence for a species, used as the shared
coordinate system everyone reports positions against. It is emphatically not
"the correct genome" — every individual differs from it at millions of
positions, which is normal variation, not error.

*Example:* `GRCh38.p14`, GenBank assembly accession `GCA_000001405.29`.
*Easy to confuse with:* **consensus sequence**. A consensus is computed from a
set of sequences you have in front of you. A reference is a curated artefact
with a version number, a release date, and a patch history.

### 12. Variant (genomic sense)

A position where a sample's sequence differs from the reference. It is a
statement about a *comparison*, not about the sample alone — change the
reference version and the variant list changes with it.

*Example:* the VCF row `22  16050075  rs141297151  A  G` says this sample has a
G where GRCh38 chromosome 22 has an A at position 16,050,075.
*Easy to confuse with:* **mutation**. In clinical genetics "mutation" carries a
connotation of pathogenicity; most variants are harmless. The field has largely
moved to "variant" plus an explicit classification.

### 13. SNP (single-nucleotide polymorphism)

A single-base difference from the reference that is common enough in a
population to be considered polymorphic — the usual informal threshold is 1%
minor allele frequency. Roughly 30–80 of them are enough to identify a person
uniquely, which is the whole reason Lecture 2 exists.

*Example:* `rs6265` in the `BDNF` gene, a common coding variant that has been
reported to influence a range of neurological phenotypes.
*Easy to confuse with:* **SNV (single-nucleotide variant)**. An SNV is any
single-base difference, however rare — a somatic tumour variant seen once is an
SNV, not a SNP. Variant callers emit SNVs; population databases catalogue SNPs.

### 14. FASTA

The plain-text format for sequences with no quality information. A `>` header
line, then one or more lines of sequence, repeated. Used for reference genomes,
gene sets, protein databases, and BLAST queries.

*Example:*
`>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1`
followed by ~499 lines of 60 bases.
*Easy to confuse with:* **FASTQ**. Similar name, different job. A FASTQ record
is always exactly four lines and carries a per-base quality string; a FASTA
record is one header plus *any* number of sequence lines and carries no
quality. A parser that assumes a fixed line count will break on the first
wrapped FASTA record it meets.

### 15. FASTQ

The raw-sequencer format: four lines per read — `@`-header, sequence,
`+`-separator, quality string. Each quality character is a Phred score encoded
as a single ASCII byte, so the sequence and quality lines are always the same
length.

*Example:* a quality character of `I` is ASCII 73; under the Phred+33 offset
used by modern Illumina that is Q40, meaning a 1-in-10,000 chance the base is
wrong.
*Easy to confuse with:* **the `+` line being a separator**. It is a *record*
line that may optionally repeat the header, which means you cannot split a
FASTQ file on `@` — quality strings legitimately contain `@` (ASCII 64, Q31).

### 16. SAM / BAM

The format for reads *after* alignment: where each read mapped, how well, and
how it differs from the reference. SAM is tab-separated text with an `@`-prefixed
header block; BAM is the same information in a binary, block-compressed
container that is roughly 5× smaller and can be indexed for random access.

*Example:* `samtools view -b aln.sam > aln.bam` converts one to the other;
`samtools index aln.bam` then builds the `.bai` so you can ask for chr22 alone.
*Easy to confuse with:* **CRAM**. CRAM is a third member of the family that
compresses further by storing only the differences from the reference — which
means you must keep the exact reference FASTA around or the file is unreadable.

### 17. VCF

Variant Call Format — a table of positions where one or more samples differ
from the reference, with a `##`-prefixed metadata header and one row per
position. Columns 1–9 describe the site; every column after that is one sample.

*Example:* a 1000 Genomes phase-3 chromosome-22 VCF carries 2,504 sample
columns and roughly a million rows.
*Easy to confuse with:* **gVCF**. A gVCF also has rows for positions where the
sample matches the reference, encoded as reference blocks. That is what makes
joint genotyping across samples possible, and it is why a gVCF is enormous
compared to the VCF from the same sample.

### 18. GFF / GTF

Tab-separated annotation: "this feature, of this type, lives from here to here
on this sequence, on this strand." Nine columns, the ninth being a bag of
key–value attributes. GTF is a stricter dialect of GFF with a fixed attribute
convention.

*Example:* a GTF line of the shape
`chr17  HAVANA  exon  1000  1072  .  -  .  gene_id "…"; transcript_id "…";`
places one 73-base exon on the minus strand of chr17.
*Easy to confuse with:* **BED**. Both describe intervals; their coordinate
systems differ. GFF/GTF is 1-based and inclusive at both ends. BED is 0-based
with an exclusive end. That same 73-base exon is `1000  1072` in GTF and
`999  1072` in BED — same interval, different arithmetic — and getting it wrong
is a classic off-by-one that silently shifts every feature by one base.

### 19. Read (sequencer sense)

One sequence record produced by a sequencing machine — a short observed stretch
of a much longer molecule. Reads are the raw material; everything else in a
pipeline is inference from a pile of them.

*Example:* a paired-end Illumina run emits reads of ~150 bp in two files, R1
and R2, where record *n* in each file is the two ends of the same DNA fragment.
*Easy to confuse with:* **fragment** (or **insert**). The fragment is the whole
physical molecule that went into the machine — often 300–500 bp. The two 150 bp
reads are its ends, and the unsequenced gap between them is the inner distance.

### 20. Coverage

How many reads overlap a given position. Reported as a mean over a region
("30× whole-genome"), it is the single most useful number for judging whether a
variant call at a position is trustworthy.

*Example:* a 30× whole-genome Illumina run on a human sample means each base is
covered by ~30 reads on average — which for a 3.2 Gb genome is about 96 Gb of
sequence, or roughly 640 million 150 bp reads.
*Easy to confuse with:* **breadth of coverage**. Depth is how many reads sit on
a position; breadth is what fraction of the target was covered at all. A
capture experiment can have superb mean depth and terrible breadth if the
probes missed a third of the exons.

## Sources I checked against

- NHGRI, *A Brief Guide to Genomics* — <https://www.genome.gov/about-genomics/fact-sheets/A-Brief-Guide-to-Genomics>
- SAM/BAM specification, samtools — <https://samtools.github.io/hts-specs/SAMv1.pdf>
- VCF v4.4 specification, samtools — <https://samtools.github.io/hts-specs/VCFv4.4.pdf>
- GFF3 specification, Sequence Ontology — <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>
```

### Why it works

Three things separate this from a stack of definitions, and each of them is one of the exercise's three required parts doing real work.

**The definitions are relational, not encyclopedic.** Every entry says what the term *is not* as much as what it is. "A gene is a locus, not a molecule." "A variant is a statement about a comparison." That framing is what lets you catch your own mistakes later: when you write `for gene in transcripts:` you will feel the friction, because you defined the two words against each other rather than in isolation.

**The examples are all falsifiable.** Not "for example, a gene" but `TP53`, `ENST00000269305`, `GCA_000001405.29`. This is the lab-notebook voice from [Lecture 1 §5](../lecture-notes/01-the-central-dogma-in-90-minutes.md) applied to your own notes. It also means that when a fact rots — an assembly gets patched, a MANE Select transcript changes — you can tell, because there is a specific thing to check. A vague glossary can never be wrong, which is exactly what makes it useless.

**The "easy to confuse with" entries are where the actual learning lives.** Notice how many of them are engineering hazards, not biology trivia: FASTA versus FASTQ is a parser bug; GTF versus BED is an off-by-one; VCF versus gVCF is a disk-space surprise; SNP versus SNV is a Week-11 cancer-genomics mistake. The exercise asks for "the closest concept this term is *not*" because in practice the two adjacent things are what you will actually mix up under time pressure. Nobody confuses a codon with a chromosome.

One deliberate structural choice: entries 14–18 (the file formats) each pair with the format most likely to be mistaken for it, rather than with a biology term. Formats are engineering objects, and their failure modes are engineering failure modes.

### What a grader is looking for

- **All 20 terms present, each with all three parts.** Missing an "easy to confuse with" is the single most common gap, because it is the hardest part.
- **Definitions between 2 and 4 sentences.** Longer than four means you drifted into mechanism — how splicing works belongs in the lecture notes, not the glossary.
- **No two entries share an example.** This is an anti-laziness rule: if `ATGC` appears under nucleotide, codon, and FASTA, you wrote one example and pasted it three times.
- **Field voice.** "Is associated with," "has been reported to influence," "is a region of." Not "is what makes you tall." Check entry 13 above: `rs6265` "has been reported to influence," which is the honest strength of that claim.
- **A "Sources I checked against" section if you consulted anything**, which you almost certainly did.

Acceptable range is wide. Your definitions will and should read nothing like ours — that is the point of the exercise. What is *not* acceptable is a definition that is merely a shorter Wikipedia sentence, or an "easy to confuse with" that names a distant concept ("a gene is easy to confuse with a protein") rather than the genuinely adjacent one.

### Common wrong turns

**Defining a gene as "a sequence of DNA that codes for a protein."** It is the definition everyone writes first and it is wrong in two directions. Plenty of genes never make protein — `MALAT1`, the rRNA genes, the tRNA genes. And the *gene* is not the coding sequence; the coding sequence is a subset of the exons of some of its transcripts. You will feel this error in Week 7, when a GTF file hands you `gene_biotype "lncRNA"` and your mental model has no slot for it.

**Writing "transcript = mRNA".** Same shape of error. Every mRNA is a transcript; not every transcript is an mRNA. Week 7's kallisto index is built over transcripts, plural and heterogeneous.

**Reusing one example across several entries.** The acceptance criteria forbid it explicitly, and the reason is diagnostic: if you can only produce one concrete example for twenty terms, you have twenty definitions and one piece of knowledge.

**Determinism creep.** "The BRCA1 gene causes breast cancer." Pathogenic BRCA1 variants are *associated with* substantially elevated lifetime risk, with penetrance well under 100%. See the don't-say/say-instead tables in [Lecture 2 §7](../lecture-notes/02-data-ethics-and-public-data-sources.md). If a sentence in your glossary would make you uncomfortable in a methods section, rewrite it now — this habit is much cheaper to build in Week 1 than to retrofit in Week 12.

### How to verify

There is no test runner for prose, so verify structurally. From your week-01 working directory:

```sh
grep -c '^### ' notes/week-01-glossary.md
```

Expected output:

```text
20
```

Then check that every entry has both required annotations — the counts must all match:

```sh
grep -c '^\*Example:\*' notes/week-01-glossary.md
grep -c '^\*Easy to confuse with:\*' notes/week-01-glossary.md
```

Expected output:

```text
20
20
```

Finally, the no-duplicate-examples rule. This prints any example line that appears more than once:

```sh
grep '^\*Example:\*' notes/week-01-glossary.md | sort | uniq -d
```

Expected output: nothing at all. Any line printed is a violation.

Then commit:

```sh
git add notes/week-01-glossary.md
git commit -m "Week 1 glossary, first pass"
```

---

## Exercise 2 — FASTA by Hand

**Assignment:** indexed at [../exercises/README.md](../exercises/README.md) as `exercise-02-fasta-by-hand.py`. **The stub file is not in the repository.** The specification below is taken verbatim from the two places the week does describe it: the [Week-1 README](../README.md) ("Parse a FASTA file in pure Python — no Biopython, no regex tricks, just `open()` and a `for` loop") and [Homework Problem 4](../homework.md) ("it may import the `parse_fasta` from your `exercise-02-fasta-by-hand.py`").

**Assumed interface:** a function `parse_fasta(handle)` that takes an open text file object and yields `(header, sequence)` tuples, where `header` is the text after `>` with surrounding whitespace stripped and `sequence` is every sequence line for that record joined with no separator.

### The complete answer

```python
"""Exercise 2 — FASTA by hand. Pure stdlib, no Biopython, no regex."""

from __future__ import annotations

import sys
from typing import Iterator, TextIO


def parse_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) for every record in an open FASTA file.

    header   -- everything after '>' on the header line, stripped
    sequence -- all sequence lines for that record, concatenated
    """
    header: str | None = None
    chunks: list[str] = []

    for raw in handle:
        line = raw.strip()
        if not line:                      # blank lines are legal padding
            continue
        if line.startswith(">"):
            if header is not None:        # flush the record we just finished
                yield header, "".join(chunks)
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("sequence data before the first '>' header line")
            chunks.append(line)

    if header is not None:                # the last record has no '>' after it
        yield header, "".join(chunks)


def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Convenience wrapper: open a path and stream its records."""
    with open(path, "r", encoding="utf-8") as handle:
        yield from parse_fasta(handle)


if __name__ == "__main__":
    path = sys.argv[1]
    for header, seq in read_fasta(path):
        record_id = header.split()[0] if header else "(empty header)"
        print(f"{record_id}\t{len(seq)}\t{seq[:30]}")
```

Test file — save as `data/good.fasta`:

```
>seq1 first test record
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
>seq2 second test record
ACGTACGTACGTACGTACGT
>seq3 third test record
GGGGCCCCGGGGCCCC
```

### Why it works

The parser is a **flush-on-header state machine**, and that structure is forced by the format rather than chosen for elegance. A FASTA record has no terminator. You do not know a record has ended until you see the *next* record's header — or until the file runs out. That gives the loop exactly two places where a record can be emitted, and both are in the code above: the `if header is not None` inside the `>` branch, and the identical check after the loop. Deleting the second one is the single most common FASTA-parser bug in existence: everything works, and you silently drop the last record of every file.

`header is None` is doing double duty. It is the "have we started a record yet?" flag *and* the guard that rejects a file whose first non-blank line is not a header. Note that it is `is None`, not falsy — a record can legitimately have an empty header (`>` alone on a line), and `if header:` would treat that as "no record started" and lose the sequence.

The function is a **generator**, not a list-builder. That is the whole reason for `open()` plus `for` rather than `handle.read().split(">")`. A GRCh38 FASTA is about 3.1 GB on disk; `read()` puts all of it in memory at once, and the split-on-`>` version briefly holds two copies. The generator holds one record at a time. Week 2's [streaming-large-FASTA challenge](../../week-02-fasta-fastq-io/challenges/challenge-01-streaming-large-fasta.md) is this idea taken further, so it is worth getting the habit here.

Three small format details the code handles that a naive version does not:

- **Sequence wrapping.** Real FASTA wraps at 60, 70, or 80 columns. The `chunks` list plus `"".join()` reassembles them. Appending to a string with `+=` in the loop is quadratic; joining a list at the end is linear.
- **Blank lines.** Legal and common between records. `if not line: continue` absorbs them.
- **Line endings.** `.strip()` removes `\r` as well as `\n`, so a file written on Windows and parsed on Linux does not end up with a carriage return glued to the last base of every line. This is not hypothetical — `\r` in a sequence is a genuinely annoying bug because it prints invisibly.

`encoding="utf-8"` is explicit rather than left to the platform default. On Windows the default is still often cp1252, and a FASTA header containing a species name with a non-ASCII character will raise `UnicodeDecodeError` on one machine and not another.

### Common wrong turns

**Splitting the whole file on `>`.**

```python
records = open(path).read().split(">")[1:]     # do not do this
```

It reads the entire file into memory, it breaks the moment a `>` appears inside a header description (rare but legal), and it leaks the file handle. It is also the version that will make Week 2 impossible on a real reference genome.

**Forgetting the final flush.** Run the parser on the three-record file above without the trailing `if header is not None:` block and you get:

```text
seq1	60	ATGTTTGTTTTTCTTGTTTTATTGCCACTA
seq2	20	ACGTACGTACGTACGTACGT
```

Two records, not three, and no error. Silent data loss is the worst failure mode there is. Always test a parser against a file whose record count you counted by hand.

**Building the sequence with `+=`.**

```python
seq = ""
for line in handle:
    seq += line.strip()          # O(n²) on a 250 Mb chromosome
```

Python strings are immutable, so each `+=` allocates a new string and copies the old one. On chr1 (~248 Mb) this turns a two-second parse into a coffee break. The `chunks` list plus one `join` is the idiom.

**Reaching for `re`.** The assignment says no regex tricks, and the reason is not purity. `line.startswith(">")` is one comparison; a compiled regex for the same job is slower and less readable. Regex is the right tool for parsing FASTA *headers* (Week 4 will pull accessions out of them), not for splitting records.

### How to verify

```sh
python exercises/exercise-02-fasta-by-hand.py data/good.fasta
```

Expected output, exactly (the separator is a literal tab):

```text
seq1	60	ATGTTTGTTTTTCTTGTTTTATTGCCACTA
seq2	20	ACGTACGTACGTACGTACGT
seq3	16	GGGGCCCCGGGGCCCC
```

Three records, and the first is 60 nt — the Homework Problem 2 sequence, so you can reuse this file there.

Then check the guard fires on a malformed file:

```sh
printf 'ACGT\n>late header\nACGT\n' > data/headerless.fasta
python exercises/exercise-02-fasta-by-hand.py data/headerless.fasta
```

Expected output:

```text
Traceback (most recent call last):
  ...
ValueError: sequence data before the first '>' header line
```

---

## Exercise 3 — Public-Data Inventory

**Assignment:** [../exercises/exercise-03-data-inventory.md](../exercises/exercise-03-data-inventory.md)

### The complete answer

A model `notes/week-01-data-inventory.md`. All eight required datasets, all nine columns, no `n/a` without an explanation.

**Read the note under "Why it works" before you copy any version number out of this table.** Several of these are release-numbered series that advance every few months. The pins below are the ones this reference answer was written against; the *exercise* is checking the landing page yourself and writing down what you actually saw on the day you looked.

| Dataset | Provider | URL | Data type | Access tier | Consent scope | Version / release | Citation requirement | C10 week |
|---|---|---|---|---|---|---|---|---|
| NCBI RefSeq | NCBI, US National Library of Medicine | <https://www.ncbi.nlm.nih.gov/refseq/> | Curated reference nucleotide and protein sequences; FASTA and GenBank flat file | Open, no account | Not donor data — curated from public submissions, so no individual consent applies | RefSeq *Homo sapiens* Annotation Release 110 on assembly `GCF_000001405.40` (GRCh38.p14) | O'Leary NA *et al.*, *Nucleic Acids Res* 44:D733–D745 (2016), doi:10.1093/nar/gkv1189 | 2 onwards |
| Ensembl reference genome (GRCh38) | EMBL-EBI, Ensembl project | <https://www.ensembl.org/Homo_sapiens/> | Reference genome FASTA plus GTF/GFF3 annotation, VEP cache, REST API | Open, no account | Not donor data — GRCh38 is a curated mosaic assembly, not one person's genome | GRCh38.p14, GenBank assembly `GCA_000001405.29`; annotation from Ensembl release 113 (October 2024) | The current Ensembl NAR database-issue paper (e.g. *Ensembl 2024*, *Nucleic Acids Res* 52:D891–D899) plus the release number | 5–8 |
| 1000 Genomes Project, phase 3 | IGSR, hosted by EMBL-EBI | <https://www.internationalgenome.org/> | Human germline variation: VCF genotypes, plus alignments and sequence data | Open, no account | Donors consented to broad public research and educational use; data released under the Fort Lauderdale Principles | Phase-3 integrated call set v5a, files dated 2013-05-02 (`…phase3_shapeit2_mvncall_integrated_v5a.20130502…`); GRCh38 liftover released 2017-05-04 | 1000 Genomes Project Consortium, *Nature* 526:68–74 (2015), doi:10.1038/nature15393 | 2, 6 |
| GTEx, public summary tier | NIH Common Fund / Broad Institute, GTEx Portal | <https://gtexportal.org/> | Tissue-level gene-expression summaries (median TPM per gene per tissue) and eQTL summary statistics | Open for the summary tier; the individual-level tier is controlled access via dbGaP | Deceased donors consented via next-of-kin authorisation to broad research use; individual-level data was *not* consented for open release, which is why it sits behind dbGaP | GTEx v8 (`phs000424.v8.p2`, 2017-06-05 data freeze; 948 donors, 54 tissues) | GTEx Consortium, *Science* 369:1318–1330 (2020), doi:10.1126/science.aaz1776, **plus** the portal version you downloaded | 7 |
| GISAID | GISAID Initiative | <https://gisaid.org/> | Viral genome sequences and metadata, chiefly SARS-CoV-2 (EpiCoV) and influenza (EpiFlu) | Registration required; free academic terms | Submitting laboratories consent to sharing under GISAID's terms; downstream users must acknowledge every originating and submitting lab | EpiCoV database, cited by the access date plus the `EPI_SET` DOI generated for your exact download | Elbe S & Buckland-Merrett G, *Global Challenges* 1:33–46 (2017), and Shu Y & McCauley J, *Eurosurveillance* 22(13):30494 (2017); **plus** the `EPI_SET` DOI and the acknowledgement table | 9 |
| UniProt | UniProt Consortium (EMBL-EBI, SIB, PIR) | <https://www.uniprot.org/> | Protein sequences and functional annotation; Swiss-Prot is manually reviewed, TrEMBL is automatic | Open, no account; CC BY 4.0 | Not donor data — protein records are curated from literature and sequence submissions | UniProt releases are dated `YYYY_MM`; pin the exact one you downloaded, e.g. `2024_06` | The UniProt Consortium, *Nucleic Acids Res* 51:D523–D531 (2023), doi:10.1093/nar/gkac1052, or the newer database-issue paper | 4 |
| dbSNP | NCBI, US National Library of Medicine | <https://www.ncbi.nlm.nih.gov/snp/> | Catalogue of short human genetic variation keyed by `rs` identifier; VCF and API | Open, no account | Aggregated submissions from many studies; individual-level genotypes are *not* in dbSNP, only variant records and allele frequencies | Build 156 (released October 2022), on GRCh38.p14 coordinates | Sherry ST *et al.*, *Nucleic Acids Res* 29:308–311 (2001), doi:10.1093/nar/29.1.308, plus the build number | 6 |
| NCBI SRA | NCBI, US National Library of Medicine | <https://www.ncbi.nlm.nih.gov/sra/> | Raw sequencing reads from published studies; FASTQ via `fasterq-dump` | Mixed. Public studies are open; human subject data is frequently controlled access via dbGaP | Varies per study and is stated on the BioProject record. Human-derived open-access runs are usually cell lines; anything from identifiable individuals is normally controlled | SRA has no global version. Pin the accession: study `PRJNA…`/`SRP…`, experiment `SRX…`, run `SRR…` | Leinonen R, Sugawara H, Shumway M, *Nucleic Acids Res* 39:D19–D21 (2011), doi:10.1093/nar/gkq1019, plus the accession and the study's own paper | 5 |

**Notes on access**

- **No account, no application:** RefSeq, Ensembl, 1000 Genomes phase 3, UniProt, dbSNP, and the GTEx public summary tier. Everything C10 asks you to do sits inside this group.
- **Registration required:** GISAID. Free, academic terms, and the account is personal — you may not share credentials or redistribute downloaded sequences. Do not register until Week 9, when you actually need it.
- **Controlled-access application required:** GTEx individual-level data (dbGaP `phs000424`), and many human SRA studies (dbGaP). Both need an institutional signing official and an approved Data Use Certification. **C10 never asks you to apply for either.** They are listed so you recognise the boundary — if a tutorial online tells you to download individual-level GTEx, it is describing a workflow that is not available to you and should not be.
- **Attribution beyond a citation:** GISAID additionally requires an acknowledgement table naming every originating and submitting laboratory, generated from your `EPI_SET` DOI. Budget time for this in Week 9; it is not optional and it is not a formality.

### Why it works

**Every version cell names a thing that can be looked up.** `GCA_000001405.29` resolves to exactly one assembly forever. `phs000424.v8.p2` resolves to exactly one GTEx freeze. `SRR…` resolves to exactly one run. "Current release" resolves to whatever is current on the day someone reads your paper, which is the failure this column exists to prevent. This is the difference the exercise's own good/bad examples are pointing at, and it is the single most portable habit in the week.

**Now the honest part, because it is the actual lesson.** Some rows above name a *series* — Ensembl release 113, dbSNP build 156, UniProt `2024_06`, RefSeq annotation release 110. Those were the current releases when this answer was written, and by the time you read it at least one of them has almost certainly advanced. That is not a defect in the answer; it is the exercise's whole point rendered concretely. A pinned version is a claim about a specific moment. Your inventory should carry a "checked on" date next to each one, and when you cite the data in your Week-12 capstone you re-check rather than trusting a note from ten weeks earlier. If you copy this table without opening a single landing page, you have produced a document that looks like due diligence and contains none.

**The consent column is written to be useful, not decorative.** Compare the RefSeq row ("not donor data — no individual consent applies") with the GTEx row ("individual-level data was *not* consented for open release, which is why it sits behind dbGaP"). The second one explains the *mechanism*: the access tier is a consequence of the consent scope, not an arbitrary bureaucratic layer. Once you see that link, the two C10 rules in [Lecture 2 §2](../lecture-notes/02-data-ethics-and-public-data-sources.md) stop being arbitrary too.

**Access tier is a spectrum, and SRA proves it.** Seven of the eight rows have one tier. SRA has both, per study, and the row says so explicitly instead of averaging it into a single misleading word. When a dataset does not fit the column, widen the answer rather than rounding it — the acceptance criterion "no row says n/a without an explanation" is really a criterion about not flattening something you did not understand.

### What a grader is looking for

- **All 8 required datasets, all 9 columns filled.** An empty cell fails; "n/a" with the next cell explaining why is fine.
- **Version identifiers are specific.** An accession, a build number, a release number, or a date. The word "current" anywhere in that column is an automatic deduction.
- **Consent scope is a sentence, not a word.** "Open" is an access tier, not a consent scope. The grader wants evidence you found and read the data-use page — or an honest "not explicitly stated on the landing page" if you looked and it was not there, which is itself a real finding.
- **Citation lines are methods-section shaped.** Author, journal, volume, pages, year, DOI. "The 1000 Genomes paper" is not a citation.
- **A "Notes on access" section** that distinguishes registration from controlled access. These are different things with different effort profiles and different reasons for existing.

Acceptable variation: extra datasets are welcome and encouraged if you have a capstone idea. Different but equally specific version pins are fine — if you cite the 1000 Genomes 30× NYGC resequencing release instead of phase 3, that is a *better* answer as long as you say which one and why.

### Common wrong turns

**Writing "Open" in the consent column.** By far the most frequent error. Access tier answers "can I download it?"; consent scope answers "what did the people who provided it agree to?" A dataset can be openly downloadable and still consented only for non-commercial research. Those are different columns because they are different facts.

**Treating GTEx as one dataset.** It is two tiers with two different consent stories, and C10 uses only the public one. A row that says "GTEx — open" is wrong about half of GTEx and will lead you to a Week-7 tutorial you cannot follow.

**Registering for GISAID in Week 1 "to get it out of the way."** The exercise explicitly says not to — "do *not* register yet." Accounts have terms you accept on the day you accept them, and you should read those terms in the week you actually need the data, not skim them nine weeks early.

**Citing a URL as the citation.** A landing-page link is not a citation requirement. Nearly every one of these datasets has a "How to cite" page naming a specific paper, and several want two things: the consortium paper *and* the version DOI. GISAID wants three: two papers, the `EPI_SET` DOI, and an acknowledgement table.

### How to verify

Structural checks first — from your week-01 working directory:

```sh
grep -c '^| ' notes/week-01-data-inventory.md
```

Expected output (one header row, one separator row, eight data rows):

```text
10
```

Then the "no vague versions" check. This looks for the words the exercise forbids:

```sh
grep -in 'current release\|latest\|n/a' notes/week-01-data-inventory.md
```

Expected output: nothing, unless the hit is inside an explicit explanation you wrote on purpose.

Then check every URL resolves. On any machine with `curl`:

```sh
grep -o 'https\?://[^ >)]*' notes/week-01-data-inventory.md | sort -u | \
  while read -r u; do printf '%s %s\n' "$(curl -o /dev/null -s -w '%{http_code}' -L "$u")" "$u"; done
```

Expected output: every line starts with `200`. Anything starting `3`, `4`, or `5` is a link to fix — broken links are an explicit deduction in the mini-project rubric, and this one-liner is worth keeping for Week 12.

Then commit:

```sh
git add notes/week-01-data-inventory.md
git commit -m "Week 1 public-data inventory"
```

---

## Stretch goals

The [Week-1 README](../README.md) lists four optional pushes. None have a single right answer; here is what each is actually for and what "done" looks like.

**Read the introductions only of three recent Nature Methods papers.** What to notice, concretely: (1) the last paragraph of a Methods-paper introduction almost always states the gap the tool fills and the claim it makes — that is the sentence the rest of the paper defends; (2) version pinning appears in the *introduction* of good tool papers, not only the methods, because the comparison baselines have to be pinned to be meaningful; (3) the data-availability signal is usually a single sentence naming an archive and an accession. Done looks like three short notes in `notes/`, one per paper, each ending with the one sentence you think the paper is defending.

**Skim the 1000 Genomes consent documents.** The thing to come away with: donors consented to *broad public research use, including education*, with no restriction on who the researcher is — which is precisely why this dataset, and not your friend's lab VCF, is the right one for coursework. Note also what is *absent*: there is no per-project re-consent mechanism, which is the operational face of "consent cannot be revoked" from [Lecture 2 §1.3](../lecture-notes/02-data-ethics-and-public-data-sources.md).

**Read the NHGRI ELSI page end to end and bookmark two trusted sources.** Two defensible picks: the NHGRI ELSI programme itself (US-centric, well maintained, links out to primary sources) and the GA4GH Framework for Responsible Sharing (international, standards-oriented, the one your future colleagues will cite in a data-access discussion). Both are in [resources.md](../resources.md). Bookmark them somewhere you will actually find them in Week 11 when cancer-genomics data-use questions come up.

**Read a chapter of *Bioinformatics Data Skills*.** Chapter 1 is the highest-value Week-1 read: it is the reproducibility argument, made before you have any code to make irreproducible. Chapter 2 (Unix data tools) is the highest-value read if you want Week 2 to be easier.

---

*Continue to [challenges.md](./challenges.md).*
