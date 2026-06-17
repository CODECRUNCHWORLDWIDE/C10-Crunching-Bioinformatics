# Lecture 1 — The Central Dogma in 90 Minutes

> **Duration:** ~2 hours of reading + hands-on.
> **Outcome:** You can explain DNA → RNA → protein on a whiteboard to a non-biologist, use the words *genome*, *chromosome*, *gene*, *transcript*, *exon*, *codon*, *nucleotide* correctly, and name the canonical bioinformatics file formats without conflating them.

If you only remember one thing from this lecture, remember this:

> **Bioinformatics is the work of treating biological molecules as strings.** A genome is a (very long) string over a four-letter alphabet. A protein is a string over a twenty-letter alphabet. Most of the field's algorithms — alignment, variant calling, phylogenetics — are string-processing algorithms that happen to be answering biological questions. Hold that frame, and the vocabulary stops being scary.

---

## 1. What bioinformatics is, and is not

Bioinformatics is the application of computational methods to **molecular biology data** — primarily sequences (DNA, RNA, protein), structures, and the high-throughput measurements that produce them. The defining feature is *scale*: a single human genome is about 3.2 billion base pairs, a single RNA-seq experiment produces tens of millions of reads, a single SARS-CoV-2 phylogeny in 2026 has millions of leaves.

Bioinformatics is **not**:

- *Computational biology* in the broad sense — that includes ecology, neuroscience modeling, mechanistic simulation. Bioinformatics is the sequence-and-omics slice.
- *Biostatistics* — there is overlap, but biostatistics is about study design and inference; bioinformatics is about pipelines and large biological datasets. Most working bioinformaticians use biostatistics as a tool.
- *Wet lab work* — we do not run pipettes. We process what the pipette-people produce.
- *Machine learning on health data* — there's a real overlap (single-cell, AlphaFold), but vanilla "ML for genomics" without grounding in the molecules tends to produce papers that don't replicate.

The bioinformatics field is **structurally open-source-first**: NCBI, Ensembl, EBI, Bioconductor, Biopython are all free, and the major journals require code and data deposition. This is why a Code Crunch track on bioinformatics is feasible — almost everything you need is free and public, which is also why the data-ethics conversation (Lecture 2) matters so much.

---

## 2. The molecules, briefly

You need three molecules in your head.

### DNA — *deoxyribonucleic acid*

Long double-stranded polymer. Each strand is a chain of **nucleotides**: one of four bases (**A**denine, **C**ytosine, **G**uanine, **T**hymine) attached to a sugar-phosphate backbone. The two strands are **complementary**: A pairs with T, C pairs with G. The two strands also run in opposite directions (anti-parallel) — one is 5' to 3', the other 3' to 5'. The numbers refer to carbon positions on the sugar; for our purposes the rule is: **sequences are conventionally written 5' to 3'**.

A DNA molecule is "just" a string over the alphabet `{A, C, G, T}`. Your genome is the concatenation of 23 such strings (the chromosomes), or 24 if you count Y.

### RNA — *ribonucleic acid*

A close cousin of DNA. Single-stranded, uses **U**racil instead of Thymine. Same alphabet otherwise: `{A, C, G, U}`. Several types: messenger RNA (**mRNA**) carries gene information out of the nucleus, ribosomal RNA (**rRNA**) and transfer RNA (**tRNA**) do the translation, plus many regulatory types we'll meet in Week 7.

### Protein

Strings over a 20-letter alphabet — the **amino acids**. Standard one-letter codes are `{A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V}`. Proteins fold into 3D structures and do most of the cell's mechanical and chemical work. AlphaFold predicts the fold; bioinformatics in 2026 still spends most of its time *before* fold-prediction — finding the gene, the transcript, the variants.

---

## 3. The central dogma

```
                 transcription              translation
        DNA   ───────────────────▶   mRNA   ───────────────────▶   protein
       (gene)                                                     (folded)
        ▲                                                            │
        │                                                            │
        └──── (no path back — proteins do not modify DNA letters) ───┘
```

Francis Crick proposed this in 1957. Read in plain English: **information flows from DNA to RNA to protein, and not back**. Modern biology has added many small exceptions (reverse transcription in retroviruses, epigenetic feedback, prion-like propagation) but the dogma is still the right mental model for the bulk of cellular biochemistry.

Two steps:

1. **Transcription** — an enzyme called RNA polymerase reads a stretch of DNA (a **gene**) and synthesizes a complementary RNA copy. In eukaryotes (including humans) the raw RNA is then **spliced** — non-coding stretches (**introns**) are cut out, leaving the coding stretches (**exons**) joined into a mature **mRNA**.
2. **Translation** — a **ribosome** reads the mRNA three letters at a time. Each three-letter **codon** specifies one **amino acid** (or a stop signal). The ribosome chains the amino acids together; the chain folds into a protein.

> **The vocabulary problem already starts here.** Outside of bioinformatics, people use "DNA" and "gene" and "genome" interchangeably. Inside the field they are different things:
>
> - **Genome** = the *entire* DNA content of one organism (~3.2 billion bp for human).
> - **Chromosome** = one of the (usually 23 in human) large DNA molecules into which the genome is partitioned.
> - **Gene** = one stretch of one chromosome, typically thousands of bases long, that codes for a product.
> - **Transcript** = one specific RNA produced from a gene; a single gene can produce multiple transcripts (isoforms) via alternative splicing.
> - **Nucleotide** = one *letter* of DNA or RNA. The smallest unit.

If you take one habit from this lecture, take this: when you read a paper or a tweet about genetics, mentally substitute the precise word every time someone says "DNA" or "gene." Most laypeople, and many press releases, get this wrong.

---

## 4. The genetic code

A **codon** is a 3-letter window over an mRNA. There are `4^3 = 64` codons. There are 20 standard amino acids plus a stop signal. So the code is **redundant** — multiple codons map to the same amino acid. This redundancy is mostly at the *third* base of the codon, the **wobble** position.

The standard codon table (most organisms — there are slight exotic variants in mitochondria and some unicellular eukaryotes):

```
       2nd base
       U          C          A          G
1st  ┌──────────┬──────────┬──────────┬──────────┐
 U   │ UUU Phe  │ UCU Ser  │ UAU Tyr  │ UGU Cys  │ U
     │ UUC Phe  │ UCC Ser  │ UAC Tyr  │ UGC Cys  │ C    3rd
     │ UUA Leu  │ UCA Ser  │ UAA STOP │ UGA STOP │ A    base
     │ UUG Leu  │ UCG Ser  │ UAG STOP │ UGG Trp  │ G
     ├──────────┼──────────┼──────────┼──────────┤
 C   │ CUU Leu  │ CCU Pro  │ CAU His  │ CGU Arg  │ U
     │ CUC Leu  │ CCC Pro  │ CAC His  │ CGC Arg  │ C
     │ CUA Leu  │ CCA Pro  │ CAA Gln  │ CGA Arg  │ A
     │ CUG Leu  │ CCG Pro  │ CAG Gln  │ CGG Arg  │ G
     ├──────────┼──────────┼──────────┼──────────┤
 A   │ AUU Ile  │ ACU Thr  │ AAU Asn  │ AGU Ser  │ U
     │ AUC Ile  │ ACC Thr  │ AAC Asn  │ AGC Ser  │ C
     │ AUA Ile  │ ACA Thr  │ AAA Lys  │ AGA Arg  │ A
     │ AUG Met* │ ACG Thr  │ AAG Lys  │ AGG Arg  │ G
     ├──────────┼──────────┼──────────┼──────────┤
 G   │ GUU Val  │ GCU Ala  │ GAU Asp  │ GGU Gly  │ U
     │ GUC Val  │ GCC Ala  │ GAC Asp  │ GGC Gly  │ C
     │ GUA Val  │ GCA Ala  │ GAA Glu  │ GGA Gly  │ A
     │ GUG Val  │ GCG Ala  │ GAG Glu  │ GGG Gly  │ G
     └──────────┴──────────┴──────────┴──────────┘
* AUG also signals "start of translation" in addition to coding for Methionine.
```

You will write a Python `translate()` function in [challenge-01](../04-challenges/challenge-01-reverse-complement.py) that uses this table. You will not memorize the table — *nobody does* — but you should be able to:

- Recognize that there are 64 codons and 20 amino acids.
- Identify the three stop codons (`UAA`, `UAG`, `UGA`) and the start codon (`AUG`).
- Explain why a single-nucleotide change at the third codon position is *often* silent.

---

## 5. Reference genomes

A **reference genome** is a community-curated example sequence for a species. It is not "the" genome of the species — every individual differs from it at millions of positions. It is the agreed-upon coordinate system everyone uses to talk about positions.

Key reference genomes you will meet in C10:

- **GRCh38** (also called hg38) — the current human reference. Genome Reference Consortium, build 38. Released 2013, patched continuously (`GRCh38.p14` was current as of late 2024).
- **T2T-CHM13** — the first truly *complete* (telomere-to-telomere) human reference, 2022. Adds ~200 Mb of sequence (mostly repetitive centromeres) that GRCh38 was missing.
- **GRCm39** — the current mouse reference.
- **NC_045512.2** — SARS-CoV-2 reference, the Wuhan-Hu-1 isolate. We use this in Week 9.

> **Cite the version.** "Mapped against GRCh38.p14 using bwa-mem2 v2.2.1" — not "mapped against the human genome." Position 123,456,789 on chromosome 22 in GRCh37 is a *different* base than position 123,456,789 in GRCh38. Mixing versions is the bioinformatics equivalent of mixing units.

---

## 6. The file formats — a tourist's overview

We will go deep on each of these in later weeks. For now: recognize them on sight.

### FASTA — plain sequences

```
>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
```

A header line beginning with `>` (the **FASTA chevron** — the visual signature of C10), then one or more lines of sequence. Used for reference genomes, gene sets, protein databases, query sequences for BLAST. Hundreds of millions of files in this format exist on NCBI alone.

### FASTQ — sequences with quality scores

```
@SRR062634.1 HWI-EAS110_103327062:6:1:1092:8469/1
GATTGTCAAGCTGAGTAAAGTCGCATCATGAACAATGGTAATTCCATATTC
+
@@CFFFFFGHHHHJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
```

Four lines per record: `@`-header, sequence, `+`-separator (sometimes repeats the header), quality string. Each quality character is the **Phred score** of the corresponding base, ASCII-encoded. Used for raw sequencer output. We parse FASTQ in Week 2.

### SAM / BAM — aligned reads

SAM is **Sequence Alignment/Map** — tab-separated text with a header section (`@HD`, `@SQ`, `@RG` lines) and one record per aligned read. Each record encodes the read, where it mapped, how well, and how it differs from the reference (the **CIGAR** string).

BAM is the binary, compressed version of SAM. Identical information, ~5× smaller, indexable with `samtools index`. We use both in Weeks 5 and 6.

### VCF — variants

**Variant Call Format**. Each row is one position where some individuals differ from the reference:

```
##fileformat=VCFv4.4
##reference=GRCh38
#CHROM  POS      ID         REF  ALT  QUAL   FILTER  INFO              FORMAT  SAMPLE1
22      16050075 rs141297151 A    G    998.0  PASS    AF=0.34;DP=42     GT:DP   0/1:42
```

Variants are how we represent the difference between an individual's genome and the reference. We go deep on VCF in Week 6.

### GFF / GTF — genome annotation

Tab-separated annotation files: "feature X (gene/exon/CDS) is located at this position on this chromosome on this strand." Used to describe gene structure. GTF is a stricter GFF variant. Critical for transcriptomics in Weeks 7–8.

### BED — genomic intervals

Bare-minimum genomic coordinates. Three required columns: chromosome, start, end. Used everywhere — peak files, region masks, capture targets.

---

## 7. The vocabulary problem

The same word can mean different things to a biologist and to an engineer working on the same dataset. A non-exhaustive list:

| Word | What an engineer might mean | What a biologist usually means |
|------|---------------------------|----------------------------|
| "DNA" | The whole genome | The DNA *polymer* (the molecule itself) |
| "Gene" | A coding sequence | A locus on a chromosome that produces one or more transcripts |
| "Expression" | An expression in code (`a + b`) | The amount of a gene's transcript being made |
| "Read" | An IO operation | One sequence record from a sequencer (~150 bp) |
| "Variant" | A version of software | A position where an individual differs from the reference |
| "Alignment" | Layout, formatting | Matching a read to its location on the reference |
| "Coverage" | Test coverage | How many reads map to a given position |
| "Pipeline" | A CI/CD pipeline | A linked series of bioinformatics tools |

This is not a joke and it is not a one-time confusion — it persists for *years* in mixed teams. The fix is the same fix as for any technical vocabulary mismatch: when you are not sure, slow down and ask "do you mean *X* in the biology sense or in the software sense?" Senior bioinformaticians do this constantly.

---

## 8. What sequencing actually produces

You will hear "Illumina reads," "Oxford Nanopore," "PacBio HiFi" — these are sequencing platforms. Each produces a different kind of raw data:

- **Illumina** — short reads (~150 bp), very accurate (>99.9% per base), cheap, ubiquitous. The default for most genomics work in 2026.
- **Oxford Nanopore (ONT)** — long reads (up to >1 Mb), portable hardware (the MinION fits in your hand), higher error rate (~95–99% per base depending on chemistry). Good for assembly and structural variants.
- **PacBio HiFi** — long reads (~15–25 kb) with very high accuracy (>99.9%). The current state of the art for *de novo* assembly.

For Week 1 you do not need to know any of this in detail. You will meet all three by Week 6.

---

## 9. What we will do over the next 12 weeks

A reminder of the arc, anchored in the vocabulary we just installed:

- **Weeks 2–3:** parse FASTA/FASTQ, write pairwise alignment from scratch.
- **Week 4:** run BLAST, build a small taxonomy classifier.
- **Weeks 5–6:** align reads to a reference, call variants, write a VCF.
- **Weeks 7–8:** quantify gene expression, run differential expression with DESeq2.
- **Week 9:** build a phylogenetic tree of SARS-CoV-2 sequences with bootstrap support.
- **Weeks 10–12:** publication-quality figures, reproducible pipelines with Snakemake, a capstone on a real biological question.

Every later week assumes you know the words we just defined. The mini-project for this week is a 1-page glossary *in your own words* and a personal data inventory. Both will come back in Week 12.

---

## 10. Self-check

Before moving on, you should be able to answer all of these out loud, without notes:

1. What are the four DNA bases? Which one is replaced in RNA, and by what?
2. What is the difference between a *gene* and a *transcript*?
3. How many codons are there? How many amino acids? Why doesn't this match?
4. Name three bioinformatics file formats and what each is for.
5. What does "GRCh38" refer to, and why does the version number matter?
6. Why is a sequencing read called a "read"?
7. Give one example of a word that means something different to an engineer and a biologist.

If any of those felt shaky, re-read the relevant section before continuing to Lecture 2.

---

*Continue to [Lecture 2 — Data Ethics and Public Data Sources](./02-data-ethics-and-public-data-sources.md).*
