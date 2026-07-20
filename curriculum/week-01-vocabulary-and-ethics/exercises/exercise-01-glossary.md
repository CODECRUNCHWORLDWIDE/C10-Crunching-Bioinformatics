# Exercise 1 — Glossary in Your Own Words

**Estimated time:** 30 minutes.

## Goal

Create a personal glossary of the Week-1 vocabulary, **in your own words**. Not Wikipedia's words. Not the lecture's words verbatim. Yours.

This is the exercise that distinguishes someone who skimmed the lecture from someone who can actually use these terms in a code review. If you cannot define a thing in your own words, you cannot reliably use it.

## Output

A single Markdown file at `notes/week-01-glossary.md` in your week-01 working directory. One entry per term, in this shape:

```markdown
### Nucleotide

A single letter of DNA or RNA. There are four DNA nucleotides (A, C, G, T)
and four RNA nucleotides (A, C, G, U — uracil replaces thymine in RNA).
A nucleotide is the smallest unit of a sequence; a whole human genome is
about 3.2 billion of them.

*Example:* the FASTA sequence `ATGC` contains four nucleotides.
*Easy to confuse with:* base pair (which is a pairing of two nucleotides
across the two strands of DNA, e.g. A-T or G-C).
```

Each entry must include:

1. **A plain-English definition** in 2–3 sentences. Pretend you are explaining to a software engineer who has never heard of biology.
2. **One concrete example** (a sequence, a number, a real gene name — whatever fits).
3. **One "easy to confuse with"** entry — the closest concept this term is *not*. This is the part that distinguishes a real glossary from a stack of definitions.

## Terms to define

Define all 20 of the following:

1. Nucleotide
2. Genome
3. Chromosome
4. Gene
5. Transcript
6. Exon
7. Intron
8. Codon
9. Amino acid
10. Protein
11. Reference genome
12. Variant (in the genomic sense)
13. SNP (single-nucleotide polymorphism)
14. FASTA
15. FASTQ
16. SAM / BAM
17. VCF
18. GFF / GTF
19. Read (in the sequencer sense)
20. Coverage

## Rules

- **Write each definition before re-checking the lecture.** If you cannot start the sentence, that is the signal that you should re-read that section. Closing the lecture and writing from memory is the *point*.
- **No copy-paste from Wikipedia, NCBI, or the lecture text.** If you find yourself reciting a sentence verbatim, paraphrase it.
- **Use the field's voice** — "is associated with," "is a region of," not "is what makes you tall." See the voice rules in [Lecture 2 §7](../lecture-notes/02-data-ethics-and-public-data-sources.md).
- **Length:** 2–4 sentences per definition. If yours is longer than 4 sentences, you are explaining mechanism — that belongs in the lecture notes, not the glossary.

## Acceptance criteria

- [ ] File exists at `notes/week-01-glossary.md`.
- [ ] All 20 terms are defined.
- [ ] Each entry has the three required parts (definition, example, easy-to-confuse-with).
- [ ] No two entries share an example.
- [ ] You did not copy from another source — and if you did consult one to check, you cited it at the bottom of the file in a "Sources I checked against" section.
- [ ] Committed to your week-01 repo with a message like `Week 1 glossary, first pass`.

## Why this matters

In Week 4 (BLAST) you will write a function called `parse_blast_xml`. In Week 6 (variants) you will write `apply_hard_filters`. If "BLAST" or "variant" is a fog-word for you, those functions will be fog too. The glossary is your one chance to install the vocabulary so it does not slow you down for the next 11 weeks.

This glossary is also a **deliverable for the Week-1 mini-project**, where it becomes a polished 1-page document. Treat this exercise as the first draft.

## When you are done

Commit the file and continue to [Exercise 2 — FASTA by hand](exercise-02-fasta-by-hand.py).
