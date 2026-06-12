# C10 · Crunching Bioinformatics — Brand Guide

> **Voice:** lab-notebook precision. Empirical, reproducible, allergic to over-claim.
> **Feel:** the methods section of a clean Nature paper — measured, restrained, citation-heavy.

Extends the family brand. C10-specific overrides only.

---

## Identity

- **Full name:** Crunching Bioinformatics
- **Program code:** C10
- **Full title in copy:** *C10 · Crunching Bioinformatics*
- **Tagline (short):** Reproducible biology, in Python.
- **Tagline (long):** A free, open-source twelve-week bioinformatics track — from FASTA parsing to a published-quality phylogenetic tree, with reproducibility as a non-negotiable.
- **Canonical URL:** `codecrunchglobal.vercel.app/course-c10-bioinformatics`
- **License:** GPL-3.0

---

## Where C10 diverges from the family palette

Inherits Ink/Parchment/Gold. Adds one **Lab Green** for "alive" / "biological data" semantics:

| Role | Name | Hex | Use |
|------|------|-----|-----|
| Accent | Lab Green | `#15803D` | The C10 mark, "successful run" indicators, FASTA chevrons |
| Accent deep | Lab Green deep | `#14532D` | Hover, eyebrows |
| Accent soft | Lab Green soft | `#BBF7D0` | Subtle row highlight in expression tables |

```css
:root {
  --lab-green:       #15803D;
  --lab-green-deep:  #14532D;
  --lab-green-soft:  #BBF7D0;
}
```

### Typography

EB Garamond display, Lora body. **JetBrains Mono for every sequence, gene symbol, accession number, organism abbreviation, and command-line invocation.** Mono is the "data" face — anything you'd otherwise wrap in `<code>` gets mono.

---

## Recurring page elements

### The FASTA-style sequence header

The visual signature of C10 is the FASTA `>` chevron rendered in Lab Green:

```
>NC_045512.2 SARS-CoV-2 isolate Wuhan-Hu-1, complete genome
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAAT
CTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACG
```

Use this header style for any biological-sequence example. The chevron is Lab Green; the rest is mono Ink.

### The "reproducibility receipt"

Every analysis exercise includes a small box certifying reproducibility:

```
┌─────────────────────────────────────────────────────────┐
│  REPRODUCIBILITY                                        │
│                                                         │
│  Data source:   1000 Genomes phase 3, chr22, 50 samples │
│  Pipeline:      snakemake v8.4.2                        │
│  Container:     ghcr.io/code-crunch-club/c10:latest      │
│  Command:       snakemake --use-conda -j 4              │
│  Wall time:     ~14 min (on M1, 16 GB)                  │
└─────────────────────────────────────────────────────────┘
```

Mono. Always five rows. The Wall-time row is non-optional — students get used to citing it.

---

## Voice rules

- **Cite the version.** "Biopython 1.83" — not "the latest Biopython."
- **Cite the dataset by accession.** "GSE52778" — not "an RNA-seq dataset."
- **Distinguish biological claim from statistical claim.** "Gene X is differentially expressed (FDR < 0.05, log2FC > 1)" — not "Gene X is important."
- **Acknowledge biological uncertainty.** "Suggests a role in" — not "proves."
- **No genetic-determinism language.** "Variant X is associated with" — not "variant X causes."
- **Methods sections are not optional.** Every analysis page includes one.

---

## Course page conventions

The course page (`course-c10-bioinformatics.html`, future) uses a slightly Lab-Green-tinted parchment with a FASTA-style sequence as the hero element. The 12-week ladder is rendered as a "lab-notebook table of contents" — each phase numbered like an experimental day.

---

*GPL-3.0. Fork freely.*
