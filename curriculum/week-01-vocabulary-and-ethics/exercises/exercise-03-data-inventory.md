# Exercise 3 — Public-Data Inventory

**Estimated time:** 30 minutes.

## Goal

Catalog the public datasets you will use across C10. By the end of this exercise you will have a single Markdown table that you can refer back to in every later week — and that will be the first thing your Week-12 capstone methods section cites.

This is *not* a creative exercise. It is a *due-diligence* exercise. The skill is reading a dataset's landing page carefully enough to know what you are actually working with.

## Output

A single Markdown file at `notes/week-01-data-inventory.md` in your week-01 working directory.

The file must contain a table with the following columns, one row per dataset:

| Column | What goes here |
|--------|----------------|
| Dataset | Name (e.g. "1000 Genomes Project, phase 3") |
| Provider | Hosting institution (e.g. "EMBL-EBI / IGSR") |
| URL | Canonical landing page |
| Data type | E.g. "human germline variation, VCF + BAM" |
| Access tier | "Open" / "Registration required" / "Controlled access" |
| Consent scope | One sentence summary of what donors consented to |
| Version / release | The specific version or release date you would cite |
| Citation requirement | The required citation in one line |
| C10 week | The week(s) you expect to use it |

## Required datasets

At minimum, your inventory must include all of the following:

1. **NCBI RefSeq** — used in Week 2 onwards for reference sequences.
2. **Ensembl reference genome (GRCh38)** — used in Weeks 5–8.
3. **1000 Genomes Project, phase 3** — used in Weeks 2 and 6.
4. **GTEx public summary tier** — used in Week 7.
5. **GISAID** — used in Week 9.
6. **UniProt** — used in Week 4.
7. **dbSNP** — used in Week 6.
8. **NCBI SRA** — used in Week 5.

You may add others you intend to use in your capstone.

## How to find each piece of information

For each dataset, visit the landing page and look for:

- **Consent scope:** usually in an "About" or "Data Use" page. If you cannot find one, this is itself information — note it as "Consent scope: not explicitly stated on landing page" and assume conservative use.
- **Version / release:** look for a release-notes link or a "current release" banner. Pin the exact version number or date.
- **Citation requirement:** every major dataset has a "How to cite" page. Some require multiple citations (e.g. GTEx wants both the consortium paper and the version DOI).
- **Access tier:** the landing page nearly always states this. If a tier requires registration, do *not* register yet — note that it does, and we will deal with it the week you actually use the data.

## Voice and citations

Write each row as a methods-section author would, not as a tweet would. Examples of the level of precision expected:

- Good: "1000 Genomes Project, phase 3 release v5a, March 2017."
- Bad: "1000 Genomes, current release."
- Good: "GRCh38.p14, GenBank assembly accession GCA_000001405.29."
- Bad: "the human reference genome."

## Acceptance criteria

- [ ] File exists at `notes/week-01-data-inventory.md`.
- [ ] All 8 required datasets are catalogued.
- [ ] Every row has all 9 columns filled in.
- [ ] No row says "n/a" without an explanation in the next column.
- [ ] At the bottom of the file, you have a short "Notes on access" section listing which datasets require registration and which require controlled-access applications.
- [ ] File is committed with a message like `Week 1 public-data inventory`.

## What this prepares you for

- **Week 2** will ask you to download a 1000 Genomes subset; you will know in advance what license terms apply.
- **Week 6** will ask you to use dbSNP; you will know in advance what version you are pinning to.
- **Week 12 capstone** has a "Data" section that is essentially this inventory, distilled to whatever subset you actually used.

The bioinformaticians who write reproducible papers do this work *up front*, not on the eve of submission. Build the habit now.

## When you are done

Commit the file. You are now done with Week-1 exercises. Move on to [Challenge 1](../challenges/challenge-01-reverse-complement.py).
