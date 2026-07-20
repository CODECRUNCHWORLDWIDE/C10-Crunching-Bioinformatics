# Mini-Project — Personal Glossary and Public-Data Inventory

> Produce two short, polished documents that you will rely on for every later week of C10: a one-page glossary of bioinformatics vocabulary **in your own words**, and a one-page inventory of the public datasets you will use across the course.

This is the only mini-project in C10 that produces no executable pipeline. The point is to install the vocabulary and the data-source map *before* you ever run a tool. Every later week assumes both.

**Estimated time:** 7 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

A new public GitHub repository named `crunch-bio-portfolio-<yourhandle>` (this is also your portfolio repo for the rest of C10). It must contain, at minimum:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md
├── LICENSE                   (GPL-3.0 or MIT — pick one and commit to it)
├── .gitignore
└── week-01/
    ├── glossary.md           ~1 page, 21 polished entries
    ├── data-inventory.md     ~1 page, ≥8 datasets in a table
    └── reflection.md         (your homework problem 6 reflection, copied in)
```

By the end you will have a clean, public repo you can point a recruiter at — and a glossary and inventory you will *actually* refer back to in Weeks 2–12.

---

## Rules

- **You may** read the Week-1 lecture notes, the NCBI/Ensembl/EBI/etc. landing pages, free chapters of *Bioinformatics Data Skills*, and any other source you cite.
- **You may NOT** copy-paste definitions or table rows from any source. Every line must be in your own words. Citing a source you consulted (with a link) is required if you used it.
- **You may NOT** include any non-public dataset in the inventory. If you are unsure whether a dataset qualifies, the answer is to leave it out.
- The repo must be **public**. The license must be permissive.

---

## Acceptance criteria

- [ ] A public GitHub repo `crunch-bio-portfolio-<yourhandle>` exists.
- [ ] `README.md` at repo root introduces the portfolio in 1–2 paragraphs and links to `week-01/`.
- [ ] `LICENSE` is present (GPL-3.0 or MIT).
- [ ] `.gitignore` excludes `__pycache__/`, `.venv/`, `.DS_Store`, `.env`.
- [ ] `week-01/glossary.md`:
  - Contains 21 entries (20 from Exercise 1 + the "vocabulary problem" entry from Homework Problem 5).
  - Each entry has a plain-English definition, an example, and an "easy to confuse with" note.
  - At least 5 entries have *two* "easy to confuse with" notes.
  - Length: 1 printed page when rendered (roughly 800–1,200 words).
  - "Sources I consulted" section at the bottom with at least 3 links.
- [ ] `week-01/data-inventory.md`:
  - Contains the 9-column table specified in [Exercise 3](../exercises/exercise-03-data-inventory.md).
  - Includes all 8 required datasets.
  - Every row's `Version / release` is a specific identifier (date, accession, or version number) — no "current."
  - Notes on access section at the bottom.
- [ ] `week-01/reflection.md` is your Homework Problem 6 file, copied in.
- [ ] Commit history has at least 5 meaningful commits (not just one big "first commit").
- [ ] The repo passes a "fresh-eyes" check: a stranger landing on the README understands what the repo is for in under 30 seconds.

---

## Suggested order of operations

Build incrementally rather than trying to write the whole thing in one sitting.

### Phase 1 — Repo skeleton (~30 min)

1. Create the repo on GitHub. Initialize with a README and a `.gitignore` (Python template).
2. Clone locally. Add a `LICENSE` file (GPL-3.0 or MIT — your call).
3. Create `week-01/` directory.
4. First commit: `Initial portfolio skeleton`.

### Phase 2 — Move in your exercise drafts (~30 min)

1. Copy your Exercise 1 glossary draft to `week-01/glossary.md`.
2. Copy your Exercise 3 data inventory to `week-01/data-inventory.md`.
3. Commit: `Week 1 exercise drafts checked in`.

These are drafts. The next phases are where you *polish them into the actual deliverable*.

### Phase 3 — Polish the glossary (~2 hours)

1. Read each entry out loud. If a definition does not sound like something you would say in a code review, rewrite it.
2. Add the "vocabulary problem" entry from Homework Problem 5.
3. For at least 5 entries, add a second "easy to confuse with" note.
4. Add a "Sources I consulted" section at the bottom — *honest* about what you looked at. Three citation links minimum.
5. Open the rendered Markdown (GitHub preview or `grip` locally). The file should be a *visually* clean one-pager. If it looks dense, break the sections with `---` rules.
6. Commit: `Glossary v2: polished, sources cited`.

### Phase 4 — Polish the data inventory (~2 hours)

1. Re-visit each dataset's landing page. Confirm the version / release number you cited is still the current one.
2. For each dataset, write the citation line *exactly as you would put it in a methods section*. This is the part that pays off in Week 12.
3. Add a "Notes on access" section at the bottom listing which datasets require registration (e.g. GISAID) and which require controlled-access applications (e.g. dbGaP for individual-level GTEx).
4. Sanity-check: every URL works (you can paste them into a browser). Broken links are an automatic deduction.
5. Commit: `Data inventory v2: rechecked versions and citations`.

### Phase 5 — Write the repo README (~1 hour)

Your `README.md` at repo root should:

- Open with a one-sentence description: "*My C10 bioinformatics portfolio. Tracks weeks 1–12 of the [Code Crunch Crunching Bioinformatics track](https://github.com/CODE-CRUNCH-CLUB/C10-Crunching-Bioinformatics).*"
- Have a "Layout" section showing the directory tree.
- Have a "Week 1" section with a one-paragraph summary plus links to the glossary and inventory.
- Have a "How to read this repo" section: link order for a stranger arriving at it.
- Be polite and present-tense; not a tweet thread.

Commit: `Repo README, week-01 section`.

### Phase 6 — Copy in your reflection (~15 min)

Copy your Homework Problem 6 reflection into `week-01/reflection.md`. Commit: `Week 1 reflection`.

### Phase 7 — Final polish (~30 min)

- Skim every file. Fix typos.
- Confirm the `LICENSE` is correct.
- Confirm `.gitignore` is complete.
- Run `git log --oneline` and confirm at least 5 commits with meaningful messages.
- Push final.

---

## Reproducibility receipt

This mini-project produces no pipeline run, so the standard 5-row reproducibility receipt is replaced with a 3-row "documentation receipt" at the bottom of each polished file:

```
┌───────────────────────────────────────────────────────────┐
│  DOCUMENTATION                                            │
│                                                           │
│  Author:        <yourhandle>                              │
│  Last revised:  YYYY-MM-DD                                │
│  Sources:       see "Sources I consulted" below           │
└───────────────────────────────────────────────────────────┘
```

Mono. Always three rows. Include one of these blocks at the bottom of `glossary.md` and `data-inventory.md`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Glossary correctness | 25% | All 21 entries correct, in your own words, with examples and "easy to confuse with" |
| Inventory completeness | 25% | All 8 datasets, all 9 columns, every version pinned to a specific identifier |
| Voice and precision | 15% | Reads like a methods section, not a tweet. No determinism language. Cites versions. |
| Repo hygiene | 15% | README, LICENSE, .gitignore, ≥5 commits, no committed venv |
| Sources cited honestly | 10% | At least 3 sources listed per polished file; no copy-paste |
| Reflection | 10% | Honest, specific, addresses all four reflection prompts |

---

## What this prepares you for

- **Week 2** assumes you can read a FASTA / FASTQ file and know what each format is for — your glossary is the cheat sheet.
- **Week 6** assumes you know what a variant, a VCF, and a reference genome are — your glossary, again.
- **Week 12 capstone** has a "Data" section that is essentially this inventory, narrowed to whatever subset you actually used.

The bioinformaticians who write reproducible papers do this work *up front*. You are doing it up front for the same reason.

---

## Submission

When done:

1. Push your repo to GitHub with a public URL.
2. Make sure `README.md` introduces the repo and links to `week-01/`.
3. Make sure every URL in `data-inventory.md` works.
4. Open Week 2 only after you have committed. The Week-2 mini-project will live in `week-02/` in the same repo.

---

## Resources

- The reference quality bar: a senior bioinformatician's portfolio (search GitHub for `bioinformatics portfolio site:github.com` and skim a few — note what the polished ones look like).
- The C10 [`branding/BRAND.md`](../../../branding/BRAND.md) — the voice and layout you are matching.
- [resources.md](../resources.md) — every source you might cite.
