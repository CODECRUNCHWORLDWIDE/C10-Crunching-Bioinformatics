# Mini-Project — BLAST-Driven Taxonomy Classifier for 20 Unknown Sequences

> Build a reproducible, batch-mode taxonomy classifier that takes ~20 unknown 16S rRNA sequences, BLASTs each against a local copy of the NCBI `16S_ribosomal_RNA` database, fetches the taxonomy lineage of every top hit via `Bio.Entrez`, applies a top-hit-with-LCA-fallback classifier, and produces a per-query classification with E-value/identity/confidence, a confusion matrix against ground-truth labels, and an honest precision-recall write-up.

This is the first C10 mini-project that produces a **methods-section-quality classifier with measured performance**, not just a single-query demonstration. By the end of it you will have a `classify_taxonomy.py` script you can point a recruiter at, a results directory with a confusion matrix and per-query confidence, and a write-up that defends each classifier choice and names the failure modes you observed.

**Estimated time:** 8 hours (split across Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-04/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 4 section)
└── week-04/
    ├── README.md                   one-page report (~800 words)
    ├── run.sh                      one-command reproduction script
    ├── env.yml                     conda environment file pinning all tool versions
    ├── data/
    │   ├── queries.fasta           20 unknown 16S rRNA sequences
    │   ├── ground_truth.tsv        per-query expected genus and species
    │   └── 16S_ribosomal_RNA.*     downloaded BLAST DB files (large, in .gitignore)
    ├── classify_taxonomy.py        the classifier pipeline
    ├── parse_blast.py              the BLAST-output reader (from Exercise 3)
    ├── fetch_taxonomy.py           the Entrez lineage fetcher (cached)
    ├── cache/
    │   ├── blast/                  per-query JSON cache of parsed hits
    │   └── taxonomy.json           taxid → lineage cache
    └── results/
        ├── classifications.tsv     one row per query, predicted vs true
        ├── confusion_matrix.tsv    predicted-genus by true-genus counts
        ├── confusion_matrix.png    heatmap
        └── per_query_confidence.png  bar plot of classification confidence
```

By the end you will have a clean, reproducible Week 4 directory you can point a recruiter at — and `classify_taxonomy.py` is the kind of pipeline that opens conversations with working bioinformaticians and metagenomics shops.

---

## The dataset

You will work with **20 unknown 16S rRNA sequences** that represent the broad diversity of cultivable bacteria plus a handful of designed-to-fail edge cases. Source the queries by fetching one published 16S reference per organism in the table below, using `Bio.Entrez.efetch(db="nuccore", id=<accession>, rettype="fasta")`. Save the concatenated FASTA as `data/queries.fasta`. Save the ground-truth labels as `data/ground_truth.tsv`.

| Query ID | NCBI accession | True genus | True species | Notes |
|----------|----------------|-----------|--------------|-------|
| `q01` | `NR_117741.1` | *Staphylococcus* | *S. aureus* | Easy. High-quality reference. |
| `q02` | `NR_074549.1` | *Escherichia* | *E. coli* | Easy. The textbook reference. |
| `q03` | `NR_119213.1` | *Bacillus* | *B. subtilis* | Multiple paralogous 16S copies — test top-N robustness. |
| `q04` | `NR_117150.1` | *Pseudomonas* | *P. aeruginosa* | Easy. |
| `q05` | `NR_118889.1` | *Mycobacterium* | *M. tuberculosis* | Slow-growing genus; sparse reference coverage. |
| `q06` | `NR_115405.1` | *Lactobacillus* | *L. acidophilus* | Recently split genus — taxonomy revision risk. |
| `q07` | `NR_113132.1` | *Clostridium* | *C. difficile* | Heavily reclassified; some references say *Clostridioides difficile* (the current accepted name). |
| `q08` | `NR_117044.1` | *Salmonella* | *S. enterica* | Close cousin to *E. coli*; test genus discrimination. |
| `q09` | `NR_113236.1` | *Streptococcus* | *S. pneumoniae* | Easy. |
| `q10` | `NR_113308.1` | *Helicobacter* | *H. pylori* | Easy. |
| `q11` | `NR_117920.1` | *Listeria* | *L. monocytogenes* | Closely related to *Bacillus*. |
| `q12` | `NR_103973.1` | *Vibrio* | *V. cholerae* | Easy. |
| `q13` | `NR_113617.1` | *Neisseria* | *N. meningitidis* | Easy. |
| `q14` | `NR_044993.1` | *Campylobacter* | *C. jejuni* | Easy. |
| `q15` | `NR_115365.1` | *Yersinia* | *Y. pestis* | Closely related to other Enterobacteriaceae; test discrimination. |
| `q16` | `NR_115991.1` | *Bordetella* | *B. pertussis* | Easy. |
| `q17` | `NR_115529.1` | *Treponema* | *T. pallidum* | Spirochete — distinctive 16S. |
| `q18` | `NR_103934.1` | *Borrelia* | *B. burgdorferi* | Another spirochete; tests cross-genus separation. |
| `q19` | `NR_117794.1` | *Haemophilus* | *H. influenzae* | Easy. |
| `q20` | `NR_044878.1` | *Rickettsia* | *R. rickettsii* | Obligate intracellular; reference is sparse. |

If any of these accessions have been retired or superseded by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt. The classifier should be robust to a few accession swaps.

### Why 16S and not whole genomes

16S rRNA is the standard bacterial taxonomy marker for three reasons: (1) it is universally present in bacteria and archaea, (2) it is small enough to sequence cheaply (~1.5 kb), and (3) the curated NCBI `16S_ribosomal_RNA` database is small (~25,000 entries) and queryable in seconds. For whole-genome identification, ANI (average nucleotide identity) is the equivalent metric but requires the full genome and a more expensive pipeline. Week 5 will introduce read-alignment-based identification at the whole-genome level.

---

## Rules

- **You may** use BLAST+ 2.15, Biopython 1.83, pandas, matplotlib, and the standard library.
- **You may** consult Lectures 1 and 2, the BLAST+ user manual, Biopython's tutorial chapter 7, and your Week 4 exercises and challenge.
- **You may NOT** copy a pre-written classifier from the internet. The point is to *build* the pipeline. Reading scikit-bio's source or DIAMOND's documentation for inspiration is fine; copy-pasting their classifier is not.
- **You must** cache all BLAST and Entrez results to disk. A second run of `bash run.sh` should not hit NCBI for any data it has already fetched.
- The repo must be **public** and the mini-project must be reproducible from `run.sh` on a fresh checkout, given the environment file. The only network access on a fresh run should be the initial database download.

---

## Acceptance criteria

- [ ] `week-04/classify_taxonomy.py` exports a function `classify(queries_fasta, db_prefix, *, evalue_cutoff=1e-50, top_n=5) -> pd.DataFrame` returning a DataFrame indexed by query ID with columns `predicted_genus`, `predicted_species`, `top_hit_accession`, `top_hit_evalue`, `top_hit_pident`, `top_hit_align_length`, `n_hits_above_cutoff`, `confidence` (a float in [0, 1] derived from the top-N consensus).
- [ ] The classifier implements a **two-stage decision rule**:
  1. If the top hit has `pident > 98.7` and `evalue < 1e-50`, emit the top hit's genus and species with `confidence = 1.0`.
  2. Otherwise, compute the majority genus among the top `N = 5` hits below `evalue_cutoff`. If at least `ceil(N / 2)` agree, emit that genus, with the species marked as "uncertain" and `confidence = top_count / N`. If no genus reaches the majority threshold, emit `predicted_genus = None` with `confidence = 0`.
- [ ] `week-04/README.md` is a one-page (≤ 1,000 word) report containing:
  - One-sentence description of the dataset, the question, and the classifier's two-stage decision rule.
  - Methods section in C10 voice: every tool pinned ("BLAST+ 2.15.0", "Biopython 1.83"), every parameter explicit ("`blastn -task blastn -evalue 1e-50 -max_target_seqs 10 -max_hsps 1`", "top-N = 5", "majority threshold = ceil(N/2) = 3").
  - Quantitative findings: per-classifier accuracy, precision, recall. Confusion matrix as an inline table. A per-query confidence histogram described in prose.
  - A failure-modes section listing every query where the classifier returned the wrong answer or declined to classify, with one-sentence diagnoses (paralog? sparse reference? taxonomy revision?).
  - A reproducibility receipt block.
- [ ] `week-04/run.sh` runs end-to-end on a fresh clone:
  - Activates the conda env from `env.yml`.
  - Downloads the `16S_ribosomal_RNA` BLAST database via `update_blastdb.pl` if not already cached on disk.
  - Fetches the 20 query sequences from NCBI via `Bio.Entrez.efetch` if not already on disk.
  - Runs `classify_taxonomy.py` to BLAST, classify, and produce results.
  - Renders `confusion_matrix.png` and `per_query_confidence.png` via matplotlib.
- [ ] `week-04/env.yml` pins every tool to an exact version:
  ```yaml
  name: c10-week-04
  channels:
    - conda-forge
    - bioconda
  dependencies:
    - python=3.11
    - numpy=1.26.4
    - biopython=1.83
    - pandas
    - matplotlib
    - blast=2.15
    - pip
  ```
- [ ] `week-04/results/classifications.tsv` is a TSV with one row per query and columns: `qseqid`, `true_genus`, `true_species`, `predicted_genus`, `predicted_species`, `confidence`, `top_hit_accession`, `top_hit_pident`, `top_hit_evalue`, `top_hit_bitscore`, `n_hits_above_cutoff`, `correct` (bool).
- [ ] `week-04/results/confusion_matrix.tsv` is a TSV with predicted-genus rows and true-genus columns, integer counts.
- [ ] `week-04/results/confusion_matrix.png` is a matplotlib heatmap with axes labelled.
- [ ] `week-04/results/per_query_confidence.png` is a bar plot of per-query confidence colored by `correct`.
- [ ] The repo passes a fresh-clone test: `git clone`, `cd week-04`, `bash run.sh` reproduces everything in `results/` (modulo a ±1 difference in counts if NCBI's database has been updated).

---

## Suggested order of operations

### Phase 1 — Environment setup (~30 min)

1. Create `week-04/env.yml` (see acceptance criteria for the exact contents).
2. `conda env create -f week-04/env.yml`. Activate it. Confirm: `python -c "import Bio; print(Bio.__version__)"` → `1.83`, and `blastn -version` → `blastn: 2.15.0+`.
3. Commit: `Week 4 env.yml`.

### Phase 2 — Download the database and queries (~30 min)

1. Write `scripts/fetch_database.sh`:
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail
   mkdir -p data
   cd data
   update_blastdb.pl --decompress 16S_ribosomal_RNA
   ```
   Run it once. The 16S_ribosomal_RNA database is small (~25 MB) so the download is fast.
2. Write `scripts/fetch_queries.py` that takes the table in this README and fetches each accession via `Bio.Entrez.efetch`, then writes `data/queries.fasta` (multi-FASTA, 20 records) and `data/ground_truth.tsv`. Set `Bio.Entrez.email`. Cache: if `data/queries.fasta` already exists, skip.
3. Commit: `Database fetched, queries staged`.

### Phase 3 — Port the parsing layer from Exercise 3 (~30 min)

1. Copy `parse_blast.py` from your Exercise 3 implementation into `week-04/`.
2. Confirm the four functions (`read_tabular`, `read_xml`, `filter_top_hit_per_query`, `attach_lineage_stub`) still pass their self-tests when invoked from this directory.
3. Replace `attach_lineage_stub` with a real implementation that calls `fetch_taxonomy.lookup_lineage(taxid)` (next phase).
4. Commit: `parse_blast.py copied and lineage call upgraded`.

### Phase 4 — Build the taxonomy fetcher (~1 hour)

1. Write `fetch_taxonomy.py` that exports two functions:
   - `accession_to_taxid(accession: str) -> int` — calls `Bio.Entrez.esummary(db="nuccore", id=accession)` and reads the `TaxId` field. Cache by accession to `cache/accession_to_taxid.json`.
   - `lookup_lineage(taxid: int) -> dict` — calls `Bio.Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")` and returns `{"scientific_name": ..., "rank": ..., "lineage": [list of dicts with rank+name]}`. Cache by taxid to `cache/taxonomy.json`.
2. Both functions must be **idempotent and offline-safe** on cache hit. Test by deleting the cache, running once (network), running again (offline — disconnect Wi-Fi to confirm).
3. Commit: `Taxonomy fetcher with disk cache`.

### Phase 5 — Wire up the classifier (~2 hours)

1. Write `classify_taxonomy.py`. Outline:
   ```python
   def classify(queries_fasta: Path, db_prefix: Path,
                *, evalue_cutoff=1e-50, top_n=5) -> pd.DataFrame:
       hits_per_query = run_blastn_batch(queries_fasta, db_prefix, evalue_cutoff)
       lineages_per_hit = {a: lookup_lineage(accession_to_taxid(a))
                           for a in distinct_accessions(hits_per_query)}
       rows = []
       for qseqid, hits in hits_per_query.items():
           rows.append(classify_one(qseqid, hits, lineages_per_hit,
                                    evalue_cutoff, top_n))
       return pd.DataFrame(rows).set_index("qseqid")

   def classify_one(qseqid, hits, lineages, evalue_cutoff, top_n) -> dict:
       above_cutoff = [h for h in hits if h["evalue"] < evalue_cutoff]
       if not above_cutoff:
           return {"qseqid": qseqid, "predicted_genus": None, ...}
       top = above_cutoff[0]  # sorted by evalue ascending
       # Stage 1: top-hit, high-confidence rule.
       if top["pident"] > 98.7 and top["evalue"] < 1e-50:
           lineage = lineages[top["sseqid"]]
           return {"qseqid": qseqid,
                   "predicted_genus": genus_from_lineage(lineage),
                   "predicted_species": species_from_lineage(lineage),
                   "confidence": 1.0, ...}
       # Stage 2: top-N majority.
       window = above_cutoff[:top_n]
       genera = [genus_from_lineage(lineages[h["sseqid"]]) for h in window]
       counts = Counter(genera)
       top_genus, top_count = counts.most_common(1)[0]
       if top_count >= math.ceil(top_n / 2):
           return {"qseqid": qseqid,
                   "predicted_genus": top_genus,
                   "predicted_species": "uncertain",
                   "confidence": top_count / top_n, ...}
       return {"qseqid": qseqid, "predicted_genus": None, "confidence": 0, ...}
   ```
2. Make sure every BLAST call is batched (single `blastn` invocation on the multi-query FASTA, not 20 separate calls) and the result is cached.
3. Commit: `classify_taxonomy.py end to end`.

### Phase 6 — Evaluate and produce the confusion matrix (~1 hour)

1. Add an `evaluate(classifications: pd.DataFrame, ground_truth: pd.DataFrame) -> dict` function that returns `{"accuracy": ..., "precision": ..., "recall": ..., "confusion_matrix": pd.DataFrame}`.
2. Definitions:
   - **Accuracy** = (correct calls + correct `None`-on-ambiguous) / total queries. For your 20 queries with non-ambiguous ground truth, this is just correct / 20.
   - **Precision** = correct calls / non-`None` predictions.
   - **Recall** = correct calls / queries with a true label (all 20 here).
3. Write the confusion matrix to `results/confusion_matrix.tsv` and render a heatmap to `results/confusion_matrix.png`.
4. Write per-query confidence bar plot to `results/per_query_confidence.png`, colored green for correct calls and red for incorrect.
5. Commit: `Evaluation + confusion matrix`.

### Phase 7 — Write the report (~2 hours)

Open `week-04/README.md`. Structure:

```
# Week 4 — BLAST-driven taxonomy classifier for 20 unknown 16S sequences

## Dataset and biological question
- 20 unknown 16S rRNA sequences, accessions listed in data/ground_truth.tsv
- Question: can a top-hit + top-N-fallback BLAST classifier achieve > 90%
  accuracy on a small, taxonomically diverse panel of bacterial 16S queries?

## Methods
[200 words. Local NCBI `16S_ribosomal_RNA` database (~25k entries) via
update_blastdb.pl. `blastn -task blastn -evalue 1e-50 -max_target_seqs 10`.
Taxonomy lookups via Bio.Entrez.efetch with disk cache. Two-stage classifier:
top-hit at pident > 98.7% else top-N=5 majority vote at majority threshold=3.
Implementation in pure Biopython 1.83 + pandas; no external classifier library.]

## Findings
- Accuracy: 18/20 = 0.900
- Precision: 18/19 = 0.947 (one query declined classification)
- Recall: 18/20 = 0.900
- Confusion matrix: [inline table]
- Per-query confidence: see results/per_query_confidence.png

## Failure modes observed
- q07 (Clostridium): predicted "Clostridioides" — taxonomy revision; database
  uses the current accepted name.
- q06 (Lactobacillus): predicted with confidence 0.6 — top hits split between
  two recently-split genera (Lactobacillus / Lactobacillaceae).
- q20 (Rickettsia): declined classification — only 2 hits above cutoff in
  the sparse reference; top-N rule rejected.

## Reproducibility
[the receipt block, below]
```

Aim for **clear, quantitative, undramatic**. C10 voice. No "the classifier looks good"; give numbers.

Commit: `Week 4 report v1`.

### Phase 8 — Reproducibility receipt and polish (~30 min)

Add the receipt block at the bottom of `week-04/README.md`:

```
+--------------------------------------------------------------------+
|  REPRODUCIBILITY                                                   |
|                                                                    |
|  Data source:   NCBI 16S_ribosomal_RNA (downloaded YYYY-MM-DD).    |
|                  20 query accessions in data/ground_truth.tsv.     |
|  Pipeline:      bash run.sh                                        |
|  Container:     conda env (env.yml pinned to BLAST+ 2.15.0,        |
|                  Biopython 1.83, Python 3.11)                      |
|  Command:       bash run.sh                                        |
|  Wall time:     ~45 s (M1, 16 GB RAM, network for first run only)  |
+--------------------------------------------------------------------+
```

Sanity-check: every file path in the README is correct, every URL works, `bash run.sh` succeeds on a fresh checkout. Commit: `Mini-project reproducibility receipt`.

### Phase 9 — Update the repo root README (~15 min)

Add a Week 4 section to your portfolio repo's top-level `README.md` linking to `week-04/README.md` with a one-paragraph summary. Commit: `Week 4 entry in portfolio README`.

---

## Rubric

| Criterion | Weight | What "great" looks like |
|----------|-------:|-------------------------|
| Correctness | 25% | Per-query predictions match the ground-truth genus for at least 17 of 20 queries. Honest accounting of the failures. |
| Reproducibility | 20% | `bash run.sh` works on a fresh clone. Versions pinned. Data cached. |
| Code quality | 15% | `classify_taxonomy.py` has module docstrings, structured return types, input validation, inline tests. Caching is correct and observable (delete cache, observe network fetch; restore cache, observe instant). |
| Quantitative report | 15% | Every claim has a number. C10 voice throughout. Failure modes named per query. |
| Confusion matrix | 10% | The matrix is correct, the heatmap is labelled, the off-diagonal cells are explained in the writeup. |
| Voice and precision | 10% | Reads like a methods section. No "fast" without seconds, no "good" without numbers. |
| Plot quality | 5% | `confusion_matrix.png` and `per_query_confidence.png` are labelled, captioned, and biologically interpretable. |

---

## What this prepares you for

- **Week 5** introduces read-alignment-based taxonomy from short reads (FastQ → BAM → coverage-based species inference). The classifier-evaluation framework (precision/recall/confusion-matrix) carries directly over.
- **Week 9** revisits classification in the phylogenetic-tree setting (tree-based placement via `pplacer` or `EPA-ng`). Your top-N consensus and LCA fallback are simpler cousins of tree-placement methods.
- **Week 11** asks you to convert this pipeline to Snakemake. The current `run.sh` is the rehearsal: every step is a `python script.py` invocation that can become a Snakemake rule.
- **Week 12 capstone** will likely include a classifier of some kind. The honest evaluation methodology you build here is the template.

---

## Common pitfalls

A short list, from instructor experience:

- **Forgetting `-max_hsps 1`.** A single (query, subject) pair can produce multiple HSPs (when the alignment has internal gaps that split it). Without `-max_hsps 1`, your "top 5 hits" can be 5 HSPs from the *same* subject, defeating the top-N consensus. Set it.
- **Hard-coded NCBI database paths.** `update_blastdb.pl` puts the database files in the current working directory by default. If `run.sh` changes directories midway through, the second BLAST call cannot find the database. Use absolute paths or `pushd`/`popd`.
- **Caching the wrong things.** Cache BLAST results per *query batch* (one cache file per `queries.fasta` checksum), not per individual query — otherwise iterating on the classifier code means re-running BLAST. Caching the parsed JSON is even better.
- **Treating `None` as a wrong call.** A classifier that *declines* to classify an ambiguous query is correctly behaving; it should not count against accuracy. Be explicit in your definitions of precision and recall. The challenge writeup walks through this; do not skip it.
- **Reporting only accuracy.** Accuracy is a misleading metric when the classifier can decline to answer. Always report all three: accuracy, precision, recall. The writeup should call out the difference.

---

## Submission

When done:

1. Confirm `bash run.sh` works on a fresh clone.
2. Confirm `week-04/README.md` renders cleanly on GitHub.
3. Confirm `results/classifications.tsv` is valid TSV with the right columns.
4. Confirm `results/confusion_matrix.png` and `results/per_query_confidence.png` are committed and are meaningful figures.
5. Push.
6. Open Week 5 — only after the report is committed and the failure-modes section names every wrong call.

---

## Resources

- [NCBI BLAST+ user manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/) — the command-line reference.
- [Bio.Blast tutorial chapter](https://biopython.org/DIST/docs/tutorial/Tutorial.pdf) — Chapter 7.
- [NCBI E-utilities help](https://www.ncbi.nlm.nih.gov/books/NBK25497/) — `efetch`, `esummary`, `elink`.
- [NCBI 16S_ribosomal_RNA database documentation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs#16S) — the curated 16S database used in this mini-project.
- [Lecture 1 — How BLAST actually works](../lecture-notes/01-how-blast-actually-works.md)
- [Lecture 2 — Running BLAST locally and via NCBI](../lecture-notes/02-running-blast-locally-and-via-ncbi.md)
- [Challenge 1 — Taxonomy classifier (the rehearsal at N=8)](../challenges/challenge-01-taxonomy-classifier.md)
- [resources.md](../resources.md) — full week resource list.
