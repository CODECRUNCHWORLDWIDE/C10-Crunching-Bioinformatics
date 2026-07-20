# Challenge 1 — Taxonomy Classifier

> **Estimated time:** 120 minutes.
> **Goal:** Implement two BLAST-based taxonomy-classifier strategies — a top-hit classifier and a top-N majority-vote classifier — on a small labelled 16S rRNA dataset. Evaluate both on accuracy and confidence. Write up which one wins on this dataset and why, and identify at least one failure mode your top-hit classifier exhibits that the top-N classifier fixes (or vice versa).

Lecture 1 §4 discussed the seed-and-extend heuristic; Lecture 2 §7 introduced the biological identity thresholds (`> 98.7%` for species, `94–98.7%` for genus). This challenge ties them together. Given a BLAST hit table for a set of unknown 16S sequences, you need to assign a taxonomic label per query — and *that assignment is a classifier design choice*, not a single right answer.

## Background — the two classifier strategies

### Strategy A: Top hit classifier

For each query, take the BLAST hit with the lowest E-value. Look up its lineage. Assign the query that lineage's genus (or species, if you are feeling brave and `pident > 98.7%`).

```
classify_top_hit(query):
    hits = blast(query)
    hits = [h for h in hits if h.evalue < cutoff]
    if not hits: return None
    top = min(hits, key=lambda h: h.evalue)
    return lineage(top.subject_accession).genus
```

**Strengths:** Simple. Fast. Always gives an answer (when a hit is below cutoff). Works perfectly when the closest database entry is the right species.

**Weaknesses:** A single misannotated database entry produces a confident wrong call. A paralog with higher identity than the ortholog gives the wrong species. A chimeric reference produces a hit that is genuinely closer to the query than any true homolog — but pointing at a non-existent organism.

### Strategy B: Top-N majority classifier

For each query, take the top `N` hits (default `N = 5`). Look up the lineage of each. Take the most common genus among those `N`. If there is a tie or no genus dominates, return `None` (declined classification).

```
classify_top_n(query, n=5):
    hits = blast(query)
    hits = sorted(filter(lambda h: h.evalue < cutoff, hits),
                  key=lambda h: h.evalue)[:n]
    if len(hits) < n: return None
    lineages = [lineage(h.subject_accession) for h in hits]
    counts = Counter(l.genus for l in lineages)
    top_genus, top_count = counts.most_common(1)[0]
    if top_count >= ceil(n / 2):
        return top_genus
    return None
```

**Strengths:** A single misannotated entry cannot tip the answer if the next four are correct. Naturally produces a confidence estimate (`top_count / N`). Declines to answer rather than guess on ambiguous queries.

**Weaknesses:** A query whose closest reference is at 99.5% identity but whose 2nd through 5th hits are at 80% identity gets classified by the *less-similar* hits as much as by the most-similar one — which can be wrong. Slower (must process N hits per query, not 1). Returns `None` more often, which is fine if you care about precision and a problem if you care about recall.

## Task

Implement both classifiers, evaluate them on the labelled dataset below, and produce a writeup.

### 1. The labelled test set

Use this set of 8 unknown queries with ground-truth genus labels. The queries are real published 16S rRNA sequences from the NCBI 16S database whose accessions we have *deliberately withheld* — your classifier should not know them.

| Query ID | True genus | True species | Notes |
|----------|-----------|--------------|-------|
| `q01` | *Staphylococcus* | *S. aureus* | High-quality reference; should be easy. |
| `q02` | *Escherichia* | *E. coli* | High-quality reference; should be easy. |
| `q03` | *Bacillus* | *B. subtilis* | Multiple paralogous 16S copies in genome — risk of confused top-N. |
| `q04` | *Pseudomonas* | *P. aeruginosa* | Should be easy. |
| `q05` | *Mycobacterium* | *M. tuberculosis* | Slow-growing genus; reference is good but database has few entries. |
| `q06` | *Lactobacillus* | *L. acidophilus* | Recently split genus (now *Lactobacillus*, *Lactobacillaceae*, etc.) — taxonomy revisions are a real failure mode. |
| `q07` | *Clostridium* (*sensu stricto*) | *C. difficile* | *Clostridium* is heavily reclassified; older references may say *Peptoclostridium*. |
| `q08` | UNKNOWN | environmental | Deliberately ambiguous — top hits span two genera. A correct classifier returns `None` rather than guessing. |

You will produce the actual sequences for these queries by selecting a real 16S rRNA reference from each genus, downloading via `Bio.Entrez.efetch`, and saving as `data/queries.fasta`. The mini-project will give you 20 such queries; the challenge is the smaller-scale rehearsal.

### 2. Implement and run

Write `challenges/taxonomy_classifier.py` in your portfolio repo. The script should:

1. Read `data/queries.fasta` (8 queries) and `data/ground_truth.tsv` (the table above as TSV).
2. BLAST each query against either:
   - A pre-built local 16S database from Exercise 2 (preferred — fast), or
   - The NCBI `16S_ribosomal_RNA` database via `Bio.Blast.NCBIWWW` (slower, but no setup).
3. Implement `classify_top_hit(blast_hits, evalue_cutoff)` and `classify_top_n(blast_hits, n, evalue_cutoff)`.
4. For each query, fetch the taxonomy lineage of each hit (cached on disk so re-runs are instant).
5. Run both classifiers on all 8 queries.
6. Compute accuracy (correct calls / total queries), precision (correct calls / non-`None` calls), and recall (correct calls / total queries that had a correct answer to give) for each classifier.
7. Print a comparison table.

### 3. Write up

In `challenges/notes/classifier-comparison.md` (200–400 words), answer:

- Which classifier scored higher on **accuracy**? On **precision**? On **recall**?
- Were the rankings different across metrics? Why?
- On which queries did the two classifiers disagree? For each disagreement, which was right and what does the disagreement tell you about the failure mode of the other?
- For query `q08` (deliberately ambiguous), what did each classifier do? Which behavior is preferable in a real metagenomics pipeline, and why?
- Pick one of the failure modes from Lecture 1 §6 (low complexity, contamination, paralog/ortholog) and propose how you would *detect* it from the BLAST hit table alone (without ground-truth labels).

## Acceptance criteria

- `python challenges/taxonomy_classifier.py` runs without crashing on the 8 queries.
- Output includes a per-query table with columns `qseqid`, `true_genus`, `predicted_top_hit`, `predicted_top_n`, `top_hit_evalue`, `top_hit_pident`, `confidence_top_n`.
- A summary block at the end reports accuracy, precision, and recall for both classifiers as three-decimal-place floats.
- `notes/classifier-comparison.md` is 200–400 words and addresses all five prompts above.
- Taxonomy lookups are cached to `cache/taxonomy.json` so the second run is offline-only.

## Hints (do not peek for at least 30 minutes)

<details>
<summary>Hint 1 — How do I structure the BLAST cache?</summary>

For each query ID, cache the parsed top-10 hits to a per-query JSON file under `cache/blast/<qseqid>.json`. The structure can be the list of dicts your Exercise 3 reader produces. Two reasons: (1) iterating on the classifier should not re-run BLAST; (2) the cache is committable so a reviewer can re-run your evaluation without hitting NCBI.

</details>

<details>
<summary>Hint 2 — How do I look up taxonomy lineages efficiently?</summary>

`Bio.Entrez.efetch(db="taxonomy", id="1280", retmode="xml")` gives you the lineage for taxid 1280. To go from an accession to a taxid in the first place, use:

```python
handle = Entrez.esummary(db="nuccore", id=accession)
summary = Entrez.read(handle)[0]
taxid = summary["TaxId"]
```

Or, if your local BLAST database was built with `-taxid` / `-taxid_map`, the `staxid` column is in the tabular output (`-outfmt "6 std staxid"`) and no extra Entrez call is needed.

Either way, cache the lineage by taxid to `cache/taxonomy.json` so you only fetch each lineage once across all queries and across all runs.

</details>

<details>
<summary>Hint 3 — How do I parse the lineage string?</summary>

The NCBI taxonomy lineage looks like:

```
cellular organisms; Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcaceae; Staphylococcus; Staphylococcus aureus
```

Semicolon-separated, ordered domain → kingdom → ... → species. The genus is typically the second-to-last token; the species is the last (and often duplicates the genus name, e.g. "Staphylococcus aureus" is the species, "Staphylococcus" is the genus). Split on `"; "` and take `parts[-2]` for genus, `parts[-1]` for species. Edge case: not every record has an organism name at the species rank — handle the missing case gracefully.

</details>

<details>
<summary>Hint 4 — Precision vs recall on a multi-class classifier?</summary>

Define them on the "did we get the right genus call?" axis:

- **Accuracy** = (number of queries where the predicted genus equals the true genus) / (total queries).
- **Precision** = (number correct, *of the queries we attempted to classify*) / (number of queries we attempted, i.e. did not return `None`).
- **Recall** = (number correct, *of the queries that had a correct answer*) / (number of queries with a true label in the database).

For query `q08` (ambiguous, true_genus = UNKNOWN), a `None` prediction counts as *correct on accuracy* (we did not make a wrong call), *not counted* in precision (we did not make a call), and *correct on recall* (there was no true label to recover). Document your accounting in the writeup.

</details>

<details>
<summary>Hint 5 — What does a high-confidence wrong call look like?</summary>

The most informative wrong call is a `top_hit` prediction that has E-value < 1e-100 (high confidence on the surface) but turns out to be wrong because the database entry it hit was misannotated or chimeric. The top-N classifier may catch this because the other 4 hits would correctly point to the true genus, outvoting the misannotated top hit. This is the failure-mode example you should target for your "detection from the hit table alone" answer in the writeup: a top hit whose pident is much higher than any of the next 4 hits' pident, and whose subject taxid disagrees with the rest. That asymmetry is a contamination/misannotation signature.

</details>

## Stretch

If you finish under time and want more:

- Implement an LCA (lowest common ancestor) classifier as a third strategy. The LCA of the top `N` hits is the deepest taxonomic node shared by all `N`. If all `N` hits agree at genus, the LCA call is the genus. If they disagree at genus but agree at family, the LCA call is the family. This produces a *graceful degradation* with ambiguity instead of a binary present/absent. Compare its precision/recall against your top-N classifier on the same 8 queries.
- Add a `bit_score_normalized` column to your hit table: `bitscore / query_length`. Plot the distribution of normalized bit scores for the top hits across the 8 queries; the distribution should be bimodal (the 7 confident queries cluster high, the ambiguous `q08` is in the middle). This is the kind of diagnostic plot a methods paper publishes.
- Read the Kraken paper (Wood & Salzberg, *Genome Biol.* 15:R46, 2014) on k-mer-based LCA classification, and write a one-paragraph comparison of Kraken's approach vs your BLAST-based LCA. Kraken trades exact alignment for a precomputed k-mer→LCA index; it is ~1000x faster than BLAST at similar accuracy in metagenomics benchmarks.

## What you should be able to do after this

- Implement two competing classifier strategies on top of a BLAST hit table without confusing the algorithm and its inputs.
- Pick the right metric for the question you are asking (accuracy for headline reporting, precision for "are my calls right?", recall for "did I get the calls I could have gotten?").
- Diagnose at least one BLAST classifier failure mode from the hit table alone, without ground-truth labels.
- Defend a classifier choice in writing with reference to specific failure modes on specific queries.

---

*Submit by committing `challenges/taxonomy_classifier.py` and `challenges/notes/classifier-comparison.md` to your portfolio repo.*
