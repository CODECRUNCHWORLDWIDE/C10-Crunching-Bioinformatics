# Week 4 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (BLAST+ 2.15.0, Biopython 1.83, Python 3.11), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Altschul, Gish, Miller, Myers, Lipman (1990)** — the original BLAST paper. *Journal of Molecular Biology* 215:403–410. Free PDF widely mirrored; the official copy is paywalled at ScienceDirect but the publisher's open-access policy permits archival:
  <https://www.sciencedirect.com/science/article/pii/S0022283605803602>
  Mirror (often easier): <https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/altschul1990.pdf>
- **Altschul et al. (1997)** — BLAST 2.0, *gapped* BLAST and PSI-BLAST, introducing the two-hit rule. *Nucleic Acids Research* 25:3389. Free full text:
  <https://academic.oup.com/nar/article/25/17/3389/1061651>
- **Karlin & Altschul (1990)** — the K-A E-value statistical framework. *PNAS* 87:2264. Free full text:
  <https://www.pnas.org/doi/10.1073/pnas.87.6.2264>
- **NCBI BLAST+ user manual** — the definitive reference for the command-line tools (`makeblastdb`, `blastn`, `blastp`, `tblastn`, `blastdbcmd`):
  <https://www.ncbi.nlm.nih.gov/books/NBK279690/>
- **Biopython Tutorial and Cookbook — Chapter 7 (BLAST)** — `Bio.Blast.NCBIWWW`, `Bio.Blast.NCBIXML`, `Bio.Blast.Applications`:
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **NCBI Entrez Programming Utilities Help (E-utils)** — the public API: `efetch`, `esearch`, `elink`, `summary`:
  <https://www.ncbi.nlm.nih.gov/books/NBK25497/>

## BLAST databases (the public ones you will query)

- **`nt`** — the full GenBank non-redundant nucleotide database. ~500 GB. Public, queryable via NCBI web BLAST or downloadable via `update_blastdb.pl --decompress nt`:
  <https://ftp.ncbi.nlm.nih.gov/blast/db/>
- **`nr`** — the non-redundant protein database. ~300 GB. Same access pattern.
- **`refseq_rna`** — the curated RefSeq RNA database. Smaller and cleaner than `nt`; preferred for transcript-level identification work.
- **`refseq_protein`** — the curated RefSeq protein database. Same idea for proteins.
- **`16S_ribosomal_RNA`** — the curated 16S rRNA database from NCBI. ~25 MB. The canonical reference for bacterial taxonomy from a 16S amplicon — and the database you will use in the Week 4 mini-project.
- **`SwissProt`** — the manually curated portion of UniProt. ~600,000 entries, much smaller than `nr`, much higher quality. Available via NCBI BLAST or directly from UniProt:
  <https://www.uniprot.org/uniprotkb?query=*&fields=accession,id,protein_name,organism_name>

## BLAST+ tool reference (the command-line surface)

Cheat-sheet of the commands you will use this week. All ship with BLAST+ 2.15.0.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `makeblastdb` | Build a database from a FASTA file | `-in`, `-dbtype` (`nucl`/`prot`), `-out`, `-parse_seqids`, `-taxid`, `-taxid_map` |
| `blastn` | Nucleotide query vs nucleotide DB | `-query`, `-db`, `-out`, `-outfmt`, `-evalue`, `-word_size`, `-task` (`megablast`/`blastn`/`dc-megablast`/`blastn-short`) |
| `blastp` | Protein query vs protein DB | `-query`, `-db`, `-out`, `-outfmt`, `-evalue`, `-matrix` (default BLOSUM62), `-word_size` (default 3) |
| `tblastn` | Protein query vs translated nucleotide DB | as `blastp` plus `-db_gencode` (default 1 = standard genetic code) |
| `blastx` | Translated nucleotide query vs protein DB | as `blastp` plus `-query_gencode` |
| `blastdbcmd` | Extract sequences/metadata from a BLAST DB | `-db`, `-entry`, `-range`, `-outfmt`, `-info` |
| `update_blastdb.pl` | Download/refresh prebuilt NCBI databases | `--decompress`, `--passive`, database name |

### BLAST tabular columns (`-outfmt 6`, default 12)

| # | Name | Meaning |
|---|------|---------|
| 1 | `qseqid` | Query sequence ID |
| 2 | `sseqid` | Subject (database hit) sequence ID |
| 3 | `pident` | Percent identity over the aligned region |
| 4 | `length` | Alignment length (including gaps) |
| 5 | `mismatch` | Number of mismatches |
| 6 | `gapopen` | Number of gap openings |
| 7 | `qstart` | Start of alignment in query (1-based) |
| 8 | `qend` | End of alignment in query (1-based, inclusive) |
| 9 | `sstart` | Start of alignment in subject |
| 10 | `send` | End of alignment in subject |
| 11 | `evalue` | E-value |
| 12 | `bitscore` | Bit score |

Add `staxid` (subject taxon ID) and `sscinames` (subject scientific name) with a custom `-outfmt "6 std staxid sscinames"`.

## NCBI Entrez E-utilities (programmatic access)

- **`esearch`** — search a database (e.g. `nuccore`, `protein`, `taxonomy`, `pubmed`) by a query string. Returns a `WebEnv` token plus a list of UIDs.
- **`efetch`** — fetch the records for a list of UIDs in a chosen format (`fasta`, `gb`, `xml`).
- **`esummary`** — get compact metadata for a list of UIDs (sequence length, organism, definition).
- **`elink`** — navigate cross-references (e.g. "for these nucleotide accessions, what are the linked taxonomy IDs?").

All four are available in Biopython 1.83 as `Bio.Entrez.esearch`, `Bio.Entrez.efetch`, `Bio.Entrez.esummary`, `Bio.Entrez.elink`. Set `Bio.Entrez.email = "you@example.com"` once at program start.

### Rate limits

- Without an API key: **3 requests/second** maximum.
- With an [NCBI API key](https://www.ncbi.nlm.nih.gov/account/settings/) (free, 30 seconds to register): **10 requests/second**.
- Bulk submissions outside US business hours are politely encouraged.
- 429 responses include `Retry-After` headers. Biopython's `Bio.Entrez` honors them automatically with exponential backoff.

## NCBI Taxonomy

- **Web interface:** <https://www.ncbi.nlm.nih.gov/taxonomy>
- **Taxonomy dump (the flat-file release):** <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/> — `taxdump.tar.gz` is the canonical hierarchical taxonomy as `nodes.dmp` (parent-child edges) and `names.dmp` (display names). ~80 MB.
- **Programmatic lookup with Biopython:**

  ```python
  from Bio import Entrez
  Entrez.email = "you@example.com"
  handle = Entrez.efetch(db="taxonomy", id="9606", retmode="xml")
  records = Entrez.read(handle)
  print(records[0]["ScientificName"])   # "Homo sapiens"
  print(records[0]["Lineage"])
  # "cellular organisms; Eukaryota; ...; Mammalia; ...; Hominidae; Homo"
  ```

- **Reverse lookup** (accession → taxon ID): use `esummary` against `nuccore` or `protein`, read the `TaxId` field.

## Tools you will install this week

- **BLAST+ 2.15.0** — `conda install -c bioconda blast=2.15` (or download the installer from <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/>). Adds `makeblastdb`, `blastn`, `blastp`, `tblastn`, `blastx`, `tblastx`, `blastdbcmd`, `update_blastdb.pl` to your PATH.
- **Biopython 1.83** — `pip install biopython==1.83`. The `Bio.Blast`, `Bio.Entrez`, and `Bio.SeqIO` modules live here.
- **pandas** — `pip install pandas`. For parsing BLAST tabular output and joining with taxonomy.
- **(optional)** **DIAMOND 2.1+** — `conda install -c bioconda diamond`. ~100x faster than BLAST+ for protein search. Used in the Stretch section of the mini-project.

## Free books and chapters

- **Biopython Tutorial and Cookbook** — Chapters 7 (BLAST), 9 (Entrez):
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **Bioinformatics Algorithms: An Active Learning Approach (Compeau & Pevzner)** — Chapter 5 is the canonical undergrad treatment of BLAST seed-and-extend with diagrams; free chapters at:
  <https://www.bioinformaticsalgorithms.org/>
- **Durbin, Eddy, Krogh, Mitchison — Biological Sequence Analysis (1998)** — Chapter 4 is the rigorous statistical treatment of database search; not free in print but most university libraries have it.
- **Mount, *Bioinformatics: Sequence and Genome Analysis* (2nd ed., 2004)** — older but still the most thorough treatment of BLAST mechanics. Public copies common in academic libraries.

## Worked examples online (free)

- **NCBI BLAST web tutorial** — a one-hour walkthrough of the web interface with real example queries:
  <https://www.ncbi.nlm.nih.gov/books/NBK1734/>
- **NCBI E-utilities Quick Start** — a one-page guide to `esearch` + `efetch`:
  <https://www.ncbi.nlm.nih.gov/books/NBK25500/>
- **Rosalind — Bioinformatics Stronghold, BLAST problems:**
  <https://rosalind.info/problems/locations/>
- **EBI BLAST tutorial** — the EBI web interface is a useful fallback when NCBI is busy:
  <https://www.ebi.ac.uk/Tools/sss/ncbiblast/>

## Open-source code to read this week

You can learn more from one hour reading the BLAST+ source than from three hours of tutorials. The BLAST+ codebase is a large C++ project; the parts worth reading are small:

- **`ncbi_blast_seqsrc.cpp`** in the BLAST+ source — the database-iterator interface. Shows how BLAST loads sequences out of a packed database file (`.nhr`/`.nin`/`.nsq` for nucleotide, `.phr`/`.pin`/`.psq` for protein) without copying them all into RAM:
  <https://github.com/ncbi/ncbi-cxx-toolkit-public>
- **Biopython `Bio/Blast/NCBIXML.py`** — the BLAST XML parser, ~600 lines of pure Python:
  <https://github.com/biopython/biopython/blob/master/Bio/Blast/NCBIXML.py>
- **Biopython `Bio/Entrez/__init__.py`** — the Entrez wrapper, including rate-limit and retry handling:
  <https://github.com/biopython/biopython/blob/master/Bio/Entrez/__init__.py>
- **DIAMOND** (Buchfink et al.) — modern, SIMD-accelerated protein aligner, ~100x faster than BLAST+ for `blastp`/`blastx`. Read `src/data/sorted_list.h` and `src/align/align.cpp` for the seed-and-extend core:
  <https://github.com/bbuchfink/diamond>

## E-value cheat sheet

A rough rule of thumb for interpreting an E-value, given a reasonably large database:

| E-value | Interpretation | Action |
|--------:|----------------|--------|
| `< 1e-50` | Almost certainly homologous | Trust the hit; check coverage and identity for the biology |
| `1e-10` to `1e-50` | Strong homology | Trust the hit; verify nothing weird about the subject |
| `1e-3` to `1e-10` | Probable homology | Look at coverage, identity, and the alignment by eye |
| `1e-3` to `1` | Possible / weak / could be chance | Need supporting evidence (synteny, structure, additional hits) |
| `> 1` | Probably chance | Do not report as a hit without strong external evidence |

These cutoffs are *not* universal. For a 16S classifier against a curated database of ~25,000 sequences, `E < 1e-50` is conservative; for a protein query against the 600M-entry `nr`, `E < 1e-3` is the published practical floor. The cutoff is part of your method, document it.

## Bit-score cheat sheet

Bit scores are normalized: a bit score of `S` means the hit is `2^S` times less likely than chance to arise in a random search. A few reference points:

| Bit score | Interpretation |
|----------:|----------------|
| `< 40` | Essentially noise; ignore. |
| `40–50` | Borderline; useful only with corroborating evidence. |
| `50–80` | Probable homology in small-to-medium databases. |
| `> 80` | Almost certainly homologous regardless of database size. |
| `> 200` | "Obvious" homolog; the kind of hit that appears in textbook screenshots. |

The advantage of the bit score over the E-value is that it does *not* depend on database size. Comparing a hit's biological strength across two searches against different databases is a job for bit scores; ranking hits *within* a single search is a job for E-values.

## BLAST family — which one when

| Question | Tool |
|----------|------|
| "I have a DNA sequence; what known DNA matches it?" | `blastn` |
| "I have a protein sequence; what known proteins match it?" | `blastp` |
| "I have a protein sequence; is it encoded somewhere in this genome?" | `tblastn` |
| "I have an unannotated DNA sequence; does it encode a known protein?" | `blastx` |
| "I have two distantly related DNAs and want to find conserved coding regions?" | `tblastx` (slow; rarely the answer) |
| "I have a short primer or ~20 bp sequence?" | `blastn -task blastn-short` (lowers `W` and adjusts heuristics) |

---

*If a link 404s, please open an issue so we can replace it.*
