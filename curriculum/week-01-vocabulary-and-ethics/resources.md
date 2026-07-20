# Week 1 — Resources

Every resource on this page is **free** and **publicly accessible**. No paywalled textbooks, no proprietary PDFs. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **NCBI Handbook** — the operational manual for the world's largest public sequence database:
  <https://www.ncbi.nlm.nih.gov/books/NBK143764/>
- **Ensembl — Getting Started**:
  <https://www.ensembl.org/info/website/tutorials/index.html>
- **EMBL-EBI Training — Bioinformatics for the terrified** (90 min, free):
  <https://www.ebi.ac.uk/training/online/courses/bioinformatics-terrified/>
- **NHGRI — A Brief Guide to Genomics**:
  <https://www.genome.gov/about-genomics/fact-sheets/A-Brief-Guide-to-Genomics>
- **1000 Genomes Project — Data portal and consent overview**:
  <https://www.internationalgenome.org/>

## File-format references (skim now, return in Weeks 2–6)

- **FASTA / FASTQ formats — NCBI quick reference**:
  <https://blast.ncbi.nlm.nih.gov/doc/blast-topics/>
- **SAM/BAM specification (samtools)**:
  <https://samtools.github.io/hts-specs/SAMv1.pdf>
- **VCF specification (samtools)**:
  <https://samtools.github.io/hts-specs/VCFv4.4.pdf>
- **GFF3 specification (Sequence Ontology)**:
  <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>

## Data ethics — required for Lecture 2

- **NHGRI — Ethical, Legal and Social Implications (ELSI) program**:
  <https://www.genome.gov/about-nhgri/Division-of-Genomics-and-Society/ELSI-Research-Program>
- **Gymrek et al., 2013 — "Identifying personal genomes by surname inference"** (*Science*). The foundational re-identification paper. Open via NCBI:
  <https://pubmed.ncbi.nlm.nih.gov/23329047/>
- **All of Us Research Program — Consent and ethics**:
  <https://allofus.nih.gov/about/protocol>
- **GA4GH (Global Alliance for Genomics and Health) — Framework for Responsible Sharing of Genomic and Health-Related Data**:
  <https://www.ga4gh.org/framework/>
- **HHS (US) — HIPAA Privacy Rule overview** (high level only):
  <https://www.hhs.gov/hipaa/for-professionals/privacy/index.html>
- **European Commission — GDPR and health data** (high level):
  <https://commission.europa.eu/law/law-topic/data-protection_en>

## Free books (chapter-level, not whole books)

- **Bioinformatics Data Skills**, Vince Buffalo — author posts free chapters on his site. Chapters 1, 2, 5 are the right Week-1/2 reads:
  <https://vincebuffalo.com/bds/>
- **Biopython Tutorial and Cookbook** (full PDF, free):
  <https://biopython.org/DIST/docs/tutorial/Tutorial.pdf>
- **Computational Genomics with R**, Akalin (open access):
  <https://compgenomr.github.io/book/>
- **Modern Statistics for Modern Biology**, Holmes and Huber (open access):
  <https://www.huber.embl.de/msmb/>

## Official tool documentation (you will bounce between these all course)

- **Biopython**: <https://biopython.org/>
- **samtools / bcftools**: <https://www.htslib.org/>
- **bwa**: <https://bio-bwa.sourceforge.net/>
- **minimap2**: <https://lh3.github.io/minimap2/>
- **Bioconductor**: <https://bioconductor.org/>
- **Snakemake**: <https://snakemake.readthedocs.io/>

## Public datasets you will meet in C10

- **NCBI GenBank / RefSeq** — annotated reference sequences: <https://www.ncbi.nlm.nih.gov/refseq/>
- **Ensembl reference genomes** — GRCh38 for human, GRCm39 for mouse: <https://www.ensembl.org/>
- **1000 Genomes Project (phase 3)** — population-scale human variation, consent-cleared: <https://www.internationalgenome.org/>
- **GTEx Portal (public summary data)** — tissue-level gene expression: <https://gtexportal.org/>
- **GISAID** — viral sequences (SARS-CoV-2, influenza). Registration required; academic terms: <https://gisaid.org/>
- **ENA (European Nucleotide Archive)**: <https://www.ebi.ac.uk/ena/>
- **UniProt** — protein sequences and annotations: <https://www.uniprot.org/>
- **dbSNP** — known human variants: <https://www.ncbi.nlm.nih.gov/snp/>
- **ClinVar** — clinically-curated variants: <https://www.ncbi.nlm.nih.gov/clinvar/>

## Videos (free, no signup)

- **Khan Academy — Central dogma of molecular biology** (12 min):
  <https://www.khanacademy.org/science/ap-biology/gene-expression-and-regulation/transcription-and-rna-processing/v/dna-transcription>
- **iBiology — short lectures from working scientists** (free, peer-reviewed):
  <https://www.ibiology.org/>
- **NCBI Insights — webinars on NCBI tools**:
  <https://www.youtube.com/@NCBINLM>

## Open-source projects to read this week

You can learn more from one hour reading other people's code than from three hours of tutorials. Pick one:

- **Biopython** — `biopython/biopython` on GitHub. Open `Bio/SeqIO/FastaIO.py`: <https://github.com/biopython/biopython>
- **samtools** — battle-tested C: <https://github.com/samtools/samtools>
- **Nextstrain** — phylogenetic surveillance for public-health viruses, end-to-end open: <https://github.com/nextstrain>

## Glossary cheat sheet

Keep this open in a tab. We expand it in [exercise-01-glossary.md](./exercises/exercise-01-glossary.md).

| Term | Plain English |
|------|---------------|
| **Nucleotide** | One letter of DNA or RNA (A, C, G, T, or U) |
| **Genome** | The entire DNA content of one organism |
| **Chromosome** | One large DNA molecule; humans have 23 pairs |
| **Gene** | A stretch of DNA that codes for a product (often a protein) |
| **Transcript** | The RNA copy of a gene |
| **Exon / Intron** | Coding / non-coding stretches inside a gene |
| **Codon** | A 3-letter "word" of RNA that maps to one amino acid |
| **Reference genome** | A community-curated example sequence (e.g. GRCh38) |
| **Variant** | A position where an individual differs from the reference |
| **FASTA / FASTQ** | Text formats for sequences (plain / with quality scores) |
| **SAM / BAM** | Text / binary format for aligned reads |
| **VCF** | Variant Call Format — list of differences from a reference |
| **GFF / GTF** | Genomic feature format — gene annotations |
| **Consent-cleared** | Donors agreed to public, research-only use |
| **Re-identification** | Linking "anonymized" data back to an individual |

---

*If a link 404s, please open an issue so we can replace it.*
