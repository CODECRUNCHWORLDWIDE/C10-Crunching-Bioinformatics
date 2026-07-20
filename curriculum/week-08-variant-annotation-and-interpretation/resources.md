# Week 8 — Resources

> **Educational and research use only.** None of the resources, tools, or databases listed below is intended as a substitute for clinical-grade variant interpretation. Annotations produced by these tools must not be used for diagnosis, prognosis, or treatment selection. The same point lives at the top of every Week 8 file, by design.

Every resource on this page is **free** and **publicly accessible**. Where we name a version (VEP 110.1, SnpEff 5.2, ClinVar VCF release 2024-09-01, gnomAD v4.1.0, cyvcf2 0.30.28, pysam 0.22.1), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **McLaren, Gil, Hunt, Riat, Ritchie, Thormann, Flicek, Cunningham (2016)** — "The Ensembl Variant Effect Predictor." The VEP paper. *Genome Biology* 17:122. Free full text:
  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4>
  Tool documentation:
  <https://www.ensembl.org/info/docs/tools/vep/index.html>
- **Cingolani, Platts, Wang, Coon, Nguyen, Wang, Land, Lu, Ruden (2012)** — "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff." The SnpEff paper. *Fly* 6:80. Free full text:
  <https://www.tandfonline.com/doi/full/10.4161/fly.19695>
  Tool documentation:
  <https://pcingola.github.io/SnpEff/>
- **Richards, Aziz, Bale, Bick, Das, Gastier-Foster, Grody, Hegde, Lyon, Spector, Voelkerding, Rehm (2015)** — "Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology." The ACMG/AMP 2015 paper. *Genetics in Medicine* 17:405. Free full text at PMC:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/>
  This is the most important paper of the week. Read it end to end; it is ~12 pages and forms the foundation of clinical variant interpretation worldwide.
- **Karczewski, Francioli, Tiao, Cummings, Alfoldi, Wang, Collins, Laricchia, Ganna, Birnbaum, Gauthier, Brand, Solomonson, Watts, Rhodes, Singer-Berk, England, Seaby, Kosmicki, Walters, Tashman, Farjoun, Banks, Poterba, Wang, Seed, Whiffin, Chong, Samocha, Pierce-Hoffman, Zappala, O'Donnell-Luria, Minikel, Weisburd, Lek, Ware, Vittal, Armean, Bergelson, Cibulskis, Connolly, Covarrubias, Donnelly, Ferriera, Gabriel, Gentry, Gupta, Jeandet, Kaplan, Llanwarne, Munshi, Novod, Petrillo, Roazen, Ruano-Rubio, Saltzman, Schleicher, Soto, Tibbetts, Tolonen, Wade, Talkowski, MacArthur, Daly, Neale (2020)** — "The mutational constraint spectrum quantified from variation in 141,456 humans." The gnomAD v2 flagship paper. *Nature* 581:434. Free full text:
  <https://www.nature.com/articles/s41586-020-2308-7>
  Database landing page:
  <https://gnomad.broadinstitute.org/>
- **Landrum, Lee, Benson, Brown, Chao, Chitipiralla, Gu, Hart, Hoffman, Jang, Karapetyan, Katz, Liu, Maddipatla, Malheiro, McDaniel, Ovetsky, Riley, Zhou, Holmes, Kattman, Maglott (2018)** — "ClinVar: improving access to variant interpretations and supporting evidence." The ClinVar paper. *Nucleic Acids Research* 46:D1062. Free full text:
  <https://academic.oup.com/nar/article/46/D1/D1062/4641904>
  Database landing page:
  <https://www.ncbi.nlm.nih.gov/clinvar/>
- **Sherry, Ward, Kholodov, Baker, Phan, Smigielski, Sirotkin (2001)** — "dbSNP: the NCBI database of genetic variation." The dbSNP paper. *Nucleic Acids Research* 29:308. Free full text:
  <https://academic.oup.com/nar/article/29/1/308/1116004>
  Database landing page:
  <https://www.ncbi.nlm.nih.gov/snp/>
- **Ng, Henikoff (2003)** — "SIFT: predicting amino acid changes that affect protein function." The SIFT paper. *Nucleic Acids Research* 31:3812. Free full text:
  <https://academic.oup.com/nar/article/31/13/3812/2904440>
- **Adzhubei, Schmidt, Peshkin, Ramensky, Gerasimova, Bork, Kondrashov, Sunyaev (2010)** — "A method and server for predicting damaging missense mutations." The PolyPhen-2 paper. *Nature Methods* 7:248. Free preprint at PMC:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2855889/>
- **Whirl-Carrillo, Huddart, Gong, Sangkuhl, Thorn, Whaley, Klein (2021)** — "An evidence-based framework for evaluating pharmacogenomics knowledge for personalized medicine." The PharmGKB paper. *Clinical Pharmacology and Therapeutics* 110:563. Free full text:
  <https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2350>
  Database landing page:
  <https://www.pharmgkb.org/>
- **Pedersen, Quinlan (2017)** — "cyvcf2: fast, flexible variant analysis with Python." The cyvcf2 paper. *Bioinformatics* 33:1867. Free full text:
  <https://academic.oup.com/bioinformatics/article/33/12/1867/2978460>
  Tool documentation:
  <https://brentp.github.io/cyvcf2/>

## Tool reference (the command-line surface)

### VEP 110.1 (Ensembl Variant Effect Predictor)

VEP is the Ensembl-maintained variant annotator. It is a Perl CLI; the conda package installs the binary plus a wrapper that fetches the cache. The cache is a static download (~25 GB for human GRCh38) and the annotation runs offline once the cache is in place.

| Flag | Purpose |
|------|---------|
| `-i input.vcf` | Input VCF |
| `-o output.vcf` | Output VCF (annotated, with CSQ INFO field) |
| `--vcf` | Output as VCF (default; alternatives: `--json`, `--tab`) |
| `--cache --dir_cache /path/to/cache` | Use a local cache (no network at runtime) |
| `--offline` | Force offline mode |
| `--species homo_sapiens --assembly GRCh38` | Required species and assembly |
| `--cache_version 110` | Pin the cache release |
| `--sift b --polyphen b` | Include SIFT and PolyPhen-2 in the output (`b`oth score and prediction) |
| `--canonical` | Flag the canonical transcript |
| `--mane` | Flag the MANE Select transcript |
| `--symbol` | Include the gene symbol |
| `--biotype` | Include the transcript biotype |
| `--af --af_gnomadg --af_gnomade` | Include 1000 Genomes, gnomAD genome, gnomAD exome allele frequencies |
| `--check_existing` | Include the rsID (Existing_variation) |
| `--fork 4` | 4 worker processes |
| `--force_overwrite` | Allow overwriting the output file |
| `--no_stats` | Skip the HTML stats page if you do not need it |

#### The canonical VEP call

```bash
vep \
    -i input.vcf \
    -o output.vep.vcf \
    --vcf \
    --cache \
    --dir_cache ${HOME}/.vep \
    --offline \
    --species homo_sapiens \
    --assembly GRCh38 \
    --cache_version 110 \
    --sift b --polyphen b \
    --canonical --mane --symbol --biotype \
    --af --af_gnomadg --af_gnomade \
    --check_existing \
    --fork 4 \
    --force_overwrite \
    --no_stats
```

For a 50-variant demo VCF, this runs in ~5 seconds on a laptop (the slow step is reading the cache; once cached in OS file cache, subsequent runs are < 1 second).

### SnpEff 5.2

SnpEff is a Java-based variant annotator. The conda package installs a wrapper script `snpEff` plus the jar.

| Command | Purpose |
|---------|---------|
| `snpEff download GRCh38.105` | Download a SnpEff database for a genome build |
| `snpEff databases | grep GRCh38` | List available human databases |
| `snpEff -v GRCh38.105 input.vcf > output.snpeff.vcf` | Annotate (default verbose; emits the ANN INFO field) |
| `snpEff -v -csvStats stats.csv GRCh38.105 input.vcf > output.vcf` | Annotate and emit per-effect summary stats |

#### The canonical SnpEff call

```bash
snpEff \
    -v \
    -csvStats stats.snpeff.csv \
    -htmlStats stats.snpeff.html \
    GRCh38.105 \
    input.vcf \
    > output.snpeff.vcf
```

Output VCF gains an `ANN` INFO field per variant, pipe-separated, with one record per overlapping transcript. The SnpEff impact column (`HIGH`, `MODERATE`, `LOW`, `MODIFIER`) maps approximately onto VEP's IMPACT but is computed independently; expect ~95% agreement and ~5% disagreement on the impact category, almost always on edge cases (splice region vs intronic, missense vs synonymous on a non-canonical transcript).

### cyvcf2 0.30.28

`cyvcf2` is a Cython binding to htslib's vcf.c. It is the fastest VCF parser in the Python ecosystem (~10x faster than `pysam.VariantFile`, ~100x faster than `PyVCF`). It exposes a `VCF` reader, a `Variant` object with the standard CHROM/POS/REF/ALT/INFO/FORMAT fields, and convenient access to per-sample genotypes as numpy arrays.

```python
from cyvcf2 import VCF, Variant

vcf = VCF("input.vep.vcf")
for variant in vcf:
    chrom: str = variant.CHROM
    pos: int = variant.POS
    ref: str = variant.REF
    alts: list[str] = variant.ALT  # list (ALT is a list because of multi-allelics)
    csq_field: str = variant.INFO.get("CSQ", "")
    # CSQ is pipe-separated; multiple records per variant separated by commas.
```

### pysam 0.22.1

`pysam` provides BAM, VCF, and FASTA access. For Week 8 we use it mostly for FASTA access (extracting the reference codon around a variant) and as a fallback for VCF when cyvcf2 is not available. The Week 5 patterns for `pysam.AlignmentFile` and `pysam.FastaFile` carry over.

### ClinVar VCF release

ClinVar publishes a VCF release every two weeks. The current release is at:

```bash
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

Each record has at minimum: CHROM, POS, REF, ALT, RS (dbSNP rsID), and the following INFO fields: `CLNSIG` (clinical significance; Benign / Likely_benign / Uncertain_significance / Likely_pathogenic / Pathogenic / Conflicting_interpretations_of_pathogenicity), `CLNREVSTAT` (review status; `practice_guideline` is 4-star, `reviewed_by_expert_panel` is 3-star, `criteria_provided,_multiple_submitters,_no_conflicts` is 2-star, `criteria_provided,_single_submitter` is 1-star, `no_assertion_criteria_provided` is 0-star), `CLNDN` (disease names), `CLNDISDB` (disease database IDs), `MC` (molecular consequence).

To query for a specific coordinate after downloading:

```bash
tabix clinvar.vcf.gz chr7:117559590-117559590
```

For programmatic JSON access without downloading:

```python
import requests
r = requests.get(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    params={"db": "clinvar", "id": "<variation_id>", "retmode": "json"},
    timeout=30,
)
data = r.json()
```

The E-utilities are free but rate-limited to 3 requests/second without an API key. Add `&api_key=<your_key>` to bump to 10 requests/second.

### gnomAD GraphQL API

gnomAD's API is a GraphQL endpoint at:

```
https://gnomad.broadinstitute.org/api
```

A minimal query for a single variant:

```python
import requests

query = """
{
  variant(variantId: "7-117559590-G-A", dataset: gnomad_r4) {
    variant_id
    rsid
    exome {
      ac
      an
      af
      populations {
        id
        ac
        an
      }
    }
    genome {
      ac
      an
      af
    }
  }
}
"""
r = requests.post(
    "https://gnomad.broadinstitute.org/api",
    json={"query": query},
    timeout=30,
)
data = r.json()["data"]["variant"]
```

The dataset identifiers as of mid-2024:

- `gnomad_r2_1` — gnomAD v2.1.1 exomes (~125,000 samples; GRCh37).
- `gnomad_r3` — gnomAD v3.1.2 genomes (~76,000 samples; GRCh38).
- `gnomad_r4` — gnomAD v4.0/v4.1 combined (807,162 samples; GRCh38). Use this by default.

Sub-populations (`populations[].id`): `afr` (African), `amr` (Admixed American / Latino), `asj` (Ashkenazi Jewish), `eas` (East Asian), `fin` (Finnish), `nfe` (Non-Finnish European), `sas` (South Asian), `oth` (Other), `mid` (Middle Eastern, new in v4).

The gnomAD API is rate-limited to ~10 requests/second per IP. The polite pattern is 200 ms sleep between requests plus a 3-attempt retry on 429 / 503.

### dbSNP via NCBI E-utilities

```python
import requests

r = requests.get(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    params={"db": "snp", "id": "121908755", "retmode": "json"},
    timeout=30,
)
data = r.json()
```

The rsID is the integer after the `rs` prefix (`rs121908755` -> `121908755`).

### PharmGKB API

PharmGKB exposes a free REST API at:

```
https://api.pharmgkb.org/v1
```

A query for a gene's annotated variants:

```python
import requests

r = requests.get(
    "https://api.pharmgkb.org/v1/data/gene",
    params={"symbol": "CYP2D6", "view": "max"},
    timeout=30,
)
data = r.json()
```

CPIC guidelines are bundled with PharmGKB under `/v1/data/guideline`. Each guideline links a gene-drug pair to a recommendation table indexed by metabolizer phenotype (ultrarapid / normal / intermediate / poor).

### ACMG criteria as code

There is no canonical free CLI implementation of the ACMG criteria; what exists is the Richards et al. 2015 decision table, several open implementations (`InterVar`, `CharGer`, `PathoMan`), and a body of literature on edge cases. The Week 8 exercises implement a minimal subset by hand. Real clinical pipelines use one of the tools above plus extensive manual curation.

| Criterion | What it means | Mechanically computable? |
|-----------|---------------|--------------------------|
| **PVS1** | Null variant (nonsense, frameshift, canonical +/-1,2 splice site) in a gene where LOF is a known disease mechanism | Yes, if you have a curated LOF-mechanism gene list |
| **PS1** | Same amino acid change as a known pathogenic variant | Yes, with a ClinVar Pathogenic table |
| **PS2** | De novo (paternity/maternity confirmed) | No (requires trio data + parentage) |
| **PS3** | Well-established functional studies | No (literature review) |
| **PS4** | Increased prevalence in affected vs controls | No (case-control study) |
| **PM1** | In a mutational hotspot / functional domain | Partial (requires curated domains) |
| **PM2** | Absent or very rare in population databases | **Yes** (gnomAD AF < 0.0001) |
| **PM3** | In trans with pathogenic for recessive | No (phasing + ClinVar) |
| **PM4** | Protein length change | Yes (VEP consequence) |
| **PM5** | Novel missense at a residue with a different pathogenic missense | Yes, with a ClinVar table |
| **PM6** | Assumed de novo (no confirmation) | No |
| **PP1** | Cosegregation with disease | No (pedigree) |
| **PP2** | Missense in a gene with low rate of benign missense | Partial (gene constraint) |
| **PP3** | Multiple lines of computational evidence | **Yes** (SIFT + PolyPhen + CADD) |
| **PP4** | Phenotype highly specific for a single-gene disorder | No (clinical) |
| **PP5** | Reputable source reports as pathogenic | Yes (ClinVar) |
| **BA1** | AF > 5% in any general population | **Yes** (gnomAD popmax > 0.05) |
| **BS1** | AF greater than expected for the disorder | Yes, if you know the disorder prevalence |
| **BS2** | Observed homozygous in a healthy adult for a recessive | Partial (gnomAD) |
| **BS3** | Well-established functional study showing no effect | No |
| **BS4** | Lack of segregation in affected family | No |
| **BP1** | Missense where only truncating cause disease | Partial (gene-specific) |
| **BP2** | In trans with pathogenic for dominant / in cis with pathogenic | No |
| **BP3** | In-frame indel in a repetitive region | Yes |
| **BP4** | Multiple lines of computational evidence of no impact | **Yes** (SIFT + PolyPhen + CADD) |
| **BP5** | Variant in a case with an alternate cause | No |
| **BP6** | Reputable source reports as benign | Yes (ClinVar) |
| **BP7** | Synonymous with no splicing impact | Yes (VEP consequence + SpliceAI) |

Of the 28 criteria, **8 are cleanly mechanically computable** (bold above): PVS1 (with a gene list), PM2, PP3, PP5, BA1, BS1 (with a disorder prevalence), BP4, BP6, plus PM4, PM5, PS1, BP3, BP7 with the right data. The remaining ~12 require manual review and cannot be mechanized in a Week-8-scoped pipeline. We will be explicit about which criteria are "computed" and which are "skipped" in the report output.

## Compute requirements

| Step | Demo (50 variants) | Real WGS (~5 M variants) |
|------|--------------------|---------------------------|
| `vep --offline` | ~5 sec | ~30 min |
| `snpEff` | ~5 sec | ~10 min |
| ClinVar tabix lookup | ~50 ms per variant | (use bulk join) |
| gnomAD API query | ~100 ms per variant, rate-limited | (use bulk download) |
| cyvcf2 parse | ~1 ms per variant | ~10 sec |
| Report build (pandas to HTML) | ~1 sec | ~30 sec on filtered subset |

For real-WGS-scale pipelines, do not query the gnomAD API per variant — download the gnomAD VCF and join on coordinate. For Week 8's demo VCFs (50-100 variants), the API is fine and the network call is the slow step.

A laptop with 16 GB RAM is enough for everything. The VEP cache (~25 GB) is the biggest disk requirement; if you do not have the space, use the Codespaces or Colab path.

## Datasets

All datasets used in Week 8 are no-cost and publicly distributed.

### The demo VCF (used by Exercise 1, 3 and the mini-project)

A 50-variant subset of a curated mixture: ~30 variants from the **GIAB NA12878 truth set** v4.2.1 (Zook et al. 2014, *Nature Biotechnology* 32:246) and ~20 hand-picked variants from **ClinVar** with known Pathogenic or Likely Pathogenic classifications across well-characterized disease genes (`BRCA1`, `BRCA2`, `CFTR`, `MLH1`, `MSH2`, `APC`, `TP53`, `FBN1`, `LDLR`, `HEXA`). The hand-picked variants are *for educational examples only* — they are real published pathogenic variants but they are not from a patient sample and the inclusion of a variant in this teaching set does not constitute clinical information about anyone.

The file is bundled with the curriculum at `data/demo.vcf.gz` (~10 KB).

### The Variant Effect Predictor cache

```bash
# Use the official VEP cache mirror. ~25 GB for human GRCh38.
vep_install \
    --AUTO cf \
    --SPECIES homo_sapiens \
    --ASSEMBLY GRCh38 \
    --CACHEDIR ${HOME}/.vep \
    --NO_HTSLIB
```

Or download manually from:

```
https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
```

(~14 GB compressed.)

### The SnpEff database

```bash
snpEff download -v GRCh38.105
```

(~2.5 GB to `~/snpEff/data/GRCh38.105/`.)

### ClinVar VCF release

```bash
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

(~200 MB.) Released every two weeks. The Week 8 mini-project pins the release date in the run-info JSON.

### gnomAD downloads (optional; the API is fine for small VCFs)

```bash
# v4.1 exomes, chr-by-chr. Each chromosome is ~1-5 GB.
curl -sLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz
```

Only download what you need; the full v4.1 exome is ~80 GB.

### PharmGKB downloads

```bash
# Star-allele definitions and CPIC guideline tables.
curl -sLO https://api.pharmgkb.org/v1/download/file/data/allele.zip
curl -sLO https://api.pharmgkb.org/v1/download/file/data/guideline.zip
```

Combined ~50 MB.

## Free-tier compute and storage

- **GitHub Codespaces** (free 60 hours/month for personal accounts): a 4-core / 16 GB instance with 32 GB persistent storage is enough for the VEP cache plus the mini-project pipeline. The codespace template in the curriculum repo has the cache pre-baked into a devcontainer image.
- **Google Colab** (free tier): also enough for the mini-project, *if you skip the local VEP cache and use Ensembl's REST API instead* (slower per variant but no disk). The Colab notebook in the mini-project directory documents this path.
- **Local laptop**: any 16 GB RAM laptop with 50 GB free is enough.

## Style guide for this week

- **Pin tool versions.** Always write "VEP 110.1", "SnpEff 5.2", "ClinVar VCF release 2024-09-01", "gnomAD v4.1.0", "cyvcf2 0.30.28" — never just "VEP" or "the latest ClinVar." Annotations drift; reproducibility requires pinning.
- **Pin reference assemblies.** Always write "GRCh38" or "GRCh37"; never "the human reference." Mixing assemblies silently shifts coordinates by hundreds of base pairs and breaks every join.
- **Cite tools by paper.** "VEP (McLaren et al. 2016)" or "ClinVar (Landrum et al. 2018)" — not just "VEP" or "ClinVar". This is standard practice in publications and what reviewers expect.
- **Report numbers, not adjectives.** "gnomAD v4 reports this variant at AF 0.0007 in non-Finnish Europeans (allele count 14 of 19,532)" is a sentence. "Rare in healthy populations" is not — at minimum, give the AF and the source.
- **Always print the database version, the run date, and the disclaimer in the report header.** Without these, you have an undated PDF of "results" that nobody can reproduce or trust.
- **Use VEP for the canonical annotation in production; SnpEff for the second opinion.** Most clinical pipelines run both and flag variants where the two disagree. This week's exercises do the same.

## Common questions

**Q. VEP and SnpEff disagree on the consequence type for one of my variants. Which is right?**

Neither is uniquely "right" — they are two implementations of the same Sequence Ontology call against the same gene model. Disagreements almost always trace to (a) different transcript-of-reference choices (one tool used the canonical transcript, the other used MANE Select), (b) edge cases at splice-region boundaries, or (c) different versions of the gene model. Investigate the disagreement on the actual transcript coordinates; do not just pick one and move on.

**Q. ClinVar reports this variant as both Pathogenic and Likely Benign by different submitters. What does that mean?**

This is called a "Conflicting interpretations of pathogenicity" record in ClinVar. ~3% of ClinVar records have it. The right action is *not* to pick the submitter you like best; it is to read all submissions, weigh by review status (4-star expert panel > 1-star single submitter), and flag the variant as "conflicting, requires manual review" in your output. Treating a conflicting record as a single classification is the most common Week 8 error.

**Q. Why does gnomAD report a different allele frequency than 1000 Genomes for the same variant?**

Different cohorts, different ascertainment criteria, different sample sizes. gnomAD v4 is ~807K samples drawn from biobanks and case studies (with severe-pediatric-disease cases excluded); 1000 Genomes phase 3 is 2,504 unrelated samples drawn from defined populations. gnomAD is the default for "is this variant common in healthy people" but you should be aware that gnomAD's sample selection is not population-representative in a strict sampling-theory sense.

**Q. Can I use this pipeline on my own VCF from 23andMe or AncestryDNA?**

You can run it, but be very careful. 23andMe and AncestryDNA produce *array-based genotypes*, not sequencing-based variant calls; the genotype quality is high for the targeted SNPs but the genome coverage is sparse (~600,000 of ~3 billion positions). Most clinical-grade interpretations require sequencing. Treating an array-based call of "Pathogenic per ClinVar" as actionable is exactly the kind of mistake the educational disclaimer at the top of every Week 8 file is meant to prevent.

**Q. What is the difference between "Likely Pathogenic" and "Pathogenic" in the ACMG framework?**

Both are in the "deleterious" tier, but they reflect different levels of evidence. Pathogenic is ≥ 1 PVS plus ≥ 2 PS plus ≥ 2 PM, or the equivalent across categories. Likely Pathogenic is the same evidence with one tier less. In practice, the clinical action on Likely Pathogenic and Pathogenic is often the same — the distinction matters most for the "ongoing surveillance" cohort where evidence may accumulate over time.

**Q. Do I need to commit the VEP cache to git?**

No. The cache is 25 GB and is a static download. Gitignore it. The run-info JSON should record the cache version (e.g. `110`), and a fresh checkout downloads the cache once. Same for the SnpEff database and the gnomAD VCFs.

**Q. The mini-project asks for an HTML report. Should I use a templating engine like Jinja?**

Yes, Jinja2 is the lab-standard choice and the Week 8 mini-project uses it. The alternative is `pandas.DataFrame.to_html()` with a static CSS file; that works too and is what the demo report template uses. Either is acceptable as long as the report renders in a browser without an internet connection.
