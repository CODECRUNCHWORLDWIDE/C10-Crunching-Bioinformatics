# Challenge 2 — Annotate PASS Variants with OncoKB and CIViC Evidence

> **Educational and research use only.** OncoKB and CIViC are research-grade clinical-interpretation databases. Their evidence levels reflect the published literature curated by the database teams. The Week 11 pipeline looks up these levels for educational purposes; it does not constitute clinical advice. Real treatment decisions require a validated clinical pipeline reviewed by a molecular pathologist and a multi-disciplinary tumor board.

## What you will do

Take the FilterMutectCalls PASS variants from Exercise 2, restrict to variants in known cancer-related genes (the COSMIC Cancer Gene Census; free at <https://cancer.sanger.ac.uk/census>), translate each variant to a protein-change notation, and look up each in **OncoKB** (Chakravarty et al. 2017, *JCO Precision Oncology* 1:1; free at <https://www.oncokb.org/>) and **CIViC** (Griffith et al. 2017, *Nature Genetics* 49:170; free at <https://civicdb.org/>). Emit a Markdown report that lists, per variant, the gene, the protein change, the OncoKB Mutation Effect and FDA evidence level, and the CIViC evidence statements with their star ratings.

This is the **interpretation** layer of the pipeline. It is not a clinical call; it is a research-grade summary of what the published literature says about each variant.

## Why

A filtered VCF of somatic SNVs is a starting point, not a conclusion. The variant-level interpretation layer answers questions like:

- Is this variant in a known cancer driver gene?
- Has this exact mutation been reported in tumors before, and at what frequency?
- Is there published evidence that this variant predicts response or resistance to a specific therapy?
- Is this variant a known hotspot (e.g. TP53 R175H, KRAS G12D, BRAF V600E)?

Without this layer, a PASS variant is just a position in the genome. With this layer, it becomes a *prioritized* hypothesis about the tumor's biology.

## Prerequisites

- Exercise 2 done; `results/ex02/filtered.vcf.gz` exists.
- A variant annotator that emits gene names and amino-acid changes — Funcotator (in GATK), VEP (Ensembl), or SnpEff. We use Funcotator for consistency with the GATK ecosystem.
- The COSMIC Cancer Gene Census TSV (free download with academic registration at <https://cancer.sanger.ac.uk/census>) at `data/cancer_gene_census.tsv`.
- The CIViC evidence release (free, no registration; download from <https://civicdb.org/releases/main>) at `data/civic_evidence.tsv`.
- An OncoKB public-tier API token (free with email registration at <https://www.oncokb.org/apiAccess>); store in `~/.oncokb_token` and chmod 600.
- Python 3.11+ with `requests` and `pandas` installed.

## Suggested approach

### Step 1 — annotate the filtered VCF with gene and protein changes (45 min)

Funcotator is GATK's annotation tool. The data sources bundle is ~30 GB; for the chr22 didactic dataset we can use a smaller GENCODE subset.

```bash
# Download the Funcotator data sources (one-time):
gatk FuncotatorDataSourceDownloader --somatic --validate-integrity \
  --extract-after-download

# Annotate:
gatk Funcotator \
  -R data/chr22_GRCh38.fasta \
  -V results/ex02/filtered.vcf.gz \
  --ref-version hg38 \
  --data-sources-path /path/to/funcotator/dataSourcesFolder \
  --output-file-format MAF \
  -O results/ch02/filtered_annotated.maf
```

The MAF (Mutation Annotation Format) output is TSV with one row per variant and columns including `Hugo_Symbol`, `HGVSp_Short` (the protein change in `p.R175H` format), `Variant_Classification`, `Variant_Type`, `Tumor_Sample_Barcode`.

If Funcotator is not available you can fall back to a smaller annotation tool. The simplest is VEP via its public REST API (no installation required) — see the alternative path below.

### Step 2 — restrict to COSMIC Cancer Gene Census genes (15 min)

The COSMIC Cancer Gene Census is a curated list of ~700 genes with established roles in cancer. Variants in CGC genes are *a priori* more likely to be drivers; variants outside CGC may still be drivers but the prior is much smaller.

```python
import pandas as pd
from pathlib import Path

def load_cgc(path: Path) -> set[str]:
    """Load the COSMIC Cancer Gene Census TSV and return the gene-symbol set."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    return set(df["Gene Symbol"].dropna().astype(str).str.upper())


def filter_maf_to_cgc(maf_path: Path, cgc: set[str]) -> pd.DataFrame:
    """Read a Funcotator MAF and filter to rows whose Hugo_Symbol is in CGC."""
    df = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False)
    return df[df["Hugo_Symbol"].astype(str).str.upper().isin(cgc)].copy()
```

On a typical tumor exome you expect 5-20 PASS variants in CGC genes; on the chr22 didactic subset you may have 0-3. The Markdown report will be short but the pattern is exactly the same on a full-exome run.

### Step 3 — look up each variant in OncoKB (30 min)

OncoKB's public-tier annotation API accepts a HUGO symbol and an HGVS protein change and returns the Mutation Effect, the highest FDA evidence level, the linked drugs, and the linked clinical trials.

```python
import requests
from pathlib import Path
import time

def load_oncokb_token() -> str:
    """Read the OncoKB token from ~/.oncokb_token."""
    token_path: Path = Path.home() / ".oncokb_token"
    if not token_path.exists():
        raise FileNotFoundError(
            "OncoKB token not found at ~/.oncokb_token. "
            "Register at https://www.oncokb.org/apiAccess to obtain one."
        )
    return token_path.read_text().strip()


def oncokb_lookup(
    hugo: str,
    protein_change: str,
    token: str,
    tumor_type: str = "BRCA",
    timeout: float = 30.0,
) -> dict:
    """Look up a variant in OncoKB. Returns the JSON response or {} on failure.

    protein_change is the HGVS short form without the 'p.' prefix, e.g. 'R175H'.
    tumor_type is the OncoKB OncoTree code (BRCA = Breast Cancer, LUAD = Lung
    Adenocarcinoma, ...).
    """
    url = "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
    params = {
        "hugoSymbol": hugo,
        "alteration": protein_change,
        "tumorType": tumor_type,
    }
    headers = {
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
    }
    try:
        resp = requests.get(url, params=params, headers=headers, timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        print(f"[oncokb] lookup failed for {hugo} {protein_change}: {exc}")
        return {}
```

For a batch of variants, rate-limit to ~1 request per second:

```python
def batch_oncokb_lookup(
    variants: list[tuple[str, str]],
    token: str,
    tumor_type: str = "BRCA",
) -> dict[tuple[str, str], dict]:
    """Look up each (hugo, protein_change) in OncoKB. Returns a dict keyed by pair."""
    out: dict[tuple[str, str], dict] = {}
    for hugo, pc in variants:
        out[(hugo, pc)] = oncokb_lookup(hugo, pc, token, tumor_type)
        time.sleep(1.0)
    return out
```

The response includes:

- `mutationEffect.knownEffect`: one of "Loss-of-function", "Gain-of-function", "Switch-of-function", "Likely Loss-of-function", "Likely Gain-of-function", "Likely Switch-of-function", "Inconclusive", "Likely Neutral", "Neutral", "Unknown".
- `oncogenic`: one of "Oncogenic", "Likely Oncogenic", "Likely Neutral", "Inconclusive", "Unknown".
- `highestSensitiveLevel`: the FDA evidence level (e.g. "LEVEL_1", "LEVEL_3B", "LEVEL_4").
- `treatments`: a list of linked therapies.

### Step 4 — look up each variant in CIViC (30 min)

CIViC is fully free and open. The simplest pattern is to download the TSV release and join in pandas:

```python
import pandas as pd
from pathlib import Path


def load_civic_tsv(path: Path) -> pd.DataFrame:
    """Load the CIViC evidence TSV release. Returns a DataFrame."""
    return pd.read_csv(path, sep="\t", low_memory=False)


def civic_lookup(
    civic: pd.DataFrame,
    hugo: str,
    protein_change: str,
) -> pd.DataFrame:
    """Return CIViC evidence rows matching (hugo, protein_change).

    The CIViC TSV's 'gene' column is the symbol; the 'variant' column has
    forms like 'R175H' or 'p.R175H' depending on the curator. We match on
    case-insensitive equality of gene and on protein-change substring.
    """
    mask = (
        (civic["gene"].astype(str).str.upper() == hugo.upper())
        & civic["variant"].astype(str).str.contains(protein_change, case=False, na=False)
    )
    return civic[mask].copy()
```

The CIViC release TSV has rows for each evidence item, with columns including `gene`, `variant`, `evidence_type` (Predictive, Prognostic, Diagnostic, Predisposing, Functional, Oncogenic), `evidence_level` (1-5 star rating), `clinical_significance` (Sensitivity/Response, Resistance, Better Outcome, Poor Outcome, etc.), `drugs`, `disease`, `citation`.

### Step 5 — render the Markdown report (30 min)

Combine the per-variant OncoKB and CIViC data into a Markdown table.

```markdown
# Challenge 2 — Clinical-Interpretation Lookup

## Inputs

- Filtered VCF: `results/ex02/filtered.vcf.gz`
- Annotation: Funcotator (GATK 4.5.0.0; data sources v1.7.20200521s)
- COSMIC Cancer Gene Census: 2024-Q1 release (~700 genes)
- OncoKB: public-tier annotation API; queried 2026-05-14
- CIViC: 2025-03 release

## PASS variants in CGC genes

| Gene  | Position       | REF>ALT | Protein   | OncoKB Effect       | OncoKB Level | CIViC evidence |
|-------|----------------|---------|-----------|---------------------|--------------|----------------|
| TP53  | chr22:23456789 | C>T     | p.R175H   | Loss-of-function    | LEVEL_3B     | 3 items (4 star) |
| BRCA2 | chr22:34567890 | G>A     | p.E2691K  | Likely Oncogenic    | LEVEL_3A     | 2 items (3 star) |
| ...

## Per-variant detail

### TP53 p.R175H (chr22:23456789 C>T)

- OncoKB Mutation Effect: Loss-of-function ("Inferred from this mutation
  disrupting the TP53 DNA-binding domain")
- OncoKB Highest Level: LEVEL_3B (predictive evidence in another tumor type)
- OncoKB Linked Therapies: APR-246 (in clinical trials)
- CIViC Evidence Item 1: Predictive, 4 stars, Sensitivity/Response to MK-1775
  in TP53-mutant ovarian cancer (Pubmed 28652570).
- CIViC Evidence Item 2: Prognostic, 3 stars, Poor outcome in TP53-mutant
  AML (Pubmed 27191167).
- CIViC Evidence Item 3: Oncogenic, 5 stars, well-characterized
  loss-of-function mutation in the DNA-binding domain.
- COSMIC frequency: 3,247 samples across pan-cancer; the most-mutated
  position in TP53; hotspot status confirmed.

### BRCA2 p.E2691K (chr22:34567890 G>A)

- OncoKB Mutation Effect: ...
- ...
```

## Acceptance criteria

- [ ] `results/ch02/filtered_annotated.maf` exists (from Funcotator or alternative annotator).
- [ ] `results/ch02/cgc_filtered.tsv` exists with the CGC-restricted subset.
- [ ] `results/ch02/oncokb_results.json` exists with the per-variant OncoKB JSON responses.
- [ ] `results/ch02/civic_results.tsv` exists with the per-variant CIViC evidence rows.
- [ ] `results/ch02/interpretation_report.md` exists with the Markdown table and per-variant detail for at least three variants (or all CGC variants if fewer than three).
- [ ] `results/ch02/run-info.json` records the Funcotator data sources version, the COSMIC CGC release date, the OncoKB API date, and the CIViC release date.

## Optional extensions

- **Tumor-type-specific OncoKB query.** Pass `tumorType=LUAD` (Lung Adenocarcinoma) or another OncoTree code. The evidence level changes with tumor type: a variant might be LEVEL_1 in breast cancer but LEVEL_3B in lung. Report the level for the most-likely tumor type given the patient context.
- **Cross-check with the COSMIC variant page.** Look up each variant in the COSMIC website and report the cancer-type frequency. Useful for distinguishing recurrent hotspot mutations from one-off observations.
- **Generate a clinical-style summary.** Add a "Top three actionable findings" section that ranks variants by OncoKB level, gene status (oncogene vs tumor suppressor), and CIViC evidence stars. This mirrors the structure of a clinical tumor board report (but is still research-grade).
- **Programmatic literature link.** For each variant, fetch the abstract of the top Pubmed citation from the CIViC evidence and include the first sentence in the report. The Entrez E-utilities API is free and unauthenticated for up to 3 requests per second.

## What this challenge does NOT do

- It does not return a treatment recommendation. The OncoKB level is a *research summary* of the published literature; turning that into a treatment decision requires (a) confirmation of the variant in an accredited clinical lab, (b) clinical correlation by a molecular pathologist, and (c) multi-disciplinary tumor board review.
- It does not validate the variant. A PASS variant in a research pipeline may still be a false positive. Clinical pipelines re-validate via orthogonal sequencing and / or PCR.
- It does not catch all relevant variants. The pipeline runs Mutect2 + FilterMutectCalls; it does not call CNVs or structural variants. A real clinical interpretation considers all variant classes.

## References

- Chakravarty D, Gao J, Phillips SM, et al. *OncoKB: A Precision Oncology Knowledge Base.* **JCO Precision Oncology** 1:1-16 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5586540/>.
- Griffith M, Spies NC, Krysiak K, et al. *CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer.* **Nature Genetics** 49:170-174 (2017). Free at PMC: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5367263/>.
- Sondka Z, Dhir NB, Carvalho-Silva D, et al. *COSMIC: a curated database of somatic variants and clinical data for cancer.* **Nucleic Acids Research** 52:D1210-D1217 (2024). Free at <https://academic.oup.com/nar/article/52/D1/D1210/7416441>.
- Sondka Z, Bamford S, Cole CG, et al. *The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers.* **Nature Reviews Cancer** 18:696-705 (2018). Closed-access at the journal; the CGC list itself is freely downloadable.
- Funcotator documentation: <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Funcotator>.
- VEP documentation: <https://www.ensembl.org/info/docs/tools/vep/index.html>.
