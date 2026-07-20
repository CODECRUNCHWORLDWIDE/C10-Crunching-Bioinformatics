# Challenge 2 — Build a pharmacogenomics-style report

> **Educational and research use only.** This challenge teaches the *mechanics* of mapping variants in pharmacogenes to drug-response recommendations from the published CPIC tables. The output is not a clinical prescription, dosing recommendation, or therapeutic decision. Pharmacogenomic testing in a clinical context requires CLIA-certified star-allele calling, a clinical pharmacist's review, and integration with the patient's clinical context. The same disclaimer that opens every Week 8 file applies.

**Estimated time:** 3 hours.
**Goal:** Take a VCF restricted to the CPIC tier-1 pharmacogenes (`CYP2D6`, `CYP2C19`, `CYP2C9`, `TPMT`, `DPYD`, `SLCO1B1`, `VKORC1`, `CYP3A5`, `UGT1A1`, `HLA-B`) and produce a per-drug recommendation report built from PharmGKB's free-tier knowledge base. The pipeline maps each variant to a known PharmGKB-curated star allele (where possible), looks up the published CPIC recommendation by phenotype, and emits a per-drug recommendation card.

This challenge is the "what other applied bioinformatics field uses these patterns" complement to the disease-variant work in Exercise 3 and Challenge 1. Pharmacogenomics has its own evidence framework (CPIC, not ACMG), its own canonical database (PharmGKB), and its own "the automated output is not the clinical decision" boundary.

---

## Background — Why this is different from ACMG

ACMG/AMP asks "is this variant causing disease?" Pharmacogenomics asks "how will this patient metabolize this drug?" The two questions differ in three important ways:

1. **The unit of analysis is a star allele, not a point variant.** A star allele like `CYP2D6*4` is a defined haplotype (combination of variants on one chromosome copy). Calling a star allele requires knowing the **phase** of multiple variants — which ones are on the same chromosome copy. Short-read sequencing struggles with this beyond a few hundred base pairs. Specialized tools (Aldy, Stargazer, Cyrius — all free) do star-allele calling from short reads; for Week 8 we use a simplified approximation: we look up *known* PharmGKB-indexed variants (by rsID) and assume the genotype indicates the star allele.

2. **The output is a phenotype + recommendation, not a classification.** Star alleles map to metabolizer phenotypes (ultrarapid / normal / intermediate / poor). The CPIC published tables map each phenotype to a per-drug recommendation: "consider alternative therapy," "reduce dose by 50%," "no change to standard dosing," etc.

3. **CPIC tier-1 evidence is much stronger than typical ACMG VUS evidence.** The CPIC tier-1 gene-drug pairs have published, peer-reviewed dosing guidelines based on clinical pharmacokinetic studies. The recommendations are conservative and standardized. They are the closest thing in clinical bioinformatics to "the textbook tells you what to do."

But — and this is the disclaimer — the recommendation tables are *for clinical pharmacists to interpret*. The bioinformatics pipeline produces the input rows; the pharmacist makes the decision. The teaching point of Challenge 2 is the same as the teaching point of Exercise 3 and Challenge 1: the automated output is not the clinical decision.

---

## Task

Build a Python module that takes a VCF restricted to the pharmacogenes and produces a per-drug recommendation report.

### Layout

```
crunch-bio-portfolio-<yourhandle>/
└── week-08/
    └── challenge-02/
        ├── README.md
        ├── env.yml
        ├── pgx_report.py             the orchestration script
        ├── pharmgkb_client.py        PharmGKB API client
        ├── data/
        │   ├── pgx_demo.vcf          a small demo VCF with ~10 pharmacogene variants
        │   ├── pgx_demo.vep.vcf      VEP-annotated version
        │   └── cpic_tier1_drugs.tsv  curated list of CPIC tier-1 gene-drug pairs
        └── results/
            ├── pgx_recommendations.csv   per-variant + per-drug recommendation table
            ├── pgx_recommendations.html  HTML report with disclaimer
            └── run-info.json             versions and date
```

### Required functions

In `pharmgkb_client.py`:

```python
from __future__ import annotations
import time
from typing import Any

import requests


PHARMGKB_BASE = "https://api.pharmgkb.org/v1"
PHARMGKB_SLEEP_SEC = 0.3  # rate-limit politeness


def fetch_pharmgkb_gene(symbol: str, timeout: float = 30.0) -> dict[str, Any]:
    """Fetch the PharmGKB gene record for a gene symbol.

    Returns the parsed JSON. Empty dict if the gene is not in PharmGKB.
    """
    ...


def fetch_pharmgkb_variant(rsid: str, timeout: float = 30.0) -> dict[str, Any]:
    """Fetch the PharmGKB variant annotation for a dbSNP rsID."""
    ...


def fetch_cpic_guideline(gene: str, drug: str,
                          timeout: float = 30.0) -> dict[str, Any]:
    """Fetch the CPIC guideline for a gene-drug pair."""
    ...
```

In `pgx_report.py`:

```python
from __future__ import annotations
from pathlib import Path
from typing import Any

import pandas as pd


def variants_to_pharmgkb_rows(annotated_vcf: Path,
                               pharmgkb_cache: dict[str, Any]) -> pd.DataFrame:
    """For each variant in the VCF, query PharmGKB by rsID.

    For variants with a PharmGKB-indexed match, emit one row per
    (gene, drug) pair. The row contains the variant coordinates, the
    PharmGKB clinical annotation level (1A / 1B / 2A / 2B / 3 / 4),
    the metabolizer phenotype assigned, and the CPIC recommendation.
    """
    ...


def build_pgx_html_report(rows: pd.DataFrame, output_path: Path,
                           run_info: dict[str, Any]) -> None:
    """Render the per-drug recommendation table as HTML.

    Must include the disclaimer, the run metadata, and a per-gene
    section with all relevant drug recommendations.
    """
    ...
```

### The demo VCF

The bundled `data/pgx_demo.vcf` contains ~10 variants in CPIC tier-1 pharmacogenes:

| Gene    | rsID         | Variant (HGVS)           | Star allele |
|---------|--------------|---------------------------|-------------|
| CYP2D6  | rs3892097    | c.1846G>A                 | CYP2D6*4    |
| CYP2D6  | rs1065852    | c.100C>T                  | CYP2D6*10   |
| CYP2C19 | rs4244285    | c.681G>A                  | CYP2C19*2   |
| CYP2C19 | rs12248560   | c.-806C>T                 | CYP2C19*17  |
| CYP2C9  | rs1799853    | c.430C>T                  | CYP2C9*2    |
| CYP2C9  | rs1057910    | c.1075A>C                 | CYP2C9*3    |
| TPMT    | rs1142345    | c.719A>G                  | TPMT*3C     |
| SLCO1B1 | rs4149056    | c.521T>C                  | SLCO1B1*5   |
| VKORC1  | rs9923231    | c.-1639G>A                | -           |
| DPYD    | rs3918290    | c.1905+1G>A               | DPYD*2A     |

All ten are well-known CPIC tier-1 variants with published recommendations. Note that one variant per position is shown here; the actual genotype (homozygous, heterozygous, or absent) plus the phase to other variants in the same gene determines the star-allele diplotype.

### The recommendation logic (simplified for Week 8)

For each variant in the input VCF:

1. Extract the rsID (from the VEP `Existing_variation` field).
2. Query PharmGKB by rsID. If no match, skip the variant.
3. From the PharmGKB response, extract the matched star allele(s) and the matched gene.
4. Look up the CPIC guideline for the matched gene plus each CPIC tier-1 drug for that gene.
5. Assign a phenotype based on the genotype (heterozygous vs homozygous):
   - Heterozygous for a loss-of-function (`*4`, `*5`, etc.) -> "Intermediate metabolizer."
   - Homozygous for a LOF -> "Poor metabolizer."
   - Heterozygous for a gain-of-function (`*17`) -> "Ultrarapid metabolizer" (CYP2C19 only).
   - Wild-type / `*1` -> "Normal metabolizer."

(Real star-allele calling does much more than this; the simplified rule is fine for Week 8's didactic purpose.)

6. Map the phenotype to the per-drug recommendation using the CPIC table.

### Acceptance criteria

- [ ] `pgx_report.py` runs end to end on the bundled `data/pgx_demo.vep.vcf` and produces `results/pgx_recommendations.csv`.
- [ ] The CSV has at least the columns: `gene, rsid, variant_hgvs, star_allele, genotype, phenotype, drug, recommendation, evidence_level, source`.
- [ ] The HTML report (`results/pgx_recommendations.html`) groups recommendations by gene, prominently features the disclaimer, and includes a per-recommendation evidence level (1A is the strongest CPIC tier).
- [ ] At least 6 of the 10 demo variants produce at least one drug recommendation. (Some pharmacogenes have specific phenotype-genotype tables that this simplified pipeline cannot resolve; those should be reported as "indeterminate phenotype, see specialist tool.")
- [ ] `results/run-info.json` records: run date, PharmGKB API version (from response headers), CPIC guidelines version (where available), and the warning that this pipeline does *not* do real star-allele calling.
- [ ] `README.md` clearly states the limits: no real star-allele calling, no phasing, no consideration of compound heterozygosity, no integration with the patient's clinical context.

### Stretch goals

- **Integrate `Aldy` or `Stargazer`** for real star-allele calling from a BAM (not just a VCF). This requires the upstream BAM file, which is not bundled with the curriculum repo; for the stretch goal, use the GIAB NA12878 BAM.
- **Add the HLA-B*57:01 detection** for abacavir contraindication. HLA typing from short-read WGS is its own subfield (HLAforest, OptiType, HLA-LA — all free); the stretch goal integrates one of these and produces the HLA-B*57:01 yes/no flag.
- **Build a one-page per-patient summary** that lists every CPIC tier-1 drug the patient takes (input as a free-text drug list), the matched gene(s), and the recommendation. This is closest to how a real pharmacy-informed pharmacogenomics report looks.
- **Add the FDA Table of Pharmacogenomic Biomarkers** (free) as a second source. The FDA list overlaps with CPIC tier-1 but has different evidence criteria and slightly different recommendations.

### Submission

Push the `challenge-02/` directory to your portfolio repo with a commit message like:

```
Week 8 challenge 2: PharmGKB-driven PGx report on demo VCF,
no real star-allele calling, disclaimers prominent.
```

Open a PR to the curriculum repo with a brief description.

---

## Validation: what does "right" look like?

A correct implementation, given the demo VCF with heterozygous CYP2D6 rs3892097 (genotype `0/1`):

- Identifies rs3892097 as CYP2D6*4.
- Assigns the heterozygous CYP2D6*4 / *1 diplotype.
- Maps the diplotype to "Intermediate metabolizer" phenotype.
- Looks up the CPIC recommendation for CYP2D6 + codeine: "Avoid use of codeine; use morphine or non-opioid analgesics instead."
- Reports the evidence level: 1A (the highest CPIC tier).

A *good* implementation labels the genotype source (the VCF GT field), the database version (PharmGKB date), and the limit ("real star-allele calling requires Aldy or Stargazer; this pipeline uses a simplified rsID lookup").

A *great* implementation produces the per-patient summary page and clearly flags the boundary between automated lookup and clinical decision. The boundary is the point of the challenge.
