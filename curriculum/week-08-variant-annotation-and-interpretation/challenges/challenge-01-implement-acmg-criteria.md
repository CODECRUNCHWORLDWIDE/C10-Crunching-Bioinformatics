# Challenge 1 — Implement the mechanically computable ACMG criteria

> **Educational and research use only.** This challenge implements a subset of the ACMG/AMP 2015 criteria in Python. The output is not a clinical classification. Clinical classifications require a clinical geneticist's review and incorporate criteria (functional studies, pedigree segregation, patient phenotype) that cannot be evaluated from a VCF alone. The same disclaimer that opens every Week 8 file applies.

**Estimated time:** 3 hours.
**Goal:** Extend the per-variant report from Exercise 3 with the ACMG criteria implementation from Lecture 3. For each variant, compute which of the mechanically computable criteria fire (PVS1, PM2, PM4, PP3, BA1, BS2, BP4, BP7), apply the Table 5 decision tree from Richards et al. 2015, and emit the resulting classification. Be explicit about which criteria were *skipped* because they cannot be evaluated mechanically.

This challenge is the bridge between "I have a variant annotation report" and "I have a variant classification (with disclaimers)." It is the closest a Week-8-scoped pipeline gets to a clinical-style output. The 20 criteria you do not implement are the reason the disclaimer exists.

---

## Background — Why this is harder than it looks

Two things make a faithful ACMG implementation harder than it sounds:

1. **The criteria are not symmetric.** PVS1 needs a curated gene list (which genes have LOF as a disease mechanism). PM2 needs a population-frequency threshold (which gnomAD release? popmax or global? the 2015 paper says "absent in controls," the 2018 SVI update softened to "< 1/100,000 in popmax"). PP3 changed in 2024 (Biesecker et al.) to require calibrated meta-predictor scores rather than just SIFT + PolyPhen agreement. A truly correct implementation would track all of these per-criterion refinements; a Week-8-scoped implementation pins to the 2015 framework with the 2018 PP5/BP6 deprecation and lists everything else as future work.

2. **The decision tree has subtle precedence rules.** A single BA1 hit makes the variant Benign regardless of any other criteria — but this is *because* gnomAD popmax > 0.05 essentially refutes the rarity assumption that pathogenicity rests on. If you also have PVS1 + PM2 + PP3 firing on the same variant, something is wrong with your evidence pipeline (PM2 requires popmax < 0.0001 and BA1 requires popmax > 0.05; they should be mutually exclusive). The classifier needs to handle this edge case cleanly — either flag the variant as having inconsistent evidence or trust the more specific criterion.

---

## Task

Build a Python module that takes the per-variant frame from Exercise 3 and adds:

- One column per mechanically computable ACMG criterion (`pvs1`, `pm2`, `pm4`, `pp3`, `ba1`, `bs2`, `bp4`, `bp7`). The value is a string: `"applied"`, `"not_applied"`, or `"not_evaluable"` (where `not_evaluable` means the input data was insufficient — e.g. SIFT score missing for a missense variant).
- A column listing the criteria *skipped* (`skipped_criteria`): a semicolon-separated list of ACMG criterion codes that this implementation does not compute (PS2, PS3, PS4, PM3, PM6, PP1, PP4, BS3, BS4, BP2, BP5). This list is the same for every row and exists so the report reader knows what was left out.
- A column with the final classification (`acmg_classification`): one of `Pathogenic`, `Likely_pathogenic`, `Uncertain_significance`, `Likely_benign`, `Benign`.
- A column listing the criteria that fired (`acmg_evidence`): a semicolon-separated string like `PVS1;PM2;PP3`.
- A column flagging inconsistent evidence (`acmg_warnings`): for example, "BA1+PM2 both fire" or "BS2 fires but classification is Likely_pathogenic."

### Layout

```
crunch-bio-portfolio-<yourhandle>/
└── week-08/
    └── challenge-01/
        ├── README.md             how-to-run + write-up
        ├── env.yml               conda env file
        ├── acmg_classifier.py    the core implementation
        ├── apply_acmg.py         orchestration script
        ├── data/
        │   ├── demo.vep.vcf      VEP-annotated demo VCF (from Exercise 1)
        │   └── lof_gene_set.tsv  curated LOF-disease-mechanism gene list
        └── results/
            ├── classified_variants.csv     per-variant classification table
            ├── classified_variants.html    HTML report with disclaimer
            └── run-info.json               versions and date
```

### Required functions

In `acmg_classifier.py`:

```python
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class ACMGClassification(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely_pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain_significance"
    LIKELY_BENIGN = "Likely_benign"
    BENIGN = "Benign"


@dataclass
class ACMGEvidence:
    pathogenic_very_strong: set[str] = field(default_factory=set)
    pathogenic_strong: set[str] = field(default_factory=set)
    pathogenic_moderate: set[str] = field(default_factory=set)
    pathogenic_supporting: set[str] = field(default_factory=set)
    benign_standalone: set[str] = field(default_factory=set)
    benign_strong: set[str] = field(default_factory=set)
    benign_supporting: set[str] = field(default_factory=set)
    skipped: set[str] = field(default_factory=set)
    warnings: list[str] = field(default_factory=list)


def evaluate_pvs1(consequence: str, gene_symbol: str,
                  lof_gene_set: frozenset[str]) -> bool: ...

def evaluate_pm2(gnomad_popmax_af: float | None) -> bool: ...

def evaluate_ba1(gnomad_popmax_af: float | None) -> bool: ...

def evaluate_pp3(sift_score: float | None,
                 polyphen_score: float | None) -> bool: ...

def evaluate_bp4(sift_score: float | None,
                 polyphen_score: float | None) -> bool: ...

def evaluate_bp7(consequence: str, in_splice_region: bool) -> bool: ...

def evaluate_pm4(consequence: str, in_repeat: bool) -> bool: ...

def evaluate_bs2(gnomad_hom_count: int | None,
                 disease_inheritance: str) -> bool: ...


def classify_variant(*, consequence: str, gene_symbol: str,
                     lof_gene_set: frozenset[str],
                     gnomad_popmax_af: float | None,
                     gnomad_hom_count: int | None,
                     sift_score: float | None,
                     polyphen_score: float | None,
                     in_repeat: bool = False,
                     in_splice_region: bool = False,
                     disease_inheritance: str = "recessive"
                     ) -> tuple[ACMGClassification, ACMGEvidence]: ...
```

The decision tree (from Richards et al. 2015 Table 5):

```text
Stand-alone Benign:
  >= 1 BA  -> Benign

Strong Benign:
  >= 2 BS  -> Benign
  1 BS + 1 BP  -> Likely benign
  >= 2 BP  -> Likely benign

Pathogenic:
  PVS1 + (>=1 PS) -> Pathogenic
  PVS1 + (>=2 PM) -> Pathogenic
  PVS1 + (1 PM + 1 PP) -> Pathogenic
  PVS1 + (>=2 PP) -> Pathogenic
  >= 2 PS -> Pathogenic
  1 PS + >=3 PM -> Pathogenic
  1 PS + 2 PM + >=2 PP -> Pathogenic
  1 PS + 1 PM + >=4 PP -> Pathogenic

Likely pathogenic:
  PVS1 + 1 PM -> Likely_pathogenic
  1 PS + (1-2 PM) -> Likely_pathogenic
  1 PS + >=2 PP -> Likely_pathogenic
  >= 3 PM -> Likely_pathogenic
  2 PM + >=2 PP -> Likely_pathogenic
  1 PM + >=4 PP -> Likely_pathogenic

Otherwise: Uncertain_significance
```

### The LOF gene set

Use the **ClinGen Haploinsufficiency Tier 3 gene list** as the curated source of "LOF is a known disease mechanism." It is published as a free TSV at:

```
https://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh38.bed
```

(Or the JSON version under `https://search.clinicalgenome.org/...`.) Pin to a release date; the list updates as new genes are curated.

For Week 8, a small curated subset is bundled with the curriculum at `data/lof_gene_set.tsv`. It contains ~50 well-characterized LOF-disease genes including BRCA1, BRCA2, TP53, MLH1, MSH2, APC, CFTR, FBN1, ATM, LDLR, NF1, RB1, PTEN, VHL, MEN1, BMPR2, SCN1A, COL3A1, RYR1.

### The orchestration script

```python
# apply_acmg.py
from __future__ import annotations
import json
import datetime as dt
from pathlib import Path
from typing import Any
import pandas as pd

from acmg_classifier import classify_variant, ACMGClassification


def load_lof_gene_set(path: Path) -> frozenset[str]:
    df = pd.read_csv(path, sep="\t")
    return frozenset(df["gene_symbol"].dropna().astype(str).tolist())


def apply_acmg_to_report(report: pd.DataFrame,
                         lof_gene_set: frozenset[str]) -> pd.DataFrame:
    rows = []
    for _, variant in report.iterrows():
        consequence: str = str(variant.get("consequence", ""))
        gene: str = str(variant.get("gene", ""))
        sift_score: float | None = (variant.get("sift_score")
                                     if pd.notna(variant.get("sift_score")) else None)
        pp_score: float | None = (variant.get("polyphen_score")
                                   if pd.notna(variant.get("polyphen_score")) else None)
        popmax: float | None = (float(variant.get("gnomad_popmax_af"))
                                 if pd.notna(variant.get("gnomad_popmax_af")) else None)
        classification, evidence = classify_variant(
            consequence=consequence,
            gene_symbol=gene,
            lof_gene_set=lof_gene_set,
            gnomad_popmax_af=popmax,
            gnomad_hom_count=None,
            sift_score=sift_score,
            polyphen_score=pp_score,
        )
        rows.append({
            **variant.to_dict(),
            "acmg_classification": classification.value,
            "acmg_evidence": ";".join(
                sorted(evidence.pathogenic_very_strong
                       | evidence.pathogenic_strong
                       | evidence.pathogenic_moderate
                       | evidence.pathogenic_supporting
                       | evidence.benign_standalone
                       | evidence.benign_strong
                       | evidence.benign_supporting)
            ),
            "skipped_criteria": ";".join(sorted(evidence.skipped)),
            "acmg_warnings": ";".join(evidence.warnings),
        })
    return pd.DataFrame(rows)
```

### Acceptance criteria

- [ ] `acmg_classifier.py` exposes a `classify_variant(...)` function and an `ACMGClassification` enum.
- [ ] The function is type-hinted on every parameter (use `Any` where the cyvcf2 stub is incomplete).
- [ ] All eight implemented criteria have a one-paragraph docstring referencing the Richards et al. 2015 criterion definition.
- [ ] `apply_acmg.py` runs end to end on the bundled `data/demo.vep.vcf` (or, in `--no-network` mode, on the bundled demo rows from Exercise 3).
- [ ] `results/classified_variants.csv` exists with at least the columns from Exercise 3 plus `acmg_classification`, `acmg_evidence`, `skipped_criteria`, `acmg_warnings`.
- [ ] `results/classified_variants.html` exists with the disclaimer, the metadata box, the per-variant table, and the `skipped_criteria` made prominent (so the reader knows which criteria were not evaluated).
- [ ] The classifier produces *consistent* classifications on the demo set:
  - BRCA1 missense at chr17:43094077 (popmax 0.00018, SIFT 0.02, PolyPhen 0.94) -> Likely_pathogenic by PM2 + PP3.
  - BRCA2 stop_gained at chr13:32398489 (popmax 0.012, common) -> Benign by BA1 (a stand-alone hit; the BA1 logic overrides PVS1 here, and the classifier emits a warning).
  - CFTR synonymous at chr7:117559590 (popmax 0.18, BP7-eligible) -> Benign by BA1; the BP7 also fires, recorded in evidence.
  - APC frameshift at chr5:112815473 (absent gnomAD, in LOF gene set) -> Likely_pathogenic by PVS1 + PM2.
  - TP53 stop_gained at chr17:7674220 (absent gnomAD, in LOF gene set) -> Likely_pathogenic by PVS1 + PM2.
  - LDLR missense at chr19:11116107 (absent gnomAD, SIFT 0.00, PolyPhen 1.00, in LOF gene set) -> Likely_pathogenic by PM2 + PP3.
- [ ] `results/run-info.json` records the date, the VEP cache version, the LOF gene list source and date.
- [ ] `README.md` documents the run, the expected outputs, and a paragraph explicitly listing which criteria the classifier *does not* implement and why.

### Stretch goals

- **Add PM5** by joining against a ClinVar table of known-pathogenic missense indexed by `gene_symbol + protein_position`. Fire PM5 when the variant is missense at a residue where a different missense is recorded as Pathogenic in ClinVar.
- **Add PVS1 refinements** from Walsh et al. 2023 (the 2023 SVI update for PVS1). The refinements consider exon location, NMD-rule status (last-exon escapes NMD), and the presence of alternative start codons.
- **Add CADD as a third PP3/BP4 predictor.** CADD score >= 20 plus SIFT < 0.05 plus PolyPhen > 0.85 is the modern PP3_Strong threshold; agreement across all three is a stronger signal than the original SIFT+PolyPhen pair.
- **Integrate SpliceAI** for splice-region variants. SpliceAI score >= 0.5 is the modern "predicted splice-altering" cutoff; it catches variants in the splice region that VEP buckets as LOW.
- **Generate a ClinGen-style summary line per variant** using their published variant interpretation summary template. The template combines the evidence and the classification into a single sentence suitable for a clinical report.

### Submission

Push the `challenge-01/` directory to your portfolio repo with a commit message like:

```
Week 8 challenge 1: ACMG classifier on demo VCF, 8 criteria implemented,
20 explicitly skipped.
```

Open a PR to the curriculum repo with a brief description.

---

## Validation: what does "right" look like?

A correct implementation produces consistent classifications on the demo set, with the BA1-vs-PM2 inconsistency on BRCA2 chr13:32398489 caught and flagged. A *good* implementation in addition:

- Uses `Optional[float]` or `float | None` consistently for fields that may be missing.
- Handles the "consequence not in the SO vocabulary" edge case (some tools emit free-text consequence strings; the classifier should not crash on those — it should log a warning and skip the criteria that need the consequence).
- Includes a test file (`test_acmg_classifier.py`) with a `pytest` suite covering each criterion in isolation plus the decision-tree edges.

A *great* implementation also documents the **boundary** between automated and manual review in its README. The boundary is the point of the challenge.
