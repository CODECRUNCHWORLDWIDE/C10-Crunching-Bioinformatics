# Lecture 3 — ACMG Classification and Pharmacogenomics

> **Educational and research use only.** The ACMG/AMP framework and the CPIC pharmacogenomics guidelines describe how clinical labs combine evidence into a classification. Implementing a subset of those rules in code does not produce a clinical classification. Clinical classifications require a clinician's review of the patient, the pedigree, and the case. The point of this lecture is to make you fluent in the *mechanics* — what each criterion means, which can be computed from a VCF + databases, and which require human judgment — not to enable you to assign a classification.

> **Duration:** ~3 hours of reading + a Python implementation session.
> **Outcome:** You can list the 28 ACMG criteria, identify the ~8 that are mechanically computable from a VCF + the databases of Lecture 2, implement a minimal ACMG classifier in Python, query PharmGKB programmatically, and explain the boundary between automated and clinician-led interpretation.

If you only remember one thing from this lecture, remember this:

> **The ACMG/AMP 2015 framework defines 28 criteria across four evidence categories. About 8 of those are mechanically computable from a VCF plus VEP / SnpEff / ClinVar / gnomAD / SIFT / PolyPhen. The remaining 20 require manual review of literature, functional studies, pedigrees, or patient phenotype. An automated pipeline that claims to "apply the ACMG framework" without naming which criteria it computes and which it skips is either overclaiming or quietly broken. The Week 8 exercises name everything they compute and refuse to claim anything they do not.**

Lecture 2 produced a per-variant frame with consequence, impact, gnomAD AF, popmax, ClinVar CLNSIG, and SIFT / PolyPhen scores. Lecture 3 maps that frame into the ACMG criteria.

---

## 1. The ACMG/AMP 2015 framework in one page

Richards et al. 2015 (*Genetics in Medicine* 17:405; free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/>) define the framework that essentially every clinical lab in the United States uses for variant classification. The 28 criteria are grouped into four categories by evidence type:

- **Population data** — frequency in healthy populations.
- **Computational and predictive data** — in-silico predictors and consequence type.
- **Functional data** — wet-lab assays and curated databases.
- **Segregation and case data** — family pedigree, case-control, de novo.

Within each category, criteria are graded by strength:

- **PVS** — Pathogenic Very Strong (one criterion: PVS1).
- **PS** — Pathogenic Strong (PS1-PS4).
- **PM** — Pathogenic Moderate (PM1-PM6).
- **PP** — Pathogenic Supporting (PP1-PP5).
- **BA** — Benign Stand-Alone (one criterion: BA1; a single BA1 hit is sufficient to classify a variant as Benign).
- **BS** — Benign Strong (BS1-BS4).
- **BP** — Benign Supporting (BP1-BP7).

The 28 criteria combine via a decision table (Table 5 of Richards et al. 2015) to produce one of five classifications:

| Pathogenic | (i) 1 PVS + (a) >= 1 PS, (b) >= 2 PM, (c) 1 PM + 1 PP, (d) >= 2 PP; (ii) >= 2 PS; (iii) 1 PS + (a) >= 3 PM, (b) 2 PM + >= 2 PP, (c) 1 PM + >= 4 PP |
| Likely Pathogenic | (i) 1 PVS + 1 PM; (ii) 1 PS + 1-2 PM; (iii) 1 PS + >= 2 PP; (iv) >= 3 PM; (v) 2 PM + >= 2 PP; (vi) 1 PM + >= 4 PP |
| Benign | (i) 1 BA1; (ii) >= 2 BS |
| Likely Benign | (i) 1 BS + 1 BP; (ii) >= 2 BP |
| Uncertain Significance (VUS) | All other combinations, including contradictory pathogenic + benign criteria |

Two important observations:

- **A single BA1 hit is sufficient to call Benign.** If a variant has gnomAD popmax > 0.05, the framework says "the variant is too common to cause a rare Mendelian disease, regardless of any other evidence." (BA1 has been refined since 2015; the 2019 SVI recommendation lowered the threshold for some specific disease contexts to 0.02 or 0.01, but the 5% default is still the standard.)
- **VUS is the default for ambiguous evidence.** Most variants — > 90% of incidental findings in a typical exome — land in VUS. VUS is not "uninformative"; it is "insufficient evidence at this time." A VUS variant should be re-examined as evidence accumulates.

---

## 2. The 28 criteria, with annotations

The table below restates each criterion in plain English, lists what data you need to apply it, and flags whether it is mechanically computable.

### Pathogenic criteria

| Criterion | Meaning | Data required | Computable? |
|-----------|---------|---------------|-------------|
| **PVS1** | Null variant (nonsense, frameshift, canonical +/- 1,2 splice site, initiation codon, single-exon or multi-exon deletion) in a gene where loss-of-function is a known disease mechanism | VEP consequence (`IMPACT=HIGH` of the appropriate type) + a curated LOF-disease-mechanism gene list (ClinGen Haploinsufficiency tier 3 is the canonical reference) | Yes (with gene list) |
| **PS1** | Same amino acid change as a previously reported pathogenic variant (different nucleotide) | ClinVar Pathogenic table indexed by HGVS protein change | Yes |
| **PS2** | De novo (paternity and maternity confirmed) in a patient with disease and no family history | Trio VCF + parentage verification | No (clinical) |
| **PS3** | Well-established in vivo or in vitro functional studies | Literature review of functional assays | No (literature) |
| **PS4** | Prevalence in affected individuals significantly increased over controls | Case-control study results | No (epidemiology) |
| **PM1** | Located in a mutational hotspot and/or critical and well-established functional domain (e.g. active site, binding site) | Curated hotspot / domain table (UniProt, Pfam, ClinGen) | Partial |
| **PM2** | Absent from controls (or extremely low frequency: gnomAD popmax < 0.0001 is the standard threshold; SVI 2018 update softened to "PM2_Supporting" for moderate rarity) | gnomAD AF | **Yes** |
| **PM3** | For recessive disorders, detected in trans with a pathogenic variant | Phasing + ClinVar | No (phasing) |
| **PM4** | Protein length change due to in-frame deletions/insertions in a non-repeat region or stop-loss | VEP consequence | Yes |
| **PM5** | Novel missense at an amino acid residue where a different missense has been seen previously as pathogenic | ClinVar Pathogenic table indexed by amino acid position | Yes |
| **PM6** | Assumed de novo (without parentage confirmation) | Trio VCF | No |
| **PP1** | Cosegregation with disease in multiple affected family members | Pedigree | No (clinical) |
| **PP2** | Missense variant in a gene that has a low rate of benign missense variation and where missense is a common disease mechanism | Gene constraint scores (gnomAD missense Z-score, ExAC obs/exp) + a curated gene list | Partial |
| **PP3** | Multiple lines of computational evidence support a deleterious effect (conservation, in-silico predictors) | SIFT + PolyPhen + CADD (or equivalent) | **Yes** |
| **PP4** | Patient phenotype highly specific for a disease with a single genetic etiology | HPO phenotype + OMIM lookup | No (clinical) |
| **PP5** | Reputable source reports as pathogenic without independent evidence | ClinVar | Yes (deprecated as of 2018 SVI update — see below) |

### Benign criteria

| Criterion | Meaning | Data required | Computable? |
|-----------|---------|---------------|-------------|
| **BA1** | Allele frequency > 5% in any general population (gnomAD popmax > 0.05) | gnomAD popmax | **Yes** |
| **BS1** | Allele frequency greater than expected for the disorder | gnomAD popmax + estimated disease prevalence | Yes (with prevalence) |
| **BS2** | Observed homozygous in a healthy adult for a recessive condition (or hemizygous for X-linked) | gnomAD homozygous count | Partial |
| **BS3** | Well-established functional studies show no damaging effect | Literature review | No (literature) |
| **BS4** | Lack of segregation in affected family members | Pedigree | No |
| **BP1** | Missense variant in a gene where only truncating variants cause disease | Gene-specific curation | Partial |
| **BP2** | Observed in trans with pathogenic for dominant disorder, or in cis with pathogenic | Phasing + ClinVar | No (phasing) |
| **BP3** | In-frame indel in a repetitive region without known function | VEP consequence + repeat annotation | Yes |
| **BP4** | Multiple lines of computational evidence suggest no impact | SIFT + PolyPhen + CADD (all suggesting tolerant) | **Yes** |
| **BP5** | Variant observed in a case with an alternate molecular basis for disease | Case-level data | No |
| **BP6** | Reputable source reports as benign without independent evidence | ClinVar | Yes (deprecated as of 2018 SVI update) |
| **BP7** | Synonymous (silent) variant where no impact on splicing is predicted and the variant is not conserved | VEP consequence (synonymous_variant) + SpliceAI / MaxEntScan | Yes |

### A note on PP5 and BP6

The 2018 ACMG SVI update (Biesecker and Harrison 2018, *Genetics in Medicine* 20:1687) recommended **discontinuing** the use of PP5 and BP6 because they are circular: "ClinVar says Pathogenic" is not independent evidence when you are trying to apply the ACMG framework that ClinVar itself uses. The 2024 PP3/BP4 update (Biesecker et al. 2024) made PP3/BP4 stricter by requiring quantitative thresholds on predictor scores. Week 8 uses the 2015 framework as the default with the 2018 PP5/BP6 deprecation but does not yet incorporate the 2024 PP3/BP4 refinements; those would require integrating CADD and a calibrated meta-predictor.

---

## 3. The mechanically computable subset

For a Week-8-scoped pipeline, the cleanly mechanically computable criteria from a VCF + the databases of Lecture 2 are:

| Criterion | Implementation rule |
|-----------|---------------------|
| **PVS1** | VEP IMPACT = HIGH (specifically: stop_gained, frameshift_variant, splice_donor_variant, splice_acceptor_variant, start_lost) AND gene is in the ClinGen Haploinsufficiency tier 3 list |
| **PS1** | VEP `Amino_acids` field matches a known Pathogenic missense at the same position in ClinVar |
| **PM2** | gnomAD popmax AF < 0.0001 (or variant absent from gnomAD entirely) |
| **PM4** | VEP consequence in (inframe_insertion, inframe_deletion, stop_lost) AND not in a repeat region |
| **PM5** | VEP missense at an amino acid position where a different missense is recorded as Pathogenic in ClinVar |
| **PP3** | SIFT score < 0.05 (deleterious) AND PolyPhen score > 0.85 (probably_damaging) |
| **BA1** | gnomAD popmax AF > 0.05 (5%) |
| **BS1** | gnomAD popmax AF > expected for the disorder (requires per-disease prevalence) |
| **BS2** | gnomAD homozygous count > 0 for a recessive disorder |
| **BP3** | VEP consequence is inframe indel AND the indel is in a repeat-annotated region (RepeatMasker) |
| **BP4** | SIFT score > 0.05 (tolerated) AND PolyPhen score < 0.5 (benign) |
| **BP7** | VEP consequence is synonymous_variant AND not in a splice region |

Eight of these are crisply mechanically computable from the databases at hand: **PVS1, PM2, PM4, PP3, BA1, BS2, BP4, BP7**. Three more — **PS1, PM5, BS1** — are computable with an additional table (ClinVar missense lookup, disease prevalence estimate). **BP3** requires RepeatMasker output. The other ~17 criteria are not mechanically computable from a VCF.

Week 8's exercises implement the eight cleanly computable criteria. The mini-project's ACMG classifier reports each criterion as `applied / not_applied / not_applicable_mechanical`, where the third category explicitly flags that the criterion was *not evaluated* because the data was not available.

---

## 4. A minimal ACMG classifier in Python

The classifier takes the per-variant frame from Lecture 2 plus a gene-list and a SpliceAI/CADD score (optional) and returns the set of criteria that fired plus the resulting classification.

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
    """The set of ACMG criteria that fired for a single variant.

    Each criterion is stored as the string code ("PVS1", "PM2", ...).
    The `skipped` set records criteria that could not be mechanically
    evaluated and which a human reviewer must consider separately.
    """
    pathogenic_very_strong: set[str] = field(default_factory=set)  # PVS
    pathogenic_strong: set[str] = field(default_factory=set)       # PS
    pathogenic_moderate: set[str] = field(default_factory=set)     # PM
    pathogenic_supporting: set[str] = field(default_factory=set)   # PP
    benign_standalone: set[str] = field(default_factory=set)       # BA
    benign_strong: set[str] = field(default_factory=set)           # BS
    benign_supporting: set[str] = field(default_factory=set)       # BP
    skipped: set[str] = field(default_factory=set)

    def all_pathogenic(self) -> list[str]:
        return sorted(
            self.pathogenic_very_strong
            | self.pathogenic_strong
            | self.pathogenic_moderate
            | self.pathogenic_supporting
        )

    def all_benign(self) -> list[str]:
        return sorted(
            self.benign_standalone
            | self.benign_strong
            | self.benign_supporting
        )


HIGH_IMPACT_LOF = frozenset({
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "transcript_ablation",
})


def evaluate_pvs1(consequence: str, gene_symbol: str,
                  lof_gene_set: frozenset[str]) -> bool:
    """PVS1: null variant in a gene where LOF is a known disease mechanism."""
    if not consequence:
        return False
    is_null = any(term in consequence for term in HIGH_IMPACT_LOF)
    return is_null and gene_symbol in lof_gene_set


def evaluate_pm2(gnomad_popmax_af: float | None) -> bool:
    """PM2: absent or extremely rare in gnomAD (popmax AF < 0.0001)."""
    if gnomad_popmax_af is None:
        # Absent from gnomAD entirely.
        return True
    return gnomad_popmax_af < 1e-4


def evaluate_ba1(gnomad_popmax_af: float | None) -> bool:
    """BA1: gnomAD popmax AF > 0.05 (5%) — stand-alone benign."""
    if gnomad_popmax_af is None:
        return False
    return gnomad_popmax_af > 0.05


def evaluate_pp3(sift_score: float | None, polyphen_score: float | None) -> bool:
    """PP3: multiple lines of computational evidence support deleterious."""
    if sift_score is None or polyphen_score is None:
        return False
    return sift_score < 0.05 and polyphen_score > 0.85


def evaluate_bp4(sift_score: float | None, polyphen_score: float | None) -> bool:
    """BP4: multiple lines of computational evidence suggest no impact."""
    if sift_score is None or polyphen_score is None:
        return False
    return sift_score > 0.05 and polyphen_score < 0.5


def evaluate_bp7(consequence: str, splice_region: bool) -> bool:
    """BP7: synonymous and not in a splice region."""
    return "synonymous_variant" in consequence and not splice_region


def evaluate_pm4(consequence: str, in_repeat: bool) -> bool:
    """PM4: in-frame indel (or stop_lost) in a non-repetitive region."""
    is_length_changing = any(term in consequence for term in
                              ("inframe_insertion", "inframe_deletion", "stop_lost"))
    return is_length_changing and not in_repeat


def evaluate_bs2(gnomad_hom_count: int | None,
                 disease_inheritance: str = "recessive") -> bool:
    """BS2: observed homozygous in a healthy adult for a recessive disorder.

    Conservative implementation: any homozygous count > 0 for a
    recessive disorder fires BS2.
    """
    if gnomad_hom_count is None:
        return False
    if disease_inheritance != "recessive":
        return False
    return gnomad_hom_count > 0


def classify_variant(consequence: str,
                     gene_symbol: str,
                     lof_gene_set: frozenset[str],
                     gnomad_popmax_af: float | None,
                     gnomad_hom_count: int | None,
                     sift_score: float | None,
                     polyphen_score: float | None,
                     in_repeat: bool = False,
                     in_splice_region: bool = False,
                     disease_inheritance: str = "recessive") -> tuple[ACMGClassification, ACMGEvidence]:
    """Apply the mechanically computable ACMG criteria.

    Returns the resulting classification and the evidence set.
    Criteria that cannot be mechanically evaluated are listed in
    `evidence.skipped` so a downstream report can flag them for
    manual review.
    """
    evidence = ACMGEvidence()

    if evaluate_pvs1(consequence, gene_symbol, lof_gene_set):
        evidence.pathogenic_very_strong.add("PVS1")
    if evaluate_pm2(gnomad_popmax_af):
        evidence.pathogenic_moderate.add("PM2")
    if evaluate_pm4(consequence, in_repeat):
        evidence.pathogenic_moderate.add("PM4")
    if evaluate_pp3(sift_score, polyphen_score):
        evidence.pathogenic_supporting.add("PP3")
    if evaluate_ba1(gnomad_popmax_af):
        evidence.benign_standalone.add("BA1")
    if evaluate_bs2(gnomad_hom_count, disease_inheritance):
        evidence.benign_strong.add("BS2")
    if evaluate_bp4(sift_score, polyphen_score):
        evidence.benign_supporting.add("BP4")
    if evaluate_bp7(consequence, in_splice_region):
        evidence.benign_supporting.add("BP7")

    # Mark the non-mechanical criteria as skipped (informational).
    for cid in ("PS2", "PS3", "PS4", "PM3", "PM6",
                 "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"):
        evidence.skipped.add(cid)

    return _resolve_classification(evidence), evidence


def _resolve_classification(ev: ACMGEvidence) -> ACMGClassification:
    """Apply Table 5 of Richards et al. 2015 to the evidence set."""
    pvs = len(ev.pathogenic_very_strong)
    ps = len(ev.pathogenic_strong)
    pm = len(ev.pathogenic_moderate)
    pp = len(ev.pathogenic_supporting)
    ba = len(ev.benign_standalone)
    bs = len(ev.benign_strong)
    bp = len(ev.benign_supporting)

    # Stand-alone benign overrides everything.
    if ba >= 1:
        return ACMGClassification.BENIGN
    # Pathogenic combinations.
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
        return ACMGClassification.PATHOGENIC
    if ps >= 2:
        return ACMGClassification.PATHOGENIC
    if ps >= 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)):
        return ACMGClassification.PATHOGENIC
    # Likely pathogenic.
    if pvs >= 1 and pm >= 1:
        return ACMGClassification.LIKELY_PATHOGENIC
    if ps >= 1 and (1 <= pm <= 2):
        return ACMGClassification.LIKELY_PATHOGENIC
    if ps >= 1 and pp >= 2:
        return ACMGClassification.LIKELY_PATHOGENIC
    if pm >= 3:
        return ACMGClassification.LIKELY_PATHOGENIC
    if pm >= 2 and pp >= 2:
        return ACMGClassification.LIKELY_PATHOGENIC
    if pm >= 1 and pp >= 4:
        return ACMGClassification.LIKELY_PATHOGENIC
    # Benign / likely benign.
    if bs >= 2:
        return ACMGClassification.BENIGN
    if bs >= 1 and bp >= 1:
        return ACMGClassification.LIKELY_BENIGN
    if bp >= 2:
        return ACMGClassification.LIKELY_BENIGN
    return ACMGClassification.UNCERTAIN_SIGNIFICANCE
```

A few notes on this implementation:

- The classifier is *mechanically applying* a published rule table. It is not making a clinical judgment. The same evidence set with manual review by a clinical geneticist may resolve differently because the clinician can incorporate PS2/PS3/PS4 (which the classifier cannot).
- The `skipped` set is the most important field of the output. A report that says "Likely Pathogenic by PVS1 + PM2 + PP3, with PS2/PS3/PS4/PP1/PP4 not evaluated" is honest. A report that says "Likely Pathogenic" without naming the missed criteria is overclaiming.
- The classifier assumes the gnomAD popmax and SIFT/PolyPhen are floats or None. None means "not available" — for example, SIFT is undefined for synonymous variants, so the PP3 and BP4 criteria simply do not fire.
- The 2018 SVI updates to PVS1 (Walsh et al. 2023 was an extensive refinement) are not implemented in this minimal classifier. The "is the gene in the LOF gene set" question is a coarse approximation of the full PVS1 decision tree. Real clinical pipelines use the ClinGen-published Sequence Variant Interpretation working group PVS1 tool (free, available at <https://www.clinicalgenome.org/working-groups/sequence-variant-interpretation/>).

---

## 5. SIFT and PolyPhen-2 in detail

The two computational missense predictors that PP3/BP4 rely on.

### SIFT (Ng and Henikoff 2003)

SIFT (**S**orting **I**ntolerant **F**rom **T**olerant) builds a multiple sequence alignment of orthologs of the protein and scores how conserved each residue is. A non-conservative substitution at a highly conserved residue gets a low SIFT score (close to 0.0, "deleterious"); a substitution at a poorly conserved residue gets a high score (close to 1.0, "tolerated"). The standard cutoff is 0.05.

| SIFT score | Prediction |
|------------|------------|
| < 0.05     | deleterious |
| >= 0.05    | tolerated |

SIFT is *sequence-only*. It does not see protein structure. The strength of SIFT is its simplicity; the limitation is that some functionally important residues (catalytic residues, binding-site residues) are not well-flagged by conservation alone.

### PolyPhen-2 (Adzhubei et al. 2010)

PolyPhen-2 (**Pol**ymorphism **Phen**otyping, version 2) extends SIFT-style conservation analysis with structural features: solvent accessibility, secondary structure type, predicted stability change, predicted impact on protein-protein interfaces. The output is a score in [0, 1] plus a categorical prediction:

| PolyPhen-2 score | Prediction |
|------------------|------------|
| 0.0-0.452        | benign |
| 0.453-0.957      | possibly_damaging |
| > 0.957          | probably_damaging |

(The cutoffs vary slightly by training set; the values above are the standard HumDiv-trained model. PolyPhen-2 also publishes a HumVar-trained model with slightly different cutoffs.)

The standard ACMG PP3 threshold for PolyPhen-2 is "probably_damaging," i.e. score > 0.85. The standard BP4 threshold is "benign," i.e. score < 0.5.

### Using SIFT + PolyPhen-2 together

PP3 fires when SIFT and PolyPhen agree that the variant is deleterious. BP4 fires when they agree that it is tolerated. If they disagree — SIFT says deleterious, PolyPhen says benign, or vice versa — neither PP3 nor BP4 fires. Disagreement happens for ~10-15% of missense variants and is a flag for manual inspection.

VEP's `--sift b --polyphen b` flags emit both scores in the CSQ output as `deleterious(0.02)` and `probably_damaging(0.94)` (the score in parentheses). Parsing the float out:

```python
import re


def parse_sift(sift_field: str) -> tuple[str | None, float | None]:
    """Parse VEP's SIFT field like 'deleterious(0.02)' or 'tolerated(0.5)'.

    Returns (prediction, score). Both None if the field is empty.
    """
    if not sift_field:
        return None, None
    match = re.match(r"([a-z_]+)\(([0-9.]+)\)", sift_field)
    if not match:
        return None, None
    return match.group(1), float(match.group(2))


def parse_polyphen(pp_field: str) -> tuple[str | None, float | None]:
    """Parse VEP's PolyPhen field like 'probably_damaging(0.94)'."""
    return parse_sift(pp_field)  # same format
```

---

## 6. PharmGKB and pharmacogenomics

The third clinical use of variant annotation — after population frequency and disease pathogenicity — is **pharmacogenomics**: matching a variant to a drug response.

**PharmGKB** (Whirl-Carrillo et al. 2021, *Clinical Pharmacology and Therapeutics* 110:563) is the Stanford-maintained free knowledge base for pharmacogenomic variants. The free tier covers the ~50 most clinically relevant pharmacogenes plus the ~30 **CPIC** (Clinical Pharmacogenetics Implementation Consortium) gene-drug pairs that have published recommendations.

The pharmacogenomics workflow is different from the ACMG framework in two important ways:

- **Star alleles, not point variants.** Pharmacogenes like `CYP2D6` and `CYP2C19` have many polymorphic positions, and the clinically relevant unit is the *combination* of variants on one haplotype, named the **star allele** (e.g. `CYP2D6*4`, `CYP2C19*17`). Calling a star allele requires knowing the **phase** of multiple variants — which ones are on the same chromosome copy. Short-read sequencing struggles with phasing across more than a few hundred base pairs; modern pharmacogenomics pipelines either use a star-allele caller (Stargazer, Aldy, Cyrius — all free) or default to "indeterminate" when phasing is ambiguous.
- **The output is a phenotype, not a classification.** Star alleles map to functional categories: **ultrarapid metabolizer**, **normal metabolizer**, **intermediate metabolizer**, **poor metabolizer**, plus a few drug-specific phenotypes. CPIC publishes a recommendation table per gene-drug pair that maps phenotype to dosing recommendation.

### The CPIC tier-1 gene-drug pairs (the ~20 with strong evidence)

| Gene     | Drugs (selected)                                  | Effect of poor metabolizer |
|----------|---------------------------------------------------|----------------------------|
| CYP2D6   | codeine, tramadol, oxycodone, fluoxetine, paroxetine | reduced analgesia / increased side effects |
| CYP2C19  | clopidogrel, omeprazole, citalopram, escitalopram | reduced efficacy / increased side effects |
| CYP2C9   | warfarin, phenytoin                               | bleeding risk / toxicity   |
| TPMT     | 6-mercaptopurine, azathioprine, thioguanine       | severe myelotoxicity       |
| DPYD     | 5-fluorouracil, capecitabine                      | severe fluoropyrimidine toxicity |
| UGT1A1   | irinotecan, atazanavir                            | hyperbilirubinemia / neutropenia |
| SLCO1B1  | simvastatin (and other statins)                   | statin-induced myopathy    |
| VKORC1   | warfarin                                          | dose sensitivity           |
| HLA-B    | abacavir (HLA-B*57:01), allopurinol (HLA-B*58:01) | severe cutaneous adverse reaction |
| HLA-A    | carbamazepine (HLA-A*31:01), abacavir             | severe cutaneous adverse reaction |
| CYP3A5   | tacrolimus                                        | dosing requirement         |
| CFTR     | ivacaftor                                         | drug responsiveness        |
| G6PD     | rasburicase, primaquine, dapsone                  | hemolytic anemia           |
| IFNL3 (IL28B) | peginterferon + ribavirin                    | hepatitis C response rate  |

The ~30 less-well-evidenced pairs add tier-2 and tier-3 genes. PharmGKB indexes all of these and links each to the published CPIC recommendation.

### Querying PharmGKB

```python
import requests


def fetch_pharmgkb_gene(symbol: str, timeout: float = 30.0) -> dict[str, Any]:
    """Fetch the PharmGKB record for a gene symbol."""
    response = requests.get(
        "https://api.pharmgkb.org/v1/data/gene",
        params={"symbol": symbol, "view": "max"},
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()


def fetch_pharmgkb_guideline(guideline_id: str, timeout: float = 30.0) -> dict[str, Any]:
    """Fetch the CPIC guideline record for a given guideline ID."""
    response = requests.get(
        f"https://api.pharmgkb.org/v1/data/guideline/{guideline_id}",
        params={"view": "max"},
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()
```

A typical pharmacogenomics pipeline:

1. Take a VCF restricted to the CPIC tier-1 pharmacogenes.
2. Call star alleles with a star-allele caller (Aldy or Stargazer).
3. Map the star-allele diplotype to the metabolizer phenotype using the CPIC-published lookup table.
4. For each drug in the CPIC tier-1 list, look up the recommendation by phenotype.
5. Emit a per-drug report: drug, phenotype, recommendation, evidence level.

For Week 8 we implement a **simplified** version that does not do star-allele calling (which is out of scope) but does the rest: take a VCF, look up known PharmGKB variants (by rsID), and report the matched gene-drug pairs and recommendations. Challenge 2 implements this.

---

## 7. The limits of automation

The whole point of this lecture is to be honest about what an automated annotation can and cannot do. To summarize:

**An automated pipeline CAN**:

- Annotate the molecular consequence on every overlapping transcript.
- Look up the population frequency in gnomAD.
- Look up the clinical assertion in ClinVar.
- Compute the SIFT and PolyPhen scores.
- Apply the ~8 mechanically computable ACMG criteria.
- Flag conflicting ClinVar interpretations.
- Look up CPIC-graded pharmacogenomic variants.
- Produce a structured, reproducible, version-pinned report.

**An automated pipeline CANNOT**:

- Apply the PS2/PS3/PS4 functional/case-level criteria (requires literature review).
- Apply the PM3/PP1 segregation/phasing criteria (requires pedigree / trio data).
- Apply the PP4 phenotype-specificity criterion (requires patient phenotype).
- Make a "clinical action" recommendation. Even given a Pathogenic classification, the action depends on the patient's clinical context, family history, and consent.
- Call star alleles confidently from short-read data without a specialized caller.
- Resolve the "this looks pathogenic in silico but is in a gene where LOF is tolerated" cases without a curated gene list.

A clinical lab's variant interpretation report goes through manual review by a board-certified clinical geneticist *for every reportable variant*. The automated annotation is the *starting point* for that review — it puts the evidence on one page so the geneticist can read across it — but the classification on the report is the geneticist's, not the pipeline's.

For Week 8's pipeline output, we will always include the disclaimer at the top of the report:

> *This report is for educational and research use only. It is not a clinical interpretation. Variant interpretation in a clinical context requires review by a CLIA-certified laboratory and a board-certified clinical geneticist. Do not act on this report.*

If you find yourself wanting to remove this disclaimer because "the pipeline is good," that is exactly the moment to leave it in.

---

## 8. What to remember

- **28 ACMG/AMP criteria; ~8 mechanically computable**. Be explicit about which you applied and which you skipped.
- **A single BA1 hit (gnomAD popmax > 0.05) is sufficient to call Benign**. PVS1 + PM2 alone is Likely Pathogenic, not Pathogenic.
- **SIFT and PolyPhen-2 are the canonical computational predictors**. Use them together for PP3/BP4; disagreement is a flag for manual review.
- **PP5 and BP6 were deprecated in 2018**. Do not implement them in new pipelines.
- **Pharmacogenomics speaks in star alleles, not point variants**. Calling star alleles requires phasing; short-read data is fragile. Use a dedicated caller or note the limitation.
- **CPIC tier-1 pharmacogenes** are the ~14 genes with strong published recommendations. PharmGKB indexes them and the CPIC guideline tables map phenotype to recommendation.
- **The automated report is not a clinical interpretation**. Every report carries the disclaimer; the clinical work is done by clinicians.

You have now finished the three Week 8 lectures. Continue to the exercises to put the patterns into code.
