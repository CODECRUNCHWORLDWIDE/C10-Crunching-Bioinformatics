# Lecture 1 — The Three Axes of Variant Interpretation

> **Educational and research use only.** This lecture describes the *mechanics* of variant annotation. It does not describe how to make a clinical interpretation. Variant interpretation in a patient-care context requires a CLIA-certified laboratory, board-certified clinical geneticists, and the full ACMG/AMP workflow with manual review. The pipelines you build this week are the same mechanical steps a clinical pipeline performs, but they are not those pipelines.

> **Duration:** ~3 hours of reading + a brief cyvcf2 sanity check.
> **Outcome:** You can describe the three orthogonal axes of variant interpretation, name the canonical free tool and database for each, run VEP and SnpEff on a small VCF, and read the resulting CSQ / ANN INFO fields.

If you only remember one thing from this lecture, remember this:

> **A variant annotation has three orthogonal axes — what the variant DOES (functional consequence: VEP, SnpEff), how OFTEN it is seen in healthy populations (population frequency: gnomAD), and what is KNOWN about it in the literature (clinical knowledge: ClinVar, dbSNP). All three axes are computed by mechanical lookup against free public databases. None of the three is a clinical interpretation. A clinical interpretation requires a clinician, a patient phenotype, a pedigree, and the ACMG/AMP framework, and it lives outside the scope of this week.**

Week 6's VCF is your input. Week 7's counts matrix is parked. Week 8's output is a per-variant report.

---

## 1. Where we are in the pipeline

The end-to-end DNA-sequencing pipeline, recapping Weeks 5 and 6, runs:

```
FASTQ -> trim (fastp) -> align (bwa mem) -> sort (samtools sort) -> markdup (samtools markdup) ->
        BAM -> call variants (bcftools / GATK HaplotypeCaller) -> VCF (raw) ->
        hard-filter (bcftools filter / GATK VariantFiltration) -> VCF (PASS subset) ->
        annotate (VEP, SnpEff) -> VCF (annotated) ->
        cross-join with ClinVar + gnomAD + dbSNP -> CSV / HTML report
```

Week 6 left off at the PASS subset. Week 8 picks up at "annotate" and finishes at "report."

The mechanics from Week 6 are unchanged. The VCF spec is the same. The PASS filter is the same. The new work is **what to do with each PASS variant** — and the answer is, mechanically: annotate, look up, and report.

---

## 2. The three axes

Every variant has three independent questions associated with it. They are independent in the sense that a variant can be common (frequent in healthy populations) AND also annotated as missense AND also recorded in ClinVar as Benign; or rare AND nonsense AND Pathogenic; or common AND nonsense AND Benign (~1% of nonsense variants are tolerated). The interpretation comes from combining all three answers, not from any single axis.

### Axis 1 — Functional consequence: what does the variant DO?

The canonical free tools: **VEP** (McLaren et al. 2016, *Genome Biology* 17:122) and **SnpEff** (Cingolani et al. 2012, *Fly* 6:80). Both walk the variant against a gene model (Ensembl, RefSeq, or GENCODE) and emit a Sequence Ontology consequence term for every overlapping transcript:

- `synonymous_variant` — same amino acid (silent at the protein level).
- `missense_variant` — different amino acid (single-residue substitution).
- `stop_gained` (also called "nonsense") — premature stop codon, truncates the protein.
- `stop_lost` — readthrough into the 3' UTR, extends the protein.
- `frameshift_variant` — indel of length not divisible by 3 in a CDS, shifts the reading frame.
- `inframe_insertion` / `inframe_deletion` — indel of length divisible by 3 in a CDS, inserts or deletes whole codons.
- `splice_donor_variant` / `splice_acceptor_variant` — variant at the canonical GT-AG splice site dinucleotide (+/- 1, 2 bp from the exon-intron boundary).
- `splice_region_variant` — variant in the 1-3 bp window around the splice site, less mechanistic than donor/acceptor.
- `intron_variant`, `5_prime_UTR_variant`, `3_prime_UTR_variant`, `intergenic_variant` — non-coding regions.

VEP and SnpEff both also report an **impact tier** that buckets the consequence terms:

- **HIGH** — frameshift, stop_gained, stop_lost, splice_donor, splice_acceptor, start_lost, transcript_ablation. These are the "null variants" of the ACMG PVS1 criterion.
- **MODERATE** — missense, inframe_indel, protein_altering_variant. Same length but altered protein sequence.
- **LOW** — synonymous, splice_region, stop_retained. Coding but expected mild effect.
- **MODIFIER** — intron, UTR, intergenic, downstream, upstream. Non-coding.

This impact tier is the most useful single field in the annotation, because most downstream filtering starts with "drop MODIFIER variants" (~80% of any whole-genome VCF) and works from there.

A subtle point: **the consequence depends on the transcript**. A single genomic variant overlaps multiple transcripts in most multi-isoform genes. VEP and SnpEff emit one consequence record per transcript, separated by commas in the CSQ/ANN INFO field. The variant might be missense in transcript A, synonymous in transcript B, and intronic in transcript C. The canonical transcript (Ensembl `--canonical`) or MANE Select (Morales et al. 2022, *Nature* 604:310) is usually what is reported as "the" consequence, but you should be aware that the "canonical" choice is a convention, not a biological fact.

### Axis 2 — Population frequency: how OFTEN is the variant seen in healthy people?

The canonical free database: **gnomAD** (Karczewski et al. 2020, *Nature* 581:434), the Broad-Institute-maintained aggregate of ~807,000 samples (combining v4 exomes and genomes) from biobanks and case studies, with severe pediatric disease cases excluded. gnomAD reports per-variant allele frequency (AF) overall and stratified by sub-population: African (afr), Admixed American / Latino (amr), Ashkenazi Jewish (asj), East Asian (eas), Finnish (fin), Non-Finnish European (nfe), South Asian (sas), Middle Eastern (mid, v4 only), and Other (oth).

The "popmax" AF — the maximum AF across sub-populations — is what most clinical pipelines pivot on. A variant with popmax > 0.05 is *common* and almost certainly not the cause of a rare Mendelian disease (the ACMG BA1 stand-alone benign criterion). A variant with popmax < 0.0001 is *rare* and supports a pathogenic interpretation (the ACMG PM2 moderate criterion).

Sub-population stratification matters. The CFTR `c.3454G>C` variant is at AF 0.07 in Finnish (common) but AF 0.0003 globally (rare). For a Finnish patient, this variant is essentially a Finnish polymorphism; for a Japanese patient, it is rare. The popmax field captures the higher of the two.

gnomAD also flags the per-variant **filter** status (`PASS`, `RF`, `AC0`, `InbreedingCoeff`, etc.) — a variant in gnomAD with a non-PASS filter has been called but flagged as low-quality, and should be discounted in the frequency estimate.

A critical caveat: **gnomAD does not contain pediatric disease cohorts**. Severe Mendelian disease alleles are *underrepresented* in gnomAD relative to the general population. This is intentional — the goal is to estimate "frequency in the apparently healthy population" — but it means that a Pathogenic variant for a recessive disease may show AF 0.001 in gnomAD (carrier frequency) even though the disease prevalence implies a much higher allele frequency. The ACMG framework accounts for this with PM2 ("absent or rare") rather than "absent" alone.

### Axis 3 — Clinical knowledge: what is KNOWN about the variant?

The canonical free databases: **ClinVar** (Landrum et al. 2018, *Nucleic Acids Research* 46:D1062), the NCBI clinical-variant database; and **dbSNP** (Sherry et al. 2001), the NCBI single-nucleotide-polymorphism database.

ClinVar records have, for each variant submission:

- **CLNSIG** (clinical significance): Benign / Likely_benign / Uncertain_significance (VUS) / Likely_pathogenic / Pathogenic / drug_response / risk_factor / association / protective / other / Conflicting_interpretations_of_pathogenicity.
- **CLNREVSTAT** (review status): a 0-to-4-star scale. 4 stars (`practice_guideline`) is the highest evidence; 3 stars (`reviewed_by_expert_panel`) is from organizations like ENIGMA (BRCA1/2) or ClinGen; 2 stars (`criteria_provided,_multiple_submitters,_no_conflicts`) is multiple labs agreeing; 1 star is a single lab assertion; 0 stars is no assertion criteria.
- **CLNDN** / **CLNDISDB**: associated disease name and ontology IDs.
- **MC**: molecular consequence (mostly redundant with VEP annotation, useful for cross-checking).

A variant with a 4-star Pathogenic assertion in ClinVar is the strongest individual piece of evidence you can get from any single source. A variant with a 1-star Pathogenic assertion is informative but not authoritative. A variant with conflicting submissions is a flag for manual review.

dbSNP is older and broader. Every common variant has a dbSNP rsID (e.g. `rs121908755`). The rsID is what the literature uses to refer to a variant; if you find a paper saying "patients with rs121908755 had increased risk of X," dbSNP is where you confirm the variant coordinates and ALT allele. dbSNP itself does not assert pathogenicity; it is a coordinate-to-rsID-to-literature index.

ClinVar and dbSNP are cross-linked: most ClinVar records carry a dbSNP RS field, and vice versa.

### The three axes are independent

A variant can be:

- **Common, missense, ClinVar Benign**: e.g. a known polymorphism. Action: drop from clinical consideration.
- **Common, synonymous, ClinVar absent**: a synonymous polymorphism. Action: drop.
- **Rare, missense, ClinVar absent, in a gene relevant to the phenotype**: a VUS by ACMG criteria. Action: flag, defer, manual review.
- **Rare, frameshift, ClinVar Pathogenic 3-star**: a likely actionable pathogenic finding. Action: report to the clinical team, *not* to the patient — the clinical team confirms in a CLIA lab and decides.
- **Common (popmax > 0.05), frameshift, ClinVar Pathogenic**: this is the conflicting case. The high frequency contradicts the pathogenic call, usually because the ClinVar submission was based on outdated or weak evidence. Action: flag the conflict, do not act on either assertion alone.

The job of the per-variant report is to put all three axes on one row and let a human read across them.

---

## 3. VEP (Variant Effect Predictor)

VEP is the Ensembl-maintained annotator. The "canonical" of the two free tools.

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

What each flag does:

- `-i / -o`: input/output VCFs. VEP can also emit JSON (`--json`) or tab-separated (`--tab`).
- `--cache --dir_cache --offline`: use a local cache, no network. The cache is ~25 GB for human GRCh38 and is the static part of the annotation; once downloaded, every annotation run is offline.
- `--species --assembly --cache_version`: pin the species, assembly, and cache release. The cache release matches an Ensembl release (110, 111, 112, ...) and drift between releases is ~1-3% per major release.
- `--sift b --polyphen b`: include SIFT and PolyPhen-2 in the output. `b` means "both score and prediction" (e.g. `SIFT=0.02 (deleterious)`).
- `--canonical`: flag the Ensembl canonical transcript.
- `--mane`: flag the MANE Select transcript (the joint Ensembl-NCBI canonical transcript, available since 2022).
- `--symbol --biotype`: include the gene symbol (e.g. `BRCA1`) and the transcript biotype (`protein_coding`, `processed_transcript`, etc.).
- `--af --af_gnomadg --af_gnomade`: include allele frequencies from 1000 Genomes, gnomAD genomes, gnomAD exomes.
- `--check_existing`: include the existing rsID.
- `--fork 4`: 4 worker processes.

The output is an extended VCF with a `CSQ` INFO field on each variant. The `CSQ` field is pipe-separated and contains one record per overlapping transcript, with the fields documented in the VCF header:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|SIFT|PolyPhen|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADg_AF">
```

A real example line:

```
chr7  117559590  .  G  A  500.0  PASS  CSQ=A|missense_variant|MODERATE|CFTR|ENSG00000001626|Transcript|ENST00000003084|protein_coding|11/27||ENST00000003084.11:c.1521_1523delCTT|ENSP00000003084.6:p.Phe508del|1521-1523|1521-1523|508|F/-|TTC/---|rs113993960&CM900037|||1|YES|ENST00000003084.11|deleterious(0)|probably_damaging(1)|0.007|...
```

That is one transcript record for the famous CFTR `p.Phe508del` variant: it is a missense (technically inframe deletion; VEP sometimes reports `inframe_deletion` instead — version-dependent) with IMPACT MODERATE on the canonical CFTR transcript, with SIFT `deleterious(0)` and PolyPhen `probably_damaging(1)`, at gnomAD AF 0.007.

To parse the CSQ field in Python, split on `,` (one per transcript), then split each on `|` (one per field), and map to the field names from the header. We do this with cyvcf2 in Lecture 2.

---

## 4. SnpEff

SnpEff is the independent annotator. The "second opinion" of the two free tools.

```bash
snpEff \
    -v \
    -csvStats stats.snpeff.csv \
    -htmlStats stats.snpeff.html \
    GRCh38.105 \
    input.vcf \
    > output.snpeff.vcf
```

SnpEff emits an `ANN` INFO field with the format documented in its header:

```
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
```

Format is very similar to VEP's CSQ but the field order and the impact strings differ slightly. SnpEff's impact column uses `HIGH` / `MODERATE` / `LOW` / `MODIFIER`, mapping closely (but not identically) to VEP's.

A real example for the same CFTR `p.Phe508del`:

```
ANN=A|missense_variant|MODERATE|CFTR|CFTR|transcript|NM_000492.4|protein_coding|11/27|c.1521_1523del|p.Phe508del|1521/4443|1521/4443|508/1480||
```

The two annotations agree on the consequence (missense or inframe_deletion, depending on tool/version), the impact (MODERATE), the gene (CFTR), and the protein change (Phe508del). They disagree on the canonical transcript ID, because VEP uses Ensembl IDs (`ENST00000003084`) and SnpEff uses RefSeq IDs (`NM_000492.4`) by default. Both are correct; they are referring to the same biological entity.

**Best practice**: run both. Where they agree, the annotation is robust. Where they disagree, you have a list of variants to inspect manually. ~95% of variants will agree on impact tier; ~5% will not, and that 5% is mostly splice-region, intron-near-splice, and the missense-vs-inframe-deletion edge.

---

## 5. The consequence taxonomy in pictures

```
Protein-coding region:

  exon 1                     exon 2                        exon 3
  ATG GTA CGA  ...  AGT GAT  |  intron  |  TAC GGA  ...  TAA
  M   V   R         S   D   |          |  Y   G         *
                                                          ^ stop

Some example variants and their consequence calls:

  Position   Change      Codon change      Protein     VEP consequence
  ----       ATG->ACG    M -> T            start lost  start_lost (HIGH)
  ----       GTA->GTC    V (silent)        V           synonymous_variant (LOW)
  ----       CGA->CGT    R (silent)        R           synonymous_variant (LOW)
  ----       AGT->AAT    S -> N            N           missense_variant (MODERATE)
  ----       GAT->TAT    D -> Y            Y           missense_variant (MODERATE)
  ----       AGT->TGA    S -> *            *           stop_gained (HIGH)
  ----       insert T at "AGT" position    frameshift  frameshift_variant (HIGH)
  ----       intron GT->AT (donor)         splice loss splice_donor_variant (HIGH)
  ----       intron 5bp into the intron    splice region splice_region_variant (LOW)
  ----       TAA->CAA    * -> Q            stop lost   stop_lost (HIGH)
```

The consequence types are defined by the Sequence Ontology and shared across VEP, SnpEff, ANNOVAR, and most other annotators. The full SO consequence vocabulary has ~50 terms; the ~15 terms above cover ~99% of variants in a typical exome.

---

## 6. The impact tier as a filter

The first thing most clinical pipelines do is drop the MODIFIER tier:

| Tier      | Variants in a typical exome | Fraction |
|-----------|-----------------------------|----------|
| HIGH      | ~30-100                    | < 0.1%   |
| MODERATE  | ~9,000-12,000              | ~20%     |
| LOW       | ~12,000-15,000             | ~25%     |
| MODIFIER  | ~25,000-30,000             | ~55%     |

For a whole exome (~50,000 variants), the MODIFIER bucket is the majority by count. Most are UTR / intron / upstream / downstream variants that are very unlikely to be clinically actionable, *with the caveat* that some non-coding variants are pathogenic (e.g. promoter / enhancer / branch-point variants in `BRCA1`, `MLH1`, others). Modern clinical pipelines do not drop the MODIFIER tier outright; they down-weight it and re-examine if the HIGH and MODERATE buckets produce no plausible findings.

For Week 8's small didactic VCFs, we keep all tiers because we want to see the annotation behave on every category.

---

## 7. The "transcript-of-reference" problem

A single genomic variant in `BRCA1` overlaps 21 annotated Ensembl transcripts, of which 11 are protein-coding and 10 are processed transcripts or retained-intron isoforms. VEP emits 21 records in the CSQ field. The consequence may differ across the 11 protein-coding transcripts because of alternative exon usage, alternative start codons, or readthrough.

Which transcript do you report?

The historical answer was "the Ensembl canonical transcript," which is the longest protein-coding transcript by default. The new answer (since 2022) is the **MANE Select** transcript (Morales et al. 2022, *Nature* 604:310): a jointly Ensembl-and-NCBI-curated single canonical transcript per gene, chosen to be the best-supported and most-conserved coding transcript. MANE Select is identical between Ensembl and NCBI annotations and is now the recommended reference for clinical reporting.

For 19,000+ human protein-coding genes, MANE Select is defined for >95% as of 2024. The remaining few percent are genes where Ensembl and NCBI annotations diverge.

VEP's `--mane` flag adds a `MANE_SELECT` field to the CSQ output; SnpEff does not yet (5.2) directly report MANE Select but its canonical transcript usually matches.

When reporting variant interpretations, use MANE Select where available; fall back to Ensembl canonical where MANE Select is undefined.

---

## 8. A Python sanity check: parsing the CSQ field with cyvcf2

cyvcf2 (Pedersen and Quinlan 2017) gives you fast Python access to a VCF. Here is the minimal pattern for extracting the CSQ field and breaking it into per-transcript records:

```python
from __future__ import annotations
from pathlib import Path
from typing import Any

from cyvcf2 import VCF


def parse_csq_field(csq_field: str, csq_header: list[str]) -> list[dict[str, str]]:
    """Split a VEP CSQ INFO string into a list of per-transcript dicts.

    Args:
        csq_field:  the raw CSQ value (e.g. "A|missense_variant|MODERATE|...|,A|intron_variant|MODIFIER|...|")
                    Multiple records are separated by commas; fields within a record by pipes.
        csq_header: ordered list of field names from the VCF header, parsed from
                    the CSQ INFO description line.

    Returns:
        list of dicts, one per overlapping transcript, mapping field name to value.
        Empty strings are preserved (VEP often leaves fields empty when irrelevant).
    """
    if not csq_field:
        return []
    records: list[dict[str, str]] = []
    for raw in csq_field.split(","):
        fields: list[str] = raw.split("|")
        if len(fields) != len(csq_header):
            # Malformed; skip with a warning. This sometimes happens at the
            # boundary of a multi-allelic site if VEP did not normalize first.
            continue
        records.append(dict(zip(csq_header, fields)))
    return records


def extract_csq_header(vcf_path: Path) -> list[str]:
    """Return the ordered CSQ field names from a VEP-annotated VCF header.

    Reads the ##INFO=<ID=CSQ,...> header line, locates the Description
    that starts with "Format: ", and splits the rest on "|".
    """
    vcf = VCF(str(vcf_path))
    for line in vcf.raw_header.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            # Description="...Format: A|B|C|..."
            marker = "Format: "
            idx = line.find(marker)
            if idx == -1:
                continue
            fmt = line[idx + len(marker):].rstrip('">')
            return fmt.split("|")
    raise ValueError(f"No CSQ INFO field in {vcf_path}")


def summarize_variant(variant: Any, csq_header: list[str]) -> dict[str, Any]:
    """Pull out a per-variant summary from a cyvcf2 Variant.

    Picks the canonical (or first) CSQ record as the headline annotation.
    """
    csq_field: str = variant.INFO.get("CSQ", "")
    records = parse_csq_field(csq_field, csq_header)
    canonical = next((r for r in records if r.get("CANONICAL") == "YES"), None)
    if canonical is None and records:
        canonical = records[0]
    if canonical is None:
        canonical = {}
    return {
        "chrom": variant.CHROM,
        "pos": int(variant.POS),
        "ref": variant.REF,
        "alt": ",".join(variant.ALT),
        "rsid": canonical.get("Existing_variation", ""),
        "gene": canonical.get("SYMBOL", ""),
        "consequence": canonical.get("Consequence", ""),
        "impact": canonical.get("IMPACT", ""),
        "hgvsc": canonical.get("HGVSc", ""),
        "hgvsp": canonical.get("HGVSp", ""),
        "sift": canonical.get("SIFT", ""),
        "polyphen": canonical.get("PolyPhen", ""),
        "gnomade_af": canonical.get("gnomADe_AF", ""),
        "gnomadg_af": canonical.get("gnomADg_AF", ""),
    }
```

Run this on a 50-variant demo VCF and you have a 13-column table of "the basic facts" per variant. The next step (Lecture 2) is to enrich each row with ClinVar's clinical assertion and gnomAD's full sub-population frequency, neither of which is in the VEP cache by default.

---

## 9. What to remember

- **Three orthogonal axes**: functional consequence, population frequency, clinical knowledge. Always report all three.
- **VEP and SnpEff are both free, both canonical, run both** for any annotation work that will be reviewed by anyone other than yourself. Disagreements between the two are signals for manual review.
- **Pin everything**: VEP cache version, SnpEff database version, ClinVar release date, gnomAD version, reference assembly. Drift is real and unfixable after the fact.
- **The consequence depends on the transcript**. Use MANE Select where defined; fall back to Ensembl canonical otherwise.
- **The impact tier (HIGH / MODERATE / LOW / MODIFIER) is the most useful single field for filtering**. Most clinical pipelines start by dropping MODIFIER and working from there.
- **Educational use only**. Every report you generate this week must include the disclaimer in the header. Annotation is mechanical; interpretation requires a clinician.

Continue to **Lecture 2 — Population and Clinical Databases** to see how to add gnomAD frequencies and ClinVar assertions to each row of your annotation, and how the rate-limiting on the public APIs shapes the implementation.
