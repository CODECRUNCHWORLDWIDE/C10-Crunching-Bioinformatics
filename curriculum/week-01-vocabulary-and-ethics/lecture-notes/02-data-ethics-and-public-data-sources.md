# Lecture 2 — Data Ethics and Public Data Sources

> **Duration:** ~2 hours of reading + a short written reflection.
> **Outcome:** You can articulate why genetic data is uniquely sensitive, explain the two C10 data-ethics rules and *why* they exist, identify at least four public consent-cleared datasets and what each is for, and recognize when language about genetics has crossed into determinism.

This lecture is non-optional. Other weeks teach you how to operate the tools; this one teaches you *when not to point them at someone*. If you skip it you will produce technically competent work that nonetheless violates the norms of the field. Worse: you may violate them on data that is genuinely sensitive to a real person who consented to one use and not another.

If you only remember one thing from this lecture, remember this:

> **Genetic data is not like other data.** It identifies the donor, identifies the donor's relatives who never consented, and cannot be revoked once shared. Treat it accordingly.

---

## 1. Why genetic data is uniquely sensitive

Three properties make genetic data different from, say, a clickstream log.

### 1.1 It identifies the donor — usually permanently

A whole-genome sequence is a unique fingerprint of one person. About 30–80 SNPs (single-nucleotide variants) are enough to uniquely identify any individual against population reference panels. The 23andMe-class genotyping arrays measure hundreds of thousands of SNPs, far more than needed for identification.

"De-identification" — stripping name, date of birth, address — does not protect a genome. The genome *is* the identifier.

### 1.2 It identifies the donor's relatives

Your siblings share ~50% of your genome. Your parents and children share ~50%. Your first cousins share ~12.5%. This means **publishing your genome publishes information about people who did not consent**.

This is not theoretical. The Golden State Killer was identified in 2018 not from his own DNA, which was not in any database, but from a third cousin's genealogy upload to GEDmatch. The cousin had not consented to law-enforcement use; the cousin had not even known what she was consenting to in the genealogy context.

### 1.3 It cannot be revoked

Once a sequence is on the internet, it is on the internet. A donor who consents in 2024 and changes their mind in 2034 has no operational way to recall the data. This makes consent-at-the-time-of-collection ethically heavier than for other data types.

### 1.4 It can be re-identified

The Gymrek et al. 2013 paper in *Science* showed that **men's surnames can often be inferred** from short tandem repeats on the Y chromosome by cross-referencing public genealogy databases. The authors took "anonymous" male genomes from the 1000 Genomes Project and recovered surnames for ~12% of them — and from there, additional public records yielded full identities.

This was a methodological wake-up call. The field's consensus since: **assume re-identification is possible** for sufficiently rich genetic data, and design consent and access accordingly.

---

## 2. The two C10 data-ethics rules

These rules apply to **every** exercise, challenge, homework problem, mini-project, and capstone in C10. They are non-negotiable.

### Rule 1 — Use public, consent-cleared datasets for all coursework

Datasets where donors explicitly consented to public research use are designed for situations exactly like ours. The canonical ones for C10:

- **1000 Genomes Project, phase 3** — 2,504 individuals from 26 populations. Open access. Consent explicitly covers research use including by students.
- **GTEx (public tier)** — tissue-level gene-expression summary data. The individual-level data lives behind dbGaP; we only use the public summary tier in C10.
- **GISAID** (viral) — SARS-CoV-2 and influenza sequences. Registration required, academic use covered, attribution required.
- **NCBI RefSeq / Ensembl reference genomes** — these are *not* data from individuals; they are community-curated reference sequences. Free to use.
- **NCBI SRA datasets explicitly marked open access** — read the consent statement on each. Most "open access" SRA datasets are from cell lines or model organisms, not from individuals with their genomes on display.

### Rule 2 — Never analyze a friend or family member's DNA for the course

Even if they "say it's fine."

This rule sounds restrictive. It exists because:

- **You cannot revoke a consent decision** they make casually today and regret in two years.
- **Their relatives have not consented** — and as Section 1 covered, your friend's genome carries information about their relatives.
- **You may discover medically actionable findings** (a BRCA1 variant, a Huntington's expansion) that you are not equipped to communicate or counsel on.
- **You are not under IRB oversight.** Academic researchers who study family genetics do so under protocols that require ethics-board review precisely because the dynamics are complex.

If you have a genuine research interest in your family's genetics, the path is: do it as a research project at a university, under an IRB-approved protocol, with a faculty advisor. Not as a Week-6 mini-project.

> **A note for adopted students or anyone interested in their own DNA.** Consumer testing (23andMe, AncestryDNA) is a separate decision — you are a single consenting adult. You can use a personal sequence for your own learning *outside* of submitted coursework. But do not commit it to a public repo, do not share it on the course tracker, and recognize that you are also publishing information about biological relatives you have not consulted.

---

## 3. IRB, GDPR, HIPAA — at a high level

We are not lawyers. The point of this section is not to make you a compliance expert — it is to make sure you know *that these regulations exist* and *who to ask* when you encounter them in real work.

### IRB — Institutional Review Board

In the US, any research involving human subjects at a university or hospital must be reviewed and approved by an IRB. The IRB checks consent, data handling, risk-benefit balance, and the right to withdraw. Bioinformatics work using consent-cleared *public* data is usually exempt; bioinformatics work using individual-level data from an institution is usually not.

Outside the US, the analog goes by different names: Research Ethics Committee (UK, Canada, Australia), Ethics Commission (continental Europe). Same role.

### HIPAA — Health Insurance Portability and Accountability Act (US)

A US federal law that governs **protected health information (PHI)** held by covered entities (hospitals, insurers, their business associates). Genetic data from a clinical context is PHI. The Privacy Rule lists 18 identifiers that must be removed for "safe harbor" de-identification — but genetic data is famously *not* one of them, which is part of why de-identification can be insufficient.

If you ever work at a US healthcare institution and someone asks "is this dataset HIPAA-covered?" — that is a question for the institution's privacy office, not for you.

### GDPR — General Data Protection Regulation (EU)

The EU's data-protection law. Treats **genetic data as a special category** of personal data requiring stronger protections than ordinary personal data. Notable principles:

- **Lawful basis** must be explicit (consent is one option; there are others).
- **Right to erasure** — though in practice, for genetic data on the internet, this is not always operationally possible.
- **Data minimization** — collect only what's needed.

If you work in Europe, or with European subjects, GDPR applies regardless of where your servers are.

### The actionable advice

You do not need to memorize any of this. You need to:

1. **Recognize the names** so you know what people are talking about in a code review.
2. **Default to public consent-cleared data** for everything in C10 — which sidesteps almost all of these regimes.
3. **Know who to ask** when you encounter the real thing: your institution's IRB / ethics committee / privacy office.

---

## 4. The public-data ecosystem

Most bioinformatics datasets you will ever use live in one of a small number of public archives. Know what each is for.

### NCBI — National Center for Biotechnology Information (US)

The largest public bioinformatics repository. Subdivisions you will use:

- **GenBank** — sequence records, donor-submitted.
- **RefSeq** — community-curated reference sequences for genes and genomes.
- **SRA — Sequence Read Archive** — raw sequencer output from published studies. Petabytes.
- **dbSNP** — known human variants.
- **ClinVar** — variants with curated clinical interpretation.
- **PubMed** — the biomedical literature index.

Free, no account required for most reads. You will use NCBI from Week 2 onward.

### Ensembl (EMBL-EBI, UK)

The European counterpart to NCBI for genome annotation. Hosts reference genomes and rich annotation for many species, plus the **Ensembl REST API** which we hit programmatically in Week 4. Also hosts the **VEP** (Variant Effect Predictor) which we use in Week 6.

### EBI — European Bioinformatics Institute

Parent organization of Ensembl, also hosts:

- **ENA — European Nucleotide Archive** — the European mirror/peer of NCBI's GenBank/SRA.
- **UniProt** — protein sequences and annotations. Used heavily in Week 4.
- **EMBL-EBI Training** — free courses.

### 1000 Genomes Project

A finished project (2008–2015) that sequenced 2,504 individuals from 26 populations at low-to-moderate coverage. The phase-3 release is the workhorse dataset for population genetics teaching. Fully open, donors consented to research and education use, distributed under the Fort Lauderdale Principles. Hosted at IGSR (International Genome Sample Resource).

### GTEx — Genotype-Tissue Expression

A US National Institutes of Health project that profiled gene expression across ~50 tissues in nearly 1,000 deceased donors. **Two access tiers**: a *public summary tier* (median expression per gene per tissue, eQTL summary statistics) and a *controlled-access tier* (individual-level data, dbGaP application required). **In C10 we use only the public tier.**

### TCGA — The Cancer Genome Atlas

Cancer genomics across 33 tumor types. Controlled access for individual-level data. Outside C10's scope, but you should know it exists; many of your future colleagues will have used it.

### GISAID

The viral-sequence archive that hosted most SARS-CoV-2 sequences during the COVID-19 pandemic. Free registration, academic terms of use, mandatory attribution. We use GISAID in Week 9 for the phylogenetics module.

### UCSC Genome Browser

A genome-visualization site at UC Santa Cruz that has been running since 2000. Hosts a different version of the same reference genomes plus thousands of annotation "tracks." Useful for sanity-checking ("does this gene really start where I think it starts?"). We touch it in Week 10.

---

## 5. The Fort Lauderdale and Bermuda Principles — why bioinformatics is open

Two foundational meetings shaped the field's open-data culture:

- **Bermuda Principles (1996)** — established at a Wellcome Trust meeting during the Human Genome Project. Result: human-genome sequence data goes into public databases **within 24 hours of generation**, before publication. Radical at the time.
- **Fort Lauderdale Principles (2003)** — extended to community-resource projects more broadly. Result: data generators get a publication window; data users get pre-publication access but must respect a citation embargo.

These principles are why a 12-week open-source bioinformatics curriculum is possible. Most of the data we will use was *deliberately* released under permissive terms by the consortia that produced it.

---

## 6. Citation and reproducibility — the operational ethics

The flip side of using free public data: **cite it properly**.

A real methods-section excerpt looks like this:

```
RNA-seq reads were aligned to the GRCh38.p14 human reference genome
using STAR v2.7.11a [Dobin et al. 2013] with default parameters.
Read counts per gene were generated using featureCounts v2.0.6
[Liao et al. 2014] against GENCODE v45 annotation [Frankish et al. 2023].
Differential expression analysis was performed with DESeq2 v1.42.0
[Love et al. 2014] in R 4.3.2.
```

Notice:

- **Every tool has a version.**
- **Every dataset has a version or accession.**
- **Every method has a citation.**

This is not optional in a real paper. It is also the practice we adopt in C10 from Week 1. The [`reproducibility receipt`](../../../assets/branding/BRAND.md) box in every mini-project enforces this — five rows, always, including wall time.

---

## 7. Language — the things you do not say

The voice of C10 (see `assets/branding/BRAND.md`) is lab-notebook precise. A specific subset of that voice is about **how we talk about genetics**.

### Things to avoid

| Don't say | Why |
|-----------|-----|
| "Gene X *causes* disease Y" | Almost no gene single-handedly causes a complex disease. Even Mendelian variants have variable penetrance. |
| "Gene for intelligence" / "gene for X behaviour" | The vast majority of behavioural traits are polygenic and gene-environment-interactive. This framing is essentially always wrong. |
| "His DNA *made* him do it" | Genetic determinism. Modern behavioral genetics does not support this for any complex behaviour. |
| "The genome decides" | Same problem. The genome *constrains and influences*; environment, development, and chance also matter. |
| "Race is genetic" | Self-identified race is a poor predictor of any specific genetic variant. Ancestry (continental genetic similarity) is a population-genetics concept and is not synonymous with social race. |
| "Anonymized genome" | See Section 1.4. Use "de-identified" and assume re-identification is possible. |

### Things to say instead

| Better | Why |
|--------|-----|
| "Gene X is associated with disease Y" | Honest about the statistical nature of the claim. |
| "Variant rs6265 has been reported to influence …" | Cites the level of evidence. |
| "Suggests a role in …" | Acknowledges uncertainty. |
| "Phenotypic variation in this trait is partly heritable, with an estimated heritability of h² ≈ 0.4" | Gives the actual statistical claim, with a number. |
| "De-identified" | Honest about what was done; does not overclaim safety. |

Reading good methods sections (the homework's Problem 1) will tune your ear. The single best self-check: re-read your own writing and ask *"would I be comfortable defending this sentence in a methods section?"* If not, rewrite.

---

## 8. The open-science norm

Bioinformatics is structurally open-source-friendly. Most major journals in the field (*Bioinformatics*, *Nature Methods*, *Genome Research*, *PLOS Comp Bio*) require:

- Code deposited in a public repository (GitHub, Zenodo).
- Data deposited in a public archive (GenBank, ENA, SRA, GEO) **before** publication.
- Container/environment files (Conda, Docker) when feasible.

This is *the field's social contract*. Your portfolio repo at the end of C10 — `crunch-bio-portfolio-<yourhandle>` — is your participation in it. Build it in public from Week 1 onward.

---

## 9. A short ethical scenario — work through it

A friend, a biology Ph.D. student, casually offers you a small VCF file from their lab's project. Two samples. They say "it's de-identified, just play with it for practice." What do you do?

The honest answer:

1. **Ask what consent the donors gave.** If "research use including by students at our institution," the lab's IRB should be able to confirm. If unsure, the answer is no.
2. **Ask whether the data is already public.** If it is — use the public version with its DOI, not the file your friend handed you. Public data has a clean provenance trail.
3. **Default to declining** if either of the above is uncertain. The cost of saying no is small (you find a public alternative — there are dozens). The cost of saying yes incorrectly is real.

This is the Week-1-of-a-bioinformatics-career version of "don't use production data in dev." Once you internalize it, it's easy.

---

## 10. Self-check

Before moving on, you should be able to answer all of these out loud, without notes:

1. Name three reasons genetic data is uniquely sensitive.
2. What are the two C10 data-ethics rules?
3. What is the difference between IRB, HIPAA, and GDPR?
4. Name four public datasets and what each is for.
5. What is the Fort Lauderdale Principle, in one sentence?
6. Re-read a sentence from a recent news article about genetics and identify one phrase that crosses into determinism.
7. What three things should a real methods sentence include?

If any felt shaky, re-read the relevant section before starting the exercises.

---

*Return to [Week 1 README](../README.md) and start the [exercises](../exercises/README.md).*
