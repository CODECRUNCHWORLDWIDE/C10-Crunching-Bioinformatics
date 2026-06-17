# Mini-Project — End-to-End Somatic Variant Calling Pipeline with Mutational Signature Decomposition

> **Educational and research use only.** This mini-project produces a research-grade somatic variant call set and mutational-signature decomposition from a simulated tumor-normal BAM pair. None of the outputs are clinically actionable. **Real treatment decisions require an accredited clinical lab operating under CAP / CLIA or equivalent oversight; this pipeline does not have the validation, version locking, or quality oversight needed to support patient care. Use only on simulated or anonymized public data.**

> **Reproducibility note.** Somatic pipelines are reproducible only if you ship the inputs (BAMs, reference, PON, germline resource), the parameters (FilterMutectCalls thresholds, contamination table, COSMIC catalog version), the tool versions, and the seed alongside the VCFs. Without the `run-info.json` the call set is an opinion; with it, the call set is a reproducible result. **The outputs of this pipeline must travel with their run-info, every time.**

Build a reproducible somatic variant-calling pipeline that takes a tumor-normal BAM pair (the bundled `data/tumor_chr22.bam` and `data/normal_chr22.bam` plus the matched references), runs **Mutect2** in tumor-normal mode, runs **GetPileupSummaries** and **CalculateContamination** on both samples, runs **FilterMutectCalls** with `--contamination-table`, parses the PASS VCF, computes the 96-class trinucleotide spectrum, decomposes against the **COSMIC v3.3 SBS catalog** using **SigProfilerAssignment**, optionally cross-checks against **Strelka2** (if installed), renders a Markdown report with the FILTER tally and the top signatures, and emits a `run-info.json` recording every parameter.

End with: the unfiltered VCF, the contamination table, the filtered VCF (with FILTER tally), the 96-class spectrum TSV, the SigProfilerAssignment activities and stats, the Markdown report, optional Strelka2 PASS-PASS overlap, and a 700-1000 word write-up that defends every parameter and names every limit.

This is the C10 mini-project that produces a **methods-section-quality somatic variant call with measured provenance**, not just a one-shot demonstration. By the end of it you will have a `cancer_pipeline.py` script and a `run.sh` wrapper, a results directory with the VCFs and the signature decomposition, and a write-up that defends every parameter and explicitly names what the pipeline cannot reliably claim.

**Estimated time:** 7 hours (split across Wednesday, Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-11/mini-project/` directory:

```text
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 11 section)
└── week-11/
    └── mini-project/
        ├── README.md               one-page report (~700-1000 words)
        ├── run.sh                  one-command reproduction script
        ├── env.yml                 conda env file pinning all tool versions
        ├── data/
        │   ├── tumor_chr22.bam
        │   ├── normal_chr22.bam
        │   ├── chr22_GRCh38.fasta + .fai + .dict
        │   ├── chr22_pon.vcf.gz + .tbi
        │   ├── chr22_gnomad.vcf.gz + .tbi
        │   ├── chr22_common_biallelic.vcf.gz + .tbi
        │   └── chr22_intervals.bed
        ├── cancer_pipeline.py      the orchestration script
        ├── starter.py              skeleton implementation with TODOs
        └── results/
            ├── unfiltered.vcf.gz                  Mutect2 raw output
            ├── tumor.pileups.table
            ├── normal.pileups.table
            ├── contamination.table
            ├── filtered.vcf.gz                    FilterMutectCalls output
            ├── filter_tally.md                    PASS / filter breakdown
            ├── pass.vcf.gz                        PASS-only subset
            ├── spectrum_96.tsv                    96-class spectrum
            ├── sigprofiler_out/                   SigProfilerAssignment dir
            ├── signature_summary.md               top signatures + cosine
            ├── strelka_out/ (optional)            Strelka2 working dir
            ├── cross_check.md (optional)          Mutect2 vs Strelka2
            ├── qc_report.md                       combined Markdown report
            └── run-info.json                      run provenance
```

By the end you will have a clean, reproducible Week 11 directory you can point a recruiter at — and `cancer_pipeline.py` is the kind of pipeline that opens conversations with working bioinformaticians and biotech / academic cancer-genomics shops, *as long as* you can speak to its limits.

---

## The dataset

You will work with the bundled `data/tumor_chr22.bam` + `data/normal_chr22.bam` pair — a 1.4 Mb chr22 subset of a public GRCh38 BAM pair with simulated somatic variants seeded at known positions. The corresponding `data/truth_chr22.vcf.gz` lists the seeded variants and can be used to compute precision and recall after FilterMutectCalls.

The choice of "simulated with known truth" rather than "real TCGA tumor" is for three reasons:

- **Privacy.** Real tumor BAMs require institutional credentialing; simulated data is freely distributable.
- **Reproducibility.** Real sequencer data has a non-reproducible random-error component; pinning the simulator seed gives byte-identical inputs across runs.
- **Pedagogical clarity.** The seeded variants are at known positions; you can verify your pipeline recovered them and explain any misses.

Expected behavior of the pipeline on this dataset:

- Mutect2 should produce ~200-300 candidate variants in 3-5 minutes.
- FilterMutectCalls should PASS 60-75% of them; the rest will be `germline`, `weak_evidence`, `panel_of_normals`, and `clustered_events`.
- The contamination estimate should be 0-2% (the simulator does not introduce contamination).
- The 96-class spectrum will be small (~100-150 SNVs); this is at the low end for stable signature decomposition. The write-up should name this limitation.
- The signature decomposition should fit with cosine >= 0.85; SBS1 and SBS5 (clock-like signatures) will dominate the synthetic spectrum.

### Optional real-data extension

If you want a more interesting signature profile, extend the pipeline to a real public dataset:

- **TCGA tumor**: select a small tumor-normal pair from <https://portal.gdc.cancer.gov/>. You will need a dbGaP authorization for the controlled-access whole-exome BAMs; consider the open-access tier first.
- **PCAWG sample**: the ICGC's free tier offers a few public-access samples; <https://dcc.icgc.org/> for the index.
- **Synthetic with stronger signature**: re-simulate the tumor with a tobacco-signature mutation burden to verify that SBS4 dominates.

---

## Rules

- **You may** use GATK 4.5.0.0 (Mutect2, FilterMutectCalls, GetPileupSummaries, CalculateContamination, optionally Funcotator), Strelka 2.9.10 (optional), samtools 1.20, bcftools 1.20, pysam 0.22.1, SigProfilerAssignment 0.1.4, SigProfilerMatrixGenerator 1.2.26, biopython 1.84, pandas 2.2.2, numpy 1.26.4, scipy 1.13.1, and the standard library.
- **You may** consult Lectures 1, 2, 3, the Mutect2 paper (Cibulskis et al. 2013), the Strelka2 paper (Kim et al. 2018), the COSMIC signatures paper (Alexandrov et al. 2020), the COSMIC database paper (Sondka et al. 2024), the OncoKB paper (Chakravarty et al. 2017), the CIViC paper (Griffith et al. 2017), and the GATK Best Practices for somatic variant discovery.
- **You may NOT** copy a pre-written somatic pipeline from the internet. Reading the GATK Best Practices article for inspiration is fine; copy-pasting a complete pipeline is not.
- **You must** pin every tool version, every input mode flag, the reference build (GRCh38), the PON source, the germline-resource VCF source and date, the contamination model, the COSMIC catalog version (v3.3), the genome build for SigProfiler, and the seed (where any tool accepts one).
- **You must** verify the BAM `@RG SM:` headers and the tumor / normal flag mapping before Mutect2 runs.
- **You should** ship a Markdown QC report alongside the run-info JSON; the JSON is for machines and the Markdown is for humans.
- **You must** commit the input BAMs (small), the references, the PON / germline / common-biallelic subsets, the unfiltered VCF, the contamination table, the filtered VCF, the spectrum TSV, the SigProfiler activities + stats, the Markdown reports, and the `run-info.json`. Gitignore the GATK intermediate logs, the SigProfiler reference-genome cache (~3 GB; reproducible), and any large reference files that exceed 100 MB.
- **The educational-research-only disclaimer is mandatory** in both the mini-project README and at the top of the Markdown reports. Removing it makes the pipeline misrepresent itself.

---

## Acceptance criteria

- [ ] `mini-project/cancer_pipeline.py` exports a function `run_cancer_pipeline(tumor_bam: Path, normal_bam: Path, reference: Path, out_dir: Path) -> Path` that runs the full pipeline and returns the path to the final `run-info.json`.
- [ ] The pipeline implements **seven stages**:
  1. **Validate inputs.** Confirm every input exists; verify BAM `@RG SM:` headers; verify reference indexes; check that the tumor and normal sample names differ.
  2. **Mutect2.** Run GATK Mutect2 in tumor-normal mode with `--panel-of-normals` and `--germline-resource`. Pin `--native-pair-hmm-threads`.
  3. **GetPileupSummaries + CalculateContamination.** One pileup per sample; matched-normal contamination estimate.
  4. **FilterMutectCalls.** Apply the standard filter set with `--contamination-table`.
  5. **Build the 96-class spectrum.** Walk PASS SNVs; pyrimidine-normalize; emit TSV.
  6. **SigProfilerAssignment.** Decompose against COSMIC v3.3 SBS; parse activities and cosine.
  7. **Render reports.** `filter_tally.md`, `signature_summary.md`, `qc_report.md`, `run-info.json`.
- [ ] `cancer_pipeline.py` produces:
  - `results/unfiltered.vcf.gz` and `.tbi`.
  - `results/{tumor,normal}.pileups.table`.
  - `results/contamination.table`.
  - `results/filtered.vcf.gz` + `.tbi`.
  - `results/filter_tally.md`.
  - `results/pass.vcf.gz` (the PASS-only subset, created with `bcftools view -f PASS`).
  - `results/spectrum_96.tsv`.
  - `results/sigprofiler_out/...`.
  - `results/signature_summary.md`.
  - `results/qc_report.md`.
  - `results/run-info.json`.
- [ ] `mini-project/README.md` is a one-page (~700-1000 word) report containing:
  - One-sentence description of the input pair, the caller, the filter pipeline, and the signature method.
  - Methods section in C10 voice: every tool pinned ("GATK 4.5.0.0", "SigProfilerAssignment 0.1.4", "SigProfilerMatrixGenerator 1.2.26", "COSMIC v3.3"), every parameter explicit, the reference build stated, the PON / germline-resource sources cited.
  - Results section in C10 voice: candidate variant count, PASS count, contamination estimate, top three signatures with fractional contributions, the cosine similarity.
  - Discussion section: 200-300 words on the limits. Tumor purity, contamination, sub-clonality, the count-vs-signal trade-off in signature decomposition, the reference-build risk, the educational-use disclaimer, the things the pipeline does *not* do (CNVs, structural variants, gene fusions).
- [ ] `run.sh` is a single bash script that, given a fresh checkout + `conda env create -f env.yml + conda activate cancer-w11 + SigProfiler GRCh38 install`, reproduces the entire pipeline from scratch in under fifteen minutes.
- [ ] The repo is **public** and at least one classmate or instructor has been added as a collaborator.
- [ ] **Most importantly**: the run-info JSON is complete, the educational-use disclaimer is prominent, and the discussion section is honest about the limits.

---

## Suggested approach (rough timeline)

### Wednesday (1 hour)

1. (15 min) `git clone`, set up `mini-project/` directory.
2. (30 min) Write `env.yml` with pinned tool versions; create the conda env (`conda env create -f env.yml`); install SigProfiler GRCh38 reference (`python -c "from SigProfilerMatrixGenerator import install; install.install('GRCh38')"`); verify each tool is on the PATH (gatk, samtools, bcftools, pysam, SigProfilerAssignment).
3. (15 min) Read this README end to end and the Challenge 1 / Challenge 2 READMEs. Sketch the seven-stage flow on paper.

### Thursday (2 hours)

1. (60 min) Implement Stages 1-2 (validate, Mutect2) in `cancer_pipeline.py`. Reuse the Exercise 1 helpers; they cover most of this.
2. (60 min) Implement Stages 3-4 (GetPileupSummaries, CalculateContamination, FilterMutectCalls). Reuse the Exercise 2 helpers.

### Friday (2 hours)

1. (60 min) Implement Stages 5-6 (96-class spectrum, SigProfilerAssignment). Reuse the Exercise 3 helpers.
2. (30 min) Implement Stage 7 (render the reports). Write the qc_report.md template; it should combine filter_tally.md and signature_summary.md plus the run-info summary.
3. (30 min) Write `run.sh` and verify end-to-end reproducibility on a fresh `conda env create`.

### Saturday (2 hours)

1. (45 min) Optional: complete Challenge 1 (Strelka2 cross-check) and inline the result into the mini-project's qc_report.md.
2. (45 min) Write the README.md write-up (methods + results + 250-word discussion).
3. (30 min) Push, add a collaborator, write the commit message with specific numbers.

---

## Tips for the write-up

- **Lead with numbers.** "Mutect2 4.5.0.0 in tumor-normal mode on the 1.4 Mb chr22 didactic BAM pair emitted 217 candidate variants (198 SNVs, 19 indels) in 3 minutes 12 seconds; CalculateContamination estimated 1.2% cross-sample contamination; FilterMutectCalls flagged 38 as germline, 12 as weak_evidence, 9 as panel_of_normals, and PASSed 158 variants (73% of candidates). The 142 PASS SNVs with valid trinucleotide context decomposed against COSMIC v3.3 SBS into SBS5 (41%, clock-like), SBS1 (34%, methylation deamination), and SBS18 (19%, reactive oxygen) with reconstructed-spectrum cosine similarity 0.93. Concordance with Strelka2 PASS-PASS was Jaccard 0.81 across 170 union variants."
- **Defend every parameter.** "We used GATK 4.5.0.0 against GRCh38 with the Broad-public PON (`1000g_pon.hg38.vcf.gz`; chr22 subset) and the gnomAD v2 allele-frequency annotation as the germline resource (`af-only-gnomad.hg38.vcf.gz`; chr22 subset). FilterMutectCalls consumed the CalculateContamination output via `--contamination-table`. We used `--native-pair-hmm-threads 4` because the chr22 dataset is small. We pinned the SigProfilerAssignment COSMIC version to v3.3 and the genome build to GRCh38 to match the BAM alignment build."
- **Name two failure modes the pipeline is susceptible to.** "The 142-SNV count is at the low end for stable signature decomposition (Alexandrov-lab recommendation: at least 100 SNVs); the reported fractional contributions are stable enough to read but the signatures with < 5% contribution are noise-dominated. The pipeline does not call CNVs; if a PASS variant is in a region with loss-of-heterozygosity, the biological interpretation (biallelic loss) would differ from the SNV-only reading."
- **Acknowledge what is missing.** "The pipeline does not implement Funcotator gene-level annotation, OncoKB / CIViC clinical-interpretation lookup, or copy-number / structural-variant calling. These are Challenge 2 / extension territory. The pipeline produces a research-grade somatic variant set and a signature decomposition; it does not produce clinical advice and is not validated for any patient-care use."
- **Carry the disclaimer.** Every Markdown output of the pipeline reproduces the educational-research-only disclaimer at the top.

---

## Stretch goals (optional)

- Run the Strelka2 cross-check inline (Challenge 1) and report the Jaccard index in the main qc_report.md.
- Run the OncoKB / CIViC interpretation pass (Challenge 2) and report the top-3 clinically-actionable variants in the main qc_report.md.
- Substitute the synthetic data for a TCGA tumor-normal pair (requires dbGaP authorization or an open-access tier sample). Re-run and observe how the signature decomposition changes; report whether SBS4 (tobacco), SBS3 (HRD), or another signature dominates.
- Add a per-AF-bin Mutect2 / Strelka2 agreement plot (Homework Optional Problem extension) to the qc_report.md.
- Implement a per-variant probability of being a true somatic call by combining the Mutect2 TLOD with the gnomAD allele-frequency annotation by hand. Compare with the FilterMutectCalls decision.

## What to commit

By the end of the mini-project your `week-11/mini-project/` should contain:

```text
mini-project/
    README.md
    env.yml
    run.sh
    cancer_pipeline.py
    starter.py
    data/
        tumor_chr22.bam + .bai
        normal_chr22.bam + .bai
        chr22_GRCh38.fasta + .fai + .dict
        chr22_pon.vcf.gz + .tbi
        chr22_gnomad.vcf.gz + .tbi
        chr22_common_biallelic.vcf.gz + .tbi
        chr22_intervals.bed
        truth_chr22.vcf.gz + .tbi               (optional, for precision/recall check)
    results/
        unfiltered.vcf.gz + .tbi + .stats
        tumor.pileups.table
        normal.pileups.table
        contamination.table
        filtered.vcf.gz + .tbi
        filter_tally.md
        pass.vcf.gz + .tbi
        spectrum_96.tsv
        sigprofiler_out/...                     (the Assignment_Solution tree)
        signature_summary.md
        qc_report.md
        run-info.json
```

Gitignore the SigProfilerMatrixGenerator reference-genome cache (`~/.SigProfilerMatrixGenerator/`; ~3 GB; reproducible from a one-line install command), the Funcotator data-sources bundle if you used it (~30 GB), and any uncompressed reference FASTA files that exceed 100 MB. The commit message for the final pipeline run should be specific, e.g. `mini-project: somatic pipeline on chr22 sim pair, 158 PASS, top sigs SBS1+SBS5+SBS18, cosine 0.93`.

---

## What this mini-project does NOT do

- **It does not call CNVs.** Copy-number changes are routinely involved in tumor biology (LOH at tumor suppressors, focal amplifications of oncogenes); a SNV-only pipeline systematically under-interprets a tumor. Lecture 1 §6 discusses this.
- **It does not call structural variants.** Translocations, large inversions, complex rearrangements. Tools like Manta, GRIDSS, or SvABA fill this gap; not in Week 11.
- **It does not produce a clinical report.** OncoKB / CIViC lookups (Challenge 2) are research-grade summaries. A clinical report requires accredited-lab validation, a molecular pathologist's signature, and a multi-disciplinary tumor board's review.
- **It does not detect contamination at very low (< 0.5%) levels.** The CalculateContamination model has a noise floor; a 0.2% contamination from a high-AF germline variant of another patient can still produce false positives at AF ~0.2%.
- **It does not handle hematologic malignancies correctly.** Leukemia patients' blood is leukemic; the matched-normal model needs a non-hematopoietic normal (buccal swab, skin). The didactic pipeline assumes the standard solid-tumor / blood-normal pair.

Be honest about all of the above in your write-up. The honesty is the deliverable.
