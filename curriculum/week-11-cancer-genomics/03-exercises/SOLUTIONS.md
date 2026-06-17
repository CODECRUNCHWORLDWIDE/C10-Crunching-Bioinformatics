# Week 11 — Exercise Solutions

> **Educational and research use only.** None of the example outputs in this file are clinically actionable. The pipelines are pedagogical implementations of the underlying methods; clinical decisions require validated, accredited assays.

Worked solutions and discussion for the three Week 11 exercises. Read each exercise's docstring first; the solution notes here explain *why* each step is the way it is, not just *what* the code does.

---

## Exercise 1 — Mutect2 via subprocess

### What it does

Takes a tumor BAM, a matched normal BAM, a reference FASTA, a panel-of-normals VCF, a germline-resource VCF, and an intervals BED; verifies the sample names in each BAM's `@RG SM:` header; runs `gatk Mutect2` in tumor-normal mode via `subprocess.run(check=True, capture_output=True)`; parses the resulting unfiltered VCF; emits a `run-info.json` recording every parameter.

### The non-obvious steps

**Sample-name verification.** The most damaging silent failure in a Mutect2 pipeline is flipping the tumor and normal BAMs at the command line. `-tumor SAMPLE_X` tells Mutect2 "treat the sample named SAMPLE_X as tumor"; the sample name is matched against the BAM's `@RG SM:` header. If you accidentally pass `-tumor NORMAL_SAMPLE -normal TUMOR_SAMPLE`, Mutect2 runs the inverse comparison and produces calls for "variants in the normal that are not in the tumor", which is almost always an empty list. The result looks defensible (Mutect2 ran, no errors) but is wrong. The `get_sample_name()` helper reads `@RG SM:` directly from the BAM and the orchestrator verifies that the two BAMs have different sample names.

Run this verification *before* you call Mutect2. The verification is cheap (a header read); a Mutect2 call wasted on flipped samples is minutes-to-hours of compute.

**The `--native-pair-hmm-threads` knob.** Mutect2's compute is dominated by the pair-HMM likelihood computation. The `--native-pair-hmm-threads` flag controls the per-process thread count for that step; 4 is a reasonable default for a laptop, 16 for a workstation. Higher counts give diminishing returns past 8 because the active-region scanner becomes the bottleneck. Do not confuse this with the JVM thread count (controlled via `--java-options`); they are distinct.

**The Java heap.** GATK4 is a JVM application. The default heap may be insufficient on real-data BAMs; for whole-genome runs, pass `--java-options "-Xmx16g"`. The chr22 didactic dataset runs comfortably in 4 GB.

**Why `check=True` and `capture_output=True`.** `check=True` raises `CalledProcessError` if GATK exits non-zero; the orchestrator catches the exception, writes a partial `run-info.json` with `skipped: true` and a reason, and returns gracefully. `capture_output=True, text=True` collects stdout and stderr as strings so we can write them to a log file rather than dumping multi-line GATK progress to the console.

### A correct invocation

```bash
python exercise-01-mutect2-via-subprocess.py \
  --reference data/chr22_GRCh38.fasta \
  --tumor-bam data/tumor_chr22.bam \
  --normal-bam data/normal_chr22.bam \
  --pon data/chr22_pon.vcf.gz \
  --germline-resource data/chr22_gnomad.vcf.gz \
  --intervals data/chr22_intervals.bed \
  --out-dir results/ex01 \
  --threads 4
```

Expected output: `results/ex01/unfiltered.vcf.gz` with ~200-300 candidate variants on the didactic dataset; `results/ex01/run-info.json` with the GATK version, the two sample names, the input file MD5s, and the candidate counts.

### Common failure modes

- **`gatk: command not found`.** GATK is not installed or not on the PATH. The orchestrator's `skip_if_no_gatk=True` mode (the default) writes a `run-info.json` with `skipped: true` and exits cleanly so the rest of the curriculum can proceed; pass `--no-skip-if-missing` to make missing tools a hard error.
- **Sample names match.** If `tumor.bam` and `normal.bam` both have `SM:SAMPLE_1`, the verification step raises a `ValueError` before Mutect2 is called. This is usually a library-prep mis-tag; re-tag the BAMs with `samtools addreplacerg -r 'ID:NEW_ID\tSM:NEW_SAMPLE'`.
- **Mutect2 fails with "stats file conflict".** The output `.stats` file already exists from a previous run. Delete the output directory and re-run.

---

## Exercise 2 — FilterMutectCalls + CalculateContamination

### What it does

Takes the unfiltered Mutect2 VCF from Exercise 1; runs `gatk GetPileupSummaries` on both BAMs against a common-biallelic VCF; runs `gatk CalculateContamination` with the matched normal pileups; runs `gatk FilterMutectCalls` with `--contamination-table`; tallies the FILTER column of the filtered VCF; emits a Markdown report and a `run-info.json`.

### The non-obvious steps

**The two-step pileup + contamination flow.** GATK's contamination model is "estimate the cross-sample contamination from the allele-frequency distribution of common biallelic germline variants in the BAM". `GetPileupSummaries` walks a list of common biallelic SNVs (typically gnomAD restricted to AF 0.01-0.2; the Broad publishes a `small_exac_common_3.hg38.vcf.gz` chr-distributed subset) and counts the reads supporting each allele at each site. `CalculateContamination` then takes the pileups and computes the maximum-likelihood contamination fraction under the matched-normal model. The two steps are separate tools because the pileup is reusable across analyses (it is the same regardless of whether you are calling variants, computing contamination, or doing tumor-only analysis).

**The `--matched` flag.** `CalculateContamination -matched normal_pileups` uses the matched-normal model: the normal sample's allele frequencies at common biallelic sites are the reference distribution; deviations in the tumor pileups indicate contamination. Without `-matched`, the contamination estimate is less precise (it uses a population-level prior).

**The `--stats` flag for FilterMutectCalls.** Mutect2 emits a `.stats` file alongside the unfiltered VCF. FilterMutectCalls needs it to know the called intervals and the per-site coverage statistics. The orchestrator passes the `.stats` path explicitly via `--stats`; if the file does not exist (e.g. Mutect2 was re-run with a different output path), pass it anyway via the appropriate flag or re-run Mutect2.

**Tallying the FILTER column.** A variant's FILTER value can be empty (treated as PASS), `["PASS"]`, or a list of filter names (each contributes one count to that filter's tally). The `tally_filters()` helper iterates the VCF records and counts; multi-filter variants increment each filter's count separately. The orchestrator then reports the count of each filter; PASS is the most informative (the call set you would use downstream).

**The Markdown summary.** The pipeline writes `filter_tally.md` with a table of filter counts and the contamination estimate. The Markdown is for humans (a glance at the report tells you whether the filter distribution is healthy); the JSON is for machines.

### A correct invocation

```bash
python exercise-02-filter-and-contamination.py \
  --reference data/chr22_GRCh38.fasta \
  --unfiltered-vcf results/ex01/unfiltered.vcf.gz \
  --tumor-bam data/tumor_chr22.bam \
  --normal-bam data/normal_chr22.bam \
  --common-biallelic data/chr22_common_biallelic.vcf.gz \
  --out-dir results/ex02
```

Expected output:
- `results/ex02/contamination.table` with a single data row.
- `results/ex02/filtered.vcf.gz` with FILTER-annotated records.
- `results/ex02/filter_tally.md` showing the distribution.
- `results/ex02/run-info.json` with the PASS count (typically 50-70% of total candidates), the contamination fraction, and the per-filter counts.

A healthy report has PASS as the largest category, `germline` second-largest, and `contamination` rare.

### Common failure modes

- **Empty contamination table.** `CalculateContamination` produces no rows when the input pileups are sparse (e.g. the chr22 subset has very few common biallelic sites covered at sufficient depth). The orchestrator's `parse_contamination_table` returns `(0.0, 0.0)` in that case; FilterMutectCalls still runs but does not consume the contamination signal.
- **High `weak_evidence` count.** Suggests the input BAM coverage is low. Check `samtools coverage tumor.bam` and the `MEAN_COVERAGE` from `gatk CollectWgsMetrics`.
- **FilterMutectCalls warns about missing INFO field.** Suggests the unfiltered VCF was produced by a different Mutect2 version; re-run Mutect2 with the version pinned in `run-info.json`.

---

## Exercise 3 — 96-class spectrum + COSMIC v3 signature decomposition

### What it does

Takes the filtered Mutect2 VCF from Exercise 2; extracts the trinucleotide context of each PASS SNV via `pysam.FastaFile.fetch`; normalizes each SNV to its pyrimidine-reference form; builds the 96-class spectrum; writes it as a SigProfiler-compatible TSV; runs `SigProfilerAssignment.cosmic_fit()` against the COSMIC v3.3 SBS catalog; parses the activities TSV; emits a Markdown summary with the top signatures and the cosine similarity.

### The non-obvious steps

**Pyrimidine normalization.** Every SNV is reported on whichever strand the reference happens to be; A>G and T>C are the same mutation, observed from opposite strands. The Alexandrov-lab convention is to normalize every SNV to the pyrimidine reference: if REF is A or G, complement REF and ALT and reverse-complement the trinucleotide context. The result is that every SNV in the spectrum has a C or T reference, which gives exactly six substitution types and the 6*16 = 96 classes.

Implementation note: when you reverse-complement the context, the left and right flanks swap (the base 5' of the original substitution is the base 3' of the complemented one). The `reverse_complement()` helper handles this; the trinucleotide context emerges in the correct order.

**The canonical 96-class order.** The Alexandrov-lab order is: six substitution types in `C>A, C>G, C>T, T>A, T>C, T>G` order; within each substitution, the 16 trinucleotide contexts in nested `ACA, ACC, ACG, ACT, CCA, CCC, ...` order (left flank then right flank, A < C < G < T). SigProfilerAssignment expects this order; if you reorder you will get nonsensical signature assignments. The `canonical_96_classes()` helper enumerates them.

**The SigProfilerAssignment.cosmic_fit() API.** The Python API is:

```python
Analyzer.cosmic_fit(
    samples=str(spectrum_tsv),
    output=str(out_dir),
    input_type="matrix",
    context_type="96",
    genome_build="GRCh38",
    cosmic_version="3.3",
    collapse_to_SBS96=True,
    make_plots=False,
)
```

It writes a directory tree under `output/`. The activities are at `output/Assignment_Solution/Activities/Assignment_Solution_Activities.txt`; the sample stats (including cosine similarity) are at `output/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt`. The orchestrator finds these via `Path.rglob()` so that minor version differences in the directory naming do not break the parse.

**The cosine similarity check.** A clean decomposition has cosine similarity >= 0.85 (good) or >= 0.95 (excellent). Below 0.85 the decomposition has failed to fit a substantial fraction of the observed spectrum; report this in the Markdown summary as a flag.

**The minimum-mutation-count check.** Below 50 SNVs the decomposition is dominated by noise; the Markdown summary inserts a warning when the count is too low. On the didactic chr22 subset, the SNV count is at the low end (~50-100); on a real whole-genome tumor it would be 1,000-20,000.

### A correct invocation

```bash
python exercise-03-trinucleotide-and-signatures.py \
  --filtered-vcf results/ex02/filtered.vcf.gz \
  --reference data/chr22_GRCh38.fasta \
  --out-dir results/ex03 \
  --sample-name TUMOR \
  --genome-build GRCh38 \
  --cosmic-version 3.3
```

Expected output:
- `results/ex03/spectrum_96.tsv` with all 96 rows.
- `results/ex03/sigprofiler_out/` with the SigProfilerAssignment outputs.
- `results/ex03/signature_summary.md` with the top signatures and the cosine.
- `results/ex03/run-info.json` with the per-signature counts and the cosine.

### Common failure modes

- **`SigProfilerAssignment is not importable`.** Not installed in the active conda env. The orchestrator's `skip_if_missing=True` writes a partial `run-info.json` with `skipped: true`. Install with `conda install -c bioconda sigprofilerassignment=0.1.4`.
- **`SigProfilerMatrixGenerator reference not found`.** The GRCh38 reference data needs to be downloaded once. Run `python -c "from SigProfilerMatrixGenerator import install; install.install('GRCh38')"`; this is a ~3 GB one-time download.
- **Cosine similarity is very low (< 0.7) and SBS39 dominates.** Likely a degeneracy between SBS3 and SBS39 on a low-count spectrum. Either run on more mutations (whole-genome rather than chr22-only) or report the SBS3 / SBS39 ambiguity explicitly.
- **The trinucleotide context fetch returns 'N'.** The variant is in a region of the reference that is masked / unspecified. Skip the variant.
- **The trinucleotide context's central base does not match the VCF REF.** The reference FASTA and the VCF are on different builds, or the VCF's REF column has a strand convention different from the FASTA. Verify the FAI matches the BAM header.

---

## Cross-exercise notes

### The exercise-to-exercise pipeline

The three exercises form a pipeline: Exercise 1 produces the unfiltered VCF, Exercise 2 produces the filtered VCF, Exercise 3 produces the signature decomposition. The mini-project wraps them as a single `bash run.sh`.

### The `run-info.json` chain

Each exercise's `run-info.json` records the inputs (with MD5s where applicable), the tool versions, the parameters, and the run date. The mini-project concatenates the three into a single combined provenance record.

### Why every step has a "graceful skip" mode

Cancer-genomics tools are large, slow, and have non-trivial install requirements. A learner without GATK installed should still be able to read the exercise, understand the pattern, and produce a valid `run-info.json` showing what was attempted. The `skip_if_missing=True` (default) flag implements this: every wrapper writes a partial `run-info.json` with `skipped: true` and a `skip_reason`. The `--no-skip-if-missing` flag makes missing tools a hard error, which is appropriate for a CI-validated run.

### Reproducibility checklist

Before you publish a somatic variant call, the combined `run-info.json` records:

- GATK version, Strelka version (if used), samtools, bcftools, pysam, SigProfilerAssignment versions.
- Reference build (GRCh38) and MD5 of the FASTA.
- PON source URL and MD5.
- Germline-resource VCF source URL and MD5.
- Tumor and normal BAM paths, their `@RG SM:` sample names, and their MD5s.
- Intervals file (or "whole genome").
- Mutect2 thread count and JVM heap.
- Contamination estimate (the value from CalculateContamination, with error).
- FilterMutectCalls parameter set.
- COSMIC catalog version (v3.3), the genome build, the top signatures, the cosine similarity.
- Run date and the host.

Without these, the somatic call cannot be reproduced.

---

## Worked example: reading the output

Suppose Exercise 1 produces `results/ex01/run-info.json` with `candidate_variants: 217`, `candidate_snvs: 198`, `candidate_indels: 19`. Exercise 2's `run-info.json` reports `n_total: 217`, `n_pass: 158`, and `filter_counts: {"PASS": 158, "germline": 38, "weak_evidence": 12, "panel_of_normals": 9, "clustered_events": 6, "strand_bias": 4, "slippage": 3, "contamination": 1}` and `contamination_fraction: 0.012`. Exercise 3's `run-info.json` reports `n_snvs_used: 142` (a few PASS SNVs were skipped because the context fetch returned N), `top_signatures: [{"signature": "SBS1", "fraction": 0.34}, {"signature": "SBS5", "fraction": 0.41}, {"signature": "SBS18", "fraction": 0.19}, ...]`, and `cosine_similarity: 0.93`.

The combined story: "Mutect2 4.5.0.0 in tumor-normal mode on the chr22 didactic BAM pair emitted 217 candidate variants (198 SNVs, 19 indels) in approximately 3 minutes; FilterMutectCalls with the contamination estimate of 1.2% flagged 38 as germline, 12 as weak_evidence, 9 as panel_of_normals, and PASSed 158 (73% of candidates). SigProfilerAssignment 0.1.4 against the COSMIC v3.3 SBS catalog (genome build GRCh38) decomposed the 142 usable PASS SNV spectrum into SBS5 (41%, clock-like), SBS1 (34%, methylation deamination), and SBS18 (19%, reactive oxygen damage) with cosine similarity 0.93. The signature attribution is well-fit; the clock-like signatures dominate, consistent with a non-hypermutator background; SBS18 contribution suggests reactive-oxygen exposure as a minor mechanism."

That is the lab-notebook-voice paragraph you want in your Methods section. Replace specific numbers with your own run's numbers.
