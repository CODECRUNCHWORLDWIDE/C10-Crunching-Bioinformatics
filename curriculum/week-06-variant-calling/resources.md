# Week 6 — Resources

Every resource on this page is **free** and **publicly accessible**. Where we name a version (bcftools 1.19, GATK 4.5.0.0, ensembl-vep 110, pysam 0.22), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Li (2011)** — the `bcftools` genotype-likelihood model. *Bioinformatics* 27:2987. Free full text:
  <https://academic.oup.com/bioinformatics/article/27/21/2987/217423>
- **Danecek et al. (2021)** — `bcftools` 1.x, the modern reference for the toolchain. *GigaScience* 10:giab008. Free full text:
  <https://academic.oup.com/gigascience/article/10/2/giab008/6137722>
- **DePristo et al. (2011)** — the GATK 1.x variant-calling framework. *Nature Genetics* 43:491. Free full text:
  <https://www.nature.com/articles/ng.806>
- **Poplin et al. (2018)** — `HaplotypeCaller`, the local-reassembly caller used in GATK 4.x. *bioRxiv*. Free preprint:
  <https://www.biorxiv.org/content/10.1101/201178v3>
- **Danecek et al. (2011)** — the VCF format paper. *Bioinformatics* 27:2156. Free full text:
  <https://academic.oup.com/bioinformatics/article/27/15/2156/402296>
- **McLaren et al. (2016)** — the Ensembl Variant Effect Predictor. *Genome Biology* 17:122. Free full text:
  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4>
- **The VCFv4.3 specification** — the canonical reference. ~40 pages, denser than SAM but worth reading once:
  <https://samtools.github.io/hts-specs/VCFv4.3.pdf>
- **The GATK Best Practices** — the Broad's canonical short-variant pipeline. Read at least the "Hard-filtering germline short variants" page:
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035890471>

## Tool reference (the command-line surface)

### bcftools 1.19

The `bcftools` suite is ~25 subcommands that together read, write, filter, normalize, intersect, and annotate VCF/BCF files. It ships with the same `htslib` library as `samtools` and follows the same conventions.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `bcftools mpileup` | Walk a BAM column by column, emit per-position genotype likelihoods to a BCF | `-f` (reference FASTA, required), `-r` (region), `-Ou` (uncompressed BCF for piping), `-a` (annotations: `AD,DP,SP`), `--max-depth` (cap on per-position read count) |
| `bcftools call` | Apply the genotype model to mpileup output, emit VCF | `-m` (multiallelic-and-rare-variant model — the default modern choice), `-c` (legacy consensus model — do not use), `-v` (output variants only, drop reference-calls), `-Oz` (output bgzipped VCF), `--ploidy` (`1` for bacteria, `2` for human) |
| `bcftools filter` | Apply hard filters by `INFO`/`FORMAT` expressions | `-s` (filter tag name, e.g. `LowQual`), `-e` (exclude expression, e.g. `'QUAL<20 \|\| DP<10'`), `--SnpGap` (mask SNPs near indels), `--IndelGap` (mask close indels) |
| `bcftools norm` | Normalize indel representation (left-align, split multiallelic) | `-f` (reference FASTA), `-m` (`+any`, `-any`, `-snps`, `-indels`), `-d` (deduplicate) |
| `bcftools view` | Read/filter VCF/BCF | `-Oz` (output bgzipped VCF), `-f PASS` (keep only PASS variants), `-i` (include by expression), `-r` (region), `-s` (sample subset) |
| `bcftools query` | Format VCF as a tab-separated table | `-f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n'` (a printf-style format string) |
| `bcftools isec` | Intersect/subtract multiple VCFs | `-p outdir` (output directory), `-n` (set operation: `=2` for "in both") |
| `bcftools stats` | Summary statistics over a VCF | — |
| `bcftools index` | Build a `.csi` (default) or `.tbi` index for a bgzipped VCF | `-t` (build `.tbi` instead of `.csi`) |
| `bcftools annotate` | Add annotations from a TSV/VCF file | `-a` (annotation file), `-c` (column mapping) |

#### The canonical BAM-to-VCF one-liner

```bash
bcftools mpileup -Ou -f ref/ecoli.fa --max-depth 250 -a 'AD,DP,SP' \
    aln/sample.markdup.bam \
| bcftools call -m -v --ploidy 1 -Oz -o calls/sample.vcf.gz - \
&& bcftools index -t calls/sample.vcf.gz
```

For *E. coli* (haploid, ~4.6 Mb at 50x coverage), this takes ~30 seconds and emits ~150-200 variants. For human (diploid, 3 Gb at 30x coverage), this takes ~2-4 hours per sample with default threads and emits ~3 million variants. The `--ploidy 1` flag is critical for bacteria; the default is 2 (diploid).

#### The canonical hard-filter step

```bash
bcftools filter -Oz -o calls/sample.filtered.vcf.gz \
    -s LowQual \
    -e 'QUAL<20 || INFO/DP<10 || INFO/MQ<40 || INFO/SP>60' \
    calls/sample.vcf.gz \
&& bcftools index -t calls/sample.filtered.vcf.gz
```

Variants matching the expression get a `FILTER` value of `LowQual`; variants that pass get `PASS`. The exact threshold values are dataset-dependent — see §"GATK hard-filter thresholds" below.

### GATK 4.5.0.0

The Genome Analysis Toolkit is a ~1 GB Java application that bundles ~50 tools for the Broad Institute's variant-calling pipeline. Only a handful are relevant to Week 6:

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `gatk HaplotypeCaller` | Call germline SNPs and indels via local haplotype reassembly | `-R` (reference FASTA), `-I` (input BAM), `-O` (output VCF), `--sample-ploidy` (1 for bacteria, 2 for human), `-ERC GVCF` (emit a gVCF for joint genotyping) |
| `gatk GenotypeGVCFs` | Joint-genotype multiple gVCFs into a single VCF | `-R`, `-V` (one or more gVCFs), `-O` |
| `gatk VariantFiltration` | Apply hard filters by expression (the GATK-native equivalent of `bcftools filter`) | `--filter-expression "QD < 2.0" --filter-name "QD2"` |
| `gatk SelectVariants` | Subset/filter variants | `--select-type-to-include SNP`, `-select 'FILTER == "PASS"'` |
| `gatk VariantsToTable` | Convert VCF to a tab-separated table | `-F CHROM -F POS -F QUAL -F QD -F FS` |
| `gatk CollectVariantCallingMetrics` | QC metrics over a VCF (Ts/Tv, het/hom ratio, count by type) | `--DBSNP <known-vcf>` |

#### A minimal GATK HaplotypeCaller call

```bash
gatk HaplotypeCaller \
    -R ref/ecoli.fa \
    -I aln/sample.markdup.bam \
    -O calls/sample.gatk.vcf.gz \
    --sample-ploidy 1
```

For an *E. coli* BAM this takes ~5 minutes (GATK is slower than bcftools but does local reassembly, which is more accurate at indels). For a human genome, this takes ~12-24 hours single-threaded; in production it is parallelized over genomic intervals with `--intervals`.

### ensembl-vep 110

The Ensembl Variant Effect Predictor takes a VCF and a species-specific transcript database (cache) and adds a `CSQ` annotation to each variant naming the consequence on every overlapping transcript.

| Command | Purpose | Most-used flags |
|---------|---------|-----------------|
| `vep` | Run VEP on a VCF | `--input_file` (input VCF), `--output_file` (output VCF or text), `--species` (e.g. `escherichia_coli_str_k_12_substr_mg1655`), `--cache` (use offline cache), `--dir_cache` (cache location), `--vcf` (output as VCF with `CSQ` in `INFO`), `--symbol` (add gene symbols), `--canonical` (mark canonical transcripts) |
| `vep_install` | Install a species cache | `--AUTO` (interactive), `--SPECIES homo_sapiens` (or another species), `--ASSEMBLY GRCh38`, `--CACHEDIR` |
| `filter_vep` | Filter VEP-annotated VCF by consequence | `-i input.vcf -filter "Consequence is missense_variant"` |

#### A minimal VEP call

```bash
vep --input_file calls/sample.filtered.vcf.gz \
    --output_file calls/sample.vep.vcf \
    --species escherichia_coli_str_k_12_substr_mg1655 \
    --cache --dir_cache vep_cache/ \
    --vcf --symbol --canonical --force_overwrite
```

For an *E. coli* VCF with 200 variants, VEP runs in ~2 seconds. For a human VCF with 3 million variants, ~30 minutes offline (much longer online via the REST API).

VEP also has an HTTP REST API at <https://rest.ensembl.org/> for one-off lookups. The REST API has a rate limit of ~15 requests/second and is not appropriate for whole-VCF annotation; use the offline cache for any production use.

### pysam 0.22 (VCF support)

In addition to BAM access (Week 5), pysam wraps htslib's VCF support:

| Class / function | Purpose |
|------------------|---------|
| `pysam.VariantFile(path, mode)` | Open a VCF/BCF file. Mode `"r"` for read (auto-detects bgzipped or plain). |
| `VariantFile.header` | The `VariantHeader` object: contigs, INFO/FORMAT/FILTER definitions, sample names. |
| `VariantFile.fetch(contig, start, stop)` | Iterate over records in a region (requires `.tbi` or `.csi` index). |
| `VariantRecord.chrom`, `.pos`, `.id`, `.ref`, `.alts`, `.qual`, `.filter`, `.info`, `.samples` | The eight mandatory VCF columns plus per-sample format data. |
| `VariantRecord.info` | Dict-like access to `INFO` fields. `rec.info["DP"]` returns an int; multi-value fields return tuples. |
| `VariantRecord.samples[name]["GT"]` | Genotype as a tuple of allele indices, e.g. `(0, 1)` for heterozygous, `(1, 1)` for homozygous-alt. |

#### A 10-line variant count

```python
import pysam

vcf = pysam.VariantFile("calls/sample.filtered.vcf.gz")
total = pass_count = snp_count = indel_count = 0
for rec in vcf:
    total += 1
    if rec.filter.keys() == ["PASS"] or not rec.filter.keys():
        pass_count += 1
    if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts):
        snp_count += 1
    else:
        indel_count += 1
print(f"{total} variants, {pass_count} PASS, {snp_count} SNPs, {indel_count} indels")
```

### cyvcf2 0.30

An alternative VCF library by Brent Pedersen, faster than pysam for whole-VCF scans because it skips the python-object construction in pysam. Same API surface, different implementation. Use cyvcf2 when iterating over a 3 M-variant human VCF; use pysam when you also need BAM access in the same script.

```python
from cyvcf2 import VCF

for variant in VCF("calls/sample.filtered.vcf.gz"):
    print(variant.CHROM, variant.POS, variant.REF, variant.ALT, variant.QUAL)
```

## Reference dataset accessions

Cited by NCBI / SRA accession so you can verify your data is the same as the curriculum's:

- **`NC_000913.3`** — *Escherichia coli* str. K-12 substr. MG1655 complete genome (4,641,652 bp). The reference for the mini-project; re-used from Week 5.
- **`SRR1770413`** — Illumina HiSeq 2500 paired-end resequencing of *E. coli* K-12 MG1655. The read set; re-used from Week 5's mini-project BAM.
- **`NC_000022.11`** — Human chromosome 22 (50,818,468 bp). Used in the homework for a small slice of human variant calling.
- **GIAB `HG001` (NA12878)** — the high-confidence variant truth set for the Genome in a Bottle benchmark sample. Used in the homework for a precision/recall comparison.
- **GFF3 annotation for `NC_000913.3`** — Bacterial gene annotation. Fetched from NCBI's RefSeq FTP. Used to build the VEP cache for the mini-project.

If any of these accessions have been retired or updated by the time you take the course, swap to the current versioned accession and note it in your reproducibility receipt.

## VCF column reference — the canonical eight

| # | Name | Type | Meaning |
|---|------|------|---------|
| 1 | `CHROM` | string | Contig name (matches `@SQ SN:` in the BAM header). |
| 2 | `POS` | int | 1-based leftmost position of the variant on the reference. |
| 3 | `ID` | string | Variant identifier (`.` if none, `rsNNN` for dbSNP variants). |
| 4 | `REF` | string | Reference allele as a substring of the reference. For SNPs, 1 base; for indels, can be multiple bases including the anchor base. |
| 5 | `ALT` | string | Comma-separated alternative alleles. For SNPs, one base each; for indels, can be multiple bases. |
| 6 | `QUAL` | float | Phred-scaled probability that the call is wrong. `30` ≈ `10^-3`, `60` ≈ `10^-6`. |
| 7 | `FILTER` | string | `PASS` if the variant survives all filters, or a semicolon-separated list of failed filter names. `.` if not yet filtered. |
| 8 | `INFO` | string | Semicolon-separated `KEY=VALUE` pairs (per-variant). |

Plus, if the VCF has samples:

| # | Name | Type | Meaning |
|---|------|------|---------|
| 9 | `FORMAT` | string | Colon-separated keys (e.g. `GT:PL:DP:AD`) that apply to every per-sample column. |
| 10+ | `<sample-name>` | string | Colon-separated values matching `FORMAT` keys. |

## INFO field reference (the most common bcftools/GATK keys)

| Key | Source | Meaning |
|-----|--------|---------|
| `DP` | bcftools, GATK | Total read depth at the position (across all samples). |
| `AF` | bcftools | Allele frequency in the called sample(s). |
| `MQ` | bcftools, GATK | RMS mapping quality of supporting reads. |
| `SP` | bcftools | Phred-scaled strand-bias p-value (Fisher's exact test). |
| `FS` | GATK | Phred-scaled Fisher's strand bias p-value. The GATK analog of `SP`. |
| `SOR` | GATK | Strand-Odds-Ratio, a strand-bias alternative robust to large allele counts. |
| `MQRankSum` | GATK | Mann-Whitney U test of mapping quality, reference allele vs alternative. |
| `ReadPosRankSum` | GATK | Mann-Whitney U test of read position, reference allele vs alternative. |
| `QD` | GATK | `QUAL / DP` — quality per unit of depth. Low values indicate low-quality calls in over-covered regions. |
| `INDEL` | bcftools | Set (no value) if the variant is an indel. |
| `MQ0` | bcftools | Number of MAPQ-0 reads at the position (an extra strand-bias-like signal). |

## FORMAT field reference

| Key | Meaning |
|-----|---------|
| `GT` | Genotype: `0/0`, `0/1`, `1/1` (diploid unphased), `0\|1` (diploid phased), `0`, `1` (haploid), `.` (missing). |
| `PL` | Phred-scaled likelihoods of the three genotypes (RR, RA, AA), comma-separated. |
| `DP` | Per-sample read depth at the position. |
| `AD` | Per-sample comma-separated allele depths: reference,alt1,alt2,... |
| `GQ` | Phred-scaled genotype quality (probability the GT call is wrong). |
| `SP` | Per-sample Phred-scaled strand bias (bcftools). |

## GATK hard-filter thresholds — the canonical recipe

The Best Practices reference is <https://gatk.broadinstitute.org/hc/en-us/articles/360035890471>. The recipe for **human germline** short variants, applied with `bcftools filter` or `gatk VariantFiltration`:

### SNPs

| Filter name | Expression | What it removes |
|-------------|------------|-----------------|
| `QD2` | `QD < 2.0` | Calls where the QUAL is low relative to depth — usually low-quality calls over-supported by many low-quality reads. |
| `FS60` | `FS > 60.0` | Calls supported only on one strand. |
| `MQ40` | `MQ < 40.0` | Calls in regions of low mean mapping quality (multimappers, low-complexity). |
| `MQRankSum-12.5` | `MQRankSum < -12.5` | Calls where reads supporting the alt allele have systematically worse mapping quality than reference reads. |
| `ReadPosRankSum-8` | `ReadPosRankSum < -8.0` | Calls where the alt allele is preferentially at one end of the supporting reads (an artifact of alignment-end errors). |
| `SOR3` | `SOR > 3.0` | An alternate strand-bias metric, robust to large allele counts. |

### Indels

| Filter name | Expression | What it removes |
|-------------|------------|-----------------|
| `QD2` | `QD < 2.0` | Same as SNPs. |
| `FS200` | `FS > 200.0` | The threshold is much looser than for SNPs because indels are intrinsically more strand-biased. |
| `ReadPosRankSum-20` | `ReadPosRankSum < -20.0` | Same idea as SNPs with a looser threshold. |
| `SOR10` | `SOR > 10.0` | Loose SOR for indels. |

### Bacterial / haploid datasets

For bacteria (and any haploid dataset), `MQRankSum` and `ReadPosRankSum` are not always computed because there is no allele-specific signal to compare. A reduced recipe:

| Filter name | Expression | What it removes |
|-------------|------------|-----------------|
| `LowQUAL` | `QUAL < 30` | Low-confidence calls overall. |
| `LowDP` | `DP < 10` | Calls supported by too few reads. |
| `HighDP` | `DP > 3 * mean_DP` | Calls in over-covered (likely repetitive) regions. |
| `LowMQ` | `MQ < 40` | Calls in low-mapping-quality regions. |
| `HighSP` | `SP > 60` | Strand-biased calls. |

This is the recipe the mini-project uses for *E. coli*.

## VEP consequence terms — the ones you must memorize

The full list (~30 terms) is at <https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html>. The ones you will see daily:

| Term | What it means | Severity |
|------|---------------|----------|
| `intergenic_variant` | Not in any annotated transcript. | Lowest (often noise). |
| `synonymous_variant` | Codon-level change with no amino-acid change. | Low. |
| `missense_variant` | Single-amino-acid change. | Medium (functional impact varies). |
| `splice_region_variant` | Within 1-3 bp of a splice site. | Medium-high (may disrupt splicing). |
| `splice_donor_variant` / `splice_acceptor_variant` | At the canonical splice site. | High (almost always disruptive). |
| `stop_gained` | Substitution creates a premature stop codon. | High (truncated protein). |
| `stop_lost` | Substitution destroys the stop codon. | High (read-through into 3' UTR). |
| `start_lost` | Substitution destroys the start codon. | High (loss of protein). |
| `frameshift_variant` | Indel of non-multiple-of-3 length. | High (downstream amino-acid sequence garbled). |
| `inframe_insertion` / `inframe_deletion` | Indel of multiple-of-3 length. | Medium (in-frame amino-acid insertion/deletion). |
| `5_prime_UTR_variant` / `3_prime_UTR_variant` | In an untranslated region of a transcript. | Low. |
| `non_coding_transcript_variant` | In a non-coding RNA. | Low. |

The "severity" column is a soft ranking used by VEP's `--most_severe` flag to pick a single consequence per variant when there are multiple overlapping transcripts.

## Free books, chapters, and tutorials

- **Heng Li's blog** — the author of `bcftools` and `samtools`. His posts often explain corner cases that the papers do not:
  <https://lh3.github.io/>
- **The bcftools tutorials** at <https://samtools.github.io/bcftools/howtos/> — official, version-pinned, kept up to date with each release. The "Variant calling" and "Filtering" tutorials are essential.
- **GATK forum** — Q&A from the people who wrote the tools. Search before posting:
  <https://gatk.broadinstitute.org/hc/en-us/community/topics>
- **The Ensembl VEP tutorial** — hands-on walkthrough with example data:
  <https://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html>
- **The 1000 Genomes Project variant-calling pipeline** — for a real-world example of the full pipeline at scale:
  <https://www.internationalgenome.org/category/variant-calling/>

## Public VCF datasets to read for practice

- **dbSNP 156** — the canonical database of known human variants. For chromosome 22 only:
  <https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/>
- **gnomAD 4.0** — population allele frequencies across ~800,000 exomes and ~70,000 genomes:
  <https://gnomad.broadinstitute.org/downloads>
- **GIAB (Genome in a Bottle) NA12878 truth VCF** — the small-variant gold standard for benchmarking:
  <https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/>
- **The 1000 Genomes phase 3 VCFs** — high-coverage variant calls across 2,504 samples:
  <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>

## Open-source code to read this week

- **`bcftools/bcftools.c`** — the dispatcher for every `bcftools` subcommand. Read the `main_call` and `main_mpileup` functions in particular:
  <https://github.com/samtools/bcftools>
- **`bcftools/mcall.c`** — the multiallelic-and-rare-variant calling model. ~2,000 lines. The math of Li 2011 + Danecek 2014 in C:
  <https://github.com/samtools/bcftools/blob/develop/mcall.c>
- **`gatk/HaplotypeCaller`** — the GATK 4 reassembly caller (Java, much more complex than bcftools but well-documented):
  <https://github.com/broadinstitute/gatk>
- **`ensembl-vep`** — Perl, surprisingly readable for a 50,000-line codebase. The `Bio::EnsEMBL::VEP::Pipeline` module is the entry point:
  <https://github.com/Ensembl/ensembl-vep>

## A note on ploidy

The default `bcftools call --ploidy` is 2 (human diploid). For bacteria you must pass `--ploidy 1`. For mitochondria, mixed populations, polyploid plants, and cancer (where the effective ploidy is tumor-purity-dependent), the answer is more nuanced — see the bcftools `--ploidy` documentation for the supported special values (`1`, `2`, `X`, `Y`).

Forgetting `--ploidy 1` for *E. coli* is the most common Week 6 bug. The result is that every variant is called as heterozygous (`GT=0/1`) instead of homozygous (`GT=1`), and the QUAL values are about half what they should be. Always set ploidy.

## Tools you will install this week

- **bcftools 1.19** — `conda install -c bioconda bcftools=1.19`. Adds `bcftools` to your PATH.
- **GATK 4.5.0.0** — `conda install -c bioconda gatk4=4.5.0.0` or download from <https://github.com/broadinstitute/gatk/releases>. Adds `gatk` to your PATH.
- **ensembl-vep 110** — `conda install -c bioconda ensembl-vep=110`. Adds `vep`, `vep_install`, `filter_vep` to your PATH.
- **cyvcf2 0.30** — `pip install cyvcf2==0.30` (optional, faster than pysam for whole-VCF scans).

### A complete environment file for the week

```yaml
name: c10-week-06
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - numpy=1.26.4
  - matplotlib
  - pandas
  - biopython=1.83
  - bwa=0.7.17
  - samtools=1.19
  - bcftools=1.19
  - pysam=0.22
  - gatk4=4.5.0.0
  - ensembl-vep=110
  - sra-tools=3.0.10
  - pip
  - pip:
      - cyvcf2==0.30
```

`conda env create -f env.yml && conda activate c10-week-06`. Total size on disk after install: ~3 GB (GATK alone is ~1 GB of Java; VEP cache for *E. coli* is small but for human is ~25 GB).

---

*If a link 404s, please open an issue so we can replace it.*
