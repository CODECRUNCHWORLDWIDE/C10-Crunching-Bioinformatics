# Exercise 1 — Annotate a small VCF with VEP

> **Educational and research use only.** This exercise teaches the mechanics of running VEP against a public demo VCF. The output is not a clinical interpretation. The same disclaimer that opens every Week 8 file applies here.

**Estimated time:** 90 minutes (45 minutes setup + cache install, 30 minutes running and reading the output, 15 minutes writing the notes).

**Goal:** Install VEP locally, download a small VEP cache, annotate the bundled `data/demo.vcf` file end to end, and read the resulting `CSQ` INFO field by hand. Build the muscle memory of "VCF in, annotated VCF out, eyes on the consequence column."

## Acceptance criteria

- `vep --help` runs and prints the VEP help page.
- A VEP cache is installed at `${HOME}/.vep/homo_sapiens_vep_110_GRCh38` (or wherever you configured `--dir_cache`).
- `data/demo.vep.vcf` exists and contains a `##INFO=<ID=CSQ,...>` header line.
- For each of the ~50 variants in the demo VCF, you can name the canonical consequence, the gene symbol, and the impact tier from the CSQ field.
- `notes/exercise-01-vep-summary.md` records:
  1. The VEP version (`vep --help | head -2`).
  2. The cache version (read from the run-info output).
  3. A 3-row table summarizing impact counts (HIGH / MODERATE / LOW / MODIFIER) across the 50 variants.
  4. Three variants picked by you, with their gene, consequence, impact, SIFT, PolyPhen, and gnomAD AF copied out by hand.

## Requirements

- `conda install -c bioconda ensembl-vep=110.1` (the conda recipe; about 5 minutes).
- ~25 GB of free disk for the human GRCh38 cache. If you do not have the space, fall back to the **Codespaces or Colab path** in the appendix.
- A working internet connection for the cache download (first run only). After the cache is on disk, VEP runs offline.

## Steps

### 1. Install VEP

```bash
conda create -n c10-week08 -c bioconda -c conda-forge \
    "python=3.11" "ensembl-vep=110.1" "snpeff=5.2" \
    "cyvcf2=0.30" "pysam=0.22" "pandas=2.2" "requests=2.32" \
    "jinja2=3.1" "numpy=1.26"
conda activate c10-week08
vep --help | head -2
```

If the last line prints `#ENSEMBL VARIANT EFFECT PREDICTOR v110.1`, VEP is on your PATH.

### 2. Install the VEP cache

VEP needs a cache (~25 GB for human GRCh38 release 110). Install it once:

```bash
mkdir -p ${HOME}/.vep
vep_install \
    --AUTO cf \
    --SPECIES homo_sapiens \
    --ASSEMBLY GRCh38 \
    --CACHEDIR ${HOME}/.vep \
    --NO_HTSLIB \
    --CACHE_VERSION 110
```

Expect ~5-15 minutes of download. The cache is split into ~270 chromosome chunks under `${HOME}/.vep/homo_sapiens_vep_110_GRCh38/`. After install, `du -sh ${HOME}/.vep` should report ~25 GB.

If the download fails or you do not have space, see the **Appendix: Codespaces and Colab paths** at the end of this exercise.

### 3. Get the demo VCF

The Week 8 demo VCF lives at `data/demo.vcf.gz` in this repository. Copy it into your working directory:

```bash
mkdir -p week-08/exercise-01
cp ../../data/demo.vcf.gz week-08/exercise-01/
cd week-08/exercise-01
gunzip -k demo.vcf.gz
head -30 demo.vcf
```

You should see:

```
##fileformat=VCFv4.2
##source=c10-week08-demo
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
...
#CHROM  POS         ID    REF  ALT  QUAL  FILTER  INFO
chr1    11790853    .     C    T    500   PASS    .
chr1    27107250    .     G    A    500   PASS    .
chr2    47403534    .     C    A    500   PASS    .
chr2    47445413    .     C    T    500   PASS    .
chr3    37001008    .     G    A    500   PASS    .
chr5    112815473   .     A    T    500   PASS    .
chr7    117559590   .     G    A    500   PASS    .   <-- (CFTR-adjacent)
chr11   108235829   .     G    A    500   PASS    .
chr13   32398489    .     G    A    500   PASS    .   <-- (BRCA2 region)
chr17   43094077    .     G    A    500   PASS    .   <-- (BRCA1 region)
...
```

The file contains ~50 variants picked to cover a range of impact tiers and known disease genes. None of the variants is from a real patient; the file is a synthetic teaching VCF.

### 4. Run VEP

The canonical command:

```bash
vep \
    -i demo.vcf \
    -o demo.vep.vcf \
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

Expect ~5-15 seconds. The output `demo.vep.vcf` is a VCF with the same variants plus a `CSQ` INFO field on each.

### 5. Read the CSQ header

```bash
grep '##INFO=<ID=CSQ' demo.vep.vcf
```

You should see one long line that ends with:

```
... Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|...|SIFT|PolyPhen|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADg_AF">
```

Note the field order. We use it to parse the CSQ field in Exercise 3.

### 6. Read a few CSQ records by hand

```bash
grep -v '^##' demo.vep.vcf | head -5
```

The first non-comment line is the column header. The next lines are variants. For the BRCA1-region variant at `chr17:43094077`:

```
chr17  43094077  .  G  A  500.0  PASS  CSQ=A|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000357654|protein_coding|10/23||ENST00000357654.9:c.4837G>A|ENSP00000350283.5:p.Glu1613Lys|4837|4837|1613|E/K|GAA/AAA|rs28897696|||-1|YES|ENST00000357654.9|deleterious(0.02)|probably_damaging(0.94)|0.0001|0|0.0001|0|0.0001|0|0.00012|0.00018,A|...
```

Parse this by hand. Split on `,` first (one record per overlapping transcript), then split each on `|`. The fields, in order: `Allele=A, Consequence=missense_variant, IMPACT=MODERATE, SYMBOL=BRCA1, Gene=ENSG00000012048, Feature_type=Transcript, Feature=ENST00000357654, BIOTYPE=protein_coding, EXON=10/23, INTRON=, HGVSc=ENST00000357654.9:c.4837G>A, HGVSp=ENSP00000350283.5:p.Glu1613Lys, cDNA_position=4837, CDS_position=4837, Protein_position=1613, Amino_acids=E/K, Codons=GAA/AAA, Existing_variation=rs28897696, DISTANCE=, STRAND=-1, FLAGS=, SYMBOL_SOURCE=, HGNC_ID=, CANONICAL=YES, MANE_SELECT=ENST00000357654.9, SIFT=deleterious(0.02), PolyPhen=probably_damaging(0.94), AF=0.0001, AFR_AF=0, AMR_AF=0.0001, EAS_AF=0, EUR_AF=0.0001, SAS_AF=0, gnomADe_AF=0.00012, gnomADg_AF=0.00018`.

Now you have a row of evidence: a known variant (rs28897696), in BRCA1, a missense at residue 1613 (E -> K), with SIFT deleterious and PolyPhen probably damaging, at gnomAD AF ~0.00012 (rare). This is exactly the kind of row the report you will build in Exercise 3 turns into one CSV line.

### 7. Build the impact summary

```bash
# Count variants by impact tier. Quick-and-dirty awk one-liner.
grep -v '^#' demo.vep.vcf | \
    awk -F'\t' '{
        n = split($8, a, ";");
        for (i = 1; i <= n; i++) {
            if (a[i] ~ /^CSQ=/) {
                csq = substr(a[i], 5);
                split(csq, b, ",");
                # take only the first (canonical-ish) record per variant.
                split(b[1], c, "|");
                impact = c[3];
                counts[impact]++;
            }
        }
    } END {
        for (k in counts) printf "%-12s %d\n", k, counts[k];
    }' | sort -k1
```

You should see something close to:

```
HIGH         3
LOW          8
MODERATE     27
MODIFIER     12
```

(Numbers may vary by ~+/-2 if the demo VCF is updated between releases.)

### 8. Pick three variants and write them up

In `notes/exercise-01-vep-summary.md`, write a table with three picked variants:

```markdown
| chrom | pos | ref | alt | gene | consequence | impact | rsid | sift | polyphen | gnomADe_AF |
|-------|-----|-----|-----|------|-------------|--------|------|------|----------|------------|
| chr17 | 43094077 | G | A | BRCA1 | missense_variant | MODERATE | rs28897696 | deleterious(0.02) | probably_damaging(0.94) | 0.00012 |
| chr7  | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| chr13 | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
```

Pick one HIGH-impact variant, one MODERATE-impact variant, and one LOW or MODIFIER variant. Read the row by hand — do not script it. The point is to look at the CSQ field with your own eyes.

### 9. Commit

```bash
git add demo.vep.vcf notes/exercise-01-vep-summary.md
git commit -m "Exercise 1: VEP-annotated demo VCF, 3 hand-curated examples"
```

## Common mistakes

- **Forgetting `--offline` or `--cache`.** Without these, VEP tries to query Ensembl's REST API and fails offline (or takes minutes per variant online). Always pass `--cache --offline` once the cache is installed.
- **Using a different cache version than the one you installed.** The `--cache_version` and the directory under `${HOME}/.vep/homo_sapiens_vep_<N>_GRCh38` must match. Mismatched versions error out with a "cache not found" message.
- **Omitting `--canonical --mane --symbol`.** Without these flags, the CSQ field has no `SYMBOL` or `CANONICAL` field, and downstream parsing in Exercise 3 fails. They cost nothing; always include them.
- **Forgetting `--sift b --polyphen b`.** The default is `--sift p --polyphen p` (prediction only, no score). For ACMG PP3/BP4 you need the scores. Use `b` for "both."
- **Running on a wrong-assembly VCF.** If your VCF has GRCh37 coordinates but you ran VEP with `--assembly GRCh38`, the annotations will be on the wrong genes. Check the VCF's `##reference=` header line before running VEP and use the matching assembly.

## Appendix: Codespaces and Colab paths (if you do not have 25 GB free)

### GitHub Codespaces

The Week 8 curriculum repo includes a `.devcontainer/devcontainer.json` that pre-installs VEP plus the cache. Launch a codespace from the repo, and the first terminal command `vep --help` works out of the box. The cache lives at `/workspaces/.vep/`.

### Google Colab

Add the following cell at the top of your Colab notebook:

```python
!conda install -c bioconda -y "ensembl-vep=110.1" 2>&1 | tail
!vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 \
    --CACHEDIR /content/.vep --NO_HTSLIB --CACHE_VERSION 110
```

The Colab free-tier disk is 50-100 GB depending on the runtime, so the cache fits. The install takes ~10 minutes on Colab's slower network.

### Ensembl REST API fallback (no cache)

If you cannot install the cache at all, VEP can query Ensembl's REST API at run time. Drop `--cache --offline --dir_cache` and add `--rest_url https://rest.ensembl.org`. This is ~10x slower per variant (~1 sec per variant network round-trip) but produces the same output. Acceptable for the 50-variant demo VCF; do not use for production pipelines.

## What you learned

- The VEP CLI flag set for offline annotation.
- The structure of the `CSQ` INFO field: pipe-separated, one record per overlapping transcript, fields documented in the header.
- The impact tier (HIGH / MODERATE / LOW / MODIFIER) and how to count it by tier.
- The "read one row by hand" discipline: do not start scripting until you can read the data.

Continue to **Exercise 2 — query ClinVar and gnomAD programmatically**.
