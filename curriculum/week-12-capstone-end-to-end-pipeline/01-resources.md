# Week 12 — Resources

> **Educational and research use only.** Every URL on this page links to a free, open-access resource. Snakemake, Nextflow, bioconda, Conda, Singularity / Apptainer, Docker, and Zenodo are all free for academic and personal use. Nothing on this page is a clinical decision-support tool; capstone outputs are educational artefacts only.

This page collects the reference papers, documentation, and primary sources for Week 12. Bookmark it; you will hit it many times this week and again next semester when you re-run your capstone.

## Primary papers

### Snakemake

- **Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., and Köster, J. (2021). Sustainable data analysis with Snakemake.** *F1000Research* 10:33. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8114187/>. The 2021 paper is the canonical Snakemake reference. It supersedes the 2012 *Bioinformatics* paper (Köster and Rahmann 2012, *Bioinformatics* 28:2520) as the citation; the 2012 paper introduces the system, the 2021 paper documents its mature state. Sections to read: §2 (the DAG execution model), §3 (the conda integration), §5 (the cluster execution and the Singularity integration).
- **Köster, J. and Rahmann, S. (2012). Snakemake — a scalable bioinformatics workflow engine.** *Bioinformatics* 28(19):2520-2522. Free at <https://academic.oup.com/bioinformatics/article/28/19/2520/290322>. The original Snakemake paper. Cite this only when you need the historical reference; the 2021 paper is the modern citation.

### Nextflow

- **Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., and Notredame, C. (2017). Nextflow enables reproducible computational workflows.** *Nature Biotechnology* 35(4):316-319. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8224876/>. The canonical Nextflow paper. The dataflow / channel model in §2-§3 is the technical differentiator from Snakemake's file-pattern model.
- **nf-core: Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., and Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines.** *Nature Biotechnology* 38(3):276-278. Free at <https://www.nature.com/articles/s41587-020-0439-x>. The community library of Nextflow pipelines, with strict reproducibility conventions. <https://nf-co.re/> hosts 80+ pipelines as of 2026.

### Conda and bioconda

- **Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J., and The Bioconda Team (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences.** *Nature Methods* 15(7):475-476. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11070151/>. The bioconda paper. The continuous-integration model and the 8000+-package channel are documented here.
- **Conda documentation:** <https://docs.conda.io/projects/conda/en/stable/>. The user guide. The relevant chapters: "Managing environments", "Building and using packages", "Configuration".
- **Mamba (drop-in faster Conda replacement):** <https://mamba.readthedocs.io/>. The C++-reimplemented dependency solver. 10-100x faster than Conda's Python solver for large environments. Most Snakemake / Nextflow users have switched to mamba via `conda install -n base mamba` and then `mamba env create -f environment.yml`. Snakemake honours the `--conda-frontend mamba` flag.

### Singularity / Apptainer

- **Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute.** *PLoS ONE* 12(5):e0177459. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426675/>. The original Singularity paper. The single-file SIF format and the rootless execution model are documented in §2-§3.
- **Apptainer documentation:** <https://apptainer.org/docs/>. Singularity forked in 2021; Apptainer is the community-maintained continuation, governed by the Linux Foundation. The `apptainer` and `singularity` CLIs are largely interchangeable for the C10-Week-12 use cases. The capstone uses `apptainer` as the canonical CLI name.
- **Singularity (Sylabs commercial fork) documentation:** <https://docs.sylabs.io/>. The Sylabs-maintained commercial Singularity branding. Compatible with Apptainer for most read-only operations on `.sif` images.

### Docker

- **Merkel, D. (2014). Docker: Lightweight Linux Containers for Consistent Development and Deployment.** *Linux Journal* 2014(239):2. Free at <https://www.linuxjournal.com/content/docker-lightweight-linux-containers-consistent-development-and-deployment>. The 2014 article that introduced Docker to the developer community.
- **Docker documentation:** <https://docs.docker.com/>. The Dockerfile reference is the relevant section for capstone work.

### Zenodo and DOI deposit

- **Sicilia, M. A., García-Barriocanal, E., and Sánchez-Alonso, S. (2017). Community curation in open dataset repositories: insights from Zenodo.** *Procedia Computer Science* 106:54-60. Free at <https://www.sciencedirect.com/science/article/pii/S187705091730159X>. The community-curation model for Zenodo.
- **Zenodo:** <https://zenodo.org/>. The CERN-hosted free data repository. Up to 50 GB per upload. Issues DOIs. Has a one-click GitHub integration that archives every release.
- **DataCite (the DOI registry for research data):** <https://datacite.org/>. Zenodo issues DOIs through DataCite.

### Reproducibility studies

- **Ioannidis, J. P. A., Allison, D. B., Ball, C. A., Coulibaly, I., Cui, X., Culhane, A. C., Falchi, M., Furlanello, C., Game, L., Jurman, G., Mangion, J., Mehta, T., Nitzberg, M., Page, G. P., Petretto, E., and van Noort, V. (2009). Repeatability of published microarray gene expression analyses.** *Nature Genetics* 41(2):149-155. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2778103/>. The 2009 reproducibility-of-microarray-analyses study. 18 papers, 2 reproduced. The empirical foundation for the reproducibility-crisis discussion.
- **Mangul, S., Mosqueiro, T., Abdill, R. J., Duong, D., Mitchell, K., Sarwal, V., Hill, B., Brito, J., Littman, R. J., Statz, B., Lam, A. K.-M., Dayama, G., Grieneisen, L., Martin, L. S., Flint, J., Eskin, E., and Blekhman, R. (2019). Challenges and recommendations to improve the installability and archival stability of omics computational tools.** *PLOS Biology* 17(6):e3000333. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6592561/>. The follow-up survey: 26 of 50 surveyed tools could not be installed from published instructions without modification.
- **Beaulieu-Jones, B. K. and Greene, C. S. (2017). Reproducibility of computational workflows is automated using continuous analysis.** *Nature Biotechnology* 35(4):342-346. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6103790/>. The continuous-analysis proposal: a CI workflow that re-runs the pipeline on every commit.

### Data sources

- **Genome in a Bottle (GIAB):** <https://www.nist.gov/programs-projects/genome-bottle>. NIST-curated reference truth-set samples (HG001 / NA12878, HG002, HG005, etc.) for variant-calling benchmarking. The HG002 chr20 subset is the standard didactic dataset for variant-calling exercises. Free.
- **NCBI GEO (Gene Expression Omnibus):** <https://www.ncbi.nlm.nih.gov/geo/>. The canonical repository for transcriptomics data. Free for download. The RNA-seq capstone track uses a small GEO series (e.g. GSE52778; 8 samples) as the input.
- **NCBI Virus:** <https://www.ncbi.nlm.nih.gov/labs/virus/>. Pre-filtered viral genome sequences with metadata. The MSA + phylogeny capstone track uses 30-100 viral genomes from here.
- **NCBI SRA (Sequence Read Archive):** <https://www.ncbi.nlm.nih.gov/sra>. The canonical repository for raw sequencing reads. The long-read assembly track downloads Nanopore reads from here.
- **GENCODE:** <https://www.gencodegenes.org/>. The canonical human (and mouse) reference annotation set. v44 is the current release as of 2026. Free.
- **Ensembl:** <https://www.ensembl.org/>. Genome browser and annotation source. Free.

## Tool documentation

### Snakemake

- **Snakemake documentation:** <https://snakemake.readthedocs.io/>. The canonical user guide. The relevant chapters: "Writing workflows" (the Snakefile language), "Distribution and reproducibility" (Conda and Singularity integration), "Cluster and cloud execution".
- **Snakemake tutorial:** <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>. The official five-step tutorial. Recommended end-to-end walkthrough for first-time users.
- **Snakemake wrappers:** <https://snakemake-wrappers.readthedocs.io/>. A community-maintained library of pre-written Snakemake rules for common bioinformatics tools. Use `wrapper: "v3.0.0/bio/salmon/quant"` to invoke a pre-written rule instead of writing the shell command yourself.

### Nextflow

- **Nextflow documentation:** <https://www.nextflow.io/docs/latest/>. The canonical user guide. The relevant chapters: "Channels", "Processes", "Configuration", "Containers".
- **Nextflow training:** <https://training.nextflow.io/>. The Seqera-maintained tutorial. Recommended for first-time Nextflow users.
- **nf-core:** <https://nf-co.re/>. The community library of curated Nextflow pipelines (80+ pipelines as of 2026; rnaseq, sarek, viralrecon, atacseq, methylseq, etc.). Each pipeline is a publication-grade artefact with conventions you can lift for your own capstone.

### Conda and bioconda

- **bioconda:** <https://bioconda.github.io/>. The community channel for bioinformatics packages. 8000+ packages. Activate with `conda config --add channels conda-forge --add channels bioconda` (order matters; conda-forge first).
- **conda-forge:** <https://conda-forge.org/>. The community channel for general-purpose packages (numpy, pandas, scipy, etc.). bioconda packages depend on conda-forge packages.
- **mamba:** <https://mamba.readthedocs.io/>. Faster Conda solver. Drop-in replacement.

### Singularity / Apptainer

- **Apptainer documentation:** <https://apptainer.org/docs/>. The canonical Apptainer user guide. The relevant chapters: "Build a Container", "Definition Files", "Use a Container".
- **Singularity Hub / Apptainer registry alternatives:** Singularity Hub was deprecated in 2021; the community now uses Docker Hub plus `apptainer pull docker://` to convert Docker images on the fly. The Apptainer registry at <https://apptainer.org/cloud> is the current canonical hosting.

### Docker

- **Docker Hub:** <https://hub.docker.com/>. The canonical Docker image registry. Pull bioconda-maintained images via `docker pull quay.io/biocontainers/<tool>:<version>` (biocontainers are auto-built from bioconda recipes).
- **BioContainers:** <https://biocontainers.pro/>. The bioconda-Docker integration. Every bioconda package gets an auto-built Docker image at `quay.io/biocontainers/<package>:<version>--<build-hash>`.

## Workflow manager comparison

| Feature | Snakemake | Nextflow |
|---------|-----------|----------|
| Syntax | Python | Groovy (Java JVM) |
| Workflow model | File-pattern DAG | Channel-based dataflow |
| Conda integration | Native (`conda:` directive per rule) | Native (`conda` directive per process) |
| Singularity integration | Native (`singularity:` directive per rule) | Native (`container` directive per process) |
| Cluster integration | Native (SLURM, SGE, LSF, PBS, Kubernetes) | Native (SLURM, SGE, LSF, PBS, Kubernetes, AWS Batch, Google Life Sciences) |
| Cloud integration | Via cluster integration | Native (AWS Batch, Google Cloud, Azure) |
| DAG rendering | `snakemake --dag \| dot -Tsvg` | `nextflow run main.nf -with-dag dag.svg` |
| Community library | Snakemake-wrappers | nf-core |
| First release | 2012 | 2013 |
| Canonical paper | Mölder et al. 2021 | Di Tommaso et al. 2017 |
| License | MIT | Apache 2.0 |
| Cost | Free | Free |

**Recommendation for the capstone:** pick whichever workflow manager you find easier to read. Snakemake's Python syntax is friendlier if you have spent the prior 11 weeks of C10 in Python. Nextflow's channel model is more flexible at the streaming-data scale but has a steeper learning curve. The capstone evaluation does not care which you pick.

## Container recipe comparison

| Feature | Singularity / Apptainer | Docker |
|---------|--------------------------|--------|
| Build file | `Singularity.def` | `Dockerfile` |
| Image format | Single-file SIF | Layered images |
| Runtime daemon | None (rootless) | Daemon |
| HPC-friendly | Yes (native; no daemon) | No (requires daemon or rootless-Docker) |
| Cloud-friendly | Via `apptainer pull docker://` | Native |
| Cost | Free | Free |
| Canonical CLI | `apptainer` or `singularity` | `docker` |

**Recommendation for the capstone:** build with `Singularity.def` and deposit `.sif` files; convert from Docker via `apptainer pull docker://` when a tool is only distributed as a Docker image.

## Project track resources

### Track 1: variant discovery (GIAB)

- **Krusche, P., Trigg, L., Boutros, P. C., Mason, C. E., De La Vega, F. M., Moore, B. L., Gonzalez-Porta, M., Eberle, M. A., Tezak, Z., Lababidi, S., Truty, R., Asimenos, G., Funke, B., Fleharty, M., Chapman, B. A., Salit, M., and Zook, J. M. (2019). Best practices for benchmarking germline small-variant calls in human genomes.** *Nature Biotechnology* 37(5):555-560. Free at <https://www.nature.com/articles/s41587-019-0054-x>. The canonical benchmarking-against-GIAB paper.
- **hap.py (Illumina, the standard variant-comparison tool):** <https://github.com/Illumina/hap.py>. Free.
- **DeepVariant:** <https://github.com/google/deepvariant>. The Google-released convolutional-network variant caller. Free.
- **GATK HaplotypeCaller:** <https://gatk.broadinstitute.org/hc/en-us>. The Broad Institute's canonical germline caller. Free.

### Track 2: RNA-seq DE (GEO)

- **Salmon (Patro et al. 2017, *Nature Methods* 14:417):** <https://salmon.readthedocs.io/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/>.
- **tximport (Soneson et al. 2015, *F1000Research* 4:1521):** <https://bioconductor.org/packages/tximport/>. Free.
- **DESeq2 (Love et al. 2014, *Genome Biology* 15:550):** <https://bioconductor.org/packages/DESeq2/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/>.
- **GSE52778 (the airway smooth muscle / dexamethasone study):** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778>. 8 samples, two conditions. The canonical didactic RNA-seq dataset.

### Track 3: MSA + phylogeny (NCBI Virus)

- **MAFFT (Katoh and Standley 2013, *Mol Biol Evol* 30:772):** <https://mafft.cbrc.jp/alignment/software/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3603318/>.
- **IQ-TREE 2 (Minh et al. 2020, *Mol Biol Evol* 37:1530):** <http://www.iqtree.org/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7182206/>.
- **trimAl (Capella-Gutiérrez et al. 2009, *Bioinformatics* 25:1972):** <http://trimal.cgenomics.org/>. Free.
- **ggtree (Yu et al. 2017, *Methods Ecol Evol* 8:28):** <https://bioconductor.org/packages/ggtree/>. Free.

### Track 4: long-read assembly (SRA)

- **Flye (Kolmogorov et al. 2019, *Nature Biotechnology* 37:540):** <https://github.com/fenderglass/Flye>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8404023/>.
- **Medaka:** <https://github.com/nanoporetech/medaka>. ONT's polisher. Free.
- **BUSCO (Manni et al. 2021, *Mol Biol Evol* 38:4647):** <https://busco.ezlabgit.io/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8476166/>.
- **Quast (Gurevich et al. 2013, *Bioinformatics* 29:1072):** <http://quast.sourceforge.net/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/>.

## Reproducibility checklists

- **The TIER protocol (Teaching Integrity in Empirical Research):** <https://www.projecttier.org/>. A reproducibility checklist for empirical social-science research; applicable to bioinformatics with minor adaptation.
- **The FAIR principles (Wilkinson et al. 2016, *Scientific Data* 3:160018):** <https://www.go-fair.org/fair-principles/>. Free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4792175/>. The Findable / Accessible / Interoperable / Reusable principles for research data. The capstone should be FAIR.
- **CITATION.cff:** <https://citation-file-format.github.io/>. The standard for citation metadata in software repositories. GitHub honours `CITATION.cff` and surfaces a "Cite this repository" button on the repo page. Add one to your capstone.

## Books

- **Buffalo, V. (2015). Bioinformatics Data Skills.** O'Reilly. The canonical bioinformatics-in-Unix textbook. Chapters 1-3 cover the reproducibility and version-control patterns; Chapter 12 introduces workflow managers (pre-Snakemake-mainstream; the chapter's `make`-based examples translate cleanly to Snakemake).
- **Snakemake-wrappers documentation as a recipe book:** <https://snakemake-wrappers.readthedocs.io/>. Treat it as a recipe collection for common bioinformatics steps.

## Free MOOCs and courses

- **Bioinformatics Specialization (UC San Diego, Pevzner and Compeau):** <https://www.coursera.org/specializations/bioinformatics>. The Compeau-Pevzner "Bioinformatics Algorithms" textbook companion course. Free to audit.
- **Galaxy Training Network:** <https://training.galaxyproject.org/>. The Galaxy project's free training library. Includes Snakemake / Nextflow integration tutorials.
- **nf-core training events:** <https://nf-co.re/events>. Free virtual training sessions for Nextflow / nf-core pipelines.

## A note on what to cite

When you write the capstone `report.md`, cite the workflow manager (Snakemake / Nextflow), the container runtime (Singularity / Apptainer), the package manager (bioconda), and every primary tool you use. The convention in computational biology is to cite tools by their *paper*, not their URL. Example citation block:

> The pipeline was implemented in Snakemake (Mölder et al. 2021) with Conda (Anaconda Inc., 2024) environments specified via the bioconda channel (Grüning et al. 2018). Tools were containerized with Apptainer (Kurtzer et al. 2017). Read alignment used BWA-MEM (Li and Durbin 2009) and the GRCh38 reference. Variant calling used DeepVariant (Poplin et al. 2018). Variant comparison against the Genome in a Bottle truth set (Zook et al. 2014) used hap.py (Krusche et al. 2019). The capstone was deposited on Zenodo (Sicilia et al. 2017; DOI: <your-DOI>).

Every citation above is from a free paper. The capstone reading list is entirely open access.
