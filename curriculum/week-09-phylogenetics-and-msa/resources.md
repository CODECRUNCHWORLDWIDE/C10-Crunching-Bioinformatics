# Week 9 — Resources

> **Reproducibility note.** Every tool, file format, and reference paper on this page is free and publicly accessible. Where we name a version (MAFFT 7.526, IQ-TREE 2.3.6, RAxML-NG 1.2.2, Biopython 1.84, ete3 3.1.3, trimAl 1.4.1), use that exact version when running locally — it pins your reproducibility. If a link breaks, please open an issue.

## Required reading (work it into your week)

- **Felsenstein, Joseph (1981)** — "Evolutionary trees from DNA sequences: a maximum likelihood approach." The maximum likelihood paper. *Journal of Molecular Evolution* 17:368. Free full text at PubMed Central:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>
  This is the most important paper of the week. Felsenstein lays out the probabilistic model and the pruning algorithm that underlies every modern ML phylogenetics tool, including RAxML and IQ-TREE. Read it end to end; it is ~20 pages and is the foundation of phylogenetic inference as it is practiced today.
- **Saitou, Naruya; Nei, Masatoshi (1987)** — "The neighbor-joining method: a new method for reconstructing phylogenetic trees." The NJ paper. *Molecular Biology and Evolution* 4:406. Free full text:
  <https://academic.oup.com/mbe/article/4/4/406/1029664>
  The classic fast distance-based tree-building algorithm. Reading it once is enough; the algorithm is the kind of thing you implement in 50 lines of Python and then never look at again because Biopython does it for you.
- **Felsenstein, Joseph (1985)** — "Confidence limits on phylogenies: an approach using the bootstrap." The bootstrap paper. *Evolution* 39:783. Free preprint scan via JSTOR / archive copies:
  <https://www.jstor.org/stable/2408678>
  The original case for the column-resampling bootstrap as a branch-support measure. Short, readable, and still the canonical reference when you put bootstrap percentages on a published tree.
- **Katoh, Kazutaka; Standley, Daron M. (2013)** — "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." The MAFFT paper. *Molecular Biology and Evolution* 30:772. Free full text:
  <https://academic.oup.com/mbe/article/30/4/772/1073398>
  Tool documentation:
  <https://mafft.cbrc.jp/alignment/software/>
- **Sievers, Fabian; Wilm, Andreas; Dineen, David; Gibson, Toby J.; Karplus, Kevin; Li, Weizhong; Lopez, Rodrigo; McWilliam, Hamish; Remmert, Michael; Soeding, Johannes; Thompson, Julie D.; Higgins, Desmond G. (2011)** — "Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega." The Clustal Omega paper. *Molecular Systems Biology* 7:539. Free full text:
  <https://www.embopress.org/doi/full/10.1038/msb.2011.75>
  Tool documentation:
  <http://www.clustal.org/omega/>
- **Edgar, Robert C. (2022)** — "Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny." The MUSCLE 5 paper. *Nature Communications* 13:6968. Free full text:
  <https://www.nature.com/articles/s41467-022-34630-w>
  Tool documentation:
  <https://drive5.com/muscle/> and <https://github.com/rcedgar/muscle>
- **Minh, Bui Quang; Schmidt, Heiko A.; Chernomor, Olga; Schrempf, Dominik; Woodhams, Michael D.; von Haeseler, Arndt; Lanfear, Robert (2020)** — "IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era." The IQ-TREE 2 paper. *Molecular Biology and Evolution* 37:1530. Free full text:
  <https://academic.oup.com/mbe/article/37/5/1530/5721363>
  Tool documentation:
  <http://www.iqtree.org/doc/>
- **Kozlov, Alexey M.; Darriba, Diego; Flouri, Tomas; Morel, Benoit; Stamatakis, Alexandros (2019)** — "RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference." The RAxML-NG paper. *Bioinformatics* 35:4453. Free full text:
  <https://academic.oup.com/bioinformatics/article/35/21/4453/5487384>
  Tool documentation:
  <https://github.com/amkozlov/raxml-ng/wiki>
- **Huerta-Cepas, Jaime; Serra, Francois; Bork, Peer (2016)** — "ETE 3: reconstruction, analysis, and visualization of phylogenomic data." The ete3 paper. *Molecular Biology and Evolution* 33:1635. Free full text:
  <https://academic.oup.com/mbe/article/33/6/1635/2579822>
  Tool documentation:
  <http://etetoolkit.org/docs/latest/>
- **Cock, Peter J. A.; Antao, Tiago; Chang, Jeffrey T.; Chapman, Brad A.; Cox, Cymon J.; Dalke, Andrew; Friedberg, Iddo; Hamelryck, Thomas; Kauff, Frank; Wilczynski, Bartek; de Hoon, Michiel J. L. (2009)** — "Biopython: freely available Python tools for computational molecular biology and bioinformatics." The Biopython paper. *Bioinformatics* 25:1422. Free full text:
  <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>
  Tutorial (the canonical reference for `Bio.AlignIO`, `Bio.Phylo`, `Bio.Phylo.TreeConstruction`):
  <https://biopython.org/docs/latest/Tutorial/index.html>
- **Capella-Gutierrez, Salvador; Silla-Martinez, Jose M.; Gabaldon, Toni (2009)** — "trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses." The trimAl paper. *Bioinformatics* 25:1972. Free full text:
  <https://academic.oup.com/bioinformatics/article/25/15/1972/213148>
  Tool documentation:
  <http://trimal.cgenomics.org/>
- **Hoang, Diep Thi; Chernomor, Olga; von Haeseler, Arndt; Minh, Bui Quang; Vinh, Le Sy (2018)** — "UFBoot2: improving the ultrafast bootstrap approximation." The UFBoot2 paper. *Molecular Biology and Evolution* 35:518. Free full text:
  <https://academic.oup.com/mbe/article/35/2/518/4565479>

## Tool reference (the command-line surface)

### MAFFT 7.526

MAFFT is the default modern multiple sequence aligner. It is a C program with a Perl front-end; the conda package installs the `mafft` binary plus the support scripts.

| Flag | Purpose |
|------|---------|
| `--auto` | Automatic algorithm selection by input size (FFT-NS-2 for small, FFT-NS-i for medium, L-INS-i for small high-accuracy). Not reproducible across MAFFT versions; pin a specific algorithm in production. |
| `--retree 2` | Two rounds of progressive alignment (the FFT-NS-2 algorithm). The fastest default. |
| `--maxiterate 1000` | Iterative refinement with up to 1,000 iterations (the FFT-NS-i algorithm). Higher quality, ~2-5x slower. |
| `--localpair --maxiterate 1000` | L-INS-i: local pairwise alignment plus iterative refinement. Highest-accuracy on small inputs (< ~200 sequences). |
| `--globalpair --maxiterate 1000` | G-INS-i: global pairwise alignment plus iterative refinement. For sequences with similar lengths end-to-end. |
| `--nuc` | Nucleotide input (auto-detected; pin in scripts). |
| `--amino` | Amino acid input (auto-detected; pin in scripts). |
| `--anysymbol` | Permit ambiguous IUPAC symbols (N, X, etc.) without warning. |
| `--quiet` | Suppress progress messages on stderr. |
| `--thread 4` | Use 4 threads. |

#### The canonical MAFFT call (Week 9 default)

```bash
mafft \
    --retree 2 \
    --maxiterate 0 \
    --nuc \
    --anysymbol \
    --quiet \
    --thread 4 \
    input.fasta \
    > aligned.fasta
```

For the demo cytochrome b panel (~10 sequences, ~1,140 columns) this runs in ~0.4 seconds. For an L-INS-i refinement add `--localpair --maxiterate 1000`; expect ~3 seconds on the same input. MAFFT writes the aligned FASTA to stdout; we redirect to a file.

**Reproducibility caveat.** MAFFT's `--auto` is *not* reproducible across versions: the same input can pick a different algorithm under MAFFT 7.490 vs 7.526 and produce a different alignment. Always pin the algorithm (`--retree 2`, `--maxiterate 1000`, or `--localpair --maxiterate 1000`) in any script that needs to be reproducible.

### Clustal Omega 1.2.4

Clustal Omega is the modern successor to ClustalW. It uses HMM profiles and scales to very large inputs.

| Flag | Purpose |
|------|---------|
| `-i input.fasta` | Input FASTA. |
| `-o aligned.fasta` | Output FASTA. |
| `--seqtype=DNA` or `--seqtype=Protein` | Sequence type. |
| `--full` | Use the full distance matrix (slower; higher quality). |
| `--iter=N` | N iterations of refinement (default 0). |
| `--threads=4` | 4 threads. |
| `--force` | Overwrite existing output. |

```bash
clustalo \
    -i input.fasta \
    -o aligned.clustalo.fasta \
    --seqtype=DNA \
    --iter=2 \
    --threads=4 \
    --force
```

Clustal Omega is comparable to MAFFT on most inputs but is noticeably faster on > 10,000-sequence inputs.

### MUSCLE 5.1

MUSCLE 5 is the 2022 rewrite with the super5 algorithm.

| Flag | Purpose |
|------|---------|
| `-align input.fasta` | Standard alignment (the PPP algorithm; high accuracy). |
| `-super5 input.fasta` | The super5 algorithm; faster on large inputs. |
| `-output aligned.fasta` | Output path. |
| `-threads 4` | 4 threads. |

```bash
muscle \
    -align input.fasta \
    -output aligned.muscle.fasta \
    -threads 4
```

MUSCLE 5 reports the highest published per-second accuracy on small-to-medium inputs (see the *Nature Communications* benchmark in Edgar 2022). It is interchangeable with MAFFT for most Week 9 purposes.

### IQ-TREE 2.3.6

IQ-TREE 2 is the modern fast maximum-likelihood phylogenetic inference tool.

| Flag | Purpose |
|------|---------|
| `-s alignment.fasta` | Input alignment. |
| `-m GTR+I+G` | Substitution model. `MFP` runs ModelFinder Plus to pick the best model automatically. |
| `-B 1000` | Ultrafast bootstrap (Hoang et al. 2018) with 1,000 replicates. |
| `-bb 1000` | Same as `-B` (older flag spelling). |
| `-alrt 1000` | SH-like approximate likelihood ratio test with 1,000 replicates (complementary support). |
| `-T 4` or `-nt 4` | 4 threads. |
| `-seed 42` | Random seed for the ML search and the bootstrap. |
| `-pre prefix` | Output file prefix. |
| `-redo` | Overwrite existing output files. |
| `-o outgroup_taxon` | Specify an outgroup for rooting. |

#### The canonical IQ-TREE 2 call (challenge)

```bash
iqtree2 \
    -s aligned.fasta \
    -m GTR+I+G \
    -B 1000 \
    -alrt 1000 \
    -T 4 \
    -seed 42 \
    -pre results/iqtree \
    -redo
```

For the demo cytochrome b panel this runs in ~15-30 seconds on a laptop. The output files are: `results/iqtree.treefile` (the ML tree in Newick), `results/iqtree.iqtree` (the human-readable log with the model, the log-likelihood, the AIC/BIC, and the per-branch support values), `results/iqtree.log` (the run log), and several intermediate files (`*.ckp.gz`, `*.bionj`, `*.mldist`) that are safe to delete.

### RAxML-NG 1.2.2

RAxML-NG is the successor to the classical RAxML. The CLI is a little different from IQ-TREE 2 but the output is the same shape: a Newick tree with bootstrap support.

```bash
raxml-ng \
    --all \
    --msa aligned.fasta \
    --model GTR+G \
    --bs-trees 1000 \
    --threads 4 \
    --seed 42 \
    --prefix results/raxmlng
```

`--all` is shorthand for "search + bootstrap + bootstrap-mapped support tree." For Week 9's small inputs, IQ-TREE 2 and RAxML-NG produce essentially the same tree. We default to IQ-TREE 2 because the CLI is slightly cleaner and the ultrafast bootstrap is faster on small alignments.

### Biopython 1.84

Biopython is the Python lab assistant. The relevant Week 9 modules:

- `Bio.SeqIO` — FASTA, GenBank, etc. read and write. (Carried over from Week 2.)
- `Bio.AlignIO` — multiple sequence alignment read and write. Formats: `fasta`, `phylip`, `clustal`, `stockholm`, `nexus`.
- `Bio.Phylo` — tree read and write. Formats: `newick`, `nexus`, `phyloxml`.
- `Bio.Phylo.TreeConstruction` — `DistanceCalculator`, `DistanceTreeConstructor` (NJ, UPGMA).
- `Bio.Phylo.Consensus` — bootstrap and consensus trees.

```python
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

aln = AlignIO.read("aligned.fasta", "fasta")
calc = DistanceCalculator("identity")  # raw p-distance; for JC/K2P we compute by hand
matrix = calc.get_distance(aln)
tree = DistanceTreeConstructor().nj(matrix)
Phylo.write(tree, "tree.nwk", "newick")
```

For nucleotide-specific JC and K2P distances, Biopython's `DistanceCalculator` does not implement the substitution-model corrections; we compute them by hand in Exercise 2.

### ete3 3.1.3

ete3 is the publication-quality tree renderer. Install via conda (`conda install -c bioconda ete3`) or pip (`pip install ete3`); the pip install also pulls `PyQt5` for the GUI viewer, which is optional for headless use.

```python
from ete3 import Tree, TreeStyle

tree = Tree("((A:0.1,B:0.2):0.3,C:0.4);", format=1)
ts = TreeStyle()
ts.show_branch_support = True
ts.show_branch_length = True
tree.render("tree.png", w=800, tree_style=ts)
```

For Week 9 we lazy-import ete3 inside any function that uses it (the import has C dependencies and can fail to resolve on some platforms; the Biopython-based code paths are the always-available fallback).

## File formats

### Newick

The de facto standard for tree exchange. Felsenstein 1986. Parenthetical with optional branch lengths and node labels:

```
((human:0.05,chimp:0.05):0.10,mouse:0.20);
```

- Parentheses define subtrees.
- Commas separate siblings.
- `:length` is the branch length (optional).
- A leaf label is just a string (no parentheses).
- An internal node label can be a bootstrap value: `((A:0.1,B:0.2)95:0.3,C:0.4);` where `95` is the bootstrap percentage on the (A,B) clade.
- The semicolon terminates the tree.

The Newick spec is informal. Read `Bio.Phylo.NewickIO.Parser` for the canonical Biopython parser; the gory details are documented at <https://evolution.genetics.washington.edu/phylip/newick_doc.html>.

### Nexus

A more structured format with blocks for taxa, characters, and trees. Maddison, Swofford, Maddison 1997, *Systematic Biology* 46:590. Free full text:
<https://academic.oup.com/sysbio/article/46/4/590/1629654>

```
#NEXUS

BEGIN TAXA;
DIMENSIONS NTAX=3;
TAXLABELS human chimp mouse;
END;

BEGIN TREES;
TREE tree1 = ((human:0.05,chimp:0.05):0.10,mouse:0.20);
END;
```

Biopython reads and writes Nexus via `Bio.Phylo.write(tree, "tree.nex", "nexus")`. The Nexus block structure preserves annotations that bare Newick cannot.

### PhyloXML

The XML-based tree format. Han and Zmasek 2009. Useful when you want to attach metadata (species, sequence, GenBank accession) to each node. Biopython supports it via `Bio.Phylo.write(tree, "tree.xml", "phyloxml")`.

## Distance models — the math you need to know

### Jukes-Cantor (JC; Jukes and Cantor 1969)

The simplest substitution model: all four nucleotides are interchangeable with equal rate. The JC-corrected distance is:

```
d_JC = -3/4 * ln(1 - 4/3 * p)
```

where `p` is the observed proportion of differing columns (the p-distance). The correction undoes the saturation effect: at high `p`, the observed difference plateaus near 0.75 and the JC formula maps it back to the true expected number of substitutions per site.

JC breaks when `p >= 0.75` (the log argument goes non-positive). In practice we clamp `p` at 0.7499 and warn the user that the distance is saturated.

### Kimura two-parameter (K2P; Kimura 1980, *Journal of Molecular Evolution* 16:111)

K2P allows transitions (A<->G, C<->T) and transversions (everything else) to have different rates. The K2P distance is:

```
d_K2P = -1/2 * ln(1 - 2*P - Q) - 1/4 * ln(1 - 2*Q)
```

where `P` is the proportion of transition differences and `Q` is the proportion of transversion differences. Same saturation caveats apply.

For most vertebrate cytochrome b panels, K2P gives distances ~5-15% larger than JC for distant sequences and effectively the same as JC for close sequences.

### Why correct at all?

Imagine two sequences that each accumulated five substitutions at the same site over evolutionary time. If both reverted to the original base, the observed difference is zero, even though the true distance is ten substitutions per site. The JC and K2P corrections statistically recover the true expected distance from the observed difference, under the assumption that the substitution process is memoryless.

## Compute requirements

| Step | Demo (10 seqs, 1,140 cols) | Real panel (100 seqs, 5,000 cols) |
|------|----------------------------|----------------------------------------|
| MAFFT `--retree 2` | ~0.4 sec | ~5 sec |
| MAFFT `L-INS-i` | ~3 sec | ~60 sec |
| Distance matrix (JC) | ~10 ms | ~1 sec |
| NJ tree | ~5 ms | ~100 ms |
| IQ-TREE 2 ML + UFBoot 1000 | ~20 sec | ~10 min |
| `Bio.Phylo.draw` to PNG | ~200 ms | ~1 sec |
| ete3 rendering | ~500 ms | ~5 sec |

A laptop with 8 GB of RAM is enough for everything in Week 9.

## Datasets

All datasets used in Week 9 are free and publicly distributed.

### The cytochrome b panel (used by Exercise 1, 2, 3 and the mini-project)

Ten vertebrate cytochrome b nucleotide sequences (full CDS, ~1,140 bp each), drawn from NCBI RefSeq:

- *Homo sapiens* (human), *Pan troglodytes* (chimpanzee), *Mus musculus* (mouse), *Rattus norvegicus* (rat), *Bos taurus* (cow), *Sus scrofa* (pig), *Gallus gallus* (chicken), *Danio rerio* (zebrafish), *Xenopus laevis* (frog), *Ornithorhynchus anatinus* (platypus).

The file is bundled at `data/cytb_vertebrates.fasta` (~12 KB). Sequence accessions are recorded in `data/cytb_vertebrates.tsv`. The platypus is the conventional outgroup for the eutherian + bird + fish + frog subset.

To regenerate from RefSeq (Week 4 patterns):

```python
from Bio import Entrez, SeqIO
Entrez.email = "you@example.com"
handle = Entrez.efetch(db="nuccore", id="NC_012920.1", rettype="fasta_cds_na", retmode="text")
# pick the cytb record and write it out.
```

### The 16S rRNA panel (homework problem 4)

Twenty 16S rRNA sequences from cultured bacteria across the major phyla (Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes, etc.), drawn from SILVA / NCBI. ~1,500 bp each. File: `data/sixteen_s.fasta` (~30 KB).

### The cytochrome oxidase I (COI) panel (mini-project)

A hand-curated panel of forty COI sequences spanning all eleven major animal phyla (~650 bp each, the standard barcode region). File: `data/coi_barcode.fasta` (~30 KB). This is the input to the mini-project pipeline.

## Free-tier compute and storage

- **GitHub Codespaces** (free 60 hours/month): a 2-core / 8 GB instance is enough for all of Week 9. The codespace template in the curriculum repo has MAFFT, IQ-TREE 2, Biopython, and ete3 pre-installed.
- **Google Colab** (free tier): also enough; `apt install mafft` and `pip install biopython ete3 iqtree` get you running in a notebook in two minutes. IQ-TREE on Colab is reasonably fast even on small inputs.
- **Local laptop**: any 8 GB RAM laptop is enough.

## Style guide for this week

- **Pin tool versions.** Always write "MAFFT 7.526", "IQ-TREE 2.3.6", "RAxML-NG 1.2.2", "Biopython 1.84", "ete3 3.1.3", "trimAl 1.4.1" — never just "MAFFT" or "the latest IQ-TREE." Alignments and trees drift across versions; reproducibility requires pinning.
- **Pin the algorithm flag.** MAFFT `--auto` is not reproducible across versions; pick one of `--retree 2`, `--maxiterate 1000`, or `--localpair --maxiterate 1000` and write it down. Same for IQ-TREE 2's `-m MFP` (use `-m GTR+I+G` or the specific best-found model after running MFP once).
- **Pin the seed.** Bootstrap resampling, ML starting trees, and NJ tie-breaking are all randomized. Every Week 9 script accepts a `--seed` argument and records the integer in the run-info JSON. The default seed across the curriculum is `42` (we use the same integer in every week for cross-week reproducibility).
- **Cite tools by paper.** "MAFFT (Katoh and Standley 2013)" or "NJ (Saitou and Nei 1987)" — not just "MAFFT" or "NJ." This is standard practice in publications and what reviewers expect.
- **Report numbers, not adjectives.** "MAFFT aligned the ten cytochrome b sequences into a 1,140-column alignment in 0.4 seconds; after trimming columns with > 50% gaps, 1,098 columns remained; the JC pairwise distance ranged from 0.014 (human-chimp) to 0.314 (platypus-chicken)" is a sentence. "The alignment was good" is not.
- **Always print the database version, the run date, the seed, and the algorithm in the report header.** Without these you have a tree that nobody can reproduce.
- **Always pick an outgroup explicitly and root the tree on it.** "rooted by midpoint" is a fallback, not a default. For the demo cytochrome b panel the canonical outgroup is the platypus.

## Common questions

**Q. MAFFT and Clustal Omega produce different alignments of the same FASTA. Which is right?**

Neither is uniquely "right" — they are two implementations of the same NP-hard problem (the simultaneous alignment of N sequences) with different heuristics. Disagreements almost always trace to indel-rich regions where the two tools place gaps differently. For most downstream uses (tree-building, conserved-column extraction), the two alignments produce nearly identical results because the unambiguous columns are aligned the same way by both. Treat persistent disagreement on a specific region as a flag for manual inspection.

**Q. My NJ tree and my ML tree disagree on one internal branch. Which do I trust?**

The ML tree, because it is based on a probabilistic model and explicitly optimizes the likelihood under that model. NJ is a fast distance-based heuristic that does not explore the tree space. That said, if your ML tree disagrees with your NJ tree by more than a couple of branches, look at the alignment columns supporting each branch; the disagreement is usually a poorly-aligned region or a long-branch attraction artefact.

**Q. What does a 70% bootstrap mean?**

It means that the branch appeared in 70% of the bootstrap-resampled trees. It is a measure of *robustness to column resampling*, not a measure of how true the branch is in the species tree. A 70% bootstrap on a branch supported by twenty alignment columns is more informative than a 90% bootstrap on a branch supported by three; always check the underlying column count when interpreting low support.

**Q. Why do I need to trim alignment columns? Doesn't that lose information?**

Poorly-aligned columns (lots of gaps, no clear homology) introduce noise into the distance calculation that biases the tree. The trimming step (`trimAl -automated1`, or our simpler "drop columns with > 50% gaps" rule) trades a small amount of information for a much cleaner signal. For most panels, trimming changes the tree topology on at most one or two internal branches.

**Q. Should I use Jukes-Cantor or Kimura two-parameter for the demo cytochrome b panel?**

Either works. JC is the simpler model and the default we use in Exercise 2. K2P is slightly more accurate for vertebrate cytochrome b because the transition / transversion ratio in this gene is far from 1.0. For the mini-project we use K2P and explain why in the write-up.

**Q. Can I use AlphaFold structural alignments instead of MAFFT?**

That is the cutting edge as of 2024-2025 (see Foldseek, FATCAT, US-align, the DALI server). Structural alignments are more reliable on distantly related sequences where sequence identity drops below ~25%. For Week 9's vertebrate cytochrome b panel, sequence identity is high enough (~75-90%) that sequence-based alignment is appropriate and structural alignment would not improve the tree. We mention this in Lecture 1 as a stretch direction.

**Q. Do I need to commit the IQ-TREE intermediate files (`*.ckp.gz`, `*.bionj`)?**

No. Gitignore them. The `*.treefile` and the `*.iqtree` log are the artefacts to keep; the rest can be regenerated from the inputs plus the seed.

**Q. The ete3 install is failing on my Mac with a PyQt5 error. What do I do?**

The PyQt5 dependency is only needed for the interactive GUI. For headless rendering (PNG / SVG to file), the conda-forge `ete3` package on macOS works without PyQt5 if you set `QT_QPA_PLATFORM=offscreen` in the environment. Alternatively, fall back to `Bio.Phylo.draw` with matplotlib, which is what the always-available code paths in Exercise 3 use.
