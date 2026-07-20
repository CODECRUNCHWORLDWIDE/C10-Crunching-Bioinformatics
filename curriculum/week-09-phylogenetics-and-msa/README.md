# Week 9 — Phylogenetics and Multiple Sequence Alignment

> **Reproducibility note.** Phylogenetic inference is famously sensitive to input data, alignment choice, model choice, and the random seed. Two runs of the same pipeline on the same FASTA can produce trees that disagree on internal branches if the seed, the model, or the alignment column-trimming step changes. The Week 9 deliverables pin every parameter (MAFFT version, MAFFT algorithm, distance model, NJ implementation, RAxML / IQ-TREE version, seed integer, bootstrap replicate count) and emit them in a `run-info.json` that travels with every tree. Without that JSON, two trees that look identical may have been built from subtly different inputs and you will not be able to reconstruct which was which.

In Week 8 you took a per-variant VCF and produced an annotation report: VEP for the functional consequence, gnomAD for the population frequency, ClinVar for the clinical knowledge, and the mechanically computable subset of the ACMG criteria for the classification. The unit of analysis was one variant in one genome. Week 9 leaves single-genome interpretation behind and asks the comparative question: **given N homologous sequences from N organisms, how related are they, and what is the order in which they diverged?** The unit of analysis is the column of a multiple sequence alignment (MSA); the output is a tree.

The toolchain is older than next-generation sequencing — phylogenetic inference predates Sanger sequencing — and the literature is dense. The Week 9 framing is deliberately narrow. We pick one alignment program (**MAFFT**, Katoh and Standley 2013, *Molecular Biology and Evolution* 30:772), one fast tree builder (**neighbor-joining** via Biopython's `Bio.Phylo.TreeConstruction`, Saitou and Nei 1987, *Molecular Biology and Evolution* 4:406), one model-based tree builder for the challenge (**IQ-TREE 2**, Minh et al. 2020, *Molecular Biology and Evolution* 37:1530; or **RAxML-NG**, Kozlov et al. 2019, *Bioinformatics* 35:4453), and one tree-rendering library (Biopython's `Bio.Phylo` for the basic case; `ete3`, Huerta-Cepas et al. 2016, *Molecular Biology and Evolution* 33:1635 for the publication-quality case). All of those are free, open source, and pip- or conda-installable. The two MSA alternates we mention are **Clustal Omega** (Sievers et al. 2011, *Molecular Systems Biology* 7:539) and **MUSCLE 5** (Edgar 2022, *Nature Communications* 13:6968); both are free, both are good, and we explain when to prefer each.

By Friday of Week 9 you will be able to take a FASTA of homologous protein or nucleotide sequences (the canonical demo is a small panel of cytochrome b sequences from ten vertebrates), run a progressive MSA via a `subprocess` call to MAFFT, parse the alignment with Biopython, compute a pairwise distance matrix under either Jukes-Cantor (Jukes and Cantor 1969) or Kimura two-parameter (Kimura 1980, *Journal of Molecular Evolution* 16:111), build a neighbor-joining tree, optionally generate bootstrap-resampled trees to estimate branch support, write the result as a **Newick** string (Felsenstein 1986, the de facto standard), and render it as a PNG or SVG. The mini-project is the whole pipeline wrapped in a CLI script, with the run parameters pinned in a `run-info.json` and the tree colour-coded by clade.

The other half of the week is **why the tree might be wrong**, and the answer has many layers: the alignment can be wrong (poorly aligned columns drive spurious branches; this is why people trim with `trimAl` or Gblocks before tree-building), the distance model can be wrong (Jukes-Cantor assumes equal base frequencies and equal substitution rates, which is rarely true), the algorithm can be wrong (NJ is fast but does not explore the tree space; ML is slower but better-justified), the bootstrap support can be misleading (high support for the wrong branch is the failure mode), and the input data can be the wrong question (a gene tree is not a species tree; a single-gene tree from cytochrome b might disagree with a species tree from concatenated genome data — see the entire field of phylogenomics). We will name these limits in Lecture 2 and in the mini-project write-up; you will repeat them whenever you publish a tree.

## Learning objectives

By the end of this week, you will be able to:

- **Describe** the progressive alignment heuristic: pairwise distance, guide tree, profile-profile alignment along the guide tree. Name MAFFT, Clustal Omega, and MUSCLE 5 as the three canonical free progressive aligners and state when each is preferred (MAFFT for speed and quality on most inputs; Clustal Omega for very large inputs; MUSCLE 5 for the most-accurate-per-second category on small inputs).
- **Run** MAFFT from the command line and from Python via `subprocess.run(..., check=True, capture_output=True)`. Use the `--auto` algorithm-choice flag, pin a seed, and verify that a second run with the same inputs and seed produces a byte-identical alignment.
- **Parse** a FASTA alignment with `Bio.SeqIO` and `Bio.AlignIO`. Iterate alignment columns and per-record gap fractions. Trim columns with > 50% gap content before any distance computation.
- **Compute** a pairwise distance matrix under Jukes-Cantor (Jukes and Cantor 1969) and under Kimura two-parameter (Kimura 1980, *Journal of Molecular Evolution* 16:111). Explain why the raw Hamming distance is biased downward for distant sequences (multiple substitutions per site) and why the JC and K2P corrections matter once observed differences exceed ~10%.
- **Build** a neighbor-joining tree (Saitou and Nei 1987, *Molecular Biology and Evolution* 4:406) using `Bio.Phylo.TreeConstruction.DistanceTreeConstructor`. Optionally build a UPGMA tree (Sokal and Michener 1958, an older method) and explain why UPGMA is correct only under a strict molecular clock.
- **Build** a maximum-likelihood tree (Felsenstein 1981, *Journal of Molecular Evolution* 17:368, free at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>) using either **IQ-TREE 2** (Minh et al. 2020) or **RAxML-NG** (Kozlov et al. 2019). Pick a substitution model (GTR+I+G for nucleotides; LG+G for proteins) and a seed; record both in the run-info JSON.
- **Estimate** bootstrap support by resampling alignment columns with replacement, rebuilding the tree per replicate, and counting how often each internal branch appears (Felsenstein 1985, *Evolution* 39:783). Use 100-1000 replicates for the mini-project. Be honest: a "70% bootstrap" branch is supported, not proven, and high support is necessary but not sufficient for correctness.
- **Write and read** trees in Newick format (Felsenstein 1986, the parenthetical de facto standard) and Nexus format (Maddison, Swofford, Maddison 1997, *Systematic Biology* 46:590). Round-trip a tree through both formats and verify topology preservation.
- **Render** a tree with `Bio.Phylo.draw_ascii` for the terminal sanity check and with `Bio.Phylo.draw` or `ete3.TreeStyle` for a PNG / SVG. Colour internal branches by bootstrap support; label leaves by species; root on a named outgroup.
- **Identify and reason about** the standard failure modes: alignment-induced artefacts, long-branch attraction (Felsenstein 1978, *Systematic Biology* 27:401), incongruence between gene trees and species trees, the GIGO problem (low-coverage sequences pull the tree), and the seed-and-version-dependence problem that the run-info JSON is meant to defend against.

## Prerequisites

This week assumes Weeks 1-8 are **done and committed**. Specifically:

- You can read and write FASTA records from Week 2 (`SeqIO.parse`, `SeqIO.write`, the `SeqRecord` object).
- You can run a CLI tool from Python with `subprocess.run(..., check=True, capture_output=True)` from Week 5 Exercise 1 and Week 8 Exercise 1.
- You can serialize a `pandas.DataFrame` to CSV / TSV and an HTML page from Weeks 4, 7, and 8.
- You have Python 3.11+, Biopython 1.84+, numpy 1.26+, pandas 2.2+, requests 2.32+, plus the conda channels `bioconda` and `conda-forge` configured.
- You will install `mafft`, `iqtree`, optionally `raxml-ng`, and the `ete3` Python package this week. The canonical conda call: `conda install -c bioconda mafft=7.526 iqtree=2.3.6 raxml-ng=1.2.2 ete3=3.1.3 biopython=1.84`. None of these tools has a paywall; all are open source.

You do not need biology beyond "homologous sequences are sequences that share a common ancestor; an alignment column hypothesizes that all the residues in the column descend from the same ancestral residue; a substitution is a residue swap at one column over evolutionary time; a tree groups sequences by how recently they shared a common ancestor; an outgroup is a sequence known to be more distantly related to the rest, used to root the tree." Lecture 1 fixes the vocabulary.

## Topics covered

- **The progressive alignment heuristic.** Pairwise distance, guide tree, profile-profile alignment. Why it is a heuristic: the order of insertions and deletions is fixed early and never revisited. MAFFT's iterative-refinement modes (`--maxiterate 1000`) walk back some of this; Clustal Omega's HMM profile mode walks back differently; MUSCLE 5's super5 algorithm walks back differently again. All three improve on naive progressive in the same direction.
- **MAFFT 7** (Katoh and Standley 2013, *Molecular Biology and Evolution* 30:772). The default modern aligner for most inputs in the 10-10,000 sequence range. Free, open source, written in C with a Perl front-end. Algorithms: `--auto` chooses by input size (FFT-NS-2 for small; FFT-NS-i for medium-quality iterative refinement; L-INS-i for highest-accuracy local alignment on small inputs). Tool docs: <https://mafft.cbrc.jp/alignment/software/>.
- **Clustal Omega** (Sievers et al. 2011, *Molecular Systems Biology* 7:539). The successor to ClustalW. Designed for very large inputs (up to ~190,000 sequences). Slower per sequence than MAFFT on small inputs but scales better. Free. Tool docs: <http://www.clustal.org/omega/>.
- **MUSCLE 5** (Edgar 2022, *Nature Communications* 13:6968). The 2022 rewrite of MUSCLE with the super5 algorithm. Highest accuracy per second on small-to-medium inputs in published benchmarks. Free. Tool docs: <https://drive5.com/muscle/> and <https://github.com/rcedgar/muscle>.
- **Distance models for nucleotides.** Raw Hamming (p-distance) is the fraction of differing columns. JC (Jukes-Cantor 1969, the equal-rates one-parameter correction) is the basic model. K2P (Kimura 1980) adds a transition / transversion ratio. F81 (Felsenstein 1981) and GTR (Tavare 1986) further loosen the equal-frequency assumption. The Week 9 exercises stick to JC and K2P; the challenge moves to GTR via IQ-TREE.
- **Distance models for proteins.** PAM (Dayhoff et al. 1978), Blosum (Henikoff and Henikoff 1992, *PNAS* 89:10915), JTT (Jones, Taylor, Thornton 1992), WAG (Whelan and Goldman 2001), LG (Le and Gascuel 2008, *Molecular Biology and Evolution* 25:1307). LG is the modern default; we use it for the challenge.
- **Tree-building methods.** UPGMA (the unweighted pair-group method with arithmetic mean) is the simplest; assumes a strict molecular clock and is generally not what you want for sequences that evolved at different rates. Neighbor-joining (Saitou and Nei 1987) is fast, gives unrooted trees, and is the standard fast method. Maximum-likelihood (Felsenstein 1981) is slower but model-based and is the standard publishable method; modern fast implementations (IQ-TREE 2, RAxML-NG) make it tractable for thousands of sequences on a laptop. Bayesian inference (MrBayes, BEAST) is the slowest and most sophisticated; out of scope for Week 9.
- **Bootstrap support** (Felsenstein 1985, *Evolution* 39:783). Resample alignment columns with replacement, rebuild the tree per replicate, count branch frequency. The published numbers next to internal branches in a tree figure are bootstrap percentages. 100 replicates is the bare minimum; 500-1000 is typical; ultrafast bootstrap (UFBoot, Hoang et al. 2018, *Molecular Biology and Evolution* 35:518) in IQ-TREE 2 is the fast modern variant.
- **Tree file formats.** Newick (Felsenstein 1986; the parenthetical de facto standard) is `((A:0.1,B:0.2):0.3,C:0.4);`. Nexus (Maddison, Swofford, Maddison 1997, *Systematic Biology* 46:590) embeds Newick inside a structured block with taxa, characters, and tree definitions; preferred for trees with annotations. PhyloXML (Han and Zmasek 2009) is the XML alternative. Biopython's `Bio.Phylo` reads and writes all three.
- **Rendering trees.** `Bio.Phylo.draw_ascii` for the terminal sanity check. `Bio.Phylo.draw` (matplotlib-backed) for a basic PNG. `ete3.TreeStyle` for the publication-quality figure with branch colouring, leaf images, and node decorations. `iTOL` (Letunic and Bork 2024, *Nucleic Acids Research* 52:W78) is the web-based alternative; free for personal use.
- **Reproducibility hygiene.** Pin tool versions (MAFFT 7.526, IQ-TREE 2.3.6, RAxML-NG 1.2.2, Biopython 1.84, ete3 3.1.3). Pin the algorithm (MAFFT `--auto` is *not* reproducible across versions; pick one — `--retree 2`, `--maxiterate 1000`, or `L-INS-i` — and write it down). Pin the seed for any randomized step (bootstrap replicate generation; ML search starting tree; NJ tie-breaking). Pin the trim threshold (we use "drop columns with > 50% gaps" throughout). Without these, two runs on the same FASTA can produce different trees and the difference is impossible to debug.
- **Standard failure modes.** Alignment-induced artefacts (a sloppy alignment produces a sloppy tree). Long-branch attraction (Felsenstein 1978; rapidly evolving sequences cluster together regardless of true relatedness). Gene-tree / species-tree incongruence (incomplete lineage sorting, horizontal gene transfer, hybridization; a single-gene tree is not the species tree). Saturation (when the observed difference plateaus near the JC asymptote; the K2P correction can no longer recover the true distance). GIGO (low-coverage or contaminated sequences pull the tree toward themselves). The mini-project write-up names which of these its tree could plausibly suffer from.

## Weekly schedule

The schedule below adds up to approximately **33 hours**. Treat it as a target. Monday's MSA lecture is the hour that decides whether the rest of the week makes interpretive sense — the tree is downstream of the alignment, and a bad alignment makes a bad tree no matter how good the tree builder is.

| Day       | Focus                                                | Lectures | Exercises | Challenges | Quiz/Read | Homework | Mini-Project | Self-Study | Daily Total |
|-----------|------------------------------------------------------|---------:|----------:|-----------:|----------:|---------:|-------------:|-----------:|------------:|
| Monday    | Progressive alignment, MAFFT, Clustal Omega, MUSCLE  |    2h    |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5.5h    |
| Tuesday   | Distance models (JC, K2P), NJ, UPGMA                 |    1.5h  |    1.5h   |     0h     |    0.5h   |   1h     |     0h       |    0.5h    |     5h      |
| Wednesday | Maximum likelihood, bootstrap, IQ-TREE / RAxML       |    1.5h  |    1.5h   |     1h     |    0.5h   |   1h     |     1h       |    0.5h    |     7h      |
| Thursday  | Newick / Nexus I/O, rendering, the mini-project      |    1h    |    2h     |     1h     |    0.5h   |   1h     |     2h       |    0.5h    |     8h      |
| Friday    | Mini-project deep work + write-up                    |    0h    |    1h     |     0h     |    0.5h   |   1h     |     2h       |    0h      |     4.5h    |
| Saturday  | Mini-project deep work                               |    0h    |    0h     |     0h     |    0h     |   1h     |     2h       |    0h      |     3h      |
| Sunday    | Quiz, review, polish                                 |    0h    |    0h     |     0h     |    0.5h   |   0h     |     0h       |    0h      |     0.5h    |
| **Total** |                                                      | **6h**   | **7.5h**  | **2h**     | **3h**    | **6h**   | **7h**       | **2h**     | **33.5h**   |

## How to navigate this week

| File | What's inside |
|------|---------------|
| [README.md](./README.md) | This overview (you are here) |
| [resources.md](./resources.md) | MAFFT, Clustal Omega, MUSCLE 5, IQ-TREE, RAxML-NG, Biopython, ete3 docs and reference papers |
| [lecture-notes/01-msa-and-progressive-alignment.md](./lecture-notes/01-msa-and-progressive-alignment.md) | The progressive alignment heuristic; MAFFT, Clustal Omega, MUSCLE 5 in one lecture |
| [lecture-notes/02-distance-models-and-tree-building.md](./lecture-notes/02-distance-models-and-tree-building.md) | Jukes-Cantor, Kimura two-parameter, UPGMA, neighbor-joining, maximum likelihood, bootstrap |
| [lecture-notes/03-newick-nexus-and-rendering.md](./lecture-notes/03-newick-nexus-and-rendering.md) | Newick and Nexus file formats; reading, writing, rendering with Biopython and ete3 |
| [exercises/exercise-01-mafft-via-subprocess.py](./exercises/exercise-01-mafft-via-subprocess.py) | Run MAFFT on a small FASTA, parse the alignment, report column statistics |
| [exercises/exercise-02-distance-and-nj-tree.py](./exercises/exercise-02-distance-and-nj-tree.py) | Compute a JC and K2P distance matrix, build an NJ tree with Biopython, write Newick |
| [exercises/exercise-03-render-tree-with-bootstrap.py](./exercises/exercise-03-render-tree-with-bootstrap.py) | Generate bootstrap replicates, annotate branches with support, render PNG / SVG |
| [exercises/SOLUTIONS.md](./exercises/SOLUTIONS.md) | Worked solutions for all three exercises |
| [challenges/challenge-01-iqtree-ml-with-bootstrap.md](./challenges/challenge-01-iqtree-ml-with-bootstrap.md) | Switch from NJ to IQ-TREE 2 maximum likelihood with ultrafast bootstrap |
| [challenges/challenge-02-newick-nexus-roundtrip.md](./challenges/challenge-02-newick-nexus-roundtrip.md) | Round-trip a tree through Newick and Nexus; prove topology preservation |
| [quiz.md](./quiz.md) | 10 multiple-choice questions on MSA, distance models, tree-building, and reproducibility |
| [homework.md](./homework.md) | Six practice problems for the week |
| [mini-project/README.md](./mini-project/README.md) | End-to-end phylogenetic pipeline: FASTA in, annotated tree out, with full provenance |

## A note on tone

C10 is written in **lab-notebook voice**. We pin versions ("MAFFT 7.526", "IQ-TREE 2.3.6", "RAxML-NG 1.2.2", "Biopython 1.84", "ete3 3.1.3"). We cite tools by their paper ("MAFFT, Katoh and Standley 2013, *Mol Biol Evol* 30:772; NJ, Saitou and Nei 1987, *Mol Biol Evol* 4:406; ML, Felsenstein 1981, *J Mol Evol* 17:368"). We say "MAFFT aligned the ten cytochrome b sequences into a 1,140-column alignment in 0.4 seconds; after trimming columns with more than 50% gaps, 1,098 columns remained; the JC pairwise distance matrix had a maximum of 0.314 between the platypus and the chicken records" rather than "MAFFT aligned the sequences." A reproducible pipeline reports the version of every tool, the algorithm flag, the seed, the date of the run, and the alignment column count. The numbers are the work.

## A note on the data size

A ten-sequence cytochrome b FASTA is ~12 KB. A 100-sequence ribosomal-protein FASTA is ~100 KB. A 1,000-sequence panel is ~1 MB. None of the Week 9 inputs requires more than a laptop with 8 GB of RAM. The IQ-TREE 2 ML search for the challenge takes ~10-30 seconds on the demo cytochrome b panel; the ultrafast bootstrap with 1,000 replicates adds another ~30 seconds. The full mini-project pipeline runs in under two minutes on a 2020-era laptop.

Commit the input FASTA, the aligned FASTA, the trimmed alignment, the distance matrix, the NJ tree (Newick), any ML tree (Newick + IQ-TREE log), the rendered PNG / SVG, and the `run-info.json`. Gitignore the IQ-TREE intermediate files (`*.ckp.gz`, `*.bionj`, `*.mldist`) — they are large and reproducible from the inputs and the seed.

## Stretch goals

If you finish early and want to push further, try any of the following:

- Read Felsenstein 1981 (the maximum-likelihood paper) end to end at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>. It is ~20 pages and is the foundation of every modern phylogenetic inference paper.
- Add `trimAl` (Capella-Gutierrez et al. 2009, *Bioinformatics* 25:1972) for principled alignment trimming with `-automated1`. Compare your NJ tree before and after trimming. The branch topology usually moves on at least one internal node.
- Reproduce a small published phylogeny end to end. The COVID-19 spike-protein trees on Nextstrain are a good target; the FASTAs are public, the tree-building parameters are documented, and your tree should match the published one to within a handful of branch swaps.
- Add a Bayesian alternative: run MrBayes on the same alignment for 100,000 generations and compare the consensus tree topology to your ML tree. MrBayes is free but the install is fussy.
- Compute a Robinson-Foulds distance (Robinson and Foulds 1981, *Mathematical Biosciences* 53:131) between your NJ tree and your ML tree. Biopython does not implement RF; `ete3` does (`tree1.compare(tree2, unrooted=True)`). An RF of 0 means topologically identical; > 0 means at least one internal branch is in disagreement. Most ten-leaf demos produce an RF of 0-4.
- Try MUSCLE 5 (Edgar 2022) on the same FASTA and compare the alignment column count and the resulting tree. The two aligners almost always agree on most columns but disagree on indel-rich regions.

## Up next

Continue to [Week 10 — Comparative and Pan-Genomics](../week-10/) once you have pushed your mini-project to GitHub. Week 10 keeps the comparative thread but moves from single-gene trees to genome-scale comparisons: ortholog finding (OrthoFinder), synteny (MCScanX), and the pan-genome (Roary, PPanGGOLiN). The MAFFT and Newick patterns from Week 9 stay parked; the iteration over many genomes is the new work.

---

*If you find errors in this material, please open an issue or send a PR. Future learners will thank you.*
