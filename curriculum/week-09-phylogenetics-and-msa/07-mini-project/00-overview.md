# Mini-Project — End-to-end phylogenetics pipeline

> **Reproducibility note.** This mini-project produces a phylogenetic tree from a FASTA input. The tree is reproducible only if you ship the inputs, the parameters, the seed, and the tool versions alongside the Newick file. Without the `run-info.json` the tree is an opinion; with it, the tree is a reproducible result. **The output of this pipeline must travel with its run-info, every time.**

Build a reproducible phylogenetic pipeline that takes a FASTA of homologous sequences (the bundled `data/coi_barcode.fasta`, 40 cytochrome oxidase I sequences spanning eleven animal phyla), aligns them with MAFFT, trims gappy columns, computes a K2P distance matrix, builds a neighbor-joining tree, generates bootstrap support by column resampling, optionally re-runs the tree under maximum likelihood with IQ-TREE 2, renders the tree as a publication-quality PNG with branch colours by support tier, and emits a `run-info.json` recording every parameter. End with: the aligned FASTA, the trimmed FASTA, the Newick tree, the rendered PNG, the run-info JSON, and a 600-800 word write-up that defends every parameter and names every limit.

This is the C10 mini-project that produces a **methods-section-quality phylogenetic figure with measured provenance**, not just a one-shot demonstration. By the end of it you will have a `phylogenetics_pipeline.py` script and a `run.sh` wrapper, a results directory with the tree and the run-info, and a write-up that defends every parameter and explicitly names what the tree cannot reliably claim.

**Estimated time:** 7 hours (split across Wednesday, Thursday, Friday, Saturday in the suggested schedule).

---

## What you will produce

In your existing portfolio repo (`crunch-bio-portfolio-<yourhandle>`), add a new `week-09/mini-project/` directory:

```
crunch-bio-portfolio-<yourhandle>/
├── README.md                       (updated, with a Week 9 section)
└── week-09/
    └── mini-project/
        ├── README.md               one-page report (~600-800 words)
        ├── run.sh                  one-command reproduction script
        ├── env.yml                 conda environment file pinning all tool versions
        ├── data/
        │   ├── coi_barcode.fasta   input FASTA (40 sequences across 11 phyla)
        │   └── coi_metadata.tsv    per-taxon metadata (phylum, common name, accession)
        ├── phylogenetics_pipeline.py    the orchestration script
        ├── starter.py                   skeleton implementation with TODOs
        └── results/
            ├── aligned.fasta            MAFFT FFT-NS-2 output
            ├── trimmed.fasta            after the > 50% gap trim
            ├── distance_matrix_k2p.tsv  K2P pairwise distances
            ├── tree_nj.nwk              NJ tree in Newick
            ├── tree_nj_bootstrap.nwk    NJ tree with bootstrap percentages
            ├── tree_ml.nwk              optional: IQ-TREE 2 ML tree
            ├── tree_final.png           rendered tree (NJ + bootstrap, coloured)
            ├── tree_final.nex           Nexus version of the final tree
            └── run-info.json            run provenance
```

By the end you will have a clean, reproducible Week 9 directory you can point a recruiter at — and `phylogenetics_pipeline.py` is the kind of pipeline that opens conversations with working bioinformaticians and biotech / academic shops, *as long as* you can speak to its limits.

---

## The dataset

You will work with the bundled `data/coi_barcode.fasta` — a hand-curated panel of 40 cytochrome oxidase subunit I (COI) sequences spanning the eleven major animal phyla:

- *Homo sapiens*, *Mus musculus*, *Bos taurus* (Chordata, vertebrate).
- *Drosophila melanogaster*, *Apis mellifera*, *Bombyx mori* (Arthropoda, insect).
- *Daphnia pulex*, *Penaeus monodon* (Arthropoda, crustacean).
- *Caenorhabditis elegans* (Nematoda).
- *Lumbricus terrestris* (Annelida).
- *Helix pomatia*, *Mytilus edulis*, *Loligo vulgaris* (Mollusca).
- *Strongylocentrotus purpuratus* (Echinodermata).
- *Hydra vulgaris*, *Aurelia aurita* (Cnidaria).
- *Amphimedon queenslandica* (Porifera).
- *Schistosoma mansoni*, *Taenia solium* (Platyhelminthes).
- ... and so on, ~40 total, sampled to cover the major branches of the animal tree.

The COI region is ~650 bp; the FASTA is ~30 KB. COI is the canonical DNA barcoding region (Hebert et al. 2003, *Proceedings of the Royal Society B* 270:313) and the sequences are publicly available from BOLD (the Barcode of Life Data System; <https://www.boldsystems.org/>) and NCBI RefSeq.

Expected coverage by the pipeline:

- All 40 sequences should align into a ~650-column alignment (insertion-rich regions of COI are rare; the trim is small).
- The NJ tree should reproduce the conventional metazoan topology: a Bilateria clade (everything but Cnidaria, Porifera), the Bilateria split into Deuterostomia (chordates + echinoderms) and Protostomia (everything else), the Protostomia further split into Ecdysozoa (arthropods + nematodes) and Lophotrochozoa (annelids + molluscs + platyhelminths).
- Bootstrap support should be high (>= 90) on the well-established deep clades (Bilateria, Deuterostomia, Ecdysozoa) and moderate (60-90) on the Lophotrochozoa internal structure.
- The ML tree (optional) should agree with the NJ tree on the deep clades and may disagree on one or two shallow nodes; this is expected and is what the discussion section addresses.

### Why this composition

COI is fast-evolving (mitochondrial, no introns, third-codon-position-saturated in deep splits). A 40-sequence COI panel across eleven phyla has plausible signal for the well-established deep clades and ambiguous signal for the deepest splits (Cnidaria vs Bilateria; the position of Porifera). The pipeline produces a tree, the bootstrap reveals which branches are well-supported, and the discussion section explains *which conclusions are robust* and *which require longer alignments or whole-genome data*.

---

## Rules

- **You may** use MAFFT 7.526, IQ-TREE 2.3.6 (optional), Biopython 1.84, ete3 3.1.3 (optional), matplotlib 3.8, and the standard library.
- **You may** consult Lectures 1, 2, 3, the MAFFT paper (Katoh and Standley 2013), the NJ paper (Saitou and Nei 1987), the ML paper (Felsenstein 1981 at <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197550/>), the bootstrap paper (Felsenstein 1985), the IQ-TREE 2 paper (Minh et al. 2020), and the Week 9 exercises and challenges.
- **You may NOT** copy a pre-written phylogenetics pipeline from the internet. The point is to *build* the pipeline. Reading the Biopython tutorial chapter 13 for inspiration is fine; copy-pasting a complete phylogenetics pipeline is not.
- **You must** pin every tool version, every algorithm flag, the trim threshold, the distance method, the bootstrap replicate count, and the seed in the `run-info.json`.
- **You must** root the tree on a named outgroup (the conventional choice for the COI animal panel is the Cnidaria + Porifera outgroup; root on Porifera if a single taxon is required).
- **You must** ladderize the tree before rendering so the leaf order is deterministic.
- **You should** ship both Newick and Nexus versions of the final tree.
- **You must** commit the input FASTA, the alignment, the trimmed alignment, the tree files, the rendered figure, and the `run-info.json`. Gitignore the IQ-TREE 2 intermediate files.

---

## Acceptance criteria

- [ ] `mini-project/phylogenetics_pipeline.py` exports a function `build_tree(fasta_path: Path, out_dir: Path, outgroup: str, seed: int = 42) -> Path` that runs the full pipeline and returns the path to the final Newick.
- [ ] The pipeline implements **seven stages**:
  1. **Validate input.** Confirm the input FASTA exists, has at least 5 records, no duplicate record IDs, and no all-N sequences.
  2. **Align with MAFFT.** Run `mafft --retree 2 --maxiterate 0 --nuc --anysymbol --quiet --thread 4 input.fasta > aligned.fasta`. Skip if the output is newer than the input.
  3. **Trim columns.** Drop columns where the gap fraction exceeds 0.5.
  4. **Compute distance matrix.** K2P distances under pairwise deletion.
  5. **Build NJ tree.** Biopython `DistanceTreeConstructor.nj`; root on the named outgroup; ladderize.
  6. **Bootstrap.** 500 replicates with the pinned seed; bipartition support attached to internal nodes.
  7. **Render.** PNG via Bio.Phylo + matplotlib (always) and ete3 (if available); plus Newick and Nexus.
- [ ] `phylogenetics_pipeline.py` produces:
  - `results/aligned.fasta` — MAFFT output.
  - `results/trimmed.fasta` — after column trim.
  - `results/distance_matrix_k2p.tsv` — square K2P matrix.
  - `results/tree_nj_bootstrap.nwk` — NJ Newick with bootstrap percentages.
  - `results/tree_final.nex` — Nexus version.
  - `results/tree_final.png` — rendered figure.
  - `results/run-info.json` — versions, dates, parameters, and the platypus / Porifera outgroup name.
- [ ] `mini-project/README.md` is a one-page (~600-800 word) report containing:
  - One-sentence description of the dataset, the alignment tool, the tree-building method, and the bootstrap protocol.
  - Methods section in C10 voice: every tool pinned ("MAFFT 7.526", "Biopython 1.84", "IQ-TREE 2.3.6"), every parameter explicit, the seed stated.
  - Results section in C10 voice: alignment column counts (before and after trim), distance matrix extremes (smallest and largest K2P), the NJ tree topology summary (which major clades are recovered), the bootstrap support distribution (median, min, max), the deepest internal branches and their support.
  - Discussion section: 150-250 words on the limits of the pipeline. What is *not* in the tree? What does the bootstrap *not* tell you? Which standard failure modes (alignment artefacts, long-branch attraction, gene-tree-vs-species-tree, saturation) is the panel plausibly susceptible to?
- [ ] `run.sh` is a single bash script that, given a fresh checkout + `conda env create -f env.yml`, reproduces the entire pipeline from scratch in under three minutes.
- [ ] The repo is **public** and at least one classmate or instructor has been added as a collaborator.
- [ ] **Most importantly**: the run-info JSON is complete and the discussion section is honest about the limits.

---

## Suggested approach (rough timeline)

### Wednesday (1 hour)

1. (15 min) `git clone`, set up `mini-project/` directory.
2. (30 min) Write `env.yml` with pinned tool versions; create conda env; verify each tool is on the PATH (MAFFT, Biopython, optionally IQ-TREE 2 and ete3).
3. (15 min) Read this README end to end and the Challenge 1 README. Sketch the seven-stage flow on paper.

### Thursday (2 hours)

1. (45 min) Implement Stages 1-3 (validate, MAFFT, trim) in `phylogenetics_pipeline.py`. Reuse the Exercise 1 helpers; they cover most of this.
2. (45 min) Implement Stages 4-5 (distance matrix, NJ tree, root, ladderize). Reuse the Exercise 2 helpers.
3. (30 min) Implement Stage 6 (bootstrap). Reuse the Exercise 3 helpers.

### Friday (2 hours)

1. (45 min) Implement Stage 7 (rendering) and the run-info JSON writer.
2. (30 min) Write `run.sh` and `env.yml`. Do a full `bash run.sh` on a fresh `conda env create -f env.yml` to verify end-to-end reproducibility.
3. (45 min) Draft the README.md write-up (methods + results + 250-word discussion).

### Saturday (2 hours)

1. (30 min) Optional: add the IQ-TREE 2 ML branch. Compare the ML topology to the NJ topology; report disagreements.
2. (30 min) Polish the rendered figure (clade colours by phylum; branch colours by bootstrap tier; a scale bar).
3. (30 min) Final edit of the write-up.
4. (30 min) Push, add a classmate as a collaborator, write the commit messages with specific numbers.

---

## Tips for the write-up

- **Lead with numbers.** "MAFFT FFT-NS-2 aligned the 40 COI sequences into a 657-column alignment in 1.2 seconds; the > 50% gap trim removed 12 columns, leaving 645. The K2P matrix had a minimum of 0.018 (D. melanogaster vs A. mellifera, both arthropod insects) and a maximum of 0.412 (Amphimedon queenslandica, a sponge, vs Strongylocentrotus purpuratus, an echinoderm). The NJ tree recovered the Bilateria clade with 98% bootstrap support, Deuterostomia with 86%, Ecdysozoa with 93%, and Lophotrochozoa with 64%."
- **Defend every parameter.** "We used MAFFT `--retree 2 --maxiterate 0` rather than L-INS-i because the demo panel is large enough (40 sequences) that L-INS-i's accuracy gain is small and its run-time cost is significant; we used K2P rather than JC because the transition / transversion ratio in COI is far from 1.0; we used 500 bootstrap replicates because the recommended floor is 500 for inferential use; the random seed is 42 by curriculum convention."
- **Name two failure modes the panel is susceptible to.** "The COI region is third-codon-saturated in the deepest splits (Bilateria vs Cnidaria) and the K2P correction is unreliable beyond p > 0.5; the inclusion of Amphimedon queenslandica (a sponge with a long terminal branch) is a plausible long-branch attraction concern."
- **Acknowledge what is missing.** "A whole-genome tree would resolve the Lophotrochozoa internal structure better than COI alone; a multi-gene concatenated alignment is the conventional next step."

---

## Stretch goals (optional)

- Replace the NJ tree with an IQ-TREE 2 ML tree under `-m MFP -B 1000 -alrt 1000 -seed 42`. Compare the topologies. The ML tree is more defensible for the deep splits.
- Add `trimAl -automated1` instead of the simple 50% gap rule. Compare the trimmed column counts and the resulting tree topology.
- Compute a Robinson-Foulds distance (Robinson and Foulds 1981) between the NJ tree and the ML tree using `ete3`'s `tree.compare(other, unrooted=True)`. Report the RF value in the write-up.
- Add a publication-quality figure: clade backgrounds coloured by phylum, leaf labels with the common name, a scale bar, the outgroup labelled, the deepest internal branch annotated with its bootstrap support.
- Reproduce a small published metazoan tree (e.g. the Dunn et al. 2008 *Nature* 452:745 broad metazoan tree) on a similar sequence set and report agreement / disagreement.

---

## What to commit

By the end of the mini-project your `week-09/mini-project/` should contain:

```
mini-project/
    README.md
    env.yml
    run.sh
    phylogenetics_pipeline.py
    starter.py
    data/{coi_barcode.fasta, coi_metadata.tsv}
    results/{aligned.fasta, trimmed.fasta, distance_matrix_k2p.tsv,
             tree_nj_bootstrap.nwk, tree_final.nex, tree_final.png,
             run-info.json}
```

Gitignore the IQ-TREE 2 intermediate files (`*.ckp.gz`, `*.bionj`, `*.mldist`) and the matplotlib font cache. The commit message for the final pipeline run should be specific, e.g. `mini-project: COI metazoan tree (40 seqs, 11 phyla, K2P NJ, 500 UFBoot), Bilateria 98% support, deep splits resolved`.
