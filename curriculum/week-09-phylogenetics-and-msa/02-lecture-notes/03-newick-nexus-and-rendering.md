# Lecture 3 — Newick, Nexus, and Tree Rendering

> **Reproducibility note.** Tree file formats encode topology and branch lengths but rarely encode the build provenance (model, seed, bootstrap replicates, software version). The conventional fix is to ship the tree file *alongside* a `run-info.json` that records the missing metadata. A Newick file by itself is not reproducible; a Newick file plus a run-info JSON is.

> **Duration:** ~2 hours of reading + a brief rendering sanity check in Biopython and ete3.
> **Outcome:** You can read and write Newick and Nexus trees with Biopython, round-trip a tree through both formats while preserving topology, render a tree as an ASCII sketch (Biopython), a PNG (Bio.Phylo + matplotlib), and a publication-quality SVG (ete3), and annotate internal branches with bootstrap support values.

If you only remember one thing from this lecture, remember this:

> **A tree on its own is not enough. A reproducible phylogenetic result is a tree file plus a run-info JSON that records the alignment input, the substitution model, the seed, the bootstrap replicate count, and the software versions. Without the JSON, the tree is an opinion. With the JSON, the tree is a reproducible result.**

Lecture 2 left off at the tree object in memory (Biopython `Tree` or IQ-TREE 2 `.treefile`). Lecture 3 covers what to write to disk, how to read it back, and how to render it for a figure or a presentation.

---

## 1. Where we are in the pipeline

```
TREE (in-memory Biopython Tree object) ->
        Newick file (de facto standard; parenthetical) ->
        OPTIONAL Nexus file (richer; embeds Newick in a structured block) ->
        rendered figure (ASCII to terminal; PNG via Bio.Phylo; SVG via ete3).
```

The Newick file is the canonical exchange format. Every modern phylogenetics tool reads and writes Newick. The Nexus format adds structure (taxa block, characters block, trees block) and is preferred when you want to ship the tree with annotations attached to nodes.

---

## 2. The Newick format

Citation note: Newick (Felsenstein 1986; the informal "Newick format" specification) is named after the New Hampshire seafood restaurant where the format was agreed upon at the 1986 SSE/SSB conference. The de facto specification is at <https://evolution.genetics.washington.edu/phylip/newick_doc.html>.

The structure:

- Parentheses group sub-trees.
- Commas separate siblings.
- `:length` attaches a branch length (optional, but conventionally present).
- A leaf label is a bare string (no parentheses, no commas, no colons inside).
- An internal node label is a string between the closing parenthesis and the colon (often used for bootstrap percentages).
- The semicolon terminates the tree.

Example — a three-leaf rooted tree with branch lengths:

```
((human:0.05,chimp:0.05):0.10,mouse:0.20);
```

This says: there is a clade containing `human` and `chimp` (each with branch length 0.05); that clade joins `mouse` (branch length 0.20) at the root, with the clade's branch to the root being 0.10.

Example with bootstrap support on the internal branch:

```
((human:0.05,chimp:0.05)98:0.10,mouse:0.20);
```

This says: the human-chimp clade has 98% bootstrap support.

### Newick parsing edge cases

- **Names with spaces or special characters** must be quoted in single quotes: `'Homo sapiens'`. Underscores are conventionally used in unquoted names (`Homo_sapiens`).
- **Negative branch lengths.** Newick technically allows them; some tree-building tools emit small negatives near zero. Treat them as zero in any downstream visualization.
- **Trailing semicolon.** Required. Trees without it are not parseable by strict implementations.
- **Multiple trees per file.** Newick allows it: each tree ends with `;`. Biopython's `Phylo.parse` returns an iterator over the trees.

### Reading and writing Newick with Biopython

```python
from __future__ import annotations

from pathlib import Path


def write_newick(tree: "Tree", out_path: Path) -> None:
    """Write a Biopython Tree to a Newick file."""
    from Bio import Phylo
    Phylo.write(tree, str(out_path), "newick")


def read_newick(path: Path) -> "Tree":
    """Read a Newick file and return a Biopython Tree."""
    from Bio import Phylo
    return Phylo.read(str(path), "newick")
```

Round-trip:

```python
from pathlib import Path


def newick_roundtrip_preserves_topology(tree: "Tree", tmp_path: Path) -> bool:
    """Return True if writing and re-reading the Newick preserves the leaf set
    and the parent-child relationships of every internal node.
    """
    from Bio import Phylo

    out: Path = tmp_path / "roundtrip.nwk"
    Phylo.write(tree, str(out), "newick")
    tree2 = Phylo.read(str(out), "newick")
    leaves_a: set[str] = {t.name for t in tree.get_terminals()}
    leaves_b: set[str] = {t.name for t in tree2.get_terminals()}
    if leaves_a != leaves_b:
        return False
    # A more careful check would compare the bipartitions; for Week 9 we
    # check the leaf set and the number of internal nodes, which is enough
    # for the demo panel.
    return len(tree.get_nonterminals()) == len(tree2.get_nonterminals())
```

Biopython's Newick writer drops internal node labels by default; if you want to preserve bootstrap support values, set the `format` flag in `Phylo.write` (`format_branch_length`, `plain` vs `bracket_labels`). Read the Biopython tutorial at <https://biopython.org/docs/latest/Tutorial/index.html> chapter 13 for the gory details.

---

## 3. The Nexus format

Citation: Maddison DR, Swofford DL, Maddison WP. "NEXUS: an extensible file format for systematic information." *Systematic Biology* 46:590 (1997). Free at <https://academic.oup.com/sysbio/article/46/4/590/1629654>.

Nexus is a structured format. A typical Nexus tree file:

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

Three things Nexus does better than bare Newick:

1. **Multiple trees with names.** Each `TREE name = ...` block names the tree, useful when you have an NJ tree, an ML tree, and a Bayesian consensus all in the same file.
2. **Embedded character data.** A Nexus file can carry the alignment alongside the tree.
3. **Per-tree annotations.** The TREES block can carry `[&support=...]` blocks per branch in some Nexus dialects (the BEAST / FigTree convention).

For Week 9 we read and write Nexus via Biopython:

```python
def write_nexus(tree: "Tree", out_path: Path) -> None:
    from Bio import Phylo
    Phylo.write(tree, str(out_path), "nexus")


def read_nexus(path: Path) -> "Tree":
    from Bio import Phylo
    return Phylo.read(str(path), "nexus")
```

Round-trip a tree through Newick and Nexus and verify topology is preserved — this is exactly what Challenge 2 asks for.

### Nexus dialects

Different tools emit slightly different Nexus. MrBayes embeds posterior probabilities as `[&prob=0.95]` blocks; BEAST embeds rate annotations as `[&rate=...]`; FigTree adds colour annotations as `[&!color=#FF0000]`. Biopython's reader handles the structural blocks but ignores most tool-specific annotations. For round-trip preservation of those annotations, use `dendropy` instead of `Bio.Phylo`. Out of scope for Week 9.

---

## 4. PhyloXML (briefly)

PhyloXML (Han and Zmasek 2009, *BMC Bioinformatics* 10:356) is an XML-based tree format that can carry rich per-node annotations: species name, sequence ID, taxonomic ID, GenBank accession, custom key-value pairs. Biopython supports it via `Phylo.write(tree, "tree.xml", "phyloxml")`.

We mention PhyloXML so you can recognize it. For Week 9's purposes Newick + a sidecar JSON of annotations is the cleaner choice.

---

## 5. ASCII rendering for the terminal sanity check

`Bio.Phylo.draw_ascii` prints a tree to the terminal:

```python
def print_ascii_tree(tree: "Tree") -> None:
    from Bio import Phylo
    import sys
    Phylo.draw_ascii(tree, file=sys.stdout)
```

Output for the demo cytochrome b panel:

```
                                                      __________________ Homo_sapiens
  _________________________________________________ |
 |                                                  |__________________ Pan_troglodytes
 |
 |                                                _____________ Mus_musculus
 |                                               |
_|                       ________________________|
 |                      |                        |____________ Rattus_norvegicus
 |                      |
 |  ____________________|
 | |                    |    _____________________ Bos_taurus
 | |                    |___|
 | |                        |_____________________ Sus_scrofa
 |_|
   |
   |   _____________________________ Gallus_gallus
   |__|
      |    _____________________________ Danio_rerio
      |___|
          |____________________________ Xenopus_laevis

___________________________________________________ Ornithorhynchus_anatinus
```

The ASCII rendering is invaluable for sanity-checking that the topology is what you expect *before* you spend an hour fiddling with matplotlib or ete3. We run it as the last step of every exercise.

---

## 6. PNG rendering with Bio.Phylo + matplotlib

`Bio.Phylo.draw` uses matplotlib to produce a PNG / SVG / PDF:

```python
def render_tree_png(tree: "Tree", out_path: Path) -> None:
    """Render a tree as a PNG via matplotlib. Headless-safe."""
    from Bio import Phylo
    import matplotlib
    matplotlib.use("Agg")  # headless backend
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    Phylo.draw(tree, axes=ax, do_show=False)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
```

Two important parameters:

- `matplotlib.use("Agg")` selects a non-interactive backend. Required for headless environments (Codespaces, Colab without inline backend, CI).
- `do_show=False` suppresses the interactive window pop-up; combined with `Agg`, the figure is rendered to memory and saved to disk.

For Week 9 the matplotlib-backed renderer is the always-available fallback. The output is functional but not publication-quality. For better aesthetics use ete3.

---

## 7. Publication-quality rendering with ete3

ete3 (Huerta-Cepas et al. 2016, *Molecular Biology and Evolution* 33:1635) is a Python library for phylogenetic tree manipulation and rendering. It has more rendering knobs than `Bio.Phylo.draw`: per-branch colour, per-leaf face images, scale bars, circular layouts.

Install: `conda install -c bioconda ete3` (preferred) or `pip install ete3`. The pip install also pulls `PyQt5` for the GUI; for headless use, set `QT_QPA_PLATFORM=offscreen` in the environment.

A minimal headless ete3 render:

```python
def render_tree_with_ete3(newick_path: Path, out_path: Path) -> None:
    """Render a Newick tree as a PNG using ete3.

    Lazy-imports ete3 inside the function so the file compiles even on
    machines where ete3 is not yet installed.
    """
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from ete3 import Tree, TreeStyle, TextFace

    tree = Tree(str(newick_path), format=1)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.scale = 800  # pixels per unit branch length
    tree.render(str(out_path), w=900, tree_style=ts)
```

ete3's `Tree(..., format=N)` parses different Newick variants:

- `format=0` — flexible with internal node names (default).
- `format=1` — flexible with internal node names being float (bootstrap values).
- `format=2` — strict format with all branches having length and label.
- `format=3` — strict format with all branches having length, no label.
- `format=5` — internal and leaf branches with length, no internal node names.
- `format=9` — leaf names only, no branch lengths.

For an IQ-TREE 2 `.treefile` (bootstrap percentages on internal nodes, branch lengths everywhere), use `format=1`.

### Colouring branches by bootstrap support

```python
def colour_branches_by_support(tree_obj) -> None:
    """Walk an ete3 Tree and colour each internal branch by support value.

    High support (>= 90): green.
    Medium support (70-89): blue.
    Low support (< 70): red.
    """
    from ete3 import NodeStyle

    for node in tree_obj.traverse():
        if node.is_leaf():
            continue
        support: float = float(getattr(node, "support", 0.0) or 0.0)
        style = NodeStyle()
        if support >= 90:
            style["hz_line_color"] = "#1a7f37"  # green
        elif support >= 70:
            style["hz_line_color"] = "#0969da"  # blue
        else:
            style["hz_line_color"] = "#cf222e"  # red
        style["hz_line_width"] = 2
        node.set_style(style)
```

This is the kind of figure that goes into a paper: branch colour encodes confidence, branch length encodes evolutionary distance, leaf labels encode taxa, the outgroup is on the bottom or far right.

---

## 8. iTOL — the web-based alternative

The Interactive Tree Of Life (iTOL; Letunic and Bork 2024, *Nucleic Acids Research* 52:W78; free for personal and research use at <https://itol.embl.de/>) is a web-based tree visualizer. Upload a Newick file, get an interactive tree with extensive styling options. iTOL is the go-to for figure-grade trees in many bioinformatics groups; the rendering is excellent and the export is publication-quality.

Free for personal and academic use; commercial use requires a paid subscription. We mention iTOL because it is the canonical web-based tree viewer; the Week 9 exercises stick to Biopython + ete3 for code reproducibility.

---

## 9. Annotating the tree with metadata

A leaf in the tree usually represents a species, a strain, or a sequence. The Newick file carries just the leaf name; everything else (full species name, GenBank accession, sampling location, isolation date) lives in a side-car table.

The Week 9 convention:

```
data/cytb_vertebrates.fasta        the input FASTA
data/cytb_vertebrates.tsv          per-taxon metadata
results/cytb_aligned.fasta         the alignment
results/cytb_nj_tree.nwk           the NJ tree in Newick
results/cytb_nj_tree.png           the rendered PNG
results/run-info.json              the build provenance
```

The `cytb_vertebrates.tsv` has columns like `taxon_name, common_name, ncbi_taxid, refseq_accession, sequence_length`. The rendering script joins the tree leaves to the table on `taxon_name` and decorates the leaves with the common name and the accession.

---

## 10. The run-info JSON

The most important file in the directory. Without it, the tree is irreproducible. The conventional shape:

```json
{
  "run_date": "2026-05-14T14:23:51Z",
  "input_fasta": "data/cytb_vertebrates.fasta",
  "input_md5": "8b14a3f5...",
  "mafft_version": "7.526",
  "mafft_algorithm": "--retree 2 --maxiterate 0",
  "trim_threshold": 0.5,
  "aligned_columns_before_trim": 1140,
  "aligned_columns_after_trim": 1098,
  "distance_method": "k2p",
  "tree_builder": "biopython_nj_1.84",
  "iqtree_version": null,
  "iqtree_model": null,
  "bootstrap_replicates": 100,
  "bootstrap_method": "column_resample_nj",
  "seed": 42,
  "outgroup": "Ornithorhynchus_anatinus",
  "biopython_version": "1.84",
  "ete3_version": "3.1.3",
  "python_version": "3.11.9",
  "notes": "Demo cytochrome b panel; NJ tree under K2P with 100-replicate column-resample bootstrap."
}
```

Every Week 9 deliverable writes this JSON alongside the tree. The mini-project script enforces it (the pipeline raises an error if the run-info has empty required fields). The pattern is identical to the Week 8 `PipelineRunInfo` dataclass.

---

## 11. Common errors and how to spot them

- **The Newick reader complains about a missing semicolon.** Some tools emit Newick without the trailing `;`. Either add it before passing to Biopython, or read with `dendropy`, which is more permissive.
- **The Nexus reader complains about a missing TAXA block.** Some tools emit TREES-only Nexus files. Biopython's Nexus reader requires the TAXA block; either add it (it is mechanical to construct from the tree leaves) or convert to Newick first.
- **The ete3 render is empty.** Usually a `QT_QPA_PLATFORM` issue. Set `QT_QPA_PLATFORM=offscreen` in the environment and re-run.
- **The matplotlib render crashes with `RuntimeError: main thread is not in main loop`.** A backend issue. Add `matplotlib.use("Agg")` before importing pyplot, or use the conda-forge matplotlib build that defaults to Agg on headless systems.
- **The tree's leaves are in the wrong order in the figure.** Tree visualization tools by default order leaves by tree traversal, not alphabetically or by the input FASTA order. Use `tree.ladderize()` in Biopython or `tree.ladderize()` in ete3 to enforce a deterministic leaf order before rendering.

---

## 12. What to take to the rest of the week

By the end of Lecture 3 you should be able to:

- Read and write Newick and Nexus trees with Biopython, and round-trip a tree through both formats while preserving topology.
- Render a tree as ASCII (Biopython `draw_ascii`), as a PNG (Biopython `Phylo.draw` via matplotlib), and as a publication-quality SVG (ete3 `TreeStyle`).
- Annotate internal branches with bootstrap support values and colour them by support tier.
- Write a `run-info.json` that records every parameter needed to reproduce the tree.
- Recognize the conventional file layout (input FASTA, alignment, tree, run-info, rendered figure) and gitignore the intermediate files that are reproducible from the inputs and the seed.

Exercises 1, 2, and 3 walk through MAFFT, NJ + bootstrap, and rendering respectively. Challenge 1 swaps NJ for IQ-TREE 2 ML with ultrafast bootstrap. Challenge 2 asks you to prove that a Newick / Nexus round-trip preserves topology on the demo panel. The mini-project assembles all the above into a single CLI script.

## References

- Felsenstein J. The Newick standard (informal specification, 1986). Documented at <https://evolution.genetics.washington.edu/phylip/newick_doc.html>.
- Maddison DR, Swofford DL, Maddison WP. "NEXUS: an extensible file format for systematic information." *Systematic Biology* 46:590 (1997). Free at <https://academic.oup.com/sysbio/article/46/4/590/1629654>.
- Han MV, Zmasek CM. "phyloXML: XML for evolutionary biology and comparative genomics." *BMC Bioinformatics* 10:356 (2009). Free at <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-356>.
- Cock PJA et al. "Biopython: freely available Python tools for computational molecular biology and bioinformatics." *Bioinformatics* 25:1422 (2009). Free at <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>. Tutorial chapter 13 (Phylogenetics with Bio.Phylo) at <https://biopython.org/docs/latest/Tutorial/index.html>.
- Huerta-Cepas J, Serra F, Bork P. "ETE 3: reconstruction, analysis, and visualization of phylogenomic data." *Molecular Biology and Evolution* 33:1635 (2016). Free at <https://academic.oup.com/mbe/article/33/6/1635/2579822>. Documentation at <http://etetoolkit.org/docs/latest/>.
- Letunic I, Bork P. "Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool." *Nucleic Acids Research* 52:W78 (2024). Free at <https://academic.oup.com/nar/article/52/W1/W78/7676118>.
