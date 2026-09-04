# Week 1 — Homework Solutions

Answers to all six Week-1 homework problems, in assignment order. Read [00-overview.md](./00-overview.md) first if you have not.

Four of these six problems are written work, and written work has no single right answer. Each of those carries a **What a grader is looking for** section describing the range that earns full marks. The model answers are one defensible path through the problem, not the path. Where a model answer names a specific paper, a specific version, or a specific number, we say where it came from so you can check it — and so you can tell when it has gone stale.

---

## Problem 1 — Read a real methods section

**Assignment:** [../homework.md#problem-1--read-a-real-methods-section](../homework.md#problem-1--read-a-real-methods-section)

### The complete answer

The paper below is open access at PLOS Genetics, published April 2025, and its Methods section is a genuinely good specimen for this exercise — not because it is flawless, but because it is *mostly* rigorous with a few real gaps you can put your finger on. That contrast is the whole assignment.

A model `notes/methods-read.md`:

```markdown
# Methods read — Liefferinckx et al. 2025

**Title.** The identification of blood-derived response eQTLs reveals complex
effects of regulatory variants on inflammatory and infectious disease risk
**Authors.** Liefferinckx C, Stern D, Perée H, Bottieau J, Mayer A, Dubussy C,
Quertinmont E, Tafciu V, Minsart C, Petrov V, Coppieters W, Karim L,
Rahmouni S, Georges M, Franchimont D
**Journal.** PLOS Genetics
**Year.** 2025 (published 10 April 2025)
**DOI.** 10.1371/journal.pgen.1011599

## Tools, with the version the authors pinned

- **PLINK v1.9** — genotype QC and data handling.
- **PLINK v2.0** — association analyses. (The authors deliberately use both;
  1.9 and 2.0 are not a simple upgrade, they have different feature sets.)
- **BCFtools v1.17** — post-imputation filtering on INFO score.
- **QTLtools v1.3.1** — the eQTL mapping itself.
- **DESeq2 v1.36.0** — expression normalisation and differential expression.

## Tools named WITHOUT a version — recorded because the absence is the finding

- **STAR** — read alignment. No version stated.
- **featureCounts** (Subread) — read quantification. No version stated.
- **UMI-tools** — UMI deduplication for the QuantSeq 3' libraries. No version.
- **R** and the **qqman** package — no R version, no package version.

## Datasets, with accession or version

- **EGA dataset `EGAD50000001320`** — the study's own RNA-seq and genotype data,
  deposited at the European Genome-phenome Archive.
- **`Homo_sapiens.GRCh38.97.gtf`** — Ensembl release 97 annotation on GRCh38,
  used as the alignment reference. Pinned exactly, which is the good habit.
- **TOPMed Imputation Reference panel**, applied via the Michigan Imputation
  Server. Panel freeze (r2 vs r3) not stated.
- **Illumina HumanOmniExpress BeadChip**, >700K SNPs — the genotyping array.
- **HapMap** — reference for population-structure checks. No release stated.
- **Reactome** — pathway enrichment. No release version stated.

## What the Methods did well, and what it left ambiguous (≈150 words)

The pinning is real where it matters most for reproducing the statistics:
PLINK 1.9 and 2.0, BCFtools v1.17, QTLtools v1.3.1 and DESeq2 v1.36.0 are all
named with versions, and the annotation is `Homo_sapiens.GRCh38.97.gtf` rather
than "GRCh38" — a specific file, not a species. The claim is stated at the right
strength throughout: these are eQTLs, associations between genotype and
transcript abundance, and the paper says "associated with" rather than
attributing a causal mechanism to any variant. That is a statistical claim about
a population, not a biological claim about a molecule, and the authors keep the
two apart.

The ambiguity sits upstream. STAR, featureCounts and UMI-tools carry no
versions, and those three determine the count matrix everything downstream is
computed from. The TOPMed panel freeze is unnamed, so the imputed genotypes are
not fully specified either. Access is honest but closed: EGAD50000001320 is a
controlled-access EGA deposit, so an independent reader can audit the method
without being able to re-run it.

## Sources

- Paper: <https://doi.org/10.1371/journal.pgen.1011599>
- EGA dataset: <https://ega-archive.org/datasets/EGAD50000001320>
```

### Why it works

**The tool list is split into two lists, and the second one is the point.** The assignment asks for "every tool the authors used, with the version they pinned." A learner who only lists the five tools that *have* versions has done half the reading. The interesting fact about this Methods section is that its version discipline is uneven: rock-solid at the statistics end, absent at the alignment end. You cannot notice that unless you write down the tools with no version next to the tools with one.

That gap is not a nitpick. `STAR` has had meaningful changes in splice handling across the 2.5 → 2.7 line, and `featureCounts` changed default behaviour for multi-overlapping reads. Both feed the count matrix that `DESeq2 v1.36.0` then analyses with perfect version fidelity. Pinning the last step of a chain and leaving the first step unpinned gives you a reproducible analysis of an unreproducible matrix.

**The claim-strength paragraph is doing the C10-voice work.** [Lecture 2 §7](../lecture-notes/02-data-ethics-and-public-data-sources.md) gives you the don't-say / say-instead tables; the assignment's acceptance criterion "distinguishes biological claim from statistical claim" is asking you to apply them to someone else's writing rather than your own. An eQTL is a clean specimen for this. "Variant *X* is an eQTL for gene *Y* in whole blood" is a statistical statement about covariation across a sample of people. It is *not* the statement "variant *X* changes the expression of gene *Y*," which is mechanistic and would need a different experiment. The paper keeps these apart; your paragraph should say that it does, and name the mechanism by which it does.

**The data-availability line is the ethics lecture landing in a real paper.** The deposit is at the European Genome-phenome Archive under a controlled-access accession. That is the correct choice for individual-level human genotype and expression data from 406 consenting donors — exactly the tier [Lecture 2 §2](../lecture-notes/02-data-ethics-and-public-data-sources.md) describes, and exactly the reason C10 works from 1000 Genomes rather than from data like this. Noting it in your read is not padding: "can I get the data?" and "did they say where the data is?" are different questions, and only the second is a Methods-quality question.

### What a grader is looking for

- **A different paper is expected and welcome.** The requirement is a recent open-access paper in one of the four named journals. If everyone answers with the same DOI, nobody read anything.
- **At least 4 tools, each with a version.** If your paper genuinely does not pin four tools, that is a finding — say so explicitly and pick a second paper rather than inventing versions. Fabricating a version number to satisfy an acceptance criterion is the single worst outcome available on this problem.
- **At least 1 dataset accession.** An accession, a release number, or a DOI. A bare URL to a lab website is not an accession.
- **The paragraph is ~150 words and does two jobs**: names something the Methods did well *and* something it left ambiguous. "It was very thorough" is not an answer; "it pins the statistical tools and not the alignment tools" is.
- **Voice.** "The authors report an association between …" not "they found the gene that causes …". The grader is reading your paragraph for exactly the phrases in the Lecture 2 §7 tables.

Acceptable range is wide on which ambiguity you pick. Missing versions, an unspecified reference build, a filtering threshold given as "standard parameters," a batch correction with no stated design matrix, and an unstated random seed are all legitimate finds. What is not acceptable is a paragraph that could have been written without opening the paper.

### Common wrong turns

**Reading the abstract and the figures instead of the Methods.** The assignment says Methods section, end to end, and it says that because the Methods is the part nobody reads and the only part that tells you whether to believe the figures. If your notes mention a result, you drifted.

**Listing tool categories instead of tools.** "They used an aligner and a differential-expression package" is a description of a genre. `STAR` and `DESeq2 v1.36.0` are tools.

**Writing "GRCh38" in the dataset list.** GRCh38 has patch releases and the annotation on top of it has its own independent release numbering. `Homo_sapiens.GRCh38.97.gtf` names one file forever; "GRCh38" names about two hundred. This is the same failure the Exercise 3 inventory is drilling — see [exercises.md](./exercises.md#exercise-3--public-data-inventory).

**Treating an unpinned tool as a fatal flaw.** It is a real gap and you should report it, but the honest framing is proportionate. This paper's statistical pipeline is more carefully pinned than most, and the correct sentence is "the version discipline is uneven," not "the paper is irreproducible."

### How to verify

There is no test runner for a reading exercise; verify structurally, then read it back out loud.

```sh
test -f notes/methods-read.md && echo present
grep -c '^- \*\*' notes/methods-read.md
```

Expected output — the file exists, and the bulleted tool and dataset entries number at least five:

```text
present
5
```

Adjust the count to your own paper; the acceptance floor is four tool entries plus one dataset entry.

Check the ~150-word paragraph is actually ~150 words. Put the paragraph in its own section and count it:

```sh
sed -n '/^## What the Methods did well/,$p' notes/methods-read.md | wc -w
```

Expected output: a number between roughly 140 and 200 (the section heading and the sources block add a few). Under 100 means you did not say enough to be wrong about anything, which is its own kind of failure.

Confirm the DOI resolves before you commit it:

```sh
curl -o /dev/null -s -w '%{http_code}\n' -L https://doi.org/10.1371/journal.pgen.1011599
```

Expected output:

```text
200
```

Then commit:

```sh
git add notes/methods-read.md
git commit -m "Read the methods section of Liefferinckx et al. 2025"
```

---

## Problem 2 — Hand-translate a sequence

**Assignment:** [../homework.md#problem-2--hand-translate-a-sequence](../homework.md#problem-2--hand-translate-a-sequence)

### The complete answer

A model `notes/hand-translation.md`. Every codon is worked by hand against the table in [Lecture 1 §4](../lecture-notes/01-the-central-dogma-in-90-minutes.md); the code cross-check is at the bottom, which is the order the acceptance criteria ask for.

````markdown
# Hand translation — first 60 nt of the SARS-CoV-2 spike CDS

Input (DNA, 5'→3', 60 nucleotides):

```
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
```

60 / 3 = 20 codons exactly, no partial codon at the end.

Transcribed to RNA (T → U), then split three at a time:

```
AUG UUU GUU UUU CUU GUU UUA UUG CCA CUA GUC UCU AGU CAG UGU GUU AAU CUU ACA ACC
```

| # | DNA codon | RNA codon | Amino acid | One-letter |
|--:|-----------|-----------|------------|:----------:|
|  1 | ATG | AUG | Methionine (also START) | M |
|  2 | TTT | UUU | Phenylalanine | F |
|  3 | GTT | GUU | Valine | V |
|  4 | TTT | UUU | Phenylalanine | F |
|  5 | CTT | CUU | Leucine | L |
|  6 | GTT | GUU | Valine | V |
|  7 | TTA | UUA | Leucine | L |
|  8 | TTG | UUG | Leucine | L |
|  9 | CCA | CCA | Proline | P |
| 10 | CTA | CUA | Leucine | L |
| 11 | GTC | GUC | Valine | V |
| 12 | TCT | UCU | Serine | S |
| 13 | AGT | AGU | Serine | S |
| 14 | CAG | CAG | Glutamine | Q |
| 15 | TGT | UGU | Cysteine | C |
| 16 | GTT | GUU | Valine | V |
| 17 | AAT | AAU | Asparagine | N |
| 18 | CTT | CUU | Leucine | L |
| 19 | ACA | ACA | Threonine | T |
| 20 | ACC | ACC | Threonine | T |

**Protein (20 residues):** `MFVFLVLLPLVSSQCVNLTT`

No stop codon appears — expected, because 60 nt is the start of a 3,822 nt
coding sequence, not a complete one. The stop is 1,254 codons further along.

## Code cross-check

```
$ python challenges/challenge-01-reverse-complement.py \
    ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
input          ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
length         60 nt
GC content     33.33%
revcomp        GGTTGTAAGATTAACACACTGACTAGAGACTAGTGGCAATAAAACAAGAAAAACAAACAT
protein        MFVFLVLLPLVSSQCVNLTT
```

Hand result and code result agree.
````

### Why it works

**The frame is fixed by the first base, and nothing else.** Translation is not a property of the string; it is a property of the string *plus an offset*. This sequence is a coding sequence starting at its own `ATG`, so frame 1 is correct and the twenty codons fall out cleanly. Slide the window by one base and you get a completely different, completely meaningless peptide — `CFFSCFIATSLLVSVLIL…` — with no warning from any tool. Every "my protein looks like garbage" bug in the next eleven weeks is a frame bug until proven otherwise.

**Transcribing before you look up is a discipline, not a requirement.** The codon table in [Lecture 1 §4](../lecture-notes/01-the-central-dogma-in-90-minutes.md) is printed in RNA — `UUU Phe`, not `TTT Phe`. You can absolutely read `TTT` off a DNA sequence and mentally map it to Phe, and after a hundred codons you will do it without thinking. On your first twenty, write out the U's. The table and your work then sit in the same alphabet, and a mismatch is visible instead of being a mental slip you cannot audit.

**Look at codons 2, 4 and 6 — and 5, 7, 8, 10 and 18.** `UUU` appears twice for phenylalanine. `GUU` appears three times for valine. Leucine arrives four separate times from three different codons (`CUU`, `UUA`, `UUG`, `CUA`). That is the redundancy of the genetic code showing up in twenty codons of real sequence, and it is why [Lecture 1 §4](../lecture-notes/01-the-central-dogma-in-90-minutes.md) asks you to be able to explain why a third-position change is often silent. Change codon 5 from `CUU` to `CUC`, `CUA` or `CUG` and the protein is unchanged. Change its *first* base to `GUU` and leucine becomes valine.

**The absent stop codon is a positive result, not a gap.** Twenty codons in and no `UAA`/`UAG`/`UGA` is exactly what a fragment from the 5' end of a long CDS should look like. If a stop *had* appeared at, say, codon 8, that would tell you either the frame was wrong or the sequence was not what you thought it was. Reading the absence of something as evidence is a habit worth building now.

**The signal peptide, if you want the biology.** `MFVFLVLLPLVSS` is the N-terminal signal peptide of the SARS-CoV-2 spike protein — a hydrophobic run (F, V, F, L, V, L, L, P, L, V) that targets the nascent chain to the secretory pathway. You are not asked to know that. It is here so you can see that the twenty letters you just derived by hand are the beginning of a real, named, much-studied molecule rather than an exercise string.

### Common wrong turns

**Splitting the codons with your eyes on a wrapped line.** Sixty characters is exactly the FASTA wrap width, so if you copy this sequence out of a FASTA file it is one line — but if you copy it out of a PDF or a slide it may not be. A single dropped or duplicated character shifts every codon after it and the tail of your protein turns to noise while the head stays correct. That asymmetry is the diagnostic: **a hand translation that is right at the start and wrong at the end is a frame shift, not twenty independent mistakes.** Count your codons before you look any of them up. Twenty groups of three, 60 characters.

**Flipping the U and the T in the wrong direction.** DNA `T` becomes RNA `U`. There is no `U` in DNA and no `T` in RNA. If you write `AUG` in the DNA column or `ATG` in the RNA column you will still get the right amino acid — but the error compounds when you build the `CODON_TABLE` dict for Challenge 1 and half your keys are in the wrong alphabet. The characteristic failure there is immediate and total:

```text
ValueError: unknown codon 'ATG' at position 0
```

**Reading the codon table with the third base as the row.** The table in the lecture is indexed first base = row block, second base = column, third base = row *within* the block. `UUA` is Leu; `AUU` is Ile. Reversing the axes turns one into the other and there is nothing in your output to tell you. If two of your twenty codons disagree with the code and eighteen agree, you almost certainly read the table transposed for those two.

**Trusting the code over your hand and stopping there.** The acceptance criteria ask for both because the *disagreement* is where the value is. If your hand answer and `translate()` differ, do not silently adopt the code's answer — [Homework Problem 2's own hint](../homework.md#problem-2--hand-translate-a-sequence) names the three candidates (codon-boundary off-by-one, flipped U/T, `CODON_TABLE` typo) and finding which one it was takes about ninety seconds. A `CODON_TABLE` typo will pass every self-test in Challenge 1 and quietly corrupt one amino acid forever.

### How to verify

The expected protein string is given to you in the assignment, so this one is checkable exactly.

```sh
python challenges/challenge-01-reverse-complement.py ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
```

Expected output, exactly:

```text
input          ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
length         60 nt
GC content     33.33%
revcomp        GGTTGTAAGATTAACACACTGACTAGAGACTAGTGGCAATAAAACAAGAAAAACAAACAT
protein        MFVFLVLLPLVSSQCVNLTT
```

Check your written table independently of the code — this pulls the one-letter column out of your Markdown table and joins it:

```sh
grep -c '^| *[0-9]' notes/hand-translation.md
```

Expected output — twenty codon rows:

```text
20
```

And confirm the protein line in your notes is the string the assignment specifies:

```sh
grep -o 'MFVFLVLLPLVSSQCVNLTT' notes/hand-translation.md | head -1
```

Expected output:

```text
MFVFLVLLPLVSSQCVNLTT
```

Then commit:

```sh
git add notes/hand-translation.md
git commit -m "Hand translation of the first 60 nt of the spike CDS"
```

---

## Problem 3 — Re-identification thought experiment

**Assignment:** [../homework.md#problem-3--re-identification-thought-experiment](../homework.md#problem-3--re-identification-thought-experiment)

### The complete answer

A model `notes/re-id-essay.md`. Four numbered paragraphs, 463 words, one in-body citation of Gymrek et al. 2013.

```markdown
# Re-identification: reading Gymrek et al. 2013

**1. What data went in.** Gymrek et al. (2013) began from material that was
already public and already considered de-identified: whole-genome sequence from
male participants in public sequencing projects. From that sequence they
profiled short tandem repeats on the Y chromosome — Y-STRs — which reduce to a
small vector of repeat counts at a handful of loci. That is a few dozen numbers,
not a genome, and it is the part of the input that carries surname information.
Alongside it they used two ordinary pieces of metadata distributed with the
samples: the donor's approximate age and their state of residence. Neither is a
protected identifier under the HIPAA Privacy Rule, which is exactly why both
were sitting in the open tier next to the sequence.

**2. What they cross-referenced.** Two recreational genetic-genealogy databases:
Ysearch and the Sorenson Molecular Genealogy Foundation collection. Hobbyist
genealogists had voluntarily deposited their own Y-STR haplotypes next to their
surnames, so querying a target haplotype returns ranked candidate surnames. A
surname is not an identity, so the authors then triangulated: surname plus age
plus state, run through free public-record search resources, narrows a US male
to a short candidate list. Every resource in the chain was free and publicly
accessible. Reviews of the genome-privacy literature report that the surname
step succeeded for roughly 12% of US males analysed, and that around fifty
individuals from the CEPH Utah collection were traced this way. The reported
institutional response was to remove the age metadata for those samples from
open access.

**3. Which sex was vulnerable, and why.** Males. The Y chromosome passes from
father to son largely intact, and in many Western naming conventions the surname
travels the same path. That shared inheritance is the entire basis of the
correlation the attack exploits — it is a fact about naming customs at least as
much as about biology. Females carry no Y chromosome, so this specific route
does not reach them directly. It is worth stating that carefully rather than
comfortably: a woman's father, brothers and sons carry haplotypes that are
informative about her family, so "not vulnerable to this attack" is not the same
claim as "not exposed."

**4. What a 2026 consent form should say.** Three things, plainly. First, that
de-identification is a procedure and not a guarantee — the sequence is itself an
identifier, and published methods recover identity by combining it with public
records the study does not control. Second, that the risk is not fixed: genealogy
databases have grown steadily since 2013, so the honest statement to a donor is
about a trend under uncertainty, not a probability. Third, that donation is
effectively irrevocable, because a file can be withdrawn from a portal but
copies already downloaded cannot be recalled. I would also say what the data does
*not* do: it is associated with health outcomes, it does not determine them.
```

### Why it works

**Paragraph 1 separates the two inputs, and the separation is the insight.** It is tempting to write "they used genomes." They did not, in the sense that matters. They used a Y-STR haplotype — a low-dimensional summary — plus *demographic metadata that was not considered identifying*. The attack is a join across two datasets, and the second dataset was released precisely because nobody thought age and state were sensitive. [Lecture 2 §1.4](../lecture-notes/02-data-ethics-and-public-data-sources.md) makes the general claim ("assume re-identification is possible"); the specific mechanism is that de-identification decisions are made per-field and attacks operate across fields.

**Paragraph 2 names the databases, which is the acceptance criterion doing real work.** Ysearch and SMGF were built by hobbyists for hobbyists. Nobody who uploaded a haplotype to trace their family tree was consenting to underwrite a re-identification method — and, exactly as in the Golden State Killer case in [Lecture 2 §1.2](../lecture-notes/02-data-ethics-and-public-data-sources.md), the people affected were not the people who uploaded. That is the recurring structure of genetic privacy failures: **the consenting party and the exposed party are different people.** An essay that gets this shape gets the lecture.

**Paragraph 3 answers the question and then refuses the easy stopping point.** "Males, because Y chromosomes and surnames are both patrilineal" is the correct answer and it is one sentence. The paragraph earns its length by being precise about what the finding does *not* say. Writing "women are safe" would be a determinism-adjacent overclaim of exactly the kind [Lecture 2 §7](../lecture-notes/02-data-ethics-and-public-data-sources.md) is training you out of — it converts "this attack has no direct route" into "no exposure," which is a stronger claim than the evidence supports.

**Paragraph 4 is the one a grader actually reads for voice.** It commits to three specific disclosures rather than gesturing at "informed consent." Note the hedging where hedging is honest — "a trend under uncertainty, not a probability" — and the flat statement where flatness is honest — "copies already downloaded cannot be recalled." That is the register of [Lecture 2 §6](../lecture-notes/02-data-ethics-and-public-data-sources.md): say the uncertain thing with a hedge and the certain thing without one. Reversing that pairing is what makes consent forms either alarmist or misleading.

The final sentence — "it is associated with health outcomes, it does not determine them" — is doing double duty. It is a real thing a consent form should say, and it is the essay demonstrating in its own last line that it can hold the associated/causal distinction.

### What a grader is looking for

- **400–500 words, four numbered paragraphs, one per question.** The model is 463. Count before you commit; this is the criterion most often failed by accident.
- **Gymrek et al. 2013 cited by author and year in the body**, not only in a footer. "Gymrek et al. (2013) began from…" satisfies it; a bibliography alone does not.
- **Question 1 names both input types.** An essay that only says "genomes" has missed that the attack is a join. Y-STRs *and* the demographic metadata.
- **Question 2 names the genealogy databases** and, ideally, notes that the surname alone was not sufficient — the triangulation step is where identity actually falls out.
- **Question 3 says males, and gives the patrilineal-inheritance reason.** Full marks need the reason. Extra credit, informally, for noticing that the reason is partly cultural rather than purely biological.
- **Question 4 is graded on specificity, not on agreeing with us.** Any three concrete disclosures work. "We would tell them about the risks" is not a disclosure. A donor-facing sentence you would actually be willing to sign your name under is.
- **Voice.** No determinism. Uncertainty cited honestly. If you quote a number — the 12%, the ~50 individuals — attribute it, because those figures come from the paper and from subsequent genome-privacy reviews rather than from thin air.

Acceptable variation on Question 4 is broad. Reasonable alternatives to our three: a plain-language statement that relatives are affected; an explicit description of the access tier the data will sit in and who can apply; a commitment to notify donors if the access tier changes; a statement about secondary use by law enforcement. Any of those, stated concretely, earns the marks.

### Common wrong turns

**Summarising the paper instead of answering the four questions.** The assignment is four numbered questions and the acceptance criteria say four numbered paragraphs. A flowing essay that covers the same ground scores worse, because the structure *is* part of what is being assessed — this is the shape a real ethics-review response takes.

**Writing "they hacked the database."** Nothing was breached. Every resource used was public and free, and that is precisely why the result mattered. An essay that frames it as an intrusion has inverted the lesson: the failure mode here is *legitimate* data used in an unanticipated combination.

**Concluding that public data is therefore unsafe to use.** The opposite conclusion is the one the field actually drew. The response to Gymrek et al. was better consent language and better-designed access tiers, not the withdrawal of public data. 1000 Genomes is still open, still the right dataset for C10, and its donors consented with this class of risk in view. If your essay ends with "so we should not use public genomic data," re-read [Lecture 2 §2](../lecture-notes/02-data-ethics-and-public-data-sources.md).

**Determinism creep in the consent-form paragraph.** "We will tell donors their genome reveals their disease risk" overclaims twice — "reveals" and the bare "risk". A defensible version names the strength of the claim: "may indicate an elevated statistical risk for some conditions, with substantial uncertainty for any individual."

**Padding to hit 400 words.** The word floor exists to stop three-sentence answers, not to reward restatement. If you are at 340 words and stuck, the missing content is almost always in Question 4 — it is the only one of the four that asks you to design rather than to report.

### How to verify

Word count first, because it is the criterion most often missed:

```sh
wc -w notes/re-id-essay.md
```

Expected output: a number between 400 and 500 (the model body is 463; the `#` title line adds a few, so aim for the middle of the range rather than the edge).

Then the structural checks:

```sh
grep -c '^\*\*[1-4]\.' notes/re-id-essay.md
grep -c 'Gymrek' notes/re-id-essay.md
```

Expected output — four numbered paragraphs, and at least one in-body citation:

```text
4
1
```

Then the voice check. This greps for the determinism vocabulary from the [Lecture 2 §7](../lecture-notes/02-data-ethics-and-public-data-sources.md) don't-say table:

```sh
grep -in 'causes\|determines\|gene for\|anonymized\|made him' notes/re-id-essay.md
```

Expected output: nothing, unless the hit is inside a sentence where you are explicitly quoting the bad phrasing in order to criticise it.

Then commit:

```sh
git add notes/re-id-essay.md
git commit -m "Re-identification essay on Gymrek et al. 2013"
```

---

## Problem 4 — Build a tiny FASTA validator

**Assignment:** [../homework.md#problem-4--build-a-tiny-fasta-validator](../homework.md#problem-4--build-a-tiny-fasta-validator)

### The complete answer

`homework/p4_fasta_validate.py`. Stdlib only, no third-party imports.

```python
#!/usr/bin/env python3
"""Homework Problem 4 — a tiny FASTA validator. Stdlib only.

Usage:
    python homework/p4_fasta_validate.py path/to/file.fasta

Prints exactly one summary line, then one WARNING line per problem found.
Exit status 0 if there were no warnings, 1 if there was at least one.
"""

from __future__ import annotations

import sys
from typing import Iterator, TextIO

ALLOWED = frozenset("ACGTN")


def parse_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) for every record. Same parser as Exercise 2."""
    header: str | None = None
    chunks: list[str] = []

    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("sequence data before the first '>' header line")
            chunks.append(line)

    if header is not None:
        yield header, "".join(chunks)


def validate(path: str) -> tuple[str, list[str]]:
    """Return (summary_line, warnings) for one FASTA file."""
    warnings: list[str] = []
    seen: set[str] = set()

    records = 0
    total_bp = 0
    gc = 0
    acgt = 0

    with open(path, "r", encoding="utf-8") as handle:
        for index, (header, seq) in enumerate(parse_fasta(handle), start=1):
            records += 1
            total_bp += len(seq)

            upper = seq.upper()
            for base in upper:
                if base in "ACGT":
                    acgt += 1
                    if base in "GC":
                        gc += 1

            if not seq:
                warnings.append(
                    f"WARNING: record {index} {header!r} has an empty sequence"
                )

            bad = sorted(set(upper) - ALLOWED)
            if bad:
                shown = ", ".join(repr(character) for character in bad)
                warnings.append(
                    f"WARNING: record {index} {header!r} has "
                    f"{len(bad)} character(s) outside {{A,C,G,T,N}}: {shown}"
                )

            if header in seen:
                warnings.append(
                    f"WARNING: record {index} repeats header {header!r}"
                )
            seen.add(header)

    mean = total_bp / records if records else 0.0
    gc_pct = 100.0 * gc / acgt if acgt else 0.0
    summary = (
        f"{records} records, {total_bp} bp total, "
        f"mean length {mean:.1f} bp, GC {gc_pct:.2f}%"
    )
    return summary, warnings


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {argv[0]} path/to/file.fasta", file=sys.stderr)
        return 2

    summary, warnings = validate(argv[1])
    print(summary)
    for warning in warnings:
        print(warning)
    return 1 if warnings else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
```

The deliberately broken test file — save as `data/broken.fasta`. It contains exactly one instance of each of the three specified problems:

```
>seq1 first test record
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC
>empty record
>seq3 third test record
GGGGCCCCXXGGGGCCCC
>seq1 first test record
ACGTACGT
```

Record 2 has no sequence. Record 3 contains `X`. Record 4 repeats record 1's header verbatim.

### Why it works

**One pass, two outputs, printing deferred to the end.** The summary needs totals that are only known after the last record, and the acceptance criteria say the summary comes first. So `validate()` accumulates warnings in a list instead of printing them where they are found, and `main()` prints the summary and then the list. This is worth doing on purpose rather than reaching for a two-pass version: the file is opened once, and on a real multi-gigabyte reference FASTA the second pass is the difference between a fast check and a slow one.

**The exit code is computed from the warning list, not tracked with a flag.** `return 1 if warnings else 0`. A separate `had_error = True` boolean set in four places is the version that eventually grows a fifth check where somebody forgets to set it, and then the script exits 0 on a file it just warned about. Deriving the status from the data means the two can never disagree.

**About "one warning per problem."** The assignment says the broken-file run should print "the summary plus one warning per problem." That sentence has two readings — one line per *kind* of problem, or one line per *instance*. This implementation emits one line per instance, and the test file is built with exactly one instance of each kind so both readings produce the same three lines. Per-instance is the more useful behaviour for a real file: a FASTA with four empty records should tell you which four, not that "some records are empty." If you chose per-kind, say so in a comment and make your test file prove it; either is defensible, an undocumented choice is not.

**`set(upper) - ALLOWED` is the whole invalid-character check.** Building a set of the distinct characters in a record and subtracting the allowed alphabet gives you the offending characters directly, without a per-character branch and without a regex. It is *O(n)* in the sequence length with a tiny constant, it reports every distinct bad character rather than only the first, and `sorted(...)` makes the output deterministic — which matters, because Python's set iteration order for strings varies between runs under hash randomisation and you do not want a test that passes four times out of five.

**Case is handled once, at the top, and only for checking.** `upper = seq.upper()` feeds both the GC count and the character check, so a soft-masked reference region (lowercase `acgt`, the repeat convention from [Exercise 1's glossary](./exercises.md#exercise-1--glossary-in-your-own-words)) is neither miscounted as bad nor excluded from the GC denominator. `len(seq)` uses the original, because case does not change length.

**The two statistics are computed over different denominators, on purpose.** `total_bp` counts every character in every record, including `N` and including invalid ones — it is a length, and lengths do not get to be opinionated. `gc_pct` divides only by the `A`/`C`/`G`/`T` count, matching the definition `GC_content` uses in [Challenge 1](./challenges.md#challenge-1--reverse-complement-gc-and-translate). Mixing those denominators would make GC% drift downward on any file with masked regions. The rule from the challenge applies here too: pick a definition, write it down, and be consistent with your own other code.

**Why this file re-declares `parse_fasta` instead of importing it.** The assignment permits either stdlib-only or importing from `exercise-02-fasta-by-hand.py`. The catch is that the filename contains hyphens, so `import exercise-02-fasta-by-hand` is a syntax error — hyphens are the subtraction operator. Importing it needs `importlib`:

```python
import importlib.util
spec = importlib.util.spec_from_file_location(
    "ex2", "exercises/exercise-02-fasta-by-hand.py")
ex2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex2)
parse_fasta = ex2.parse_fasta
```

That works and it is worth knowing. We chose the self-contained version so the script runs from any working directory with no path assumptions, which is what you want from a validator you will point at arbitrary files. If you do import, note that the loader path above is relative to your *current directory*, not to the script — a foot-gun the standalone version simply does not have.

### Common wrong turns

**Printing warnings as you find them.** You get the warnings before the summary, which fails the "prints the summary plus one warning per problem" ordering, and there is no clean fix without buffering — which is what the list already is.

**Using `sys.exit(1)` inside the record loop.** It exits on the *first* problem, so a file with all three problems reports one. The acceptance criteria want one warning per problem; you must finish the scan.

**Forgetting that `exit code` and `printed output` are separate channels.** A common near-miss is a script that prints everything correctly and always exits 0, because `main()` returns `None`. Check it explicitly — `echo $?` after every run. `None` is falsy but `raise SystemExit(None)` exits 0, so this failure is completely invisible in the output.

**Checking `if seq == ""` after `.strip()` has already run on each line.** That part is fine. The subtle version is checking `if not header` for the duplicate test — a record can legitimately have an empty header (`>` alone), and `not header` is `True` for it. Use `header in seen`, which treats `""` as a perfectly good key.

**Counting GC with `seq.count("G") + seq.count("C")`.** `str.count` is case-sensitive, so lowercase soft-masked regions contribute zero and your GC% comes out implausibly low on exactly the repeat-rich sequence where you would expect it high. Same trap as in [Challenge 1](./challenges.md#challenge-1--reverse-complement-gc-and-translate).

**Treating `N` as an invalid character.** `N` is in `ALLOWED` because the assignment says the allowed set is `{A, C, G, T, N}`, and because real reference FASTA is full of it. A validator that warns on every `N` in GRCh38 emits tens of millions of warnings and is useless.

**Dividing by zero on an empty file.** `records` is 0, so `total_bp / records` raises `ZeroDivisionError`. Both guards are in the code above. Note the deliberate scope decision: an empty file is not one of the three specified problems, so it prints `0 records, 0 bp total, mean length 0.0 bp, GC 0.00%` and exits 0. If you want it to warn, that is a defensible extension — document it, because it changes the exit code contract.

### How to verify

**On a valid file.** Use the `data/good.fasta` from [Exercise 2](./exercises.md#exercise-2--fasta-by-hand):

```sh
python homework/p4_fasta_validate.py data/good.fasta
echo $?
```

Expected output, exactly — one summary line, zero warnings, exit 0:

```text
3 records, 96 bp total, mean length 32.0 bp, GC 47.92%
0
```

The arithmetic is checkable by hand: 60 + 20 + 16 = 96 bp over 3 records is a mean of 32.0. GC is 20 (in `seq1`) + 10 (in `seq2`) + 16 (in `seq3`) = 46 over 96 unambiguous bases = 47.9166…%, which rounds to 47.92.

**On the broken file.**

```sh
python homework/p4_fasta_validate.py data/broken.fasta
echo $?
```

Expected output, exactly — summary plus three warnings, exit 1:

```text
4 records, 86 bp total, mean length 21.5 bp, GC 47.62%
WARNING: record 2 'empty record' has an empty sequence
WARNING: record 3 'seq3 third test record' has 1 character(s) outside {A,C,G,T,N}: 'X'
WARNING: record 4 repeats header 'seq1 first test record'
0
```

Note the trailing `0` there is **wrong** if you see it — it should be `1`. This is the exit-code trap called out above, and the reason the verification prints `$?` at all. A correct run ends:

```text
WARNING: record 4 repeats header 'seq1 first test record'
1
```

**On a bad invocation:**

```sh
python homework/p4_fasta_validate.py
echo $?
```

Expected output — usage on stderr, exit 2 (distinct from both 0 and 1, so a calling script can tell "you used it wrong" from "the file is bad"):

```text
usage: homework/p4_fasta_validate.py path/to/file.fasta
2
```

Then commit:

```sh
git add homework/p4_fasta_validate.py data/broken.fasta
git commit -m "Homework 4: tiny FASTA validator with summary and warnings"
```

---

## Problem 5 — Refine your glossary

**Assignment:** [../homework.md#problem-5--refine-your-glossary](../homework.md#problem-5--refine-your-glossary)

### The complete answer

This problem has three deliverables — the read-aloud edit, five second "easy to confuse with" notes, and entry 21 — and one meta-requirement that is easy to miss: **the diff must show an edit, not a rewrite.** Work in the same file, commit on top of the v1 commit.

**Entry 21, the model answer:**

```markdown
### 21. The vocabulary problem

The failure mode where a biologist and an engineer use the same English word for
two different things and both leave the conversation confident they agreed. It is
not a beginner's confusion that wears off — it persists for years in mixed teams,
because both parties are using their word correctly within their own field.

*Example:* "coverage." I say a module has 80% coverage and mean the fraction of
lines a test suite executes. A biologist says a sample has 30× coverage and means
the average number of reads stacked over each base. Same word, unrelated units,
and neither of us is wrong.
*Easy to confuse with:* **jargon**. Jargon is a word you have not learned yet, and
the fix is to look it up. The vocabulary problem is a word you *have* learned,
used against a different definition, and the fix is to ask "do you mean that in
the biology sense or the software sense?" — out loud, in the meeting, every time
you are not certain.
*Easy to confuse with:* **a naming disagreement**. A naming disagreement is
resolvable by picking one name. The vocabulary problem is not, because both
communities are much larger than your team and neither is going to renumber its
literature for you. You manage it; you do not fix it.
```

**Five second "easy to confuse with" notes**, added to entries that already have one. These are chosen the way the assignment asks — because they are mistakes people actually make, not because they are the tidiest pairings:

```markdown
### 2. Genome
*Easy to confuse with:* **exome** — the protein-coding subset only, ~1–2% of the
human genome.
*Easy to confuse with:* **karyotype**. A karyotype is the chromosome *count and
gross structure* of a cell — 46,XX — with no sequence in it at all. "We looked at
the genome" and "we looked at the karyotype" are different experiments run on
different instruments, and the second one cannot see a SNP.

### 4. Gene
*Easy to confuse with:* **ORF (open reading frame)** — a stretch between a start
codon and an in-frame stop.
*Easy to confuse with:* **locus**. A locus is any named position or region on a
chromosome; a gene is one kind of locus. GWAS hits are reported per *locus*
precisely because the associated region often contains several genes and the
causal one is unknown. Writing "the gene associated with X" when the paper said
"the locus associated with X" silently invents a result.

### 12. Variant (genomic sense)
*Easy to confuse with:* **mutation** — which carries a pathogenicity connotation
the word "variant" deliberately drops.
*Easy to confuse with:* **allele**. A variant is the difference from the
reference; an allele is one of the specific sequences present at that position.
A biallelic site has a REF allele and an ALT allele and *one* variant record.
Counting alleles when you meant variants doubles your numbers.

### 14. FASTA
*Easy to confuse with:* **FASTQ** — similar name, four-line records, quality
scores.
*Easy to confuse with:* **the FASTA *program***. FASTA was a sequence-search tool
by Pearson and Lipman before it was a file format, and it still exists. "We ran
FASTA" and "we wrote a FASTA" are different sentences, and search results for
the format are polluted by the tool.

### 20. Coverage
*Easy to confuse with:* **breadth of coverage** — what fraction of the target was
covered at all, as opposed to how deeply.
*Easy to confuse with:* **read count**. Coverage is normalised by length; read
count is not. Two samples with identical read counts have different coverage if
their target sizes differ, which is why exome and genome runs are never compared
on read count.
```

The read-aloud pass, in practice, is where the word budget gets recovered. Example — entry 5 from the draft in [exercises.md](./exercises.md#exercise-1--glossary-in-your-own-words), before and after:

```diff
 ### 5. Transcript

-One specific RNA molecule produced from a gene. A gene is a region; a
-transcript is a product. Alternative splicing means one gene routinely produces
-several different transcripts (isoforms), which is why transcript-level and
-gene-level counts differ in RNA-seq.
+One specific RNA molecule produced from a gene: a gene is a region, a transcript
+is a product. Alternative splicing means one gene routinely produces several
+isoforms, which is why transcript-level and gene-level counts differ in RNA-seq.
```

Twelve words shorter, same content, and it now reads the way you would actually say it.

Commit message, as the assignment specifies:

```sh
git commit -am "Glossary v2: tightened wording, added vocabulary problem entry"
```

### Why it works

**Entry 21 is defined by its failure mode, not by its topic.** "The vocabulary problem is when words mean different things in different fields" is true, useless, and would be the entry most people write. The version above is a definition you could *act on*: it names who is affected (mixed teams), why it persists (both parties are locally correct), and what you do about it (ask, out loud, every time). [Lecture 1 §7](../lecture-notes/01-the-central-dogma-in-90-minutes.md) makes the same move — its table is a list of concrete word-pairs, not an abstract principle.

**The two "easy to confuse with" notes on entry 21 are load-bearing.** Distinguishing the vocabulary problem from *jargon* is the distinction that changes your behaviour: jargon has a lookup fix, this does not. Distinguishing it from *a naming disagreement* is the one that stops you from trying to solve it. Both are conclusions you can only reach by holding the concept against a near neighbour — which is the entire pedagogical argument for the "easy to confuse with" field, applied recursively to the entry about vocabulary itself.

**The five second-notes are chosen for the error, not the symmetry.** Look at what each one prevents. Gene/locus stops you misreporting a GWAS result. Variant/allele stops you double-counting. FASTA-format/FASTA-program stops you from an hour of bad search results. Genome/karyotype stops you from claiming an instrument saw something it cannot see. Coverage/read-count stops an invalid cross-experiment comparison. Every one is a mistake with a downstream consequence you can name — which is the standard the assignment sets when it says the second note should come from something "someone in the C10 community Slack actually got it wrong" about.

**The edit-don't-rewrite rule has a real reason behind it.** [The assignment's hint](../homework.md#problem-5--refine-your-glossary) puts it plainly — "the act of editing is the learning." Mechanically: if you delete v1 and write v2 from scratch, you write it from the same memory that produced v1, so you reproduce the same gaps. Editing forces you to *read* your own sentence, which is a different cognitive act from producing one. The diff is the evidence that you did the reading, which is why it is an acceptance criterion.

**Polishing means getting shorter.** This is counter-intuitive and it matters for the mini-project, where the glossary has a hard length target of roughly 800–1,200 words for 21 entries. The v1 draft in [exercises.md](./exercises.md#exercise-1--glossary-in-your-own-words) is well over that with only 20 entries. So v2 adds an entry and five notes *and still has to come in shorter* — which is only possible by cutting. Read-aloud is the tool: the sentences that are hard to say aloud are almost always the ones carrying the padding.

### What a grader is looking for

- **Entry 21 exists and is your own definition**, not a paraphrase of the lecture's table. The lecture gives examples; the entry should give a definition plus a behaviour.
- **Five entries carry a second "easy to confuse with."** Not five *new* entries — five *existing* entries with a second note appended.
- **The diff is visible in history.** `git log` should show a v1 commit and a v2 commit against the same path. A single commit adding the whole file fails this criterion even if the content is perfect.
- **v2 is measurably different from v1 across many entries**, not only the five you were told to touch. The read-aloud pass should leave fingerprints everywhere.
- **The commit message names what changed.** `Glossary v2: tightened wording, added vocabulary problem entry` — the assignment gives you this one verbatim, so use it or something equally specific. `update` is not a commit message.

Acceptable variation: which five entries get a second note is entirely yours, and your picks should be *your* confusions rather than ours. If you genuinely have never mixed up variant and allele, do not add that note — add the one for the pair you actually got wrong this week. The grader is checking that the notes are specific and consequential, not that they match this list.

### Common wrong turns

**"I re-read it and it was fine."** Then you did not read it aloud. Reading aloud is not a metaphor here; it is a mechanical trick that surfaces clause-stacking and hedging that the eye skips. Twenty entries takes about six minutes.

**Deleting the file and starting over.** Fails the diff criterion, and — more importantly — throws away the comparison that was the entire point. If you have already done it, recover: `git show HEAD~1:notes/week-01-glossary.md` prints v1, and you can diff by hand and describe the changes in your commit message.

**Adding a second "easy to confuse with" that is a restatement of the first.** "Easy to confuse with FASTQ. Also easy to confuse with the FASTQ format." That is one note written twice. The second note has to name a *different* neighbour with a *different* failure mode.

**Making entry 21 about biology.** The vocabulary problem is not a biology concept — it is a communication failure mode that happens to be acute in this field. An entry that defines it as "when biology terms are hard" has missed [Lecture 1 §7](../lecture-notes/01-the-central-dogma-in-90-minutes.md) entirely; the table there is symmetric, and the engineer is wrong exactly as often as the biologist.

**Letting v2 grow.** Adding an entry and five notes while lengthening every definition produces a two-page glossary that fails the mini-project's length criterion two days later. Add content, cut prose, come out shorter.

### How to verify

Structure first — 21 entries, and at least 26 confusion notes (21 entries × 1, plus 5 second notes):

```sh
grep -c '^### ' notes/week-01-glossary.md
grep -c '^\*Easy to confuse with:\*' notes/week-01-glossary.md
```

Expected output:

```text
21
26
```

Check entry 21 is present by name:

```sh
grep -n '^### 21\.' notes/week-01-glossary.md
```

Expected output (the line number will differ):

```text
243:### 21. The vocabulary problem
```

Then prove the file was *edited*, not replaced. This is the criterion with teeth:

```sh
git log --oneline -- notes/week-01-glossary.md
```

Expected output — at least two commits touching the same path:

```text
a1b2c3d Glossary v2: tightened wording, added vocabulary problem entry
e4f5a6b Week 1 glossary, first pass
```

And look at the shape of the change:

```sh
git diff HEAD~1 --stat -- notes/week-01-glossary.md
```

Expected output — additions *and* deletions. If the deletions column is 0, you appended without editing:

```text
 notes/week-01-glossary.md | 61 ++++++++++++++++++++++++++++++++-------------
 1 file changed, 44 insertions(+), 17 deletions(-)
```

Finally, the length check that the mini-project will enforce in two days:

```sh
wc -w notes/week-01-glossary.md
```

Expected output: a number in the 800–1,200 range. Over 1,400 and the polish pass is not finished.

---

## Problem 6 — Mini reflection essay

**Assignment:** [../homework.md#problem-6--mini-reflection-essay](../homework.md#problem-6--mini-reflection-essay)

### The complete answer

A model `notes/week-01-reflection.md`. Four numbered paragraphs, 348 words. This is a *model* in the narrowest sense — it demonstrates the shape and the specificity a good reflection has. The content should be yours and will not resemble it.

```markdown
# Week 1 reflection

**1. Easiest and hardest.** The central dogma was easiest, and I think that is
because it is a pipeline. DNA in, mRNA out, protein out — that is a shape I
already had. Vocabulary was hardest, and not for the reason I expected. I could
recite the definitions after one read. What I could not do was keep gene and
transcript apart under any pressure at all; the first time I tried to explain
RNA-seq quantification to myself I used them interchangeably in one sentence and
did not notice until I re-read it. Ethics sat in the middle: the rules are short,
the reasoning behind Rule 2 took a second pass.

**2. What I had wrong.** Two things. I believed a gene was, by definition, a
thing that codes for a protein — so I had no mental slot for the rRNA genes or
for lncRNAs, which turn out to be a large fraction of what is annotated. And I
believed the reference genome was somebody's actual genome. Learning that GRCh38
is a curated mosaic, and that individuals differ from it at millions of
positions as a matter of course, reframed what "variant" means. A variant is a
statement about a comparison. I had been reading it as a statement about a
person.

**3. The dataset I want to work with.** GISAID, for Week 9. The phylogenetics
module is the first place in the course where the output is a claim about
history — which lineage descends from which — rather than a summary statistic,
and I want to see how much confidence a bootstrap value actually buys. The
mandatory acknowledgement table is also the first time I will have to do
attribution properly, at scale, rather than dropping in a citation.

**4. What I want next.** How annotation gets made. Every week from here consumes
a GTF as if it were ground truth, and I now know it is a curated artefact with a
release number, produced by a pipeline plus human curators. I would like to know
what that pipeline is, how disagreements between RefSeq and Ensembl arise, and
how anybody decides a transcript is real.
```

### Why it works

**Every paragraph names a specific thing.** Not "vocabulary was hard" but "I used gene and transcript interchangeably in one sentence and did not notice until I re-read it." Not "I want to use a dataset" but GISAID, in Week 9, for a stated reason. A reflection made of specifics is re-readable in Week 12; a reflection made of generalities is not, and re-reading it later is the only purpose it has.

**Paragraph 2 answers the actual question, which is about being wrong.** "Did anything you previously believed turn out to be off?" invites the answer "no, it all made sense," and that answer is almost never true and always useless. The model names two concrete prior beliefs and what replaced them. Note the second one — reference-genome-as-a-person's-genome — is the exact misconception [Lecture 1 §5](../lecture-notes/01-the-central-dogma-in-90-minutes.md) and quiz Q10 are both aimed at, which is a small signal that the misconception is common enough to be worth a lecture section and a quiz question.

**Paragraph 4 is forward-looking and specific enough to act on.** "How annotation gets made" is a real gap that Week 1 genuinely does not cover — the week hands you GTF as a format without saying where its contents come from. A gap you can name is a gap you can close; "I want to learn more about bioinformatics" is not.

**The register is honest without being performative.** No self-flagellation, no triumph. This is the lab-notebook voice from the [Week-1 README](../README.md) applied to yourself: state what happened, state what you concluded, hedge where you are unsure. That is the same discipline as the Problem 3 essay, pointed inward.

### What a grader is looking for

The assignment says it plainly: *this is for you, not for a grade.* The honest description of what is being checked:

- **300–400 words.** The model is 348. This is the only hard criterion.
- **Four paragraphs, one per numbered question.** All four addressed.
- **Committed.**

Beyond that, the only meaningful standard is one you apply yourself: **would future-you, reading this after Week 12, learn anything from it?** If it is four paragraphs of "it was interesting and I learned a lot," the answer is no, and you have spent thirty minutes producing nothing. If it names two things you had wrong, the answer is yes.

There is no wrong answer to question 1 — if ethics was hardest for you and vocabulary was trivial, say so. The reflection is not being assessed against ours.

### Common wrong turns

**Writing what you think the course wants to hear.** "This week completely changed how I think about biology." Maybe it did; if so, say which part and how. If not, the sentence is worse than nothing, because it makes the document useless to you later.

**Answering question 2 with "no."** Occasionally true, usually a sign of not looking. Two candidates almost everyone gets wrong before Week 1: that a gene is by definition protein-coding, and that the reference genome is a real individual's sequence. Check whether you believed either.

**Skipping question 4 because you do not know what you do not know.** Question 4 is the most useful of the four precisely because it is hard. If nothing comes, look at the [Week-1 README topic list](../README.md) and find the topic that got one bullet where you wanted a section.

**Missing the word count in the low direction.** 300 words is about four solid paragraphs. If you are at 180, you answered the questions but did not give a reason for any of your answers — the reasons are where the reflection actually lives.

**Not committing it.** It is a graded acceptance criterion, and the mini-project copies this exact file into `week-01/reflection.md`. Uncommitted, it is not there when you need it in two days.

### How to verify

```sh
wc -w notes/week-01-reflection.md
grep -c '^\*\*[1-4]\.' notes/week-01-reflection.md
```

Expected output — word count in range, four numbered paragraphs:

```text
348 notes/week-01-reflection.md
4
```

Then commit:

```sh
git add notes/week-01-reflection.md
git commit -m "Week 1 reflection"
```

And confirm all six homework problems produced a commit, which is what the assignment preamble asked for:

```sh
git log --oneline | head -8
```

Expected output — something with this shape, one or more commits per problem:

```text
7f8e9d0 Week 1 reflection
6d7c8b9 Glossary v2: tightened wording, added vocabulary problem entry
5c6b7a8 Homework 4: tiny FASTA validator with summary and warnings
4b5a6c7 Re-identification essay on Gymrek et al. 2013
3a4b5c6 Hand translation of the first 60 nt of the spike CDS
2f3a4b5 Read the methods section of Liefferinckx et al. 2025
1e2f3a4 Week 1 public-data inventory
0d1e2f3 Week 1 glossary, first pass
```

---

*Continue to [mini-project.md](./mini-project.md).*
