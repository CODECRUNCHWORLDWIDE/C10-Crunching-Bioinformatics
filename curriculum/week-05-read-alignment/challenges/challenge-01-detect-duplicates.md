# Challenge 1 — Detect duplicates by hand

> **Estimated time:** 90 minutes.
> **Goal:** Implement a positional PCR-duplicate detector in pure Python (with `pysam` for BAM access) and compare its duplicate rate to the rate `samtools markdup` reports on the same BAM. Identify and explain at least two specific cases where your detector and `samtools markdup` disagree.

Lecture 2 §6 introduced duplicate marking as "samtools markdup detects duplicates by 5' alignment position and read orientation, marks them with flag `0x400`, and leaves them in the BAM." The one-sentence summary is correct but conceals about ten edge cases that the C source code (`bam_markdup.c`, ~700 lines) handles. This challenge has you implement the *core* of that algorithm — the 5'-position grouping — and then *measure your discrepancy* against samtools' production implementation. The point is not to match samtools' answer exactly. The point is to know what your detector misses and why.

## Background — the algorithm in three rules

Given a coordinate-sorted, mate-fixed BAM, **PCR duplicates** are read pairs that started life as the same original DNA fragment. After alignment, they appear at the same 5' reference position with the same strand orientation, **for both mates of the pair**. Concretely:

For each paired-end read pair `(R1, R2)` with primary alignments:

1. Compute the 5' reference position of `R1`:
   - If `R1` is on the forward strand (`0x10` FLAG bit unset), the 5' position is `read.reference_start` (0-based) minus any leading soft clip in the CIGAR.
   - If `R1` is on the reverse strand (`0x10` set), the 5' position is `read.reference_end` (one past the last aligned base) plus any trailing soft clip in the CIGAR. (For reverse-strand reads, the "5' end" of the original DNA molecule is the *right* end of the alignment.)
2. Compute the 5' reference position of `R2` the same way.
3. Build a duplicate key: `(contig, R1_5p, R1_strand, R2_5p, R2_strand)`. Group read pairs by this key. Any group of size ≥ 2 is a duplicate cluster; the highest-quality pair in the cluster is the *representative*, every other pair in the cluster is a *duplicate* and gets the `0x400` flag.

The "highest-quality" heuristic samtools uses is the sum of base qualities at the 5' end (a proxy for sequencing run quality at the position-determining bases). For this challenge, you can use the simpler heuristic of "first pair encountered in sorted order is the representative, rest are duplicates" — your duplicate *rate* will be the same; only *which specific pair gets the flag* differs.

```
detect_pcr_duplicates(bam_path):
    seen = {}                # key -> first pair encountered
    duplicates = []          # list of (qname, R1_or_R2) that are dups
    for read_pair in paired_reads_sorted(bam_path):
        key = duplicate_key(read_pair.r1, read_pair.r2)
        if key in seen:
            duplicates.append(read_pair.qname)
        else:
            seen[key] = read_pair.qname
    return len(duplicates), len(seen)
```

## Task

Implement the detector, run it on a real BAM, and write up the comparison.

### 1. The input BAM

Use the sorted+indexed BAM from **Exercise 1** (`exercises/aln/lambda.sorted.bam`). It has ~2000 reads from 1000 simulated read pairs against the lambda reference. The duplicate rate from simulation should be ~0% (wgsim does not produce duplicates), so this is a sanity-check dataset where both your detector and `samtools markdup` should report a near-zero duplicate rate.

For a more realistic test, also run on the **mini-project's** larger BAM (the SRR1770413 alignment), where the duplicate rate is typically 10–20%.

### 2. Implement and run

Write `challenges/detect_duplicates.py` in your portfolio repo. The script should:

1. Read a sorted+indexed BAM via `pysam.AlignmentFile`.
2. Iterate over primary, non-secondary, non-supplementary, mapped, paired reads.
3. Group reads by `QNAME` to assemble `(R1, R2)` pairs.
4. For each pair, compute the 5' position of each mate (accounting for soft clipping and strand).
5. Build the duplicate key `(contig, R1_5p, R1_strand, R2_5p, R2_strand)`.
6. Group pairs by key. Any group of size ≥ 2 is a duplicate cluster.
7. Report:
   - Total mapped read pairs.
   - Number of duplicate read pairs (pairs that are not the representative of their cluster).
   - Duplicate rate as a fraction.

### 3. Compare to `samtools markdup`

Run `samtools markdup -s` on the same BAM (after fixmate, as in Lecture 2 §6.1). Capture the `-s` summary output. Extract its reported duplicate count and rate.

Produce a comparison block:

```
                            | Your detector | samtools markdup |
----------------------------|---------------|------------------|
Total read pairs            |               |                  |
Duplicate read pairs        |               |                  |
Duplicate rate              |               |                  |
Difference (pairs)          |               |                  |
Difference (rate, pp)       |               |                  |
```

The differences are not bugs in your code (necessarily); they are the cases samtools handles that your simplified detector does not.

### 4. Write up

In `challenges/notes/duplicate-comparison.md` (200–400 words), answer:

- What was the duplicate rate from your hand-rolled detector? From `samtools markdup`?
- How many pairs were called duplicates by one tool but not the other? Look at the specific qnames in the disagreement set.
- For at least **two** specific qnames where your detector and samtools disagree, examine the raw alignment records (`samtools view aln.bam | grep <qname>`) and explain *why* samtools called it differently. Likely causes:
  - Soft-clipping correction: your detector may forget to add leading soft clip to the 5' position, off-by-one-vs-samtools.
  - Supplementary alignments: samtools considers supplementary alignment positions, your detector does not.
  - Unmapped mate: samtools handles "one mate mapped, one unmapped" pairs (which it usually does not call duplicates); your detector may skip them entirely.
  - Tie-breaking: samtools picks the representative by base quality at 5' end; you picked the first-seen. Different representative does *not* change duplicate counts but *does* change which qnames are flagged.
- Based on what you find, would you trust your detector for production use? Why or why not? In one sentence, what would have to change in your detector to match samtools' behavior?

## Acceptance criteria

- `python challenges/detect_duplicates.py <bam>` runs without crashing on any sorted+indexed BAM.
- Output includes the total read pairs, duplicate pairs, and duplicate rate as three numbers.
- A comparison block (table or list) shows your numbers next to `samtools markdup -s` numbers.
- `notes/duplicate-comparison.md` is 200–400 words and addresses all four prompts above with specific qnames.
- The script is < 200 lines (it should not be a re-implementation of samtools; just the core position-based grouping).

## Hints (do not peek for at least 20 minutes)

<details>
<summary>Hint 1 — How do I compute the 5' position with soft clipping?</summary>

For a forward-strand read (`0x10` flag unset), the 5' end on the reference is `read.reference_start - leading_soft_clip`, where `leading_soft_clip` is the length of the first CIGAR `S` operation if the first op is `S` (else 0). Negative positions are legal (the unclipped 5' end may fall *before* the reference's start).

For a reverse-strand read (`0x10` set), the 5' end is `read.reference_end + trailing_soft_clip`, where `trailing_soft_clip` is the length of the last CIGAR `S` operation if it is `S` (else 0).

In pysam:

```python
def five_prime_position(read):
    cigar = read.cigartuples or []
    if read.is_reverse:
        trailing_s = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
        return read.reference_end + trailing_s
    else:
        leading_s = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
        return read.reference_start - leading_s
```

(pysam's CIGAR op code 4 is `S` for soft clip.)

</details>

<details>
<summary>Hint 2 — How do I assemble (R1, R2) pairs from the BAM?</summary>

A coordinate-sorted BAM does not interleave mates — they can be far apart. The clean approach is to iterate twice, but for a 2000-read BAM you can do it in one pass and keep a dict:

```python
pending = {}  # qname -> first mate seen
pairs = {}    # qname -> (r1, r2) once both seen
for read in af.fetch(until_eof=True):
    if not read.is_paired or read.is_secondary or read.is_supplementary:
        continue
    if read.is_unmapped or read.mate_is_unmapped:
        continue
    if read.query_name in pending:
        first = pending.pop(read.query_name)
        # Order them: r1 is the read with is_read1=True.
        if first.is_read1:
            r1, r2 = first, read
        else:
            r1, r2 = read, first
        pairs[read.query_name] = (r1, r2)
    else:
        pending[read.query_name] = read
```

At end of file, `pending` holds reads whose mate did not appear (orphan reads — drop them or treat them as singletons; they are rare in healthy BAMs).

</details>

<details>
<summary>Hint 3 — How do I parse the samtools markdup -s output?</summary>

`samtools markdup -s` writes a stats block to stderr (not stdout). Capture stderr in subprocess:

```python
result = subprocess.run(
    ["samtools", "markdup", "-s", in_bam, out_bam],
    capture_output=True, text=True, check=True,
)
stats_text = result.stderr
```

The relevant lines look like:

```
COMMAND: samtools markdup -s ...
READ: 2000
WRITTEN: 2000
EXCLUDED: 0
EXAMINED: 2000
PAIRED: 2000
SINGLE: 0
DUPLICATE PAIR: 12
DUPLICATE SINGLE: 0
DUPLICATE PAIR OPTICAL: 0
DUPLICATE SINGLE OPTICAL: 0
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 12
DUPLICATE TOTAL: 12
ESTIMATED_LIBRARY_SIZE: ...
```

Parse with simple `for line in stats_text.splitlines(): if line.startswith("DUPLICATE PAIR:"): ...`.

The pair-level duplicate rate is `DUPLICATE PAIR / PAIRED / 2` (samtools counts both mates of a duplicate pair, so dividing by 2 gives pair-level rate).

</details>

<details>
<summary>Hint 4 — Why does my detector say zero duplicates on the wgsim BAM?</summary>

Because there are zero duplicates. `wgsim` generates reads from random positions across the reference; the probability of two pairs starting at exactly the same position is `(read_count / ref_length)^2`, which for 1000 pairs against 48.5 kb lambda is `(1000/48500)^2 ≈ 4 × 10^-4`. You will see 0–2 duplicates by chance. This is the *correct* answer for the wgsim dataset.

If you want to see your detector fire, generate a duplicate-rich dataset:

```bash
wgsim -N 100 -1 150 -2 150 -e 0.001 ref.fa r1a.fq r2a.fq
# Make 5 copies of the same data:
for i in 1 2 3 4; do
    cat r1a.fq >> reads/dup_R1.fq
    cat r2a.fq >> reads/dup_R2.fq
done
```

This produces a dataset with exactly 80% duplicates. Run your detector on the resulting BAM; you should report ~80% duplicate rate.

</details>

<details>
<summary>Hint 5 — Should I count the representative as a duplicate?</summary>

No. In a cluster of `k ≥ 2` pairs, only `k - 1` are duplicates (the cluster representative is not). Total duplicate count = `sum(size - 1 for cluster in clusters if size >= 2)` = `total_pairs - len(clusters)`. Get this off-by-one wrong and your duplicate rate will be 1/N higher than samtools' for every cluster.

</details>

## Stretch

If you finish under time and want more:

- Add **optical-duplicate** detection. Optical duplicates are positionally close on the flow cell, detected by parsing the qname (Illumina qnames look like `instrument:run:flowcell:lane:tile:x:y` — the last two fields are the flow-cell coordinates). Pairs within a 100-pixel Manhattan distance on the same tile are optical duplicates. Compare your optical-dup rate to `samtools markdup -d 100 -s`.
- Compute a **library complexity estimate** from your duplicate-cluster size distribution: if you see clusters of size 2, 3, 4, ..., fit a Lander-Waterman model to estimate the original library complexity. Picard's `ESTIMATED_LIBRARY_SIZE` is the same idea. See Wendl et al. 2009 for the math.
- Run your detector on a real human-genome BAM (e.g. the HG00096 sample from 1000 Genomes — just a 1 Mb region for tractability). Compare to the `samtools markdup` answer. Real BAMs expose edge cases that wgsim never produces.
- Read `bam_markdup.c` in the samtools source and identify the function that computes the 5'-end-with-soft-clip-correction. Compare it to your Python implementation.

## What you should be able to do after this

- Detect PCR duplicates from a BAM without any library help beyond a SAM reader.
- Explain in one paragraph what `samtools markdup` does, what your detector does, and where they disagree.
- Defend a duplicate rate in writing with reference to specific qnames and specific edge cases.
- Distinguish PCR duplicates (positional) from sequence-identical reads (which can be biologically real for very high-coverage targeted sequencing).

---

*Submit by committing `challenges/detect_duplicates.py` and `challenges/notes/duplicate-comparison.md` to your portfolio repo.*
