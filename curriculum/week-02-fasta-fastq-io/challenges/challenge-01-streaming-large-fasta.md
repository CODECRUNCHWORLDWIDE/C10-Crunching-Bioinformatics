# Challenge 1 — Stream a Large FASTA

> **Estimated time:** 75 minutes.
> **Goal:** Compute three summary statistics over a multi-gigabyte gzipped FASTA file **in a single pass**, with **bounded memory** (peak RSS ≤ 150 MB regardless of input size), using only `Bio.SeqIO.parse` (or the lower-level `SimpleFastaParser`) and the Python standard library.

Every Week-2 exercise so far has loaded its FASTA into a list. That works because the inputs were small. It will not work when you point your script at GRCh38 (3.0 GB), the human transcriptome (1.6 GB), or a SARS-CoV-2 assembly catalog (4 GB and growing). This challenge builds the IO muscle you will lean on for the rest of C10.

## Background

The canonical pattern is:

```python
import gzip
from Bio import SeqIO

with gzip.open("genome.fasta.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # process exactly one record, accumulate scalars, discard the record
        pass
```

The discipline: **never hold more than one record's worth of sequence in memory at a time**. Accumulators must be scalars (counts, sums) or fixed-size objects (histograms, top-K heaps). The moment you write `records = list(...)` or `all_seqs = []` followed by `all_seqs.append(record.seq)`, you have broken the pattern.

If the input file has 100 million records and the average record is 200 bp, the peak memory of a *correct* streaming solution is ~1 KB regardless. Memory should be O(1) in record count.

## Task

Write `stream_fasta_stats.py` that:

1. Accepts a file path on the command line. The path may end in `.fasta`, `.fa`, `.fasta.gz`, or `.fa.gz` — your script must handle gzipped input transparently.
2. Streams the file end to end, **without ever materializing more than one record at a time**.
3. Computes and prints, on stdout:
   - **`n_records`** — total record count.
   - **`total_bp`** — sum of all sequence lengths (count of `A`, `C`, `G`, `T`, `N`, and any other letter; treat the sequence as opaque length).
   - **`gc_percent`** — overall GC% computed across only `A`, `C`, `G`, `T` (case-insensitive). Exclude `N` from both numerator and denominator.
   - **`length_histogram`** — a small histogram of sequence lengths binned into the buckets:
     `[0, 50)`, `[50, 100)`, `[100, 500)`, `[500, 1000)`, `[1000, 10000)`, `[10000, 100000)`, `[100000, +infinity)`.
   - **`longest_id`** — the `record.id` of the single longest record. Tie-break: first one seen.
4. Prints output as a JSON object so it is machine-parseable.

Example output for a 12-record toy file:

```
{
  "n_records": 12,
  "total_bp": 1842031,
  "gc_percent": 41.62,
  "length_histogram": {
    "[0,50)": 0,
    "[50,100)": 1,
    "[100,500)": 3,
    "[500,1000)": 4,
    "[1000,10000)": 3,
    "[10000,100000)": 1,
    "[100000,+inf)": 0
  },
  "longest_id": "chr22:fragment_47"
}
```

## Acceptance criteria

- `python stream_fasta_stats.py <path>` runs against both gzipped and uncompressed inputs.
- The implementation **does not** call `list()`, `dict()`, `to_dict`, or any other constructor that materializes the full record set.
- Peak RSS is independent of file size. Verify on a small file (1 MB) and a medium file (300 MB+) and confirm RSS is comparable. On Linux: `/usr/bin/time -v python stream_fasta_stats.py <path>`; on macOS: `/usr/bin/time -l python stream_fasta_stats.py <path>`. Report the peak in MB in your `notes/`.
- Output is valid JSON. `python -c "import json,sys; json.load(open('out.json'))"` parses it cleanly.
- The script handles an empty file (zero records) without crashing — return all counters at 0 and `longest_id: null`.

## Acceptance test data

For self-testing without downloading multi-gigabyte data:

1. **Toy input.** Use the SARS-CoV-2 reference genome (NC_045512.2). Download:

   ```bash
   curl -L 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta' \
       -o sars2.fasta
   ```

   Expected: one record, ~29,903 bp, GC ~38%. Confirm.

2. **Medium input.** The human transcriptome from Ensembl, release 113:

   ```bash
   curl -L 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' \
       -o transcripts.fa.gz
   ```

   ~110,000 records, ~250 MB compressed, ~1.6 GB uncompressed. Your script must process it *without* uncompressing to disk.

3. **Tiny synthetic input.** Build a 5-record FASTA in `/tmp` for unit testing. (Hint: write your own generator and pipe it into the script with shell redirection or test on a temp file.)

Compare your script's output to `seqkit stats <path>` on the same input — your `total_bp` and `n_records` must agree to the digit.

## Hints (do not peek for at least 30 minutes)

<details>
<summary>Hint 1 — How do I open both .gz and uncompressed paths transparently?</summary>

```python
import gzip
from pathlib import Path

def open_maybe_gz(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")
```

Use it in a `with` block. The returned object is a text-mode file handle in either case.
</details>

<details>
<summary>Hint 2 — What goes in the inner loop?</summary>

A scalar accumulator pattern. For example:

```python
n = 0
total_bp = 0
ac = gc = at = 0  # counts of A+T, G+C
buckets = {...}    # the histogram
longest_id = None
longest_len = -1

for record in SeqIO.parse(handle, "fasta"):
    n += 1
    seq = str(record.seq).upper()
    L = len(seq)
    total_bp += L
    # GC counting; do it inline rather than building a Counter (faster).
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")
    gc += g + c
    at += a + t
    # histogram bucket lookup
    bucket = bucket_for(L)
    buckets[bucket] += 1
    if L > longest_len:
        longest_len = L
        longest_id = record.id
```

This loop is O(record-size) in time and O(1) in extra memory.
</details>

<details>
<summary>Hint 3 — How do I confirm I am not leaking memory?</summary>

Use `tracemalloc`:

```python
import tracemalloc
tracemalloc.start()
# ... run the loop ...
current, peak = tracemalloc.get_traced_memory()
print(f"Python peak: {peak/1024/1024:.1f} MB", file=sys.stderr)
```

Or measure RSS from outside the process with `/usr/bin/time -l` (macOS) / `-v` (Linux). The RSS number is what actually matters; `tracemalloc` only counts Python allocations.

If your peak is growing with input size, look for: building a list of records, storing whole sequences for later, using `Counter` over the concatenated sequence, or any code path that calls `to_dict`.
</details>

<details>
<summary>Hint 4 — Edge cases I should handle.</summary>

- Empty file → emit valid JSON with zeros and `null` for `longest_id`.
- A record with an empty sequence → length 0, lands in `[0,50)`, contributes nothing to GC.
- Mixed case → upper-case the sequence once per record (don't iterate character by character in Python).
- Non-ACGTN characters (rare but possible: IUPAC ambiguity codes `R`, `Y`, `K`, etc.) → count toward `total_bp` but not toward A/C/G/T.
</details>

## Stretch

If you finish under time and want more:

- Make the script accept a list of paths and aggregate across them, producing one JSON object.
- Add a `--by-record` flag that emits one JSON line per record (length, GC%, id) instead of a single aggregate. This is the streaming JSON-lines pattern used in production logging.
- Profile your script with `python -m cProfile` on the medium input. Identify the bottleneck (it will almost certainly be `str(record.seq)` and `.count()` calls). Try dropping to `SimpleFastaParser` and benchmark again.

## What you should be able to do after this

- Process production-scale FASTA without your laptop swapping.
- Know when to use `SeqIO.parse` (streaming) vs `SeqIO.index` (random access) vs `SeqIO.to_dict` (small in-memory).
- Read code that uses these patterns in samtools, seqkit, or Biopython itself — you will recognize the shape.

---

*Submit by committing `stream_fasta_stats.py` plus a `notes/streaming-rss.md` recording the peak RSS you measured on the toy and medium inputs.*
