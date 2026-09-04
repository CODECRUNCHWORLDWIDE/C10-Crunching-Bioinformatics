# Week 1 — Challenge Solutions

Read this **after** you have made your own version run. Read [00-overview.md](./00-overview.md) first if you have not.

---

## Challenge 1 — Reverse complement, GC, and translate

**Assignment:** indexed at [../challenges/README.md](../challenges/README.md) as `challenge-01-reverse-complement.py`. **The stub file is not in the repository.** The specification below is taken from the two places the week describes it: the [Week-1 README](../README.md) ("`reverse_complement`, `GC_content`, `translate` from scratch") and [challenges/README.md](../challenges/README.md) ("No Biopython… Run your solution against the provided test cases at the bottom of the file").

**Assumed interface:** a module exporting `reverse_complement(seq)`, `GC_content(seq)`, and `translate(seq)`, with self-tests under `if __name__ == "__main__":`. [Homework Problem 2](../homework.md) additionally requires that running the file on a sequence prints its translation, so the `__main__` block takes an optional sequence argument.

The function name `GC_content` is not PEP 8 — a function should be `gc_content` — but it is the name the week's own index specifies, so that is the name we export. The module also binds `gc_content = GC_content` so both spellings work; if a later week's code calls one and your file defines the other, that alias is why nothing breaks.

### The complete answer

```python
"""Challenge 1 — reverse_complement, GC_content, translate. No Biopython."""

from __future__ import annotations

import sys

# Standard genetic code, keyed on RNA codons exactly as printed in Lecture 1 §4.
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
}


def reverse_complement(seq: str) -> str:
    """Reverse the sequence and complement each base. Case is preserved."""
    out = []
    for base in reversed(seq):
        try:
            out.append(COMPLEMENT[base])
        except KeyError:
            raise ValueError(f"not a DNA base: {base!r}") from None
    return "".join(out)


def GC_content(seq: str) -> float:
    """Percent G+C over unambiguous bases. N is excluded from the denominator."""
    counted = 0
    gc = 0
    for base in seq.upper():
        if base in "ACGT":
            counted += 1
            if base in "GC":
                gc += 1
        elif base == "N":
            continue
        else:
            raise ValueError(f"not a DNA base: {base!r}")
    if counted == 0:
        return 0.0
    return 100.0 * gc / counted


gc_content = GC_content   # PEP 8 alias; same function


def transcribe(dna: str) -> str:
    """DNA to RNA: T becomes U. Uppercases, because CODON_TABLE is uppercase."""
    return dna.upper().replace("T", "U")


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate a DNA sequence in frame 1. Stops are '*'.

    A trailing partial codon (1 or 2 leftover bases) is silently dropped,
    which is what every standard tool does.
    """
    rna = transcribe(seq)
    peptide = []
    for i in range(0, len(rna) - len(rna) % 3, 3):
        codon = rna[i:i + 3]
        try:
            aa = CODON_TABLE[codon]
        except KeyError:
            raise ValueError(f"unknown codon {codon!r} at position {i}") from None
        if aa == "*" and to_stop:
            break
        peptide.append(aa)
    return "".join(peptide)


# First 60 nt of the SARS-CoV-2 spike CDS — the Homework Problem 2 sequence.
SPIKE_60 = "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACC"


def _self_test() -> None:
    # reverse_complement
    assert reverse_complement("ATGCGT") == "ACGCAT"          # quiz Q9
    assert reverse_complement("") == ""                       # empty is legal
    assert reverse_complement("AAAA") == "TTTT"               # homopolymer
    assert reverse_complement("ATGCNN") == "NNGCAT"           # N maps to N
    assert reverse_complement(reverse_complement("ACGTTGCA")) == "ACGTTGCA"

    # GC_content
    assert GC_content("GGCC") == 100.0
    assert GC_content("ATAT") == 0.0
    assert GC_content("ATGC") == 50.0
    assert GC_content("ATGCNN") == 50.0                       # N not in denominator
    assert GC_content("") == 0.0                              # no division by zero

    # translate
    assert translate("ATGGCCTAA") == "MA*"
    assert translate("ATGGCCTAA", to_stop=True) == "MA"
    assert translate("ATGGCCTAAG") == "MA*"                   # partial codon dropped
    assert translate(SPIKE_60) == "MFVFLVLLPLVSSQCVNLTT"      # homework 2

    print("all self-tests passed")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        s = sys.argv[1].strip().upper()
        print(f"input          {s}")
        print(f"length         {len(s)} nt")
        print(f"GC content     {GC_content(s):.2f}%")
        print(f"revcomp        {reverse_complement(s)}")
        print(f"protein        {translate(s)}")
    else:
        _self_test()
```

### Why it works

**`reverse_complement` does both operations in one pass, and the order does not matter.** Iterating `reversed(seq)` while complementing each base is the same result as complementing the whole string and then reversing it — the quiz's Q9 explanation makes the same point. The reason the operations commute is that complementing is *positionwise* (each base's image depends only on that base) while reversing is *positional* (it permutes indices without touching values). Independent transforms commute. Doing them together is one traversal and one allocation instead of two of each.

The biology under it: DNA is anti-parallel, and the convention from [Lecture 1 §2](../lecture-notes/01-the-central-dogma-in-90-minutes.md) is that sequences are written 5'→3'. If the top strand reads `ATGCGT` 5'→3', the bottom strand pairs base-for-base as `TACGCA` — but running the *other* direction. To write the bottom strand in the same 5'→3' convention you have to read it backwards, which is exactly the reverse of the complement. That is the entire reason this function exists and why it is the first thing every sequence library implements.

`N` is in the complement table because real reference FASTA is full of it — GRCh38 has tens of millions of `N` bases in centromeres and telomeres. A function that raises on `N` cannot be pointed at a real chromosome. Lowercase is in the table because reference genomes use soft-masking: repeat regions are written in lowercase, and losing the case throws away information that Week 4's BLAST and Week 5's aligners both consume.

**`GC_content` excludes `N` from the denominator, and that choice is the interesting part.** There are three defensible definitions of GC% on a sequence with ambiguity, and they disagree: G+C over the full length, G+C over unambiguous bases only, or G+C with ambiguity codes counted fractionally. We take the middle one, because it answers the question people actually mean — "of the bases we can see, what fraction are G or C?" — and because it does not silently report a lower GC% for a sequence that happens to have a masked region in it. What matters is not which definition you pick but that you **write it down**; the docstring states it, and Week 10's assembly QC will hand you a tool that made a different choice.

Returning `0.0` for an all-`N` or empty sequence rather than raising is the pragmatic call for a summary statistic: you will run this over every record in a FASTA file, and one empty record should not abort the run. The alternative — returning `None` or `float("nan")` — is arguably more honest but forces every caller to special-case it.

**`translate` transcribes first, then looks up.** That is one extra pass over the string, and it is deliberate. The codon table in [Lecture 1 §4](../lecture-notes/01-the-central-dogma-in-90-minutes.md) is printed in RNA (`UUU Phe`, not `TTT Phe`), so keying the dict on RNA codons means you can check the code against the lecture line by line without mentally substituting 64 times. Transcription is `T → U` and nothing else, so the mapping is exact and reversible. When you are typing 64 dictionary entries by hand, matching the source table character for character is worth one `O(n)` pass.

The loop bound `len(rna) - len(rna) % 3` is the whole partial-codon story in one expression. `SPIKE_60` is 60 nucleotides, `60 % 3 == 0`, so it yields 20 codons. Give it 61 and the bound becomes 60, so the trailing base is dropped rather than producing a `KeyError` on a two-character "codon". That is what Biopython, EMBOSS, and every other tool do; the difference is that Biopython emits `BiopythonWarning: Partial codon, len(sequence) not a multiple of three` while ours is silent. Ours is documented, which is the minimum bar.

`to_stop=False` by default so that `'*'` appears in the output. That is the right default for a teaching function: you want to *see* that `UAA` is a stop, not have it silently truncate your protein. `translate("ATGGCCTAA")` returning `"MA*"` tells you where the reading frame ended; returning `"MA"` does not.

### Comparing against Biopython, once you are done

[challenges/README.md](../challenges/README.md) asks you to compare your implementation to Biopython afterwards. Do it — but expect these differences, all of which are Biopython being more general rather than you being wrong:

| Behaviour | This implementation | `Bio.Seq.Seq` |
|---|---|---|
| Ambiguity codes beyond `N` (`R`, `Y`, `S`, `W`, `K`, `M`, …) | `ValueError` | Complemented correctly per IUPAC |
| Non-standard genetic codes | Standard table only | `translate(table=2)` etc., via `Bio.Data.CodonTable` |
| Partial trailing codon | Dropped silently | Dropped with a `BiopythonWarning` |
| Stop character | `*` | `*` — same convention |
| `to_stop` | Supported, same meaning | Supported, same meaning |

The one difference worth internalising is the first row. IUPAC defines a full ambiguity alphabet, and its complements are not identities: `R` (A or G) complements to `Y` (C or T), not to `R`. If you ever extend `COMPLEMENT`, that is the trap — you cannot map an ambiguity code to itself the way you can with `N`.

### Common wrong turns

**Reversing with `seq[::-1]` after building the complement with `str.replace`.**

```python
def reverse_complement(seq):                 # WRONG
    seq = seq.replace("A", "T")
    seq = seq.replace("T", "A")
    seq = seq.replace("G", "C")
    seq = seq.replace("C", "G")
    return seq[::-1]
```

This is the classic. `"ATGCGT"` returns `"ACACAT"` instead of `"ACGCAT"`. The first `replace` turns every `A` into `T`, and then the second `replace` turns *those brand-new `T`s* back into `A`. Sequential replacement is not a simultaneous substitution. `seq[::-1]` itself is fine — the bug is upstream. If you want the slice-based version, do the substitution simultaneously with `str.maketrans`:

```python
_TBL = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def reverse_complement(seq):
    return seq.translate(_TBL)[::-1]         # correct, and faster
```

`str.translate` builds the whole output in C in one pass and is the fastest pure-Python option by a wide margin. It is worth knowing after you have written the explicit loop once — which is the point of the challenge.

**Building the codon table in DNA and then feeding it RNA (or the reverse).** Symptom: a `KeyError` on the very first codon.

```text
ValueError: unknown codon 'ATG' at position 0
```

That message means your table is keyed on RNA and you skipped `transcribe`, or it is keyed on DNA and you called `transcribe` anyway. Pick one representation, put it in the docstring, and make the tests cover it.

**Iterating with `range(0, len(seq), 3)` instead of subtracting the remainder.** On a sequence whose length is not a multiple of three you get:

```text
ValueError: unknown codon 'AG' at position 60
```

A two-character key is not in a 64-entry table. The fix is the loop bound, not a `try`/`except` around the lookup — swallowing the `KeyError` would also hide a genuine typo in your table.

**Counting GC with `seq.count("G") + seq.count("C")` on a mixed-case sequence.** `str.count` is case-sensitive, so a soft-masked reference region contributes zero. Symptom: GC% that is implausibly low on exactly the repeat-rich regions where you would expect it to be high. Uppercase first, as the code above does.

**A typo in one of the 64 table entries.** This is the one that will actually cost you an hour, because 63 codons work and one silently produces the wrong amino acid. Two defences, both cheap: check the table's length (`len(CODON_TABLE) == 64`), and check that the values contain exactly 21 distinct symbols (20 amino acids plus `*`). Both are one-liners in the verification section below.

### How to verify

**Self-tests.** From your week-01 working directory:

```sh
python challenges/challenge-01-reverse-complement.py
```

Expected output, exactly:

```text
all self-tests passed
```

**On the Homework Problem 2 sequence.** This is the invocation Homework Problem 2's acceptance criteria ask for:

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

`MFVFLVLLPLVSSQCVNLTT` is the string Homework Problem 2 tells you to expect, so this run is simultaneously the challenge's test and the homework's cross-check.

**Table integrity.** Catch a typo in the 64 entries before it costs you an afternoon:

```sh
python -c "import importlib.util,sys; s=importlib.util.spec_from_file_location('c','challenges/challenge-01-reverse-complement.py'); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); print(len(m.CODON_TABLE), len(set(m.CODON_TABLE.values())), sorted(set(''.join(m.CODON_TABLE))))"
```

Expected output:

```text
64 21 ['A', 'C', 'G', 'U']
```

Sixty-four codons, twenty-one distinct output symbols (20 amino acids plus `*`), and every key built from `A`, `C`, `G`, `U` only — no stray `T` from a copy-paste out of a DNA table.

---

*Continue to [homework.md](./homework.md).*
