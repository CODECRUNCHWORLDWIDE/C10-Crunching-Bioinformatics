"""
Challenge 1 — Reverse complement, GC content, and translate.

Goal: Implement three classic operations on DNA sequences from scratch.
No Biopython. No third-party packages. Standard library only.

Estimated time: 60 minutes.

Why we do this: these three operations show up in every bioinformatics
program ever written. Biopython has them; samtools has them; your own
code in Week 6 will have them. Writing them once, by hand, on a quiet
afternoon means you will *understand* the corner cases (what happens
with N's? with lowercase? with U's?) for the rest of your career.

Acceptance criteria:
- `python challenge-01-reverse-complement.py` runs without crashing.
- All three functions below pass the `assert` checks at the bottom.
- Your implementation does NOT import Biopython, regex, or any non-stdlib package.

What you learn:
- DNA complementarity at the level of a `str.maketrans` table.
- GC content as a fraction over A/C/G/T (excluding ambiguity codes).
- The genetic code as a Python dict, and how `translate()` walks a sequence in steps of 3.
- The reading-frame question — translate from position 0, 1, or 2? Forward strand or reverse?

TO COMPLETE: implement `reverse_complement`, `gc_content`, and `translate`.
Do not look at the hint until you have tried at least 20 minutes.
"""

from __future__ import annotations


# The standard genetic code (NCBI codon table 1).
# Use this as the source of truth. Stop codons map to '*'.
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Complementarity rules:
        A <-> T,  C <-> G,  N <-> N

    The input may be mixed-case; the output should be uppercase.
    The input may contain N's (any unknown base); they map to N.
    Any character not in {A, C, G, T, N} (case-insensitive) is an error
    — raise ValueError with a helpful message.

    Examples:
        >>> reverse_complement("ACGT")
        'ACGT'
        >>> reverse_complement("AAAACCCGGT")
        'ACCGGGTTTT'
        >>> reverse_complement("NNNATGC")
        'GCATNNN'
    """
    # TODO: implement this function. A clean solution is ~5 lines using
    # str.maketrans and slicing.
    raise NotImplementedError("Implement reverse_complement")


def gc_content(seq: str) -> float:
    """Return the GC content of a DNA sequence as a percentage (0.0–100.0).

    GC content is defined as (G + C) / (A + C + G + T), counted
    case-insensitively. N's and other ambiguity codes are excluded
    from BOTH numerator and denominator — they neither help nor hurt.

    An empty sequence, or a sequence containing only N's, returns 0.0.

    Examples:
        >>> gc_content("ATGC")
        50.0
        >>> gc_content("GGGG")
        100.0
        >>> gc_content("AAAA")
        0.0
        >>> gc_content("NNNN")
        0.0
        >>> gc_content("")
        0.0
    """
    # TODO: implement this function. Iterate once over the sequence;
    # count A/C/G/T explicitly. Do NOT use Bio.SeqUtils.gc_fraction.
    raise NotImplementedError("Implement gc_content")


def translate(seq: str, frame: int = 0) -> str:
    """Translate a DNA sequence into a protein sequence.

    Reads the sequence in codons (3 nucleotides at a time) starting at
    position `frame` (0, 1, or 2). Maps each codon to its amino acid
    using the standard genetic code (CODON_TABLE above). Stop codons
    are represented as '*'.

    A trailing 1 or 2 nucleotides that do not complete a codon are
    silently dropped.

    Any codon containing an N (or any non-ACGT character) translates
    to 'X' (unknown amino acid) — do not raise.

    Examples:
        >>> translate("ATGGCCTGA")
        'MA*'
        >>> translate("ATGGCCTGAA")
        'MA*'
        >>> translate("ATGGCC", frame=0)
        'MA'
        >>> translate("AATGGCC", frame=1)
        'MA'
        >>> translate("ATGNNN")
        'MX'
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1, or 2; got {frame}")
    # TODO: implement this function. Slice `seq[frame:]`, iterate in
    # steps of 3, look each codon up in CODON_TABLE; codons containing
    # any non-standard base map to 'X'.
    raise NotImplementedError("Implement translate")


# ----------------------------------------------------------------------
# Self-test. Run this file directly with
#   python challenge-01-reverse-complement.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # reverse_complement
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"
    assert reverse_complement("a") == "T"
    assert reverse_complement("NNNATGC") == "GCATNNN"
    assert reverse_complement("") == ""
    try:
        reverse_complement("ACGU")  # RNA, not DNA — should raise
    except ValueError:
        pass
    else:
        raise AssertionError("reverse_complement must reject non-DNA characters")

    # gc_content
    assert gc_content("ATGC") == 50.0
    assert gc_content("GGGG") == 100.0
    assert gc_content("AAAA") == 0.0
    assert gc_content("NNNN") == 0.0
    assert gc_content("") == 0.0
    assert abs(gc_content("ACGTACGTNN") - 50.0) < 1e-9
    assert abs(gc_content("aaaagggg") - 50.0) < 1e-9  # case-insensitive

    # translate
    assert translate("ATGGCCTGA") == "MA*"
    assert translate("ATGGCCTGAA") == "MA*"           # trailing partial codon dropped
    assert translate("ATGGCC") == "MA"
    assert translate("AATGGCC", frame=1) == "MA"      # frame shift
    assert translate("ATGNNN") == "MX"                # N-containing codon
    assert translate("") == ""

    # The SARS-CoV-2 spike protein starts MFVFLVLLPLVSSQCV...
    # The first 51 nt of the spike CDS encode the first 17 amino acids:
    spike_first_51_nt = "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAAT"
    assert translate(spike_first_51_nt) == "MFVFLVLLPLVSSQCVN", (
        f"got {translate(spike_first_51_nt)!r}"
    )

    print("All assertions passed.")
    print("Now open Bio/Seq.py in the Biopython source and compare.")
