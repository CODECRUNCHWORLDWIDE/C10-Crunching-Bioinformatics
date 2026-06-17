"""
Exercise 2 - Parse SAM by hand.

Goal: read a small SAM file as plain text (no pysam, no samtools), split
each alignment line into its eleven mandatory columns, decode the FLAG
field bit by bit, parse the CIGAR string into (length, op) operations,
and compute the implied query length and reference span. The point is
to demystify the SAM format: it really is just tab-separated text with
a clearly-defined column order, and you should be able to read one
without any library help.

Estimated time: 40 minutes. Pure Python, no BWA or samtools required.

Acceptance criteria:
- `python exercise-02-parse-sam-by-hand.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You implemented six functions: `parse_header`, `parse_alignment`,
  `decode_flag`, `parse_cigar`, `query_length_from_cigar`, and
  `reference_span_from_cigar`.

Requirements:
    Python 3.11+ (no extra packages).

What you learn:
- The eleven mandatory SAM columns by name and meaning.
- How to decode a 12-bit FLAG field by hand using bitwise AND.
- How to parse a CIGAR string into (length, op) tuples with a regex.
- Which CIGAR ops consume query bases (M, I, S, =, X) and which
  consume reference bases (M, D, N, =, X). The other two (H, P) consume
  neither.
- Why understanding the format pays off when pysam returns weird tuples:
  it is the same data, packed differently.

TO COMPLETE: implement the six functions below. Run the file; all
assertions must pass.

Tool versions assumed:
- Python 3.11+ (only the standard library is used)
"""

from __future__ import annotations

import re
from typing import NamedTuple


# A minimal but realistic SAM file. Two header lines plus six alignment
# records: a forward-strand uniquely-mapped primary, its mate (reverse-
# strand), a soft-clipped read, an unmapped read, a multimapper (MAPQ 0),
# and a duplicate. Together these cover every flag bit we care about.
EXAMPLE_SAM = """\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:NC_001416.1\tLN:48502
read001\t99\tNC_001416.1\t1001\t60\t150M\t=\t1351\t500\tACGTACGTACGTACGT\tFFFFFFFFFFFFFFFF\tNM:i:0\tMD:Z:150\tAS:i:150\tXS:i:0\tRG:Z:lambda
read001\t147\tNC_001416.1\t1351\t60\t150M\t=\t1001\t-500\tTGCATGCATGCATGCA\tFFFFFFFFFFFFFFFF\tNM:i:0\tMD:Z:150\tAS:i:150\tXS:i:0\tRG:Z:lambda
read002\t99\tNC_001416.1\t2001\t60\t10S140M\t=\t2351\t500\tACGTACGTACGTACGT\tFFFFFFFFFFFFFFFF\tNM:i:0\tMD:Z:140\tAS:i:140\tXS:i:0\tRG:Z:lambda
read003\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGT\tFFFFFFFFFFFFFFFF\tRG:Z:lambda
read004\t99\tNC_001416.1\t3001\t0\t150M\t=\t3351\t500\tACGTACGTACGTACGT\tFFFFFFFFFFFFFFFF\tNM:i:0\tMD:Z:150\tAS:i:150\tXS:i:150\tRG:Z:lambda
read005\t1123\tNC_001416.1\t4001\t60\t130M2I18M\t=\t4351\t500\tACGTACGTACGTACGT\tFFFFFFFFFFFFFFFF\tNM:i:2\tMD:Z:148\tAS:i:142\tXS:i:0\tRG:Z:lambda
"""


# SAM flag bit constants - the canonical 12-bit field.
FLAG_PAIRED = 0x1
FLAG_PROPER_PAIR = 0x2
FLAG_UNMAP = 0x4
FLAG_MUNMAP = 0x8
FLAG_REVERSE = 0x10
FLAG_MREVERSE = 0x20
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUP = 0x400
FLAG_SUPPLEMENTARY = 0x800


class HeaderLine(NamedTuple):
    """A parsed SAM header line. The `tag` is the two-char type after
    the '@' (e.g. 'HD', 'SQ', 'PG'). The `fields` is a dict mapping
    each key (e.g. 'SN', 'LN', 'VN') to its value as a string."""
    tag: str
    fields: dict


class Alignment(NamedTuple):
    """A parsed SAM alignment record (the eleven mandatory columns plus
    a list of optional tag strings)."""
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    tags: list[str]


# A regex that picks out (length, op) pairs in a CIGAR string. The op
# alphabet is M, I, D, N, S, H, P, =, X. The length is a positive integer.
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# Sets of CIGAR ops by what they consume. Memorize these.
CIGAR_CONSUMES_QUERY = set("MIS=X")
CIGAR_CONSUMES_REF = set("MDN=X")


def parse_header(line: str) -> HeaderLine:
    """Parse a single SAM header line (starts with '@').

    Header lines look like:
        @HD\\tVN:1.6\\tSO:coordinate
        @SQ\\tSN:NC_001416.1\\tLN:48502
        @PG\\tID:bwa\\tPN:bwa\\tVN:0.7.17\\tCL:bwa mem -t 4 ...

    The first token is '@<TAG>'; the rest are tab-separated 'KEY:VALUE'
    pairs.

    Returns:
        HeaderLine(tag='HD', fields={'VN': '1.6', 'SO': 'coordinate'})
    """
    assert line.startswith("@"), f"header line must start with @, got {line!r}"
    # TODO: strip whitespace, split on tabs, peel off '@<TAG>', then
    # split each remaining token on the first ':' into (key, value).
    raise NotImplementedError("Parse the header into HeaderLine")


def parse_alignment(line: str) -> Alignment:
    """Parse a single SAM alignment line (no leading '@').

    The line has at least 11 tab-separated columns; any extra columns
    are optional tag fields like 'NM:i:0', 'MD:Z:150', 'RG:Z:lambda'.

    Numeric columns: FLAG (col 2), POS (col 4), MAPQ (col 5),
    PNEXT (col 8), TLEN (col 9). The rest are strings.

    Returns an Alignment namedtuple with all eleven fields plus a list
    of extra tag strings.
    """
    # TODO: split on tabs. Verify at least 11 columns. Cast the five
    # int columns. Slice columns 12+ into the tags list.
    raise NotImplementedError("Parse the alignment line into Alignment")


def decode_flag(flag: int) -> dict[str, bool]:
    """Decode a SAM FLAG integer into a dict of boolean flags.

    Returns a dict with twelve keys:
        is_paired, is_proper_pair, is_unmapped, is_mate_unmapped,
        is_reverse, is_mate_reverse, is_read1, is_read2,
        is_secondary, is_qcfail, is_duplicate, is_supplementary.

    Use bitwise AND against the FLAG_* constants above. E.g.,
        is_paired = bool(flag & FLAG_PAIRED)
    """
    # TODO: build and return the twelve-key dict by AND-ing flag with
    # each FLAG_* constant and casting to bool.
    raise NotImplementedError("Decode the FLAG bit field")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """Split a CIGAR string into (length, op) tuples.

    Examples:
        parse_cigar('150M') -> [(150, 'M')]
        parse_cigar('10S140M') -> [(10, 'S'), (140, 'M')]
        parse_cigar('36M2I12M1D80M') -> [(36, 'M'), (2, 'I'), (12, 'M'),
                                          (1, 'D'), (80, 'M')]

    Unmapped reads have CIGAR '*'. Return an empty list in that case.
    """
    if cigar == "*":
        return []
    # TODO: use the CIGAR_RE regex to find all (length, op) pairs.
    # Cast length to int.
    raise NotImplementedError("Parse the CIGAR string")


def query_length_from_cigar(cigar: str) -> int:
    """Total number of bases this alignment consumes from the query.

    M, I, S, =, X consume query bases. D, N, H, P do not.

    Examples:
        query_length_from_cigar('150M') -> 150
        query_length_from_cigar('10S140M') -> 150
        query_length_from_cigar('36M2I12M1D80M') -> 36+2+12+0+80 = 130
        query_length_from_cigar('*') -> 0  (unmapped)
    """
    # TODO: sum length over the (length, op) pairs where op is in
    # CIGAR_CONSUMES_QUERY.
    raise NotImplementedError("Compute query length from CIGAR")


def reference_span_from_cigar(cigar: str) -> int:
    """Total number of bases this alignment spans on the reference.

    M, D, N, =, X consume reference bases. I, S, H, P do not.

    Examples:
        reference_span_from_cigar('150M') -> 150
        reference_span_from_cigar('10S140M') -> 140 (soft clip
                                                    does not consume ref)
        reference_span_from_cigar('36M2I12M1D80M') -> 36+12+1+80 = 129
        reference_span_from_cigar('*') -> 0  (unmapped)
    """
    # TODO: sum length over the (length, op) pairs where op is in
    # CIGAR_CONSUMES_REF.
    raise NotImplementedError("Compute reference span from CIGAR")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-02-parse-sam-by-hand.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # --- Test parse_cigar and its derived helpers. --------------------
    assert parse_cigar("150M") == [(150, "M")]
    assert parse_cigar("10S140M") == [(10, "S"), (140, "M")]
    assert parse_cigar("36M2I12M1D80M") == [
        (36, "M"), (2, "I"), (12, "M"), (1, "D"), (80, "M"),
    ]
    assert parse_cigar("*") == []

    assert query_length_from_cigar("150M") == 150
    assert query_length_from_cigar("10S140M") == 150
    assert query_length_from_cigar("36M2I12M1D80M") == 130
    assert query_length_from_cigar("130M2I18M") == 150
    assert query_length_from_cigar("*") == 0

    assert reference_span_from_cigar("150M") == 150
    assert reference_span_from_cigar("10S140M") == 140
    assert reference_span_from_cigar("36M2I12M1D80M") == 129
    assert reference_span_from_cigar("130M2I18M") == 148
    assert reference_span_from_cigar("*") == 0

    print("[exercise-02] CIGAR parsing: PASS")

    # --- Test decode_flag. --------------------------------------------
    # FLAG 99 = 0x40 + 0x20 + 0x2 + 0x1
    #        = READ1 + MREVERSE + PROPER_PAIR + PAIRED
    decoded = decode_flag(99)
    assert decoded["is_paired"] is True
    assert decoded["is_proper_pair"] is True
    assert decoded["is_unmapped"] is False
    assert decoded["is_reverse"] is False
    assert decoded["is_mate_reverse"] is True
    assert decoded["is_read1"] is True
    assert decoded["is_read2"] is False
    assert decoded["is_duplicate"] is False

    # FLAG 147 = 0x80 + 0x10 + 0x2 + 0x1
    #         = READ2 + REVERSE + PROPER_PAIR + PAIRED
    decoded = decode_flag(147)
    assert decoded["is_paired"] is True
    assert decoded["is_proper_pair"] is True
    assert decoded["is_reverse"] is True
    assert decoded["is_read1"] is False
    assert decoded["is_read2"] is True
    assert decoded["is_duplicate"] is False

    # FLAG 77 = 0x40 + 0x8 + 0x4 + 0x1
    #        = READ1 + MUNMAP + UNMAP + PAIRED (unmapped pair)
    decoded = decode_flag(77)
    assert decoded["is_paired"] is True
    assert decoded["is_unmapped"] is True
    assert decoded["is_mate_unmapped"] is True
    assert decoded["is_read1"] is True

    # FLAG 1123 = 0x400 + 0x40 + 0x2 + 0x1
    #          = DUP + READ1 + PROPER_PAIR + PAIRED
    # Wait: 1123 = 1024 + 64 + 32 + 2 + 1
    #             = DUP + READ1 + MREVERSE + PROPER_PAIR + PAIRED
    decoded = decode_flag(1123)
    assert decoded["is_duplicate"] is True
    assert decoded["is_paired"] is True
    assert decoded["is_proper_pair"] is True
    assert decoded["is_mate_reverse"] is True
    assert decoded["is_read1"] is True

    print("[exercise-02] FLAG decoding: PASS")

    # --- Test parse_header and parse_alignment on the example SAM. ----
    headers = []
    alignments = []
    for line in EXAMPLE_SAM.splitlines():
        if not line:
            continue
        if line.startswith("@"):
            headers.append(parse_header(line))
        else:
            alignments.append(parse_alignment(line))

    assert len(headers) == 2, f"expected 2 headers, got {len(headers)}"
    assert headers[0].tag == "HD"
    assert headers[0].fields["VN"] == "1.6"
    assert headers[0].fields["SO"] == "coordinate"
    assert headers[1].tag == "SQ"
    assert headers[1].fields["SN"] == "NC_001416.1"
    assert headers[1].fields["LN"] == "48502"

    assert len(alignments) == 6, f"expected 6 alignments, got {len(alignments)}"

    # Read 1: primary forward, FLAG 99, 150M.
    a = alignments[0]
    assert a.qname == "read001"
    assert a.flag == 99
    assert a.rname == "NC_001416.1"
    assert a.pos == 1001
    assert a.mapq == 60
    assert a.cigar == "150M"
    assert a.rnext == "="
    assert a.pnext == 1351
    assert a.tlen == 500
    assert a.tags[0].startswith("NM:i:")

    # Read 4: unmapped, FLAG 77.
    a = alignments[3]
    assert a.qname == "read003"
    assert a.flag == 77
    assert a.rname == "*"
    assert a.pos == 0
    assert a.mapq == 0
    assert a.cigar == "*"
    assert decode_flag(a.flag)["is_unmapped"] is True

    # Read 5: multimapper MAPQ 0, FLAG 99.
    a = alignments[4]
    assert a.qname == "read004"
    assert a.mapq == 0

    # Read 6: duplicate with insertion. FLAG 1123 = DUP + READ1 + MREVERSE +
    # PROPER_PAIR + PAIRED.
    a = alignments[5]
    assert a.qname == "read005"
    assert a.flag == 1123
    assert decode_flag(a.flag)["is_duplicate"] is True
    assert a.cigar == "130M2I18M"
    assert query_length_from_cigar(a.cigar) == 150
    assert reference_span_from_cigar(a.cigar) == 148

    print("[exercise-02] SAM parsing: PASS")

    # --- Summary. -----------------------------------------------------
    print()
    print("[exercise-02] Decoded 6 alignment records and 2 header lines.")
    print("[exercise-02] Reference contig: NC_001416.1 (48,502 bp).")
    print("[exercise-02] Records: 1 forward primary, 1 reverse mate,")
    print("[exercise-02]          1 soft-clipped, 1 unmapped,")
    print("[exercise-02]          1 multimapper (MAPQ 0), 1 duplicate.")
    print()
    print("[exercise-02] All assertions passed.")
    print("[exercise-02] You can now read a SAM file column by column")
    print("[exercise-02] and decode a FLAG bit by bit. Continue to")
    print("[exercise-02] exercise-03 (coverage plot).")
