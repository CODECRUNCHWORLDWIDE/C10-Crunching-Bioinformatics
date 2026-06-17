"""
Exercise 2 - Parse a VCF by hand.

Goal: read a small VCF file as plain text (no pysam, no bcftools), split
each variant line into its eight mandatory columns plus FORMAT and per-
sample columns, decode the INFO field into a dict of KEY -> VALUE, decode
the FORMAT and per-sample column into a dict of KEY -> VALUE, and answer
basic questions: how many SNPs, how many indels, how many PASS, what is
the mean QUAL across PASS variants. The point is to demystify VCF: it
really is just tab-separated text with a clearly-defined column order,
and you should be able to read one without any library help.

Estimated time: 45 minutes. Pure Python, no bcftools required.

Acceptance criteria:
- `python exercise-02-vcf-parse.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You implemented six functions: `parse_header_line`, `parse_variant_line`,
  `parse_info`, `parse_sample`, `is_snp`, and `is_indel`.

Requirements:
    Python 3.11+ (no extra packages).

What you learn:
- The eight mandatory VCF columns by name and meaning.
- How to decode the INFO field (semicolon-separated KEY=VALUE pairs).
- How to decode the FORMAT and per-sample columns (colon-separated
  KEY:KEY:... and VALUE:VALUE:... in parallel).
- How to distinguish SNPs (REF and all ALTs of length 1) from indels
  (any allele of length > 1).
- Which INFO keys come from bcftools (DP, MQ, SP, AF, INDEL) vs from
  GATK (QD, FS, MQ, MQRankSum, ReadPosRankSum, SOR).
- Why understanding the format pays off when pysam.VariantFile returns
  a weird object: it is the same data, packed differently.

TO COMPLETE: implement the six functions below. Run the file; all
assertions must pass.

Tool versions assumed:
- Python 3.11+ (only the standard library is used).
"""

from __future__ import annotations

from typing import NamedTuple


# A minimal but realistic VCF excerpt. The header lines define the
# metadata; the variant records cover a SNP, two indels (an insertion
# and a deletion), a PASS variant, a LowQual variant, and a multiallelic
# (already split by bcftools norm -m -any).
EXAMPLE_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Failed bcftools filter expression">
##contig=<ID=NC_001416.1,length=48502>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##INFO=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele frequency">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Per-sample depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tlambda_sample
NC_001416.1\t1001\t.\tA\tG\t227.0\tPASS\tDP=42;MQ=60;SP=2;AF=1\tGT:PL:DP:AD\t1:255,0:42:0,42
NC_001416.1\t2003\trs0001\tT\tC\t198.0\tPASS\tDP=38;MQ=60;SP=3;AF=1\tGT:PL:DP:AD\t1:255,0:38:0,38
NC_001416.1\t5021\t.\tA\tAGT\t156.0\tPASS\tINDEL;DP=30;MQ=60;SP=8;AF=1\tGT:PL:DP:AD\t1:255,0:30:0,30
NC_001416.1\t12450\t.\tCAA\tC\t102.0\tPASS\tINDEL;DP=25;MQ=58;SP=12;AF=1\tGT:PL:DP:AD\t1:255,0:25:0,25
NC_001416.1\t20100\t.\tG\tA\t22.0\tLowQual\tDP=4;MQ=60;SP=1;AF=1\tGT:PL:DP:AD\t1:50,0:4:0,4
NC_001416.1\t30500\t.\tT\tA\t15.0\tLowQual\tDP=8;MQ=22;SP=70;AF=1\tGT:PL:DP:AD\t1:40,0:8:0,8
NC_001416.1\t40000\t.\tA\tG\t245.0\tPASS\tDP=45;MQ=60;SP=4;AF=1\tGT:PL:DP:AD\t1:255,0:45:0,45
"""


# The eight mandatory VCF columns plus FORMAT.
VCF_COLUMNS = [
    "CHROM",  # 0
    "POS",    # 1 (int)
    "ID",     # 2
    "REF",    # 3
    "ALT",    # 4 (comma-separated; we keep as a string here, split later)
    "QUAL",   # 5 (float; '.' if missing)
    "FILTER", # 6 (semicolon-separated; '.' if missing)
    "INFO",   # 7 (semicolon-separated KEY=VALUE pairs)
    "FORMAT", # 8 (colon-separated keys; only if samples are present)
]


class VariantRecord(NamedTuple):
    """A parsed VCF variant line."""
    chrom: str
    pos: int
    id: str
    ref: str
    alts: list[str]
    qual: float | None
    filter: list[str]
    info: dict[str, str]
    format: list[str]
    samples: dict[str, dict[str, str]]


def parse_header_line(line: str) -> tuple[str, str] | tuple[str, dict]:
    """Parse a single VCF header line (starts with '##' or '#').

    Two kinds of header line:

    1. `##KEY=VALUE` — simple key/value (e.g. `##fileformat=VCFv4.2`).
       Returns ('KEY', 'VALUE') as a 2-tuple.

    2. `##KEY=<ID=...,Number=...,Type=...,Description="...">` — structured.
       Returns ('KEY', {'ID': ..., 'Number': ..., 'Type': ..., 'Description': ...})
       as a 2-tuple of (str, dict).

    3. `#CHROM\\tPOS\\tID\\t...` — the column header line.
       Returns ('COLUMNS', [list, of, column, names]) as a 2-tuple of
       (str, list).

    Hints:
    - Strip leading '##' and trailing whitespace.
    - For the structured form, peel off the leading '<' and trailing '>',
      then split on commas — but watch out for commas inside quoted
      Description fields. A simple regex `([A-Za-z]+)=("[^"]*"|[^,]*)`
      works.
    - For '#CHROM\\t...', strip the leading '#' and split on tabs.
    """
    assert line.startswith("#"), f"header line must start with #, got {line!r}"
    # TODO: detect which of the three forms this is and parse accordingly.
    raise NotImplementedError("Parse the VCF header line")


def parse_variant_line(line: str, sample_names: list[str]) -> VariantRecord:
    """Parse a single VCF variant line (no leading '#').

    The line has at least 8 tab-separated mandatory columns; if there
    are samples, there is also a FORMAT column (9) and one column per
    sample (10+).

    Numeric columns: POS (int), QUAL (float, '.' means None).

    Returns a VariantRecord namedtuple with all the parsed fields.

    Hint: split on tabs. Cast POS to int. Cast QUAL to float (or None
    if '.'). Split ALT on commas to get a list. Split FILTER on
    semicolons (or wrap in a list if a single value). Parse INFO with
    parse_info. Parse FORMAT into a list of keys. Parse each per-sample
    column with parse_sample.
    """
    # TODO: split on tabs. Validate column count.
    # TODO: cast POS to int. Cast QUAL to float or None.
    # TODO: split ALT on commas.
    # TODO: parse FILTER, INFO, FORMAT, and per-sample columns.
    # TODO: return a VariantRecord.
    raise NotImplementedError("Parse the variant line into VariantRecord")


def parse_info(info: str) -> dict[str, str]:
    """Parse a VCF INFO column into a dict of KEY -> VALUE.

    INFO is a semicolon-separated list of KEY=VALUE pairs. Some entries
    are flags with no value (e.g. 'INDEL') — represent these as
    KEY -> '' (empty string).

    Examples:
        parse_info('DP=42;MQ=60;SP=2;AF=1')
            -> {'DP': '42', 'MQ': '60', 'SP': '2', 'AF': '1'}
        parse_info('INDEL;DP=30;MQ=60')
            -> {'INDEL': '', 'DP': '30', 'MQ': '60'}
        parse_info('.')
            -> {}  (missing/empty INFO)

    All values are returned as strings; type conversion is the caller's
    job (different fields have different types per the header definition).
    """
    if info == "." or info == "":
        return {}
    # TODO: split on ';'. For each piece, split on '=' to get (key, value).
    # TODO: if no '=' (it is a flag), use '' as the value.
    raise NotImplementedError("Parse the INFO column")


def parse_sample(format_keys: list[str], sample_col: str) -> dict[str, str]:
    """Parse a per-sample column into a dict of KEY -> VALUE.

    Both the FORMAT column and the per-sample column are colon-separated.
    The two lists are parallel: format_keys[i] is the name of the i-th
    sample value.

    Examples:
        parse_sample(['GT', 'PL', 'DP', 'AD'], '1:255,0:42:0,42')
            -> {'GT': '1', 'PL': '255,0', 'DP': '42', 'AD': '0,42'}

    If the sample column has fewer values than format_keys (which is
    allowed when trailing fields are missing), pad with '.':
        parse_sample(['GT', 'PL', 'DP'], '0/1:50,0')
            -> {'GT': '0/1', 'PL': '50,0', 'DP': '.'}
    """
    # TODO: split sample_col on ':'.
    # TODO: zip with format_keys; pad with '.' if the sample is short.
    raise NotImplementedError("Parse the per-sample column")


def is_snp(ref: str, alts: list[str]) -> bool:
    """Return True if every allele (REF and all ALTs) is a single base.

    Examples:
        is_snp('A', ['G']) -> True
        is_snp('A', ['G', 'T']) -> True   (multiallelic SNP)
        is_snp('AG', ['A']) -> False      (deletion)
        is_snp('A', ['AG']) -> False      (insertion)
        is_snp('AG', ['CT']) -> False     (MNP, treated as not-SNP here)
    """
    # TODO: return True iff len(ref) == 1 and all(len(a) == 1 for a in alts).
    raise NotImplementedError("Decide if a variant is a SNP")


def is_indel(ref: str, alts: list[str]) -> bool:
    """Return True if the REF and any ALT differ in length.

    Examples:
        is_indel('A', ['G']) -> False     (SNP)
        is_indel('AG', ['A']) -> True     (deletion)
        is_indel('A', ['AG']) -> True     (insertion)
        is_indel('A', ['G', 'AG']) -> True  (mixed multiallelic; indel-ish)

    Note: per the VCF spec, an indel always has a 1-bp anchor base
    shared between REF and ALT. We do not check that here; we just
    check the length difference.
    """
    # TODO: return True iff any(len(a) != len(ref) for a in alts).
    raise NotImplementedError("Decide if a variant is an indel")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-02-vcf-parse.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # --- Test parse_info. ---------------------------------------------
    assert parse_info("DP=42;MQ=60") == {"DP": "42", "MQ": "60"}
    assert parse_info("DP=42;MQ=60;SP=2;AF=1") == {
        "DP": "42", "MQ": "60", "SP": "2", "AF": "1",
    }
    assert parse_info("INDEL;DP=30;MQ=60") == {
        "INDEL": "", "DP": "30", "MQ": "60",
    }
    assert parse_info(".") == {}
    assert parse_info("") == {}

    print("[exercise-02] parse_info: PASS")

    # --- Test parse_sample. -------------------------------------------
    assert parse_sample(
        ["GT", "PL", "DP", "AD"], "1:255,0:42:0,42"
    ) == {"GT": "1", "PL": "255,0", "DP": "42", "AD": "0,42"}
    assert parse_sample(
        ["GT", "PL", "DP"], "0/1:50,0"
    ) == {"GT": "0/1", "PL": "50,0", "DP": "."}
    assert parse_sample(
        ["GT"], "."
    ) == {"GT": "."}

    print("[exercise-02] parse_sample: PASS")

    # --- Test is_snp / is_indel. --------------------------------------
    assert is_snp("A", ["G"]) is True
    assert is_snp("A", ["G", "T"]) is True
    assert is_snp("AG", ["A"]) is False
    assert is_snp("A", ["AG"]) is False

    assert is_indel("A", ["G"]) is False
    assert is_indel("AG", ["A"]) is True
    assert is_indel("A", ["AG"]) is True
    assert is_indel("A", ["G", "AG"]) is True

    print("[exercise-02] is_snp / is_indel: PASS")

    # --- Test parse_header_line and parse_variant_line. ---------------
    header_records = []
    column_header = None
    variants = []
    for line in EXAMPLE_VCF.splitlines():
        if not line:
            continue
        if line.startswith("##"):
            header_records.append(parse_header_line(line))
        elif line.startswith("#"):
            # The column header line.
            parsed = parse_header_line(line)
            assert parsed[0] == "COLUMNS", (
                f"expected COLUMNS, got {parsed[0]}"
            )
            column_header = parsed[1]
        else:
            # A variant line.
            sample_names = column_header[9:] if column_header else []
            variants.append(parse_variant_line(line, sample_names))

    # We expect ~10 header records (fileformat, 2x FILTER, 1x contig,
    # 5x INFO, 4x FORMAT) and one COLUMNS line.
    assert column_header is not None
    assert column_header[0] == "CHROM"
    assert column_header[1] == "POS"
    assert column_header[8] == "FORMAT"
    assert column_header[9] == "lambda_sample"

    # We expect 7 variants.
    assert len(variants) == 7, f"expected 7 variants, got {len(variants)}"

    # First variant: SNP, PASS, QUAL 227.
    v = variants[0]
    assert v.chrom == "NC_001416.1"
    assert v.pos == 1001
    assert v.ref == "A"
    assert v.alts == ["G"]
    assert v.qual == 227.0
    assert v.filter == ["PASS"]
    assert v.info["DP"] == "42"
    assert v.info["MQ"] == "60"
    assert v.format == ["GT", "PL", "DP", "AD"]
    assert v.samples["lambda_sample"]["GT"] == "1"
    assert v.samples["lambda_sample"]["AD"] == "0,42"
    assert is_snp(v.ref, v.alts) is True
    assert is_indel(v.ref, v.alts) is False

    # Third variant: insertion, PASS.
    v = variants[2]
    assert v.pos == 5021
    assert v.ref == "A"
    assert v.alts == ["AGT"]
    assert "INDEL" in v.info
    assert is_indel(v.ref, v.alts) is True

    # Fourth variant: deletion, PASS.
    v = variants[3]
    assert v.pos == 12450
    assert v.ref == "CAA"
    assert v.alts == ["C"]
    assert "INDEL" in v.info
    assert is_indel(v.ref, v.alts) is True

    # Fifth variant: LowQual (depth too low).
    v = variants[4]
    assert v.pos == 20100
    assert v.filter == ["LowQual"]
    assert int(v.info["DP"]) < 10

    # Sixth variant: LowQual (high strand bias and low MQ).
    v = variants[5]
    assert v.pos == 30500
    assert v.filter == ["LowQual"]
    assert int(v.info["SP"]) > 60
    assert int(v.info["MQ"]) < 40

    print("[exercise-02] parse_header_line / parse_variant_line: PASS")

    # --- Higher-level summary. ----------------------------------------
    n_total = len(variants)
    n_pass = sum(1 for v in variants if v.filter == ["PASS"])
    n_low_qual = sum(1 for v in variants if v.filter == ["LowQual"])
    n_snp = sum(1 for v in variants if is_snp(v.ref, v.alts))
    n_indel = sum(1 for v in variants if is_indel(v.ref, v.alts))

    assert n_total == 7
    assert n_pass == 5
    assert n_low_qual == 2
    assert n_snp == 5
    assert n_indel == 2

    pass_quals = [v.qual for v in variants if v.filter == ["PASS"]]
    mean_pass_qual = sum(pass_quals) / len(pass_quals)
    assert 150 < mean_pass_qual < 250, (
        f"mean PASS QUAL {mean_pass_qual:.1f} outside expected range"
    )

    print()
    print("[exercise-02] Decoded 7 variant records and 12+ header lines.")
    print(f"[exercise-02] Reference contig: NC_001416.1 (48,502 bp).")
    print(f"[exercise-02] {n_total} variants total: {n_snp} SNPs, "
          f"{n_indel} indels.")
    print(f"[exercise-02] {n_pass} PASS, {n_low_qual} LowQual.")
    print(f"[exercise-02] Mean QUAL across PASS variants: "
          f"{mean_pass_qual:.1f}")
    print()
    print("[exercise-02] All assertions passed.")
    print("[exercise-02] You can now read a VCF file column by column")
    print("[exercise-02] and decode INFO and FORMAT/sample fields by hand.")
    print("[exercise-02] Continue to exercise-03 (VEP annotation).")
