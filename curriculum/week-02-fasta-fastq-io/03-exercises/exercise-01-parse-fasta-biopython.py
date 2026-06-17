"""
Exercise 1 - Parse FASTA with Biopython.

Goal: re-parse the same kind of FASTA you parsed by hand in Week 1, this
time with Bio.SeqIO. Use the official API. Use SeqRecord. Use SeqIO.write
to round-trip the file. Notice how much shorter the code is, and how
much more it knows about the underlying biology than your Week-1
parser did.

Estimated time: 35 minutes.

Acceptance criteria:
- `python exercise-01-parse-fasta-biopython.py` runs without crashing.
- All `assert` checks at the bottom pass.
- You wrote a TEMPORARY FASTA file, parsed it back with SeqIO, computed
  per-record stats, wrote the records back out wrapped at 60 columns,
  and confirmed the round-trip preserved the data.
- You used `Bio.SeqIO.parse`, `Bio.SeqIO.read`, `Bio.SeqIO.to_dict`, and
  `Bio.SeqIO.write` at least once each.

Required: biopython 1.83 or compatible.
    python -m pip install biopython==1.83

What you learn:
- The four-function SeqIO API and when to use which.
- The shape of a SeqRecord (id, description, seq, annotations).
- How to round-trip a FASTA without losing anything.
- Why `SeqIO.read` (singular) is a fail-loud helper for one-record files.

TO COMPLETE: implement the three functions below. Run the file; all
assertions must pass.
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# A small inline FASTA we will write to disk and then round-trip.
# Three records, different sizes, mixed case, one record with an N.
TEST_FASTA = """>seq1 first record
ATGCGTACGT
>seq2 second record, multi-line
atgcGCGCAAAA
AACCGGTT
>seq3 third record with N's
NNNNACGTNNNN
"""

# Just the first record on its own - we use this with SeqIO.read.
SINGLE_RECORD_FASTA = """>spike_fragment NC_045512.2 partial SARS-CoV-2 spike CDS
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAAT
"""


def parse_with_biopython(fasta_path: Path) -> list[SeqRecord]:
    """Parse a FASTA file with Bio.SeqIO.parse and return a list of SeqRecords.

    Use Bio.SeqIO.parse - the streaming generator API. Materialize the
    iterator into a list before returning (the caller wants random access).

    Returns:
        list of SeqRecord, in file order.
    """
    # TODO: implement using SeqIO.parse.
    # Hint: one line of real work plus a `return`.
    raise NotImplementedError("Implement parse_with_biopython")


def fasta_to_dict(fasta_path: Path) -> dict[str, SeqRecord]:
    """Load a FASTA file into a dict keyed by record.id.

    Use Bio.SeqIO.to_dict (in-memory index, fine for small files).
    """
    # TODO: implement using SeqIO.to_dict(SeqIO.parse(...)).
    raise NotImplementedError("Implement fasta_to_dict")


def write_records(records: list[SeqRecord], out_path: Path) -> int:
    """Write records to disk as FASTA at the given path. Return record count.

    Use Bio.SeqIO.write. The default line-wrap is 60 columns - leave it
    at the default; the round-trip assertions below depend on it.
    """
    # TODO: implement using SeqIO.write.
    raise NotImplementedError("Implement write_records")


# ----------------------------------------------------------------------
# Self-test.
# Run with:  python exercise-01-parse-fasta-biopython.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        input_fasta = tmp_dir / "input.fasta"
        output_fasta = tmp_dir / "output.fasta"
        single_fasta = tmp_dir / "single.fasta"

        # Write test inputs to disk so we are exercising the real file path.
        input_fasta.write_text(TEST_FASTA)
        single_fasta.write_text(SINGLE_RECORD_FASTA)

        # 1) parse_with_biopython.
        records = parse_with_biopython(input_fasta)
        assert len(records) == 3, f"expected 3 records, got {len(records)}"
        assert records[0].id == "seq1"
        assert records[1].id == "seq2"
        assert records[2].id == "seq3"
        # Biopython preserves case in FASTA sequences.
        assert str(records[1].seq) == "atgcGCGCAAAAAACCGGTT", (
            f"seq2 multi-line join wrong: {records[1].seq!r}"
        )

        # 2) fasta_to_dict.
        idx = fasta_to_dict(input_fasta)
        assert set(idx.keys()) == {"seq1", "seq2", "seq3"}
        assert str(idx["seq3"].seq) == "NNNNACGTNNNN"

        # 3) write_records round-trip. Default line-wrap is 60 cols, so all
        # three of our test sequences (each < 60 nt) come out on one line.
        n_written = write_records(records, output_fasta)
        assert n_written == 3, f"write_records should return 3, got {n_written}"

        round_tripped = parse_with_biopython(output_fasta)
        assert len(round_tripped) == 3
        for original, after in zip(records, round_tripped):
            assert original.id == after.id, (original.id, after.id)
            assert str(original.seq) == str(after.seq), (
                f"round-trip mismatch on {original.id}: "
                f"{original.seq!r} vs {after.seq!r}"
            )

        # 4) SeqIO.read on a single-record file - must succeed.
        spike = SeqIO.read(single_fasta, "fasta")
        assert spike.id == "spike_fragment"
        assert str(spike.seq).startswith("ATGTTTGTTTTT")

        # 5) SeqIO.read on a MULTI-record file - must raise ValueError.
        try:
            SeqIO.read(input_fasta, "fasta")
        except ValueError:
            pass
        else:
            raise AssertionError(
                "SeqIO.read on a multi-record FASTA must raise ValueError"
            )

        # 6) Show off SeqRecord features: build one from scratch and round-trip.
        built = SeqRecord(
            Seq("ATGCATGCATGCATGCATGC"),
            id="synthetic_01",
            description="synthetic_01 a record I built in Python",
        )
        buf = io.StringIO()
        SeqIO.write([built], buf, "fasta")
        emitted = buf.getvalue()
        assert emitted.startswith(">synthetic_01 "), emitted[:40]
        assert "ATGCATGCATGCATGCATGC" in emitted

    print("All assertions passed. You have used the SeqIO four-function API.")
    print("Move on to Exercise 2 - FASTQ quality plotting.")
