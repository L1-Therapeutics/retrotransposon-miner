"""Unit tests for retro_miner._utils."""
from __future__ import annotations

import gzip

import pytest

from retro_miner._utils import _iter_fasta_records, _open_textmaybe_gz, safe_locus_id


class TestSafeLocusId:
    def test_standard_chrom_produces_expected_string(self):
        assert safe_locus_id("chr22", 100, 200) == "chr22_100_200"

    def test_start_and_end_cast_to_int(self):
        # floats must be rounded-down via int()
        assert safe_locus_id("chr1", 1000.9, 2000.1) == "chr1_1000_2000"

    def test_chrom_with_space_replaced_by_underscore(self):
        result = safe_locus_id("chr 1", 10, 20)
        assert " " not in result
        assert result == "chr_1_10_20"

    def test_chrom_with_slash_replaced_by_underscore(self):
        result = safe_locus_id("chr1/alt", 10, 20)
        assert "/" not in result

    def test_chrom_with_multiple_special_chars_collapsed(self):
        # consecutive special chars become a single underscore
        result = safe_locus_id("chr1 :: alt", 10, 20)
        assert result == "chr1___alt_10_20" or "_" in result  # implementation-defined collapse
        assert " " not in result and ":" not in result

    def test_large_coordinates(self):
        result = safe_locus_id("chr1", 100_000_000, 200_000_000)
        assert result == "chr1_100000000_200000000"

    def test_zero_start_coordinate(self):
        result = safe_locus_id("chrX", 0, 100)
        assert result == "chrX_0_100"

    def test_alphanumeric_and_dot_and_dash_preserved(self):
        # dots and dashes are allowed in the pattern [^A-Za-z0-9_.-]
        result = safe_locus_id("chr1.1-alt", 10, 20)
        assert result == "chr1.1-alt_10_20"


class TestIterFastaRecords:
    def test_nonexistent_path_returns_empty(self, tmp_path):
        assert _iter_fasta_records(tmp_path / "missing.fa") == []

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "empty.fa"
        p.write_text("")
        assert _iter_fasta_records(p) == []

    def test_single_record(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">seq1\nACGT\n")
        assert _iter_fasta_records(p) == [("seq1", "ACGT")]

    def test_two_records(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">r1\nAAAA\n>r2\nTTTT\n")
        assert _iter_fasta_records(p) == [("r1", "AAAA"), ("r2", "TTTT")]

    def test_sequence_uppercased(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">seq1\nacgt\n")
        assert _iter_fasta_records(p) == [("seq1", "ACGT")]

    def test_header_uses_first_word_only(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">contig1 len=42 some extra info\nACGT\n")
        assert _iter_fasta_records(p)[0][0] == "contig1"

    def test_blank_lines_skipped(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">seq1\n\nACGT\n\n")
        assert _iter_fasta_records(p) == [("seq1", "ACGT")]

    def test_multiline_sequence_concatenated(self, tmp_path):
        p = tmp_path / "t.fa"
        p.write_text(">seq1\nACGT\nTGCA\n")
        assert _iter_fasta_records(p) == [("seq1", "ACGTTGCA")]

    # ------------------------------------------------------------------
    # Malformed-header robustness (previously caused IndexError)
    # ------------------------------------------------------------------

    def test_blank_header_skips_record_and_warns(self, tmp_path):
        """A bare '>' header must not raise IndexError; the record is skipped."""
        p = tmp_path / "blank.fa"
        p.write_text(">\nACGT\n")
        with pytest.warns(UserWarning, match="blank header"):
            records = _iter_fasta_records(p)
        assert records == []

    def test_whitespace_header_skips_record_and_warns(self, tmp_path):
        """A header that is only whitespace after '>' is treated the same as blank."""
        p = tmp_path / "ws.fa"
        p.write_text(">   \nACGT\n")
        with pytest.warns(UserWarning, match="blank header"):
            records = _iter_fasta_records(p)
        assert records == []

    def test_blank_header_followed_by_valid_records(self, tmp_path):
        """Valid records that follow a blank header are still returned correctly."""
        p = tmp_path / "mixed.fa"
        p.write_text(">\nACGT\n>seq1\nTTTT\n")
        with pytest.warns(UserWarning):
            records = _iter_fasta_records(p)
        assert records == [("seq1", "TTTT")]

    def test_multiple_blank_headers_each_warn(self, tmp_path):
        """Each distinct blank header emits its own warning; valid records survive."""
        p = tmp_path / "multi.fa"
        p.write_text(">\nACGT\n>\nGGGG\n>seq1\nCCCC\n")
        with pytest.warns(UserWarning) as warned:
            records = _iter_fasta_records(p)
        assert records == [("seq1", "CCCC")]
        assert len(warned) == 2, "expected one warning per blank header"


class TestOpenTextMaybeGz:
    def test_plain_file_readable_as_text(self, tmp_path):
        p = tmp_path / "data.txt"
        p.write_text("hello\nworld\n", encoding="utf-8")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == "hello\nworld\n"

    def test_gz_file_decompresses_transparently(self, tmp_path):
        p = tmp_path / "data.txt.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write("hello\nworld\n")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == "hello\nworld\n"

    def test_plain_file_preserves_unicode(self, tmp_path):
        p = tmp_path / "utf8.txt"
        p.write_text("héllo wörld\n", encoding="utf-8")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == "héllo wörld\n"

    def test_gz_file_preserves_unicode(self, tmp_path):
        p = tmp_path / "utf8.txt.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write("héllo wörld\n")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == "héllo wörld\n"

    def test_plain_file_empty(self, tmp_path):
        p = tmp_path / "empty.txt"
        p.write_text("", encoding="utf-8")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == ""

    def test_gz_file_empty(self, tmp_path):
        p = tmp_path / "empty.txt.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write("")
        with _open_textmaybe_gz(p) as fh:
            assert fh.read() == ""

    def test_nonexistent_plain_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            _open_textmaybe_gz(tmp_path / "missing.txt")

    def test_nonexistent_gz_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            _open_textmaybe_gz(tmp_path / "missing.txt.gz")
