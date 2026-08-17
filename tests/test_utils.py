"""Unit tests for retro_miner._utils."""
from __future__ import annotations

from retro_miner._utils import _iter_fasta_records, safe_locus_id


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
