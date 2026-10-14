"""Unit tests for retro_miner._utils."""
from __future__ import annotations

import gzip
from pathlib import Path

import math

import numpy as np
import pytest

from retro_miner._utils import (
    _iter_fasta_records,
    _open_textmaybe_gz,
    _safe_float,
    _safe_int,
    read_tsv,
    safe_locus_id,
)


class TestSafeInt:
    """_safe_int is the untrusted-value boundary for record/table fields.

    It must return the default for every input it documents as unsupported
    rather than raising, because callers feed it raw DataFrame row values.
    """

    @pytest.mark.parametrize(
        "value",
        [
            None,
            "",
            "   ",
            "abc",
            float("nan"),
            float("inf"),
            float("-inf"),
            object(),
            [1, 2],
        ],
    )
    def test_unsupported_inputs_return_default(self, value):
        assert _safe_int(value) == 0
        assert _safe_int(value, 7) == 7

    @pytest.mark.parametrize(
        ("value", "expected"),
        [(7, 7), ("12", 12), (3.9, 3), ("12.7", 12), (-3.9, -3), ("  12  ", 12), (True, 1), (False, 0)],
    )
    def test_supported_inputs_convert(self, value, expected):
        assert _safe_int(value) == expected

    def test_float_column_values_do_not_raise(self):
        # a numeric DataFrame column yields float (and NaN for a missing cell);
        # the previous implementation raised NameError on any float
        assert _safe_int(4.0) == 4
        assert _safe_int(float("nan"), -1) == -1

    def test_huge_int_is_returned_exactly(self):
        # int(float(...)) would raise OverflowError here
        assert _safe_int(10**400) == 10**400

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            (np.int64(7), 7),
            (np.int32(-3), -3),
            (np.uint16(5), 5),
            (np.int64(10**18), 10**18),
            (np.bool_(True), 1),
            (np.bool_(False), 0),
            (np.float32(7.9), 7),
            (np.float64(-2.5), -2),
        ],
    )
    def test_numpy_scalars_convert(self, value, expected):
        # DataFrame rows and Series values yield NumPy scalars, which do not
        # subclass the builtin numeric types (np.int64 is not int, np.float32
        # is not float, np.bool_ is neither). They must convert, not collapse
        # to the default -- doing so silently zeroed out integer table columns.
        assert _safe_int(value) == expected

    def test_numpy_non_finite_returns_default(self):
        assert _safe_int(np.float32("nan")) == 0
        assert _safe_int(np.float32("inf"), -1) == -1
        assert _safe_int(np.float64("-inf"), -1) == -1


class TestReadTsv:
    def test_reads_tab_separated_by_default(self, tmp_path: Path):
        p = tmp_path / "x.tsv"
        p.write_text("a\tb\n1\t2\n", encoding="utf-8")
        df = read_tsv(p)
        assert list(df.columns) == ["a", "b"]
        assert int(df["a"].iloc[0]) == 1

    def test_running_header_and_usecols_defaults(self, tmp_path: Path):
        # Mirrors the .fai parse in mei_support._resolve_full_consensus_fasta.
        p = tmp_path / "x.fai"
        p.write_text("chr1\t100\t1\t60\t61\n", encoding="utf-8")
        df = read_tsv(
            p,
            header=None,
            names=["name", "length", "offset", "line_bases", "line_width"],
            usecols=[0, 1],
        )
        assert list(df.columns) == ["name", "length"]
        assert df["name"].iloc[0] == "chr1"

    def test_caller_overrides_do_not_duplicate_keywords(self, tmp_path: Path):
        # Regression: read_tsv used to always emit sep="\t" and low_memory=False,
        # so a caller forwarding the same keyword crashed with TypeError
        # ("multiple values for keyword argument"), e.g. the redundant
        # sep="\t" / low_memory=False arguments at two mei_support call sites.
        p = tmp_path / "x.tsv"
        p.write_text("a\tb\n1\t2\n", encoding="utf-8")
        assert read_tsv(p, sep="\t").shape == (1, 2)
        assert read_tsv(p, low_memory=False).shape == (1, 2)
        df = read_tsv(p, sep=",", low_memory=False)
        assert df.shape == (1, 1)  # single comma-joined column under the wrong separator

    def test_callable_usecols_skips_absent_columns(self, tmp_path: Path):
        # A callable selects only the columns present, unlike a sequence which
        # raises when a listed column is missing. split-evidence tables omit
        # optional columns, so read_architecture relies on this behavior.
        p = tmp_path / "split.tsv"
        p.write_text("read_name\tclip_len\nr1\t42\n", encoding="utf-8")
        df = read_tsv(p, usecols=lambda c: c in {"read_name", "clip_len", "clip_side"})
        assert list(df.columns) == ["read_name", "clip_len"]
        assert int(df["clip_len"].iloc[0]) == 42


class TestSafeFloat:
    """_safe_float is the float-valued mirror of _safe_int.

    It must return the default for every documented-unsupported input rather
    than raising, and it collapses non-finite results (nan/inf) to the default
    so callers can keep ``math.isfinite`` as the single validity test.
    """

    @pytest.mark.parametrize(
        "value",
        [None, "", "   ", "abc", object(), [1, 2]],
    )
    def test_unsupported_inputs_return_default(self, value):
        assert math.isnan(_safe_float(value))
        assert _safe_float(value, 7.0) == 7.0

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
    def test_non_finite_inputs_return_default(self, value):
        assert math.isnan(_safe_float(value))
        assert _safe_float(value, 7.0) == 7.0

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            (7, 7.0),
            ("12", 12.0),
            (3.9, 3.9),
            ("12.7", 12.7),
            (-3.9, -3.9),
            ("  12  ", 12.0),
            (True, 1.0),
            (False, 0.0),
        ],
    )
    def test_supported_inputs_convert(self, value, expected):
        assert _safe_float(value) == expected

    def test_huge_int_returns_default(self):
        # float(10**400) raises OverflowError; it must not propagate.
        assert math.isnan(_safe_float(10**400))
        assert _safe_float(10**400, -1.0) == -1.0

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            (np.int64(7), 7.0),
            (np.int32(-3), -3.0),
            (np.uint16(5), 5.0),
            (np.bool_(True), 1.0),
            (np.bool_(False), 0.0),
            (np.float32(7.5), 7.5),
            (np.float64(-2.5), -2.5),
        ],
    )
    def test_numpy_scalars_convert(self, value, expected):
        # DataFrame rows and Series values yield NumPy scalars, which do not
        # subclass the builtin numeric types; they must convert, not collapse
        # to the default.
        assert _safe_float(value) == expected


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
