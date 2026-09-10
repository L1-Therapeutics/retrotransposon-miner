"""Unit tests for pure-logic helpers in evidence_extract.py.

None of these tests require a real BAM file or pysam; they exercise only the
string/sequence functions that sit below the BAM-scanning layer.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from retro_miner.evidence_extract import (
    ExtractionSummary,
    _DISCORDANT_DEDUP_KEYS,
    _SPLIT_DEDUP_KEYS,
    _clip_to_poly_at_region,
    _drop_duplicate_evidence_rows,
    _longest_poly_at_span,
    _normalize_regions,
    _poly_at_breakpoint_proximal_stats,
    _poly_at_stats,
    _soft_clip_query_seq,
    _write_discordant_table,
    _write_split_table,
    clone_extraction_summary,
    clone_sample_evidence_tables,
    same_alignment_path,
)


# ---------------------------------------------------------------------------
# _normalize_regions
# ---------------------------------------------------------------------------


def test_normalize_regions_single_string() -> None:
    assert _normalize_regions("chr22") == ["chr22"]


def test_normalize_regions_list_passthrough() -> None:
    assert _normalize_regions(["chr1", "chr22"]) == ["chr1", "chr22"]


def test_normalize_regions_strips_whitespace() -> None:
    assert _normalize_regions(["  chr22  ", " chr1"]) == ["chr22", "chr1"]


def test_normalize_regions_drops_blank_entries() -> None:
    assert _normalize_regions(["chr22", "", "chr1"]) == ["chr22", "chr1"]


def test_normalize_regions_empty_list_raises() -> None:
    with pytest.raises(ValueError, match="No valid regions"):
        _normalize_regions([])


def test_normalize_regions_all_blank_raises() -> None:
    with pytest.raises(ValueError, match="No valid regions"):
        _normalize_regions(["  ", ""])


# ---------------------------------------------------------------------------
# _soft_clip_query_seq
# ---------------------------------------------------------------------------


def test_soft_clip_left() -> None:
    assert _soft_clip_query_seq("ABCDE", "L", 3) == "ABC"


def test_soft_clip_right() -> None:
    assert _soft_clip_query_seq("ABCDE", "R", 3) == "CDE"


def test_soft_clip_zero_len_returns_empty() -> None:
    assert _soft_clip_query_seq("ABCDE", "L", 0) == ""


def test_soft_clip_longer_than_seq_returns_empty() -> None:
    # clip_len > len(seq) → empty
    assert _soft_clip_query_seq("ABC", "L", 10) == ""


def test_soft_clip_unknown_side_returns_empty() -> None:
    assert _soft_clip_query_seq("ABCDE", "X", 3) == ""


def test_soft_clip_empty_seq_returns_empty() -> None:
    assert _soft_clip_query_seq("", "L", 3) == ""


# ---------------------------------------------------------------------------
# _poly_at_stats
# ---------------------------------------------------------------------------


def test_poly_at_stats_empty_seq() -> None:
    run, frac, base = _poly_at_stats("")
    assert run == 0 and frac == 0.0 and base == ""


def test_poly_at_stats_no_at_content() -> None:
    run, frac, base = _poly_at_stats("GCGCGCGC")
    assert run == 0 and frac == 0.0 and base == ""


def test_poly_at_stats_pure_polya() -> None:
    run, frac, base = _poly_at_stats("AAAA")
    assert run == 4 and frac == pytest.approx(1.0) and base == "A"


def test_poly_at_stats_pure_polyt() -> None:
    run, frac, base = _poly_at_stats("TTTTTT")
    assert run == 6 and frac == pytest.approx(1.0) and base == "T"


def test_poly_at_stats_mostly_a() -> None:
    # AAAGAAA → 6A/7 total; best run = 3
    run, frac, base = _poly_at_stats("AAAGAAA")
    assert base == "A"
    assert run == 3
    assert frac == pytest.approx(6 / 7)


def test_poly_at_stats_mixed_at_a_wins_tie() -> None:
    # ATAT → n_a == n_t → A wins by '>=' rule; best run of A = 1
    run, frac, base = _poly_at_stats("ATAT")
    assert base == "A"
    assert run == 1
    assert frac == pytest.approx(0.5)


def test_poly_at_stats_t_dominates() -> None:
    # TTTAA → n_t=3 > n_a=2 → T; best T run = 3
    run, frac, base = _poly_at_stats("TTTAA")
    assert base == "T"
    assert run == 3
    assert frac == pytest.approx(3 / 5)


# ---------------------------------------------------------------------------
# _longest_poly_at_span
# ---------------------------------------------------------------------------


_POLY_A_30 = "A" * 30
_POLY_T_30 = "T" * 30


def test_longest_poly_at_span_too_short() -> None:
    # Sequence shorter than default min_len=25 → all zeros
    length, frac, base, span = _longest_poly_at_span("A" * 20)
    assert length == 0 and frac == 0.0 and base == "" and span == ""


def test_longest_poly_at_span_pure_polya() -> None:
    length, frac, base, span = _longest_poly_at_span(_POLY_A_30)
    assert length == 30
    assert frac == pytest.approx(1.0)
    assert base == "A"
    assert span == _POLY_A_30


def test_longest_poly_at_span_pure_polyt() -> None:
    length, frac, base, span = _longest_poly_at_span(_POLY_T_30)
    assert base == "T"
    assert length == 30


def test_longest_poly_at_span_below_purity_threshold() -> None:
    # 50% A 50% G — no 90%-pure window of >= 25 bases
    seq = "AG" * 40  # 80 chars, 50% A
    length, frac, base, span = _longest_poly_at_span(seq)
    assert length == 0


def test_longest_poly_at_span_high_purity_polya() -> None:
    # ~91% A over 33 bases (3 non-A interspersed)
    seq = ("A" * 10 + "C" + "A" * 10 + "G" + "A" * 10 + "C")
    # min_frac default is 0.90; a window of 31 chars has 30A/31 = 96.8% A
    length, frac, base, span = _longest_poly_at_span(seq)
    assert base == "A"
    assert length >= 25
    assert frac >= 0.90


# ---------------------------------------------------------------------------
# _clip_to_poly_at_region  (thin wrapper)
# ---------------------------------------------------------------------------


def test_clip_to_poly_at_region_returns_span() -> None:
    result = _clip_to_poly_at_region("A" * 30)
    assert result == "A" * 30


def test_clip_to_poly_at_region_returns_empty_when_below_threshold() -> None:
    result = _clip_to_poly_at_region("GCGC" * 20)  # no A/T content
    assert result == ""


# ---------------------------------------------------------------------------
# _poly_at_breakpoint_proximal_stats
# ---------------------------------------------------------------------------


def test_poly_at_bp_proximal_left_dominates() -> None:
    # Left window = poly-A; right window = all C (no A/T)
    seq = "A" * 10 + "C" * 10
    run, frac, base, side = _poly_at_breakpoint_proximal_stats(seq, window_bases=10)
    assert side == "L" and base == "A" and run == 10


def test_poly_at_bp_proximal_right_dominates() -> None:
    seq = "C" * 10 + "T" * 10
    run, frac, base, side = _poly_at_breakpoint_proximal_stats(seq, window_bases=10)
    assert side == "R" and base == "T" and run == 10


def test_poly_at_bp_proximal_empty_seq() -> None:
    run, frac, base, side = _poly_at_breakpoint_proximal_stats("", window_bases=10)
    assert run == 0 and frac == 0.0


def test_poly_at_bp_proximal_no_at_content() -> None:
    run, frac, base, side = _poly_at_breakpoint_proximal_stats("GCGCGCGCGC", window_bases=5)
    assert run == 0 and frac == 0.0


def test_poly_at_bp_proximal_both_windows_equal_left_wins() -> None:
    # Both ends have the same poly-A stats; left wins via '>=' in comparison
    seq = "A" * 5 + "GCGC" + "A" * 5
    run, frac, base, side = _poly_at_breakpoint_proximal_stats(seq, window_bases=5)
    assert side == "L"  # tie broken by left-first '>=' comparison


def test_drop_duplicate_discordant_rows_keeps_first() -> None:
    row = {
        "read_name": "q1",
        "chrom": "chr22",
        "pos": 100,
        "is_read1": False,
        "mate_chrom": "chr1",
        "mate_pos": 5000,
        "is_reverse": True,
        "mate_is_reverse": True,
        "is_proper_pair": False,
        "mapq": 60,
    }
    df = _drop_duplicate_evidence_rows(pd.DataFrame([row, dict(row), dict(row, mapq=1)]), _DISCORDANT_DEDUP_KEYS)
    assert len(df) == 1
    assert int(df.iloc[0]["mapq"]) == 60


def test_write_discordant_table_drops_exact_duplicate_rows(tmp_path: Path) -> None:
    row = {
        "sample": "s",
        "chrom": "chr22",
        "pos": 100,
        "is_read1": False,
        "mate_chrom": "chr1",
        "mate_pos": 5000,
        "is_reverse": True,
        "mate_is_reverse": True,
        "is_proper_pair": False,
        "read_name": "q1",
        "discordant_reasons": "interchrom,improper_pair",
    }
    df = _write_discordant_table([row, row], tmp_path, "s")
    assert len(df) == 1
    written = pd.read_csv(tmp_path / "discordant_evidence.s.tsv", sep="\t")
    assert len(written) == 1


def test_write_split_table_drops_exact_duplicate_rows(tmp_path: Path) -> None:
    row = {
        "sample": "s",
        "chrom": "chr22",
        "pos": 100,
        "clip_side": "L",
        "mate_chrom": "chr1",
        "mate_pos": 5000,
        "is_reverse": False,
        "read_name": "q1",
        "clip_len": 25,
    }
    other = dict(row, clip_side="R")
    df = _write_split_table([row, row, other], tmp_path, "s")
    assert len(df) == 2
    assert set(df["clip_side"]) == {"L", "R"}


def test_drop_duplicate_split_rows_keeps_distinct_clip_sides() -> None:
    left = {
        "read_name": "q1",
        "chrom": "chr22",
        "pos": 100,
        "clip_side": "L",
        "mate_chrom": "chr1",
        "mate_pos": 5000,
        "is_reverse": False,
    }
    df = _drop_duplicate_evidence_rows(pd.DataFrame([left, left, dict(left, clip_side="R")]), _SPLIT_DEDUP_KEYS)
    assert len(df) == 2


def test_same_alignment_path_resolves_duplicates(tmp_path: Path) -> None:
    bam = tmp_path / "sample.cram"
    bam.write_bytes(b"cram")
    assert same_alignment_path(bam, tmp_path / "sample.cram")
    assert not same_alignment_path(bam, tmp_path / "other.cram")
    assert not same_alignment_path(bam, None)


def test_clone_sample_evidence_tables_rewrites_sample(tmp_path: Path) -> None:
    src = pd.DataFrame({"read_name": ["r1"], "sample": ["control"], "chrom": ["chr21"]})
    src.to_csv(tmp_path / "split_evidence.control.tsv", sep="\t", index=False)
    src.to_parquet(tmp_path / "split_evidence.control.parquet", index=False)
    written = clone_sample_evidence_tables(tmp_path, src_sample="control", dst_sample="disease")
    assert any(path.endswith("split_evidence.disease.tsv") for path in written)
    cloned = pd.read_csv(tmp_path / "split_evidence.disease.tsv", sep="\t")
    assert list(cloned["sample"]) == ["disease"]
    assert list(cloned["read_name"]) == ["r1"]


def test_clone_extraction_summary_renames_sample() -> None:
    src = ExtractionSummary(sample="control", total_reads_scanned=10, passing_reads=8, split_evidence_rows=2)
    dst = clone_extraction_summary(src, "disease")
    assert dst.sample == "disease"
    assert dst.total_reads_scanned == 10
