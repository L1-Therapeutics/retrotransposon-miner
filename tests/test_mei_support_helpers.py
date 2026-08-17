"""Unit tests for pure-logic helpers in mei_support.py.

All functions tested here are purely computational (string/int arithmetic) and
require no BAM files, pysam calls, or external tools.
"""

from __future__ import annotations

import pytest

from retro_miner.mei_support import (
    _cigar_alignment_spans,
    _cigar_query_len,
    _clip_shannon_entropy,
    _homopolymer_at_run,
    _is_informative_split_clip,
    trim_mei_consensus_terminal_polya,
)


# ---------------------------------------------------------------------------
# _cigar_alignment_spans  →  (query_aligned_bp, ref_span_bp, alnlen)
#
#  M/=/X  consume both query and reference
#  I      consumes query only
#  D/N    consumes reference only
#  S/H    consume neither (not counted here)
# ---------------------------------------------------------------------------


def test_cigar_spans_pure_match() -> None:
    q, r, a = _cigar_alignment_spans("100M")
    assert q == 100 and r == 100 and a == 100


def test_cigar_spans_match_insertion_match() -> None:
    # 50M 10I 40M:  q=50+10+40=100  r=50+40=90  a=50+10+40=100
    q, r, a = _cigar_alignment_spans("50M10I40M")
    assert q == 100 and r == 90 and a == 100


def test_cigar_spans_match_deletion_match() -> None:
    # 50M 10D 40M:  q=50+40=90  r=50+10+40=100  a=50+10+40=100
    q, r, a = _cigar_alignment_spans("50M10D40M")
    assert q == 90 and r == 100 and a == 100


def test_cigar_spans_soft_clip_excluded() -> None:
    # Soft clips do not consume reference or contribute to alignment spans
    q, r, a = _cigar_alignment_spans("10S90M")
    assert q == 90 and r == 90 and a == 90


def test_cigar_spans_empty_string() -> None:
    assert _cigar_alignment_spans("") == (0, 0, 0)


def test_cigar_spans_none_like_empty() -> None:
    # None coerced to empty string internally
    assert _cigar_alignment_spans(None) == (0, 0, 0)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# _cigar_query_len  →  full query length (M + I + =/X + S + H)
# ---------------------------------------------------------------------------


def test_cigar_query_len_pure_match() -> None:
    assert _cigar_query_len("100M") == 100


def test_cigar_query_len_with_soft_clip() -> None:
    # 10S + 90M = 100 query bases
    assert _cigar_query_len("10S90M") == 100


def test_cigar_query_len_with_hard_clip() -> None:
    # 5H + 95M = 100
    assert _cigar_query_len("5H95M") == 100


def test_cigar_query_len_insertion_adds() -> None:
    # 50M + 10I + 40M = 100
    assert _cigar_query_len("50M10I40M") == 100


def test_cigar_query_len_deletion_excluded() -> None:
    # Deletions consume reference, not query
    assert _cigar_query_len("50M10D40M") == 90


def test_cigar_query_len_empty() -> None:
    assert _cigar_query_len("") == 0


# ---------------------------------------------------------------------------
# trim_mei_consensus_terminal_polya
# ---------------------------------------------------------------------------


def test_trim_trailing_polya_long_enough() -> None:
    # 4 trailing As ≥ default min_run → stripped
    result = trim_mei_consensus_terminal_polya("ACGTAAAA", min_run=4)
    assert result == "ACGT"


def test_trim_trailing_polya_too_short() -> None:
    # Only 3 trailing As < min_run=4 → not stripped
    result = trim_mei_consensus_terminal_polya("ACGTAAA", min_run=4)
    assert result == "ACGTAAA"


def test_trim_leading_polyt_long_enough() -> None:
    result = trim_mei_consensus_terminal_polya("TTTACGT", min_run=3)
    assert result == "ACGT"


def test_trim_leading_polyt_too_short() -> None:
    result = trim_mei_consensus_terminal_polya("TTACGT", min_run=3)
    assert result == "TTACGT"


def test_trim_both_ends() -> None:
    # Leading polyT and trailing polyA both stripped
    result = trim_mei_consensus_terminal_polya("TTTACGTAAAA", min_run=3)
    assert result == "ACGT"


def test_trim_empty_seq() -> None:
    assert trim_mei_consensus_terminal_polya("") == ""


def test_trim_converts_u_to_t() -> None:
    # U is treated as T; leading UUUU should be stripped like TTTT
    result = trim_mei_consensus_terminal_polya("UUUUACGT", min_run=4)
    assert result == "ACGT"


def test_trim_all_a_becomes_empty() -> None:
    result = trim_mei_consensus_terminal_polya("AAAA", min_run=4)
    assert result == ""


# ---------------------------------------------------------------------------
# _homopolymer_at_run  →  longest consecutive A-only or T-only run
# ---------------------------------------------------------------------------


def test_homopolymer_pure_a() -> None:
    assert _homopolymer_at_run("AAAAA") == 5


def test_homopolymer_pure_t() -> None:
    assert _homopolymer_at_run("TTTTTT") == 6


def test_homopolymer_a_longer() -> None:
    assert _homopolymer_at_run("TTTTAAAAAAAGGGGG") == 7


def test_homopolymer_t_longer() -> None:
    assert _homopolymer_at_run("TTTTTTAAAAGGGGG") == 6


def test_homopolymer_no_at_content() -> None:
    assert _homopolymer_at_run("GCGCGCGC") == 0


def test_homopolymer_alternating_at() -> None:
    # Each A and T run length is 1
    assert _homopolymer_at_run("ATATATAT") == 1


def test_homopolymer_empty() -> None:
    assert _homopolymer_at_run("") == 0


# ---------------------------------------------------------------------------
# _clip_shannon_entropy
# ---------------------------------------------------------------------------


def test_entropy_uniform() -> None:
    # Equal counts of A, C, G, T → max entropy = 2.0 bits
    assert _clip_shannon_entropy("ACGT") == pytest.approx(2.0)


def test_entropy_single_base_is_zero() -> None:
    assert _clip_shannon_entropy("AAAA") == pytest.approx(0.0)


def test_entropy_empty_is_zero() -> None:
    assert _clip_shannon_entropy("") == pytest.approx(0.0)


def test_entropy_two_equal_bases() -> None:
    # 50% A, 50% C → H = 1.0 bit
    assert _clip_shannon_entropy("ACACAC") == pytest.approx(1.0)


def test_entropy_increases_with_diversity() -> None:
    polya = _clip_shannon_entropy("AAAAAAAAAA")
    mixed = _clip_shannon_entropy("ACGTACGTAC")
    assert mixed > polya


# ---------------------------------------------------------------------------
# _is_informative_split_clip
# ---------------------------------------------------------------------------


def test_informative_clip_too_short() -> None:
    # Default min_len=20; 10 chars is too short
    assert _is_informative_split_clip("ACGTACGTAC") is False


def test_informative_clip_polya_fails_non_at_fraction() -> None:
    # All A → 0% non-A/T (C or G) → fails fraction threshold
    assert _is_informative_split_clip("A" * 30) is False


def test_informative_clip_diverse_seq_passes() -> None:
    # 40 chars, equal ACGT → high entropy, 50% CG
    assert _is_informative_split_clip("ACGT" * 10) is True


def test_informative_clip_borderline_non_at_fails() -> None:
    # 22 A + 3 CG in 25 chars: non_at_fraction = 3/25 = 0.12 < 0.15
    seq = "A" * 22 + "CGT"
    assert _is_informative_split_clip(seq) is False
