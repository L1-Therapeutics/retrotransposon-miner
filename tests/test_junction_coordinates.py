"""Split-read junction positions must name real aligned bases.

pysam's ``reference_start``/``reference_end`` are 0-based, and ``reference_end``
is *exclusive*. The 1-based inclusive coordinates of the aligned block are
therefore ``reference_start + 1`` (first base) and ``reference_end`` (last base).

A previous fix added ``+ 1`` to the end, which pushed every right-clip junction
one base past the last aligned base. Nothing caught it: the alignment-end
position is only consumed by the split-row writer, which had no coverage.
"""

from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from retro_miner.evidence_extract import _longest_soft_clip_from_read  # noqa: E402


def _read(cigar: str, ref_start: int, query_seq: str) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment()
    seg.query_name = "synthetic"
    seg.reference_id = 0
    seg.reference_start = ref_start
    seg.cigarstring = cigar
    seg.query_sequence = query_seq
    seg.query_qualities = pysam.qualitystring_to_array("I" * len(query_seq))
    return seg


def test_reference_end_is_zero_based_exclusive() -> None:
    """Guard the pysam convention the rest of the file depends on."""
    seg = _read("90M", 1000, "A" * 90)
    assert seg.reference_start == 1000
    assert seg.reference_end == 1090
    # 1-based inclusive span of the aligned block is 1001..1090.
    assert seg.reference_end - seg.reference_start == 90


def test_right_clip_junction_is_last_aligned_base() -> None:
    seg = _read("90M30S", 1000, "A" * 90 + "C" * 30)
    side, clip_len, pos, _clip_seq = _longest_soft_clip_from_read(seg)
    assert side == "R"
    assert clip_len == 30
    # Last aligned base, 1-based inclusive. 1091 would be one past the read.
    assert pos == 1090


def test_left_clip_junction_is_first_aligned_base() -> None:
    seg = _read("30S90M", 2000, "C" * 30 + "A" * 90)
    side, clip_len, pos, _clip_seq = _longest_soft_clip_from_read(seg)
    assert side == "L"
    assert clip_len == 30
    assert pos == 2001


def test_both_junction_ends_stay_inside_the_aligned_block() -> None:
    """A split read's two junction ends are the first and last aligned bases."""
    seq = "A" * 60 + "C" * 25
    right = _read("60M25S", 5000, seq)
    left = _read("25S60M", 5000, seq)
    first_aligned = left.reference_start + 1
    last_aligned = right.reference_end
    assert last_aligned - first_aligned + 1 == 60

    r_side, _r_len, r_pos, _ = _longest_soft_clip_from_read(right)
    l_side, _l_len, l_pos, _ = _longest_soft_clip_from_read(left)
    assert (r_side, r_pos) == ("R", last_aligned)
    assert (l_side, l_pos) == ("L", first_aligned)
