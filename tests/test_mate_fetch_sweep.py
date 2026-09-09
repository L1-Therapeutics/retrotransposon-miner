"""Tile-coalesced mate fetch must match per-window bam.fetch results."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from retro_miner.mei_support import (
    _MATE_FETCH_QUERY_BP,
    _MATE_FETCH_TILE_BP,
    _fetch_mates_per_window,
    _fetch_mates_swept,
    _merge_mate_fetch_windows,
)


def _have_pysam() -> bool:
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return shutil.which("samtools") is not None or True


class _FetchRecorder:
    """Wrap a BAM so tests can count ``fetch`` spans without changing results."""

    def __init__(self, bam):
        self._bam = bam
        self.fetches: list[tuple[str, int, int]] = []

    def fetch(self, contig, start=None, end=None, **kwargs):
        self.fetches.append((str(contig), int(start), int(end)))
        return self._bam.fetch(contig, start, end, **kwargs)

    def __getattr__(self, name):
        return getattr(self._bam, name)


def test_merge_mate_windows_same_tile_one_tight_union() -> None:
    windows = [
        ("chr22", 100, 600, "a"),
        ("chr22", 8_000, 8_500, "b"),
        ("chr1", 10, 510, "c"),
    ]
    merged = _merge_mate_fetch_windows(windows)
    assert len(merged) == 2
    chr22 = next(m for m in merged if m[0] == "chr22")
    assert chr22[1] == 100
    assert chr22[2] == 8_500
    assert chr22[2] - chr22[1] < _MATE_FETCH_TILE_BP
    assert {q for q, _, _ in chr22[3]} == {"a", "b"}
    chr1 = next(m for m in merged if m[0] == "chr1")
    assert (chr1[1], chr1[2]) == (10, 510)
    assert chr1[2] - chr1[1] == _MATE_FETCH_QUERY_BP


def test_merge_mate_windows_singleton_is_500bp_not_16kb() -> None:
    windows = [("chr22", 1_000, 1_500, "a")]
    merged = _merge_mate_fetch_windows(windows)
    assert merged == [("chr22", 1_000, 1_500, [("a", 1_000, 1_500)])]
    assert merged[0][2] - merged[0][1] == _MATE_FETCH_QUERY_BP


def test_merge_mate_windows_50kb_apart_stay_two_fetches() -> None:
    windows = [
        ("chr22", 1_000, 1_500, "a"),
        ("chr22", 51_000, 51_500, "b"),
    ]
    merged = _merge_mate_fetch_windows(windows)
    assert len(merged) == 2
    assert {(m[1], m[2]) for m in merged} == {(1_000, 1_500), (51_000, 51_500)}
    assert all(m[2] - m[1] == _MATE_FETCH_QUERY_BP for m in merged)


def test_merge_mate_windows_does_not_glue_adjacent_tiles() -> None:
    # 16_000 is tile 0; 17_000 is tile 1. ~1 kb apart, previously gap-merged.
    windows = [
        ("chr22", 16_000, 16_500, "a"),
        ("chr22", 17_000, 17_500, "b"),
    ]
    merged = _merge_mate_fetch_windows(windows)
    assert len(merged) == 2
    assert {(m[1], m[2]) for m in merged} == {(16_000, 16_500), (17_000, 17_500)}


def _write_mate_bam(tmp_path: Path, records: list[tuple[str, str, int, str]]) -> Path:
    import pysam

    raw = tmp_path / "raw.bam"
    bam_path = tmp_path / "mates.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr22", "LN": 1_000_000}, {"SN": "chr1", "LN": 1_000_000}],
    }

    def _write(bam, qname: str, chrom: str, pos: int, seq: str) -> None:
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.query_sequence = seq
        a.flag = 1 | 64
        a.reference_id = bam.get_tid(chrom)
        a.reference_start = pos
        a.next_reference_id = bam.get_tid(chrom)
        a.next_reference_start = pos + 200
        a.mapping_quality = 60
        a.cigarstring = f"{len(seq)}M"
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        bam.write(a)

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        for qname, chrom, pos, seq in records:
            _write(bam, qname, chrom, pos, seq)

    pysam.sort("-o", str(bam_path), str(raw))
    pysam.index(str(bam_path))
    return bam_path


@pytest.mark.skipif(not _have_pysam(), reason="pysam required")
def test_sweep_matches_per_window_fetch(tmp_path: Path) -> None:
    import pysam

    bam_path = _write_mate_bam(
        tmp_path,
        [
            ("near1", "chr22", 1000, "A" * 40),
            ("near2", "chr22", 1400, "C" * 40),
            ("far", "chr22", 500_000, "G" * 40),
            ("other", "chr1", 50, "T" * 40),
        ],
    )
    windows = [
        ("chr22", 1000, 1500, "near1"),
        ("chr22", 1400, 1900, "near2"),
        ("chr22", 500_000, 500_500, "far"),
        ("chr1", 50, 550, "other"),
        ("chr22", 10, 510, "missing"),
    ]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        per_window = _fetch_mates_per_window(bam, windows)
        swept = _fetch_mates_swept(bam, windows)
    assert set(per_window) == set(swept)
    for qname, tup in per_window.items():
        assert swept[qname] == tup
    assert "missing" not in swept
    assert swept["near1"][0] == "A" * 40
    assert swept["other"][0] == "T" * 40


@pytest.mark.skipif(not _have_pysam(), reason="pysam required")
def test_same_tile_one_fetch_finds_both_mates(tmp_path: Path) -> None:
    import pysam

    bam_path = _write_mate_bam(
        tmp_path,
        [
            ("a", "chr22", 100, "A" * 40),
            ("b", "chr22", 8_000, "C" * 40),
        ],
    )
    windows = [
        ("chr22", 100, 600, "a"),
        ("chr22", 8_000, 8_500, "b"),
    ]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        recorder = _FetchRecorder(bam)
        swept = _fetch_mates_swept(recorder, windows)
        per_window = _fetch_mates_per_window(bam, windows)
    assert len(recorder.fetches) == 1
    assert recorder.fetches[0] == ("chr22", 100, 8_500)
    assert recorder.fetches[0][2] - recorder.fetches[0][1] < _MATE_FETCH_TILE_BP
    assert set(swept) == {"a", "b"}
    assert swept == per_window


@pytest.mark.skipif(not _have_pysam(), reason="pysam required")
def test_50kb_apart_two_fetches_not_chained(tmp_path: Path) -> None:
    import pysam

    bam_path = _write_mate_bam(
        tmp_path,
        [
            ("a", "chr22", 1_000, "A" * 40),
            ("b", "chr22", 51_000, "C" * 40),
        ],
    )
    windows = [
        ("chr22", 1_000, 1_500, "a"),
        ("chr22", 51_000, 51_500, "b"),
    ]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        recorder = _FetchRecorder(bam)
        swept = _fetch_mates_swept(recorder, windows)
        per_window = _fetch_mates_per_window(bam, windows)
    assert len(recorder.fetches) == 2
    assert set(recorder.fetches) == {("chr22", 1_000, 1_500), ("chr22", 51_000, 51_500)}
    assert all(end - start == _MATE_FETCH_QUERY_BP for _, start, end in recorder.fetches)
    assert set(swept) == {"a", "b"}
    assert swept == per_window
