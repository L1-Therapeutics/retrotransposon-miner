"""Sorted mate-window sweep must match per-window bam.fetch results."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from retro_miner.mei_support import _fetch_mates_per_window, _fetch_mates_swept, _merge_mate_fetch_windows


def _have_pysam() -> bool:
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return shutil.which("samtools") is not None or True


def test_merge_mate_windows_joins_nearby_and_splits_far() -> None:
    windows = [
        ("chr22", 100, 600, "a"),
        ("chr22", 700, 1200, "b"),
        ("chr22", 200_000, 200_500, "c"),
        ("chr1", 10, 510, "d"),
    ]
    merged = _merge_mate_fetch_windows(windows, max_gap_bp=2_000)
    assert len(merged) == 3
    chroms = [m[0] for m in merged]
    assert chroms.count("chr22") == 2
    assert "chr1" in chroms
    nearby = next(m for m in merged if m[0] == "chr22" and m[1] == 100)
    assert {q for q, _, _ in nearby[3]} == {"a", "b"}


@pytest.mark.skipif(not _have_pysam(), reason="pysam required")
def test_sweep_matches_per_window_fetch(tmp_path: Path) -> None:
    import pysam

    raw = tmp_path / "raw.bam"
    bam_path = tmp_path / "mates.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr22", "LN": 1_000_000}, {"SN": "chr1", "LN": 1_000_000}],
    }
    seqs = {
        "near1": "A" * 40,
        "near2": "C" * 40,
        "far": "G" * 40,
        "other": "T" * 40,
    }

    def _write(bam, qname: str, chrom: str, pos: int) -> None:
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.query_sequence = seqs[qname]
        a.flag = 1 | 64
        a.reference_id = bam.get_tid(chrom)
        a.reference_start = pos
        a.next_reference_id = bam.get_tid(chrom)
        a.next_reference_start = pos + 200
        a.mapping_quality = 60
        a.cigarstring = f"{len(seqs[qname])}M"
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seqs[qname]))
        bam.write(a)

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        _write(bam, "near1", "chr22", 1000)
        _write(bam, "near2", "chr22", 1400)
        _write(bam, "far", "chr22", 500_000)
        _write(bam, "other", "chr1", 50)

    pysam.sort("-o", str(bam_path), str(raw))
    pysam.index(str(bam_path))

    windows = [
        ("chr22", 1000, 1500, "near1"),
        ("chr22", 1400, 1900, "near2"),
        ("chr22", 500_000, 500_500, "far"),
        ("chr1", 50, 550, "other"),
        ("chr22", 10, 510, "missing"),
    ]
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        per_window = _fetch_mates_per_window(bam, windows)
        swept = _fetch_mates_swept(bam, windows, max_gap_bp=2_000)
    assert set(per_window) == set(swept)
    for qname, tup in per_window.items():
        assert swept[qname] == tup
    assert "missing" not in swept
    assert swept["near1"][0] == "A" * 40
    assert swept["other"][0] == "T" * 40
