"""Clip-alignment and indel-scan caches reload when the evidence keys match."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    ClipAlignmentSummary,
    _load_indel_evidence_cache,
    _load_mei_remap_cache,
    _write_indel_evidence_cache,
    _write_mei_remap_cache,
)


def _split_row(read_name: str = "r1") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr18"],
            "window_start": [10],
            "window_end": [20],
            "read_name": [read_name],
            "mei_hit": [True],
        }
    )


def _disc_row() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr18"],
            "window_start": [10],
            "window_end": [20],
            "read_name": ["d1"],
            "mei_hit": [True],
            "mate_mei_hit": [False],
        }
    )


def test_remap_cache_hits_same_evidence_and_misses_a_new_read(tmp_path):
    split = _split_row()
    disc = _disc_row()
    result = {
        "split_hits": split,
        "split_summary": ClipAlignmentSummary(sample="control", clip_count=4, paf_hits=1),
        "disc_hits": disc,
        "disc_summary": ClipAlignmentSummary(sample="control", clip_count=3, paf_hits=1),
        "disc_mate_summary": ClipAlignmentSummary(sample="control_mate", clip_count=3, paf_hits=0),
    }
    _write_mei_remap_cache(tmp_path, "germline", split, disc, result)

    loaded = _load_mei_remap_cache(tmp_path, "germline", split, disc, sample="control")
    assert loaded is not None
    assert loaded["split_summary"].paf_hits == 1
    assert loaded["split_summary"].clip_count == 4
    assert loaded["disc_mate_summary"].sample == "control_mate"
    assert list(loaded["split_hits"]["read_name"]) == ["r1"]

    changed = _split_row("r9")
    assert _load_mei_remap_cache(tmp_path, "germline", changed, disc, sample="control") is None


def test_indel_cache_hits_same_windows_and_misses_a_shifted_window(tmp_path):
    candidates = pd.DataFrame({"chrom": ["chr18"], "window_start": [10], "window_end": [40]})
    indels = pd.DataFrame(
        {
            "sample": ["control"],
            "chrom": ["chr18"],
            "window_start": [10],
            "window_end": [40],
            "read_name": ["r1"],
            "evidence_type": ["indel"],
            "indel_len": [20],
        }
    )
    _write_indel_evidence_cache(indels, tmp_path, "germline", candidates)

    loaded = _load_indel_evidence_cache(tmp_path, "germline", candidates)
    assert loaded is not None
    assert len(loaded) == 1
    assert int(loaded.iloc[0]["indel_len"]) == 20

    shifted = candidates.copy()
    shifted["window_end"] = [50]
    assert _load_indel_evidence_cache(tmp_path, "germline", shifted) is None
