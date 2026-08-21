"""Parity tests for vectorized midpoint collection in read_architecture.

These tests cover the behavior previously implemented via work.itertuples()
inside _infer_orientation._collect():

* rescue rows are excluded by both `polya_rescue` and `mei_hit_source`
* SR clip-side / anchor-side logic partitions left vs right flank mids
* DPE genomic-position / anchor-side logic partitions left vs right flank mids
"""

from __future__ import annotations

import pandas as pd

from retro_miner.read_architecture import _infer_orientation


def test_infer_orientation_sr_ignores_polya_rescue_rows() -> None:
    """SR midpoint collection excludes polyA-rescue rows before median calling."""
    detail = pd.DataFrame(
        [
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 100,
                "mei_end": 120,
                "clip_side": "R",
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "",
                "window_start": 1000,
                "window_end": 1100,
            },
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 420,
                "mei_end": 440,
                "clip_side": "L",
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "",
                "window_start": 1000,
                "window_end": 1100,
            },
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 900,
                "mei_end": 920,
                "clip_side": "R",
                "anchor_side": "",
                "polya_rescue": True,
                "mei_hit_source": "",
                "window_start": 1000,
                "window_end": 1100,
            },
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 50,
                "mei_end": 70,
                "clip_side": "L",
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "polya_rescue",
                "window_start": 1000,
                "window_end": 1100,
            },
        ]
    )

    got = _infer_orientation(
        detail,
        table_orientation="+",
        window_start=1000,
        window_end=1100,
    )

    assert got == "+"


def test_infer_orientation_sr_uses_anchor_side_when_clip_side_missing() -> None:
    """SR collection preserves the old clip_side-or-anchor_side fallback behavior."""
    detail = pd.DataFrame(
        [
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 100,
                "mei_end": 120,
                "clip_side": "",
                "anchor_side": "R",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
            {
                "evidence_type": "SR",
                "mei_hit": True,
                "mei_start": 420,
                "mei_end": 440,
                "clip_side": "",
                "anchor_side": "L",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
        ]
    )

    got = _infer_orientation(detail, table_orientation="-")

    assert got == "+"


def test_infer_orientation_dpe_uses_genomic_flank_positions() -> None:
    """DPE midpoint collection partitions left/right using genomic_pos vs breakpoint flanks."""
    detail = pd.DataFrame(
        [
            {
                "evidence_type": "DPE",
                "mate_mei_hit": True,
                "mate_mei_start": 110,
                "mate_mei_end": 130,
                "genomic_pos": 100,
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
            {
                "evidence_type": "DPE",
                "mate_mei_hit": True,
                "mate_mei_start": 430,
                "mate_mei_end": 450,
                "genomic_pos": 900,
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
            {
                "evidence_type": "DPE",
                "mate_mei_hit": True,
                "mate_mei_start": 950,
                "mate_mei_end": 970,
                "genomic_pos": 50,
                "anchor_side": "",
                "polya_rescue": False,
                "mei_hit_source": "polya_rescue",
            },
        ]
    )

    got = _infer_orientation(
        detail,
        table_orientation="-",
        breakpoint_left=400,
        breakpoint_right=600,
        breakpoint=500,
    )

    assert got == "+"


def test_infer_orientation_dpe_uses_anchor_side_when_breakpoints_missing() -> None:
    """DPE collection preserves the old anchor_side fallback when genomic breakpoints are unavailable."""
    detail = pd.DataFrame(
        [
            {
                "evidence_type": "DPE",
                "mate_mei_hit": True,
                "mate_mei_start": 700,
                "mate_mei_end": 720,
                "genomic_pos": 0,
                "anchor_side": "L",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
            {
                "evidence_type": "DPE",
                "mate_mei_hit": True,
                "mate_mei_start": 200,
                "mate_mei_end": 220,
                "genomic_pos": 0,
                "anchor_side": "R",
                "polya_rescue": False,
                "mei_hit_source": "",
            },
        ]
    )

    got = _infer_orientation(detail, table_orientation="+")

    assert got == "-"
