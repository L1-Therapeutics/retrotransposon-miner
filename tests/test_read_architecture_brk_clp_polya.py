"""Read-architecture plots include MEI and polyA supporting splits (no BRK_CLP)."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _build_supporting_reads_detail_table,
    _split_polya_member_mask,
)
from retro_miner.read_architecture import (
    LocusLayout,
    _filter_detail_for_locus,
    _pair_from_sr_row,
    _pair_segments,
)


def _layout() -> LocusLayout:
    return LocusLayout(
        chrom="chr22",
        window_start=100,
        window_end=200,
        breakpoint=150,
        breakpoint_left=145,
        breakpoint_right=155,
        mei_5p=1,
        mei_3p=300,
        flank_bp=200,
        mei_span_bp=300,
        orientation="+",
        polya_zone_bp=40,
    )


def test_split_polya_mask():
    split = pd.DataFrame(
        [
            {"chrom": "chr22", "window_start": 100, "window_end": 200, "read_name": "p1", "clip_side": "R", "pos": 152, "clip_len": 40, "poly_tail_rescued": True, "clip_poly_at_run": 20},
            {"chrom": "chr22", "window_start": 100, "window_end": 200, "read_name": "n1", "clip_side": "L", "pos": 148, "clip_len": 30},
        ]
    )
    polya = _split_polya_member_mask(split)
    assert bool(polya.loc[split.read_name.eq("p1").idxmax()])
    assert not bool(polya.loc[split.read_name.eq("n1").idxmax()])


def test_detail_emits_polya_without_mei_hit_but_not_bare_pileup():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "read_name": "noise",
                "clip_side": "L",
                "pos": 148,
                "clip_len": 30,
                "mei_hit": False,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "read_name": "poly1",
                "clip_side": "R",
                "pos": 152,
                "clip_len": 40,
                "mei_hit": False,
                "poly_tail_rescued": True,
                "clip_poly_at_run": 20,
                "clip_poly_base": "A",
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "read_name": "mei1",
                "clip_side": "L",
                "pos": 149,
                "clip_len": 50,
                "mei_hit": True,
                "target": "AluY#SINE/Alu",
                "target_start": 10,
                "target_end": 60,
                "target_strand": "+",
            },
        ]
    )
    detail = _build_supporting_reads_detail_table(
        split_hits=split,
        discordant_hits=pd.DataFrame(),
        discordant_mate_hits=pd.DataFrame(),
        sample="disease",
    )
    names = set(detail.read_name.astype(str))
    assert names == {"poly1", "mei1"}


def test_sr_polya_requires_run_or_rescue_not_bare_base():
    layout = _layout()
    bare = pd.Series(
        {
            "evidence_type": "SR",
            "read_name": "bare",
            "genomic_pos": 148,
            "anchor_side": "L",
            "clip_len": 30,
            "mei_hit": False,
            "polya_rescue": False,
            "poly_tail_rescued": False,
            "clip_poly_at_run": 2,
            "clip_poly_base": "T",
            "polya_base": "T",
            "mei_start": 0,
            "mei_end": 0,
        }
    )
    from retro_miner.read_architecture import _sr_is_polya_clip

    assert _sr_is_polya_clip(bare, layout) is False
    assert _pair_from_sr_row(bare, layout) is None


def test_pair_from_sr_builds_polya_without_mei():
    layout = _layout()
    poly = pd.Series(
        {
            "evidence_type": "SR",
            "read_name": "poly1",
            "genomic_pos": 152,
            "anchor_side": "R",
            "clip_len": 40,
            "mei_hit": False,
            "polya_rescue": True,
            "poly_tail_rescued": True,
            "clip_poly_at_run": 20,
            "polya_base": "A",
            "mei_start": 0,
            "mei_end": 0,
        }
    )
    poly_pair = _pair_from_sr_row(poly, layout)
    assert poly_pair is not None
    assert poly_pair["remote_kind"] == "sr_polya"
    assert poly_pair["rescue_kind"] == "polya"


def test_filter_and_pair_segments_keep_mei_polya():
    layout = _layout()
    detail = pd.DataFrame(
        [
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "mei1",
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "anchor_side": "L",
                "genomic_pos": 148,
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 30,
                "mei_start": 10,
                "mei_end": 40,
            },
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "poly1",
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "anchor_side": "R",
                "genomic_pos": 152,
                "mei_hit": False,
                "mate_mei_hit": False,
                "polya_rescue": True,
                "poly_tail_rescued": True,
                "clip_len": 40,
                "clip_poly_at_run": 20,
                "polya_base": "A",
                "mei_start": 0,
                "mei_end": 0,
            },
        ]
    )
    filtered = _filter_detail_for_locus(
        detail,
        chrom="chr22",
        window_start=100,
        window_end=200,
        sample="disease",
        breakpoint=150,
    )
    assert len(filtered) == 2
    pairs, stats = _pair_segments(filtered, layout, max_pairs=50)
    kinds = {p["remote_kind"] for p in pairs}
    assert "sr_mei" in kinds
    assert "sr_polya" in kinds
    assert stats["sr_plotted"] == 2
