"""VNTR-like soft-clips leave MEI-SR and count as VNTR_MAPPED."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _annotate_vntr_like_split_clips,
    _sva_vntr_like_score,
    _split_mei_support_eligible_mask,
)


def _ccctct(n: int = 36) -> str:
    return ("CCCTCT" * ((n // 6) + 2))[:n]


def test_ccctct_score_positive_at_clip_min_len():
    assert _sva_vntr_like_score(_ccctct(20), min_seq_len=20) >= 0.35
    assert _sva_vntr_like_score("ACGTACGTACGTACGTACGT", min_seq_len=20) < 0.35


def test_sva_mapped_vntr_clip_demoted_from_mei_sr():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1100,
                "read_name": "vntr_hit",
                "clip_side": "R",
                "pos": 1050,
                "clip_len": 30,
                "clip_seq": _ccctct(30),
                "mei_hit": True,
                "family": "SVA",
                "target": "SVA_A#Retroposon/SVA",
                "mei_score": 0.9,
            },
            {
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1100,
                "read_name": "sva_body",
                "clip_side": "R",
                "pos": 1050,
                "clip_len": 40,
                "clip_seq": "ACGT" * 10,
                "mei_hit": True,
                "family": "SVA",
                "target": "SVA_A#Retroposon/SVA",
                "mei_score": 0.9,
            },
        ]
    )
    out = _annotate_vntr_like_split_clips(split)
    elig = _split_mei_support_eligible_mask(out)
    by_name = out.set_index("read_name")
    vntr = by_name.loc["vntr_hit"]
    body = by_name.loc["sva_body"]
    assert bool(vntr.vntr_rescue)
    assert not bool(vntr.mei_hit)
    assert vntr.mei_hit_source == "vntr_rescue"
    assert not bool(elig.loc[out["read_name"].eq("vntr_hit")].iloc[0])
    assert not bool(body.vntr_rescue)
    assert bool(body.mei_hit)
    assert bool(elig.loc[out["read_name"].eq("sva_body")].iloc[0])


def test_nonhit_vntr_clip_rescued_at_sva_locus_via_dpe():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 2000,
                "window_end": 2100,
                "read_name": "vntr_miss",
                "clip_side": "R",
                "pos": 2050,
                "clip_len": 24,
                "clip_seq": _ccctct(24),
                "mei_hit": False,
                "family": "",
                "target": "",
                "mei_score": 0.0,
            }
        ]
    )
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 2000,
                "window_end": 2100,
                "read_name": f"dpe{i}",
                "pos": 2040,
                "mei_hit": True,
                "mate_mei_hit": True,
                "family": "SVA",
                "target": "SVA_A#Retroposon/SVA",
                "mei_hit_source": "minimap",
                "vntr_rescue": False,
                "polya_rescue": False,
            }
            for i in range(3)
        ]
    )
    out = _annotate_vntr_like_split_clips(split, discordant_df=disc)
    row = out.iloc[0]
    assert bool(row.vntr_rescue)
    assert not bool(row.mei_hit)
    assert not bool(_split_mei_support_eligible_mask(out).iloc[0])


def test_vntr_clip_not_rescued_without_sva_support():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "orphan_vntr",
                "clip_side": "L",
                "pos": 15,
                "clip_len": 24,
                "clip_seq": _ccctct(24),
                "mei_hit": False,
                "family": "",
                "target": "",
            }
        ]
    )
    out = _annotate_vntr_like_split_clips(split)
    assert not bool(out.iloc[0].vntr_rescue)
