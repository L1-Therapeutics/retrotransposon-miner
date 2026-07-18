"""Short MEI clips count toward SR only when consistent with a strict seed."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _add_candidate_support_info_fields,
    _annotate_short_mei_seed_rescue,
    _build_supporting_reads_detail_table,
    _split_mei_support_eligible_mask,
)


def _empty_mei() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["chrom", "window_start", "window_end", "read_name", "mei_hit"]
    )


def test_short_mei_rescued_when_consistent_with_strict_seed():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "strict1",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 20,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 20,
                "alnlen": 20,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_ok",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_wrong_strand",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "-",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_far",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 200,
                "target_end": 212,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_no_seed_side",
                "clip_side": "R",
                "pos": 112,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
        ]
    )
    out = _annotate_short_mei_seed_rescue(split)
    rescued = set(out.loc[out.short_mei_seed_rescued, "read_name"].astype(str))
    assert rescued == {"short_ok"}
    eligible = set(out.loc[_split_mei_support_eligible_mask(out), "read_name"].astype(str))
    assert eligible == {"strict1", "short_ok"}


def test_short_mei_without_strict_seed_not_eligible():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "only_short",
                "clip_side": "L",
                "pos": 15,
                "clip_len": 15,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 15,
                "alnlen": 15,
            }
        ]
    )
    out = _annotate_short_mei_seed_rescue(split)
    assert not bool(out.iloc[0].short_mei_seed_rescued)
    assert not bool(_split_mei_support_eligible_mask(out).iloc[0])


def _dpe_mei_mates(*, family="ALU", strand="+", mate_start=200, mate_end=250, n=2):
    """Two same-locus L-side DPE mates landing mid-MEI (typical insert-size offset)."""
    rows = []
    for i in range(n):
        rows.append(
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": f"dpe{i}",
                "pos": 108,
                "mei_hit": False,
                "mate_mei_hit": True,
                "mate_mei_start": mate_start,
                "mate_mei_end": mate_end,
                "mate_mei_strand": strand,
                "mate_mei_family": family,
            }
        )
    return pd.DataFrame(rows)


def test_short_mei_rescued_by_dpe_seed_when_no_strict_sr():
    """DPE MEI seeds short clips that sit junction-proximal (no coord overlap required)."""
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_5p",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_wrong_strand",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "-",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_far",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                # Far past DPE cluster along MEI axis (>500 bp gap).
                "target_start": 900,
                "target_end": 912,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_distal_3p",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                # Same side/family but distal to DPE (toward MEI 3' vs +/L entry).
                "target_start": 400,
                "target_end": 412,
                "alnlen": 12,
            },
        ]
    )
    disc = _dpe_mei_mates(mate_start=200, mate_end=250)
    out = _annotate_short_mei_seed_rescue(split, discordant_df=disc)
    rescued = set(out.loc[out.short_mei_seed_rescued, "read_name"].astype(str))
    assert rescued == {"short_5p"}
    eligible = set(out.loc[_split_mei_support_eligible_mask(out), "read_name"].astype(str))
    assert eligible == {"short_5p"}


def test_dpe_seed_ignored_when_strict_sr_seed_exists_on_side():
    """SR≥20 takes priority; DPE must not rescue shorts that fail the SR overlap gate."""
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "strict1",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 22,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 22,
                "alnlen": 22,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_nonoverlap",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                # Would pass DPE proximal gate (near mates at 200–250) but not SR overlap.
                "target_start": 200,
                "target_end": 212,
                "alnlen": 12,
            },
        ]
    )
    disc = _dpe_mei_mates(mate_start=200, mate_end=250)
    out = _annotate_short_mei_seed_rescue(split, discordant_df=disc)
    rescued = set(out.loc[out.short_mei_seed_rescued, "read_name"].astype(str))
    assert rescued == set()
    eligible = set(out.loc[_split_mei_support_eligible_mask(out), "read_name"].astype(str))
    assert eligible == {"strict1"}


def test_support_string_has_no_brk_clp_and_counts_rescued_short():
    candidates = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "insertion_breakpoint_pos": 110,
                "consensus_mei_family": "ALU",
            }
        ]
    )
    split_mei = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "strict1",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 22,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 22,
                "alnlen": 22,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "short_ok",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 12,
                "mei_hit": True,
                "family": "ALU",
                "target_strand": "+",
                "target_start": 1,
                "target_end": 12,
                "alnlen": 12,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "gata_noise",
                "clip_side": "L",
                "pos": 110,
                "clip_len": 25,
                "mei_hit": False,
            },
        ]
    )
    split_mei = _annotate_short_mei_seed_rescue(split_mei)
    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_mei,
        split_control=pd.DataFrame(columns=split_mei.columns),
        discordant_disease=pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "pos"]
        ),
        discordant_control=pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "pos"]
        ),
        split_disease_mei=split_mei.loc[_split_mei_support_eligible_mask(split_mei)].copy(),
        split_control_mei=_empty_mei(),
        discordant_disease_mei=pd.DataFrame(
            [
                {
                    "chrom": "chr22",
                    "window_start": 100,
                    "window_end": 120,
                    "read_name": "dMei1",
                    "pos": 111,
                    "mei_hit": True,
                }
            ]
        ),
        discordant_control_mei=_empty_mei(),
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "BRK_CLP" not in support
    assert "SR_L=2" in support  # strict + rescued short; gata_noise excluded
    assert "MEI_MAPPED=" in support


def test_detail_omits_non_mei_pileup_clips():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "read_name": "noise1",
                "clip_side": "L",
                "pos": 148,
                "clip_len": 30,
                "mei_hit": False,
            },
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 200,
                "read_name": "noise2",
                "clip_side": "L",
                "pos": 148,
                "clip_len": 28,
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
        ]
    )
    detail = _build_supporting_reads_detail_table(
        split_hits=split,
        discordant_hits=pd.DataFrame(),
        discordant_mate_hits=pd.DataFrame(),
        sample="disease",
    )
    names = set(detail.read_name.astype(str))
    assert names == {"poly1"}
