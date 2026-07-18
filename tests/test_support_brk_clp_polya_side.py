"""Supporting-reads strings: no BRK_CLP; polyA_side still reported."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import _add_candidate_support_info_fields


def _empty_mei() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["chrom", "window_start", "window_end", "read_name", "mei_hit"]
    )


def test_untagged_softclips_do_not_count_as_sr():
    """Soft-clips without MEI mapping must not inflate SR."""
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
    split_disease = pd.DataFrame(
        [
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": f"dL{i}", "clip_side": "L", "pos": 108}
            for i in range(5)
        ]
        + [
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": "dR1", "clip_side": "R", "pos": 112},
        ]
    )
    disc_disease_mei = pd.DataFrame(
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
    )
    disc_control_mei = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "cMei1",
                "pos": 111,
                "mei_hit": True,
            }
        ]
    )

    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_disease,
        split_control=pd.DataFrame(columns=split_disease.columns),
        discordant_disease=disc_disease_mei.copy(),
        discordant_control=disc_control_mei.copy(),
        split_disease_mei=_empty_mei(),
        split_control_mei=_empty_mei(),
        discordant_disease_mei=disc_disease_mei,
        discordant_control_mei=disc_control_mei,
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "BRK_CLP" not in support
    assert "SR_L=0" in support
    assert "SR_R=0" in support


def test_polya_clip_with_mei_hit_counts_polya_only_not_mei_mapped():
    """PolyA soft-clips must not also inflate MEI_MAPPED."""
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
    # Pretend a polyA clip also remapped to Alu 3' (legacy double-count path).
    polya_mei = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "poly1",
                "clip_side": "L",
                "pos": 108,
                "clip_len": 40,
                "clip_poly_at_run": 30,
                "poly_tail_rescued": True,
                "mei_hit": True,
                "target": "AluSq4#SINE/Alu",
                "target_start": 272,
                "target_end": 311,
            }
        ]
    )
    real_mei = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 100,
                "window_end": 120,
                "read_name": "mei1",
                "clip_side": "R",
                "pos": 112,
                "clip_len": 40,
                "clip_poly_at_run": 0,
                "poly_tail_rescued": False,
                "mei_hit": True,
                "target": "AluSq4#SINE/Alu",
                "target_start": 1,
                "target_end": 40,
            }
        ]
    )
    from retro_miner.mei_support import (
        _demote_polya_split_mei_hits,
        _split_mei_support_eligible_mask,
    )

    split_all = pd.concat([polya_mei, real_mei], ignore_index=True)
    split_all = _demote_polya_split_mei_hits(split_all)
    split_mei = split_all.loc[_split_mei_support_eligible_mask(split_all)].copy()
    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_all,
        split_control=pd.DataFrame(columns=split_all.columns),
        discordant_disease=pd.DataFrame(columns=["chrom", "window_start", "window_end", "read_name"]),
        discordant_control=pd.DataFrame(columns=["chrom", "window_start", "window_end", "read_name"]),
        split_disease_mei=split_mei,
        split_control_mei=_empty_mei(),
        discordant_disease_mei=_empty_mei(),
        discordant_control_mei=_empty_mei(),
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "polyA_MAPPED=1" in support
    assert "MEI_MAPPED=1" in support
    assert "SR_L=0" in support  # polyA side not SR
    assert "SR_R=1" in support


def test_polya_side_majority_from_split_clips():
    candidates = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 50,
                "window_end": 60,
                "insertion_breakpoint_pos": 55,
                "consensus_mei_family": "LINE1",
            }
        ]
    )
    split_disease = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 50,
                "window_end": 60,
                "read_name": f"pR{i}",
                "clip_side": "R",
                "pos": 56,
                "clip_poly_at_run": 20,
                "poly_tail_rescued": True,
            }
            for i in range(3)
        ]
        + [
            {
                "chrom": "chr1",
                "window_start": 50,
                "window_end": 60,
                "read_name": "pL1",
                "clip_side": "L",
                "pos": 54,
                "clip_poly_at_run": 12,
                "poly_tail_rescued": True,
            }
        ]
    )
    disc_mei = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 50,
                "window_end": 60,
                "read_name": "mei1",
                "pos": 55,
                "mei_hit": True,
            }
        ]
    )
    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_disease,
        split_control=pd.DataFrame(columns=split_disease.columns),
        discordant_disease=disc_mei,
        discordant_control=pd.DataFrame(columns=disc_mei.columns),
        split_disease_mei=_empty_mei(),
        split_control_mei=_empty_mei(),
        discordant_disease_mei=disc_mei,
        discordant_control_mei=_empty_mei(),
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "polyA_side=R" in support
    assert "polyA_MAPPED=4" in support
    assert "BRK_CLP" not in support
    # PolyA is not SR.
    assert "SR_L=0" in support
    assert "SR_R=0" in support


def test_indel_evidence_does_not_count_as_sr():
    """CIGAR indels must not inflate SR; only MEI-mapped splits do."""
    candidates = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 200,
                "window_end": 220,
                "insertion_breakpoint_pos": 210,
                "consensus_mei_family": "LINE1",
            }
        ]
    )
    split_disease = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 200,
                "window_end": 220,
                "read_name": "mei_sr1",
                "clip_side": "L",
                "pos": 209,
                "clip_len": 25,
                "mei_hit": True,
            }
        ]
    )
    indel_disease = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 200,
                "window_end": 220,
                "read_name": f"indel{i}",
                "clip_side": "R",
                "pos": 211,
                "clip_len": 40,
                "evidence_type": "indel",
                "indel_type": "I",
                "indel_len": 40,
            }
            for i in range(10)
        ]
    )
    disc_mei = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 200,
                "window_end": 220,
                "read_name": "dMei1",
                "pos": 211,
                "mei_hit": True,
            }
        ]
    )
    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_disease,
        split_control=pd.DataFrame(columns=split_disease.columns),
        discordant_disease=disc_mei,
        discordant_control=pd.DataFrame(columns=disc_mei.columns),
        split_disease_mei=split_disease.copy(),
        split_control_mei=_empty_mei(),
        discordant_disease_mei=disc_mei,
        discordant_control_mei=_empty_mei(),
        indel_disease=indel_disease,
        indel_control=pd.DataFrame(columns=indel_disease.columns),
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "SR_L=1" in support
    assert "SR_R=0" in support
    assert "MEI_MAPPED=" in support
