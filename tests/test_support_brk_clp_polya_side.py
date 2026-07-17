"""BRK_CLP pileup counts and polyA_side in supporting-reads strings."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import _add_candidate_support_info_fields


def _empty_mei() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["chrom", "window_start", "window_end", "read_name", "mei_hit"]
    )


def test_brk_clp_counts_modal_pileup_and_skips_strict_gate():
    """Clips near the mode count as BRK_CLP even when SR is strict-filtered."""
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
    # Disease: five L clips at mode 108 (±2), one outlier at 130.
    # Control has MEI so disease can enter strict mode; disease MEI is weak
    # (only disc MEI ≤2, no split MEI / poly split) → SR uses strict reads only.
    split_disease = pd.DataFrame(
        [
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": f"dL{i}", "clip_side": "L", "pos": 108}
            for i in range(5)
        ]
        + [
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": "dL_out", "clip_side": "L", "pos": 130},
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": "dR1", "clip_side": "R", "pos": 112},
            {"chrom": "chr22", "window_start": 100, "window_end": 120, "read_name": "dR2", "clip_side": "R", "pos": 113},
        ]
    )
    # Only one disease disc MEI read → weak_mei_only_discordant; control has MEI.
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
    disc_disease = disc_disease_mei.copy()
    disc_control = disc_control_mei.copy()

    out = _add_candidate_support_info_fields(
        candidates,
        split_disease=split_disease,
        split_control=pd.DataFrame(columns=split_disease.columns),
        discordant_disease=disc_disease,
        discordant_control=disc_control,
        split_disease_mei=_empty_mei(),
        split_control_mei=_empty_mei(),
        discordant_disease_mei=disc_disease_mei,
        discordant_control_mei=disc_control_mei,
    )
    support = str(out.iloc[0]["disease_supporting_reads"])
    assert "BRK_CLP_L=5" in support  # modal pileup; outlier excluded
    assert "BRK_CLP_R=2" in support
    # Strict mode: disease SR should not credit untagged clips.
    assert "SR_L=0" in support
    assert "SR_R=0" in support


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
    # Give the locus MEI support so counts are emitted.
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
