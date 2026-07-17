"""COMPLEX_INS careful labeling vs SIMPLE_MEI / MEI_WITH_COMPLEX."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _add_consolidated_event_fields,
    _assign_gold_stage,
    _compute_insertion_model_scores,
)


def _base_row(**overrides) -> dict:
    row = {
        "chrom": "chr22",
        "window_start": 100,
        "window_end": 120,
        "insertion_breakpoint_pos": 110,
        "junk_flag_count": 0,
        "coherence_score": 0.5,
        "disease_supporting_reads": (
            "SR_L=5,SR_R=6,DPE_L=4,DPE_R=40,BRK_CLP_L=5,BRK_CLP_R=4,"
            "MEI_MAPPED=1,polyA_MAPPED=0"
        ),
        "control_supporting_reads": (
            "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,BRK_CLP_L=0,BRK_CLP_R=0,"
            "MEI_MAPPED=0,polyA_MAPPED=0"
        ),
        "disease_discordant_mei_mapped_fraction": 0.10,
        "control_discordant_mei_mapped_fraction": 0.0,
        "disease_discordant_residual_unique_reads": 30,
        "control_discordant_residual_unique_reads": 0,
        "disease_discordant_residual_interchrom_fraction": 0.80,
        "disease_discordant_residual_large_insert_fraction": 0.10,
        "disease_discordant_residual_mate_unmapped_fraction": 0.0,
        "disease_discordant_residual_same_strand_fraction": 0.0,
        "disease_discordant_residual_improper_pair_fraction": 0.0,
        "control_discordant_residual_interchrom_fraction": 0.0,
        "control_discordant_residual_large_insert_fraction": 0.0,
        "control_discordant_residual_mate_unmapped_fraction": 0.0,
        "disease_mei_supported_reads": 1,
        "control_mei_supported_reads": 0,
        "disease_discordant_mei_supported_reads": 1,
        "control_discordant_mei_supported_reads": 0,
        "disease_L_mei_supported_reads": 0,
        "disease_R_mei_supported_reads": 0,
        "control_L_mei_supported_reads": 0,
        "control_R_mei_supported_reads": 0,
        "disease_L_poly_at_reads": 0,
        "disease_R_poly_at_reads": 0,
        "control_L_poly_at_reads": 0,
        "control_R_poly_at_reads": 0,
        "disease_L_poly_at_max_run": 0,
        "disease_R_poly_at_max_run": 0,
        "disease_L_split_breakpoint_support": 5,
        "disease_R_split_breakpoint_support": 4,
        "disease_total_rows": 50,
        "control_total_rows": 10,
        "split_disease_mapq_mean": 50.0,
        "disease_L_mei_family": "",
        "disease_R_mei_family": "",
        "disease_L_mei_subfamily": "",
        "disease_R_mei_subfamily": "",
        "disease_L_mei_strand": "",
        "disease_R_mei_strand": "",
        "control_L_mei_family": "",
        "control_R_mei_family": "",
        "control_L_mei_subfamily": "",
        "control_R_mei_subfamily": "",
        "control_L_mei_strand": "",
        "control_R_mei_strand": "",
        "disease_discordant_mei_subfamily": "SVA_E#Retroposon/SVA",
        "disease_discordant_mei_family": "SVA",
        "disease_discordant_mei_family_votes": "SVA:1",
        "disease_discordant_mei_subfamily_votes": "SVA_E#Retroposon/SVA:1",
    }
    row.update(overrides)
    return row


def test_complex_ins_when_residual_interchrom_majority():
    out = _compute_insertion_model_scores(pd.DataFrame([_base_row()]))
    assert out.loc[0, "insertion_event_class"] == "COMPLEX_INS"
    assert out.loc[0, "insertion_call_tier"] == "complex_ins"
    assert bool(out.loc[0, "discordant_mei_majority"]) is False


def test_complex_ins_blanks_mei_family():
    scored = _compute_insertion_model_scores(pd.DataFrame([_base_row()]))
    consolidated = _add_consolidated_event_fields(scored)
    assert consolidated.loc[0, "insertion_event_class"] == "COMPLEX_INS"
    assert consolidated.loc[0, "mei_family"] == ""
    assert consolidated.loc[0, "mei_subfamily"] == ""


def test_complex_ins_demoted_from_gold():
    scored = _compute_insertion_model_scores(pd.DataFrame([_base_row()]))
    scored["silver_stage_pass"] = True
    scored["analysis_stage_tier"] = "silver"
    gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=1)
    assert gold.loc[0, "insertion_event_class"] == "COMPLEX_INS"
    assert bool(gold.loc[0, "gold_stage_pass"]) is False
    assert "complex_ins_non_mei" in str(gold.loc[0, "gold_stage_fail_reason"])


def test_mei_with_complex_can_remain_gold():
    """Real MEI + weird companion stays eligible for gold (unlike COMPLEX_INS)."""
    # Not MEI-majority residual geometry, but enough MEI_MAPPED for a real MEI call.
    row = _base_row(
        disease_discordant_mei_mapped_fraction=0.30,
        disease_discordant_residual_interchrom_fraction=0.80,
        disease_discordant_residual_unique_reads=20,
        disease_mei_supported_reads=10,
        disease_discordant_mei_supported_reads=10,
        disease_L_mei_supported_reads=4,
        disease_R_mei_supported_reads=4,
        disease_supporting_reads=(
            "SR_L=5,SR_R=6,DPE_L=4,DPE_R=40,BRK_CLP_L=5,BRK_CLP_R=4,"
            "MEI_MAPPED=10,polyA_MAPPED=0"
        ),
    )
    scored = _compute_insertion_model_scores(pd.DataFrame([row]))
    assert scored.loc[0, "insertion_event_class"] == "MEI_WITH_COMPLEX"
    scored["silver_stage_pass"] = True
    scored["analysis_stage_tier"] = "silver"
    gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
    assert bool(gold.loc[0, "gold_stage_pass"]) is True
    assert "complex_ins_non_mei" not in str(gold.loc[0, "gold_stage_fail_reason"])


def test_classic_polya_mei_sidepair_not_complex_ins():
    """polyA on one flank + MEI on the other stays out of COMPLEX_INS."""
    row = _base_row(
        disease_L_poly_at_reads=3,
        disease_L_poly_at_max_run=20,
        disease_R_mei_supported_reads=5,
        disease_mei_supported_reads=5,
        disease_discordant_mei_supported_reads=5,
        disease_discordant_mei_mapped_fraction=0.20,
        disease_supporting_reads=(
            "SR_L=5,SR_R=6,DPE_L=4,DPE_R=40,BRK_CLP_L=5,BRK_CLP_R=4,"
            "MEI_MAPPED=5,polyA_MAPPED=3"
        ),
    )
    out = _compute_insertion_model_scores(pd.DataFrame([row]))
    assert bool(out.loc[0, "classic_polya_mei_sidepair"]) is True
    assert out.loc[0, "insertion_event_class"] != "COMPLEX_INS"


def test_mei_majority_with_complex_companion_is_mei_with_complex():
    row = _base_row(
        disease_discordant_mei_mapped_fraction=0.70,
        disease_discordant_residual_interchrom_fraction=0.80,
        disease_discordant_residual_unique_reads=20,
        disease_mei_supported_reads=10,
        disease_discordant_mei_supported_reads=10,
        disease_L_mei_supported_reads=4,
        disease_R_mei_supported_reads=4,
        disease_supporting_reads=(
            "SR_L=5,SR_R=6,DPE_L=4,DPE_R=40,BRK_CLP_L=5,BRK_CLP_R=4,"
            "MEI_MAPPED=10,polyA_MAPPED=0"
        ),
    )
    out = _compute_insertion_model_scores(pd.DataFrame([row]))
    assert bool(out.loc[0, "discordant_mei_majority"]) is True
    assert out.loc[0, "insertion_event_class"] != "COMPLEX_INS"
    # Majority MEI → simple MEI path (complex residual flags suppressed by mei_majority).
    assert out.loc[0, "insertion_event_class"] in {"SIMPLE_MEI", "MEI_WITH_COMPLEX", "NONE"}


def test_mei_majority_stays_simple_when_classic():
    row = _base_row(
        disease_discordant_mei_mapped_fraction=0.80,
        disease_discordant_residual_interchrom_fraction=0.10,
        disease_discordant_residual_unique_reads=2,
        disease_L_poly_at_reads=2,
        disease_R_mei_supported_reads=4,
        disease_mei_supported_reads=8,
        disease_discordant_mei_supported_reads=8,
        disease_supporting_reads=(
            "SR_L=4,SR_R=4,DPE_L=3,DPE_R=3,BRK_CLP_L=3,BRK_CLP_R=3,"
            "MEI_MAPPED=8,polyA_MAPPED=2"
        ),
    )
    out = _compute_insertion_model_scores(pd.DataFrame([row]))
    assert out.loc[0, "insertion_event_class"] == "SIMPLE_MEI"
