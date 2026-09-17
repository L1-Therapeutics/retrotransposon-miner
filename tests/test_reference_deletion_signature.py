"""Reference-deletion DPE clustering and gold-stage safeguards."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _add_reference_deletion_signature,
    _aggregate_reference_deletion_dpe_metrics,
    _assign_gold_stage,
)


def _deletion_pairs(
    *,
    count: int = 20,
    gap_bp: int = 10_000,
    mate_step: int = 20,
    mate_chrom: str = "chr22",
) -> pd.DataFrame:
    rows = []
    for idx in range(count):
        pos = 100 + idx
        rows.append(
            {
                "chrom": "chr22",
                "window_start": 90,
                "window_end": 150,
                "pos": pos,
                "mate_chrom": mate_chrom,
                "mate_pos": pos + gap_bp + idx * mate_step,
                "read_name": f"pair_{idx}",
                "is_reverse": False,
                "mate_is_reverse": True,
                "discordant_reasons": "improper_pair,large_insert",
                "mei_hit": True,
            }
        )
    return pd.DataFrame(rows)


def _signature_row(**overrides: object) -> pd.DataFrame:
    row: dict[str, object] = {
        "disease_reference_deletion_cluster_reads": 20,
        "disease_reference_deletion_cluster_fraction": 1.0,
        "disease_reference_deletion_cluster_width_bp": 400,
        "disease_reference_deletion_span_bp": 10_000,
        "disease_reference_deletion_anchor_side": "L",
        "disease_L_mei_supported_reads": 20,
        "disease_R_mei_supported_reads": 2,
        "disease_family_agreement": 1,
        "disease_strand_agreement": 1,
        "control_reference_deletion_cluster_reads": 0,
        "control_reference_deletion_cluster_fraction": 0.0,
        "control_reference_deletion_cluster_width_bp": 0,
        "control_reference_deletion_span_bp": 0,
        "control_reference_deletion_anchor_side": "",
        "event_clip_overlap_consistency": 0.0,
        "classic_polya_mei_sidepair": False,
    }
    row.update(overrides)
    return pd.DataFrame([row])


def test_clusters_ten_kb_same_chromosome_fr_pairs() -> None:
    out = _aggregate_reference_deletion_dpe_metrics(_deletion_pairs(), "disease")
    assert int(out.loc[0, "disease_reference_deletion_cluster_reads"]) == 20
    assert float(out.loc[0, "disease_reference_deletion_cluster_fraction"]) == 1.0
    assert int(out.loc[0, "disease_reference_deletion_cluster_width_bp"]) <= 1000
    assert float(out.loc[0, "disease_reference_deletion_span_bp"]) >= 10_000
    assert out.loc[0, "disease_reference_deletion_anchor_side"] == "L"


def test_clusters_hundred_kb_deletion_without_upper_span_limit() -> None:
    out = _aggregate_reference_deletion_dpe_metrics(
        _deletion_pairs(gap_bp=100_000),
        "disease",
    )
    assert int(out.loc[0, "disease_reference_deletion_cluster_reads"]) == 20
    assert float(out.loc[0, "disease_reference_deletion_span_bp"]) >= 100_000


def test_diffuse_or_interchromosomal_mates_do_not_form_deletion_cluster() -> None:
    diffuse = _aggregate_reference_deletion_dpe_metrics(
        _deletion_pairs(count=8, mate_step=1500),
        "disease",
    )
    assert int(diffuse.loc[0, "disease_reference_deletion_cluster_reads"]) == 1
    assert float(diffuse.loc[0, "disease_reference_deletion_cluster_fraction"]) < 0.70

    interchrom = _aggregate_reference_deletion_dpe_metrics(
        _deletion_pairs(mate_chrom="chr2"),
        "disease",
    )
    assert interchrom.empty


def test_one_sided_cluster_sets_reference_deletion_signature() -> None:
    out = _add_reference_deletion_signature(_signature_row())
    assert bool(out.loc[0, "disease_reference_deletion_signature"])
    assert bool(out.loc[0, "reference_deletion_signature"])
    assert int(out.loc[0, "disease_reference_deletion_opposite_mei_reads"]) == 2
    assert int(out.loc[0, "disease_reference_deletion_allowed_opposite_reads"]) == 4


def test_independent_bilateral_junction_protects_real_insertion() -> None:
    out = _add_reference_deletion_signature(
        _signature_row(
            disease_R_mei_supported_reads=3,
            event_clip_overlap_consistency=0.50,
        )
    )
    assert bool(out.loc[0, "reference_deletion_independent_insertion_junction"])
    assert not bool(out.loc[0, "reference_deletion_signature"])


def test_classic_polya_sidepair_protects_real_insertion() -> None:
    out = _add_reference_deletion_signature(
        _signature_row(classic_polya_mei_sidepair=True)
    )
    assert bool(out.loc[0, "reference_deletion_independent_insertion_junction"])
    assert not bool(out.loc[0, "reference_deletion_signature"])


def test_reference_deletion_is_downgraded_from_gold() -> None:
    row = _add_reference_deletion_signature(_signature_row()).iloc[0].to_dict()
    row.update(
        {
            "silver_stage_pass": True,
            "analysis_stage_tier": "silver",
            "disease_mei_mapped": 20,
            "control_mei_mapped": 0,
            "insertion_event_class": "SIMPLE_MEI",
            "stage_fail_reason": "",
        }
    )
    out = _assign_gold_stage(pd.DataFrame([row]), empirical_stage=False)
    assert not bool(out.loc[0, "gold_stage_pass"])
    assert out.loc[0, "analysis_stage_tier"] == "silver"
    assert "reference_deletion_signature" in out.loc[0, "gold_stage_fail_reason"]

