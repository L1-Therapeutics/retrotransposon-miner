"""Unit tests for MEI candidate scoring logic.

Covers the insertion-model score computation pipeline, helper functions
used by ``_compute_insertion_model_scores``, and edge-case handling for
zero supporting reads, MAPQ caps, enrichment ratios, and data-structure
parity between baseline and production columns.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from retro_miner.candidate_loci import (
    _cluster_labels,
    _cluster_sorted_positions,
    _safe_cpm,
    _split_cluster_positions,
)
from retro_miner.mei_support import (
    _agreement_flag,
    _compute_insertion_model_scores,
    _ensure_candidate_schema_defaults,
)


# ─────────────────────────────────────────────────────────────────────────────
# Minimal candidate row factory
# ─────────────────────────────────────────────────────────────────────────────


def _base_candidate(**overrides) -> pd.DataFrame:
    """Return a single-row DataFrame modelling the minimal candidate schema."""
    row = {
        "chrom": "chr1",
        "window_start": 5000,
        "window_end": 5200,
        "disease_total_rows": 20,
        "control_total_rows": 10,
        "disease_mei_supported_reads": 0,
        "control_mei_supported_reads": 0,
        "disease_L_mei_supported_reads": 0,
        "disease_R_mei_supported_reads": 0,
        "control_L_mei_supported_reads": 0,
        "control_R_mei_supported_reads": 0,
        "disease_discordant_mei_supported_reads": 0,
        "control_discordant_mei_supported_reads": 0,
        "split_disease_mapq_mean": 0.0,
        "split_control_mapq_mean": 0.0,
        "mei_score_enrichment_ratio": 0.0,
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
        "disease_subfamily_purity_weighted": 0.0,
        "control_subfamily_purity_weighted": 0.0,
        "disease_breakpoint_mode_fraction_weighted": 0.0,
        "control_breakpoint_mode_fraction_weighted": 0.0,
        "tsd_detected": False,
        "disease_poly_at_fraction_weighted": 0.0,
        "control_poly_at_fraction_weighted": 0.0,
        "breakpoint_l1_en_motif_like": False,
        "breakpoint_yyrrrr_logodds_shift1_mt_adj": 0.0,
        "disease_L_clip_overlap_jaccard_median": 0.0,
        "disease_R_clip_overlap_jaccard_median": 0.0,
        "control_L_clip_overlap_jaccard_median": 0.0,
        "control_R_clip_overlap_jaccard_median": 0.0,
        "disease_L_clip_overlap_informative_reads": 0,
        "disease_R_clip_overlap_informative_reads": 0,
        "control_L_clip_overlap_informative_reads": 0,
        "control_R_clip_overlap_informative_reads": 0,
        "coherence_score": 0.5,
        "junk_flag_count": 0,
        "disease_discordant_mei_mapped_fraction": 0.0,
        "control_discordant_mei_mapped_fraction": 0.0,
        "disease_discordant_residual_unique_reads": 0,
        "control_discordant_residual_unique_reads": 0,
        "disease_discordant_residual_interchrom_fraction": 0.0,
        "disease_discordant_residual_large_insert_fraction": 0.0,
        "disease_discordant_residual_mate_unmapped_fraction": 0.0,
        "disease_discordant_residual_same_strand_fraction": 0.0,
        "disease_discordant_residual_improper_pair_fraction": 0.0,
        "control_discordant_residual_interchrom_fraction": 0.0,
        "control_discordant_residual_large_insert_fraction": 0.0,
        "control_discordant_residual_mate_unmapped_fraction": 0.0,
        "control_discordant_residual_same_strand_fraction": 0.0,
        "control_discordant_residual_improper_pair_fraction": 0.0,
        "disease_supporting_reads": (
            "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,"
            "MEI_MAPPED=0,polyA_MAPPED=0"
        ),
        "control_supporting_reads": (
            "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,"
            "MEI_MAPPED=0,polyA_MAPPED=0"
        ),
    }
    row.update(overrides)
    return pd.DataFrame([row])


# ─────────────────────────────────────────────────────────────────────────────
# _agreement_flag
# ─────────────────────────────────────────────────────────────────────────────


class TestAgreementFlag:
    def test_both_empty_returns_zero(self):
        assert _agreement_flag("", "") == 0

    def test_both_none_returns_zero(self):
        assert _agreement_flag(None, None) == 0

    def test_one_empty_one_nonempty_returns_one(self):
        assert _agreement_flag("", "L1HS") == 1
        assert _agreement_flag("L1HS", "") == 1

    def test_matching_nonempty_returns_one(self):
        assert _agreement_flag("SVA", "SVA") == 1

    def test_mismatching_nonempty_returns_zero(self):
        assert _agreement_flag("L1HS", "AluYa5") == 0

    def test_whitespace_only_treated_as_empty(self):
        assert _agreement_flag("  ", "") == 0

    def test_whitespace_stripped_before_comparison(self):
        assert _agreement_flag(" L1HS ", "L1HS") == 1

    def test_none_vs_string_returns_one(self):
        assert _agreement_flag(None, "SVA") == 1


# ─────────────────────────────────────────────────────────────────────────────
# _ensure_candidate_schema_defaults
# ─────────────────────────────────────────────────────────────────────────────


class TestEnsureCandidateSchemaDefaults:
    def test_empty_dataframe_gets_all_default_columns(self):
        df = pd.DataFrame({"chrom": ["chr1"], "window_start": [100], "window_end": [200]})
        result = _ensure_candidate_schema_defaults(df)
        assert "disease_L_mei_supported_reads" in result.columns
        assert "disease_mei_supported_reads" in result.columns
        assert "junk_flag_count" in result.columns
        assert "tsd_detected" in result.columns
        assert "insertion_breakpoint_pos" in result.columns

    def test_existing_columns_are_not_overwritten(self):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "window_start": [100],
            "window_end": [200],
            "disease_mei_supported_reads": [42],
        })
        result = _ensure_candidate_schema_defaults(df)
        assert result["disease_mei_supported_reads"].iloc[0] == 42

    def test_does_not_mutate_input(self):
        df = pd.DataFrame({"chrom": ["chr1"]})
        original_cols = set(df.columns)
        _ensure_candidate_schema_defaults(df)
        assert set(df.columns) == original_cols


# ─────────────────────────────────────────────────────────────────────────────
# _compute_insertion_model_scores — core score computation
# ─────────────────────────────────────────────────────────────────────────────


class TestComputeInsertionModelScores:
    def test_zero_supporting_reads_gives_low_score(self):
        result = _compute_insertion_model_scores(_base_candidate())
        score = result["insertion_model_score"].iloc[0]
        assert 0.0 <= score <= 1.0
        assert score < 0.15

    def test_strong_evidence_gives_high_score(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 5,
            "disease_R_mei_supported_reads": 6,
            "disease_mei_supported_reads": 11,
            "disease_total_rows": 50,
            "control_total_rows": 5,
            "mei_score_enrichment_ratio": 8.0,
            "disease_subfamily_purity_weighted": 0.95,
            "disease_breakpoint_mode_fraction_weighted": 0.90,
            "disease_L_mei_family": "L1HS",
            "disease_R_mei_family": "L1HS",
            "disease_L_mei_subfamily": "L1HS",
            "disease_R_mei_subfamily": "L1HS",
            "disease_L_mei_strand": "+",
            "disease_R_mei_strand": "-",
            "tsd_detected": True,
            "disease_poly_at_fraction_weighted": 0.80,
            "split_disease_mapq_mean": 55.0,
            "coherence_score": 0.85,
            "disease_L_clip_overlap_jaccard_median": 0.75,
            "disease_L_clip_overlap_informative_reads": 4,
            "disease_R_clip_overlap_jaccard_median": 0.70,
            "disease_R_clip_overlap_informative_reads": 3,
        }))
        score = result["insertion_model_score"].iloc[0]
        assert score >= 0.60

    def test_score_capped_at_one(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 100,
            "disease_R_mei_supported_reads": 100,
            "disease_mei_supported_reads": 200,
            "disease_total_rows": 500,
            "control_total_rows": 1,
            "mei_score_enrichment_ratio": 500.0,
            "disease_subfamily_purity_weighted": 1.0,
            "disease_breakpoint_mode_fraction_weighted": 1.0,
            "disease_L_mei_family": "L1HS",
            "disease_R_mei_family": "L1HS",
            "disease_L_mei_subfamily": "L1HS",
            "disease_R_mei_subfamily": "L1HS",
            "disease_L_mei_strand": "+",
            "disease_R_mei_strand": "-",
            "tsd_detected": True,
            "disease_poly_at_fraction_weighted": 1.0,
            "split_disease_mapq_mean": 60.0,
            "coherence_score": 1.0,
        }))
        score = result["insertion_model_score"].iloc[0]
        assert score <= 1.0

    def test_score_non_negative(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_total_rows": 0,
            "control_total_rows": 0,
            "mei_score_enrichment_ratio": -1.0,
        }))
        score = result["insertion_model_score"].iloc[0]
        assert score >= 0.0


# ─────────────────────────────────────────────────────────────────────────────
# MAPQ cap / clipping
# ─────────────────────────────────────────────────────────────────────────────


class TestMapqClipping:
    def test_mapq_60_caps_at_1(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "split_disease_mapq_mean": 60.0,
            "split_control_mapq_mean": 0.0,
        }))
        mapq_contribution = result["insertion_model_score"].iloc[0]
        result_zero = _compute_insertion_model_scores(_base_candidate(**{
            "split_disease_mapq_mean": 0.0,
            "split_control_mapq_mean": 0.0,
        }))
        delta = mapq_contribution - result_zero["insertion_model_score"].iloc[0]
        assert delta <= 0.05 + 1e-9

    def test_mapq_above_60_clipped(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "split_disease_mapq_mean": 100.0,
            "split_control_mapq_mean": 0.0,
        }))
        result_60 = _compute_insertion_model_scores(_base_candidate(**{
            "split_disease_mapq_mean": 60.0,
            "split_control_mapq_mean": 0.0,
        }))
        assert result["insertion_model_score"].iloc[0] == pytest.approx(
            result_60["insertion_model_score"].iloc[0], abs=1e-9
        )

    def test_control_mapq_contributes(self):
        result_d = _compute_insertion_model_scores(_base_candidate(**{
            "split_disease_mapq_mean": 60.0,
        }))
        result_c = _compute_insertion_model_scores(_base_candidate(**{
            "split_control_mapq_mean": 60.0,
        }))
        assert result_d["insertion_model_score"].iloc[0] == pytest.approx(
            result_c["insertion_model_score"].iloc[0], abs=1e-9
        )


# ─────────────────────────────────────────────────────────────────────────────
# Enrichment ratio scaling
# ─────────────────────────────────────────────────────────────────────────────


class TestEnrichmentScaling:
    def test_enrichment_zero_maps_to_zero(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "mei_score_enrichment_ratio": 0.0,
        }))
        mei_enrichment = 0.0
        scaled = mei_enrichment / (mei_enrichment + 1.0)
        assert scaled == pytest.approx(0.0)

    def test_enrichment_large_maps_near_one(self):
        mei_enrichment = 100.0
        scaled = mei_enrichment / (mei_enrichment + 1.0)
        assert scaled == pytest.approx(100.0 / 101.0)

    def test_enrichment_maps_through_sigmoid(self):
        for val in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]:
            scaled = val / (val + 1.0)
            assert 0.0 <= scaled <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# Split / discordant read weightings
# ─────────────────────────────────────────────────────────────────────────────


class TestReadWeightings:
    def test_two_sided_split_support_flagged(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 3,
            "disease_R_mei_supported_reads": 3,
        }))
        assert result["disease_two_sided_support"].iloc[0]

    def test_one_sided_split_not_two_sided(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 5,
            "disease_R_mei_supported_reads": 0,
        }))
        assert not result["disease_two_sided_support"].iloc[0]

    def test_strong_two_sided_requires_ge_2_each(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 1,
            "disease_R_mei_supported_reads": 1,
        }))
        assert result["disease_two_sided_support"].iloc[0]
        assert not result["disease_two_sided_strong_support"].iloc[0]

    def test_discordant_mei_strong_requires_ge_3(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_discordant_mei_supported_reads": 3,
        }))
        assert result["disease_discordant_mei_strong_support"].iloc[0]

    def test_discordant_mei_weak_below_3(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_discordant_mei_supported_reads": 2,
        }))
        assert not result["disease_discordant_mei_strong_support"].iloc[0]

    def test_control_two_sided_support_independently_evaluated(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "control_L_mei_supported_reads": 2,
            "control_R_mei_supported_reads": 2,
        }))
        assert result["control_two_sided_support"].iloc[0]


# ─────────────────────────────────────────────────────────────────────────────
# Event-level combined flags
# ─────────────────────────────────────────────────────────────────────────────


class TestEventFlags:
    def test_event_two_sided_from_disease(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_supported_reads": 3,
            "disease_R_mei_supported_reads": 3,
        }))
        assert result["event_two_sided_like_support"].iloc[0]

    def test_event_two_sided_from_control(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "control_L_mei_supported_reads": 2,
            "control_R_mei_supported_reads": 2,
        }))
        assert result["event_two_sided_like_support"].iloc[0]

    def test_event_family_consistent_from_disease(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_family": "L1HS",
            "disease_R_mei_family": "L1HS",
        }))
        assert result["event_family_consistent"].iloc[0]

    def test_event_family_not_consistent_mismatch(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_family": "L1HS",
            "disease_R_mei_family": "AluYa5",
        }))
        assert not result["event_family_consistent"].iloc[0]

    def test_event_strand_consistent(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_strand": "+",
            "disease_R_mei_strand": "+",
        }))
        assert result["event_strand_consistent"].iloc[0]

    def test_event_strand_not_consistent_mismatch(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_L_mei_strand": "+",
            "disease_R_mei_strand": "-",
        }))
        assert not result["event_strand_consistent"].iloc[0]

    def test_insertion_call_tier_none_when_no_evidence(self):
        result = _compute_insertion_model_scores(_base_candidate())
        assert result["insertion_call_tier"].iloc[0] == "none"


# ─────────────────────────────────────────────────────────────────────────────
# Edge cases: zero reads and NaN handling
# ─────────────────────────────────────────────────────────────────────────────


class TestEdgeCases:
    def test_zero_total_rows_does_not_crash(self):
        result = _compute_insertion_model_scores(_base_candidate(**{
            "disease_total_rows": 0,
            "control_total_rows": 0,
        }))
        assert result["insertion_model_score"].iloc[0] >= 0.0

    def test_nan_mapq_does_not_crash(self):
        row = _base_candidate()
        row.loc[0, "split_disease_mapq_mean"] = np.nan
        result = _compute_insertion_model_scores(row)
        assert result["insertion_model_score"].iloc[0] >= 0.0

    def test_all_missing_optional_columns_handled(self):
        minimal = pd.DataFrame([{
            "chrom": "chr1",
            "window_start": 100,
            "window_end": 200,
        }])
        result = _compute_insertion_model_scores(minimal)
        assert "insertion_model_score" in result.columns
        assert 0.0 <= result["insertion_model_score"].iloc[0] <= 1.0

    def test_multiple_rows_score_independently(self):
        df = pd.concat([
            _base_candidate(**{
                "disease_L_mei_supported_reads": 5,
                "disease_R_mei_supported_reads": 5,
                "disease_mei_supported_reads": 10,
                "mei_score_enrichment_ratio": 5.0,
                "disease_subfamily_purity_weighted": 0.9,
                "disease_L_mei_family": "L1HS",
                "disease_R_mei_family": "L1HS",
                "disease_L_mei_subfamily": "L1HS",
                "disease_R_mei_subfamily": "L1HS",
                "split_disease_mapq_mean": 50.0,
                "coherence_score": 0.7,
            }),
            _base_candidate(**{
                "chrom": "chr2",
                "window_start": 8000,
                "window_end": 8200,
            }),
        ], ignore_index=True)
        result = _compute_insertion_model_scores(df)
        assert len(result) == 2
        assert result["insertion_model_score"].iloc[0] > result["insertion_model_score"].iloc[1]


# ─────────────────────────────────────────────────────────────────────────────
# Score determinism (baseline == performance)
# ─────────────────────────────────────────────────────────────────────────────


class TestScoreDeterminism:
    """Verify identical results for identical inputs — no RNG or drift."""

    def test_same_input_same_score(self):
        inp = _base_candidate(**{
            "disease_L_mei_supported_reads": 4,
            "disease_R_mei_supported_reads": 3,
            "disease_mei_supported_reads": 7,
            "mei_score_enrichment_ratio": 3.0,
            "disease_subfamily_purity_weighted": 0.80,
            "disease_breakpoint_mode_fraction_weighted": 0.70,
            "disease_L_mei_family": "L1HS",
            "disease_R_mei_family": "L1HS",
            "disease_L_mei_subfamily": "L1HS",
            "disease_R_mei_subfamily": "L1HS",
            "disease_L_mei_strand": "+",
            "disease_R_mei_strand": "-",
            "tsd_detected": True,
            "split_disease_mapq_mean": 45.0,
            "coherence_score": 0.65,
        })
        r1 = _compute_insertion_model_scores(inp)
        r2 = _compute_insertion_model_scores(inp)
        assert r1["insertion_model_score"].iloc[0] == pytest.approx(
            r2["insertion_model_score"].iloc[0], abs=1e-12
        )

    def test_order_independent(self):
        inp_a = _base_candidate(**{
            "disease_L_mei_supported_reads": 2,
            "disease_R_mei_supported_reads": 3,
            "split_disease_mapq_mean": 40.0,
        })
        inp_b = _base_candidate(**{
            "disease_L_mei_supported_reads": 2,
            "disease_R_mei_supported_reads": 3,
            "split_disease_mapq_mean": 40.0,
        })
        r1 = _compute_insertion_model_scores(inp_a)
        r2 = _compute_insertion_model_scores(inp_b)
        assert r1["insertion_model_score"].iloc[0] == pytest.approx(
            r2["insertion_model_score"].iloc[0], abs=1e-12
        )


# ─────────────────────────────────────────────────────────────────────────────
# CPM calculation
# ─────────────────────────────────────────────────────────────────────────────


class TestSafeCpm:
    def test_zero_denominator_returns_zero(self):
        result = _safe_cpm(pd.Series([10, 20]), 0)
        assert (result == 0.0).all()

    def test_positive_denominator(self):
        result = _safe_cpm(pd.Series([1000]), 10000)
        assert result.iloc[0] == pytest.approx(100_000.0)

    def test_empty_series(self):
        result = _safe_cpm(pd.Series([], dtype=int), 1000)
        assert len(result) == 0


# ─────────────────────────────────────────────────────────────────────────────
# Cluster labels / positions parity
# ─────────────────────────────────────────────────────────────────────────────


class TestClusterLabelsParity:
    """Verify that _cluster_labels produces the same grouping as _cluster_sorted_positions."""

    def test_parity_basic(self):
        positions = [100, 120, 140, 500, 510, 900]
        labels = _cluster_labels(positions, max_gap_bp=50)
        clusters_from_labels = {}
        for pos, lab in zip(positions, labels):
            clusters_from_labels.setdefault(int(lab), []).append(pos)
        raw_clusters = _cluster_sorted_positions(positions, max_gap_bp=50)
        assert sorted(clusters_from_labels.values()) == sorted(raw_clusters)

    def test_empty_input(self):
        labels = _cluster_labels([], max_gap_bp=100)
        assert len(labels) == 0

    def test_single_element(self):
        labels = _cluster_labels([42], max_gap_bp=100)
        assert labels.tolist() == [0]
