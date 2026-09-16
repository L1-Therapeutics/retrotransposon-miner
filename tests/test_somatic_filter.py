"""Tests for the somatic vs. germline / PON odds-ratio filter.

Verifies that Fisher exact p-value discrimination separates genuine low-VAF
somatic L1 insertion calls from germline alleles and panel-of-normals
sequencing artifacts, per the Evrony/iskow somatic VAF rationale.
"""

import pytest

from retro_miner.somatic_filter import (
    CONTROL_PANEL_ARTIFACT,
    GERMLINE_PASS,
    SOMATIC_CANDIDATE,
    SomaticFilterResult,
    evaluate_somatic_odds,
)


class TestEvaluateSomaticOdds:
    def test_genuine_somatic_call_at_low_vaf(self):
        result = evaluate_somatic_odds(k_alt=5, k_ref=35, c_alt=0, c_ref=120)

        assert isinstance(result, SomaticFilterResult)
        assert result.is_somatic is True
        assert result.classification == SOMATIC_CANDIDATE
        assert result.p_value < 0.001
        assert result.odds_ratio > 10

    def test_germline_allele_present_in_control_panel_flagged_as_artifact(self):
        result = evaluate_somatic_odds(k_alt=15, k_ref=15, c_alt=12, c_ref=15)

        assert result.is_somatic is False
        assert result.classification == CONTROL_PANEL_ARTIFACT
        assert result.p_value == pytest.approx(0.7921, abs=1e-4)
        assert result.odds_ratio == pytest.approx(1.25)

    def test_control_alt_reads_take_precedence_over_somatic_signal(self):
        result = evaluate_somatic_odds(k_alt=8, k_ref=32, c_alt=3, c_ref=30)

        assert result.classification == CONTROL_PANEL_ARTIFACT
        assert result.is_somatic is False

    def test_low_alt_support_not_somatic_with_clean_panel(self):
        result = evaluate_somatic_odds(k_alt=2, k_ref=38, c_alt=0, c_ref=40)

        assert result.classification == GERMLINE_PASS
        assert result.is_somatic is False

    def test_single_control_alt_read_is_inconclusive(self):
        result = evaluate_somatic_odds(k_alt=10, k_ref=30, c_alt=1, c_ref=30)

        assert result.classification == GERMLINE_PASS
        assert result.is_somatic is False

    def test_default_control_parameters_match_explicit_clean_panel(self):
        via_defaults = evaluate_somatic_odds(k_alt=5, k_ref=35)
        via_explicit = evaluate_somatic_odds(k_alt=5, k_ref=35, c_alt=0, c_ref=30)

        assert via_defaults == via_explicit
        assert via_defaults.classification == GERMLINE_PASS
        assert via_defaults.is_somatic is False

    def test_empty_depth_table_is_germline_pass(self):
        result = evaluate_somatic_odds(k_alt=0, k_ref=0, c_alt=0, c_ref=0)

        assert result.classification == GERMLINE_PASS
        assert result.is_somatic is False
        assert result.p_value == 1.0
