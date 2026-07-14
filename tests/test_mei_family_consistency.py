"""Tests for MEI family/subfamily consensus consistency.

Covers the fix for the case where consensus_mei_family could contradict
consensus_mei_subfamily (see issue #11): the two columns were selected by
independent procedures with no cross-check.
"""

import pandas as pd
import pytest

from retro_miner.mei_support import (
    _add_consolidated_event_fields,
    _choose_event_family,
    _choose_event_subfamily,
    _normalize_mei_family_token,
)


class TestNormalizeMeiFamilyToken:
    @pytest.mark.parametrize(
        "token,expected",
        [
            ("AluYb8#SINE/Alu", "ALU"),
            ("AluYa5#SINE/Alu", "ALU"),
            ("FAM#SINE/Alu", "ALU"),
            ("SVA_D#Retroposon/SVA", "SVA"),
            ("SVA_E#Retroposon/SVA", "SVA"),
            ("L1HS_5end#LINE/L1", "LINE1"),
            ("L1MB8_3end#LINE/L1", "LINE1"),
        ],
    )
    def test_known_subfamilies_map_to_expected_family(self, token, expected):
        assert _normalize_mei_family_token(token) == expected

    def test_matching_is_case_insensitive(self):
        assert _normalize_mei_family_token("aluyb8#sine/alu") == "ALU"

    @pytest.mark.parametrize("token", ["", None, "NOT_A_RETROTRANSPOSON"])
    def test_unrecognized_tokens_return_empty(self, token):
        assert _normalize_mei_family_token(token) == ""


class TestFamilySubfamilyConsistency:
    """The regression this fix addresses.

    The README example table contained a locus at chr22:23853715 annotated
    family=ALU with subfamily=SVA_E#Retroposon/SVA. The subfamily is chosen by
    weighted read support while the family was chosen by scanning candidate
    columns in a fixed order, so the two could land on different classes.
    """

    def _row(self, **overrides):
        base = {
            "disease_L_mei_family": "",
            "disease_R_mei_family": "",
            "disease_discordant_mei_family": "",
            "control_L_mei_family": "",
            "control_R_mei_family": "",
            "control_discordant_mei_family": "",
            "disease_L_mei_subfamily": "",
            "disease_R_mei_subfamily": "",
            "disease_discordant_mei_subfamily": "",
            "control_L_mei_subfamily": "",
            "control_R_mei_subfamily": "",
            "control_discordant_mei_subfamily": "",
            "disease_discordant_mei_supported_reads": 0,
            "control_discordant_mei_supported_reads": 0,
            "disease_left_supported_reads": 0,
            "disease_right_supported_reads": 0,
            "control_left_supported_reads": 0,
            "control_right_supported_reads": 0,
            "g1k_melt_insertion_subfamily": "",
            "g1k_melt_id": "",
            "lr_svan_subfamily": "",
            "lr_svan_id": "",
            "known_mei_polymorphism_family": "",
            "known_mei_polymorphism_subfamily": "",
        }
        base.update(overrides)
        return base

    def test_family_agrees_with_subfamily_when_sources_conflict(self):
        """A row whose family columns say ALU but whose winning subfamily is SVA
        must be reported as SVA."""
        row = self._row(
            disease_L_mei_family="ALU",
            disease_discordant_mei_subfamily="SVA_E#Retroposon/SVA",
            disease_discordant_mei_supported_reads=40,
        )
        out = _add_consolidated_event_fields(pd.DataFrame([row]))

        subfamily = out.loc[0, "mei_subfamily"]
        family = out.loc[0, "mei_family"]

        assert subfamily == "SVA_E#Retroposon/SVA"
        assert family == "SVA"
        assert family == _normalize_mei_family_token(subfamily)

    @pytest.mark.parametrize(
        "subfamily,expected_family",
        [
            ("AluYb8#SINE/Alu", "ALU"),
            ("SVA_E#Retroposon/SVA", "SVA"),
            ("L1HS_5end#LINE/L1", "LINE1"),
        ],
    )
    def test_family_is_always_derivable_from_subfamily(self, subfamily, expected_family):
        row = self._row(
            disease_discordant_mei_subfamily=subfamily,
            disease_discordant_mei_supported_reads=10,
        )
        out = _add_consolidated_event_fields(pd.DataFrame([row]))

        assert out.loc[0, "mei_family"] == expected_family

    def test_family_falls_back_to_column_scan_when_subfamily_is_empty(self):
        """If no subfamily is resolvable, the existing family precedence scan
        should still supply a value."""
        row = self._row(disease_L_mei_family="ALU")
        out = _add_consolidated_event_fields(pd.DataFrame([row]))

        assert out.loc[0, "mei_subfamily"] == ""
        assert out.loc[0, "mei_family"] == "ALU"

    def test_no_family_and_no_subfamily_yields_empty(self):
        out = _add_consolidated_event_fields(pd.DataFrame([self._row()]))

        assert out.loc[0, "mei_subfamily"] == ""
        assert out.loc[0, "mei_family"] == ""

