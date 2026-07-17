"""Family-first MEI identity voting and DPE mate-distance gate.

Covers the README/gold failure mode where ALU majority support was split across
Alu subfamilies while a single SVA/LINE1 label won raw subfamily vote, and where
same-chr DPE mates into adjacent reference MEIs polluted identity.
"""

from __future__ import annotations

import pandas as pd
import pytest

from retro_miner.mei_support import (
    _DPE_MEI_IDENTITY_MIN_SAME_CHR_BP,
    _add_consolidated_event_fields,
    _aggregate_discordant_mei_metrics,
    _choose_event_family_and_subfamily,
    _discordant_mate_ok_for_mei_identity,
    _discordant_rows_for_mei_identity,
    _discordant_rows_for_mei_mapped_support,
    _normalize_mei_family_token,
    _top_family_then_subfamily,
)


class TestNormalizeMeiFamilyToken:
    @pytest.mark.parametrize(
        "token,expected",
        [
            ("AluYb8#SINE/Alu", "ALU"),
            ("SVA_E#Retroposon/SVA", "SVA"),
            ("L1HS_5end#LINE/L1", "LINE1"),
            ("L1M3_orf2#LINE/L1", "LINE1"),
        ],
    )
    def test_known_tokens(self, token, expected):
        assert _normalize_mei_family_token(token) == expected


class TestFamilyFirstVoting:
    def test_split_alu_subfamilies_beat_single_sva(self):
        """ALU 6+5 beats SVA 7; subfamily must be AluY, not SVA_D."""
        row = pd.Series(
            {
                "disease_discordant_mei_subfamily": "AluY#SINE/Alu",
                "disease_discordant_mei_supported_reads": 6,
                "disease_L_mei_subfamily": "AluS#SINE/Alu",
                "disease_L_mei_supported_reads": 5,
                "control_discordant_mei_subfamily": "SVA_D#Retroposon/SVA",
                "control_discordant_mei_supported_reads": 7,
                "disease_R_mei_subfamily": "",
                "disease_R_mei_supported_reads": 0,
                "control_L_mei_subfamily": "",
                "control_L_mei_supported_reads": 0,
                "control_R_mei_subfamily": "",
                "control_R_mei_supported_reads": 0,
            }
        )
        family, subfamily = _choose_event_family_and_subfamily(row)
        assert family == "ALU"
        assert subfamily == "AluY#SINE/Alu"

    def test_pools_disease_and_control_without_precedence(self):
        """Control-only Alu majority should win over disease SVA minority."""
        row = pd.Series(
            {
                "disease_discordant_mei_subfamily": "SVA_E#Retroposon/SVA",
                "disease_discordant_mei_supported_reads": 4,
                "control_discordant_mei_subfamily": "AluYb8#SINE/Alu",
                "control_discordant_mei_supported_reads": 10,
                "disease_L_mei_subfamily": "",
                "disease_L_mei_supported_reads": 0,
                "disease_R_mei_subfamily": "",
                "disease_R_mei_supported_reads": 0,
                "control_L_mei_subfamily": "",
                "control_L_mei_supported_reads": 0,
                "control_R_mei_subfamily": "",
                "control_R_mei_supported_reads": 0,
            }
        )
        family, subfamily = _choose_event_family_and_subfamily(row)
        assert family == "ALU"
        assert subfamily == "AluYb8#SINE/Alu"

    def test_sentinel_style_conflict_keeps_alu(self):
        """chr22:23853715-style: pool ALU vs SVA across samples, not winner×total."""
        row = {
            "disease_discordant_mei_family": "SVA",
            "disease_discordant_mei_subfamily": "SVA_E#Retroposon/SVA",
            "disease_discordant_mei_supported_reads": 51,
            "disease_discordant_mei_family_votes": "ALU:36,SVA:15",
            "disease_discordant_mei_subfamily_votes": "AluSc8#SINE/Alu:17,SVA_E#Retroposon/SVA:15",
            "control_discordant_mei_family": "ALU",
            "control_discordant_mei_subfamily": "AluYe5#SINE/Alu",
            "control_discordant_mei_supported_reads": 26,
            "control_discordant_mei_family_votes": "ALU:25,SVA:1",
            "control_discordant_mei_subfamily_votes": "AluYe5#SINE/Alu:10,AluY#SINE/Alu:15",
            "disease_L_mei_subfamily": "",
            "disease_L_mei_supported_reads": 0,
            "disease_R_mei_subfamily": "",
            "disease_R_mei_supported_reads": 0,
            "control_L_mei_subfamily": "",
            "control_L_mei_supported_reads": 0,
            "control_R_mei_subfamily": "",
            "control_R_mei_supported_reads": 0,
        }
        # Pooled family votes: ALU 36+25=61 vs SVA 15+1=16 -> ALU.
        out = _add_consolidated_event_fields(pd.DataFrame([row]))
        assert out.loc[0, "mei_family"] == "ALU"
        assert _normalize_mei_family_token(out.loc[0, "mei_subfamily"]) == "ALU"

    def test_pools_family_votes_not_winner_times_total(self):
        """Disease winner SVA×51 must not beat pooled ALU when vote maps exist."""
        row = pd.Series(
            {
                "disease_discordant_mei_subfamily": "SVA_E#Retroposon/SVA",
                "disease_discordant_mei_supported_reads": 51,
                "disease_discordant_mei_family_votes": "ALU:36,SVA:15",
                "disease_discordant_mei_subfamily_votes": "AluY#SINE/Alu:36,SVA_E#Retroposon/SVA:15",
                "control_discordant_mei_subfamily": "AluYe5#SINE/Alu",
                "control_discordant_mei_supported_reads": 26,
                "control_discordant_mei_family_votes": "ALU:25,SVA:1",
                "control_discordant_mei_subfamily_votes": "AluYe5#SINE/Alu:25,SVA_E#Retroposon/SVA:1",
                "disease_L_mei_subfamily": "",
                "disease_L_mei_supported_reads": 0,
                "disease_R_mei_subfamily": "",
                "disease_R_mei_supported_reads": 0,
                "control_L_mei_subfamily": "",
                "control_L_mei_supported_reads": 0,
                "control_R_mei_subfamily": "",
                "control_R_mei_supported_reads": 0,
            }
        )
        family, subfamily = _choose_event_family_and_subfamily(row)
        assert family == "ALU"
        assert _normalize_mei_family_token(subfamily) == "ALU"

    def test_control_l1_weight_does_not_override_disease_alu_majority(self):
        """chr22:28994190-style: control L1MB2 wt=7 must not beat pooled ALU."""
        row = pd.Series(
            {
                "disease_discordant_mei_subfamily": "AluYh3#SINE/Alu",
                "disease_discordant_mei_supported_reads": 5,
                "control_discordant_mei_subfamily": "L1MB2_3end#LINE/L1",
                "control_discordant_mei_supported_reads": 7,
                "disease_L_mei_subfamily": "AluYh7#SINE/Alu",
                "disease_L_mei_supported_reads": 20,
                "disease_R_mei_subfamily": "",
                "disease_R_mei_supported_reads": 0,
                "control_L_mei_subfamily": "",
                "control_L_mei_supported_reads": 0,
                "control_R_mei_subfamily": "",
                "control_R_mei_supported_reads": 0,
            }
        )
        family, subfamily = _choose_event_family_and_subfamily(row)
        assert family == "ALU"
        assert "Alu" in subfamily or "ALU" in subfamily.upper()


class TestDpeMateDistanceGate:
    def test_threshold_is_1000bp(self):
        assert _DPE_MEI_IDENTITY_MIN_SAME_CHR_BP == 1000

    def test_interchrom_ok(self):
        df = pd.DataFrame(
            [{"chrom": "chr22", "pos": 1000, "mate_chrom": "chr2", "mate_pos": 5000}]
        )
        assert bool(_discordant_mate_ok_for_mei_identity(df).iloc[0])

    def test_same_chr_within_1kb_rejected(self):
        df = pd.DataFrame(
            [{"chrom": "chr22", "pos": 40539002, "mate_chrom": "chr22", "mate_pos": 40538200}]
        )
        assert not bool(_discordant_mate_ok_for_mei_identity(df).iloc[0])

    def test_same_chr_beyond_1kb_ok(self):
        df = pd.DataFrame(
            [{"chrom": "chr22", "pos": 1000, "mate_chrom": "chr22", "mate_pos": 2500}]
        )
        assert bool(_discordant_mate_ok_for_mei_identity(df).iloc[0])

    def test_identity_rows_prefer_mate_labels_but_keep_anchor_only(self):
        """Distance-OK anchor-only hits still vote; mate labels preferred when present."""
        df = pd.DataFrame(
            [
                {
                    "chrom": "chr22",
                    "pos": 1000,
                    "mate_chrom": "chr2",
                    "mate_pos": 5000,
                    "mate_mei_hit": False,
                    "mei_hit": True,
                    "target": "AluY#SINE/Alu",
                    "family": "ALU",
                    "mei_score": 0.9,
                    "read_name": "anchor_only_alu",
                },
                {
                    "chrom": "chr22",
                    "pos": 1000,
                    "mate_chrom": "chr2",
                    "mate_pos": 6000,
                    "mate_mei_hit": True,
                    "mei_hit": True,
                    "target": "SVA_E#Retroposon/SVA",
                    "family": "SVA",
                    "mate_mei_target": "AluY#SINE/Alu",
                    "mate_mei_family": "ALU",
                    "mei_score": 0.8,
                    "read_name": "mate_overrides",
                },
            ]
        )
        out = _discordant_rows_for_mei_identity(df)
        assert len(out) == 2
        by_name = out.set_index("read_name")
        assert by_name.loc["anchor_only_alu", "family"] == "ALU"
        assert by_name.loc["mate_overrides", "family"] == "ALU"
        assert by_name.loc["mate_overrides", "target"] == "AluY#SINE/Alu"

    def test_aggregate_ignores_adjacent_ref_sva_mates(self):
        """Adjacent same-chr SVA mates must not set discordant family/subfamily."""
        rows = []
        # 10 local SVA mates into adjacent reference element
        for i in range(10):
            rows.append(
                {
                    "chrom": "chr22",
                    "window_start": 40538991,
                    "window_end": 40539012,
                    "pos": 40539000 + i,
                    "mate_chrom": "chr22",
                    "mate_pos": 40538200 + i,
                    "read_name": f"sva_local_{i}",
                    "mei_hit": True,
                    "mate_mei_hit": True,
                    "family": "SVA",
                    "target": "SVA_F#Retroposon/SVA",
                    "mate_mei_family": "SVA",
                    "mate_mei_target": "SVA_F#Retroposon/SVA",
                    "target_strand": "+",
                    "target_start": 1000,
                    "target_end": 1100,
                    "mei_score": 1.0,
                    "template_len": 500,
                }
            )
        # 3 far/interchrom ALU mates
        for i in range(3):
            rows.append(
                {
                    "chrom": "chr22",
                    "window_start": 40538991,
                    "window_end": 40539012,
                    "pos": 40539000 + i,
                    "mate_chrom": "chr10",
                    "mate_pos": 100000 + i,
                    "read_name": f"alu_far_{i}",
                    "mei_hit": True,
                    "mate_mei_hit": True,
                    "family": "ALU",
                    "target": "AluY#SINE/Alu",
                    "mate_mei_family": "ALU",
                    "mate_mei_target": "AluY#SINE/Alu",
                    "target_strand": "+",
                    "target_start": 50,
                    "target_end": 150,
                    "mei_score": 1.0,
                    "template_len": 500,
                }
            )
        agg = _aggregate_discordant_mei_metrics(pd.DataFrame(rows), "disease")
        assert len(agg) == 1
        assert agg.iloc[0]["disease_discordant_mei_family"] == "ALU"
        assert "Alu" in str(agg.iloc[0]["disease_discordant_mei_subfamily"])
        # Raw discordant MEI support still includes all mei_hits (local + far);
        # MEI_MAPPED uses the mate-geometry gate separately.
        assert int(agg.iloc[0]["disease_discordant_mei_supported_reads"]) == 13

    def test_mei_mapped_excludes_nearby_same_chr_mates(self):
        """MEI_MAPPED must not count same-chr mates within 1 kb (nearby ref MEIs)."""
        df = pd.DataFrame(
            [
                {
                    "chrom": "chr22",
                    "pos": 40539002,
                    "mate_chrom": "chr22",
                    "mate_pos": 40538200,
                    "mei_hit": True,
                    "mate_mei_hit": True,
                    "read_name": "local_ref_mei",
                },
                {
                    "chrom": "chr22",
                    "pos": 40539002,
                    "mate_chrom": "chr10",
                    "mate_pos": 100000,
                    "mei_hit": True,
                    "mate_mei_hit": True,
                    "read_name": "interchrom_ok",
                },
                {
                    "chrom": "chr22",
                    "pos": 40539002,
                    "mate_chrom": "chr22",
                    "mate_pos": 40541000,
                    "mei_hit": True,
                    "mate_mei_hit": True,
                    "read_name": "same_chr_far_ok",
                },
            ]
        )
        out = _discordant_rows_for_mei_mapped_support(df)
        assert set(out["read_name"]) == {"interchrom_ok", "same_chr_far_ok"}


class TestTopFamilyThenSubfamily:
    def test_helper_pools_family_before_label(self):
        df = pd.DataFrame(
            [
                {"chrom": "chr22", "window_start": 1, "window_end": 2, "family": "ALU", "target": "AluY", "mei_score": 6},
                {"chrom": "chr22", "window_start": 1, "window_end": 2, "family": "ALU", "target": "AluS", "mei_score": 5},
                {"chrom": "chr22", "window_start": 1, "window_end": 2, "family": "SVA", "target": "SVA_D", "mei_score": 7},
            ]
        )
        fam, sub = _top_family_then_subfamily(
            df,
            group_cols=["chrom", "window_start", "window_end"],
            family_col="family",
            subfamily_col="target",
            score_col="mei_score",
        )
        assert fam.iloc[0]["family"] == "ALU"
        assert sub.iloc[0]["target"] == "AluY"
