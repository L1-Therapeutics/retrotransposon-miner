"""Family-consistent panel→full MEI coordinate projection.

Guards the regression where Alu tip intervals were unioned with L1 body
projections (or used alone when L1 fragments were unmapped), producing bogus
“twin-truncation near L1 5′” spans on an Alu-length axis.
"""

from __future__ import annotations

import pandas as pd
import pytest

from retro_miner.mei_support import (
    FragmentToFullMap,
    _build_gold_review_table,
    _full_axis_union_from_panel_fragments,
    _infer_twin_priming_inversion_from_read_evidence,
    _project_panel_coords_to_full,
)
from retro_miner.vcf_export import export_vcf_from_tsv


def _entry(
    frag: str,
    full: str,
    *,
    frag_start: int,
    frag_end: int,
    full_start: int,
    full_end: int,
    family: str,
    frag_len: int = 0,
    full_len: int = 0,
) -> FragmentToFullMap:
    return FragmentToFullMap(
        fragment_name=frag,
        fragment_length=frag_len or frag_end,
        full_name=full,
        full_length=full_len or full_end,
        fragment_aln_start=frag_start,
        fragment_aln_end=frag_end,
        full_aln_start=full_start,
        full_aln_end=full_end,
        strand="+",
        family=family,
    )


@pytest.fixture
def mixed_frag_map() -> dict[str, FragmentToFullMap]:
    """Map with L1 and Alu axes (as in real mei_fragment_to_full_coords.tsv)."""
    entries = [
        _entry(
            "L1MC_orf2",
            "L1HS_full",
            frag_start=1,
            frag_end=1000,
            full_start=3000,
            full_end=3999,
            family="LINE1",
            full_len=6049,
        ),
        _entry(
            "AluSc5",
            "AluSc5_full",
            frag_start=1,
            frag_end=312,
            full_start=1,
            full_end=312,
            family="ALU",
            full_len=312,
        ),
        _entry(
            "L1MCa_5end",
            "L1MA1_full",
            frag_start=791,
            frag_end=2434,
            full_start=96,
            full_end=1691,
            family="LINE1",
            full_len=5544,
        ),
    ]
    out: dict[str, FragmentToFullMap] = {}
    for entry in entries:
        out[entry.fragment_name] = entry
    return out


class TestTwinPrimingInversionReadEvidence:
    def test_gold_review_builder_exports_detected_segment_evidence(self, tmp_path):
        row = {
            "chrom": "chr22",
            "window_start": 100000,
            "window_end": 100000,
            "insertion_breakpoint_pos": 100000,
            "consensus_mei_family": "LINE1",
            "mei_family": "LINE1",
            "consensus_insertion_orientation": "+",
            "disease_L_mei_strand": "-",
            "disease_L_mei_start": 100,
            "disease_L_mei_end": 500,
            "disease_L_mei_supported_reads": 4,
            "disease_L_mei_subfamily": "L1HS_full#LINE/L1",
            "disease_R_mei_strand": "+",
            "disease_R_mei_start": 501,
            "disease_R_mei_end": 6018,
            "disease_R_mei_supported_reads": 5,
            "disease_R_mei_subfamily": "L1HS_full#LINE/L1",
            "disease_discordant_mei_two_sided_support": False,
            "control_discordant_mei_two_sided_support": False,
            "coherence_score": 0.8,
            "insertion_model_score": 0.8,
            "read_support_heuristic_score": 1.0,
            "mei_score_enrichment_ratio": 1.0,
        }
        detail = pd.DataFrame(
            [
                {
                    "sample": "disease",
                    "evidence_type": "SR",
                    "read_name": "read_inv_1",
                    "chrom": "chr22",
                    "window_start": 100000,
                    "window_end": 100000,
                    "anchor_side": "L",
                    "mei_start": 100,
                    "mei_end": 500,
                    "mei_target": "L1HS_full#LINE/L1",
                    "mei_strand": "-",
                    "mei_hit": True,
                },
                {
                    "sample": "disease",
                    "evidence_type": "SR",
                    "read_name": "read_dir_1",
                    "chrom": "chr22",
                    "window_start": 100000,
                    "window_end": 100000,
                    "anchor_side": "R",
                    "mei_start": 501,
                    "mei_end": 6018,
                    "mei_target": "L1HS_full#LINE/L1",
                    "mei_strand": "+",
                    "mei_hit": True,
                },
            ]
        )
        review = _build_gold_review_table(
            pd.DataFrame([row]),
            empirical_stage=False,
            supporting_reads_detail=detail,
        )
        assert review.loc[0, "twin_priming_inv"]
        assert tuple(review.loc[0, ["inv_5p_start", "inv_5p_end", "noninv_5p_start", "noninv_3p_end"]]) == (
            100,
            500,
            501,
            6018,
        )
        assert float(review.loc[0, "twin_priming_confidence_score"]) >= 0.75
        assert int(review.loc[0, "twin_priming_junction_gap_bp"]) == 0
        assert review.loc[0, "twin_priming_confidence_tier"] == "high"
        audit = str(review.loc[0, "twin_priming_audit_trail"])
        assert "inv_arm:L(-)[100-500|n:4]" in audit
        assert "direct_arm:R(+)[501-6018|n:5]" in audit
        assert "inv_reads:read_inv_1(-:100-500)" in audit
        assert "direct_reads:read_dir_1(+:501-6018)" in audit

        # Exercise the real gold-review TSV → VCF export path using annotations
        # inferred by the builder rather than hand-authored twin-priming columns.
        gold_tsv = tmp_path / "candidate_loci.mei.gold_review.tsv"
        gold_vcf = tmp_path / "candidate_loci.mei.gold_review.vcf"
        review.to_csv(gold_tsv, sep="\t", index=False)
        assert export_vcf_from_tsv(gold_tsv, gold_vcf) == 1
        record_lines = [
            line for line in gold_vcf.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        assert len(record_lines) == 1
        info = dict(
            part.split("=", 1) if "=" in part else (part, True)
            for part in record_lines[0].split("\t")[7].split(";")
        )
        assert info["TWIN_PRIMING_INV"] is True
        assert info["INV_5P_START"] == "100"
        assert info["INV_5P_END"] == "500"
        assert info["NONINV_5P_START"] == "501"
        assert info["NONINV_3P_END"] == "6018"
        assert info["TWIN_PRIMING_GAP_BP"] == "0"
        assert float(info["TWIN_PRIMING_CONF"]) == pytest.approx(0.846, abs=1e-3)
        assert info["TWIN_PRIMING_TIER"] == "high"
        assert "inv_reads:read_inv_1(-:100-500)" in info["TWIN_PRIMING_AUDIT"]
        assert "direct_reads:read_dir_1(+:501-6018)" in info["TWIN_PRIMING_AUDIT"]

    def test_opposite_strand_split_arms_emit_inversion_segments(self):
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_insertion_orientation": "+",
                "disease_L_mei_strand": "-",
                "disease_L_mei_start": 100,
                "disease_L_mei_end": 500,
                "disease_L_mei_supported_reads": 4,
                "disease_L_mei_subfamily": "L1HS_full#LINE/L1",
                "disease_R_mei_strand": "+",
                "disease_R_mei_start": 501,
                "disease_R_mei_end": 6018,
                "disease_R_mei_supported_reads": 5,
                "disease_R_mei_subfamily": "L1HS_full#LINE/L1",
            }
        )
        result = _infer_twin_priming_inversion_from_read_evidence(row)
        assert result["twin_priming_inv"] is True
        assert (
            result["inv_5p_start"],
            result["inv_5p_end"],
            result["noninv_5p_start"],
            result["noninv_3p_end"],
        ) == (100, 500, 501, 6018)
        assert result["twin_priming_confidence_score"] == pytest.approx(0.846, abs=1e-3)
        assert result["twin_priming_confidence_tier"] == "high"
        assert "sample:disease|axis:L1HS_full|inv_arm:L(-)[100-500|n:4]|direct_arm:R(+)[501-6018|n:5]|junction_gap_bp:0|conf:0.846(high)" == result["twin_priming_audit_trail"]

    def test_confidence_tiers_and_read_audit_trail_scale_with_evidence(self):
        weak_row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_insertion_orientation": "+",
                "disease_L_mei_strand": "-",
                "disease_L_mei_start": 100,
                "disease_L_mei_end": 120,
                "disease_L_mei_supported_reads": 2,
                "disease_L_mei_subfamily": "L1HS_full#LINE/L1",
                "disease_L_mei_subfamily_purity": 0.5,
                "disease_L_mei_breakpoint_mode_fraction": 0.5,
                "disease_R_mei_strand": "+",
                "disease_R_mei_start": 900,
                "disease_R_mei_end": 920,
                "disease_R_mei_supported_reads": 2,
                "disease_R_mei_subfamily": "L1HS_full#LINE/L1",
                "disease_R_mei_subfamily_purity": 0.5,
                "disease_R_mei_breakpoint_mode_fraction": 0.5,
            }
        )
        weak_res = _infer_twin_priming_inversion_from_read_evidence(weak_row)
        assert weak_res["twin_priming_inv"] is True
        assert 0.0 < float(weak_res["twin_priming_confidence_score"]) < 0.50
        assert weak_res["twin_priming_confidence_tier"] == "low"
        assert "junction_gap_bp:779" in str(weak_res["twin_priming_audit_trail"])

    def test_no_call_for_same_strand_or_non_mei_family(self):
        row = pd.Series(
            {
                "consensus_mei_family": "ALU",
                "consensus_insertion_orientation": "+",
                "disease_L_mei_strand": "-",
                "disease_L_mei_start": 100,
                "disease_L_mei_end": 500,
                "disease_L_mei_supported_reads": 4,
                "disease_L_mei_subfamily": "L1HS_full#LINE/L1",
                "disease_R_mei_strand": "+",
                "disease_R_mei_start": 501,
                "disease_R_mei_end": 6018,
                "disease_R_mei_supported_reads": 5,
                "disease_R_mei_subfamily": "L1HS_full#LINE/L1",
            }
        )
        result = _infer_twin_priming_inversion_from_read_evidence(row)
        assert result["twin_priming_inv"] is False
        assert result["inv_5p_start"] == -1

    def test_panel_fragment_arms_project_before_orientation_comparison(self, mixed_frag_map):
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_insertion_orientation": "+",
                "disease_L_mei_strand": "-",
                "disease_L_mei_start": 10,
                "disease_L_mei_end": 210,
                "disease_L_mei_supported_reads": 3,
                "disease_L_mei_subfamily": "L1MC_orf2",
                "disease_R_mei_strand": "+",
                "disease_R_mei_start": 300,
                "disease_R_mei_end": 600,
                "disease_R_mei_supported_reads": 4,
                "disease_R_mei_subfamily": "L1MC_orf2",
            }
        )
        result = _infer_twin_priming_inversion_from_read_evidence(row, fragment_to_full_map=mixed_frag_map)
        assert result["twin_priming_inv"] is True
        assert (result["inv_5p_start"], result["inv_5p_end"]) == (3009, 3209)
        assert (result["noninv_5p_start"], result["noninv_3p_end"]) == (3299, 3599)

    def test_antisense_poly_a_orientation_selects_the_3p_direct_arm(self):
        row = pd.Series(
            {
                "consensus_mei_family": "SVA",
                "consensus_insertion_orientation": "-",
                "disease_L_mei_strand": "+",
                "disease_L_mei_start": 1,
                "disease_L_mei_end": 430,
                "disease_L_mei_supported_reads": 5,
                "disease_L_mei_subfamily": "SVA_F_full#Retroposon/SVA",
                "disease_R_mei_strand": "-",
                "disease_R_mei_start": 431,
                "disease_R_mei_end": 1366,
                "disease_R_mei_supported_reads": 6,
                "disease_R_mei_subfamily": "SVA_F_full#Retroposon/SVA",
            }
        )
        result = _infer_twin_priming_inversion_from_read_evidence(row)
        assert result["twin_priming_inv"] is True
        assert (result["inv_5p_start"], result["inv_5p_end"]) == (1, 430)
        assert (result["noninv_5p_start"], result["noninv_3p_end"]) == (431, 1366)


class TestProjectPanelCoordsToFull:
    def test_linear_plus_strand(self, mixed_frag_map):
        got = _project_panel_coords_to_full(100, 200, "L1MC_orf2", mixed_frag_map)
        assert got == (3099, 3199, "L1HS_full")


class TestFullAxisUnionFamilySafe:
    def test_rejects_alu_when_call_is_line1(self, mixed_frag_map):
        """L1MC_orf2 body + Alu tip must not smash into 157–3971-style unions."""
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "disease_L_mei_subfamily": "L1MC_orf2",
                "disease_L_mei_start": 120,
                "disease_L_mei_end": 971,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 157,
                "disease_R_mei_end": 186,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union is not None
        assert union == (3119, 3970)

    def test_unmapped_l1_does_not_fall_back_to_alu_axis(self, mixed_frag_map):
        """Empty L1 map must not project AluSc5 onto a fake L1 full span."""
        # Drop L1MCa mapping so only Alu can project.
        frag_map = {
            k: v for k, v in mixed_frag_map.items() if not k.startswith("L1MCa")
        }
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_start": 136,
                "disease_L_mei_end": 1009,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 136,
                "disease_R_mei_end": 309,
            }
        )
        assert _full_axis_union_from_panel_fragments(row, frag_map) is None

    def test_mapped_l1_5end_uses_l1_axis_only(self, mixed_frag_map):
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "disease_L_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_start": 800,
                "disease_L_mei_end": 1200,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 10,
                "disease_R_mei_end": 50,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union is not None
        # 800–1200 on L1MCa_5end → offset from frag 791 → full 96+(800-791) .. 96+(1200-791)
        assert union == (105, 505)

    def test_alu_call_keeps_alu_axis(self, mixed_frag_map):
        row = pd.Series(
            {
                "consensus_mei_family": "ALU",
                "disease_L_mei_subfamily": "AluSc5",
                "disease_L_mei_start": 10,
                "disease_L_mei_end": 100,
                "disease_R_mei_subfamily": "L1MC_orf2",
                "disease_R_mei_start": 120,
                "disease_R_mei_end": 200,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union == (10, 100)

    def test_prefers_consensus_subfamily_full_axis(self):
        """L1HS_5end + L1HS_3end must union on L1HS_full, not a sibling L1PA* axis."""
        frag_map = {
            "L1HS_5end": _entry(
                "L1HS_5end",
                "L1PA2_full",
                frag_start=1,
                frag_end=2136,
                full_start=1,
                full_end=2135,
                family="LINE1",
                full_len=6045,
            ),
            "L1HS_3end": _entry(
                "L1HS_3end",
                "L1HS_full",
                frag_start=1,
                frag_end=902,
                full_start=5130,
                full_end=6031,
                family="LINE1",
                full_len=6049,
            ),
        }
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_mei_subfamily": "L1HS_5end#LINE/L1",
                "disease_L_mei_subfamily": "L1HS_5end",
                "disease_L_mei_start": 1,
                "disease_L_mei_end": 840,
                "disease_R_mei_subfamily": "L1HS_3end",
                "disease_R_mei_start": 483,
                "disease_R_mei_end": 889,
            }
        )
        # Only the L1HS_full axis matches the consensus subfamily preference.
        assert _full_axis_union_from_panel_fragments(row, frag_map) == (5612, 6018)
