"""Unit tests for pure-logic helpers in read_architecture.py.

No BAM files, pysam, matplotlib, or real data required.
"""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import pytest

from retro_miner.read_architecture import (
    LocusLayout,
    _bar_maps_to_inverted_5p_arm,
    _clip_span,
    _connector_endpoints,
    _draw_read_pair,
    _inverted_arm_arrow_label,
    _layout_from_row,
    _mei_axis_ticks,
    _pair_segments,
    _parse_support_counts,
    _safe_plot_stem,
    _support_score,
    _target_family,
    _twin_priming_summary_label,
    generate_gold_read_architecture_plots,
    plot_locus_architecture,
)


# ---------------------------------------------------------------------------
# _target_family  →  SVA | ALU | LINE1 | ERV | OTHER
# ---------------------------------------------------------------------------


def test_target_family_sva() -> None:
    assert _target_family("SVA_F#Retroposon/SVA") == "SVA"


def test_target_family_alu() -> None:
    assert _target_family("AluYa5#SINE/Alu") == "ALU"


def test_target_family_line1_keyword() -> None:
    assert _target_family("LINE1") == "LINE1"


def test_target_family_line1_slash() -> None:
    assert _target_family("L1HS#LINE/L1") == "LINE1"


def test_target_family_herv() -> None:
    assert _target_family("HERV-K") == "ERV"


def test_target_family_erv_substring() -> None:
    assert _target_family("MLT2A1#DNA/ERV") == "ERV"


def test_target_family_other() -> None:
    assert _target_family("MIR#DNA") == "OTHER"


def test_target_family_empty() -> None:
    assert _target_family("") == "OTHER"


def test_target_family_sva_wins_over_alu_prefix() -> None:
    # SVA is checked first; a string containing both SVA and ALU -> SVA
    assert _target_family("SVA_ALU_HYBRID") == "SVA"


# ---------------------------------------------------------------------------
# _parse_support_counts  →  dict with SR_L, SR_R, DPE_L, DPE_R, MEI_MAPPED,
#                           polyA_MAPPED, VNTR_MAPPED
# ---------------------------------------------------------------------------


_FULL_SUPPORT = (
    "SR_L=5,SR_R=3,DPE_L=2,DPE_R=1,MEI_MAPPED=10,polyA_MAPPED=4,VNTR_MAPPED=7"
)


def test_parse_support_counts_full() -> None:
    d = _parse_support_counts(_FULL_SUPPORT)
    assert d == {
        "SR_L": 5,
        "SR_R": 3,
        "DPE_L": 2,
        "DPE_R": 1,
        "MEI_MAPPED": 10,
        "polyA_MAPPED": 4,
        "VNTR_MAPPED": 7,
    }


def test_parse_support_counts_empty_string() -> None:
    d = _parse_support_counts("")
    assert all(v == 0 for v in d.values())


def test_parse_support_counts_none() -> None:
    d = _parse_support_counts(None)
    assert all(v == 0 for v in d.values())


def test_parse_support_counts_nan() -> None:
    d = _parse_support_counts(float("nan"))
    assert all(v == 0 for v in d.values())


def test_parse_support_counts_partial() -> None:
    d = _parse_support_counts("SR_L=9")
    assert d["SR_L"] == 9
    assert d["MEI_MAPPED"] == 0


# ---------------------------------------------------------------------------
# _support_score  →  (MEI_MAPPED, polyA_MAPPED, VNTR_MAPPED, SR+DPE flank)
# ---------------------------------------------------------------------------


def test_support_score_full() -> None:
    mei, polya, vntr, flank = _support_score(_FULL_SUPPORT)
    assert mei == 10
    assert polya == 4
    assert vntr == 7
    assert flank == 5 + 3 + 2 + 1  # SR_L + SR_R + DPE_L + DPE_R


def test_support_score_empty() -> None:
    assert _support_score("") == (0, 0, 0, 0)


# ---------------------------------------------------------------------------
# _clip_span  →  clamp x0/x1 to [x_min, x_max], swap if inverted, enforce
#               min width of 1.0
# ---------------------------------------------------------------------------


def test_clip_span_normal() -> None:
    x0, x1 = _clip_span(10.0, 20.0, 0.0, 100.0)
    assert x0 == pytest.approx(10.0)
    assert x1 == pytest.approx(20.0)


def test_clip_span_clamped() -> None:
    x0, x1 = _clip_span(-5.0, 150.0, 0.0, 100.0)
    assert x0 == pytest.approx(0.0)
    assert x1 == pytest.approx(100.0)


def test_clip_span_inverted_swaps() -> None:
    x0, x1 = _clip_span(20.0, 10.0, 0.0, 100.0)
    assert x0 <= x1
    assert x0 == pytest.approx(10.0) and x1 == pytest.approx(20.0)


def test_clip_span_zero_width_gets_min_width() -> None:
    x0, x1 = _clip_span(50.0, 50.0, 0.0, 100.0)
    assert x1 - x0 == pytest.approx(1.0)


def test_clip_span_small_width_gets_padded() -> None:
    # Width 0.1 < 1.0 → padded to 1.0
    x0, x1 = _clip_span(50.0, 50.1, 0.0, 100.0)
    assert x1 - x0 == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# _connector_endpoints  →  nearest edges of two bars
# ---------------------------------------------------------------------------


def test_connector_a_left_of_b() -> None:
    # a ends at 5, b starts at 10 → return (5, 10)
    left, right = _connector_endpoints((0.0, 5.0), (10.0, 15.0))
    assert left == pytest.approx(5.0)
    assert right == pytest.approx(10.0)


def test_connector_b_left_of_a() -> None:
    # b ends at 5, a starts at 10 → return (10, 5)
    edge_a, edge_b = _connector_endpoints((10.0, 15.0), (0.0, 5.0))
    assert edge_a == pytest.approx(10.0)
    assert edge_b == pytest.approx(5.0)


def test_connector_overlapping_returns_midpoints() -> None:
    # Overlapping → midpoints of each span
    ma, mb = _connector_endpoints((0.0, 10.0), (5.0, 15.0))
    assert ma == pytest.approx(5.0)   # midpoint of (0,10)
    assert mb == pytest.approx(10.0)  # midpoint of (5,15)


# ---------------------------------------------------------------------------
# _safe_plot_stem  →  filesystem-safe filename stem
# ---------------------------------------------------------------------------


def test_safe_plot_stem_no_rank() -> None:
    result = _safe_plot_stem("chr22", 100, 200, sample="disease")
    assert result == "read_arch_disease_chr22_100_200"


def test_safe_plot_stem_with_rank() -> None:
    result = _safe_plot_stem("chr22", 100, 200, sample="disease", rank=3)
    assert result == "rank003_read_arch_disease_chr22_100_200"


def test_safe_plot_stem_sanitises_special_chars() -> None:
    # Spaces and ! in chrom should become underscores
    result = _safe_plot_stem("chr 1!!", 0, 50, sample="ctrl")
    assert " " not in result
    assert "!" not in result
    assert "chr_1" in result


# ---------------------------------------------------------------------------
# Twin-priming 5' inversion layout, axis ticks, and plot rendering
# ---------------------------------------------------------------------------


def test_twin_priming_layout_and_summary_forward_and_reverse() -> None:
    row_fwd = pd.Series(
        {
            "chrom": "chr22",
            "window_start": 1000,
            "window_end": 1010,
            "consensus_insertion_breakpoint_pos": 1005,
            "consensus_mei_family": "LINE1",
            "consensus_mei_subfamily": "L1HS_full#LINE/L1",
            "consensus_insertion_orientation": "+",
            "consensus_insertion_mei_5p_coord_full": 100,
            "consensus_insertion_mei_3p_coord_full": 1100,
            "consensus_insertion_mei_span_full": 1001,
            "twin_priming_inv": True,
            "inv_5p_start": 100,
            "inv_5p_end": 400,
            "noninv_5p_start": 450,
            "noninv_3p_end": 1100,
            "twin_priming_confidence_score": 0.846,
            "twin_priming_confidence_tier": "high",
        }
    )
    layout_fwd = _layout_from_row(row_fwd, sample="disease", detail=None, flank_bp=200)
    assert layout_fwd.has_twin_priming_segments is True
    assert layout_fwd.twin_priming_junction_gap_bp == 49
    assert layout_fwd.inv_5p_span_x == pytest.approx((200.0, 500.0))
    assert layout_fwd.noninv_3p_span_x == pytest.approx((550.0, 1200.0))
    assert layout_fwd.twin_priming_junction_span_x == pytest.approx((500.0, 550.0))
    summary = _twin_priming_summary_label(layout_fwd)
    assert "twin_priming_5p_inv" in summary
    assert "inv_5p=100-400" in summary
    assert "noninv_3p=450-1100" in summary
    assert "junction_gap=49bp" in summary
    assert "conf=0.846(high)" in summary
    ticks, labels = _mei_axis_ticks(layout_fwd)
    assert labels == ["100", "400", "450", "1100"]
    assert ticks == pytest.approx([200.0, 500.0, 550.0, 1201.0])

    # Reverse-oriented insertion: 3' end is on the left, 5' inverted segment is on the right.
    row_rev = row_fwd.copy()
    row_rev["consensus_insertion_orientation"] = "-"
    layout_rev = _layout_from_row(row_rev, sample="disease", detail=None, flank_bp=200)
    assert layout_rev.reverse_oriented is True
    assert layout_rev.inv_5p_span_x == pytest.approx((900.0, 1200.0))
    assert layout_rev.noninv_3p_span_x == pytest.approx((200.0, 850.0))
    ticks_rev, labels_rev = _mei_axis_ticks(layout_rev)
    assert labels_rev == ["1100", "450", "400", "100"]


def test_plot_and_batch_index_include_twin_priming_annotations(tmp_path: Path) -> None:
    gold_df = pd.DataFrame(
        [
            {
                "analysis_stage_tier": "gold",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "discovery_window_start": 1000,
                "discovery_window_end": 1010,
                "consensus_insertion_breakpoint_pos": 1005,
                "sample_status_label": "disease_only",
                "disease_supporting_reads": "SR_L=4,SR_R=5,DPE_L=0,DPE_R=0,MEI_MAPPED=9,polyA_MAPPED=2",
                "control_supporting_reads": "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0",
                "consensus_mei_family": "LINE1",
                "consensus_mei_subfamily": "L1HS_full#LINE/L1",
                "consensus_insertion_orientation": "+",
                "consensus_poly_at_min_bp": 25,
                "consensus_insertion_mei_5p_coord_full": 100,
                "consensus_insertion_mei_3p_coord_full": 1100,
                "consensus_insertion_mei_span_full": 1001,
                "twin_priming_inv": True,
                "inv_5p_start": 100,
                "inv_5p_end": 400,
                "noninv_5p_start": 450,
                "noninv_3p_end": 1100,
                "twin_priming_confidence_score": 0.846,
                "twin_priming_confidence_tier": "high",
                "twin_priming_audit_trail": "sample:disease|axis:L1HS_full|inv_arm:L(-)[100-400|n:4]|direct_arm:R(+)[450-1100|n:5]|junction_gap_bp:49|conf:0.846(high)",
            }
        ]
    )
    detail_df = pd.DataFrame(
        [
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "sr_inv_1",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "R",
                "genomic_pos": 1005,
                "mei_start": 100,
                "mei_end": 220,
                "mei_target": "L1HS_full#LINE/L1",
                "mei_strand": "-",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            },
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "sr_dir_1",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L",
                "genomic_pos": 1005,
                "mei_start": 450,
                "mei_end": 600,
                "mei_target": "L1HS_full#LINE/L1",
                "mei_strand": "+",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            },
        ]
    )
    out_dir = tmp_path / "read_architecture"
    index = generate_gold_read_architecture_plots(
        gold_df,
        supporting_reads_detail=detail_df,
        out_dir=out_dir,
    )
    assert len(index) == 1
    assert index.loc[0, "status"] == "ok"
    assert bool(index.loc[0, "twin_priming_inv"]) is True
    assert int(index.loc[0, "inv_5p_start"]) == 100
    assert int(index.loc[0, "inv_5p_end"]) == 400
    assert int(index.loc[0, "noninv_5p_start"]) == 450
    assert int(index.loc[0, "noninv_3p_end"]) == 1100
    assert float(index.loc[0, "twin_priming_confidence_score"]) == pytest.approx(0.846)
    assert index.loc[0, "twin_priming_confidence_tier"] == "high"
    png_path = Path(str(index.loc[0, "png"]))
    assert png_path.exists() and png_path.stat().st_size > 0


def test_inverted_5p_arm_read_bars_get_hatching_and_strand_arrows() -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    layout = LocusLayout(
        chrom="chr22",
        window_start=1000,
        window_end=1010,
        breakpoint=1005,
        breakpoint_left=1000,
        breakpoint_right=1010,
        mei_5p=100,
        mei_3p=1100,
        flank_bp=200,
        mei_span_bp=1001,
        orientation="+",
        polya_zone_bp=25,
        twin_priming_inv=True,
        inv_5p_start=100,
        inv_5p_end=400,
        noninv_5p_start=450,
        noninv_3p_end=1100,
        twin_priming_confidence_score=0.846,
        twin_priming_confidence_tier="high",
    )
    assert _bar_maps_to_inverted_5p_arm(120, 260, "-", layout) is True
    assert _bar_maps_to_inverted_5p_arm(120, 260, "+", layout) is False
    assert _bar_maps_to_inverted_5p_arm(500, 700, "-", layout) is False
    assert _inverted_arm_arrow_label("-", layout) == "←-"

    detail_df = pd.DataFrame(
        [
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "sr_inv_1",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "R",
                "genomic_pos": 1005,
                "mei_start": 120,
                "mei_end": 240,
                "mei_target": "L1HS_full#LINE/L1",
                "mei_strand": "-",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            },
            {
                "sample": "disease",
                "evidence_type": "DPE",
                "read_name": "dpe_inv_1",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L",
                "genomic_pos": 900,
                "mate_chrom": "chr22",
                "mate_genomic_pos": 1005,
                "mei_start": 0,
                "mei_end": 0,
                "mate_mei_start": 150,
                "mate_mei_end": 300,
                "mate_mei_target": "L1HS_full#LINE/L1",
                "mate_mei_strand": "-",
                "mei_hit": False,
                "mate_mei_hit": True,
                "polya_rescue": False,
            },
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "sr_dir_1",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L",
                "genomic_pos": 1005,
                "mei_start": 500,
                "mei_end": 650,
                "mei_target": "L1HS_full#LINE/L1",
                "mei_strand": "+",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            },
        ]
    )
    pairs, stats = _pair_segments(detail_df, layout, max_pairs=50)
    assert stats["inv_5p_plotted"] == 2
    by_name = {str(p["read_name"]): p for p in pairs}
    assert by_name["sr_inv_1"]["inv_5p_arm"] is True
    assert by_name["dpe_inv_1"]["inv_5p_arm"] is True
    assert by_name["sr_dir_1"]["inv_5p_arm"] is False

    fig, ax = plt.subplots()
    try:
        for yi, p in enumerate(pairs):
            _draw_read_pair(ax, yi=float(yi), pair=p, layout=layout, read_width_bp=150.0)
        hatched_rects = [
            patch
            for patch in ax.patches
            if isinstance(patch, Rectangle) and patch.get_hatch() == "///"
        ]
        assert len(hatched_rects) == 2
        arrow_texts = [t.get_text() for t in ax.texts if "←-" in t.get_text()]
        assert len(arrow_texts) == 2
    finally:
        plt.close(fig)


def test_pair_sampling_prioritizes_inverted_sr_and_dpe_under_cap() -> None:
    layout = LocusLayout(
        chrom="chr22",
        window_start=1000,
        window_end=1010,
        breakpoint=1005,
        breakpoint_left=1000,
        breakpoint_right=1010,
        mei_5p=100,
        mei_3p=1100,
        flank_bp=200,
        mei_span_bp=1001,
        orientation="+",
        twin_priming_inv=True,
        inv_5p_start=100,
        inv_5p_end=400,
        noninv_5p_start=450,
        noninv_3p_end=1100,
    )
    rows = []
    for i in range(24):
        inv = i in {3, 17}
        rows.append(
            {
                "evidence_type": "SR",
                "read_name": f"sr_{i}",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "R" if i % 2 == 0 else "L",
                "genomic_pos": 1005,
                "mei_start": 120 if inv else 500,
                "mei_end": 240 if inv else 650,
                "mei_strand": "-" if inv else "+",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            }
        )
    for i in range(20):
        inv = i in {6, 18}
        rows.append(
            {
                "evidence_type": "DPE",
                "read_name": f"dpe_{i}",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L" if i % 2 == 0 else "R",
                "genomic_pos": 900 if i % 2 == 0 else 1100,
                "mate_chrom": "chr22",
                "mate_genomic_pos": 1005,
                "mei_start": 0,
                "mei_end": 0,
                "mate_mei_start": 150 if inv else 500,
                "mate_mei_end": 300 if inv else 650,
                "mate_mei_strand": "-" if inv else "+",
                "mei_hit": False,
                "mate_mei_hit": True,
                "polya_rescue": False,
            }
        )

    pairs, stats = _pair_segments(
        pd.DataFrame(rows),
        layout,
        max_pairs=8,
        rng=random.Random(2),
    )

    names = {str(pair["read_name"]) for pair in pairs}
    assert len(pairs) == 8
    assert stats["pairs_before_cap"] == 44
    assert stats["pairs_capped"] == 36
    assert stats["inv_5p_plotted"] == 4
    assert {"sr_3", "sr_17", "dpe_6", "dpe_18"} <= names
    assert sum(pair["evidence_type"] == "SR" for pair in pairs) >= 4
    assert sum(pair["evidence_type"] == "DPE" for pair in pairs) >= 2


def test_inverted_pairs_over_cap_are_sampled_with_sr_dpe_coverage() -> None:
    layout = LocusLayout(
        chrom="chr22",
        window_start=1000,
        window_end=1010,
        breakpoint=1005,
        breakpoint_left=1000,
        breakpoint_right=1010,
        mei_5p=100,
        mei_3p=1100,
        flank_bp=200,
        mei_span_bp=1001,
        orientation="+",
        twin_priming_inv=True,
        inv_5p_start=100,
        inv_5p_end=400,
        noninv_5p_start=450,
        noninv_3p_end=1100,
    )
    rows = []
    for i in range(12):
        rows.append(
            {
                "evidence_type": "SR",
                "read_name": f"sr_{i}",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L" if i % 2 else "R",
                "genomic_pos": 1005,
                "mei_start": 120,
                "mei_end": 240,
                "mei_strand": "-",
                "mei_hit": True,
                "mate_mei_hit": False,
                "polya_rescue": False,
                "clip_len": 60,
            }
        )
    for i in range(8):
        rows.append(
            {
                "evidence_type": "DPE",
                "read_name": f"dpe_{i}",
                "chrom": "chr22",
                "window_start": 1000,
                "window_end": 1010,
                "anchor_side": "L" if i % 2 else "R",
                "genomic_pos": 900 if i % 2 else 1100,
                "mate_chrom": "chr22",
                "mate_genomic_pos": 1005,
                "mei_start": 0,
                "mei_end": 0,
                "mate_mei_start": 150,
                "mate_mei_end": 300,
                "mate_mei_strand": "-",
                "mei_hit": False,
                "mate_mei_hit": True,
                "polya_rescue": False,
            }
        )

    pairs, stats = _pair_segments(pd.DataFrame(rows), layout, max_pairs=6)

    assert len(pairs) == 6
    assert stats["inv_5p_plotted"] == 6
    assert sum(pair["evidence_type"] == "SR" for pair in pairs) >= 1
    assert sum(pair["evidence_type"] == "DPE" for pair in pairs) >= 1
    assert all(pair["inv_5p_arm"] for pair in pairs)
