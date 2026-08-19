"""Unit tests for pure-logic helpers in read_architecture.py.

Pure-logic helper tests do not require BAM files, pysam, or real data.
The TestPlotFigureCleanup class patches matplotlib internals to verify
that plt.close(fig) is guaranteed even when fig.savefig raises.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from retro_miner.read_architecture import (
    _clip_span,
    _connector_endpoints,
    _parse_support_counts,
    _safe_plot_stem,
    _support_score,
    _target_family,
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
    # Spaces and ! in chrom should become underscores; a digest is embedded
    # because sanitisation changed the value.
    result = _safe_plot_stem("chr 1!!", 0, 50, sample="ctrl")
    assert " " not in result
    assert "!" not in result
    assert "chr_1" in result


def test_safe_plot_stem_sanitised_chroms_do_not_collide() -> None:
    """Distinct chrom strings that sanitise to the same form must produce distinct stems."""
    stem_slash = _safe_plot_stem("chr1/a", 10, 20, sample="disease")
    stem_under = _safe_plot_stem("chr1_a", 10, 20, sample="disease")
    stem_colon = _safe_plot_stem("chr1:a", 10, 20, sample="disease")
    # chr1_a needs no sanitisation; it must keep the digest-free format.
    assert stem_under == "read_arch_disease_chr1_a_10_20"
    # All three distinct inputs must produce distinct stems.
    assert len({stem_slash, stem_under, stem_colon}) == 3


# ---------------------------------------------------------------------------
# plot_locus_architecture — figure cleanup on exception
# ---------------------------------------------------------------------------

class TestPlotFigureCleanup:
    """plt.close(fig) is called even when fig.savefig raises an exception.

    All expensive helpers (BAM loading, layout computation, alignment) are
    patched out so the test exercises only the try/finally cleanup guarantee.
    """

    def _minimal_pair_stats(self) -> dict:
        return {
            "pairs_shown": 0,
            "pairs_before_cap": 0,
            "sr_plotted": 0,
            "dpe_plotted": 0,
            "polya_rescue_plotted": 0,
            "vntr_rescue_plotted": 0,
            "sr_skipped": 0,
            "dpe_skipped": 0,
            "detail_rows": 0,
        }

    def _minimal_layout(self) -> MagicMock:
        m = MagicMock()
        # Keep numeric attrs as real numbers to avoid arithmetic surprises.
        m.total_width = 1000
        m.polya_zone_bp = 0  # falsy → polya block skipped
        m.insertion_left_x = 100
        m.insertion_right_x = 900
        m.mei_region_start_x = 300
        m.mei_region_end_x = 700
        m.reverse_oriented = False
        m.flank_bp = 200
        m.breakpoint = 500
        m.insert_size_estimates = []  # falsy → insert_hint skipped
        return m

    def _minimal_row(self) -> pd.Series:
        return pd.Series({
            "chrom": "chr1",
            "window_start": 1000,
            "window_end": 1200,
            "discovery_window_start": 1000,
            "discovery_window_end": 1200,
            "consensus_mei_family": "ALU",
            "consensus_mei_subfamily": "AluYa5",
            "analysis_stage_tier": "gold",
            "disease_supporting_reads": 3,
        })

    def test_figure_closed_on_savefig_exception(self, tmp_path: Path) -> None:
        """plt.close(fig) is called even when fig.savefig raises RuntimeError."""
        # Minimal supporting_reads_detail file so the FileNotFoundError guard passes.
        detail_tsv = tmp_path / "detail.mei.tsv"
        detail_tsv.write_text("chrom\tpos\tread_name\n", encoding="utf-8")

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_fig.savefig.side_effect = RuntimeError("simulated disk-full")

        closed_figs: list = []
        layout = self._minimal_layout()
        pair_stats = self._minimal_pair_stats()
        row = self._minimal_row()

        patches = {
            "retro_miner.read_architecture._build_read_table_for_locus": MagicMock(
                return_value=pd.DataFrame()
            ),
            "retro_miner.read_architecture._layout_from_row": MagicMock(return_value=layout),
            "retro_miner.read_architecture._pair_segments": MagicMock(
                return_value=([], pair_stats)
            ),
            "retro_miner.read_architecture._breakpoint_interval": MagicMock(
                return_value=(1100, 1110, 1105)
            ),
            "retro_miner.read_architecture._auto_flank_bp": MagicMock(return_value=200),
            "retro_miner.read_architecture._choose_sample": MagicMock(return_value="disease"),
            "retro_miner.read_architecture._sample_status_label": MagicMock(return_value="shared"),
            "retro_miner.read_architecture._mei_axis_ticks": MagicMock(return_value=([], [])),
            "retro_miner.read_architecture.blended_transform_factory": MagicMock(
                return_value=MagicMock()
            ),
            "matplotlib.pyplot.subplots": MagicMock(return_value=(mock_fig, mock_ax)),
            "matplotlib.pyplot.close": MagicMock(side_effect=closed_figs.append),
        }

        import contextlib
        with contextlib.ExitStack() as stack:
            for target, mock_obj in patches.items():
                stack.enter_context(patch(target, mock_obj))

            with pytest.raises(RuntimeError, match="disk-full"):
                plot_locus_architecture(
                    chrom="chr1",
                    pos=1105,
                    out_png=tmp_path / "out.png",
                    # gold_tsv required even when row is provided (cache=None path).
                    gold_tsv=tmp_path / "fake_gold.tsv",
                    row=row,
                    supporting_reads_detail=detail_tsv,
                )

        assert len(closed_figs) == 1, "plt.close(fig) was not called after savefig raised"
        assert closed_figs[0] is mock_fig

    def test_figure_closed_on_successful_save(self, tmp_path: Path) -> None:
        """plt.close(fig) is also called on the normal (no-exception) code path."""
        detail_tsv = tmp_path / "detail.mei.tsv"
        detail_tsv.write_text("chrom\tpos\tread_name\n", encoding="utf-8")

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        # savefig succeeds (no side_effect).

        closed_figs: list = []
        layout = self._minimal_layout()
        pair_stats = self._minimal_pair_stats()
        row = self._minimal_row()
        out_png = tmp_path / "out.png"

        patches = {
            "retro_miner.read_architecture._build_read_table_for_locus": MagicMock(
                return_value=pd.DataFrame()
            ),
            "retro_miner.read_architecture._layout_from_row": MagicMock(return_value=layout),
            "retro_miner.read_architecture._pair_segments": MagicMock(
                return_value=([], pair_stats)
            ),
            "retro_miner.read_architecture._breakpoint_interval": MagicMock(
                return_value=(1100, 1110, 1105)
            ),
            "retro_miner.read_architecture._auto_flank_bp": MagicMock(return_value=200),
            "retro_miner.read_architecture._choose_sample": MagicMock(return_value="disease"),
            "retro_miner.read_architecture._sample_status_label": MagicMock(return_value="shared"),
            "retro_miner.read_architecture._mei_axis_ticks": MagicMock(return_value=([], [])),
            "retro_miner.read_architecture.blended_transform_factory": MagicMock(
                return_value=MagicMock()
            ),
            "matplotlib.pyplot.subplots": MagicMock(return_value=(mock_fig, mock_ax)),
            "matplotlib.pyplot.close": MagicMock(side_effect=closed_figs.append),
        }

        import contextlib
        with contextlib.ExitStack() as stack:
            for target, mock_obj in patches.items():
                stack.enter_context(patch(target, mock_obj))

            returned_png, _detail = plot_locus_architecture(
                chrom="chr1",
                pos=1105,
                out_png=out_png,
                gold_tsv=tmp_path / "fake_gold.tsv",
                row=row,
                supporting_reads_detail=detail_tsv,
            )

        assert returned_png == out_png
        assert len(closed_figs) == 1, "plt.close(fig) was not called on successful save"
        assert closed_figs[0] is mock_fig
