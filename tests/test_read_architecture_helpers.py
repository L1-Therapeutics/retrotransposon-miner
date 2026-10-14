"""Unit tests for pure-logic helpers in read_architecture.py.

No BAM files, pysam, matplotlib, or real data required.
"""

from __future__ import annotations

import pandas as pd
import pytest

from retro_miner.read_architecture import (
    _clip_span,
    _connector_endpoints,
    _enrich_detail_with_clip_lens,
    _load_split_evidence,
    _numeric_series,
    _parse_support_counts,
    _row_int,
    _safe_plot_stem,
    _support_score,
    _target_family,
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
# _row_int  →  int64 DataFrame columns must not collapse to 0
# ---------------------------------------------------------------------------


def test_row_int_reads_int64_coordinates() -> None:
    # row[col] on an int64 Series yields np.int64, which does not subclass
    # Python int; regressions here silently zeroed out every integer table
    # column (breakpoints, span estimates, window coordinates).
    row = pd.Series(
        {
            "consensus_insertion_breakpoint_pos": 12345678,
            "window_start": 100,
            "window_end": 200,
            "insertion_model_score": 0.875,
        }
    )
    assert _row_int(row, "consensus_insertion_breakpoint_pos") == 12345678
    assert _row_int(row, "window_start") == 100
    assert _row_int(row, "window_end") == 200
    assert _row_int(row, "missing_column", default=-1) == -1


def test_row_int_missing_or_na_returns_default() -> None:
    row = pd.Series({"window_start": pd.NA}, dtype="object")
    assert _row_int(row, "window_start", default=5) == 5


# ---------------------------------------------------------------------------
# _numeric_series  →  always a Series, even when the column is absent
# ---------------------------------------------------------------------------


def test_numeric_series_present_column_coerces_non_numeric() -> None:
    frame = pd.DataFrame({"x": ["1", "oops", None]})
    series = _numeric_series(frame, "x")
    assert series.iloc[0] == 1.0
    assert pd.isna(series.iloc[1])
    assert pd.isna(series.iloc[2])


def test_numeric_series_absent_column_returns_constant_series() -> None:
    # Regression: DataFrame.get(col, 0) yields a bare float for a missing
    # column, so the old pd.to_numeric(...).fillna(0) chain raised
    # AttributeError instead of degrading gracefully.
    frame = pd.DataFrame({"a": [10, 20, 30]})
    series = _numeric_series(frame, "missing")
    assert list(series) == [0.0, 0.0, 0.0]
    assert list(series.index) == [0, 1, 2]


def test_numeric_series_absent_column_honours_default() -> None:
    frame = pd.DataFrame({"a": [1, 2]})
    assert list(_numeric_series(frame, "missing", default=-1.0)) == [-1.0, -1.0]


# ---------------------------------------------------------------------------
# _enrich_detail_with_clip_lens  →  clip_len fallback without the span columns
# ---------------------------------------------------------------------------


def test_enrich_detail_with_missing_span_columns_does_not_crash() -> None:
    # Regression: when neither mei_end nor mei_start is present and there is no
    # split evidence, the fallback span calculation used to raise
    # AttributeError ('float' object has no attribute 'fillna').
    detail = pd.DataFrame({"read_name": ["r1"], "evidence_type": ["SR"]})
    out = _enrich_detail_with_clip_lens(detail, gold_tsv=None, sample="disease")
    assert "clip_len" in out.columns
    assert out["clip_len"].notna().all()


def test_enrich_detail_fallback_span_uses_mei_coordinates() -> None:
    detail = pd.DataFrame(
        {
            "read_name": ["sr1", "dpe1"],
            "evidence_type": ["SR", "DPE"],
            "mei_start": [100, 100],
            "mei_end": [150, 150],
        }
    )
    out = _enrich_detail_with_clip_lens(detail, gold_tsv=None, sample="disease")
    # SR rows take the mei_start..mei_end span; non-SR rows keep the 0 default.
    assert out.loc[out["read_name"] == "sr1", "clip_len"].iloc[0] == 51
    assert out.loc[out["read_name"] == "dpe1", "clip_len"].iloc[0] == 0


# ---------------------------------------------------------------------------
# _load_split_evidence  →  only the present split-evidence columns
# ---------------------------------------------------------------------------


def test_load_split_evidence_keeps_present_columns(tmp_path) -> None:
    path = tmp_path / "split_evidence.disease.tsv"
    path.write_text("read_name\tclip_len\tclip_side\nr1\t111\tR\n", encoding="utf-8")
    df = _load_split_evidence(path)
    assert list(df.columns) == ["read_name", "clip_len", "clip_side"]
    assert int(df["clip_len"].iloc[0]) == 111


def test_load_split_evidence_ignores_unexpected_columns(tmp_path) -> None:
    path = tmp_path / "split_evidence.control.tsv"
    path.write_text("read_name\tsecret\tpos\nr2\tx\t99\n", encoding="utf-8")
    df = _load_split_evidence(path)
    assert list(df.columns) == ["read_name", "pos"]

