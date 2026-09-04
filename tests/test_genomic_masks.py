"""Unit tests for genomic mask and BED intersection logic.

Covers ``_parse_interval_parts``, ``_normalize_track_to_bed``,
``_annotate_junk_flags``, ``_write_candidate_windows_bed``, and
boundary/edge-case handling for interval coordinates.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from retro_miner.candidate_loci import (
    _annotate_junk_flags,
    _distance_to_closed_interval,
    _merge_overlapping_intervals,
    _normalize_track_to_bed,
    _parse_interval_parts,
    _write_candidate_windows_bed,
)


# ─────────────────────────────────────────────────────────────────────────────
# _parse_interval_parts
# ─────────────────────────────────────────────────────────────────────────────


class TestParseIntervalParts:
    def test_standard_bed_three_columns(self):
        result = _parse_interval_parts(["chr1", "100", "200"])
        assert result == ("chr1", 100, 200)

    def test_standard_bed_four_columns(self):
        result = _parse_interval_parts(["chr22", "5000", "6000", "name"])
        assert result == ("chr22", 5000, 6000)

    def test_ucsc_bin_format_four_columns(self):
        result = _parse_interval_parts(["100", "chr1", "200", "300"])
        assert result == ("chr1", 200, 300)

    def test_non_chrom_first_column_returns_none(self):
        assert _parse_interval_parts(["100", "200", "300"]) is None

    def test_too_few_columns_returns_none(self):
        assert _parse_interval_parts(["chr1", "100"]) is None

    def test_non_numeric_start_returns_none(self):
        assert _parse_interval_parts(["chr1", "abc", "200"]) is None

    def test_non_numeric_end_returns_none(self):
        assert _parse_interval_parts(["chr1", "100", "xyz"]) is None

    def test_end_equal_start_returns_none(self):
        assert _parse_interval_parts(["chr1", "100", "100"]) is None

    def test_end_less_than_start_returns_none(self):
        assert _parse_interval_parts(["chr1", "200", "100"]) is None

    def test_empty_list_returns_none(self):
        assert _parse_interval_parts([]) is None

    def test_single_element_returns_none(self):
        assert _parse_interval_parts(["chr1"]) is None

    def test_chrom_x_and_y(self):
        assert _parse_interval_parts(["chrX", "0", "1000"]) == ("chrX", 0, 1000)
        assert _parse_interval_parts(["chrY", "50", "150"]) == ("chrY", 50, 150)

    def test_mitochondrial(self):
        assert _parse_interval_parts(["chrM", "0", "16569"]) == ("chrM", 0, 16569)


# ─────────────────────────────────────────────────────────────────────────────
# _normalize_track_to_bed
# ─────────────────────────────────────────────────────────────────────────────


class TestNormalizeTrackToBed:
    def test_standard_bed_passthrough(self, tmp_path):
        inp = tmp_path / "input.bed"
        out = tmp_path / "output.bed"
        inp.write_text("chr1\t100\t200\nchr2\t300\t400\n")
        _normalize_track_to_bed(inp, out)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2
        assert lines[0] == "chr1\t100\t200"

    def test_comment_lines_skipped(self, tmp_path):
        inp = tmp_path / "input.bed"
        out = tmp_path / "output.bed"
        inp.write_text("# header\nchr1\t100\t200\n# another\nchr2\t300\t400\n")
        _normalize_track_to_bed(inp, out)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2

    def test_empty_file(self, tmp_path):
        inp = tmp_path / "input.bed"
        out = tmp_path / "output.bed"
        inp.write_text("")
        _normalize_track_to_bed(inp, out)
        assert out.read_text() == ""

    def test_non_bed_lines_skipped(self, tmp_path):
        inp = tmp_path / "input.bed"
        out = tmp_path / "output.bed"
        inp.write_text("chr1\n100\n200\n")
        _normalize_track_to_bed(inp, out)
        assert out.read_text().strip() == ""

    def test_invalid_coordinates_skipped(self, tmp_path):
        inp = tmp_path / "input.bed"
        out = tmp_path / "output.bed"
        inp.write_text("chr1\t100\t200\nchr1\t200\t100\nchr1\tabc\t200\n")
        _normalize_track_to_bed(inp, out)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 1
        assert lines[0] == "chr1\t100\t200"

    def test_mappability_threshold_filters_high_scores(self, tmp_path):
        inp = tmp_path / "input.bedgraph"
        out = tmp_path / "output.bed"
        # Standard bedgraph: chrom start end score
        inp.write_text("chr1\t100\t200\t0.3\nchr1\t200\t300\t0.8\n")
        _normalize_track_to_bed(inp, out, mappability_threshold=0.5)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 1
        assert lines[0] == "chr1\t100\t200"

    def test_mappability_threshold_keeps_low_scores(self, tmp_path):
        inp = tmp_path / "input.bedgraph"
        out = tmp_path / "output.bed"
        inp.write_text("chr1\t100\t200\t0.2\nchr1\t200\t300\t0.1\n")
        _normalize_track_to_bed(inp, out, mappability_threshold=0.5)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2

    def test_bedgraph_with_leading_bin_column(self, tmp_path):
        inp = tmp_path / "input.bedgraph"
        out = tmp_path / "output.bed"
        # UCSC bin format: bin chrom start end score
        inp.write_text("100\tchr1\t100\t200\t0.3\n")
        _normalize_track_to_bed(inp, out, mappability_threshold=0.5)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 1
        assert lines[0] == "chr1\t100\t200"

    def test_mappability_threshold_none_writes_all(self, tmp_path):
        inp = tmp_path / "input.bedgraph"
        out = tmp_path / "output.bed"
        inp.write_text("chr1\t100\t200\t0.9\n")
        _normalize_track_to_bed(inp, out, mappability_threshold=None)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 1


# ─────────────────────────────────────────────────────────────────────────────
# _merge_overlapping_intervals (candidate_loci)
# ─────────────────────────────────────────────────────────────────────────────


class TestMergeOverlappingIntervals:
    def test_empty_input(self):
        assert _merge_overlapping_intervals([], max_span_bp=1000) == []

    def test_single_interval(self):
        result = _merge_overlapping_intervals([(100, 200)], max_span_bp=1000)
        assert result == [(100, 200)]

    def test_two_non_overlapping(self):
        result = _merge_overlapping_intervals([(100, 200), (300, 400)], max_span_bp=1000)
        assert result == [(100, 200), (300, 400)]

    def test_two_overlapping(self):
        result = _merge_overlapping_intervals([(100, 200), (150, 300)], max_span_bp=1000)
        assert result == [(100, 300)]

    def test_touching_at_boundary(self):
        result = _merge_overlapping_intervals([(100, 200), (200, 300)], max_span_bp=1000)
        assert result == [(100, 300)]

    def test_span_cap_prevents_merge(self):
        result = _merge_overlapping_intervals([(100, 300), (250, 600)], max_span_bp=200)
        assert len(result) == 2

    def test_three_way_chain(self):
        result = _merge_overlapping_intervals(
            [(100, 200), (180, 280), (260, 350)], max_span_bp=1000
        )
        assert result == [(100, 350)]

    def test_unsorted_input_sorted(self):
        result = _merge_overlapping_intervals([(300, 400), (100, 200)], max_span_bp=1000)
        assert result == [(100, 200), (300, 400)]

    def test_exact_boundary_no_overlap(self):
        result = _merge_overlapping_intervals([(100, 200), (201, 300)], max_span_bp=1000)
        assert result == [(100, 200), (201, 300)]


# ─────────────────────────────────────────────────────────────────────────────
# _distance_to_closed_interval — genomic-mask boundary cases
# ─────────────────────────────────────────────────────────────────────────────


class TestDistanceGenomicBoundary:
    def test_zero_width_interval(self):
        assert _distance_to_closed_interval(50, 100, 100) == 50

    def test_pos_at_start_of_zero_width(self):
        assert _distance_to_closed_interval(100, 100, 100) == 0

    def test_very_large_coordinates(self):
        assert _distance_to_closed_interval(250_000_000, 100_000_000, 200_000_000) == 50_000_000

    def test_distance_just_outside(self):
        assert _distance_to_closed_interval(99, 100, 100) == 1
        assert _distance_to_closed_interval(101, 100, 100) == 1


# ─────────────────────────────────────────────────────────────────────────────
# _write_candidate_windows_bed
# ─────────────────────────────────────────────────────────────────────────────


class TestWriteCandidateWindowsBed:
    def test_writes_bed_format(self, tmp_path):
        df = pd.DataFrame([
            {"chrom": "chr1", "window_start": 101, "window_end": 200, "row_id": 0},
            {"chrom": "chr2", "window_start": 301, "window_end": 400, "row_id": 1},
        ])
        out = tmp_path / "out.bed"
        _write_candidate_windows_bed(df, out)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2
        parts = lines[0].split("\t")
        assert parts[0] == "chr1"
        assert parts[1] == "100"
        assert parts[2] == "200"
        assert parts[3] == "0"

    def test_empty_dataframe(self, tmp_path):
        df = pd.DataFrame(columns=["chrom", "window_start", "window_end", "row_id"])
        out = tmp_path / "out.bed"
        _write_candidate_windows_bed(df, out)
        assert out.read_text() == ""


# ─────────────────────────────────────────────────────────────────────────────
# _annotate_junk_flags — full integration with mocked bedtools
# ─────────────────────────────────────────────────────────────────────────────


def _junk_loci() -> pd.DataFrame:
    return pd.DataFrame([
        {"chrom": "chr1", "window_start": 100, "window_end": 200},
        {"chrom": "chr1", "window_start": 500, "window_end": 600},
        {"chrom": "chr2", "window_start": 1000, "window_end": 1200},
    ])


def _empty_discordant() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["chrom", "pos", "read_name", "mapq", "discordant_reasons", "template_len",
                 "window_start", "window_end", "mate_chrom", "mate_pos"]
    )


class TestAnnotateJunkFlags:
    def test_all_flags_initialized_false(self):
        loci = _junk_loci()
        disc = _empty_discordant()
        result = _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=None,
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=None,
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=None,
            gap_bed=None,
            gap_min_fraction=0.1,
            encode_blacklist_bed=None,
            encode_blacklist_min_fraction=0.1,
        )
        assert "flag_segdup" in result.columns
        assert "flag_low_mappability" in result.columns
        assert "flag_outside_giab_highconf" in result.columns
        assert "flag_gap_region" in result.columns
        assert "flag_encode_blacklist" in result.columns
        assert "flag_mate_in_segdup" in result.columns
        assert "flag_mate_low_mappability" in result.columns
        assert "flag_mate_in_gap" in result.columns
        assert "flag_mate_in_blacklist" in result.columns
        assert "junk_flag_count" in result.columns
        assert "mate_junk_flag_count" in result.columns

    def test_no_tracks_all_flags_zero(self):
        loci = _junk_loci()
        disc = _empty_discordant()
        result = _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=None,
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=None,
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=None,
            gap_bed=None,
            gap_min_fraction=0.1,
            encode_blacklist_bed=None,
            encode_blacklist_min_fraction=0.1,
        )
        assert (result["junk_flag_count"] == 0).all()
        assert (result["mate_junk_flag_count"] == 0).all()

    def test_empty_loci_returns_expected_columns(self):
        loci = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
        disc = _empty_discordant()
        result = _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=None,
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=None,
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=None,
            gap_bed=None,
            gap_min_fraction=0.1,
            encode_blacklist_bed=None,
            encode_blacklist_min_fraction=0.1,
        )
        assert len(result) == 0
        assert "junk_flag_count" in result.columns

    def test_row_count_preserved(self):
        loci = _junk_loci()
        disc = _empty_discordant()
        result = _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=None,
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=None,
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=None,
            gap_bed=None,
            gap_min_fraction=0.1,
            encode_blacklist_bed=None,
            encode_blacklist_min_fraction=0.1,
        )
        assert len(result) == 3

    def test_segdup_bed_hits_flagged(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        segdup_bed = tmp_path / "segdup.bed"
        segdup_bed.write_text("chr1\t90\t210\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids_with_fraction") as mock:
            mock.return_value = {0}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=segdup_bed,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=None,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=None,
                gap_bed=None,
                gap_min_fraction=0.1,
                encode_blacklist_bed=None,
                encode_blacklist_min_fraction=0.1,
            )
        assert result.loc[result.index[0], "flag_segdup"]
        assert not result.loc[result.index[1], "flag_segdup"]

    def test_gap_bed_hits_flagged(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        gap_bed = tmp_path / "gap.bed"
        gap_bed.write_text("chr2\t990\t1210\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids_with_fraction") as mock:
            mock.return_value = {2}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=None,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=None,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=None,
                gap_bed=gap_bed,
                gap_min_fraction=0.1,
                encode_blacklist_bed=None,
                encode_blacklist_min_fraction=0.1,
            )
        assert result.loc[result.index[2], "flag_gap_region"]
        assert not result.loc[result.index[0], "flag_gap_region"]

    def test_encode_blacklist_bed_hits_flagged(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        bl_bed = tmp_path / "blacklist.bed"
        bl_bed.write_text("chr1\t490\t610\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids_with_fraction") as mock:
            mock.return_value = {1}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=None,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=None,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=None,
                gap_bed=None,
                gap_min_fraction=0.1,
                encode_blacklist_bed=bl_bed,
                encode_blacklist_min_fraction=0.1,
            )
        assert result.loc[result.index[1], "flag_encode_blacklist"]

    def test_giab_highconf_outside_flagged(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        giab_bed = tmp_path / "giab.bed"
        giab_bed.write_text("chr1\t90\t210\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids") as mock:
            mock.return_value = {0}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=None,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=None,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=giab_bed,
                gap_bed=None,
                gap_min_fraction=0.1,
                encode_blacklist_bed=None,
                encode_blacklist_min_fraction=0.1,
            )
        assert not result.loc[result.index[0], "flag_outside_giab_highconf"]
        assert result.loc[result.index[1], "flag_outside_giab_highconf"]

    def test_multiple_flags_aggregate(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        segdup_bed = tmp_path / "segdup.bed"
        segdup_bed.write_text("chr1\t90\t210\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids_with_fraction") as mock:
            mock.return_value = {0}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=segdup_bed,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=None,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=None,
                gap_bed=None,
                gap_min_fraction=0.1,
                encode_blacklist_bed=None,
                encode_blacklist_min_fraction=0.1,
            )
        assert result.loc[result.index[0], "junk_flag_count"] >= 1
        assert result.loc[result.index[1], "junk_flag_count"] == 0

    def test_nonexistent_bed_file_treated_as_none(self):
        loci = _junk_loci()
        disc = _empty_discordant()
        result = _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=Path("/nonexistent/segdup.bed"),
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=Path("/nonexistent/map.bedgraph"),
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=Path("/nonexistent/giab.bed"),
            gap_bed=Path("/nonexistent/gap.bed"),
            gap_min_fraction=0.1,
            encode_blacklist_bed=Path("/nonexistent/blacklist.bed"),
            encode_blacklist_min_fraction=0.1,
        )
        assert (result["junk_flag_count"] == 0).all()

    def test_does_not_mutate_input(self):
        loci = _junk_loci()
        original_cols = set(loci.columns)
        disc = _empty_discordant()
        _annotate_junk_flags(
            loci=loci,
            discordant_disease=disc,
            discordant_control=disc,
            segdup_bed=None,
            segdup_min_fraction=0.1,
            low_mappability_bedgraph=None,
            low_mappability_threshold=0.5,
            low_mappability_min_fraction=0.5,
            giab_highconf_bed=None,
            gap_bed=None,
            gap_min_fraction=0.1,
            encode_blacklist_bed=None,
            encode_blacklist_min_fraction=0.1,
        )
        assert set(loci.columns) == original_cols

    def test_mappability_bedgraph_hits_flagged(self, tmp_path):
        loci = _junk_loci()
        disc = _empty_discordant()
        map_bed = tmp_path / "mappability.bedgraph"
        map_bed.write_text("chr1\t490\t610\n")
        with patch("retro_miner.candidate_loci._get_overlapping_row_ids_with_fraction") as mock:
            mock.return_value = {1}
            result = _annotate_junk_flags(
                loci=loci,
                discordant_disease=disc,
                discordant_control=disc,
                segdup_bed=None,
                segdup_min_fraction=0.1,
                low_mappability_bedgraph=map_bed,
                low_mappability_threshold=0.5,
                low_mappability_min_fraction=0.5,
                giab_highconf_bed=None,
                gap_bed=None,
                gap_min_fraction=0.1,
                encode_blacklist_bed=None,
                encode_blacklist_min_fraction=0.1,
            )
        assert result.loc[result.index[1], "flag_low_mappability"]
        assert not result.loc[result.index[0], "flag_low_mappability"]
