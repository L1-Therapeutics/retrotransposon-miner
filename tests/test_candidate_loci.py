"""Unit tests for pure-logic helpers in candidate_loci.

These tests exercise the clustering and interval-merging algorithms at the
heart of the locus-building pipeline.  No BAM files or real evidence are
required; every function under test is a pure (or near-pure) data
transformation.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from retro_miner.candidate_loci import (
    _cluster_sorted_positions,
    _distance_to_closed_interval,
    _merge_overlapping_loci,
    _read_passing_counts,
    _split_cluster_positions,
)


# ─────────────────────────────────────────────────────────────────────────────
# _cluster_sorted_positions
# ─────────────────────────────────────────────────────────────────────────────

class TestClusterSortedPositions:
    def test_empty_input_returns_empty(self):
        assert _cluster_sorted_positions([], max_gap_bp=100) == []

    def test_single_position_returns_one_cluster(self):
        assert _cluster_sorted_positions([500], max_gap_bp=100) == [[500]]

    def test_all_within_gap_returns_one_cluster(self):
        result = _cluster_sorted_positions([100, 150, 200], max_gap_bp=100)
        assert result == [[100, 150, 200]]

    def test_positions_exactly_at_gap_boundary_stay_in_same_cluster(self):
        # gap between 100 and 200 is exactly max_gap_bp=100 → same cluster
        result = _cluster_sorted_positions([100, 200], max_gap_bp=100)
        assert result == [[100, 200]]

    def test_positions_one_over_gap_split_into_separate_clusters(self):
        # gap between 100 and 202 is 102 > max_gap_bp=100 → separate
        result = _cluster_sorted_positions([100, 202], max_gap_bp=100)
        assert result == [[100], [202]]

    def test_multiple_clusters(self):
        result = _cluster_sorted_positions(
            [100, 120, 140, 500, 510, 900], max_gap_bp=50
        )
        assert result == [[100, 120, 140], [500, 510], [900]]

    def test_negative_max_gap_treated_as_zero(self):
        # max_gap_bp < 0 → treated as 0; only equal positions cluster together
        result = _cluster_sorted_positions([100, 100, 200], max_gap_bp=-5)
        assert result == [[100, 100], [200]]

    def test_adjacent_positions_with_zero_gap(self):
        # positions 100, 101 differ by 1; with max_gap_bp=0 they must split
        result = _cluster_sorted_positions([100, 101], max_gap_bp=0)
        assert result == [[100], [101]]

    def test_duplicate_positions_stay_in_same_cluster(self):
        result = _cluster_sorted_positions([100, 100, 100], max_gap_bp=0)
        assert result == [[100, 100, 100]]


# ─────────────────────────────────────────────────────────────────────────────
# _split_cluster_positions
# ─────────────────────────────────────────────────────────────────────────────

class TestSplitClusterPositions:
    def test_empty_input_returns_empty(self):
        assert _split_cluster_positions([], valley_gap_bp=100, max_locus_span_bp=1000) == []

    def test_single_position_returns_one_cluster(self):
        result = _split_cluster_positions([42], valley_gap_bp=100, max_locus_span_bp=1000)
        assert result == [[42]]

    def test_unsorted_input_is_sorted_automatically(self):
        # gaps are 100 bp each; valley_gap_bp=200 keeps all in one cluster
        result = _split_cluster_positions(
            [300, 100, 200], valley_gap_bp=200, max_locus_span_bp=1000
        )
        assert result == [[100, 200, 300]]

    def test_valley_gap_splits_well_separated_positions(self):
        result = _split_cluster_positions(
            [100, 150, 600, 650], valley_gap_bp=100, max_locus_span_bp=2000
        )
        # gap 150→600 = 450 > 100 → split
        assert result == [[100, 150], [600, 650]]

    def test_span_cap_forces_split_even_within_valley_gap(self):
        # positions 100 and 300 differ by 200, valley_gap=500, max_span=150
        # span_if_added = 300 - 100 + 1 = 201 > 150 → forced split
        result = _split_cluster_positions(
            [100, 300], valley_gap_bp=500, max_locus_span_bp=150
        )
        assert result == [[100], [300]]

    def test_span_exactly_at_cap_stays_in_one_cluster(self):
        # span 100→200 = 101 bp; cap = 101 → should stay together
        result = _split_cluster_positions(
            [100, 200], valley_gap_bp=500, max_locus_span_bp=101
        )
        assert result == [[100, 200]]

    def test_zero_or_negative_max_locus_span_forces_one_per_cluster(self):
        # max_locus_span_bp <= 0 is clamped to 1
        result = _split_cluster_positions(
            [100, 200, 300], valley_gap_bp=0, max_locus_span_bp=0
        )
        # every pair exceeds span=1, so each gets its own cluster
        assert result == [[100], [200], [300]]


# ─────────────────────────────────────────────────────────────────────────────
# _distance_to_closed_interval
# ─────────────────────────────────────────────────────────────────────────────

class TestDistanceToClosedInterval:
    def test_pos_inside_interval_returns_zero(self):
        assert _distance_to_closed_interval(150, 100, 200) == 0

    def test_pos_at_left_endpoint_returns_zero(self):
        assert _distance_to_closed_interval(100, 100, 200) == 0

    def test_pos_at_right_endpoint_returns_zero(self):
        assert _distance_to_closed_interval(200, 100, 200) == 0

    def test_pos_left_of_interval(self):
        assert _distance_to_closed_interval(50, 100, 200) == 50

    def test_pos_right_of_interval(self):
        assert _distance_to_closed_interval(250, 100, 200) == 50

    def test_pos_one_left_of_left_endpoint(self):
        assert _distance_to_closed_interval(99, 100, 200) == 1

    def test_pos_one_right_of_right_endpoint(self):
        assert _distance_to_closed_interval(201, 100, 200) == 1


# ─────────────────────────────────────────────────────────────────────────────
# _merge_overlapping_loci
# ─────────────────────────────────────────────────────────────────────────────

def _loci(rows: list[tuple[str, int, int]]) -> pd.DataFrame:
    """Helper: build a minimal loci DataFrame from (chrom, start, end) tuples."""
    return pd.DataFrame(rows, columns=["chrom", "window_start", "window_end"])


class TestMergeOverlappingLoci:
    def test_empty_dataframe_returns_empty(self):
        result = _merge_overlapping_loci(pd.DataFrame(), max_locus_span_bp=2000)
        assert result.empty

    def test_none_returns_empty(self):
        result = _merge_overlapping_loci(None, max_locus_span_bp=2000)
        assert result.empty

    def test_single_row_unchanged(self):
        df = _loci([("chr1", 100, 300)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=2000)
        assert len(result) == 1
        assert result.iloc[0]["window_start"] == 100
        assert result.iloc[0]["window_end"] == 300

    def test_non_overlapping_stays_separate(self):
        # gap of 100 bp between 300 and 401 — no overlap
        df = _loci([("chr1", 100, 300), ("chr1", 401, 600)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=5000)
        assert len(result) == 2

    def test_overlapping_intervals_merge(self):
        df = _loci([("chr1", 100, 300), ("chr1", 250, 450)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=5000)
        assert len(result) == 1
        r = result.iloc[0]
        assert r["window_start"] == 100
        assert r["window_end"] == 450

    def test_touching_intervals_merge(self):
        # start of second = end of first → condition start <= cur_end is True
        df = _loci([("chr1", 100, 200), ("chr1", 200, 300)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=5000)
        assert len(result) == 1
        r = result.iloc[0]
        assert r["window_start"] == 100
        assert r["window_end"] == 300

    def test_span_cap_prevents_merger(self):
        # merged span would be 400 - 100 + 1 = 301 > cap=200
        df = _loci([("chr1", 100, 250), ("chr1", 200, 400)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=200)
        assert len(result) == 2

    def test_different_chromosomes_stay_separate(self):
        df = _loci([("chr1", 100, 300), ("chr2", 100, 300)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=5000)
        assert len(result) == 2
        chroms = set(result["chrom"])
        assert chroms == {"chr1", "chr2"}

    def test_three_way_merge(self):
        df = _loci([("chr1", 100, 200), ("chr1", 180, 280), ("chr1", 260, 350)])
        result = _merge_overlapping_loci(df, max_locus_span_bp=5000)
        assert len(result) == 1
        assert result.iloc[0]["window_start"] == 100
        assert result.iloc[0]["window_end"] == 350


# ─────────────────────────────────────────────────────────────────────────────
# _read_passing_counts
# ─────────────────────────────────────────────────────────────────────────────

class TestReadPassingCounts:
    def test_missing_file_raises_with_helpful_message(self):
        with pytest.raises(FileNotFoundError) as exc_info:
            _read_passing_counts(Path("/nonexistent/split_evidence.summary.tsv"))
        msg = str(exc_info.value)
        assert "Split evidence summary not found" in msg
        assert "rtm extract-split-evidence" in msg

    def test_reads_sample_passing_counts(self, tmp_path):
        summary = tmp_path / "split_evidence.summary.tsv"
        summary.write_text(
            "sample\ttotal_reads_scanned\tpassing_reads\tsplit_evidence_rows\n"
            "disease\t10000\t1234\t88\n"
            "control\t8000\t987\t62\n"
        )
        result = _read_passing_counts(summary)
        assert result == {"disease": 1234, "control": 987}

    def test_returns_int_not_float(self, tmp_path):
        summary = tmp_path / "split_evidence.summary.tsv"
        summary.write_text("sample\tpassing_reads\ndisease\t500\n")
        result = _read_passing_counts(summary)
        assert isinstance(result["disease"], int)
