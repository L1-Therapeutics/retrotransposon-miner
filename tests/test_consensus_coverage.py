"""Tests for scripts/consensus_coverage.py.

Covers the three pure functions that carry the LINE-1 coverage result:
``read_intervals`` (target filtering across the read and mate column pairs),
``coverage_profile`` (1-based depth with clamping at both ends) and
``largest_gap`` (interior dead zones). ``profile_locus`` is exercised only for
its missing-input contract; its happy path needs a locus fixture.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "consensus_coverage.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("consensus_coverage", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["consensus_coverage"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


cc = _load_module()


def _detail(rows):
    cols = ["mei_target", "mei_start", "mei_end",
            "mate_mei_target", "mate_mei_start", "mate_mei_end"]
    return pd.DataFrame(rows, columns=cols)


def test_read_intervals_pools_read_and_mate_columns():
    detail = _detail([
        ["L1HS_full", 10, 20, "L1HS_full", 30, 40],
        ["L1HS_full", 50, 60, None, None, None],
    ])
    assert sorted(cc.read_intervals(detail, "L1HS_full")) == [(10, 20), (30, 40), (50, 60)]


def test_read_intervals_excludes_other_targets():
    """The whole point of the metric: an L1PA3 read must not count toward L1HS."""
    detail = _detail([
        ["L1HS_full", 10, 20, None, None, None],
        ["L1PA3_full", 100, 200, "L1PA3_full", 300, 400],
    ])
    assert cc.read_intervals(detail, "L1HS_full") == [(10, 20)]
    assert sorted(cc.read_intervals(detail, "L1PA3_full")) == [(100, 200), (300, 400)]


def test_read_intervals_drops_incomplete_and_inverted_rows():
    detail = _detail([
        ["L1HS_full", None, 20, None, None, None],   # missing start
        ["L1HS_full", 80, 70, None, None, None],     # end before start
        ["L1HS_full", 10, 20, None, None, None],     # the only keeper
    ])
    assert cc.read_intervals(detail, "L1HS_full") == [(10, 20)]


def test_read_intervals_tolerates_absent_mate_columns():
    detail = pd.DataFrame([["L1HS_full", 10, 20]],
                          columns=["mei_target", "mei_start", "mei_end"])
    assert cc.read_intervals(detail, "L1HS_full") == [(10, 20)]


def test_coverage_profile_is_one_based_and_inclusive():
    depth = cc.coverage_profile([(2, 4)], length=6)
    assert depth.tolist() == [0, 1, 1, 1, 0, 0]


def test_coverage_profile_accumulates_overlaps():
    depth = cc.coverage_profile([(1, 3), (2, 5)], length=5)
    assert depth.tolist() == [1, 2, 2, 1, 1]


def test_coverage_profile_clamps_to_the_element():
    """Reads running past either consensus end must not inflate covered bases."""
    depth = cc.coverage_profile([(-5, 2), (4, 999)], length=5)
    assert depth.tolist() == [1, 1, 0, 1, 1]
    assert depth.size == 5


def test_largest_gap_is_zero_for_contiguous_coverage():
    assert cc.largest_gap(np.array([1, 1, 1, 1])) == 0


def test_largest_gap_measures_the_interior_dead_zone():
    # covered at 0 and 4 -> three uncovered bases between them
    assert cc.largest_gap(np.array([1, 0, 0, 0, 1])) == 3


def test_largest_gap_needs_two_covered_bases():
    assert cc.largest_gap(np.array([0, 0, 0])) == 0
    assert cc.largest_gap(np.array([0, 1, 0])) == 0


def test_profile_locus_returns_none_when_inputs_are_missing(tmp_path):
    assert cc.profile_locus(tmp_path) is None
