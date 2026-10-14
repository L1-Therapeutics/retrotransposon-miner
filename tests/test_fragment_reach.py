"""Tests for scripts/fragment_reach.py.

``summarise`` is the function the reach argument rests on, so these pin its
arithmetic: the reported ceiling is the 99th-percentile fragment minus the
anchor allowance, expressed as a fraction of the element. BAM reading is not
exercised here -- the root conftest substitutes a MagicMock when pysam is
unavailable, so a BAM test would assert against a mock rather than htslib.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "fragment_reach.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("fragment_reach", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fragment_reach"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


fr = _load_module()


def test_summarise_reports_the_requested_percentiles():
    out = fr.summarise(np.arange(1, 101), element_length=6018, anchor_bp=0)
    assert out["n_fragments"] == 100
    assert out["median"] == pytest.approx(50.5)
    assert out["p99"] == pytest.approx(99.01)
    assert out["max"] == 100


def test_reach_cap_subtracts_the_anchor_from_the_p99():
    lengths = np.full(100, 700)
    out = fr.summarise(lengths, element_length=6018, anchor_bp=30)
    assert out["p99"] == pytest.approx(700)
    assert out["reach_cap_bp"] == pytest.approx(670)


def test_reachable_fraction_is_reach_over_element_length():
    out = fr.summarise(np.full(100, 726), element_length=6018, anchor_bp=30)
    assert out["reachable_fraction"] == pytest.approx(696 / 6018)
    assert out["reachable_fraction"] == pytest.approx(0.1157, abs=1e-4)


def test_a_larger_anchor_allowance_shortens_the_reach():
    """Reach is monotonically decreasing in the anchor: the sensitivity sweep
    reported to the team depends on this direction."""
    lengths = np.full(100, 726)
    fracs = [fr.summarise(lengths, 6018, a)["reachable_fraction"] for a in (0, 30, 100)]
    assert fracs == sorted(fracs, reverse=True)


def test_an_alu_is_fully_reachable_where_a_line1_is_not():
    """Same library, two element lengths -- the contrast the analysis turns on."""
    lengths = np.full(100, 726)
    alu = fr.summarise(lengths, element_length=281, anchor_bp=30)
    line1 = fr.summarise(lengths, element_length=6018, anchor_bp=30)
    assert alu["reachable_fraction"] > 1.0
    assert line1["reachable_fraction"] < 0.2


def test_summarise_rejects_an_empty_library():
    with pytest.raises(SystemExit):
        fr.summarise(np.array([]), element_length=6018, anchor_bp=30)


def test_zero_element_length_yields_nan_rather_than_dividing_by_zero():
    out = fr.summarise(np.full(10, 400), element_length=0, anchor_bp=30)
    assert np.isnan(out["reachable_fraction"])
