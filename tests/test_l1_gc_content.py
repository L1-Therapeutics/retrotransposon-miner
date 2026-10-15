"""Tests for l1_gc_content pure functions (no reference/callset needed)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
import l1_gc_content as m  # noqa: E402


def test_g4hunter_base_scores_signed_runs_capped():
    # G-run of length 6 caps at 4; C-run of length 2 is -2; A/T are 0.
    seq = "ATGGGGGGCCAA"
    b = m.g4hunter_base_scores(seq)
    assert b[0] == 0 and b[1] == 0            # A, T
    assert list(b[2:8]) == [4, 4, 4, 4, 4, 4]  # 6 G's -> capped 4
    assert list(b[8:10]) == [-2, -2]           # 2 C's
    assert list(b[10:12]) == [0, 0]            # A, A
    assert len(b) == len(seq)


def test_g4hunter_symmetric_in_strand():
    # A G-run and the C-run of its reverse complement have equal |score|.
    seq = "AAAAAAAAAAAAGGGGGGGGCCCCCCCCAAAAAAAAAAAA"
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert len(seq) >= m.G4_WINDOW
    assert m.max_abs_g4hunter(seq) == pytest.approx(m.max_abs_g4hunter(rc))


def test_max_abs_g4hunter_short_is_nan():
    assert np.isnan(m.max_abs_g4hunter("ACGT"))


def test_max_abs_g4hunter_g_rich_beats_mixed():
    g_rich = "GGGGGGGGGGGGGGGGGGGGGGGGG"   # 25 G's -> window score 4.0
    mixed = "ACGTACGTACGTACGTACGTACGTAA"
    assert m.max_abs_g4hunter(g_rich) == pytest.approx(4.0)
    assert m.max_abs_g4hunter(g_rich) > m.max_abs_g4hunter(mixed)


def test_window_features_fractions_and_n_rejection():
    f = m.window_features("AAACCCGGGTTT")
    assert f["gc_frac"] == pytest.approx(0.5)
    assert f["at_frac"] == pytest.approx(0.5)
    assert f["g_frac"] == pytest.approx(0.25)
    assert f["c_frac"] == pytest.approx(0.25)
    assert m.window_features("NNNNNNNNNNACGT") is None   # >50% N
    assert m.window_features("") is None


def test_parse_subfamily_strips_end_and_family_tag():
    assert m.parse_subfamily("L1MEc_5end#LINE/L1") == "L1MEc"
    assert m.parse_subfamily("L1HS#LINE/L1") == "L1HS"
    assert m.parse_subfamily("AluYb8#SINE/Alu") == "AluYb8"
    assert m.parse_subfamily(None) is None


def test_l1_age_class_prefixes():
    assert m.l1_age_class("L1HS") == "primate/young"
    assert m.l1_age_class("L1PA3") == "primate/young"
    assert m.l1_age_class("L1MC1") == "intermediate"
    assert m.l1_age_class("L1M4c") == "old"
    assert m.l1_age_class("AluY") == "other"
    assert m.l1_age_class(None) == "unknown"
