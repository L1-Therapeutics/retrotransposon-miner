"""Tests for scripts/line1_span_support.py.

The cohort loader needs real VCFs, so these cover the three pure functions
downstream of it: the consensus-length proxy (max MEI_3P per subfamily,
sentinels excluded), the per-call annotation it feeds, and the summary counts.
The proxy is a stand-in for the Dfam panel and these tests pin that contract,
not the biological length.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "line1_span_support.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("line1_span_support", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["line1_span_support"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


ls = _load_module()


def _calls(rows):
    frame = pd.DataFrame(rows)
    frame["sentinel"] = (frame["MEI_SPAN"] <= 0) | (frame["MEI_5P"] < 0)
    return frame


def test_consensus_length_is_the_max_three_prime_per_subfamily():
    calls = _calls([
        {"sample": "A", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 100, "MEI_5P": 1, "MEI_3P": 5000},
        {"sample": "B", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 200, "MEI_5P": 1, "MEI_3P": 6018},
        {"sample": "C", "subfamily": "AluY", "MEIFAMILY": "ALU",
         "MEI_SPAN": 50, "MEI_5P": 1, "MEI_3P": 281},
    ])
    lengths = ls.estimate_consensus_lengths(calls)
    assert lengths["L1HS_full"] == 6018
    assert lengths["AluY"] == 281


def test_sentinels_do_not_contribute_to_the_consensus_proxy():
    """A sentinel's coordinates are placeholders; letting one set the axis
    length would silently rescale every span fraction in that subfamily."""
    calls = _calls([
        {"sample": "A", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 100, "MEI_5P": 1, "MEI_3P": 6018},
        {"sample": "B", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 0, "MEI_5P": -1, "MEI_3P": 99999},
    ])
    assert ls.estimate_consensus_lengths(calls)["L1HS_full"] == 6018


def test_annotate_drops_sentinels_and_flags_calls_ending_at_consensus():
    calls = _calls([
        {"sample": "A", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 6018, "MEI_5P": 1, "MEI_3P": 6018},
        {"sample": "B", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 3000, "MEI_5P": 1, "MEI_3P": 5000},
        {"sample": "C", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 0, "MEI_5P": -1, "MEI_3P": 10},
    ])
    out = ls.annotate(calls)
    assert len(out) == 2                      # sentinel removed
    assert out.ends_at_consensus.tolist() == [True, False]
    assert out.span_fraction.tolist() == [pytest.approx(1.0), pytest.approx(3000 / 6018)]


def test_summarise_counts_only_full_length_line1_at_or_above_the_threshold():
    calls = _calls([
        {"sample": "A", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 6018, "MEI_5P": 1, "MEI_3P": 6018},
        {"sample": "B", "subfamily": "L1HS_full", "MEIFAMILY": "LINE1",
         "MEI_SPAN": 400, "MEI_5P": 1, "MEI_3P": 6018},
        {"sample": "A", "subfamily": "AluY", "MEIFAMILY": "ALU",
         "MEI_SPAN": 281, "MEI_5P": 1, "MEI_3P": 281},
    ])
    out = ls.summarise(calls, ls.annotate(calls), min_full_span=5500)
    assert out["genomes"] == 2
    assert out["calls_total"] == 3
    assert out["line1_usable"] == 2
    assert out["line1_apparently_full"] == 1
    assert out["line1_full_ending_at_consensus"] == 1
    assert out["line1_full_per_genome"] == {"A": 1}


def test_an_alu_is_never_counted_as_a_full_length_line1():
    calls = _calls([
        {"sample": "A", "subfamily": "AluY", "MEIFAMILY": "ALU",
         "MEI_SPAN": 9000, "MEI_5P": 1, "MEI_3P": 9000},
    ])
    out = ls.summarise(calls, ls.annotate(calls), min_full_span=5500)
    assert out["line1_usable"] == 0
    assert out["line1_apparently_full"] == 0
