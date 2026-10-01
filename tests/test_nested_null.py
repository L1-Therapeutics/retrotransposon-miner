"""Self-insertion null models for the nested-orientation analysis."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_nested_orientation.py"


def _mod():
    spec = importlib.util.spec_from_file_location("analyze_nested_orientation", _SCRIPT)
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


mod = _mod()


def test_merge_intervals_closes_gaps_and_overlaps():
    assert mod._merge_intervals([(10, 20), (15, 25), (30, 40)]) == [(10, 25), (30, 40)]
    assert mod._merge_intervals([(5, 10), (10, 15)]) == [(5, 15)]
    assert mod._merge_intervals([]) == []


def test_merge_intervals_is_unioned_not_summed():
    """Overlapping rmsk intervals must not double-count coverage."""
    merged = mod._merge_intervals([(0, 100), (50, 150)])
    starts = np.array([0], dtype=np.int64)
    ends = np.array([200], dtype=np.int64)
    cov = mod._coverage_per_bin(starts, ends, merged)
    assert int(cov[0]) == 150


def test_coverage_per_bin_matches_brute_force():
    rng = np.random.default_rng(7)
    ivs = []
    for _ in range(200):
        a = int(rng.integers(0, 5000))
        b = int(rng.integers(0, 5000))
        if b > a:
            ivs.append((a, b))
    merged = mod._merge_intervals(ivs)
    starts = np.arange(0, 5000, 137, dtype=np.int64)
    ends = np.minimum(starts + 137, 5000)
    fast = mod._coverage_per_bin(starts, ends, merged)
    slow = [
        sum(max(0, min(e, int(b)) - max(s, int(a))) for s, e in merged)
        for a, b in zip(starts, ends)
    ]
    assert fast.tolist() == slow


def test_coverage_per_bin_handles_no_intervals():
    starts = np.arange(0, 1000, 100, dtype=np.int64)
    ends = starts + 100
    assert mod._coverage_per_bin(starts, ends, []).tolist() == [0] * 10


def test_chunk_bounds_covers_every_bin_once():
    for n_bins, bin_size, chunk in [(10, 200, 1000), (195805, 200, 20_000_000), (0, 200, 1000)]:
        bounds = mod._chunk_bounds(n_bins, bin_size, chunk)
        if n_bins == 0:
            assert bounds == []
            continue
        assert bounds[0][0] == 0
        assert bounds[-1][1] == n_bins
        for (_, end), (start, _) in zip(bounds, bounds[1:]):
            assert end == start, "chunks must be contiguous with no gap or overlap"


def test_selfinsertion_null_uses_per_call_chromosome():
    """Expected count must be the sum of per-call probabilities.

    Chromosome A has all the MEI content; chromosome B has none. Calls are
    split across both, so the expected count is driven by the A calls only.
    """
    sizes = {"chrA": 1000, "chrB": 1000}
    content = {("ALU", "chrA"): 500}
    df = pd.DataFrame(
        {
            "family": ["ALU"] * 4,
            "chrom": ["chrA", "chrA", "chrB", "chrB"],
            "nested_same_class_orientation": [
                "nested_sense",
                "nested_sense",
                "unnested",
                "unnested",
            ],
        }
    )
    out = mod.selfinsertion_null(df, content, sizes)
    row = out[out["family"] == "ALU"].iloc[0]
    # Only the two chrA calls can nest: 2 * 0.5 = 1 expected nested.
    assert row["expected_nested_count"] == pytest.approx(1.0)
    assert row["expected_nested_fraction_uniform"] == pytest.approx(0.25)


def test_selfinsertion_null_skips_families_with_no_content():
    sizes = {"chrA": 1000}
    df = pd.DataFrame(
        {
            "family": ["SVA", "ALU"],
            "chrom": ["chrA", "chrA"],
            "nested_same_class_orientation": ["nested_sense", "nested_sense"],
        }
    )
    out = mod.selfinsertion_null(df, {("ALU", "chrA"): 100}, sizes)
    assert list(out["family"]) == ["ALU"]


def test_selfinsertion_null_counts_unknown_as_nested():
    sizes = {"chrA": 1000}
    df = pd.DataFrame(
        {
            "family": ["ALU", "ALU"],
            "chrom": ["chrA", "chrA"],
            "nested_same_class_orientation": ["nested_unknown", "unnested"],
        }
    )
    out = mod.selfinsertion_null(df, {("ALU", "chrA"): 500}, sizes)
    # nested_unknown is nested, just not in either orientation class.
    assert int(out.iloc[0]["n_nested"]) == 1


def test_gc_matched_null_excludes_bins_outside_tolerance():
    """A GC-matched control pool must not contain wildly different bins."""
    bins = pd.DataFrame(
        {
            "chrom": ["chrA"] * 4,
            "start": [0, 100, 200, 300],
            "gc": [0.30, 0.31, 0.70, 0.71],
            "span": [100, 100, 100, 100],
            "ALU": [50, 50, 100, 100],
            "LINE1": [0, 0, 0, 0],
            "SVA": [0, 0, 0, 0],
        }
    )
    df = pd.DataFrame(
        {
            "family": ["ALU"],
            "chrom": ["chrA"],
            "gc": [0.305],
            "nested_same_class_orientation": ["nested_sense"],
        }
    )
    out = mod.gc_matched_null(bins, df, n_matched=1000, gc_tolerance=0.05, seed=1)
    row = out[out["family"] == "ALU"].iloc[0]
    # Only the two ~0.30 bins match, each 50% ALU-covered.
    assert row["expected_nested_fraction_gc_matched"] == pytest.approx(0.5)


def test_gc_matched_null_falls_back_when_no_bin_matches():
    bins = pd.DataFrame(
        {
            "chrom": ["chrA"] * 2,
            "start": [0, 100],
            "gc": [0.10, 0.12],
            "span": [100, 100],
            "ALU": [10, 90],
            "LINE1": [0, 0],
            "SVA": [0, 0],
        }
    )
    df = pd.DataFrame(
        {
            "family": ["ALU"],
            "chrom": ["chrA"],
            "gc": [0.90],  # nothing remotely like this
            "nested_same_class_orientation": ["unnested"],
        }
    )
    out = mod.gc_matched_null(bins, df, n_matched=1000, gc_tolerance=0.01, seed=1)
    assert out.iloc[0]["expected_nested_fraction_gc_matched"] == pytest.approx(0.5)


def test_genome_sizes_parses_two_column_file(tmp_path: Path):
    p = tmp_path / "sizes"
    p.write_text("chr1\t100\nchr2\t200\n", encoding="utf-8")
    assert mod.genome_sizes(p) == {"chr1": 100, "chr2": 200}


def test_genome_sizes_rejects_empty_file(tmp_path: Path):
    p = tmp_path / "sizes"
    p.write_text("", encoding="utf-8")
    with pytest.raises(ValueError):
        mod.genome_sizes(p)
