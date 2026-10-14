"""Tests for genome_gc_content pure functions (no real reference needed).

A tiny in-memory 2bit stand-in lets the binning, masking and GC maths be checked
on hand-computed sequences; the real hg38 run is the integration test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genome_gc_content as m  # noqa: E402


class FakeSeq:
    def __init__(self, s: str):
        self.s = s

    def __len__(self) -> int:
        return len(self.s)

    def __getitem__(self, k):
        return self.s[k]


class Fake2bit:
    def __init__(self, seqs: dict[str, str]):
        self.seqs = seqs

    def __contains__(self, c: str) -> bool:
        return c in self.seqs

    def __getitem__(self, c: str) -> FakeSeq:
        return FakeSeq(self.seqs[c])


def _toy():
    """40 bp, four 10 bp bins: all GC, all AT, half GC, mostly N."""
    seq = "GGGGGGGGGG" + "AAAAAAAAAA" + "ACGTACGTAC" + "NNNNNACGTN"
    return Fake2bit({"chrT": seq}), seq


def test_count_bases_per_bin_is_case_insensitive():
    tb = Fake2bit({"chrT": "aCgTnNaCgT"})
    counts = m.count_bases_per_bin(tb, "chrT", 10, bin_bp=10)
    assert counts["A"][0] == 2 and counts["C"][0] == 2
    assert counts["G"][0] == 2 and counts["T"][0] == 2
    assert counts["N"][0] == 2


def test_bin_table_gc_fraction_and_validity():
    tb, _ = _toy()
    b = m.bin_table(tb, "chrT", 40, bin_bp=10, gaps=None)
    assert list(b.gc_frac.round(3)) == [1.0, 0.0, 0.5, 0.5]
    # last bin is only 4 ACGT of 10 -> below the 50% rule
    assert list(b.valid) == [True, True, True, False]


def test_gap_mask_marks_overlapping_bins():
    mask = m.gap_mask("chrT", nbins=4, bin_bp=10, gaps={"chrT": [[5, 25]]})
    assert list(mask) == [True, True, True, False]


def test_sum_scope_counts_only_valid_over_acgt():
    tb, _ = _toy()
    b = m.bin_table(tb, "chrT", 40, bin_bp=10, gaps=None)
    s = m.sum_scope(b, "chrT")
    assert s["bp_valid"] == 30          # three valid bins of 10 bp
    assert s["gc_bp"] == 15             # 10 + 0 + 5
    assert s["gc_frac"] == pytest.approx(0.5)


def test_gc_histogram_mass_sums_to_one():
    tb, _ = _toy()
    b = m.bin_table(tb, "chrT", 40, bin_bp=10, gaps=None)
    h = m.gc_histogram(b, step=0.1)
    assert h.bp_fraction.sum() == pytest.approx(1.0)
    # the all-GC bin and the all-AT bin each hold one third of the valid bp
    assert h.bp_fraction.iloc[0] == pytest.approx(1 / 3)
    assert h.bp_fraction.iloc[9] == pytest.approx(1 / 3)


def test_attach_and_density_and_summary():
    tb, _ = _toy()
    bins = m.bin_table(tb, "chrT", 40, bin_bp=10, gaps=None)
    sites = pd.DataFrame({
        "chrom": ["chrT", "chrT", "chrT", "chrT"],
        "pos": [3, 15, 25, 37],                     # bins 0,1,2,3(->invalid)
        "family": ["LINE1", "LINE1", "LINE1", "ALU"],
    })
    labeled = m.attach_gc_bin(sites, bins, bin_bp=10, step=0.1)
    assert labeled.loc[0, "gc_frac"] == pytest.approx(1.0)
    assert labeled.loc[1, "gc_frac"] == pytest.approx(0.0)
    assert pd.isna(labeled.loc[3, "gc_frac"])       # masked bin -> NaN

    dens = m.mei_density_by_gc(labeled, bins, step=0.1)
    # GC bin 0.0 holds 10 bp with one LINE1 -> 100,000 per Mb
    row = dens.iloc[0]
    assert row.gc_frac_low == pytest.approx(0.0)
    assert row.LINE1_n == 1
    assert row.LINE1_per_mb == pytest.approx(1e5)

    summ = m.mei_gc_summary(labeled, bins, band=(0.4, 0.6))
    l1 = summ[summ.family == "LINE1"].iloc[0]
    assert l1.n == 3
    assert l1.mean_local_gc == pytest.approx((1.0 + 0.0 + 0.5) / 3)
    assert l1.genome_gc == pytest.approx(0.5)
    # 1 of 3 LINE1 in the band vs 1 of 3 valid bins of genome bp -> 1.0
    assert l1.band_enrichment == pytest.approx(1.0)


def test_robust_peak_ignores_thin_bins():
    dens = pd.DataFrame({
        "gc_frac_low": [0.10, 0.30, 0.50],
        "bp_valid": [1000, 10_000_000, 10_000_000],
        "LINE1_per_mb": [999.0, 2.0, 1.0],          # huge but only 1 kb of sequence
        "ALU_per_mb": [0.0, 1.0, 2.0],
        "SVA_per_mb": [0.0, 0.5, 0.5],
    })
    peak = m.robust_peak(dens, min_bp=5_000_000)
    assert peak["LINE1"] == 0.30                     # 0.10 bin excluded
    assert peak["ALU"] == 0.50
