"""Tests for scripts/exon_depletion_null.py.

Pure-function coverage of the pieces the exon-avoidance result rests on: local GC with N handling, gap
rejection, GC binning and the matched expectation (checked against a null where GC fully determines the
exon probability), feature labelling, and the null sampler on an in-memory "genome" (seeded, gap-free).
The 2bit reader itself is not exercised; a dict of strings stands in for it.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "exon_depletion_null.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("exon_depletion_null", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["exon_depletion_null"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


edn = _load_module()


def test_local_gc_ignores_case_and_n():
    assert edn.local_gc("GGCC") == 1.0 and edn.local_gc("aatt") == 0.0
    assert edn.local_gc("ggAA") == pytest.approx(0.5)
    assert edn.local_gc("GGNNNNNNNN") != edn.local_gc("GGNNNNNNNN")      # more than half N -> NaN
    assert edn.local_gc("GCATNN") == pytest.approx(0.5)                   # N excluded from the denominator


def test_in_gap_uses_half_open_zero_based_intervals():
    gaps = np.array([[100, 200], [500, 510]])
    assert [edn.in_gap(gaps, p) for p in (99, 100, 199, 200, 505, 510)] == [False, True, True, False, True, False]
    assert edn.in_gap(np.zeros((0, 2), dtype=int), 5) is False


def test_window_gcs_returns_one_value_per_window_and_clips_at_chromosome_ends():
    seq = "A" * 400 + "G" * 400                                            # AT-only left half, GC-only right half
    near_left = edn.window_gcs(seq, 5)
    centred = edn.window_gcs(seq, 400)
    assert len(near_left) == 2 and near_left == [0.0, 0.0]               # window clipped at position 0
    assert centred[0] == pytest.approx(0.5, abs=0.01) and centred[1] == pytest.approx(0.5, abs=0.02)


def test_matched_expectation_equals_the_binned_rate_when_gc_determines_exons():
    gc = np.linspace(0.2, 0.8, 2000)
    null = pd.DataFrame({"gc250": gc, "exon": gc > 0.65, "genic": True})   # GC fully decides exon status
    sites = pd.Series([0.30] * 10 + [0.75] * 10)
    exp, var = edn.matched_expectation(null, sites, "gc250", n_bins=10)
    assert exp == pytest.approx(10.0, abs=1.0)                            # only the 10 GC-rich sites are expected exonic
    assert 0 <= var <= exp


def test_matched_expectation_with_uninformative_gc_equals_the_unmatched_rate():
    rng = np.random.default_rng(1)
    null = pd.DataFrame({"gc250": rng.uniform(0.3, 0.6, 20000), "exon": rng.random(20000) < 0.09})
    sites = pd.Series(rng.uniform(0.3, 0.6, 500))
    exp, _ = edn.matched_expectation(null, sites, "gc250", n_bins=10)
    assert exp / 500 == pytest.approx(0.09, abs=0.02)


def test_label_positions_marks_exons_and_gene_bodies():
    df = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1", "chr2"], "pos": [150, 300, 5000, 100]})
    out = edn.label_positions(df, {"chr1": np.array([[100, 200]])}, {"chr1": np.array([[100, 500]])})
    assert out.exon.tolist() == [True, False, False, False]
    assert out.genic.tolist() == [True, True, False, False]


def test_sample_null_is_seeded_skips_gaps_and_covers_the_requested_chromosomes():
    rng_seq = np.random.default_rng(3)
    genome = {c: "".join(rng_seq.choice(list("ACGT"), 5000)) for c in ("chr1", "chr2")}
    gaps = {"chr1": np.array([[0, 2500]]), "chr2": np.zeros((0, 2), dtype=np.int64)}
    callable_bp = {"chr1": 2500, "chr2": 5000}
    a = edn.sample_null(genome, ["chr1", "chr2"], callable_bp, gaps, 400, np.random.default_rng(7))
    b = edn.sample_null(genome, ["chr1", "chr2"], callable_bp, gaps, 400, np.random.default_rng(7))
    assert a.equals(b)                                                    # same seed, same positions
    assert (a[a.chrom == "chr1"].pos - 1 >= 2500).all()                   # nothing drawn from the chr1 gap
    assert set(a.chrom) == {"chr1", "chr2"} and a.gc250.between(0, 1).all()
