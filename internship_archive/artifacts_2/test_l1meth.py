"""Tests for the L1 locus-methylation feasibility study."""

import numpy as np
import pytest

from l1meth import catalog as C
from l1meth import compare as K
from l1meth import panel as P
from l1meth import simulate as S


# --------------------------------------------------------------- catalogue

def test_revcomp():
    assert C.revcomp("ACGT") == "ACGT"
    assert C.revcomp("AACCGGTT") == "AACCGGTT"
    assert C.revcomp("ATGC") == "GCAT"


def test_utr_window_plus_strand_starts_at_element_start():
    # untruncated '+' element: promoter window is the first UTR_LEN bp
    assert C.utr_window(1000, 7000, "+", 0) == (1000, 1000 + C.UTR_LEN)


def test_utr_window_minus_strand_is_at_genomic_end():
    # consensus runs right-to-left for '-' strand: promoter is the LAST bp
    assert C.utr_window(1000, 7000, "-", 0) == (7000 - C.UTR_LEN, 7000)


def test_utr_window_respects_truncation():
    # a copy starting at consensus 400 has only 500 bp of promoter left
    lo, hi = C.utr_window(1000, 7000, "+", 400)
    assert hi - lo == C.UTR_LEN - 400


def test_utr_window_empty_when_promoter_absent():
    assert C.utr_window(1000, 7000, "+", C.UTR_LEN + 1) == (0, 0)


def test_promoter_cpg_order_flips_with_strand():
    seq = "A" * 1000 + "CG" + "T" * 300 + "CG" + "A" * 5000
    plus = C.promoter_cpgs(seq, 1000, 6000, "+", 0)
    assert plus == [1000, 1302]           # 5'->3' along a '+' element
    minus = C.promoter_cpgs(seq, 0, 1900, "-", 0)
    assert minus == sorted(minus, reverse=True)   # reversed for '-'


def test_promoter_cpg_capped_at_window():
    seq = "CG" * 2000
    hits = C.promoter_cpgs(seq, 0, 6000, "+", 0)
    assert len(hits) == C.N_PROMOTER_CPG


def _locus(**kw):
    base = dict(chrom="chr22", start=0, end=6000, strand="+", subfamily="L1HS",
                divergence=0.01, cons_start=0, flank_bp=500, cons_len=6200,
                cpg_genomic=[1, 2, 3, 4, 5])
    base.update(kw)
    return C.L1Locus(**base)


def test_consensus_length_uses_strand_correct_field():
    """Regression: UCSC swaps repStart/repLeft for '-' rows. Applying the '+'
    formula to a '-' row adds a coordinate to a length -- on chr22 that turned
    an 892 bp model into values up to 11,799 bp."""
    # '+' row: end 600, 292 bp of consensus remaining -> 892
    assert C.consensus_length(1, 600, -292, "+") == 892
    # '-' row with the same model: remainder in repStart, start in repLeft
    assert C.consensus_length(-292, 600, 5023, "-") == 892
    # the naive formula would have given 600 + 5023 = 5623 for that row
    assert 600 + abs(5023) != C.consensus_length(-292, 600, 5023, "-")


def test_subfamily_consensus_length_takes_the_mode_not_the_median():
    """Rows disagree even after the strand fix; the mode is robust, the median
    lands on a real value only by luck."""
    rows = ([dict(repName="L1ME4b", repStart=1, repEnd=600, repLeft=-292, strand="+")] * 5
            + [dict(repName="L1ME4b", repStart=1, repEnd=900, repLeft=-6982, strand="+")] * 2
            + [dict(repName="L1HS", repStart=1, repEnd=6000, repLeft=-155, strand="+")] * 3)
    out = C.subfamily_consensus_length(rows)
    assert out["L1ME4b"] == 892      # modal, not the mean (2418) or max (7882)
    assert out["L1HS"] == 6155


def test_short_consensus_subfamily_is_not_addressable():
    """L1ME4b's consensus model (892 bp) is shorter than the 900 bp promoter
    window, so a naive `cons_start < 900` test passes for its copies and floods
    the catalogue with degraded fragments that have no promoter at all."""
    bad = _locus(subfamily="L1ME4b", cons_len=892)
    assert not bad.full_length_consensus
    assert not bad.has_utr
    assert not bad.addressable(primate_only=False)


def test_non_primate_excluded_by_default_but_reachable():
    old = _locus(subfamily="L1MEc", cons_len=6400)
    assert old.has_utr                      # full-length consensus, retains 5' end
    assert not old.addressable()            # ...but not the human-type promoter
    assert old.addressable(primate_only=False)


def test_addressable_requires_flank_and_cpgs():
    assert _locus().addressable()
    assert not _locus(flank_bp=10).addressable()
    assert not _locus(cpg_genomic=[1]).addressable()
    assert not _locus(cons_start=5000).addressable()


def test_build_catalog_roundtrip():
    seq = "A" * 500 + "CG" * 40 + "T" * 6000
    rows = [dict(genoStart=500, genoEnd=6000, strand="+", repName="L1HS",
                 milliDiv=11, cons_start=0, flank_bp=400, cons_len=6200)]
    out = C.build_catalog(rows, seq)
    assert len(out) == 1 and out[0].subfamily == "L1HS"
    assert out[0].divergence == pytest.approx(0.011)
    assert out[0].n_cpg > 0


# ---------------------------------------------------------------- simulate

def test_simulate_shapes_and_labels():
    cfg = S.SimConfig(n_loci=20, n_normal=15, n_tumour=25, seed=1)
    X, y, sig = S.simulate(cfg)
    assert X.shape == (40, 20)
    assert y.sum() == 25 and len(y) == 40
    assert sig.shape == (20,)


def test_zero_signal_fraction_marks_no_loci():
    _, _, sig = S.simulate(S.SimConfig(signal_fraction=0.0, seed=2))
    assert sig.sum() == 0


def test_signal_loci_are_more_hypomethylated_in_tumour():
    cfg = S.SimConfig(n_loci=40, signal_fraction=0.25, locus_effect=30,
                      coverage=200, seed=3)
    X, y, sig = S.simulate(cfg)
    drop = X[y == 0].mean(axis=0) - X[y == 1].mean(axis=0)
    assert drop[sig].mean() > drop[~sig].mean()


def test_methylation_stays_in_range():
    X, _, _ = S.simulate(S.SimConfig(seed=4))
    assert X.min() >= 0.0 and X.max() <= 100.0


def test_lower_coverage_is_noisier():
    spreads = []
    for cov in (5, 100):
        X, y, _ = S.simulate(S.SimConfig(coverage=cov, signal_fraction=0.0, seed=5))
        spreads.append(X[y == 0].std(axis=0).mean())
    assert spreads[0] > spreads[1]


# ------------------------------------------------------------------ panel

def test_youngest_strategy_picks_lowest_age():
    age = np.array([0.9, 0.1, 0.5, 0.2])
    idx = P.select_panel(age, 2, strategy="youngest")
    assert set(idx) == {1, 3}


def test_unknown_strategy_rejected():
    with pytest.raises(ValueError):
        P.select_panel(np.zeros(5), 2, strategy="nope")


def test_age_informed_beats_blind_selection_only_under_enrichment():
    cfg = P.PoolConfig(signal_fraction=0.10, age_enrichment=8.0)
    blind = P.panel_hit_rate(cfg, 11, "random", n_rep=200)
    informed = P.panel_hit_rate(cfg, 11, "youngest", n_rep=200)
    assert informed["mean_signal_loci"] > blind["mean_signal_loci"]
    assert informed["p_any_signal"] > blind["p_any_signal"]


def test_hit_rate_rises_with_panel_size():
    cfg = P.PoolConfig(signal_fraction=0.05)
    small = P.panel_hit_rate(cfg, 5, "random", n_rep=200)["p_any_signal"]
    large = P.panel_hit_rate(cfg, 60, "random", n_rep=200)["p_any_signal"]
    assert large > small


def test_default_is_unenriched_so_selection_gains_nothing():
    """The default reflects the 450K evidence for carcinomas: young L1PAs are
    not preferentially hypomethylated there, so age-informed selection has no
    edge. Guards against the optimistic default silently coming back."""
    assert P.PoolConfig().age_enrichment == 1.0
    cfg = P.PoolConfig(signal_fraction=0.10)
    blind = P.panel_hit_rate(cfg, 20, "random", n_rep=300)["mean_signal_loci"]
    informed = P.panel_hit_rate(cfg, 20, "youngest", n_rep=300)["mean_signal_loci"]
    assert abs(blind - informed) < 0.6


def test_enrichment_is_what_drives_the_selection_advantage():
    """Only under an assumed enrichment -- defensible for germ-cell/embryonal
    contexts, not carcinomas -- does informed selection pull ahead."""
    germ = P.PoolConfig(signal_fraction=0.10, age_enrichment=8.0)
    blind = P.panel_hit_rate(germ, 11, "random", n_rep=300)["p_any_signal"]
    informed = P.panel_hit_rate(germ, 11, "youngest", n_rep=300)["p_any_signal"]
    assert informed > blind + 0.15


# ---------------------------------------------------------------- compare

def test_compare_returns_valid_aucs():
    X, y, _ = S.simulate(S.SimConfig(n_loci=15, seed=6))
    r = K.compare(X, y, seed=6)
    assert 0.0 <= r["auc_global"] <= 1.0
    assert 0.0 <= r["auc_locus"] <= 1.0
    assert r["n_loci"] == 15


def test_locus_loses_when_no_locus_specific_signal():
    """The published null, reproduced: with only a sample-level shift, the
    extra per-locus parameters cost accuracy rather than adding any."""
    g, lo = [], []
    for sd in range(8):
        X, y, _ = S.simulate(S.SimConfig(n_loci=35, signal_fraction=0.0, seed=sd))
        r = K.compare(X, y, seed=sd)
        g.append(r["auc_global"]); lo.append(r["auc_locus"])
    assert np.mean(lo) < np.mean(g)


def test_locus_wins_with_concentrated_signal():
    g, lo = [], []
    for sd in range(8):
        X, y, _ = S.simulate(S.SimConfig(n_loci=35, signal_fraction=0.20, seed=sd))
        r = K.compare(X, y, seed=sd)
        g.append(r["auc_global"]); lo.append(r["auc_locus"])
    assert np.mean(lo) > np.mean(g)
