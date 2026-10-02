"""Tests for the repo-local l1meth analysis package.

``l1meth/`` sits outside ``src/`` and so is not part of the installed
distribution, and it had no test coverage at all. These encode the invariants
the module docstrings assert -- strand symmetry, coordinate order, and the
signal-mass properties of the two generative models -- rather than golden
values, so they keep holding if the calibration constants are retuned.
"""

from __future__ import annotations

import numpy as np
import pytest

from l1meth.catalog import (
    MIN_CONSENSUS_LEN,
    N_PROMOTER_CPG,
    PRIMATE_PREFIXES,
    UTR_LEN,
    L1Locus,
    build_catalog,
    consensus_length,
    promoter_cpgs,
    revcomp,
    subfamily_consensus_length,
    utr_window,
)
from l1meth.panel import PoolConfig, make_pool, panel_hit_rate, select_panel
from l1meth.simulate import SimConfig, simulate


# ---------------------------------------------------------------------------
# revcomp
# ---------------------------------------------------------------------------


def test_revcomp_is_an_involution() -> None:
    for seq in ("ACGT", "GGCCAATT", "ACGTNN", ""):
        assert revcomp(revcomp(seq)) == seq


def test_revcomp_complements_each_base() -> None:
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAC") == "GTTT"
    assert revcomp("acgtn") == "nacgt"


# ---------------------------------------------------------------------------
# consensus_length -- the strand trap
# ---------------------------------------------------------------------------


def test_consensus_length_is_strand_symmetric() -> None:
    """The same element on '-' must not report a different model length.

    UCSC stores (repStart, repEnd, repLeft) in a strand-dependent order: for
    '-' rows repStart carries the negative remainder and repLeft the consensus
    start. Both layouts must yield one model length.
    """
    rep_start, rep_end, remainder = 5_000, 5_800, 4_200
    plus = consensus_length(rep_start, rep_end, -remainder, "+")
    minus = consensus_length(-remainder, rep_end, rep_start, "-")
    assert plus == minus == rep_end + remainder


def test_consensus_length_ignores_remainder_sign() -> None:
    assert consensus_length(100, 200, -50, "+") == 250
    assert consensus_length(-50, 200, 100, "-") == 250


def test_subfamily_consensus_length_takes_the_mode() -> None:
    rows = [
        {"repStart": 0, "repEnd": 600, "repLeft": -5200, "strand": "+", "repName": "L1PA2"},
        {"repStart": 0, "repEnd": 600, "repLeft": -5200, "strand": "+", "repName": "L1PA2"},
        {"repStart": 0, "repEnd": 599, "repLeft": -100, "strand": "+", "repName": "L1PA2"},
        {"repStart": 0, "repEnd": 892, "repLeft": 0, "strand": "-", "repName": "L1ME4b"},
    ]
    got = subfamily_consensus_length(rows)
    assert got["L1PA2"] == 5800
    assert got["L1ME4b"] == 892


# ---------------------------------------------------------------------------
# utr_window / promoter_cpgs -- coordinate order
# ---------------------------------------------------------------------------


def test_utr_window_plus_strand_runs_from_element_start() -> None:
    w0, w1 = utr_window(1_000, 7_000, "+", 0)
    assert (w0, w1) == (1_000, 1_000 + UTR_LEN)


def test_utr_window_minus_strand_runs_to_element_end() -> None:
    w0, w1 = utr_window(1_000, 7_000, "-", 0)
    assert (w1, w0) == (7_000, 7_000 - UTR_LEN)


def test_utr_window_is_empty_when_no_utr_retained() -> None:
    assert utr_window(1_000, 7_000, "+", UTR_LEN) == (0, 0)
    assert utr_window(1_000, 7_000, "-", UTR_LEN + 5) == (0, 0)


def test_utr_window_is_truncated_by_cons_start() -> None:
    span = UTR_LEN - 100
    assert utr_window(1_000, 7_000, "+", 100) == (1_000, 1_000 + span)
    assert utr_window(1_000, 7_000, "-", 100) == (7_000 - span, 7_000)


def test_utr_window_never_leaves_the_element() -> None:
    """A short element must be clamped, not allowed to run off either end."""
    for strand in ("+", "-"):
        for length in (50, 200, 900, 901):
            w0, w1 = utr_window(1_000, 1_000 + length, strand, 0)
            assert 1_000 <= w0 < w1 <= 1_000 + length


def _chr_with(motif_at: int, motif: str, length: int = 4_000) -> str:
    """A chromosome-scale string with ``motif`` written at ``motif_at``.

    ``promoter_cpgs`` slices ``seq`` by genomic coordinate, so the sequence
    must be as long as the coordinates it is asked about.
    """
    return "T" * motif_at + motif + "T" * (length - motif_at - len(motif))


def test_promoter_cpgs_are_ordered_five_to_three_along_the_element() -> None:
    """Each strand reports the 15 CpGs nearest *its own* 5' end, in element order.

    On a long element the '+' and '-' 5'UTR windows sit at opposite ends of the
    copy, so they do not contain the same CpGs -- only the traversal direction
    is shared.
    """
    start, end = 1_000, 2_000
    seq = _chr_with(start, "CG" * 500, length=4_000)
    plus = promoter_cpgs(seq, start, end, "+", 0)
    minus = promoter_cpgs(seq, start, end, "-", 0)

    assert len(plus) == len(minus) == N_PROMOTER_CPG
    # '+' : ascending genomic, anchored at the element's left edge.
    assert plus == [start + 2 * i for i in range(N_PROMOTER_CPG)]
    # '-' : descending genomic, anchored at the element's right edge.
    assert minus == [end - 2 - 2 * i for i in range(N_PROMOTER_CPG)]
    assert plus == sorted(plus)
    assert minus == sorted(minus, reverse=True)
    # Each window is the 900 bp 5'UTR measured from that strand's 5' end.
    assert min(plus) >= start and max(plus) < start + UTR_LEN
    assert max(minus) < end and min(minus) > end - UTR_LEN


def test_promoter_cpgs_caps_at_the_promoter_window() -> None:
    seq = "CG" * 600
    got = promoter_cpgs(seq, 0, 1_200, "+", 0)
    assert len(got) == N_PROMOTER_CPG


def test_promoter_cpgs_empty_when_no_window() -> None:
    assert promoter_cpgs("C" * 1_000, 0, 1_000, "+", UTR_LEN) == []


def test_promoter_cpgs_finds_cpg_on_either_strand() -> None:
    """CpG is palindromic: a '-' copy's CpGs are CG on the forward strand too."""
    start, end = 2_000, 3_000
    seq = _chr_with(2_100, "CG", length=3_500)
    got = promoter_cpgs(seq, start, end, "-", 0)
    assert got == [2_100]  # ordering flips, the position does not


# ---------------------------------------------------------------------------
# L1Locus gates
# ---------------------------------------------------------------------------


def _locus(**kw) -> L1Locus:
    base = dict(
        chrom="chr22", start=1_000, end=7_000, strand="+", subfamily="L1PA2",
        divergence=0.02, cons_start=0, flank_bp=300,
        cons_len=6_000, cpg_genomic=[1_000, 1_010, 1_020],
    )
    base.update(kw)
    return L1Locus(**base)


def test_has_utr_requires_full_length_consensus() -> None:
    """A fragment-derived subfamily has no promoter to speak of."""
    assert _locus(cons_len=MIN_CONSENSUS_LEN).has_utr
    assert not _locus(cons_len=899).has_utr


def test_has_utr_requires_the_copy_to_reach_the_promoter() -> None:
    assert _locus(cons_start=0).has_utr
    assert _locus(cons_start=UTR_LEN - 1).has_utr
    assert not _locus(cons_start=UTR_LEN).has_utr
    assert not _locus(cons_start=-1).has_utr


def test_primate_lineage_prefixes() -> None:
    for sub in PRIMATE_PREFIXES:
        assert _locus(subfamily=sub).primate_lineage
    assert not _locus(subfamily="L1ME4b").primate_lineage


def test_addressable_requires_primate_flank_and_cpg() -> None:
    assert _locus().addressable()
    assert not _locus(flank_bp=10).addressable()
    assert not _locus(cpg_genomic=[1_000, 1_010]).addressable()
    # ...but those relax under their own knobs.
    assert _locus(flank_bp=10).addressable(min_flank=5)
    assert _locus(cpg_genomic=[1_000, 1_010]).addressable(min_cpg=2)


def test_addressable_separates_the_lineage_gate_from_the_length_gate() -> None:
    """The two independent reasons a copy is unusable, isolated.

    Dropping the primate filter admits an ancient-but-full-length L1M copy; it
    must not admit L1ME4b, whose 892 bp consensus is shorter than the 900 bp
    promoter window itself, so "the first 900 bp" is the whole model.
    """
    old_full_length = _locus(subfamily="L1M2a", cons_len=6_000)
    assert not old_full_length.addressable()               # lineage gate
    assert old_full_length.addressable(primate_only=False)  # fine on length

    truncated = _locus(subfamily="L1ME4b", cons_len=892)  # its real model length
    assert not truncated.addressable()
    assert not truncated.addressable(primate_only=False)   # length gate still holds


def test_addressable_never_bypasses_the_utr_gate() -> None:
    """Relaxing every downstream knob must not resurrect a promoterless copy."""
    assert not _locus(cons_len=500, flank_bp=10_000, cpg_genomic=list(range(500))).addressable(
        min_flank=0, min_cpg=0, primate_only=False
    )


def test_build_catalog_populates_cpgs_only_for_utr_copies() -> None:
    # CpG-dense reference so the promoter's 15-CpG window is actually populated.
    seq = "CG" * 450 + "A" * 6_550
    rows = [
        {"genoStart": 0, "genoEnd": 6_000, "strand": "+", "repName": "L1PA2",
         "milliDiv": 20, "cons_start": 0, "flank_bp": 300, "cons_len": 6_000},
        {"genoStart": 20_000, "genoEnd": 26_000, "strand": "+", "repName": "L1ME4b",
         "milliDiv": 700, "cons_start": 0, "flank_bp": 300, "cons_len": 892},
    ]
    loci = build_catalog(rows, seq, chrom="chr22")
    assert len(loci) == 2
    assert loci[0].n_cpg == N_PROMOTER_CPG
    assert loci[0].addressable()
    # A fragment subfamily has no promoter to query, so it is never queried.
    assert loci[1].cpg_genomic == []
    assert not loci[1].has_utr
    assert not loci[1].addressable()
    assert loci[0].divergence == pytest.approx(0.020)
    assert loci[0].locus_id == "chr22:0-6000:+:L1PA2"


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


def test_simulate_shapes_and_ranges() -> None:
    cfg = SimConfig(n_loci=35, n_normal=20, n_tumour=20, seed=1)
    X, y, is_signal = simulate(cfg)
    assert X.shape == (40, 35)
    assert y.shape == (40,)
    assert is_signal.shape == (35,)
    assert set(np.unique(y)) == {0, 1}
    assert (y[: cfg.n_normal] == 0).all()
    assert (y[cfg.n_normal:] == 1).all()
    assert ((X >= 0) & (X <= 100)).all()


def test_simulate_signal_count_matches_the_requested_fraction() -> None:
    cfg = SimConfig(n_loci=100, signal_fraction=0.15, seed=3)
    _X, _y, is_signal = simulate(cfg)
    assert int(is_signal.sum()) == 15


def test_simulate_is_deterministic_for_a_seed() -> None:
    cfg = SimConfig(seed=7)
    X1, y1, s1 = simulate(cfg)
    X2, y2, s2 = simulate(cfg)
    assert np.array_equal(X1, X2)
    assert np.array_equal(y1, y2)
    assert np.array_equal(s1, s2)


def test_simulate_tumours_are_hypomethylated_versus_normals() -> None:
    """The whole point of the model: tumour < normal."""
    cfg = SimConfig(n_loci=200, n_normal=80, n_tumour=80, seed=5)
    X, y, _ = simulate(cfg)
    assert X[y == 1].mean() < X[y == 0].mean()


def test_simulate_signal_loci_carry_the_locus_specific_drop() -> None:
    cfg = SimConfig(
        n_loci=200, n_normal=200, n_tumour=200,
        signal_fraction=0.5, locus_effect=25.0, seed=11,
    )
    X, y, is_signal = simulate(cfg)
    drop = X[y == 0].mean(axis=0) - X[y == 1].mean(axis=0)
    assert drop[is_signal].mean() > drop[~is_signal].mean() + 10.0


# ---------------------------------------------------------------------------
# panel
# ---------------------------------------------------------------------------


def test_select_panel_youngest_picks_the_lowest_age() -> None:
    age = np.array([0.9, 0.1, 0.5, 0.2])
    idx = select_panel(age, 2, strategy="youngest")
    assert sorted(idx) == [1, 3]
    assert np.all(np.diff(age[idx]) >= 0)


def test_select_panel_random_is_distinct_and_sized() -> None:
    age = np.linspace(0, 1, 50)
    idx = select_panel(age, 10, strategy="random", seed=0)
    assert len(idx) == 10
    assert len(set(idx.tolist())) == 10
    assert idx.max() < 50


def test_select_panel_rejects_an_unknown_strategy() -> None:
    with pytest.raises(ValueError, match="unknown strategy"):
        select_panel(np.linspace(0, 1, 5), 2, strategy="oldest")


def test_make_pool_age_enrichment_one_carries_no_age_bias() -> None:
    """At age_enrichment=1.0 the evidence-supported default, signal is
    uncorrelated with age."""
    age, sig = make_pool(PoolConfig(n_pool=4000, signal_fraction=0.10,
                                   age_enrichment=1.0, seed=0))
    young, old = sig[age < 0.5], sig[age >= 0.5]
    assert abs(young.mean() - old.mean()) < 0.03


def test_make_pool_enrichment_raises_signal_rate_in_young_loci() -> None:
    age, sig = make_pool(PoolConfig(n_pool=4000, signal_fraction=0.10,
                                   age_enrichment=8.0, seed=0))
    assert sig[age < 0.5].mean() > 3 * sig[age >= 0.5].mean()


def test_make_pool_signal_weight_decays_with_age() -> None:
    """Per-locus probability must fall monotonically with age at enrichment>1."""
    cfg = PoolConfig(n_pool=30_000, signal_fraction=0.5, age_enrichment=8.0, seed=0)
    age, sig = make_pool(cfg)
    bins = np.array_split(np.argsort(age), 4)
    rates = [sig[b].mean() for b in bins]
    assert rates == sorted(rates, reverse=True), f"not age-decaying: {rates}"


def test_make_pool_preserves_the_requested_signal_mass() -> None:
    """Enrichment reweights *which* loci carry signal, not how many overall."""
    for enrichment in (1.0, 4.0, 8.0):
        _age, sig = make_pool(PoolConfig(n_pool=20_000, signal_fraction=0.10,
                                        age_enrichment=enrichment, seed=2))
        assert abs(sig.mean() - 0.10) < 0.015, f"enrichment={enrichment}"


def test_panel_hit_rate_agrees_with_the_signal_mass_it_simulates() -> None:
    """The reported hit rate must match the pool's actual signal density."""
    cfg = PoolConfig(n_pool=2_000, signal_fraction=0.10, age_enrichment=1.0)
    res = panel_hit_rate(cfg, k=50, strategy="random", n_rep=120)
    density = cfg.signal_fraction
    expected = 1.0 - (1.0 - density) ** 50
    assert res["p_any_signal"] == pytest.approx(expected, abs=0.06)
    assert res["mean_signal_loci"] == pytest.approx(density * 50, rel=0.15)
    assert res["k"] == 50
    assert res["strategy"] == "random"
