"""Tests for the Phase 1 position-133 enrichment test.

The null model is the part worth testing hard. A window test that resamples in
one coordinate frame and counts in another looks perfectly reasonable and
produces a confidently wrong enrichment, so the frame is pinned explicitly
here.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "test_position_enrichment.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("test_position_enrichment", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


m = _load_module()


# --------------------------------------------------------------------------
# Window counting
# --------------------------------------------------------------------------


def test_observed_count_is_half_open():
    assert m.observed_count([127, 128, 138, 139], 128, 139) == 2
    assert m.observed_count([], 128, 139) == 0


def test_expected_uniform_matches_a_known_case():
    # One element of length 100, window [25, 35) has 10 admissible positions.
    assert m.expected_uniform([0], [100], 25, 35) == pytest.approx(0.10)


def test_expected_uniform_respects_element_length():
    # A short element has proportionally fewer chances to land in the window.
    short = m.expected_uniform([0], [100], 0, 10)
    long_ = m.expected_uniform([0], [300], 0, 10)
    assert short == pytest.approx(0.10)
    assert long_ == pytest.approx(10 / 300)
    assert long_ < short


def test_expected_uniform_clips_window_to_element():
    # Window wider than the element can never exceed the element's own length.
    assert m.expected_uniform([0], [20], 0, 100) == pytest.approx(1.0)


# --------------------------------------------------------------------------
# One definition of "window opportunity", shared by every null
# --------------------------------------------------------------------------


def test_window_opportunity_counts_the_window_inside_the_element():
    # Window [128, 139) is 11 wide and the element is 300 long: 11 positions.
    assert m.window_opportunity(offset=0, span=300, lo=128, hi=139) == 11


def test_window_opportunity_clips_to_the_element():
    # Element is 20 long, window 100 wide: the whole element, and no more.
    assert m.window_opportunity(offset=0, span=20, lo=0, hi=100) == 20
    # Element ends before the window starts: nothing at all.
    assert m.window_opportunity(offset=0, span=5, lo=128, hi=139) == 0


def test_window_opportunity_is_offset_by_the_call():
    # `offset` shifts the window relative to the element. An event 134 bases in
    # leaves only the first 5 of the window's 11 positions inside [0, span).
    assert m.window_opportunity(offset=134, span=300, lo=128, hi=139) == 5
    # Past the window entirely: the element has no positions left in it.
    assert m.window_opportunity(offset=200, span=300, lo=128, hi=139) == 0


def test_a_degenerate_window_or_element_gives_no_opportunity():
    assert m.window_opportunity(offset=0, span=300, lo=139, hi=139) == 0
    assert m.window_opportunity(offset=0, span=0, lo=0, hi=10) == 0


def test_a_weighted_null_at_unit_weights_is_exactly_the_uniform_null():
    """This is the invariant that pins Null B and Null C to Null A.

    They are Null A under a different weight, nothing else. Null B stopped being
    that when it dropped the clip, and published an expectation 2x Null A's
    while the median host mappable fraction was 1.0 -- a gap with nothing to do
    with mappability, reported as if it were the track's doing.
    """
    offsets = [0, 30, 128, 140, 290]
    spans = [300, 300, 150, 150, 12]
    for lo, hi in ((128, 139), (0, 11), (133, 134)):
        assert m.expected_weighted(
            [1.0] * len(offsets), offsets, spans, lo, hi
        ) == pytest.approx(m.expected_uniform(offsets, spans, lo, hi))


def test_a_weighted_null_cannot_exceed_one_expected_event_per_event():
    """The clip is what makes this true, and Null B's old formula broke it.

    Null B computed `(hi - lo) * f / span`, so a 100-wide window over a 20-long
    element contributed 5 expected events from one call.
    """
    offsets = [0] * 10
    spans = [20] * 10
    weighted = m.expected_weighted([1.0] * 10, offsets, spans, 0, 100)
    assert weighted == pytest.approx(10.0)
    # The unclipped formula this replaces, for the record.
    unclipped = sum((100 - 0) * 1.0 / 20 for _ in range(10))
    assert unclipped == pytest.approx(50.0)


def test_a_weighted_null_scales_linearly_in_the_weight():
    offsets, spans = [0, 50], [300, 300]
    unit = m.expected_weighted([1.0, 1.0], offsets, spans, 128, 139)
    half = m.expected_weighted([0.5, 0.5], offsets, spans, 128, 139)
    assert half == pytest.approx(unit * 0.5)


def test_expected_uniform_agrees_with_the_shared_helper_element_by_element():
    """The refactor that introduced `window_opportunity` must not have moved A."""
    offsets = [0, 7, 50, 128, 133, 139, 145, 290]
    spans = [300, 300, 12, 150, 150, 20, 20, 300]
    for lo, hi in ((128, 139), (0, 100), (133, 134), (10, 11)):
        assert m.expected_uniform(offsets, spans, lo, hi) == pytest.approx(
            sum(
                m.window_opportunity(o, s, lo, hi) / s
                for o, s in zip(offsets, spans)
                if s > 0
            )
        )


# --------------------------------------------------------------------------
# Coordinate frame: the bug this module exists to prevent
# --------------------------------------------------------------------------


def test_consensus_opportunity_prefers_aligned_span():
    assert m.consensus_opportunity(280, 300) == 280
    assert m.consensus_opportunity(None, 300) == 300
    assert m.consensus_opportunity(0, 300) == 300


def test_host_opportunity_uses_host_length():
    assert m.host_opportunity(0) == 1
    assert m.host_opportunity(301) == 301


def test_resample_stays_inside_each_element():
    """Every null draw must fall within its own element's extent."""
    import numpy as np

    rng = np.random.default_rng(7)
    opps = [50, 120, 300, 275]
    for _ in range(200):
        draws = m.resample_offsets([0, 0, 0, 0], opps, rng)
        for d, span in zip(draws, opps):
            assert 0 <= d < span


# --------------------------------------------------------------------------
# Null count distribution
# --------------------------------------------------------------------------


def test_null_window_counts_mean_matches_expectation():
    import numpy as np

    rng = np.random.default_rng(11)
    opps = [300] * 400
    counts = m._null_window_counts(opps, 128, 139, 4000, rng)
    assert len(counts) == 4000
    observed = m.expected_uniform([0] * 400, opps, 128, 139)
    assert counts.mean() == pytest.approx(observed, rel=0.10)


def test_null_window_counts_edge_and_centre_windows_are_equivalent():
    """An edge window and a central window of the same width are equal.

    Every element's own coordinate frame starts at 0, so a 10 bp window at
    [0,10) and one at [145,155) are both 10 of 300 positions. An earlier
    draft of this test assumed the edge window was scarcer, which would have
    been true only if positions were measured from a shared frame.
    """
    import numpy as np

    rng = np.random.default_rng(13)
    opps = [300] * 800
    edge = m._null_window_counts(opps, 0, 10, 4000, rng).mean()
    centre = m._null_window_counts(opps, 145, 155, 4000, rng).mean()
    assert edge == pytest.approx(centre, rel=0.05)
    assert edge == pytest.approx(800 * 10 / 300, rel=0.05)


def test_null_bin_maxima_shape_and_range():
    import numpy as np

    rng = np.random.default_rng(17)
    maxima = m._null_bin_maxima([300] * 200, 500, rng)
    assert maxima.shape == (500,)
    assert maxima.min() >= 1
    assert maxima.max() <= 200


# --------------------------------------------------------------------------
# Binning
# --------------------------------------------------------------------------


def test_bin_index_maps_offsets():
    assert m.bin_index(0) == 0
    assert m.bin_index(19) == 0
    assert m.bin_index(20) == 1
    assert m.bin_index(133) == 6
    assert m.bin_index(319) == 15


def test_bin_index_rejects_out_of_range_and_negative():
    assert m.bin_index(-1) is None
    assert m.bin_index(400) is None


# --------------------------------------------------------------------------
# Max-statistic FWER
# --------------------------------------------------------------------------


def test_max_statistic_fwer_detects_a_real_peak():
    import numpy as np

    rng = np.random.default_rng(19)
    offsets = [130] * 300 + [rng.integers(0, 320) for _ in range(300)]
    opps = [320] * len(offsets)
    res = m.max_statistic_fwer(offsets, opps, 2000, rng)
    assert res["n_bins"] == 16
    assert res["null_max_bin_p"] < 0.01


def test_max_statistic_fwer_gives_no_signal_for_uniform_data():
    import numpy as np

    rng = np.random.default_rng(23)
    rng2 = np.random.default_rng(29)
    offsets = [int(rng2.integers(0, 320)) for _ in range(600)]
    opps = [320] * 600
    res = m.max_statistic_fwer(offsets, opps, 2000, rng)
    assert res["null_max_bin_p"] > 0.05


# --------------------------------------------------------------------------
# Event loading
# --------------------------------------------------------------------------


def _row(**kw):
    base = {
        "chrom": "chr1",
        "pos": "1000",
        "nested_in_alu_host": "1",
        "host_len": "300",
        "host_start0": "900",
        "host_end0": "1200",
        "host_strand": "+",
        "host_name": "AluY",
        "consensus_mapping": "unambiguous",
        "consensus_offset": "133",
        "consensus_match_name": "AluY",
        "consensus_span_bp": "295",
        "host_offset_5p_0based": "138",
        "offset_drift_bp": "-5",
        "site_allele_count": "3",
        "insert_strand": "+",
        "perc_resolved": "99.0",
        "not_canonical": ".",
    }
    base.update(kw)
    return base


def test_load_nested_rows_keeps_near_full_hosts():
    out = m.load_nested_rows([_row()], require_unambiguous=True)
    assert len(out) == 1
    assert out[0]["offset"] == 133
    assert out[0]["consensus_span_bp"] == 295


def test_load_nested_rows_excludes_short_hosts():
    out = m.load_nested_rows([_row(host_len="150")])
    assert out == []


def test_load_nested_rows_excludes_ambiguous_by_default():
    out = m.load_nested_rows([_row(consensus_mapping="ambiguous_indel")])
    assert out == []
    out = m.load_nested_rows(
        [_row(consensus_mapping="ambiguous_indel")], require_unambiguous=False
    )
    assert len(out) == 1


def test_load_nested_rows_excludes_unnested():
    assert m.load_nested_rows([_row(nested_in_alu_host="0")]) == []


def test_load_nested_rows_can_use_host_coordinate():
    out = m.load_nested_rows([_row()], coordinate="host_offset_5p_0based")
    assert out[0]["offset"] == 138


def test_load_nested_rows_drops_unusable_coordinates():
    assert m.load_nested_rows([_row(consensus_offset="")]) == []


# --------------------------------------------------------------------------
# Column access: one discipline, so a missing column can never empty the cohort
# --------------------------------------------------------------------------


def test_the_required_and_optional_column_sets_are_disjoint():
    """A column cannot be load-bearing for membership and merely informative."""
    assert not set(m.REQUIRED_COHORT_COLUMNS) & set(m.OPTIONAL_COHORT_COLUMNS)


def test_the_declared_column_sets_cover_everything_load_nested_rows_reads():
    """The declarations are the contract, so nothing may be read outside them.

    A column read by `_cell` but absent from all three categories would be
    neither validated nor documented, and the partition would be a fiction. The
    third category is the coordinate columns, which are required only when the
    caller selects that frame -- `load_nested_rows` validates whichever one was
    asked for, so they belong to no fixed list.
    """
    declared = (
        set(m.REQUIRED_COHORT_COLUMNS)
        | set(m.OPTIONAL_COHORT_COLUMNS)
        | {"consensus_offset", "host_offset_5p_0based"}
    )
    assert {"chrom", "pos", "host_start0", "host_end0"} <= declared
    assert set(_row()) - declared == set(), (
        f"the fixture carries column(s) {sorted(set(_row()) - declared)} that "
        "`load_nested_rows` neither declares nor reads"
    )


@pytest.mark.parametrize("column", m.REQUIRED_COHORT_COLUMNS)
def test_a_missing_required_column_is_a_schema_error(column):
    """The gate columns used to fail open here, and that was the defect.

    `nested_in_alu_host` and the other three clauses were read with
    `row.get(...)`, so a column that went missing made every row fail the gate,
    emptied the cohort, and produced a report of zeros reading as "no enrichment
    here". The identity columns raised `KeyError` instead -- three behaviours in
    one function. All of them now raise `CohortSchemaError`.
    """
    row = {k: v for k, v in _row().items() if k != column}
    with pytest.raises(m.CohortSchemaError) as excinfo:
        m.load_nested_rows([row])
    message = str(excinfo.value)
    assert column in message
    assert "chrom" in message  # the header is reported, so the fix is obvious


def test_the_schema_error_is_distinguishable_from_an_empty_cohort():
    """An empty analysis set is a data outcome; a missing column is not.

    Collapsing the two is what let a schema break present as a clean null, so the
    distinction has to survive into the exception type.
    """
    assert m.load_nested_rows([_row(nested_in_alu_host="0")]) == []
    assert m.load_nested_rows([]) == []


@pytest.mark.parametrize("column", m.OPTIONAL_COHORT_COLUMNS)
def test_a_missing_optional_column_degrades_and_keeps_the_row(column):
    """Evidence columns are not membership, so absence must not exclude a row."""
    row = {k: v for k, v in _row().items() if k != column}
    out = m.load_nested_rows([row])
    assert len(out) == 1, f"dropping {column} changed the analysis-set size"


def test_a_degraded_optional_column_reads_as_empty_not_as_its_old_value():
    out = m.load_nested_rows([{k: v for k, v in _row().items() if k != "host_strand"}])
    assert out[0]["host_strand"] == ""


@pytest.mark.parametrize("column", ("pos", "host_start0", "host_end0"))
def test_a_present_but_unparseable_identity_is_loud(column):
    """Present-but-garbage must not quietly leave the denominator.

    These three are parsed with `int()` rather than `_to_int`, so an unusable
    value raises instead of excluding the row. A row that vanishes from the
    denominator for a reason nobody reported is the failure this whole change is
    about.
    """
    with pytest.raises(ValueError):
        m.load_nested_rows([_row(**{column: "not-a-number"})])


# --------------------------------------------------------------------------
# Estimand: sites and carriers are never merged
# --------------------------------------------------------------------------


def test_carrier_weighted_counts_reports_both():
    events = [
        {"offset": 133, "allele_count": 1},
        {"offset": 133, "allele_count": 400},
        {"offset": 10, "allele_count": 5},
    ]
    res = m.carrier_weighted_counts(events, 128, 139)
    assert res["n_sites_total"] == 3
    assert res["n_sites_in_window"] == 2
    assert res["n_carriers_total"] == 406
    assert res["n_carriers_in_window"] == 401
    # Two values in the window, so the median is the upper of the two.
    assert res["median_allele_count_in_window"] == 400


def test_leave_one_host_out_reports_influence():
    # Two events share host "a", one in host "b", one out of window.
    hosts = [("a",), ("a",), ("b",), ("b",)]
    offsets = [133, 133, 133, 10]
    res = m.leave_one_host_out(hosts, offsets, 128, 139)
    assert res["observed"] == 3
    assert res["n_hosts"] == 2
    assert res["max_events_in_one_host"] == 2
    assert res["observed_minus_worst_host"] == 1


def test_stratify_splits_by_orientation():
    events = [
        # sense: insert strand matches host strand
        {"insert_strand": "+", "host_strand": "+", "consensus_match_name": "AluY", "offset": 133},
        {"insert_strand": "-", "host_strand": "-", "consensus_match_name": "AluY", "offset": 133},
        # antisense: insert strand opposes host strand
        {"insert_strand": "-", "host_strand": "+", "consensus_match_name": "AluSx", "offset": 133},
        {"insert_strand": "+", "host_strand": "-", "consensus_match_name": "AluSx", "offset": 50},
    ]
    res = m.stratify(events, 128, 139)
    assert res["sense"]["n_sites"] == 2
    assert res["sense"]["n_in_window"] == 2
    assert res["antisense"]["n_sites"] == 2
    assert res["antisense"]["n_in_window"] == 1
    assert res["by_host_subfamily"]["AluY"]["n_sites"] == 2
    assert res["host_plus_strand"]["n_sites"] == 2
    assert res["host_minus_strand"]["n_sites"] == 2


def test_stratify_reports_unresolved_orientation_separately():
    events = [
        {"insert_strand": ".", "host_strand": "+", "consensus_match_name": "AluY", "offset": 133},
    ]
    res = m.stratify(events, 128, 139)
    # An unresolvable insert strand is folded into neither sense nor antisense.
    assert res["unresolved_orientation"]["n_sites"] == 1
    assert res["sense"]["n_sites"] == 0
    assert res["antisense"]["n_sites"] == 0
