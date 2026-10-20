"""Tests for the Phase 2a boundary-perturbation falsification.

The point of this module is to avoid two symmetric mistakes: declaring every
terminal peak an artefact, and declaring none of them one. The tests below pin
both the mechanical-baseline logic that makes the perturbation fair, and the
fact that a genuinely artefactual region would still be struck.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "test_boundary_artifact.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("test_boundary_artifact", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


m = _load_module()


def _event(offset, host_len=300, strand="+", tsd_len=15, tsd_seq="ACGTTGCAACGTTGCA", poly_a=20):
    return {
        "offset": offset,
        "host_len": host_len,
        "host_strand": strand,
        "insert_strand": "+",
        "tsd_len": tsd_len,
        "tsd_seq": tsd_seq,
        "polyA_len": poly_a,
        "conformation": "FOR+POLYA",
        "perc_resolved": 99.0,
        "not_canonical": ".",
    }


# --------------------------------------------------------------------------
# Boundary perturbation
# --------------------------------------------------------------------------


def _row(pos=1000, start0=900, end0=1200, strand="+", host_len="300"):
    return {
        "pos": str(pos),
        "host_start0": str(start0),
        "host_end0": str(end0),
        "host_strand": strand,
        "host_len": host_len,
        "nested_in_alu_host": "1",
        "consensus_mapping": "unambiguous",
        "consensus_offset": "133",
    }


def test_shifted_offset_plus_strand():
    # pos 1000 -> 0-based 999, host starts at 900, so offset 99.
    assert m.shifted_offset(_row(), 0, 0) == 99
    # Shrinking the left edge by 5 pushes the host start to 905.
    assert m.shifted_offset(_row(), 5, -5) == 94


def test_shifted_offset_minus_strand_swaps_the_ends():
    """A - strand host's genomic-left end is its 3' end, so the shifts swap."""
    row = _row(pos=1000, start0=900, end0=1200, strand="-")
    # 0-based 999; for a - strand host offset = end-1-pos0 = 1200-1-999 = 200
    assert m.shifted_offset(row, 0, 0) == 200
    # shift left by 5 and right by -5: end becomes 1195
    assert m.shifted_offset(row, 5, -5) == 1195 - 1 - 999


def test_shifted_offset_returns_none_when_event_falls_outside():
    # Shrink by 200: the host no longer contains the breakpoint.
    assert m.shifted_offset(_row(), 200, -200) is None


def test_shifted_offset_handles_missing_fields():
    row = _row()
    row["host_start0"] = ""
    assert m.shifted_offset(row, 0, 0) is None


# --------------------------------------------------------------------------
# The mechanical baseline: the reason the perturbation is fair
# --------------------------------------------------------------------------


def test_mechanical_shift_counts_edge_events_that_must_leave():
    """Events within the shift of the 5' end are pushed outside a + host."""
    events = [_event(2), _event(150), _event(298)]
    # shift 5: only offset 2 is inside [0,5) on a + strand host
    assert m.mechanical_shift_expectation(events, (0, 20), 5) == 1


def test_mechanical_shift_respects_host_strand():
    """On a - strand host the same genomic end is the element's 3' end."""
    plus = [_event(2, strand="+")]
    minus = [_event(2, strand="-")]
    assert m.mechanical_shift_expectation(plus, (0, 20), 5) == 1
    assert m.mechanical_shift_expectation(minus, (0, 20), 5) == 0


def test_mechanical_shift_ignores_events_outside_the_region():
    events = [_event(150)]
    assert m.mechanical_shift_expectation(events, (0, 20), 5) == 0


def test_stability_uses_the_mechanical_tolerance():
    """A large swing is only a strike once it exceeds what geometry forces.

    An edge region must not be struck for moving by the amount the
    perturbation is guaranteed to move it.
    """
    # 100 events at the 5' end, none mechanical, so any movement is a strike.
    stable_perturb = {
        "as_reported": {"five_prime": 100, "three_prime": 0, "linker": 0, "n_sites": 0, "profile": []},
        "shrink_5bp_each_side": {"five_prime": 100, "three_prime": 0, "linker": 0, "n_sites": 0, "profile": []},
        "shrink_2bp_each_side": {"five_prime": 100, "three_prime": 0, "linker": 0, "n_sites": 0, "profile": []},
        "expand_2bp_each_side": {"five_prime": 100, "three_prime": 0, "linker": 0, "n_sites": 0, "profile": []},
    }
    events = [_event(10 + i) for i in range(100)]
    st = m.peak_stability(stable_perturb, events)
    # Mechanical expectation for a 5bp shift with offsets 10..109 is zero.
    assert st["five_prime"]["mechanical_expectation"] == 0
    assert st["five_prime"]["stable"] is True

    # Now make it swing beyond the mechanical allowance.
    swinging = {
        k: dict(v) for k, v in stable_perturb.items()
    }
    swinging["shrink_5bp_each_side"]["five_prime"] = 60
    st2 = m.peak_stability(swinging, events)
    assert st2["five_prime"]["observed_swing"] == 40
    assert st2["five_prime"]["mechanical_expectation"] == 0
    assert st2["five_prime"]["stable"] is False


# --------------------------------------------------------------------------
# Junction-evidence triage
# --------------------------------------------------------------------------


def test_sentinel_tsd_is_its_own_class():
    assert m.classify_junction_evidence(_event(10, tsd_len=41)) == "sentinel_tsd"
    assert m.classify_junction_evidence(_event(10, tsd_len=82)) == "sentinel_tsd"


def test_absent_tsd_is_its_own_class():
    assert m.classify_junction_evidence(_event(10, tsd_len=0, tsd_seq="")) == "tsd_absent"


def test_poly_at_tsd_is_flagged():
    assert m.classify_junction_evidence(_event(10, tsd_seq="AAAAAAAAAA")) == "poly_at_tsd"


def test_tsd_outside_de_novo_range_is_flagged_not_dropped():
    """Short or long TSDs are described, never hard-thresholded away."""
    assert m.classify_junction_evidence(_event(10, tsd_len=3, tsd_seq="ACG")) == "tsd_outside_de_novo_range"
    assert m.classify_junction_evidence(_event(10, tsd_len=60, tsd_seq="ACG")) == "tsd_outside_de_novo_range"


def test_tsd_in_range_without_polya_is_its_own_class():
    assert m.classify_junction_evidence(_event(10, poly_a=0)) == "tsd_in_range_no_polya"


def test_full_triad_is_informative():
    assert m.classify_junction_evidence(_event(10)) == "triad_informative"


def test_evidence_profile_counts_every_class():
    events = [
        _event(5, tsd_len=41),        # sentinel
        _event(6, tsd_len=0, tsd_seq=""),  # absent
        _event(7),                    # informative
        _event(300),                  # outside region
    ]
    prof = m.evidence_profile(events, (0, 20), (41, 82))
    assert prof["n"] == 3
    assert prof["class_counts"]["sentinel_tsd"] == 1
    assert prof["class_counts"]["tsd_absent"] == 1
    assert prof["class_counts"]["triad_informative"] == 1
    assert prof["informative_fraction"] == 1 / 3


# --------------------------------------------------------------------------
# Verdicts: must be capable of striking
# --------------------------------------------------------------------------


def _stability(stable=True):
    return {"stable": stable, "as_reported": 100, "min": 100, "max": 100,
            "observed_swing": 0, "mechanical_expectation": 0, "tolerance": 10,
            "relative_swing": 0.0}


def test_verdict_strikes_an_unstable_region():
    v = m.verdict_for_region(
        "three_prime", _stability(stable=False),
        {"enrichment": 2.0}, {"informative_fraction": 0.5},
        {"informative_fraction": 0.5},
    )
    assert v["verdict"] == "struck_as_artefact"
    assert "region_count_moves_under_boundary_perturbation" in v["strike_conditions_met"]


def test_verdict_strikes_a_region_at_chance():
    v = m.verdict_for_region(
        "three_prime", _stability(stable=True),
        {"enrichment": 0.9}, {"informative_fraction": 0.5},
        {"informative_fraction": 0.5},
    )
    assert v["verdict"] == "struck_as_artefact"
    assert "region_enrichment_not_above_uniform_after_its_own_denominator" in v["strike_conditions_met"]


def test_verdict_strikes_a_region_with_worse_evidence_than_control():
    v = m.verdict_for_region(
        "three_prime", _stability(stable=True),
        {"enrichment": 2.0}, {"informative_fraction": 0.1},
        {"informative_fraction": 0.5},
    )
    assert v["verdict"] == "struck_as_artefact"
    assert "region_junction_evidence_worse_than_the_mid_host_control" in v["strike_conditions_met"]


def test_verdict_survives_a_clean_region():
    v = m.verdict_for_region(
        "linker", _stability(stable=True),
        {"enrichment": 2.0}, {"informative_fraction": 0.5},
        {"informative_fraction": 0.4},
    )
    assert v["verdict"] == "survives_artefact_tests"
    assert v["strike_conditions_met"] == []


def test_verdict_does_not_favour_either_terminal_end():
    """The same inputs must give the same verdict for 5' and 3'."""
    args = (_stability(stable=True), {"enrichment": 1.5},
            {"informative_fraction": 0.4}, {"informative_fraction": 0.4})
    five = m.verdict_for_region("five_prime", *args)
    three = m.verdict_for_region("three_prime", *args)
    assert five["verdict"] == three["verdict"]
    assert five["strike_conditions_met"] == three["strike_conditions_met"]


# --------------------------------------------------------------------------
# Enrichment denominators
# --------------------------------------------------------------------------


def test_region_enrichment_uses_region_width():
    # 30 hosts, a 60 bp region in a 300 bp element -> 6 expected.
    res = m.region_enrichment(30, (240, 300), 30)
    assert res["expected_uniform"] == 6.0
    assert res["enrichment"] == 5.0


def test_region_count_is_half_open():
    assert m.region_count([0, 19, 20], (0, 20)) == 2
