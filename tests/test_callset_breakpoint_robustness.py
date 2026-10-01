"""Regression tests for rtm_breakpoint / _as_int robustness.

Two defects fixed together:

1. ``int(float("inf"))`` and ``float(10**400)`` raise ``OverflowError``, which
   is neither ``TypeError`` nor ``ValueError``. The handlers only caught the
   latter two, so an ``inf`` / ``1e400`` cell in ``window_start`` or
   ``window_end`` crashed with an unhandled ``OverflowError``. ``nan`` was
   always safe because it raises ``ValueError``.

2. A window-midpoint fallback after the column loop was unreachable: the loop
   already returns for any ``window_start`` that parses positive, and
   ``_as_int`` applies the identical conversion, so its ``start > 0`` guard
   could never hold.
"""

from __future__ import annotations

import itertools

import pytest

from retro_miner.callset_eval import _as_int, rtm_breakpoint

MISSING = object()

INFINITE = ["inf", "-inf", "1e400", "-1e400", 10**400, -(10**400)]
NON_FINITE = ["nan", "-nan"]
NON_NUMERIC = ["", "  ", "abc", "12abc", None, MISSING]


@pytest.mark.parametrize("value", INFINITE)
def test_as_int_handles_infinity(value):
    """infinity must fall back to the default, not raise OverflowError."""
    assert _as_int(value, -1) == -1


def test_as_int_handles_huge_int():
    assert _as_int(10**400, -1) == -1
    assert _as_int(-(10**400), -1) == -1


@pytest.mark.parametrize("value", NON_FINITE)
def test_as_int_still_handles_nan(value):
    assert _as_int(value, -1) == -1


@pytest.mark.parametrize("value", INFINITE)
def test_rtm_breakpoint_does_not_crash_on_infinite_window_start(value):
    assert rtm_breakpoint({"window_start": value}) == 0


@pytest.mark.parametrize("value", INFINITE)
def test_rtm_breakpoint_does_not_crash_on_infinite_window_end(value):
    assert rtm_breakpoint({"window_start": 0, "window_end": value}) == 0


def test_rtm_breakpoint_infinite_end_with_real_start_is_ignored():
    """An infinite window_end must not influence the returned breakpoint."""
    assert rtm_breakpoint({"window_start": 100, "window_end": "1e400"}) == 100


def test_rtm_breakpoint_prefers_consensus_then_insertion_then_window_start():
    row = {
        "consensus_insertion_breakpoint_pos": 10,
        "insertion_breakpoint_pos": 20,
        "window_start": 30,
    }
    assert rtm_breakpoint(row) == 10
    row.pop("consensus_insertion_breakpoint_pos")
    assert rtm_breakpoint(row) == 20
    row.pop("insertion_breakpoint_pos")
    assert rtm_breakpoint(row) == 30


@pytest.mark.parametrize("start", [0, -5, "0", None, "", "nan", "abc"])
def test_rtm_breakpoint_returns_zero_when_no_usable_field(start):
    assert rtm_breakpoint({"window_start": start, "window_end": 999}) == 0


def test_rtm_breakpoint_truncates_fractional_values():
    assert rtm_breakpoint({"window_start": "100.9"}) == 100


def test_rtm_breakpoint_never_uses_window_midpoint():
    """The deleted midpoint fallback was unreachable, so a valid start is
    returned as-is and never averaged with window_end."""
    assert rtm_breakpoint({"window_start": 100, "window_end": 200}) == 100
    assert rtm_breakpoint({"window_start": 100, "window_end": 100_000}) == 100


def test_matches_loop_only_reference_over_input_space():
    """Behaviour must equal 'first positive of the three columns, else 0'."""
    values = [
        MISSING, None, "", "  ", "abc", 0, "0", -5, "1", "1", 42, "3.9",
        "1e3", "  77  ", "nan", True, False, "inf", "1e400",
    ]
    cols = [
        "consensus_insertion_breakpoint_pos",
        "insertion_breakpoint_pos",
        "window_start",
        "window_end",
    ]

    def reference(row):
        for col in cols[:3]:
            raw = row.get(col)
            try:
                pos = int(float(raw))
            except (TypeError, ValueError, OverflowError):
                continue
            if pos > 0:
                return pos
        return 0

    for combo in itertools.product(values, repeat=4):
        row = {c: v for c, v in zip(cols, combo) if v is not MISSING}
        assert rtm_breakpoint(dict(row)) == reference(row), row
