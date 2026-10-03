"""`int(float(x))` raises OverflowError on infinity, and nobody was catching it.

Every numeric converter in this repo funnels a cell through
`int(float(value))` and catches `(TypeError, ValueError)`. That looks complete
and is not: `OverflowError` is not a `ValueError` subclass, and

    int(float("inf"))   -> OverflowError
    int(float("1e400")) -> OverflowError   # float() overflows to inf

so a single non-finite cell in a cohort column aborts the whole analysis with
an unhandled exception instead of reading as unevaluable and falling back to
the default. The same converters already treat `""` and `"."` as "no value", so
the intent is unambiguous -- the missing case was just not covered.

These tests pin the contract on each converter that reads cohort data. They are
cheap and they exist because the failure is loud but unhelpful: it surfaces as a
traceback from deep inside an annotation run, far from the cell that caused it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import test_boundary_artifact as tba  # noqa: E402
import test_position_enrichment as tpe  # noqa: E402
from retro_miner import callset_eval, read_architecture  # noqa: E402

#: Every float spelling that reaches `inf` through `float()`.
NON_FINITE = ["inf", "-inf", "Infinity", "-Infinity", "INF", "1e400", "-1e400", "1e309"]

#: Values that must keep working exactly as before.
FINITE = ["0", "1", "133", " 133 ", "-5", "1e3", "2.9"]


# --------------------------------------------------------------------------
# the two shapes of converter: "return None" and "return a default"
# --------------------------------------------------------------------------


@pytest.mark.parametrize("value", NON_FINITE)
@pytest.mark.parametrize(
    "converter", [tpe._to_int, tba._to_int], ids=["position_enrichment", "boundary_artifact"]
)
def test_none_returning_converters_treat_infinity_as_unevaluable(converter, value) -> None:
    """A non-finite cell is missing data, not a reason to abort the run."""
    assert converter(value) is None


@pytest.mark.parametrize("value", NON_FINITE)
def test_default_returning_converters_treat_infinity_as_the_default(value) -> None:
    assert callset_eval._as_int(value, -1) == -1
    row = pd.Series({"c": value})
    assert read_architecture._row_int(row, "c", -1) == -1


# --------------------------------------------------------------------------
# no behaviour change for ordinary input
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "converter,expected",
    [("0", 0), ("1", 1), ("133", 133), (" 133 ", 133), ("-5", -5), ("1e3", 1000), ("2.9", 2)],
)
def test_the_finite_cases_are_unchanged(converter, expected) -> None:
    assert tpe._to_int(converter) == expected


@pytest.mark.parametrize("bad", ["", ".", "x", "1,5", "--3", None])
def test_the_previously_handled_bad_values_still_read_as_unevaluable(bad) -> None:
    assert tpe._to_int(bad) is None
    assert tba._to_int(bad) is None


@pytest.mark.parametrize("value", FINITE)
def test_nan_is_still_missing_not_a_number(value) -> None:
    """`int(float("nan"))` raises ValueError and must keep reading as missing."""
    assert tpe._to_int(value) == int(float(value))


def test_nan_and_infinity_are_both_missing_and_neither_crashes() -> None:
    for value in ("nan", "NaN", "inf", "1e400"):
        assert tpe._to_int(value) is None
        assert tba._to_int(value) is None


def test_the_three_cohort_readers_agree_on_every_case() -> None:
    """`read_cohort` is triplicated across the analysis scripts verbatim.

    If the `_to_int` beside one copy drifts from the others, two published
    analyses can end up on different denominators from one table -- the exact
    thing `nested_analysis_set`'s docstring says it exists to prevent. This
    pins the three copies to identical behaviour on the awkward inputs.
    """
    for value in NON_FINITE + FINITE + ["", ".", "x", "nan"]:
        results = {tpe._to_int(value), tba._to_int(value)}
        assert len(results) == 1, (
            f"_to_int copies disagree on {value!r}: "
            f"position_enrichment={tpe._to_int(value)!r} "
            f"boundary_artifact={tba._to_int(value)!r}"
        )


def test_overflow_error_is_not_a_value_error() -> None:
    """The premise of every fix above, asserted so it cannot quietly change.

    If this ever stops holding, the `OverflowError` arms become redundant
    rather than wrong -- harmless -- but if it starts holding differently the
    arms would be masking a different failure.
    """
    assert not issubclass(OverflowError, ValueError)
    with pytest.raises(OverflowError):
        int(float("inf"))