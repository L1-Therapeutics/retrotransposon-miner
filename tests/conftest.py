"""Test-harness guards for the optional ``bedtools`` binary.

The nested-insertion annotation feature (``_annotate_nested_retrotransposon``)
shells out to ``bedtools``, so its tests are marked
``skipif(shutil.which("bedtools") is None)``. That gate is invisible in a
passing run: on a machine without bedtools the entire feature goes untested and
the suite still reports green. Two mechanisms make that visible:

* **Opt-in strict mode.** With ``RTM_REQUIRE_BEDTOOLS=1`` the run aborts at
  collection time (pytest exit code 4) when bedtools is missing, instead of
  skipping. Unset, behaviour is byte-for-byte what it was before.
* **End-of-run notice.** Whenever one or more tests skipped with a reason
  mentioning bedtools, the terminal summary prints how many skipped and that the
  nested-insertion annotation feature is unverified by this run.

Both read skip reasons off the report objects pytest already produces, so they
keep working for gated tests added later without editing this file.
"""

from __future__ import annotations

import os
import shutil
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from _pytest.config import Config
    from _pytest.terminal import TerminalReporter

# The `pytester` fixture, used by tests/test_bedtools_gate.py.
pytest_plugins = ["pytester"]

#: Set this environment variable to "1" to make a missing bedtools fatal.
REQUIRE_ENV_VAR = "RTM_REQUIRE_BEDTOOLS"

#: The external binary whose absence gates the nested-insertion tests.
REQUIRED_BINARY = "bedtools"

#: Substring that marks a skip reason as belonging to the bedtools gate.
GATE_MARKER = "bedtools"

#: Stable phrase in the end-of-run notice, asserted on by the tests.
NOTICE_TITLE = "bedtools-gated tests skipped"


def _strict_mode_enabled() -> bool:
    """True when the user asked for a missing bedtools to fail the run."""
    return os.environ.get(REQUIRE_ENV_VAR, "").strip() == "1"


def _skip_reason(report: pytest.TestReport) -> str:
    """The human-readable skip reason for a report, as a plain string."""
    longrepr = report.longrepr
    if isinstance(longrepr, tuple) and len(longrepr) == 3:
        # (path, lineno, "Skipped: <reason>") — the skipif/call form.
        return str(longrepr[2])
    return str(longrepr)


def _is_gate_skip(report: pytest.TestReport) -> bool:
    """True for any skip whose *reason* blames the bedtools gate.

    Deliberately reads the reason rather than re-deriving the gated set, so a
    twelfth gated test is covered without touching this file.
    """
    return report.outcome == "skipped" and GATE_MARKER in _skip_reason(report).lower()


def _usage_error(binary: str) -> pytest.UsageError:
    return pytest.UsageError(
        f"{REQUIRE_ENV_VAR}=1 is set but the {binary!r} binary was not found on PATH.\n"
        f"The nested-insertion annotation feature is exercised only by {binary}-gated tests\n"
        "(tests/test_nested_rmsk_bedtools.py). A run that skipped them reports a false pass,\n"
        f"which is what {REQUIRE_ENV_VAR}=1 exists to prevent.\n"
        f"  Fix: install {binary} (e.g. `brew install {binary}`, `apt-get install {binary}`),\n"
        f"  or unset {REQUIRE_ENV_VAR} to let {binary}-gated tests skip again."
    )


class _BedtoolsGatePlugin:
    """Per-session state for the two guards. One instance per pytest session."""

    def __init__(self) -> None:
        self.gate_skips: list[str] = []

    def pytest_collection_modifyitems(self, config: Config, items: list[pytest.Item]) -> None:
        # Strict mode only. An unset RTM_REQUIRE_BEDTOOLS leaves collection alone.
        if not _strict_mode_enabled():
            return
        if shutil.which(REQUIRED_BINARY) is not None:
            return
        # UsageError => pytest exits 4 with "ERROR: ..." and blames no test.
        raise _usage_error(REQUIRED_BINARY)

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if _is_gate_skip(report):
            self.gate_skips.append(report.nodeid)

    def pytest_terminal_summary(self, terminalreporter: TerminalReporter) -> None:
        count = len(self.gate_skips)
        if not count:
            return
        binary = REQUIRED_BINARY
        plural = "test" if count == 1 else "tests"
        write = terminalreporter.write_line
        terminalreporter.write_sep("=", NOTICE_TITLE.upper())
        write("")
        write(f"  {count} {plural} skipped because the {binary} binary is not on PATH.")
        write("")
        write("  The nested-insertion annotation feature (_annotate_nested_retrotransposon)")
        write("  is UNVERIFIED by this run. A green suite is NOT evidence it works.")
        write("")
        write(f"  Fix: install {binary} (e.g. `brew install {binary}`) and re-run. To turn a")
        write(f"  missing {binary} into a hard error instead of a silent skip, set")
        write(f"  {REQUIRE_ENV_VAR}=1.")
        write("")


def pytest_configure(config: Config) -> None:
    config.pluginmanager.register(_BedtoolsGatePlugin(), "rtm-bedtools-gate")
