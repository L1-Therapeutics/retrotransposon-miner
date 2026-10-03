"""Tests for the bedtools gate visibility machinery in tests/conftest.py.

Each test builds a throwaway suite in a pytester tmp dir and copies the real
``tests/conftest.py`` into it, so the guard is exercised as shipped. Absence of
bedtools is simulated by pointing PATH at an empty directory, never by
uninstalling anything, so these tests behave identically on a laptop with
bedtools and on a fresh CI box without it.
"""

from __future__ import annotations

import re
import shutil
import stat
from pathlib import Path

import pytest

# The conftest under test, copied verbatim into each throwaway suite.
REPO_CONFTEST = Path(__file__).parent / "conftest.py"

REQUIRE_ENV_VAR = "RTM_REQUIRE_BEDTOOLS"
NOTICE_TITLE = "BEDTOOLS-GATED TESTS SKIPPED"

# Two bedtools-gated tests, mirroring the skipif pattern used across
# tests/test_nested_rmsk_bedtools.py.
GATED_SUITE = """
import shutil

import pytest

gate = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")


@gate
def test_gated_first():
    pass


@gate
def test_gated_second():
    pass
"""

# A suite with no bedtools-gated tests at all.
UNGATED_SUITE = """
def test_plain_first():
    pass


def test_plain_second():
    pass
"""

# A suite that skips, but never for bedtools reasons.
UNRELATED_SKIP_SUITE = """
import pytest


@pytest.mark.skip(reason="waiting on a fixture from the sequencing core")
def test_skipped_for_another_reason():
    pass


def test_plain():
    pass
"""


def _write_suite(pytester: pytest.Pytester, test_source: str) -> None:
    """Lay down a minimal suite guarded by the real tests/conftest.py."""
    (pytester.path / "conftest.py").write_text(REPO_CONFTEST.read_text(encoding="utf-8"), encoding="utf-8")
    (pytester.path / "test_probe_suite.py").write_text(test_source, encoding="utf-8")


def _hide_bedtools(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Simulate bedtools being absent, and assert the simulation took."""
    empty_bin = tmp_path / "bin-with-no-bedtools"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    assert shutil.which("bedtools") is None, "PATH simulation failed to hide bedtools"


def _provide_bedtools(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Simulate bedtools being present, and assert the simulation took."""
    bin_dir = tmp_path / "bin-with-bedtools"
    bin_dir.mkdir()
    stub = bin_dir / "bedtools"
    stub.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", str(bin_dir))
    assert shutil.which("bedtools") is not None, "PATH simulation failed to expose bedtools"


def _run(pytester: pytest.Pytester) -> tuple[int, str, str]:
    """Run the throwaway suite in a subprocess. Returns (exit code, stdout, stdout+stderr).

    Unrelated auto-loaded plugins are disabled: they are slow to import and
    cannot affect the gate, and the gate is the only thing under test here.
    """
    result = pytester.runpytest_subprocess(
        "-p",
        "no:cacheprovider",
        "-p",
        "no:hypothesis",
        "-p",
        "no:timeout",
        "-p",
        "no:anyio",
        "-rs",
    )
    stdout = result.stdout.str()
    return int(result.ret), stdout, stdout + result.stderr.str()


def test_unset_env_var_preserves_skip_and_success(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """Default behaviour: no env var, no bedtools -> skip, and the run passes."""
    _write_suite(pytester, GATED_SUITE)
    _hide_bedtools(monkeypatch, tmp_path)
    monkeypatch.delenv(REQUIRE_ENV_VAR, raising=False)

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert re.search(r"\b2 skipped\b", out), f"expected both gated tests to skip\n{both}"
    assert "bedtools not on PATH" in out


def test_strict_mode_fails_loudly_when_bedtools_absent(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """RTM_REQUIRE_BEDTOOLS=1 turns the silent skip into a usage error."""
    _write_suite(pytester, GATED_SUITE)
    _hide_bedtools(monkeypatch, tmp_path)
    monkeypatch.setenv(REQUIRE_ENV_VAR, "1")

    ret, _out, both = _run(pytester)

    assert ret == pytest.ExitCode.USAGE_ERROR, f"expected usage error, got {ret}\n{both}"
    assert "bedtools" in both
    assert REQUIRE_ENV_VAR in both


def test_strict_mode_passes_when_bedtools_present(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """Strict mode must not fire when the binary is there."""
    _write_suite(pytester, GATED_SUITE)
    _provide_bedtools(monkeypatch, tmp_path)
    monkeypatch.setenv(REQUIRE_ENV_VAR, "1")

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert re.search(r"\b2 passed\b", out), f"gated tests should have run\n{both}"


def test_terminal_summary_prints_notice_on_gate_skip(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """A gate skip is announced, with the count, the stakes, and the fix."""
    _write_suite(pytester, GATED_SUITE)
    _hide_bedtools(monkeypatch, tmp_path)
    monkeypatch.delenv(REQUIRE_ENV_VAR, raising=False)

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert NOTICE_TITLE in out, f"notice header missing\n{both}"
    assert re.search(r"2 tests skipped because the bedtools binary is not on PATH", out), both
    assert "UNVERIFIED" in out
    assert "_annotate_nested_retrotransposon" in out
    assert "brew install bedtools" in out
    assert f"{REQUIRE_ENV_VAR}=1" in out


def test_no_notice_when_gated_tests_ran(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """bedtools present -> gated tests execute -> nothing to announce."""
    _write_suite(pytester, GATED_SUITE)
    _provide_bedtools(monkeypatch, tmp_path)
    monkeypatch.delenv(REQUIRE_ENV_VAR, raising=False)

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert NOTICE_TITLE not in both, f"notice printed with no gate skip\n{both}"


def test_no_notice_when_no_gated_tests_exist(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """A suite with nothing gated must not produce the banner either."""
    _write_suite(pytester, UNGATED_SUITE)
    _hide_bedtools(monkeypatch, tmp_path)
    monkeypatch.delenv(REQUIRE_ENV_VAR, raising=False)

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert NOTICE_TITLE not in both, f"notice printed with no gate skip\n{both}"


def test_no_notice_for_unrelated_skips(pytester: pytest.Pytester, monkeypatch, tmp_path):
    """Only bedtools-reasoned skips trigger the banner, not skips in general."""
    _write_suite(pytester, UNRELATED_SKIP_SUITE)
    _hide_bedtools(monkeypatch, tmp_path)
    monkeypatch.delenv(REQUIRE_ENV_VAR, raising=False)

    ret, out, both = _run(pytester)

    assert ret == pytest.ExitCode.OK, f"expected success, got {ret}\n{both}"
    assert re.search(r"\b1 skipped\b", out), f"the unrelated skip should still be reported\n{both}"
    assert NOTICE_TITLE not in both, f"banner fired on a non-bedtools skip\n{both}"
