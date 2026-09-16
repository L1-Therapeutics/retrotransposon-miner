"""Focused tests for scripts/validate_environment.sh.

The real shell script is executed against fake tool binaries placed first on
PATH. The fakes use `#!<real python>` as their interpreter (never `env
python3`) so a fake `python3` on PATH cannot re-dispatch to itself. No real
third-party binaries, credentials, or network are used, and no subprocess is
ever launched with shell=True.
"""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "validate_environment.sh"

REQUIRED_BINS = [
    "git",
    "jupyter",
    "samtools",
    "bedtools",
    "minimap2",
    "bwa",
    "bwa-mem2",
    "bcftools",
    "liftOver",
    "bigBedToBed",
    "bigWigToBedGraph",
]
OPTIONAL_BINS = ["igv", "spades.py", "java", "Xvfb", "xvfb-run"]


class Runner:
    def __init__(self, root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        bindir = root / "bin"
        bindir.mkdir()
        self.bindir = bindir
        self.home = root / "home"
        self.home.mkdir()
        self.tmpdir = root / "tmpdir"
        self.tmpdir.mkdir()

        shebang = f"#!{sys.executable}"
        for name in REQUIRED_BINS:
            self._fake(name, shebang, f'print("{name} 0.1")')
        self._fake("bwa", shebang, 'import sys\nprint("usage: bwa <command>", file=sys.stderr)\nsys.exit(1)')
        self._fake("bwa-mem2", shebang, 'print("bwa-mem2 2.2.1")')
        self._fake("liftOver", shebang, 'import sys\nprint("liftOver: no args", file=sys.stderr)\nsys.exit(255)')
        self._fake_python(shebang, 'print("Python modules OK")')

        monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
        monkeypatch.setenv("HOME", str(self.home))
        monkeypatch.setenv("PYTHON_BIN", "python3")
        monkeypatch.delenv("DISPLAY", raising=False)

        self.returncode: int | None = None
        self.stdout = ""
        self.stderr = ""

    def _fake(self, name: str, shebang: str, body: str) -> None:
        exe = self.bindir / name
        exe.write_text(shebang + "\n" + body + "\n", encoding="utf-8")
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def _fake_python(self, shebang: str, body: str) -> None:
        self._fake("python3", shebang, body)

    def remove(self, name: str) -> None:
        (self.bindir / name).unlink()

    def isolated_path(self) -> str:
        realdir = self.tmpdir.parent / "real"
        realdir.mkdir()
        for tool in ("bash", "awk", "uname"):
            resolved = shutil.which(tool)
            assert resolved is not None, f"missing real {tool} for symlink"
            (realdir / tool).symlink_to(Path(resolved))
        return str(self.bindir) + os.pathsep + str(realdir)

    def run(self, **env_overrides: str) -> None:
        env = dict(os.environ)
        env.update(env_overrides)
        env["TMPDIR"] = str(self.tmpdir)

        assert shutil.which("python3", path=env["PATH"]) == str(self.bindir / "python3"), (
            "the fake python3 must be what PATH resolves, never the real interpreter"
        )
        assert shutil.which("samtools", path=env["PATH"]) == str(self.bindir / "samtools"), (
            "the fake samtools must be what PATH resolves, never the real binary"
        )

        proc = subprocess.run(
            [str(SCRIPT)],
            capture_output=True,
            text=True,
            env=env,
            cwd=str(REPO_ROOT),
        )
        self.returncode = proc.returncode
        self.stdout = proc.stdout
        self.stderr = proc.stderr

    def assert_no_temp_leftovers(self) -> None:
        assert list(self.tmpdir.iterdir()) == [], "temp files leaked into TMPDIR"


@pytest.fixture()
def runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Runner:
    return Runner(tmp_path, monkeypatch)


def test_all_tools_present_exits_zero(runner: Runner) -> None:
    runner.run()

    assert runner.returncode == 0
    assert "Validating retrotransposon-miner environment..." in runner.stdout
    assert "All required tools detected." in runner.stdout
    assert "samtools: samtools 0.1" in runner.stdout
    assert "bwa-mem2: bwa-mem2 2.2.1" in runner.stdout
    assert "bwa: usage: bwa <command>" in runner.stdout
    assert "Environment validation complete." in runner.stdout
    assert "ERROR:" not in runner.stderr
    runner.assert_no_temp_leftovers()


def test_report_section_tolerates_nonzero_version_exits(runner: Runner) -> None:
    runner._fake("samtools", f"#!{sys.executable}", 'import sys\nprint("samtools explode", file=sys.stderr)\nsys.exit(1)')
    runner._fake("bcftools", f"#!{sys.executable}", 'import sys\nsys.exit(9)')

    runner.run()

    # Informational version prints must never fail a successful validation.
    assert runner.returncode == 0
    assert "All required tools detected." in runner.stdout
    assert "samtools: samtools explode" in runner.stdout
    assert "Environment validation complete." in runner.stdout


def test_lift_over_nonzero_usage_exit_is_not_fatal(runner: Runner) -> None:
    runner.run()

    # liftOver invoked with no args prints usage and exits 255; the validation
    # must still succeed and report the tool was found.
    assert runner.returncode == 0
    assert "liftOver: liftOver: no args" in runner.stdout
    assert "Environment validation complete." in runner.stdout


def test_missing_required_binary_fails(runner: Runner) -> None:
    runner.remove("bwa-mem2")

    runner.run(PATH=runner.isolated_path())

    assert runner.returncode == 1
    assert "ERROR: missing binary: bwa-mem2" in runner.stderr
    assert "Environment validation complete." not in runner.stdout


def test_missing_python_module_fails(runner: Runner) -> None:
    runner._fake_python(
        f"#!{sys.executable}",
        'import sys\nprint("ERROR: missing python modules: matplotlib, click", file=sys.stderr)\nsys.exit(1)',
    )

    runner.run()

    assert runner.returncode == 1
    assert "ERROR: missing python modules: matplotlib, click" in runner.stderr


def test_missing_python_bin_fails(runner: Runner) -> None:
    runner.run(PYTHON_BIN="python-does-not-exist-xyz")

    assert runner.returncode == 1
    assert "python-does-not-exist-xyz not found in PATH." in runner.stderr


def test_optional_binary_missing_is_warn_only(runner: Runner) -> None:
    runner.run(PATH=runner.isolated_path())

    assert runner.returncode == 0
    for name in OPTIONAL_BINS:
        assert f"WARN: optional binary not found: {name}" in runner.stderr
    assert "ERROR:" not in runner.stderr


def test_fake_python_is_the_resolved_interpreter(runner: Runner) -> None:
    assert shutil.which("python3", path=os.environ["PATH"]) == str(runner.bindir / "python3")
    assert shutil.which("samtools", path=os.environ["PATH"]) == str(runner.bindir / "samtools")
