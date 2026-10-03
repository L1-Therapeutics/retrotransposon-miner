"""Guards the pytest `pythonpath` setting that pins imports to this checkout.

This repo is checked out as several git worktrees sharing one .git, and an
editable install of one of them puts a *sibling* worktree's src/ on sys.path.
Without `pythonpath = ["src"]` the suite can import another branch's code and
still pass, with no error and no warning.
"""

from __future__ import annotations

import importlib.util
import tomllib
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_pyproject_pins_pythonpath_to_src():
    config = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    ini_options = config["tool"]["pytest"]["ini_options"]
    pythonpath = ini_options.get("pythonpath", [])
    assert "src" in pythonpath, (
        "tool.pytest.ini_options.pythonpath must contain 'src' so pytest "
        "prepends this checkout's src/ to sys.path"
    )


def test_retro_miner_resolves_inside_this_checkout():
    """The suite must exercise THIS checkout's src, never a sibling worktree's.

    Skips when retro_miner is not installed at all: bare pytest runs get the
    package from pythonpath = ["src"], so there is nothing to police unless
    some sys.path entry (e.g. an editable install) could supply a copy.
    """
    spec = importlib.util.find_spec("retro_miner")
    if spec is None or spec.origin is None:
        pytest.skip("retro_miner not installed; pythonpath = ['src'] supplies it")
    module_path = Path(spec.origin).resolve()
    assert module_path.is_relative_to(REPO_ROOT / "src"), (
        f"retro_miner resolved to {module_path}, which is outside "
        f"{REPO_ROOT / 'src'} -- a sibling worktree's editable install won on sys.path"
    )
