"""Guards the pytest `pythonpath` setting that pins imports to this checkout.

This repo is checked out as several git worktrees sharing one .git, and an
editable install of one of them puts a *sibling* worktree's src/ on sys.path.
Without `pythonpath = ["src"]` the suite can import another branch's code and
still pass, with no error and no warning.
"""

import tomllib
from pathlib import Path

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
    # This is the assertion that would have caught the real problem: it is the
    # one that fails, silently-but-loudly-here, when sys.path resolves
    # retro_miner to a sibling worktree instead of the tree these tests live in.
    import retro_miner

    module_path = Path(retro_miner.__file__).resolve()
    assert module_path.is_relative_to(REPO_ROOT / "src"), (
        f"retro_miner resolved to {module_path}, which is outside "
        f"{REPO_ROOT / 'src'} -- a sibling worktree's editable install won on sys.path"
    )
