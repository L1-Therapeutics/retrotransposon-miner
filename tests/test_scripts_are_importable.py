"""Every script in `scripts/` must import cleanly against this repo's own code.

`pyproject.toml` sets `pythonpath = ["src"]`, so these scripts import the
vendored `src/retro_miner/` in this repository rather than the sibling
`retrotransposon-miner` package. Nothing asserted that, and the two copies have
diverged: this one carries seven modules the sibling lacks, and the sibling
carries two this one lacks. A cross-repo import can therefore fail outright
when the wrong copy lands first on the path, and nothing in the suite noticed.

The check runs in-process. Spawning one interpreter per script costs ~143 s
because each pays for pandas and pysam; importing them all in one process
costs ~26 s and catches the same failure class -- a missing name or a module
that no longer exists raises at import time, long before any work happens.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
VENDORED_SRC = REPO / "src"
VENDORED = VENDORED_SRC / "retro_miner"

SCRIPT_NAMES = sorted(path.stem for path in SCRIPTS.glob("*.py"))


def test_there_are_scripts_to_check():
    """A silently empty glob would make this file pass without checking anything."""
    assert len(SCRIPT_NAMES) > 20, SCRIPT_NAMES


@pytest.mark.skipif(not VENDORED.is_dir(), reason="no vendored retro_miner in this repo")
def test_every_script_imports_without_error():
    for entry in (str(VENDORED_SRC), str(SCRIPTS)):
        if entry not in sys.path:
            sys.path.insert(0, entry)

    failures = []
    for stem in SCRIPT_NAMES:
        path = SCRIPTS / f"{stem}.py"
        name = f"_import_smoke_{stem}"
        try:
            spec = importlib.util.spec_from_file_location(name, path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
        except Exception as error:  # noqa: BLE001 - the point is to report any of them
            failures.append(f"{path.name}: {type(error).__name__}: {error}")
        finally:
            sys.modules.pop(name, None)

    assert not failures, "scripts failed to import:\n  " + "\n  ".join(failures)


@pytest.mark.skipif(not VENDORED.is_dir(), reason="no vendored retro_miner in this repo")
def test_scripts_resolve_retro_miner_from_this_repository():
    """Pin which copy of the package the scripts actually get.

    Nothing else here inserts a path: this relies on the `pythonpath = ["src"]`
    in pyproject.toml, so it fails if that configuration stops taking effect
    and a sibling checkout silently becomes the code under test.
    """
    import retro_miner

    resolved = Path(retro_miner.__file__).resolve()
    assert resolved.is_relative_to(VENDORED_SRC.resolve()), (
        f"retro_miner resolved to {resolved}, outside this repository")
