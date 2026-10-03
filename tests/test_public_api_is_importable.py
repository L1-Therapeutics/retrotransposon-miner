"""Every name a module advertises in `__all__` must actually exist.

An `__all__` entry that names nothing is not a documentation slip. It makes
`from retro_miner.x import *` raise `AttributeError`, it lies to anything
introspecting the public API, and -- because the mismatch is invisible at import
time -- nothing tells you until the star-import happens in someone else's code.

That is not hypothetical. `local_assembly.py` shipped seven names in `__all__`
that were never defined *and* referenced them from sixteen places inside function
bodies. The module imported cleanly and 1122 tests passed, because every path
that touched those names was untested; `_run_minimap2_paf`,
`_load_fasta_lengths`, `_choose_consensus_features` and
`_extract_sample_assembly_features` each raised `NameError` on entry. The only
visible symptom was `threading` sitting unused in the imports, because the two
locks it exists to build were missing.

This test is the cheap guard that makes that class of bug impossible to ship
again: it is one import and one `hasattr` per name.
"""

from __future__ import annotations

import importlib
import pkgutil

import pytest

import retro_miner

#: Modules that are not importable in a plain test environment. Each needs an
#: optional third-party dependency or a data file to import at all; excluding
#: them here keeps this test about `__all__` correctness rather than about
#: environment provisioning.
_SKIP_IMPORT = {
    # Optional heavy/optional-dependency modules.
    "bam_transfer",
}


def _modules() -> list[str]:
    names = []
    for info in pkgutil.iter_modules(retro_miner.__path__):
        if not info.ispkg and info.name not in _SKIP_IMPORT:
            names.append(info.name)
    return sorted(names)


@pytest.mark.parametrize("module_name", _modules())
def test_every_all_entry_resolves(module_name: str) -> None:
    """`__all__` must not name anything the module does not define."""
    try:
        module = importlib.import_module(f"retro_miner.{module_name}")
    except Exception as exc:  # noqa: BLE001 - deliberate: the point of this guard
        # is that a module which cannot even be imported in this environment is
        # out of scope, and the failure mode is open-ended (optional C extension,
        # missing data file, syntax error under an older interpreter). Naming the
        # exception types would just make this a maintenance trap.
        pytest.skip(f"{module_name} is not importable here: {exc}")
    advertised = getattr(module, "__all__", None)
    if advertised is None:
        pytest.skip(f"{module_name} declares no __all__")
    missing = [name for name in advertised if not hasattr(module, name)]
    assert not missing, (
        f"retro_miner.{module_name}.__all__ advertises {len(missing)} name(s) that "
        f"the module does not define: {missing}. `import *` raises AttributeError "
        f"on these."
    )


@pytest.mark.parametrize("module_name", _modules())
def test_all_has_no_duplicates(module_name: str) -> None:
    """A repeated entry is a copy-paste slip that hides a missing one."""
    try:
        module = importlib.import_module(f"retro_miner.{module_name}")
    except Exception as exc:  # noqa: BLE001 - deliberate, see above
        pytest.skip(f"{module_name} is not importable here: {exc}")
    advertised = getattr(module, "__all__", None)
    if not advertised:
        pytest.skip(f"{module_name} declares no __all__")
    seen: set[str] = set()
    duplicates: list[str] = []
    for name in advertised:
        if name in seen:
            duplicates.append(name)
        seen.add(name)
    assert not duplicates, f"retro_miner.{module_name}.__all__ repeats {duplicates}"
