"""Root conftest: pandas future flags and optional pysam import fallback.

pysam (and by extension htslib) requires a system-level libhts shared
library that is not available on all CI / notebook compute environments.
When pysam cannot be imported, inject a MagicMock so pure-Python unit tests
can still be collected. When pysam is installed, use the real library.
"""

from __future__ import annotations

import sys
from unittest.mock import MagicMock

import pandas as _pd

# Pandas 2.2 emits FutureWarning when fillna() silently downcasts object
# Series. Opt into the 2.2+ behavior so those warnings are not test noise.
# Pandas 3 already uses that behavior and deprecates the option (Pandas4Warning).
_pd_major_minor = tuple(int(x) for x in _pd.__version__.split(".")[:2])
if _pd_major_minor < (3, 0):
    _pd.set_option("future.no_silent_downcasting", True)

try:
    import pysam  # noqa: F401
except ImportError:
    _mock = MagicMock()
    sys.modules["pysam"] = _mock
    for _sub in ("libchtslib", "libcsamfile", "libcbcf", "libcfaidx"):
        sys.modules[f"pysam.{_sub}"] = MagicMock()
