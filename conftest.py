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

# Opt into pandas 2.2+ future behavior for fillna/ffill/bfill on object-dtype
# Series.  Every callsite in this codebase already performs an explicit
# .astype() after fillna, so the new non-downcasting behavior is safe and
# this eliminates FutureWarnings during test runs.
_pd.set_option("future.no_silent_downcasting", True)

try:
    import pysam  # noqa: F401
except ImportError:
    _mock = MagicMock()
    sys.modules["pysam"] = _mock
    for _sub in ("libchtslib", "libcsamfile", "libcbcf", "libcfaidx"):
        sys.modules[f"pysam.{_sub}"] = MagicMock()
