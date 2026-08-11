"""Root conftest: stub out heavy C-extension / system-binary dependencies.

pysam (and by extension htslib) requires a system-level libhts shared
library that is not available on all CI / notebook compute environments.
The unit tests in this suite are pure-Python logic tests that import from
modules which happen to import pysam at the top level, but no test
actually calls into pysam itself.

By injecting a MagicMock into sys.modules *before* pytest begins
collecting test files, all ``import pysam`` statements resolve to the
mock and the C extension is never loaded.  On environments where pysam
is already present (i.e. already in sys.modules from a successful
earlier import), the real library is used instead.
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

if "pysam" not in sys.modules:
    _mock = MagicMock()
    sys.modules["pysam"] = _mock
    for _sub in ("libchtslib", "libcsamfile", "libcbcf", "libcfaidx"):
        sys.modules[f"pysam.{_sub}"] = MagicMock()
