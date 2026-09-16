#!/usr/bin/env python3
"""CLI wrapper for schematic MEI read-architecture plots.

Implementation lives in ``retro_miner.read_architecture`` so annotate can reuse
the same batch/cached path.
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Repo-root sys.path resolution for standalone execution
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from retro_miner.read_architecture import main

if __name__ == "__main__":
    main()
