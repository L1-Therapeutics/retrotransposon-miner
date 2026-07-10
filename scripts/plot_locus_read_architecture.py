#!/usr/bin/env python3
"""CLI wrapper for schematic MEI read-architecture plots.

Implementation lives in ``retro_miner.read_architecture`` so annotate can reuse
the same batch/cached path.
"""

from __future__ import annotations

from retro_miner.read_architecture import main


if __name__ == "__main__":
    main()
