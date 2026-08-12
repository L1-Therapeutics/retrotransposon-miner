"""Shared internal utilities for retro_miner.

This module intentionally has no intra-package imports so it can be used
by any module in the package without risk of circular imports.
"""
from __future__ import annotations

import re


def safe_locus_id(chrom: str, start: int, end: int) -> str:
    """Return a filesystem-safe locus identifier string.

    Special characters in *chrom* (e.g. spaces, slashes) are replaced with
    underscores so the result can be used as a directory or file-name stem.
    """
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    return f"{chrom_safe}_{int(start)}_{int(end)}"
