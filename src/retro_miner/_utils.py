"""Shared internal utilities for retro_miner.

This module intentionally has no intra-package imports so it can be used
by any module in the package without risk of circular imports.
"""
from __future__ import annotations

import gzip
import math
import numbers
import re
import warnings
import pandas as pd
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, IO


def safe_locus_id(chrom: str, start: int, end: int) -> str:
    """Return a filesystem-safe locus identifier string.

    Special characters in *chrom* (e.g. spaces, slashes) are replaced with
    underscores so the result can be used as a directory or file-name stem.
    """
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    return f"{chrom_safe}_{int(start)}_{int(end)}"


def read_tsv(
    path: Path,
    *,
    usecols: Sequence[str] | Callable[[str], bool] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Read a TSV file with sensible defaults.

    Wraps :func:`pandas.read_csv` with ``sep="\t"`` and ``low_memory=False``
    unless the caller already supplied them. Additional keyword arguments are
    forwarded unchanged. When ``usecols`` is a sequence, every listed column
    must exist; a callable selects only the columns present in the file.
    """
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("low_memory", False)
    if usecols is not None:
        kwargs.setdefault("usecols", usecols)
    return pd.read_csv(path, **kwargs)


def _safe_int(value: object, default: int = 0) -> int:
    """Safely convert *value* to int, returning *default* on failure.

    - int values are returned directly (no float round-trip, so huge ints
      like ``10**400`` are preserved exactly).
    - float values are converted if finite, otherwise *default*.
    - strings are stripped and parsed via ``float()`` first.
    - NumPy/pandas numeric scalars (``np.int64``, ``np.float32``,
      ``np.bool_``, ...) are converted even though they do not subclass the
      builtin ``int``/``float``; table rows and Series values produce these.
    - everything else returns *default*.
    """
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return default if not math.isfinite(value) else int(value)
    if isinstance(value, str):
        try:
            number = float(value.strip())
        except ValueError:
            return default
        return default if not math.isfinite(number) else int(number)
    # NumPy/pandas scalars do not always subclass the builtin numeric types:
    # np.int64 is a numbers.Integral but not an int, np.float32 is a
    # numbers.Real but not a float, and np.bool_ is neither (it still defines
    # __int__). Without these checks integer table columns were read as 0.
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        number = float(value)
        return default if not math.isfinite(number) else int(number)
    to_int = getattr(value, "__int__", None)
    if to_int is None:
        return default
    try:
        return int(to_int())
    except (TypeError, ValueError, OverflowError):
        return default


def _open_textmaybe_gz(path: Path) -> IO[str]:
    """Return an open text file handle for *path*, decompressing .gz transparently.

    The returned handle is itself a context manager::

        with _open_textmaybe_gz(path) as handle:
            for line in handle:
                ...
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _poly_at_stats(seq: str) -> tuple[int, float, str]:
    """PolyA/T stats: longest dominant-base run, dominant-base fraction, base.

    Purity is ``max(n_A, n_T) / length`` — mostly A **or** mostly T — not
    combined A+T. Mixed AT sequence scores ~0.5 and fails typical thresholds.
    """
    s = (seq or "").upper()
    if not s:
        return (0, 0.0, "")
    n_a = s.count("A")
    n_t = s.count("T")
    if n_a <= 0 and n_t <= 0:
        return (0, 0.0, "")
    if n_a >= n_t:
        base = "A"
        n_dom = n_a
    else:
        base = "T"
        n_dom = n_t
    frac = float(n_dom) / float(len(s))
    best = 0
    cur = 0
    for ch in s:
        if ch == base:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return (int(best), float(frac), base)


def _longest_poly_at_span(
    seq: str,
    *,
    min_frac: float = 0.90,
    min_len: int = 25,
) -> tuple[int, float, str, str]:
    """Longest substring that is mostly polyA **or** mostly polyT.

    For each base in {A,T}, two-pointer search for the longest window with
    ``count(base) / window_len ≥ min_frac``. Returns
    ``(length, purity, base, span_seq)``. Length 0 if none found.

    If the span covers essentially the whole read (≥140 bp or within 2 bp of
    read length), that is a full-read polyA/T (few mismatches allowed).
    """
    s = "".join(ch for ch in (seq or "").upper() if ch in {"A", "C", "G", "T"})
    n = len(s)
    if n < int(min_len):
        return (0, 0.0, "", "")
    best_len = 0
    best_frac = 0.0
    best_base = ""
    best_ij = (0, 0)
    thr = float(min_frac)
    for base in ("A", "T"):
        left = 0
        n_base = 0
        for right in range(n):
            if s[right] == base:
                n_base += 1
            while left <= right and (n_base / float(right - left + 1)) < thr:
                if s[left] == base:
                    n_base -= 1
                left += 1
            cur_len = right - left + 1
            if cur_len >= int(min_len) and n_base > 0:
                frac = float(n_base) / float(cur_len)
                if cur_len > best_len or (cur_len == best_len and frac > best_frac):
                    best_len = cur_len
                    best_frac = frac
                    best_base = base
                    best_ij = (left, right + 1)
    if best_len <= 0:
        return (0, 0.0, "", "")
    span = s[best_ij[0] : best_ij[1]]
    if best_len >= 140 or best_len >= n - 2:
        # Promote to the whole read, then recompute purity over that span —
        # keeping the matched window's purity would overstate a widened span.
        span = s
        best_len = n
        best_frac = float(span.count(best_base)) / float(n)
    return (int(best_len), float(best_frac), best_base, span)


def _iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into a list of (name, sequence) tuples.

    Returns an empty list if *path* does not exist.  Sequences are uppercased
    and blank lines are skipped.  Only the first word of each header line is
    used as the record name.

    Records with a blank or whitespace-only header (e.g. a bare ``>``) are
    **skipped** and a :class:`UserWarning` is emitted so the caller can detect
    the malformed input.  This prevents an :exc:`IndexError` that would
    otherwise crash the parser.
    """
    if not path.exists():
        return []
    out: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq_parts)))
                raw_header = line[1:].strip()
                if not raw_header:
                    warnings.warn(
                        f"Skipping FASTA record with blank header in {path!s}",
                        UserWarning,
                        stacklevel=2,
                    )
                    name = ""
                    seq_parts = []
                else:
                    name = raw_header.split()[0]
                    seq_parts = []
            else:
                seq_parts.append(line.upper())
    if name:
        out.append((name, "".join(seq_parts)))
    return out
