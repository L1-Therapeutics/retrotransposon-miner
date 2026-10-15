"""Shared internal utilities for retro_miner.

This module intentionally has no intra-package imports so it can be used
by any module in the package without risk of circular imports.
"""

import subprocess


def run_subprocess(
    cmd: list[str],
    *,
    label: str = "command",
    check: bool = True,
    timeout: float | None = None,
    capture_output: bool = True,
    text: bool = True,
    env: dict[str, str] | None = None,
    timeout_default: float | None = 3600.0,
) -> subprocess.CompletedProcess[str]:
    """Run a subprocess with consistent error handling and logging.

    Parameters
    ----------
    cmd : list[str]
        Command and arguments to run.
    label : str
        Human-readable label for error messages (e.g. "bwa mem", "bedtools intersect").
    check : bool
        If True, raise :exc:`RuntimeError` on non-zero exit. If False, return the
        completed process and let the caller check ``returncode``.
    timeout : float | None
        Optional timeout in seconds. Pass ``None`` for no timeout, or use
        ``timeout_default`` for a sensible default.
    capture_output : bool
        If True, capture stdout and stderr.
    text : bool
        If True, return stdout/stderr as strings instead of bytes.
    env : dict[str, str] | None
        Optional environment variables to pass to the subprocess.
    timeout_default : float | None
        Default timeout in seconds when ``timeout`` is ``None``. Pass ``None``
        to allow no timeout. Defaults to 1 hour for long-running bioinformatics tools.

    Returns
    -------
    subprocess.CompletedProcess[str]
        The completed process. When ``check=False``, the caller should check
        ``proc.returncode`` and inspect ``proc.stderr`` on failure.

    Raises
    ------
    RuntimeError
        When ``check=True`` and the command exits non-zero, or when the command
        times out.
    FileNotFoundError
        When the executable is not found on PATH.
    """
    effective_timeout = timeout if timeout is not None else timeout_default
    try:
        proc = subprocess.run(
            cmd,
            check=False,
            capture_output=capture_output,
            text=text,
            timeout=effective_timeout,
            env=env,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"{label} not found: {cmd[0]!r} — is it on PATH?"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"{label} timed out after {effective_timeout:.1f}s: {' '.join(cmd)}"
        ) from exc

    if check and proc.returncode != 0:
        err = (proc.stderr or "").strip()
        out = (proc.stdout or "").strip()
        msg = f"{label} failed (exit={proc.returncode}): {' '.join(cmd)}"
        if err:
            msg += f"\nstderr:\n{err}"
        if out:
            truncated = out if len(out) <= 2000 else out[:2000] + "\n... (truncated)"
            msg += f"\nstdout (truncated):\n{truncated}"
        raise RuntimeError(msg)

    return proc


def run_subprocess_ignore_error(
    cmd: list[str],
    *,
    label: str = "command",
    timeout: float | None = None,
) -> subprocess.CompletedProcess[str] | None:
    """Run a subprocess, returning None on any failure.

    Use for best-effort operations where failure is acceptable (e.g. indexing
    a FASTA file that may already be indexed).

    Parameters
    ----------
    cmd : list[str]
        Command and arguments to run.
    label : str
        Human-readable label for log messages.
    timeout : float | None
        Optional timeout in seconds.

    Returns
    -------
    subprocess.CompletedProcess[str] | None
        The completed process on success, or ``None`` on failure.
    """
    try:
        return run_subprocess(cmd, label=label, check=False, timeout=timeout)
    except (RuntimeError, FileNotFoundError, subprocess.TimeoutExpired):
        return None

from __future__ import annotations

import gzip
import re
import warnings
from pathlib import Path
from typing import IO


def safe_locus_id(chrom: str, start: int, end: int) -> str:
    """Return a filesystem-safe locus identifier string.

    Special characters in *chrom* (e.g. spaces, slashes) are replaced with
    underscores so the result can be used as a directory or file-name stem.
    """
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    return f"{chrom_safe}_{int(start)}_{int(end)}"


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
        best_len = n
        best_frac = float(span.count(best_base)) / float(len(span)) if span else best_frac
        span = s
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
