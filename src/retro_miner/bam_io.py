"""Open BAM or CRAM with the correct pysam mode.

HTSlib rejects CRAM when the handle is opened as ``rb``. Use ``rc`` for
``.cram`` and pass ``reference_filename`` (or ``RTM_ALIGNMENT_REFERENCE``).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pysam


def alignment_path_is_cram(path: str | Path) -> bool:
    text = str(path).strip()
    parsed = urlparse(text)
    name = Path(parsed.path).name.lower() if parsed.scheme else Path(text.split("?", 1)[0]).name.lower()
    return name.endswith(".cram")


def alignment_open_mode(path: str | Path) -> str:
    return "rc" if alignment_path_is_cram(path) else "rb"


def resolve_alignment_reference(explicit: str | Path | None = None) -> str | None:
    if explicit:
        text = str(explicit).strip()
        return text or None
    env = (os.environ.get("RTM_ALIGNMENT_REFERENCE") or "").strip()
    return env or None


def open_alignment(
    path: str | Path,
    *,
    reference_filename: str | Path | None = None,
    threads: int = 1,
) -> pysam.AlignmentFile:
    """Open a BAM or CRAM. BAM callers can omit the reference."""
    kwargs: dict[str, Any] = {"threads": int(threads)}
    ref = resolve_alignment_reference(reference_filename)
    if ref:
        kwargs["reference_filename"] = ref
    return pysam.AlignmentFile(str(path), alignment_open_mode(path), **kwargs)
