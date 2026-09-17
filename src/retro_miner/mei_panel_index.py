"""One-time classic ``bwa index`` for MEI remap panels.

Disease and control remaps share the same Dfam/full-consensus FASTA. Index those
panels during ``download_public_data.py`` postprocess (and cache on S3). Annotate
calls the same helpers if the index is missing, **serially**, before the
disease∥control thread pool.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

PUBLIC_MEI_REMAP_FASTA_RELPATHS = (
    "retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta",
    "retrotransposon_db/full_consensus/mei_full_canonical.panel.fa",
    "retrotransposon_db/full_consensus/mei_full_canonical.ucsc_repeatbrowser.fa",
)


def nopolya_sidecar_path(src: Path) -> Path:
    src = Path(src)
    return src.with_name(f"{src.stem}.nopolya{src.suffix}")


def has_bwa_index(fasta: Path) -> bool:
    fasta = Path(fasta)
    return all(Path(f"{fasta}{suffix}").exists() for suffix in BWA_INDEX_SUFFIXES)


def drop_bwa_index(fasta: Path) -> None:
    fasta = Path(fasta)
    for suffix in BWA_INDEX_SUFFIXES:
        idx = Path(f"{fasta}{suffix}")
        if idx.exists():
            idx.unlink()


def ensure_polya_trimmed_mei_fasta(mei_fasta: Path) -> Path:
    """Return a polyA-trimmed sidecar, refreshing when the source FASTA is newer."""
    from retro_miner.mei_support import _write_polya_trimmed_fasta

    src = Path(mei_fasta)
    if not src.exists():
        raise FileNotFoundError(f"MEI FASTA not found: {src}")
    dst = nopolya_sidecar_path(src)
    src_mtime = src.stat().st_mtime
    needs = (not dst.exists()) or (dst.stat().st_mtime < src_mtime) or (dst.stat().st_size <= 0)
    if needs:
        _write_polya_trimmed_fasta(src, dst)
        drop_bwa_index(dst)
    return dst


def ensure_bwa_index(fasta: Path, *, force: bool = False) -> dict[str, Any]:
    """Build a classic ``bwa index`` next to ``fasta`` when missing."""
    fasta = Path(fasta)
    if not fasta.exists() or fasta.stat().st_size <= 0:
        raise FileNotFoundError(f"FASTA not found or empty: {fasta}")
    if not force and has_bwa_index(fasta):
        return {"status": "skipped_exists", "path": str(fasta)}
    if shutil.which("bwa") is None:
        raise RuntimeError("bwa not found on PATH; cannot index MEI remap FASTA")
    if force:
        drop_bwa_index(fasta)
    subprocess.run(
        ["bwa", "index", str(fasta)],
        check=True,
        capture_output=True,
        text=True,
    )
    return {"status": "indexed", "path": str(fasta)}


def ensure_mei_remap_bwa_index(mei_fasta: Path, *, force: bool = False) -> dict[str, Any]:
    """Ensure the polyA-trimmed sidecar and its bwa index exist.

    Used by ``download_public_data.py`` postprocess and by annotate (once, before
    parallel disease/control remaps).
    """
    remap_fa = ensure_polya_trimmed_mei_fasta(mei_fasta)
    result = ensure_bwa_index(remap_fa, force=force)
    result["source_fasta"] = str(Path(mei_fasta))
    result["nopolya_fasta"] = str(remap_fa)
    return result


def index_public_mei_remap_fastas(outdir: Path, *, force: bool = False) -> dict[str, Any]:
    """Index every MEI remap FASTA under a public-data outdir, if present."""
    indexed: list[dict[str, Any]] = []
    skipped_missing: list[str] = []
    for rel in PUBLIC_MEI_REMAP_FASTA_RELPATHS:
        path = Path(outdir) / rel
        if not path.exists() or path.stat().st_size <= 0:
            skipped_missing.append(rel)
            continue
        indexed.append(ensure_mei_remap_bwa_index(path, force=force))
    status = "ok" if indexed else "skipped_missing_mei_fasta"
    return {
        "status": status,
        "indexed": indexed,
        "skipped_missing": skipped_missing,
    }
