"""One-time classic ``bwa index`` for MEI remap panels.

Disease and control remaps share the same Dfam/full-consensus FASTA. Index those
panels during ``download_public_data.py`` postprocess (and cache on S3). Annotate
calls the same helpers if the index is missing, **serially**, before the
disease∥control thread pool.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

PUBLIC_MEI_REMAP_FASTA_RELPATHS = (
    "retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta",
    "retrotransposon_db/full_consensus/mei_full_canonical.panel.fa",
    "retrotransposon_db/full_consensus/mei_full_canonical.ucsc_repeatbrowser.fa",
)


@dataclass(frozen=True)
class ConsensusLength:
    """A consensus length together with the axis it is measured on.

    ``axis`` is ``"full"`` for a full-length canonical record and
    ``"fragment"`` for a Dfam consensus fragment; the two are not
    interchangeable, so the number never travels without it.  ``basis`` says how
    the record was chosen: the subfamily's own full record, the family canonical
    full record, or the fragment itself.

    The choice mirrors the projection the annotate step performs
    (``_build_mei_fragment_to_full_coord_map``): a fragment with no full record of
    its own is projected onto its family canonical, so its coordinates are on
    that axis.  Checked against the cohort rather than assumed: ``L1HS_5end``
    calls reach ``MEI_3P`` 6018 (the ``L1HS_full`` length, not their own 2136),
    ``L1P3_orf2`` reaches 6128 with no ``L1P3_full`` in the panel, and
    ``L1PA8_3end`` reaches 6470 against its own ``L1PA8_full`` of 6470.
    """

    bp: int
    target: str
    axis: str
    panel: str
    basis: str
    mei_class: str = ""


#: Human MEI consensus axis lengths, for the analyses that used to estimate
#: them from their own data because no panel was on disk.
#:
#: Two panels matter, and they are different axes:
#:   * the Dfam fragment panel the remap aligns reads to
#:     (``L1HS_5end`` = 2136 bp),
#:   * the full-length canonical panel derived from it
#:     (``L1HS_full`` = 6018 bp).
#: The caller's ``MEI_5P``/``MEI_3P`` coordinates are on the *full* axis: a call
#: labelled ``L1HS_5end`` reports ``MEI_3P`` up to 6018 against its 2136 bp
#: fragment.  A length read off the wrong axis is how an analysis concludes a
#: span exceeds its own consensus, so the axis is carried with the number.
FRAGMENT_SUFFIXES = ("_5end", "_3end", "_orf1", "_orf2", "_short")
FULL_SUFFIX = "_full"

#: Fallback record per family, matching the projection map's own choice.
FAMILY_CANONICAL_FULL = {
    "LINE1": "L1HS_full",
    "ALU": "AluY_full",
    "SVA": "SVA_F_full",
}

CLASS_TO_FAMILY = {
    "LINE/L1": "LINE1",
    "LINE/L1-Tx1": "LINE1",
    "SINE/Alu": "ALU",
    "Retroposon/SVA": "SVA",
}


def family_of_label(label: str, mei_class: str = "") -> str:
    """Family of a panel label, from its class when known, else from its name."""
    if mei_class in CLASS_TO_FAMILY:
        return CLASS_TO_FAMILY[mei_class]
    name = (label or "").strip()
    if name.startswith("Alu"):
        return "ALU"
    if name.startswith("SVA"):
        return "SVA"
    if name.startswith(("L1", "HAL1")):
        return "LINE1"
    return ""

PANEL_ENV_VARS = {
    "fragment": "RTM_DFAM_MEI_PANEL_FASTA",
    "full": "RTM_MEI_FULL_PANEL_FASTA",
}

PANEL_KINDS = (
    ("fragment", PUBLIC_MEI_REMAP_FASTA_RELPATHS[0]),
    ("full", PUBLIC_MEI_REMAP_FASTA_RELPATHS[1]),
)


def full_axis_name(label: str) -> str:
    """``L1HS_5end`` and ``L1HS`` both name the ``L1HS_full`` consensus."""
    base = (label or "").strip()
    for suffix in FRAGMENT_SUFFIXES + (FULL_SUFFIX,):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return f"{base}{FULL_SUFFIX}"


def read_panel_labels(path: str | Path) -> dict[str, tuple[str, int]]:
    """Map each panel record name to ``(class/family, sequence_bp)``.

    Headers are ``>name#class/family`` (optionally with ``|source=...`` appended,
    as the derived full panel writes them); the sequence length is counted from
    the bases rather than taken on trust.
    """
    labels: dict[str, tuple[str, int]] = {}
    name: str | None = None
    mei_class = ""
    length = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    labels[name] = (mei_class, length)
                head = line[1:].strip().split("|", 1)[0]
                name, _, mei_class = head.partition("#")
                length = 0
            elif name:
                length += len(line.strip())
    if name:
        labels[name] = (mei_class, length)
    return labels


def find_panel_paths(
    explicit: "Sequence[str | Path]" = (),
    *,
    public_data_dir: str | Path | None = None,
) -> list[tuple[str, Path]]:
    """Locate the MEI panels; return ``[(kind, path)]`` for those that exist.

    Search order per panel: an explicit path, its environment variable, then
    ``RTM_PUBLIC_DATA_DIR`` or ``public_data_dir`` or the working directory, each
    joined with the relpath this package already publishes.
    """
    explicit = [Path(path) for path in explicit]
    roots = []
    if public_data_dir:
        roots.append(Path(public_data_dir))
    env_root = os.environ.get("RTM_PUBLIC_DATA_DIR")
    if env_root:
        roots.append(Path(env_root))
    roots.append(Path.cwd())

    found: list[tuple[str, Path]] = []
    for kind, relpath in PANEL_KINDS:
        candidates: list[Path] = []
        if kind == "fragment" and explicit:
            candidates.extend(explicit)
        env_path = os.environ.get(PANEL_ENV_VARS[kind])
        if env_path:
            candidates.append(Path(env_path))
        candidates.extend(root / relpath for root in roots)
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_size > 0:
                found.append((kind, candidate))
                break
    return found


def resolve_consensus_lengths(
    labels: "Iterable[str]",
    panels: "Sequence[tuple[str, Path]]" = (),
) -> dict[str, ConsensusLength | None]:
    """Length on the consensus axis for each label, or ``None`` when unknown.

    The full axis is tried first, by ``full_axis_name(label)`` and then by the
    label itself; the fragment axis is used only for a label that exists nowhere
    on the full axis.  ``None`` means *this panel does not know*, which is a
    different statement from a length of zero -- callers should keep their
    data-derived fallback and say which labels used it.
    """
    by_kind: dict[str, dict[str, tuple[str, int]]] = {}
    for kind, path in panels:
        by_kind.setdefault(kind, {}).update(read_panel_labels(path))

    full_panel = by_kind.get("full", {})
    fragment_panel = by_kind.get("fragment", {})

    def _full(name: str, basis: str) -> ConsensusLength | None:
        if name not in full_panel:
            return None
        mei_class, bp = full_panel[name]
        return ConsensusLength(bp=bp, target=name, axis="full", panel="full",
                               basis=basis, mei_class=mei_class)

    resolved: dict[str, ConsensusLength | None] = {}
    for label in labels:
        mei_class = fragment_panel.get(label, ("", 0))[0]
        hit = (_full(full_axis_name(label), "base_full")
               or _full(label, "base_full")
               or _full(FAMILY_CANONICAL_FULL.get(family_of_label(label, mei_class), ""),
                        "family_canonical"))
        if hit is None and label in fragment_panel:
            frag_class, bp = fragment_panel[label]
            hit = ConsensusLength(bp=bp, target=label, axis="fragment", panel="fragment",
                                  basis="fragment", mei_class=frag_class)
        resolved[label] = hit
    return resolved


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
