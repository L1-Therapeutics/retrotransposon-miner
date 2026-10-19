#!/usr/bin/env python3
"""Build ``dfam_human_mei_l1_alu_sva.fasta`` from Dfam's RepeatMasker library.

Why this exists
---------------
``download_public_data.py`` prepares the human MEI panel in
``_prepare_dfam_mei_library``: download the FamDB archive
(``dfam40.0.h5.gz``), decompress it, ``git clone`` the FamDB repository, run
``famdb.py families -ad 9606 --curated``, then filter to L1/Alu/SVA.  That path
needs ``h5py``/``pyhmmer`` and a multi-gigabyte download, so on a machine
without them the step fails and every consumer of the panel records the file as
*absent*:

* ``src/retro_miner/mei_panel_index.py`` lists it as the primary bwa remap
  target (and skips it silently when missing),
* ``download_public_data.py`` builds the full-consensus panel and the
  fragment-to-full coordinate map from it, both skipping when it is absent,
* ``span_work/line1_span_support.py`` estimates each subfamily's consensus
  length as the largest ``MEI_3P`` seen in the cohort instead of reading it,
* ``span_work/consensus_coverage.py`` bounds the consensus axis from below by
  the furthest read coordinate instead of knowing where the axis ends.

This script produces the same file from a source that is one HTTP request
instead of a pipeline: ``Dfam-RepeatMasker.lib.gz`` from Dfam's release
directory (15.8 MB, and Dfam publishes its md5 next to it).  It is the
RepeatMasker-format library built from the same Dfam release, with headers
``>name#class/family @clade [flags]``, so it carries the same fragment names the
pipeline already emits (``L1HS_5end``, ``L1PA2_3end``, ``AluYa5``, ``SVA_F``,
``L1M1_orf2``, ...) and the same sequences.

What "human" means here
-----------------------
The RepeatMasker library is pan-clade: of its 26,292 records, 2,454 are
L1/Alu/SVA but most belong to other mammals and fish.  A record is kept when
its clade tag names a node on the human lineage (Homo_sapiens up through
Eutheria, plus ``Theria_mammals``, which is the tag Dfam uses for the L1M*
families that do occur in hg38).  The filter is checked against hg38 rather than
assumed: every one of the 192 distinct L1/Alu/SVA ``repName`` values in this
project's hg38 RepeatMasker track has a human-lineage tag, so no name hg38
carries is dropped by the allowlist.

Output format
-------------
``>{name}#{class}/{family}`` with the sequence wrapped at 60 columns and
terminal poly(A) trimmed by the same function the downloader's subset filter
uses, so this file is a drop-in replacement for the famdb-derived one in every
consumer above.

Usage
-----
    python scripts/fetch_dfam_mei_panel.py --out-dir /path/to/public
    python scripts/fetch_dfam_mei_panel.py --verify /path/to/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import time
import urllib.request
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

# The downloader's own poly(A) trim, reused rather than reimplemented: the panel
# this script writes has to sit on the same coordinate axis as the panel the
# downloader would have written.
from download_public_data import _trim_mei_consensus_terminal_polya  # noqa: E402

DEFAULT_LIBRARY_URL = "https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz"
PANEL_RELPATH = Path("retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta")
PANEL_MANIFEST_SUFFIX = ".manifest.json"

#: Classes kept.  ``Retroposon/SVA`` is matched by name below, exactly as
#: ``_export_human_mei_subset`` does, so the two producers agree on membership.
KEPT_CLASSES = ("LINE/L1", "SINE/Alu")
SVA_CLASS = "Retroposon/SVA"

#: Clade tags naming a node on the human lineage, root-most last.  ``Eutheria``
#: is the floor that keeps ``HAL1b``; ``Theria_mammals`` is spelled the way Dfam
#: spells it (the L1M* families hg38 contains carry that tag).
HUMAN_LINEAGE_CLADES = (
    "Homo_sapiens",
    "Homininae",
    "Hominidae",
    "Hominoidea",
    "Catarrhini",
    "Simiiformes",
    "Haplorrhini",
    "Primates",
    "Euarchonta",
    "Euarchontoglires",
    "Boreoeutheria",
    "Eutheria",
    "Theria_mammals",
    "Mammalia",
    "Amniota",
    "Vertebrata_vertebrates",
)

HEADER_RE = re.compile(r"^>(?P<name>\S+?)#(?P<cls>[^ @\[]+)(?P<rest>.*)$")
CLADE_RE = re.compile(r"@(\S+)")

__all__ = [
    "DEFAULT_LIBRARY_URL",
    "PANEL_RELPATH",
    "HUMAN_LINEAGE_CLADES",
    "PanelRecord",
    "parse_repeatmasker_library",
    "select_human_mei",
    "write_panel",
    "build_panel",
    "main",
]


@dataclass(frozen=True)
class PanelRecord:
    """One library record: a Dfam consensus or consensus fragment."""

    name: str
    cls: str
    clade: str
    seq: str

    @property
    def header(self) -> str:
        return f">{self.name}#{self.cls}"

    @property
    def trimmed(self) -> str:
        return _trim_mei_consensus_terminal_polya(self.seq)


def _read_text(path: Path) -> str:
    """Read plain or gzipped text (the library ships gzipped)."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8", errors="replace")


def parse_repeatmasker_library(text: str) -> list[PanelRecord]:
    """Parse RepeatMasker library text into records, in file order.

    Header form is ``>name#class/family @clade [fields]``.  A header without a
    ``#class/family`` part is dropped: every consumer of this panel splits the
    header on ``#``, so a record that cannot be split is not usable.
    """
    records: list[PanelRecord] = []
    name = cls = clade = None
    parts: list[str] = []

    def flush() -> None:
        if name is not None and cls is not None:
            clade_match = CLADE_RE.search(clade or "")
            records.append(
                PanelRecord(name=name, cls=cls,
                            clade=clade_match.group(1) if clade_match else "",
                            seq="".join(parts).upper().replace("U", "T"))
            )

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            match = HEADER_RE.match(line.strip())
            name = match.group("name") if match else None
            cls = match.group("cls") if match else None
            clade = match.group("rest") if match else ""
            parts = []
        elif name is not None:
            parts.append(line.strip())
    flush()
    return records


def select_human_mei(
    records: list[PanelRecord],
    *,
    clades: tuple[str, ...] = HUMAN_LINEAGE_CLADES,
    kept_classes: tuple[str, ...] = KEPT_CLASSES,
) -> tuple[list[PanelRecord], dict[str, object]]:
    """Keep the L1/Alu/SVA records whose clade is on the human lineage.

    Returns the selection (sorted by class then name, so the output hash is
    reproducible) and the counts needed to audit what was dropped.
    """
    allowed = set(clades)
    kept: list[PanelRecord] = []
    meis = [record for record in records
            if record.cls in kept_classes or record.cls == SVA_CLASS]
    dropped_clades: Counter[str] = Counter()
    for record in meis:
        if record.clade in allowed:
            kept.append(record)
        else:
            dropped_clades[record.clade or "(no clade)"] += 1
    kept.sort(key=lambda record: (record.cls, record.name))
    stats: dict[str, object] = {
        "library_records": len(records),
        "mei_records_in_library": len(meis),
        "kept": len(kept),
        "kept_by_class": dict(sorted(Counter(r.cls for r in kept).items())),
        "kept_by_clade": dict(sorted(Counter(r.clade for r in kept).items())),
        "dropped_non_human_clades": len(dropped_clades),
        "dropped_clade_examples": dict(dropped_clades.most_common(5)),
    }
    return kept, stats


def write_panel(records: list[PanelRecord], out_path: Path, *, width: int = 60) -> dict[str, object]:
    """Write the panel FASTA; return its sha256 and the per-record lengths."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, object]] = []
    with out_path.open("w", encoding="utf-8") as handle:
        for record in records:
            body = record.trimmed
            if not body:
                continue
            handle.write(record.header + "\n")
            for start in range(0, len(body), width):
                handle.write(body[start : start + width] + "\n")
            entries.append({
                "name": record.name,
                "class": record.cls,
                "clade": record.clade,
                "library_bp": len(record.seq),
                "panel_bp": len(body),
                "polya_trimmed_bp": len(record.seq) - len(body),
                "seq_sha256": hashlib.sha256(body.encode()).hexdigest(),
            })
    digest = hashlib.sha256(out_path.read_bytes()).hexdigest()
    return {"fasta_sha256": digest, "fasta_bytes": out_path.stat().st_size, "entries": entries}


def _download(url: str, dest: Path, *, timeout: int = 600, force: bool = False) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=timeout) as response, tmp.open("wb") as handle:
        while chunk := response.read(1 << 20):
            handle.write(chunk)
    tmp.replace(dest)


def _published_md5(url: str, *, timeout: int = 60) -> str:
    """Dfam publishes ``<file>.md5`` as ``<hex>  <filename>``."""
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode().split()[0].strip().lower()


def _md5(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324  (Dfam publishes md5; this checks their file)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_panel(
    out_path: Path,
    *,
    source: Path | None = None,
    url: str = DEFAULT_LIBRARY_URL,
    force: bool = False,
    clades: tuple[str, ...] = HUMAN_LINEAGE_CLADES,
    timeout: int = 600,
) -> dict[str, object]:
    """Download (or reuse) the library, select the human MEI records, write the panel.

    The library's md5 is checked against the one Dfam publishes before anything
    is read from it, and both hashes land in the manifest: the input hash pins
    which release of the source was used, the output hash pins the panel every
    consumer reads.
    """
    started = time.time()
    manifest_path = out_path.with_name(out_path.name + PANEL_MANIFEST_SUFFIX)
    published = None
    try:
        published = _published_md5(f"{url}.md5")
    except Exception as exc:  # noqa: BLE001  (a missing .md5 must not stop the build)
        published_md5_error = f"{type(exc).__name__}: {exc}"
    else:
        published_md5_error = ""

    if source is None:
        source = out_path.parent / Path(url).name
        _download(url, source, timeout=timeout, force=force)
    source_md5 = _md5(source)
    if published and source_md5 != published:
        raise SystemExit(
            f"{source}: md5 {source_md5} does not match the published {published}; "
            "refusing to build a panel from an unverified library"
        )

    records = parse_repeatmasker_library(_read_text(source))
    kept, selection = select_human_mei(records, clades=clades)
    written = write_panel(kept, out_path)

    manifest: dict[str, object] = {
        "panel": str(out_path),
        **written,
        **selection,
        "source_url": url,
        "source_path": str(source),
        "source_bytes": source.stat().st_size,
        "source_md5": source_md5,
        "source_md5_published": published or "",
        "source_md5_verified": bool(published) and source_md5 == published,
        "source_md5_error": published_md5_error,
        "human_lineage_clades": list(clades),
        "kept_classes": list(KEPT_CLASSES),
        "header_format": ">{name}#{class}/{family}",
        "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "built_by": "scripts/fetch_dfam_mei_panel.py",
        "build_seconds": round(time.time() - started, 1),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
    manifest["manifest"] = str(manifest_path)
    return manifest


def build_full_panel(out_root: Path, *, timeout: int = 600, force: bool = False) -> dict[str, object]:
    """Derive the full-length canonical panel the coordinate axis actually uses.

    The caller's ``MEI_5P``/``MEI_3P`` coordinates are *not* on the fragment axis
    this script writes: a ``L1HS_5end`` call reports ``MEI_3P`` up to 6018 while
    its fragment consensus is 2136 bp.  They are on the full-length axis, which
    the downloader derives from this panel plus UCSC ``hg38reps.fa`` (one 1.7 MB
    file).  Calling that derivation here rather than reimplementing it keeps both
    producers of the panel in agreement; the axis lengths the span and coverage
    scripts read come from its output.
    """
    from download_public_data import _build_mei_full_consensus_panel

    return _build_mei_full_consensus_panel(
        outdir=out_root, ds_map={}, timeout_sec=timeout, force=force
    )


def verify_panel(out_path: Path) -> int:
    """Re-check an existing panel against its manifest: sizes and hashes."""
    manifest_path = out_path.with_name(out_path.name + PANEL_MANIFEST_SUFFIX)
    if not manifest_path.exists():
        print(f"no manifest at {manifest_path}", file=sys.stderr)
        return 2
    manifest = json.loads(manifest_path.read_text())
    actual = hashlib.sha256(out_path.read_bytes()).hexdigest()
    expected = str(manifest.get("fasta_sha256", ""))
    ok = actual == expected
    print(f"{out_path}")
    print(f"  sha256 {actual}")
    print(f"  pinned {expected}  -> {'MATCH' if ok else 'MISMATCH'}")
    print(f"  records {manifest.get('kept')} | bytes {out_path.stat().st_size} "
          f"| source md5 verified {manifest.get('source_md5_verified')}")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("."),
                        help="public-data root; the panel lands at "
                             f"{PANEL_RELPATH} under it")
    parser.add_argument("--out", type=Path, default=None,
                        help="explicit panel path (overrides --out-dir)")
    parser.add_argument("--source", type=Path, default=None,
                        help="use this local copy of Dfam-RepeatMasker.lib.gz "
                             "instead of downloading it")
    parser.add_argument("--url", default=DEFAULT_LIBRARY_URL)
    parser.add_argument("--force", action="store_true",
                        help="re-download even when the library is already present")
    parser.add_argument("--verify", type=Path, default=None,
                        help="verify an existing panel against its manifest and exit")
    parser.add_argument("--no-full-panel", action="store_true",
                        help="skip deriving the full-length canonical panel (the axis "
                             "the caller's MEI_5P/MEI_3P coordinates live on)")
    args = parser.parse_args(argv)

    if args.verify is not None:
        return verify_panel(args.verify)

    out_path = args.out or (args.out_dir / PANEL_RELPATH)
    manifest = build_panel(out_path, source=args.source, url=args.url, force=args.force)
    print(f"panel {manifest['panel']}")
    print(f"  library records {manifest['library_records']} | MEI in library "
          f"{manifest['mei_records_in_library']} | kept {manifest['kept']} "
          f"{manifest['kept_by_class']}")
    print(f"  kept by clade {manifest['kept_by_clade']}")
    print(f"  dropped (non-human clades) {manifest['dropped_non_human_clades']} "
          f"e.g. {manifest['dropped_clade_examples']}")
    print(f"  source md5 {manifest['source_md5']} published "
          f"{manifest['source_md5_published'] or '(unavailable)'} "
          f"verified {manifest['source_md5_verified']}")
    print(f"  panel sha256 {manifest['fasta_sha256']} ({manifest['fasta_bytes']} bytes)")
    print(f"  manifest {manifest['manifest']}")

    if not args.no_full_panel:
        full = build_full_panel(args.out_dir, timeout=600, force=args.force)
        print(f"\nfull-consensus panel {full.get('panel_fasta')}")
        print(f"  status {full.get('status')} | records {full.get('selected_total')} "
              f"(L1 {full.get('selected_line1')}, Alu {full.get('selected_alu')}, "
              f"SVA {full.get('selected_sva')}) | samtools faidx {full.get('panel_faidx')}")
        print("  this is the axis the MEI_5P/MEI_3P coordinates are on; the fragment "
              "panel above is what the remap aligns to")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
