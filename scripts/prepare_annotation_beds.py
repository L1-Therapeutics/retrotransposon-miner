#!/usr/bin/env python3
"""Convert downloaded UCSC annotation tables into BED tracks.

Materializes the BED tracks consumed by the discovery/annotation pipeline from
the raw UCSC table dumps under the hg38 annotation directory:

* ``segdup/genomicSuperDups.txt.gz`` -> ``segdup/genomicSuperDups.bed``
  4 columns: ``chrom  chromStart  chromEnd  otherChrom``
* ``repeats/rmsk.txt.gz`` -> ``repeats/rmsk.bed``
  6 columns: ``genoName  genoStart  genoEnd  repName  repClass  repFamily``

Column indices are 0-based and match the UCSC table schemas (and the layout the
nested-rmsk annotation already expects in ``mei_support._write_rmsk_mei_bed``).
Existing non-empty target BEDs are skipped unless ``--force`` is given.
"""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

DEFAULT_OUTDIR = Path.home() / "retrotransposon-workdir" / "data" / "public" / "annotation" / "hg38"


@dataclass(frozen=True)
class BedSpec:
    """One UCSC ``.txt.gz`` -> ``.bed`` conversion job.

    ``columns`` is a tuple of ``(output_name, source_column_index)`` pairs in
    output order. By convention the 2nd and 3rd entries are the 0-based
    start/end coordinates used for row validation.
    """

    name: str
    src_rel: str
    dst_rel: str
    columns: tuple[tuple[str, int], ...]


SEG_DUP_SPEC = BedSpec(
    name="segdup",
    src_rel="segdup/genomicSuperDups.txt.gz",
    dst_rel="segdup/genomicSuperDups.bed",
    columns=(
        ("chrom", 1),
        ("chromStart", 2),
        ("chromEnd", 3),
        ("otherChrom", 7),
    ),
)

RMSK_SPEC = BedSpec(
    name="rmsk",
    src_rel="repeats/rmsk.txt.gz",
    dst_rel="repeats/rmsk.bed",
    columns=(
        ("genoName", 5),
        ("genoStart", 6),
        ("genoEnd", 7),
        ("repName", 10),
        ("repClass", 11),
        ("repFamily", 12),
    ),
)

SPECS: tuple[BedSpec, ...] = (SEG_DUP_SPEC, RMSK_SPEC)


def _convert_one(outdir: Path, spec: BedSpec, force: bool) -> dict[str, object]:
    """Convert a single spec's source table to a BED file. Returns a status dict."""
    src = outdir / spec.src_rel
    dst = outdir / spec.dst_rel

    if dst.exists() and dst.stat().st_size > 0 and not force:
        return {"name": spec.name, "status": "skipped_exists", "path": str(dst), "bytes": dst.stat().st_size}
    if not src.exists():
        return {"name": spec.name, "status": "missing_src", "path": str(src)}

    max_idx = max(idx for _, idx in spec.columns)
    start_idx = spec.columns[1][1]
    end_idx = spec.columns[2][1]

    dst.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    skipped = 0
    with gzip.open(src, "rt", encoding="utf-8") as hin, dst.open("w", encoding="utf-8") as hout:
        for line in hin:
            line = line.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max_idx:
                skipped += 1
                continue
            try:
                start_i = int(parts[start_idx])
                end_i = int(parts[end_idx])
            except ValueError:
                skipped += 1
                continue
            if end_i <= start_i:
                skipped += 1
                continue
            fields = [parts[idx] for _, idx in spec.columns]
            if not fields[0]:
                skipped += 1
                continue
            hout.write("\t".join(fields) + "\n")
            rows += 1

    return {"name": spec.name, "status": "converted", "path": str(dst), "rows": rows, "skipped": skipped}


def _print_result(result: dict[str, object]) -> None:
    name = str(result["name"])
    status = str(result["status"])
    path = str(result["path"])
    if status == "skipped_exists":
        print(f"[{name}] skipped (exists, non-empty): {path} ({result['bytes']} bytes)")
    elif status == "missing_src":
        print(f"[{name}] WARNING: source not found: {path}")
    else:
        print(f"[{name}] converted {path} (rows={result['rows']} skipped={result['skipped']})")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert UCSC segdup/rmsk annotation tables to BED tracks.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="hg38 annotation directory (default: %(default)s)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-convert even if a non-empty target BED already exists",
    )
    args = parser.parse_args()

    outdir = args.outdir
    missing = False
    for spec in SPECS:
        result = _convert_one(outdir, spec, args.force)
        _print_result(result)
        if result["status"] == "missing_src":
            missing = True

    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
