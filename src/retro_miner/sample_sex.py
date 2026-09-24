"""Decide whether chrY should run from whole-chromosome index counts.

``samtools idxstats`` reads the BAM index, not the reads. Sex is the chrX
mapped-reads-per-base divided by the same value on chr1. Males are near 0.5
and females are near 1.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Halfway between the male cluster near 0.5 and the female cluster near 1.
MALE_X_TO_AUTOSOME_MAX = 0.75


def _canonical_chrom(name: str) -> str:
    if name == "*":
        return name
    if name.startswith("chr"):
        return name
    if name in {"X", "Y", "M", "MT"} or name.isdigit():
        return f"chr{name}"
    return name


def mapped_per_base(idxstats_text: str) -> dict[str, float]:
    """Map chromosome name to mapped reads divided by chromosome length."""
    density: dict[str, float] = {}
    for line in idxstats_text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chrom = _canonical_chrom(parts[0])
        try:
            length = int(parts[1])
            mapped = int(parts[2])
        except ValueError:
            continue
        if length <= 0 or chrom in density:
            continue
        density[chrom] = mapped / length
    return density


def x_to_autosome_ratio(idxstats_text: str, autosome: str = "chr1") -> float | None:
    density = mapped_per_base(idxstats_text)
    auto = density.get(autosome)
    chrom_x = density.get("chrX")
    if auto is None or chrom_x is None or auto <= 0:
        return None
    return chrom_x / auto


def classify_sex(idxstats_text: str, autosome: str = "chr1") -> tuple[str, float | None]:
    """Return ``(male|female|unknown, ratio)``."""
    ratio = x_to_autosome_ratio(idxstats_text, autosome=autosome)
    if ratio is None:
        return "unknown", None
    if ratio < MALE_X_TO_AUTOSOME_MAX:
        return "male", ratio
    return "female", ratio


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        print("usage: python -m retro_miner.sample_sex IDXSTATS_TSV", file=sys.stderr)
        return 2
    text = Path(args[0]).read_text(encoding="utf-8")
    label, ratio = classify_sex(text)
    if ratio is None:
        print(label)
    else:
        print(f"{label} {ratio:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
