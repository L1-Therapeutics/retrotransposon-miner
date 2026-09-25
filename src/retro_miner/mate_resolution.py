"""Resolve full-genome BAMs/CRAMs for interchromosomal discordant mates.

Chr-sliced test BAMs still record off-chromosome ``mate_chrom``/``mate_pos``
on the anchor, but the mate itself is not in the slice. Annotate then silently
drops those sequences unless ``--disease-mate-bam`` points at the WGS file.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Chr-sliced public test BAMs whose WGS object lives under bam_stage / S3 full/.
_SLICE_BASENAME_TO_FULL: dict[str, tuple[str, str]] = {
    "hg00100.shortread.chr22.hg38.bam": (
        "HG00100.final.cram",
        "hg00100_shortread_highcov_cram",
    ),
}

# Fail annotate when this many off-chrom discordants still lack mate_seq.
MIN_EMPTY_INTERCHROM_MATES_TO_FAIL = 10


@dataclass(frozen=True)
class InterchromMateSeqStats:
    off_chrom: int
    off_empty: int
    off_with_seq: int


def _workdir(explicit: str | Path | None = None) -> Path:
    if explicit:
        return Path(explicit).expanduser()
    env = (os.environ.get("RTM_WORKDIR") or "").strip()
    if env:
        return Path(env).expanduser()
    return Path.home() / "retrotransposon-workdir"


def _public_data_dir(workdir: Path) -> Path:
    env = (os.environ.get("RTM_PUBLIC_DATA_DIR") or "").strip()
    if env:
        return Path(env).expanduser()
    return workdir / "data" / "public"


def _bam_stage_dir(workdir: Path) -> Path:
    env = (os.environ.get("RTM_BAM_STAGE_DIR") or "").strip()
    if env:
        return Path(env).expanduser()
    return workdir / "data" / "bam_stage"


def infer_full_genome_mate_bam(
    scan_bam: str | Path | None,
    *,
    workdir: str | Path | None = None,
) -> Path | None:
    """Return a local WGS BAM/CRAM for a known chromosome-sliced scan BAM."""
    env = (os.environ.get("RTM_DISEASE_MATE_BAM") or os.environ.get("RTM_MATE_BAM") or "").strip()
    if env:
        path = Path(env).expanduser()
        if path.is_file():
            return path
    if scan_bam is None:
        return None
    name = Path(scan_bam).name
    mapping = _SLICE_BASENAME_TO_FULL.get(name)
    if mapping is None:
        return None
    full_name, dataset_id = mapping
    root = _workdir(workdir)
    candidates = (
        _bam_stage_dir(root) / full_name,
        _public_data_dir(root) / "test_data" / "full" / dataset_id / full_name,
    )
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def interchrom_mate_seq_stats(discordant_df: pd.DataFrame) -> InterchromMateSeqStats:
    if discordant_df is None or discordant_df.empty:
        return InterchromMateSeqStats(0, 0, 0)
    if "mate_chrom" not in discordant_df.columns or "chrom" not in discordant_df.columns:
        return InterchromMateSeqStats(0, 0, 0)
    chrom = discordant_df["chrom"].fillna("").astype(str)
    mate_chrom = discordant_df["mate_chrom"].fillna("").astype(str)
    off = (mate_chrom != "") & (mate_chrom != "*") & (mate_chrom != chrom)
    if "mate_seq" in discordant_df.columns:
        seq_len = discordant_df["mate_seq"].fillna("").astype(str).str.len()
    else:
        seq_len = pd.Series(0, index=discordant_df.index)
    off_n = int(off.sum())
    off_empty = int((off & seq_len.eq(0)).sum())
    return InterchromMateSeqStats(
        off_chrom=off_n,
        off_empty=off_empty,
        off_with_seq=off_n - off_empty,
    )


def missing_interchrom_mates_message(
    stats: InterchromMateSeqStats,
    *,
    mate_bam: str | Path | None = None,
) -> str:
    mate = str(mate_bam) if mate_bam else "(scan BAM / none)"
    return (
        f"{stats.off_empty}/{stats.off_chrom} interchromosomal discordant mates have empty "
        f"mate_seq after fetch from {mate}. The scan BAM is a chromosome slice, so those "
        "mates live on other chromosomes. Pass --disease-mate-bam / --control-mate-bam "
        "to the full-genome BAM or CRAM (for HG00100: bam_stage/HG00100.final.cram or "
        "s3://l1tx-data/public/test_data/full/hg00100_shortread_highcov_cram/), or rebuild "
        "the chr slice with include_discordant_mates. "
        "Override with --allow-missing-interchrom-mates only if this is intentional."
    )


def require_interchrom_mate_sequences(
    discordant_df: pd.DataFrame,
    *,
    mate_bam: str | Path | None = None,
    allow_missing: bool = False,
    min_empty: int = MIN_EMPTY_INTERCHROM_MATES_TO_FAIL,
) -> InterchromMateSeqStats:
    """Raise if off-chrom discordants still have no mate sequence after fetch."""
    import click

    stats = interchrom_mate_seq_stats(discordant_df)
    if allow_missing or stats.off_empty < int(min_empty):
        return stats
    raise click.ClickException(missing_interchrom_mates_message(stats, mate_bam=mate_bam))
