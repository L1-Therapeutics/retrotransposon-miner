"""Spatial Index-Guided Interval Streaming Engine for BAM Evidence Fetching.

Provides O(log N) spatial index queries against BAM/CSI-indexed alignment files
to extract per-locus evidence reads without full-file sequential iteration.

Literature Anchors: Li et al. (2009) SAMtools BAI/CSI / Cameron et al. (2017) GRIDSS.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .bam_io import open_alignment


@dataclass(frozen=True)
class LocusEvidence:
    """Container for evidence reads extracted from a single genomic locus interval.

    Attributes:
        chrom: Reference contig name.
        start: 0-based interval start.
        end: 0-based interval end (exclusive).
        spanning_pairs: Count of concordant properly-paired reference reads
            fully spanning the interval (used as VAF denominator).
        soft_clips: Soft-clipped query sequences anchored at the locus.
        discordant_mates: Mate pairs mapping to distinct reference contigs.
        raw_reads: All raw pysam read records returned by the indexed fetch,
            preserved for downstream consumer processing.
    """

    chrom: str
    start: int
    end: int
    spanning_pairs: int = 0
    soft_clips: list[str] = field(default_factory=list)
    discordant_mates: list[str] = field(default_factory=list)
    raw_reads: list[pysam.AlignedSegment] = field(default_factory=list)


def has_bam_index(bam_path: Path) -> bool:
    """Return True if a BAM index file exists alongside *bam_path*."""
    bai = Path(str(bam_path) + ".bai")
    if bai.exists():
        return True
    csi = Path(str(bam_path) + ".csi")
    return csi.exists()


def _classify_read(read: pysam.AlignedSegment, start: int, end: int, min_mapq: int = 20) -> tuple[bool, bool, bool]:
    """Classify a single read into spanning / soft-clip / discordant categories.

    Returns:
        Tuple of ``(is_spanning, is_soft_clip, is_discordant)`` booleans.
    """
    if read.is_unmapped or read.mapping_quality < min_mapq:
        return False, False, False

    ref_start = int(read.reference_start) if read.reference_start is not None else -1
    ref_end = int(read.reference_end) if read.reference_end is not None else -1

    spans = ref_start >= 0 and ref_end >= 0 and ref_start <= end and ref_end >= start

    has_soft_clip = False
    if read.cigartuples is not None:
        has_soft_clip = any(op == 4 for op, _ in read.cigartuples)

    is_discordant = False
    if read.is_paired and not read.mate_is_unmapped:
        if read.reference_id != read.next_reference_id:
            is_discordant = True
        elif abs(read.template_length) >= 1000:
            is_discordant = True

    return bool(spans), has_soft_clip, is_discordant


def fetch_locus_spanning_pairs(
    bam_path: Path,
    chrom: str,
    start: int,
    end: int,
    min_mapq: int = 20,
    max_reads: int = 10000,
) -> LocusEvidence:
    """Fetch evidence reads for a single candidate locus using spatial index.

    Opens *bam_path* with its accompanying index (``.bai`` or ``.csi``) and
    performs an ``HTSlib``-native interval fetch.  Reads are classified into
    spanning concordant pairs, soft-clipped split anchors, and discordant
    mate pairs.

    Args:
        bam_path: Path to sorted BAM file.  An index file with the same
            basename and ``.bai`` or ``.csi`` extension must exist.
        chrom: Reference contig name.
        start: 0-based start coordinate of the candidate locus.
        end: 0-based end coordinate of the candidate locus (exclusive).
        min_mapq: Minimum mapping quality for included reads.
        max_reads: Safety cap on returned raw reads to bound memory.

    Returns:
        LocusEvidence with classified read counts and sequences.

    Raises:
        FileNotFoundError: If *bam_path* or its index does not exist.
        ValueError: If the index file cannot be opened by pysam.
    """
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")
    if not has_bam_index(bam_path):
        raise FileNotFoundError(
            f"BAM index not found for {bam_path}; expected {bam_path}.bai or {bam_path}.csi"
        )

    spanning = 0
    soft_clips: list[str] = []
    discordant: list[str] = []
    raw: list[pysam.AlignedSegment] = []

    with open_alignment(bam_path) as bam:
        try:
            iterator = bam.fetch(chrom, max(0, start), max(start + 1, end))
        except (ValueError, KeyError) as exc:
            raise ValueError(f"Failed indexed fetch for {chrom}:{start}-{end}: {exc}") from exc

        for read in iterator:
            if len(raw) >= max_reads:
                break
            raw.append(read)

            is_spanning, is_soft_clip, is_discordant = _classify_read(read, start, end, min_mapq)
            if is_spanning:
                spanning += 1
            if is_soft_clip:
                seq = read.query_sequence or ""
                if read.cigartuples is not None:
                    first_op, first_len = read.cigartuples[0]
                    last_op, last_len = read.cigartuples[-1]
                    if first_op == 4 and first_len > 0:
                        soft_clips.append(seq[:first_len])
                    elif last_op == 4 and last_len > 0:
                        soft_clips.append(seq[-last_len:])
            if is_discordant:
                discordant.append(read.query_name or "")

    return LocusEvidence(
        chrom=chrom,
        start=start,
        end=end,
        spanning_pairs=spanning,
        soft_clips=soft_clips,
        discordant_mates=discordant,
        raw_reads=raw,
    )
