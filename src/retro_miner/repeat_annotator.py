"""Spatial Interval Tree Annotator for Reference Repeat Overlap.

Provides O(log N) nested containment list (NClist) interval tree queries
against RepeatMasker BED annotations to flag candidate MEI loci that overlap
pre-existing reference repeats.

Literature Anchors: Szak et al. (2002) Genome Res / Li (2011) SAMtools.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from intervaltree import IntervalTree


@dataclass(frozen=True)
class RepeatAnnotation:
    """Single RepeatMasker annotation record.

    Attributes:
        chrom: Reference contig name.
        start: 0-based start coordinate.
        end: 0-based exclusive end coordinate.
        name: Repeat consensus name (e.g. ``AluY``, ``L1HS``).
        repeat_class: RepeatMasker class (e.g. ``SINE``, ``LINE``, ``SVA``).
        family: Repeat family (e.g. ``Alu``, ``LINE1``).
        strand: Strand (``+``, ``-``, or ``.``).
        divergence: Divergence score from consensus (percent or raw).
        length: Length in bp.
    """

    chrom: str
    start: int
    end: int
    name: str
    repeat_class: str
    family: str
    strand: str = "."
    divergence: float = 0.0
    length: int = 0


class RepeatIntervalTree:
    """O(log N) spatial index for RepeatMasker reference annotations.

    Wraps a per-chromosome ``IntervalTree`` with RepeatMasker-specific
    parsing and query helpers.  Build once, query many candidate loci.

    Attributes:
        trees: Mapping of chromosome name to ``IntervalTree``.
        annotations: Flat list of loaded :class:`RepeatAnnotation` records.
    """

    def __init__(self) -> None:
        self.trees: dict[str, IntervalTree] = {}
        self.annotations: list[RepeatAnnotation] = []

    def build_from_bed(self, bed_path: str | Path) -> None:
        """Load RepeatMasker BED regions into the interval tree.

        Supports standard BED3/BED6/BED9 formats plus UCSC rmsk-derived
        BEDs with extra repeat metadata columns.

        Args:
            bed_path: Path to BED file (plain or ``.gz``).
        """
        bed_path = Path(bed_path)
        opener = gzip.open if str(bed_path).endswith(".gz") else open
        with opener(bed_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                if end <= start:
                    continue

                name = parts[3] if len(parts) > 3 else "."
                strand = "."
                repeat_class = ""
                family = ""
                divergence = 0.0
                length = end - start

                if len(parts) >= 6:
                    strand = parts[5]
                if len(parts) >= 7:
                    repeat_class = parts[6]
                if len(parts) >= 8:
                    family = parts[7]
                if len(parts) >= 9:
                    try:
                        divergence = float(parts[8])
                    except ValueError:
                        divergence = 0.0

                ann = RepeatAnnotation(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    repeat_class=repeat_class,
                    family=family,
                    strand=strand,
                    divergence=divergence,
                    length=length,
                )
                self.annotations.append(ann)
                tree = self.trees.setdefault(chrom, IntervalTree())
                tree.addi(start, end, ann)

    def query_overlap(
        self,
        chrom: str,
        pos: int,
        buffer_bp: int = 50,
    ) -> list[dict[str, Any]]:
        """Query overlapping repeats for a candidate position in O(log N).

        Args:
            chrom: Reference contig name.
            pos: Candidate breakpoint coordinate.
            buffer_bp: Expansion around *pos* used as interval query bounds.

        Returns:
            List of dicts with overlap metadata.  Empty list when no
            overlapping repeats are found.
        """
        tree = self.trees.get(chrom)
        if tree is None:
            return []

        query_start = max(0, pos - buffer_bp)
        query_end = pos + buffer_bp
        hits = list(tree.at(pos))
        if not hits:
            return []

        results: list[dict[str, Any]] = []
        for iv in hits:
            ann = iv.data
            overlap_bp = max(0, min(ann.end, query_end) - max(ann.start, query_start))
            if overlap_bp == 0 and not (ann.start <= pos < ann.end):
                continue
            if overlap_bp == 0:
                overlap_bp = 1
            results.append(
                {
                    "chrom": ann.chrom,
                    "start": ann.start,
                    "end": ann.end,
                    "name": ann.name,
                    "repeat_class": ann.repeat_class,
                    "family": ann.family,
                    "strand": ann.strand,
                    "divergence": ann.divergence,
                    "length": ann.length,
                    "overlap_bp": overlap_bp,
                    "overlap_fraction": overlap_bp / max(1, ann.length),
                }
            )
        return results

    def __len__(self) -> int:
        return len(self.annotations)
