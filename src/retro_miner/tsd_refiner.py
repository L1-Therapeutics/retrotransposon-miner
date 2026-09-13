"""Micro-homology aware Target Site Duplication (TSD) refiner.

Single-base precision resolution of the duplicated target site flanking
Target-Primed Reverse Transcription (TPRT) insertions (LINE-1, Alu, SVA).

The refiner ingests a genomic locus (``chrom``/``pos``), the supporting
soft-clipped reads, and a reference genome FASTA.  It returns a
:class:`TSDResult` describing the resolved TSD sequence, its length, a
confidence score, and the status of the insertion's 3' poly(A)/poly(T)
tail.

Algorithm
---------
1. Partition soft-clipped reads into the 5' (left-anchored) and 3'
   (right-anchored) arms using the clip-to-junction geometry relative to
   the insertion coordinate.
2. Detect the 3' poly(A)/poly(T) tail by Shannon entropy
   (``H = -sum p_i * log2(p_i)``) discrimination (``H < 1.2``) over the
   insert-facing homopolymer run of each clip.  ``A``-rich poly(A) tails
   and ``T``-rich poly(T) tails are both recognised.
3. Solve the two breakpoints from the median anchor coordinate of each
   arm: the 5' junction (left, 0-based exclusive) and the 3' junction
   (right, 0-based inclusive).
4. Find the duplicated target-site block by exact k-mer matching between
   the left clip's genomic portion and the right clip's reference flank
   (maximal flush match straddling the resolved junction).
5. Resolve micro-homology.  If the overlap is a tandem repeat of a
   shorter unit (a common motif such as ``AGCT`` shared between the donor
   end and the target site), collapse to the minimal period.  If a
   low-entropy poly(A)/poly(T) run abuts a breakpoint, the tail is placed
   *downstream* of the resolved TSD rather than inside it.
6. Reject pure poly(A)/poly(T) runs as artifact TSDs (mirrors
   ``_poly_at_artifact_tsd_mask`` in ``mei_support``) so a homopolymer
   flank is never reported as a genuine duplication.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

MIN_TSD_LEN = 4
MAX_TSD_LEN = 40
MIN_CLIP_LEN = 5
POLY_TAIL_MIN_RUN = 5
POLY_TAIL_ENTROPY_THRESHOLD = 1.2

CONFIDENCE_TSD_POLYA = 0.95
CONFIDENCE_TSD = 0.85
CONFIDENCE_MICROHOMOLOGY = 0.80
CONFIDENCE_BLUNT = 0.60
CONFIDENCE_DELETION = 0.60
CONFIDENCE_SHORT_OVERLAP = 0.30
CONFIDENCE_SINGLE_ARM = 0.20
CONFIDENCE_NONE = 0.0


@dataclass(frozen=True)
class TSDResult:
    """Resolved TSD and poly(A)/poly(T) tail status for one insertion locus."""

    chrom: str
    pos: int
    tsd_seq: str
    tsd_length: int
    tsd_confidence_score: float
    polyA_tail_detected: bool
    entropy: float
    method: str
    microhomology_seq: str = ""
    poly_tail_base: str = ""
    poly_tail_length: int = 0
    left_breakpoint: int = 0
    right_breakpoint: int = 0
    n_left_anchored: int = 0
    n_right_anchored: int = 0


def calculate_sequence_entropy(seq: str) -> float:
    """Shannon entropy ``H = -sum p_i * log2(p_i)`` of a nucleotide sequence."""
    if not seq:
        return 0.0
    seq_upper = seq.upper()
    length = len(seq_upper)
    counts: dict[str, int] = {}
    for char in seq_upper:
        counts[char] = counts.get(char, 0) + 1

    entropy = 0.0
    for count in counts.values():
        p = count / length
        entropy -= p * math.log2(p)
    return round(entropy, 4)


def is_poly_a_or_t(seq: str, max_entropy: float = POLY_TAIL_ENTROPY_THRESHOLD) -> bool:
    """True if *seq* is a low-entropy poly(A) or poly(T) fragment (``H < 1.2``)."""
    if len(seq) < 5:
        return False
    if calculate_sequence_entropy(seq) >= max_entropy:
        return False
    seq_upper = seq.upper()
    a_or_t_fraction = max(seq_upper.count("A"), seq_upper.count("T")) / len(seq_upper)
    return a_or_t_fraction >= 0.75


def find_longest_common_substring(s1: str, s2: str) -> str:
    """Longest exact matching substring shared by two sequences."""
    if not s1 or not s2:
        return ""
    m = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    longest_len = 0
    longest_end = 0

    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            if s1[i - 1].upper() == s2[j - 1].upper():
                m[i][j] = m[i - 1][j - 1] + 1
                if m[i][j] > longest_len:
                    longest_len = m[i][j]
                    longest_end = i

    return s1[longest_end - longest_len : longest_end].upper()


def _edge_homopolymer_run(seq: str, at_start: bool) -> tuple[str, int]:
    """Length and base of the homopolymer run at the start or end of *seq*."""
    if not seq:
        return "", 0
    ordered = seq if at_start else seq[::-1]
    base = ordered[0]
    run = 1
    for ch in ordered[1:]:
        if ch == base:
            run += 1
        else:
            break
    return base, run


def _minimal_period(seq: str) -> tuple[str, int]:
    """Return ``(minimal_unit, repeat_count)`` if *seq* is an exact tandem repeat.

    ``"AGCTAGCT"`` collapses to ``("AGCT", 2)``; a non-repetitive sequence
    returns ``(seq, 1)``.
    """
    n = len(seq)
    for unit_len in range(1, n // 2 + 1):
        if n % unit_len != 0:
            continue
        unit = seq[:unit_len]
        if unit * (n // unit_len) == seq:
            return unit, n // unit_len
    return seq, 1


def _max_flush_match(left: str, right: str) -> int:
    """Maximal ``k`` such that ``left[-k:] == right[:k]``.

    This is the exact junction-flush k-mer match between the left clip's
    genomic portion (*left*) and the right clip's reference flank (*right*).
    """
    k_max = min(len(left), len(right))
    for k in range(k_max, 0, -1):
        if left[-k:].upper() == right[:k].upper():
            return k
    return 0


def _trim_leading_poly_run(tsd_seq: str, poly_tail_base: str) -> str:
    """Trim a leading poly(A)/poly(T) run off *tsd_seq* (tail-adjacent micro-homology).

    When the target site shares a homopolymer run with the insertion's 3'
    poly(A)/poly(T) tail, the breakpoint is ambiguous.  Favor the breakpoint
    that places the tail downstream of the TSD by removing the run.
    """
    if not tsd_seq or poly_tail_base not in {"A", "T"}:
        return tsd_seq
    base, run = _edge_homopolymer_run(tsd_seq, at_start=True)
    if base == poly_tail_base and 2 <= run < len(tsd_seq):
        if len(tsd_seq) - run >= MIN_TSD_LEN and not is_poly_a_or_t(tsd_seq[run:]):
            return tsd_seq[run:]
    return tsd_seq


def _detect_poly_tail(clip_seqs: Sequence[str]) -> tuple[str, int]:
    """Detect the dominant low-entropy poly(A)/poly(T) run in *clip_seqs*.

    The insert-facing end of a soft clip carries the 3' tail, so the
    homopolymer run at either edge of each clip is inspected.  Returns
    ``(base, longest_run)`` or ``("", 0)`` when no tail is found.
    """
    best_base = ""
    best_run = 0
    for clip in clip_seqs:
        if len(clip) < POLY_TAIL_MIN_RUN:
            continue
        for at_start in (True, False):
            base, run = _edge_homopolymer_run(clip, at_start)
            if base in {"A", "T"} and run >= POLY_TAIL_MIN_RUN and run > best_run:
                if is_poly_a_or_t(base * run):
                    best_base = base
                    best_run = run
    return best_base, best_run


def _median_int(values: Sequence[int]) -> int:
    return int(round(statistics.median(values)))


def _resolve_overlap(
    overlap: str,
    poly_tail_base: str,
    min_tsd_len: int,
) -> tuple[str, str, str, bool]:
    """Resolve micro-homology in *overlap* -> ``(tsd, micro_seq, method, resolved)``."""
    unit, _reps = _minimal_period(overlap)
    collapsed = len(unit) >= min_tsd_len and len(unit) < len(overlap)
    micro_seq = ""
    if collapsed:
        micro_seq = overlap[len(unit):]
        overlap = unit

    trimmed = False
    trimmed_tsd = _trim_leading_poly_run(overlap, poly_tail_base)
    if trimmed_tsd != overlap:
        trimmed = True
        removed_run = overlap[: len(overlap) - len(trimmed_tsd)]
        if micro_seq:
            micro_seq = removed_run + micro_seq
        else:
            micro_seq = removed_run
        overlap = trimmed_tsd

    resolved = collapsed or trimmed
    if resolved:
        return overlap, micro_seq, "microhomology", True
    return overlap, micro_seq, "", False


def refine_tsd_boundaries(
    chrom: str,
    pos: int,
    soft_clips: Sequence[pysam.AlignedSegment],
    ref_genome_fasta: str | Path | pysam.FastaFile,
    min_tsd_len: int = MIN_TSD_LEN,
    max_tsd_len: int = MAX_TSD_LEN,
    clip_min_len: int = MIN_CLIP_LEN,
    flank_window_bp: int = 80,
) -> TSDResult:
    """Refine the TSD boundaries for one insertion locus.

    Args:
        chrom: Reference chromosome/contig name.
        pos: 1-based insertion breakpoint estimate (left/5' side).
        soft_clips: Supporting reads carrying soft-clipped insert sequence.
        ref_genome_fasta: Path to a reference FASTA (with .fai) or an open
            :class:`pysam.FastaFile` (which will not be closed here).
        min_tsd_len: Smallest biological TSD length considered real.
        max_tsd_len: Largest TSD length considered a single block.
        clip_min_len: Minimum soft-clip length used as evidence.
        flank_window_bp: Reference window flanking each junction used for the
            exact k-mer resolution.

    Returns:
        :class:`TSDResult` with resolved TSD, confidence, and tail status.
    """
    if min_tsd_len < 1:
        min_tsd_len = 1
    if max_tsd_len < min_tsd_len:
        max_tsd_len = min_tsd_len
    if clip_min_len < 1:
        clip_min_len = 1
    if flank_window_bp < max_tsd_len:
        flank_window_bp = max_tsd_len

    pos0 = int(pos) - 1
    slack = max(max_tsd_len, 50) + 10

    left_anchors: list[int] = []
    right_anchors: list[int] = []
    clip_seqs: list[str] = []

    for read in soft_clips:
        if read is None or read.is_unmapped or read.query_sequence is None:
            continue
        if read.cigartuples is None:
            continue
        qs = read.query_sequence
        first_op, last_op = read.cigartuples[0], read.cigartuples[-1]
        leading_clip = first_op[1] if first_op[0] == 4 else 0
        trailing_clip = last_op[1] if last_op[0] == 4 else 0

        if leading_clip >= clip_min_len and read.reference_start is not None:
            if read.reference_start >= pos0 - slack:
                right_anchors.append(int(read.reference_start))
                clip_seqs.append(qs[:leading_clip].upper())

        if trailing_clip >= clip_min_len and read.reference_end is not None:
            if read.reference_end <= pos0 + slack:
                left_anchors.append(int(read.reference_end))
                clip_seqs.append(qs[len(qs) - trailing_clip :].upper())

    if not left_anchors and not right_anchors:
        return TSDResult(
            chrom=chrom,
            pos=int(pos),
            tsd_seq="",
            tsd_length=0,
            tsd_confidence_score=CONFIDENCE_NONE,
            polyA_tail_detected=False,
            entropy=0.0,
            method="no_soft_clips",
        )

    poly_tail_base, poly_tail_length = _detect_poly_tail(clip_seqs)
    polyA_tail_detected = poly_tail_base != ""
    avg_entropy = round(
        sum(calculate_sequence_entropy(c) for c in clip_seqs) / len(clip_seqs),
        4,
    )

    close_reference = isinstance(ref_genome_fasta, str) or isinstance(
        ref_genome_fasta, Path
    )
    try:
        if close_reference:
            ref = pysam.FastaFile(str(ref_genome_fasta))
        else:
            ref = ref_genome_fasta
        try:
            fetched = ref.fetch(str(chrom))
            if isinstance(fetched, bytes):
                fetched = fetched.decode("ascii")
            chrom_seq = fetched.upper()
        except (KeyError, ValueError):
            chrom_seq = ""
    finally:
        if close_reference:
            ref.close()

    if not chrom_seq:
        return TSDResult(
            chrom=chrom,
            pos=int(pos),
            tsd_seq="",
            tsd_length=0,
            tsd_confidence_score=CONFIDENCE_NONE,
            polyA_tail_detected=polyA_tail_detected,
            entropy=avg_entropy,
            method="no_reference_contig",
            poly_tail_base=poly_tail_base,
            poly_tail_length=poly_tail_length,
        )

    left_bp = _median_int(left_anchors) if left_anchors else None
    right_bp = _median_int(right_anchors) if right_anchors else None

    tsd_seq = ""
    tsd_length = 0
    confidence = CONFIDENCE_NONE
    method = "none"
    microhomology_seq = ""

    if left_bp is not None and right_bp is not None:
        left_bp_val = left_bp
        right_bp_val = right_bp
        overlap = left_bp_val - right_bp_val
        if overlap < 0:
            method = "deletion"
            confidence = CONFIDENCE_DELETION
        elif overlap == 0:
            method = "blunt"
            confidence = CONFIDENCE_BLUNT
        elif overlap < min_tsd_len:
            candidate = chrom_seq[right_bp_val:left_bp_val]
            method = "short_overlap"
            confidence = CONFIDENCE_SHORT_OVERLAP
            tsd_seq = candidate
            tsd_length = len(candidate)
        else:
            left_genomic = chrom_seq[
                max(0, left_bp_val - flank_window_bp):left_bp_val
            ]
            right_flank = chrom_seq[
                right_bp_val : right_bp_val + flank_window_bp
            ]
            flush_len = _max_flush_match(left_genomic, right_flank)
            if min_tsd_len <= flush_len <= max_tsd_len and flush_len >= overlap:
                solved_end = right_bp_val + flush_len
            else:
                solved_end = left_bp_val
            candidate = chrom_seq[right_bp_val:solved_end]
            if len(candidate) > max_tsd_len:
                candidate = chrom_seq[right_bp_val:left_bp_val]

            tsd_seq, microhomology_seq, micro_method, resolved = _resolve_overlap(
                candidate, poly_tail_base, min_tsd_len
            )
            if resolved:
                method = micro_method
                confidence = CONFIDENCE_MICROHOMOLOGY
            else:
                if is_poly_a_or_t(tsd_seq):
                    tsd_seq = ""
                    method = "polyA_artifact"
                    confidence = CONFIDENCE_NONE
                elif len(tsd_seq) < min_tsd_len:
                    method = "short_overlap"
                    confidence = CONFIDENCE_SHORT_OVERLAP
                else:
                    method = "tsd_flush_match"
                    confidence = (
                        CONFIDENCE_TSD_POLYA if polyA_tail_detected else CONFIDENCE_TSD
                    )
            tsd_length = len(tsd_seq)
    elif left_bp is not None:
        method = "left_arm_only"
        confidence = CONFIDENCE_SINGLE_ARM
        left_bp_val = left_bp
        right_bp_val = 0
    elif right_bp is not None:
        method = "right_arm_only"
        confidence = CONFIDENCE_SINGLE_ARM
        left_bp_val = 0
        right_bp_val = right_bp
    else:
        left_bp_val = 0
        right_bp_val = 0

    return TSDResult(
        chrom=chrom,
        pos=int(pos),
        tsd_seq=tsd_seq,
        tsd_length=tsd_length,
        tsd_confidence_score=confidence,
        polyA_tail_detected=polyA_tail_detected,
        entropy=avg_entropy,
        method=method,
        microhomology_seq=microhomology_seq,
        poly_tail_base=poly_tail_base,
        poly_tail_length=poly_tail_length,
        left_breakpoint=left_bp_val,
        right_breakpoint=right_bp_val,
        n_left_anchored=len(left_anchors),
        n_right_anchored=len(right_anchors),
    )