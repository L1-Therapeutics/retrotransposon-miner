"""L1 structural completeness and copy-number depth estimator.

Distinguishes full-length, potentially autonomous LINE-1 elements (~6,019 bp)
from 5'-truncated non-autonomous fragments produced during target-primed
reverse transcription (TPRT).

Active L1 elements are bicistronic: a 5' UTR RNA-polymerase-II internal
promoter, a ~1,042 bp ORF1, a ~3,966 bp ORF2, and a short 3' UTR.  More than
95% of genomic L1 insertions are 5'-truncated as a result of low
retrotranscriptase processivity and are therefore non-autonomous
(Beck et al., 2010, Annual Review of Genomics and Human Genetics 11:109-127).

Unmapped mate alignments are mapped onto the 6,019 bp L1HS consensus sequence
(1-based inclusive coordinates ``1..6,019``).  The 5' boundary of the element
is the 5'-most (minimum) covered position, the total functional span is
``3' boundary - 5' boundary + 1``, and the presence of each functional domain
is decided by the fraction of that domain covered by mate alignment intervals.
"""

from __future__ import annotations

from dataclasses import dataclass

L1HS_CONSENSUS_LENGTH = 6019

# Canonical L1HS domain boundaries (approximate; curated against a
# RepeatMasker / GenBank L1HS annotation before production use).
_L1HS_5PRIME_UTR = (1, 910)
_L1HS_ORF1 = (911, 1952)
_L1HS_ORF2 = (1953, 5918)
_L1HS_3PRIME_UTR = (5919, 6019)

# A domain counts as present when at least this fraction is covered by mates.
_DOMAIN_COVERAGE_THRESHOLD = 0.5

# An element is treated as full length when its 5' boundary lies within this
# many bp of the consensus 5' terminus.
_FULL_LENGTH_FIVE_PRIME_TOLERANCE = 100


@dataclass(frozen=True)
class L1CompletenessProfile:
    """Structural completeness / autonomy profile of an L1 element.

    Attributes:
        estimated_length: Total functional span (bp) covered by the mapped
            mates, ``3' boundary - 5' boundary + 1``.
        has_5prime_utr: At least ``_DOMAIN_COVERAGE_THRESHOLD`` of the 5' UTR
            is covered by mates.
        has_orf1: At least ``_DOMAIN_COVERAGE_THRESHOLD`` of ORF1 is covered.
        has_orf2: At least ``_DOMAIN_COVERAGE_THRESHOLD`` of ORF2 is covered.
        is_potentially_autonomous: The 5' boundary is within
            ``_FULL_LENGTH_FIVE_PRIME_TOLERANCE`` bp of the consensus start
            and the 5' UTR, ORF1 and ORF2 are all substantially covered.
    """

    estimated_length: int
    has_5prime_utr: bool
    has_orf1: bool
    has_orf2: bool
    is_potentially_autonomous: bool


def _merge_contiguous(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge inclusive intervals that overlap or touch, preserving order."""
    if not intervals:
        return []
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _covered_fraction(
    domain: tuple[int, int],
    intervals: list[tuple[int, int]],
) -> float:
    """Fraction of *domain* (1-based inclusive) covered by *intervals*."""
    domain_start, domain_end = domain
    clipped: list[tuple[int, int]] = []
    for start, end in intervals:
        if end < domain_start or start > domain_end:
            continue
        clipped.append((max(start, domain_start), min(end, domain_end)))
    covered = sum(
        end - start + 1 for start, end in _merge_contiguous(clipped)
    )
    return covered / (domain_end - domain_start + 1)


def profile_l1_completeness(
    unmapped_mate_alignments: list[tuple[int, int]],
) -> L1CompletenessProfile:
    """Estimate the structural completeness of an L1 element from mate alignments.

    Args:
        unmapped_mate_alignments: ``(start, end)`` intervals (1-based inclusive,
            either orientation accepted) of unmapped mate alignments projected
            onto the 6,019 bp L1HS consensus sequence.  Coordinates are clamped
            to ``1..6019``.

    Returns:
        An :class:`L1CompletenessProfile` describing the covered functional
        span and predicted autonomy.
    """
    intervals: list[tuple[int, int]] = []
    for start, end in unmapped_mate_alignments or []:
        lo, hi = sorted((int(start), int(end)))
        lo = max(1, min(lo, L1HS_CONSENSUS_LENGTH))
        hi = max(1, min(hi, L1HS_CONSENSUS_LENGTH))
        intervals.append((lo, hi))

    if not intervals:
        return L1CompletenessProfile(0, False, False, False, False)

    five_prime_boundary = min(start for start, _end in intervals)
    three_prime_boundary = max(end for _start, end in intervals)
    estimated_length = three_prime_boundary - five_prime_boundary + 1

    has_5prime_utr = (
        _covered_fraction(_L1HS_5PRIME_UTR, intervals) >= _DOMAIN_COVERAGE_THRESHOLD
    )
    has_orf1 = _covered_fraction(_L1HS_ORF1, intervals) >= _DOMAIN_COVERAGE_THRESHOLD
    has_orf2 = _covered_fraction(_L1HS_ORF2, intervals) >= _DOMAIN_COVERAGE_THRESHOLD
    is_potentially_autonomous = (
        five_prime_boundary <= _FULL_LENGTH_FIVE_PRIME_TOLERANCE
        and has_5prime_utr
        and has_orf1
        and has_orf2
    )
    return L1CompletenessProfile(
        estimated_length=estimated_length,
        has_5prime_utr=has_5prime_utr,
        has_orf1=has_orf1,
        has_orf2=has_orf2,
        is_potentially_autonomous=is_potentially_autonomous,
    )