"""Tests for L1 structural completeness / autonomy profiling.

Synthetic mate-alignment scenarios:
1. Reads spanning positions 1-6019 of the L1HS consensus -> full length,
   potentially autonomous.
2. 5'-truncated reads spanning positions 4500-6019 -> non-autonomous (5' UTR
   and ORF1 absent / uncovered).
3. Empty input -> zero-length, non-autonomous.
4. Partial 5' UTR and partial ORF2 coverage -> not autonomous.
5. Interleaved full-length + truncated reads -> 5' boundary at position 1.
6. Reversed tuples and out-of-range coordinates are normalised.
"""

from retro_miner.completeness_profiler import (
    L1CompletenessProfile,
    profile_l1_completeness,
)

FULL_LENGTH_INTERVALS = [(1, 2500), (2501, 4700), (4701, 6019)]
ORF2_ONLY_INTERVALS = [(4500, 5210), (5211, 6019)]


def test_full_length_reads_are_autonomous():
    profile = profile_l1_completeness(FULL_LENGTH_INTERVALS)
    assert isinstance(profile, L1CompletenessProfile)
    assert profile.estimated_length == 6019
    assert profile.has_5prime_utr is True
    assert profile.has_orf1 is True
    assert profile.has_orf2 is True
    assert profile.is_potentially_autonomous is True


def test_single_full_length_span_is_autonomous():
    profile = profile_l1_completeness([(1, 6019)])
    assert profile.estimated_length == 6019
    assert profile.is_potentially_autonomous is True


def test_five_prime_truncated_reads_not_autonomous():
    profile = profile_l1_completeness(ORF2_ONLY_INTERVALS)
    assert profile.estimated_length == 1520
    assert profile.has_5prime_utr is False
    assert profile.has_orf1 is False
    assert profile.has_orf2 is False
    assert profile.is_potentially_autonomous is False


def test_partial_domain_coverage_not_autonomous():
    profile = profile_l1_completeness([(750, 3000)])
    assert profile.has_5prime_utr is False
    assert profile.has_orf2 is False
    assert profile.is_potentially_autonomous is False


def test_orf2_present_without_orfs_5prime_is_not_autonomous():
    profile = profile_l1_completeness([(3000, 6019)])
    assert profile.estimated_length == 3020
    assert profile.has_5prime_utr is False
    assert profile.has_orf1 is False
    assert profile.has_orf2 is True
    assert profile.is_potentially_autonomous is False


def test_full_length_plus_truncated_reads_still_autonomous():
    profile = profile_l1_completeness(
        [(1, 6019)] + ORF2_ONLY_INTERVALS
    )
    assert profile.estimated_length == 6019
    assert profile.is_potentially_autonomous is True


def test_empty_input_is_not_autonomous():
    profile = profile_l1_completeness([])
    assert profile.estimated_length == 0
    assert profile.has_5prime_utr is False
    assert profile.has_orf1 is False
    assert profile.has_orf2 is False
    assert profile.is_potentially_autonomous is False


def test_reversed_and_out_of_range_coordinates_are_normalised():
    profile = profile_l1_completeness([(6019, 1), (0, 6100)])
    assert profile.estimated_length == 6019
    assert profile.is_potentially_autonomous is True