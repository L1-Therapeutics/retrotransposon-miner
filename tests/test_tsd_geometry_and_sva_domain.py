"""Unit tests for TSD (target-site duplication) geometry resolution and
SVA multi-domain classification scoring.

All functions under test are module-level helpers in ``mei_support`` that
require no BAM files, reference FASTA, or external tools.

TSD geometry coverage
---------------------
* Canonical TSD (2-30 bp)              — classic happy path
* 1-bp micro-TSD                       — edge of canonical range
* 0-bp blunt insertion (right < left)  — common in 5'-truncated LINE-1
* Target-site deletion (tsd_len < 0)   — reference bases removed at insertion
* ±2 bp rescue for near-miss pairs     — including blunt + deletion rescues
* Higher support wins tie-breaks
* Empty candidates                     — must return sentinel (0, 0, 0, '')
* tsd_len > 30                         — outside accepted range, rejected

SVA domain scoring coverage
----------------------------
* VNTR hexamer score on pure CCCTCT sequence
* SINE-R k-mer score on HERV-K K-box motif-bearing sequence
* Combined score is weighted blend (VNTR-dominant)
* Alu-like sequence scores near-zero on both SINE-R and combined score
* Short sequences below minimum length return 0.0
"""

from __future__ import annotations

import pytest

from retro_miner.mei_support import (
    _resolve_tsd_pair_from_candidates,
    _sva_combined_domain_score,
    _sva_siner_domain_score,
    _sva_vntr_like_score,
)


# ===========================================================================
# Helpers
# ===========================================================================

def _make_candidate(
    left: int,
    right: int,
    support: int = 5,
    source: str = "tsd_disease",
) -> tuple[int, int, int, str]:
    """Convenience factory matching (left, right, support, source) convention."""
    return (left, right, support, source)


# ===========================================================================
# _resolve_tsd_pair_from_candidates — empty / sentinel
# ===========================================================================


def test_empty_candidates_returns_sentinel() -> None:
    left, right, length, source = _resolve_tsd_pair_from_candidates([])
    assert left == 0 and right == 0 and length == 0 and source == ""


# ===========================================================================
# Canonical TSD (2-30 bp) — the happy path
# ===========================================================================


def test_canonical_tsd_5bp() -> None:
    # right - left + 1 = 5 → classic canonical TSD
    left, right, length, source = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1000, 1004)]
    )
    assert left == 1000
    assert right == 1004
    assert length == 5
    assert source == "tsd_disease"


def test_canonical_tsd_exact_maximum_30bp() -> None:
    # TSD of exactly 30 bp must be accepted.
    left, right, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(500, 529)]
    )
    assert length == 30


def test_canonical_tsd_31bp_rescued_to_30() -> None:
    # 31 bp fails the strict pass, but the rescue pass applies a -1 shift on
    # the right coordinate (dr=-1, shift_penalty=1) to yield tsd_len=30.
    # This is the intended behaviour: small measurement errors are recoverable.
    left, right, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(500, 530)]  # strict: 31 bp; rescue: 30 bp with dr=-1
    )
    assert length == 30


def test_1bp_micro_tsd_accepted() -> None:
    # tsd_len = 1 (right == left) is within the new extended range
    left, right, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(200, 200)]
    )
    assert length == 1


# ===========================================================================
# 0-bp blunt insertion  (right == left - 1  →  tsd_len == 0)
# ===========================================================================


def test_blunt_insertion_0bp_accepted() -> None:
    """0-bp blunt insertion: right == left - 1, tsd_len == 0.

    This geometry is biologically real for 5'-truncated LINE-1 events that
    insert without duplicating a target site.  The old guard ``2 <= tsd_len``
    silently dropped these; the new guard ``0 <= tsd_len`` accepts them.
    """
    # left=1001, right=1000 → tsd_len = 1000 - 1001 + 1 = 0
    left, right, length, source = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1001, 1000)]
    )
    assert length == 0
    assert left == 1001
    assert right == 1000
    assert source == "tsd_disease"


def test_blunt_insertion_preferred_over_deletion_when_both_present() -> None:
    """When both a 0-bp blunt and a target-site deletion are candidates,
    the blunt insertion (longer non-negative TSD) is preferred."""
    blunt = _make_candidate(1001, 1000, support=3)   # tsd_len = 0
    deletion = _make_candidate(1005, 1000, support=3)  # tsd_len = -4
    _, _, length, _ = _resolve_tsd_pair_from_candidates([blunt, deletion])
    assert length == 0  # blunt wins over deletion at equal support


def test_canonical_preferred_over_blunt_at_equal_support() -> None:
    """Canonical TSDs outrank 0-bp blunt insertions in the sort key."""
    blunt = _make_candidate(1001, 1000, support=5)   # tsd_len = 0
    canonical = _make_candidate(1001, 1005, support=5)  # tsd_len = 5
    _, _, length, _ = _resolve_tsd_pair_from_candidates([blunt, canonical])
    assert length == 5


def test_blunt_preferred_over_canonical_when_support_higher() -> None:
    """Higher read support overrides the canonical-TSD length preference."""
    blunt = _make_candidate(1001, 1000, support=20)   # tsd_len = 0
    canonical = _make_candidate(1001, 1005, support=3)  # tsd_len = 5
    _, _, length, source = _resolve_tsd_pair_from_candidates([blunt, canonical])
    assert length == 0
    assert source == "tsd_disease"


# ===========================================================================
# Target-site deletion  (tsd_len < 0)
# ===========================================================================


def test_target_site_deletion_minus3_accepted() -> None:
    """A -3 bp target-site deletion is within _TSD_MAX_DELETION_BP=10."""
    # left=1005, right=1001 → tsd_len = 1001 - 1005 + 1 = -3
    left, right, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1005, 1001)]
    )
    assert length == -3
    assert left == 1005
    assert right == 1001


def test_target_site_deletion_minus10_boundary_accepted() -> None:
    """Exactly -_TSD_MAX_DELETION_BP=-10 should be accepted (boundary)."""
    # tsd_len = right - left + 1 = -10  →  right = left - 11
    left, right, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1011, 1000)]
    )
    assert length == -10


def test_target_site_deletion_11bp_rescued_to_boundary() -> None:
    """A deletion of -11 bp fails the strict pass but is recoverable by the
    rescue pass: a shift of dl=-1 (or dr=+1) yields tsd_len=-10, which sits
    exactly at _TSD_MAX_DELETION_BP=10 and is accepted.  The rescue pass
    intentionally applies ±2 bp coordinate shifts to tolerate small
    measurement errors in breakpoint mode estimation."""
    # left=1012, right=1000: strict tsd_len = 1000 - 1012 + 1 = -11 (rejected)
    # rescue with dl=-1: ll=1011, rr=1000, tsd_len = -10 (accepted at boundary)
    _, _, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1012, 1000)]
    )
    assert length == -10


def test_target_site_deletion_far_beyond_max_returns_sentinel() -> None:
    """A deletion so large that even the rescue ±2 shift cannot bring it within
    [-_TSD_MAX_DELETION_BP, 30] returns the (0, 0, 0, '') sentinel.

    With _TSD_MAX_DELETION_BP=10 and strict rejection at tsd_len < -10, the
    rescue pass can shift by at most +4 (dl=-2, dr=+2).  A -15 bp deletion
    becomes at best -11 under rescue, still outside the window.
    """
    # left=1016, right=1000: strict tsd_len = -15; best rescue: -11 (still rejected)
    _, _, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1016, 1000)]
    )
    assert length == 0


# ===========================================================================
# ±2 bp rescue pass
# ===========================================================================


def test_rescue_near_miss_canonical_recovered() -> None:
    """A pair that misses canonical range by 1 bp is rescued with a ±1 shift."""
    # tsd_len = 31 bp (1 over limit), shift dr=-1 gives tsd_len=30 → accepted
    _, _, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(500, 530)]  # strict: 31 bp, rescue: 30 bp
    )
    assert length == 30


def test_rescue_blunt_insertion_near_miss() -> None:
    """A near-blunt pair (tsd_len = -1, just below 0) can be rescued to 0."""
    # left=1002, right=1000 → tsd_len = -1; shift dl=-1 gives ll=1001, tsd_len=0
    _, _, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(1002, 1000)]
    )
    # Rescue can produce 0-bp (shift) or keep -1 (within -10 bound);
    # the strict pass accepts -1, so result must be -1 (strict preferred).
    assert length == -1


def test_rescue_returns_sentinel_when_no_shift_helps() -> None:
    """When even ±2 bp shifts cannot produce a valid TSD, return sentinel."""
    # tsd_len = 35; no ±2 shift can bring it within [−10, 30]
    _, _, length, _ = _resolve_tsd_pair_from_candidates(
        [_make_candidate(100, 134)]  # tsd_len = 35
    )
    assert length == 0


# ===========================================================================
# Source (disease vs control) and multi-candidate tie-breaks
# ===========================================================================


def test_rescue_pass_disease_source_preferred_over_control_at_equal_support() -> None:
    """In the rescue pass, 'tsd_disease' has sample_priority=0 vs
    'tsd_control' sample_priority=1.  When shift_penalty and support are
    equal, disease source wins.  The strict pass does not encode source
    priority; this test uses candidates that require rescue (tsd_len=31)."""
    # Both need rescue (31 bp > 30 bp limit); equal support; disease should win
    disease = _make_candidate(500, 530, support=5, source="tsd_disease")  # tsd_len=31
    control = _make_candidate(500, 530, support=5, source="tsd_control")  # tsd_len=31
    # With dr=-1 both yield tsd_len=30; tie-break: sample_priority disease=0 < control=1
    _, _, length, source = _resolve_tsd_pair_from_candidates([control, disease])
    assert length == 30
    assert source == "tsd_disease"


def test_higher_support_wins_regardless_of_source() -> None:
    """Higher read support always outranks disease-priority."""
    low_disease = _make_candidate(1000, 1005, support=2, source="tsd_disease")
    high_control = _make_candidate(1000, 1005, support=10, source="tsd_control")
    _, _, _, source = _resolve_tsd_pair_from_candidates([low_disease, high_control])
    assert source == "tsd_control"


# ===========================================================================
# _sva_vntr_like_score
# ===========================================================================


def test_vntr_score_pure_ccctct_repeats_high() -> None:
    """A string of CCCTCT hexamers must score >= 0.35 (threshold for rescue)."""
    seq = "CCCTCT" * 12  # 72 bp of pure VNTR
    score = _sva_vntr_like_score(seq)
    assert score >= 0.35


def test_vntr_score_random_dna_low() -> None:
    """Random non-VNTR sequence scores near zero."""
    seq = "ACGTACGTGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"  # 48 bp
    score = _sva_vntr_like_score(seq)
    assert score < 0.20


def test_vntr_score_short_seq_returns_zero() -> None:
    """Sequences below _VNTR_MIN_SEQ_LEN (40 bp) return 0.0."""
    assert _sva_vntr_like_score("CCCTCT" * 5) == pytest.approx(0.0)  # 30 bp


def test_vntr_score_empty_returns_zero() -> None:
    assert _sva_vntr_like_score("") == pytest.approx(0.0)


# ===========================================================================
# _sva_siner_domain_score
# ===========================================================================


def test_siner_score_kbox_motif_positive() -> None:
    """A sequence containing the TTGCAAACCAA K-box motif must score > 0."""
    # Embed K-box inside random context to avoid length issues.
    seq = "ACGTACGTACGTACGTTTGCAAACCAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
    score = _sva_siner_domain_score(seq)
    assert score > 0.0


def test_siner_score_multiple_motifs_higher_than_one() -> None:
    """More SINE-R k-mer hits increase the score."""
    one_hit = "A" * 30 + "TTGCAAACCAA" + "G" * 30
    two_hits = "A" * 20 + "TTGCAAACCAA" + "C" * 10 + "AACGCAAACCAA" + "G" * 10
    assert _sva_siner_domain_score(two_hits) >= _sva_siner_domain_score(one_hit)


def test_siner_score_alu_consensus_near_zero() -> None:
    """Alu AluY consensus lacks SINE-R k-mers; score must be 0.0."""
    # Representative AluY body (no HERV-K K-box motifs).
    alu_fragment = (
        "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAG"
        "GTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTA"
    )
    assert _sva_siner_domain_score(alu_fragment) == pytest.approx(0.0)


def test_siner_score_short_seq_returns_zero() -> None:
    """Sequences below _SINER_MIN_SEQ_LEN (30 bp) return 0.0."""
    assert _sva_siner_domain_score("TTGCAAACCAA") == pytest.approx(0.0)  # 11 bp


def test_siner_score_empty_returns_zero() -> None:
    assert _sva_siner_domain_score("") == pytest.approx(0.0)


# ===========================================================================
# _sva_combined_domain_score
# ===========================================================================


def test_combined_score_is_weighted_blend() -> None:
    """Combined score lies between the two component scores (weighted blend)."""
    seq = "CCCTCT" * 12  # pure VNTR, no SINE-R
    v = _sva_vntr_like_score(seq)
    r = _sva_siner_domain_score(seq, min_seq_len=len(seq))
    combined = _sva_combined_domain_score(seq)
    # With default weights (0.65 VNTR + 0.35 SINE-R), combined = 0.65*v + 0.35*r
    expected = 0.65 * v + 0.35 * r
    assert combined == pytest.approx(expected, abs=1e-9)


def test_combined_score_higher_for_sva_than_alu() -> None:
    """SVA-like sequence (VNTR + SINE-R) outscores Alu-like sequence."""
    sva_seq = "CCCTCT" * 8 + "TTGCAAACCAA" * 2 + "CCCTCT" * 3  # VNTR + K-box
    alu_seq = (
        "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAG"
        "GTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTA"
    )
    assert _sva_combined_domain_score(sva_seq) > _sva_combined_domain_score(alu_seq)


def test_combined_score_alu_sequence_below_rescue_threshold() -> None:
    """Alu consensus must fall below the 0.35 VNTR-rescue threshold."""
    alu_seq = (
        "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAG"
        "GTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTA"
    )
    assert _sva_combined_domain_score(alu_seq) < 0.35


def test_combined_score_empty_returns_zero() -> None:
    assert _sva_combined_domain_score("") == pytest.approx(0.0)


def test_combined_custom_weights_respected() -> None:
    """Passing vntr_weight=1.0, siner_weight=0.0 degrades to pure VNTR score."""
    seq = "CCCTCT" * 12
    pure_vntr = _sva_vntr_like_score(seq)
    weighted = _sva_combined_domain_score(seq, vntr_weight=1.0, siner_weight=0.0)
    assert weighted == pytest.approx(pure_vntr, abs=1e-9)
