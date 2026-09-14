"""Scientific unit tests for the L1 EN cleavage-motif PWM scorer (Jurka 1997).

Covers the canonical ``5'-TTTT/AA-3'`` consensus, a degenerate nick-adjacent
variant, and non-canonical genomic context, plus window extraction, padding and
dataclass semantics.
"""

import math

import pytest

from retro_miner.cleavage_motif import (
    CANONICAL_TPRT_THRESHOLD,
    CLASS_CANONICAL_TPRT,
    CLASS_DEGENERATE_TPRT,
    CLASS_EN_INDEPENDENT,
    DEGENERATE_TPRT_THRESHOLD,
    CleavageMotifResult,
    score_cleavage_motif,
)

CANONICAL_MOTIF = "TTTTAAAA"
DEGENERATE_MOTIF = "TTTAAAAC"
EN_INDEPENDENT_MOTIF = "GCAGCTAG"


def test_canonical_motif_is_canonical_tprt():
    res = score_cleavage_motif(CANONICAL_MOTIF)
    assert res.is_canonical_tprt is True
    assert res.cleavage_class == CLASS_CANONICAL_TPRT
    assert res.pwm_score >= 4.0
    assert res.pwm_score >= CANONICAL_TPRT_THRESHOLD
    assert res.breakpoint_motif == CANONICAL_MOTIF


def test_degenerate_motif_is_degenerate_tprt():
    res = score_cleavage_motif(DEGENERATE_MOTIF)
    assert res.is_canonical_tprt is False
    assert res.cleavage_class == CLASS_DEGENERATE_TPRT
    assert DEGENERATE_TPRT_THRESHOLD <= res.pwm_score < CANONICAL_TPRT_THRESHOLD
    assert res.breakpoint_motif == DEGENERATE_MOTIF


def test_non_canonical_context_is_en_independent():
    res = score_cleavage_motif(EN_INDEPENDENT_MOTIF)
    assert res.is_canonical_tprt is False
    assert res.cleavage_class == CLASS_EN_INDEPENDENT
    assert res.pwm_score < DEGENERATE_TPRT_THRESHOLD
    assert res.breakpoint_motif == EN_INDEPENDENT_MOTIF


def test_log_odds_score_relative_to_uniform_background():
    res = score_cleavage_motif(CANONICAL_MOTIF)

    def log_odds(frequency: float) -> float:
        return math.log2(frequency / 0.25)

    expected = (
        3 * log_odds(0.70)  # positions 1-3: P(T)
        + log_odds(0.84)  # position 4: nick-adjacent P(T)
        + 2 * log_odds(0.78)  # positions 5-6: P(A)
        + log_odds(0.34)  # position 7: weak P(A)
        + log_odds(0.32)  # position 8: weak P(A)
    )
    assert res.pwm_score == pytest.approx(expected, abs=1e-3)
    assert res.pwm_score > 0.0


def test_window_extraction_from_flanked_reference():
    flank_upstream = "AACCGG"
    flank_downstream = "GGAATT"
    ref = flank_upstream + "TTTTAAAA" + flank_downstream
    res = score_cleavage_motif(ref, breakpoint_offset=len(flank_upstream) + 4)
    assert res.breakpoint_motif == CANONICAL_MOTIF
    assert res.is_canonical_tprt is True


def test_breakpoint_offset_controls_window():
    ref = "AACC" + "TTTTAAAA" + "GG"
    res = score_cleavage_motif(ref, breakpoint_offset=4)
    assert res.breakpoint_motif == "AACCTTTT"
    res_shifted = score_cleavage_motif(ref, breakpoint_offset=8)
    assert res_shifted.breakpoint_motif == CANONICAL_MOTIF
    assert res_shifted.is_canonical_tprt is True


def test_short_reference_is_padded_with_neutral_n():
    res = score_cleavage_motif(CANONICAL_MOTIF)
    assert "N" not in res.breakpoint_motif
    short = score_cleavage_motif("TTTT")
    assert len(short.breakpoint_motif) == 8
    assert short.breakpoint_motif == "TTTTNNNN"
    expected = 3 * math.log2(0.70 / 0.25) + math.log2(0.84 / 0.25)
    assert short.pwm_score == pytest.approx(expected, abs=1e-3)


def test_lowercase_reference_is_normalized():
    res = score_cleavage_motif("ttttaaaa")
    assert res.breakpoint_motif == CANONICAL_MOTIF
    assert res.is_canonical_tprt is True
    assert res.pwm_score == score_cleavage_motif(CANONICAL_MOTIF).pwm_score


def test_iupac_ambiguous_bases_score_neutral():
    res = score_cleavage_motif("TTTRAAA?", breakpoint_offset=4)
    assert res.breakpoint_motif == "TTTRAAA?"
    assert res.pwm_score == score_cleavage_motif("TTTNAAAN").pwm_score


def test_boundary_threshold_zero_for_empty_reference():
    res = score_cleavage_motif("", breakpoint_offset=0)
    assert res.breakpoint_motif == "NNNNNNNN"
    assert res.cleavage_class == CLASS_EN_INDEPENDENT
    assert res.pwm_score == 0.0


def test_result_is_frozen_dataclass():
    res = score_cleavage_motif(CANONICAL_MOTIF)
    assert isinstance(res, CleavageMotifResult)
    with pytest.raises(AttributeError):
        res.breakpoint_motif = "CHANGED"


def test_canonical_ranked_above_degenerate():
    canonical = score_cleavage_motif(CANONICAL_MOTIF).pwm_score
    degenerate = score_cleavage_motif(DEGENERATE_MOTIF).pwm_score
    assert canonical > degenerate