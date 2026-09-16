"""Unit tests for EN-independent (DSB repair) MEI integration classification."""

from retro_miner.enb_independent import (
    EN_MOTIF_SCORE_THRESHOLD,
    classify_en_independent_event,
)


def test_zero_tsd_with_25bp_target_deletion_is_en_independent():
    call = classify_en_independent_event(tsd_length=0, en_score=0.0, span_deletion_bp=25)

    assert call.is_en_independent is True
    assert call.deletion_size == 25
    assert call.cleavage_motif_score == 0.0
    assert call.confidence > 0.0


def test_no_target_deletion_is_not_en_independent():
    call = classify_en_independent_event(tsd_length=0, en_score=0.0, span_deletion_bp=0)

    assert call.is_en_independent is False
    assert call.deletion_size == 0
    assert call.confidence == 0.0


def test_present_tsd_is_not_en_independent():
    call = classify_en_independent_event(tsd_length=7, en_score=0.0, span_deletion_bp=200)

    assert call.is_en_independent is False
    assert call.deletion_size == 200
    assert call.confidence == 0.0


def test_canonical_en_motif_is_not_en_independent():
    call = classify_en_independent_event(
        tsd_length=0,
        en_score=EN_MOTIF_SCORE_THRESHOLD,
        span_deletion_bp=500,
    )

    assert call.is_en_independent is False
    assert call.confidence == 0.0


def test_negative_inputs_are_clamped():
    call = classify_en_independent_event(
        tsd_length=-1, en_score=-0.5, span_deletion_bp=-10
    )

    assert call.is_en_independent is False
    assert call.deletion_size == 0
    assert call.cleavage_motif_score == 0.0


def test_en_score_above_threshold_is_neutralized_by_deletion_presence():
    call = classify_en_independent_event(tsd_length=0, en_score=0.9, span_deletion_bp=25)

    assert call.is_en_independent is False
    assert call.cleavage_motif_score == 0.9