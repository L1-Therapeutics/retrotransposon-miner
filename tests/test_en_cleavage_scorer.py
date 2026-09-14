"""Unit tests for L1 EN cleavage site preference scoring."""

from retro_miner.en_cleavage_scorer import (
    CANONICAL_EN_SITE_THRESHOLD,
    CLASS_CANONICAL_EN,
    CLASS_EN_INDEPENDENT,
    EN_CONSENSUS,
    score_en_cleavage_site,
)


def test_canonical_tittt_aa_motif():
    context = "GTTTTAAACC"
    result = score_en_cleavage_site(context)

    assert result.is_canonical_en_site is True
    assert result.raw_score >= 0.90
    assert result.motif == EN_CONSENSUS
    assert result.cleavage_class == CLASS_CANONICAL_EN


def test_noncanonical_gcgg_cc_motif():
    context = "AGCGGCCCTA"
    result = score_en_cleavage_site(context)

    assert result.is_canonical_en_site is False
    assert result.raw_score < 0.30
    assert result.cleavage_class != CLASS_CANONICAL_EN


def test_canonical_threshold_is_exposed():
    assert CANONICAL_EN_SITE_THRESHOLD == 0.70


def test_empty_context_is_en_independent():
    result = score_en_cleavage_site("")

    assert result.is_canonical_en_site is False
    assert result.raw_score == 0.0
    assert result.cleavage_class == CLASS_EN_INDEPENDENT


def test_lowercase_context_is_upper_cased():
    lower = score_en_cleavage_site("gttttaaacc")
    upper = score_en_cleavage_site("GTTTTAAACC")

    assert lower.raw_score == upper.raw_score
    assert lower.motif == EN_CONSENSUS


def test_degenerate_ttttag_scores_above_threshold():
    context = "GTTTTAGCCC"
    result = score_en_cleavage_site(context)

    assert result.motif == "TTTTAG"
    assert result.raw_score >= CANONICAL_EN_SITE_THRESHOLD
    assert result.is_canonical_en_site is True