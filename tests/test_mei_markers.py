"""Tests for literature-anchored structural markers (mei_markers.py).

The constants encode numbers read from the primary papers; the tests pin
those numbers so a future edit cannot silently drift them.
"""

from __future__ import annotations

import pytest

from retro_miner.mei_markers import (
    EN_MOTIF_WEIGHTS,
    en_motif_score,
    en_motif_site_rank,
    en_motif_site_weight,
    ALU_HOTSPOT_POSITIONS,
    EN_MOTIF_7MERS_97PCT,
    EN_MOTIF_7MERS_TOTAL,
    EN_NICK_OFFSET,
    INTERNAL_PRIMING_MIN_TSD,
    MH_MAX_LEN,
    MH_TRUNCATED_L1_COUNTS,
    MH_TRUNCATED_L1_RATE,
    RouteEvidence,
    TSD_MAX_BP,
    TSD_MIN_BP,
    TRIAD_FRACTION,
    classify_route,
    en_motif_opportunity,
    flag_complex_event,
    internal_priming_signature,
    levy_hotspot_likelihood,
    microhomology_at_5p_junction,
    revcomp,
    scan_en_motif_strand,
    tsd_len_in_range,
    ComplexEventKind,
    flasch_weights_provenance,
)


# ---------------------------------------------------------------------------
# Literature constants
# ---------------------------------------------------------------------------

def test_flasch_motif_constants():
    assert EN_MOTIF_7MERS_TOTAL == 12_288
    assert EN_MOTIF_7MERS_97PCT == 743
    assert pytest.approx(743 / 12_288, abs=0.005) == 0.06


def test_en_nick_offset():
    # Empirical consensus TTTTT/AA: five T's 5' of the nick, two A's 3'.
    assert EN_NICK_OFFSET == 5


def test_wagstaff_tsd_range_and_triad():
    assert (TSD_MIN_BP, TSD_MAX_BP) == (5, 27)
    assert TRIAD_FRACTION == pytest.approx(0.96)


def test_zingler_microhomology_rates():
    assert MH_TRUNCATED_L1_COUNTS == (1_503, 2_427)
    assert MH_TRUNCATED_L1_RATE == pytest.approx(1_503 / 2_427, abs=0.005)
    assert MH_MAX_LEN == 12


def test_internal_priming_min_tsd():
    assert INTERNAL_PRIMING_MIN_TSD == 6


def test_tsd_len_in_range():
    assert tsd_len_in_range(5)
    assert tsd_len_in_range(27)
    assert not tsd_len_in_range(4)
    assert not tsd_len_in_range(28)
    assert not tsd_len_in_range(None)


# ---------------------------------------------------------------------------
# Motif scanning
# ---------------------------------------------------------------------------

def test_flasch_table_paper_constraints():
    """The embedded tables must satisfy the paper's own constraints.

    Position 6 can never be T (Fig. 2E / STAR method V constraint), all
    sites are 7-nt ACGT, weights are descending, and the canonical site's
    rank/weight match the paper text exactly.
    """
    corrected = EN_MOTIF_WEIGHTS["corrected"]
    uncorrected = EN_MOTIF_WEIGHTS["uncorrected"]
    for table in (corrected, uncorrected):
        for site, weight in table.items():
            assert len(site) == 7 and set(site) <= set("ACGT")
            assert site[5] != "T", site  # V constraint at position 6
            assert 0 < weight <= 1.0
        weights = list(table.values())
        assert weights == sorted(weights, reverse=True)
    assert len(corrected) == 21 and len(uncorrected) == 20
    # Paper text: TTTTT/AA is rank 21 under the corrected model.
    assert en_motif_site_rank("TTTTTAA", corrected=True) == 21
    assert en_motif_site_weight("TTTTTAA", corrected=True) == pytest.approx(0.2186)
    # ...and rank 1 (weight 1.0) under the uncorrected model.
    assert en_motif_site_rank("TTTTTAA", corrected=False) == 1
    assert en_motif_site_weight("TTTTTAA", corrected=False) == 1.0
    # Top of the corrected table is a C-containing site, not the canonical.
    assert en_motif_site_rank("TTTTCAA", corrected=True) == 1
    assert en_motif_site_weight("TTTTCAA", corrected=True) == 1.0


def test_en_motif_score_tiers():
    # Table tier: TTTTCAA is the top corrected site.
    weight, tier = en_motif_score("GATTTTCAAGC", corrected=True)
    assert tier == "table"
    assert weight == pytest.approx(1.0)
    # Fallback tier: hexamer present but no published 7-mer site.
    weight, tier = en_motif_score("GGTTTTAACC", corrected=True)
    assert tier == "fallback" and weight is None
    # None tier.
    weight, tier = en_motif_score("GGCCGGCCGGCC", corrected=True)
    assert tier == "none" and weight is None
    # Empty input.
    assert en_motif_score("") == (None, "none")


def test_en_motif_score_model_choice_changes_verdict():
    # TTTTTAA: top of the uncorrected table, bottom of the corrected one.
    wu, tu = en_motif_score("TTTTTAA", corrected=False)
    wc, tc = en_motif_score("TTTTTAA", corrected=True)
    assert (tu, wu) == ("table", 1.0)
    assert tc == "table"
    assert wc == pytest.approx(0.2186)


def test_flasch_provenance_roundtrip():
    prov = flasch_weights_provenance()
    assert prov["sha256"].startswith("b9cc818ca1c25ab6")
    assert prov["pmcid"] == "PMC6558663"
    assert "Flasch" in str(prov["source_paper"])
    assert " mmc4" in str(prov["source_file"]) or "mmc4" in str(prov["source_file"])
    assert prov["retrieved"] == "2026-09-30"


def test_scan_en_motif_finds_canonical():
    # Canonical 5'-TTTT/AA-3' written as contiguous TTTTAA on the strand.
    seq = "GATTTTAAGCC"
    hits = scan_en_motif_strand(seq)
    assert hits == [2]  # TTTTAA at 0-based 2
    assert seq[2:8] == "TTTTAA"


def test_scan_en_motif_multiple_and_overlapping():
    # Two non-overlapping copies; TTTTAA cannot self-overlap (T... vs ...A).
    assert scan_en_motif_strand("TTTTAATTTTAA") == [0, 6]
    assert scan_en_motif_strand("TTTTTTAA") == [2]  # only the exact 6-mer
    assert scan_en_motif_strand("CCATTTTT") == []  # only 4 Ts before AA side


def test_scan_en_motif_requires_full_hexamer():
    assert scan_en_motif_strand("TTTTA") == []  # 5 bp, too short
    assert scan_en_motif_strand("") == []


def test_en_motif_opportunity():
    assert en_motif_opportunity("TTTTAA") == pytest.approx(1.0)
    # 2 matches in 7 windows of a 12-mer.
    assert en_motif_opportunity("TTTTAATTTTAA") == pytest.approx(2 / 7)
    assert en_motif_opportunity("GGCCGGCCGGC") == 0.0
    assert en_motif_opportunity("") == 0.0


def test_revcomp_roundtrip():
    assert revcomp("TTTTAAA") == "TTTAAAA"
    assert revcomp(revcomp("ACGTN")) == "ACGTN"


# ---------------------------------------------------------------------------
# 5' junction microhomology
# ---------------------------------------------------------------------------

def test_microhomology_counts_overlap_at_junction():
    # 4-bp overlap between flank end and element start.
    assert microhomology_at_5p_junction("ACGTACGT", "ACGTAAAA") == 4
    assert microhomology_at_5p_junction("ACGTACGT", "TTTTAAAA") == 1  # only the final T


def test_microhomology_prefers_longest():
    flank = "AAAAACGT"
    element = "ACGTACGT"  # both 4-bp ACGT overlap
    assert microhomology_at_5p_junction(flank, element) == 4


def test_microhomology_caps_at_max_len():
    flank = "AAAAAAAAAAAAAAAA"
    element = "AAAAAAAAAAAAAAAA"
    assert microhomology_at_5p_junction(flank, element) == MH_MAX_LEN


def test_microhomology_empty_inputs():
    assert microhomology_at_5p_junction("", "AAA") == 0
    assert microhomology_at_5p_junction("AAA", "") == 0


# ---------------------------------------------------------------------------
# Internal priming (Srikanta signature)
# ---------------------------------------------------------------------------

def test_internal_priming_full_signature_is_suspect():
    sig = internal_priming_signature(
        tsd_len_bp=8,
        poly_a_tail_present=False,
        element_3p_truncated=True,
        en_motif_present=False,
    )
    assert sig.is_suspect
    assert any("internal-priming signature complete" in r for r in sig.reasons)


def test_internal_priming_partial_signature_is_not_suspect():
    # 3 of 4 features — individually common, joint architecture is the test.
    sig = internal_priming_signature(
        tsd_len_bp=8,
        poly_a_tail_present=False,
        element_3p_truncated=True,
        en_motif_present=True,  # canonical EN site present -> not the signature
    )
    assert not sig.is_suspect
    assert any("incomplete" in r for r in sig.reasons)


def test_internal_priming_short_tsd_blocks():
    sig = internal_priming_signature(
        tsd_len_bp=4,  # < 6 bp: ambiguous with tail polyA
        poly_a_tail_present=False,
        element_3p_truncated=True,
        en_motif_present=False,
    )
    assert not sig.is_suspect


def test_internal_priming_none_blocks():
    sig = internal_priming_signature(
        tsd_len_bp=None,
        poly_a_tail_present=False,
        element_3p_truncated=True,
        en_motif_present=False,
    )
    assert not sig.is_suspect


# ---------------------------------------------------------------------------
# Complex-event flags
# ---------------------------------------------------------------------------

def test_flag_complex_inversion_and_deletion():
    flags = flag_complex_event(inserted_len_bp=900, expected_len_bp=6000, inversion_fraction=0.25)
    assert ComplexEventKind.INVERSION in flags.kinds
    assert ComplexEventKind.DELETION in flags.kinds  # 900 < 0.3 * 6000


def test_flag_complex_negative_inputs_stay_clean():
    flags = flag_complex_event(inserted_len_bp=6000, expected_len_bp=6000, inversion_fraction=0.0)
    assert flags.kinds == frozenset()


def test_flag_complex_en_and_repeat_negative_pair():
    only_en = flag_complex_event(en_site_negative=True, repeat_negative=False)
    assert ComplexEventKind.EN_AND_REPEAT_NEGATIVE not in only_en.kinds
    both = flag_complex_event(en_site_negative=True, repeat_negative=True)
    assert ComplexEventKind.EN_AND_REPEAT_NEGATIVE in both.kinds


def test_flag_chimeric_and_recombination():
    flags = flag_complex_event(chimeric=True, recombination_like=True)
    assert ComplexEventKind.CHIMERIC in flags.kinds
    assert ComplexEventKind.RECOMBINATION_LIKE in flags.kinds


# ---------------------------------------------------------------------------
# Alu hotspot scoring
# ---------------------------------------------------------------------------

def test_hotspot_likelihood_inside_and_outside():
    assert levy_hotspot_likelihood(130) == 1.0  # linker 118-140
    assert levy_hotspot_likelihood(290) == 1.0  # 3p tail 280-300
    assert levy_hotspot_likelihood(200) == 0.0
    assert levy_hotspot_likelihood(160) == 0.0  # > 10 bp falloff past linker


def test_hotspot_falloff_is_linear():
    edge = levy_hotspot_likelihood(145)  # 5 bp past the linker end, falloff 10
    assert edge == pytest.approx(0.4)
    assert levy_hotspot_likelihood(290.5) == 1.0  # fractional pos inside interval


def test_hotspot_table_sources():
    names = {hs.name for hs in ALU_HOTSPOT_POSITIONS}
    assert names == {"a_rich_linker", "3p_a_tail"}
    for hs in ALU_HOTSPOT_POSITIONS:
        assert hs.source  # every named hotspot carries its provenance


# ---------------------------------------------------------------------------
# Route classification
# ---------------------------------------------------------------------------

def test_route_clean_tprt():
    ev = RouteEvidence(
        family="ALU",
        tsd_len_bp=14,
        en_motif_present=True,
        poly_a_tail_present=True,
        element_5p_truncated=False,
        microhomology_bp=0,
    )
    result = classify_route(ev)
    assert result.route == "clean_tprt"
    assert any("triad" in r for r in result.reasons)


def test_route_l1_truncated_with_microhomology():
    ev = RouteEvidence(
        family="LINE1",
        element_5p_truncated=True,
        microhomology_bp=5,
        tsd_len_bp=16,
        poly_a_tail_present=True,
    )
    result = classify_route(ev)
    assert result.route == "l1_truncated_route"
    assert any("Zingler" in r for r in result.reasons)


def test_route_l1_truncated_no_tail():
    ev = RouteEvidence(
        family="LINE1",
        element_5p_truncated=True,
        microhomology_bp=0,
        poly_a_tail_present=False,
    )
    assert classify_route(ev).route == "l1_truncated_route"


def test_route_truncated_alu_without_mh_stays_unknown():
    # Zingler: Alu is NOT microhomology-dependent, so a truncated Alu without
    # microhomology must not be stamped l1_truncated_route.
    ev = RouteEvidence(
        family="ALU",
        element_5p_truncated=True,
        microhomology_bp=0,
        poly_a_tail_present=True,
    )
    assert classify_route(ev).route == "unknown"


def test_route_internal_priming_wins_priority():
    ev = RouteEvidence(
        family="LINE1",
        tsd_len_bp=10,
        poly_a_tail_present=False,
        element_3p_truncated=True,
        en_motif_present=False,
        element_5p_truncated=True,
    )
    assert classify_route(ev).route == "internal_priming_suspect"


def test_route_complex_flags_beat_tprt():
    ev = RouteEvidence(
        family="ALU",
        tsd_len_bp=14,
        en_motif_present=True,
        poly_a_tail_present=True,
        complex_kinds=frozenset({ComplexEventKind.INVERSION}),
    )
    assert classify_route(ev).route == "complex"


def test_route_unknown_when_evidence_missing():
    result = classify_route(RouteEvidence(family="ALU"))
    assert result.route == "unknown"


def test_route_evidence_to_dict_stable_kinds():
    ev = RouteEvidence(
        family="LINE1",
        complex_kinds=frozenset({ComplexEventKind.INVERSION, ComplexEventKind.DELETION}),
    )
    d = ev.to_dict()
    assert d["complex_kinds"] == "deletion,inversion"  # sorted, stable


def test_route_classification_rejects_bad_route():
    from retro_miner.mei_markers import RouteClassification

    with pytest.raises(ValueError):
        RouteClassification(route="not_a_route", reasons=())
