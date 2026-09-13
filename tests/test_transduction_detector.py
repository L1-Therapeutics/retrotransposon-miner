"""Biological test suite for L1 3' transduction detection (Tubio et al. 2015).

Synthetic scenarios exercised:
1. Intact L1HS 3' UTR + poly(A) tail + 45 bp ``GATTACA...`` unique segment
   downstream -> ``3_PRIME_PARTNERED``, high confidence, exact extraction.
2. Poly(A) tail extension without transduction -> ``has_transduction == False``.
3. Poly(A) tail followed by low-complexity repeat -> gated out (entropy < 1.5).
4. ``5'``-truncated L1 body (no 3' UTR signal) + poly(A) + transduction ->
   ``3_PRIME_ORPHAN``.
5. Reverse-strand contig: minus-strand polyadenylation signal (``TTTATT``) +
   poly(T) tail -> ``3_PRIME_PARTNERED``.
6. Truncated downstream (< ``min_transduction_len``) -> ``NONE``.
7. High-complexity downstream directly abutting the L1 (no poly(A) tail) ->
   retained as ``3_PRIME_PARTNERED`` at reduced confidence.
"""

from retro_miner.transduction_detector import (
    TransductionResult,
    detect_3prime_transduction,
)

L1_AT_PADDING = (
    "AATTTTGTACAAAGTTACTTCCTCATCCCTCCCTTTCTTCTCATTCCTCTTCTCCAATCT"
)

L1HS_3UTR_CORE = "TGAGGGGCCTGTGCTTCTCCTTCAAGGAAATAAAATTAACCTTGAGCT"
L1HS_3UTR_BODY = L1_AT_PADDING + L1HS_3UTR_CORE

L1HS_3UTR_BODY_MINUS = L1_AT_PADDING + "TGAGGGGCCTGTGCTTCTCCTTCAAGGAAATTTTATTTAACCTTGAGCT"

L1_TRUNCATED_BODY = (
    "TAGGCCATGCCAGTACTCACATAGACACACACACAAAAGTACACACACACACACAGGATC"
    "CCTCAGCACTTTTCCCACCTTAGACATAG"
)

TRANSDUCED_45 = "GATTACA" * 6 + "GAT"
POLY_A_TAIL = "AAAAAA"
POLY_T_TAIL = "TTTTTT"


def _partnered_contig(tail: str = POLY_A_TAIL) -> str:
    return L1HS_3UTR_BODY + tail + TRANSDUCED_45


def _orphan_contig(tail: str = POLY_A_TAIL) -> str:
    return L1_TRUNCATED_BODY + tail + TRANSDUCED_45


def test_partnered_transduction_detected_exactly():
    contig = _partnered_contig()
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.transduction_length == 45
    assert res.confidence_score >= 0.90


def test_partnered_transduction_with_trailing_flank():
    flank = "ACGTACGTAACGTCGATCTGATC"
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + TRANSDUCED_45 + flank
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq.startswith(TRANSDUCED_45)
    assert res.transduction_length >= 45
    assert res.confidence_score >= 0.90


def test_pure_polyA_extension_is_not_transduction():
    contig = L1HS_3UTR_BODY + "A" * 40
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"
    assert res.confidence_score == 0.0
    assert res.transduction_seq == ""


def test_low_complexity_downstream_gated_out():
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + "AG" * 30
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_orphan_transduction_5prime_truncated_l1():
    contig = _orphan_contig()
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1_TRUNCATED_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_ORPHAN"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.transduction_length == 45
    assert res.confidence_score == 0.85


def test_reverse_strand_polyT_partnered_transduction():
    contig = L1HS_3UTR_BODY_MINUS + POLY_T_TAIL + TRANSDUCED_45
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY_MINUS))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.confidence_score == 0.95


def test_short_downstream_below_min_length_is_none():
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + TRANSDUCED_45[:7]
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_no_poly_tail_low_confidence_partnered():
    contig = L1HS_3UTR_BODY + TRANSDUCED_45
    res = detect_3prime_transduction(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.confidence_score == 0.80


def test_min_length_param_respected():
    contig = _partnered_contig()
    res = detect_3prime_transduction(
        contig,
        mei_alignment_end=len(L1HS_3UTR_BODY),
        min_transduction_len=100,
    )
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_empty_contig_no_transduction():
    res = detect_3prime_transduction("", 0)
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"
    assert res.transduction_length == 0
    assert res.confidence_score == 0.0


def test_result_is_frozen_dataclass():
    res = detect_3prime_transduction(_partnered_contig(), len(L1HS_3UTR_BODY))
    assert isinstance(res, TransductionResult)
    try:
        res.transduction_seq = "X"
        raised = False
    except Exception:
        raised = True
    assert raised