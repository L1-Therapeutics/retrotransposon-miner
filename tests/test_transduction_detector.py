"""Biological test suite for L1 3' transduction detection.

Contig-level classification (``detect_3prime_transduction_from_contig``),
synthetic scenarios exercised:
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

Read-level donor profiling (``detect_3prime_transduction``), synthetic
scenarios exercised:
8. Chimeric read of ~100 bp L1HS 3' UTR + poly(A) tail + 40 bp unique genomic
   sequence -> unique non-MEI tail detected and mapped back to its donor locus.
9. Pure poly(A)/poly(T) reads -> filtered (no transduction).
10. MEI-consensus-only reads -> stripped and filtered.
11. Alu-masked reads (non-MEI clips input) -> non-MEI tail still recovered.
12. Reverse-orientation reads -> donor localised on the matching strand.
"""

from retro_miner.transduction_detector import (
    MEI_CONSENSUS_FRAGMENTS,
    TransductionResult,
    detect_3prime_transduction,
    detect_3prime_transduction_from_contig,
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

# Unique 40 bp "transduction fingerprint" used by the read-level profiler
# tests.  Avoids any >=5 bp homopolymer runs and begins/ends in G/C so the
# synthetic poly(A) tract is unambiguously the read-through boundary.
TRANSDUCED_40_UNIQUE = "CCTCCAAACTTATACGATGCTCACTGTCGTACCTAAACGC"
TRANSDUCED_40_ABSENT = "CTCCGTCGAGCAGAAGCTTGTTTGACAGTTCGCGAGCCTG"

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _partnered_contig(tail: str = POLY_A_TAIL) -> str:
    return L1HS_3UTR_BODY + tail + TRANSDUCED_45


def _orphan_contig(tail: str = POLY_A_TAIL) -> str:
    return L1_TRUNCATED_BODY + tail + TRANSDUCED_45


def test_partnered_transduction_detected_exactly():
    contig = _partnered_contig()
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.transduction_length == 45
    assert res.confidence_score >= 0.90


def test_partnered_transduction_with_trailing_flank():
    flank = "ACGTACGTAACGTCGATCTGATC"
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + TRANSDUCED_45 + flank
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq.startswith(TRANSDUCED_45)
    assert res.transduction_length >= 45
    assert res.confidence_score >= 0.90


def test_pure_polyA_extension_is_not_transduction():
    contig = L1HS_3UTR_BODY + "A" * 40
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"
    assert res.confidence_score == 0.0
    assert res.transduction_seq == ""


def test_low_complexity_downstream_gated_out():
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + "AG" * 30
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_orphan_transduction_5prime_truncated_l1():
    contig = _orphan_contig()
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1_TRUNCATED_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_ORPHAN"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.transduction_length == 45
    assert res.confidence_score == 0.85


def test_reverse_strand_polyT_partnered_transduction():
    contig = L1HS_3UTR_BODY_MINUS + POLY_T_TAIL + TRANSDUCED_45
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY_MINUS))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.transduction_seq == TRANSDUCED_45
    assert res.confidence_score == 0.95


def test_short_downstream_below_min_length_is_none():
    contig = L1HS_3UTR_BODY + POLY_A_TAIL + TRANSDUCED_45[:7]
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_no_poly_tail_low_confidence_partnered():
    contig = L1HS_3UTR_BODY + TRANSDUCED_45
    res = detect_3prime_transduction_from_contig(contig, mei_alignment_end=len(L1HS_3UTR_BODY))
    assert res.has_transduction is True
    assert res.transduction_type == "3_PRIME_PARTNERED"
    assert res.confidence_score == 0.80


def test_min_length_param_respected():
    contig = _partnered_contig()
    res = detect_3prime_transduction_from_contig(
        contig,
        mei_alignment_end=len(L1HS_3UTR_BODY),
        min_transduction_len=100,
    )
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"


def test_empty_contig_no_transduction():
    res = detect_3prime_transduction_from_contig("", 0)
    assert res.has_transduction is False
    assert res.transduction_type == "NONE"
    assert res.transduction_length == 0
    assert res.confidence_score == 0.0


def test_result_is_frozen_dataclass():
    res = detect_3prime_transduction_from_contig(_partnered_contig(), len(L1HS_3UTR_BODY))
    assert isinstance(res, TransductionResult)
    try:
        res.transduction_seq = "X"
        raised = False
    except Exception:
        raised = True
    assert raised


def test_auto_mei_alignment_end_detects_polyA_junction():
    contig = L1HS_3UTR_BODY + "A" * 8 + TRANSDUCED_40_UNIQUE
    res = detect_3prime_transduction_from_contig(contig)
    assert res.has_transduction is True
    assert res.transduction_seq == TRANSDUCED_40_UNIQUE


def _write_synthetic_reference(tmp_path) -> str:
    """Write a small deterministic reference with two donor loci.

    The 40 bp unique tail is placed once in the plus orientation at
    ``chr22:150-190`` and once in reverse-complement orientation at
    ``chr22:400-440`` (modeling a minus-strand donor locus); the rest is a
    deterministic pseudo-random filler.
    """
    import random

    ref = tmp_path / "synthetic_ref.fa"
    rng = random.Random(11)
    seq = [rng.choice("ACGT") for _ in range(700)]
    seq[150:190] = list(TRANSDUCED_40_UNIQUE)
    seq[400:440] = list(_rc(TRANSDUCED_40_UNIQUE))
    ref.write_text(">chr22\n{}\n".format("".join(seq)), encoding="utf-8")
    return str(ref)


def test_synthetic_chimeric_read_with_l1hs_utr_and_unique_tail(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    read = L1HS_3UTR_BODY + "A" * 8 + TRANSDUCED_40_UNIQUE
    call = detect_3prime_transduction([read], [], ref_fasta)
    assert call.has_transduction is True
    assert call.tail_sequence == TRANSDUCED_40_UNIQUE
    assert call.tail_length == 40
    assert call.entropy > 1.8
    assert call.poly_a_trimmed_bp == 8
    assert call.source_read_count == 1
    assert call.best_donor_locus is not None
    assert call.best_donor_locus.chrom == "chr22"
    assert call.best_donor_locus.donor_start == 150
    assert call.best_donor_locus.donor_end == 190
    assert call.best_donor_locus.strand == "+"
    assert call.best_donor_locus.identity == 1.0
    loci = {(locus.chrom, locus.donor_start, locus.strand) for locus in call.donor_loci}
    assert ("chr22", 150, "+") in loci
    assert ("chr22", 400, "-") in loci


def test_pure_polyA_reads_are_filtered(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    call = detect_3prime_transduction(["A" * 60, "T" * 60], ["AAAAATTTTAAAA"], ref_fasta)
    assert call.has_transduction is False
    assert call.tail_sequence == ""
    assert call.donor_loci == ()


def test_mei_consensus_only_reads_are_filtered(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    l1_only = L1HS_3UTR_BODY + "A" * 8
    call = detect_3prime_transduction([l1_only], [], ref_fasta)
    assert call.has_transduction is False


def test_alu_masked_read_still_yields_tail(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    alu_head = MEI_CONSENSUS_FRAGMENTS[1]
    read = alu_head + "A" * 8 + TRANSDUCED_40_UNIQUE
    call = detect_3prime_transduction([], [read], ref_fasta)
    assert call.has_transduction is True
    assert call.tail_sequence == TRANSDUCED_40_UNIQUE
    assert call.best_donor_locus is not None
    assert call.best_donor_locus.donor_start == 150


def test_reverse_strand_tail_maps_to_minus_donor(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    read = L1HS_3UTR_BODY + "A" * 7 + _rc(TRANSDUCED_40_UNIQUE)
    call = detect_3prime_transduction([read], [], ref_fasta)
    assert call.has_transduction is True
    assert call.tail_sequence == _rc(TRANSDUCED_40_UNIQUE)
    assert call.best_donor_locus is not None
    loci = {(locus.chrom, locus.donor_start, locus.strand) for locus in call.donor_loci}
    assert ("chr22", 400, "+") in loci
    assert ("chr22", 150, "-") in loci


def test_read_stored_in_reverse_orientation_is_detected(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    read = _rc(L1HS_3UTR_BODY + "A" * 8 + TRANSDUCED_40_UNIQUE)
    call = detect_3prime_transduction([read], [], ref_fasta)
    assert call.has_transduction is True
    assert call.tail_sequence == TRANSDUCED_40_UNIQUE
    assert call.best_donor_locus is not None
    assert call.best_donor_locus.donor_start == 150
    assert call.best_donor_locus.strand == "+"


def test_unique_tail_absent_from_reference_still_reports_candidate(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    read = L1HS_3UTR_BODY + "A" * 8 + TRANSDUCED_40_ABSENT
    call = detect_3prime_transduction([read], [], ref_fasta)
    assert call.has_transduction is True
    assert call.tail_sequence == TRANSDUCED_40_ABSENT
    assert call.tail_length == 40
    assert call.donor_loci == ()


def test_empty_read_inputs_return_no_call(tmp_path):
    ref_fasta = _write_synthetic_reference(tmp_path)
    call = detect_3prime_transduction([], [], ref_fasta)
    assert call.has_transduction is False
    assert call.best_donor_locus is None