"""Scientific test suite for the micro-homology aware TSD refiner.

Synthetic TPRT (Target-Primed Reverse Transcription) events are simulated
against a tiny reference genome FASTA with soft-clipped supporting reads
built as ``pysam.AlignedSegment`` objects:

* Canonical 15-bp TSD integration (``5'-ATTGCA...ATTGCA-3'`` footprint with
  a ``3'`` poly(A) tail).
* Zero-TSD blunt integration and insertion-site deletion.
* Micro-homology resolution where a shared ``AGCT`` motif makes the
  breakpoint ambiguous.
* Poly(A)/poly(T) artifact rejection so a homopolymer flank is never
  reported as a real duplication.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pysam

from retro_miner.tsd_refiner import (
    TSDResult,
    calculate_sequence_entropy,
    is_poly_a_or_t,
    refine_tsd_boundaries,
)

# Repeating/flank motifs chosen to avoid accidental junction matches.
M1 = "GTCATCGGTACGAACCTGTA"
M2 = "CTAGGCATTACGTACCGAGT"
INSERT5_HEAD = "CGTTACCTGAGTACGATCAGGACTGCATCGATCGTACCAA"  # 40 bp non-poly 5' head
TSD15 = "ATTGCAGCTAGCTAG"


# ─────────────────────────────────────────────────────────────────────────────
# Fixture helpers
# ─────────────────────────────────────────────────────────────────────────────


def write_reference(path: Path, chrom: str, sequence: str) -> None:
    path.write_text(f">{chrom}\n" + "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)) + "\n")
    pysam.faidx(str(path))


def build_read(
    query_seq: str,
    cigar: list[tuple[int, int]],
    reference_start: int,
    reference_id: int = 0,
) -> pysam.AlignedSegment:
    """Construct an AlignedSegment without setting the read-only reference_end."""
    read = pysam.AlignedSegment()
    read.query_name = "tsd_test"
    read.query_sequence = query_seq
    read.flag = 0
    read.reference_id = reference_id
    read.reference_start = reference_start
    read.mapping_quality = 60
    read.cigartuples = cigar
    # reference_end is computed by pysam from reference_start + cigar.
    return read


def build_tprt_reads(
    ref_seq: str,
    left_bp: int,
    right_bp: int,
    poly_a_run: int = 20,
    insert_head: str = INSERT5_HEAD,
) -> tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """Build a left-anchored read (5' junction) and right-anchored read (3' junction)."""
    assert left_bp >= 1 and right_bp >= 0

    left_query = ref_seq[:left_bp] + insert_head
    left_read = build_read(
        left_query,
        [(0, left_bp), (4, len(insert_head))],
        reference_start=0,
    )

    aligned_len = 60
    right_ref_start = right_bp
    assert right_ref_start + aligned_len <= len(ref_seq)
    right_query = "A" * poly_a_run + ref_seq[right_ref_start : right_ref_start + aligned_len]
    right_read = build_read(
        right_query,
        [(4, poly_a_run), (0, aligned_len)],
        reference_start=right_ref_start,
    )
    return left_read, right_read


# ─────────────────────────────────────────────────────────────────────────────
# Shannon entropy and poly(A)/poly(T) discrimination
# ─────────────────────────────────────────────────────────────────────────────

def test_shannon_entropy_homopolymer_vs_genomic():
    assert calculate_sequence_entropy("AAAAAAAAAAAAAAAA") == 0.0
    assert calculate_sequence_entropy("TTTTTTTTTT") == 0.0
    assert calculate_sequence_entropy("ACGTACGTACGTACGT") == 2.0
    assert calculate_sequence_entropy("") == 0.0


def test_poly_a_tail_detection():
    assert is_poly_a_or_t("AAAAAAAAAAAAAAAT") is True
    assert is_poly_a_or_t("TTTTTTTTTTTTTTTG") is True
    assert is_poly_a_or_t("CAGTAGCTAGCTAGCT") is False
    assert is_poly_a_or_t("AAAA") is False  # too short for a tail


def test_entropy_discriminates_poly_from_genomic():
    # Poly(A) with a single wobble has low entropy; balanced 4-mer has H = 2.
    assert calculate_sequence_entropy("AAAAAAAAAAAAAAAA") < 1.2
    assert calculate_sequence_entropy("ACGTACGTACGTACGT") >= 1.2


# ─────────────────────────────────────────────────────────────────────────────
# Canonical 15-bp TSD TPRT event
# ─────────────────────────────────────────────────────────────────────────────

def test_canonical_15bp_tprt_tsd_refinement(tmp_path):
    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + len(TSD15)  # 75 (0-based exclusive 5' junction)
    right_bp = len(M1) * 3              # 60 (0-based inclusive 3' junction)
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=20)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert isinstance(res, TSDResult)
    assert res.tsd_seq == TSD15
    assert res.tsd_length == 15
    assert res.polyA_tail_detected is True
    assert res.poly_tail_base == "A"
    assert res.poly_tail_length >= 20
    assert res.tsd_confidence_score >= 0.90
    assert res.method == "tsd_flush_match"
    assert res.left_breakpoint == left_bp
    assert res.right_breakpoint == right_bp
    assert res.n_left_anchored == 1
    assert res.n_right_anchored == 1


def test_tsd_confidence_lower_without_poly_tail(tmp_path):
    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + len(TSD15)
    right_bp = len(M1) * 3
    pos = left_bp + 1
    # Right read uses a non-poly insert head instead of a poly(A) tail.
    right_query = INSERT5_HEAD + ref_seq[right_bp : right_bp + 60]
    right_read = build_read(
        right_query,
        [(4, len(INSERT5_HEAD)), (0, 60)],
        reference_start=right_bp,
    )
    left_query = ref_seq[:left_bp] + INSERT5_HEAD
    left_read = build_read(
        left_query,
        [(0, left_bp), (4, len(INSERT5_HEAD))],
        reference_start=0,
    )

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == TSD15
    assert res.tsd_length == 15
    assert res.polyA_tail_detected is False
    assert res.tsd_confidence_score == 0.85


# ─────────────────────────────────────────────────────────────────────────────
# Zero-TSD blunt integration and deletion at the insertion site
# ─────────────────────────────────────────────────────────────────────────────

def test_blunt_integration_zero_tsd(tmp_path):
    ref_seq = M1 * 3 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3   # 60
    right_bp = len(M1) * 3  # 60 -> both junctions coincide, no TSD
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=16)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == ""
    assert res.tsd_length == 0
    assert res.method == "blunt"
    assert res.polyA_tail_detected is True
    assert res.left_breakpoint == right_bp
    assert res.right_breakpoint == left_bp


def test_deletion_at_insertion_site(tmp_path):
    ref_seq = M1 * 3 + M2 * 3 + M1 * 3  # padded so the downstream read fits
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3   # 60
    right_bp = left_bp + 15  # 75 -> 15 bp of reference deleted between junctions
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=16)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == ""
    assert res.tsd_length == 0
    assert res.method == "deletion"
    assert res.polyA_tail_detected is True


# ─────────────────────────────────────────────────────────────────────────────
# Micro-homology resolution: shared AGCT motif
# ─────────────────────────────────────────────────────────────────────────────

def test_microhomology_agct_overlap_resolution(tmp_path):
    # The AGCT motif is shared between the donor end and the target site, so
    # the naive breakpoint overlap is a tandem "AGCTAGCT" that collapses to
    # the 4-bp minimal period.
    ref_seq = M1 * 3 + "AGCTAGCT" + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + 8
    right_bp = len(M1) * 3
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=16)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == "AGCT"
    assert res.tsd_length == 4
    assert res.method == "microhomology"
    assert res.microhomology_seq == "AGCT"
    assert res.tsd_confidence_score == 0.80
    assert res.polyA_tail_detected is True


def test_poly_tail_microhomology_keeps_tail_downstream(tmp_path):
    # A homopolymer A-run shared between the target and the poly(A) tail must
    # be resolved by placing the tail downstream: the leading run is trimmed
    # off the resolved TSD (which starts with a non-A base here).
    tsd_no_a_start = "GCATGCAAGTCGATCGA"
    ref_seq = M1 * 3 + "AAA" + tsd_no_a_start + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + 3 + len(tsd_no_a_start)
    right_bp = len(M1) * 3
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=20)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == tsd_no_a_start
    assert res.tsd_length == len(tsd_no_a_start)
    assert res.polyA_tail_detected is True
    assert res.method == "microhomology"
    assert res.microhomology_seq == "AAA"


def test_poly_a_homopolymer_flank_is_artifact(tmp_path):
    # A pure poly(A) overlap is a tail artifact, never a TSD.
    ref_seq = M1 * 3 + "A" * 16 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + 16
    right_bp = len(M1) * 3
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=16)

    res = refine_tsd_boundaries("chr1", pos, [left_read, right_read], fasta)

    assert res.tsd_seq == ""
    assert res.tsd_length == 0
    assert res.method == "polyA_artifact"
    assert res.polyA_tail_detected is True


# ─────────────────────────────────────────────────────────────────────────────
# Degenerate and single-arm cases
# ─────────────────────────────────────────────────────────────────────────────

def test_no_soft_clips_returns_empty_tsd(tmp_path):
    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    res = refine_tsd_boundaries("chr1", 76, [], fasta)

    assert res.method == "no_soft_clips"
    assert res.tsd_seq == ""
    assert res.tsd_length == 0
    assert res.polyA_tail_detected is False


def test_single_left_arm_is_low_confidence(tmp_path):
    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + len(TSD15)
    pos = left_bp + 1
    left_query = ref_seq[:left_bp] + INSERT5_HEAD
    left_read = build_read(
        left_query,
        [(0, left_bp), (4, len(INSERT5_HEAD))],
        reference_start=0,
    )

    res = refine_tsd_boundaries("chr1", pos, [left_read], fasta)

    assert res.method == "left_arm_only"
    assert res.tsd_seq == ""
    assert res.tsd_length == 0
    assert res.tsd_confidence_score == 0.20


def test_missing_reference_contig_returns_empty_tsd(tmp_path):
    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + len(TSD15)
    right_bp = len(M1) * 3
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=20)

    res = refine_tsd_boundaries("chr2", pos, [left_read, right_read], fasta)

    assert res.method == "no_reference_contig"
    assert res.tsd_seq == ""
    assert res.tsd_length == 0


# ─────────────────────────────────────────────────────────────────────────────
# Exact k-mer flush match primitive
# ─────────────────────────────────────────────────────────────────────────────

def test_max_flush_match_detects_shared_kmer():
    from retro_miner.tsd_refiner import _max_flush_match

    left = "CGTATGATCGCGAT"
    right = "GCGATTTAGGCCTA"
    assert _max_flush_match(left, right) == 5

    assert _max_flush_match("AAAACCCC", "GGGGTTTT") == 0
    assert _max_flush_match("CCCCAGCT", "AGCTAGCTAA") == 4


# ─────────────────────────────────────────────────────────────────────────────
# candidate_loci pipeline integration
# ─────────────────────────────────────────────────────────────────────────────

def test_candidate_loci_tsd_annotation_integration(tmp_path) -> None:
    from retro_miner.candidate_loci import _annotate_tsd_refinement

    ref_seq = M1 * 3 + TSD15 + M2 * 3
    fasta = tmp_path / "ref.fa"
    write_reference(fasta, "chr1", ref_seq)

    left_bp = len(M1) * 3 + len(TSD15)
    right_bp = len(M1) * 3
    pos = left_bp + 1
    left_read, right_read = build_tprt_reads(ref_seq, left_bp, right_bp, poly_a_run=20)

    bam_path = tmp_path / "support.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(ref_seq)}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        out.write(left_read)
        out.write(right_read)
    pysam.index(str(bam_path))

    loci = pd.DataFrame([{"chrom": "chr1", "window_start": 1, "window_end": 200}])
    split_disease = pd.DataFrame(
        [{"chrom": "chr1", "pos": pos, "read_name": "x", "mapq": 60, "clip_len": 20, "has_sa": False,
          "window_start": 1, "window_end": 200}]
    )
    split_control = pd.DataFrame(
        columns=["chrom", "pos", "read_name", "mapq", "clip_len", "has_sa", "window_start", "window_end"]
    )

    out = _annotate_tsd_refinement(
        loci=loci,
        split_disease=split_disease,
        split_control=split_control,
        bam_path=Path(bam_path),
        reference_fasta=Path(fasta),
        tsd_refine_flank_bp=60,
    )

    assert "tsd_seq" in out.columns
    assert "tsd_length" in out.columns
    assert "tsd_confidence_score" in out.columns
    assert "polyA_tail_detected" in out.columns
    assert out.iloc[0]["tsd_seq"] == TSD15
    assert out.iloc[0]["tsd_length"] == 15
    assert out.iloc[0]["tsd_confidence_score"] >= 0.90
    assert bool(out.iloc[0]["polyA_tail_detected"]) is True


def test_candidate_loci_tsd_defaults_without_bam(tmp_path) -> None:
    from retro_miner.candidate_loci import _annotate_tsd_refinement

    loci = pd.DataFrame([{"chrom": "chr1", "window_start": 1, "window_end": 200}])
    split_disease = pd.DataFrame(
        [{"chrom": "chr1", "pos": 76, "read_name": "x", "mapq": 60, "clip_len": 20, "has_sa": False,
          "window_start": 1, "window_end": 200}]
    )
    out = _annotate_tsd_refinement(
        loci=loci,
        split_disease=split_disease,
        split_control=pd.DataFrame(),
        bam_path=None,
        reference_fasta=Path(tmp_path / "ref.fa"),
        tsd_refine_flank_bp=60,
    )

    assert out.iloc[0]["tsd_seq"] == ""
    assert out.iloc[0]["tsd_length"] == 0
    assert out.iloc[0]["tsd_confidence_score"] == 0.0
    assert bool(out.iloc[0]["polyA_tail_detected"]) is False