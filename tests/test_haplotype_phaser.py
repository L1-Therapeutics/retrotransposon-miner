"""Unit tests for read-spanning haplotype linkage & phase estimation."""

import pysam

from retro_miner.haplotype_phaser import (
    PHASE_BLOCK_H1,
    PHASE_BLOCK_H2,
    PHASE_UNPHASED,
    PhaseLinkageResult,
    phase_mei_locus,
)
from retro_miner.vcf_export import write_mei_vcf

CHROM = "chr1"
BREAKPOINT = 5000
SNP_POS = 5100
SNP_REF = "G"
SNP_ALT = "A"


def _mei_read(seq_index: int) -> pysam.AlignedSegment:
    seq = "T" * 20 + "A" * 300 + "T" * 130
    bases = list(seq)
    bases[420] = SNP_ALT  # base at BREAKPOINT-aligned SNP_POS
    read = pysam.AlignedSegment()
    read.query_name = f"mei_read_{seq_index}"
    read.query_sequence = "".join(bases)
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 4980
    read.mapping_quality = 60
    read.cigar = ((0, 20), (1, 300), (0, 130))  # 300 bp MEI insertion
    return read


def _ref_read(seq_index: int) -> pysam.AlignedSegment:
    seq = "T" * 150
    bases = list(seq)
    bases[120] = SNP_REF  # base at BREAKPOINT-aligned SNP_POS
    read = pysam.AlignedSegment()
    read.query_name = f"ref_read_{seq_index}"
    read.query_sequence = "".join(bases)
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 4980
    read.mapping_quality = 60
    read.cigar = ((0, 150),)
    return read


def _write_synth_bam(path, mei_reads, ref_reads):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 20000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for read in mei_reads + ref_reads:
            out.write(read)
    sorted_path = path.with_name(path.stem + ".sorted.bam")
    pysam.sort("-o", str(sorted_path), str(path))
    pysam.index(str(sorted_path))
    return sorted_path


def test_mei_reads_linked_to_snp_alt_haplotype_h1(tmp_path):
    bam_path = _write_synth_bam(
        tmp_path / "h1_linked.bam",
        [_mei_read(i) for i in range(10)],
        [_ref_read(i) for i in range(10)],
    )

    result = phase_mei_locus(
        bam_path=str(bam_path),
        chrom=CHROM,
        pos=BREAKPOINT,
        het_snps=[(SNP_POS, SNP_REF, SNP_ALT)],
        window_bp=1000,
    )

    assert result.haplotype_assigned == PHASE_BLOCK_H1
    assert result.linked_snps_count == 1
    assert result.phase_confidence >= 13.0
    assert result.phase_block_id == f"{CHROM}:4000-6000"


def test_reference_reads_linked_to_snp_ref_haplotype_h2(tmp_path):
    def _mei_carrying_ref(seq_index):
        seq = "T" * 20 + "A" * 300 + "T" * 130
        bases = list(seq)
        bases[420] = SNP_REF
        read = pysam.AlignedSegment()
        read.query_name = f"mei_ref_snp_{seq_index}"
        read.query_sequence = "".join(bases)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 4980
        read.mapping_quality = 60
        read.cigar = ((0, 20), (1, 300), (0, 130))
        return read

    def _ref_carrying_alt(seq_index):
        seq = "T" * 150
        bases = list(seq)
        bases[120] = SNP_ALT
        read = pysam.AlignedSegment()
        read.query_name = f"ref_alt_snp_{seq_index}"
        read.query_sequence = "".join(bases)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 4980
        read.mapping_quality = 60
        read.cigar = ((0, 150),)
        return read

    bam_path = _write_synth_bam(
        tmp_path / "h2_linked.bam",
        [_mei_carrying_ref(i) for i in range(10)],
        [_ref_carrying_alt(i) for i in range(10)],
    )

    result = phase_mei_locus(
        bam_path=str(bam_path),
        chrom=CHROM,
        pos=BREAKPOINT,
        het_snps=[(SNP_POS, SNP_REF, SNP_ALT)],
        window_bp=1000,
    )

    assert result.haplotype_assigned == PHASE_BLOCK_H2
    assert result.linked_snps_count == 1
    assert result.phase_confidence >= 13.0


def test_no_het_snps_returns_unphased(tmp_path):
    bam_path = _write_synth_bam(
        tmp_path / "no_snps.bam",
        [_mei_read(i) for i in range(10)],
        [_ref_read(i) for i in range(10)],
    )

    result = phase_mei_locus(
        bam_path=str(bam_path),
        chrom=CHROM,
        pos=BREAKPOINT,
        het_snps=[],
        window_bp=1000,
    )

    assert result.haplotype_assigned == PHASE_UNPHASED
    assert result.linked_snps_count == 0
    assert result.phase_confidence == 0.0


def test_missing_bam_returns_unphased():
    result = phase_mei_locus(
        bam_path="/nonexistent/path/reads.bam",
        chrom=CHROM,
        pos=BREAKPOINT,
        het_snps=[(SNP_POS, SNP_REF, SNP_ALT)],
        window_bp=1000,
    )

    assert result.haplotype_assigned == PHASE_UNPHASED
    assert result.phase_confidence == 0.0


def test_phase_fields_emitted_in_vcf(tmp_path):
    record = {
        "chrom": CHROM,
        "pos": BREAKPOINT,
        "ref_base": "N",
        "family": "L1",
        "k_alt": 10,
        "k_ref": 8,
        "supporting_reads": 10,
        "tsd_seq": "TTTAAA",
        "poly_a_detected": 1,
        "subfamily": "L1HS",
        "mei_llr": 3.2,
        "transduction_type": "NONE",
        "phase_linkage": PhaseLinkageResult(
            haplotype_assigned=PHASE_BLOCK_H1,
            linked_snps_count=1,
            phase_confidence=42.0,
            phase_block_id="chr1:4000-6000",
        ),
    }
    out_path = write_mei_vcf([record], str(tmp_path / "out.vcf"))

    text = out_path.read_text()
    assert "##FORMAT=<ID=HP" in text
    assert "##FORMAT=<ID=PQ" in text
    assert "##INFO=<ID=HAPLOTYPE" in text
    assert "##INFO=<ID=PHASE_BLOCK" in text
    assert "GT:GQ:VAF:AD:HP:PQ" in text
    assert "HAPLOTYPE=H1" in text
    assert "PHASE_BLOCK=chr1:4000-6000" in text
    assert ":H1:42.00" in text
