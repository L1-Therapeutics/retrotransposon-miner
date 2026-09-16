import tempfile
from pathlib import Path

import pysam

from retro_miner.pipeline_engine import (
    process_candidate_locus,
    run_scientific_mei_pipeline,
)


def test_process_candidate_locus_refinement():
    left_clips = ["CGATCGATTGCAGCTAGCTAG"]
    right_clips = ["ATTGCAGCTAGCTAGAAAAAAAAAAAAAAAA"]

    res = process_candidate_locus(
        chrom="chr1",
        pos=50000,
        left_clips=left_clips,
        right_clips=right_clips,
        k_alt=15,
        k_ref=15,
        family_hint="L1",
    )

    assert res["chrom"] == "chr1"
    assert res["pos"] == 50000
    assert res["genotype"] == "0/1"
    assert res["vaf"] == 0.50
    assert res["genotype_quality"] >= 50.0
    assert res["poly_a_detected"] is True
    assert res["tsd_seq"] == "ATTGCAGCTAGCTAG"


def test_run_scientific_mei_pipeline_e2e_vcf():
    raw_loci = [
        {
            "chrom": "chr1",
            "pos": 100000,
            "family": "L1",
            "left_clips": ["ATTGCAGCTAGCTAG"],
            "right_clips": ["ATTGCAGCTAGCTAGAAAAAAAAAAAAAAAA"],
            "k_alt": 20,
            "k_ref": 2,
        }
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_out = Path(tmpdir) / "output_refined.vcf"
        out_path = run_scientific_mei_pipeline(raw_loci, vcf_out, sample_name="NA12878")

        assert out_path.exists()
        content = out_path.read_text()

        # Verify VCF header & variant lines
        assert "##fileformat=VCFv4.3" in content
        assert "chr1\t100000\tMEI_chr1_100000" in content
        assert "POLYA=1" in content
        assert "1/1:" in content or "0/1:" in content


def test_process_candidate_locus_phases_against_het_snps(tmp_path):
    chrom = "chr1"
    breakpoint = 5000
    snp_pos = 5100
    snp_ref = "G"
    snp_alt = "A"

    def _mei_read(i):
        seq = "T" * 20 + "A" * 300 + "T" * 130
        bases = list(seq)
        bases[420] = snp_alt
        read = pysam.AlignedSegment()
        read.query_name = f"mei_read_{i}"
        read.query_sequence = "".join(bases)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 4980
        read.mapping_quality = 60
        read.cigar = ((0, 20), (1, 300), (0, 130))
        return read

    def _ref_read(i):
        seq = "T" * 150
        bases = list(seq)
        bases[120] = snp_ref
        read = pysam.AlignedSegment()
        read.query_name = f"ref_read_{i}"
        read.query_sequence = "".join(bases)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 4980
        read.mapping_quality = 60
        read.cigar = ((0, 150),)
        return read

    bam_path = tmp_path / "phased.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 20000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for read in [_mei_read(i) for i in range(10)] + [_ref_read(i) for i in range(10)]:
            out.write(read)
    sorted_path = bam_path.with_name(bam_path.stem + ".sorted.bam")
    pysam.sort("-o", str(sorted_path), str(bam_path))
    pysam.index(str(sorted_path))

    res = process_candidate_locus(
        chrom=chrom,
        pos=breakpoint,
        left_clips=["ATTGCAGCTAGCTAG"],
        right_clips=["ATTGCAGCTAGCTAGAAAAA"],
        k_alt=15,
        k_ref=15,
        family_hint="L1",
        bam_path=sorted_path,
        locus_start=breakpoint - 20,
        locus_end=breakpoint + 20,
        het_snps=[(snp_pos, snp_ref, snp_alt)],
        phase_window_bp=1000,
    )

    assert res["phase_linkage"].haplotype_assigned == "H1"
    assert res["phase_linkage"].phase_confidence >= 13.0
    assert res["phase_linkage"].phase_block_id == f"{chrom}:4000-6000"
