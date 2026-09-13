import pytest
import tempfile
from pathlib import Path
from retro_miner.vcf_export import write_mei_vcf

def test_vcf_export_header_and_body_formatting():
    candidate_records = [
        {
            "chrom": "chr1",
            "pos": 15000,
            "ref_base": "A",
            "family": "L1",
            "subfamily": "L1HS",
            "genotype": "0/1",
            "genotype_quality": 75.0,
            "vaf": 0.48,
            "k_ref": 12,
            "k_alt": 11,
            "tsd_seq": "ATTGCAG",
            "tsd_len": 7,
            "poly_a_detected": True,
            "mei_llr": 4.25,
        }
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "test_output.vcf"
        write_mei_vcf(candidate_records, vcf_file, sample_name="NA12878")

        assert vcf_file.exists()
        content = vcf_file.read_text()

        # Check VCF specification headers
        assert "##fileformat=VCFv4.3" in content
        assert "##ALT=<ID=INS:MEI:L1HS" in content
        assert "#CHROM\tPOS\tID\tREF\tALT" in content
        assert "NA12878" in content

        # Check variant line contents
        assert "chr1\t15000\tMEI_chr1_15000\tA\t<INS:MEI:L1HS>" in content
        assert "PASS" in content
        assert "TSD=ATTGCAG" in content
        assert "POLYA=1" in content
        assert "0/1:75.0:0.48:12,11" in content

def test_low_quality_filter_assignment():
    candidate_records = [
        {
            "chrom": "chr2",
            "pos": 20000,
            "genotype": "0/1",
            "genotype_quality": 12.5,  # Below PASS threshold
            "vaf": 0.15,
            "k_ref": 20,
            "k_alt": 3,
        }
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "low_qual.vcf"
        write_mei_vcf(candidate_records, vcf_file)

        content = vcf_file.read_text()
        assert "LowQual" in content
