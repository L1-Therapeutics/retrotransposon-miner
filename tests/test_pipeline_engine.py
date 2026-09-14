import pytest
import tempfile
from pathlib import Path
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
            "right_clips": ["ATTGCAGCTAGCTAGAAAAA"],
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
