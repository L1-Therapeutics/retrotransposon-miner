import pytest
import pysam
import tempfile
from pathlib import Path
from scripts.generate_synthetic_bam import generate_synthetic_bam
from scripts.compare_callset_parity import compare_parity
from benchmarks.benchmark_one_pass_extraction import profile_extraction

"""
Integration test verifying end-to-end pipeline components:
1. Synthetic BAM generation
2. One-pass extraction profiling
3. Callset parity validation across VCFs
"""

def test_end_to_end_synthetic_parity_pipeline():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "pipeline_input.bam"
        vcf_a = Path(tmpdir) / "calls_baseline.vcf"
        vcf_b = Path(tmpdir) / "calls_candidate.vcf"

        # 1. Generate synthetic alignment dataset
        generate_synthetic_bam(bam_path, num_reads=300)
        assert bam_path.exists()

        # 2. Run extraction performance profiling
        metrics = profile_extraction(bam_path)
        assert metrics["read_count"] > 0
        assert metrics["peak_memory_mb"] > 0.0

        # 3. Verify VCF parity comparator on mock outputs
        vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        vcf_record = "chr1\t10500\t.\tA\t<INS:MEI:L1HS>\t60\tPASS\t.\n"

        with open(vcf_a, "w") as f:
            f.write(vcf_header + vcf_record)
        with open(vcf_b, "w") as f:
            f.write(vcf_header + vcf_record)

        parity = compare_parity(vcf_a, vcf_b)
        assert parity["exact_match"] is True
        assert parity["concordance_pct"] == 100.0
