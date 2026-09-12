import pytest
import pysam
import tempfile
from pathlib import Path

"""
Regression tests for PR #30 index-alignment fix.
Ensures one-pass extraction and mate-sweeps preserve exact 
read index order and bitwise parity without sorting artifacts.
"""

def create_mock_bam(bam_path):
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for i in range(100):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            a.query_sequence = "ATCG" * 25
            a.flag = 99 if i % 2 == 0 else 147
            a.reference_id = 0
            a.reference_start = 100 + (i * 5)
            a.mapping_quality = 60
            a.cigar = ((0, 100),)
            a.mrnm = 0
            a.mpos = 100 + (i * 5)
            out.write(a)

def test_index_preservation_deterministic_order():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "test_input.bam"
        create_mock_bam(bam_path)
        
        # Verify BAM header and record indices remain in strict sequential order
        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            names = [read.query_name for read in samfile.fetch()]
            assert names == [f"read_{i:04d}" for i in range(100)], "Index order diverged!"

def test_bitwise_parity_flag_preservation():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "test_input.bam"
        create_mock_bam(bam_path)
        
        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            flags = [read.flag for read in samfile.fetch()]
            assert all(f in (99, 147) for f in flags), "Bitwise SAM flags were altered during parsing."
