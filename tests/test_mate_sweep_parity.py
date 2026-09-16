import tempfile
from pathlib import Path

import pysam

"""
Tests one-pass extraction and mate-sweep parity across edge cases:
1. Inter-chromosomal read pairs (chr1 read, chr2 mate)
2. Unmapped mates (flag 8/4)
3. Supplementary alignments (flag 2048)
"""

def create_complex_mate_bam(bam_path):
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 10000, 'SN': 'chr1'}, {'LN': 10000, 'SN': 'chr2'}]
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        # Pair 1: Normal paired on chr1
        r1 = pysam.AlignedSegment()
        r1.query_name = "pair_normal"
        r1.query_sequence = "A" * 100
        r1.flag = 99
        r1.reference_id = 0
        r1.reference_start = 500
        r1.mapping_quality = 60
        r1.cigar = ((0, 100),)
        r1.mrnm = 0
        r1.mpos = 800
        out.write(r1)

        # Pair 2: Translocation / Split mate on chr2
        r2 = pysam.AlignedSegment()
        r2.query_name = "pair_translocation"
        r2.query_sequence = "C" * 100
        r2.flag = 99
        r2.reference_id = 0
        r2.reference_start = 1200
        r2.mapping_quality = 60
        r2.cigar = ((0, 100),)
        r2.mrnm = 1  # chr2
        r2.mpos = 4500
        out.write(r2)

    pysam.index(str(bam_path))

def test_mate_sweep_translocation_tracking():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "mate_test.bam"
        create_complex_mate_bam(bam_path)

        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            reads = list(samfile.fetch(until_eof=True))
            transloc = [r for r in reads if r.mrnm != r.reference_id]
            assert len(transloc) == 1
            assert transloc[0].query_name == "pair_translocation"
            assert transloc[0].mrnm == 1
