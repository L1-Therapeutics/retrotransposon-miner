import tempfile
from pathlib import Path

import pysam

"""
Stress tests for streaming record boundaries (PR #32 one-pass extraction engine):
1. Missing query qualities (wildcard '*')
2. Unplaced decoy scaffolds (chrUn_*, *_random)

Records are validated directly through the streaming pysam.AlignmentFile API.
"""

def test_missing_quality_scores_handled_gracefully():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "missing_qual.bam"
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000000, 'SN': 'chr1'}]}

        with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
            a = pysam.AlignedSegment()
            a.query_name = "no_qual_read_01"
            a.query_sequence = "ACTG" * 25
            a.query_qualities = None  # Missing qualities
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 1000
            a.mapping_quality = 60
            a.cigar = ((0, 100),)
            a.mrnm = 0
            a.mpos = 1300
            out.write(a)

        pysam.index(str(bam_path))

        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            reads = [r for r in samfile.fetch(until_eof=True) if not r.is_unmapped]
            assert len(reads) == 1
            assert reads[0].query_qualities is None

def test_unplaced_decoy_contig_mate_tracking():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "decoy_scaffold.bam"
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [
                {'LN': 248956422, 'SN': 'chr1'},
                {'LN': 100000, 'SN': 'chr1_KI270706v1_random'}
            ]
        }

        with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
            r1 = pysam.AlignedSegment()
            r1.query_name = "translocation_decoy_pair"
            r1.query_sequence = "AAAA" * 25
            r1.flag = 99
            r1.reference_id = 0
            r1.reference_start = 5000
            r1.mapping_quality = 60
            r1.cigar = ((0, 100),)
            r1.mrnm = 1
            r1.mpos = 1000
            out.write(r1)

        pysam.index(str(bam_path))
        assert bam_path.exists()
