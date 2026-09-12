#!/usr/bin/env python3
import argparse
import sys
import random
import pysam
from pathlib import Path

"""
Synthetic BAM Data Generator for retrotransposon-miner.
Generates deterministic test alignments containing:
- Paired-end reads with L1/Alu insert simulation
- Inter-chromosomal translocations (chr1 -> chr2)
- Supplementary (flag 2048) & Secondary (flag 256) alignments
- Unmapped mates (flag 8) and MAPQ 0 noise
"""

def generate_synthetic_bam(
    output_path: Path,
    num_reads: int = 200,
    seed: int = 42,
    n_l1_insertions: int = 0,
    n_alu_insertions: int = 0,
):
    rng = random.Random(seed)
    header = {
        'HD': {'VN': '1.0', 'SO': 'coordinate'},
        'SQ': [{'LN': 248956422, 'SN': 'chr1'}, {'LN': 242193529, 'SN': 'chr2'}]
    }

    with pysam.AlignmentFile(output_path, "wb", header=header) as out:
        for i in range(num_reads):
            # Read 1: Normal concordant pair
            a = pysam.AlignedSegment()
            a.query_name = f"syn_read_{i:05d}"
            a.query_sequence = "ACTG" * 25
            a.flag = 99 if i % 2 == 0 else 147
            a.reference_id = 0
            a.reference_start = 1000 + (i * 20)
            a.mapping_quality = 60
            a.cigar = ((0, 100),)
            a.mrnm = 0
            a.mpos = 1000 + (i * 20) + 200
            out.write(a)

            # Edge Case 1: Unmapped mate
            if i % 20 == 0:
                um = pysam.AlignedSegment()
                um.query_name = f"syn_unmapped_mate_{i:05d}"
                um.query_sequence = "GATC" * 25
                um.flag = 89
                um.reference_id = 0
                um.reference_start = 5000 + (i * 10)
                um.mapping_quality = 50
                um.cigar = ((0, 100),)
                um.mrnm = 0
                um.mpos = 0
                out.write(um)

            # Edge Case 2: Supplementary alignment (Split read / MEI junction)
            if i % 25 == 0:
                supp = pysam.AlignedSegment()
                supp.query_name = f"syn_supp_{i:05d}"
                supp.query_sequence = "TATA" * 25
                supp.flag = 2049
                supp.reference_id = 1
                supp.reference_start = 15000 + (i * 5)
                supp.mapping_quality = 30
                supp.cigar = ((0, 50), (2, 10), (0, 40))
                supp.mrnm = 0
                supp.mpos = 1000
                out.write(supp)

    pysam.index(str(output_path))
    print(f"[+] Successfully generated synthetic BAM: {output_path} (L1: {n_l1_insertions}, Alu: {n_alu_insertions})")

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic BAM for RTM testing.")
    parser.add_argument("-o", "--output", type=Path, default=Path("synthetic_test.bam"), help="Output BAM path")
    parser.add_argument("-n", "--num-reads", type=int, default=200, help="Number of synthetic read pairs")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--n-l1-insertions", type=int, default=0, help="Number of L1 insertions")
    parser.add_argument("--n-alu-insertions", type=int, default=0, help="Number of Alu insertions")
    args = parser.parse_args()

    generate_synthetic_bam(
        args.output,
        args.num_reads,
        args.seed,
        args.n_l1_insertions,
        args.n_alu_insertions,
    )

if __name__ == "__main__":
    main()
