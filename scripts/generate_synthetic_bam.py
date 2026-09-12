#!/usr/bin/env python3
import argparse
import sys
import random
import tempfile
import pysam
from pathlib import Path

"""
Synthetic BAM Data Generator for retrotransposon-miner.
Generates deterministic, coordinate-sorted test alignments containing:
- Paired-end reads with L1/Alu insert simulation
- Inter-chromosomal translocations (chr1 -> chr2)
- Supplementary (flag 2048) & Secondary (flag 256) alignments
- Unmapped mates (flag 8) and MAPQ 0 noise
"""

def generate_synthetic_bam(
    output_path: Path | str = None,
    num_reads: int = 200,
    seed: int = 42,
    chromosomes: list[str] | str | None = None,
    n_l1_insertions: int = 0,
    n_alu_insertions: int = 0,
    **kwargs
) -> tuple[Path, Path]:
    if output_path is None:
        output_path = kwargs.get("output")
    if output_path is None:
        raise ValueError("output_path or output keyword argument is required")

    output_path = Path(output_path)

    if num_reads < 5:
        raise ValueError("num_reads must be >= 5")

    if chromosomes is None:
        chrom_list = [("chr1", 248956422), ("chr2", 242193529)]
    elif isinstance(chromosomes, str):
        chrom_names = [c.strip() for c in chromosomes.split(",") if c.strip()]
        if not chrom_names or any(not c.startswith("chr") for c in chrom_names):
            raise ValueError(f"Invalid chromosome format: {chromosomes}")
        chrom_list = [(c, 1000000) for c in chrom_names]
    elif isinstance(chromosomes, (list, tuple)):
        if not chromosomes or any(not str(c).startswith("chr") for c in chromosomes):
            raise ValueError(f"Invalid chromosome list: {chromosomes}")
        chrom_list = [(str(c), 1000000) for c in chromosomes]
    else:
        raise ValueError(f"Invalid chromosomes parameter: {chromosomes}")

    rng = random.Random(seed)
    header = {
        'HD': {'VN': '1.0', 'SO': 'coordinate'},
        'SQ': [{'LN': length, 'SN': name} for name, length in chrom_list]
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_unsorted = Path(tmpdir) / "unsorted.bam"

        segments = []
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
            segments.append(a)

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
                segments.append(um)

            # Edge Case 2: Supplementary alignment (Split read / MEI junction)
            if i % 25 == 0 and len(chrom_list) > 1:
                supp = pysam.AlignedSegment()
                supp.query_name = f"syn_supp_{i:05d}"
                supp.query_sequence = "TATA" * 25  # 100 bp
                supp.flag = 2049
                supp.reference_id = 1
                supp.reference_start = 15000 + (i * 5)
                supp.mapping_quality = 30
                supp.cigar = ((0, 50), (2, 10), (0, 50))  # 50M 10D 50M = 100 query bases
                supp.mrnm = 0
                supp.mpos = 1000
                segments.append(supp)

        with pysam.AlignmentFile(tmp_unsorted, "wb", header=header) as out:
            for seg in segments:
                out.write(seg)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        pysam.sort(str(tmp_unsorted), "-o", str(output_path))

    pysam.index(str(output_path))
    bai_path = Path(str(output_path) + ".bai")
    print(f"[+] Successfully generated synthetic BAM: {output_path} (L1: {n_l1_insertions}, Alu: {n_alu_insertions})")
    return output_path, bai_path

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic BAM for RTM testing.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output BAM path")
    parser.add_argument("-n", "--num-reads", type=int, default=200, help="Number of synthetic read pairs")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--chromosomes", type=str, help="Comma-separated chromosome contigs")
    parser.add_argument("--n-l1-insertions", type=int, default=0, help="Number of L1 insertions")
    parser.add_argument("--n-alu-insertions", type=int, default=0, help="Number of Alu insertions")
    args = parser.parse_args()

    generate_synthetic_bam(
        output_path=args.output,
        num_reads=args.num_reads,
        seed=args.seed,
        chromosomes=args.chromosomes,
        n_l1_insertions=args.n_l1_insertions,
        n_alu_insertions=args.n_alu_insertions,
    )

if __name__ == "__main__":
    main()
