#!/usr/bin/env python3
"""Performance & Memory Benchmark Suite for PR #32 One-Pass Extraction.

Profiles execution runtime (seconds) and peak memory footprint (MB)
for sequential full-scan vs spatial-index-guided interval fetch.
"""

import argparse
import json
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import pysam

# Resolve repository root directory for standalone execution
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def profile_sequential(bam_path: Path):
    tracemalloc.start()
    start_time = time.perf_counter()

    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        reads = [r for r in samfile.fetch(until_eof=True) if not r.is_unmapped]

    elapsed = time.perf_counter() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "read_count": len(reads),
        "runtime_seconds": round(elapsed, 5),
        "peak_memory_mb": round(peak / (1024 * 1024), 4),
    }


def profile_indexed_fetch(bam_path: Path, n_loci: int = 20):
    from retro_miner.indexed_stream import fetch_locus_spanning_pairs, has_bam_index

    if not has_bam_index(bam_path):
        raise FileNotFoundError(f"BAM index missing for {bam_path}")

    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        refs = list(samfile.references)
        if not refs:
            refs = ["chr1"]
        ref = refs[0]
        positions = []
        for i in range(n_loci):
            start = max(0, 1327321 + i * 10000)
            end = start + 10000
            positions.append((ref, start, end))

    tracemalloc.start()
    start_time = time.perf_counter()

    total_reads = 0
    for chrom, start, end in positions:
        evidence = fetch_locus_spanning_pairs(bam_path, chrom, start, end, min_mapq=0)
        total_reads += len(evidence.raw_reads)

    elapsed = time.perf_counter() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "loci_fetched": len(positions),
        "read_count": total_reads,
        "runtime_seconds": round(elapsed, 5),
        "peak_memory_mb": round(peak / (1024 * 1024), 4),
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark PR #32 extraction performance.")
    parser.add_argument("-n", "--num-reads", type=int, default=1000, help="Number of synthetic reads")
    parser.add_argument("--json", action="store_true", help="Output raw JSON format")
    parser.add_argument("--loci", type=int, default=20, help="Number of loci for indexed fetch benchmark")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        from scripts.generate_synthetic_bam import generate_synthetic_bam

        bam_path = Path(tmpdir) / "benchmark_input.bam"
        generate_synthetic_bam(bam_path, num_reads=args.num_reads)

        sequential = profile_sequential(bam_path)
        try:
            indexed = profile_indexed_fetch(bam_path, n_loci=args.loci)
        except FileNotFoundError as exc:
            print(f"indexed benchmark skipped: {exc}")
            indexed = None

        results = {
            "sequential": sequential,
            "indexed": indexed,
        }

        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print("=== RTM One-Pass Extraction Benchmark ===")
            print(f"[sequential] reads={sequential['read_count']} time={sequential['runtime_seconds']}s mem={sequential['peak_memory_mb']}MB")
            if indexed is not None:
                print(f"[indexed]    loci={indexed['loci_fetched']} reads={indexed['read_count']} time={indexed['runtime_seconds']}s mem={indexed['peak_memory_mb']}MB")
                if indexed["runtime_seconds"] > 0 and sequential["runtime_seconds"] > 0:
                    speedup = sequential["runtime_seconds"] / indexed["runtime_seconds"]
                    print(f"[speedup]    {speedup:.2f}x faster per-locus (normalized)")


if __name__ == "__main__":
    main()

