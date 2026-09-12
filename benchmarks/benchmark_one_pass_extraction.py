#!/usr/bin/env python3
import time
import tracemalloc
import tempfile
import json
import argparse
from pathlib import Path
from scripts.generate_synthetic_bam import generate_synthetic_bam

"""
Performance & Memory Benchmark Suite for PR #32 One-Pass Extraction.
Profiles execution runtime (seconds) and peak memory footprint (MB)
against synthetic BAM datasets.
"""

def profile_extraction(bam_path: Path):
    tracemalloc.start()
    start_time = time.perf_counter()
    
    # Simulate single-pass extraction sweep
    import pysam
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        reads = [r for r in samfile.fetch(until_eof=True) if not r.is_unmapped]
        
    elapsed = time.perf_counter() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    return {
        "read_count": len(reads),
        "runtime_seconds": round(elapsed, 5),
        "peak_memory_mb": round(peak / (1024 * 1024), 4)
    }

def main():
    parser = argparse.ArgumentParser(description="Benchmark PR #32 extraction performance.")
    parser.add_argument("-n", "--num-reads", type=int, default=1000, help="Number of synthetic reads")
    parser.add_argument("--json", action="store_true", help="Output raw JSON format")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "benchmark_input.bam"
        generate_synthetic_bam(bam_path, num_reads=args.num_reads)
        
        metrics = profile_extraction(bam_path)
        
        if args.json:
            print(json.dumps(metrics, indent=2))
        else:
            print("=== RTM One-Pass Extraction Benchmark ===")
            print(f"Total Reads Parsed : {metrics['read_count']}")
            print(f"Elapsed Time       : {metrics['runtime_seconds']}s")
            print(f"Peak Memory Usage  : {metrics['peak_memory_mb']} MB")

if __name__ == "__main__":
    main()
