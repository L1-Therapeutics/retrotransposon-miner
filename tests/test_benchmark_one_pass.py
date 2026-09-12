import pytest
import tempfile
from pathlib import Path
from benchmarks.benchmark_one_pass_extraction import profile_extraction
from scripts.generate_synthetic_bam import generate_synthetic_bam

def test_benchmark_metrics_collection():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "test_bench.bam"
        generate_synthetic_bam(bam_path, num_reads=100)
        
        metrics = profile_extraction(bam_path)
        
        assert metrics["read_count"] > 0
        assert metrics["runtime_seconds"] >= 0.0
        assert metrics["peak_memory_mb"] > 0.0
