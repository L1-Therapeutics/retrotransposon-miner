import tempfile
from pathlib import Path

from benchmarks.benchmark_one_pass_extraction import profile_indexed_fetch, profile_sequential
from scripts.generate_synthetic_bam import generate_synthetic_bam


def test_benchmark_metrics_collection():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "test_bench.bam"
        generate_synthetic_bam(bam_path, num_reads=100)

        metrics = profile_sequential(bam_path)

        assert metrics["read_count"] > 0
        assert metrics["runtime_seconds"] >= 0.0
        assert metrics["peak_memory_mb"] > 0.0


def test_indexed_fetch_benchmark_metrics_collection():
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "test_bench.bam"
        generate_synthetic_bam(bam_path, num_reads=200)
        import pysam
        pysam.index(str(bam_path))

        metrics = profile_indexed_fetch(bam_path, n_loci=5)

        assert metrics["loci_fetched"] == 5
        assert metrics["runtime_seconds"] >= 0.0
        assert metrics["peak_memory_mb"] > 0.0
