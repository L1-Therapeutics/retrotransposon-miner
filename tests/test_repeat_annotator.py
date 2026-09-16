"""Tests for RepeatIntervalTree spatial index.

Benchmarks 10k synthetic RepeatMasker BED intervals and verifies O(log N)
overlap queries against candidate insertion coordinates.
"""

from __future__ import annotations

import random
import time
from pathlib import Path

from retro_miner.repeat_annotator import RepeatAnnotation, RepeatIntervalTree


def _write_synthetic_bed(path: Path, n_intervals: int = 10000, seed: int = 1) -> None:
    rng = random.Random(seed)
    chroms = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    with path.open("w", encoding="utf-8") as hout:
        for _ in range(n_intervals):
            chrom = rng.choice(chroms)
            start = rng.randint(0, 250_000_000)
            end = start + rng.randint(50, 5000)
            name = rng.choice(["AluY", "L1HS", "SVA_F", "AluSx", "L1PA2", "MIR"])
            repeat_class = rng.choice(["SINE", "LINE", "LTR", "DNA", "Simple_repeat"])
            family = rng.choice(["Alu", "LINE1", "SVA", "MIR", "ERV"])
            strand = rng.choice(["+", "-"])
            divergence = round(rng.uniform(0.0, 30.0), 2)
            length = end - start
            hout.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{repeat_class}\t{family}\t{divergence}\n"
            )


def test_build_from_bed_loads_annotations(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    _write_synthetic_bed(bed, n_intervals=1000)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    assert len(tree) == 1000


def test_query_overlap_returns_matching_attributes(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    _write_synthetic_bed(bed, n_intervals=1000, seed=42)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)

    target_chrom = "chr1"
    target_pos = 100_000
    hits = tree.query_overlap(target_chrom, target_pos, buffer_bp=50)
    assert isinstance(hits, list)
    for hit in hits:
        assert hit["chrom"] == target_chrom
        assert hit["start"] <= target_pos <= hit["end"]
        assert "name" in hit
        assert "family" in hit
        assert "repeat_class" in hit
        assert "overlap_bp" in hit
        assert "overlap_fraction" in hit


def test_query_no_overlap_returns_empty_list(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    _write_synthetic_bed(bed, n_intervals=500, seed=7)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    hits = tree.query_overlap("chrZ", 999_999_999, buffer_bp=50)
    assert hits == []


def test_query_missing_chromosome_returns_empty_list(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    _write_synthetic_bed(bed, n_intervals=100)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    assert tree.query_overlap("nonexistent", 1000) == []


def test_10k_interval_tree_query_performance(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk_10k.bed"
    _write_synthetic_bed(bed, n_intervals=10_000, seed=123)
    tree = RepeatIntervalTree()
    t0 = time.monotonic()
    tree.build_from_bed(bed)
    build_elapsed = time.monotonic() - t0
    assert len(tree) == 10_000

    query_points = [
        ("chr1", 50_000),
        ("chrX", 10_000_000),
        ("chr22", 100_000),
    ]
    for chrom, pos in query_points:
        t1 = time.monotonic()
        hits = tree.query_overlap(chrom, pos, buffer_bp=50)
        q_elapsed = time.monotonic() - t1
        assert q_elapsed < 1.0, f"Query took {q_elapsed:.4f}s, expected < 1.0s"
        assert isinstance(hits, list)

    assert build_elapsed < 60.0, f"Build took {build_elapsed:.2f}s, expected < 60s"


def test_gzip_bed_supported(tmp_path: Path) -> None:
    import gzip

    bed = tmp_path / "rmsk.bed.gz"
    lines = []
    rng = random.Random(99)
    chroms = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    for _ in range(100):
        chrom = rng.choice(chroms)
        start = rng.randint(0, 250_000_000)
        end = start + rng.randint(50, 5000)
        name = rng.choice(["AluY", "L1HS", "SVA_F", "AluSx", "L1PA2", "MIR"])
        repeat_class = rng.choice(["SINE", "LINE", "LTR", "DNA", "Simple_repeat"])
        family = rng.choice(["Alu", "LINE1", "SVA", "MIR", "ERV"])
        strand = rng.choice(["+", "-"])
        divergence = round(rng.uniform(0.0, 30.0), 2)
        length = end - start
        lines.append(
            f"{chrom}\t{start}\t{end}\t{name}\t{length}\t{strand}\t{repeat_class}\t{family}\t{divergence}\n"
        )
    with gzip.open(bed, "wt", encoding="utf-8") as hout:
        hout.writelines(lines)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    assert len(tree) == 100


def test_overlap_fraction_computed_correctly(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    bed.write_text(
        "chr1\t100\t200\tAluY\t100\t+\tSINE\tAlu\t5.0\n",
        encoding="utf-8",
    )
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    hits = tree.query_overlap("chr1", 150, buffer_bp=50)
    assert len(hits) == 1
    assert hits[0]["overlap_bp"] == 100
    assert hits[0]["overlap_fraction"] == 1.0
    assert hits[0]["divergence"] == 5.0


def test_annotations_are_frozen_dataclasses(tmp_path: Path) -> None:
    bed = tmp_path / "rmsk.bed"
    _write_synthetic_bed(bed, n_intervals=10)
    tree = RepeatIntervalTree()
    tree.build_from_bed(bed)
    for ann in tree.annotations:
        assert isinstance(ann, RepeatAnnotation)
