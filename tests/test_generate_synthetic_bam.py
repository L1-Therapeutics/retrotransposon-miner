"""Tests for scripts/generate_synthetic_bam.py.

Validates the output structure of the synthetic BAM generator:

  - the BAM and its auto-generated ``.bai`` index exist,
  - total read count matches ``num_reads`` exactly,
  - the header declares the requested hg38 contigs,
  - alignment edge cases are present (secondary 256, supplementary 2048,
    unmapped mate 8, MAPQ 0),
  - the file is coordinate-sorted (indexable),
  - output is reproducible for a fixed seed, and
  - the argparse CLI end-to-end writes a valid indexed BAM.

The script is loaded via importlib so running it does not create a package
dependency for the test suite.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "generate_synthetic_bam.py"


def _load_generate_synthetic_bam():
    spec = importlib.util.spec_from_file_location("generate_synthetic_bam", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def gen():
    return _load_generate_synthetic_bam()


def _all_reads(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        return list(sf.fetch(until_eof=True))


def test_bam_and_index_created(tmp_path, gen) -> None:
    bam, bai = gen.generate_synthetic_bam(tmp_path / "syn.bam", num_reads=200, seed=7)
    assert bam.exists() and bam.stat().st_size > 0
    assert bai.exists() and bai.stat().st_size > 0


def test_total_read_count_matches_requested(tmp_path, gen) -> None:
    num_reads = 237
    bam, _ = gen.generate_synthetic_bam(tmp_path / "count.bam", num_reads=num_reads, seed=3)
    assert len(_all_reads(bam)) == num_reads


def test_header_declares_hg38_contigs(tmp_path, gen) -> None:
    chroms = ["chr1", "chr22", "chrX", "chrM"]
    bam, _ = gen.generate_synthetic_bam(
        tmp_path / "contigs.bam", num_reads=150, seed=11, chromosomes=chroms
    )
    with pysam.AlignmentFile(bam, "rb") as sf:
        assert [c["SN"] for c in sf.header["SQ"]] == chroms
        assert all(sf.get_reference_length(c) == gen.HG38_CHROM_LENGTHS[c] for c in chroms)
    assert gen.HG38_CHROM_LENGTHS["chrM"] == 16569


def test_edge_case_flags_present(tmp_path, gen) -> None:
    bam, _ = gen.generate_synthetic_bam(
        tmp_path / "edge.bam",
        num_reads=500,
        seed=1,
        n_l1_insertions=1,
        n_alu_insertions=1,
    )
    reads = _all_reads(bam)
    assert any(r.flag & 256 for r in reads), "expected at least one secondary alignment (256)"
    assert any(r.flag & 2048 for r in reads), "expected at least one supplementary alignment (2048)"
    assert any(r.flag & 8 for r in reads), "expected at least one unmapped mate (flag 8)"
    mapped = [r for r in reads if not r.is_unmapped]
    assert any(r.mapping_quality == 0 for r in mapped), "expected at least one MAPQ 0 read"
    assert any(r.is_secondary for r in reads)
    assert any(r.is_supplementary for r in reads)


def test_reads_are_coordinate_sorted(tmp_path, gen) -> None:
    bam, _ = gen.generate_synthetic_bam(tmp_path / "sorted.bam", num_reads=400, seed=13)
    prev: tuple[int, int] | None = None
    for r in _all_reads(bam):
        key = (r.reference_id, r.reference_start)
        if prev is not None:
            assert key >= prev, "reads out of coordinate order"
        prev = key


def test_deterministic_given_same_seed(tmp_path, gen) -> None:
    a = tmp_path / "a.bam"
    b = tmp_path / "b.bam"
    gen.generate_synthetic_bam(a, num_reads=300, seed=99)
    gen.generate_synthetic_bam(b, num_reads=300, seed=99)
    assert a.read_bytes() == b.read_bytes()


def test_different_seed_changes_reads(tmp_path, gen) -> None:
    a = tmp_path / "a.bam"
    b = tmp_path / "b.bam"
    gen.generate_synthetic_bam(a, num_reads=300, seed=99)
    gen.generate_synthetic_bam(b, num_reads=300, seed=101)
    assert a.read_bytes() != b.read_bytes()


def test_invalid_chromosome_rejected(tmp_path, gen) -> None:
    with pytest.raises(ValueError, match="not an hg38 contig"):
        gen.generate_synthetic_bam(
            tmp_path / "bad.bam", num_reads=200, seed=5, chromosomes=["chr22", "chrZZ"]
        )


def test_too_few_reads_rejected(tmp_path, gen) -> None:
    with pytest.raises(ValueError, match="below minimum"):
        gen.generate_synthetic_bam(tmp_path / "tiny.bam", num_reads=3, seed=5)


def test_cli_writes_indexed_bam(tmp_path) -> None:
    out = tmp_path / "cli.bam"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--output",
            str(out),
            "--num-reads",
            "150",
            "--seed",
            "5",
            "--chromosomes",
            "chr22,chrX",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert out.exists() and out.stat().st_size > 0
    assert Path(str(out) + ".bai").exists()
    assert "150 reads" in result.stdout


def test_cli_requires_output(tmp_path) -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--num-reads", "50", "--seed", "1"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "--output" in result.stderr
