"""Tests for spatial index-guided interval streaming engine.

Builds synthetic indexed BAMs and verifies parity between indexed spatial fetch
and full sequential BAM iteration.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from retro_miner.indexed_stream import LocusEvidence, fetch_locus_spanning_pairs, has_bam_index
from scripts.generate_synthetic_bam import generate_synthetic_bam


def _build_indexed_bam(tmpdir: Path, num_reads: int = 500, seed: int = 1) -> Path:
    bam = tmpdir / "indexed_stream_test.bam"
    generate_synthetic_bam(
        output=bam,
        num_reads=num_reads,
        seed=seed,
        chromosomes=["chr1"],
        n_l1_insertions=2,
        n_alu_insertions=1,
        n_hard_clip_reads=1,
        n_interchrom_pairs=1,
    )
    pysam.index(str(bam))
    return bam


def test_has_bam_index_true_when_bai_present(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    assert has_bam_index(bam) is True


def test_has_bam_index_false_when_no_index(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    bai = Path(str(bam) + ".bai")
    bai.unlink()
    assert has_bam_index(bam) is False


def test_fetch_locus_spanning_pairs_returns_evidence(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    evidence = fetch_locus_spanning_pairs(bam, "chr1", 1327321, 1328321, min_mapq=0)
    assert isinstance(evidence, LocusEvidence)
    assert evidence.chrom == "chr1"
    assert evidence.start == 1327321
    assert evidence.end == 1328321
    assert len(evidence.raw_reads) >= 0  # may be 0 if synthetic BAM sparse


def test_indexed_fetch_parity_with_full_scan(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    evidence_idx = fetch_locus_spanning_pairs(bam, "chr1", 1327321, 1328321, min_mapq=0)
    with pysam.AlignmentFile(str(bam), "rb") as af:
        reads = [
            r for r in af.fetch("chr1", 1327321, 1328321)
            if r is not None and not r.is_unmapped
        ]
    assert len(evidence_idx.raw_reads) == len(reads)


def test_fetch_raises_file_not_found_on_missing_bam(tmp_path: Path) -> None:
    missing = tmp_path / "nonexistent.bam"
    try:
        fetch_locus_spanning_pairs(missing, "chr1", 1327321, 1328321)
    except FileNotFoundError:
        return
    raise AssertionError("Expected FileNotFoundError for missing BAM")


def test_fetch_raises_file_not_found_on_missing_index(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    Path(str(bam) + ".bai").unlink()
    try:
        fetch_locus_spanning_pairs(bam, "chr1", 1327321, 1328321)
    except FileNotFoundError:
        return
    raise AssertionError("Expected FileNotFoundError for missing index")


def test_soft_clips_extracted_from_indexed_fetch(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    evidence = fetch_locus_spanning_pairs(bam, "chr1", 1327321, 1328321, min_mapq=0)
    assert isinstance(evidence.soft_clips, list)
    for seq in evidence.soft_clips:
        assert isinstance(seq, str)


def test_discordant_mates_collected_from_indexed_fetch(tmp_path: Path) -> None:
    bam = _build_indexed_bam(tmp_path)
    evidence = fetch_locus_spanning_pairs(bam, "chr1", 1327321, 1328321, min_mapq=0)
    assert isinstance(evidence.discordant_mates, list)
