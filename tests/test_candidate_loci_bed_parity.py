"""Unit tests for vectorized BED output parity in candidate_loci.py."""

from pathlib import Path
import pandas as pd
from retro_miner.candidate_loci import _write_candidate_windows_bed


def test_write_candidate_windows_bed_vectorized_parity(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "chrom": ["chr1", "chr22", "chrX"],
        "window_start": [100, 5000, 123456],
        "window_end": [200, 5500, 123999],
        "row_id": [1, 42, 999],
    })

    out_file = tmp_path / "candidates.bed"
    _write_candidate_windows_bed(df, out_file)

    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    assert lines[0] == "chr1\t99\t200\t1"
    assert lines[1] == "chr22\t4999\t5500\t42"
    assert lines[2] == "chrX\t123455\t123999\t999"


def test_write_candidate_windows_bed_empty(tmp_path: Path) -> None:
    df = pd.DataFrame(columns=["chrom", "window_start", "window_end", "row_id"])
    out_file = tmp_path / "empty.bed"
    _write_candidate_windows_bed(df, out_file)
    assert out_file.read_text(encoding="utf-8") == ""


def test_write_candidate_windows_bed_type_robustness(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "chrom": ["chr1", "chr2", "chr3"],
        "window_start": ["100", 5000.0, "100"],
        "window_end": [200, "5500", 1000],
        "row_id": [1, 42, 999],
    })
    out_file = tmp_path / "robust.bed"
    _write_candidate_windows_bed(df, out_file)

    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    assert lines[0] == "chr1\t99\t200\t1"
    assert lines[1] == "chr2\t4999\t5500\t42"


def test_write_candidate_windows_bed_na_handling(tmp_path: Path) -> None:
    df = pd.DataFrame({
        "chrom": ["chr1", "chr2"],
        "window_start": [100, None],
        "window_end": [200, None],
        "row_id": [1, 2],
    })
    out_file = tmp_path / "na.bed"
    _write_candidate_windows_bed(df, out_file)

    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    assert lines[0] == "chr1\t99\t200\t1"
    assert lines[1] == "chr2\t0\t0\t2"