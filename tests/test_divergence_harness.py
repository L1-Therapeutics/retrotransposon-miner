"""Validation harness for NM-derived divergence vs rmsk milliDiv."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from retro_miner.mei_support import _write_rmsk_mei_bed, panel_divergence_from_nm


def test_panel_divergence_from_nm_normal():
    assert panel_divergence_from_nm(5, 100) == pytest.approx(0.05)
    assert panel_divergence_from_nm(0, 50) == 0.0
    assert panel_divergence_from_nm(1, 1) == 1.0


def test_panel_divergence_from_nm_invalid():
    assert panel_divergence_from_nm(0, 0) is None
    assert panel_divergence_from_nm(-1, 10) is None
    assert panel_divergence_from_nm(5, -1) is None


def test_write_rmsk_mei_bed_milli_div_roundtrip(tmp_path: Path):
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "0\tfoo\t42\tbaz\tqux\tchr22\t100000\t100200\t.\t+\tL1HS\tLINE\tLINE1\n",
        encoding="utf-8",
    )
    out_bed = tmp_path / "out.bed"
    n = _write_rmsk_mei_bed(rmsk, out_bed)
    assert n == 1
    lines = out_bed.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    parts = lines[0].split("\t")
    assert parts[9] == "42"


def test_write_rmsk_mei_bed_fallback_emits_minus_one(tmp_path: Path):
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n",
        encoding="utf-8",
    )
    out_bed = tmp_path / "out.bed"
    n = _write_rmsk_mei_bed(rmsk, out_bed)
    assert n == 1
    lines = out_bed.read_text(encoding="utf-8").splitlines()
    parts = lines[0].split("\t")
    assert parts[9] == "-1"


def _make_synthetic_fixture(tmp_path: Path) -> tuple[Path, Path]:
    calls = tmp_path / "calls.tsv"
    calls.write_text(
        "chrom\twindow_start\twindow_end\tconsensus_mei_family\tnm_mean\taln_len_mean\n"
        "chr22\t100050\t100250\tLINE1\t5\t100\n"
        "chr22\t200100\t200300\tLINE1\t20\t100\n"
        "chr22\t300050\t300250\tALU\t2\t100\n",
        encoding="utf-8",
    )
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "0\tfoo\t50\tbaz\tqux\tchr22\t100000\t100500\t.\t+\tL1HS\tLINE\tLINE1\n"
        "0\tfoo\t100\tbaz\tqux\tchr22\t200000\t200500\t.\t+\tL1HS\tLINE\tLINE1\n"
        "0\tfoo\t10\tbaz\tqux\tchr22\t300000\t300500\t.\t+\tAluSx\tSINE\tAlu\n",
        encoding="utf-8",
    )
    return calls, rmsk


def test_validate_divergence_script_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import os
    calls, rmsk = _make_synthetic_fixture(tmp_path)
    out_tsv = tmp_path / "pairs.tsv"
    script = Path(__file__).resolve().parent.parent / "scripts" / "validate_divergence_vs_rmsk.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent / "src")
    proc = subprocess.run(
        [sys.executable, str(script), "--calls", str(calls), "--rmsk", str(rmsk),
         "--nm-col", "nm_mean", "--alnlen-col", "aln_len_mean",
         "--family-col", "consensus_mei_family", "--out", str(out_tsv)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    assert out_tsv.exists()
    pairs = pd.read_csv(out_tsv, sep="\t")
    assert len(pairs) == 3
    assert {"nm_divergence", "rmsk_milli_div_fraction", "family"}.issubset(pairs.columns)


def test_validate_divergence_perfect_agreement(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import os
    calls = tmp_path / "calls.tsv"
    calls.write_text(
        "chrom\twindow_start\twindow_end\tconsensus_mei_family\tnm_mean\taln_len_mean\n"
        "chr22\t100050\t100250\tLINE1\t50\t100\n"
        "chr22\t200100\t200300\tLINE1\t100\t100\n",
        encoding="utf-8",
    )
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "0\tfoo\t50\tbaz\tqux\tchr22\t100000\t100500\t.\t+\tL1HS\tLINE\tLINE1\n"
        "0\tfoo\t100\tbaz\tqux\tchr22\t200000\t200500\t.\t+\tL1HS\tLINE\tLINE1\n",
        encoding="utf-8",
    )
    out_tsv = tmp_path / "pairs.tsv"
    script = Path(__file__).resolve().parent.parent / "scripts" / "validate_divergence_vs_rmsk.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent / "src")
    proc = subprocess.run(
        [sys.executable, str(script), "--calls", str(calls), "--rmsk", str(rmsk),
         "--nm-col", "nm_mean", "--alnlen-col", "aln_len_mean",
         "--family-col", "consensus_mei_family", "--out", str(out_tsv)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    pairs = pd.read_csv(out_tsv, sep="\t")
    assert len(pairs) == 2
    r = np.corrcoef(pairs["nm_divergence"], pairs["rmsk_milli_div_fraction"])[0, 1]
    assert abs(r - 1.0) < 1e-9
