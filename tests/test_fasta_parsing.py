"""Regression tests: blank/whitespace-only FASTA headers must not raise IndexError."""

from __future__ import annotations

from pathlib import Path

from retro_miner.igv_plots import _iter_fasta_records as _igv_iter_fasta_records
from retro_miner.local_assembly import _iter_fasta_records as _asm_iter_fasta_records

_BLANK_HEADER_FASTA = ">\nACGT\n>   \nTGCA"


def _write(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "blank_headers.fasta"
    path.write_text(content, encoding="utf-8")
    return path


def test_igv_iter_fasta_records_handles_blank_header(tmp_path: Path):
    path = _write(tmp_path, _BLANK_HEADER_FASTA)
    records = _igv_iter_fasta_records(path)
    assert records == [("unnamed_record", "ACGT"), ("unnamed_record", "TGCA")]


def test_local_assembly_iter_fasta_records_handles_blank_header(tmp_path: Path):
    path = _write(tmp_path, _BLANK_HEADER_FASTA)
    records = _asm_iter_fasta_records(path)
    assert records == [("unnamed_record", "ACGT"), ("unnamed_record", "TGCA")]


def test_igv_iter_fasta_records_still_parses_named_headers(tmp_path: Path):
    path = _write(tmp_path, ">seq1\nACGT\n>seq2\nTGCA")
    records = _igv_iter_fasta_records(path)
    assert records == [("seq1", "ACGT"), ("seq2", "TGCA")]
