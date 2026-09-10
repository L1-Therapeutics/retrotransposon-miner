from __future__ import annotations

from retro_miner.bam_io import (
    alignment_open_mode,
    alignment_path_is_cram,
    resolve_alignment_reference,
)


def test_alignment_open_mode_bam_vs_cram() -> None:
    assert alignment_path_is_cram("HG00100.final.cram") is True
    assert alignment_path_is_cram("/tmp/HG00100.final.bam") is False
    assert alignment_open_mode("s3://bucket/x.cram") == "rc"
    assert alignment_open_mode("https://example.com/x.cram") == "rc"
    assert alignment_open_mode("/data/x.bam") == "rb"


def test_resolve_alignment_reference_env(monkeypatch) -> None:
    monkeypatch.delenv("RTM_ALIGNMENT_REFERENCE", raising=False)
    assert resolve_alignment_reference(None) is None
    monkeypatch.setenv("RTM_ALIGNMENT_REFERENCE", "/ref/hg38.fa")
    assert resolve_alignment_reference(None) == "/ref/hg38.fa"
    assert resolve_alignment_reference("/explicit.fa") == "/explicit.fa"
