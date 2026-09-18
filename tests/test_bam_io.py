from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from retro_miner.bam_io import (
    alignment_open_mode,
    alignment_path_is_cram,
    bind_alignment_reference,
    resolve_alignment_reference,
)
from retro_miner.mei_support import _collect_indel_breakpoint_evidence


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


def test_bind_alignment_reference_exports_env(monkeypatch) -> None:
    monkeypatch.delenv("RTM_ALIGNMENT_REFERENCE", raising=False)
    assert bind_alignment_reference("/ref/hg38.fa") == "/ref/hg38.fa"
    assert os.environ["RTM_ALIGNMENT_REFERENCE"] == "/ref/hg38.fa"


class _EmptyBam:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def fetch(self, *args, **kwargs):
        return iter(())


def test_indel_collection_passes_reference_into_open_alignment(monkeypatch) -> None:
    seen: dict[str, object] = {}

    def fake_open(path, *, reference_filename=None, **kwargs):
        seen["path"] = path
        seen["reference_filename"] = reference_filename
        return _EmptyBam()

    monkeypatch.setattr("retro_miner.mei_support.open_alignment", fake_open)
    loci = pd.DataFrame({"chrom": ["chr22"], "window_start": [100], "window_end": [200]})
    _collect_indel_breakpoint_evidence(
        Path("/tmp/HG00100.final.cram"),
        loci,
        sample="disease",
        reference_fasta=Path("/ref/hg38.fa"),
    )
    assert seen["reference_filename"] == Path("/ref/hg38.fa")
