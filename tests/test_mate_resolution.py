from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from click import ClickException

from retro_miner.mate_resolution import (
    infer_full_genome_mate_bam,
    interchrom_mate_seq_stats,
    require_interchrom_mate_sequences,
)


def test_infer_full_genome_mate_bam_prefers_bam_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    work = tmp_path / "workdir"
    stage = work / "data" / "bam_stage"
    stage.mkdir(parents=True)
    cram = stage / "HG00100.final.cram"
    cram.write_bytes(b"CRAM")
    slice_bam = tmp_path / "hg00100.shortread.chr22.hg38.bam"
    slice_bam.write_bytes(b"BAM")
    monkeypatch.delenv("RTM_DISEASE_MATE_BAM", raising=False)
    monkeypatch.delenv("RTM_MATE_BAM", raising=False)
    monkeypatch.delenv("RTM_BAM_STAGE_DIR", raising=False)
    monkeypatch.delenv("RTM_PUBLIC_DATA_DIR", raising=False)
    found = infer_full_genome_mate_bam(slice_bam, workdir=work)
    assert found == cram


def test_infer_full_genome_mate_bam_env_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    explicit = tmp_path / "other.cram"
    explicit.write_bytes(b"CRAM")
    monkeypatch.setenv("RTM_DISEASE_MATE_BAM", str(explicit))
    found = infer_full_genome_mate_bam(tmp_path / "hg00100.shortread.chr22.hg38.bam")
    assert found == explicit


def test_require_interchrom_mate_sequences_fails_when_off_chrom_empty() -> None:
    df = pd.DataFrame(
        {
            "chrom": ["chr22"] * 12,
            "mate_chrom": ["chr2"] * 12,
            "mate_seq": [""] * 12,
        }
    )
    with pytest.raises(ClickException, match="interchromosomal discordant mates"):
        require_interchrom_mate_sequences(df, mate_bam="/slice.bam")


def test_require_interchrom_mate_sequences_ok_when_filled() -> None:
    df = pd.DataFrame(
        {
            "chrom": ["chr22", "chr22"],
            "mate_chrom": ["chr2", "chr22"],
            "mate_seq": ["ACGT" * 10, ""],
        }
    )
    stats = require_interchrom_mate_sequences(df)
    assert stats.off_chrom == 1
    assert stats.off_empty == 0


def test_interchrom_stats_ignore_same_chrom() -> None:
    df = pd.DataFrame(
        {
            "chrom": ["chr22"] * 20,
            "mate_chrom": ["chr22"] * 20,
            "mate_seq": [""] * 20,
        }
    )
    stats = interchrom_mate_seq_stats(df)
    assert stats.off_chrom == 0
    require_interchrom_mate_sequences(df)
