"""CLI regression tests for numeric option boundary validation."""

from __future__ import annotations

from click.testing import CliRunner

from retro_miner.cli import cli


def test_extract_split_evidence_rejects_negative_min_mapq():
    result = CliRunner().invoke(cli, ["extract-split-evidence", "--min-mapq", "-1"])
    assert result.exit_code == 2
    assert "Invalid value for '--min-mapq'" in result.output


def test_extract_split_evidence_rejects_zero_min_clip_len():
    result = CliRunner().invoke(cli, ["extract-split-evidence", "--min-clip-len", "0"])
    assert result.exit_code == 2
    assert "Invalid value for '--min-clip-len'" in result.output


def test_build_candidate_loci_rejects_zero_window_size():
    result = CliRunner().invoke(cli, ["build-candidate-loci", "--window-size", "0"])
    assert result.exit_code == 2
    assert "Invalid value for '--window-size'" in result.output
