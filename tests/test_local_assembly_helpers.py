"""Unit tests for subprocess-timeout handling in local_assembly.py.

All tests use mocked subprocess.run so that neither SPAdes nor minimap2
need to be installed on the test machine.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from retro_miner.local_assembly import _run_minimap2_paf, _run_spades


# ---------------------------------------------------------------------------
# _run_spades — subprocess timeout handling
# ---------------------------------------------------------------------------

class TestRunSpadesTimeout:
    """_run_spades returns a non-zero code and informative message on timeout."""

    def test_timeout_returns_nonzero_exit_code(self, tmp_path: Path) -> None:
        """subprocess.TimeoutExpired is caught; return code must be non-zero."""
        fastq = tmp_path / "reads.fastq.gz"
        fastq.write_bytes(b"")
        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd=["spades.py"], timeout=600),
        ):
            rc, _msg = _run_spades("spades.py", fastq, tmp_path / "out", threads=1, memory_gb=1)
        assert rc != 0

    def test_timeout_message_mentions_spades(self, tmp_path: Path) -> None:
        """The returned message names SPAdes so callers can surface a useful log entry."""
        fastq = tmp_path / "reads.fastq.gz"
        fastq.write_bytes(b"")
        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd=["spades.py"], timeout=600),
        ):
            _rc, msg = _run_spades("spades.py", fastq, tmp_path / "out", threads=1, memory_gb=1)
        assert "spades" in msg.lower()

    def test_timeout_message_mentions_seconds(self, tmp_path: Path) -> None:
        """The returned message includes the timeout threshold for easy diagnosis."""
        fastq = tmp_path / "reads.fastq.gz"
        fastq.write_bytes(b"")
        with patch(
            "subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd=["spades.py"], timeout=600),
        ):
            _rc, msg = _run_spades("spades.py", fastq, tmp_path / "out", threads=1, memory_gb=1)
        assert "600" in msg

    def test_success_path_still_returns_returncode(self, tmp_path: Path) -> None:
        """A non-timeout run still returns (returncode, combined_output)."""
        fastq = tmp_path / "reads.fastq.gz"
        fastq.write_bytes(b"")
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "SPAdes finished"
        mock_proc.stderr = ""
        with patch("subprocess.run", return_value=mock_proc):
            rc, msg = _run_spades("spades.py", fastq, tmp_path / "out", threads=1, memory_gb=1)
        assert rc == 0
        assert "SPAdes finished" in msg


# ---------------------------------------------------------------------------
# _run_minimap2_paf — subprocess timeout handling
# ---------------------------------------------------------------------------

class TestRunMinimap2PafTimeout:
    """_run_minimap2_paf returns [] when the alignment subprocess times out."""

    def _make_run_side_effect(self, *, timeout_on_align: bool) -> object:
        """Return a subprocess.run side_effect:
        - index-build call ('-d' in cmd): returns returncode=1 (failure, no index cached)
        - main alignment call: raises TimeoutExpired when timeout_on_align is True
        """
        def _side_effect(cmd, **kwargs):
            if "-d" in cmd:
                r = MagicMock()
                r.returncode = 1
                return r
            if timeout_on_align:
                raise subprocess.TimeoutExpired(cmd=cmd, timeout=120)
            r = MagicMock()
            r.returncode = 0
            r.stdout = ""
            return r

        return _side_effect

    def test_alignment_timeout_raises_runtime_error(self, tmp_path: Path) -> None:
        """TimeoutExpired on the alignment subprocess is re-raised as RuntimeError."""
        query = tmp_path / "query.fa"
        target = tmp_path / "target.fa"
        query.write_text(">q1\nACGT\n", encoding="utf-8")
        target.write_text(">t1\nACGT\n", encoding="utf-8")
        with patch("subprocess.run", side_effect=self._make_run_side_effect(timeout_on_align=True)):
            with pytest.raises(RuntimeError, match="minimap2 alignment timed out"):
                _run_minimap2_paf(query, target, preset="asm5", threads=1)

    def test_alignment_empty_stdout_returns_empty_list(self, tmp_path: Path) -> None:
        """An alignment run with no output rows correctly returns []."""
        query = tmp_path / "query.fa"
        target = tmp_path / "target.fa"
        query.write_text(">q1\nACGT\n", encoding="utf-8")
        target.write_text(">t1\nACGT\n", encoding="utf-8")
        with patch("subprocess.run", side_effect=self._make_run_side_effect(timeout_on_align=False)):
            result = _run_minimap2_paf(query, target, preset="asm5", threads=1)
        assert result == []

    def test_missing_input_files_returns_empty_list(self, tmp_path: Path) -> None:
        """If either input file does not exist, [] is returned without calling subprocess."""
        query = tmp_path / "missing_query.fa"
        target = tmp_path / "missing_target.fa"
        result = _run_minimap2_paf(query, target)
        assert result == []
