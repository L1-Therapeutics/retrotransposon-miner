"""Security and correctness tests for _build_assembly_contig_track (igv_plots).

Verifies that:
  - subprocess calls never use shell=True;
  - paths with spaces and shell metacharacters are passed verbatim as list elements;
  - minimap2 or samtools failure returns None (skip behavior);
  - a missing BAM after exit-0 returns None;
  - a partial BAM is cleaned up on failure;
  - a successful pipeline returns the BAM path.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from retro_miner._utils import safe_locus_id
from retro_miner.igv_plots import _build_assembly_contig_track, _validate_igv_chrom, build_igv_batch_script


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_CHROM = "chr1"
_START = 1000
_END = 2000
_CONTIG_ID = "NODE_1"


def _make_variants(
    chrom: str = _CHROM,
    start: int = _START,
    end: int = _END,
    contig_id: str = _CONTIG_ID,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": [chrom],
            "window_start": [start],
            "window_end": [end],
            "discovery_window_start": [start],
            "discovery_window_end": [end],
            "assembly_best_contig_id": [contig_id],
        }
    )


def _create_assembly_cache(
    cache_dir: Path,
    chrom: str = _CHROM,
    start: int = _START,
    end: int = _END,
    contig_id: str = _CONTIG_ID,
) -> None:
    """Create the minimal assembly-cache directory tree expected by the function."""
    locus_id = safe_locus_id(chrom, start, end)
    locus_dir = cache_dir / locus_id
    spades_dir = locus_dir / "disease.spades.pad250"
    spades_dir.mkdir(parents=True)
    (locus_dir / "assembly_manifest.json").write_text(
        '{"interval": {"pad_bp": 250}}', encoding="utf-8"
    )
    (spades_dir / "contigs.fasta").write_text(
        f">{contig_id}\nACGTACGTACGT\n", encoding="utf-8"
    )


def _fake_minimap2(returncode: int = 0) -> MagicMock:
    """Return a minimal mock representing a minimap2 Popen process."""
    proc = MagicMock()
    proc.stdout = MagicMock()  # supports .close() without error
    proc.stderr = MagicMock()
    proc.stderr.read.return_value = b"" if returncode == 0 else b"minimap2 error\n"
    proc.wait.return_value = returncode
    return proc


def _fake_samtools(returncode: int = 0, bam_path: Path | None = None) -> MagicMock:
    """Return a minimal mock representing a samtools sort Popen process."""
    proc = MagicMock()
    proc.returncode = returncode

    def _communicate(*args, **kwargs):
        if returncode == 0 and bam_path is not None:
            bam_path.write_bytes(b"MOCK BAM DATA")
        stderr = b"" if returncode == 0 else b"samtools error\n"
        return (b"", stderr)

    proc.communicate.side_effect = _communicate
    return proc


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def track_setup(tmp_path: Path):
    """Minimal valid inputs to reach the subprocess section of the function."""
    cache_dir = tmp_path / "cache"
    _create_assembly_cache(cache_dir)
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    snap = tmp_path / "snap"
    snap.mkdir()
    return {
        "variants": _make_variants(),
        "assembly_cache_dir": cache_dir,
        "reference_fasta": ref,
        "snapshot_dir": snap,
        "bam_path": snap / "assembly_selected_contigs.bam",
    }


# ---------------------------------------------------------------------------
# Helpers for building call-counting Popen side-effects
# ---------------------------------------------------------------------------


def _popen_side_effect(minimap2_rc: int = 0, samtools_rc: int = 0, bam_path: Path | None = None):
    """Return a side_effect callable for patching subprocess.Popen.

    The first call is assumed to be minimap2; the second samtools sort.
    """
    calls = [0]

    def _side(args, **kwargs):
        calls[0] += 1
        if calls[0] == 1:
            return _fake_minimap2(minimap2_rc)
        return _fake_samtools(samtools_rc, bam_path)

    return _side


# ---------------------------------------------------------------------------
# Security: no shell invocation
# ---------------------------------------------------------------------------


class TestNoShellInvocation:
    """subprocess.Popen must never be called with shell=True."""

    def test_popen_args_are_lists_not_shell_strings(self, track_setup):
        """Every subprocess.Popen call must receive a list, not a shell command string."""
        bam_path = track_setup["bam_path"]
        captured: list[dict] = []

        def _spy(args, **kwargs):
            captured.append({"args": list(args), "kwargs": dict(kwargs)})
            if len(captured) == 1:
                return _fake_minimap2(0)
            # Do NOT create the BAM — function returns None at alignment check,
            # never reaching the subprocess.run indexing step.
            return _fake_samtools(0)

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch("retro_miner.igv_plots.subprocess.Popen", side_effect=_spy):
                _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert captured, "subprocess.Popen was never called"
        for call_info in captured:
            assert call_info["kwargs"].get("shell") is not True, (
                f"Popen called with shell=True: {call_info}"
            )
            assert isinstance(call_info["args"], list), (
                f"Popen args are not a list (possible shell string): {call_info['args']!r}"
            )

    def test_subprocess_run_not_used_with_shell_true(self, track_setup):
        """subprocess.run must not be called with shell=True anywhere in the alignment flow."""
        bam_path = track_setup["bam_path"]
        shell_calls: list[dict] = []

        def _spy_run(args, **kwargs):
            if kwargs.get("shell") is True:
                shell_calls.append({"args": args, "kwargs": kwargs})
            result = MagicMock()
            result.returncode = 0
            result.stderr = ""
            result.stdout = ""
            # Simulate samtools index creating the .bai file
            if isinstance(args, list) and len(args) >= 2 and args[1] == "index":
                Path(f"{bam_path}.bai").write_bytes(b"FAKE INDEX")
            return result

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                side_effect=_popen_side_effect(0, 0, bam_path),
            ):
                with patch("retro_miner.igv_plots.subprocess.run", side_effect=_spy_run):
                    _build_assembly_contig_track(
                        track_setup["variants"],
                        assembly_cache_dir=track_setup["assembly_cache_dir"],
                        reference_fasta=track_setup["reference_fasta"],
                        snapshot_dir=track_setup["snapshot_dir"],
                    )

        assert not shell_calls, (
            f"subprocess.run was called with shell=True: {shell_calls}"
        )

    def test_minimap2_is_first_popen_arg(self, track_setup):
        """The first subprocess.Popen call must invoke minimap2 as a list element."""
        bam_path = track_setup["bam_path"]
        captured: list[list] = []

        def _spy(args, **kwargs):
            captured.append(list(args))
            if len(captured) == 1:
                return _fake_minimap2(0)
            # No BAM created — function exits at alignment check before indexing step.
            return _fake_samtools(0)

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch("retro_miner.igv_plots.subprocess.Popen", side_effect=_spy):
                _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert captured, "subprocess.Popen was never called"
        assert captured[0][0] == "minimap2", (
            f"First Popen call is not minimap2: {captured[0]}"
        )
        # The samtools call must be separate
        assert len(captured) >= 2 and captured[1][0] == "samtools", (
            f"Second Popen call is not samtools: {captured[1] if len(captured) >= 2 else '(missing)'}"
        )
        # No element of either call should be a shell pipe character
        for call_args in captured:
            for arg in call_args:
                assert "|" not in arg, f"Pipe character in argument: {arg!r}"


# ---------------------------------------------------------------------------
# Security: paths with special characters
# ---------------------------------------------------------------------------


class TestPathArgumentHandling:
    """Paths containing spaces or shell metacharacters must be passed verbatim."""

    def test_reference_path_with_spaces_is_single_argument(self, tmp_path: Path):
        """A reference FASTA at a path containing spaces must appear as one list element."""
        ref_dir = tmp_path / "ref data dir"  # space in directory name
        ref_dir.mkdir()
        ref = ref_dir / "reference genome.fa"  # space in filename
        ref.write_text(">chr1\nACGT\n", encoding="utf-8")

        cache_dir = tmp_path / "cache"
        _create_assembly_cache(cache_dir)
        snap = tmp_path / "snap"
        snap.mkdir()
        bam_path = snap / "assembly_selected_contigs.bam"

        captured_minimap2_args: list[str] = []
        call_count = [0]

        def _spy(args, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                captured_minimap2_args.extend(list(args))
                return _fake_minimap2(0)
            # No BAM created — function exits at alignment check before indexing step.
            return _fake_samtools(0)

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch("retro_miner.igv_plots.subprocess.Popen", side_effect=_spy):
                _build_assembly_contig_track(
                    _make_variants(),
                    assembly_cache_dir=cache_dir,
                    reference_fasta=ref,
                    snapshot_dir=snap,
                )

        assert captured_minimap2_args, "minimap2 Popen was never called"
        expected_ref_path = str(ref.resolve())
        assert expected_ref_path in captured_minimap2_args, (
            f"Reference path with spaces not found as a single argument.\n"
            f"Expected: {expected_ref_path!r}\n"
            f"Got args: {captured_minimap2_args}"
        )

    def test_no_arg_contains_shell_pipe_or_semicolon(self, tmp_path: Path):
        """No argument passed to minimap2 or samtools should contain | or ;."""
        cache_dir = tmp_path / "cache"
        _create_assembly_cache(cache_dir)
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\nACGT\n", encoding="utf-8")
        snap = tmp_path / "snap"
        snap.mkdir()
        bam_path = snap / "assembly_selected_contigs.bam"

        all_args: list[str] = []
        call_count = [0]

        def _spy(args, **kwargs):
            call_count[0] += 1
            all_args.extend(str(a) for a in args)
            if call_count[0] == 1:
                return _fake_minimap2(0)
            # No BAM created — function exits at alignment check before indexing step.
            return _fake_samtools(0)

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch("retro_miner.igv_plots.subprocess.Popen", side_effect=_spy):
                _build_assembly_contig_track(
                    _make_variants(),
                    assembly_cache_dir=cache_dir,
                    reference_fasta=ref,
                    snapshot_dir=snap,
                )

        assert all_args, "subprocess.Popen was never called"
        for arg in all_args:
            assert "|" not in arg, f"Pipe found in subprocess argument: {arg!r}"
            assert ";" not in arg, f"Semicolon found in subprocess argument: {arg!r}"
            assert "`" not in arg, f"Backtick found in subprocess argument: {arg!r}"


# ---------------------------------------------------------------------------
# Failure behavior: skip and return None
# ---------------------------------------------------------------------------


class TestFailureBehavior:
    """The function must return None (skip) on any process failure."""

    def test_minimap2_nonzero_exit_returns_none(self, track_setup):
        """minimap2 exiting non-zero must cause the function to return None."""
        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                side_effect=_popen_side_effect(minimap2_rc=1, samtools_rc=0, bam_path=track_setup["bam_path"]),
            ):
                result = _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert result is None

    def test_samtools_nonzero_exit_returns_none(self, track_setup):
        """samtools exiting non-zero must cause the function to return None."""
        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                # samtools fails and does NOT create the BAM
                side_effect=_popen_side_effect(minimap2_rc=0, samtools_rc=1, bam_path=None),
            ):
                result = _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert result is None

    def test_both_exit_zero_but_no_bam_returns_none(self, track_setup):
        """If both processes exit 0 but no BAM is written to disk, return None."""
        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                # bam_path=None means samtools does not create the file
                side_effect=_popen_side_effect(minimap2_rc=0, samtools_rc=0, bam_path=None),
            ):
                result = _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert result is None

    def test_partial_bam_deleted_on_samtools_failure(self, track_setup):
        """A partial BAM already on disk must be removed when the pipeline fails."""
        bam_path = track_setup["bam_path"]
        # Simulate a partial write left from a previous or current run
        bam_path.write_bytes(b"PARTIAL DATA")

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                side_effect=_popen_side_effect(minimap2_rc=0, samtools_rc=1, bam_path=None),
            ):
                result = _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert result is None
        assert not bam_path.exists(), "Partial BAM must be removed on pipeline failure"

    def test_minimap2_missing_returns_none_before_subprocess(self, track_setup):
        """If minimap2 is not in PATH, the function returns None without calling Popen."""

        def _no_minimap2(name):
            return None if name == "minimap2" else "/usr/bin/samtools"

        with patch("retro_miner.igv_plots.shutil.which", side_effect=_no_minimap2):
            with patch("retro_miner.igv_plots.subprocess.Popen") as mock_popen:
                result = _build_assembly_contig_track(
                    track_setup["variants"],
                    assembly_cache_dir=track_setup["assembly_cache_dir"],
                    reference_fasta=track_setup["reference_fasta"],
                    snapshot_dir=track_setup["snapshot_dir"],
                )

        assert result is None
        mock_popen.assert_not_called()


# ---------------------------------------------------------------------------
# Success behavior
# ---------------------------------------------------------------------------


class TestSuccessBehavior:
    """Verify the function returns the BAM path on a fully successful pipeline run."""

    def test_success_returns_bam_path(self, track_setup):
        """On success the function returns the BAM path (after samtools index)."""
        bam_path = track_setup["bam_path"]
        bai_path = Path(f"{bam_path}.bai")

        def _spy_run(args, **kwargs):
            result = MagicMock()
            result.returncode = 0
            result.stderr = ""
            result.stdout = ""
            if isinstance(args, list) and len(args) >= 2 and args[1] == "index":
                bai_path.write_bytes(b"FAKE INDEX")
            return result

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                side_effect=_popen_side_effect(0, 0, bam_path),
            ):
                with patch("retro_miner.igv_plots.subprocess.run", side_effect=_spy_run):
                    result = _build_assembly_contig_track(
                        track_setup["variants"],
                        assembly_cache_dir=track_setup["assembly_cache_dir"],
                        reference_fasta=track_setup["reference_fasta"],
                        snapshot_dir=track_setup["snapshot_dir"],
                    )

        assert result == bam_path
        assert bam_path.exists()

    def test_success_requires_index_to_exist(self, track_setup):
        """If samtools index exits 0 but creates no index file, the function returns None."""
        bam_path = track_setup["bam_path"]

        def _spy_run(args, **kwargs):
            result = MagicMock()
            result.returncode = 0
            result.stderr = ""
            result.stdout = ""
            # Intentionally do NOT create the .bai file
            return result

        with patch("retro_miner.igv_plots.shutil.which", return_value="/usr/bin/fake"):
            with patch(
                "retro_miner.igv_plots.subprocess.Popen",
                side_effect=_popen_side_effect(0, 0, bam_path),
            ):
                with patch("retro_miner.igv_plots.subprocess.run", side_effect=_spy_run):
                    result = _build_assembly_contig_track(
                        track_setup["variants"],
                        assembly_cache_dir=track_setup["assembly_cache_dir"],
                        reference_fasta=track_setup["reference_fasta"],
                        snapshot_dir=track_setup["snapshot_dir"],
                    )

        assert result is None


# ---------------------------------------------------------------------------
# IGV batch-command serialization: control-character injection
# ---------------------------------------------------------------------------


def _make_batch_variants(*chroms: str) -> pd.DataFrame:
    """Minimal variants DataFrame for build_igv_batch_script tests."""
    return pd.DataFrame(
        {
            "chrom": list(chroms),
            "window_start": [1000] * len(chroms),
            "window_end": [2000] * len(chroms),
            "discovery_window_start": [1000] * len(chroms),
            "discovery_window_end": [2000] * len(chroms),
            "assembly_best_contig_id": [""] * len(chroms),
        }
    )


@pytest.fixture()
def batch_setup(tmp_path: Path):
    """Minimal file tree for build_igv_batch_script (no real BAM/FASTA data needed)."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr22\nACGT\n", encoding="utf-8")

    disease_bam = tmp_path / "disease.bam"
    disease_bam.write_bytes(b"")
    Path(f"{disease_bam}.bai").write_bytes(b"")  # _resolve_bam_index looks for {bam}.bai first

    control_bam = tmp_path / "control.bam"
    control_bam.write_bytes(b"")
    Path(f"{control_bam}.bai").write_bytes(b"")

    snap = tmp_path / "snap"
    snap.mkdir()

    return {
        "reference_fasta": ref,
        "disease_bam": disease_bam,
        "control_bam": control_bam,
        "snapshot_dir": snap,
    }


class TestIgvBatchValidation:
    """goto batch commands must not contain control characters from chromosome data."""

    def _call(self, batch_setup, variants):
        """Call build_igv_batch_script with _estimate_panel_height mocked out."""
        with patch("retro_miner.igv_plots._estimate_panel_height", return_value=250):
            return build_igv_batch_script(variants, **batch_setup)

    def test_newline_in_chromosome_skips_row(self, batch_setup):
        """A chromosome containing \\n must not appear in the batch output."""
        batch = self._call(batch_setup, _make_batch_variants("chr22\nINJECT"))
        assert "INJECT" not in batch
        assert "goto" not in batch  # the only row was skipped

    def test_carriage_return_in_chromosome_skips_row(self, batch_setup):
        """A chromosome containing \\r must not appear in the batch output."""
        batch = self._call(batch_setup, _make_batch_variants("chr22\rINJECT"))
        assert "INJECT" not in batch
        assert "goto" not in batch

    def test_tab_in_chromosome_skips_row(self, batch_setup):
        """A chromosome containing \\t must not appear in the batch output."""
        batch = self._call(batch_setup, _make_batch_variants("chr22\tINJECT"))
        assert "INJECT" not in batch
        assert "goto" not in batch

    def test_valid_chromosome_produces_goto(self, batch_setup):
        """A canonical chromosome name must produce a properly formed goto command."""
        batch = self._call(batch_setup, _make_batch_variants("chr22"))
        assert "goto chr22:" in batch

    def test_bad_chromosome_skipped_good_chromosome_kept(self, batch_setup):
        """Rows with invalid chromosomes are skipped; valid rows are preserved."""
        batch = self._call(batch_setup, _make_batch_variants("chr22\nINJECT", "chr1"))
        assert "INJECT" not in batch
        assert "goto chr1:" in batch

    def test_no_goto_line_contains_embedded_control_chars(self, batch_setup):
        """No goto line in valid batch output may contain embedded control characters."""
        batch = self._call(batch_setup, _make_batch_variants("chr22", "chr1"))
        for line in batch.splitlines():
            if line.startswith("goto"):
                assert "\n" not in line
                assert "\r" not in line
                assert "\t" not in line

    def test_validate_igv_chrom_raises_on_newline(self):
        """_validate_igv_chrom raises ValueError for a newline-containing chromosome."""
        with pytest.raises(ValueError, match="forbidden character"):
            _validate_igv_chrom("chr22\nBAD")

    def test_validate_igv_chrom_raises_on_tab(self):
        """_validate_igv_chrom raises ValueError for a tab-containing chromosome."""
        with pytest.raises(ValueError, match="forbidden character"):
            _validate_igv_chrom("chr22\tBAD")

    def test_validate_igv_chrom_passes_for_canonical_names(self):
        """_validate_igv_chrom does not raise for standard chromosome names."""
        for chrom in ("chr1", "chr22", "chrX", "chrY", "chrM"):
            _validate_igv_chrom(chrom)  # must not raise
