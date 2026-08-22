"""Tests for build_single_bam_igv_batch in scripts/make_locus_zoom_gif.py.

Verifies that:
  - chromosome values containing newline, carriage-return, or tab raise ValueError
    *before* any batch content is written to disk;
  - canonical chromosome names are accepted and produce a valid batch file.

The script is loaded via importlib so its top-level imports (PIL, matplotlib)
do not create a package dependency for the test suite.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "make_locus_zoom_gif.py"


def _load_zoom_module():
    """Load make_locus_zoom_gif.py as a named module.

    The module must be registered in sys.modules *before* exec_module runs
    because @dataclass looks up ``sys.modules.get(cls.__module__)`` during
    class creation to resolve type annotations.  Without this registration
    the lookup returns None and dataclass raises AttributeError.
    """
    import sys as _sys
    spec = importlib.util.spec_from_file_location("make_locus_zoom_gif", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    _sys.modules["make_locus_zoom_gif"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def zoom():
    """Loaded make_locus_zoom_gif module (loaded once per test module)."""
    return _load_zoom_module()


@pytest.fixture()
def batch_setup(tmp_path: Path):
    """Minimal file tree for build_single_bam_igv_batch (no real BAM data needed)."""
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr22\nACGT\n", encoding="utf-8")

    bam = tmp_path / "test.bam"
    bam.write_bytes(b"")
    # _resolve_bam_index looks for {bam}.bai first
    Path(f"{bam}.bai").write_bytes(b"")

    snap = tmp_path / "snap"
    snap.mkdir()

    return {"reference_fasta": ref, "bam_path": bam, "snapshot_dir": snap}


class TestBuildSingleBamIgvBatchValidation:
    """build_single_bam_igv_batch must reject chromosomes with IGV batch-forbidden chars."""

    def test_newline_in_chrom_raises_value_error(self, zoom, batch_setup):
        """A newline in the chromosome name must raise ValueError."""
        with pytest.raises(ValueError, match="forbidden character"):
            zoom.build_single_bam_igv_batch(
                reference_fasta=batch_setup["reference_fasta"],
                bam_path=batch_setup["bam_path"],
                snapshot_dir=batch_setup["snapshot_dir"],
                windows=[("chr1\nBAD:100", 1000, 2000)],
                panel_heights=[380],
            )

    def test_carriage_return_in_chrom_raises_value_error(self, zoom, batch_setup):
        """A carriage-return in the chromosome name must raise ValueError."""
        with pytest.raises(ValueError, match="forbidden character"):
            zoom.build_single_bam_igv_batch(
                reference_fasta=batch_setup["reference_fasta"],
                bam_path=batch_setup["bam_path"],
                snapshot_dir=batch_setup["snapshot_dir"],
                windows=[("chr1\rBAD:100", 1000, 2000)],
                panel_heights=[380],
            )

    def test_tab_in_chrom_raises_value_error(self, zoom, batch_setup):
        """A tab in the chromosome name must raise ValueError."""
        with pytest.raises(ValueError, match="forbidden character"):
            zoom.build_single_bam_igv_batch(
                reference_fasta=batch_setup["reference_fasta"],
                bam_path=batch_setup["bam_path"],
                snapshot_dir=batch_setup["snapshot_dir"],
                windows=[("chr1\tBAD", 1000, 2000)],
                panel_heights=[380],
            )

    def test_batch_file_not_written_on_invalid_chrom(self, zoom, batch_setup):
        """No batch file must exist after a validation failure."""
        with pytest.raises(ValueError):
            zoom.build_single_bam_igv_batch(
                reference_fasta=batch_setup["reference_fasta"],
                bam_path=batch_setup["bam_path"],
                snapshot_dir=batch_setup["snapshot_dir"],
                windows=[("chr1\nBAD", 1000, 2000)],
                panel_heights=[380],
            )
        assert not (batch_setup["snapshot_dir"] / "igv_zoom_batch.txt").exists()

    def test_canonical_chrom_does_not_raise(self, zoom, batch_setup):
        """Standard chromosome names must produce a batch file without error."""
        batch_path = zoom.build_single_bam_igv_batch(
            reference_fasta=batch_setup["reference_fasta"],
            bam_path=batch_setup["bam_path"],
            snapshot_dir=batch_setup["snapshot_dir"],
            windows=[("chr22", 1000, 2000)],
            panel_heights=[380],
        )
        assert batch_path.exists()

    def test_valid_batch_contains_correct_goto_line(self, zoom, batch_setup):
        """A canonical chromosome must produce a well-formed goto line."""
        batch_path = zoom.build_single_bam_igv_batch(
            reference_fasta=batch_setup["reference_fasta"],
            bam_path=batch_setup["bam_path"],
            snapshot_dir=batch_setup["snapshot_dir"],
            windows=[("chr22", 1000, 2000)],
            panel_heights=[380],
        )
        content = batch_path.read_text(encoding="utf-8")
        assert "goto chr22:1000-2000\n" in content

    def test_second_window_invalid_chrom_prevents_batch_write(self, zoom, batch_setup):
        """Validation failure in the second window must also prevent the batch file."""
        with pytest.raises(ValueError):
            zoom.build_single_bam_igv_batch(
                reference_fasta=batch_setup["reference_fasta"],
                bam_path=batch_setup["bam_path"],
                snapshot_dir=batch_setup["snapshot_dir"],
                windows=[("chr22", 1000, 2000), ("chr1\nBAD", 3000, 4000)],
                panel_heights=[380, 380],
            )
        assert not (batch_setup["snapshot_dir"] / "igv_zoom_batch.txt").exists()


class TestBuildSingleBamIgvBatchPathQuoting:
    """Paths in batch header lines are double-quoted to support spaces in filenames.

    IGV >=2.x supports double-quoted paths in batch scripts; without quoting
    a space in any path would cause IGV to misparse the command.
    """

    def _call(self, zoom, batch_setup):
        return zoom.build_single_bam_igv_batch(
            reference_fasta=batch_setup["reference_fasta"],
            bam_path=batch_setup["bam_path"],
            snapshot_dir=batch_setup["snapshot_dir"],
            windows=[("chr22", 1000, 2000)],
            panel_heights=[380],
        )

    def test_genome_line_is_quoted(self, zoom, batch_setup):
        """The 'genome' batch line wraps the reference FASTA path in double quotes."""
        content = self._call(zoom, batch_setup).read_text(encoding="utf-8")
        ref = str(batch_setup["reference_fasta"].resolve())
        assert f'genome "{ref}"' in content

    def test_snapshotdirectory_line_is_quoted(self, zoom, batch_setup):
        """The 'snapshotDirectory' batch line wraps the directory path in double quotes."""
        content = self._call(zoom, batch_setup).read_text(encoding="utf-8")
        snap = str(batch_setup["snapshot_dir"].resolve())
        assert f'snapshotDirectory "{snap}"' in content

    def test_load_bam_path_is_quoted(self, zoom, batch_setup):
        """The 'load' batch line wraps the BAM and index paths in double quotes."""
        content = self._call(zoom, batch_setup).read_text(encoding="utf-8")
        bam = str(batch_setup["bam_path"].resolve())
        assert f'load "{bam}"' in content

    def test_space_in_genome_path_survives_untruncated(self, zoom, tmp_path):
        """A genome path containing spaces is quoted and the full path appears in the batch."""
        space_dir = tmp_path / "ref dir with spaces"
        space_dir.mkdir()
        ref = space_dir / "ref genome.fa"
        ref.write_text(">chr22\nACGT\n", encoding="utf-8")
        bam = tmp_path / "test.bam"
        bam.write_bytes(b"")
        Path(f"{bam}.bai").write_bytes(b"")
        snap = tmp_path / "snap"
        snap.mkdir()
        batch_path = zoom.build_single_bam_igv_batch(
            reference_fasta=ref,
            bam_path=bam,
            snapshot_dir=snap,
            windows=[("chr22", 1000, 2000)],
            panel_heights=[380],
        )
        content = batch_path.read_text(encoding="utf-8")
        ref_str = str(ref.resolve())
        # Full space-containing path must appear quoted; without quoting IGV
        # would only see the text up to the first space.
        assert f'genome "{ref_str}"' in content
