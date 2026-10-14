"""Unit tests for pure-logic helpers in local_assembly.

These cover the manifest-cache helpers that decide whether an assembly can be
reused.  No BAM/FASTA files are required: ``_manifest_interval`` normalises a
parsed ``assembly_manifest.json`` mapping, and ``_has_assembled_cache`` scans
for the SPAdes output directory named after the manifest pad.
"""
from __future__ import annotations

from pathlib import Path

from retro_miner.local_assembly import _has_assembled_cache, _manifest_interval


# ─────────────────────────────────────────────────────────────────────────────
# _manifest_interval
# ─────────────────────────────────────────────────────────────────────────────

class TestManifestInterval:
    def test_returns_mapping_when_present(self):
        manifest = {"interval": {"pad_bp": 5000}}
        assert _manifest_interval(manifest) == {"pad_bp": 5000}

    def test_missing_key_returns_empty(self):
        assert _manifest_interval({"status": "assembled"}) == {}

    def test_non_mapping_returns_empty(self):
        # A hand-edited/corrupt manifest may have `interval` as a scalar.
        assert _manifest_interval({"interval": "chr1:1-2"}) == {}
        assert _manifest_interval({"interval": None}) == {}
        assert _manifest_interval({"interval": ["pad"]}) == {}


# ─────────────────────────────────────────────────────────────────────────────
# _has_assembled_cache
# ─────────────────────────────────────────────────────────────────────────────

def _touch_contigs(locus_dir: Path, sample: str, pad_bp: int) -> None:
    spades_dir = locus_dir / f"{sample}.spades.pad{pad_bp}"
    spades_dir.mkdir(parents=True, exist_ok=True)
    (spades_dir / "contigs.fasta").write_text(">ctg1\nACGT\n", encoding="utf-8")


class TestHasAssembledCache:
    def test_none_manifest_is_not_a_hit(self, tmp_path: Path):
        assert _has_assembled_cache(tmp_path, None, interval_pad_bp=5000) is False

    def test_non_assembled_status_is_not_a_cache_hit(self, tmp_path: Path):
        manifest = {"status": "failed"}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is False

    def test_status_prefix_match(self):
        # `startswith("assembled")` intentionally accepts the longer status
        # string the caller round-trips through the manifest.
        manifest = {"status": "assembled_retry"}
        assert _has_assembled_cache(Path("/nonexistent"), manifest, interval_pad_bp=5000) is False

    def test_missing_contigs_returns_false(self, tmp_path: Path):
        manifest = {"status": "assembled", "interval": {"pad_bp": 5000}}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is False

    def test_finds_disease_contigs_at_manifest_pad(self, tmp_path):
        _touch_contigs(tmp_path, "disease", pad_bp=2500)
        manifest = {"status": "assembled", "interval": {"pad_bp": 2500}}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is True

    def test_finds_control_contigs_at_manifest_pad(self, tmp_path):
        _touch_contigs(tmp_path, "control", pad_bp=2500)
        manifest = {"status": "assembled", "interval": {"pad_bp": 2500}}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is True

    def test_ignores_contigs_at_other_pad(self, tmp_path):
        _touch_contigs(tmp_path, "disease", pad_bp=5000)
        manifest = {"status": "assembled", "interval": {"pad_bp": 2500}}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is False

    def test_malformed_interval_falls_back_to_caller_pad(self, tmp_path):
        _touch_contigs(tmp_path, "disease", pad_bp=5000)
        manifest = {"status": "assembled", "interval": "not-a-mapping"}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is True

    def test_missing_interval_falls_back_to_caller_pad(self, tmp_path):
        _touch_contigs(tmp_path, "control", pad_bp=5000)
        manifest = {"status": "assembled"}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is True

    def test_non_numeric_pad_falls_back_to_caller_pad(self, tmp_path):
        _touch_contigs(tmp_path, "disease", pad_bp=5000)
        manifest = {"status": "assembled", "interval": {"pad_bp": "oops"}}
        assert _has_assembled_cache(tmp_path, manifest, interval_pad_bp=5000) is True
