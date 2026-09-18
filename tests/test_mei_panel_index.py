"""MEI remap panel bwa index is built once (download) and reused at annotate."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from retro_miner.mei_panel_index import (
    PUBLIC_MEI_REMAP_FASTA_RELPATHS,
    ensure_bwa_index,
    ensure_mei_remap_bwa_index,
    has_bwa_index,
    index_public_mei_remap_fastas,
    nopolya_sidecar_path,
)


def _write_fasta(path: Path, name: str = "AluY#SINE/Alu", seq: str = "ACGT" * 20) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{seq}\n", encoding="utf-8")
    return path


def _touch_index(fasta: Path) -> None:
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{fasta}{suffix}").write_text("idx", encoding="utf-8")


def test_has_bwa_index_false_until_all_sidecars_exist(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "panel.fa")
    assert has_bwa_index(fasta) is False
    Path(f"{fasta}.bwt").write_text("x", encoding="utf-8")
    Path(f"{fasta}.sa").write_text("x", encoding="utf-8")
    assert has_bwa_index(fasta) is False
    _touch_index(fasta)
    assert has_bwa_index(fasta) is True


def test_ensure_bwa_index_skips_when_complete(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "panel.fa")
    _touch_index(fasta)
    with patch("retro_miner.mei_panel_index.subprocess.run") as run:
        result = ensure_bwa_index(fasta)
    assert result["status"] == "skipped_exists"
    run.assert_not_called()


def test_ensure_bwa_index_runs_bwa_when_missing(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "panel.fa")

    def _fake_bwa(cmd, **_kwargs):
        assert cmd[:2] == ["bwa", "index"]
        _touch_index(Path(cmd[2]))

    with (
        patch("retro_miner.mei_panel_index.shutil.which", return_value="/usr/bin/bwa"),
        patch("retro_miner.mei_panel_index.subprocess.run", side_effect=_fake_bwa) as run,
    ):
        result = ensure_bwa_index(fasta)
    assert result["status"] == "indexed"
    run.assert_called_once()
    assert has_bwa_index(fasta)


def test_ensure_mei_remap_bwa_index_writes_nopolya_then_indexes(tmp_path: Path) -> None:
    src = _write_fasta(tmp_path / "dfam_human_mei_l1_alu_sva.fasta", seq="ACGT" + ("A" * 20))
    sidecar = nopolya_sidecar_path(src)

    def _fake_bwa(cmd, **_kwargs):
        _touch_index(Path(cmd[2]))

    with (
        patch("retro_miner.mei_panel_index.shutil.which", return_value="/usr/bin/bwa"),
        patch("retro_miner.mei_panel_index.subprocess.run", side_effect=_fake_bwa),
    ):
        result = ensure_mei_remap_bwa_index(src)
    assert sidecar.exists()
    assert result["nopolya_fasta"] == str(sidecar)
    assert has_bwa_index(sidecar)
    with patch("retro_miner.mei_panel_index.subprocess.run") as run:
        again = ensure_mei_remap_bwa_index(src)
    assert again["status"] == "skipped_exists"
    run.assert_not_called()


def test_index_public_mei_remap_fastas_indexes_present_dfam_only(tmp_path: Path) -> None:
    dfam_rel = PUBLIC_MEI_REMAP_FASTA_RELPATHS[0]
    _write_fasta(tmp_path / dfam_rel)

    def _fake_bwa(cmd, **_kwargs):
        _touch_index(Path(cmd[2]))

    with (
        patch("retro_miner.mei_panel_index.shutil.which", return_value="/usr/bin/bwa"),
        patch("retro_miner.mei_panel_index.subprocess.run", side_effect=_fake_bwa),
    ):
        result = index_public_mei_remap_fastas(tmp_path)
    assert result["status"] == "ok"
    assert len(result["indexed"]) == 1
    assert dfam_rel not in result["skipped_missing"]
    assert PUBLIC_MEI_REMAP_FASTA_RELPATHS[1] in result["skipped_missing"]


def test_download_script_wires_postprocess_step() -> None:
    src = Path(__file__).resolve().parents[1] / "scripts" / "download_public_data.py"
    text = src.read_text(encoding="utf-8")
    assert "from retro_miner.mei_panel_index import index_public_mei_remap_fastas" in text
    assert '"index_mei_remap_bwa"' in text
    assert "index_public_mei_remap_fastas(outdir, force=force)" in text


def test_annotate_indexes_panel_before_any_remaps() -> None:
    src = Path(__file__).resolve().parents[1] / "src" / "retro_miner" / "mei_support.py"
    text = src.read_text(encoding="utf-8")
    marker = "ensuring MEI panel bwa index before remaps"
    assert marker in text
    after = text.split(marker, 1)[1]
    before_pool, _sep, _rest = after.partition("ThreadPoolExecutor(max_workers=2)")
    assert "ensure_mei_remap_bwa_index(Path(mei_fasta))" in before_pool
