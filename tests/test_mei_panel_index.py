"""MEI remap panel bwa index is built once (download) and reused at annotate."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from retro_miner.mei_panel_index import (
    PUBLIC_MEI_REMAP_FASTA_RELPATHS,
    ConsensusLength,
    ensure_bwa_index,
    ensure_mei_remap_bwa_index,
    find_panel_paths,
    full_axis_name,
    has_bwa_index,
    index_public_mei_remap_fastas,
    nopolya_sidecar_path,
    read_panel_labels,
    resolve_consensus_lengths,
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


# --- consensus-axis lengths from the panels ---------------------------------------

def _write_panel(path: Path, records: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()), encoding="utf-8")
    return path


def _panels(tmp_path: Path) -> list[tuple[str, Path]]:
    fragment = _write_panel(tmp_path / "retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta", {
        "L1HS_5end#LINE/L1": "ACGT" * 10,
        "L1M1_orf2#LINE/L1": "ACGT" * 8,
        "AluYa5#SINE/Alu": "ACGT" * 6,
    })
    full = _write_panel(tmp_path / "retrotransposon_db/full_consensus/mei_full_canonical.panel.fa", {
        "L1HS_full#LINE/L1|source=L1HS": "ACGT" * 100,
        "L1PA2_full#LINE/L1|source=L1PA2": "ACGT" * 90,
        "AluY_full#SINE/Alu|source=AluY": "ACGT" * 20,
        "SVA_F_full#Retroposon/SVA|source=SVA_F": "ACGT" * 30,
    })
    return [("fragment", fragment), ("full", full)]


def test_read_panel_labels_counts_bases_and_splits_attributes(tmp_path: Path) -> None:
    labels = read_panel_labels(_write_panel(tmp_path / "p.fa", {
        "L1HS_5end#LINE/L1": "ACGT" * 3,
        "L1PA2_full#LINE/L1|source=L1PA2": "AC" * 7,
    }))
    assert labels["L1HS_5end"] == ("LINE/L1", 12)
    # The ``|source=`` suffix the derived panel writes is not part of the class.
    assert labels["L1PA2_full"] == ("LINE/L1", 14)


def test_full_axis_name_maps_every_fragment_suffix() -> None:
    for label in ("L1HS_5end", "L1HS_3end", "L1M1_orf2", "L1MC_orf1", "AluY_short", "AluY"):
        assert full_axis_name(label).endswith("_full")
    assert full_axis_name("L1HS_5end") == "L1HS_full"
    assert full_axis_name("L1HS_full") == "L1HS_full"   # idempotent


def test_a_fragment_label_takes_the_full_axis_not_its_own_length(tmp_path: Path) -> None:
    # The caller's coordinates are on the full axis: L1HS_5end calls reach MEI_3P 6018
    # while their fragment consensus is 2136 bp.  Reading the fragment length here
    # would make a legitimate full-length call look like it overruns its consensus.
    resolved = resolve_consensus_lengths(["L1HS_5end"], _panels(tmp_path))
    hit = resolved["L1HS_5end"]
    assert hit is not None and hit.bp == 400 and hit.target == "L1HS_full"
    assert hit.axis == "full" and hit.basis == "base_full"


def test_a_label_without_a_full_record_projects_onto_the_family_canonical(tmp_path: Path) -> None:
    # No L1M1_full in the panel, and L1M1_orf2 calls reach 5677 -- past both their own
    # 3294 bp fragment and every other L1 axis except the family canonical, which is
    # what the annotate projection uses for such a fragment.
    resolved = resolve_consensus_lengths(["L1M1_orf2", "L1MB7_3end"], _panels(tmp_path))
    assert resolved["L1M1_orf2"] == ConsensusLength(
        bp=400, target="L1HS_full", axis="full", panel="full", basis="family_canonical",
        mei_class="LINE/L1")
    # A fragment canonical also carries the class through from the fragment panel.
    assert resolved["L1MB7_3end"] is not None
    assert resolved["L1MB7_3end"].basis == "family_canonical"  # type: ignore[union-attr]


def test_a_label_the_panels_do_not_know_resolves_to_nothing(tmp_path: Path) -> None:
    # None is "unknown", which callers must keep apart from a length of zero: the
    # scripts fall back to their data-derived estimate and say so.
    resolved = resolve_consensus_lengths(["NotARepeat_x9"], _panels(tmp_path))
    assert resolved["NotARepeat_x9"] is None
    assert resolve_consensus_lengths([], _panels(tmp_path)) == {}


def test_family_canonical_covers_alu_and_sva_too(tmp_path: Path) -> None:
    resolved = resolve_consensus_lengths(["AluJb", "SVA_A"], _panels(tmp_path))
    assert resolved["AluJb"] is not None and resolved["AluJb"].target == "AluY_full"
    assert resolved["SVA_A"] is not None and resolved["SVA_A"].target == "SVA_F_full"


def test_find_panel_paths_honours_env_then_explicit_then_root(tmp_path: Path) -> None:
    panels = _panels(tmp_path)
    fragment_relpath = PUBLIC_MEI_REMAP_FASTA_RELPATHS[0]
    found = find_panel_paths(public_data_dir=tmp_path)
    assert [kind for kind, _path in found] == ["fragment", "full"]

    # An explicit path wins for the fragment panel; the env var covers the other.
    other = _write_panel(tmp_path / "explicit.fa", {"L1HS_5end#LINE/L1": "ACGT"})
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setenv("RTM_MEI_FULL_PANEL_FASTA", str(panels[1][1]))
        found = find_panel_paths([other], public_data_dir=tmp_path)
        assert found[0][1] == other
        assert found[1][1] == panels[1][1]
        assert (tmp_path / fragment_relpath).exists()
    finally:
        monkeypatch.undo()


def test_find_panel_paths_returns_nothing_when_no_panel_exists(tmp_path: Path) -> None:
    assert find_panel_paths(public_data_dir=tmp_path / "empty") == []
    assert resolve_consensus_lengths(["L1HS_5end"], [])["L1HS_5end"] is None


def test_download_script_falls_back_to_the_repeatmasker_library() -> None:
    src = Path(__file__).resolve().parents[1] / "scripts" / "download_public_data.py"
    text = src.read_text(encoding="utf-8")
    # The h5/famdb route stays first; both of its failure points must reach the builder
    # that needs neither, or the panel stays absent exactly as before.
    assert "from fetch_dfam_mei_panel import PANEL_RELPATH, build_panel" in text
    assert text.count("_prepare_dfam_mei_subset_without_famdb(") >= 3
    assert 'reason="missing_archive"' in text
    assert "reason=tool_status" in text


def test_span_and_coverage_scripts_read_the_panel_axis() -> None:
    scripts = Path(__file__).resolve().parents[1] / "scripts" / "analysis"
    span = (scripts / "line1_span_support.py").read_text(encoding="utf-8")
    coverage = (scripts / "consensus_coverage.py").read_text(encoding="utf-8")
    for text in (span, coverage):
        assert "from retro_miner.mei_panel_index import" in text
        assert "find_panel_paths" in text
        assert "resolve_consensus_lengths" in text
    # The data-derived estimate survives only as a named fallback.
    assert "PROXY_SOURCE = \"proxy_max_MEI_3P\"" in span
    assert "OBSERVED_SOURCE = \"data_observed_max_end\"" in coverage
    assert "consensus_source" in span and "axis_source" in coverage
