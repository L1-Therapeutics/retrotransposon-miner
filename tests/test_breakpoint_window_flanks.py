"""Gold flanks use the inferred breakpoint, not the discovery-bin midpoint."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import (
    _GOLD_MIN_FLANK_MEI_READS,
    _assign_gold_stage,
    _inferred_flank_breakpoint_series,
    annotate_flanks_on_breakpoint_windows,
)

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"

# Five chr19 unique-1KG sites that looked two-sided in architecture plots
# but dropped gold when flanks were split at the discovery-window midpoint.
REAL_LOCI = (
    "SvimAsm00155226",
    "SvimAsm00155806",
    "L1_umary_LINE1_2883",
    "SvimAsm00159362",
    "SvimAsm00159431",
)
# 001 DPE-only is still two-sided around this bin midpoint; the other four
# replay the one-sided DPE split that dropped gold.
MIDPOINT_DROP_LOCI = tuple(c for c in REAL_LOCI if c != "SvimAsm00155226")


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _retag(df: pd.DataFrame, chrom: str, window_start: int, window_end: int) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    work = df.copy()
    work["chrom"] = chrom
    work["window_start"] = int(window_start)
    work["window_end"] = int(window_end)
    return work


def _evidence_frames(catalog_id: str, chrom: str, window_start: int, window_end: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail_path = _fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv"
    if detail_path.exists():
        detail = pd.read_csv(detail_path, sep="\t", low_memory=False)
        ev = detail["evidence_type"].fillna("").astype(str).str.upper()
        split = _retag(detail.loc[ev.eq("SR")].copy(), chrom, window_start, window_end)
        disc = _retag(detail.loc[ev.eq("DPE")].copy(), chrom, window_start, window_end)
        if not split.empty or not disc.empty:
            return split, disc
    split = pd.read_parquet(_fixture_dir(catalog_id) / "split_evidence.disease.parquet")
    disc = pd.read_parquet(_fixture_dir(catalog_id) / "discordant_evidence.disease.parquet")
    return _retag(split, chrom, window_start, window_end), _retag(disc, chrom, window_start, window_end)


def _score_locus(catalog_id: str, *, use_discovery_midpoint: bool) -> pd.DataFrame:
    manifest = _load_manifest(catalog_id)
    locus = pd.read_csv(_fixture_dir(catalog_id) / "gold_locus.tsv", sep="\t")
    assert not locus.empty
    chrom = str(manifest["chrom"])
    window_start = int(manifest["discovery_window_start"])
    window_end = int(manifest["discovery_window_end"])
    drop = [c for c in locus.columns if "flank_" in c]
    locus = locus.drop(columns=drop, errors="ignore")
    locus["chrom"] = chrom
    locus["window_start"] = window_start
    locus["window_end"] = window_end
    # Stacked gold_locus rows are silver; isolate the flank gate.
    locus["silver_stage_pass"] = True
    locus["junk_flag_count"] = 0
    locus["insertion_event_class"] = "SIMPLE_MEI"
    if "disease_mei_mapped" in locus.columns:
        locus["disease_mei_mapped"] = pd.to_numeric(locus["disease_mei_mapped"], errors="coerce").fillna(3).clip(lower=3)
    if "control_mei_mapped" in locus.columns:
        locus["control_mei_mapped"] = pd.to_numeric(locus["control_mei_mapped"], errors="coerce").fillna(3).clip(lower=3)
    if use_discovery_midpoint:
        locus["insertion_breakpoint_pos"] = -1
        locus["consensus_insertion_breakpoint_pos"] = -1
        for col in (
            "insertion_breakpoint_interval_start",
            "insertion_breakpoint_interval_end",
            "consensus_breakpoint_interval_start",
            "consensus_breakpoint_interval_end",
        ):
            if col in locus.columns:
                locus[col] = -1
    else:
        bp = int(manifest["expected_breakpoint"])
        locus["insertion_breakpoint_pos"] = bp
        locus["consensus_insertion_breakpoint_pos"] = bp
    split, disc = _evidence_frames(catalog_id, chrom, window_start, window_end)
    scored = annotate_flanks_on_breakpoint_windows(
        locus,
        split_disease=split,
        split_control=split,
        discordant_disease=disc,
        discordant_control=disc,
    )
    return _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)


def _named_detail_read(catalog_id: str, breakpoint: int) -> str:
    detail = pd.read_csv(_fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv", sep="\t", low_memory=False)
    names = detail["read_name"].fillna("").astype(str)
    detail = detail.loc[names.str.len() > 0].copy()
    pos = pd.to_numeric(
        detail["genomic_pos"] if "genomic_pos" in detail.columns else detail["pos"],
        errors="coerce",
    )
    detail["_dist"] = (pos - int(breakpoint)).abs()
    pick = detail.sort_values("_dist").iloc[0]
    return str(pick["read_name"])


def _score_dpe_only(catalog_id: str, *, use_discovery_midpoint: bool) -> pd.DataFrame:
    """Replay the annotate-time DPE split (SR clip geometry does not use the BP)."""
    manifest = _load_manifest(catalog_id)
    scored = _score_locus(catalog_id, use_discovery_midpoint=use_discovery_midpoint)
    # Recompute flanks from DPE only so the discovery-midpoint bug is isolated.
    chrom = str(manifest["chrom"])
    window_start = int(manifest["discovery_window_start"])
    window_end = int(manifest["discovery_window_end"])
    _split, disc = _evidence_frames(catalog_id, chrom, window_start, window_end)
    locus = scored.drop(columns=[c for c in scored.columns if "flank_" in c], errors="ignore")
    return annotate_flanks_on_breakpoint_windows(
        locus,
        split_disease=pd.DataFrame(),
        split_control=pd.DataFrame(),
        discordant_disease=disc,
        discordant_control=disc,
    )


def test_inferred_breakpoint_prefers_junction_not_discovery_mid():
    manifest = _load_manifest("SvimAsm00155226")
    row = pd.DataFrame(
        [
            {
                "window_start": int(manifest["discovery_window_start"]),
                "window_end": int(manifest["discovery_window_end"]),
                "insertion_breakpoint_pos": -1,
                "consensus_insertion_breakpoint_pos": int(manifest["expected_breakpoint"]),
            }
        ]
    )
    inferred = int(_inferred_flank_breakpoint_series(row).iloc[0])
    mid = (int(manifest["discovery_window_start"]) + int(manifest["discovery_window_end"])) // 2
    assert inferred == int(manifest["expected_breakpoint"])
    row["consensus_insertion_breakpoint_pos"] = -1
    assert int(_inferred_flank_breakpoint_series(row).iloc[0]) == mid


@pytest.mark.parametrize("catalog_id", MIDPOINT_DROP_LOCI)
def test_discovery_midpoint_still_drops_the_five_1kg_sites(catalog_id: str):
    gold = _assign_gold_stage(
        _score_dpe_only(catalog_id, use_discovery_midpoint=True),
        empirical_stage=False,
        min_mei_mapped=3,
    )
    assert bool(gold.loc[0, "gold_stage_pass"]) is False
    assert "one_sided_or_inconsistent_flank_support" in str(gold.loc[0, "gold_stage_fail_reason"])
    left = int(gold.loc[0, "disease_left_flank_mei_reads"])
    right = int(gold.loc[0, "disease_right_flank_mei_reads"])
    assert min(left, right) < _GOLD_MIN_FLANK_MEI_READS


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_inferred_breakpoint_flanks_restore_gold(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    gold = _score_locus(catalog_id, use_discovery_midpoint=False)
    assert bool(manifest.get("expect_gold")) is True
    assert bool(gold.loc[0, "gold_stage_pass"]) is True
    left = int(gold.loc[0, "disease_left_flank_mei_reads"])
    right = int(gold.loc[0, "disease_right_flank_mei_reads"])
    assert min(left, right) >= _GOLD_MIN_FLANK_MEI_READS


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_fixture_bam_contains_named_detail_read(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / str(manifest["bam"])
    assert bam_path.exists(), f"missing BAM snippet {bam_path}"
    read_name = _named_detail_read(catalog_id, int(manifest["expected_breakpoint"]))
    bam = pysam.AlignmentFile(str(bam_path))
    try:
        names = {read.query_name for read in bam.fetch()}
    finally:
        bam.close()
    assert read_name in names
