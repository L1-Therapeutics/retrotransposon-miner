"""Segdup flag uses the inferred breakpoint window, not the discovery cluster."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.candidate_loci import annotate_segdup_on_breakpoint_windows
from retro_miner.mei_support import _assign_bronze_silver_stages, _assign_gold_stage

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"
CATALOG_ID = "nssv14077938"


def _fixture_dir() -> Path:
    return FIXTURE_ROOT / CATALOG_ID


def _load_manifest() -> dict:
    return json.loads((_fixture_dir() / "manifest.json").read_text())


def _named_detail_read(breakpoint: int) -> str:
    detail = pd.read_csv(_fixture_dir() / "supporting_reads_detail.mei.tsv", sep="\t", low_memory=False)
    names = detail["read_name"].fillna("").astype(str)
    detail = detail.loc[names.str.len() > 0].copy()
    pos = pd.to_numeric(
        detail["genomic_pos"] if "genomic_pos" in detail.columns else detail["pos"],
        errors="coerce",
    )
    detail["_dist"] = (pos - int(breakpoint)).abs()
    pick = detail.sort_values("_dist").iloc[0]
    return str(pick["read_name"])


def _score_locus(*, use_discovery_window: bool) -> pd.DataFrame:
    manifest = _load_manifest()
    locus = pd.read_csv(_fixture_dir() / "gold_locus.tsv", sep="\t")
    assert not locus.empty
    drop = [col for col in locus.columns if "flank_" in col]
    locus = locus.drop(columns=drop, errors="ignore")
    if use_discovery_window:
        locus["insertion_breakpoint_interval_start"] = int(manifest["discovery_window_start"])
        locus["insertion_breakpoint_interval_end"] = int(manifest["discovery_window_end"])
        locus["insertion_breakpoint_pos"] = -1
    flagged = annotate_segdup_on_breakpoint_windows(
        locus,
        segdup_bed=_fixture_dir() / manifest["segdup_bed"],
        min_fraction=0.1,
    )
    silver = _assign_bronze_silver_stages(flagged)
    return _assign_gold_stage(silver, empirical_stage=False, min_mei_mapped=3)


@pytest.mark.skipif(not (_fixture_dir() / "gold_locus.tsv").exists(), reason="missing 006 fixture")
def test_nssv14077938_discovery_window_still_hits_the_old_segdup_clip():
    gold = _score_locus(use_discovery_window=True)
    assert bool(gold.loc[0, "flag_segdup"]) is True
    assert bool(gold.loc[0, "silver_stage_pass"]) is False
    assert "junk_region_flagged" in str(gold.loc[0, "silver_stage_fail_reason"])
    assert bool(gold.loc[0, "gold_stage_pass"]) is False


@pytest.mark.skipif(not (_fixture_dir() / "gold_locus.tsv").exists(), reason="missing 006 fixture")
def test_nssv14077938_breakpoint_window_is_segdup_clean_and_gold():
    manifest = _load_manifest()
    gold = _score_locus(use_discovery_window=False)
    assert str(gold.loc[0, "segdup_query_source"]) == "breakpoint_interval"
    assert bool(gold.loc[0, "flag_segdup"]) is False
    assert bool(gold.loc[0, "silver_stage_pass"]) is True
    assert bool(gold.loc[0, "gold_stage_pass"]) is True
    assert bool(manifest.get("expect_gold")) is True
    assert bool(manifest.get("expect_flag_segdup")) is False


@pytest.mark.skipif(not (_fixture_dir() / "gold_locus.tsv").exists(), reason="missing 006 fixture")
def test_nssv14077938_bam_has_named_support_read():
    manifest = _load_manifest()
    bam_path = _fixture_dir() / manifest["bam"]
    read_name = _named_detail_read(int(manifest["expected_breakpoint"]))
    names = {aln.query_name for aln in pysam.AlignmentFile(str(bam_path), "rb")}
    assert read_name in names
