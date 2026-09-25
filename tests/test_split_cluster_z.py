"""Split-cluster binomial z on chr18 sentinels, with a gold drop below z = -2."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import _SPLIT_CLUSTER_Z_CUTOFF, _assign_gold_stage, _build_gold_review_table
from retro_miner.split_cluster_z import (
    BINOMIAL_Z_COL,
    CLUSTERED_READS_COL,
    WINDOW_READS_COL,
    annotate_split_cluster_binomial_z,
    binomial_z,
)

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"
# Read-weighted cluster rate among chr18 silver calls with at least 8 splits.
CHR18_SILVER_RATE = 0.2353
SENTINELS = (
    "rank001_chr18_30220890",
    "rank019_chr18_57632912",
    "rank109_chr18_31960576",
)


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _load_call(catalog_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    gold = pd.read_csv(_fixture_dir(catalog_id) / "gold_locus.tsv", sep="\t")
    gold["silver_stage_pass"] = True
    splits = pd.read_parquet(_fixture_dir(catalog_id) / "split_evidence.disease.parquet")
    return gold, splits


def test_rank19_split_reads_do_not_cluster() -> None:
    gold, splits = _load_call("rank019_chr18_57632912")
    scored = annotate_split_cluster_binomial_z(gold, splits)
    assert int(scored[WINDOW_READS_COL].iloc[0]) == 43
    assert int(scored[CLUSTERED_READS_COL].iloc[0]) == 0
    assert binomial_z(0, 43, CHR18_SILVER_RATE) < -3.0


def test_rank1_split_reads_cluster() -> None:
    gold, splits = _load_call("rank001_chr18_30220890")
    scored = annotate_split_cluster_binomial_z(gold, splits)
    assert int(scored[WINDOW_READS_COL].iloc[0]) == 41
    assert int(scored[CLUSTERED_READS_COL].iloc[0]) == 40
    assert binomial_z(40, 41, CHR18_SILVER_RATE) > 0.0


def test_rank109_is_below_minus_2_at_chr18_silver_rate() -> None:
    gold, splits = _load_call("rank109_chr18_31960576")
    scored = annotate_split_cluster_binomial_z(gold, splits)
    n = int(scored[WINDOW_READS_COL].iloc[0])
    k = int(scored[CLUSTERED_READS_COL].iloc[0])
    assert n == 38
    assert k == 3
    assert binomial_z(k, n, CHR18_SILVER_RATE) < _SPLIT_CLUSTER_Z_CUTOFF


def test_binomial_z_depends_on_split_count() -> None:
    assert binomial_z(0, 60, CHR18_SILVER_RATE) < binomial_z(0, 8, CHR18_SILVER_RATE)
    assert binomial_z(0, 8, CHR18_SILVER_RATE) > _SPLIT_CLUSTER_Z_CUTOFF
    assert binomial_z(0, 43, CHR18_SILVER_RATE) < _SPLIT_CLUSTER_Z_CUTOFF


def test_z_below_minus_2_drops_gold_and_keeps_rank1() -> None:
    frames = []
    splits = []
    for catalog_id in SENTINELS:
        gold, split = _load_call(catalog_id)
        frames.append(gold)
        splits.append(split)
    scored = annotate_split_cluster_binomial_z(
        pd.concat(frames, ignore_index=True),
        pd.concat(splits, ignore_index=True),
    )
    scored["disease_local_bam_peak_depth"] = 40.0
    scored["control_local_bam_peak_depth"] = 40.0
    out = _assign_gold_stage(scored, empirical_stage=False)
    by_window = {int(row.discovery_window_start): row for row in out.itertuples(index=False)}
    rank1 = by_window[30220369]
    rank19 = by_window[57632449]
    rank109 = by_window[31960184]
    assert bool(rank1.gold_stage_pass)
    assert float(rank1.split_cluster_binomial_z) > _SPLIT_CLUSTER_Z_CUTOFF
    assert not bool(rank19.gold_stage_pass)
    assert "split_cluster_low_z" in str(rank19.gold_stage_fail_reason)
    assert not bool(rank109.gold_stage_pass)
    assert "split_cluster_low_z" in str(rank109.gold_stage_fail_reason)
    review = _build_gold_review_table(out, empirical_stage=False)
    assert BINOMIAL_Z_COL in review.columns
    review_gold = {
        int(row.discovery_window_start)
        for row in review.itertuples(index=False)
        if str(row.analysis_stage_tier).lower() == "gold"
    }
    assert 30220369 in review_gold
    assert 57632449 not in review_gold
    assert 31960184 not in review_gold


@pytest.mark.parametrize("catalog_id", SENTINELS)
def test_fixture_bam_contains_named_detail_read(catalog_id: str) -> None:
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / str(manifest["bam"])
    detail = pd.read_csv(_fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv", sep="\t")
    names = detail["read_name"].fillna("").astype(str)
    detail = detail.loc[names.str.len() > 0]
    read_name = str(detail.iloc[0]["read_name"])
    bam = pysam.AlignmentFile(str(bam_path))
    try:
        bam_names = {read.query_name for read in bam.fetch()}
    finally:
        bam.close()
    assert read_name in bam_names
