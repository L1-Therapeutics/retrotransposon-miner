"""Split-cluster binomial z on the chr18 rank 19 and rank 1 sentinels."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import _build_gold_review_table
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


def test_binomial_z_depends_on_split_count() -> None:
    assert binomial_z(0, 60, CHR18_SILVER_RATE) < binomial_z(0, 8, CHR18_SILVER_RATE)
    assert binomial_z(0, 8, CHR18_SILVER_RATE) > -3.0
    assert binomial_z(0, 60, CHR18_SILVER_RATE) < -3.0


def test_score_is_written_on_gold_review() -> None:
    frames = []
    splits = []
    for catalog_id in ("rank019_chr18_57632912", "rank001_chr18_30220890"):
        gold, split = _load_call(catalog_id)
        frames.append(gold)
        splits.append(split)
    scored = annotate_split_cluster_binomial_z(pd.concat(frames, ignore_index=True), pd.concat(splits, ignore_index=True))
    review = _build_gold_review_table(scored, empirical_stage=False)
    assert BINOMIAL_Z_COL in review.columns
    assert WINDOW_READS_COL in review.columns
    by_window = {
        int(row.discovery_window_start): float(row.split_cluster_binomial_z)
        for row in review.itertuples(index=False)
    }
    assert by_window[57632449] < -3.0
    assert by_window[30220369] > by_window[57632449]


@pytest.mark.parametrize(
    "catalog_id",
    ["rank019_chr18_57632912", "rank001_chr18_30220890"],
)
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
