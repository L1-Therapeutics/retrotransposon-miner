"""Family-level MEI overlap pile is a gold-table training feature, not a cutoff."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import (
    _assign_gold_stage,
    _build_gold_review_table,
    annotate_mei_overlap_piles,
)

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"

# chr19 sentinels: scattered rank 072 vs real Alu/SVA/L1 piles.
REAL_LOCI = (
    ("rank072_chr19_877655", 3),
    ("rank021_chr19_23466909", 17),
    ("SvimAsm00157599", 17),
    ("rank020_chr19_3128000", 18),
    ("rank010_chr19_3989864", 24),
)


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _load_locus(catalog_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    gold = pd.read_csv(_fixture_dir(catalog_id) / "gold_locus.tsv", sep="\t")
    detail = pd.read_csv(_fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv", sep="\t")
    return gold, detail


@pytest.mark.parametrize("catalog_id,expected_overlap", REAL_LOCI)
def test_overlap_pile_is_written_on_gold_review(catalog_id: str, expected_overlap: int) -> None:
    gold, detail = _load_locus(catalog_id)
    scored = annotate_mei_overlap_piles(gold, detail)
    assert "mei_consensus_overlap_reads" in scored.columns
    assert int(scored["mei_consensus_overlap_reads"].iloc[0]) == expected_overlap
    review = _build_gold_review_table(scored, empirical_stage=False)
    assert "mei_consensus_overlap_reads" in review.columns
    assert int(review["mei_consensus_overlap_reads"].iloc[0]) == expected_overlap


def test_low_overlap_does_not_fail_gold() -> None:
    gold, detail = _load_locus("rank072_chr19_877655")
    scored = annotate_mei_overlap_piles(gold, detail)
    assert int(scored["mei_consensus_overlap_reads"].iloc[0]) == 3
    scored["mei_consensus_overlap_reads"] = 1
    scored["silver_stage_pass"] = True
    out = _assign_gold_stage(scored, empirical_stage=False)
    reasons = out["gold_stage_fail_reason"].fillna("").astype(str).iloc[0]
    assert "scattered_mei_no_overlap_pile" not in reasons
    assert bool(out["gold_stage_pass"].iloc[0])


@pytest.mark.parametrize("catalog_id,expected_overlap", REAL_LOCI)
def test_fixture_bam_contains_named_detail_read(catalog_id: str, expected_overlap: int) -> None:
    del expected_overlap
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / str(manifest["bam"])
    if not bam_path.exists():
        pytest.skip(f"missing BAM snippet {bam_path.name}")
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
