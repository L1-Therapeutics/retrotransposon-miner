"""MEI span min/max is restricted to the locus winning family."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import _aggregate_detail_mei_extents, _normalize_mei_family_token
from retro_miner.read_architecture import _clustered_coord_extent, _mei_coords_from_detail

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"

# chr19 gold-review ranks whose published span mixed Alu/SVA with L1 coords.
REAL_LOCI = (
    "rank010_chr19_3989864",
    "SvimAsm00157599",
    "rank020_chr19_3128000",
    "rank021_chr19_23466909",
    "rank072_chr19_877655",
)


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _load_detail(catalog_id: str) -> pd.DataFrame:
    path = _fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv"
    detail = pd.read_csv(path, sep="\t", low_memory=False)
    manifest = _load_manifest(catalog_id)
    chrom = str(manifest["chrom"])
    window_start = int(manifest["discovery_window_start"])
    window_end = int(manifest["discovery_window_end"])
    detail["chrom"] = chrom
    if "window_start" not in detail.columns:
        detail["window_start"] = window_start
        detail["window_end"] = window_end
    else:
        tagged = detail.loc[
            detail["chrom"].astype(str).eq(chrom)
            & pd.to_numeric(detail["window_start"], errors="coerce").eq(window_start)
        ]
        if not tagged.empty:
            detail = tagged
        else:
            detail["window_start"] = window_start
            detail["window_end"] = window_end
    return detail


def _row_family(detail: pd.DataFrame) -> pd.Series:
    mei = detail["mei_target"].fillna("").astype(str).map(_normalize_mei_family_token) if "mei_target" in detail.columns else pd.Series("", index=detail.index)
    mate = (
        detail["mate_mei_target"].fillna("").astype(str).map(_normalize_mei_family_token)
        if "mate_mei_target" in detail.columns
        else pd.Series("", index=detail.index)
    )
    return mei.where(mei.ne(""), mate)


def _unfiltered_extent(detail: pd.DataFrame) -> tuple[int, int]:
    hit = detail.loc[detail["mei_hit"].fillna(False).astype(bool)].copy()
    lo = pd.to_numeric(hit["mei_start"], errors="coerce")
    hi = pd.to_numeric(hit["mei_end"], errors="coerce")
    lo = lo.where(lo.gt(0))
    hi = hi.where(hi.gt(0))
    return int(lo.min()), int(hi.max())


def _family_extent(detail: pd.DataFrame, family: str) -> tuple[int, int]:
    hit = detail.loc[detail["mei_hit"].fillna(False).astype(bool)].copy()
    fam = _row_family(hit)
    hit = hit.loc[fam.eq(family)]
    lo = pd.to_numeric(hit["mei_start"], errors="coerce")
    hi = pd.to_numeric(hit["mei_end"], errors="coerce")
    lo = lo.where(lo.gt(0))
    hi = hi.where(hi.gt(0))
    return int(lo.min()), int(hi.max())


def _family_robust_extent(detail: pd.DataFrame, family: str) -> tuple[int, int]:
    hit = detail.loc[detail["mei_hit"].fillna(False).astype(bool)].copy()
    fam = _row_family(hit)
    hit = hit.loc[fam.eq(family)]
    names = hit["read_name"] if "read_name" in hit.columns else None
    lo, hi = _clustered_coord_extent(
        pd.to_numeric(hit["mei_start"], errors="coerce"),
        pd.to_numeric(hit["mei_end"], errors="coerce"),
        names,
    )
    return int(lo), int(hi)


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_unfiltered_minmax_replays_published_mixed_span(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    detail = _load_detail(catalog_id)
    lo, hi = _unfiltered_extent(detail)
    assert lo == int(manifest["old_5p"])
    assert hi == int(manifest["old_3p"])
    assert (hi - lo + 1) == int(manifest["old_span"])


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_aggregate_extents_use_winning_family_only(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    detail = _load_detail(catalog_id)
    family = str(manifest["family"])
    want_lo, want_hi = _family_robust_extent(detail, family)
    extents = _aggregate_detail_mei_extents(detail)
    assert not extents.empty
    row = extents.iloc[0]
    got_lo = int(pd.to_numeric(row["detail_mei_start_min"], errors="coerce"))
    got_hi = int(pd.to_numeric(row["detail_mei_end_max"], errors="coerce"))
    assert got_lo == want_lo
    assert got_hi == want_hi
    mixed_lo, mixed_hi = _unfiltered_extent(detail)
    if mixed_hi > want_hi or mixed_lo < want_lo:
        assert (got_hi - got_lo + 1) < (mixed_hi - mixed_lo + 1)


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_fixture_bam_contains_named_detail_read(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / str(manifest["bam"])
    assert bam_path.exists(), f"missing BAM snippet {bam_path}"
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


def test_rank018_plot_span_is_sva_not_l1_end():
    """Plot axis must use SVA 16–1366, not Tukey-on-mixed 16–3620."""
    catalog_id = "SvimAsm00157599"
    detail = _load_detail(catalog_id)
    plot = _mei_coords_from_detail(detail, family="SVA")
    assert plot == (16, 1366)
    slim = detail.drop(
        columns=[c for c in ("mei_target", "mate_mei_target", "family", "target") if c in detail.columns]
    )
    mixed = _mei_coords_from_detail(slim, family="SVA")
    assert mixed == (16, 1366)


def test_rank021_singleton_1520_does_not_set_span():
    """One DPE read at 1498–1520 must not stretch the L1MA1 pile at 551–573."""
    catalog_id = "rank021_chr19_23466909"
    detail = _load_detail(catalog_id)
    hit = detail.loc[detail["mei_hit"].fillna(False).astype(bool)].copy()
    hit["mei_end"] = pd.to_numeric(hit["mei_end"], errors="coerce")
    outlier = hit.loc[hit["mei_end"].ge(1000)]
    assert not outlier.empty
    assert outlier["read_name"].nunique() == 1
    lo, hi = _unfiltered_extent(detail)
    assert hi == 1520
    extents = _aggregate_detail_mei_extents(detail)
    got_hi = int(pd.to_numeric(extents.iloc[0]["detail_mei_end_max"], errors="coerce"))
    got_lo = int(pd.to_numeric(extents.iloc[0]["detail_mei_start_min"], errors="coerce"))
    assert got_hi <= 576
    assert got_lo >= 494
    assert (got_hi - got_lo + 1) < 200
    plot = _mei_coords_from_detail(detail, family="LINE1")
    assert plot is not None
    assert plot[1] == got_hi
    assert plot[0] == got_lo


def test_rank072_scattered_l1_is_not_full_length():
    """Random ~20 bp L1 seeds must not draw a 6 kb LINE-1 axis."""
    catalog_id = "rank072_chr19_877655"
    detail = _load_detail(catalog_id)
    lo, hi = _unfiltered_extent(detail)
    assert hi >= 6000
    assert (hi - lo + 1) > 5000
    plot = _mei_coords_from_detail(detail, family="LINE1")
    assert plot is not None
    assert (plot[1] - plot[0] + 1) < 200
    extents = _aggregate_detail_mei_extents(detail)
    got_lo = int(pd.to_numeric(extents.iloc[0]["detail_mei_start_min"], errors="coerce"))
    got_hi = int(pd.to_numeric(extents.iloc[0]["detail_mei_end_max"], errors="coerce"))
    assert (got_hi - got_lo + 1) < 200
