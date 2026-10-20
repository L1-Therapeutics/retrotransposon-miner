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
    """Return the winning MEI family for each detail row.

    Prefers ``mei_target`` over ``mate_mei_target`` when both are present.
    """
    mei = detail["mei_target"].fillna("").astype(str).map(_normalize_mei_family_token) if "mei_target" in detail.columns else pd.Series("", index=detail.index)
    mate = (
        detail["mate_mei_target"].fillna("").astype(str).map(_normalize_mei_family_token)
        if "mate_mei_target" in detail.columns
        else pd.Series("", index=detail.index)
    )
    return mei.where(mei.ne(""), mate)


def _hit_frame(detail: pd.DataFrame) -> pd.DataFrame:
    """Filter detail to mei_hit rows with positive coordinates."""
    hit = detail.loc[detail["mei_hit"].fillna(False).astype(bool)].copy()
    for col in ("mei_start", "mei_end"):
        hit[col] = pd.to_numeric(hit[col], errors="coerce")
    return hit.loc[hit["mei_start"].gt(0) & hit["mei_end"].gt(0)]


def _unfiltered_extent(detail: pd.DataFrame) -> tuple[int, int]:
    hit = _hit_frame(detail)
    return int(hit["mei_start"].min()), int(hit["mei_end"].max())


def _family_extent(detail: pd.DataFrame, family: str) -> tuple[int, int]:
    hit = _hit_frame(detail)
    fam = _row_family(hit)
    hit = hit.loc[fam.eq(family)]
    if hit.empty:
        raise ValueError(f"no {family} mei_hit rows in detail")
    return int(hit["mei_start"].min()), int(hit["mei_end"].max())


def _family_robust_extent(detail: pd.DataFrame, family: str) -> tuple[int, int]:
    hit = _hit_frame(detail)
    fam = _row_family(hit)
    hit = hit.loc[fam.eq(family)]
    if hit.empty:
        raise ValueError(f"no {family} mei_hit rows in detail")
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
    assert lo == int(manifest["old_5p"]), (
        f"{catalog_id}: unfiltered min {lo} != manifest old_5p {manifest['old_5p']}"
    )
    assert hi == int(manifest["old_3p"]), (
        f"{catalog_id}: unfiltered max {hi} != manifest old_3p {manifest['old_3p']}"
    )
    assert (hi - lo + 1) == int(manifest["old_span"]), (
        f"{catalog_id}: unfiltered span {hi - lo + 1} != manifest old_span {manifest['old_span']}"
    )


@pytest.mark.parametrize("catalog_id", REAL_LOCI)
def test_aggregate_extents_use_winning_family_only(catalog_id: str):
    manifest = _load_manifest(catalog_id)
    detail = _load_detail(catalog_id)
    family = str(manifest["family"])
    want_lo, want_hi = _family_robust_extent(detail, family)
    extents = _aggregate_detail_mei_extents(detail)
    assert not extents.empty, f"{catalog_id}: _aggregate_detail_mei_extents returned empty"
    row = extents.iloc[0]
    got_lo = int(pd.to_numeric(row["detail_mei_start_min"], errors="coerce"))
    got_hi = int(pd.to_numeric(row["detail_mei_end_max"], errors="coerce"))
    assert got_lo == want_lo, (
        f"{catalog_id}: robust min {got_lo} != expected {want_lo} for family {family}"
    )
    assert got_hi == want_hi, (
        f"{catalog_id}: robust max {got_hi} != expected {want_hi} for family {family}"
    )
    mixed_lo, mixed_hi = _unfiltered_extent(detail)
    if mixed_hi > want_hi or mixed_lo < want_lo:
        assert (got_hi - got_lo + 1) < (mixed_hi - mixed_lo + 1), (
            f"{catalog_id}: family-restricted span {got_hi - got_lo + 1} not smaller than mixed span {mixed_hi - mixed_lo + 1}"
        )


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
    hit = _hit_frame(detail)
    hit["mei_end"] = pd.to_numeric(hit["mei_end"], errors="coerce")
    outlier = hit.loc[hit["mei_end"].ge(1000)]
    assert not outlier.empty, f"{catalog_id}: expected outlier mei_end >= 1000, none found"
    assert outlier["read_name"].nunique() == 1, (
        f"{catalog_id}: expected 1 outlier read, found {outlier['read_name'].nunique()}"
    )
    lo, hi = _unfiltered_extent(detail)
    assert hi == 1520, f"{catalog_id}: unfiltered max {hi} != 1520"
    extents = _aggregate_detail_mei_extents(detail)
    got_hi = int(pd.to_numeric(extents.iloc[0]["detail_mei_end_max"], errors="coerce"))
    got_lo = int(pd.to_numeric(extents.iloc[0]["detail_mei_start_min"], errors="coerce"))
    assert got_hi <= 576, f"{catalog_id}: robust max {got_hi} > 576"
    assert got_lo >= 494, f"{catalog_id}: robust min {got_lo} < 494"
    assert (got_hi - got_lo + 1) < 200, (
        f"{catalog_id}: robust span {got_hi - got_lo + 1} >= 200"
    )
    plot = _mei_coords_from_detail(detail, family="LINE1")
    assert plot is not None, f"{catalog_id}: LINE1 family coords are None"
    assert plot[1] == got_hi, (
        f"{catalog_id}: LINE1 max {plot[1]} != robust max {got_hi}"
    )
    assert plot[0] == got_lo, (
        f"{catalog_id}: LINE1 min {plot[0]} != robust min {got_lo}"
    )


def test_rank072_scattered_l1_is_not_full_length():
    """Random ~20 bp L1 seeds must not draw a 6 kb LINE-1 axis."""
    catalog_id = "rank072_chr19_877655"
    detail = _load_detail(catalog_id)
    lo, hi = _unfiltered_extent(detail)
    assert hi >= 6000, f"{catalog_id}: unfiltered max {hi} < 6000"
    assert (hi - lo + 1) > 5000, (
        f"{catalog_id}: unfiltered span {hi - lo + 1} <= 5000"
    )
    plot = _mei_coords_from_detail(detail, family="LINE1")
    assert plot is not None, f"{catalog_id}: LINE1 family coords are None"
    assert (plot[1] - plot[0] + 1) < 1300, (
        f"{catalog_id}: LINE1 clustered span {plot[1] - plot[0] + 1} >= 1300"
    )
    extents = _aggregate_detail_mei_extents(detail)
    got_lo = int(pd.to_numeric(extents.iloc[0]["detail_mei_start_min"], errors="coerce"))
    got_hi = int(pd.to_numeric(extents.iloc[0]["detail_mei_end_max"], errors="coerce"))
    assert (got_hi - got_lo + 1) < 1300, (
        f"{catalog_id}: aggregated span {got_hi - got_lo + 1} >= 1300"
    )
