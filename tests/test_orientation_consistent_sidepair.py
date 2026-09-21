"""Gold two-sided flank rule on real HG00100 chr22 sentinel fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import (
    _GOLD_MIN_FLANK_MEI_READS,
    _assign_gold_stage,
    _genomic_flank_evidence_table,
)

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"

# Named gold-review sentinels used to design the rule. expect_gold is the
# replay of _assign_gold_stage on the WGS-mate annotate tables.
REAL_LOCI = (
    ("rank012_chr22_23935321", False),
    ("ALU_umary_ALU_12520", False),
    ("rank022_chr22_31226464", True),
    ("rank023_chr22_50495209", False),
    ("rank026_chr22_50351027", False),
    ("rank042_chr22_35735283", False),
    ("nssv14073986", True),
)

_FLANK_COLS = (
    "left_flank_mei_reads",
    "right_flank_mei_reads",
    "left_flank_polya_reads",
    "right_flank_polya_reads",
)


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _filter_or_tag_window(df: pd.DataFrame, chrom: str, window_start: int, window_end: int) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    work = df.copy()
    if "chrom" in work.columns and "window_start" in work.columns and "window_end" in work.columns:
        tagged = work.loc[
            work["chrom"].astype(str).eq(chrom)
            & pd.to_numeric(work["window_start"], errors="coerce").eq(window_start)
            & pd.to_numeric(work["window_end"], errors="coerce").eq(window_end)
        ]
        if not tagged.empty:
            return tagged.copy()
    work["chrom"] = chrom
    work["window_start"] = window_start
    work["window_end"] = window_end
    return work


def _load_extract(catalog_id: str, name: str, chrom: str, window_start: int, window_end: int) -> pd.DataFrame:
    path = _fixture_dir(catalog_id) / name
    if not path.exists():
        return pd.DataFrame()
    return _filter_or_tag_window(pd.read_parquet(path), chrom, window_start, window_end)


def _score_fixture_locus(catalog_id: str) -> tuple[dict, pd.DataFrame]:
    manifest = _load_manifest(catalog_id)
    gold_path = _fixture_dir(catalog_id) / "gold_locus.tsv"
    assert gold_path.exists(), f"missing gold_locus.tsv for {catalog_id}"
    locus = pd.read_csv(gold_path, sep="\t")
    assert not locus.empty, f"empty gold_locus.tsv for {catalog_id}"
    chrom = str(manifest["chrom"])
    window_start = int(manifest["discovery_window_start"])
    window_end = int(manifest["discovery_window_end"])
    bp = int(manifest["expected_breakpoint"])
    if "insertion_breakpoint_pos" in locus.columns:
        published = pd.to_numeric(locus["insertion_breakpoint_pos"], errors="coerce").dropna()
        if not published.empty and int(published.iloc[0]) > 0:
            bp = int(published.iloc[0])
    scored = locus.copy()
    scored["chrom"] = chrom
    scored["window_start"] = window_start
    scored["window_end"] = window_end
    scored["insertion_breakpoint_pos"] = bp
    drop = [f"{prefix}_{col}" for prefix in ("disease", "control") for col in _FLANK_COLS]
    scored = scored.drop(columns=[c for c in drop if c in scored.columns])
    breakpoints = scored.loc[:, ["chrom", "window_start", "window_end", "insertion_breakpoint_pos"]].drop_duplicates()
    for prefix, split_name, disc_name in (
        ("disease", "split_evidence.disease.parquet", "discordant_evidence.disease.parquet"),
        ("control", "split_evidence.control.parquet", "discordant_evidence.control.parquet"),
    ):
        split = _load_extract(catalog_id, split_name, chrom, window_start, window_end)
        disc = _load_extract(catalog_id, disc_name, chrom, window_start, window_end)
        flanks = _genomic_flank_evidence_table(split, disc, breakpoints, prefix=prefix)
        if flanks.empty:
            for col in (f"{prefix}_{c}" for c in _FLANK_COLS):
                scored[col] = 0
            continue
        scored = scored.merge(flanks, on=["chrom", "window_start", "window_end"], how="left")
        for col in (f"{prefix}_{c}" for c in _FLANK_COLS):
            scored[col] = pd.to_numeric(scored[col], errors="coerce").fillna(0).astype(int)
    return manifest, scored


def _named_detail_read(catalog_id: str, breakpoint: int) -> tuple[str, int]:
    detail_path = _fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv"
    detail = pd.read_csv(detail_path, sep="\t", low_memory=False)
    assert not detail.empty, f"empty supporting_reads_detail for {catalog_id}"
    names = detail["read_name"].fillna("").astype(str)
    detail = detail.loc[names.str.len() > 0].copy()
    pos = pd.to_numeric(
        detail["genomic_pos"] if "genomic_pos" in detail.columns else detail["pos"],
        errors="coerce",
    )
    detail["_dist"] = (pos - int(breakpoint)).abs()
    preferred = detail
    if "mei_hit" in detail.columns and bool(detail["mei_hit"].fillna(False).astype(bool).any()):
        preferred = detail.loc[detail["mei_hit"].fillna(False).astype(bool)]
    elif "evidence_type" in detail.columns:
        ev = detail["evidence_type"].fillna("").astype(str)
        hit = ev.isin(["SR", "DPE", "polyA", "POLYA"])
        if bool(hit.any()):
            preferred = detail.loc[hit]
    pick = preferred.sort_values("_dist").iloc[0]
    read_name = str(pick["read_name"])
    read_pos = int(pd.to_numeric(pick.get("genomic_pos", pick.get("pos")), errors="coerce") or breakpoint)
    return read_name, read_pos


@pytest.mark.parametrize("catalog_id,expect_gold", REAL_LOCI)
def test_real_locus_two_sided_gold_matches_chr22_replay(catalog_id: str, expect_gold: bool):
    """HG00100 chr22 extract tables, not invented SR/DPE piles."""
    manifest, scored = _score_fixture_locus(catalog_id)
    assert bool(manifest.get("expect_gold")) is expect_gold
    gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
    assert bool(gold.loc[0, "gold_stage_pass"]) is expect_gold
    if not expect_gold:
        assert "one_sided_or_inconsistent_flank_support" in str(gold.loc[0, "gold_stage_fail_reason"])


@pytest.mark.parametrize("catalog_id,expect_gold", REAL_LOCI)
def test_real_locus_old_gold_without_flank_gate(catalog_id: str, expect_gold: bool):
    """These sentinels were gold before the flank gate; dropping the columns restores that."""
    _manifest, scored = _score_fixture_locus(catalog_id)
    drop = [f"{prefix}_{col}" for prefix in ("disease", "control") for col in _FLANK_COLS]
    legacy = scored.drop(columns=[c for c in drop if c in scored.columns])
    old = _assign_gold_stage(legacy, empirical_stage=False, min_mei_mapped=3)
    assert bool(old.loc[0, "gold_stage_pass"]) is True
    if not expect_gold:
        new = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
        assert bool(new.loc[0, "gold_stage_pass"]) is False


@pytest.mark.parametrize("catalog_id,_expect_gold", REAL_LOCI)
def test_real_locus_bam_contains_named_detail_read(catalog_id: str, _expect_gold: bool):
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / str(manifest["bam"])
    assert bam_path.exists(), f"missing BAM snippet {bam_path}"
    read_name, read_pos = _named_detail_read(catalog_id, int(manifest["expected_breakpoint"]))
    bam = pysam.AlignmentFile(str(bam_path))
    try:
        names = {read.query_name for read in bam.fetch(str(manifest["chrom"]), max(0, read_pos - 80), read_pos + 80)}
        if read_name not in names:
            names = {read.query_name for read in bam.fetch()}
    finally:
        bam.close()
    assert read_name in names


def test_one_sided_sentinels_lack_multiple_mei_on_both_flanks():
    """Ranks 12/23/26/42: architecture piles sit on one genomic flank."""
    for catalog_id in (
        "rank012_chr22_23935321",
        "rank023_chr22_50495209",
        "rank026_chr22_50351027",
        "rank042_chr22_35735283",
    ):
        _manifest, scored = _score_fixture_locus(catalog_id)
        left = int(scored.loc[0, "disease_left_flank_mei_reads"])
        right = int(scored.loc[0, "disease_right_flank_mei_reads"])
        assert max(left, right) >= _GOLD_MIN_FLANK_MEI_READS
        assert min(left, right) < _GOLD_MIN_FLANK_MEI_READS


def test_alu_umary_parked_dpe_do_not_create_opposite_flank():
    """Rank 13 / ALU_umary_ALU_12520: token DPE_L/R looks two-sided; parked DPE do not."""
    catalog_id = "ALU_umary_ALU_12520"
    _manifest, scored = _score_fixture_locus(catalog_id)
    support = str(scored.loc[0, "disease_supporting_reads"])
    assert "DPE_L=" in support and "DPE_R=" in support
    left = int(scored.loc[0, "disease_left_flank_mei_reads"])
    right = int(scored.loc[0, "disease_right_flank_mei_reads"])
    assert min(left, right) < _GOLD_MIN_FLANK_MEI_READS
    gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
    assert bool(gold.loc[0, "gold_stage_pass"]) is False


def test_keep_controls_have_orientation_consistent_two_sided_support():
    """nssv14073986 and rank 22 stay gold; both genomic flanks carry MEI or MEI+polyA."""
    for catalog_id in ("nssv14073986", "rank022_chr22_31226464"):
        _manifest, scored = _score_fixture_locus(catalog_id)
        gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
        assert bool(gold.loc[0, "gold_stage_pass"]) is True
        left = int(scored.loc[0, "disease_left_flank_mei_reads"])
        right = int(scored.loc[0, "disease_right_flank_mei_reads"])
        l_poly = int(scored.loc[0, "disease_left_flank_polya_reads"])
        r_poly = int(scored.loc[0, "disease_right_flank_polya_reads"])
        two_sided_mei = left >= _GOLD_MIN_FLANK_MEI_READS and right >= _GOLD_MIN_FLANK_MEI_READS
        mei_polya = (left >= _GOLD_MIN_FLANK_MEI_READS and r_poly >= 1) or (
            right >= _GOLD_MIN_FLANK_MEI_READS and l_poly >= 1
        )
        assert two_sided_mei or mei_polya
