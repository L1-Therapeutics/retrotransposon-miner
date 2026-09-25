"""COMPLEX_INS keep/drop on real HG00100 chr19 sentinel fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from retro_miner.mei_support import _assign_gold_stage, _compute_insertion_model_scores

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "loci"

# 002: TSD-looking 1KG miss. MEI=12 / DPE=64 = 0.1875, rescued by MEI >= 8.
# A remaining COMPLEX_INS silver is extracted as the stay-drop control.
REAL_LOCI = (
    ("SvimAsm00156362", False, True),
)


def _fixture_dir(catalog_id: str) -> Path:
    return FIXTURE_ROOT / catalog_id


def _load_manifest(catalog_id: str) -> dict:
    return json.loads((_fixture_dir(catalog_id) / "manifest.json").read_text())


def _named_detail_read(catalog_id: str, breakpoint: int) -> tuple[str, int]:
    detail = pd.read_csv(_fixture_dir(catalog_id) / "supporting_reads_detail.mei.tsv", sep="\t", low_memory=False)
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
    pick = preferred.sort_values("_dist").iloc[0]
    read_name = str(pick["read_name"])
    read_pos = int(pd.to_numeric(pick.get("genomic_pos", pick.get("pos")), errors="coerce") or breakpoint)
    return read_name, read_pos


def _score_complex_ins_only(catalog_id: str) -> tuple[dict, pd.DataFrame]:
    manifest = _load_manifest(catalog_id)
    locus = pd.read_csv(_fixture_dir(catalog_id) / "gold_locus.tsv", sep="\t")
    assert not locus.empty
    # Isolate COMPLEX_INS from the later two-sided flank gate.
    drop = [
        col
        for col in locus.columns
        if col.endswith("_left_flank_mei_reads")
        or col.endswith("_right_flank_mei_reads")
        or col.endswith("_left_flank_polya_reads")
        or col.endswith("_right_flank_polya_reads")
    ]
    locus = locus.drop(columns=drop)
    scored = _compute_insertion_model_scores(locus)
    scored["silver_stage_pass"] = True
    scored["analysis_stage_tier"] = "silver"
    gold = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
    return manifest, gold


def _discover_loci() -> list[tuple[str, bool, bool]]:
    found: list[tuple[str, bool, bool]] = []
    for catalog_id, expect_complex, expect_gold in REAL_LOCI:
        if (_fixture_dir(catalog_id) / "gold_locus.tsv").exists():
            found.append((catalog_id, expect_complex, expect_gold))
    extra = FIXTURE_ROOT / "chr19_complex_ins_stay"
    if (extra / "gold_locus.tsv").exists() and extra.name not in {row[0] for row in found}:
        man = json.loads((extra / "manifest.json").read_text())
        found.append((extra.name, bool(man.get("expect_complex_ins", True)), bool(man.get("expect_gold", False))))
    return found


@pytest.mark.parametrize("catalog_id,expect_complex,expect_gold", _discover_loci())
def test_real_locus_complex_ins_rule(catalog_id: str, expect_complex: bool, expect_gold: bool):
    manifest, gold = _score_complex_ins_only(catalog_id)
    klass = str(gold.loc[0, "insertion_event_class"])
    if expect_complex:
        assert klass == "COMPLEX_INS"
        assert bool(gold.loc[0, "gold_stage_pass"]) is False
        assert "complex_ins_non_mei" in str(gold.loc[0, "gold_stage_fail_reason"])
    else:
        assert klass != "COMPLEX_INS"
        assert "complex_ins_non_mei" not in str(gold.loc[0, "gold_stage_fail_reason"])
        assert bool(gold.loc[0, "gold_stage_pass"]) is expect_gold
        assert bool(manifest.get("expect_complex_ins", True)) is False


@pytest.mark.parametrize("catalog_id,expect_complex,expect_gold", _discover_loci())
def test_real_locus_bam_has_named_support_read(catalog_id: str, expect_complex: bool, expect_gold: bool):
    del expect_complex, expect_gold
    manifest = _load_manifest(catalog_id)
    bam_path = _fixture_dir(catalog_id) / manifest["bam"]
    read_name, _read_pos = _named_detail_read(catalog_id, int(manifest["expected_breakpoint"]))
    names = {aln.query_name for aln in pysam.AlignmentFile(str(bam_path), "rb")}
    assert read_name in names
