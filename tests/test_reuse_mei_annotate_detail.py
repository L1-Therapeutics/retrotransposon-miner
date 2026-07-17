"""Reuse prior supporting_reads_detail to hydrate MEI hits without remapping."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _build_supporting_reads_detail_table,
    _hydrate_sample_mei_hits_from_detail,
)


def test_hydrate_marks_detail_dpe_hits_and_keeps_non_hits():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "s1",
                "clip_side": "L",
                "pos": 12,
            }
        ]
    )
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "d_hit",
                "pos": 15,
                "discordant_reasons": "interchrom",
            },
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "d_miss",
                "pos": 16,
                "discordant_reasons": "interchrom",
            },
        ]
    )
    detail = pd.DataFrame(
        [
            {
                "sample": "disease",
                "evidence_type": "DPE",
                "read_name": "d_hit",
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "mei_target": "AluY#SINE/Alu",
                "mei_start": 1,
                "mei_end": 50,
                "mei_strand": "+",
                "mate_mei_target": "",
                "mate_mei_start": 0,
                "mate_mei_end": 0,
                "mate_mei_strand": "",
                "mei_hit": True,
                "mate_mei_hit": False,
                "vntr_rescue": False,
                "polya_rescue": False,
                "mei_hit_source": "",
            },
            {
                "sample": "disease",
                "evidence_type": "SR",
                "read_name": "s1",
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "mei_target": "L1HS_full#LINE/L1",
                "mei_start": 100,
                "mei_end": 150,
                "mei_strand": "-",
                "mei_hit": True,
                "mate_mei_hit": False,
            },
        ]
    )
    out = _hydrate_sample_mei_hits_from_detail(
        sample="disease",
        split_df=split,
        discordant_df=disc,
        detail=detail,
    )
    assert out["split_summary"].paf_hits == 1
    assert out["disc_summary"].paf_hits == 1
    assert bool(out["split_hits"].loc[0, "mei_hit"]) is True
    assert out["split_hits"].loc[0, "target"] == "L1HS_full#LINE/L1"
    assert out["split_hits"].loc[0, "target_strand"] == "-"
    hit = out["disc_hits"].loc[out["disc_hits"]["read_name"] == "d_hit"].iloc[0]
    miss = out["disc_hits"].loc[out["disc_hits"]["read_name"] == "d_miss"].iloc[0]
    assert bool(hit["mei_hit"]) is True
    assert hit["target"] == "AluY#SINE/Alu"
    assert hit["target_strand"] == "+"
    assert bool(miss["mei_hit"]) is False
    assert len(out["disc_hits"]) == 2


def test_build_detail_persists_strands_and_hydrate_roundtrips():
    split = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "s1",
                "clip_side": "L",
                "pos": 12,
                "mei_hit": True,
                "target": "AluYb8#SINE/Alu",
                "target_start": 10,
                "target_end": 40,
                "target_strand": "+",
            }
        ]
    )
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 10,
                "window_end": 20,
                "read_name": "d1",
                "pos": 11,
                "mate_chrom": "chr2",
                "mate_pos": 99,
                "mei_hit": False,
                "mate_mei_hit": True,
                "vntr_rescue": False,
                "polya_rescue": False,
                "target": "",
                "target_start": 0,
                "target_end": 0,
                "target_strand": "",
                "mate_mei_target": "L1HS_3end#LINE/L1",
                "mate_mei_start": 500,
                "mate_mei_end": 700,
                "mate_mei_strand": "-",
            }
        ]
    )
    detail = _build_supporting_reads_detail_table(
        split_hits=split,
        discordant_hits=disc,
        discordant_mate_hits=pd.DataFrame(),
        sample="disease",
    )
    assert set(detail["mei_strand"]).issuperset({"+", ""})
    sr = detail.loc[detail["evidence_type"].eq("SR")].iloc[0]
    dpe = detail.loc[detail["evidence_type"].eq("DPE")].iloc[0]
    assert sr["mei_strand"] == "+"
    assert dpe["mate_mei_strand"] == "-"

    hydrated = _hydrate_sample_mei_hits_from_detail(
        sample="disease",
        split_df=split.drop(columns=["target", "target_start", "target_end", "target_strand", "mei_hit"]),
        discordant_df=disc.drop(
            columns=[
                "target",
                "target_start",
                "target_end",
                "target_strand",
                "mate_mei_target",
                "mate_mei_start",
                "mate_mei_end",
                "mate_mei_strand",
                "mei_hit",
                "mate_mei_hit",
                "vntr_rescue",
                "polya_rescue",
            ]
        ),
        detail=detail,
    )
    assert hydrated["split_hits"].loc[0, "target_strand"] == "+"
    hit = hydrated["disc_hits"].iloc[0]
    assert hit["mate_mei_strand"] == "-"
    assert hit["mate_mei_target"] == "L1HS_3end#LINE/L1"
