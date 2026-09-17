"""Pure-logic tests for sample callset overlap (no VCF/BAM required)."""

from __future__ import annotations

from retro_miner.callset_eval import (
    Variant,
    is_carrier_alleles,
    is_carrier_gt_string,
    is_mei_like_text,
    normalize_strand,
    strand_from_meinfo,
    label_rtm_calls,
    match_variants,
    normalize_chrom,
    normalize_mei_family,
    overlap_metrics,
    padded_interval,
    rtm_breakpoint,
    rtm_family,
    suggest_heuristic_cutoff,
    sweep_cutoffs,
)


def test_normalize_chrom_adds_prefix() -> None:
    assert normalize_chrom("22") == "chr22"
    assert normalize_chrom("chr22") == "chr22"


def test_normalize_mei_family_aliases() -> None:
    assert normalize_mei_family("<INS:ME:ALU>") == "ALU"
    assert normalize_mei_family("L1HS") == "LINE1"
    assert normalize_mei_family("SVA_E") == "SVA"
    assert normalize_mei_family("VNTR") == ""


def test_is_mei_like_and_carrier() -> None:
    assert is_mei_like_text("ALU_umary_ALU_1", "<INS:ME:ALU>")
    assert not is_mei_like_text("DEL_pindel_1", "<DEL>")
    assert is_carrier_alleles((0, 1))
    assert is_carrier_alleles((1, 1))
    assert not is_carrier_alleles((0, 0))
    assert not is_carrier_alleles((None, None))
    assert is_carrier_gt_string("0|1")
    assert is_carrier_gt_string("1/1")
    assert not is_carrier_gt_string("0|0")
    assert not is_carrier_gt_string("./.")


def test_overlap_with_pad_recovers_nearby_truth() -> None:
    truth = [Variant("chr22", 1000, 1000, "t1", "ALU", "melt")]
    calls = [Variant("chr22", 1150, 1150, "c1", "ALU", "rtm")]
    hits = match_variants(truth, calls, pad_bp=200)
    assert hits["t1"][0].variant_id == "c1"
    miss = match_variants(truth, calls, pad_bp=50)
    assert miss["t1"] == []


def test_family_filter_optional() -> None:
    truth = [Variant("chr22", 1000, 1000, "t1", "ALU", "ont")]
    calls = [Variant("chr22", 1000, 1000, "c1", "SVA", "rtm")]
    assert match_variants(truth, calls, pad_bp=0, require_family=False)["t1"]
    assert match_variants(truth, calls, pad_bp=0, require_family=True)["t1"] == []


def test_recall_and_novel() -> None:
    truth = [
        Variant("chr22", 1000, 1000, "t1", "ALU", "ont"),
        Variant("chr22", 5000, 5000, "t2", "SVA", "ont"),
    ]
    calls = [
        Variant("chr22", 1010, 1010, "c1", "ALU", "rtm"),
        Variant("chr22", 9000, 9000, "c2", "ALU", "rtm"),
    ]
    metrics = overlap_metrics(truth, calls, pad_bp=50)
    assert metrics["truth_recovered"] == 1
    assert metrics["recall"] == 0.5
    assert metrics["novel_calls"] == 1


def test_label_and_suggest_cutoff() -> None:
    melt = [Variant("chr22", 1000, 1000, "m1", "ALU", "melt")]
    ont = [Variant("chr22", 1000, 1000, "o1", "ALU", "ont")]
    calls = [
        Variant(
            "chr22",
            1000,
            1000,
            "keep",
            "ALU",
            "rtm",
            extra={"insertion_call_tier": "high_conf_two_sided", "insertion_model_score": 0.8},
        ),
        Variant(
            "chr22",
            8000,
            8000,
            "noise",
            "ALU",
            "rtm",
            extra={"insertion_call_tier": "none", "insertion_model_score": 0.1},
        ),
    ]
    labeled = label_rtm_calls(calls, melt, ont, pad_bp=50)
    assert labeled[0]["overlap_ont"] is True
    assert labeled[1]["novel"] is True
    sweep = sweep_cutoffs(labeled, melt, ont, pad_bp=50)
    suggestion = suggest_heuristic_cutoff(sweep, min_ont_recall=0.8)
    assert suggestion is not None
    assert suggestion["ont_recall"] >= 0.8
    assert suggestion["novel_calls"] <= 1


def test_rtm_breakpoint_and_family_from_row() -> None:
    assert rtm_breakpoint({"consensus_insertion_breakpoint_pos": 123}) == 123
    assert rtm_breakpoint({"window_start": 100, "window_end": 200}) == 100
    assert rtm_family({"consensus_mei_family": "AluYa5"}) == "ALU"


def test_strand_from_melt_meinfo() -> None:
    assert normalize_strand("+") == "+"
    assert normalize_strand("minus") == "-"
    assert strand_from_meinfo("SVA,48,1315,-") == "-"
    assert strand_from_meinfo("ALU,1,280,+") == "+"
    assert strand_from_meinfo("ALU") == ""


def test_padded_interval_is_half_open_positive() -> None:
    start, end = padded_interval(10, 10, 5)
    assert start == 5
    assert end == 15
