"""Pure-logic tests for sample callset overlap (no VCF/BAM required)."""

from __future__ import annotations

from types import SimpleNamespace

import retro_miner.callset_eval as callset_eval

from retro_miner.callset_eval import (
    Variant,
    catalog_overlap,
    filter_variants_to_regions,
    is_carrier_alleles,
    is_carrier_gt_string,
    is_insertion_svtype,
    is_mei_like_text,
    is_melt_mei_insertion,
    normalize_strand,
    strand_from_meinfo,
    label_rtm_calls,
    match_variants,
    merge_unique_events,
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
    assert normalize_mei_family("BI_GS_DEL1_B1_P4126_65") == ""


def test_only_mei_insertions_are_truth_records() -> None:
    assert is_insertion_svtype("INS", "<INS:ME:ALU>")
    assert is_insertion_svtype("ALU", "<INS:ME:ALU>")
    assert not is_insertion_svtype("DEL", "<DEL>")
    assert is_melt_mei_insertion(
        "ALU_umary_ALU_12446", "<INS:ME:ALU>", "ALU", "ALU,1,280,+"
    )
    assert not is_melt_mei_insertion("DEL_ALU_1", "<DEL>", "DEL")
    assert not is_melt_mei_insertion("BI_GS_DEL1_B1_P4126_65", "<DEL>", "DEL")


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


def test_shared_callable_filter_and_unique_events(tmp_path) -> None:
    callable_bed = tmp_path / "callable.bed"
    callable_bed.write_text("chr22\t900\t6000\n", encoding="utf-8")
    junk_bed = tmp_path / "junk.bed"
    junk_bed.write_text("chr22\t4900\t5100\n", encoding="utf-8")
    melt = [
        Variant("chr22", 1000, 1000, "m1", "ALU", "melt"),
        Variant("chr22", 5000, 5000, "m_junk", "SVA", "melt"),
    ]
    ont = [
        Variant("22", 1050, 1050, "o1", "ALU", "ont"),
        Variant("chr22", 3000, 3000, "o2", "LINE1", "ont"),
    ]
    region_args = {"include_beds": [callable_bed], "exclude_beds": [junk_bed]}
    melt = filter_variants_to_regions(melt, **region_args)
    ont = filter_variants_to_regions(ont, **region_args)
    assert [variant.variant_id for variant in melt] == ["m1"]
    assert [variant.variant_id for variant in ont] == ["o1", "o2"]

    overlap = catalog_overlap(melt, ont, pad_bp=200, require_family=True)
    assert [variant.variant_id for variant in overlap["shared"]] == ["m1"]
    assert [variant.variant_id for variant in overlap["right_only"]] == ["o2"]

    events = merge_unique_events([melt, ont], pad_bp=200, require_family=True)
    assert len(events) == 2
    assert events[0].source == "melt,ont"
    assert events[0].extra["source_ids"] == {"melt": ["m1"], "ont": ["o1"]}
    assert events[1].variant_id == "ont:o2"


def test_whole_genome_query_omits_region(monkeypatch, tmp_path) -> None:
    commands = []
    monkeypatch.setattr(callset_eval.shutil, "which", lambda _: "/usr/bin/bcftools")

    def fake_run(command, **_kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="chr1\t100\nchr22\t200\n", stderr="")

    monkeypatch.setattr(callset_eval.subprocess, "run", fake_run)
    rows = callset_eval._query_vcf_rows(
        tmp_path / "input.bcf",
        chrom=None,
        sample="HG00100",
        fmt="%CHROM\t%POS\n",
    )
    assert rows == [["chr1", "100"], ["chr22", "200"]]
    assert "-r" not in commands[0]
    assert commands[0][commands[0].index("-s") + 1] == "HG00100"


def test_query_retries_chromosome_alias_after_empty_success(monkeypatch, tmp_path) -> None:
    commands = []
    monkeypatch.setattr(callset_eval.shutil, "which", lambda _: "/usr/bin/bcftools")

    def fake_run(command, **_kwargs):
        commands.append(command)
        stdout = "" if command[command.index("-r") + 1] == "chr22" else "22\t100\n"
        return SimpleNamespace(returncode=0, stdout=stdout, stderr="")

    monkeypatch.setattr(callset_eval.subprocess, "run", fake_run)
    rows = callset_eval._query_vcf_rows(
        tmp_path / "input.vcf.gz",
        chrom="chr22",
        sample="HG00100",
        fmt="%CHROM\t%POS\n",
    )
    assert rows == [["22", "100"]]
    assert [command[command.index("-r") + 1] for command in commands] == ["chr22", "22"]
