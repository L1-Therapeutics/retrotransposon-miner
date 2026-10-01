"""Gold-table training labels from one genome-wide ranked list per sample."""

from __future__ import annotations

import pandas as pd

from retro_miner.gold_classifier import (
    CATEGORICAL_FEATURES,
    NUMERIC_FEATURES,
    bernoulli_changepoint,
    label_gold_table,
    rank_gold_like_review,
    sample_matched_negatives,
)


def _support(mei_mapped: int) -> str:
    return f"SR_L=1,SR_R=1,DPE_L=0,DPE_R=0,MEI_MAPPED={mei_mapped},polyA_MAPPED=0,VNTR_MAPPED=0"


def _gold_row(family: str, known: bool, mei_mapped: int, chrom: str) -> dict:
    return {
        "chrom": chrom,
        "window_start": 1,
        "window_end": 2,
        "consensus_insertion_breakpoint_pos": 10,
        "consensus_mei_family": family,
        "known_mei_polymorphism": known,
        "analysis_stage_tier": "gold",
        "consensus_poly_at_min_bp": 10,
        "consensus_tsd_len_estimate": 8,
        "split_cluster_reads": 4,
        "split_cluster_window_reads": 8,
        "disease_left_flank_mei_reads": 3,
        "disease_right_flank_mei_reads": 2,
        "disease_supporting_reads": _support(mei_mapped),
        "control_supporting_reads": _support(mei_mapped),
        "read_support_heuristic_score": mei_mapped / 10.0,
        "insertion_model_score": 0.5,
    }


def test_changepoint_keeps_the_dense_prefix() -> None:
    known = [True] * 8 + [False] * 40 + [True]
    assert bernoulli_changepoint(known) == 8
    assert bernoulli_changepoint([False, False, False]) == 0


def test_genome_wide_rank_follows_support_not_chromosome_or_catalog() -> None:
    # Weaker catalog hit on chr1 is listed first. Review priority should still
    # place the stronger non-catalog chr22 call above it.
    frame = pd.DataFrame(
        [
            _gold_row("Alu", True, mei_mapped=2, chrom="chr1"),
            _gold_row("Alu", False, mei_mapped=40, chrom="chr22"),
        ]
    )
    ranked = rank_gold_like_review(frame)
    assert ranked["chrom"].tolist() == ["chr22", "chr1"]
    assert ranked["known_mei_polymorphism"].tolist() == [False, True]


def test_bottom_known_is_excluded_and_negatives_match_family() -> None:
    rows = [_gold_row("Alu", True, mei_mapped=50 - i, chrom="chr1") for i in range(3)]
    rows.extend(_gold_row("Alu", False, mei_mapped=20, chrom="chr2") for _ in range(4))
    rows.append(_gold_row("Alu", True, mei_mapped=1, chrom="chr22"))
    rows.extend(_gold_row("L1", False, mei_mapped=3, chrom="chr3") for _ in range(3))
    labeled = label_gold_table(rank_gold_like_review(pd.DataFrame(rows)))
    labeled["sample"] = "HG00100"
    positives = labeled.loc[labeled["train_role"].eq("positive")]
    assert len(positives) == 3
    assert set(positives["consensus_mei_family"]) == {"Alu"}
    outliers = labeled.loc[labeled["train_role"].eq("excluded_outlier_known")]
    assert len(outliers) == 1
    assert outliers["chrom"].tolist() == ["chr22"]
    training = sample_matched_negatives(labeled, seed=1)
    alu = training.loc[training["mei_family"].eq("Alu")]
    assert int(alu["train_label"].sum()) == int((alu["train_label"] == 0).sum())
    assert "excluded_outlier_known" not in set(training["train_role"])


def test_in_cluster_non_known_stays_unlabeled() -> None:
    known = [True, False, True] + [False] * 30
    labeled = label_gold_table(pd.DataFrame({"known_mei_polymorphism": known, "consensus_mei_family": "Alu"}))
    in_cluster = labeled.loc[labeled["in_top_cluster"]]
    non_known = in_cluster.loc[~in_cluster["known_mei_polymorphism"]]
    assert len(non_known) == 1
    assert non_known["train_role"].tolist() == ["unlabeled"]


def test_rank_coordinates_and_catalog_fields_are_not_features() -> None:
    banned = {
        "chrom",
        "window_start",
        "consensus_insertion_breakpoint_pos",
        "known_mei_polymorphism",
        "known_mei_polymorphism_id",
        "read_support_heuristic_score",
        "insertion_model_score",
    }
    used = set(NUMERIC_FEATURES) | set(CATEGORICAL_FEATURES)
    assert banned.isdisjoint(used)
