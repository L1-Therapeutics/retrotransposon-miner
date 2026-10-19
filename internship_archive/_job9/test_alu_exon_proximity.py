import numpy as np
import pandas as pd
import pytest

from alu_exon_proximity import (
    build,
    covered,
    merge_intervals,
    nearest_exon,
    summarise,
    window_test,
    zone_of,
)

EX = np.array([[1000, 1100], [5000, 5200]], dtype=np.int64)


def test_merge_intervals_sorts_then_merges():
    a = np.array([[100, 110], [10, 20], [20, 30]], dtype=np.int64)
    assert merge_intervals(a).tolist() == [[10, 30], [100, 110]]


def test_covered_detects_membership_and_gaps():
    iv = np.array([[10, 20], [30, 40]], dtype=np.int64)
    got = covered(iv, np.array([5, 10, 19, 20, 35, 45]))
    assert got.tolist() == [False, True, True, False, True, False]


def test_covered_on_empty_intervals_is_all_false():
    assert covered(np.empty((0, 2), np.int64), np.array([1, 2])).tolist() == [False, False]


def test_element_entirely_left_of_an_exon_measures_the_gap():
    d, s = nearest_exon(EX, np.array([800]), np.array([900]))
    assert d.tolist() == [100] and s.tolist() == [-1]


def test_element_entirely_right_of_an_exon_measures_the_gap():
    d, s = nearest_exon(EX, np.array([1200]), np.array([1300]))
    assert d.tolist() == [100] and s.tolist() == [1]


def test_element_overlapping_an_exon_has_distance_zero():
    d, s = nearest_exon(EX, np.array([1050]), np.array([1150]))
    assert d.tolist() == [0] and s.tolist() == [0]


def test_element_containing_an_exon_has_distance_zero():
    d, s = nearest_exon(EX, np.array([900]), np.array([1200]))
    assert d.tolist() == [0] and s.tolist() == [0]


def test_element_between_two_exons_takes_the_closer_one():
    # 1150-1250 is 50 past exon one, 3750 short of exon two
    d, s = nearest_exon(EX, np.array([1150]), np.array([1250]))
    assert d.tolist() == [50] and s.tolist() == [1]


def test_nearest_exon_with_no_exons_returns_the_sentinel():
    d, s = nearest_exon(np.empty((0, 2), np.int64), np.array([100]), np.array([200]))
    assert d.tolist() == [-1] and s.tolist() == [0]


def _repeats(start, end, strand):
    return pd.DataFrame(
        {"chrom": ["chr1"], "start": [start], "end": [end], "strand": [strand], "milli_div": [50]}
    )


PLUS = {"chr1": np.array([[0, 9000]], dtype=np.int64)}
MINUS = {"chr1": np.array([[0, 9000]], dtype=np.int64)}
EXD = {"chr1": EX}


def test_left_of_exon_is_upstream_on_a_plus_strand_gene():
    out = build(_repeats(800, 900, "+"), PLUS, {}, EXD)
    assert out.iloc[0].side == "upstream" and out.iloc[0].sense


def test_left_of_exon_is_downstream_on_a_minus_strand_gene():
    out = build(_repeats(800, 900, "+"), {}, MINUS, EXD)
    assert out.iloc[0].side == "downstream"


def test_right_of_exon_is_upstream_on_a_minus_strand_gene():
    out = build(_repeats(1200, 1300, "+"), {}, MINUS, EXD)
    assert out.iloc[0].side == "upstream"


def test_sense_is_element_strand_matching_gene_strand():
    assert not build(_repeats(800, 900, "-"), PLUS, {}, EXD).iloc[0].sense
    assert build(_repeats(800, 900, "-"), {}, MINUS, EXD).iloc[0].sense


def test_a_copy_inside_genes_on_both_strands_is_excluded():
    assert len(build(_repeats(800, 900, "+"), PLUS, MINUS, EXD)) == 0


def test_a_copy_outside_any_gene_is_excluded():
    plus = {"chr1": np.array([[20000, 30000]], dtype=np.int64)}
    assert len(build(_repeats(800, 900, "+"), plus, {}, EXD)) == 0


def _zoned(dist, side, sense):
    return pd.DataFrame({"dist_to_exon": dist, "side": side, "sense": sense})


def test_zone_labels_separate_overlap_near_and_far():
    df = _zoned([0, 100, 100, 9000], ["overlap", "upstream", "downstream", "upstream"], [True] * 4)
    assert zone_of(df).tolist() == [
        "overlaps_exon",
        "within_150_upstream",
        "within_150_downstream",
        "intronic_far",
    ]


def test_window_test_recovers_a_planted_depletion():
    df = _zoned(
        [50] * 100 + [9000] * 1000,
        ["upstream"] * 1100,
        [False] * 10 + [True] * 90 + [False] * 550 + [True] * 450,
    )
    r = window_test(df, "upstream")
    assert r["antisense_frac_near"] == pytest.approx(0.10)
    assert r["odds_ratio"] < 1 and r["p"] < 1e-6


def test_window_test_is_null_when_the_fractions_match():
    df = _zoned([50] * 100 + [9000] * 100, ["upstream"] * 200, ([False] * 50 + [True] * 50) * 2)
    r = window_test(df, "upstream")
    assert r["odds_ratio"] == pytest.approx(1.0) and r["p"] > 0.9


def test_summarise_accounts_for_every_copy():
    df = _zoned([0, 50, 9000], ["overlap", "upstream", "upstream"], [True, False, True])
    df["milli_div"] = [10, 20, 30]
    assert summarise(df).n.sum() == 3
