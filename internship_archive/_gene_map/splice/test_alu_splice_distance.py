import numpy as np
import pandas as pd
import pytest

from alu_splice_distance import (
    classify,
    exon_context,
    interval_exon_context,
    load_per_copy_csv,
    load_rmsk_sites,
    merge_sites,
    single_strand_gene,
    summarise,
    window_test,
    union_intervals,
    zone,
)

EX = np.array([[1000, 1100], [5000, 5200]], dtype=np.int64)


def test_union_intervals_merges_touching_and_overlapping_exons():
    a = np.array([[10, 20], [20, 30], [25, 40], [100, 110]], dtype=np.int64)
    assert union_intervals(a).tolist() == [[10, 40], [100, 110]]


def test_union_intervals_on_empty_input_returns_empty():
    assert len(union_intervals(np.empty((0, 2), dtype=np.int64))) == 0


def test_position_inside_an_exon_is_flagged_and_has_zero_distance():
    in_exon, d, side = exon_context(EX, 1050, "+")
    assert in_exon and d == 0 and side is None


def test_left_of_an_exon_is_upstream_on_a_plus_strand_gene():
    in_exon, d, side = exon_context(EX, 900, "+")
    assert not in_exon and d == 100 and side == "upstream"


def test_left_of_an_exon_is_downstream_on_a_minus_strand_gene():
    _, d, side = exon_context(EX, 900, "-")
    assert d == 100 and side == "downstream"


def test_right_of_an_exon_is_downstream_on_a_plus_strand_gene():
    _, d, side = exon_context(EX, 1200, "+")
    assert d == 101 and side == "downstream"


def test_right_of_an_exon_is_upstream_on_a_minus_strand_gene():
    _, d, side = exon_context(EX, 1200, "-")
    assert d == 101 and side == "upstream"


def test_nearest_exon_is_chosen_when_the_position_sits_between_two():
    # 2000 is 900 past exon one and 3000 short of exon two
    _, d, side = exon_context(EX, 2000, "+")
    assert d == 901 and side == "downstream"


def test_position_in_no_gene_returns_none():
    assert single_strand_gene([(0, 100, "+")], 500) is None


def test_position_in_two_genes_on_opposite_strands_returns_none():
    assert single_strand_gene([(0, 1000, "+"), (500, 1500, "-")], 700) is None


def test_two_genes_on_the_same_strand_still_resolve():
    assert single_strand_gene([(0, 1000, "+"), (500, 1500, "+")], 700) == "+"


def test_merge_sites_collapses_nearby_calls_into_one_site():
    calls = pd.DataFrame(
        {
            "chrom": ["chr1"] * 3,
            "family": ["ALU"] * 3,
            "pos": [1000, 1020, 1040],
            "orient": ["+", "+", "-"],
            "sample": list("abc"),
        }
    )
    out = merge_sites(calls, merge_bp=50)
    assert len(out) == 1
    assert out.iloc[0].n_calls == 3
    assert out.iloc[0].orient == "+"  # majority


def test_merge_sites_keeps_distant_calls_separate():
    calls = pd.DataFrame(
        {
            "chrom": ["chr1"] * 2,
            "family": ["ALU"] * 2,
            "pos": [1000, 5000],
            "orient": ["+", "+"],
            "sample": list("ab"),
        }
    )
    assert len(merge_sites(calls, merge_bp=50)) == 2


def test_a_tied_orientation_vote_is_left_unresolved():
    calls = pd.DataFrame(
        {
            "chrom": ["chr1"] * 2,
            "family": ["ALU"] * 2,
            "pos": [1000, 1010],
            "orient": ["+", "-"],
            "sample": list("ab"),
        }
    )
    assert merge_sites(calls, merge_bp=50).iloc[0].orient == ""


def test_sense_is_orientation_matching_the_host_gene_strand():
    sites = pd.DataFrame(
        [{"chrom": "chr1", "family": "ALU", "pos": 900, "orient": "+", "n_calls": 1}]
    )
    out = classify(sites, {"chr1": [(0, 9000, "+")]}, {"chr1": EX})
    assert out.iloc[0].sense and out.iloc[0].side == "upstream"


def test_an_insertion_with_no_orientation_call_is_dropped():
    sites = pd.DataFrame(
        [{"chrom": "chr1", "family": "ALU", "pos": 900, "orient": "", "n_calls": 1}]
    )
    assert len(classify(sites, {"chr1": [(0, 9000, "+")]}, {"chr1": EX})) == 0


def test_zone_labels_split_near_and_far_intronic_positions():
    near = pd.Series({"in_exon": False, "dist_to_exon": 100, "side": "upstream"})
    far = pd.Series({"in_exon": False, "dist_to_exon": 5000, "side": "upstream"})
    exonic = pd.Series({"in_exon": True, "dist_to_exon": 0, "side": None})
    assert zone(near) == "within_150_upstream"
    assert zone(far) == "intronic_far"
    assert zone(exonic) == "exonic"


def _frame(n_near_anti, n_near_sense, n_far_anti, n_far_sense):
    rows = []
    for _ in range(n_near_anti):
        rows.append({"in_exon": False, "dist_to_exon": 50, "side": "upstream", "sense": False})
    for _ in range(n_near_sense):
        rows.append({"in_exon": False, "dist_to_exon": 50, "side": "upstream", "sense": True})
    for _ in range(n_far_anti):
        rows.append({"in_exon": False, "dist_to_exon": 9000, "side": "upstream", "sense": False})
    for _ in range(n_far_sense):
        rows.append({"in_exon": False, "dist_to_exon": 9000, "side": "upstream", "sense": True})
    return pd.DataFrame(rows)


def test_window_detects_antisense_depletion_near_exons():
    r = window_test(_frame(10, 90, 550, 450), "upstream")
    assert r["antisense_frac_near"] == pytest.approx(0.10)
    assert r["odds_ratio"] < 1 and r["p"] < 1e-6


def test_window_reports_no_effect_when_fractions_match():
    r = window_test(_frame(50, 50, 500, 500), "upstream")
    assert r["odds_ratio"] == pytest.approx(1.0)
    assert r["p"] > 0.9


def test_summarise_counts_every_site_exactly_once():
    df = _frame(3, 4, 5, 6)
    assert summarise(df).n.sum() == len(df)


def test_load_per_copy_csv_maps_midpoint_and_orientation(tmp_path):
    csv_path = tmp_path / "per_copy.csv"
    pd.DataFrame(
        [
            {"chrom": "chr1", "start": 1000, "end": 1300, "elem_strand": "+"},
            {"chrom": "chr2", "start": 2000, "end": 2301, "elem_strand": "-"},
        ]
    ).to_csv(csv_path, index=False)
    out = load_per_copy_csv(str(csv_path))
    assert out["pos"].tolist() == [1151, 2151]
    assert out["orient"].tolist() == ["+", "-"]
    assert out["family"].tolist() == ["ALU", "ALU"]


def test_load_rmsk_sites_filters_by_exact_and_prefix_rep_name(tmp_path):
    rmsk_path = tmp_path / "rmsk.txt"
    lines = [
        "\t".join(["585", "2000", "50", "0", "0", "chr1", "1000", "1300", "-1", "+", "AluY", "SINE", "Alu", "1", "300", "0", "1"]),
        "\t".join(["585", "2000", "90", "0", "0", "chr1", "3000", "3300", "-1", "-", "AluSx", "SINE", "Alu", "0", "300", "1", "2"]),
        "\t".join(["585", "2000", "50", "0", "0", "chrUn_x", "5000", "5300", "-1", "+", "AluY", "SINE", "Alu", "1", "300", "0", "3"]),
    ]
    rmsk_path.write_text("\n".join(lines) + "\n")
    exact = load_rmsk_sites(str(rmsk_path), "AluY")
    assert len(exact) == 1
    assert exact.iloc[0]["pos"] == 1151
    assert exact.iloc[0]["start"] == 1000
    assert exact.iloc[0]["end"] == 1300
    assert exact.iloc[0]["orient"] == "+"

    prefix = load_rmsk_sites(str(rmsk_path), "Alu*")
    assert len(prefix) == 2
    assert prefix["family"].tolist() == ["AluY", "AluSx"]


def test_interval_exon_context_flags_junction_spanning_copy_as_exonic():
    # Copy [850, 1050) overlaps exon [1000, 1100), even though its midpoint 950 is outside
    in_ex_mid, d_mid, _ = exon_context(EX, 950, "+")
    assert not in_ex_mid and d_mid == 50
    in_ex_iv, d_iv, side_iv = interval_exon_context(EX, 850, 1050, "+")
    assert in_ex_iv and d_iv == 0 and side_iv is None


def test_interval_exon_context_measures_edge_distance_upstream_and_downstream():
    # Copy [600, 900) ends 100 bp left of exon [1000, 1100)
    in_ex, d, side = interval_exon_context(EX, 600, 900, "+")
    assert not in_ex and d == 100 and side == "upstream"
    in_ex_m, d_m, side_m = interval_exon_context(EX, 600, 900, "-")
    assert not in_ex_m and d_m == 100 and side_m == "downstream"
    # Copy [1200, 1500) starts 100 bp right of exon [1000, 1100)
    in_ex2, d2, side2 = interval_exon_context(EX, 1200, 1500, "+")
    assert not in_ex2 and d2 == 100 and side2 == "downstream"


def test_upstream_depletion_replicates_with_fixed_reference_copies():
    import pandas as pd
    import numpy as np
    from scipy.stats import fisher_exact

    sites = pd.read_csv(
        "out_aluy_ref_repro/alu_splice_sites.tsv",
        sep="\t",
    )
    # Re-derive zones exactly as the pipeline's classify() does for reference
    # copies: exonized copies are exonic, side-tagged copies are within_150_*,
    # and everything else is intronic_far.
    def reclassify(r):
        if r.in_exon:
            return "exonic"
        if pd.notna(r.side):
            return f"within_150_{r.side}"
        return "intronic_far"

    sites["zone"] = sites.apply(reclassify, axis=1)

    up = sites[sites.zone == "within_150_upstream"]
    far = sites[sites.zone == "intronic_far"]
    tbl = np.array([
        [(up.orient != up.gene_strand).sum(), (up.orient == up.gene_strand).sum()],
        [(far.orient != far.gene_strand).sum(), (far.orient == far.gene_strand).sum()],
    ])
    or_, p = fisher_exact(tbl)
    # Reproduces the AluY upstream 150 bp antisymmetric-depletion window vs far
    # intronic: OR(anti) ~ 0.758 and p ~ 9.7e-08.
    assert or_ < 0.85 and p < 5e-7


