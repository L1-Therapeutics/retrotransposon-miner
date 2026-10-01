"""Tests for Phase 2b (scripts/score_tprt_hallmarks.py).

The load is placed on the failure modes this phase actually hit while being
written, rather than on line coverage:

  * `unevaluable` must never be laundered into `absent`;
  * the permutation null must preserve each stratum's margins -- the first
    implementation subtracted the permuted event count from the *observed*
    reference count instead of the stratum total, which silently distorted
    every p-value in this phase;
  * the Mantel-Haenszel risk difference must equal the crude difference when
    there is only one stratum, which it did not until the N_i weight was added;
  * matching must never read a hallmark column;
  * a sentinel TSD length must never enter a TSD-length statement.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import random
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "score_tprt_hallmarks", REPO_ROOT / "scripts" / "score_tprt_hallmarks.py"
)
tprt = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(tprt)

SENTINELS = frozenset({41, 82})


# --------------------------------------------------------------------------
# fakes
# --------------------------------------------------------------------------


class FakeReference:
    """Minimal stand-in for the pyfaidx wrapper: one contig, 0-based half-open."""

    def __init__(self, chrom: str, seq: str) -> None:
        self.chrom = chrom
        self.seq = seq

    def fetch(self, chrom: str, start0: int, end0: int) -> str | None:
        if chrom != self.chrom or start0 < 0 or end0 <= start0 or end0 > len(self.seq):
            return None
        return self.seq[start0:end0]


def make_row(**overrides) -> dict[str, str]:
    row = {
        "chrom": "chr1",
        "pos": "1000",
        "tsd_len": "",
        "tsd_seq": "",
        "polya_len": "",
        "polya_seq": "",
        "insert_strand": "+",
        "host_name": "AluY",
        "host_len": "300",
        "consensus_offset": "133",
        "consensus_mapping": "unambiguous",
        "nested_in_alu_host": "1",
        "conformation": "FOR+POLYA",
        "site_allele_freq": "0.0005",
        "perc_resolved": "99.9",
        "alignment_identity": "0.9",
        "offset_drift_bp": "0",
        "rt_len": "281",
        "host_strand": "+",
    }
    row.update(overrides)
    return row


# --------------------------------------------------------------------------
# primitive helpers
# --------------------------------------------------------------------------


def test_hamming_counts_substitutions_and_penalises_length_mismatch():
    assert tprt._hamming("ACGT", "ACGT") == 0
    assert tprt._hamming("ACGT", "ACGA") == 1
    # unequal lengths cannot be aligned base-for-base; the penalty is the
    # longer length so the caller never reads a partial match as a good one.
    assert tprt._hamming("ACG", "ACGT") == 4


def test_revcomp_round_trips():
    assert tprt.revcomp("ACGT") == "ACGT"
    assert tprt.revcomp("AAAG") == "CTTT"
    assert tprt.revcomp(tprt.revcomp("GGCCGGGCGCGGTGG")) == "GGCCGGGCGCGGTGG"


# --------------------------------------------------------------------------
# TSD: unevaluable is not absent
# --------------------------------------------------------------------------


def test_missing_tsd_is_unevaluable_not_absent():
    out = tprt.tsd_reference_check(make_row(), FakeReference("chr1", "A" * 2000), SENTINELS)
    assert out["state"] == "unevaluable"
    assert out["reason"] == "no_caller_tsd"


@pytest.mark.parametrize("length", [41, 82])
def test_sentinel_tsd_length_is_unevaluable_not_inconsistent(length):
    seq = "A" * length
    ref = FakeReference("chr1", "N" * 999 + seq + "N" * 1500)
    out = tprt.tsd_reference_check(
        make_row(tsd_len=str(length), tsd_seq=seq), ref, SENTINELS
    )
    assert out["state"] == "unevaluable"
    assert out["reason"] == "sentinel_tsd_length"


def test_tsd_len_disagreeing_with_tsd_seq_is_unevaluable():
    out = tprt.tsd_reference_check(
        make_row(tsd_len="12", tsd_seq="AAAA"), FakeReference("chr1", "A" * 2000), SENTINELS
    )
    assert out["state"] == "unevaluable"
    assert out["reason"] == "tsd_seq_length_disagrees_with_tsd_len"


def test_unreadable_locus_is_unevaluable():
    out = tprt.tsd_reference_check(
        make_row(tsd_len="4", tsd_seq="AAAA"), FakeReference("chr2", "A" * 2000), SENTINELS
    )
    assert out["state"] == "unevaluable"
    assert out["reason"] == "insertion_point_unreadable"


def test_exact_tsd_match_is_reference_verifiable():
    seq = "GTCAG"
    ref = FakeReference("chr1", "N" * 999 + seq + "N" * 1500)
    out = tprt.tsd_reference_check(
        make_row(tsd_len="5", tsd_seq=seq), ref, SENTINELS
    )
    assert out["state"] == "reference_verifiable"
    assert out["exact"] is True
    assert out["n_mismatches"] == 0


def test_tsd_within_tolerance_is_verifiable_but_not_exact():
    seq = "GTCAG"
    ref = FakeReference("chr1", "N" * 999 + "GTCAT" + "N" * 1500)
    out = tprt.tsd_reference_check(make_row(tsd_len="5", tsd_seq=seq), ref, SENTINELS)
    assert out["state"] == "reference_verifiable"
    assert out["exact"] is False
    assert out["n_mismatches"] == 1


def test_tsd_beyond_tolerance_is_inconsistent_not_unevaluable():
    seq = "GTCAG"
    ref = FakeReference("chr1", "N" * 999 + "AAAAA" + "N" * 1500)
    out = tprt.tsd_reference_check(make_row(tsd_len="5", tsd_seq=seq), ref, SENTINELS)
    assert out["state"] == "reference_inconsistent"
    assert out["n_mismatches"] == 4  # GTCAG vs AAAAA differs at four positions


def test_mismatch_offsets_are_measured_from_the_three_prime_end():
    """The terminal-mismatch concentration is the diagnostic, so pin the axis."""
    seq = "AAAAAGGGG"  # 9 bases; the final base differs from the reference
    ref = FakeReference("chr1", "N" * 999 + "AAAAAGGGT" + "N" * 1500)
    out = tprt.tsd_reference_check(make_row(tsd_len="9", tsd_seq=seq), ref, SENTINELS)
    assert out["mismatch_offsets_from_3prime_end"] == [0]


def test_tsd_without_a_reference_is_unevaluable():
    out = tprt.tsd_reference_check(
        make_row(tsd_len="4", tsd_seq="AAAA"), None, SENTINELS
    )
    assert out["state"] == "unevaluable"
    assert out["reason"] == "no_reference_available"


def test_tsd_length_class_is_categorical_and_separates_sentinels():
    assert tprt.tsd_length_class(12, SENTINELS) == "de_novo_range_5_27"
    assert tprt.tsd_length_class(40, SENTINELS) == "longer_than_de_novo_28_plus"
    assert tprt.tsd_length_class(41, SENTINELS) == "sentinel_pileup"
    assert tprt.tsd_length_class(None, SENTINELS) == "unevaluable_no_length"
    # Gate 0d's range endpoints are inclusive and nothing else leaks in.
    assert tprt.tsd_length_class(5, SENTINELS) == "de_novo_range_5_27"
    assert tprt.tsd_length_class(27, SENTINELS) == "de_novo_range_5_27"
    assert tprt.tsd_length_class(28, SENTINELS) == "longer_than_de_novo_28_plus"


def test_tsd_structure_excludes_sentinels_from_the_length_distribution():
    rows = [{"tsd_len": v, "offset_drift_bp": "0"} for v in (10, 12, 41, 82, 20)]
    out = tprt.tsd_structure(rows, SENTINELS)
    assert out["with_sentinels"]["n"] == 5
    assert out["sentinel_excluded"]["n"] == 3
    assert out["sentinel_excluded"]["max"] == 20
    assert out["n_sentinel_rows_excluded"] == 2


# --------------------------------------------------------------------------
# poly(A) and the L1 endonuclease pentamer
# --------------------------------------------------------------------------


def test_missing_polya_field_is_unevaluable_not_absent():
    out = tprt.polya_state(make_row())
    assert out["state"] == "unevaluable"
    assert out["reason"] == "no_caller_polya_field"


def test_polya_length_disagreeing_with_sequence_is_unevaluable():
    out = tprt.polya_state(make_row(polya_seq="AAAA", polya_len=""))
    assert out["state"] == "unevaluable"


def test_short_polya_is_absent_and_long_polya_is_present():
    assert tprt.polya_state(make_row(polya_seq="AAAA", polya_len="4"))["state"] == "absent"
    assert tprt.polya_state(make_row(polya_seq="A" * 20, polya_len="20"))["state"] == "present"


def test_l1_en_pentamer_is_reported_on_both_sides():
    # cut = pos - 1 = 999; the 5' side is ref[994:999]
    ref = FakeReference("chr1", "N" * 994 + "TTTTA" + "N" * 1500)
    out = tprt.l1_en_nick_state(make_row(), ref)
    assert out["state"] == "tttt_a_five_prime_only"
    assert out["five_prime_side"] is True
    assert out["three_prime_side"] is False


def test_l1_en_pentamer_on_the_three_prime_side():
    ref = FakeReference("chr1", "N" * 999 + "TTTTA" + "N" * 1500)
    out = tprt.l1_en_nick_state(make_row(), ref)
    assert out["state"] == "tttt_a_three_prime_only"
    assert out["five_prime_side"] is False


def test_l1_en_absent_is_not_unevaluable():
    ref = FakeReference("chr1", "G" * 2000)
    assert tprt.l1_en_nick_state(make_row(), ref)["state"] == "no_tttt_a_at_nick"


def test_l1_en_without_a_reference_is_unevaluable():
    assert tprt.l1_en_nick_state(make_row(), None)["state"] == "unevaluable"


# --------------------------------------------------------------------------
# EN motif strand compatibility
# --------------------------------------------------------------------------


MOTIF = "GGCCGGGCGCGGTGG"
#: cut = pos - 1 = 999. The search window is ref[cut - flank - len, cut + flank + len),
#: so a left-side motif ending exactly at the nick occupies ref[999 - len : 999].
MOTIF_LEFT_START = 999 - len(MOTIF)


def test_plus_insert_expects_the_motif_on_the_left_of_the_nick():
    ref = FakeReference("chr1", "N" * MOTIF_LEFT_START + MOTIF + "N" * 1500)
    out = tprt.en_motif_state(make_row(insert_strand="+"), ref, MOTIF)
    assert out["expected_side"] == "left"
    assert out["state"] == "motif_compatible"


def test_minus_insert_expects_the_reverse_complement_on_the_right():
    ref = FakeReference("chr1", "N" * 999 + tprt.revcomp(MOTIF) + "N" * 1500)
    out = tprt.en_motif_state(make_row(insert_strand="-"), ref, MOTIF)
    assert out["expected_side"] == "right"
    assert out["state"] == "motif_compatible"


def test_motif_on_the_wrong_side_is_incompatible_not_absent():
    """The distinction the whole strand-compatibility rule exists to make."""
    ref = FakeReference("chr1", "N" * MOTIF_LEFT_START + MOTIF + "N" * 1500)
    out = tprt.en_motif_state(make_row(insert_strand="-"), ref, MOTIF)
    assert out["state"] == "motif_incompatible"
    assert out["hit_side"] == "left"


def test_en_motif_tolerates_a_single_mismatch():
    one_off = MOTIF[:-1] + "A"
    ref = FakeReference("chr1", "N" * MOTIF_LEFT_START + one_off + "N" * 1500)
    out = tprt.en_motif_state(make_row(insert_strand="+"), ref, MOTIF)
    assert out["state"] == "motif_compatible"
    assert out["mismatches"] == 1


def test_unresolved_insert_strand_makes_the_en_probe_unevaluable():
    ref = FakeReference("chr1", "G" * 2000)
    out = tprt.en_motif_state(make_row(insert_strand="."), ref, MOTIF)
    assert out["state"] == "unevaluable"
    assert out["reason"] == "insert_strand_unresolved"


# --------------------------------------------------------------------------
# triad assembly
# --------------------------------------------------------------------------


def _components(tsd_state, polya_state, en_state):
    return {
        "tsd": {"state": tsd_state},
        "polya": {"state": polya_state},
        "en": {"state": en_state},
    }


def test_triad_counts_unevaluable_in_neither_numerator_nor_denominator():
    out = tprt.score_triad(
        _components("reference_verifiable", "unevaluable", "motif_compatible")
    )
    assert out["n_evaluable"] == 2
    assert out["n_observed"] == 2
    assert out["triad_complete"] is False
    assert out["hallmarks_per_evaluable"] == 1.0


def test_triad_complete_requires_all_three_evaluable():
    out = tprt.score_triad(
        _components("reference_verifiable", "present", "motif_absent")
    )
    assert out["triad_complete"] is True
    assert out["triad_complete_observed"] is False
    assert out["n_observed"] == 2


def test_triad_never_scores_a_refuted_tsd_as_observed():
    out = tprt.score_triad(
        _components("reference_inconsistent", "present", "motif_compatible")
    )
    assert out["n_evaluable"] == 3
    assert out["n_observed"] == 2


def test_fully_unevaluable_triad_has_no_rate():
    out = tprt.score_triad(
        _components("unevaluable", "unevaluable", "unevaluable")
    )
    assert out["n_evaluable"] == 0
    assert out["hallmarks_per_evaluable"] is None


# --------------------------------------------------------------------------
# Mantel-Haenszel weighting
# --------------------------------------------------------------------------


def test_single_stratum_risk_difference_equals_the_crude_difference():
    tables = [(30, 70, 10, 90)]
    assert tprt._pooled_risk_difference(tables) == pytest.approx(0.30 - 0.10)


def test_equal_strata_average_the_stratum_differences():
    # Two strata, both 100 rows, stratum risks 0.4 vs 0.1 and 0.0 vs 0.5.
    tables = [(40, 60, 10, 90), (0, 100, 50, 50)]
    assert tprt._pooled_risk_difference(tables) == pytest.approx((0.3 + -0.5) / 2)


def test_single_stratum_odds_ratio_is_the_ordinary_one():
    assert tprt._pooled_odds_ratio([(30, 70, 10, 90)]) == pytest.approx((30 * 90) / (70 * 10))


def test_strata_missing_an_arm_are_skipped():
    assert tprt._pooled_risk_difference([(5, 5, 0, 0), (10, 90, 20, 80)]) == pytest.approx(
        0.10 - 0.20
    )
    assert tprt._pooled_risk_difference([(5, 5, 0, 0)]) is None


# --------------------------------------------------------------------------
# permutation null
# --------------------------------------------------------------------------


def _strata_from_flags(flags: list[tuple[int, int, int, int]]) -> list[dict]:
    """Build strata from (a, b, c, d) tables so a test can state an exact null."""
    strata = []
    for i, (a, b, c, d) in enumerate(flags):
        exposed = [{"v": True}] * a + [{"v": False}] * b
        reference = [{"v": True}] * c + [{"v": False}] * d
        strata.append({"key": [str(i)], "exposed": exposed, "reference": reference})
    return strata


def test_permutation_preserves_each_stratum_margin():
    """Regression: the null used to subtract from the observed reference count.

    Every replicate must leave the stratum's event totals exactly where they
    started, or the null is not the null.
    """
    strata = _strata_from_flags([(3, 5, 20, 80), (1, 4, 2, 8)])
    res = tprt.stratified_tests(
        strata, {"v": lambda r: r["v"]}, permutations=50, seed=3
    )["v"]
    for table in res["strata"]:
        assert table["a"] + table["b"] + table["c"] + table["d"] == (
            108 if table["key"] == ["0"] else 15
        )


def test_permutation_preserves_group_sizes_with_unevaluable_rows():
    strata = [
        {
            "key": ["0"],
            "exposed": [{"v": True}, {"v": None}, {"v": False}],
            "reference": [{"v": True}, {"v": None}, {"v": True}, {"v": None}],
        }
    ]
    res = tprt.stratified_tests(strata, {"v": lambda r: r["v"]}, permutations=10, seed=1)["v"]
    obs = res["observed"]
    # unevaluable rows are excluded from both numerators and both denominators
    assert obs["exposed_events"] == 1 and obs["exposed_nonevents"] == 1
    assert obs["reference_events"] == 2 and obs["reference_nonevents"] == 0
    assert obs["exposed_unevaluable"] == 1
    assert obs["reference_unevaluable"] == 2
    assert obs["exposed_rate"] == pytest.approx(0.5)
    assert obs["reference_rate"] == pytest.approx(1.0)


def test_no_difference_gives_a_null_centred_on_zero():
    """A null shifted away from zero would make every p-value meaningless."""
    rng = random.Random(11)
    strata = []
    for s in range(6):
        exposed = [{"v": rng.random() < 0.4} for _ in range(10)]
        reference = [{"v": rng.random() < 0.4} for _ in range(40)]
        strata.append({"key": [str(s)], "exposed": exposed, "reference": reference})
    res = tprt.stratified_tests(strata, {"v": lambda r: r["v"]}, permutations=2000, seed=9)
    assert res["v"]["permutation_p"] > 0.05
    # the estimate itself should sit near the middle of the null
    assert abs(res["v"]["observed"]["risk_difference"]) < 0.2


def test_permutation_p_matches_an_exhaustive_enumeration():
    """Monte-Carlo inference has to agree with the exact answer it approximates."""
    rng = random.Random(7)
    strata = []
    for s, (ne, nr, p) in enumerate([(2, 3, 0.5), (1, 2, 0.0)]):
        strata.append(
            {
                "key": [str(s)],
                "exposed": [{"v": rng.random() < p} for _ in range(ne)],
                "reference": [{"v": rng.random() < p} for _ in range(nr)],
            }
        )
    event = lambda r: r["v"]  # noqa: E731
    res = tprt.stratified_tests(strata, {"v": event}, permutations=5000, seed=5)["v"]
    observed = res["observed"]["risk_difference"]

    exact = []
    for combo in itertools.product(
        *[
            itertools.permutations(range(len(s["exposed"]) + len(s["reference"])))
            for s in strata
        ]
    ):
        tables = []
        for s, perm in zip(strata, combo):
            pool = s["exposed"] + s["reference"]
            ne = len(s["exposed"])
            drawn = [pool[i] for i in perm[:ne]]
            a = sum(1 for x in drawn if x["v"] is True)
            b = sum(1 for x in drawn if x["v"] is False)
            c = sum(1 for x in pool if x["v"] is True) - a
            d = sum(1 for x in pool if x["v"] is False) - b
            tables.append((a, b, c, d))
        value = tprt._pooled_risk_difference(tables)
        if value is not None:
            exact.append(value)
    extreme = sum(1 for v in exact if abs(v) >= abs(observed) - 1e-12)
    exact_p = extreme / len(exact)
    # 5000 Monte-Carlo draws against 720 exact permutations: allow sampling slack.
    assert res["permutation_p"] == pytest.approx(exact_p, abs=0.05)


def test_bh_adjustment_is_monotone_and_bounded():
    ps = [0.01, 0.04, 0.03, None, 0.5]
    adj = tprt._benjamini_hochberg(ps)
    assert adj[3] is None
    finite = [a for a in adj if a is not None]
    assert all(0.0 <= a <= 1.0 for a in finite)
    # BH never decreases p
    for original, adjusted in zip(ps, adj):
        if original is not None:
            assert adjusted >= original - 1e-12
    # ordering by original p is preserved
    assert adj[0] <= adj[2] <= adj[1] <= adj[4]


# --------------------------------------------------------------------------
# matching rules
# --------------------------------------------------------------------------


def test_no_matching_covariate_reads_a_hallmark_column():
    """Matching on an outcome would make the comparison circular."""
    forbidden = ("_tsd", "_polya", "_en", "_l1_en", "_triad", "tsd_", "polya_")
    for comparison, fns in tprt.COVARIATE_SETS.items():
        for name, fn in fns.items():
            source = fn.__doc__ or ""
            assert fn.__code__.co_names is not None
            for const in fn.__code__.co_consts:
                if isinstance(const, str):
                    source += const
            for token in forbidden:
                assert token not in source, f"{comparison}.{name} reads {token}"


def test_ladder_defines_every_level_used_by_the_covariate_sets():
    for comparison, ladder in tprt.LADDERS.items():
        available = set(tprt.COVARIATE_SETS[comparison])
        for covariates in ladder.values():
            assert set(covariates) <= available


def test_coarsening_reduces_or_keeps_the_covariate_set():
    for ladder in tprt.LADDERS.values():
        sizes = [len(v) for v in ladder.values()]
        assert sizes == sorted(sizes, reverse=True)
        assert sizes[-1] == 0


def test_coarsened_match_falls_back_when_the_finest_rung_has_no_support():
    exposed = [make_row(host_name=f"AluZ{i % 40}") for i in range(40)]
    reference = [make_row(host_name=f"AluZ{i % 40}") for i in range(40)]
    match = tprt.coarsened_match(exposed, reference, "peak_vs_rest_of_host")
    # no stratum reaches the reference-support floor, so the ladder must coarsen
    assert match["level"] != "full"
    assert match["exposed_matched"] == 40  # the floor rung keeps everything


def test_coarsened_match_keeps_the_finest_supported_rung():
    exposed = [make_row(host_name="AluY") for _ in range(10)]
    reference = [make_row(host_name="AluY") for _ in range(200)]
    match = tprt.coarsened_match(exposed, reference, "peak_vs_rest_of_host")
    assert match["level"] == "full"
    assert match["exposed_matched"] == 10


def test_coarsened_match_reports_the_rows_it_had_to_drop():
    exposed = [make_row(host_name="AluY") for _ in range(10)] + [
        make_row(host_name="AluRare")
    ]
    reference = [make_row(host_name="AluY") for _ in range(200)]
    match = tprt.coarsened_match(exposed, reference, "peak_vs_rest_of_host")
    assert match["exposed_dropped_unsupported"] == 1
    assert sum(a["exposed_dropped"] for a in match["attempts"]) >= 1


def test_unmatched_comparison_is_flagged_as_unmatched():
    exposed = [make_row(host_name=f"AluZ{i % 40}") for i in range(40)]
    reference = [make_row(host_name=f"AluZ{i % 40}") for i in range(40)]
    out = tprt.run_comparison("peak_vs_rest_of_host", exposed, reference, 1, 20)
    assert out["matched"] is False
    assert "UNMATCHED" in out["interpretation"]


def test_balance_flags_a_covariate_that_never_varies():
    strata = [
        {
            "key": ["AluY", "high"],
            "exposed": [make_row(alignment_identity="0.9") for _ in range(5)],
            "reference": [make_row(alignment_identity="0.9") for _ in range(30)],
        }
    ]
    table = tprt.balance_table(strata, "peak_vs_rest_of_host", "full")
    by_name = {b["covariate"]: b for b in table}
    assert by_name["alignment_confidence"]["max_abs_smd"] == pytest.approx(0.0)
    # A covariate with a single level has nothing to balance; the table says so
    # instead of reporting a confident SMD of zero.
    assert by_name["alignment_confidence"]["levels"][0]["level_value"] == "high"
    assert by_name["alignment_confidence"]["estimable"] is True


# --------------------------------------------------------------------------
# analysis set and covariates
# --------------------------------------------------------------------------


def test_nested_analysis_set_matches_the_phase1_filter():
    rows = [
        make_row(nested_in_alu_host="1", host_len="300", consensus_mapping="unambiguous"),
        make_row(nested_in_alu_host="1", host_len="120", consensus_mapping="unambiguous"),
        make_row(nested_in_alu_host="1", host_len="300", consensus_mapping="ambiguous_indel"),
        make_row(nested_in_alu_host="0", host_len="300", consensus_mapping="unambiguous"),
        make_row(nested_in_alu_host="1", host_len="300", consensus_mapping="unambiguous",
                 consensus_offset=""),
    ]
    assert len(tprt.nested_analysis_set(rows)) == 1


def test_peak_and_control_windows_do_not_overlap():
    peak, control = tprt.split_by_offset(
        [
            make_row(consensus_offset=str(o))
            for o in (127, 128, 133, 138, 139, 159, 160, 259, 260)
        ]
    )
    assert [r["consensus_offset"] for r in peak] == ["128", "133", "138"]
    assert [r["consensus_offset"] for r in control] == ["160", "259"]


def test_confidence_class_is_ordered_and_has_an_unknown_bucket():
    assert tprt._confidence_class(make_row(alignment_identity="0.4")) == "low"
    assert tprt._confidence_class(make_row(alignment_identity="0.6")) == "mid"
    assert tprt._confidence_class(make_row(alignment_identity="0.95")) == "high"
    assert tprt._confidence_class(make_row(alignment_identity="")) == "unknown"


def test_af_and_support_bins_partition_their_ranges():
    assert tprt._af_bin(make_row(site_allele_freq="0.0001")) == "lt_0.001"
    assert tprt._af_bin(make_row(site_allele_freq="0.005")) == "0.001_0.01"
    assert tprt._af_bin(make_row(site_allele_freq="0.05")) == "0.01_0.1"
    assert tprt._af_bin(make_row(site_allele_freq="0.9")) == "ge_0.1"
    assert tprt._support_bin(make_row(perc_resolved="100")) == "ge_99"
    assert tprt._support_bin(make_row(perc_resolved="95")) == "90_99"
    assert tprt._support_bin(make_row(perc_resolved="50")) == "lt_90"


def test_truncation_class_uses_the_conformation_prefix():
    assert tprt._truncation_class(make_row(conformation="TRUN+FOR+POLYA")) == "truncated"
    assert tprt._truncation_class(make_row(conformation="FOR+POLYA")) == "full_or_other"


def test_pure_at_tsd_fraction_ignores_sentinels_and_bad_lengths():
    rows = [
        {"tsd_len": "4", "tsd_seq": "AAAA"},
        {"tsd_len": "4", "tsd_seq": "AAAA"},  # duplicated reference row still counts
        {"tsd_len": "4", "tsd_seq": "GGCG"},
        {"tsd_len": "41", "tsd_seq": "A" * 41},  # sentinel: excluded
        {"tsd_len": "4", "tsd_seq": "AAA"},  # length disagreement: excluded
    ]
    out = tprt.pure_at_tsd_fraction(rows, SENTINELS)
    assert out["n_evaluable"] == 3
    assert out["n_pure_at"] == 2
    assert out["fraction"] == pytest.approx(2 / 3)


def test_en_motif_proxy_is_labelled_a_proxy():
    motif, meta = tprt.derive_en_motif_proxy(tprt.DEFAULT_ALU_CONSENSUS)
    assert len(motif) == tprt.EN_MOTIF_LENGTH_BP
    assert meta["status"] == "consensus_derived_proxy_not_curated_en_consensus"
    assert meta["n_alu_consensus_entries"] > 10
    assert 0.0 < meta["per_column_support_min"] <= 1.0


def test_sentinels_come_from_the_manifest_when_present(tmp_path):
    manifest = tmp_path / "m.json"
    manifest.write_text(json.dumps({"gate_0d_tsd_sentinels": {"sentinel_values": {"41": {}, "82": {}}}}))
    assert tprt.load_sentinels(manifest) == ([41, 82], "manifest")
    assert tprt.load_sentinels(tmp_path / "missing.json") == (
        list(tprt.FALLBACK_SENTINELS),
        "fallback",
    )


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def _synthetic_cohort() -> list[dict[str, str]]:
    rng = random.Random(4)
    rows = []
    for i in range(120):
        nested = i < 40
        offset = 133 if nested and i % 3 == 0 else 200
        host_len = "300" if nested else ""
        tsd = rng.choice(["12", "20", "41", "82"])
        rows.append(
            make_row(
                pos=str(2000 + i * 7),
                nested_in_alu_host="1" if nested else "0",
                host_len=host_len,
                consensus_offset=str(offset) if nested else "",
                consensus_mapping="unambiguous" if nested else "",
                host_name="AluY" if nested else "",
                tsd_len=tsd,
                tsd_seq="A" * int(tsd),
                polya_len="30",
                polya_seq="A" * 30,
                site_allele_freq="0.0005",
                insert_strand="+",
            )
        )
    return rows


def test_main_runs_end_to_end_on_a_synthetic_cohort(tmp_path, monkeypatch, capsys):
    cohort = tmp_path / "cohort.csv"
    fields = list(make_row())
    with cohort.open("w", newline="") as fh:
        fh.write(",".join(fields) + "\n")
        for row in _synthetic_cohort():
            fh.write(",".join(str(row.get(f, "")) for f in fields) + "\n")

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"gate_0d_tsd_sentinels": {"sentinel_values": {"41": {}, "82": {}}}}))

    outdir = tmp_path / "out"
    code = tprt.main(
        [
            "--cohort", str(cohort),
            "--manifest", str(manifest),
            "--outdir", str(outdir),
            "--permutations", "50",
            "--no-reference",
        ]
    )
    assert code == 0
    report = json.loads((outdir / "tprt_hallmarks.json").read_text())
    assert report["phase"] == "2b"
    assert report["tsd_sentinels"] == [41, 82]
    # no reference was supplied, so every reference-dependent component is
    # unevaluable -- and must not have been promoted to "absent".
    for name in ("analysis_set", "peak", "mid_host_control", "unnested"):
        assert set(report["component_summary"][name]["en_states"]) == {"unevaluable"}
        assert set(report["component_summary"][name]["l1_en_states"]) == {"unevaluable"}
        assert set(report["component_summary"][name]["tsd_states"]) == {"unevaluable"}
    text = (outdir / "tprt_hallmarks.md").read_text()
    assert "Phase 2b" in text
    # The plan reserves "validated" for cross-method evidence that does not exist.
    assert "validated" not in text.lower()
    assert "no precision, recall" in text.lower()


def test_report_states_the_alu_composition_confound(tmp_path, monkeypatch):
    cohort = tmp_path / "cohort.csv"
    fields = list(make_row())
    with cohort.open("w", newline="") as fh:
        fh.write(",".join(fields) + "\n")
        for row in _synthetic_cohort():
            fh.write(",".join(str(row.get(f, "")) for f in fields) + "\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"gate_0d_tsd_sentinels": {"sentinel_values": {}}}))
    outdir = tmp_path / "out"
    tprt.main(
        [
            "--cohort", str(cohort),
            "--manifest", str(manifest),
            "--outdir", str(outdir),
            "--permutations", "20",
            "--no-reference",
        ]
    )
    report = json.loads((outdir / "tprt_hallmarks.json").read_text())
    assert any("confound" in d.lower() for d in report["disclosures"])
    assert any(d["status"] == "unavailable" for d in report["deviations"])


# --------------------------------------------------------------------------
# the shipped cohort, marked slow
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_shipped_cohort_columns_are_all_present():
    cohort = tprt.DEFAULT_COHORT
    if not cohort.exists():
        pytest.skip(f"{cohort} not available")
    rows = tprt.read_cohort(cohort)
    assert rows
    required = {
        "chrom", "pos", "nested_in_alu_host", "host_len", "consensus_offset",
        "consensus_mapping", "tsd_len", "tsd_seq", "polya_len", "polya_seq",
        "insert_strand", "host_name", "conformation", "site_allele_freq",
        "perc_resolved", "alignment_identity",
    }
    assert required <= set(rows[0])