"""Tests for the same-host recurrence test.

The definitions here are pre-specified, so the tests are mostly about the
classification's precedence and its refusals: a missing TSD must not be read as a
TSD of length zero, absent data must not mask evidence that is already present,
and no pair may be adjudicated in code.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import nested_multi_sample_common as common  # noqa: E402
import recurrence_test as rt  # noqa: E402


def _call(
    subfamily="AluYb9", tsd="ACGTACGTACGT", orientation="+", sample="S1", pos=1000
) -> dict:
    return {
        "pos": pos,
        "id": "x",
        "family": "ALU",
        "subfamily": subfamily,
        "orientation": orientation,
        "nested_raw": "nested",
        "nested_state": common.LEGACY_NESTED_UNLABELED,
        "tsd": tsd,
        "tsd_len": len(tsd),
        "tsd_reported": bool(tsd),
        "mei_span": 282,
        "known_mei": "",
        "call_tier": "",
    }


# --------------------------------------------------------------------------
# TSD preprocessing
# --------------------------------------------------------------------------


def test_a_poly_a_tail_is_stripped_before_any_comparison():
    assert rt.strip_terminal_homopolymer("GGGCCCAAAA") == "GGGCCC"
    assert rt.strip_terminal_homopolymer("GGGGTTTT") == "GGGG"


def test_a_sequence_ending_in_a_non_homopolymer_base_is_unchanged():
    assert rt.strip_terminal_homopolymer("ACGTACG") == "ACGTACG"
    assert rt.strip_terminal_homopolymer("ACGTACC") == "ACGTACC"


def test_a_poly_t_tail_is_also_stripped():
    """T is the reverse-complement poly-A tail and carries the same information."""
    assert rt.strip_terminal_homopolymer("ACGTACGTTT") == "ACGTACG"


def test_an_empty_sequence_stays_empty():
    assert rt.strip_terminal_homopolymer("") == ""


def test_cores_shorter_than_the_threshold_are_not_comparable():
    assert rt.tsd_cores_comparable("AAA", "AAA")[0] is False


def test_cores_at_the_threshold_are_comparable():
    comparable, a, b = rt.tsd_cores_comparable("ACGTAC", "ACGTAG")
    assert comparable is True
    assert a == "ACGTAC" and b == "ACGTAG"


# --------------------------------------------------------------------------
# classification: the branches
# --------------------------------------------------------------------------


def test_matching_subfamily_tsd_and_core_is_one_inherited_event():
    out = rt.classify_pair(_call(), _call())
    assert out["verdict"] == rt.IBD
    assert out["tsd_sequence_discordant"] is False


def test_a_different_child_subfamily_is_independent_recurrence():
    out = rt.classify_pair(_call(subfamily="AluYb9"), _call(subfamily="AluYb8"))
    assert out["verdict"] == rt.CANDIDATE
    assert "subfamily differs" in out["reason"]


def test_a_tsd_length_difference_beyond_tolerance_is_independent_recurrence():
    # Identical cores, tails differing by 4 bp, so length is the only discriminator.
    a = _call(tsd="ACGTACGTACGTGC")        # 15 bp
    b = _call(tsd="ACGTACGTACGTGCAAAA")    # 19 bp
    out = rt.classify_pair(a, b)
    assert out["verdict"] == rt.CANDIDATE
    assert out["tsd_len_abs_difference"] == 4
    assert "TSD length differs" in out["reason"]


def test_a_tsd_length_difference_inside_tolerance_is_not_a_discordance():
    # Same core, tails differing by 2 bp: the tail is stripped, so the cores match.
    a = _call(tsd="ACGTACGTACGTGC")
    b = _call(tsd="ACGTACGTACGTGCAA")
    out = rt.classify_pair(a, b)
    assert out["tsd_len_abs_difference"] == 2
    assert out["verdict"] == rt.IBD


def test_a_differing_tsd_core_is_a_discordance_even_at_equal_length():
    a = _call(tsd="ACGTACGTACGTGC")
    b = _call(tsd="ACGTACGTACGTGG")
    out = rt.classify_pair(a, b)
    assert out["verdict"] == rt.CANDIDATE
    assert out["tsd_sequence_discordant"] is True


def test_differing_only_in_the_shared_tail_is_not_a_discordance():
    """A poly-A tail is shared by construction and carries no identity."""
    a = _call(tsd="ACGTACGTAA")
    b = _call(tsd="ACGTACGTAAAA")
    out = rt.classify_pair(a, b)
    assert out["tsd_sequence_discordant"] is False
    assert out["verdict"] == rt.IBD


# --------------------------------------------------------------------------
# classification: the refusals
# --------------------------------------------------------------------------


def test_a_missing_tsd_is_unevaluable_not_a_zero_length_tsd():
    out = rt.classify_pair(_call(tsd=""), _call(tsd="ACGTACGT"))
    assert out["verdict"] == rt.UNEVALUABLE
    assert out["tsd_reported_both"] is False
    assert out["tsd_len_abs_difference"] is None
    assert "missing TSD" in out["reason"] or "reports no TSD" in out["reason"]


def test_two_missing_tsds_are_still_unevaluable_not_concordant():
    out = rt.classify_pair(_call(tsd=""), _call(tsd=""))
    assert out["verdict"] == rt.UNEVALUABLE


def test_a_subfamily_discordance_is_not_masked_by_a_missing_tsd():
    """Precedence regression.

    The definition ORs its discordance conditions, so a differing child
    subfamily is sufficient on its own. Testing 'is a TSD missing?' first lets
    absent data hide positive evidence that is already present. This exact case
    was misclassified as unevaluable in an earlier revision.
    """
    out = rt.classify_pair(
        _call(subfamily="AluYk12", tsd=""), _call(subfamily="AluSq2", tsd="")
    )
    assert out["verdict"] == rt.CANDIDATE
    assert "sufficient" in out["reason"]


def test_a_short_tsd_core_is_unevaluable_rather_than_concordant():
    out = rt.classify_pair(_call(tsd="AAA"), _call(tsd="AAA"))
    assert out["verdict"] == rt.UNEVALUABLE
    assert rt.TSD_CORE_MIN_BP == 4


# --------------------------------------------------------------------------
# pair discovery
# --------------------------------------------------------------------------


def _site(site_id, pos, host_start=1000, host_end=1300, host="AluSx1",
          family="ALU", offset=None, carriers=("S1",)):
    return {
        "site_id": site_id, "chrom": "chr1", "pos": pos, "family": family,
        "orientation": "+", "host_name": host, "host_start0": host_start,
        "host_end0": host_end, "host_strand": "+", "host_len": host_end - host_start,
        "host_offset": offset if offset is not None else pos - host_start,
        "n_carriers": len(carriers), "carriers": list(carriers), "private": False,
        "calls": [_call(sample=s, pos=pos) for s in carriers],
    }


def test_only_hosts_with_more_than_one_event_contribute():
    d = rt.find_same_host_pairs(
        [_site("a", 1100), _site("b", 1200, host_start=5000, host_end=5300, host="AluY")],
        10,
    )
    assert d["hosts_total"] == 2
    assert d["hosts_with_more_than_one_nested_event"] == 0
    assert d["pairs_within_window"] == 0


def test_a_pair_inside_the_window_is_found():
    d = rt.find_same_host_pairs([_site("a", 1100), _site("b", 1105)], 10)
    assert d["pairs_within_window"] == 1
    assert d["pairs"][0]["pos_delta_bp"] == 5


def test_a_pair_outside_the_window_is_not_found():
    d = rt.find_same_host_pairs([_site("a", 1100), _site("b", 1200)], 10)
    assert d["pairs_within_window"] == 0


def test_the_widening_window_finds_more_pairs():
    sites = [_site("a", 1100), _site("b", 1115)]
    assert rt.find_same_host_pairs(sites, 5)["pairs_within_window"] == 0
    assert rt.find_same_host_pairs(sites, 20)["pairs_within_window"] == 1


def test_a_host_is_keyed_by_interval_not_by_the_inserted_element():
    same_host = [_site("a", 1100, family="ALU"), _site("b", 1105, family="SVA")]
    assert rt.find_same_host_pairs(same_host, 10)["hosts_with_more_than_one_nested_event"] == 1


def test_discovery_breaks_once_the_separation_exceeds_the_window():
    sites = [_site("a", 1100), _site("b", 1105), _site("c", 1200), _site("d", 1205)]
    d = rt.find_same_host_pairs(sites, 10)
    # (a,b) and (c,d) only; a naive all-pairs scan would also emit (a,c) etc.
    assert d["pairs_within_window"] == 2


# --------------------------------------------------------------------------
# chance expectations
# --------------------------------------------------------------------------


def test_the_positional_chance_rate_uses_two_w_plus_one_over_host_length():
    sites = [_site("a", 1100), _site("b", 1200)]
    out = rt.positional_chance_rate(
        rt.find_same_host_pairs(sites, 10), sites, 10
    )
    assert out["expected_pairs_under_chance"] == pytest.approx(21 / 300)


def test_a_short_host_has_a_higher_positional_chance_rate():
    short = [_site("a", 1100, host_start=1000, host_end=1050),
             _site("b", 1020, host_start=1000, host_end=1050)]
    long = [_site("c", 1100, host_start=1000, host_end=1300),
            _site("d", 1200, host_start=1000, host_end=1300)]
    a = rt.positional_chance_rate(rt.find_same_host_pairs(short, 10), short, 10)
    b = rt.positional_chance_rate(rt.find_same_host_pairs(long, 10), long, 10)
    assert a["expected_pairs_under_chance"] > b["expected_pairs_under_chance"]


def test_the_discordance_rate_is_estimated_from_the_observed_distribution():
    sites = [_site("a", 1100, carriers=("S1", "S2", "S3"))]
    sites[0]["calls"] = [
        _call(subfamily="AluYb9", sample="S1"),
        _call(subfamily="AluYb8", sample="S2"),
        _call(subfamily="AluYb7", sample="S3"),
    ]
    out = rt.discordance_chance(sites, 400, seed=1)["Alu"]
    assert out["status"] == "estimated"
    # three distinct subfamilies: two distinct calls can never agree
    assert out["p_subfamily_discordant"] == pytest.approx(1.0)


def test_the_discordance_rate_does_not_let_a_call_pair_with_itself():
    """Two independent insertions are two distinct calls.

    Sampling with replacement and allowing i == j would let a call agree with
    itself, biasing the discordance rate down by 1/n and flattering the
    definition.
    """
    sites = [_site("a", 1100, carriers=("S1", "S2"))]
    sites[0]["calls"] = [
        _call(subfamily="AluYb9", sample="S1"),
        _call(subfamily="AluYb9", sample="S2"),
    ]
    out = rt.discordance_chance(sites, 400, seed=1)["Alu"]
    # the only distinct pair shares a subfamily, so nothing can be discordant
    assert out["p_subfamily_discordant"] == pytest.approx(0.0)


def test_the_ibd_affirmation_rate_is_measured_not_assumed():
    sites = [_site("a", 1100, carriers=("S1", "S2", "S3"))]
    sites[0]["calls"] = [
        _call(subfamily="AluYb9", tsd="ACGTACGT", sample="S1"),
        _call(subfamily="AluYb8", tsd="ACGTACGT", sample="S2"),
        _call(subfamily="AluYb7", tsd="ACGTACGT", sample="S3"),
    ]
    out = rt.ibd_affirmation_chance(sites, 200, seed=1)["Alu"]
    assert out["status"] == "estimated"
    # distinct subfamilies cannot satisfy the IBD definition
    assert out["p_would_be_called_ibd_under_chance"] == pytest.approx(0.0)


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


HEADER = [
    "site_id", "chrom", "representative_pos", "representative_sample", "family",
    "insertion_orientation", "host_name", "host_start0", "host_end0", "host_strand",
    "host_len", "host_offset_5p_0based", "host_selection_rule", "n_carriers",
    "samples", "private", "private_to_sample", "window_bp", "tsd_overlap_within_site",
    "site_allele_count", "site_allele_freq", "pooled_site_matched",
]

VCF_HEADER = "\n".join(
    [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=NESTED,Number=1,Type=String,Description="nested state">',
        '##INFO=<ID=ORIENT,Number=1,Type=String,Description="orientation">',
        '##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="family">',
        '##INFO=<ID=MEISUBFAMILY,Number=1,Type=String,Description="subfamily">',
        '##INFO=<ID=TSD,Number=1,Type=String,Description="tsd">',
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
)


def _fixture(tmp_path: Path):
    """Two sites in one host within the window, discordant, plus singleton hosts."""
    callset = tmp_path / "callsets"
    callset.mkdir()
    vcf = [VCF_HEADER]
    rows = []

    def add(site_id, pos, start, subfamily, tsd):
        rows.append(
            {
                "site_id": site_id, "chrom": "chr1", "representative_pos": str(pos),
                "representative_sample": "S1", "family": "ALU",
                "insertion_orientation": "+", "host_name": "AluSx1",
                "host_start0": str(start), "host_end0": str(start + 300),
                "host_strand": "+", "host_len": "300",
                "host_offset_5p_0based": str(pos - start),
                "host_selection_rule": "longest_span", "n_carriers": "1",
                "samples": "S1", "private": "True", "private_to_sample": "S1",
                "window_bp": "10", "tsd_overlap_within_site": "0",
                "site_allele_count": "1", "site_allele_freq": "0.2",
                "pooled_site_matched": "False",
            }
        )
        vcf.append(
            f"chr1\t{pos}\t{site_id}\tN\t<INS>\t.\t.\tNESTED=nested;ORIENT=+;"
            f"MEIFAMILY=ALU;MEISUBFAMILY={subfamily}#SINE/Alu;TSD={tsd}\tGT\t0/1"
        )

    # the candidate pair: same host, 5 bp apart, different child subfamily
    add("US00001", 1_000_130, 1_000_000, "AluYb8", "ACGTACGTACGTAC")
    add("US00002", 1_000_135, 1_000_000, "AluSq2", "ACGTACGTACGTAC")
    # a lone site in another host
    add("US00003", 2_000_130, 2_000_000, "AluYb8", "ACGTACGTACGTAC")
    # a site with no TSD anywhere, for the unevaluable branch
    add("US00004", 3_000_100, 3_000_000, "AluYb8", "")
    add("US00005", 3_000_104, 3_000_000, "AluYb8", "")

    (callset / "S1.vcf").write_text("\n".join(vcf) + "\n")
    us = tmp_path / "unique_sites.csv"
    with us.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return us, callset


def _run(tmp_path, draws=200):
    us, callset = _fixture(tmp_path)
    out = tmp_path / "out"
    rc = rt.main(
        [
            "--unique-sites", str(us),
            "--callset-dir", str(callset),
            "--outdir", str(out),
            "--chance-draws", str(draws),
        ]
    )
    assert rc == 0
    return json.loads((out / "recurrence.json").read_text()), out


def test_main_runs_end_to_end_on_a_synthetic_cohort(tmp_path):
    report, out = _run(tmp_path)
    assert report["analysis"] == "same_host_recurrence_ibd_vs_independent"
    assert (out / "recurrence.md").exists()
    assert (out / "recurrence_pairs.csv").exists()


def test_the_planted_candidate_pair_is_recovered(tmp_path):
    report, _ = _run(tmp_path)
    assert report["any_candidate_recurrent"] == 1
    assert report["verdict_summary"][rt.CANDIDATE]["n_pairs"] == 1


def test_the_planted_no_tsd_pair_is_unevaluable_not_ibd(tmp_path):
    report, _ = _run(tmp_path)
    assert report["verdict_summary"][rt.UNEVALUABLE]["n_pairs"] == 1
    assert report["verdict_summary"][rt.IBD]["n_pairs"] == 0


def test_no_pair_is_adjudicated_in_the_output(tmp_path):
    report, out = _run(tmp_path)
    for row in report["candidate_evidence_lines"]:
        assert row["adjudicated"].startswith("no")
    text = (out / "recurrence.md").read_text()
    assert "does not make it" in text or "deliberately does not" in text


def test_candidate_evidence_lines_carry_the_full_record(tmp_path):
    report, _ = _run(tmp_path)
    row = report["candidate_evidence_lines"][0]
    for key in (
        "chrom", "pos_a", "pos_b", "pos_delta_bp", "host_name", "host_start0",
        "host_end0", "decisive_child_subfamily_a", "decisive_child_subfamily_b",
        "decisive_sample_a", "decisive_sample_b", "decisive_reason", "carriers_a",
        "carriers_b",
    ):
        assert key in row, key
    assert row["decisive_child_subfamily_a"] != row["decisive_child_subfamily_b"]


def test_the_window_sweep_is_reported_in_full(tmp_path):
    report, _ = _run(tmp_path)
    assert set(report["sweep"]) == {"10", "5", "20"}
    for window, block in report["sweep"].items():
        assert "chance_expected_pairs" in block
        assert "verdicts" in block


def test_the_headline_states_the_chance_comparison_not_just_a_count(tmp_path):
    _, out = _run(tmp_path)
    text = (out / "recurrence.md").read_text()
    assert "does not clear its own chance expectation" in text
    assert "first test of its kind" in text


def test_the_ibd_affirmation_caveat_is_reported(tmp_path):
    report, out = _run(tmp_path)
    assert "ibd_affirmation_chance" in report
    text = (out / "recurrence.md").read_text()
    assert "not evidence that no shared events exist" in text


def test_the_denominator_caveat_rejects_comparison_to_pooled_svAN(tmp_path):
    report, out = _run(tmp_path)
    assert "26/2,559" in out.joinpath("recurrence.md").read_text()
    assert "not comparable" in out.joinpath("recurrence.md").read_text()
    for block in report["per_family"].values():
        assert "hosts_total" in block
        assert "carrier_observations" in block
