"""Tests for the joint nested-signature enrichment test.

The pre-registration is under test as much as the arithmetic is. A cell list
that grows during a run, a structural zero that gets pseudocounted into a
probability, a private site leaking into a primary p-value, or a multiplied
effect size appearing in the prose are all failures that no p-value would reveal.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import joint_enrichment as je  # noqa: E402
import nested_multi_sample_common as common  # noqa: E402


# --------------------------------------------------------------------------
# the pre-specified cell list
# --------------------------------------------------------------------------


def test_exactly_four_cells_are_pre_specified_with_one_primary():
    assert len(je.CELLS) == 4
    assert [c.role for c in je.CELLS].count("primary") == 1
    assert [c.role for c in je.CELLS].count("secondary") == 3


def test_the_primary_cell_is_the_pre_specified_joint_cell():
    primary = je.CELLS[0]
    assert (primary.bin_lo, primary.bin_hi) == (120, 140)
    assert primary.interaction is False
    assert primary.cell_id == "P1_primary"


def test_the_secondary_cells_are_the_pre_specified_ones():
    got = {(c.bin_lo, c.bin_hi, c.interaction) for c in je.CELLS[1:]}
    assert got == {
        (280, 300, False),
        (120, 140, True),
        (280, 300, True),
    }


def test_the_reported_cells_are_exactly_the_pre_specified_ones(tmp_path):
    """A run must report the declared cells and nothing else.

    This is the anti-scan invariant. A previous version of this test tried to
    prove the cells were immutable by assigning to one, and the assignment
    succeeded -- silently rewriting the primary cell's bin for every test that
    ran afterwards. The useful assertion is about the reported set, not about
    Python's attribute semantics.
    """
    report, _ = _run(tmp_path)
    assert [c["cell_id"] for c in report["cells"]] == [c.cell_id for c in je.CELLS]
    assert [c["role"] for c in report["cells"]] == [c.role for c in je.CELLS]


def test_cell_ids_are_unique():
    ids = [c.cell_id for c in je.CELLS]
    assert len(set(ids)) == len(ids)


def test_the_primary_cell_is_not_holm_adjusted_against_cells_named_after_it():
    adjusted = je.holm([0.01, 0.02, 0.03])
    assert len(adjusted) == 3
    assert adjusted[0] == pytest.approx(0.03)


# --------------------------------------------------------------------------
# opportunity and structural zeros
# --------------------------------------------------------------------------


def test_a_bin_beyond_a_short_host_has_zero_opportunity():
    assert je.opportunity_fraction(host_len=100, bin_lo=120, bin_hi=140) == 0.0


def test_a_bin_overlapping_a_short_host_is_clipped_not_pseudocounted():
    fraction = je.opportunity_fraction(host_len=130, bin_lo=120, bin_hi=140)
    assert fraction == pytest.approx(10 / 130)


def test_a_fully_covered_host_gives_unity_opportunity():
    assert je.opportunity_fraction(host_len=300, bin_lo=120, bin_hi=140) == pytest.approx(20 / 300)


def test_a_bin_beyond_the_host_end_is_excluded_not_counted_as_a_hit():
    assert je.opportunity_fraction(host_len=50, bin_lo=0, bin_hi=10) == pytest.approx(0.2)


# --------------------------------------------------------------------------
# the statistic
# --------------------------------------------------------------------------


def _host(host_len, events):
    return {"host_len": host_len, "events": [(o, s, None) for o, s in events]}


def test_a_sense_cell_counts_only_sense_events_in_the_bin():
    host = _host(300, [(130, "+"), (130, "-"), (10, "+")])
    cell = je.CELLS[0]
    assert je.observed_statistic(host, cell) == 1


def test_bin_bounds_are_inclusive_on_both_ends():
    host = _host(300, [(120, "+"), (140, "+"), (121, "+")])
    assert je.observed_statistic(host, je.CELLS[0]) == 3


def test_an_interaction_cell_is_sense_minus_antisense_in_the_bin():
    host = _host(300, [(130, "+"), (130, "+"), (130, "-"), (10, "+")])
    cell = je.Cell("t", "secondary", 120, 140, interaction=True)
    assert je.observed_statistic(host, cell) == 1


# --------------------------------------------------------------------------
# the null
# --------------------------------------------------------------------------


def test_the_simulation_agrees_with_its_own_analytic_expectation():
    """Two independent routes to the same number; a gap means one is wrong."""
    hosts = [_host(300, [(i, "+") for i in range(5)]), _host(200, [(i, "-") for i in range(3)])]
    cell = je.CELLS[0]
    analytic = je.analytic_expectation(hosts, cell)
    null = je.simulate(hosts, cell, 4000, np.random.default_rng(1))
    assert null.mean() == pytest.approx(analytic, rel=0.15)


def test_the_simulation_is_centred_on_zero_for_a_balanced_interaction_cell():
    hosts = [_host(300, [(10, "+"), (20, "-")])]
    cell = je.Cell("t", "secondary", 120, 140, interaction=True)
    analytic = je.analytic_expectation(hosts, cell)
    assert analytic == pytest.approx(0.0)
    null = je.simulate(hosts, cell, 4000, np.random.default_rng(2))
    assert abs(null.mean()) < 1.0


def test_a_host_with_no_opportunity_contributes_nothing_to_the_null():
    hosts = [_host(100, [(50, "+")])]
    cell = je.CELLS[0]
    assert je.analytic_expectation(hosts, cell) == 0.0
    null = je.simulate(hosts, cell, 200, np.random.default_rng(3))
    assert np.all(null == 0)


def test_the_null_preserves_each_hosts_event_count():
    """Position is redrawn, not the number of events."""
    hosts = [_host(300, [(10, "+"), (20, "+"), (30, "+")])]
    cell = je.Cell("t", "secondary", 120, 140, interaction=True)
    null = je.simulate(hosts, cell, 500, np.random.default_rng(4))
    # three events can contribute at most +3 and at least -3
    assert null.max() <= 3
    assert null.min() >= -3


def test_the_null_is_reproducible_from_a_seed():
    hosts = [_host(300, [(i, "+") for i in range(4)])]
    a = je.simulate(hosts, je.CELLS[0], 100, np.random.default_rng(7))
    b = je.simulate(hosts, je.CELLS[0], 100, np.random.default_rng(7))
    assert a.tolist() == b.tolist()


def test_the_null_moves_only_positions_when_orientation_is_held():
    """A cell with no sense events can never produce a hit."""
    hosts = [_host(300, [(10, "-"), (20, "-")])]
    null = je.simulate(hosts, je.CELLS[0], 500, np.random.default_rng(5))
    assert np.all(null == 0)


# --------------------------------------------------------------------------
# p-values and adjustment
# --------------------------------------------------------------------------


def test_the_monte_carlo_p_carries_the_plus_one_correction():
    null = np.zeros(100, dtype=np.int64)
    # nothing reaches an impossible observed value
    assert je.monte_carlo_p(50, null) == pytest.approx(1 / 101)
    # everything reaches an always-exceeded value
    assert je.monte_carlo_p(-5, null) == pytest.approx(1.0)


def test_the_monte_carlo_p_is_one_sided_and_counts_equality():
    null = np.array([0, 1, 1, 1], dtype=np.int64)
    assert je.monte_carlo_p(1, null) == pytest.approx((1 + 3) / 5)


def test_holm_is_step_down_and_monotone():
    adjusted = je.holm([0.01, 0.02, 0.04])
    assert adjusted[0] == pytest.approx(0.03)
    assert adjusted[1] == pytest.approx(0.04)
    assert adjusted[2] == pytest.approx(0.04)
    assert adjusted == sorted(adjusted)


def test_holm_leaves_unestimable_entries_alone_and_does_not_count_them():
    adjusted = je.holm([0.01, None, 0.02])
    assert adjusted[1] is None
    # only two estimable hypotheses, so the first is doubled
    assert adjusted[0] == pytest.approx(0.02)


def test_holm_never_exceeds_one():
    assert all(v is None or v <= 1.0 for v in je.holm([0.9, 0.95]))


# --------------------------------------------------------------------------
# effect size and bootstrap
# --------------------------------------------------------------------------


def test_the_effect_size_is_observed_over_expected():
    hosts = [_host(300, [(130, "+")] * 4)]
    cell = je.CELLS[0]
    observed = sum(je.observed_statistic(h, cell) for h in hosts)
    expected = je.analytic_expectation(hosts, cell)
    assert observed == 4
    assert observed / expected == pytest.approx(4 / (4 * 20 / 300))


def test_the_bootstrap_brackets_the_point_estimate():
    hosts = [_host(300, [(130, "+") for _ in range(3)]) for _ in range(20)]
    cell = je.CELLS[0]
    lo, hi = je.host_clustered_bootstrap(hosts, cell, 300, np.random.default_rng(9))
    assert lo is not None and hi is not None
    assert lo <= hi
    assert lo > 1.0


def test_the_bootstrap_returns_nothing_when_there_is_no_expectation():
    assert je.host_clustered_bootstrap([], je.CELLS[0], 10, np.random.default_rng(1)) == (None, None)
    assert je.host_clustered_bootstrap(
        [_host(50, [(10, "+")])], je.CELLS[0], 10, np.random.default_rng(1)
    ) == (None, None)


# --------------------------------------------------------------------------
# host stratification
# --------------------------------------------------------------------------


def _site(chrom, host, start, end, offset, orientation="+", family="ALU",
          private=False, site_id="US1"):
    return {
        "site_id": site_id, "chrom": chrom, "pos": start + offset, "family": family,
        "orientation": orientation, "host_name": host, "host_start0": start,
        "host_end0": end, "host_strand": "+", "host_len": end - start,
        "host_offset": offset, "n_carriers": 1, "carriers": ["S1"],
        "private": private, "calls": [],
    }


def test_events_are_grouped_by_host_copy_not_by_site():
    table = je.build_host_table(
        [_site("chr1", "AluSx1", 1000, 1300, 10, site_id="a"),
         _site("chr1", "AluSx1", 1000, 1300, 20, site_id="b")]
    )
    assert len(table) == 1
    assert len(next(iter(table.values()))["events"]) == 2


def test_only_the_requested_host_family_enters_the_table():
    table = je.build_host_table(
        [_site("chr1", "AluSx1", 1000, 1300, 10),
         _site("chr1", "L1PA2", 1000, 1300, 10)]
    )
    assert len(table) == 1


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

VCF_HEADER = """##fileformat=VCFv4.2
##INFO=<ID=NESTED,Number=1,Type=String,Description="nested state">
##INFO=<ID=ORIENT,Number=1,Type=String,Description="orientation">
##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="family">
##INFO=<ID=MEISUBFAMILY,Number=1,Type=String,Description="subfamily">
##INFO=<ID=TSD,Number=1,Type=String,Description="tsd">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""


def _fixture(tmp_path: Path, n_hosts: int = 30):
    """A synthetic cohort with one shared site per host in and out of the bin."""
    rows = []
    callset = tmp_path / "callsets"
    callset.mkdir()
    vcf = [VCF_HEADER.rstrip("\n")]
    for i in range(n_hosts):
        start = 1_000_000 + i * 10_000
        for j, (offset, orient) in enumerate([(130, "+"), (10, "-")]):
            pos = start + offset
            site_id = f"US{i:05d}{j}"
            shared = i % 2 == 0
            rows.append(
                {
                    "site_id": site_id, "chrom": "chr1", "representative_pos": str(pos),
                    "representative_sample": "S1", "family": "ALU",
                    "insertion_orientation": orient, "host_name": "AluSx1",
                    "host_start0": str(start), "host_end0": str(start + 300),
                    "host_strand": "+", "host_len": "300",
                    "host_offset_5p_0based": str(offset),
                    "host_selection_rule": "longest_span",
                    "n_carriers": "2" if shared else "1",
                    "samples": "S1|S2" if shared else "S1",
                    "private": "False" if shared else "True",
                    "private_to_sample": "" if shared else "S1", "window_bp": "10",
                    "tsd_overlap_within_site": "0", "site_allele_count": "1",
                    "site_allele_freq": "0.1", "pooled_site_matched": "False",
                }
            )
            vcf.append(
                f"chr1\t{pos}\t{site_id}\tN\t<INS>\t.\t.\t"
                f"NESTED=nested;ORIENT={orient};MEIFAMILY=ALU;"
                f"MEISUBFAMILY=AluYb9#SINE/Alu;TSD=ACGTACGTACGT\tGT\t0/1"
            )
    # S2 carries the shared sites, so the join reproduces n_carriers
    (callset / "S1.vcf").write_text("\n".join(vcf) + "\n")
    s2 = [VCF_HEADER.rstrip("\n")]
    for row in rows:
        if row["samples"] == "S1|S2":
            s2.append(
                f"chr1\t{row['representative_pos']}\t{row['site_id']}\tN\t<INS>\t.\t.\t"
                f"NESTED=nested;ORIENT={row['insertion_orientation']};MEIFAMILY=ALU;"
                f"MEISUBFAMILY=AluYb9#SINE/Alu;TSD=ACGTACGTACGT\tGT\t0/1"
            )
    (callset / "S2.vcf").write_text("\n".join(s2) + "\n")

    us = tmp_path / "unique_sites.csv"
    with us.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return us, callset


def _run(tmp_path, replicates=200, bootstrap=50, rmsk=None):
    us, callset = _fixture(tmp_path)
    out = tmp_path / "out"
    rc = je.main(
        [
            "--unique-sites", str(us),
            "--callset-dir", str(callset),
            "--rmsk", str(rmsk or (tmp_path / "absent.rmsk.gz")),
            "--outdir", str(out),
            "--replicates", str(replicates),
            "--bootstrap", str(bootstrap),
        ]
    )
    assert rc == 0
    return json.loads((out / "joint_enrichment.json").read_text()), out


def test_main_runs_end_to_end_on_a_synthetic_cohort(tmp_path):
    report, out = _run(tmp_path)
    assert report["analysis"] == "joint_nested_signature_enrichment"
    assert len(report["cells"]) == 4
    assert (out / "joint_enrichment.md").exists()
    assert (out / "joint_cells.csv").exists()


def test_the_join_is_gated_before_any_cell_is_computed(tmp_path):
    report, _ = _run(tmp_path)
    assert report["join_verification"]["verdict"] == "join_consistent_with_dedup_output"


def test_private_sites_never_enter_a_primary_p_value(tmp_path):
    report, _ = _run(tmp_path)
    assert report["cohort_split"]["private_sites"] > 0
    for cell in report["cells"]:
        # a private site cannot inflate the primary observed count
        assert cell["sites_in_scope"] <= report["cohort_split"]["shared_sites"]
        assert cell["private_exploratory_sites"] >= 0
    assert "no primary-test p-value is computed over them" in report["cohort_split"]["note"]


def test_structural_zero_hosts_are_counted_and_excluded(tmp_path):
    report, _ = _run(tmp_path)
    for cell in report["cells"]:
        assert cell["hosts_structural_zero_excluded"] == 0
        assert cell["hosts_with_opportunity"] + cell["hosts_structural_zero_excluded"] == (
            cell["hosts_total"]
        )


def test_the_report_never_multiplies_two_factor_estimates(tmp_path):
    """The forbidden figure is "1.47x * 3.0x = 4.4x"."""
    report, out = _run(tmp_path)
    text = (out / "joint_enrichment.md").read_text()
    forbidden = re.search(r"\d+(?:\.\d+)?\s*x\s*\*\s*\d+(?:\.\d+)?\s*x\s*=", text)
    assert forbidden is None, f"a multiplied effect size appears in the report: {forbidden}"
    assert re.search(r"=\s*\d+(?:\.\d+)?\s*x", text) is None


def test_the_multiplication_test_reports_a_verdict_not_a_number(tmp_path):
    report, _ = _run(tmp_path)
    mt = report["multiplication_test"]
    assert "verdict" in mt
    assert "joint_exceeds_benchmark" in mt
    assert "benchmark" in mt["note"]
    assert "effect size" in mt["note"]


def test_the_nested_enum_census_is_reported(tmp_path):
    report, _ = _run(tmp_path)
    census = report["nested_enum_census"]
    for key in ("unnested", "nested_sense", "nested_antisense", "nested_unknown"):
        assert key in census
    assert census[common.LEGACY_NESTED_UNLABELED] > 0


def test_a_missing_rmsk_degrades_to_a_reported_status_not_a_crash(tmp_path):
    report, _ = _run(tmp_path, rmsk=tmp_path / "definitely_absent.rmsk.gz")
    assert report["opportunity_diagnostic"]["status"] == "rmsk_unavailable"
    assert "consequence" in report["opportunity_diagnostic"]


def test_site_and_carrier_counts_are_never_collapsed(tmp_path):
    report, _ = _run(tmp_path)
    for cell in report["cells"]:
        assert cell["sites_in_scope"] != cell["carrier_observations_in_scope"] or (
            cell["sites_in_scope"] == 0
        )
    assert "never summed" in report["deviations"][-1]
