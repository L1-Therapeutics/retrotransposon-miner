"""Tests for the Phase 3 genotype-concordance scorer.

The recurring theme here is that every quantity Phase 3 reports is a *negative*
result or a bounded one: no missingness channel, a comparator that cannot see
the locus, a recurrence axis that is mostly singletons. Tests therefore spend
most of their effort on the places where such a phase can quietly lie -- by
collapsing an unevaluable row into an absent one, by mis-aligning genotypes
against the sample list, by splitting the cohort in a way that depends on the
genotypes, or by emitting a metric the gates forbid.
"""

from __future__ import annotations

import ast
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import score_genotype_concordance as p3  # noqa: E402
import score_tprt_hallmarks as tprt  # noqa: E402


# --------------------------------------------------------------------------
# genotype parsing
# --------------------------------------------------------------------------


def test_a_non_ref_allele_makes_a_carrier():
    tokens = ["0/0", "0/1", "1/1", "1|0", "0|1"]
    mask = p3.parse_genotype_block(tokens)
    assert mask.tolist() == [False, True, True, True, True]


def test_a_missing_genotype_is_neither_carrier_nor_demonstrated_negative():
    tokens = ["./.", ".", "", "0/0", "0/1"]
    mask = p3.parse_genotype_block(tokens)
    assert mask.tolist() == [False, False, False, False, True]


def test_multi_allelic_and_haploid_tokens_parse():
    tokens = ["0/2", "2", "0|2", "0/0"]
    assert p3.parse_genotype_block(tokens).tolist() == [True, True, True, False]


def test_the_query_format_expands_one_field_per_sample():
    """A bare `[%GT]` concatenates genotypes with no separator.

    That failure is silent: it yields a one-field record rather than an error, and
    the array it produces is short enough to mis-align against the sample list.
    The separator therefore has to be inside the brackets.
    """
    fmt = p3.GENOTYPE_QUERY_FORMAT
    assert fmt == "%CHROM\\t%POS[\\t%GT]\\n"
    bracket = fmt[fmt.index("[") + 1 : fmt.index("]")]
    assert bracket.startswith("\\t"), "separator must be inside the brackets"
    assert "%GT" in bracket


def test_a_record_with_the_wrong_sample_count_is_refused_not_truncated(tmp_path, monkeypatch):
    """A short genotype array must be named, never silently used."""
    bcf = tmp_path / "fake.bcf"
    bcf.write_text("")
    monkeypatch.setattr(p3, "sample_names", lambda path: ["S1", "S2", "S3"])
    monkeypatch.setattr(
        p3,
        "_bcftools",
        lambda args: "chr1\t100\t0/0\t0/1\n",  # only two of three samples
    )
    _, carriers, counters = p3.read_site_genotypes(bcf, [("chr1", 100)])
    assert carriers == {}, "a short record must not be used"
    assert counters["n_records_with_unexpected_sample_count"] == 1
    assert counters["records_with_unexpected_sample_count"] == ["chr1:100"]


def test_well_formed_records_are_kept_and_counted(tmp_path, monkeypatch):
    bcf = tmp_path / "fake.bcf"
    bcf.write_text("")
    monkeypatch.setattr(p3, "sample_names", lambda path: ["S1", "S2", "S3"])
    monkeypatch.setattr(p3, "_bcftools", lambda args: "chr1\t100\t0/0\t0/1\t1/1\n")
    _, carriers, counters = p3.read_site_genotypes(bcf, [("chr1", 100)])
    assert carriers[("chr1", 100)].tolist() == [False, True, True]
    assert counters["n_missing_genotypes"] == 0
    assert counters["n_records_with_unexpected_sample_count"] == 0


def test_a_record_outside_the_requested_sites_is_ignored(tmp_path, monkeypatch):
    """A +/-1 bp region window returns neighbours; they are not analysis sites."""
    bcf = tmp_path / "fake.bcf"
    bcf.write_text("")
    monkeypatch.setattr(p3, "sample_names", lambda path: ["S1"])
    monkeypatch.setattr(p3, "_bcftools", lambda args: "chr1\t100\t1/1\nchr1\t101\t0/1\n")
    _, carriers, _ = p3.read_site_genotypes(bcf, [("chr1", 100)])
    assert set(carriers) == {("chr1", 100)}


# --------------------------------------------------------------------------
# contig naming
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "left,right",
    [("chr1", "1"), ("chrX", "X"), ("1", "1"), ("chrM", "M"), ("GL000009v2", "GL000009v2")],
)
def test_contig_styles_agree(left, right):
    assert p3.normalise_contig(left) == p3.normalise_contig(right)


def test_a_contig_named_chr_only_at_the_front_is_not_mangled():
    assert p3.normalise_contig("chrchr1") == "chr1"


# --------------------------------------------------------------------------
# the half split
# --------------------------------------------------------------------------


def test_the_split_is_balanced_to_within_one_genome():
    halves = p3.split_halves([f"S{i:03d}" for i in range(908)])
    assert halves["n_half_a"] == 454
    assert halves["n_half_b"] == 454
    assert halves["n_samples"] == 908


def test_the_split_does_not_depend_on_the_input_order():
    names = [f"S{i:03d}" for i in range(50)]
    shuffled = list(reversed(names))
    a = p3.split_halves(names)
    b = p3.split_halves(shuffled)
    assert a["order"] == b["order"]
    assert a["index_a"].tolist() == b["index_a"].tolist()
    assert a["index_b"].tolist() == b["index_b"].tolist()


def test_the_two_halves_partition_the_cohort_without_overlap():
    halves = p3.split_halves([f"S{i:03d}" for i in range(37)])
    a = set(halves["index_a"].tolist())
    b = set(halves["index_b"].tolist())
    assert a & b == set()
    assert a | b == set(range(37))


# --------------------------------------------------------------------------
# recurrence
# --------------------------------------------------------------------------


def _mask(*carriers: int, n: int = 8) -> np.ndarray:
    m = np.zeros(n, dtype=bool)
    for i in carriers:
        m[i] = True
    return m


def test_a_site_needs_a_carrier_in_each_half_to_be_cross_half_supported():
    halves = p3.split_halves([f"S{i:02d}" for i in range(8)])
    keys = [("chr1", 1), ("chr1", 2), ("chr1", 3), ("chr1", 4)]
    carriers = {
        # one carrier in each half -> supported
        ("chr1", 1): _mask(0, 5),
        # both carriers in half A -> not supported
        ("chr1", 2): _mask(0, 2),
        # both carriers in half B -> not supported
        ("chr1", 3): _mask(5, 7),
        # a single carrier cannot be cross-half by arithmetic
        ("chr1", 4): _mask(3),
    }
    rec = p3.site_recurrence(keys, carriers, halves)
    assert rec[("chr1", 1)]["cross_half_supported"] is True
    assert rec[("chr1", 2)]["cross_half_supported"] is False
    assert rec[("chr1", 3)]["cross_half_supported"] is False
    assert rec[("chr1", 4)]["cross_half_supported"] is False
    assert rec[("chr1", 4)]["private_to_one_genome"] is True


def test_a_site_with_no_genotype_record_is_unevaluable_not_absent():
    halves = p3.split_halves([f"S{i:02d}" for i in range(4)])
    rec = p3.site_recurrence([("chr1", 9)], {}, halves)
    entry = rec[("chr1", 9)]
    assert entry["state"] == "unevaluable_no_genotype_record"
    assert entry["cross_half_supported"] is None
    assert entry["carriers_total"] is None, "an unevaluable site has no carrier count"
    assert entry["private_to_one_genome"] is None


def test_carrier_counts_in_the_two_halves_sum_to_the_total():
    halves = p3.split_halves([f"S{i:02d}" for i in range(10)])
    rec = p3.site_recurrence(
        [("chr1", 1)], {("chr1", 1): _mask(0, 1, 2, 7, 9, n=10)}, halves
    )[("chr1", 1)]
    assert rec["carriers_half_a"] + rec["carriers_half_b"] == rec["carriers_total"]
    assert rec["carriers_total"] == 5


# --------------------------------------------------------------------------
# accumulation
# --------------------------------------------------------------------------


def test_unique_sites_saturate_and_never_exceed_the_number_of_sites():
    keys = [("chr1", i) for i in range(10)]
    carriers = {
        ("chr1", 0): _mask(0, 1, n=4),
        ("chr1", 1): _mask(1, n=4),
        ("chr1", 2): _mask(3, n=4),
    }
    curve = p3.accumulation_curve(keys, carriers, [f"S{i}" for i in range(4)])
    counts = [pt["unique_sites"] for pt in curve["checkpoints"]]
    assert counts == sorted(counts), "unique sites cannot decrease"
    assert curve["final_unique_sites"] == 3
    assert max(counts) <= len(keys)


def test_a_site_carried_twice_is_one_site_and_two_observations():
    keys = [("chr1", 0), ("chr1", 1)]
    carriers = {("chr1", 0): _mask(0, 1, n=4), ("chr1", 1): _mask(2, n=4)}
    curve = p3.accumulation_curve(keys, carriers, [f"S{i}" for i in range(4)])
    assert curve["final_unique_sites"] == 2
    assert curve["final_carrier_observations"] == 3
    assert curve["final_carrier_observations"] > curve["final_unique_sites"]
    assert curve["mean_carriers_per_unique_site"] == pytest.approx(1.5)


def test_saturation_genome_counts_are_monotone_and_bounded():
    keys = [("chr1", i) for i in range(20)]
    carriers = {("chr1", i): _mask(i % 4, n=6) for i in range(20)}
    curve = p3.accumulation_curve(keys, carriers, [f"S{i}" for i in range(6)])
    reached = [
        v
        for v in curve["genomes_to_reach_fraction_of_final_unique_sites"].values()
        if v is not None
    ]
    assert reached == sorted(reached)
    assert all(1 <= v <= 6 for v in reached)


def test_an_unreached_saturation_target_is_reported_as_not_reached():
    """A fraction of zero sites is trivially reached; the guard must not lie."""
    curve = p3.accumulation_curve([], {}, ["S0", "S1"])
    assert curve["final_unique_sites"] == 0
    assert curve["mean_carriers_per_unique_site"] is None


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------


def test_fisher_exact_is_one_when_the_two_arms_are_identical():
    assert p3.fisher_exact_two_sided(10, 10, 10, 10) == pytest.approx(1.0)


def test_fisher_exact_detects_a_planted_absence():
    p = p3.fisher_exact_two_sided(2, 98, 50, 50)
    assert p is not None and p < 1e-6


def test_fisher_exact_is_symmetric_in_its_arms():
    a = p3.fisher_exact_two_sided(5, 45, 20, 30)
    b = p3.fisher_exact_two_sided(20, 30, 5, 45)
    assert a == pytest.approx(b)


def test_fisher_exact_declines_a_degenerate_table():
    assert p3.fisher_exact_two_sided(0, 0, 0, 0) is None
    assert p3.fisher_exact_two_sided(0, 10, 0, 10) is None


# --------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------


def test_gate_0c_reports_non_independence_when_the_focus_sample_is_inside():
    gate = p3.gate_0c_independence(["HG03086", "HG00096"])
    assert gate["verdict"] == "not_independent"
    assert gate["tier"] == "restricted_concordance_only"
    assert gate["focus_sample_in_genotype_cohort"] is True


def test_gate_0c_flips_when_the_focus_sample_is_absent():
    gate = p3.gate_0c_independence(["HG00096", "HG00099"])
    assert gate["verdict"] == "independent_of_pooled_samples"
    assert gate["focus_sample_in_genotype_cohort"] is False


def test_gate_0c_tier_does_not_depend_on_the_verdict():
    """Even an independent cohort stays in the restricted tier under plan s5."""
    assert (
        p3.gate_0c_independence(["HG00096"])["tier"] == "restricted_concordance_only"
    )


def test_gate_3b_flags_a_callset_with_no_missingness_channel():
    gate = p3.gate_3b_genotype_channel({"n_missing_genotypes": 0}, 100, 908)
    assert gate["verdict"] == "no_missingness_channel_zero_zero_uninterpretable"
    assert gate["missingness_channel_present"] is False
    assert "NOT EXECUTABLE" in gate["consequence"]


def test_gate_3b_detects_a_missingness_channel_when_one_exists():
    gate = p3.gate_3b_genotype_channel({"n_missing_genotypes": 17}, 100, 908)
    assert gate["verdict"] == "missingness_channel_present"
    assert gate["missingness_channel_present"] is True
    assert gate["genotypes_examined"] == 90_800


def test_gate_3c_reports_column_order_divergence_as_a_by_name_join():
    genotypes = ["S1", "S2", "S3", "S4"]
    cross = ["S1", "S3", "S2", "S4"]
    gate = p3.gate_3c_sample_identity(genotypes, cross, Path("x.vcf.gz"))
    assert gate["verdict"] == "samples_joined_by_name"
    assert gate["column_order_identical"] is False
    assert gate["join_is_by_name"] is True
    assert gate["n_samples_shared"] == 4


def test_gate_3c_degrades_honestly_without_a_comparator():
    gate = p3.gate_3c_sample_identity(["S1"], None, None)
    assert gate["verdict"] == "cross_method_comparator_unavailable"
    assert gate["consequence"]


def test_gate_3d_declines_to_measure_without_a_comparator():
    gate = p3.gate_3d_cross_method_floor(None, [("chr1", 1)], [("chr1", 2)])
    assert gate["verdict"] == "cross_method_comparator_unavailable"
    assert "nested_arm" not in gate


# --------------------------------------------------------------------------
# the matching registration must not disturb Phase 2b
# --------------------------------------------------------------------------


def test_registering_the_phase3_comparison_leaves_phase2b_definitions_alone():
    before_covariates = {
        k: dict(v) for k, v in tprt.COVARIATE_SETS.items() if k != p3.COMPARISON
    }
    before_ladders = {k: dict(v) for k, v in tprt.LADDERS.items() if k != p3.COMPARISON}
    p3.register_comparison()
    assert {
        k: dict(v) for k, v in tprt.COVARIATE_SETS.items() if k != p3.COMPARISON
    } == before_covariates
    assert {k: dict(v) for k, v in tprt.LADDERS.items() if k != p3.COMPARISON} == before_ladders


def test_phase2b_comparisons_are_still_exactly_the_pre_registered_ones():
    assert set(tprt.COVARIATE_SETS["peak_vs_rest_of_host"]) == {
        "host_subfamily",
        "alignment_confidence",
        "local_opportunity",
    }
    assert set(tprt.COVARIATE_SETS["nested_vs_unnested"]) == {
        "truncation",
        "allele_frequency",
        "local_opportunity",
    }
    assert list(tprt.LADDERS["peak_vs_rest_of_host"])[-1] == "no_strata"


def test_the_phase3_ladder_excludes_allele_frequency():
    """Carrier count is algebraically tied to the outcome, so it cannot be a covariate."""
    p3.register_comparison()
    for level, covariates in tprt.LADDERS[p3.COMPARISON].items():
        assert "allele_frequency" not in covariates, level
        assert "allele_freq" not in covariates, level


# --------------------------------------------------------------------------
# the comparison itself
# --------------------------------------------------------------------------


def _row(chrom: str, pos: int, **extra: str) -> dict[str, str]:
    row = {
        "chrom": chrom,
        "pos": str(pos),
        "family": "ALU",
        "nested_in_alu_host": "1",
        "consensus_offset": "133",
        "consensus_mapping": "unambiguous",
        "consensus_match_name": "AluSx1",
        "host_len": "300",
        "host_name": "AluSx1",
        "host_strand": "+",
        "host_start0": "1000000",
        "host_end0": "1000300",
        "alignment_identity": "99.0",
        "perc_resolved": "99.5",
        "not_canonical": "0",
        "site_allele_freq": "0.01",
        "tsd_len": "",
        "tsd_seq": "",
        "polya_len": "",
        "polya_seq": "",
    }
    row.update(extra)
    return row


def test_an_empty_exposure_arm_is_reported_not_silently_dropped():
    p3.register_comparison()
    rows = [
        dict(_row("chr1", 1), _tsd_length_class="de_novo_range_5_27", _cross_half_supported="1"),
        dict(_row("chr1", 2), _tsd_length_class="de_novo_range_5_27", _cross_half_supported="0"),
    ]
    result = p3.run_recurrence_comparison(rows, seed=1, permutations=20)
    block = result["exposures"]["tsd_in_de_novo_range_5_27"]
    assert block["status"] == "not_estimable_empty_arm"
    assert block["exposed_input"] == 2
    assert block["reference_input"] == 0


def test_recurrence_outcomes_treat_a_missing_genotype_as_unevaluable():
    assert p3.event_recurrence_cross_half({"_cross_half_supported": ""}) is None
    assert p3.event_recurrence_cross_half({"_cross_half_supported": "1"}) is True
    assert p3.event_recurrence_cross_half({"_cross_half_supported": "0"}) is False
    assert p3.event_recurrence_multi_genome({"_carriers_total": ""}) is None
    assert p3.event_recurrence_multi_genome({"_carriers_total": "1"}) is False
    assert p3.event_recurrence_multi_genome({"_carriers_total": "2"}) is True


def test_every_pre_declared_hallmark_is_tested_as_an_exposure():
    p3.register_comparison()
    rows = [
        dict(
            _row("chr1", i),
            _tsd_state="reference_inconsistent",
            _tsd_length_class="de_novo_range_5_27" if i % 2 else "longer_than_de_novo_28_plus",
            _polya_state="present",
            _en_state="compatible",
            _l1_en_state="absent",
            _triad_all_observed="1",
            _cross_half_supported="1" if i % 3 else "0",
            _carriers_total=str(i % 5),
        )
        for i in range(60)
    ]
    result = p3.run_recurrence_comparison(rows, seed=7, permutations=50)
    assert set(result["exposures"]) == set(tprt.EVENTS)
    for label, block in result["exposures"].items():
        assert block["status"] in {"estimated", "not_estimable_empty_arm"}, label


def test_bh_correction_is_applied_across_exposures():
    p3.register_comparison()
    rows = [
        dict(
            _row("chr1", i),
            _tsd_state="reference_inconsistent",
            _tsd_length_class="de_novo_range_5_27" if i % 2 else "longer_than_de_novo_28_plus",
            _polya_state="present",
            _cross_half_supported="1" if i % 3 else "0",
        )
        for i in range(80)
    ]
    result = p3.run_recurrence_comparison(rows, seed=3, permutations=100)
    for block in result["exposures"].values():
        if block["status"] != "estimated":
            continue
        for res in block["outcomes"].values():
            raw = res["permutation_p"]
            adj = res["permutation_p_bh_across_exposures"]
            if raw is None or adj is None:
                continue
            assert adj >= raw - 1e-12, "BH adjustment can only raise a p-value"


# --------------------------------------------------------------------------
# prohibited language
# --------------------------------------------------------------------------


def test_the_guard_catches_a_planted_prohibited_key():
    payload = {"comparison": {"recall_of_calls": 0.9}}
    assert p3.forbidden_metric_keys(payload) == ["/comparison/recall_of_calls"]


def test_the_guard_reaches_into_lists():
    payload = {"outcomes": [{"name": "x"}, {"precision": 1}]}
    assert p3.forbidden_metric_keys(payload) == ["/outcomes[1]/precision"]


def test_the_guard_allows_prose_that_disclaims_a_metric():
    """A disclaimer is allowed to name what is not being claimed."""
    payload = {
        "note": "this is not a precision or recall result",
        "verdict": "genotype_concordance_only",
    }
    assert p3.forbidden_metric_keys(payload) == []


def test_the_guard_covers_every_prohibited_term():
    for term in p3.FORBIDDEN_METRIC_TERMS:
        assert p3.forbidden_metric_keys({f"x_{term}_y": 1})


def test_no_identifier_in_the_module_is_named_after_a_prohibited_metric():
    """Identifiers are checked separately from prose, via the parsed AST."""
    tree = ast.parse((SCRIPTS / "score_genotype_concordance.py").read_text())
    banned = set(p3.FORBIDDEN_METRIC_TERMS)
    offenders = []
    for node in ast.walk(tree):
        names = []
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            names.append(node.id)
        elif isinstance(node, ast.keyword) and node.arg:
            names.append(node.arg)
        for name in names:
            lowered = name.lower()
            for term in banned:
                if term in lowered and "guard" not in lowered:
                    offenders.append(name)
    assert offenders == [], f"prohibited metric naming in identifiers: {offenders}"


def test_the_only_permitted_subject_term_is_genotype_concordance():
    """The phase's subject label is fixed, so the report cannot drift vocabulary."""
    tree = ast.parse((SCRIPTS / "score_genotype_concordance.py").read_text())
    subjects = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "SUBJECT":
                    subjects.append(node.value.value)
    assert subjects == ["genotype_concordance"]
    assert p3.SUBJECT == "genotype_concordance"


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------


def _synthetic_cohort(n: int = 120) -> list[dict[str, str]]:
    rows = []
    for i in range(n):
        rows.append(
            dict(
                _row(f"chr{(i % 22) + 1}", 1_000_000 + i * 37),
                _tsd_state="reference_inconsistent",
                _tsd_exact="0",
                _tsd_length_class=(
                    "de_novo_range_5_27" if i % 3 == 0 else "longer_than_de_novo_28_plus"
                ),
                _polya_state="present" if i % 4 else "absent",
                _en_state="compatible" if i % 5 == 0 else "incompatible",
                _l1_en_state="absent",
                _triad_all_observed="1" if i % 6 == 0 else "0",
                _in_peak_window="1" if i % 2 == 0 else "0",
            )
        )
    return rows


def _install_fake_bcftools(monkeypatch, n_samples: int = 12, n_carriers_every: int = 3):
    samples = [f"HG{i:05d}" for i in range(n_samples)]

    def fake_bcftools(args):
        if args[0] == "query" and args[1] == "-l":
            return "\n".join(samples) + "\n"
        if args[0] == "query" and args[1] == "-f" and "%POS\n" in args[2]:
            return ""
        if args[0] == "query" and args[1] == "-f" and "GT" in args[2]:
            lines = []
            for i in range(0, 40):
                gts = []
                for j in range(n_samples):
                    gts.append("0/1" if (i % n_carriers_every == 0 and j < 2) else "0/0")
                lines.append(f"chr1\t{1_000_000 + i * 37}\t" + "\t".join(gts))
            return "\n".join(lines) + "\n"
        if args[0] == "view" and args[1] == "-h":
            return "##contig=<ID=chr1,length=248956422>\n"
        raise AssertionError(f"unexpected bcftools call: {args[:3]}")

    monkeypatch.setattr(p3, "_bcftools", fake_bcftools)
    monkeypatch.setattr(p3, "sample_names", lambda path: samples)
    return samples


def _write_cohort(path: Path, rows) -> None:
    """Write only the cohort-visible columns; the `_`-prefixed keys are added later."""
    fields = list(_row("chr1", 1).keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def test_main_runs_end_to_end_on_a_synthetic_cohort(tmp_path, monkeypatch):
    _install_fake_bcftools(monkeypatch)
    cohort = tmp_path / "cohort.csv"
    _write_cohort(cohort, _synthetic_cohort())

    rc = p3.main(
        [
            "--cohort", str(cohort),
            "--manifest", str(tmp_path / "missing_manifest.json"),
            "--alu-consensus", str(tprt.DEFAULT_ALU_CONSENSUS),
            "--genotype-bcf", str(tmp_path / "geno.bcf"),
            "--site-bcf", str(tmp_path / "geno.bcf"),
            "--cross-method-bcf", str(tmp_path / "nope.vcf.gz"),
            "--outdir", str(tmp_path / "out"),
            "--permutations", "50",
            "--no-reference",
        ]
    )
    assert rc == 0
    report = json.loads((tmp_path / "out" / "genotype_concordance.json").read_text())
    assert report["phase"] == "3"
    assert report["subject"] == "genotype_concordance"
    assert set(report["gates"]) == {"0c", "3a", "3b", "3c", "3d"}
    assert report["gates"]["3d"]["verdict"] == "cross_method_comparator_unavailable"
    assert report["recurrence"]["analysis_sites"] == 120
    assert report["accumulation"]["final_unique_sites"] > 0
    assert (tmp_path / "out" / "genotype_concordance.md").exists()
    assert (tmp_path / "out" / "phase3_site_genotypes.csv").exists()


def test_a_degenerate_run_still_writes_a_report(tmp_path, monkeypatch):
    """An empty cohort must produce a readable report, not a traceback.

    A run with no evaluable sites is the run most likely to be silently
    discarded, so the writer is required to survive it.
    """
    _install_fake_bcftools(monkeypatch)
    cohort = tmp_path / "empty.csv"
    fields = list(_row("chr1", 1).keys())
    with cohort.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
    rc = p3.main(
        [
            "--cohort", str(cohort),
            "--manifest", str(tmp_path / "m.json"),
            "--alu-consensus", str(tprt.DEFAULT_ALU_CONSENSUS),
            "--genotype-bcf", str(tmp_path / "geno.bcf"),
            "--site-bcf", str(tmp_path / "geno.bcf"),
            "--cross-method-bcf", str(tmp_path / "nope.vcf.gz"),
            "--outdir", str(tmp_path / "out"),
            "--permutations", "20",
            "--no-reference",
        ]
    )
    assert rc == 0
    report = json.loads((tmp_path / "out" / "genotype_concordance.json").read_text())
    assert report["recurrence"]["analysis_sites"] == 0
    assert report["accumulation"]["mean_carriers_per_unique_site"] is None
    assert "n/a" in (tmp_path / "out" / "genotype_concordance.md").read_text()


def test_the_end_to_end_report_contains_no_prohibited_key(tmp_path, monkeypatch):
    _install_fake_bcftools(monkeypatch)
    cohort = tmp_path / "cohort.csv"
    _write_cohort(cohort, _synthetic_cohort())
    p3.main(
        [
            "--cohort", str(cohort),
            "--manifest", str(tmp_path / "m.json"),
            "--alu-consensus", str(tprt.DEFAULT_ALU_CONSENSUS),
            "--genotype-bcf", str(tmp_path / "geno.bcf"),
            "--site-bcf", str(tmp_path / "geno.bcf"),
            "--cross-method-bcf", str(tmp_path / "nope.vcf.gz"),
            "--outdir", str(tmp_path / "out"),
            "--permutations", "50",
            "--no-reference",
        ]
    )
    report = json.loads((tmp_path / "out" / "genotype_concordance.json").read_text())
    assert p3.forbidden_metric_keys(report) == []


def test_the_site_table_lists_carrier_names_so_the_counts_can_be_audited(
    tmp_path, monkeypatch
):
    samples = _install_fake_bcftools(monkeypatch)
    annotated = [
        dict(_row("chr1", 1_000_000), _carriers_total="2", _carriers_half_a="1",
             _carriers_half_b="1", _cross_half_supported="1", _in_peak_window="1"),
    ]
    mask = np.zeros(len(samples), dtype=bool)
    mask[0] = True
    mask[5] = True
    halves = p3.split_halves(samples)
    path = tmp_path / "sites.csv"
    written = p3.write_site_table(
        path, annotated, samples, halves, {("chr1", 1_000_000): mask}
    )
    assert written == 1
    rows = list(csv.DictReader(path.open()))
    assert rows[0]["carriers_total"] == "2"
    assert rows[0]["cross_half_supported"] == "1"
    assert rows[0]["carriers_half_a_samples"] == samples[0]
    assert rows[0]["carriers_half_b_samples"] == samples[5]


def test_the_shipped_cohort_has_the_columns_phase3_reads():
    """Slow: touches the real cohort file if it is present."""
    cohort = tprt.DEFAULT_COHORT
    if not cohort.exists():
        pytest.skip(f"{cohort} not built")
    with cohort.open() as fh:
        header = next(csv.reader(fh))
    required = {"chrom", "pos", "family", "nested_in_alu_host", "consensus_offset"}
    assert required <= set(header)


def test_bcftools_is_available_for_the_real_runs():
    """Phase 3's production path needs bcftools; fail loudly in CI if it is gone."""
    if p3.DEFAULT_GENOTYPE_BCF.exists():
        assert subprocess.run(["bcftools", "--version"], capture_output=True).returncode == 0
    else:
        pytest.skip("genotype BCF not present")


def test_the_prohibited_terms_list_matches_the_gate_0c_tier_language():
    """Every term the plan forbids for this tier must be in the guard list."""
    for term in ("precision", "recall", "accuracy", "validated"):
        assert term in p3.FORBIDDEN_METRIC_TERMS
    gate = p3.gate_0c_independence(["HG03086"])
    assert gate["tier"] == "restricted_concordance_only"
    assert "no detection-performance metric is licensed" in gate["consequence"]


# --------------------------------------------------------------------------
# Gate 3d overlap: a record's length decides how far back to look
# --------------------------------------------------------------------------


def _gate_3d(monkeypatch, rows, nested, control, tolerance=50):
    """Drive gate 3d off a canned bcftools table: (chrom, pos, svtype, svlen)."""
    text = "".join("\t".join(str(f) for f in row) + "\n" for row in rows)
    monkeypatch.setattr(p3, "_bcftools", lambda args: text)
    return p3.gate_3d_cross_method_floor(
        Path("canned.bcf"), nested, control, tolerance=tolerance
    )


def test_a_record_longer_than_300bp_that_spans_the_site_is_still_counted(
    monkeypatch,
):
    """The old prefilter assumed no record exceeded 300 bp.

    It bisected on start over `[pos - (tol + 300), pos + tol + 300]`, so a 2 kb
    SVA whose span starts 1 kb to the left of the site was never examined even
    though it plainly covers it. Gate 3d reports a *coverage* measurement, so
    every miss biases it downward -- in the gate whose job is to establish that
    the comparator cannot see these sites, a downward bias makes the conclusion
    look better founded than it is.
    """
    pos, tol, length = 1_000_000, 50, 2000
    start = pos - tol - 1000
    gate = _gate_3d(
        monkeypatch,
        [("chr1", start, "SVA", length)],
        nested=[("chr1", pos)],
        control=[("chr1", 5_000_000)],
        tolerance=tol,
    )
    assert gate["nested_arm"]["covered"] == 1
    assert gate["nested_arm"]["fraction_covered"] == pytest.approx(1.0)


def test_a_site_outside_every_record_span_is_not_covered(monkeypatch):
    """The fix must not over-match: widening the scan cannot cover everything."""
    gate = _gate_3d(
        monkeypatch,
        [("chr1", 1_000, "ALU", 300)],
        nested=[("chr1", 9_000_000)],
        control=[("chr1", 1_000)],
        tolerance=50,
    )
    assert gate["nested_arm"]["covered"] == 0
    assert gate["control_arm_non_nested_alu"]["covered"] == 1


def test_coverage_respects_the_proximity_tolerance(monkeypatch):
    """A site 49 bp past the end is covered; 51 bp past is not."""
    start, tol, length = 1_000, 50, 300
    end = start + length
    rows = [("chr1", start, "ALU", length)]
    inside = _gate_3d(monkeypatch, rows, [("chr1", end + tol - 1)], [], tolerance=tol)
    outside = _gate_3d(monkeypatch, rows, [("chr1", end + tol + 1)], [], tolerance=tol)
    assert inside["nested_arm"]["covered"] == 1
    assert outside["nested_arm"]["covered"] == 0


def test_a_missing_svlen_is_treated_as_a_point_not_as_unbounded(monkeypatch):
    """SVLEN "." means unknown length; it must not become a whole-chromosome span."""
    gate = _gate_3d(
        monkeypatch,
        [("chr1", 1_000, "INS", ".")],
        nested=[("chr1", 1_040)],
        control=[("chr1", 1_000_500)],
        tolerance=50,
    )
    assert gate["nested_arm"]["covered"] == 1      # 40 bp from the point
    assert gate["control_arm_non_nested_alu"]["covered"] == 0


# --------------------------------------------------------------------------
# Sample order: masks must live in the order the halves are defined on
# --------------------------------------------------------------------------


def test_an_already_sorted_sample_list_is_left_alone():
    """The shipped 908-sample BCF is stored sorted, so this must be a no-op."""
    samples = [f"S{i:03d}" for i in range(10)]
    mask = np.array([i % 2 == 0 for i in range(10)])
    out_samples, out_carriers, info = p3.carriers_in_sorted_order(
        samples, {("chr1", 1): mask}
    )
    assert out_samples == samples
    assert info["reordered"] is False
    assert np.array_equal(out_carriers[("chr1", 1)], mask)


def test_an_unsorted_sample_list_permutes_the_masks_to_match():
    """Half membership is by position in the *sorted* order, so masks must be.

    `read_site_genotypes` builds masks in the BCF's own sample order, while
    `split_halves` assigns halves by position in the sorted order and
    `accumulation_curve` walks the sorted order while indexing masks by column.
    Without the permutation, an unsorted cohort attributes every carrier to the
    wrong half -- silently, with no gate to catch it.
    """
    samples = ["S002", "S000", "S003", "S001"]
    # File order: index 0 -> S002, 1 -> S000, 2 -> S003, 3 -> S001.
    file_order_mask = np.array([True, False, False, True])
    sorted_samples, carriers, info = p3.carriers_in_sorted_order(
        samples, {("chr1", 7): file_order_mask}
    )
    assert info["reordered"] is True
    assert sorted_samples == ["S000", "S001", "S002", "S003"]
    # Sorted order is S000, S001, S002, S003 -> file indices 1, 3, 0, 2.
    assert np.array_equal(
        carriers[("chr1", 7)], np.array([False, True, True, False])
    )


def test_permuted_masks_land_in_the_half_the_sample_belongs_to():
    """End to end: the carrier must count toward the half holding its name."""
    samples = ["S002", "S000", "S003", "S001"]
    only_S002_carries = np.array([True, False, False, False])  # index 0 == S002
    sorted_samples, carriers, _ = p3.carriers_in_sorted_order(
        samples, {("chr1", 1): only_S002_carries}
    )
    halves = p3.split_halves(sorted_samples)
    rec = p3.site_recurrence([("chr1", 1)], carriers, halves)[("chr1", 1)]
    # S002 sits at sorted index 2, which is half A.
    assert rec["carriers_total"] == 1
    assert rec["carriers_half_a"] == 1
    assert rec["carriers_half_b"] == 0
    assert rec["cross_half_supported"] is False


def test_the_accumulation_curve_and_the_halves_agree_on_one_sample_order():
    """`order` and the mask columns must be the same ordering, by construction."""
    samples = ["S002", "S000", "S003", "S001"]
    mask = np.array([True, False, True, False])
    sorted_samples, carriers, _ = p3.carriers_in_sorted_order(
        samples, {("chr1", i): mask for i in range(4)}
    )
    halves = p3.split_halves(sorted_samples)
    curve = p3.accumulation_curve(
        [("chr1", i) for i in range(4)], carriers, halves["order"]
    )
    assert curve["genome_order"] == "sorted_sample_name"
    assert curve["n_genomes"] == len(sorted_samples)
    assert curve["final_unique_sites"] == 4
