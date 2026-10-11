"""Tests for scripts/gene_orientation_summary.py.

Pure-function coverage of the pieces the sense/antisense result rests on: the Wilson interval, site
grouping (single linkage per chromosome and family), per-site collapse (a site whose member calls
disagree on ORIENT has no orientation), the register/neighbour diagnostics, the sense-fraction table
and the LINE-1 strand-asymmetry test. ``load_calls`` is exercised on a two-record VCF written to tmp_path.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "gene_orientation_summary.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("gene_orientation_summary", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gene_orientation_summary"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


gos = _load_module()


def _calls(rows):
    base = {"sample": "S1", "chrom": "chr1", "MEIFAMILY": "ALU", "ORIENT": "+", "CSQ": "intron_variant",
            "gene_l": [], "geneid_l": [], "strand_l": [], "has_updown": False, "strands": frozenset()}
    out = []
    for r in rows:
        row = {**base, **r}
        row["strands"] = frozenset(row["strand_l"])
        out.append(row)
    return pd.DataFrame(out)


def test_wilson_matches_a_known_value_and_edge_cases():
    lo, hi = gos.wilson(50, 100)
    assert lo == pytest.approx(0.4038, abs=1e-3) and hi == pytest.approx(0.5962, abs=1e-3)
    assert gos.wilson(0, 10)[0] == pytest.approx(0.0, abs=1e-12)
    assert all(np.isnan(v) for v in gos.wilson(0, 0))


def test_split_list_handles_missing_and_empty_values():
    assert gos.split_list("A,B,,C") == ["A", "B", "C"]
    assert gos.split_list(None) == [] and gos.split_list(float("nan")) == [] and gos.split_list("") == []


def test_assign_sites_is_single_linkage_per_chromosome_and_family():
    calls = _calls([
        {"pos": 100}, {"pos": 140}, {"pos": 190},          # chain: each gap <= 50 -> one site
        {"pos": 300},                                       # gap 110 -> new site
        {"pos": 100, "chrom": "chr2"},                      # other chromosome
        {"pos": 100, "MEIFAMILY": "LINE1"},                 # other family at the same coordinate
    ])
    out = gos.assign_sites(calls, merge_bp=50)
    assert out.site.nunique() == 4
    chain = out[(out.chrom == "chr1") & (out.MEIFAMILY == "ALU") & (out.pos <= 190)]
    assert chain.site.nunique() == 1


def test_assign_sites_ignores_families_outside_the_three_mei_classes():
    calls = _calls([{"pos": 100}, {"pos": 120, "MEIFAMILY": "OTHER"}])
    assert len(gos.assign_sites(calls, 50)) == 1


def test_collapse_sites_drops_orientation_when_members_disagree_and_counts_samples():
    calls = gos.assign_sites(_calls([
        {"pos": 100, "sample": "A", "ORIENT": "+", "geneid_l": ["G1"], "strand_l": ["+"]},
        {"pos": 110, "sample": "B", "ORIENT": "-", "geneid_l": ["G1"], "strand_l": ["+"]},
        {"pos": 500, "sample": "A", "ORIENT": "+", "geneid_l": ["G2"], "strand_l": ["-"]},
        {"pos": 510, "sample": "B", "ORIENT": "+", "geneid_l": ["G2"], "strand_l": ["-"]},
    ]), 50)
    sites = gos.collapse_sites(calls).sort_values("pos").reset_index(drop=True)
    assert sites.n_samples.tolist() == [2, 2]
    assert pd.isna(sites.orient[0]) and sites.orient[1] == "+"
    assert sites.geneid.tolist() == ["G1", "G2"]


def test_inspect_annotation_counts_register_and_neighbour_problems():
    calls = _calls([
        {"gene_l": ["A", "B"], "geneid_l": ["g1", "g2", "g3"], "strand_l": ["+", "+", "-"], "has_updown": True},
        {"gene_l": ["A"], "geneid_l": ["g1"], "strand_l": ["+"]},
        {"gene_l": [], "geneid_l": [], "strand_l": [], "CSQ": "intergenic_variant"},
        {"gene_l": ["A", "B"], "geneid_l": ["g1", "g2"], "strand_l": ["+", "-"], "has_updown": False},
    ])
    out = gos.inspect_annotation(calls)
    assert out["gene_vs_geneid_length_mismatch"] == 1
    assert out["strand_vs_geneid_length_mismatch"] == 0
    assert out["intronic_calls"] == 3
    assert out["intronic_mixed_strand"] == 2
    assert out["intronic_mixed_strand_with_updown_term"] == 1
    assert out["intronic_mixed_strand_without_updown_term"] == 1


def _sites(rows):
    base = {"family": "ALU", "chrom": "chr1", "n_samples": 1, "orient": "+", "csq": "intron_variant",
            "has_updown": False, "strands": frozenset({"+"})}
    return pd.DataFrame([{**base, **r} for r in rows])


def test_analysed_sites_keeps_only_unambiguous_intronic_sites_and_defines_sense():
    sites = _sites([
        {"orient": "+", "strands": frozenset({"+"})},                    # sense
        {"orient": "-", "strands": frozenset({"+"})},                    # antisense
        {"orient": None, "strands": frozenset({"+"})},                   # orientation undefined
        {"orient": "+", "strands": frozenset({"+", "-"})},               # mixed strand
        {"orient": "+", "strands": frozenset({"+"}), "csq": "intergenic_variant"},
        {"orient": "-", "strands": frozenset({"-"}), "n_samples": 3},    # sense, shared
    ])
    out = gos.analysed_sites(sites)
    assert len(out) == 3
    assert out.sense.tolist() == [True, False, True]
    assert out.sharing.tolist() == ["private", "private", "shared"]


def test_sense_table_counts_match_the_input_for_every_stratum():
    sites = _sites([{"orient": "+", "strands": frozenset({"+"})}] * 3
                   + [{"orient": "-", "strands": frozenset({"+"})}] * 7)
    table = gos.sense_table(gos.analysed_sites(sites))
    row = table[(table.set == "A_single_strand_intronic") & (table.family == "ALU")
                & (table.gene_strand == "all") & (table.sharing == "all")].iloc[0]
    assert (row.n_sites, row.n_sense) == (10, 3)
    assert row.sense_frac == pytest.approx(0.3)
    assert row["p_vs_0.5"] == pytest.approx(0.3438, abs=1e-3)


def test_strand_asymmetry_detects_a_planted_difference_and_is_seed_reproducible():
    rows = ([{"family": "LINE1", "orient": "+", "strands": frozenset({"+"})}] * 10    # + genes: 10 sense
            + [{"family": "LINE1", "orient": "-", "strands": frozenset({"+"})}] * 30  # + genes: 30 antisense
            + [{"family": "LINE1", "orient": "-", "strands": frozenset({"-"})}] * 20  # - genes: 20 sense
            + [{"family": "LINE1", "orient": "+", "strands": frozenset({"-"})}] * 20)  # - genes: 20 antisense
    a = gos.analysed_sites(_sites(rows))
    first = gos.strand_asymmetry(a, "LINE1", n_perm=200, seed=1)
    again = gos.strand_asymmetry(a, "LINE1", n_perm=200, seed=1)
    assert first["plus_sense"] == 10 and first["plus_n"] == 40
    assert first["minus_sense"] == 20 and first["minus_n"] == 40
    assert first["diff_plus_minus"] == pytest.approx(0.25 - 0.5)
    assert first["fisher_p"] < 0.05
    assert first["perm_p"] == again["perm_p"]


def test_strand_asymmetry_survives_a_single_chromosome_callset():
    """Leave-one-chromosome-out has nothing left to compare when all sites share a chromosome."""
    rows = ([{"family": "LINE1", "orient": "+", "strands": frozenset({"+"})}] * 5
            + [{"family": "LINE1", "orient": "-", "strands": frozenset({"-"})}] * 5)
    out = gos.strand_asymmetry(gos.analysed_sites(_sites(rows)), "LINE1", n_perm=20, seed=1)
    assert out["loo_n_chrom"] == 0 and "loo_min" not in out
    assert "perm_p" in out


def test_strand_asymmetry_leave_one_chromosome_out_range_spans_the_per_chromosome_drops():
    def block(chrom, plus_sense, plus_anti, minus_sense, minus_anti):
        mk = lambda o, s: {"family": "LINE1", "chrom": chrom, "orient": o, "strands": frozenset({s})}
        return ([mk("+", "+")] * plus_sense + [mk("-", "+")] * plus_anti
                + [mk("-", "-")] * minus_sense + [mk("+", "-")] * minus_anti)
    rows = block("chr1", 2, 8, 5, 5) + block("chr2", 2, 8, 5, 5) + block("chr3", 5, 5, 5, 5)
    out = gos.strand_asymmetry(gos.analysed_sites(_sites(rows)), "LINE1", n_perm=20, seed=1)
    assert out["loo_n_chrom"] == 3
    assert out["loo_min"] <= out["diff_plus_minus"] <= out["loo_max"]


def test_strand_asymmetry_returns_early_when_one_gene_strand_is_missing():
    a = gos.analysed_sites(_sites([{"family": "SVA", "orient": "+", "strands": frozenset({"+"})}] * 3))
    out = gos.strand_asymmetry(a, "SVA", n_perm=10, seed=1)
    assert "perm_p" not in out and out["minus_n"] == 0


def test_load_calls_parses_gene_fields_and_flags_neighbour_terms(tmp_path):
    vcf = tmp_path / "S1_classifier.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tN\t<INS:ME:ALU>\t.\tPASS\tMEIFAMILY=ALU;ORIENT=+;CSQ=intron_variant;"
        "CSQ_TERMS=intron_variant,upstream_gene_variant;GENE=A;GENEID=g1,g2;GENE_STRAND=+,-\n"
        "chr1\t900\t.\tN\t<INS:ME:ALU>\t.\tPASS\tMEIFAMILY=ALU;ORIENT=-;CSQ=intergenic_variant\n")
    out = gos.load_calls(tmp_path)
    assert len(out) == 2 and set(out["sample"]) == {"S1"}
    first = out.iloc[0]
    assert first.gene_l == ["A"] and first.geneid_l == ["g1", "g2"] and first.strand_l == ["+", "-"]
    assert bool(first.has_updown) and first.strands == frozenset({"+", "-"})
    assert out.iloc[1].geneid_l == [] and not bool(out.iloc[1].has_updown)


def test_load_calls_rejects_an_empty_directory(tmp_path):
    with pytest.raises(FileNotFoundError):
        gos.load_calls(tmp_path)
