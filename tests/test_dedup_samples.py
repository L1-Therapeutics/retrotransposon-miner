"""Tests for the five-callset unique-site MEI analysis."""
from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "dedup_samples.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("dedup_samples", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


m = _load_module()


def call(sample: str, pos: int, *, family="ALU", orient="+", host=None, nested="nested_sense"):
    c = m.Call(
        sample=sample, source_index=pos, chrom="chr1", pos=pos, family=family,
        orientation=orient, raw_nested=nested, info={}, sample_column=sample,
        same_family_host=host, match_host=host, nested_state=nested,
    )
    return c


def test_checksum_gate_compares_both_bytes_before_selecting(tmp_path):
    for sample in m.SAMPLES:
        (tmp_path / f"{sample}.vcf").write_bytes(b"identical")
        (tmp_path / f"{sample}_classifier_ge_0.997.genes.snpeff.vcf").write_bytes(b"identical")
    pairs = m.input_paths(tmp_path, tmp_path)
    assert list(pairs) == list(m.SAMPLES)
    (tmp_path / "HG03172_classifier_ge_0.997.genes.snpeff.vcf").write_bytes(b"changed")
    with pytest.raises(m.InputGateError, match="SHA-256 MISMATCH.*STOP"):
        m.input_paths(tmp_path, tmp_path)


def test_family_normalization_and_four_state_nesting_are_separate():
    assert m.normalize_family("AluY#SINE/Alu") == "ALU"
    assert m.normalize_family("L1HS#LINE/L1") == "LINE1"
    assert m.normalize_family("SVA_F#Retroposon/SVA") == "SVA"
    host = m.Host("chr1", 100, 400, "+", "AluY", "ALU")
    assert call("HG03086", 201, orient="+", host=host).nested_state == "nested_sense"
    assert call("HG03086", 201, orient="-", host=host, nested="nested_antisense").nested_state == "nested_antisense"
    assert call("HG03086", 201, orient=".", host=host, nested="nested_unknown").nested_state == "nested_unknown"
    assert call("HG03086", 201, host=None, nested="unnested").nested_state == "unnested"


def test_family_chromosome_orientation_host_and_window_define_site():
    host = m.Host("chr1", 100, 400, "+", "AluY", "ALU")
    base = call("HG03086", 200, host=host)
    match = call("HG01474", 210, host=host)
    other_orientation = call("HG01566", 200, orient="-", host=host)
    other_family = call("HG03172", 200, family="LINE1", host=host)
    other_host = call("NA18498", 200, host=m.Host("chr1", 100, 400, "+", "AluSx", "ALU"))
    sites = m.anchored_groups([base, match, other_orientation, other_family, other_host], 10, allow_same_sample=False)
    carriers = sorted(len(site.samples) for site in sites)
    assert carriers == [1, 1, 1, 2]
    assert sum(site.samples == ("HG03086", "HG01474") for site in sites) == 1


def test_no_transitive_chaining_and_private_sites_remain():
    host = m.Host("chr1", 100, 500, "+", "AluY", "ALU")
    calls = [call("HG03086", 100, host=host), call("HG01474", 110, host=host), call("HG01566", 120, host=host)]
    sites = m.anchored_groups(calls, 10, allow_same_sample=False)
    assert sorted(len(site.samples) for site in sites) == [1, 2]
    assert sum(site.anchor.pos == 120 and site.samples == ("HG01566",) for site in sites) == 1


def test_exact_and_window_sweep_boundaries():
    host = m.Host("chr1", 100, 500, "+", "AluY", "ALU")
    calls = [call("HG03086", 100, host=host), call("HG01474", 105, host=host), call("HG01566", 120, host=host)]
    assert len(m.anchored_groups(calls, 0, allow_same_sample=False)) == 3
    assert len(m.anchored_groups(calls, 5, allow_same_sample=False)) == 2
    assert len(m.anchored_groups(calls, 10, allow_same_sample=False)) == 2
    assert len(m.anchored_groups(calls, 20, allow_same_sample=False)) == 1


def test_within_sample_duplicate_records_collapse_but_remain_provenanced():
    host = m.Host("chr1", 100, 500, "+", "AluY", "ALU")
    a, b = call("HG03086", 200, host=host), call("HG03086", 205, host=host)
    a.source_ids, b.source_ids = ["a"], ["b"]
    sites, collapsed = m.unique_sites([a, b, call("HG01474", 203, host=host)], 10)
    assert collapsed == 1
    assert len(sites) == 1 and sites[0].samples == ("HG03086", "HG01474")
    assert sum(member.source_count for member in sites[0].members) == 3


def test_repeated_dedup_passes_do_not_compound_provenance_totals():
    """`run` dedups the same Call objects once per window, so the within-sample
    collapse must not write its totals back onto the caller's objects.

    Regression: the in-place version re-summed its own inflated predecessor, so a
    three-record site reported 3 -> 4 -> 5 source records and grew its id list
    "a|b|c" -> "a|b|b|c" -> "a|b|b|b|c".
    """
    host = m.Host("chr1", 100, 500, "+", "AluY", "ALU")
    a, b = call("HG03086", 200, host=host), call("HG03086", 205, host=host)
    c = call("HG01474", 203, host=host)
    a.source_ids, b.source_ids, c.source_ids = ["a"], ["b"], ["c"]
    calls = [a, b, c]
    for window in (m.PRIMARY_WINDOW, *m.SWEEP_WINDOWS):
        sites, _ = m.unique_sites(calls, window)
        if window != m.PRIMARY_WINDOW:
            continue
        assert sum(member.source_count for member in sites[0].members) == 3
        ids = [i for member in sites[0].members for i in member.source_ids]
        assert sorted(ids) == ["a", "b", "c"]
    # the caller's own records are left exactly as parsed
    assert [x.source_count for x in calls] == [1, 1, 1]
    assert [x.source_ids for x in calls] == [["a"], ["b"], ["c"]]


def test_matching_and_profile_never_read_gene_fields_or_genotypes():
    host = m.Host("chr1", 100, 400, "+", "AluY", "ALU")
    a, b = call("HG03086", 200, host=host), call("HG01474", 205, host=host)
    a.info.update({"GENE": "GENE_A", "GT": "0/0", "GQ": "99"})
    b.info.update({"GENE": "GENE_B", "GT": "1/1", "GQ": "60"})
    sites = m.anchored_groups([a, b], 10, allow_same_sample=False)
    assert len(sites) == 1
    assert len(sites[0].samples) == 2
    assert m.profile_events([], sites, "ALU", "unique") == [(99, 300)]


def test_gc_matched_probabilities_use_local_gc_and_per_family_coverage():
    host = m.Host("chr1", 100, 400, "+", "AluY", "ALU")
    sites = [m.Site(call("HG03086", 150, family="LINE1", host=host), []),
             m.Site(call("HG01474", 250, family="LINE1", host=host), [])]
    sites[0].members = [sites[0].anchor]
    sites[1].members = [sites[1].anchor]
    bundle = {
        "gc_by_site": {("chr1", 0): 0.30, ("chr1", 1): 0.80},
        "frame": {"gc": np.asarray([0.30, 0.31, 0.80]),
                  "coverage": {"LINE1": np.asarray([0.1, 0.3, 0.9]), "ALU": np.zeros(3), "SVA": np.zeros(3)},
                  "sorted_order": np.asarray([0, 1, 2])},
    }
    probabilities = m.gc_probabilities(sites, "LINE1", bundle, np.random.default_rng(3))
    assert probabilities.tolist() == pytest.approx([0.2, 0.9])


def test_bernoulli_null_uses_reproducible_resamples():
    p = np.asarray([0.25, 0.5, 0.75])
    a = m.resample_nested(p, 1000, np.random.default_rng(44))
    b = m.resample_nested(p, 1000, np.random.default_rng(44))
    assert len(a) == 1000
    assert np.array_equal(a, b)
    assert abs(float(a.mean()) - 1.5) < 0.12


def test_position_profile_reports_the_two_windows_and_uniform_expectation():
    events = [(130, 300)] * 9 + [(290, 300)] * 11 + [(30, 300)] * 20
    rows = m.position_profile(events, "S1", "ALU", "per_sample")
    linker = next(row for row in rows if row["bin_start_bp"] == 120)
    tail = next(row for row in rows if row["bin_start_bp"] == 280)
    assert linker["observed"] == 9 and tail["observed"] == 11
    assert linker["expected_uniform"] == pytest.approx(40 * 20 / 300)
    assert linker["enrichment"] > 1 and tail["enrichment"] > 1


def test_fast_fasta_reader_handles_wrapped_sequence():
    raw = b"ACGTAC\nGTACGT\n"
    # .fai tuple: sequence length, first-base byte offset, bases/line, bytes/line.
    got = m.fasta_sequence_chunk(io.BytesIO(raw), (12, 0, 6, 7), 4, 10)
    assert got == b"ACGTAC"


def _load_common():
    spec = importlib.util.spec_from_file_location(
        "nested_multi_sample_common", SCRIPT.parent / "nested_multi_sample_common.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_dedup_table_is_loadable_by_the_phase4_consumers(tmp_path):
    """Producer and consumer must agree on the columns, not just on the filename.

    Both scripts write and read `nested_analysis/results_multi_sample/unique_sites.csv`.
    This one wrote host geometry only under `same_family_host_*` names, while the
    Phase 4 loader requires `host_name`/`host_start0`/`host_end0`/`host_len`, so
    `joint_enrichment.py` and `recurrence_test.py` both died with a SchemaError
    on the real file while every unit test here still passed. Nothing exercised
    the seam between the two modules, which is exactly why it went unnoticed.

    The contract is asserted directly: every column the loader requires must be
    present in a table this script produces.
    """
    common = _load_common()
    host = m.Host("chr1", 900, 1200, "+", "AluY", "ALU")
    sites = m.anchored_groups(
        [call("HG03086", 1000, host=host), call("HG01474", 1005, host=host)],
        10,
        allow_same_sample=False,
    )
    frame = m.build_sites_table(sites, 10)
    missing = [
        column
        for column in common.REQUIRED_UNIQUE_SITE_COLUMNS
        if column not in frame.columns
    ]
    assert missing == [], f"dedup output is missing loader-required columns: {missing}"

    # And it must survive the real loader, not just a column-presence check.
    out = tmp_path / "unique_sites.csv"
    frame.to_csv(out, index=False)
    rows, report = common.load_unique_sites(out)
    assert report["rows_in_scope"] == 1
    assert rows[0]["host_name"] == "AluY"
    assert rows[0]["host_start0"] == 900
    assert rows[0]["host_end0"] == 1200
    assert rows[0]["host_len"] == 300
    # Carriers are emitted in the fixed SAMPLES order, not alphabetical order.
    assert rows[0]["carriers"] == ["HG03086", "HG01474"]


def test_the_host_geometry_columns_agree_with_the_same_family_columns(tmp_path):
    """The two spellings must describe the same host, not drift apart."""
    host = m.Host("chr1", 900, 1200, "-", "AluY", "ALU")
    sites = m.anchored_groups([call("HG03086", 1000, host=host)], 10,
                              allow_same_sample=True)
    row = m.site_row(sites[0], "US00001", 10)
    assert row["host_name"] == row["same_family_host_name"] == "AluY"
    assert row["host_strand"] == row["same_family_host_strand"] == "-"
    assert row["host_len"] == row["same_family_host_length"] == 300
    assert row["host_end0"] - row["host_start0"] == row["host_len"]


def test_a_hostless_site_keeps_empty_geometry_rather_than_a_fake_zero():
    """Unnested calls have no host; zero would look like a real 0 bp element."""
    sites = m.anchored_groups([call("HG03086", 1000, host=None)], 10,
                              allow_same_sample=True)
    row = m.site_row(sites[0], "US00001", 10)
    assert row["host_name"] == "" and row["host_len"] == ""
    assert row["host_start0"] == "" and row["host_end0"] == ""
    assert row["host_strand"] == "" and row["host_offset_5p_0based"] == ""
    # `same_family_host_length` stayed blank before this fix too; a 0 would be
    # indistinguishable from a real zero-length element to any downstream reader.


def test_the_misleading_any_mei_host_column_is_not_emitted():
    """`any_MEI_host_id_for_matching` is gone and must not come back.

    The name promises a host drawn from any MEI family, but the frozen matching
    key in `matching_audit.md` is the *same-family* host, so the column could
    only ever hold a byte-for-byte duplicate of `same_family_host_id`. It had no
    consumer and the authoritative producer never emitted it. Populating it from
    the all-family sweep would have named the wrong host -- the longest element
    of any family, which for an Alu inside an L1 is the L1 -- so the column is
    removed rather than filled.
    """
    host = m.Host("chr1", 900, 1200, "+", "AluY", "ALU")
    sites = m.anchored_groups([call("HG03086", 1000, host=host)], 10,
                              allow_same_sample=True)
    row = m.site_row(sites[0], "US00001", 10)
    assert "any_MEI_host_id_for_matching" not in row


def test_the_emitted_schema_has_no_field_whose_name_contradicts_its_value():
    """Both host vocabularies stay; only the false name is gone.

    The `same_family_host_*` spelling is what the matching audit reads and the
    plain `host_*` spelling is what `nested_multi_sample_common.load_unique_sites`
    requires. Both describe the same host, and the test above this one already
    pins that they agree. What must not survive is a third spelling claiming a
    different host.
    """
    host = m.Host("chr1", 900, 1200, "+", "AluY", "ALU")
    sites = m.anchored_groups([call("HG03086", 1000, host=host)], 10,
                              allow_same_sample=True)
    row = m.site_row(sites[0], "US00001", 10)
    for column in ("host_name", "host_start0", "host_end0", "host_strand",
                   "host_len", "same_family_host_id", "same_family_host_name",
                   "same_family_host_strand", "same_family_host_length"):
        assert column in row
    assert not [c for c in row if c.startswith("any_")]
