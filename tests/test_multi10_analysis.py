"""Tests for the ten-genome dedup core (the frozen matching rules).

These pin the rules the task specification freezes:

* match key = chromosome + family + orientation + host element + breakpoint
  within the window;
* no transitive chaining -- every additional member must match the site anchor;
* exactly one call per sample per site, so one sample contributes one carrier;
* the primary +/-10 bp window and the exact/+/-5/+/-20 sweep.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sys.path.insert(0, str(ROOT / "scripts"))
engine = _load("dedup_samples", "scripts/dedup_samples.py")
m = _load("multi10_analysis", "scripts/multi10_analysis.py")

S1, S2, S3 = m.SAMPLES[0], m.SAMPLES[1], m.SAMPLES[2]


def host(chrom="chr1", start=100, end=400, strand="+", name="AluY", family="ALU"):
    return engine.Host(chrom, start, end, strand, name, family)


def call(sample, pos, *, family="ALU", orient="+", h=None, nested="nested_sense", sid=None):
    return m.Call(
        sample=sample,
        source_index=pos,
        chrom="chr1",
        pos=pos,
        family=family,
        orientation=orient,
        raw_nested="nested",
        info={"MEIFAMILY": family, "ORIENT": orient, "NESTED": "nested"},
        same_family_host=h,
        match_host=h,
        nested_state=nested,
        source_ids=[sid or f"{sample}:{pos}"],
    )


def test_calls_carry_info_for_known_mei_exclusions():
    # mei_reference_opportunity.known_exclusions reads call.info['KNOWNMEI']
    c = call(S1, 150)
    c.info["KNOWNMEI"] = "True"
    assert c.info["KNOWNMEI"].lower() == "true"
    assert c.pos0 == 149


def test_sweep_and_primary_window_match_the_spec():
    assert m.PRIMARY_WINDOW == 10
    assert m.SWEEP_WINDOWS == (0, 5, 10, 20)
    assert len(m.SAMPLES) == 10
    assert m.ORIGINAL_FIVE == m.SAMPLES[:5]


def test_host_containment_uses_start_lt_pos_le_end():
    h = host(start=100, end=200)
    # VCF POS is 1-based; a host contains POS when start0 < POS <= end0
    assert not (h.start0 < 100 <= h.end0)
    assert h.start0 < 101 <= h.end0
    assert h.start0 < 200 <= h.end0
    assert not (h.start0 < 201 <= h.end0)


def test_match_key_uses_family_orientation_and_host():
    h = host()
    a = call(S1, 150, h=h)
    b = call(S2, 152, h=h)
    assert a.match_key == b.match_key
    # different family, orientation or host must not collide
    assert call(S1, 150, family="LINE1", h=host(family="LINE1")).match_key != a.match_key
    assert call(S1, 150, orient="-", h=h).match_key != a.match_key
    assert call(S1, 150, h=host(name="AluSx")).match_key != a.match_key
    # an unnested call gets the explicit no-host key
    assert call(S1, 150, h=None, nested="unnested").match_key != a.match_key


def test_no_transitive_chaining():
    # A@100 and B@110 pair (distance 10).  C@120 is within 10 of B but 20 of A,
    # so C must NOT join the site.
    h = host(start=0, end=500)
    calls = [call(S1, 100, h=h), call(S2, 110, h=h), call(S3, 120, h=h)]
    sites = m.anchored_groups(calls, 10, allow_same_sample=False)
    sizes = sorted(len(s.members) for s in sites)
    assert sizes == [1, 2]
    assert sorted(len(s.samples) for s in sites) == [1, 2]


def test_one_sample_contributes_one_carrier():
    h = host(start=0, end=500)
    calls = [call(S1, 100, h=h), call(S1, 108, h=h), call(S2, 104, h=h)]
    # unique_sites collapses S1's two records to one representative FIRST, so
    # S1 cannot appear twice in a site and contributes exactly one carrier.
    sites, collapsed = m.unique_sites(calls, 10)
    assert collapsed == 1
    assert len(sites) == 1
    site = sites[0]
    assert site.samples == (S1, S2)
    assert sum(1 for member in site.members if member.sample == S1) == 1
    # the collapsed-away record is retained as provenance, not deleted
    assert sum(member.source_count for member in site.members) == 3


def test_anchored_groups_alone_refuses_a_second_same_sample_call():
    # without the within-sample collapse, a second S1 call in the window must
    # not join the site S1 already belongs to
    h = host(start=0, end=500)
    calls = [call(S1, 100, h=h), call(S1, 108, h=h), call(S2, 104, h=h)]
    sites = m.anchored_groups(calls, 10, allow_same_sample=False)
    shared = [s for s in sites if len(s.samples) >= 2]
    assert len(shared) == 1
    assert sum(1 for member in shared[0].members if member.sample == S1) == 1


def test_collapse_within_sample_keeps_provenance_and_is_idempotent():
    h = host(start=0, end=500)
    a, b = call(S1, 100, h=h, sid="a"), call(S1, 105, h=h, sid="b")
    c = call(S2, 103, h=h, sid="c")
    calls = [a, b, c]
    sites, collapsed = m.unique_sites(calls, 10)
    assert collapsed == 1
    assert len(sites) == 1 and sites[0].samples == (S1, S2)
    assert sum(member.source_count for member in sites[0].members) == 3
    ids = sorted(i for member in sites[0].members for i in member.source_ids)
    assert ids == ["a", "b", "c"]
    # repeated passes over the same Call objects must not compound the totals.
    # Window 0 groups nothing (no two records share a coordinate), so only the
    # windows that actually form the three-record site are checked.
    for window in (5, 10, 20):
        again, _ = m.unique_sites(calls, window)
        assert sum(member.source_count for member in again[0].members) == 3
        ids = sorted(i for member in again[0].members for i in member.source_ids)
        assert ids == ["a", "b", "c"]
    # the exact-match sweep still separates the records
    exact, _ = m.unique_sites(calls, 0)
    assert len(exact) == 3
    # and the caller's records are untouched
    assert [x.source_count for x in calls] == [1, 1, 1]
    assert [x.source_ids for x in calls] == [["a"], ["b"], ["c"]]


def test_window_sweep_boundaries():
    h = host(start=0, end=500)
    calls = [call(S1, 100, h=h), call(S2, 105, h=h), call(S3, 120, h=h)]
    assert len(m.anchored_groups(calls, 0, allow_same_sample=False)) == 3
    assert len(m.anchored_groups(calls, 5, allow_same_sample=False)) == 2
    assert len(m.anchored_groups(calls, 10, allow_same_sample=False)) == 2
    assert len(m.anchored_groups(calls, 20, allow_same_sample=False)) == 1


def test_layers_split_private_and_shared():
    h = host(start=0, end=500)
    calls = [
        call(S1, 100, h=h),
        call(S2, 104, h=h),  # shares with S1
        call(S3, 300, h=h),  # private
    ]
    sites, _ = m.unique_sites(calls, 10)
    private = [s for s in sites if len(s.samples) == 1]
    shared = [s for s in sites if len(s.samples) >= 2]
    assert len(private) == 1 and private[0].samples == (S3,)
    assert len(shared) == 1 and shared[0].samples == (S1, S2)
    frame = m.build_sites_table(sites, 10)
    # sites come out ordered by position, so carrier counts follow position
    assert list(frame["representative_pos"]) == [100, 300]
    assert list(frame["n_carriers"]) == [2, 1]
    assert frame["presence_bitmap"].str.len().eq(len(m.SAMPLES)).all()
    assert set(frame.loc[frame["n_carriers"] == 1, "private_to_sample"]) == {S3}


def test_unknown_orientation_is_never_folded_into_antisense():
    h = host(start=0, end=500)
    a = call(S1, 100, orient=".", h=h, nested="nested_unknown")
    b = call(S2, 104, orient="+", h=h, nested="nested_sense")
    c = call(S3, 104, orient="-", h=h, nested="nested_antisense")
    states = {x.nested_state for x in (a, b, c)}
    assert states == {"nested_unknown", "nested_sense", "nested_antisense"}
    assert "nested_unknown" not in {"nested_antisense"}


def test_chrY_is_never_assigned_a_host():
    # assign_hosts only considers primary chromosomes
    y_call = m.Call(
        sample=S1, source_index=0, chrom="chrY", pos=1000, family="ALU",
        orientation="+", raw_nested="nested",
    )
    m.assign_hosts([y_call], {("chrY", "ALU"): [host(chrom="chrY", start=900, end=1100)]})
    assert y_call.same_family_host is None
    assert y_call.nested_state == "unnested"


def test_longest_host_wins_then_leftmost():
    short = host(start=100, end=200, name="short")
    long = host(start=100, end=400, name="long")
    left = host(start=50, end=400, name="left")
    c = call(S1, 150, h=None)
    m.assign_hosts([c], {("chr1", "ALU"): [short, long]})
    assert c.same_family_host.name == "long"
    c2 = call(S1, 150, h=None)
    m.assign_hosts([c2], {("chr1", "ALU"): [long, left]})
    assert c2.same_family_host.name == "left"
