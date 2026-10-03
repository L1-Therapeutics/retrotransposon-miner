"""Reproduce the published ten-genome numbers from the real inputs.

The unit tests in `test_analyze_ten_genome_mei.py` pin the *rules* on synthetic
fixtures. This file checks the other half: that the pipeline still turns the ten
real callsets into the numbers that were actually published. A rule can be
pinned and still be applied wrongly, and `analysis_metrics.json` recorded hashes
of the code but never of the inputs, so nothing would have caught that.

Cheap checks (published CSVs and the metrics block) run by default. The checks
that must parse the 156 MB `rmsk.txt.gz` take about five minutes and are opt-in:

    RTM_REPRODUCE_SLOW=1 pytest tests/test_published_ten_genome_reproduction.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
WORKSPACE = REPO.parent
CALLSETS = WORKSPACE / "nested_analysis" / "callsets"
RESULTS = WORKSPACE / "nested_analysis" / "results_multi_sample"
RMSK = WORKSPACE / "nested_analysis" / "data" / "rmsk.txt.gz"
MANIFEST = CALLSETS / "manifest.sha256"
METRICS = RESULTS / "analysis_metrics.json"
SITES_CSV = RESULTS / "unique_sites.csv"
PRIVATE_CSV = RESULTS / "private_sites.csv"

SLOW = os.environ.get("RTM_REPRODUCE_SLOW") == "1"
requires_slow = pytest.mark.skipif(
    not SLOW,
    reason="set RTM_REPRODUCE_SLOW=1 to re-derive sites from rmsk.txt.gz (~5 min)",
)


def _require(*paths):
    missing = [p for p in paths if not p.exists()]
    if missing:
        pytest.skip(f"published artifacts or inputs absent on this host: "
                    f"{[p.name for p in missing]}")


if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import analyze_ten_genome_mei as ten  # noqa: E402


@pytest.fixture(scope="module")
def engine():
    return ten.load_engine(str(SCRIPTS))


@pytest.fixture(scope="module")
def published():
    _require(METRICS, SITES_CSV, PRIVATE_CSV, MANIFEST, CALLSETS)
    return json.loads(METRICS.read_text())


@pytest.fixture(scope="module")
def site_rows():
    _require(SITES_CSV)
    with SITES_CSV.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


@pytest.fixture(scope="module")
def ingestion(engine):
    """Re-read the ten callsets. ~5 s; no reference data needed."""
    _require(CALLSETS)
    calls, qc = ten.read_calls(engine, CALLSETS)
    return calls, qc


@pytest.fixture(scope="module")
def matched(engine, ingestion):
    """Re-derive sites with hosts assigned. Opt-in: needs rmsk."""
    _require(RMSK)
    calls, _ = ingestion
    all_calls = [c for sample in ten.SAMPLES for c in calls[sample]]
    hosts = engine.read_rmsk(RMSK)
    engine.assign_hosts(all_calls, hosts)
    return all_calls


# ---------------------------------------------------------------------------
# ingestion
# ---------------------------------------------------------------------------
def test_ingestion_reproduces_the_published_qc_block(ingestion, published):
    """Every per-sample count, chromosome census, flag and hash must be identical."""
    _, qc = ingestion
    assert qc == published["qc"]


def test_published_callset_hashes_still_match_the_ingestion_manifest(published):
    """Keeps re-verifying the hashes the analysis recorded against the manifest.

    Section 2 established these once. Nothing re-checked them afterwards, so a
    silently replaced callset would still have been analysed and reported.
    """
    recorded = {entry["sample"]: entry["sha256"] for entry in published["qc"]}
    manifest = {
        line.split()[1].rsplit("/", 1)[-1].removesuffix(".vcf"): line.split()[0]
        for line in MANIFEST.read_text().splitlines() if line.strip()
    }
    assert set(recorded) <= set(manifest)
    for sample, digest in recorded.items():
        assert digest == manifest[sample], f"{sample} callset changed since the run"


def test_ingestion_excludes_chr22_mei_slice_and_keeps_every_genome(ingestion):
    """chr22_mei.vcf is not a cohort member, and no callset loses a chromosome."""
    calls, _ = ingestion
    assert set(calls) == set(ten.SAMPLES)
    assert "chr22_mei" not in calls
    for sample in ten.SAMPLES:
        assert {c.chrom for c in calls[sample]} == set(engine_chroms())


def engine_chroms():
    return [f"chr{i}" for i in range(1, 23)] + ["chrX"]


# ---------------------------------------------------------------------------
# published table
# ---------------------------------------------------------------------------
def test_published_table_row_count_matches_the_reported_union(site_rows, published):
    union = sum(row["union"] for row in published["ten"]["layers"])
    assert len(site_rows) == union


def test_carrier_count_agrees_with_the_ten_bit_presence_bitmap(site_rows):
    """`n_carriers` and the bitmap are two encodings of one fact."""
    for row in site_rows:
        assert int(row["n_carriers"]) == row["presence_bitmap"].count("1")
        assert len(row["presence_bitmap"]) == 10
        # the writer emits Python bool repr, not 0/1
        assert (row["private"] == "True") == (row["n_carriers"] == "1")
        assert (row["private_to_sample"] != "") == (row["n_carriers"] == "1")


def test_every_genome_carries_at_least_one_unique_site(site_rows):
    for sample in ten.SAMPLES:
        assert any(row[f"present_{sample}"] == "1" for row in site_rows), sample


def test_private_table_is_exactly_the_carrier_count_one_rows(site_rows, published):
    _require(PRIVATE_CSV)
    with PRIVATE_CSV.open(encoding="utf-8") as handle:
        private = list(csv.DictReader(handle))
    expected = [row for row in site_rows if row["n_carriers"] == "1"]
    assert len(private) == len(expected)
    assert sum(row["private"] for row in published["ten"]["layers"]) == len(expected)
    assert all(row["n_carriers"] == "1" for row in private)


def test_published_nested_states_never_fold_unknown_into_antisense(site_rows):
    """`nested_unknown` must survive as its own state in the shipped table.

    Folding it into `nested_antisense` would manufacture a polarity the data
    does not have, and nothing downstream could tell.
    """
    seen = Counter(row["same_family_nested_state"] for row in site_rows)
    assert set(seen) <= set(ten.STATES) | {""}
    assert seen["nested_antisense"], "no antisense sites at all is itself suspicious"


def test_published_rows_carry_the_opportunity_caveat_verbatim(site_rows):
    for row in site_rows:
        assert row["opportunity_caveat"] == ten.CAVEAT


# ---------------------------------------------------------------------------
# matching stage (opt-in: parses the 156 MB rmsk)
# ---------------------------------------------------------------------------
@requires_slow
def test_window_sweep_reproduces_exactly(matched, engine, published):
    got = []
    for window in (0, 5, 10, 20):
        sites, collapsed = ten.dedup(engine, matched, window)
        got.append({"window": window, "union": len(sites),
                    "private": sum(len(s.samples) == 1 for s in sites),
                    "shared": sum(len(s.samples) >= 2 for s in sites),
                    "within_sample_collapsed": collapsed})
    assert got == published["sweep"]


@requires_slow
def test_ten_genome_layers_reproduce_exactly(matched, engine, published):
    sites, _ = ten.dedup(engine, matched)
    assert ten.layer_counts(sites) == published["ten"]["layers"]


@requires_slow
def test_five_genome_layers_reproduce_exactly(engine, ingestion, published):
    calls, _ = ingestion
    sites, _ = ten.dedup(engine, [c for s in ten.ORIGINAL for c in calls[s]])
    assert ten.layer_counts(sites) == published["five"]["layers"]


@requires_slow
def test_five_genome_registered_counts_reproduce(engine, ingestion, published):
    """The four counts the run gates on, re-derived rather than trusted."""
    calls, _ = ingestion
    sites, _ = ten.dedup(engine, [c for s in ten.ORIGINAL for c in calls[s]])
    nested_l1 = sum(1 for s in sites if s.representative.family == "LINE1"
                    and s.representative.nested_state.startswith("nested"))
    assert (len(sites),
            sum(len(s.samples) == 1 for s in sites),
            sum(len(s.samples) >= 2 for s in sites),
            nested_l1) == (4998, 3051, 1947, 231)
