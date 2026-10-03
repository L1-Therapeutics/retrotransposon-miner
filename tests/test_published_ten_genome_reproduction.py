"""Reproduce the published ten-genome numbers from the real inputs.

The unit tests in `test_analyze_ten_genome_mei.py` pin the *rules* on synthetic
fixtures. This file checks the other half: that the pipeline still turns the ten
real callsets into the numbers that were actually published. A rule can be
pinned and still be applied wrongly, and `analysis_metrics.json` recorded hashes
of the code but never of the inputs, so nothing would have caught that.

Cheap checks (published CSVs, the metrics block, and the rendered report) run by
default. The checks that must parse the 156 MB `rmsk.txt.gz` take about five
minutes and are opt-in:

    RTM_REPRODUCE_SLOW=1 pytest tests/test_published_ten_genome_reproduction.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter, defaultdict
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
SUMMARY_MD = RESULTS / "analysis_summary.md"
ACCUMULATION_CSV = RESULTS / "accumulation_curve.csv"
PER_SAMPLE_CSV = RESULTS / "per_sample_calls.csv"
PROFILES_CSV = RESULTS / "position_profiles.csv"
PROJECTION_CSV = RESULTS / "l1_consensus_projection.csv"
HOST_EVIDENCE_CSV = RESULTS / "host_parent_evidence.csv"

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
import mei_reference_opportunity as opportunity  # noqa: E402


def _phase4_status_path():
    return RESULTS / "phase4_status.json"


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


# ---------------------------------------------------------------------------
# the written report against the machine-readable record
#
# `analysis_summary.md` is the artifact a person actually reads, and it is
# rendered separately from `analysis_metrics.json`. A number that drifts between
# them would be invisible to every other check here, so each row is rebuilt from
# the metrics through the same formatter the writer uses and required to appear
# verbatim.
# ---------------------------------------------------------------------------
def report_text():
    _require(SUMMARY_MD)
    return SUMMARY_MD.read_text(encoding="utf-8")


def row(cells):
    return "| " + " | ".join(cells) + " |"


def test_qc_table_is_rendered_from_the_published_qc_block(published):
    text = report_text()
    for entry in published["qc"]:
        assert row([entry["sample"], ten.fmt(entry["calls"]), ten.fmt(entry["ALU"]),
                    ten.fmt(entry["LINE1"]), ten.fmt(entry["SVA"]), ten.fmt(entry["chrY"]),
                    "PASS", entry["schema"], entry["note"]]) in text


def test_delta_table_is_rendered_from_the_published_delta(published):
    text = report_text()
    for entry in published["delta"]:
        assert row([entry["headline"], ten.fmt(entry["five"]), ten.fmt(entry["ten"]),
                    ten.fmt(entry["percent_change"]), entry["reason"]]) in text


def test_enrichment_table_is_rendered_from_the_published_enrichment(published):
    text = report_text()
    primary = published["ten"]["enrichment"]["private"]["primary"]
    tier2 = published["ten"]["enrichment"]["private"]["tier2"]
    for observed, sensitivity in zip(primary, tier2):
        assert observed["family"] == sensitivity["family"]
        assert row([observed["family"], ten.fmt(observed["nested"]), ten.fmt(observed["expected"]),
                    ten.fmt(observed["enrichment"]), ten.fmt(observed["p"]),
                    ten.fmt(sensitivity["enrichment"])]) in text


def test_embedded_provenance_block_is_the_published_provenance(published):
    """The report quotes its own provenance; it must not drift from the record."""
    text = report_text()
    blocks = text.split("```json")
    assert len(blocks) > 1, "report carries no provenance block"
    quoted = json.loads(blocks[-1].split("```")[0])
    assert quoted == published["provenance"]


def test_report_states_every_standing_caveat(published):
    """A caveat dropped from the prose is the one most likely to be missed.

    The report is generated from the same data, so a rule can be enforced in the
    pipeline and still go unmentioned in the text a reader relies on.
    """
    text = report_text()
    required = [
        "GT/GQ are never parsed",
        "Gene/snpEff fields are not analysis inputs or matching keys",
        "chrY counts appear in QC only",
        "chr22_mei.vcf is an excluded HG03086 slice",
        "unknown never becomes antisense",
        "no transitive chaining, one carrier per sample",
        "HG03086 raw 261 is reproduced",
        ten.CAVEAT,
        opportunity.TIER2,
        "NOT carriage denominators",
    ]
    missing = [phrase for phrase in required if phrase not in text]
    assert not missing, f"standing caveats missing from the report: {missing}"


def test_report_carries_the_opportunity_caveat_on_every_enrichment_claim(published):
    """Each of the headline answers must be conditioned, not a bare number."""
    text = report_text()
    answers = text.split("## William’s three answers")[1].split("## Standing caveats")[0]
    for number in (1, 2, 3):
        line = next(line for line in answers.splitlines()
                    if line.strip().startswith(f"{number}."))
        if "enrichment" in line or "x versus" in line:
            assert ten.CAVEAT in line, f"answer {number} reports a number unconditioned"


# ---------------------------------------------------------------------------
# the derived artifacts
#
# Five further CSVs carry the bulk of the scientific output. Each is checked
# against the metrics block or against the identity that defines it, so a
# number cannot quietly stop meaning what its column name says.
# ---------------------------------------------------------------------------
def read_csv(path):
    _require(path)
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_accumulation_curve_lands_on_the_published_union_counts(published):
    """The curve is built by re-running dedup on growing cohort prefixes.

    Its fifth and tenth rows must therefore reproduce the five-genome and
    ten-genome unions the run gated on, which makes it an independent witness
    to both.
    """
    rows = read_csv(ACCUMULATION_CSV)
    assert len(rows) == len(ten.SAMPLES)
    assert [row["sample_added"] for row in rows] == list(ten.SAMPLES)
    assert [int(row["samples_included"]) for row in rows] == list(range(1, 11))

    counts = [int(row["unique_sites"]) for row in rows]
    assert counts == sorted(counts), "adding a genome must never remove a site"
    five = sum(entry["union"] for entry in published["five"]["layers"])
    ten_union = sum(entry["union"] for entry in published["ten"]["layers"])
    assert counts[4] == five == 4998
    assert counts[9] == ten_union == 7132
    for row in rows:
        assert int(row["ALU"]) + int(row["LINE1"]) + int(row["SVA"]) == int(row["unique_sites"])


def test_per_sample_rows_are_each_private_to_their_own_genome():
    """A per-genome table collapses within-sample records, so every row is n=1."""
    rows = read_csv(PER_SAMPLE_CSV)
    assert rows
    for row in rows:
        index = ten.SAMPLES.index(row["sample"])
        assert row["n_carriers"] == "1"
        assert row["presence_bitmap"].count("1") == 1
        assert row["presence_bitmap"][index] == "1"
        assert row["private_to_sample"] == row["sample"]


def test_l1_projection_status_census_matches_the_metrics_block(published):
    """The CSV and the diagnostics block are written at different points."""
    rows = read_csv(PROJECTION_CSV)
    census = Counter(row["status"] for row in rows)
    assert dict(census) == published["diagnostics"]["opportunity"]["projection_status"]
    assert census["eligible"] > 0, "no L1 reached the consensus-projected eligible state"


def profile_groups():
    groups = defaultdict(list)
    for row in read_csv(PROFILES_CSV):
        groups[(row["cohort_layer"], row["family"], row["method"])].append(row)
    return groups


def test_position_profile_expected_follows_the_opportunity_null():
    """`expected` is the group's observed count rescaled by that bin's share of
    reference opportunity -- recomputed here rather than taken on trust."""
    for key, rows in profile_groups().items():
        observed_total = sum(float(row["observed"]) for row in rows)
        opportunity_total = sum(float(row["opportunity_bp"]) for row in rows)
        for row in rows:
            share = float(row["opportunity_bp"]) / opportunity_total if opportunity_total else 0.0
            assert float(row["expected"]) == pytest.approx(
                observed_total * share, rel=1e-6, abs=1e-6), key


def test_position_profile_enrichment_is_observed_over_expected():
    for row in read_csv(PROFILES_CSV):
        expected = float(row["expected"])
        assert float(row["bin_start"]) < float(row["bin_end"])
        if expected > 0:
            assert float(row["enrichment"]) == pytest.approx(
                float(row["observed"]) / expected, rel=1e-6)
        else:
            # A zero-opportunity bin is left at zero, never pseudocounted, so
            # it cannot carry an enrichment.
            assert float(row["observed"]) == 0
            assert float(row["enrichment"]) == 0


def test_no_profile_bin_reports_events_where_opportunity_is_zero():
    """Events can only be counted against callable reference."""
    for row in read_csv(PROFILES_CSV):
        if float(row["opportunity_bp"]) == 0:
            assert float(row["observed"]) == 0, row
            assert row["zero_opportunity_excluded"] in {"True", "False"}


def test_host_parent_evidence_carriage_counts_are_self_consistent():
    rows = read_csv(HOST_EVIDENCE_CSV)
    assert rows
    for row in rows:
        samples = row["samples"].split(",") if row["samples"] else []
        assert int(row["n_genomes"]) == len(samples)
        assert len(set(samples)) == len(samples)
        assert all(sample in ten.SAMPLES for sample in samples)
        assert samples == [s for s in ten.SAMPLES if s in set(samples)], "not in cohort order"
        assert row["window"] in {"linker", "tail", "any L1-host position"}


def test_every_derived_artifact_carries_the_opportunity_caveat():
    """Each is quotable on its own, so each must carry the qualification."""
    for path in (ACCUMULATION_CSV, PER_SAMPLE_CSV, PROFILES_CSV, PROJECTION_CSV,
                 HOST_EVIDENCE_CSV, SITES_CSV, PRIVATE_CSV):
        rows = read_csv(path)
        assert rows, path.name
        assert {row["opportunity_caveat"] for row in rows} == {ten.CAVEAT}, path.name


def test_the_shipped_report_is_exactly_what_the_current_code_renders(tmp_path, published):
    """Re-render the report from the published artifacts and demand equality.

    Every other check here compares selected fields. This one rebuilds the whole
    document from the machine-readable record and compares it byte for byte, so
    no wording, ordering or rounding in the published summary can drift away
    from the numbers behind it unnoticed.

    `write_report` takes the dependent Phase 4 result as an argument; it is
    reconstructed from `phase4_status.json`, minus the two bookkeeping blocks
    that gate records rather than report content.
    """
    status = json.loads(_phase4_status_path().read_text())
    dependent = {k: v for k, v in status.items() if k not in ("provenance", "load")}

    carriage = [
        {"family": row["family"], "window": row["window"], "host_id": row["host_id"],
         "samples": row["samples"], "n_genomes": int(row["n_genomes"]),
         "opportunity_caveat": row["opportunity_caveat"]}
        for row in read_csv(HOST_EVIDENCE_CSV)
    ]

    ten.write_report(tmp_path, published["five"], published["ten"], published["qc"],
                     published["delta"], published["provenance"], published["diagnostics"],
                     carriage, dependent)

    rendered = (tmp_path / "analysis_summary.md").read_text(encoding="utf-8")
    shipped = SUMMARY_MD.read_text(encoding="utf-8")
    assert rendered == shipped, "the shipped report is not what this code renders"
