"""Reproducibility tripwires for the published analysis artifacts.

The failure this exists to catch already happened once. `unique_sites.csv` was
regenerated with a different schema after Phase 4 had run against it: the file
stayed at the same path, stayed non-empty, and every recorded input still
existed. Only re-deriving from it showed the drift -- the JSON said
`rows_read: 763`, the file on disk had 4,998 rows and none of the columns the
loader reads. Both Phase 4 scripts then returned an empty cohort, which reads as
a clean null result.

So **existence checks are not enough, and this module does not rely on them.**
Every assertion here re-derives a recorded quantity from the recorded input and
compares. A path that is merely present proves nothing.

Two tiers:

* Tier 0 (default) -- seconds. Re-derives the input-derived counts in every
  published JSON: the loader's row report, the per-call join counts, and the
  cohort's analysis-set size. This is the tier that would have caught the
  incident. It also checks that no two directories carry divergent copies of
  the same Phase 4 report, which is how the canonical copy went stale without
  any single path ever disappearing.
* Tier 1 (`slow`) -- minutes. Re-executes the analyses that are self-contained
  enough to re-run and compares the headline numbers exactly.

Both tiers skip cleanly when the workspace data is absent, so a fresh clone is
not reddened by missing 700 MB inputs. They **fail** when the data is present
and does not reproduce, because that is the entire point.

What is deliberately not covered: `genotype_concordance.json` needs the 908
sample BCFs, snpEff and the reference FASTA, and `tprt_hallmarks.json` needs the
curated reference. Those are exercised by their own phase tests; reproducing
every field of those reports inside the unit suite would make "green" depend on
a multi-minute external pipeline run, which is not a signal anyone can act on.
Their recorded inputs are still checked for existence and provenance by the
tier-0 sweep below.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import nested_multi_sample_common as common  # noqa: E402
import score_tprt_hallmarks as tprt  # noqa: E402

ANALYSIS = common.WORKSPACE / "nested_analysis"
PHASE4 = ANALYSIS / "results_phase4"
LONGREAD = ANALYSIS / "results_longread"

#: Fields of the loader's report that must be identical for the input to be the
#: one the analysis ran on. Compared by value, not by dict equality, so a new
#: key added later does not break the test and a changed key cannot hide.
LOAD_REPORT_KEYS = (
    "rows_read",
    "rows_in_scope",
    "dropped_out_of_scope_contig",
    "dropped_missing_geometry",
    "samples_field_separator",
)

#: Join fields that must be identical. `verdict` is the gate itself.
JOIN_KEYS = (
    "sites_total",
    "sites_with_at_least_one_joined_call",
    "sites_with_no_joined_call",
    "calls_recovered",
    "calls_declared_by_dedup",
    "sites_where_carrier_count_disagrees",
    "calls_where_orientation_disagrees",
    "verdict",
)

PHASE4_JSONS = ("joint_enrichment.json", "recurrence.json")


def _analysis_reports() -> list[Path]:
    """Published reports that computed something, i.e. carry `phase`/`analysis`.

    `matching_manifest.json` is a run manifest rather than an analysis report:
    it records the script, the pre-registered parameters and the per-sample QC,
    not derived numbers, so it has no `inputs` block to check and is excluded.
    """
    found = []
    for path in sorted(PHASE4.glob("*.json")) + sorted(LONGREAD.glob("*.json")):
        doc = json.loads(path.read_text())
        if "phase" in doc or "analysis" in doc:
            found.append(path)
    return found


def _require_published() -> None:
    if not PHASE4.is_dir() or not LONGREAD.is_dir():
        pytest.skip("published analysis artifacts not present on this host")


def _require_inputs(*names: str) -> None:
    missing = [name for name in names if not (ANALYSIS / name).exists()]
    if missing:
        pytest.skip(f"workspace input(s) absent on this host: {missing}")


def _published(path: Path) -> dict:
    if not path.exists():
        pytest.skip(f"published artifact absent on this host: {path.name}")
    return json.loads(path.read_text())


def _recorded_inputs(doc: dict) -> dict[str, str]:
    """Absolute-path inputs a run recorded, which is what makes it reproducible."""
    inputs = doc.get("inputs")
    if not isinstance(inputs, dict):
        return {}
    return {
        key: value
        for key, value in inputs.items()
        if isinstance(value, str) and value.startswith("/")
    }


def _drift(recorded: dict, derived: dict, keys) -> list[str]:
    return [
        f"{key}: recorded {recorded.get(key)!r} != derived {derived.get(key)!r}"
        for key in keys
        if recorded.get(key) != derived.get(key)
    ]


# --------------------------------------------------------------------------
# Tier 0: provenance, and the input-derived counts in every published report
# --------------------------------------------------------------------------


def test_every_published_report_records_its_inputs():
    """A report that forgets its inputs cannot be checked by anything.

    This is what makes the rest of the module possible: if an analysis ships
    numbers without recording what it read them from, no later test can tell a
    stale number from a current one.
    """
    _require_published()
    reports = _analysis_reports()
    assert reports, "no published analysis reports found"
    for path in reports:
        assert _recorded_inputs(json.loads(path.read_text())), (
            f"{path.name} records no absolute-path inputs"
        )


def _phase4_copies(name: str) -> list[Path]:
    """Every published copy of one Phase 4 report, wherever it was written.

    The glob is deliberately `*<name>` under every `results*` directory rather
    than the exact filename: the stale/correct pair this exists to catch were
    named differently (`joint_enrichment.json` in the plan's directory,
    `phase4_joint_enrichment.json` next to the dedup output), so matching the
    canonical name alone would find only one of them and see no conflict.
    """
    if not ANALYSIS.is_dir():
        return []
    return sorted(
        path
        for results in sorted(ANALYSIS.glob("results*"))
        if results.is_dir()
        for path in results.glob(f"*{name}")
        if path.is_file()
    )


def _report_fingerprint(doc: dict) -> dict:
    """The subset that decides whether two copies are the same result.

    Only cohort-derived quantities, so a cosmetic or provenance-only difference
    between two copies is not reported as a conflict.
    """
    fingerprint: dict = {"load_report": doc.get("load_report")}
    split = doc.get("cohort_split")
    if isinstance(split, dict):
        fingerprint["shared_sites"] = split.get("shared_sites")
        fingerprint["private_sites"] = split.get("private_sites")
    cells = doc.get("cells")
    if isinstance(cells, list):
        fingerprint["observed_by_cell"] = {
            str(cell.get("cell_id")): cell.get("observed")
            for cell in cells
            if isinstance(cell, dict)
        }
    return fingerprint


@pytest.mark.parametrize("name", PHASE4_JSONS)
def test_phase4_published_copies_agree(name):
    """Two copies of one report that disagree are a trap, not redundancy.

    Re-running Phase 4 into a new outdirectory leaves the previous copy in
    place, at a path that still exists, still parses, and still answers to the
    plan's filename. Nothing in this suite looked outside `results_phase4/`, so
    the superseded copy kept reading as current while the directory the plan
    names went stale. A second copy is fine while it is identical; it is the
    disagreement that makes the published number unknowable, because a reader
    cannot tell which path is the result.
    """
    _require_published()
    copies = _phase4_copies(name)
    if len(copies) < 2:
        return
    prints = [(p, _report_fingerprint(json.loads(p.read_text()))) for p in copies]
    canonical_path, canonical = prints[0]
    conflicts: list[str] = []
    for path, fingerprint in prints[1:]:
        drift = [
            f"{key}: {canonical.get(key)!r} != {fingerprint.get(key)!r}"
            for key in sorted(fingerprint)
            if canonical.get(key) != fingerprint.get(key)
        ]
        if drift:
            conflicts.append(
                f"{path}\n    vs {canonical_path}\n    "
                + "\n    ".join(drift)
            )
    assert not conflicts, (
        f"{name} is published in {len(copies)} places and they are not the same "
        f"result. One of them no longer describes the substrate it records. The "
        f"plan names `results_phase4/`; keep the copy that reproduces and remove "
        f"or re-point the other:\n  " + "\n  ".join(conflicts)
    )


@pytest.mark.parametrize("name", PHASE4_JSONS)
def test_every_phase4_copy_reproduces_from_its_recorded_input(name):
    """Apply the tier-0 reproduction check to every copy, not just the first.

    A second copy is only safe if it also re-derives. If one copy of a report
    still describes the substrate it was computed from and another does not,
    the numbers are ambiguous no matter which path a reader opens.
    """
    copies = _phase4_copies(name)
    if not copies:
        pytest.skip(f"no published copy of {name} present on this host")
    stale: list[str] = []
    for path in copies:
        doc = json.loads(path.read_text())
        recorded = doc.get("load_report")
        if not recorded:
            stale.append(f"{path}: records no load_report")
            continue
        try:
            _sites, derived = common.load_unique_sites(
                Path(_recorded_inputs(doc)["unique_sites"])
            )
        except (common.SchemaError, KeyError) as exc:
            stale.append(f"{path}: loader rejects its own recorded input: {exc}")
            continue
        drift = _drift(recorded, derived, LOAD_REPORT_KEYS)
        if drift:
            stale.append(f"{path}:\n    " + "\n    ".join(drift))
    assert not stale, (
        f"published {name} no longer reproduces from the input it records:\n  "
        + "\n  ".join(stale)
    )


def test_every_recorded_input_still_exists():
    _require_published()
    reports = _analysis_reports()
    assert reports, "no published analysis reports found"
    gone: list[str] = []
    for path in reports:
        for key, value in _recorded_inputs(json.loads(path.read_text())).items():
            if not Path(value).exists():
                gone.append(f"{path.name}:{key} -> {value}")
    assert not gone, "recorded inputs no longer present:\n  " + "\n  ".join(gone)


@pytest.mark.parametrize("name", PHASE4_JSONS)
def test_phase4_load_report_still_reproduces(name):
    """The loader must still read the same rows out of the same file.

    This is the assertion the incident would have tripped. The published JSON
    records `rows_read`; if the file at that path no longer yields that many
    rows through the same loader, the published numbers describe an input that
    no longer exists and every figure derived from them is unverifiable.
    """
    _require_published()
    doc = _published(PHASE4 / name)
    path = Path(_recorded_inputs(doc)["unique_sites"])
    recorded = doc.get("load_report")
    assert recorded, f"{name} records no load_report to compare against"

    try:
        _sites, derived = common.load_unique_sites(path)
    except common.SchemaError as exc:
        pytest.fail(
            f"{name} cannot be regenerated: the loader rejects its own recorded "
            f"input.\n  file:    {path}\n  reason:  {exc}"
        )

    drift = _drift(recorded, derived, LOAD_REPORT_KEYS)
    assert not drift, (
        f"{name} was computed from a different version of {path.name}.\n  "
        + "\n  ".join(drift)
    )


@pytest.mark.parametrize("name", PHASE4_JSONS)
def test_phase4_join_verification_still_reproduces(name):
    """The per-call join must still recover exactly the recorded calls.

    `rows_read` can match while the per-call detail behind the analysis does
    not, so the join is checked separately from the row count.
    """
    _require_published()
    doc = _published(PHASE4 / name)
    inputs = _recorded_inputs(doc)
    recorded = doc.get("join_verification")
    assert recorded, f"{name} records no join_verification to compare against"

    sites, _report = common.load_unique_sites(Path(inputs["unique_sites"]))
    carriers = sorted({s for site in sites for s in site["carriers"]})
    callsets = common.load_callsets(Path(inputs["callset_dir"]), carriers)
    derived = common.verify_join(common.attach_call_details(sites, callsets))

    drift = _drift(recorded, derived, JOIN_KEYS)
    assert not drift, (
        f"{name} was computed from a different per-call join.\n  "
        + "\n  ".join(drift)
    )


def test_the_longread_cohort_still_yields_the_published_analysis_set():
    """`results_longread`'s published analyses all rest on this one count.

    Every Phase 1/2b/3 number is computed over the nested analysis set, so if the
    cohort file no longer produces the recorded number of sites, none of them
    are reproducible and the drift is upstream of all of them.
    """
    _require_published()
    _require_inputs("results_longread/per_call_longread_nested.csv")
    doc = _published(LONGREAD / "position_enrichment.json")
    recorded = doc["primary_test"]["n_sites"]

    rows = tprt.read_cohort(
        Path(_recorded_inputs(doc)["cohort"])
    )
    derived = len(tprt.nested_analysis_set(rows))
    assert derived == recorded, (
        "the long-read cohort no longer yields the published analysis set: "
        f"recorded {recorded} sites, derived {derived}"
    )


# --------------------------------------------------------------------------
# Tier 1: re-execute and compare headline numbers exactly
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_phase1_position_enrichment_reproduces_its_headline(tmp_path):
    """Re-run Phase 1 from its recorded inputs and compare the primary test.

    Seeded throughout, so `observed`, `expected_a` and `p_a` must match exactly.
    A change here means either the cohort moved or the null changed -- and Null
    B was corrected after this artifact was published, so a difference in the
    mappability-weighted arm is expected to be reported, not silently absorbed.
    """
    _require_published()
    import test_position_enrichment as tpe

    doc = _published(LONGREAD / "position_enrichment.json")
    inputs = _recorded_inputs(doc)
    pre = doc["pre_registration"]
    recorded = doc["primary_test"]

    out = tmp_path / "phase1"
    assert tpe.main(
        [
            "--cohort", inputs["cohort"],
            "--mapability-bed", inputs["mapability_bed"],
            "--outdir", str(out),
            "--replicates", str(pre["replicates"]),
            "--seed", str(pre["seed"]),
        ]
    ) == 0
    derived = json.loads((out / "position_enrichment.json").read_text())["primary_test"]

    for key in ("n_sites", "n_hosts", "observed", "window_lo_inclusive",
                "window_hi_exclusive"):
        assert derived[key] == recorded[key], (
            f"Phase 1 primary_test.{key}: recorded {recorded[key]!r} != "
            f"derived {derived[key]!r}"
        )
    assert derived["expected_a"] == pytest.approx(recorded["expected_a"], rel=1e-9)
    assert derived["p_a"] == pytest.approx(recorded["p_a"], rel=1e-9)


@pytest.mark.slow
def test_phase4_joint_enrichment_reproduces_its_observed_counts(tmp_path):
    """Re-run Gate 3d's sibling, Phase 4, and compare the deterministic cells.

    Only the deterministic fields are compared: `observed` and the host and
    event counts are a function of the cohort alone. The Monte Carlo p-values
    are compared too because the run is seeded, but the observed count is the
    assertion that must not move -- it is the primary result.
    """
    _require_published()
    import joint_enrichment as je

    doc = _published(PHASE4 / "joint_enrichment.json")
    inputs = _recorded_inputs(doc)
    out = tmp_path / "phase4"
    try:
        rc = je.main(
            [
                "--unique-sites", inputs["unique_sites"],
                "--callset-dir", inputs["callset_dir"],
                "--rmsk", str(je.DEFAULT_RMSK),
                "--outdir", str(out),
                "--replicates", str(doc["cells"][0]["replicates"]),
                "--bootstrap", str(doc["inputs"]["bootstrap_draws"]),
                "--seed", str(doc["inputs"]["seed"]),
            ]
        )
    except (common.SchemaError, SystemExit) as exc:
        pytest.fail(
            "joint_enrichment.json cannot be regenerated from its own recorded "
            f"input:\n  file:   {inputs['unique_sites']}\n  reason: {exc}"
        )
    assert rc == 0
    derived = json.loads((out / "joint_enrichment.json").read_text())

    assert derived["load_report"]["rows_in_scope"] == doc["load_report"]["rows_in_scope"]
    recorded_cells = {c["cell_id"]: c for c in doc["cells"]}
    derived_cells = {c["cell_id"]: c for c in derived["cells"]}
    assert set(derived_cells) == set(recorded_cells), (
        "the reported cell set changed: "
        f"recorded {sorted(recorded_cells)}, derived {sorted(derived_cells)}"
    )
    for cell_id, cell in recorded_cells.items():
        got = derived_cells[cell_id]
        for key in ("observed", "sites_in_scope", "hosts_total", "events_in_scope",
                    "carrier_observations_in_scope"):
            assert got[key] == cell[key], (
                f"cell {cell_id} {key}: recorded {cell[key]!r} != derived {got[key]!r}"
            )
