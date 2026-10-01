"""Generated-report drift guard backed by the single JSON correction registry."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = REPO_ROOT / "artifacts"
sys.path.insert(0, str(ARTIFACTS))
from retraction_guard import find_retracted_strings, load_registry, validate_sources  # noqa: E402


def _sources() -> dict[Path, str]:
    report = ARTIFACTS / "selfins_report.md"
    literature = ARTIFACTS / "nested_orientation_literature.csv"
    panel = ARTIFACTS / "panel.py"
    catalog = ARTIFACTS / "l1_locus_catalog_chr22.csv"
    mechanism = ARTIFACTS / "alu_mechanism_literature.csv"
    return {
        report: report.read_text(encoding="utf-8"),
        literature: literature.read_text(encoding="utf-8"),
        panel: panel.read_text(encoding="utf-8"),
        catalog: catalog.read_text(encoding="utf-8"),
        mechanism: mechanism.read_text(encoding="utf-8"),
    }


def test_generated_report_and_literature_have_corrected_claims_and_no_retractions():
    sources = _sources()
    registry = load_registry()
    errors = [error for path, text in sources.items() for error in find_retracted_strings(text, path, registry)]
    assert not errors, "\n".join(errors)
    validate_sources(sources)



def test_guard_rejects_fixture_containing_retracted_phrase_with_file_and_line():
    registry = load_registry()
    fixture = Path(__file__).parent / "fixtures" / "retraction" / "stale_claim.md"
    text = fixture.read_text(encoding="utf-8")
    errors = find_retracted_strings(text, fixture, registry)
    assert errors == [f"{fixture}:2: retracted claim '−0.8 pp'"]
    with pytest.raises(ValueError, match=r"stale_claim\.md:2: retracted claim"):
        validate_sources({**_sources(), fixture: text})


def test_derived_ttc28_replacement_is_explicit_and_required(tmp_path):
    sources = _sources()
    registry_path = ARTIFACTS / "retraction_replacements.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    ttc28 = next(entry for entry in registry if entry["id"] == "ttc28-locus-count")
    corrected = ttc28["replacement"].format(observed=7, denominator=15)
    assert corrected in ttc28["description"]
    validate_sources(sources)

    ttc28["description"] = ttc28["description"].replace(corrected, "catalog count correction")
    broken_registry = tmp_path / "registry.json"
    broken_registry.write_text(json.dumps(registry), encoding="utf-8")
    with pytest.raises(ValueError, match=r"ttc28-locus-count: corrected claim missing"):
        validate_sources(sources, broken_registry)


def test_registry_is_json_and_covers_required_claim_families():
    raw = json.loads((ARTIFACTS / "retraction_replacements.json").read_text(encoding="utf-8"))
    ids = {entry["id"] for entry in raw}
    assert {
        "panel-delta",
        "ttc28-locus-count",
        "age-enrichment-default",
        "levy-category-ratio",
        "conley-jordan-authors",
        "tajaddod-authors",
    } <= ids
