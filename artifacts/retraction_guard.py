"""Single-source retraction checks shared by report generation and pytest."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping

REGISTRY_PATH = Path(__file__).with_name("retraction_replacements.json")


def load_registry(path: Path = REGISTRY_PATH) -> list[dict[str, object]]:
    """Read and minimally validate the machine-readable correction registry."""
    registry = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(registry, list) or not registry:
        raise ValueError(f"{path}: expected a non-empty JSON list")
    for entry in registry:
        if not isinstance(entry.get("retracted"), list) or not entry["retracted"]:
            raise ValueError(f"{path}: each entry needs non-empty retracted strings")
        if not isinstance(entry.get("id"), str) or not entry["id"]:
            raise ValueError(f"{path}: each entry needs an id")
        if not isinstance(entry.get("replacement"), str) or not entry["replacement"]:
            raise ValueError(f"{path}: each entry needs a replacement string")
        if not isinstance(entry.get("required_files"), list):
            raise ValueError(f"{path}: required_files must be a list")
        if not entry["required_files"] and not entry.get("derived_assertion"):
            raise ValueError(f"{path}: each entry needs a required file or derived assertion")
    return registry


def find_retracted_strings(text: str, path: Path, registry: list[dict[str, object]]) -> list[str]:
    """Return one actionable diagnostic per stale phrase occurrence."""
    errors = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        for entry in registry:
            for stale in entry["retracted"]:
                if stale in line:
                    errors.append(f"{path}:{line_number}: retracted claim {stale!r}")
    return errors


def validate_sources(
    sources: Mapping[Path, str], registry_path: Path = REGISTRY_PATH
) -> None:
    """Fail on stale text or when an approved replacement is absent from its target."""
    registry = load_registry(registry_path)
    errors = [error for path, text in sources.items() for error in find_retracted_strings(text, path, registry)]
    for entry in registry:
        replacement = str(entry["replacement"])
        assertion = entry.get("derived_assertion")
        if assertion:
            target = str(assertion["source"])
            source_path = next((path for path in sources if path.as_posix().endswith(target)), None)
            text = sources.get(source_path) if source_path is not None else None
            if text is None:
                errors.append(f"{target}: derived correction source was not scanned")
            else:
                rows = list(csv.DictReader(text.splitlines()))
                filters = assertion["filters"]
                observed = sum(
                    row.get(assertion["column"]) == assertion["value"]
                    and all(row.get(key) == value for key, value in filters.items())
                    for row in rows
                )
                denominator = sum(
                    all(row.get(key) == value for key, value in assertion["denominator_filters"].items())
                    for row in rows
                )
                if observed != int(assertion["expected"]):
                    errors.append(f"{target}: corrected catalog count expected {assertion['expected']}, found {observed}")
                if denominator != int(assertion["denominator"]):
                    errors.append(
                        f"{target}: corrected denominator expected {assertion['denominator']}, found {denominator}"
                    )
                replacement = replacement.format(observed=observed, denominator=denominator)
                description = str(entry.get("description", ""))
                if replacement not in description:
                    errors.append(f"{entry['id']}: corrected claim missing from registry description: {replacement!r}")
        for required_file in entry["required_files"]:
            target = str(required_file)
            source_path = next((path for path in sources if path.as_posix().endswith(target)), None)
            text = sources.get(source_path) if source_path is not None else None
            if text is None:
                errors.append(f"{target}: required correction target was not scanned")
            elif replacement not in text:
                errors.append(f"{target}: corrected claim missing: {replacement!r}")
    if errors:
        raise ValueError("\n".join(errors))
