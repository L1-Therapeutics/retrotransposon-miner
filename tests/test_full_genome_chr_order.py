"""--chr all starts sex chromosomes so chrX fills a first-wave slot."""

from __future__ import annotations

import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "scripts" / "run_candidate_discovery_and_annotation.sh"


def _resolve(chr_arg: str) -> list[str]:
    text = WRAPPER.read_text(encoding="utf-8")
    start = text.index("normalize_chr_token() {")
    end = text.index("\nset_reference_build_defaults()", start)
    functions = text[start:end]
    script = f"""
set -euo pipefail
REGION=chr22
{functions}
resolve_chr_list {chr_arg!r}
"""
    out = subprocess.run(
        ["bash", "-c", script],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in out.stdout.splitlines() if line.strip()]


def test_chr_all_starts_with_sex_chromosomes_then_autosomes():
    got = _resolve("all")
    assert got[:4] == ["chrX", "chrY", "chr1", "chr2"]
    assert got[15] == "chr14"
    assert got[-1] == "chr22"
    assert len(got) == 24
    assert got.index("chrX") < 16


def test_explicit_chr_list_keeps_caller_order():
    assert _resolve("chr22,chr1,X") == ["chr22", "chr1", "chrX"]
