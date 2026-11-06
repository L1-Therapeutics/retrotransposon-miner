"""--chr all starts sex chromosomes so chrX fills a first-wave slot."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "scripts" / "run_candidate_discovery_and_annotation.sh"


def _resolve(chr_arg: str, reference_build: str = "") -> list[str]:
    text = WRAPPER.read_text(encoding="utf-8")
    start = text.index("normalize_chr_token() {")
    end = text.index("\nset_reference_build_defaults()", start)
    functions = text[start:end]
    script = f"""
set -euo pipefail
REGION=chr22
REFERENCE_BUILD={reference_build}
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


def test_mm39_chr_all_is_mouse_chromosomes():
    got = _resolve("all", "mm39")
    assert got[:4] == ["chrX", "chrY", "chr1", "chr2"]
    assert got[-1] == "chr19"
    assert "chr22" not in got
    assert len(got) == 21


def test_mm39_rejects_human_chr22():
    with pytest.raises(subprocess.CalledProcessError):
        _resolve("22", "mm39")


def test_mm39_defaults_use_mouse_annotation_files():
    text = WRAPPER.read_text(encoding="utf-8")
    start = text.index("\n    mm39)")
    end = text.index("\n    *)", start)
    block = text[start:end]
    for path in (
        "reference/mm39/mm39.fa",
        "annotation/mm39/repeats/rmsk.txt.gz",
        "annotation/mm39/masks/gap.txt.gz",
        "annotation/mm39/blacklist/blacklist.mm39.bed",
        "annotation/mm39/segdup/genomicSuperDups.bed",
        "annotation/mm39/mappability/k100.Umap.MultiTrackMappability.bedGraph",
        "annotation/mm39/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed",
        "annotation/mm39/junk/junk_exclusion_merged.bed",
        "annotation/mm39/gencode/gencode.vM39.annotation.gtf.gz",
        "retrotransposon_db/dfam/dfam_mouse_mei_l1_sine_erv.fasta",
        "polymorphism/mm39/mgp/mgpv5.SV_insertions.retroseq.bed",
    ):
        assert path in block


def test_explicit_chr_list_keeps_caller_order():
    assert _resolve("chr22,chr1,X") == ["chr22", "chr1", "chrX"]
