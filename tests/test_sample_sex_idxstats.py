"""chrY is skipped for females using samtools idxstats on real BAM indexes.

HG03172.header.bam.bai is the .bai from that sample's whole-genome BAM.
HG00100.header.bam.bai carries the mapped counts from idxstats of the
HG00100 whole-genome BAM. That original .bai was deleted with the alignment;
the lengths are from the HG00100 CRAM header.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from retro_miner.sample_sex import classify_sex

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_candidate_discovery_and_annotation.sh"
FIXTURES = ROOT / "tests" / "fixtures" / "idxstats"

MALE_BAM = FIXTURES / "HG03172" / "HG03172.header.bam"
FEMALE_BAM = FIXTURES / "HG00100" / "HG00100.header.bam"

# Mapped counts from the real idxstats runs.
MALE_MAPPED = {"chr1": 61552506, "chrX": 18661761, "chrY": 6251340}
FEMALE_MAPPED = {"chr1": 59072342, "chrX": 34861419, "chrY": 1269068}


def _require_samtools() -> str:
    samtools = shutil.which("samtools")
    if samtools is None:
        pytest.skip("samtools is required to read the fixture .bai files")
    return samtools


def _idxstats(bam: Path) -> str:
    samtools = _require_samtools()
    result = subprocess.run(
        [samtools, "idxstats", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def _mapped(idxstats_text: str) -> dict[str, int]:
    counts = {}
    for line in idxstats_text.splitlines():
        chrom, _length, mapped, *_rest = line.split("\t")
        if chrom in {"chr1", "chrX", "chrY"}:
            counts[chrom] = int(mapped)
    return counts


def _shell_chroms(bam: Path) -> list[str]:
    _require_samtools()
    text = SCRIPT.read_text(encoding="utf-8")
    start = text.index("# BEGIN drop_chry_when_female")
    end = text.index("# END drop_chry_when_female") + len("# END drop_chry_when_female")
    function = text[start:end]
    script = f"""
set -euo pipefail
{function}
is_remote_alignment() {{ return 1; }}
run_python_module() {{ python -m "$@"; }}
CHR_LIST=(chrX chrY chr1)
DISEASE_BAM={bam}
drop_chry_when_female
printf '%s\\n' "${{CHR_LIST[@]}}"
"""
    result = subprocess.run(
        ["bash", "-c", script],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={
            "PATH": __import__("os").environ["PATH"],
            "PYTHONPATH": str(ROOT / "src"),
        },
    )
    chroms = [
        line
        for line in result.stdout.splitlines()
        if line.startswith("chr")
    ]
    return result.stdout, chroms


def test_hg03172_male_index_keeps_chry():
    stats = _idxstats(MALE_BAM)
    assert _mapped(stats) == MALE_MAPPED
    label, ratio = classify_sex(stats)
    assert label == "male"
    assert ratio is not None and ratio < 0.75
    log, chroms = _shell_chroms(MALE_BAM)
    assert "male; keeping chrY" in log
    assert chroms == ["chrX", "chrY", "chr1"]


def test_hg00100_female_index_skips_chry():
    stats = _idxstats(FEMALE_BAM)
    assert _mapped(stats) == FEMALE_MAPPED
    label, ratio = classify_sex(stats)
    assert label == "female"
    assert ratio is not None and ratio >= 0.75
    log, chroms = _shell_chroms(FEMALE_BAM)
    assert "female; skipping chrY" in log
    assert chroms == ["chrX", "chr1"]


def test_missing_chrx_keeps_the_decision_unknown():
    label, ratio = classify_sex("chr1\t248956422\t100\t0\n")
    assert label == "unknown"
    assert ratio is None
