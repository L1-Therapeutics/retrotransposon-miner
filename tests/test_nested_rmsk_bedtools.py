"""Nested rmsk annotation via bedtools."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from retro_miner.mei_support import (
    _annotate_nested_retrotransposon,
    _bed_field,
    _write_bed_row,
)


def test_bed_field_never_emits_empty_or_tabs():
    assert _bed_field("") == "."
    assert _bed_field(None) == "."
    assert _bed_field(float("nan")) == "."
    assert _bed_field("ALU") == "ALU"
    assert _bed_field("+\tbad") == "+ bad"
    assert _bed_field("", default="x") == "x"


def test_write_bed_row_avoids_trailing_empty_fields(tmp_path: Path):
    path = tmp_path / "x.bed"
    with path.open("w", encoding="utf-8") as hout:
        _write_bed_row(hout, ["chr22", 1, 2, 0, 0, ".", "ALU", ""])
        _write_bed_row(hout, ["chr22", 3, 4, 1, 0, ".", "", "+"])
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "chr22\t1\t2\t0\t0\t.\tALU\t."
    assert lines[1] == "chr22\t3\t4\t1\t0\t.\t.\t+"
    assert all(not ln.endswith("\t") for ln in lines)


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_nested_rmsk_accepts_missing_orientation(tmp_path: Path):
    """Empty orientation used to emit trailing tabs and crash bedtools."""
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n",
        encoding="utf-8",
    )
    cand = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 90,
                "window_end": 110,
                "insertion_breakpoint_pos": 150,
                # no orientation / consensus fields on purpose
                "disease_L_mei_subfamily": "AluY",
                "disease_L_mei_supported_reads": 5,
            },
            {
                "chrom": "chr22",
                "window_start": 5000,
                "window_end": 5200,
                "insertion_breakpoint_pos": 5100,
                "disease_L_mei_subfamily": "AluY",
                "disease_L_mei_supported_reads": 3,
            },
        ]
    )
    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert bool(out.loc[0, "nested_repeat_overlap"]) is True
    assert out.loc[0, "nested_mei_family"] == "ALU"
    assert bool(out.loc[0, "nested_same_orientation"]) is False
    assert bool(out.loc[1, "nested_repeat_overlap"]) is False


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_nested_rmsk_bedtools_flags_same_family_hits(tmp_path: Path):
    # BED-like parser expects: chrom start end . . strand repName repClass repFamily
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n"
        "chr22\t1000\t1500\t.\t0\t-\tL1HS\tLINE\tL1\n",
        encoding="utf-8",
    )

    cand = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 90,
                "window_end": 110,
                "insertion_breakpoint_pos": 150,  # inside Alu
                "consensus_insertion_orientation": "+",
                "disease_L_mei_subfamily": "AluY",
                "disease_L_mei_supported_reads": 5,
            },
            {
                "chrom": "chr22",
                "window_start": 900,
                "window_end": 1100,
                "insertion_breakpoint_pos": 1200,  # inside L1
                "consensus_insertion_orientation": "-",
                "disease_L_mei_subfamily": "L1HS",
                "disease_L_mei_supported_reads": 4,
            },
            {
                "chrom": "chr22",
                "window_start": 5000,
                "window_end": 5200,
                "insertion_breakpoint_pos": 5100,  # unnested
                "consensus_insertion_orientation": "+",
                "disease_L_mei_subfamily": "AluY",
                "disease_L_mei_supported_reads": 3,
            },
        ]
    )

    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert bool(out.loc[0, "nested_repeat_overlap"]) is True
    assert out.loc[0, "nested_repeat_name"] == "AluY"
    assert out.loc[0, "nested_mei_family"] == "ALU"
    assert bool(out.loc[1, "nested_repeat_overlap"]) is True
    assert out.loc[1, "nested_mei_family"] == "LINE1"
    assert bool(out.loc[2, "nested_repeat_overlap"]) is False
    assert out.loc[2, "nested_same_class_orientation"] == "unnested"


# ------------------------------------------------------- three-valued nesting vocabulary
# "unnested" means no same-family element overlaps the breakpoint. An insertion
# inside a same-family element in the opposite orientation is still nested, so it
# gets its own value rather than borrowing "unnested" from the no-overlap case.


def _candidate(chrom: str, pos: int, subfamily: str, orientation: str) -> dict:
    return {
        "chrom": chrom,
        "window_start": pos - 60,
        "window_end": pos + 60,
        "insertion_breakpoint_pos": pos,
        "consensus_insertion_orientation": orientation,
        "disease_L_mei_subfamily": subfamily,
        "disease_L_mei_supported_reads": 5,
    }


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_opposite_orientation_nesting_is_not_reported_as_unnested(tmp_path: Path):
    """An opposite-orientation call is nested in the element; only the strand differs."""
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text("chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n", encoding="utf-8")
    cand = pd.DataFrame([_candidate("chr22", 150, "AluY", "-")])  # insertion - inside a + AluY

    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert bool(out.loc[0, "nested_repeat_overlap"]) is True
    assert bool(out.loc[0, "nested_same_class"]) is True
    assert bool(out.loc[0, "nested_same_orientation"]) is False
    assert out.loc[0, "nested_same_class_orientation"] == "nested_opposite_orientation"
    # The exact contradiction this fixes: an overlapping call was labelled unnested.
    assert out.loc[0, "nested_same_class_orientation"] != "unnested"


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_same_orientation_nesting_is_labelled_nested_same_orientation(tmp_path: Path):
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text("chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n", encoding="utf-8")
    cand = pd.DataFrame([_candidate("chr22", 150, "AluY", "+")])  # insertion + inside a + AluY

    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert bool(out.loc[0, "nested_repeat_overlap"]) is True
    assert bool(out.loc[0, "nested_same_orientation"]) is True
    assert out.loc[0, "nested_same_class_orientation"] == "nested_same_orientation"


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_no_overlapping_element_stays_unnested(tmp_path: Path):
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text("chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n", encoding="utf-8")
    cand = pd.DataFrame([_candidate("chr22", 9000, "AluY", "-")])  # far from the element

    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert bool(out.loc[0, "nested_repeat_overlap"]) is False
    assert bool(out.loc[0, "nested_same_class"]) is False
    assert out.loc[0, "nested_same_class_orientation"] == "unnested"


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not on PATH")
def test_opposite_and_same_orientation_nesting_are_distinguished(tmp_path: Path):
    """The two overlapping cases must not collapse onto the same value."""
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text(
        "chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n"
        "chr22\t1000\t1500\t.\t0\t-\tAluYb\tSINE\tAlu\n",
        encoding="utf-8",
    )
    cand = pd.DataFrame(
        [
            _candidate("chr22", 150, "AluY", "-"),    # opposite to the + element
            _candidate("chr22", 1200, "AluYb", "-"),  # same as the - element
        ]
    )

    out = _annotate_nested_retrotransposon(cand, rmsk)
    assert out.loc[0, "nested_same_class_orientation"] == "nested_opposite_orientation"
    assert out.loc[1, "nested_same_class_orientation"] == "nested_same_orientation"
    assert out.loc[0, "nested_same_class_orientation"] != out.loc[1, "nested_same_class_orientation"]


def test_nested_rmsk_requires_bedtools(monkeypatch, tmp_path: Path):
    monkeypatch.setattr("retro_miner.mei_support.shutil.which", lambda _name: None)
    rmsk = tmp_path / "rmsk.txt"
    rmsk.write_text("chr22\t100\t400\t.\t0\t+\tAluY\tSINE\tAlu\n", encoding="utf-8")
    cand = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "window_start": 90,
                "window_end": 110,
                "insertion_breakpoint_pos": 150,
                "disease_L_mei_subfamily": "AluY",
                "disease_L_mei_supported_reads": 5,
            }
        ]
    )
    with pytest.raises(RuntimeError, match="requires bedtools"):
        _annotate_nested_retrotransposon(cand, rmsk)
