"""PolyA trim + length-aware MEI alignment gate."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _apply_mei_align_quality_gate,
    _demote_polya_split_mei_hits,
    _mei_align_quality_ok,
    _split_mei_support_eligible_mask,
    _trim_poly_at_from_clip,
    trim_mei_consensus_terminal_polya,
)


def test_trim_polyA_keeps_right_tip_for_left_clip() -> None:
    # L soft-clip: polyA often abuts the junction; tip is 3' of the polyA run.
    tip = "GGCCGGGCGCGGTGGCTCAC"
    seq = ("A" * 40) + tip
    assert _trim_poly_at_from_clip(seq, "L") == tip


def test_trim_polyT_keeps_left_tip_for_right_clip() -> None:
    tip = "GGCCGGGCGCGGTGGCTCAC"
    seq = tip + ("T" * 40)
    assert _trim_poly_at_from_clip(seq, "R") == tip


def test_trim_noop_without_long_homopolymer() -> None:
    seq = "GGCCGGGCGCGGTGGCTCACGCCTGTAATC"
    assert _trim_poly_at_from_clip(seq, "L") == seq


def test_trim_pure_polya_returns_empty() -> None:
    # Must not fall back to the untrimmed A-run (that remaps as MEI_MAPPED).
    assert _trim_poly_at_from_clip("A" * 40, "L") == ""
    assert _trim_poly_at_from_clip("T" * 40, "R") == ""


def test_polya_split_not_mei_mapped_eligible() -> None:
    df = pd.DataFrame(
        [
            {
                "mei_hit": True,
                "clip_len": 40,
                "clip_poly_at_run": 25,
                "poly_tail_rescued": True,
            },
            {
                "mei_hit": True,
                "clip_len": 40,
                "clip_poly_at_run": 0,
                "poly_tail_rescued": False,
            },
        ]
    )
    ok = _split_mei_support_eligible_mask(df)
    assert bool(ok.iloc[0]) is False
    assert bool(ok.iloc[1]) is True


def test_demote_polya_clears_mei_hit() -> None:
    df = pd.DataFrame(
        [
            {
                "mei_hit": True,
                "mei_hit_coord": True,
                "target": "AluY#SINE/Alu",
                "family": "ALU",
                "target_start": 272,
                "target_end": 311,
                "clip_poly_at_run": 30,
                "poly_tail_rescued": True,
            }
        ]
    )
    out = _demote_polya_split_mei_hits(df)
    assert bool(out.loc[0, "mei_hit"]) is False
    assert bool(out.loc[0, "mei_hit_coord"]) is False
    assert str(out.loc[0, "target"]) == ""
    assert int(out.loc[0, "target_start"]) == 0


def test_trim_mei_consensus_strips_terminal_polya() -> None:
    body = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTT"
    assert trim_mei_consensus_terminal_polya(body + ("A" * 40)) == body
    assert trim_mei_consensus_terminal_polya(("T" * 20) + body) == body
    assert trim_mei_consensus_terminal_polya(body + "AAAAAAA") == body + "AAAAAAA"  # <8 kept


def test_write_polya_trimmed_fasta_preserves_single_gt(tmp_path) -> None:
    from retro_miner.mei_support import _write_polya_trimmed_fasta

    src = tmp_path / "mei.fa"
    dst = tmp_path / "mei.nopolya.fa"
    src.write_text(">AluY#SINE/Alu\n" + ("ACGT" * 10) + ("A" * 20) + "\n", encoding="utf-8")
    assert _write_polya_trimmed_fasta(src, dst) == 1
    first = dst.read_text(encoding="utf-8").splitlines()[0]
    assert first == ">AluY#SINE/Alu"
    assert not first.startswith(">>")


def test_short_query_requires_high_qcov() -> None:
    df = pd.DataFrame(
        {
            "target": ["AluY", "AluY"],
            "qcov": [0.95, 0.50],
            "pid": [0.95, 0.95],
            "alnlen": [19, 10],
        }
    )
    qlen = pd.Series([20, 20])
    ok = _mei_align_quality_ok(df, query_len=qlen)
    assert bool(ok.iloc[0]) is True
    assert bool(ok.iloc[1]) is False


def test_short_query_qcov_implies_min_aln() -> None:
    # qcov≥0.80 on len≥12 already covers ≳10 bp; no separate alnlen floor.
    df = pd.DataFrame(
        {
            "target": ["AluY"],
            "qcov": [0.80],
            "pid": [0.95],
            "alnlen": [10],
        }
    )
    assert bool(_mei_align_quality_ok(df, query_len=pd.Series([12])).iloc[0]) is True


def test_long_query_allows_tip_alignment() -> None:
    # 100 bp clip with only a 25 bp high-identity tip must still pass.
    df = pd.DataFrame(
        {
            "target": ["L1HS"],
            "qcov": [0.25],
            "pid": [0.96],
            "alnlen": [25],
        }
    )
    ok = _mei_align_quality_ok(df, query_len=pd.Series([100]))
    assert bool(ok.iloc[0]) is True


def test_gate_clears_weak_hits() -> None:
    df = pd.DataFrame(
        {
            "target": ["AluY", "AluY"],
            "family": ["ALU", "ALU"],
            "target_strand": ["+", "+"],
            "target_start": [1, 1],
            "target_end": [20, 10],
            "target_len": [300, 300],
            "alnlen": [19, 10],
            "mapq": [0, 0],
            "pid": [0.95, 0.95],
            "qcov": [0.95, 0.40],
            "mei_score": [0.9, 0.5],
        }
    )
    gated = _apply_mei_align_quality_gate(df, query_len=pd.Series([20, 20]))
    assert gated.loc[0, "target"] == "AluY"
    assert gated.loc[1, "target"] == ""
    assert float(gated.loc[1, "qcov"]) == 0.0
