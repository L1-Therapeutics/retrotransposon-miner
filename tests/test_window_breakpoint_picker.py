"""In-window MEI breakpoint pile selection: split > TSD > polyA > DPE."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _choose_window_breakpoints,
    _derive_breakpoint_interval_fields,
    _WINDOW_BREAKPOINT_PILE_GAP_BP,
)


WINDOW = {
    "chrom": "chr22",
    "window_start": 49878612,
    "window_end": 49880399,
}

LEFT_PILE = 49879145
RIGHT_PILE = 49879732  # nssv14073986 / sentinel catalog site


def _candidate(**extra) -> pd.DataFrame:
    row = {
        **WINDOW,
        "tsd_left_breakpoint": 0,
        "tsd_right_breakpoint": 0,
        "tsd_len_estimate": 0,
        "tsd_detected": False,
        "tsd_evidence_source": "",
    }
    row.update(extra)
    return pd.DataFrame([row])


def _split_row(pos: int, name: str, *, mei: bool = False, polya: bool = False) -> dict[str, object]:
    return {
        **WINDOW,
        "pos": pos,
        "read_name": name,
        "clip_side": "L",
        "clip_len": 40,
        "mei_hit": mei,
        "clip_poly_at_run": 12 if polya else 0,
        "poly_tail_rescued": polya,
    }


def _dpe_row(pos: int, name: str) -> dict[str, object]:
    return {
        **WINDOW,
        "pos": pos,
        "read_name": name,
        "mei_hit": True,
        "mate_mei_hit": True,
    }


def _pick(candidates, *, split=None, disc=None):
    return _choose_window_breakpoints(
        candidates,
        split_frames=[split if split is not None else pd.DataFrame()],
        discordant_frames=[disc if disc is not None else pd.DataFrame()],
    )


def test_sentinel_split_pile_beats_heavier_dpe_pile():
    """chr22:49879732 split MEI wins over the nearby DPE-heavy left pile."""
    split = pd.DataFrame(
        [_split_row(RIGHT_PILE, f"sr{i}", mei=True) for i in range(4)]
    )
    disc = pd.DataFrame(
        [_dpe_row(LEFT_PILE, f"dpe{i}") for i in range(17)]
    )
    picked = _pick(_candidate(), split=split, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == RIGHT_PILE
    assert picked.iloc[0]["breakpoint_evidence_source"] == "split_mei"


def test_tsd_pile_beats_distant_dpe_when_no_split():
    tsd_left, tsd_right = 200, 215
    candidates = _candidate(
        chrom="chr1",
        window_start=100,
        window_end=1000,
        tsd_left_breakpoint=tsd_left,
        tsd_right_breakpoint=tsd_right,
        tsd_len_estimate=16,
        tsd_detected=True,
        tsd_evidence_source="tsd_disease",
    )
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 100,
                "window_end": 1000,
                "pos": 800,
                "read_name": f"dpe{i}",
                "mei_hit": True,
            }
            for i in range(20)
        ]
    )
    picked = _pick(candidates, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == (tsd_left + tsd_right) // 2
    assert picked.iloc[0]["breakpoint_evidence_source"] == "tsd_disease"


def test_polya_pile_beats_distant_dpe_when_no_split_or_tsd():
    split = pd.DataFrame(
        [_split_row(500, f"pa{i}", polya=True) for i in range(3)]
    )
    split["chrom"] = "chr1"
    split["window_start"] = 100
    split["window_end"] = 1000
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 100,
                "window_end": 1000,
                "pos": 900,
                "read_name": f"dpe{i}",
                "mei_hit": True,
            }
            for i in range(10)
        ]
    )
    picked = _pick(
        _candidate(chrom="chr1", window_start=100, window_end=1000),
        split=split,
        disc=disc,
    )
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == 500
    assert picked.iloc[0]["breakpoint_evidence_source"] == "polyA"


def test_dpe_only_picks_the_pile_with_more_mei_reads():
    disc = pd.DataFrame(
        [_dpe_row(LEFT_PILE, f"l{i}") for i in range(17)]
        + [_dpe_row(RIGHT_PILE, f"r{i}") for i in range(4)]
    )
    picked = _pick(_candidate(), disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == LEFT_PILE
    assert picked.iloc[0]["breakpoint_evidence_source"] == "dpe_mei"


def test_derive_interval_stays_on_chosen_pile_not_distant_mode():
    df = pd.DataFrame(
        [
            {
                **WINDOW,
                "insertion_breakpoint_pos": RIGHT_PILE,
                "tsd_left_breakpoint": 0,
                "tsd_right_breakpoint": 0,
                "dpe_gap_left": 0,
                "dpe_gap_right": 0,
                "disease_L_mei_breakpoint_mode": LEFT_PILE,
                "disease_R_mei_breakpoint_mode": RIGHT_PILE,
                "disease_L_split_breakpoint_mode": LEFT_PILE,
                "disease_R_split_breakpoint_mode": RIGHT_PILE,
                "control_L_mei_breakpoint_mode": 0,
                "control_R_mei_breakpoint_mode": 0,
                "control_L_split_breakpoint_mode": 0,
                "control_R_split_breakpoint_mode": 0,
                "dpe_soft_clip_mode": 0,
            }
        ]
    )
    out = _derive_breakpoint_interval_fields(
        df,
        breakpoint_pos_col="insertion_breakpoint_pos",
        output_prefix="insertion_",
    )
    lo = int(out.iloc[0]["insertion_breakpoint_interval_start"])
    hi = int(out.iloc[0]["insertion_breakpoint_interval_end"])
    bp = int(out.iloc[0]["insertion_breakpoint_pos"])
    assert bp == RIGHT_PILE
    assert abs(lo - RIGHT_PILE) <= _WINDOW_BREAKPOINT_PILE_GAP_BP
    assert abs(hi - RIGHT_PILE) <= _WINDOW_BREAKPOINT_PILE_GAP_BP
    assert (hi - lo) < abs(RIGHT_PILE - LEFT_PILE)
