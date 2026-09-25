"""In-window MEI breakpoint pile selection: weighted split > TSD > polyA > DPE."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pysam

from retro_miner.mei_support import (
    _choose_window_breakpoints,
    _derive_breakpoint_interval_fields,
    _score_breakpoint_pile,
    _WINDOW_BREAKPOINT_PILE_GAP_BP,
)

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "loci" / "nssv14073986"
SENTINEL_POLYA_READ = "A00297:35:HFKWWDSXX:2:1678:17761:4006"


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


def test_weights_keep_split_ahead_of_similar_dpe():
    assert _score_breakpoint_pile(n_split=4, n_polya=0, n_dpe=0, n_clip=0, has_tsd=False) > (
        _score_breakpoint_pile(n_split=0, n_polya=0, n_dpe=17, n_clip=0, has_tsd=False)
    )
    assert _score_breakpoint_pile(n_split=0, n_polya=2, n_dpe=0, n_clip=0, has_tsd=False) < (
        _score_breakpoint_pile(n_split=0, n_polya=0, n_dpe=44, n_clip=0, has_tsd=False)
    )


def test_sentinel_split_pile_beats_heavier_dpe_pile():
    """chr22:49879732 split MEI still wins over a nearby DPE pile of similar size."""
    split = pd.DataFrame(
        [_split_row(RIGHT_PILE, f"sr{i}", mei=True) for i in range(4)]
    )
    disc = pd.DataFrame(
        [_dpe_row(LEFT_PILE, f"dpe{i}") for i in range(17)]
    )
    picked = _pick(_candidate(), split=split, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == RIGHT_PILE
    assert picked.iloc[0]["breakpoint_evidence_source"] == "split_mei"


def test_few_polya_reads_do_not_beat_heavier_dpe_pile():
    """Far DPE (beyond smear radius) can still form their own pile."""
    split = pd.DataFrame(
        [_split_row(500, f"pa{i}", polya=True) for i in range(3)]
    )
    split["chrom"] = "chr1"
    split["window_start"] = 100
    split["window_end"] = 2000
    disc = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 100,
                "window_end": 2000,
                "pos": 1500,
                "read_name": f"dpe{i}",
                "mei_hit": True,
            }
            for i in range(10)
        ]
    )
    picked = _pick(
        _candidate(chrom="chr1", window_start=100, window_end=2000),
        split=split,
        disc=disc,
    )
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == 1500
    assert picked.iloc[0]["breakpoint_evidence_source"] == "dpe_mei"


def test_dpe_smear_does_not_replace_junction_breakpoint():
    """nssv14073986: DPE anchors smear ~200 bp left of the polyA/TSD junction."""
    split = pd.DataFrame(
        [_split_row(LEFT_PILE - 100, "pa_left", polya=True), _split_row(RIGHT_PILE, "pa_true", polya=True)]
        + [_split_row(RIGHT_PILE, f"clip{i}") for i in range(4)]
    )
    smear = [49879437, 49879441, 49879477, 49879519, 49879521, 49879527, 49879527, 49879533, 49879540]
    disc = pd.DataFrame(
        [_dpe_row(pos, f"dpe{i}") for i, pos in enumerate(smear * 3)]
    )
    picked = _pick(_candidate(), split=split, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == RIGHT_PILE
    assert picked.iloc[0]["breakpoint_evidence_source"] == "polyA"


def test_sentinel_left_polya_loses_to_catalog_insertion_pile():
    """nssv14073986: two left polyA reads must not beat catalog DPE+polyA+clips."""
    split = pd.DataFrame(
        [_split_row(LEFT_PILE, f"pa{i}", polya=True) for i in range(2)]
        + [_split_row(RIGHT_PILE, "pa_true", polya=True)]
        + [_split_row(RIGHT_PILE, f"clip{i}") for i in range(4)]
    )
    disc = pd.DataFrame(
        [_dpe_row(RIGHT_PILE, f"dpe{i}") for i in range(44)]
    )
    picked = _pick(
        _candidate(
            tsd_left_breakpoint=RIGHT_PILE,
            tsd_right_breakpoint=RIGHT_PILE + 15,
            tsd_len_estimate=16,
            tsd_detected=True,
            tsd_evidence_source="tsd_disease",
        ),
        split=split,
        disc=disc,
    )
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == RIGHT_PILE + 7
    assert picked.iloc[0]["breakpoint_evidence_source"] == "tsd_disease"


def test_few_left_split_mei_lose_to_catalog_dpe_pile():
    split = pd.DataFrame(
        [_split_row(LEFT_PILE, f"sr{i}", mei=True) for i in range(2)]
    )
    disc = pd.DataFrame(
        [_dpe_row(RIGHT_PILE, f"dpe{i}") for i in range(44)]
    )
    picked = _pick(_candidate(), split=split, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == RIGHT_PILE
    assert picked.iloc[0]["breakpoint_evidence_source"] == "dpe_mei"


def test_tsd_bonus_breaks_near_tie_but_not_a_veto():
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
    disc_near = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": 100,
                "window_end": 1000,
                "pos": 800,
                "read_name": f"dpe{i}",
                "mei_hit": True,
            }
            for i in range(17)
        ]
        + [
            {
                "chrom": "chr1",
                "window_start": 100,
                "window_end": 1000,
                "pos": 208,
                "read_name": f"local{i}",
                "mei_hit": True,
            }
            for i in range(14)
        ]
    )
    picked = _pick(candidates, disc=disc_near)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == (tsd_left + tsd_right) // 2
    assert picked.iloc[0]["breakpoint_evidence_source"] == "tsd_disease"

    disc_heavy = pd.DataFrame(
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
    picked_heavy = _pick(candidates, disc=disc_heavy)
    assert int(picked_heavy.iloc[0]["insertion_breakpoint_pos"]) == 800
    assert picked_heavy.iloc[0]["breakpoint_evidence_source"] == "dpe_mei"


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


def _attach_discovery_window(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["chrom"] = WINDOW["chrom"]
    out["window_start"] = WINDOW["window_start"]
    out["window_end"] = WINDOW["window_end"]
    return out


def test_nssv14073986_real_extract_publishes_catalog_breakpoint():
    """HG00100 chr22 extract for nssv14073986, not a cartoon pile."""
    manifest = json.loads((FIXTURE_DIR / "manifest.json").read_text())
    cand = pd.read_csv(FIXTURE_DIR / "candidate.tsv", sep="\t")
    split = pd.concat(
        [
            _attach_discovery_window(pd.read_parquet(FIXTURE_DIR / "split_evidence.disease.parquet")),
            _attach_discovery_window(pd.read_parquet(FIXTURE_DIR / "split_evidence.control.parquet")),
        ],
        ignore_index=True,
    )
    disc = pd.concat(
        [
            _attach_discovery_window(pd.read_parquet(FIXTURE_DIR / "discordant_evidence.disease.parquet")),
            _attach_discovery_window(pd.read_parquet(FIXTURE_DIR / "discordant_evidence.control.parquet")),
        ],
        ignore_index=True,
    )
    picked = _pick(cand, split=split, disc=disc)
    assert int(picked.iloc[0]["insertion_breakpoint_pos"]) == int(manifest["expected_breakpoint"])
    assert picked.iloc[0]["breakpoint_evidence_source"] == "polyA"


def test_nssv14073986_bam_snippet_contains_catalog_polya_read():
    bam_path = FIXTURE_DIR / "HG00100.chr22_49877812_49881199.bam"
    bam = pysam.AlignmentFile(str(bam_path))
    try:
        names = {read.query_name for read in bam.fetch("chr22", 49879720, 49879740)}
    finally:
        bam.close()
    assert SENTINEL_POLYA_READ in names
