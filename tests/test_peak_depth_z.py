"""Peak-depth z-score is stored on gold and drops calls at z >= 2."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from retro_miner.mei_support import _assign_gold_stage, _build_gold_review_table

FIXTURE = (
    Path(__file__).resolve().parent
    / "fixtures"
    / "loci"
    / "rank010_chr19_3989864"
    / "gold_locus.tsv"
)


def _pileup_table() -> pd.DataFrame:
    base = pd.read_csv(FIXTURE, sep="\t").iloc[[0]].copy()
    base["silver_stage_pass"] = True
    base["known_mei_polymorphism"] = False
    rows = []
    for peak in [40.0] * 30 + [4000.0, 4000.0]:
        row = base.copy()
        row["disease_local_bam_peak_depth"] = peak
        row["control_local_bam_peak_depth"] = peak
        rows.append(row)
    out = pd.concat(rows, ignore_index=True)
    out.loc[out.index[-1], "known_mei_polymorphism"] = True
    return out


def test_peak_depth_z_is_on_gold_review_and_drops_at_two() -> None:
    scored = _pileup_table()
    out = _assign_gold_stage(scored, empirical_stage=False, min_mei_mapped=3)
    z = pd.to_numeric(out["local_bam_peak_depth_z"], errors="coerce")
    assert z.iloc[0] < 2.0
    assert bool(out["gold_stage_pass"].iloc[0])
    assert z.iloc[-2] >= 2.0
    assert not bool(out["gold_stage_pass"].iloc[-2])
    assert "depth_pileup_artifact" in str(out["gold_stage_fail_reason"].iloc[-2])
    assert z.iloc[-1] >= 2.0
    assert not bool(out["gold_stage_pass"].iloc[-1])
    assert "depth_pileup_artifact" in str(out["gold_stage_fail_reason"].iloc[-1])
    review = _build_gold_review_table(out, empirical_stage=False)
    assert "local_bam_peak_depth" in review.columns
    assert "local_bam_peak_depth_z" in review.columns
    assert float(review["local_bam_peak_depth_z"].iloc[-2]) >= 2.0
