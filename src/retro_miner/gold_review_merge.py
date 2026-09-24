"""Concatenate per-chromosome gold-review tables and re-rank them as one genome.

Used at the end of ``--chr all``. The sort is ``_prioritize_mei_candidates``,
the same order a single chromosome already writes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd

from retro_miner.mei_support import _prioritize_mei_candidates

# Columns the review sort treats as booleans. CSV re-reads them as strings,
# and bool("False") is True, which would scramble stage and catalog flags.
_BOOL_COLUMNS = (
    "gold_stage_pass",
    "silver_stage_pass",
    "known_mei_polymorphism",
    "consensus_tsd_detected",
    "tsd_detected",
    "consensus_poly_at_supported",
    "tsd_poly_at_filter_applied",
    "poly_at_supported",
    "two_sided_support",
)


def _parse_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "t", "yes"})


def merge_gold_review_tables(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Stack gold-review tables and sort them like one chromosome."""
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(list(frames), ignore_index=True)
    for col in _BOOL_COLUMNS:
        if col in merged.columns:
            merged[col] = _parse_bool(merged[col])
    ranked = _prioritize_mei_candidates(merged, stage_first=True)
    for col in ranked.columns:
        if ranked[col].dtype == bool:
            ranked[col] = ranked[col].astype(int)
    return ranked


def write_merged_gold_review(inputs: Sequence[Path], output: Path) -> int:
    """Read per-chromosome TSVs, re-rank, and write one genome-wide table."""
    frames: list[pd.DataFrame] = []
    for path in inputs:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        frame = frame.loc[:, ~frame.columns.duplicated()]
        for col in frame.columns:
            frame[col] = frame[col].where(frame[col] != "", other=pd.NA)
        frames.append(frame)
    ranked = merge_gold_review_tables(frames)
    if ranked.empty:
        return 0
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(output, sep="\t", index=False)
    return int(len(ranked))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Re-rank per-chromosome gold reviews as one genome.")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args(argv)
    missing = [str(path) for path in args.inputs if not path.is_file()]
    if missing:
        raise SystemExit(f"gold review table not found: {missing[0]}")
    n_rows = write_merged_gold_review(args.inputs, args.output)
    print(f"[gold-review-merge] rows={n_rows} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
