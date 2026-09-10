#!/usr/bin/env python3
"""Compare swept vs per-window mate fetch, including annotate's next-step table."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from retro_miner.mei_support import (
    _discordant_mate_mei_query_seq,
    _fetch_discordant_mate_sequences,
    _fetch_mates_per_window,
    _fetch_mates_swept,
)

_NEXT_STEP_COLS = [
    "read_name",
    "mate_chrom",
    "mate_pos",
    "mate_seq",
    "mate_ref_start",
    "mate_ref_end",
    "mate_soft_clip_side",
    "mate_soft_clip_len",
    "mate_soft_clip_seq",
    "mei_query_seq",
]


def _with_mei_query(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["mei_query_seq"] = [
        _discordant_mate_mei_query_seq(row) for row in out.itertuples(index=False)
    ]
    return out


def _next_step_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in _NEXT_STEP_COLS if c in df.columns]
    return df.loc[:, cols].reset_index(drop=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--bam", type=Path, required=True)
    p.add_argument("--discordant-tsv", type=Path, required=True)
    p.add_argument(
        "--nrows",
        type=int,
        default=5000,
        help="Discordant rows to read. 0 = entire table.",
    )
    args = p.parse_args()

    read_kw: dict = {"sep": "\t"}
    if args.nrows > 0:
        read_kw["nrows"] = args.nrows
    df = pd.read_csv(args.discordant_tsv, **read_kw)
    print(f"rows={len(df)} bam={args.bam}")

    t0 = time.monotonic()
    per_window = _fetch_discordant_mate_sequences(
        df, args.bam, fetch_fn=_fetch_mates_per_window
    )
    per_s = time.monotonic() - t0

    t1 = time.monotonic()
    swept = _fetch_discordant_mate_sequences(df, args.bam, fetch_fn=_fetch_mates_swept)
    sweep_s = time.monotonic() - t1

    per_next = _next_step_table(_with_mei_query(per_window))
    sweep_next = _next_step_table(_with_mei_query(swept))

    if list(per_next.columns) != list(sweep_next.columns):
        raise SystemExit(f"COLUMN_MISMATCH {list(per_next.columns)} vs {list(sweep_next.columns)}")
    if len(per_next) != len(sweep_next):
        raise SystemExit(f"ROWCOUNT_MISMATCH {len(per_next)} vs {len(sweep_next)}")

    mismatches = []
    for col in per_next.columns:
        left = per_next[col]
        right = sweep_next[col]
        if left.dtype != right.dtype:
            left = left.astype(str)
            right = right.astype(str)
        bad = left.ne(right) & ~(left.isna() & right.isna())
        n_bad = int(bad.sum())
        if n_bad:
            mismatches.append(f"{col}={n_bad}")
    if mismatches:
        raise SystemExit(f"MISMATCH columns: {', '.join(mismatches)}")

    fetched = int(sweep_next["mate_seq"].fillna("").astype(str).str.len().gt(0).sum())
    next_queries = int(sweep_next["mei_query_seq"].fillna("").astype(str).str.len().gt(0).sum())
    print(
        f"MATCH next_step_rows={len(sweep_next)} mate_seq_filled={fetched} "
        f"mei_query_nonempty={next_queries} per_window_s={per_s:.1f} "
        f"sweep_s={sweep_s:.1f} speedup={per_s / sweep_s:.1f}x"
    )


if __name__ == "__main__":
    main()
