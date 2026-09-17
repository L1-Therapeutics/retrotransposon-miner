#!/usr/bin/env python3
"""Compare two annotate-through-gold table pairs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _tier_counts(df: pd.DataFrame) -> str:
    for col in ("review_tier", "tier", "call_tier", "gold_tier"):
        if col in df.columns:
            vc = df[col].value_counts(dropna=False)
            return f"{col} " + ", ".join(f"{k}={v}" for k, v in vc.items())
    return "no_tier_column"


def _compare(left: pd.DataFrame, right: pd.DataFrame, label: str) -> None:
    print(f"=== {label} ===")
    print(f"left_rows={len(left)} right_rows={len(right)}")
    print(f"left_cols={len(left.columns)} right_cols={len(right.columns)}")
    if list(left.columns) != list(right.columns):
        only_l = [c for c in left.columns if c not in right.columns]
        only_r = [c for c in right.columns if c not in left.columns]
        print(f"COLUMN_MISMATCH only_left={only_l} only_right={only_r}")
        cols = [c for c in left.columns if c in right.columns]
        left = left.loc[:, cols]
        right = right.loc[:, cols]
    if len(left) != len(right):
        raise SystemExit(f"ROWCOUNT_MISMATCH {label}: {len(left)} vs {len(right)}")

    mismatches = []
    for col in left.columns:
        bad = left[col].fillna("").astype(str).ne(right[col].fillna("").astype(str))
        n_bad = int(bad.sum())
        if n_bad:
            mismatches.append((col, n_bad))
    if mismatches:
        print("MISMATCH " + ", ".join(f"{c}={n}" for c, n in mismatches))
        raise SystemExit(f"MISMATCH {label}")
    print(f"MATCH {label} rows={len(left)}")
    print(_tier_counts(left))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--baseline-dir", type=Path, required=True)
    p.add_argument("--new-dir", type=Path, required=True)
    args = p.parse_args()
    for name in ("candidate_loci.mei.tsv", "candidate_loci.mei.gold_review.tsv"):
        _compare(_load(args.baseline_dir / name), _load(args.new_dir / name), name)


if __name__ == "__main__":
    main()
