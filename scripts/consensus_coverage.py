"""What part of a mobile element's consensus do the supporting reads actually reach?

Each locus fixture under tests/fixtures/loci/ ships supporting_reads_detail.mei.tsv,
which carries per-read coordinates on the consensus axis of whatever panel target the
read aligned to (mei_start / mei_end), plus the same for the mate (mate_mei_start /
mate_mei_end). Those columns are the read-level evidence that a reported span either
does or does not rest on.

For the subfamily named in gold_locus.tsv this builds a coverage profile over the
consensus axis and reports the covered extent, the covered fraction, and the largest
internal gap. Reads aligning to a different subfamily are reported separately rather
than merged, because mixing targets is how a span can exceed its own consensus length.

Usage:
    python consensus_coverage.py <repo-root> [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

READ_COLS = ("mei_target", "mei_start", "mei_end")
MATE_COLS = ("mate_mei_target", "mate_mei_start", "mate_mei_end")


def read_intervals(detail: pd.DataFrame, target: str) -> list[tuple[int, int]]:
    """Consensus intervals contributed by reads and by mates aligning to `target`."""
    spans: list[tuple[int, int]] = []
    for tcol, scol, ecol in (READ_COLS, MATE_COLS):
        if tcol not in detail.columns:
            continue
        hit = detail[detail[tcol].fillna("") == target].dropna(subset=[scol, ecol])
        spans.extend((int(s), int(e)) for s, e in zip(hit[scol], hit[ecol]) if e >= s)
    return spans


def coverage_profile(spans: list[tuple[int, int]], length: int) -> np.ndarray:
    depth = np.zeros(max(length, 1) + 1, dtype=int)
    for start, end in spans:
        depth[max(start, 1):min(end, length) + 1] += 1
    return depth[1:]


def largest_gap(depth: np.ndarray) -> int:
    covered = np.flatnonzero(depth > 0)
    if covered.size < 2:
        return 0
    return int(np.max(np.diff(covered)) - 1)


def profile_locus(locus_dir: str | Path) -> dict | None:
    locus_dir = Path(locus_dir)
    detail_path = locus_dir / "supporting_reads_detail.mei.tsv"
    gold_path = locus_dir / "gold_locus.tsv"
    manifest_path = locus_dir / "manifest.json"
    if not (detail_path.exists() and gold_path.exists()):
        return None
    detail = pd.read_csv(detail_path, sep="\t")
    gold = pd.read_csv(gold_path, sep="\t")
    target = str(gold.get("consensus_mei_subfamily", pd.Series([""]))[0])
    if not target or target == "nan":
        return None

    spans = read_intervals(detail, target)
    if not spans:
        return None
    # The Dfam panel is not shipped, so the consensus length of THIS target is bounded
    # from below by the furthest coordinate any read reached on it. Coordinates from
    # other subfamilies are deliberately excluded: mixing targets is the bug this script
    # exists to expose, not a source of axis length.
    length = int(max(end for _, end in spans))
    depth = coverage_profile(spans, length)
    covered = np.flatnonzero(depth > 0)

    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    return {
        "locus": manifest.get("catalog_id", locus_dir.name),
        "subfamily": target.split("#")[0],
        "family": ("LINE1" if target.startswith("L1") else
                   "ALU" if target.startswith("Alu") else "SVA"),
        "n_intervals": len(spans),
        "cov_lo": int(covered.min()) + 1,
        "cov_hi": int(covered.max()) + 1,
        "cov_bp": int((depth > 0).sum()),
        "largest_internal_gap": largest_gap(depth),
        # assumption-free: consensus positions upstream of the first covered base that
        # no read reached. Needs no consensus length at all.
        "unobserved_5prime_bp": int(covered.min()),
        "axis_len": length,
        "reported_span": int(gold.get("consensus_insertion_mei_span", pd.Series([0]))[0]),
        "reported_span_full": int(gold.get("consensus_insertion_mei_span_full", pd.Series([0]))[0]),
        "depth": depth,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("repo_root")
    ap.add_argument("--out-dir", default="coverage_out")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = [r for r in (profile_locus(d) for d in
                           sorted(glob.glob(os.path.join(args.repo_root,
                                                         "tests/fixtures/loci/*")))) if r]
    if not results:
        raise SystemExit("no locus fixture carried both supporting_reads_detail.mei.tsv "
                         "and gold_locus.tsv")
    table = pd.DataFrame([{k: v for k, v in r.items() if k != "depth"} for r in results])
    table["cov_pct_of_axis"] = (100 * table.cov_bp / table.axis_len).round(1)
    table.to_csv(out / "consensus_coverage.tsv", sep="\t", index=False)
    np.savez(out / "coverage_profiles.npz",
             **{r["locus"]: r["depth"] for r in results})
    print(table.to_string(index=False))
    print(f"\nwrote {out}/consensus_coverage.tsv and {out}/coverage_profiles.npz")


if __name__ == "__main__":
    main()
