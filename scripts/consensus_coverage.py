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

The axis length used to come from the data: the furthest coordinate any read reached,
because the panel was not on disk. That made the profile's upper bound a property of
the reads being measured, so "the reads reach 100% of the consensus" was unanswerable
by construction. Axis lengths are now read from the panels
(`scripts/fetch_dfam_mei_panel.py`; hashes in `retrotransposon_db/dfam/*.manifest.json`),
on the same full-length axis the coordinates live on: the fixture's target
`L1HS_full#LINE/L1` is 6018 bp, which is also the largest mei_end any read reports.
`axis_source` says which of the two supplied the number, and every locus whose reads
reach past the panel axis is flagged rather than silently trimmed.
Reads are never counted on the wrong axis without saying so.

Usage:
    python consensus_coverage.py <repo-root> [--out-dir DIR] [--public-data-dir DIR]
                                             [--mei-panel FASTA]...
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_PACKAGE_ROOTS = tuple(
    candidate
    for parent in Path(__file__).resolve().parents
    for candidate in (parent / "src", parent / "retrotransposon-miner" / "src")
)
for _root in _PACKAGE_ROOTS:
    if (_root / "retro_miner" / "__init__.py").exists():
        if str(_root) not in sys.path:
            sys.path.insert(0, str(_root))
        break

from retro_miner.mei_panel_index import (  # noqa: E402
    find_panel_paths,
    resolve_consensus_lengths,
)

READ_COLS = ("mei_target", "mei_start", "mei_end")
MATE_COLS = ("mate_mei_target", "mate_mei_start", "mate_mei_end")

OBSERVED_SOURCE = "data_observed_max_end"


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


def profile_locus(locus_dir: str | Path,
                  panel_lengths: dict[str, object] | None = None) -> dict | None:
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
    observed_end = int(max(end for _, end in spans))
    hit = (panel_lengths or {}).get(target.split("#")[0])
    if hit is not None:
        length = int(hit.bp)          # type: ignore[attr-defined]
        axis_source = f"panel:{hit.basis}"   # type: ignore[attr-defined]
    else:
        # No panel record for this target: the old data-derived lower bound, named
        # as such so a reader knows the upper end of the axis is not measured.
        length = observed_end
        axis_source = OBSERVED_SOURCE
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
        "axis_source": axis_source,
        "axis_observed_max_end": observed_end,
        # Reads reaching past the panel axis would mean the target and the axis disagree,
        # which is the failure this script exists to expose, so it is flagged.
        "reads_beyond_axis_bp": max(observed_end - length, 0),
        "reported_span": int(gold.get("consensus_insertion_mei_span", pd.Series([0]))[0]),
        "reported_span_full": int(gold.get("consensus_insertion_mei_span_full", pd.Series([0]))[0]),
        "depth": depth,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("repo_root")
    ap.add_argument("--out-dir", default="coverage_out")
    ap.add_argument("--public-data-dir", default=None,
                    help="root holding retrotransposon_db/{dfam,full_consensus}; defaults to "
                         "$RTM_PUBLIC_DATA_DIR or the working directory")
    ap.add_argument("--mei-panel", action="append", default=[], type=Path,
                    help="explicit Dfam MEI panel FASTA (repeatable); overrides discovery")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    locus_dirs = sorted(glob.glob(os.path.join(args.repo_root, "tests/fixtures/loci/*")))
    targets: set[str] = set()
    for locus_dir in locus_dirs:
        gold_path = os.path.join(locus_dir, "gold_locus.tsv")
        if not os.path.exists(gold_path):
            continue
        gold = pd.read_csv(gold_path, sep="\t")
        target = str(gold.get("consensus_mei_subfamily", pd.Series([""]))[0])
        if target and target != "nan":
            targets.add(target.split("#")[0])

    panels = find_panel_paths(list(args.mei_panel), public_data_dir=args.public_data_dir)
    panel_lengths = resolve_consensus_lengths(sorted(targets), panels) if panels else {}

    results = [r for r in (profile_locus(d, panel_lengths) for d in locus_dirs) if r]
    if not results:
        raise SystemExit("no locus fixture carried both supporting_reads_detail.mei.tsv "
                         "and gold_locus.tsv")
    table = pd.DataFrame([{k: v for k, v in r.items() if k != "depth"} for r in results])
    table["cov_pct_of_axis"] = (100 * table.cov_bp / table.axis_len).round(1)
    table.to_csv(out / "consensus_coverage.tsv", sep="\t", index=False)
    np.savez(out / "coverage_profiles.npz",
             **{r["locus"]: r["depth"] for r in results})
    print(table.to_string(index=False))
    print(f"\npanels: {[str(path) for _kind, path in panels] or '(none found)'}")
    print(f"axis resolved from a panel for {int((table.axis_source != OBSERVED_SOURCE).sum())}"
          f"/{len(table)} loci"
          + (f"; the rest fell back to the observed read extent: "
             f"{sorted(set(table.loc[table.axis_source == OBSERVED_SOURCE, 'subfamily']))}"
             if (table.axis_source == OBSERVED_SOURCE).any() else ""))
    beyond = table[table.reads_beyond_axis_bp > 0]
    if len(beyond):
        print(f"loci whose reads reach past the panel axis: {len(beyond)} "
              f"{beyond[['locus', 'subfamily', 'axis_len', 'axis_observed_max_end']].to_dict('records')}")
    print(f"\nwrote {out}/consensus_coverage.tsv and {out}/coverage_profiles.npz")


if __name__ == "__main__":
    main()
