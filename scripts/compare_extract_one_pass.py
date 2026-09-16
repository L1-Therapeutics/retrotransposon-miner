#!/usr/bin/env python3
"""Compare one-pass extract tables to an existing two-pass outdir."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Repo-root sys.path resolution for standalone execution
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from retro_miner.evidence_extract import extract_split_and_discordant_evidence


def _load(path: Path) -> pd.DataFrame:
    """Load and sort an evidence TSV for comparison."""
    df = pd.read_csv(path, sep="\t")
    return df.sort_values(list(df.columns), kind="mergesort").reset_index(drop=True)


def main() -> None:
    """Compare 1-pass extraction output against a 2-pass baseline directory."""
    parser = argparse.ArgumentParser(description="Compare one-pass vs two-pass extract tables.")
    parser.add_argument("--bam", type=Path, required=True, help="Input BAM path")
    parser.add_argument("--sample", default="disease", help="Sample name")
    parser.add_argument("--region", default="chr22", help="Region string")
    parser.add_argument("--baseline-dir", type=Path, required=True, help="2-pass baseline output dir")
    parser.add_argument("--outdir", type=Path, required=True, help="1-pass output dir")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()
    split_s, disc_s = extract_split_and_discordant_evidence(
        bam_path=args.bam,
        sample_name=args.sample,
        outdir=args.outdir,
        regions=args.region,
        min_mapq=20,
        min_mapq_discordant=20,
        fetch_mate_seq=False,
    )
    elapsed = time.monotonic() - t0
    print(
        f"one_pass_elapsed_s={elapsed:.1f} split_rows={split_s.split_evidence_rows} "
        f"disc_rows={disc_s.discordant_evidence_rows} "
        f"insert_threshold={disc_s.insert_size_threshold} "
        f"weak_only={disc_s.weak_only_discordant_filtered_rows}"
    )

    for kind in ("split_evidence", "discordant_evidence"):
        base = _load(args.baseline_dir / f"{kind}.{args.sample}.tsv")
        new = _load(args.outdir / f"{kind}.{args.sample}.tsv")
        pd.testing.assert_frame_equal(base, new, check_dtype=False, obj=kind)
        print(f"MATCH {kind} rows={len(new)}")


if __name__ == "__main__":
    main()
