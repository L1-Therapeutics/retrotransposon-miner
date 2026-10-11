"""How much of a reported LINE-1 span is supported by read alignment?

Brandler, 2026-10-09: "we expect some line-1 insertions to be full length but this
should be based on read alignment and not imputed."

Background. In `_sample_insertion_span_and_orientation` (mei_support.py) there are
branches that set the far end of the consensus footprint to the subfamily's consensus
length when one side has a strong anchor and the other has a strong poly-A run:

    elif l_strong and r_poly_strong:
        tlen = max(l_target_len, r_target_len)
        raw_end = max(tlen, max(l_start, l_end))

A span produced that way ends exactly at the consensus length without a read ever
reaching it. A genuinely full-length insertion ends there too, so "ends at the
consensus length" is necessary but NOT sufficient evidence of imputation. This script
measures how large the ambiguous set is; it cannot resolve it from VCF fields alone.

CAVEAT on consensus lengths. The Dfam panel FASTA is not available locally, so the
consensus length of each subfamily is ESTIMATED as the largest MEI_3P observed for that
subfamily across the whole cohort. With thousands of calls this is a good estimate for
common subfamilies and a poor one for rare ones. Replace with the real panel lengths
(`retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta`) before quoting anything.

Usage:
    python line1_span_support.py <dir-of-vcfs> [--out-dir DIR] [--min-full-span 5500]
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path

import pandas as pd

SUPPORT_KEYS = ("SR_L", "SR_R", "DPE_L", "DPE_R", "MEI_MAPPED", "polyA_MAPPED")
NUMERIC_INFO = ("MEI_SPAN", "MEI_5P", "MEI_3P", "POLYA_MIN_BP")


def load_cohort(vcf_dir: str | Path) -> pd.DataFrame:
    """One row per call across every *.vcf in the directory, with SUPPORT parsed out."""
    from retro_miner.gene_annotation import read_vcf

    paths = sorted(glob.glob(os.path.join(str(vcf_dir), "*.vcf")))
    if not paths:
        raise FileNotFoundError(f"no *.vcf in {vcf_dir}")
    frames = []
    for path in paths:
        _, records = read_vcf(path)
        frame = pd.DataFrame([{"chrom": r.chrom, "pos": r.pos, **r.info} for r in records])
        frame["sample"] = os.path.basename(path).split("_")[0]
        frames.append(frame)
    calls = pd.concat(frames, ignore_index=True)
    for col in NUMERIC_INFO:
        calls[col] = pd.to_numeric(calls.get(col), errors="coerce")
    support = calls.get("SUPPORT", pd.Series("", index=calls.index)).fillna("")
    for key in SUPPORT_KEYS:
        calls[key] = pd.to_numeric(
            support.str.extract(rf"(?:^|\|){re.escape(key)}:(\d+)")[0], errors="coerce").fillna(0)
    calls["polyA_side"] = support.str.extract(r"polyA_side:(\w)")[0]
    calls["subfamily"] = calls["MEISUBFAMILY"].fillna("").str.split("#").str[0]
    calls["sentinel"] = (calls["MEI_SPAN"] <= 0) | (calls["MEI_5P"] < 0)
    return calls


def estimate_consensus_lengths(calls: pd.DataFrame) -> pd.Series:
    """Largest MEI_3P seen per subfamily. A PROXY for the Dfam consensus length."""
    usable = calls[~calls.sentinel]
    return usable.groupby("subfamily")["MEI_3P"].max().rename("consensus_len")


def annotate(calls: pd.DataFrame) -> pd.DataFrame:
    usable = calls[~calls.sentinel].join(estimate_consensus_lengths(calls), on="subfamily")
    usable["ends_at_consensus"] = usable["MEI_3P"] == usable["consensus_len"]
    usable["span_fraction"] = usable["MEI_SPAN"] / usable["consensus_len"]
    return usable


def summarise(calls: pd.DataFrame, usable: pd.DataFrame, min_full_span: int) -> dict:
    line1 = usable[usable.MEIFAMILY == "LINE1"]
    full = line1[line1.MEI_SPAN >= min_full_span]
    return {
        "genomes": int(calls["sample"].nunique()),
        "calls_total": len(calls),
        "calls_by_family": calls.MEIFAMILY.value_counts().to_dict(),
        "sentinels_by_family": calls[calls.sentinel].MEIFAMILY.value_counts().to_dict(),
        "line1_usable": len(line1),
        "line1_apparently_full": len(full),
        "line1_full_ending_at_consensus": int(full.ends_at_consensus.sum()),
        "line1_full_per_genome": full.groupby("sample").size().to_dict(),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("--out-dir", default="line1_span_out")
    ap.add_argument("--min-full-span", type=int, default=5500,
                    help="span at or above which a LINE-1 call is treated as apparently full length")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    calls = load_cohort(args.vcf_dir)
    usable = annotate(calls)
    stats = summarise(calls, usable, args.min_full_span)

    line1 = usable[usable.MEIFAMILY == "LINE1"]
    bands = [(20, 99), (100, 999), (1000, args.min_full_span - 1), (args.min_full_span, 10 ** 9)]
    rows = []
    for lo, hi in bands:
        sel = line1[line1.MEI_SPAN.between(lo, hi)]
        if len(sel):
            rows.append({"band": f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}",
                         "n": len(sel),
                         "ends_at_consensus": int(sel.ends_at_consensus.sum()),
                         "pct": round(100 * sel.ends_at_consensus.mean(), 1)})
    pd.DataFrame(rows).to_csv(out / "line1_span_bands.tsv", sep="\t", index=False)
    usable[usable.MEIFAMILY == "LINE1"].to_csv(out / "line1_calls.tsv", sep="\t", index=False)

    print(f"genomes {stats['genomes']} | calls {stats['calls_total']} | {stats['calls_by_family']}")
    print(f"sentinels {stats['sentinels_by_family']}")
    print(f"\nLINE-1 with a usable footprint: {stats['line1_usable']}")
    print(f"apparently full length (span >= {args.min_full_span}): {stats['line1_apparently_full']}")
    print(f"  of which the 3' end sits exactly at the estimated consensus length: "
          f"{stats['line1_full_ending_at_consensus']}")
    print("\nby span band:")
    for row in rows:
        print(f"  {row['band']:>10s}  n={row['n']:5d}  ends at consensus {row['ends_at_consensus']:5d} "
              f"({row['pct']}%)")
    print("\nper genome, apparently full-length LINE-1:")
    for sample, n in sorted(stats["line1_full_per_genome"].items()):
        print(f"  {sample}  {n}")
    print(f"\nwrote {out}/line1_span_bands.tsv and {out}/line1_calls.tsv")
    print("NOTE: consensus lengths are estimated from the data, not read from the Dfam panel.")


if __name__ == "__main__":
    main()
