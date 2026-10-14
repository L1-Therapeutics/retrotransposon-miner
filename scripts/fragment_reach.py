"""How far can a paired-end library reach into a mobile element insertion?

Brandler, 2026-10-09: LINE-1 full-length spans "should be based on read alignment and
not imputed." Whether that is achievable is bounded by physics, not by code: a read pair
anchored in unique genomic sequence beside the insertion can only reach as far into the
inserted sequence as its fragment is long.

This measures the fragment-length distribution directly from a locus fixture BAM and
converts it into a reach ceiling, so the limit is a measured number rather than a
remembered rule of thumb.

Definitions, stated because they change the answer:

  template_length (TLEN)  signed outer distance between the leftmost and rightmost
                          mapped bases of a pair. Only read from properly-paired,
                          primary, non-duplicate records, where it is meaningful.
  reach                   the furthest a mate can sit from the junction while its
                          partner still anchors in unique flanking sequence. Bounded
                          above by the fragment length; smaller once the anchor needs
                          a usable number of bases outside the element.

Usage:
    python fragment_reach.py <locus-fixture-dir> [--element-length 6018] [--anchor-bp 30]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import numpy as np
import pysam


def fragment_lengths(bam_path: str | Path) -> np.ndarray:
    """|TLEN| for properly-paired primary records, deduplicated by read name.

    Both mates of a pair report the same |TLEN|, so counting each record would
    double-count every fragment and halve the standard error.
    """
    seen: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_unmapped or read.mate_is_unmapped:
                continue
            tlen = abs(int(read.template_length or 0))
            if tlen > 0:
                seen.setdefault(read.query_name, tlen)
    return np.array(sorted(seen.values()), dtype=int)


def summarise(lengths: np.ndarray, element_length: int, anchor_bp: int) -> dict:
    if lengths.size == 0:
        raise SystemExit("no properly-paired records with a usable TLEN in this BAM")
    pct = {p: float(np.percentile(lengths, p)) for p in (50, 90, 99)}
    reach_cap = pct[99] - anchor_bp          # mate reach once the anchor is accounted for
    return {
        "n_fragments": int(lengths.size),
        "median": pct[50],
        "p90": pct[90],
        "p99": pct[99],
        "max": int(lengths.max()),
        "reach_cap_bp": reach_cap,
        "element_length": element_length,
        "reachable_fraction": reach_cap / element_length if element_length else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("locus_dir", help="a tests/fixtures/loci/<id> directory")
    ap.add_argument("--element-length", type=int, default=6018,
                    help="consensus length of the element at this locus (L1HS is 6018)")
    ap.add_argument("--anchor-bp", type=int, default=30,
                    help="bases the anchoring mate must keep in unique flanking sequence")
    args = ap.parse_args()

    bams = sorted(glob.glob(os.path.join(args.locus_dir, "*.bam")))
    if not bams:
        raise SystemExit(f"no BAM in {args.locus_dir}")
    lengths = fragment_lengths(bams[0])
    stats = summarise(lengths, args.element_length, args.anchor_bp)

    print(f"BAM: {os.path.basename(bams[0])}")
    print(f"properly-paired fragments: {stats['n_fragments']}")
    print(f"fragment length  median {stats['median']:.0f}  p90 {stats['p90']:.0f}  "
          f"p99 {stats['p99']:.0f}  max {stats['max']}")
    print(f"\nreach ceiling, allowing {args.anchor_bp} bp of anchor in unique sequence:")
    print(f"  a mate can sit at most {stats['reach_cap_bp']:.0f} bp inside the insertion")
    print(f"  that is {100 * stats['reachable_fraction']:.1f}% of a {args.element_length} bp element")
    print(f"  leaving {args.element_length - stats['reach_cap_bp']:.0f} bp unreachable "
          f"from either junction by any pair in this library")
    print("\n" + json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
