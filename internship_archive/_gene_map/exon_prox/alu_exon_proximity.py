"""Are fixed reference Alu copies antisense-depleted near exons?

Lev-Maor 2008 reported that antisense Alus are depleted within 150 bp upstream
of exons -- the position from which an antisense right arm can supply a 3' splice
site and become an exon (Gal-Mark 2008). Our polymorphic callset cannot test this:
only 55 unique sites fall in that window, giving 80% power against nothing smaller
than a 34% relative drop.

Fixed reference copies can. There are ~115k AluY copies in hg38 and ~67k sit
inside single-strand genes, which puts thousands in the window.

Usage:
  python alu_exon_proximity.py <rmsk.txt.gz> <gencode.gtf.gz> --out-dir out
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

WINDOW_BP = 150

# UCSC rmsk column positions (0-based), file has no header
RM_MILLIDIV, RM_CHROM, RM_START, RM_END, RM_STRAND, RM_NAME, RM_FAMILY = 2, 5, 6, 7, 9, 10, 12


def load_repeats(path: str, rep_name: str) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or f[RM_NAME] != rep_name:
                continue
            rows.append(
                (f[RM_CHROM], int(f[RM_START]), int(f[RM_END]), f[RM_STRAND], int(f[RM_MILLIDIV]))
            )
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "milli_div"])


def merge_intervals(a: np.ndarray) -> np.ndarray:
    if len(a) == 0:
        return np.empty((0, 2), dtype=np.int64)
    a = a[np.argsort(a[:, 0], kind="stable")]
    out = [a[0].copy()]
    for s, e in a[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append(np.array([s, e]))
    return np.array(out, dtype=np.int64)


def load_annotation(path: str) -> tuple[dict, dict, dict]:
    """Merged + gene, - gene, and exon intervals per chromosome (0-based half-open)."""
    plus, minus, exon = defaultdict(list), defaultdict(list), defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            chrom, s, e, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if f[2] == "exon":
                exon[chrom].append((s, e))
            elif strand == "+":
                plus[chrom].append((s, e))
            elif strand == "-":
                minus[chrom].append((s, e))
    to_iv = lambda d: {c: merge_intervals(np.array(v, dtype=np.int64)) for c, v in d.items()}
    return to_iv(plus), to_iv(minus), to_iv(exon)


def covered(iv: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Boolean: is each pos inside any [start, end) interval of iv."""
    if len(iv) == 0:
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(iv[:, 0], pos, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos), dtype=bool)
    out[ok] = pos[ok] < iv[i[ok], 1]
    return out


def nearest_exon(iv: np.ndarray, start: np.ndarray, end: np.ndarray):
    """Distance from each [start, end) element to the nearest exon, and which
    genomic side the element lies on (-1 left of exon, +1 right, 0 overlapping)."""
    n = len(start)
    if len(iv) == 0:
        return np.full(n, -1, dtype=np.int64), np.zeros(n, dtype=np.int8)
    es, ee = iv[:, 0], iv[:, 1]
    dist = np.full(n, np.iinfo(np.int64).max, dtype=np.int64)
    side = np.zeros(n, dtype=np.int8)
    # candidate exon just left of, and just right of, the element start
    j = np.searchsorted(es, start, side="right") - 1
    for off in (0, 1):
        k = j + off
        valid = (k >= 0) & (k < len(iv))
        if not valid.any():
            continue
        idx = np.where(valid)[0]
        kk = k[idx]
        # gap between the two intervals; negative means they overlap
        gap_left = es[kk] - end[idx]        # element entirely left of exon
        gap_right = start[idx] - ee[kk]     # element entirely right of exon
        d = np.maximum(gap_left, gap_right)
        s = np.where(gap_left >= gap_right, -1, 1).astype(np.int8)
        d = np.where(d < 0, 0, d)
        s = np.where(d == 0, 0, s).astype(np.int8)
        better = d < dist[idx]
        dist[idx[better]] = d[better]
        side[idx[better]] = s[better]
    return dist, side


def build(repeats: pd.DataFrame, plus: dict, minus: dict, exons: dict) -> pd.DataFrame:
    out = []
    for chrom, g in repeats.groupby("chrom", sort=False):
        if chrom not in plus and chrom not in minus:
            continue
        g = g.sort_values("start")
        start = g.start.values
        end = g.end.values
        mid = (start + end) // 2
        p = covered(plus.get(chrom, np.empty((0, 2), np.int64)), mid)
        m = covered(minus.get(chrom, np.empty((0, 2), np.int64)), mid)
        single = p ^ m
        if not single.any():
            continue
        gene_strand = np.where(p[single], "+", "-")
        d, s = nearest_exon(exons.get(chrom, np.empty((0, 2), np.int64)), start[single], end[single])
        # genomic side -> transcription-direction side
        gs = gene_strand
        side = np.where(d == 0, "overlap",
               np.where(((s == -1) & (gs == "+")) | ((s == 1) & (gs == "-")),
                        "upstream", "downstream"))
        out.append(
            pd.DataFrame(
                {
                    "chrom": chrom,
                    "start": start[single],
                    "end": end[single],
                    "strand": g.strand.values[single],
                    "milli_div": g.milli_div.values[single],
                    "gene_strand": gs,
                    "dist_to_exon": d,
                    "side": side,
                }
            )
        )
    cols = ["chrom", "start", "end", "strand", "milli_div", "gene_strand", "dist_to_exon", "side"]
    df = pd.concat(out, ignore_index=True) if out else pd.DataFrame({c: [] for c in cols})
    df["sense"] = df.strand.values == df.gene_strand.values
    return df


def zone_of(df: pd.DataFrame) -> pd.Series:
    z = pd.Series("intronic_far", index=df.index, dtype=object)
    z[df.side == "overlap"] = "overlaps_exon"
    near = (df.side != "overlap") & (df.dist_to_exon <= WINDOW_BP)
    z[near & (df.side == "upstream")] = f"within_{WINDOW_BP}_upstream"
    z[near & (df.side == "downstream")] = f"within_{WINDOW_BP}_downstream"
    return z


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    d = df.assign(zone=zone_of(df))
    rows = []
    for z, g in d.groupby("zone"):
        anti = int((~g.sense).sum())
        rows.append(
            {
                "zone": z,
                "n": len(g),
                "antisense": anti,
                "antisense_frac": anti / len(g),
                "median_milli_div": float(g.milli_div.median()),
            }
        )
    return pd.DataFrame(rows).sort_values("zone").reset_index(drop=True)


def window_test(df: pd.DataFrame, side: str) -> dict:
    d = df.assign(zone=zone_of(df))
    near = d[d.zone == f"within_{WINDOW_BP}_{side}"]
    far = d[d.zone == "intronic_far"]
    t = [
        [int((~near.sense).sum()), int(near.sense.sum())],
        [int((~far.sense).sum()), int(far.sense.sum())],
    ]
    odds, p = fisher_exact(t)
    return {
        "side": side,
        "n_near": len(near),
        "antisense_frac_near": t[0][0] / len(near) if len(near) else np.nan,
        "n_far": len(far),
        "antisense_frac_far": t[1][0] / len(far) if len(far) else np.nan,
        "odds_ratio": float(odds),
        "p": float(p),
    }


def by_divergence(df: pd.DataFrame, side: str, bins: int = 4) -> pd.DataFrame:
    d = df.assign(zone=zone_of(df))
    d = d[d.zone.isin([f"within_{WINDOW_BP}_{side}", "intronic_far"])].copy()
    d["div_bin"] = pd.qcut(d.milli_div, bins, duplicates="drop")
    rows = []
    for b, g in d.groupby("div_bin", observed=True):
        r = window_test(g, side)
        r["div_bin"] = str(b)
        rows.append(r)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rmsk")
    ap.add_argument("gencode")
    ap.add_argument("--rep-name", default="AluY")
    ap.add_argument("--out-dir", default="exonprox_out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)

    rep = load_repeats(a.rmsk, a.rep_name)
    plus, minus, exons = load_annotation(a.gencode)
    df = build(rep, plus, minus, exons)
    df.to_csv(os.path.join(a.out_dir, f"{a.rep_name}_exon_proximity.tsv.gz"),
              sep="\t", index=False, compression="gzip")

    summary = summarise(df)
    tests = pd.DataFrame([window_test(df, "upstream"), window_test(df, "downstream")])
    strat = by_divergence(df, "upstream")
    summary.to_csv(os.path.join(a.out_dir, "zones.tsv"), sep="\t", index=False)
    tests.to_csv(os.path.join(a.out_dir, "tests.tsv"), sep="\t", index=False)
    strat.to_csv(os.path.join(a.out_dir, "by_divergence.tsv"), sep="\t", index=False)
    json.dump(
        {
            "rep_name": a.rep_name,
            "copies_total": int(len(rep)),
            "copies_in_single_strand_genes": int(len(df)),
            "overall_sense_fraction": round(float(df.sense.mean()), 4),
        },
        open(os.path.join(a.out_dir, "summary.json"), "w"),
        indent=1,
    )
    print(summary.to_string(index=False))
    print()
    print(tests.to_string(index=False))
    print()
    print(strat.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
