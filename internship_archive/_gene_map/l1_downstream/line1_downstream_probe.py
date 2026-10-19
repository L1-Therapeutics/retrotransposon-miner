"""Is the antisense-LINE-1 enrichment downstream of exons a gene-end effect?

Found while testing something else: pooled L1PA3-L1PA7 copies are ENRICHED for
antisense within 150 bp downstream of an exon (0.683 vs 0.594, OR 1.47,
p=0.0004, n=391) -- opposite in direction and position to the Alu result.

The obvious confounder is gene position. Lavi & Carmel 2018 report that Alu
orientation bias REVERSES at gene ends (77.7% sense there against 54.9%
antisense in introns generally), so a signal concentrated near the 3' end of
genes would be that known effect, not a new one. This stratifies the same test
by where in the gene the element sits and by which exon it is near.

Usage:
  python line1_downstream_probe.py <rmsk.txt.gz> <gencode.gtf.gz> --out-dir out
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
from collections import defaultdict

import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# alu_exon_proximity.py must sit beside this file. Add our own directory rather
# than relying on the caller's cwd, so the script runs from anywhere.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from alu_exon_proximity import load_repeats, merge_intervals, nearest_exon  # noqa: E402

WINDOW_BP = 150


def load_genes_and_exons(path: str):
    """Per-chromosome gene records with their own exon lists.

    genes[chrom] -> list of (start, end, strand, gene_id)
    gene_exons[gene_id] -> merged np.ndarray of that gene's exons
    """
    genes = defaultdict(list)
    gex = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            chrom, s, e, strand, attr = f[0], int(f[3]) - 1, int(f[4]), f[6], f[8]
            gid = attr.split('gene_id "', 1)[1].split('"', 1)[0] if 'gene_id "' in attr else None
            if gid is None:
                continue
            if f[2] == "gene":
                if strand in ("+", "-"):
                    genes[chrom].append((s, e, strand, gid))
            else:
                gex[gid].append((s, e))
    gene_exons = {g: merge_intervals(np.array(sorted(v), dtype=np.int64)) for g, v in gex.items()}
    return dict(genes), gene_exons


def containing_gene(gene_list, pos: int):
    """The single gene covering pos, or None when zero or several strands cover it."""
    hits = [g for g in gene_list if g[0] <= pos < g[1]]
    if not hits:
        return None
    if len({g[2] for g in hits}) != 1:
        return None
    # when several genes on one strand overlap, take the longest (the host locus)
    return max(hits, key=lambda g: g[1] - g[0])


def exon_rank(exons: np.ndarray, idx: int, strand: str) -> str:
    """Is exon `idx` the first, last, or an internal exon in transcription order."""
    n = len(exons)
    if n == 1:
        return "only"
    first, last = (0, n - 1) if strand == "+" else (n - 1, 0)
    if idx == first:
        return "first"
    if idx == last:
        return "last"
    return "internal"


def nearest_exon_index(exons: np.ndarray, start: int, end: int) -> int:
    best_i, best_d = -1, None
    for i, (s, e) in enumerate(exons):
        d = max(s - end, start - e)
        d = 0 if d < 0 else int(d)
        if best_d is None or d < best_d:
            best_d, best_i = d, i
    return best_i


def build(repeats: pd.DataFrame, genes: dict, gene_exons: dict) -> pd.DataFrame:
    out = []
    for r in repeats.itertuples(index=False):
        g = containing_gene(genes.get(r.chrom, []), (r.start + r.end) // 2)
        if g is None:
            continue
        gs, ge, gstrand, gid = g
        ex = gene_exons.get(gid)
        if ex is None or len(ex) == 0:
            continue
        d, s = nearest_exon(ex, np.array([r.start]), np.array([r.end]))
        dist, gside = int(d[0]), int(s[0])
        if dist == 0:
            continue
        side = ("upstream" if ((gside == -1) == (gstrand == "+")) else "downstream")
        i = nearest_exon_index(ex, r.start, r.end)
        # position along the gene in transcription direction, 0 = TSS, 1 = TES
        mid = (r.start + r.end) / 2
        frac = (mid - gs) / max(ge - gs, 1)
        if gstrand == "-":
            frac = 1 - frac
        out.append(
            {
                "chrom": r.chrom, "start": r.start, "end": r.end,
                "strand": r.strand, "milli_div": r.milli_div,
                "gene_id": gid, "gene_strand": gstrand, "gene_len": ge - gs,
                "dist_to_exon": dist, "side": side,
                "exon_rank": exon_rank(ex, i, gstrand), "n_exons": len(ex),
                "gene_frac": round(float(frac), 4),
                "dist_to_tes": int(ge - mid if gstrand == "+" else mid - gs),
                "sense": r.strand == gstrand,
            }
        )
    return pd.DataFrame(out)


def test(df: pd.DataFrame, side: str) -> dict:
    near = df[(df.side == side) & (df.dist_to_exon <= WINDOW_BP)]
    far = df[df.dist_to_exon > WINDOW_BP]
    t = [[int((~near.sense).sum()), int(near.sense.sum())],
         [int((~far.sense).sum()), int(far.sense.sum())]]
    odds, p = fisher_exact(t)
    return {"side": side, "n_near": len(near),
            "anti_near": round(t[0][0] / len(near), 4) if len(near) else np.nan,
            "n_far": len(far),
            "anti_far": round(t[1][0] / len(far), 4) if len(far) else np.nan,
            "OR": round(float(odds), 3), "p": float(p)}


def stratified(df: pd.DataFrame, side: str, by: str) -> pd.DataFrame:
    rows = []
    for k, g in df.groupby(by, observed=True):
        r = test(g, side)
        r[by] = str(k)
        rows.append(r)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rmsk")
    ap.add_argument("gencode")
    ap.add_argument("--rep-name", default="L1PA3,L1PA4,L1PA5,L1PA6,L1PA7")
    ap.add_argument("--out-dir", default="l1probe_out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)

    rep = load_repeats(a.rmsk, a.rep_name)
    genes, gene_exons = load_genes_and_exons(a.gencode)
    df = build(rep, genes, gene_exons)
    df.to_csv(os.path.join(a.out_dir, "copies.tsv.gz"), sep="\t", index=False, compression="gzip")

    overall = pd.DataFrame([test(df, "downstream"), test(df, "upstream")])
    by_rank = stratified(df, "downstream", "exon_rank")
    df = df.copy()
    df["gene_third"] = pd.cut(df.gene_frac, [-0.01, 1 / 3, 2 / 3, 1.01],
                              labels=["first_third", "middle", "last_third"])
    by_third = stratified(df, "downstream", "gene_third")
    df["tes_bin"] = pd.cut(df.dist_to_tes, [-1, 5000, 20000, 50000, 10 ** 9],
                           labels=["<=5kb", "5-20kb", "20-50kb", ">50kb"])
    by_tes = stratified(df, "downstream", "tes_bin")

    for name, t in [("overall", overall), ("by_exon_rank", by_rank),
                    ("by_gene_third", by_third), ("by_tes_distance", by_tes)]:
        t.to_csv(os.path.join(a.out_dir, f"{name}.tsv"), sep="\t", index=False)
        print(f"== {name} ==")
        print(t.to_string(index=False))
        print()
    json.dump({"rep_name": a.rep_name, "n_copies": int(len(df))},
              open(os.path.join(a.out_dir, "summary.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
