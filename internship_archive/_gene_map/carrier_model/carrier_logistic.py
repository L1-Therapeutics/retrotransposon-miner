"""Logistic model for the poly(A) carrier asymmetry, replacing 20-stratum Mantel-Haenszel.

Two weaknesses of the stratified estimate this addresses:

1. Mantel-Haenszel bins divergence and length into 20 cells, so within-cell
   residual confounding is invisible. Here both enter as continuous terms, with
   a quadratic on divergence because its relationship to orientation need not be
   linear over a 0-290 milliDiv range.
2. The cohort was "copies inside single-strand genes", which includes copies
   overlapping an exon. Exonic Alu is a different population -- some of it is
   exonized Alu, the very thing selection is acting on -- so the model is also
   fitted on intronic copies only.

A pooled fit with a carrier x subfamily interaction tests directly whether the
effect really differs between subfamilies, rather than eyeballing three
separate confidence intervals.

Usage:
  python carrier_logistic.py <gencode.gtf.gz> <AluY.csv.gz> <AluSx.csv.gz> <AluJo.csv.gz> --out-dir out
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.api as sm

EXON_RE = re.compile(r"\t(exon)\t")


def load_exons(path: str) -> dict:
    """Merged exon intervals per chromosome, 0-based half-open."""
    rows = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            rows[f[0]].append((int(f[3]) - 1, int(f[4])))
    out = {}
    for c, v in rows.items():
        a = np.array(sorted(v), dtype=np.int64)
        merged = [a[0].copy()]
        for s, e in a[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append(np.array([s, e]))
        out[c] = np.array(merged, dtype=np.int64)
    return out


def flag_exonic(df: pd.DataFrame, exons: dict) -> np.ndarray:
    """True when the copy overlaps any exon at all."""
    hit = np.zeros(len(df), dtype=bool)
    for chrom, g in df.groupby("chrom", sort=False):
        iv = exons.get(chrom)
        if iv is None or len(iv) == 0:
            continue
        starts, ends = iv[:, 0], iv[:, 1]
        # side="left" so an exon starting exactly where the copy ends is NOT a hit:
        # intervals are half-open, so [900,1000) and [1000,1200) do not overlap.
        i = np.searchsorted(starts, g.end.values, side="left") - 1
        ok = i >= 0
        overlap = np.zeros(len(g), dtype=bool)
        overlap[ok] = ends[i[ok]] > g.start.values[ok]
        hit[df.index.get_indexer(g.index)] = overlap
    return hit


def prepare(df: pd.DataFrame, rep: str, exonic: np.ndarray) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "sense": df.sense.astype(int).values,
            "carrier_any": (df.n_hex_any > 0).astype(int).values,
            "carrier_body": (df.n_hex_body > 0).astype(int).values,
            "div10": df.milli_div.values / 10.0,
            "len10": df.length_bp.values / 10.0,
            "exonic": exonic.astype(int),
            "rep": rep,
        }
    )
    out["div10_sq"] = out.div10 ** 2
    return out


def fit(d: pd.DataFrame, carrier_col: str, extra: list[str] | None = None):
    """Logistic regression of sense on carrier status plus covariates."""
    cols = [carrier_col, "div10", "div10_sq", "len10"] + (extra or [])
    X = sm.add_constant(d[cols].astype(float), has_constant="add")
    model = sm.GLM(d["sense"].values, X, family=sm.families.Binomial())
    return model.fit()


def summarise(res, term: str, label: str, n: int) -> dict:
    ci = res.conf_int()
    lo, hi = ci.loc[term]
    return {
        "model": label,
        "n": int(n),
        "term": term,
        "OR": round(float(np.exp(res.params[term])), 4),
        "ci_lo": round(float(np.exp(lo)), 4),
        "ci_hi": round(float(np.exp(hi)), 4),
        "p": float(res.pvalues[term]),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("gencode")
    ap.add_argument("per_copy", nargs="+", help="alu_per_copy.csv.gz files, one per subfamily")
    ap.add_argument("--names", default="AluY,AluSx,AluJo")
    ap.add_argument("--out-dir", default="logit_out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    names = [x.strip() for x in a.names.split(",")]
    if len(names) != len(a.per_copy):
        raise SystemExit("--names must have one entry per per_copy file")

    exons = load_exons(a.gencode)
    frames = []
    for rep, path in zip(names, a.per_copy):
        df = pd.read_csv(path).reset_index(drop=True)
        frames.append(prepare(df, rep, flag_exonic(df, exons)))
    allc = pd.concat(frames, ignore_index=True)
    allc.to_csv(os.path.join(a.out_dir, "model_input.csv.gz"), index=False, compression="gzip")

    rows = []
    for rep in names:
        d = allc[allc.rep == rep]
        di = d[d.exonic == 0]
        for carrier in ("carrier_any", "carrier_body"):
            rows.append(summarise(fit(d, carrier), carrier, f"{rep} all copies", len(d)))
            rows.append(summarise(fit(di, carrier), carrier, f"{rep} intronic only", len(di)))
    per_rep = pd.DataFrame(rows)

    # pooled fit with a carrier x subfamily interaction
    pooled_rows = []
    for scope, d in (("all copies", allc), ("intronic only", allc[allc.exonic == 0])):
        d = d.copy()
        for rep in names[1:]:
            d[f"is_{rep}"] = (d.rep == rep).astype(int)
            d[f"carrier_x_{rep}"] = d.carrier_any * d[f"is_{rep}"]
        extra = [f"is_{r}" for r in names[1:]] + [f"carrier_x_{r}" for r in names[1:]]
        res = fit(d, "carrier_any", extra)
        for term in ["carrier_any"] + [f"carrier_x_{r}" for r in names[1:]]:
            pooled_rows.append(summarise(res, term, f"pooled, {scope}", len(d)))
    pooled = pd.DataFrame(pooled_rows)

    exonic_share = allc.groupby("rep").exonic.mean().round(4).to_dict()
    per_rep.to_csv(os.path.join(a.out_dir, "per_subfamily.tsv"), sep="\t", index=False)
    pooled.to_csv(os.path.join(a.out_dir, "pooled_interaction.tsv"), sep="\t", index=False)
    json.dump({"exonic_share_by_subfamily": exonic_share,
               "n_total": int(len(allc)),
               "n_intronic": int((allc.exonic == 0).sum())},
              open(os.path.join(a.out_dir, "summary.json"), "w"), indent=1)
    print(per_rep.to_string(index=False))
    print()
    print(pooled.to_string(index=False))
    print()
    print("exonic share by subfamily:", exonic_share)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
