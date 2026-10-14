"""Exon avoidance of MEI insertions against a uniform null and a GC-matched, within-gene null.

Inputs: a directory of annotated MEI VCFs, a GENCODE GTF, the hg38 2bit file, and two JSON files
(hg38 chromosome sizes {chrom: bp} and assembly gaps {chrom: [[start, end], ...]}, UCSC coordinates).

Method (each unique site counts once; sites by single linkage, --merge-bp):
  * a site is exonic if it lies inside the union of GENCODE exons, genic if inside a gene body
    (an insertion after POS lies inside a feature when start <= POS < end);
  * UNIFORM null: probability that a position drawn uniformly from the non-gap genome is exonic, weighted
    by each family's chromosome distribution (exact interval arithmetic);
  * MATCHED null: --n-null random non-gap positions (seeded) with their local GC; the expected exonic count
    among a family's GENIC sites is the sum, over those sites, of P(exon | genic, GC bin of the site)
    from the null, in --n-bins quantile bins of GC. Done for +/-250 bp and +/-50 bp windows.

Usage:
    python exon_depletion_null.py VCF_DIR GTF HG38.2bit SIZES.json GAPS.json --out-dir OUT
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_gtf_check import inside, union_intervals  # noqa: E402
from gene_orientation_summary import FAMILIES, assign_sites, collapse_sites, load_calls, wilson  # noqa: E402
from gene_overlap_table import load_gene_model  # noqa: E402

WINDOWS = ((250, "gc250"), (50, "gc50"))


def local_gc(seq: str) -> float:
    """GC fraction of A/C/G/T bases; NaN when more than half the window is N."""
    s = seq.upper()
    at = s.count("A") + s.count("T")
    gc = s.count("G") + s.count("C")
    return gc / (at + gc) if at + gc >= 0.5 * len(s) and at + gc > 0 else float("nan")


def window_gcs(chrom_seq, pos1: int) -> list[float]:
    """Local GC for each window in WINDOWS around a 1-based position."""
    p0, length = pos1 - 1, len(chrom_seq)
    return [local_gc(str(chrom_seq[max(0, p0 - h):min(length, p0 + h)])) for h, _ in WINDOWS]


def in_gap(gap_iv: np.ndarray, p0: int) -> bool:
    if len(gap_iv) == 0:
        return False
    i = int(np.searchsorted(gap_iv[:, 0], p0, side="right")) - 1
    return i >= 0 and gap_iv[i, 0] <= p0 < gap_iv[i, 1]


def bin_edges(values: pd.Series, n_bins: int) -> np.ndarray:
    e = np.unique(np.quantile(values, np.linspace(0, 1, n_bins + 1)))
    return e[1:-1]                                     # interior cut points for np.digitize


def matched_expectation(null_genic: pd.DataFrame, site_gc: pd.Series, col: str, n_bins: int) -> tuple[float, float]:
    """Expected exonic count (and its variance) among genic sites, given their GC and the genic null."""
    g = null_genic.dropna(subset=[col])
    cuts = bin_edges(g[col], n_bins)
    p_bin = g.assign(b=np.digitize(g[col], cuts)).groupby("b").exon.mean()
    p = p_bin.reindex(np.digitize(site_gc.dropna(), cuts)).values
    return float(np.nansum(p)), float(np.nansum(p * (1 - p)))


def sample_null(tb, chroms, callable_bp, gaps, n, rng) -> pd.DataFrame:
    weights = np.array([callable_bp[c] for c in chroms], float)
    weights /= weights.sum()
    picks = rng.choice(len(chroms), size=n, p=weights)
    rows = []
    for k in range(n):
        c = chroms[picks[k]]
        seq = tb[c]
        p1 = int(rng.integers(300, len(seq) - 300))
        if in_gap(gaps[c], p1 - 1):
            continue
        rows.append((c, p1, *window_gcs(seq, p1)))
    return pd.DataFrame(rows, columns=["chrom", "pos"] + [name for _, name in WINDOWS])


def label_positions(df: pd.DataFrame, exons: dict, bodies: dict) -> pd.DataFrame:
    df = df.copy()
    df["exon"], df["genic"] = False, False
    for c, ix in df.groupby("chrom").groups.items():
        pos = df.loc[ix, "pos"].values
        if c in exons:
            df.loc[ix, "exon"] = inside(exons[c], pos)
        if c in bodies:
            df.loc[ix, "genic"] = inside(bodies[c], pos)
    return df


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("vcf_dir", "gtf", "twobit", "sizes", "gaps"):
        ap.add_argument(a)
    ap.add_argument("--out-dir", default="exon_null_out")
    ap.add_argument("--merge-bp", type=int, default=50)
    ap.add_argument("--n-null", type=int, default=200_000)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args(argv)
    import twobitreader

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tb = twobitreader.TwoBitFile(args.twobit)
    sizes, gap_json = json.load(open(args.sizes)), json.load(open(args.gaps))
    chroms = [c for c in sizes if c != "chrM" and c in tb.keys()]
    gaps = {c: np.array(gap_json.get(c, []), dtype=np.int64).reshape(-1, 2) for c in chroms}
    callable_bp = {c: sizes[c] - int(sum(e - s for s, e in gap_json.get(c, []))) for c in chroms}

    genes, exons = load_gene_model(args.gtf)
    bodies = {c: union_intervals(x) for c, x in genes.groupby("chrom")}
    sites = collapse_sites(assign_sites(load_calls(args.vcf_dir), args.merge_bp))
    sites = label_positions(sites[["site", "family", "chrom", "pos"]].astype({"pos": int}), exons, bodies)
    sites["chrom"] = sites.chrom.astype(str)

    # UNIFORM null: exact interval arithmetic per chromosome
    p_exon = {c: sum(e - s for s, e in exons.get(c, np.zeros((0, 2), int))) / callable_bp[c] for c in chroms}
    p_genic = {c: sum(e - s for s, e in bodies.get(c, np.zeros((0, 2), int))) / callable_bp[c] for c in chroms}
    rng = np.random.default_rng(args.seed)
    null = label_positions(sample_null(tb, chroms, callable_bp, gaps, args.n_null, rng), exons, bodies)
    null_genic = null[null.genic]
    for h, name in WINDOWS:
        sites[name] = [window_gcs(tb[r.chrom], int(r.pos))[[w for w, _ in WINDOWS].index(h)] for r in sites.itertuples()]

    rows = []
    for fam in (*FAMILIES, "ALL"):
        x = sites if fam == "ALL" else sites[sites.family == fam]
        cnt = x.chrom.value_counts()
        n = len(x)
        for label, key in (("exonic", "exon"), ("genic", "genic")):
            exp_frac = sum(cnt[c] * (p_exon if key == "exon" else p_genic)[c] for c in cnt.index if c in p_exon) / n
            k = int(x[key].sum())
            lo, hi = wilson(k, n)
            rows.append({"null": "uniform", "family": fam, "class": label, "sites": n, "observed": k,
                         "expected": exp_frac * n, "obs_over_exp": k / (exp_frac * n),
                         "ci_lo": lo / exp_frac, "ci_hi": hi / exp_frac})
        gx = x[x.genic]
        for h, name in WINDOWS:
            exp, var = matched_expectation(null_genic, gx[name], name, args.n_bins)
            k, ng = int(gx.exon.sum()), len(gx.dropna(subset=[name]))
            lo, hi = wilson(k, ng)
            rows.append({"null": f"gc_matched_within_genes_{h}bp", "family": fam, "class": "exonic|genic",
                         "sites": ng, "observed": k, "expected": exp, "obs_over_exp": k / exp,
                         "ci_lo": lo * ng / exp, "ci_hi": hi * ng / exp})
    table = pd.DataFrame(rows)
    table.to_csv(out / "exon_depletion_null.tsv", sep="\t", index=False)
    summary = {"null_positions": int(len(null)), "null_genic_fraction": float(null.genic.mean()),
               "null_exon_fraction": float(null.exon.mean()), "seed": args.seed, "n_bins": args.n_bins,
               "median_gc250_by_family": sites.groupby("family").gc250.median().to_dict(),
               "median_gc250_null_genic": float(null_genic.gc250.median())}
    (out / "exon_depletion_summary.json").write_text(json.dumps(summary, indent=2))
    print(table.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
