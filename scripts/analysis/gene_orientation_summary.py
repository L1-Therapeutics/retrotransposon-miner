"""Gene-annotation inspection and sense/antisense summary for annotated MEI VCFs.

Reads every ``*.vcf`` in a directory (snpEff-annotated retrotransposon-miner
output with GENE/GENEID/GENE_STRAND/CSQ/CSQ_TERMS) and reports:

1. Register check: calls whose GENE and GENEID lists differ in length, i.e.
   where zipping GENE against GENE_STRAND would assign the wrong strand.
2. Neighbour-gene leakage: intronic calls whose CSQ_TERMS also contain an
   upstream/downstream term, so GENEID includes genes the breakpoint does not
   overlap, and how many mixed-strand gene lists that explains.
3. Sense fraction of intronic insertions per unique site (single linkage,
   ``--merge-bp``, per chromosome and family), by family, gene strand and
   private-vs-shared status, with Wilson intervals and exact binomial tests
   against 0.5.
4. LINE-1 gene-strand asymmetry: Fisher test, leave-one-chromosome-out range
   and a within-chromosome permutation test.

A site counts once however many genomes carry it (sites shared between genomes
are not independent observations). Sites whose member calls disagree on ORIENT
are excluded from strand tests and counted in the funnel.

Usage:
    python gene_orientation_summary.py VCF_DIR --out-dir OUT
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
from scipy.stats import binomtest, fisher_exact

FAMILIES = ("ALU", "LINE1", "SVA")
UPDOWN = {"upstream_gene_variant", "downstream_gene_variant"}


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (centre - half, centre + half)


def split_list(value) -> list[str]:
    return [t for t in str(value).split(",") if t] if isinstance(value, str) and value else []


def load_calls(vcf_dir: str | Path) -> pd.DataFrame:
    from retro_miner.gene_annotation import read_vcf

    paths = sorted(glob.glob(os.path.join(str(vcf_dir), "*.vcf")))
    if not paths:
        raise FileNotFoundError(f"no *.vcf in {vcf_dir}")
    rows = []
    for path in paths:
        _, records = read_vcf(path)
        sample = os.path.basename(path).split("_")[0]
        for r in records:
            rows.append({"sample": sample, "chrom": r.chrom, "pos": r.pos, **r.info})
    df = pd.DataFrame(rows)
    for col in ("GENE", "GENEID", "GENE_STRAND", "CSQ", "CSQ_TERMS", "ORIENT", "MEIFAMILY"):
        if col not in df:
            df[col] = np.nan
    df["gene_l"] = df.GENE.map(split_list)
    df["geneid_l"] = df.GENEID.map(split_list)
    df["strand_l"] = df.GENE_STRAND.map(split_list)
    df["has_updown"] = df.CSQ_TERMS.map(lambda x: bool(set(split_list(x)) & UPDOWN))
    df["strands"] = df.strand_l.map(frozenset)
    return df


def assign_sites(calls: pd.DataFrame, merge_bp: int) -> pd.DataFrame:
    d = calls[calls.MEIFAMILY.isin(FAMILIES)].sort_values(["chrom", "MEIFAMILY", "pos"]).reset_index(drop=True)
    new = (d.chrom != d.chrom.shift()) | (d.MEIFAMILY != d.MEIFAMILY.shift()) | (d.pos.diff() > merge_bp)
    d["site"] = new.cumsum()
    return d


def collapse_sites(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for site, x in d.groupby("site"):
        first = x.iloc[0]
        orients = x.ORIENT.dropna().unique()
        rows.append({
            "site": site, "family": first.MEIFAMILY, "chrom": first.chrom, "pos": first.pos,
            "n_samples": x["sample"].nunique(),
            "orient": orients[0] if len(orients) == 1 else None,
            "csq": first.CSQ, "has_updown": first.has_updown, "strands": first.strands,
            "geneid": first.geneid_l[0] if first.geneid_l else None,
        })
    return pd.DataFrame(rows)


def inspect_annotation(d: pd.DataFrame) -> dict:
    n_g, n_id, n_s = d.gene_l.map(len), d.geneid_l.map(len), d.strand_l.map(len)
    intr = d[d.CSQ == "intron_variant"]
    mixed = intr[intr.strands.map(len) == 2]
    return {
        "calls": int(len(d)),
        "calls_without_gene": int((n_id == 0).sum()),
        "gene_vs_geneid_length_mismatch": int((n_g != n_id).sum()),
        "strand_vs_geneid_length_mismatch": int((n_s != n_id).sum()),
        "intronic_calls": int(len(intr)),
        "intronic_with_updown_term": int(intr.has_updown.sum()),
        "intronic_mixed_strand": int(len(mixed)),
        "intronic_mixed_strand_with_updown_term": int(mixed.has_updown.sum()),
        "intronic_mixed_strand_without_updown_term": int((~mixed.has_updown).sum()),
        "updown_only_calls_carrying_gene_ids": int(
            ((d.CSQ.isin(UPDOWN)) & (n_id > 0)).sum()),
    }


def analysed_sites(sites: pd.DataFrame) -> pd.DataFrame:
    s = sites[(sites.csq == "intron_variant") & sites.orient.notna() & (sites.strands.map(len) == 1)].copy()
    s["gene_strand"] = s.strands.map(lambda x: next(iter(x)))
    s["sense"] = s.orient == s.gene_strand
    s["sharing"] = np.where(s.n_samples == 1, "private", "shared")
    return s


def sense_table(a: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sets = {"A_single_strand_intronic": a, "B_no_neighbour_term": a[~a.has_updown]}
    for set_name, base in sets.items():
        for fam in FAMILIES:
            f = base[base.family == fam]
            for gs in ("all", "+", "-"):
                g = f if gs == "all" else f[f.gene_strand == gs]
                for sh in ("all", "private", "shared"):
                    x = g if sh == "all" else g[g.sharing == sh]
                    n, k = len(x), int(x.sense.sum())
                    lo, hi = wilson(k, n)
                    rows.append({"set": set_name, "family": fam, "gene_strand": gs, "sharing": sh,
                                 "n_sites": n, "n_sense": k,
                                 "sense_frac": k / n if n else float("nan"),
                                 "ci_lo": lo, "ci_hi": hi,
                                 "p_vs_0.5": binomtest(k, n, 0.5).pvalue if n else float("nan")})
    return pd.DataFrame(rows)


def strand_asymmetry(a: pd.DataFrame, family: str, n_perm: int, seed: int) -> dict:
    x = a[a.family == family]
    pos, neg = x[x.gene_strand == "+"], x[x.gene_strand == "-"]
    ka, na, kb, nb = int(pos.sense.sum()), len(pos), int(neg.sense.sum()), len(neg)
    out = {"family": family, "plus_n": na, "plus_sense": ka, "minus_n": nb, "minus_sense": kb,
           "fisher_p": float(fisher_exact([[ka, na - ka], [kb, nb - kb]])[1])}
    if na == 0 or nb == 0:
        return out
    out["diff_plus_minus"] = ka / na - kb / nb
    loo = []
    for c in sorted(x.chrom.unique()):
        y = x[x.chrom != c]
        yp, yn = y[y.gene_strand == "+"].sense, y[y.gene_strand == "-"].sense
        if len(yp) and len(yn):
            loo.append(float(yp.mean() - yn.mean()))
    out["loo_n_chrom"] = len(loo)
    if loo:                      # empty when every site is on one chromosome (e.g. a chr22-only callset)
        out["loo_min"], out["loo_max"] = min(loo), max(loo)
    rng = np.random.default_rng(seed)
    gs, sn, ch = x.gene_strand.values, x.sense.values.astype(float), x.chrom.values
    groups = [np.where(ch == c)[0] for c in np.unique(ch)]
    obs, hits = abs(out["diff_plus_minus"]), 0
    for _ in range(n_perm):
        g2 = gs.copy()
        for ix in groups:
            g2[ix] = rng.permutation(gs[ix])
        if abs(sn[g2 == "+"].mean() - sn[g2 == "-"].mean()) >= obs:
            hits += 1
    out["perm_p"], out["n_perm"] = (hits + 1) / (n_perm + 1), n_perm
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("--out-dir", default="gene_orientation_out")
    ap.add_argument("--merge-bp", type=int, default=50)
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    calls = assign_sites(load_calls(args.vcf_dir), args.merge_bp)
    sites = collapse_sites(calls)
    a = analysed_sites(sites)

    report = {
        "genomes": int(calls["sample"].nunique()), "merge_bp": args.merge_bp,
        "unique_sites": int(len(sites)),
        "sites_by_family": sites.family.value_counts().to_dict(),
        "sites_with_discordant_orient": int(sites.orient.isna().sum()),
        "annotation": inspect_annotation(calls),
        "analysed_sites_by_family": a.family.value_counts().to_dict(),
        "line1_strand_asymmetry": strand_asymmetry(a, "LINE1", args.n_perm, args.seed),
        "alu_strand_asymmetry": strand_asymmetry(a, "ALU", args.n_perm, args.seed),
    }
    sense_table(a).to_csv(out / "sense_fraction.tsv", sep="\t", index=False)
    (out / "gene_annotation_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
