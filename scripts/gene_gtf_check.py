"""Independent re-annotation of MEI sites against a GENCODE GTF, and orientation cross-check.

Companion to gene_orientation_summary.py (imports its helpers; keep both files in one folder).

Reports, for the unique sites of a directory of annotated MEI VCFs:
1. Database fingerprint: fraction of the callset's Ensembl gene IDs present in the GTF.
2. Class concordance: GTF interval overlap (gene body, exon union) against snpEff CSQ.
3. Strand resolution: intronic sites whose snpEff gene list is mixed-strand but whose
   gene-body overlap is single-strand, and strand agreement where both are single-strand.
4. Sense fraction on the GTF-resolved set, by family and gene strand.
5. Orientation cross-check: ORIENT (split-read strand) against the poly(A) side
   (polyA_side:L = left-clipped = tail on the right = '+'), and the LINE-1 gene-strand
   asymmetry under each orientation measure.

Usage:
    python gene_gtf_check.py VCF_DIR GENCODE.annotation.gtf.gz --out-dir OUT
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_orientation_summary import FAMILIES, assign_sites, collapse_sites, load_calls, wilson  # noqa: E402


def load_gtf(path: str) -> pd.DataFrame:
    cols = ["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"]
    g = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols, dtype={"chrom": str, "attr": str},
                    usecols=[0, 2, 3, 4, 6, 8], compression="infer")
    g = g[g.feat.isin(("gene", "exon"))].copy()
    g["gene_id"] = g.attr.str.extract(r'gene_id "([^"]+)"')[0].str.split(".").str[0]
    return g.drop(columns="attr")


def union_intervals(x: pd.DataFrame) -> np.ndarray:
    out: list[list[int]] = []
    for s, e in x.sort_values("start")[["start", "end"]].values:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return np.array(out)


def inside(iv: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """True where an insertion after `pos` lies inside a feature (feature covers pos and pos+1)."""
    i = np.searchsorted(iv[:, 0], pos, side="right") - 1
    ok = i >= 0
    i = np.clip(i, 0, len(iv) - 1)
    return ok & (iv[i, 0] <= pos) & (pos < iv[i, 1])


def overlay_gtf(sites: pd.DataFrame, gtf: pd.DataFrame) -> pd.DataFrame:
    genes = gtf[gtf.feat == "gene"]
    exons = {c: union_intervals(x) for c, x in gtf[gtf.feat == "exon"].groupby("chrom")}
    bodies = {c: x for c, x in genes.groupby("chrom")}
    rows = []
    for chrom, x in sites.groupby("chrom"):
        pos = x.pos.values
        body = bodies.get(chrom)
        if body is None:
            rows += [(s, 0, frozenset(), False) for s in x.site]
            continue
        hit = (body.start.values[None, :] <= pos[:, None]) & (pos[:, None] < body.end.values[None, :])
        in_exon = inside(exons[chrom], pos) if chrom in exons else np.zeros(len(pos), bool)
        for k, s in enumerate(x.site.values):
            idx = np.where(hit[k])[0]
            rows.append((s, len(idx), frozenset(body.strand.values[idx]), bool(in_exon[k])))
    out = sites.merge(pd.DataFrame(rows, columns=["site", "n_gene", "g_strands", "g_exon"]), on="site")
    out["gtf_class"] = np.where(out.n_gene == 0, "intergenic", np.where(out.g_exon, "exonic", "intronic"))
    return out


def sense_summary(z: pd.DataFrame, col: str) -> dict:
    out = {}
    for fam in FAMILIES:
        x = z[(z.family == fam) & z[col].notna()]
        sense = x[col] == x.gene_strand
        k, n = int(sense.sum()), len(x)
        a, b = x[x.gene_strand == "+"], x[x.gene_strand == "-"]
        ka, kb = int((a[col] == a.gene_strand).sum()), int((b[col] == b.gene_strand).sum())
        lo, hi = wilson(k, n)
        out[fam] = {"n": n, "n_sense": k, "sense_frac": k / n if n else None, "ci": [lo, hi],
                    "p_vs_0.5": binomtest(k, n, 0.5).pvalue if n else None,
                    "plus_genes": [ka, len(a)], "minus_genes": [kb, len(b)],
                    "asymmetry_fisher_p": float(fisher_exact([[ka, len(a) - ka], [kb, len(b) - kb]])[1])
                    if len(a) and len(b) else None}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("gtf")
    ap.add_argument("--out-dir", default="gene_gtf_out")
    ap.add_argument("--merge-bp", type=int, default=50)
    args = ap.parse_args(argv)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    calls = assign_sites(load_calls(args.vcf_dir), args.merge_bp)
    calls["polyA_side"] = calls.SUPPORT.fillna("").str.extract(r"polyA_side:(\w+)")[0]
    gtf = load_gtf(args.gtf)
    report: dict = {"gtf": Path(args.gtf).name}

    prod = {g for ids in calls.geneid_l for g in ids}
    ids = set(gtf[gtf.feat == "gene"].gene_id)
    report["gene_ids_in_callset"] = len(prod)
    report["gene_ids_found_in_gtf"] = len(prod & ids)

    sites = collapse_sites(calls)
    s2 = overlay_gtf(sites, gtf)
    snp = s2.csq.map(lambda c: "intron" if c == "intron_variant" else "intergenic" if c == "intergenic_variant"
                     else "up/downstream" if "stream" in c else "exon/UTR")
    report["class_concordance"] = pd.crosstab(s2.gtf_class, snp).to_dict()

    intr = s2[s2.gtf_class == "intronic"].copy()
    intr["snp_single"] = intr.strands.map(len) == 1
    intr["gtf_single"] = intr.g_strands.map(len) == 1
    both = intr[intr.snp_single & intr.gtf_single]
    report["strand_resolution"] = {
        "intronic_sites": len(intr),
        "both_single_strand": len(both),
        "strand_agrees": int((both.strands.map(lambda s: next(iter(s))) == both.g_strands.map(lambda s: next(iter(s)))).sum()),
        "snpeff_mixed_but_gtf_single": int(((~intr.snp_single) & intr.gtf_single).sum()),
        "snpeff_single_but_gtf_mixed": int((intr.snp_single & ~intr.gtf_single).sum()),
    }

    z = intr[intr.gtf_single & intr.orient.notna()].copy()
    z["gene_strand"] = z.g_strands.map(lambda s: next(iter(s)))
    site_pa = (calls.groupby("site").polyA_side.apply(lambda s: s.map({"L": "+", "R": "-"}).dropna().unique())
               .map(lambda u: u[0] if len(u) == 1 else None).rename("pa_orient"))
    z = z.join(site_pa, on="site")
    report["orientation_agreement_calls"] = {
        fam: {"agreement": float((x.ORIENT == np.where(x.polyA_side == "L", "+", "-")).mean()), "n": len(x)}
        for fam, x in calls[calls.polyA_side.isin(["L", "R"]) & calls.ORIENT.isin(["+", "-"])].groupby("MEIFAMILY")}
    report["analysed_sites"] = {"n": len(z), **z.family.value_counts().to_dict()}
    report["sense_by_ORIENT"] = sense_summary(z, "orient")
    report["sense_by_polyA"] = sense_summary(z, "pa_orient")
    agree = z[z.orient == z.pa_orient]
    report["sense_by_ORIENT_where_measures_agree"] = sense_summary(agree, "orient")

    (out / "gene_gtf_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
