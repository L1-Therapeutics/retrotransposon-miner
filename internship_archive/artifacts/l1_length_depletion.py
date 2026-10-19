"""Length-dependence of the intronic LINE1 sense depletion.

Motivation (Finding 22). Han, Szak & Boeke 2004 (Nature 429:268) show L1 blocks
transcriptional elongation, and Ustyugova, Lebedev & Sverdlov 2006 (Genetica
128:261) show the expression-reducing effect is confined to *nearly full-sized*
L1 insertions in introns. Both imply a length-dependent effect, so the
sense-orientation depletion that Task 1 measured for intronic LINE1 should be
stronger for near-full-length elements than for 5'-truncated ones.

This script tests that prediction on the sample MEI call sets that carry both
element length (MEI_SPAN) and gene context (GENE_STRAND, CSQ). It is an
exploratory follow-up, not a pre-registered test: n is small and the direction
was suggested by the very literature it is testing.

Design
------
* Unit of analysis: unique intronic insertion site (de-duplicated by CHROM/POS
  across genomes, because MEI_SPAN and the reference gene strand are properties
  of the site, not of the sample).
* Keep single-gene, unambiguous-orientation LINE1 calls with CSQ intron_variant.
* orientation: sense if ORIENT == GENE_STRAND, else antisense (same as Task 1).
* Prediction: sense fraction is lower for near-full-length LINE1
  (MEI_SPAN >= 5000 bp, vs the ~6.0 kb full-length consensus) than for shorter.
* Null for comparison: ALU (a family Han's elongation mechanism is not about and
  whose span barely varies), reported the same way.

Pure-Python statistics (no scipy available in this environment).
"""

from __future__ import annotations

import csv
import glob
import math
import os
from collections import defaultdict

FULL_LENGTH_BP = 5000  # "nearly full-sized" cutoff (consensus L1 ~6.0 kb)


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def parse_info(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def _int(info: dict[str, str], key: str) -> int:
    try:
        return int(info.get(key, "0"))
    except ValueError:
        return 0


def load_sites(paths: list[str]) -> dict[str, dict[tuple[str, int], dict]]:
    """family -> {(chrom, pos): record} for clean single-gene intronic calls."""
    sites: dict[str, dict[tuple[str, int], dict]] = {
        "LINE1": {},
        "ALU": {},
    }
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                info = parse_info(f[7])
                fam = info.get("MEIFAMILY")
                if fam not in sites:
                    continue
                genes = [g for g in (info.get("GENE", "")).split(",") if g not in ("", ".")]
                strands = [s for s in (info.get("GENE_STRAND", "")).split(",") if s in ("+", "-")]
                terms = (info.get("CSQ_TERMS") or info.get("CSQ") or "").split(",")
                if len(genes) != 1 or len(strands) != 1:
                    continue  # single-gene, unambiguous strand only
                if "intron_variant" not in terms:
                    continue
                orient = info.get("ORIENT")
                if orient not in ("+", "-"):
                    continue
                span = _int(info, "MEI_SPAN")
                if span <= 0:
                    continue
                rec = {
                    "chrom": f[0],
                    "pos": int(f[1]),
                    "family": fam,
                    "subfamily": info.get("MEISUBFAMILY", ""),
                    "orient": orient,
                    "gene_strand": strands[0],
                    "span": span,
                    "mei5p": _int(info, "MEI_5P"),
                    "sense": orient == strands[0],
                }
                sites[fam].setdefault((f[0], int(f[1])), rec)
    return sites


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #
def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (centre - half, centre + half)


def _hypergeom_pmf(a: int, r: int, n1: int, n: int) -> float:
    # P(X = a) for X ~ Hypergeometric(N=n, K=r, n=n1)
    return math.comb(r, a) * math.comb(n - r, n1 - a) / math.comb(n, n1)


def fisher_exact(a: int, b: int, c: int, d: int) -> dict[str, float]:
    """2x2 [[a,b],[c,d]] = [[sense_long, anti_long],[sense_short, anti_short]].

    Returns two-sided p and the one-sided p for the predicted direction
    (long has a *lower* sense fraction than short).
    """
    n1, n2 = a + b, c + d
    r, n = a + c, a + b + c + d
    lo, hi = max(0, n1 - (n - r)), min(n1, r)

    def pmf(x: int) -> float:
        return _hypergeom_pmf(x, r, n1, n)

    obs = pmf(a)
    two_sided = sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= obs + 1e-12)
    # predicted direction: sense_long lower than expected -> small a
    one_sided = sum(pmf(x) for x in range(lo, a + 1))
    return {"p_two_sided": min(1.0, two_sided), "p_less": min(1.0, one_sided)}


def cochran_armitage(counts: list[tuple[float, int, int]]) -> dict[str, float]:
    """counts = [(score, n_i, sense_i), ...]; tests monotone trend in sense fraction."""
    n = sum(c[1] for c in counts)
    r = sum(c[2] for c in counts)
    if n == 0 or r == 0 or r == n:
        return {"z": float("nan"), "p_two_sided": float("nan"), "p_less": float("nan")}
    p_bar = r / n
    t = sum(x * (ri - ni * p_bar) for x, ni, ri in counts)
    s_x2 = sum(ni * x * x for x, ni, _ in counts)
    s_x = sum(ni * x for x, ni, _ in counts)
    var = p_bar * (1 - p_bar) * (s_x2 - s_x * s_x / n)
    if var <= 0:
        return {"z": float("nan"), "p_two_sided": float("nan"), "p_less": float("nan")}
    z = t / math.sqrt(var)
    two = math.erfc(abs(z) / math.sqrt(2))
    # p_less: one-sided for a *decreasing* trend in sense with increasing score
    less = 0.5 * math.erfc(-z / math.sqrt(2))  # P(Z <= z)
    return {"z": z, "p_two_sided": two, "p_less": less}


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    paths = sorted(
        p
        for p in glob.glob(os.path.join(root, "*classifier_ge_0.997*snpeff*.vcf"))
        if " (" not in os.path.basename(p)  # drop Finder " 2"/" (1)" duplicates
    )
    print(f"VCFs used ({len(paths)}):")
    for p in paths:
        print("  ", os.path.basename(p))

    sites = load_sites(paths)
    rows: list[dict] = []

    for fam in ("LINE1", "ALU"):
        recs = list(sites[fam].values())
        recs.sort(key=lambda r: (r["chrom"], r["pos"]))
        n = len(recs)
        k = sum(1 for r in recs if r["sense"])
        lo, hi = wilson(k, n)
        print(f"\n=== {fam}: {n} unique intronic sites, sense {k} ({k/n:.3f}) "
              f"95% CI [{lo:.3f}, {hi:.3f}] ===")
        rows.append({"family": fam, "stratum": "all", "n": n, "sense": k,
                     "sense_frac": k / n, "ci_low": lo, "ci_high": hi})

        long_ = [r for r in recs if r["span"] >= FULL_LENGTH_BP]
        short = [r for r in recs if r["span"] < FULL_LENGTH_BP]
        a = sum(1 for r in long_ if r["sense"])
        b = len(long_) - a
        c = sum(1 for r in short if r["sense"])
        d = len(short) - c
        print(f"  near-full-length (>= {FULL_LENGTH_BP} bp): n={len(long_)} sense={a} "
              f"({a/len(long_) if long_ else float('nan'):.3f})")
        print(f"  shorter                     : n={len(short)} sense={c} "
              f"({c/len(short) if short else float('nan'):.3f})")
        if long_ and short:
            fe = fisher_exact(a, b, c, d)
            print(f"  Fisher 2x2 [[{a},{b}],[{c},{d}]]: "
                  f"p(two-sided)={fe['p_two_sided']:.4f}  "
                  f"p(long less sense)={fe['p_less']:.4f}")
            for name, sub in (("near_full_length", long_), ("shorter", short)):
                kk = sum(1 for r in sub if r["sense"])
                l, h = wilson(kk, len(sub))
                rows.append({"family": fam, "stratum": name, "n": len(sub), "sense": kk,
                             "sense_frac": kk / len(sub), "ci_low": l, "ci_high": h,
                             "p_fisher_two_sided": fe["p_two_sided"],
                             "p_fisher_less": fe["p_less"]})

        # ordered trend across span bins (equal-width, line1-relevant)
        bins = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000),
                (4000, 5000), (5000, 10**9)]
        counts = []
        for i, (lo_b, hi_b) in enumerate(bins):
            sub = [r for r in recs if lo_b <= r["span"] < hi_b]
            if sub:
                sn = sum(1 for r in sub if r["sense"])
                counts.append((float(i), len(sub), sn))
        if len(counts) >= 3:
            ca = cochran_armitage(counts)
            print(f"  Cochran-Armitage trend over {len(counts)} span bins: "
                  f"z={ca['z']:.3f}  p(two-sided)={ca['p_two_sided']:.4f}  "
                  f"p(decreasing)={ca['p_less']:.4f}")
            rows.append({"family": fam, "stratum": "trend", "n": len(recs),
                         "sense": k, "sense_frac": k / n,
                         "z_trend": ca["z"], "p_trend_two_sided": ca["p_two_sided"],
                         "p_trend_decreasing": ca["p_less"]})

        # per gene strand (LINE1 depletion is known to sit on + strand genes)
        for strand in ("+", "-"):
            sub = [r for r in recs if r["gene_strand"] == strand]
            if not sub:
                continue
            sn = sum(1 for r in sub if r["sense"])
            l, h = wilson(sn, len(sub))
            sub_long = [r for r in sub if r["span"] >= FULL_LENGTH_BP]
            sub_short = [r for r in sub if r["span"] < FULL_LENGTH_BP]
            fe2 = None
            if sub_long and sub_short:
                aa = sum(1 for r in sub_long if r["sense"])
                bb = len(sub_long) - aa
                cc = sum(1 for r in sub_short if r["sense"])
                dd = len(sub_short) - cc
                fe2 = fisher_exact(aa, bb, cc, dd)
            print(f"  gene strand {strand}: n={len(sub)} sense={sn} "
                  f"({sn/len(sub):.3f})  near-full n={len(sub_long)}"
                  + (f"  Fisher p_less={fe2['p_less']:.4f}" if fe2 else ""))
            rows.append({"family": fam, "stratum": f"gene_strand_{strand}", "n": len(sub),
                         "sense": sn, "sense_frac": sn / len(sub), "ci_low": l, "ci_high": h,
                         **({"p_fisher_less": fe2["p_less"]} if fe2 else {})})

    fields = ["family", "stratum", "n", "sense", "sense_frac", "ci_low", "ci_high",
              "p_fisher_two_sided", "p_fisher_less", "z_trend", "p_trend_two_sided",
              "p_trend_decreasing"]
    out = os.path.join(here, "l1_length_depletion.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
