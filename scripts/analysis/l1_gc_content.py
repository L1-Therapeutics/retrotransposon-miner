"""GC / G-quadruplex composition at MEI insertion sites, genome-wide.

Question this answers
---------------------
Do mobile-element insertions land in GC-poor (AT-rich) DNA, and does LINE-1
differ from Alu/SVA? This re-runs Finding 17 (artifacts/l1_literature_tracker.md),
which tested the "salt -> guanine -> G4" hypothesis on 30 chr22 calls and was
left UNRESOLVED for LINE-1: only 7 of those 30 were L1, and the tracker says
"the L1-specific question is UNRESOLVED, not answered. Re-run on the
genome-wide callset." This script is that re-run.

Method
------
  * Calls are read from a directory of retrotransposon-miner VCFs, collapsed to
    unique sites by single linkage (per chrom x family, --merge-bp) so a site
    carried by several genomes counts once.
  * For every site the reference window of +/- --window bp (default 150, i.e.
    301 bp, matching Finding 17) is read from the hg38 2bit and scored:
    G / C / GC / AT fraction over A/C/G/T, and max |G4Hunter| over 25-nt windows
    on both strands.
  * Background is drawn uniformly from the non-gap genome (chromosome-weighted
    by non-gap length), keeping only N-free windows with the same A/C/G/T
    requirement, so the null has the same length and quality filter.
  * G4Hunter follows the implementation used in Finding 17:
    per-base score = signed homopolymer G/C run length capped at 4, A/T = 0;
    window score = mean over 25 nt; reported as the max absolute score across
    both strands and both window positions.
  * The target-site duplication (TSD) carried in each call is summarised as an
    internal positive control: the L1 endonuclease nicks an AT-rich motif, so
    L1 TSDs should be AT-rich irrespective of the reference background.

Usage:
    python l1_gc_content.py VCF_DIR HG38.2bit SIZES.json GAPS.json --out-dir OUT
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_orientation_summary import assign_sites, collapse_sites, load_calls  # noqa: E402

G4_WINDOW = 25
G4_RUN_CAP = 4
G4_THRESHOLD = 1.2
FEATURES = ("g_frac", "c_frac", "gc_frac", "at_frac", "max_abs_g4")

# LINE-1 age classes, by subfamily prefix (RepeatMasker naming).
L1_AGE = (
    (("L1HS", "L1PA"), "primate/young"),
    (("L1PB", "L1MC", "L1MB"), "intermediate"),
    (("L1M",), "old"),
)


def g4hunter_base_scores(seq: str) -> np.ndarray:
    """Finding 17's G4Hunter: signed homopolymer G/C run length, capped at 4."""
    out = np.zeros(len(seq))
    for m in re.finditer(r"G+|C+", seq):
        run = min(len(m.group()), G4_RUN_CAP)
        out[m.start():m.end()] = run if m.group()[0] == "G" else -run
    return out


def max_abs_g4hunter(seq: str) -> float:
    """Max |G4Hunter| over 25-nt windows; both strands are covered because the
    reverse complement swaps G-runs for C-runs, which the base scores give as
    the negated signal already."""
    b = g4hunter_base_scores(seq)
    if len(b) < G4_WINDOW:
        return float("nan")
    k = np.convolve(b, np.ones(G4_WINDOW) / G4_WINDOW, "valid")
    return float(np.max(np.abs(k))) if len(k) else float("nan")


def window_features(seq: str) -> dict | None:
    """Composition of one reference window; None when it is mostly N."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return None
    a = s.count("A")
    c = s.count("C")
    g = s.count("G")
    t = s.count("T")
    acgt = a + c + g + t
    if acgt < 0.5 * n or acgt == 0:          # same N rule as exon_depletion_null.local_gc
        return None
    return {
        "g_frac": g / acgt,
        "c_frac": c / acgt,
        "gc_frac": (g + c) / acgt,
        "at_frac": (a + t) / acgt,
        "max_abs_g4": max_abs_g4hunter(s),
    }


def fetch(tb, chrom: str, pos1: int, half: int) -> str | None:
    """Reference window [pos-half, pos+half] (0-based half-open, 2*half+1 bp)."""
    if chrom not in tb:
        return None
    length = len(tb[chrom])
    start = pos1 - 1 - half
    end = pos1 + half
    if start < 0 or end > length:
        return None
    return str(tb[chrom][start:end]).upper()


def non_gap_lengths(sizes: dict, gaps: dict) -> tuple[dict, float]:
    out = {}
    for chrom, bp in sizes.items():
        if chrom not in sizes:
            continue
        covered = sum(e - s for s, e in gaps.get(chrom, []))
        out[chrom] = max(0, bp - covered)
    return out, float(sum(out.values()))


def in_gap(gaps_chrom: list, p0: int) -> bool:
    for s, e in gaps_chrom:
        if s <= p0 < e:
            return True
    return False


def sample_background(tb, sizes, gaps, half, n_bg, seed):
    """Uniform non-gap, N-free reference windows, chromosome-weighted."""
    rng = np.random.default_rng(seed)
    weights, _ = non_gap_lengths(sizes, gaps)
    chroms = [c for c in sizes if c in tb and weights[c] > 0]
    w = np.array([weights[c] for c in chroms], dtype=float)
    w /= w.sum()
    feats, tries = [], 0
    max_tries = n_bg * 40
    while len(feats) < n_bg and tries < max_tries:
        tries += 1
        chrom = chroms[int(rng.choice(len(chroms), p=w))]
        p0 = int(rng.integers(0, sizes[chrom]))
        if in_gap(gaps.get(chrom, []), p0):
            continue
        seq = fetch(tb, chrom, p0 + 1, half)
        if seq is None:
            continue
        f = window_features(seq)
        if f is not None:
            feats.append({"chrom": chrom, **f})
    return pd.DataFrame(feats)


def parse_subfamily(value) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    token = value.split("#", 1)[0]
    return re.sub(r"_(5|3)end$", "", token) or None


def l1_age_class(sub: str | None) -> str:
    if not sub:
        return "unknown"
    for prefixes, label in L1_AGE:
        if sub.startswith(prefixes):
            return label
    return "other"


def tsd_composition(df: pd.DataFrame) -> pd.DataFrame:
    """AT/GC fraction of the target-site duplication, by family (positive control)."""
    rows = []
    for fam, x in df.groupby("MEIFAMILY"):
        vals = []
        for tsd in x.TSD.dropna():
            s = str(tsd).upper()
            acgt = sum(s.count(b) for b in "ACGT")
            if acgt:
                vals.append((s.count("A") + s.count("T")) / acgt)
        if vals:
            rows.append({"family": fam, "n_tsd": len(vals),
                         "tsd_at_frac_median": float(np.median(vals))})
    return pd.DataFrame(rows)


def compare(site_feats: pd.DataFrame, bg: pd.DataFrame, label: str) -> pd.DataFrame:
    from scipy.stats import fisher_exact, mannwhitneyu

    rows = []
    n_site, n_bg = len(site_feats), len(bg)
    for feat in FEATURES:
        a = site_feats[feat].dropna().values
        b = bg[feat].dropna().values
        if len(a) == 0 or len(b) == 0:
            continue
        two = mannwhitneyu(a, b, alternative="two-sided")
        less = mannwhitneyu(a, b, alternative="less")
        row = {
            "set": label, "feature": feat, "n_sites": len(a), "n_bg": len(b),
            "site_median": float(np.median(a)), "bg_median": float(np.median(b)),
            "delta": float(np.median(a) - np.median(b)),
            "p_two_sided": float(two.pvalue), "p_site_less": float(less.pvalue),
        }
        if feat == "max_abs_g4":
            sa = int((a >= G4_THRESHOLD).sum())
            sb = int((b >= G4_THRESHOLD).sum())
            row["site_frac_ge_thr"] = sa / len(a)
            row["bg_frac_ge_thr"] = sb / len(b)
            row["fisher_or"] = float(fisher_exact(
                [[sa, len(a) - sa], [sb, len(b) - sb]])[0])
            row["fisher_p"] = float(fisher_exact(
                [[sa, len(a) - sa], [sb, len(b) - sb]])[1])
        rows.append(row)
    return pd.DataFrame(rows)


def make_figure(site_feats: pd.DataFrame, bg: pd.DataFrame, out: Path) -> None:
    """GC-fraction and max|G4Hunter| distributions: each family vs background."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                                  # pragma: no cover
        print(f"figure skipped: {exc}")
        return
    fams = [f for f in ("LINE1", "ALU", "SVA") if (site_feats.family == f).any()]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    colors = {"LINE1": "tab:red", "ALU": "tab:blue", "SVA": "tab:green"}
    for ax, feat, title in zip(axes, ("gc_frac", "max_abs_g4"),
                               ("GC fraction", "max |G4Hunter| (25 nt)")):
        data = [bg[feat].dropna().values] + [site_feats[site_feats.family == f][feat].dropna().values
                                             for f in fams]
        labels = ["background"] + fams
        bp = ax.boxplot(data, showfliers=False, patch_artist=True)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
        ax.axhline(float(bg[feat].median()), color="0.5", ls="--", lw=1)
        for patch, label in zip(bp["boxes"], labels):
            patch.set_facecolor("0.85" if label == "background" else colors[label])
        ax.set_title(title)
    fig.suptitle("MEI insertion-site composition vs genome background (hg38, 301 bp)")
    fig.tight_layout()
    fig.savefig(out / "l1_gc_content.png", dpi=130)
    print(f"wrote {out / 'l1_gc_content.png'}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("twobit")
    ap.add_argument("sizes_json")
    ap.add_argument("gaps_json")
    ap.add_argument("--out-dir", default="l1_gc_content_out")
    ap.add_argument("--window", type=int, default=150, help="half-window bp (301 bp default)")
    ap.add_argument("--merge-bp", type=int, default=50)
    ap.add_argument("--n-bg", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--families", default="LINE1,ALU,SVA")
    ap.add_argument("--figure", action="store_true", help="save a GC/G4 boxplot PNG")
    args = ap.parse_args(argv)

    import twobitreader

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sizes = json.loads(Path(args.sizes_json).read_text())
    gaps = json.loads(Path(args.gaps_json).read_text())
    families = [f.strip() for f in args.families.split(",") if f.strip()]

    calls = load_calls(args.vcf_dir)
    d = assign_sites(calls, args.merge_bp)
    sites = collapse_sites(d)
    sites["subfamily"] = sites.site.map(d.groupby("site").MEISUBFAMILY.first())
    sites = sites[sites.family.isin(families)].copy()

    tb = twobitreader.TwoBitFile(args.twobit)

    # Per-site features.
    recs = []
    for r in sites.itertuples():
        seq = fetch(tb, r.chrom, int(r.pos), args.window)
        f = window_features(seq) if seq else None
        if f is None:
            continue
        sub = parse_subfamily(r.subfamily)
        recs.append({"site": r.site, "family": r.family, "chrom": r.chrom,
                     "pos": r.pos, "n_samples": r.n_samples, "orient": r.orient,
                     "subfamily": sub, "l1_age": l1_age_class(sub) if r.family == "LINE1" else None,
                     **f})
    site_feats = pd.DataFrame(recs)

    # Background.
    bg = sample_background(tb, sizes, gaps, args.window, args.n_bg, args.seed)
    bg.to_csv(out / "background_features.csv", index=False)

    # Comparisons: each family vs background, and each family vs LINE1.
    summaries = [compare(site_feats[site_feats.family == fam], bg, fam) for fam in families
                 if (site_feats.family == fam).any()]
    l1 = site_feats[site_feats.family == "LINE1"]
    for fam in families:
        if fam != "LINE1" and (site_feats.family == fam).any():
            summaries.append(compare(site_feats[site_feats.family == fam], l1, f"{fam}_vs_LINE1"))
    summary = pd.concat(summaries, ignore_index=True)

    # LINE-1 subfamily age gradient.
    age_rows = []
    for age, x in l1.groupby("l1_age"):
        age_rows.append({"l1_age": age, "n": len(x),
                         "gc_frac_median": float(x.gc_frac.median()),
                         "max_abs_g4_median": float(x.max_abs_g4.median())})
    age_df = pd.DataFrame(age_rows).sort_values("l1_age")

    tsd_df = tsd_composition(calls[calls.MEIFAMILY.isin(families)])

    site_feats.to_csv(out / "site_features.csv", index=False)
    summary.to_csv(out / "summary.csv", index=False)
    age_df.to_csv(out / "l1_age_gradient.csv", index=False)
    tsd_df.to_csv(out / "tsd_composition.csv", index=False)

    report = {
        "vcf_dir": str(args.vcf_dir), "window_bp": 2 * args.window + 1,
        "merge_bp": args.merge_bp, "n_bg": int(len(bg)),
        "n_genomes": int(calls["sample"].nunique()),
        "n_calls": int(len(calls)),
        "sites_by_family": sites.family.value_counts().to_dict(),
        "scored_sites_by_family": site_feats.family.value_counts().to_dict(),
        "background_gc_median": float(bg.gc_frac.median()),
        "tsd": tsd_df.to_dict(orient="records"),
        "l1_age_gradient": age_df.to_dict(orient="records"),
        "summary": summary.to_dict(orient="records"),
    }
    (out / "l1_gc_content_report.json").write_text(json.dumps(report, indent=2, default=str))
    if args.figure:
        make_figure(site_feats, bg, out)

    pd.set_option("display.width", 200)
    print("=== family vs genome background ===")
    print(summary[summary.set.isin(families)].to_string(index=False))
    print("\n=== LINE-1 age gradient ===")
    print(age_df.to_string(index=False))
    print("\n=== TSD (positive control) ===")
    print(tsd_df.to_string(index=False))
    print(f"\nbackground GC median = {bg.gc_frac.median():.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
