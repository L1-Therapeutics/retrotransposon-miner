"""GC content of the genome, measured from the reference sequence itself.

Question
--------
What is the GC content of the genome(s) we are working with, and where does the
MEI callset sit relative to it? The umbrella task is "L1 vs GC content". Finding
17 (artifacts/l1_literature_tracker.md) tested GC at 30 chr22 insertion sites;
Finding 19 quotes Medstrand et al. 2002, who bin the genome by GC fraction and
report that young L1 peak in the 38-42% GC band. Both of those need a genome-wide
*GC-content map* as the denominator, not a sample of insertion sites.

This script builds that map from the hg38 2bit and then reads the MEI calls out
of the VCFs on top of it.

What it computes
----------------
  1. Base composition per --bin-bp window for every requested chromosome
     (A/C/G/T/N counts), so GC is measured, never inferred.
  2. Genome-wide and per-chromosome GC fraction over A/C/G/T (N excluded), plus
     the distribution of window GC.
  3. MEI (LINE1 / ALU / SVA) insertion density - insertions per Mb of *equal-GC*
     sequence - as a function of window GC. This is the genome-wide version of
     Medstrand's GC-fraction histogram and needs no per-sample genome because the
     calls are already genome-wide coordinates on hg38.

Note on the VCFs
----------------
The per-sample VCFs are symbolic MEI calls (REF = N) and carry no sequence, so
they cannot give a per-sample genome GC. Every one of the ten samples is human,
so the reference hg38 GC measured here is the estimate for each of them; the
script reports the genome sequence it actually read so that is auditable.

Usage:
    python genome_gc_content.py HG38.2bit SIZES.json [--gaps GAPS.json]
        [--vcf-dir vcf] [--bin-bp 1000] [--out-dir genome_gc_out] [--figure]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_orientation_summary import assign_sites, collapse_sites, load_calls  # noqa: E402

BASES = {"A": 65, "C": 67, "G": 71, "T": 84, "N": 78}
FAMILIES = ("LINE1", "ALU", "SVA")


def count_bases_per_bin(tb, chrom: str, length: int, bin_bp: int,
                        chunk_bins: int = 8000) -> dict:
    """Per-bin A/C/G/T/N counts for one chromosome, read in bin-aligned chunks.

    A whole-chromosome numpy array would be ~0.25 GB for chr1 and several such
    allocations per chromosome; chunking at ~8 Mb keeps the peak tiny and the
    chunk boundaries land exactly on bin edges so no bin is split.
    """
    nbins = length // bin_bp
    out = {b: np.zeros(nbins, dtype=np.int64) for b in BASES}
    usable = nbins * bin_bp
    step = bin_bp * chunk_bins
    for start in range(0, usable, step):
        stop = min(start + step, usable)
        seq = str(tb[chrom][start:stop]).upper()           # 2bit is soft-masked
        if len(seq) != stop - start:                       # short read: stop cleanly
            break
        arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).reshape(-1, bin_bp)
        b0 = start // bin_bp
        b1 = b0 + arr.shape[0]
        for b, code in BASES.items():
            out[b][b0:b1] = (arr == code).sum(axis=1)
    return out


def gap_mask(chrom: str, nbins: int, bin_bp: int, gaps: dict | None) -> np.ndarray:
    """True where a bin overlaps an assembly gap (masked out of every stat)."""
    mask = np.zeros(nbins, dtype=bool)
    for s, e in (gaps or {}).get(chrom, []):
        lo = max(0, s // bin_bp)
        hi = min(nbins, -(-e // bin_bp))                   # ceil
        mask[lo:hi] = True
    return mask


def bin_table(tb, chrom: str, length: int, bin_bp: int, gaps: dict | None) -> pd.DataFrame:
    """One row per window: coordinates, base counts, GC fraction, validity."""
    counts = count_bases_per_bin(tb, chrom, length, bin_bp)
    n = len(counts["A"])
    a, c, g, t, nn = (counts[k] for k in "ACGTN")
    acgt = a + c + g + t
    gc = c + g
    valid = (acgt >= 0.5 * bin_bp) & ~gap_mask(chrom, n, bin_bp, gaps)
    return pd.DataFrame({
        "chrom": chrom,
        "bin": np.arange(n),
        "start": np.arange(n) * bin_bp,
        "end": np.minimum((np.arange(n) + 1) * bin_bp, length),
        "a": a, "c": c, "g": g, "t": t, "n": nn,
        "acgt": acgt,
        "gc": gc,
        "gc_frac": np.where(acgt > 0, gc / np.maximum(acgt, 1), np.nan),
        "valid": valid,
    })


def sum_scope(bins: pd.DataFrame, label: str) -> dict:
    v = bins[bins.valid]
    acgt = int(v.acgt.sum())
    gc = int(v.gc.sum())
    return {
        "scope": label,
        "bins_total": int(len(bins)),
        "bins_valid": int(len(v)),
        "bp_valid": acgt,
        "gc_bp": gc,
        "gc_frac": (gc / acgt) if acgt else float("nan"),
        "gc_frac_median_window": float(v.gc_frac.median()) if len(v) else float("nan"),
    }


def gc_histogram(bins: pd.DataFrame, step: float = 0.01) -> pd.DataFrame:
    """Base-pair mass in each GC bin (the genome GC landscape)."""
    v = bins[bins.valid]
    edges = np.arange(0.0, 1.0 + step, step)
    idx = np.clip((v.gc_frac.values / step).astype(int), 0, len(edges) - 2)
    bp = np.bincount(idx, weights=v.acgt.values, minlength=len(edges) - 1)
    nbins = np.bincount(idx, minlength=len(edges) - 1)
    total = bp.sum() or 1
    return pd.DataFrame({
        "gc_frac_low": edges[:-1],
        "gc_frac_high": edges[1:],
        "bp": bp.astype(np.int64),
        "n_windows": nbins.astype(np.int64),
        "bp_fraction": bp / total,
    })


def attach_gc_bin(sites: pd.DataFrame, bins: pd.DataFrame, bin_bp: int,
                  step: float = 0.01) -> pd.DataFrame:
    """Give each insertion the GC fraction of the window it falls in (NaN if masked)."""
    v = bins[bins.valid].copy()
    v["gc_bin"] = np.clip((v.gc_frac.values / step).astype(int), 0, int(1 / step) - 1)
    keys = pd.MultiIndex.from_arrays([v.chrom.values, v.bin.values])
    gc_lookup = pd.Series(v.gc_frac.values, index=keys)
    bin_lookup = pd.Series(v.gc_bin.values, index=keys)
    site_key = pd.MultiIndex.from_arrays(
        [sites.chrom.values, (sites.pos.values // bin_bp).astype(np.int64)])
    return sites.assign(gc_frac=site_key.map(gc_lookup).values,
                        gc_bin=site_key.map(bin_lookup).values)


def mei_density_by_gc(labeled: pd.DataFrame, bins: pd.DataFrame,
                      step: float = 0.01) -> pd.DataFrame:
    """Insertions per Mb of equal-GC sequence, per family.

    Density, not raw count: the valid genome is strongly bimodal in GC, so a raw
    count per GC bin would just restate where the GC-rich CpG islands are.
    """
    v = bins[bins.valid].copy()
    v["gc_bin"] = np.clip((v.gc_frac.values / step).astype(int), 0, int(1 / step) - 1)
    mass = v.groupby("gc_bin").acgt.sum()
    rows = []
    for idx, gc_low in enumerate(np.arange(0.0, 1.0, step)):
        bp = int(mass.get(idx, 0))
        row = {"gc_frac_low": round(float(gc_low), 2), "bp_valid": bp}
        for fam in FAMILIES:
            k = int(((labeled.family == fam) & (labeled.gc_bin == idx)).sum())
            row[f"{fam}_n"] = k
            row[f"{fam}_per_mb"] = (k / bp * 1e6) if bp else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def mei_gc_summary(labeled: pd.DataFrame, bins: pd.DataFrame,
                   band: tuple[float, float] = (0.38, 0.42)) -> pd.DataFrame:
    """Per-family local GC at insertion sites vs the genome, incl. a Medstrand band.

    Finding 19 records Medstrand et al. 2002 as reporting that young L1 peak in
    the 38-42% GC fraction. ``band_enrichment`` is (fraction of insertions in the
    band) / (fraction of the genome's bp in the band), so 1.0 means the family
    sits in that GC band exactly as often as sheer sequence availability allows.
    """
    genome_gc = float(bins.loc[bins.valid, "gc"].sum() / bins.loc[bins.valid, "acgt"].sum())
    in_band = bins.valid & bins.gc_frac.between(band[0], band[1], inclusive="left")
    genome_band_bp = float(bins.loc[in_band, "acgt"].sum() / bins.loc[bins.valid, "acgt"].sum())
    rows = []
    for fam in FAMILIES:
        x = labeled[(labeled.family == fam) & labeled.gc_frac.notna()]
        if not len(x):
            continue
        obs_band = float(x.gc_frac.between(band[0], band[1], inclusive="left").mean())
        rows.append({
            "family": fam,
            "n": int(len(x)),
            "mean_local_gc": float(x.gc_frac.mean()),
            "median_local_gc": float(x.gc_frac.median()),
            "genome_gc": genome_gc,
            "delta_vs_genome": float(x.gc_frac.mean() - genome_gc),
            f"frac_in_{band[0]:.2f}_{band[1]:.2f}": obs_band,
            "genome_frac_in_band": genome_band_bp,
            "band_enrichment": (obs_band / genome_band_bp) if genome_band_bp else float("nan"),
        })
    return pd.DataFrame(rows)


def robust_peak(dens: pd.DataFrame, min_bp: int = 5_000_000) -> dict:
    """GC bin with maximum density among bins carrying enough sequence to trust."""
    good = dens[dens.bp_valid >= min_bp]
    out = {}
    for fam in FAMILIES:
        col = f"{fam}_per_mb"
        if len(good) and good[col].notna().any():
            out[fam] = round(float(good.loc[good[col].idxmax(), "gc_frac_low"]), 2)
    return out


def make_figure(hist: pd.DataFrame, dens: pd.DataFrame, out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                                # pragma: no cover
        print(f"figure skipped: {exc}")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    mid = (hist.gc_frac_low + hist.gc_frac_high) / 2
    axes[0].fill_between(mid, hist.bp_fraction, color="0.6")
    axes[0].set_xlabel("window GC fraction")
    axes[0].set_ylabel("fraction of valid genome bp")
    axes[0].set_title("Genome GC landscape (hg38, 1 kb windows)")
    for fam, color in zip(FAMILIES, ("tab:red", "tab:blue", "tab:green")):
        axes[1].plot(dens.gc_frac_low, dens[f"{fam}_per_mb"], marker="o", ms=3,
                     label=fam, color=color)
    axes[1].set_xlabel("window GC fraction")
    axes[1].set_ylabel("insertions per Mb")
    axes[1].legend()
    axes[1].set_title("MEI density vs local GC (10 genomes, unique sites)")
    fig.tight_layout()
    fig.savefig(out / "genome_gc_content.png", dpi=130)
    print(f"wrote {out / 'genome_gc_content.png'}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("twobit")
    ap.add_argument("sizes_json")
    ap.add_argument("--gaps", default=None, help="assembly-gap intervals JSON")
    ap.add_argument("--vcf-dir", default=None, help="MEI VCF dir for the density overlay")
    ap.add_argument("--bin-bp", type=int, default=1000)
    ap.add_argument("--merge-bp", type=int, default=50)
    ap.add_argument("--chroms", default=None, help="comma list; default all in sizes")
    ap.add_argument("--out-dir", default="genome_gc_out")
    ap.add_argument("--figure", action="store_true")
    args = ap.parse_args(argv)

    import twobitreader

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sizes = json.loads(Path(args.sizes_json).read_text())
    gaps = json.loads(Path(args.gaps).read_text()) if args.gaps else None

    tb = twobitreader.TwoBitFile(args.twobit)
    chroms = ([c.strip() for c in args.chroms.split(",")] if args.chroms
              else [c for c in sizes if c in tb])

    frames = []
    for chrom in chroms:
        length = sizes[chrom]
        if chrom not in tb:
            print(f"skip {chrom}: not in 2bit")
            continue
        frames.append(bin_table(tb, chrom, length, args.bin_bp, gaps))
        print(f"binned {chrom} ({length:,} bp)")
    bins = pd.concat(frames, ignore_index=True)
    if _has_parquet():
        bins.to_parquet(out / "gc_bins.parquet")

    per_chrom = pd.DataFrame([sum_scope(bins[bins.chrom == c], c) for c in chroms])
    genome_scope = sum_scope(bins, "genome")

    hist = gc_histogram(bins)

    summary = pd.DataFrame([genome_scope])
    summary.to_csv(out / "genome_gc_summary.csv", index=False)
    per_chrom.to_csv(out / "gc_per_chromosome.csv", index=False)
    hist.to_csv(out / "gc_landscape.csv", index=False)

    dens = None
    gc_summary = None
    if args.vcf_dir:
        calls = load_calls(args.vcf_dir)
        sites = collapse_sites(assign_sites(calls, args.merge_bp))
        sites = sites[sites.family.isin(FAMILIES)].copy()
        labeled = attach_gc_bin(sites, bins, args.bin_bp)
        labeled.to_csv(out / "mei_site_gc.csv", index=False)
        dens = mei_density_by_gc(labeled, bins)
        dens.to_csv(out / "mei_density_by_gc.csv", index=False)
        gc_summary = mei_gc_summary(labeled, bins)
        gc_summary.to_csv(out / "mei_gc_summary.csv", index=False)

    report = {
        "twobit": str(args.twobit),
        "bin_bp": args.bin_bp,
        "n_chroms": len(chroms),
        "genome": genome_scope,
        "per_chromosome": per_chrom.to_dict(orient="records"),
    }
    if dens is not None:
        report["mei_density_peak_gc_bin"] = robust_peak(dens)
        report["mei_gc_summary"] = gc_summary.to_dict(orient="records")
    (out / "genome_gc_report.json").write_text(json.dumps(report, indent=2, default=str))
    if args.figure:
        make_figure(hist, dens if dens is not None else pd.DataFrame(
            {"gc_frac_low": hist.gc_frac_low, **{f"{f}_per_mb": np.nan for f in FAMILIES}}), out)

    pd.set_option("display.width", 200)
    print("\n=== genome GC ===")
    print(pd.DataFrame([genome_scope]).to_string(index=False))
    print("\n=== per chromosome ===")
    print(per_chrom.to_string(index=False))
    if dens is not None:
        print("\n=== MEI local GC vs genome ===")
        print(gc_summary.to_string(index=False))
        print("\nMEI density peak GC bin (>=5 Mb):", report.get("mei_density_peak_gc_bin"))
    return 0


def _has_parquet() -> bool:
    try:
        import pyarrow  # noqa: F401
        return True
    except Exception:
        return False


if __name__ == "__main__":
    sys.exit(main())
