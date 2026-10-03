#!/usr/bin/env python3
"""Validate NM-derived panel divergence against rmsk milliDiv.

Pairs each candidate call that overlaps a reference element of the same
family with (a) the divergence derived from panel-alignment NM values and
(b) the reference element's ``milliDiv / 1000``.  Writes a TSV of paired
values and prints Pearson/Spearman correlation plus median absolute
difference.

  PYTHONPATH=src python scripts/validate_divergence_vs_rmsk.py \
    --calls candidate_loci.tsv \
    --rmsk rmsk.txt \
    --nm-col nm_mean \
    --alnlen-col aln_len_mean \
    --family-col consensus_mei_family \
    --out divergence_pairs.tsv \
    --plot
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from retro_miner.mei_support import _write_rmsk_mei_bed, panel_divergence_from_nm


def _norm_chrom(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        return "chr" + text[3:]
    return f"chr{text}"


def _load_calls(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"chrom", "window_start", "window_end"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{path} missing columns {sorted(missing)}")
    out = df.loc[:, list(df.columns)].copy()
    out["chrom"] = out["chrom"].map(_norm_chrom)
    for col in ("window_start", "window_end"):
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)
    return out


def _write_calls_bed(df: pd.DataFrame, bed_path: Path, nm_col: str, alnlen_col: str, family_col: str) -> None:
    rows = []
    for rec in df.itertuples(index=False):
        chrom = _norm_chrom(getattr(rec, "chrom", ""))
        ws = int(getattr(rec, "window_start", 0))
        we = int(getattr(rec, "window_end", 0))
        if not chrom or ws <= 0 or we <= 0:
            continue
        start0 = ws - 1
        family = str(getattr(rec, family_col, "") or "").strip()
        nm = pd.to_numeric(getattr(rec, nm_col, np.nan), errors="coerce")
        aln_len = pd.to_numeric(getattr(rec, alnlen_col, np.nan), errors="coerce")
        nm = int(nm) if pd.notna(nm) else -1
        aln_len = int(aln_len) if pd.notna(aln_len) else -1
        rows.append(f"{chrom}\t{start0}\t{we}\t{family}\t{nm}\t{aln_len}\n")
    bed_path.write_text("".join(rows), encoding="utf-8")


def _run_bedtools_intersect(a: Path, b: Path) -> str:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found on PATH; install bedtools to run this script")
    proc = subprocess.run(
        ["bedtools", "intersect", "-a", str(a), "-b", str(b), "-wa", "-wb", "-f", "0.5", "-r", "-e"],
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout


def _parse_intersect(
    raw: str,
    a_nm_col: str = "nm",
    a_alnlen_col: str = "aln_len",
    a_family_col: str = "family",
) -> pd.DataFrame:
    if not raw.strip():
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                "call_family",
                "nm",
                "aln_len",
                "nm_divergence",
                "rmsk_milli_div",
                "rmsk_milli_div_fraction",
            ]
        )
    cols = [
        "a_chrom", "a_start0", "a_end0",
        "call_family", "nm", "aln_len",
        "b_chrom", "b_start0", "b_end0",
        "b_repName", "b_length", "b_strand",
        "b_repClass", "b_repFamily", "b_normFamily", "b_milliDiv",
    ]
    df = pd.read_csv(
        pd.io.common.StringIO(raw),
        sep="\t",
        header=None,
        names=cols,
        dtype={"b_milliDiv": object},
    )
    df["nm"] = pd.to_numeric(df["nm"], errors="coerce").fillna(-1).astype(int)
    df["aln_len"] = pd.to_numeric(df["aln_len"], errors="coerce").fillna(-1).astype(int)
    df["b_milliDiv"] = pd.to_numeric(df["b_milliDiv"], errors="coerce").fillna(-1).astype(int)
    df["call_family"] = df["call_family"].fillna("").astype(str)
    df["b_normFamily"] = df["b_normFamily"].fillna("").astype(str)
    df = df.loc[df["call_family"].ne("") & df["b_normFamily"].ne("")].copy()
    df = df.loc[df["call_family"].str.upper() == df["b_normFamily"].str.upper()].copy()
    df["nm_divergence"] = df.apply(
        lambda r: panel_divergence_from_nm(int(r["nm"]), int(r["aln_len"])), axis=1
    )
    df = df.loc[df["nm_divergence"].notna()].copy()
    df["rmsk_milli_div"] = df["b_milliDiv"].where(df["b_milliDiv"] >= 0, np.nan)
    df["rmsk_milli_div_fraction"] = df["rmsk_milli_div"] / 1000.0
    df["window_start"] = df["a_start0"] + 1
    df["window_end"] = df["a_end0"]
    df = df.rename(columns={"a_chrom": "chrom", "call_family": "family"})
    out = df.loc[:, ["chrom", "window_start", "window_end", "family", "nm_divergence", "rmsk_milli_div", "rmsk_milli_div_fraction"]].copy()
    return out.reset_index(drop=True)


def _pearson(x: np.ndarray, y: np.ndarray) -> float | None:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    x = x[mask]
    y = y[mask]
    n = len(x)
    sx = x.sum()
    sy = y.sum()
    sx2 = (x ** 2).sum()
    sy2 = (y ** 2).sum()
    sxy = (x * y).sum()
    num = n * sxy - sx * sy
    den = np.sqrt(max(0.0, (n * sx2 - sx ** 2) * (n * sy2 - sy ** 2)))
    if den == 0:
        return None
    return float(num / den)


def _spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    try:
        from scipy.stats import spearmanr

        r, _ = spearmanr(x[mask], y[mask])
        return float(r)
    except Exception:
        return None


def _median_abs_diff(x: np.ndarray, y: np.ndarray) -> float | None:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() == 0:
        return None
    diffs = np.abs(x[mask] - y[mask])
    return float(np.median(diffs))


def _write_plot(pairs: pd.DataFrame, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(pairs["rmsk_milli_div_fraction"], pairs["nm_divergence"], alpha=0.6, s=20, edgecolors="none")
    lo = min(pairs["rmsk_milli_div_fraction"].min(), pairs["nm_divergence"].min())
    hi = max(pairs["rmsk_milli_div_fraction"].max(), pairs["nm_divergence"].max())
    pad = (hi - lo) * 0.05 or 0.01
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "r--", lw=1, label="y=x")
    ax.set_xlabel("rmsk milliDiv / 1000")
    ax.set_ylabel("NM-derived divergence (nm / aln_len)")
    ax.set_title("Panel divergence vs rmsk milliDiv")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--calls", type=Path, required=True, help="Candidate loci TSV")
    p.add_argument("--rmsk", type=Path, required=True, help="UCSC RepeatMasker rmsk table")
    p.add_argument("--nm-col", default="nm_mean", help="NM column in calls TSV (default: nm_mean)")
    p.add_argument("--alnlen-col", default="aln_len_mean", help="Alignment length column in calls TSV (default: aln_len_mean)")
    p.add_argument("--family-col", default="consensus_mei_family", help="Family column in calls TSV (default: consensus_mei_family)")
    p.add_argument("--out", type=Path, required=True, help="Output paired TSV")
    p.add_argument("--plot", action="store_true", help="Write scatter plot PNG alongside --out")
    args = p.parse_args()

    calls = _load_calls(args.calls)
    if args.nm_col not in calls.columns:
        raise SystemExit(f"calls TSV missing NM column: {args.nm_col}")
    if args.alnlen_col not in calls.columns:
        raise SystemExit(f"calls TSV missing alignment-length column: {args.alnlen_col}")
    if args.family_col not in calls.columns:
        raise SystemExit(f"calls TSV missing family column: {args.family_col}")

    with tempfile.TemporaryDirectory(prefix="rtm_divergence_") as tmpdir:
        tmp = Path(tmpdir)
        calls_bed = tmp / "calls.bed"
        rmsk_bed = tmp / "rmsk.bed"
        _write_calls_bed(calls, calls_bed, args.nm_col, args.alnlen_col, args.family_col)
        n_rmsk = _write_rmsk_mei_bed(args.rmsk, rmsk_bed)
        if n_rmsk == 0:
            raise SystemExit(f"No MEI entries written from {args.rmsk}")
        raw_intersect = _run_bedtools_intersect(calls_bed, rmsk_bed)
        pairs = _parse_intersect(raw_intersect)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.out, sep="\t", index=False)

    n = len(pairs)
    x = pairs["nm_divergence"].to_numpy()
    y = pairs["rmsk_milli_div_fraction"].to_numpy()
    r_pearson = _pearson(x, y)
    r_spearman = _spearman(x, y)
    mad = _median_abs_diff(x, y)

    print(f"matched\t{n}")
    if r_pearson is not None:
        print(f"pearson_r\t{r_pearson:.4f}")
    else:
        print("pearson_r\tNA")
    if r_spearman is not None:
        print(f"spearman_r\t{r_spearman:.4f}")
    else:
        print("spearman_r\tNA")
    if mad is not None:
        print(f"median_abs_diff\t{mad:.6f}")
    else:
        print("median_abs_diff\tNA")

    if args.plot and n > 0:
        plot_path = args.out.with_suffix(".png")
        _write_plot(pairs, plot_path)
        print(f"plot\t{plot_path}")


if __name__ == "__main__":
    main()
