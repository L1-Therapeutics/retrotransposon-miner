#!/usr/bin/env python3
"""Window-overlap recall of candidate loci vs every known MELT/ONT MEI.

Reads a loci TSV (chrom, window_start, window_end) and one or more truth
TSVs (chrom, pos [, id, family, source]). Reports every site, not just one
example locus.

  PYTHONPATH=src python scripts/eval_loci_vs_known_mei.py \
    --loci candidate_loci.tsv \
    --truth melt_ins.tsv ont_ins.tsv \
    --pad 200 \
    --out recall.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _norm_chrom(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        return "chr" + text[3:]
    return f"chr{text}"


def load_loci(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    need = {"chrom", "window_start", "window_end"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"{path} missing columns {sorted(missing)}")
    out = df.loc[:, ["chrom", "window_start", "window_end"]].copy()
    out["chrom"] = out["chrom"].map(_norm_chrom)
    out["window_start"] = pd.to_numeric(out["window_start"], errors="coerce").fillna(0).astype(int)
    out["window_end"] = pd.to_numeric(out["window_end"], errors="coerce").fillna(0).astype(int)
    return out


def load_truth(path: Path, source_fallback: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    pos_col = "pos" if "pos" in df.columns else ("start" if "start" in df.columns else "")
    if "chrom" not in df.columns or not pos_col:
        raise SystemExit(f"{path} needs chrom and pos/start")
    out = pd.DataFrame(
        {
            "chrom": df["chrom"].map(_norm_chrom),
            "pos": pd.to_numeric(df[pos_col], errors="coerce").fillna(0).astype(int),
            "id": df["id"] if "id" in df.columns else df.get("variant_id", ""),
            "family": df["family"] if "family" in df.columns else "",
            "source": df["source"] if "source" in df.columns else source_fallback,
        }
    )
    return out.loc[out["pos"] > 0].copy()


def annotate_hits(truth: pd.DataFrame, loci: pd.DataFrame, pad: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for rec in truth.itertuples(index=False):
        chrom = str(rec.chrom)
        pos = int(rec.pos)
        hits = loci.loc[
            (loci["chrom"] == chrom)
            & (loci["window_start"] - int(pad) <= pos)
            & (loci["window_end"] + int(pad) >= pos)
        ]
        if hits.empty:
            rows.append(
                {
                    "source": rec.source,
                    "id": rec.id,
                    "family": rec.family,
                    "chrom": chrom,
                    "pos": pos,
                    "recovered": False,
                    "n_windows": 0,
                    "window_start": 0,
                    "window_end": 0,
                    "window_span": 0,
                }
            )
            continue
        best = hits.assign(span=hits["window_end"] - hits["window_start"] + 1).sort_values(
            ["span", "window_start"], kind="mergesort"
        ).iloc[0]
        rows.append(
            {
                "source": rec.source,
                "id": rec.id,
                "family": rec.family,
                "chrom": chrom,
                "pos": pos,
                "recovered": True,
                "n_windows": int(len(hits)),
                "window_start": int(best["window_start"]),
                "window_end": int(best["window_end"]),
                "window_span": int(best["span"]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--loci", type=Path, required=True)
    p.add_argument("--truth", type=Path, nargs="+", required=True)
    p.add_argument("--pad", type=int, default=200)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    loci = load_loci(args.loci)
    truth_frames = [
        load_truth(path, source_fallback=path.stem) for path in args.truth
    ]
    truth = pd.concat(truth_frames, ignore_index=True)
    detail = annotate_hits(truth, loci, pad=int(args.pad))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(args.out, sep="\t", index=False)

    print(f"loci\t{len(loci)}")
    print(f"truth_sites\t{len(detail)}")
    print(f"recovered\t{int(detail['recovered'].sum())}/{len(detail)}")
    for source, grp in detail.groupby("source", sort=False):
        n = len(grp)
        hit = int(grp["recovered"].sum())
        print(f"{source}\t{hit}/{n}")
    missed = detail.loc[~detail["recovered"]]
    if not missed.empty:
        print("missed_sites")
        for rec in missed.itertuples(index=False):
            print(f"  {rec.source}\t{rec.id}\t{rec.chrom}:{rec.pos}")
    wide = detail.loc[detail["recovered"] & (detail["window_span"] > 1000)]
    if not wide.empty:
        print("wide_recovered_windows_gt_1kb")
        for rec in wide.itertuples(index=False):
            print(
                f"  {rec.source}\t{rec.id}\t{rec.chrom}:{rec.pos}\t"
                f"{rec.window_start}-{rec.window_end}\tspan={rec.window_span}"
            )


if __name__ == "__main__":
    main()
