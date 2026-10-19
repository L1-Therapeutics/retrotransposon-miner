"""Where intronic Alu insertions sit relative to splice sites, by orientation.

Exonization needs an ANTISENSE intronic Alu whose right arm supplies a 3' splice
site to the exon downstream of it (Gal-Mark 2008; Lev-Maor 2008). Lev-Maor
reported that fixed antisense Alus are DEPLETED within 150 bp upstream of exons.
This asks the same question of our polymorphic calls, which are young and have
had far less time under selection.

Usage:
  python alu_splice_distance.py <vcf_dir> <gencode.gtf.gz> --out-dir out
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

INFO_RE = re.compile(r"([^=;]+)=([^;]*)")
WINDOW_BP = 150


MAIN_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def parse_info(s: str) -> dict:
    return {m.group(1): m.group(2) for m in INFO_RE.finditer(s)}


def load_per_copy_csv(path: str) -> pd.DataFrame:
    """Load precomputed per-copy Alu table (e.g. alu_per_copy.csv.gz).

    Preserves 0-based ``[start, end)`` interval boundaries alongside the
    midpoint ``(start + end) // 2 + 1`` so ``classify`` can use exact
    interval-to-exon distance for reference copies while remaining backward
    compatible with point-coordinate VCF calls.
    """
    raw = pd.read_csv(path)
    start = raw["start"].astype(np.int64)
    end = raw["end"].astype(np.int64)
    pos = (start + end) // 2 + 1
    return pd.DataFrame(
        {
            "chrom": raw["chrom"].astype(str),
            "family": "ALU",
            "pos": pos.astype(np.int64),
            "start": start,
            "end": end,
            "orient": raw["elem_strand"].astype(str),
            "n_calls": 1,
        }
    )


def load_rmsk_sites(rmsk_path: str, rep_name: str = "AluY", chroms: set[str] = MAIN_CHROMS) -> pd.DataFrame:
    """Load fixed reference RepeatMasker copies for ``rep_name`` (or prefix ``Alu*``).

    If ``rep_name`` ends with ``*`` (e.g. ``Alu*``), matches any repeat whose
    name starts with ``rep_name[:-1]``; otherwise requires an exact match.
    """
    prefix = rep_name[:-1] if rep_name.endswith("*") else None
    rows = []
    opener = gzip.open if str(rmsk_path).endswith(".gz") else open
    with opener(rmsk_path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16 or f[5] not in chroms:
                continue
            rname = f[10]
            if prefix is not None:
                if not rname.startswith(prefix):
                    continue
            elif rname != rep_name:
                continue
            start, end, strand = int(f[6]), int(f[7]), f[9]
            rows.append(
                {
                    "chrom": f[5],
                    "family": rname,
                    "pos": (start + end) // 2 + 1,
                    "start": start,
                    "end": end,
                    "orient": strand,
                    "n_calls": 1,
                }
            )
    return pd.DataFrame(rows)


def load_calls(vcf_dir: str) -> pd.DataFrame:
    """chrom, pos, family, orient, sample for every non-header VCF record."""
    paths = sorted(glob.glob(os.path.join(vcf_dir, "*.vcf")))
    if not paths:
        raise FileNotFoundError(f"no *.vcf in {vcf_dir}")
    rows = []
    for path in paths:
        sample = os.path.basename(path).split("_")[0]
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                info = parse_info(f[7])
                rows.append(
                    {
                        "sample": sample,
                        "chrom": f[0],
                        "pos": int(f[1]),
                        "family": info.get("MEIFAMILY", ""),
                        "orient": info.get("ORIENT", ""),
                    }
                )
    return pd.DataFrame(rows)


def merge_sites(calls: pd.DataFrame, merge_bp: int = 50) -> pd.DataFrame:
    """Collapse calls seen in several genomes into unique sites.

    Sites are merged within a family on one chromosome when consecutive
    positions are within merge_bp. Orientation is the majority call.
    """
    out = []
    for (chrom, fam), g in calls.groupby(["chrom", "family"], sort=False):
        g = g.sort_values("pos")
        start = prev = None
        bucket: list = []
        for pos, orient in zip(g.pos.values, g.orient.values):
            if start is None or pos - prev > merge_bp:
                if bucket:
                    out.append(_collapse(chrom, fam, bucket))
                bucket = []
                start = pos
            bucket.append((pos, orient))
            prev = pos
        if bucket:
            out.append(_collapse(chrom, fam, bucket))
    return pd.DataFrame(out)


def _collapse(chrom: str, fam: str, bucket: list) -> dict:
    positions = [p for p, _ in bucket]
    orients = [o for _, o in bucket if o in ("+", "-")]
    orient = ""
    if orients:
        plus = orients.count("+")
        minus = orients.count("-")
        orient = "+" if plus > minus else ("-" if minus > plus else "")
    return {
        "chrom": chrom,
        "family": fam,
        "pos": int(np.median(positions)),
        "orient": orient,
        "n_calls": len(bucket),
    }


def load_gencode(path: str) -> tuple[dict, dict]:
    """Return per-chromosome gene records and sorted exon-boundary arrays.

    genes[chrom] -> list of (start, end, strand)
    exons[chrom] -> sorted np.ndarray of (start, end) rows
    Coordinates are 0-based half-open.
    """
    genes = defaultdict(list)
    exon_rows = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("gene", "exon"):
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if feat == "gene":
                genes[chrom].append((start, end, strand))
            else:
                exon_rows[chrom].append((start, end))
    exons = {}
    for chrom, rows in exon_rows.items():
        a = np.array(sorted(rows), dtype=np.int64)
        exons[chrom] = union_intervals(a)
    return dict(genes), exons


def union_intervals(a: np.ndarray) -> np.ndarray:
    """Merge overlapping [start, end) rows; input must be sorted by start."""
    if len(a) == 0:
        return a
    out = [a[0].copy()]
    for s, e in a[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append(np.array([s, e]))
    return np.array(out, dtype=np.int64)


def single_strand_gene(gene_list, pos: int):
    """Strand of the gene(s) containing pos, or None when 0 or mixed strands."""
    hits = {st for s, e, st in gene_list if s <= pos < e}
    if len(hits) != 1:
        return None
    return hits.pop()


def exon_context(exons: np.ndarray, pos: int, strand: str):
    """Distance to the nearest exon and whether pos is 5' or 3' of it.

    Returns (in_exon, distance_bp, side) where side is 'upstream' when pos lies
    on the 5' side of the nearest exon in transcription direction.
    """
    if len(exons) == 0:
        return False, None, None
    starts = exons[:, 0]
    i = int(np.searchsorted(starts, pos, side="right")) - 1
    if 0 <= i < len(exons) and exons[i][0] <= pos < exons[i][1]:
        return True, 0, None
    best_d, best_side = None, None
    for j in (i, i + 1):
        if not (0 <= j < len(exons)):
            continue
        s, e = exons[j]
        if pos < s:
            d, genomic_side = s - pos, "left"
        else:
            d, genomic_side = pos - (e - 1), "right"
        if best_d is None or d < best_d:
            best_d = int(d)
            # on a + strand gene, being to the LEFT of an exon is upstream of it
            if strand == "+":
                best_side = "upstream" if genomic_side == "left" else "downstream"
            else:
                best_side = "downstream" if genomic_side == "left" else "upstream"
    return False, best_d, best_side


def interval_exon_context(exons: np.ndarray, start: int, end: int, strand: str):
    """Edge-to-exon distance for a 0-based half-open reference interval [start, end).

    Unlike a polymorphic VCF call (which is a 0-bp insertion point in the
    reference), a fixed reference Alu occupies ~300 bp of reference coordinates.
    Using the midpoint misclassifies junction-spanning exonized Alus (which are
    ~74-86% antisense per Lavi & Carmel 2018) as 'intronic within 150 bp'.
    Returns (in_exon, distance_bp, side) using the interval's nearest edge.
    """
    if len(exons) == 0:
        return False, None, None
    i = int(np.searchsorted(exons[:, 1], start + 1, side="left"))
    if i < len(exons) and exons[i][0] < end:
        return True, 0, None
    best_d, best_side = None, None
    for j in (i - 1, i):
        if not (0 <= j < len(exons)):
            continue
        s, e = exons[j]
        if end <= s:
            d, genomic_side = s - end, "left"
        else:
            d, genomic_side = start - e, "right"
        if best_d is None or d < best_d:
            best_d = int(d)
            if strand == "+":
                best_side = "upstream" if genomic_side == "left" else "downstream"
            else:
                best_side = "downstream" if genomic_side == "left" else "upstream"
    return False, best_d, best_side


def classify(sites: pd.DataFrame, genes: dict, exons: dict) -> pd.DataFrame:
    has_interval = "start" in sites.columns and "end" in sites.columns
    recs = []
    for r in sites.itertuples(index=False):
        gl = genes.get(r.chrom)
        if not gl:
            continue
        strand = single_strand_gene(gl, r.pos)
        if strand is None or r.orient not in ("+", "-"):
            continue
        ex_arr = exons.get(r.chrom, np.empty((0, 2), np.int64))
        if has_interval and pd.notna(r.start) and pd.notna(r.end):
            in_exon, dist, side = interval_exon_context(ex_arr, int(r.start), int(r.end), strand)
        else:
            in_exon, dist, side = exon_context(ex_arr, r.pos, strand)
        recs.append(
            {
                "chrom": r.chrom,
                "pos": r.pos,
                "family": r.family,
                "orient": r.orient,
                "gene_strand": strand,
                "sense": r.orient == strand,
                "in_exon": in_exon,
                "dist_to_exon": dist,
                "side": side,
            }
        )
    return pd.DataFrame(recs)


def zone(row) -> str:
    if row.in_exon:
        return "exonic"
    if row.dist_to_exon is None:
        return "no_exon"
    if row.dist_to_exon <= WINDOW_BP:
        return f"within_{WINDOW_BP}_{row.side}"
    return "intronic_far"


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["zone"] = df.apply(zone, axis=1)
    rows = []
    for z, g in df.groupby("zone"):
        n = len(g)
        anti = int((~g.sense).sum())
        rows.append({"zone": z, "n": n, "antisense": anti, "antisense_frac": anti / n if n else np.nan})
    return pd.DataFrame(rows).sort_values("zone").reset_index(drop=True)


def window_test(df: pd.DataFrame, side: str) -> dict:
    """Fisher test: antisense fraction in the near-exon window vs far intronic."""
    df = df.copy()
    df["zone"] = df.apply(zone, axis=1)
    near = df[df.zone == f"within_{WINDOW_BP}_{side}"]
    far = df[df.zone == "intronic_far"]
    table = [
        [int((~near.sense).sum()), int(near.sense.sum())],
        [int((~far.sense).sum()), int(far.sense.sum())],
    ]
    odds, p = fisher_exact(table)
    return {
        "side": side,
        "n_near": len(near),
        "antisense_near": table[0][0],
        "antisense_frac_near": table[0][0] / len(near) if len(near) else np.nan,
        "n_far": len(far),
        "antisense_frac_far": table[1][0] / len(far) if len(far) else np.nan,
        "odds_ratio": odds,
        "p": p,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf_dir", nargs="?", default=None, help="Directory of VCF files (optional when --per-copy-csv or --rmsk is used)")
    ap.add_argument("gencode", help="GENCODE annotation GTF(.gz)")
    ap.add_argument("--family", default="ALU")
    ap.add_argument("--per-copy-csv", default=None, help="Precomputed per-copy CSV(.gz) from alu_polya_asymmetry.py")
    ap.add_argument("--rmsk", default=None, help="UCSC RepeatMasker rmsk.txt(.gz) for fixed reference copies")
    ap.add_argument("--rep-name", default="AluY", help="RepeatMasker subfamily (e.g. AluY, AluSx, AluJo, or Alu*) when --rmsk is used")
    ap.add_argument("--out-dir", default="splice_out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)

    if a.per_copy_csv:
        sites = load_per_copy_csv(a.per_copy_csv)
    elif a.rmsk:
        sites = load_rmsk_sites(a.rmsk, a.rep_name)
    elif a.vcf_dir:
        calls = load_calls(a.vcf_dir)
        sites = merge_sites(calls[calls.family == a.family])
    else:
        raise ValueError("Provide vcf_dir, --per-copy-csv, or --rmsk")
    genes, exons = load_gencode(a.gencode)
    df = classify(sites, genes, exons)
    df.to_csv(os.path.join(a.out_dir, "alu_splice_sites.tsv"), sep="\t", index=False)

    summary = summarise(df)
    summary.to_csv(os.path.join(a.out_dir, "alu_splice_zones.tsv"), sep="\t", index=False)
    tests = pd.DataFrame([window_test(df, "upstream"), window_test(df, "downstream")])
    tests.to_csv(os.path.join(a.out_dir, "alu_splice_tests.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False))
    print()
    print(tests.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
