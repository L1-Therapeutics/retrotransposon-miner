"""Sense vs antisense nesting of MEIs inside same-family reference elements.

Answers William's two questions on one callset:

1. Do MEIs insert into their own family's reference copies more often than
   chance predicts (overall self-insertion enrichment)?
2. Do they prefer the SAME orientation as the reference element they land in,
   relative to insertions nested in the OPPOSITE orientation?

Nesting is recomputed here with the four-valued annotator on
``fix/nested-orientation-semantics`` rather than read from an existing column,
because every VCF/TSV written before this branch carries only the binary
``nested``/``unnested`` pair and therefore cannot express the antisense class.

The annotator needs only breakpoint coordinates and an rmsk table, so this
runs on a laptop: no BAM, no AWS, no pipeline re-run.

``nested_unknown`` is excluded from both the sense and antisense counts. It
means "nested, but orientation unresolvable", so folding it into antisense
would inflate antisense and bias exactly the comparison being made here.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from retro_miner.mei_support import (  # noqa: E402
    _annotate_nested_retrotransposon,
    _normalize_mei_family_token,
)

FAMILY_ORDER = ["ALU", "LINE1", "SVA"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", required=True, type=Path)
    p.add_argument("--rmsk", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-confidence", default="high", choices=["high", "high,medium", "all"])
    p.add_argument("--knownmei-only", action="store_true")
    return p.parse_args()


def load_candidates(path: Path, min_conf: str, knownmei_only: bool) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = df.rename(columns={"chrom": "chrom"})

    df["family"] = df["consensus_mei_family"].fillna("").map(_normalize_mei_family_token)
    df["orientation"] = df["consensus_insertion_orientation"].fillna("").str.strip()
    df["bp"] = pd.to_numeric(df["consensus_insertion_breakpoint_pos"], errors="coerce")
    df["ws"] = pd.to_numeric(df["window_start"], errors="coerce")
    df["we"] = pd.to_numeric(df["window_end"], errors="coerce")
    df["subfamily"] = df.get("mei_subfamily", pd.Series(dtype=str)).fillna("")

    keep = pd.Series(True, index=df.index)
    if min_conf != "all":
        allowed = set(min_conf.split(","))
        keep &= df["consensus_breakpoint_confidence_tier"].fillna("").isin(allowed)
    if knownmei_only:
        keep &= df.get("known_mei_polymorphism", pd.Series("0", index=df.index)).astype(str) == "1"
    df = df.loc[keep].copy()

    df["insertion_breakpoint_pos"] = df["bp"].fillna(0).astype(int)
    df["consensus_insertion_orientation"] = df["orientation"]
    df["mei_family_resolved"] = df["family"]
    # `_choose_event_family` does not read consensus_mei_family / family. On this
    # callset known_mei_polymorphism_family is empty for ~98% of high-conf rows,
    # so without this mirror almost every breakpoint is treated as family-less
    # and nesting collapses to unnested. Prefer the consensus family for the
    # annotator (same token the counts_by_family table uses).
    df["known_mei_polymorphism_family"] = df["family"]
    return df


def require_bedtools() -> None:
    import shutil

    if shutil.which("bedtools") is None:
        raise SystemExit(
            "bedtools not found on PATH. Install with one of:\n"
            "  brew install bedtools\n"
            "  conda activate rtm-miner   # environment.yml includes bedtools>=2.31\n"
            "Do not silently skip: nested annotation requires bedtools."
        )


def annotate(df: pd.DataFrame, rmsk: Path) -> pd.DataFrame:
    require_bedtools()
    out = _annotate_nested_retrotransposon(df, rmsk_table_path=rmsk)
    return out


def counts_by_family(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fam in FAMILY_ORDER:
        sub = df[df["family"] == fam]
        vc = sub["nested_same_class_orientation"].value_counts()
        n_sense = int(vc.get("nested_sense", 0))
        n_anti = int(vc.get("nested_antisense", 0))
        n_unk = int(vc.get("nested_unknown", 0))
        n_res = n_sense + n_anti
        rows.append(
            {
                "family": fam,
                "n_total": int(len(sub)),
                "n_nested": n_sense + n_anti + n_unk,
                "n_sense": n_sense,
                "n_antisense": n_anti,
                "n_unknown": n_unk,
                "n_resolved": n_res,
                "sense_fraction": (n_sense / n_res) if n_res else np.nan,
                "antisense_fraction": (n_anti / n_res) if n_res else np.nan,
            }
        )
    return pd.DataFrame(rows)


def holm(pvals: list[float]) -> list[float]:
    finite = [(i, p) for i, p in enumerate(pvals) if np.isfinite(p)]
    adj = [np.nan] * len(pvals)
    order = sorted(finite, key=lambda t: t[1])
    m = len(order)
    running = 0.0
    for rank, (idx, p) in enumerate(order):
        val = min(1.0, (m - rank) * p)
        running = max(running, val)
        adj[idx] = running
    return adj


def sense_vs_antisense(counts: pd.DataFrame) -> pd.DataFrame:
    from scipy import stats

    p_raw, rows = [], []
    for _, r in counts.iterrows():
        n_res, n_sense = int(r["n_resolved"]), int(r["n_sense"])
        n_anti = int(r["n_antisense"])
        if n_res > 0:
            p_one = float(stats.binomtest(n_sense, n_res, 0.5, alternative="greater").pvalue)
            p_two = float(stats.binomtest(n_sense, n_res, 0.5, alternative="two-sided").pvalue)
            status = "ok"
            note = ""
            if n_anti == 0:
                note = "antisense count is 0; not a null biological result (underpowered for depletion)"
            elif n_sense == 0:
                note = "sense count is 0; not a null biological result (underpowered for enrichment)"
        else:
            p_one = p_two = np.nan
            status = "underpowered"
            note = "n_resolved==0 (nested_unknown excluded; no sense/antisense events)"
        p_raw.append(p_one)
        rows.append(
            {
                "family": r["family"],
                "n_sense": n_sense,
                "n_antisense": n_anti,
                "n_resolved": n_res,
                "sense_fraction": r["sense_fraction"],
                "fold_vs_chance_sense": (r["sense_fraction"] / 0.5) if n_res else np.nan,
                "fold_vs_chance_antisense": (r["antisense_fraction"] / 0.5) if n_res else np.nan,
                "p_one_sided_greater": p_one,
                "p_two_sided": p_two,
                "status": status,
                "note": note,
            }
        )
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(p_raw)
    return out


def mei_base_content(rmsk: Path, chroms: set[str]) -> dict[tuple[str, str], int]:
    """Per-family, per-chromosome MEI base content, for the self-insertion null.

    Column layout is UCSC rmsk with leading bin:
    5 chrom, 6 start0, 7 end0, 9 strand, 10 repName, 11 repClass, 12 repFamily.
    Family is resolved exactly as the annotator does, by normalising
    repName + repClass + repFamily together, so the null covers the same
    intervals the annotator was able to nest into.
    """
    content: dict[tuple[str, str], int] = {}
    with gzip_open(rmsk) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13 or not parts[5].startswith("chr"):
                continue
            chrom = parts[5]
            if chrom not in chroms:
                continue
            try:
                start0, end0 = int(parts[6]), int(parts[7])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            fam = _normalize_mei_family_token(f"{parts[10]} {parts[11]} {parts[12]}")
            if not fam:
                continue
            key = (fam, chrom)
            content[key] = content.get(key, 0) + (end0 - start0)
    return content


def gzip_open(path: Path):
    import gzip

    return gzip.open(path, "rt", encoding="utf-8", errors="replace")


def genome_sizes(path: Path) -> dict[str, int]:
    """Parse a two-column UCSC chrom.sizes file into {chrom: length}."""
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError:
                continue
    if not sizes:
        raise ValueError(f"no chrom sizes parsed from {path}")
    return sizes


def selfinsertion_null(
    df: pd.DataFrame,
    content: dict[tuple[str, str], int],
    chrom_sizes: dict[str, int],
) -> pd.DataFrame:
    """Expected nested-in-same-family fraction if calls landed at random.

    Per call, the chance of landing inside a same-family element is that
    family's base content divided by the length of the chromosome the call
    sits on. Each call therefore gets its own null probability, and the
    expected nested count is the sum of those rather than n * one global rate,
    because the calls are not spread evenly across chromosomes.

    The uniform null is the conservative choice for Alu, where GC-matching
    partly removes the signal, and the optimistic choice for LINE-1/SVA.
    Each row records which null produced it. A GC-matched null is NOT
    computed: it needs per-window GC content, which rmsk does not carry.
    """
    from scipy import stats

    rows = []
    for fam in FAMILY_ORDER:
        sub = df[df["family"] == fam]
        if sub.empty:
            continue
        usable = {
            c: content.get((fam, c), 0) / s
            for c, s in chrom_sizes.items()
            if s > 0 and content.get((fam, c), 0) > 0
        }
        if not usable:
            continue
        # Per-call null probability, using that call's own chromosome.
        p_call = np.array(
            [usable.get(str(ch), 0.0) for ch in sub["chrom"].astype(str)], dtype=float
        )
        exp = float(p_call.sum())
        n = len(sub)
        n_nested = int((sub["nested_same_class_orientation"] != "unnested").sum())
        obs = n_nested / n
        p = float(stats.binomtest(n_nested, n, exp / n, alternative="greater").pvalue)
        rows.append(
            {
                "family": fam,
                "n_calls": n,
                "n_nested": n_nested,
                "observed_nested_fraction": obs,
                "expected_nested_fraction_uniform": exp / n,
                "expected_nested_count": exp,
                "enrichment_uniform": (obs / (exp / n)) if exp > 0 else np.nan,
                "p_greater_uniform": p,
                "null": "uniform",
                "n_chroms_in_null": len(usable),
            }
        )
    return pd.DataFrame(rows)


def write_report(outdir: Path, counts: pd.DataFrame, sva: pd.DataFrame, nulls: pd.DataFrame, meta: dict) -> None:
    def fmt(v, spec=".4g"):
        return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, spec)

    lines = [
        "# Nested sense vs antisense orientation",
        "",
        "Every number below is read back from the CSVs in this directory at write time.",
        "",
        "## Inputs",
        "",
        f"- candidates: `{meta['candidates']}`",
        f"- rmsk: `{meta['rmsk']}`",
        f"- filter: confidence tier in {{{meta['min_conf']}}}, knownmei_only={meta['knownmei_only']}",
        f"- calls analysed: {meta['n_calls']}",
        "",
        "## Counts by family",
        "",
        "| family | n_total | n_sense | n_antisense | n_unknown | n_resolved |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in counts.iterrows():
        lines.append(
            f"| {r['family']} | {int(r['n_total'])} | {int(r['n_sense'])} | "
            f"{int(r['n_antisense'])} | {int(r['n_unknown'])} | {int(r['n_resolved'])} |"
        )
    lines += [
        "",
        "`nested_unknown` is excluded from both nested columns, so n_resolved != n_nested whenever unknown > 0.",
        "",
        "## Sense vs antisense (H0: 50/50 among resolved nested calls)",
        "",
        "| family | n_resolved | sense_frac | fold_sense | fold_antisense | p_one | p_holm |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in sva.iterrows():
        lines.append(
            f"| {r['family']} | {int(r['n_resolved'])} | {fmt(r['sense_fraction'])} | "
            f"{fmt(r['fold_vs_chance_sense'])} | {fmt(r['fold_vs_chance_antisense'])} | "
            f"{fmt(r['p_one_sided_greater'])} | {fmt(r['p_holm'])} |"
        )
    lines += [
        "",
        "fold = fraction / 0.5, so 1.0 is exactly chance and >1.0 favours that class.",
        "",
        "## Overall self-insertion vs chance",
        "",
        "| family | n_calls | n_nested | obs_frac | exp_frac | enrichment | p |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in nulls.iterrows():
        lines.append(
            f"| {r['family']} | {int(r['n_calls'])} | {int(r['n_nested'])} | "
            f"{fmt(r['observed_nested_fraction'])} | {fmt(r['expected_nested_fraction_uniform'])} | "
            f"{fmt(r['enrichment_uniform'])} | {fmt(r['p_greater_uniform'])} |"
        )
    lines += [
        "",
        "Null is uniform-random placement. GC-matched null is NOT computed here:",
        "it needs per-window GC content, which the rmsk table does not carry.",
        "",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = load_candidates(args.candidates, args.min_confidence, args.knownmei_only)
    if df.empty:
        print("no calls after filtering", file=sys.stderr)
        return 1
    print(f"loaded {len(df)} calls", flush=True)

    ann = annotate(df, args.rmsk)
    print(f"annotated: {ann['nested_same_class_orientation'].value_counts().to_dict()}", flush=True)

    counts = counts_by_family(ann)
    sva = sense_vs_antisense(counts)

    chroms = set(ann["chrom"].astype(str)) - {""}
    gsz = args.outdir / "hg38.chrom.sizes"
    if not gsz.exists():
        urllib.request.urlretrieve(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", gsz
        )
    chrom_sizes = genome_sizes(gsz)
    print(f"chrom sizes for {len(chrom_sizes)} chromosomes", flush=True)
    print("scanning rmsk for MEI base content", flush=True)
    content = mei_base_content(args.rmsk, chroms)
    print(f"MEI content for {len(content)} (family, chrom) pairs", flush=True)
    nulls = selfinsertion_null(ann, content, chrom_sizes)

    ann.to_csv(args.outdir / "per_call_annotated.csv", index=False)
    counts.to_csv(args.outdir / "counts_by_family.csv", index=False)
    sva.to_csv(args.outdir / "sense_vs_antisense.csv", index=False)
    nulls.to_csv(args.outdir / "selfinsertion_enrichment.csv", index=False)

    write_report(
        args.outdir,
        counts,
        sva,
        nulls,
        {
            "candidates": args.candidates,
            "rmsk": args.rmsk,
            "min_conf": args.min_confidence,
            "knownmei_only": args.knownmei_only,
            "n_calls": len(ann),
        },
    )
    print(counts.to_string(index=False))
    print()
    print(sva.to_string(index=False))
    print()
    print(nulls.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())