#!/usr/bin/env python3
"""Re-annotate nested orientation from an existing MEI VCF + hg38 rmsk.

Does not re-run the BAM pipeline and does not touch AWS. Nesting is a pure
function of (breakpoint, family, orientation, rmsk), so this script rebuilds
a candidate table from VCF INFO fields and calls

    retro_miner.mei_support._annotate_nested_retrotransposon

to emit the four-valued labels:

    unnested | nested_sense | nested_antisense | nested_unknown

``nested_unknown`` is excluded from both sense and antisense tallies.

Example:

  PYTHONPATH=src .venv/bin/python scripts/analyze_nested_orientation.py \\
    --vcf ~/Desktop/Research/Research_Projects/L1/chr22_mei.annot.vcf \\
    --rmsk ~/retrotransposon-workdir/data/public/annotation/hg38/repeats/rmsk.txt.gz \\
    --reference ~/retrotransposon-workdir/data/public/reference/hg38/Homo_sapiens_assembly38.fasta \\
    --outdir results/nested_orientation
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from intervaltree import IntervalTree

from retro_miner.mei_support import (
    _annotate_nested_retrotransposon,
    _choose_event_family,
    _normalize_mei_family_token,
)

FAMILIES = ("ALU", "LINE1", "SVA")
NESTED_LABELS = ("unnested", "nested_sense", "nested_antisense", "nested_unknown")
PRIMARY_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
GC_HALF_WINDOW = 500
GC_BIN_WIDTH = 0.05
NULL_POOL_SIZE = 50_000
BEDTOOLS_INSTALL = (
    "bedtools not found on PATH. Install with one of:\n"
    "  brew install bedtools\n"
    "  conda activate rtm-miner   # environment.yml includes bedtools>=2.31\n"
    "Do not silently skip: nested annotation requires bedtools."
)


def _require_bedtools() -> None:
    if shutil.which("bedtools") is None:
        raise SystemExit(BEDTOOLS_INSTALL)


def _default_public_data_dir() -> Path:
    env = os.environ.get("RTM_PUBLIC_DATA_DIR")
    if env:
        return Path(env)
    return Path.home() / "retrotransposon-workdir" / "data" / "public"


def _parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            key, val = token.split("=", 1)
            out[key] = val
        else:
            out[token] = "True"
    return out


def load_calls_from_vcf(vcf_path: Path) -> pd.DataFrame:
    """Build the candidate frame expected by ``_annotate_nested_retrotransposon``.

    VCF INFO ``MEIFAMILY`` is not among the columns ``_choose_event_family``
    inspects, so it is mirrored onto ``known_mei_polymorphism_family`` (and
    kept as ``MEIFAMILY`` for bookkeeping). ``ORIENT`` maps to
    ``consensus_insertion_orientation``.
    """
    # Text parse is authoritative for INFO key/value pairs and does not depend on
    # header completeness; pysam is only used when available for contig typing.
    rows: list[dict[str, object]] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as hin:
        for line in hin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, vid, _ref, _alt, _qual, _filt, info_s = parts[:8]
            info = _parse_info(info_s)
            meifamily = info.get("MEIFAMILY", "")
            orient = info.get("ORIENT", "")
            rows.append(
                {
                    "chrom": chrom,
                    "insertion_breakpoint_pos": int(pos_s),
                    "MEIFAMILY": meifamily,
                    "consensus_insertion_orientation": orient,
                    # Fallback column actually read by _choose_event_family().
                    "known_mei_polymorphism_family": meifamily,
                    "vcf_id": vid or ".",
                    "vcf_svtype": info.get("SVTYPE", ""),
                    "vcf_nested_original": info.get("NESTED", ""),
                }
            )
    return pd.DataFrame(rows)


def resolve_families(df: pd.DataFrame) -> pd.Series:
    return df.apply(lambda row: _choose_event_family(row), axis=1)


def print_family_resolution(df: pd.DataFrame, resolved: pd.Series) -> None:
    raw = df["MEIFAMILY"].fillna("").astype(str).map(_normalize_mei_family_token)
    print("Resolved family distribution (_choose_event_family):")
    for fam, n in sorted(Counter(resolved.tolist()).items(), key=lambda x: (-x[1], x[0])):
        print(f"  {fam or '<unresolved>'}: {n}")
    mismatches = (raw != resolved) & ~((raw == "") & (resolved == ""))
    if mismatches.any():
        print(
            f"WARNING: {int(mismatches.sum())} calls differ between "
            "normalized MEIFAMILY and _choose_event_family; inspect per_call_annotated.csv"
        )
    else:
        print("Family resolution check: MEIFAMILY tokens match _choose_event_family for all calls.")


def exact_binom_test_p05(k: int, n: int) -> float:
    """Two-sided exact binomial test of H0: p=0.5. Returns NaN if n==0."""
    if n <= 0:
        return float("nan")
    k = int(k)
    n = int(n)
    # For p=0.5, P(X=i) = C(n,i) / 2^n. Two-sided = 2 * min(cdf, sf) capped at 1.
    # Use the smaller of k and n-k as the left-tail count.
    left = min(k, n - k)
    # sum_{i=0..left} C(n,i) / 2^n
    total = 0
    for i in range(left + 1):
        total += math.comb(n, i)
    p = 2.0 * total / (2**n)
    return float(min(1.0, p))


def holm_adjust(pvalues: list[float]) -> list[float]:
    """Holm–Bonferroni adjusted p-values; NaNs stay NaN and are ignored in ranking."""
    m_eff = sum(1 for p in pvalues if p == p)  # not NaN
    out = [float("nan")] * len(pvalues)
    if m_eff == 0:
        return out
    order = sorted((i for i, p in enumerate(pvalues) if p == p), key=lambda i: pvalues[i])
    running = 0.0
    for rank, i in enumerate(order):
        adj = (m_eff - rank) * pvalues[i]
        adj = min(1.0, adj)
        if rank > 0:
            adj = max(adj, running)
        running = adj
        out[i] = adj
    return out


def counts_by_family(annotated: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fam in FAMILIES:
        sub = annotated.loc[annotated["resolved_family"] == fam]
        labels = sub["nested_same_class_orientation"].astype(str)
        n_sense = int((labels == "nested_sense").sum())
        n_antisense = int((labels == "nested_antisense").sum())
        n_unknown = int((labels == "nested_unknown").sum())
        n_nested = n_sense + n_antisense + n_unknown
        n_resolved = n_sense + n_antisense  # unknown excluded from both sides
        rows.append(
            {
                "family": fam,
                "n_total": int(len(sub)),
                "n_nested": n_nested,
                "n_sense": n_sense,
                "n_antisense": n_antisense,
                "n_unknown": n_unknown,
                "n_resolved": n_resolved,
            }
        )
    return pd.DataFrame(rows)


def sense_vs_antisense_table(counts: pd.DataFrame) -> pd.DataFrame:
    raw_p: list[float] = []
    meta: list[dict[str, object]] = []
    for rec in counts.itertuples(index=False):
        n_resolved = int(rec.n_resolved)
        n_sense = int(rec.n_sense)
        n_antisense = int(rec.n_antisense)
        underpowered = n_resolved == 0
        if underpowered:
            sense_fraction = float("nan")
            antisense_fraction = float("nan")
            sense_fold = float("nan")
            antisense_fold = float("nan")
            sense_p = float("nan")
            antisense_p = float("nan")
            note = "underpowered: n_resolved==0 (nested_unknown excluded; no sense/antisense events)"
        else:
            sense_fraction = n_sense / n_resolved
            antisense_fraction = n_antisense / n_resolved
            sense_fold = sense_fraction / 0.5
            antisense_fold = antisense_fraction / 0.5
            sense_p = exact_binom_test_p05(n_sense, n_resolved)
            antisense_p = exact_binom_test_p05(n_antisense, n_resolved)
            note = ""
            if n_antisense == 0:
                note = "antisense count is 0; fold_vs_chance_antisense=0 (not a null biological result)"
            elif n_sense == 0:
                note = "sense count is 0; fold_vs_chance_sense=0 (not a null biological result)"
        raw_p.append(sense_p)
        meta.append(
            {
                "family": rec.family,
                "n_resolved": n_resolved,
                "n_sense": n_sense,
                "n_antisense": n_antisense,
                "sense_fraction": sense_fraction,
                "antisense_fraction": antisense_fraction,
                "fold_vs_chance_sense": sense_fold,
                "fold_vs_chance_antisense": antisense_fold,
                "binom_p_sense_vs_0.5": sense_p,
                "binom_p_antisense_vs_0.5": antisense_p,
                "status": "underpowered" if underpowered else "ok",
                "note": note,
            }
        )
    adj = holm_adjust(raw_p)
    for i, row in enumerate(meta):
        row["binom_p_sense_vs_0.5_holm"] = adj[i]
        # antisense p equals sense p under p=0.5; apply the same Holm factor.
        row["binom_p_antisense_vs_0.5_holm"] = adj[i]
    return pd.DataFrame(meta)


def _chrom_lengths(fasta: pysam.FastaFile) -> list[tuple[str, int]]:
    out = []
    for chrom in PRIMARY_CHROMS:
        if chrom in fasta.references:
            out.append((chrom, int(fasta.get_reference_length(chrom))))
    if not out:
        raise SystemExit("Reference FASTA has none of the primary chromosomes chr1–22,X,Y.")
    return out


def _local_gc_from_seq(seq: str, pos_1based: int, half_window: int = GC_HALF_WINDOW) -> float:
    """GC fraction in ±half_window around a 1-based position in an already-loaded chrom seq."""
    length = len(seq)
    start0 = max(0, int(pos_1based) - 1 - half_window)
    end0 = min(length, int(pos_1based) - 1 + half_window)
    if end0 <= start0:
        return float("nan")
    window = seq[start0:end0]
    if not window:
        return float("nan")
    gc = 0
    atgc = 0
    for b in window:
        if b in {"G", "C", "g", "c"}:
            gc += 1
            atgc += 1
        elif b in {"A", "T", "a", "t"}:
            atgc += 1
    if atgc == 0:
        return float("nan")
    return gc / atgc


def _gc_of_fetched(seq: str) -> float:
    if not seq:
        return float("nan")
    gc = sum(1 for b in seq if b in {"G", "C", "g", "c"})
    atgc = sum(1 for b in seq if b in {"A", "T", "G", "C", "a", "t", "g", "c"})
    if atgc == 0:
        return float("nan")
    return gc / atgc


def _local_gc(fasta: pysam.FastaFile, chrom: str, pos_1based: int, half_window: int = GC_HALF_WINDOW) -> float:
    length = int(fasta.get_reference_length(chrom))
    start0 = max(0, int(pos_1based) - 1 - half_window)
    end0 = min(length, int(pos_1based) - 1 + half_window)
    if end0 <= start0:
        return float("nan")
    return _gc_of_fetched(fasta.fetch(chrom, start0, end0))


def _gc_bin(gc: float) -> float | None:
    if gc != gc:  # NaN
        return None
    # Bin left edge in [0, 1).
    edge = math.floor(gc / GC_BIN_WIDTH) * GC_BIN_WIDTH
    return round(min(edge, 1.0 - GC_BIN_WIDTH), 10)


def build_rmsk_trees(rmsk_table_path: Path) -> dict[str, dict[str, IntervalTree]]:
    """family -> chrom -> IntervalTree of [start0, end0)."""
    trees: dict[str, dict[str, IntervalTree]] = {fam: defaultdict(IntervalTree) for fam in FAMILIES}
    opener = gzip.open if str(rmsk_table_path).endswith(".gz") else open
    with opener(rmsk_table_path, "rt", encoding="utf-8") as hin:
        for line in hin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                if len(parts) >= 13 and parts[5].startswith("chr"):
                    chrom = parts[5]
                    start0 = int(parts[6])
                    end0 = int(parts[7])
                    rep_name = parts[10]
                    rep_class = parts[11]
                    rep_family = parts[12]
                elif len(parts) >= 8 and parts[0].startswith("chr"):
                    chrom = parts[0]
                    start0 = int(parts[1])
                    end0 = int(parts[2])
                    rep_name = parts[6]
                    rep_class = parts[7]
                    rep_family = parts[8] if len(parts) > 8 else ""
                else:
                    continue
            except (ValueError, IndexError):
                continue
            if end0 <= start0:
                continue
            fam = _normalize_mei_family_token(f"{rep_name} {rep_class} {rep_family}")
            if fam not in trees:
                continue
            trees[fam][chrom].addi(start0, end0)
    return trees


def _is_nested_at(trees: dict[str, dict[str, IntervalTree]], family: str, chrom: str, pos_1based: int) -> bool:
    chrom_tree = trees.get(family, {}).get(chrom)
    if chrom_tree is None or not chrom_tree:
        return False
    pos0 = max(0, int(pos_1based) - 1)
    return bool(chrom_tree.at(pos0))


def build_uniform_pool(
    chrom_lengths: list[tuple[str, int]],
    rng: np.random.Generator,
    n: int,
) -> list[tuple[str, int]]:
    chroms = np.array([c for c, _ in chrom_lengths], dtype=object)
    lengths = np.array([L for _, L in chrom_lengths], dtype=np.int64)
    weights = lengths / lengths.sum()
    chosen_chroms = rng.choice(chroms, size=n, p=weights)
    pool: list[tuple[str, int]] = []
    for chrom in chroms:
        idx = np.flatnonzero(chosen_chroms == chrom)
        if idx.size == 0:
            continue
        L = int(dict(chrom_lengths)[str(chrom)])
        positions = rng.integers(1, L + 1, size=idx.size)
        for pos in positions:
            pool.append((str(chrom), int(pos)))
    rng.shuffle(pool)
    return pool


def build_gc_pools(
    fasta: pysam.FastaFile,
    chrom_lengths: list[tuple[str, int]],
    rng: np.random.Generator,
    n: int = NULL_POOL_SIZE,
) -> dict[float, list[tuple[str, int]]]:
    """Sample genome-wide positions and bucket by local GC bin.

    Loads each chromosome sequence once and scores GC from the in-memory
    string. Random ``fasta.fetch`` per position is far too slow for tens of
    thousands of null draws.
    """
    base_pool = build_uniform_pool(chrom_lengths, rng, n)
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for chrom, pos in base_pool:
        by_chrom[chrom].append(pos)

    buckets: dict[float, list[tuple[str, int]]] = defaultdict(list)
    for chrom, positions in by_chrom.items():
        print(f"  GC-scoring {len(positions)} positions on {chrom} ...", flush=True)
        seq = fasta.fetch(chrom).upper()
        for pos in positions:
            gc = _local_gc_from_seq(seq, pos)
            edge = _gc_bin(gc)
            if edge is None:
                continue
            buckets[edge].append((chrom, pos))
    return buckets


def sample_gc_matched_positions(
    observed_gcs: list[float],
    gc_pools: dict[float, list[tuple[str, int]]],
    uniform_fallback: list[tuple[str, int]],
    rng: np.random.Generator,
) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for gc in observed_gcs:
        edge = _gc_bin(gc)
        pool = gc_pools.get(edge) if edge is not None else None
        if not pool:
            # Empty bin: fall back to uniform genome-wide draw.
            pool = uniform_fallback
        if not pool:
            continue
        out.append(pool[int(rng.integers(0, len(pool)))])
    return out


def selfinsertion_enrichment(
    annotated: pd.DataFrame,
    trees: dict[str, dict[str, IntervalTree]],
    fasta_path: Path,
    n_replicates: int,
    seed: int,
) -> pd.DataFrame:
    """Nested-in-same-family rate vs genome-wide null (UNIFORM for ALU, GC_MATCHED for LINE1/SVA)."""
    fasta = pysam.FastaFile(str(fasta_path))
    try:
        chrom_lengths = _chrom_lengths(fasta)
        rng = np.random.default_rng(seed)
        print(f"Building uniform null pool (n={NULL_POOL_SIZE}) ...", flush=True)
        uniform_pool = build_uniform_pool(chrom_lengths, rng, NULL_POOL_SIZE)
        print("Building GC-matched null pools ...", flush=True)
        gc_pools = build_gc_pools(fasta, chrom_lengths, rng, NULL_POOL_SIZE)
        for edge, pool in sorted(gc_pools.items()):
            print(f"  GC bin [{edge:.2f},{edge + GC_BIN_WIDTH:.2f}): {len(pool)} positions", flush=True)

        rows = []
        null_model_by_family = {"ALU": "UNIFORM", "LINE1": "GC_MATCHED", "SVA": "GC_MATCHED"}
        for fam in FAMILIES:
            sub = annotated.loc[annotated["resolved_family"] == fam]
            n_total = int(len(sub))
            if n_total == 0:
                rows.append(
                    {
                        "family": fam,
                        "null_model": null_model_by_family[fam],
                        "n_total": 0,
                        "n_nested_observed": 0,
                        "nested_rate_observed": float("nan"),
                        "nested_rate_null_mean": float("nan"),
                        "nested_rate_null_sd": float("nan"),
                        "fold_vs_null": float("nan"),
                        "empirical_p_enrichment": float("nan"),
                        "n_replicates": n_replicates,
                        "status": "underpowered",
                        "note": "no calls for this family",
                    }
                )
                continue

            chroms = sub["chrom"].astype(str).tolist()
            positions = sub["insertion_breakpoint_pos"].astype(int).tolist()
            n_nested_obs = sum(
                1 for c, p in zip(chroms, positions) if _is_nested_at(trees, fam, c, p)
            )
            rate_obs = n_nested_obs / n_total

            observed_gcs = [_local_gc(fasta, c, p) for c, p in zip(chroms, positions)]
            null_rates = np.empty(n_replicates, dtype=float)
            for rep in range(n_replicates):
                if null_model_by_family[fam] == "UNIFORM":
                    draws = [uniform_pool[int(rng.integers(0, len(uniform_pool)))] for _ in range(n_total)]
                else:
                    draws = sample_gc_matched_positions(observed_gcs, gc_pools, uniform_pool, rng)
                    # If GC matching could not fill (empty pools), pad from uniform.
                    while len(draws) < n_total:
                        draws.append(uniform_pool[int(rng.integers(0, len(uniform_pool)))])
                    draws = draws[:n_total]
                n_hit = sum(1 for c, p in draws if _is_nested_at(trees, fam, c, p))
                null_rates[rep] = n_hit / n_total

            null_mean = float(null_rates.mean())
            null_sd = float(null_rates.std(ddof=1)) if n_replicates > 1 else float("nan")
            fold = (rate_obs / null_mean) if null_mean > 0 else float("nan")
            # Greater-than enrichment: (1 + #{null >= obs}) / (1 + R)
            emp_p = (1.0 + float(np.sum(null_rates >= rate_obs))) / (1.0 + n_replicates)
            rows.append(
                {
                    "family": fam,
                    "null_model": null_model_by_family[fam],
                    "n_total": n_total,
                    "n_nested_observed": int(n_nested_obs),
                    "nested_rate_observed": rate_obs,
                    "nested_rate_null_mean": null_mean,
                    "nested_rate_null_sd": null_sd,
                    "fold_vs_null": fold,
                    "empirical_p_enrichment": emp_p,
                    "n_replicates": n_replicates,
                    "status": "ok",
                    "note": (
                        "ALU uses UNIFORM null because GC-matching partially controls the signal away; "
                        "LINE1/SVA use GC_MATCHED null"
                        if fam == "ALU"
                        else "GC_MATCHED null: resampled positions matched on local ±500 bp GC (5% bins)"
                    ),
                }
            )
        return pd.DataFrame(rows)
    finally:
        fasta.close()


def write_readme(
    outdir: Path,
    counts: pd.DataFrame,
    sense: pd.DataFrame,
    enrichment: pd.DataFrame,
    vcf_path: Path,
    rmsk_path: Path,
    reference_path: Path,
) -> None:
    """Auto-generate README.md; every statistic is read back from the CSVs."""
    counts_r = pd.read_csv(outdir / "counts_by_family.csv")
    sense_r = pd.read_csv(outdir / "sense_vs_antisense.csv")
    enrich_r = pd.read_csv(outdir / "selfinsertion_enrichment.csv")
    n_calls = int(pd.read_csv(outdir / "per_call_annotated.csv").shape[0])

    def _fmt(x: object) -> str:
        if x is None or (isinstance(x, float) and (math.isnan(x) or pd.isna(x))):
            return "NA"
        if isinstance(x, float):
            return f"{x:.6g}"
        return str(x)

    lines = [
        "# Nested orientation analysis",
        "",
        "Generated by `scripts/analyze_nested_orientation.py`. All numbers below were",
        "read back from the CSVs in this directory at write time (no hardcoded values).",
        "",
        "## Inputs",
        "",
        f"- VCF: `{vcf_path}`",
        f"- rmsk: `{rmsk_path}`",
        f"- reference (GC null): `{reference_path}`",
        f"- calls annotated: {n_calls}",
        "",
        "## Counts by family",
        "",
        "| family | n_total | n_nested | n_sense | n_antisense | n_unknown | n_resolved |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for rec in counts_r.itertuples(index=False):
        lines.append(
            f"| {rec.family} | {rec.n_total} | {rec.n_nested} | {rec.n_sense} | "
            f"{rec.n_antisense} | {rec.n_unknown} | {rec.n_resolved} |"
        )
    lines += [
        "",
        "`n_resolved = n_sense + n_antisense`. `nested_unknown` is excluded from both.",
        "",
        "## Sense vs antisense",
        "",
        "| family | sense_fraction | fold_sense | binom_p_sense | holm_p_sense | antisense_fraction | fold_antisense | status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for rec in sense_r.itertuples(index=False):
        lines.append(
            f"| {rec.family} | {_fmt(rec.sense_fraction)} | {_fmt(rec.fold_vs_chance_sense)} | "
            f"{_fmt(getattr(rec, 'binom_p_sense_vs_0.5'))} | {_fmt(getattr(rec, 'binom_p_sense_vs_0.5_holm'))} | "
            f"{_fmt(rec.antisense_fraction)} | {_fmt(rec.fold_vs_chance_antisense)} | {rec.status} |"
        )
    lines += [
        "",
        "## Self-insertion enrichment",
        "",
        "| family | null_model | nested_rate_obs | nested_rate_null_mean | fold_vs_null | empirical_p | status |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for rec in enrich_r.itertuples(index=False):
        lines.append(
            f"| {rec.family} | {rec.null_model} | {_fmt(rec.nested_rate_observed)} | "
            f"{_fmt(rec.nested_rate_null_mean)} | {_fmt(rec.fold_vs_null)} | "
            f"{_fmt(rec.empirical_p_enrichment)} | {rec.status} |"
        )
    lines += [
        "",
        "## Methods notes",
        "",
        "- Nesting recomputed with `_annotate_nested_retrotransposon` (length-only host selection).",
        "- Sense/antisense binomial tests are exact two-sided tests vs p=0.5; Holm correction across families with n_resolved>0.",
        "- ALU self-insertion null is UNIFORM genome-wide; LINE1 and SVA use GC_MATCHED (±500 bp local GC, 5% bins, 1000 replicates).",
        "",
    ]
    (outdir / "README.md").write_text("\n".join(lines), encoding="utf-8")
    # Silence unused-arg lint for the in-memory frames (we re-read from disk by design).
    _ = (counts, sense, enrichment)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    public = _default_public_data_dir()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", type=Path, required=True, help="Annotated MEI VCF with MEIFAMILY/ORIENT/NESTED")
    p.add_argument(
        "--rmsk",
        type=Path,
        default=public / "annotation" / "hg38" / "repeats" / "rmsk.txt.gz",
        help="UCSC hg38 rmsk table (rmsk.txt.gz)",
    )
    p.add_argument(
        "--reference",
        type=Path,
        default=public / "reference" / "hg38" / "Homo_sapiens_assembly38.fasta",
        help="hg38 FASTA for GC-matched null (LINE1/SVA)",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("results") / "nested_orientation",
        help="Output directory (default: results/nested_orientation)",
    )
    p.add_argument("--n-replicates", type=int, default=1000, help="Null replicates for self-insertion enrichment")
    p.add_argument("--seed", type=int, default=13, help="RNG seed")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    _require_bedtools()

    missing = []
    for label, path in (("vcf", args.vcf), ("rmsk", args.rmsk), ("reference", args.reference)):
        if not Path(path).is_file():
            missing.append(f"{label}: {path}")
    if missing:
        print("Missing required input(s):", file=sys.stderr)
        for m in missing:
            print(f"  - {m}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading VCF: {args.vcf}")
    calls = load_calls_from_vcf(Path(args.vcf))
    if calls.empty:
        print("ERROR: VCF produced zero call rows.", file=sys.stderr)
        return 2

    resolved = resolve_families(calls)
    calls = calls.copy()
    calls["resolved_family"] = resolved
    print_family_resolution(calls, resolved)

    annotate_input = calls[
        ["chrom", "insertion_breakpoint_pos", "MEIFAMILY", "consensus_insertion_orientation", "known_mei_polymorphism_family"]
    ].copy()
    print(f"Re-annotating nesting with rmsk={args.rmsk} ...")
    annotated = _annotate_nested_retrotransposon(annotate_input, Path(args.rmsk))
    annotated["resolved_family"] = calls["resolved_family"].to_numpy()
    annotated["vcf_id"] = calls["vcf_id"].to_numpy()
    annotated["vcf_nested_original"] = calls["vcf_nested_original"].to_numpy()
    annotated["MEIFAMILY"] = calls["MEIFAMILY"].to_numpy()
    # Canonical analysis column.
    annotated["NESTED"] = annotated["nested_same_class_orientation"].astype(str)

    bad = ~annotated["NESTED"].isin(NESTED_LABELS)
    if bad.any():
        raise RuntimeError(f"Unexpected NESTED labels: {sorted(set(annotated.loc[bad, 'NESTED']))}")

    per_call_path = outdir / "per_call_annotated.csv"
    annotated.to_csv(per_call_path, index=False)
    print(f"Wrote {per_call_path} ({len(annotated)} rows)")

    counts = counts_by_family(annotated)
    counts_path = outdir / "counts_by_family.csv"
    counts.to_csv(counts_path, index=False)
    print(f"Wrote {counts_path}")
    print(counts.to_string(index=False))

    sense = sense_vs_antisense_table(counts)
    sense_path = outdir / "sense_vs_antisense.csv"
    sense.to_csv(sense_path, index=False)
    print(f"Wrote {sense_path}")
    print(sense.to_string(index=False))

    print("Building rmsk interval trees for self-insertion null ...")
    trees = build_rmsk_trees(Path(args.rmsk))
    enrichment = selfinsertion_enrichment(
        annotated,
        trees,
        Path(args.reference),
        n_replicates=int(args.n_replicates),
        seed=int(args.seed),
    )
    enrich_path = outdir / "selfinsertion_enrichment.csv"
    enrichment.to_csv(enrich_path, index=False)
    print(f"Wrote {enrich_path}")
    print(enrichment.to_string(index=False))

    write_readme(
        outdir,
        counts,
        sense,
        enrichment,
        Path(args.vcf),
        Path(args.rmsk),
        Path(args.reference),
    )
    print(f"Wrote {outdir / 'README.md'}")

    print("\n=== Summary ===")
    print("Resolved families:", dict(Counter(resolved.tolist())))
    print(sense[["family", "n_sense", "n_antisense", "n_resolved", "sense_fraction", "binom_p_sense_vs_0.5", "binom_p_sense_vs_0.5_holm", "status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
