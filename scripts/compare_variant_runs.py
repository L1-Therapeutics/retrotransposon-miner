#!/usr/bin/env python3
"""Compare baseline vs refactored/vectorized annotate runs across ALL candidate loci.

Reads two ``candidate_loci.mei.gold_review.tsv`` files (one per run) and performs
an all-by-all spatial outer join on ``chrom`` + ``consensus_insertion_breakpoint_pos``
within ``--window-bp`` base pairs using a fast per-chromosome 1D range search
(``np.searchsorted`` over sorted breakpoint coordinates). Every candidate row in
each TSV is covered (the full ranked candidate list, not just top-ranked subsets),
and every baseline candidate is matched to its nearest perf candidate within the
window tolerance. Loci entering/leaving the top ``--top-rank`` tier are flagged.

Outputs written under ``--outdir``:

* ``variant_discordance_summary.csv`` - one row per candidate locus (every baseline
  row, plus any perf-only rows) carrying rank shift, gold/high-confidence status
  changes, metric deltas, and 1000G overlap flags.
* ``1000g_overlap_venn.png`` - 2-way Venn diagram of the ``1000g_mei_overlap == True``
  locus sets for the two runs.
* ``rank_decay_1000g_density.png`` - rolling (window 50) 1000G overlap density
  versus review rank for both runs.

Column names are resolved with fallbacks so both canonical names
(e.g. ``total_split_reads``) and pipeline-native names (e.g.
``disease_split_mei_supported_reads``) are accepted.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

CHROM_COL = "chrom"
#: Accepted names for the insertion breakpoint coordinate (new post-refactor names
#: first, then legacy pre-refactor names) so baseline and perf gold-review TSVs
#: with different schemas can be compared.
BP_COL_FALLBACKS = ("consensus_insertion_breakpoint_pos", "insertion_breakpoint_pos", "breakpoint_pos")
FAMILY_COL = "consensus_mei_family"
FAMILY_COL_FALLBACKS = ("consensus_mei_family", "mei_family", "known_mei_polymorphism_family")
ROLLING_WINDOW = 50

#: Canonical metric name -> accepted column names (first present wins).
METRIC_FALLBACKS: dict[str, tuple[str, ...]] = {
    "total_split_reads": (
        "total_split_reads",
        "split_reads_total",
        "mei_split_reads",
        "disease_split_mei_supported_reads",
        "split_mei_supported_reads",
    ),
    "total_discordant_pairs": (
        "total_discordant_pairs",
        "discordant_pairs_total",
        "disease_discordant_mei_supported_reads",
        "discordant_mei_supported_reads",
    ),
    "polyA_tail_len": (
        "polyA_tail_len",
        "polya_tail_len",
        "polya_tail_length",
        "consensus_poly_at_min_bp",
        "poly_at_max_run",
    ),
    "score": (
        "score",
        "insertion_model_score",
        "coherence_score",
        "read_support_heuristic_score",
    ),
}

#: Metrics that should be summed across disease/control when both sides exist.
PAIRED_TOTALS: dict[str, tuple[str, str]] = {
    "total_split_reads": ("disease_split_mei_supported_reads", "control_split_mei_supported_reads"),
    "total_discordant_pairs": (
        "disease_discordant_mei_supported_reads",
        "control_discordant_mei_supported_reads",
    ),
}

TRUTHY = {"true", "1", "yes", "t"}


@dataclass
class Run:
    """Normalized view of one annotate run's gold_review TSV."""

    label: str
    path: Path
    chrom: pd.Series  # str
    bp: pd.Series  # int64, 0 when missing
    rank: pd.Series  # int64, 1-based review rank
    is_gold: pd.Series  # bool
    is_high_conf: pd.Series  # bool (gold + high-confidence breakpoint/assembly)
    is_1000g: pd.Series  # bool
    family: pd.Series  # str
    metrics: dict[str, pd.Series] = field(default_factory=dict)  # float, may be NaN

    def __len__(self) -> int:
        return len(self.chrom)


def _first_col(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in df.columns:
            return name
    return None


def _coerce_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.astype(bool)
    return series.map(str).str.strip().str.lower().isin(TRUTHY)


def _numeric(df: pd.DataFrame, names: tuple[str, ...]) -> pd.Series | None:
    col = _first_col(df, names)
    if col is None:
        return None
    return pd.to_numeric(df[col], errors="coerce")


def _resolve_metric(df: pd.DataFrame, name: str) -> pd.Series | None:
    if name in PAIRED_TOTALS and all(col in df.columns for col in PAIRED_TOTALS[name]):
        left, right = PAIRED_TOTALS[name]
        return (
            pd.to_numeric(df[left], errors="coerce").fillna(0.0)
            + pd.to_numeric(df[right], errors="coerce").fillna(0.0)
        )
    return _numeric(df, METRIC_FALLBACKS[name])


def _rank_series(df: pd.DataFrame) -> pd.Series:
    col = _first_col(df, ("rank", "review_rank"))
    if col is not None:
        ranked = pd.to_numeric(df[col], errors="coerce")
        if ranked.notna().any():
            return ranked.fillna(np.nan).astype("float64")
    return pd.Series(np.arange(1, len(df) + 1, dtype=np.int64), index=df.index)


def _gold_mask(df: pd.DataFrame) -> pd.Series:
    tier_col = _first_col(df, ("analysis_stage_tier", "stage_tier", "tier"))
    if tier_col is not None:
        tier = df[tier_col].fillna("").astype(str).str.strip().str.lower()
        gold = tier.eq("gold")
    else:
        gold = pd.Series(False, index=df.index)
    pass_col = _first_col(df, ("gold_stage_pass", "is_gold", "gold"))
    if pass_col is not None:
        gold = gold | _coerce_bool(df[pass_col])
    return gold


def _high_conf_mask(df: pd.DataFrame, gold: pd.Series) -> pd.Series:
    tier_col = _first_col(
        df,
        (
            "consensus_breakpoint_confidence_tier",
            "breakpoint_confidence_tier",
            "insertion_breakpoint_confidence_tier",
        ),
    )
    if tier_col is not None:
        high = df[tier_col].fillna("").astype(str).str.strip().str.lower().eq("high")
    else:
        high = pd.Series(False, index=df.index)
    conf = _numeric(df, ("assembly_confidence_score",))
    if conf is not None:
        high = high | (conf.fillna(0.0) >= 0.5)
    return gold & high


def _1000g_mask(df: pd.DataFrame) -> pd.Series:
    col = _first_col(df, ("1000g_mei_overlap", "1000_genomes_mei_overlap", "g1k_mei_overlap"))
    if col is not None:
        return _coerce_bool(df[col])
    id_col = _first_col(df, ("g1k_melt_id", "g1k_melt_region_id"))
    if id_col is not None:
        hit = df[id_col].fillna("").astype(str).str.strip().ne("")
    else:
        hit = pd.Series(False, index=df.index)
    src_col = _first_col(df, ("known_mei_polymorphism_source",))
    if src_col is not None:
        hit = hit | df[src_col].fillna("").astype(str).str.contains("melt_1kg", case=False, na=False)
    bool_col = _first_col(df, ("known_mei_polymorphism", "is_known_mei_polymorphism"))
    if bool_col is not None:
        hit = hit | _coerce_bool(df[bool_col])
    return hit


def load_run(path: Path, label: str) -> Run:
    if not path.exists():
        raise SystemExit(f"{label} TSV not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if CHROM_COL not in df.columns:
        raise SystemExit(f"{path}: missing required column: {CHROM_COL}")
    bp_col = _first_col(df, BP_COL_FALLBACKS)
    if bp_col is None:
        raise SystemExit(
            f"{path}: missing breakpoint column (tried: {', '.join(BP_COL_FALLBACKS)})"
        )

    chrom = df[CHROM_COL].fillna("").astype(str).str.strip()
    bp = pd.to_numeric(df[bp_col], errors="coerce").fillna(0).astype("int64")
    gold = _gold_mask(df)
    metrics: dict[str, pd.Series] = {}
    for name in METRIC_FALLBACKS:
        series = _resolve_metric(df, name)
        metrics[name] = series if series is not None else pd.Series(np.nan, index=df.index)
    fam_col = _first_col(df, FAMILY_COL_FALLBACKS)
    family = df[fam_col].fillna("").astype(str).str.strip() if fam_col is not None else pd.Series("", index=df.index)
    return Run(
        label=label,
        path=path,
        chrom=chrom,
        bp=bp,
        rank=_rank_series(df),
        is_gold=gold,
        is_high_conf=_high_conf_mask(df, gold),
        is_1000g=_1000g_mask(df),
        family=family,
        metrics=metrics,
    )


def spatial_pairs(base: Run, perf: Run, window_bp: int) -> pd.DataFrame:
    """All (baseline_row, perf_row) pairs on the same chrom within ``window_bp``.

    Fast 1D range search: per-chromosome ``np.searchsorted`` over the sorted perf
    breakpoint coordinates (O((n+m) log m) total, no exact-coordinate joins).
    """
    b_chrom = base.chrom.to_numpy()
    p_chrom = perf.chrom.to_numpy()
    b_bp = base.bp.to_numpy()
    p_bp = perf.bp.to_numpy()
    b_bp_orig = b_bp
    p_bp_orig = p_bp
    # Rows without a resolved breakpoint (placeholder coordinate <= 0) cannot be
    # matched spatially: they would otherwise explode the window join (every
    # unresolved row within ``window_bp`` of every other on the same chrom).
    # Keep original file-order indices so downstream row lookups stay correct.
    b_keep = np.nonzero(b_bp > 0)[0]
    p_keep = np.nonzero(p_bp > 0)[0]
    b_chrom = b_chrom[b_keep]
    b_bp = b_bp[b_keep]
    p_chrom = p_chrom[p_keep]
    p_bp = p_bp[p_keep]
    if b_chrom.size == 0 or p_chrom.size == 0:
        return pd.DataFrame(columns=["b_row", "p_row", "bp_delta"])
    b_parts: list[np.ndarray] = []
    p_parts: list[np.ndarray] = []
    d_parts: list[np.ndarray] = []
    for chrom_name in np.unique(b_chrom):
        p_idx_all = np.where(p_chrom == chrom_name)[0]
        if p_idx_all.size == 0:
            continue
        order = np.argsort(p_bp[p_idx_all], kind="stable")
        p_bp_sorted = p_bp[p_idx_all][order]
        p_idx_sorted = p_keep[p_idx_all][order]
        b_idx = np.where(b_chrom == chrom_name)[0]
        b_orig = b_keep[b_idx]
        b_bp_c = b_bp[b_idx]
        lo = np.searchsorted(p_bp_sorted, b_bp_c - window_bp, side="left")
        hi = np.searchsorted(p_bp_sorted, b_bp_c + window_bp, side="right")
        counts = hi - lo
        total = int(counts.sum())
        if total == 0:
            continue
        block_start = np.repeat(lo, counts)
        block_offsets = np.arange(total) - np.repeat(
            np.concatenate(([0], np.cumsum(counts)))[:-1], counts
        )
        j_rep = block_start + block_offsets
        b_rep = np.repeat(b_orig, counts)
        p_rep = p_idx_sorted[j_rep]
        b_parts.append(b_rep)
        p_parts.append(p_rep)
        d_parts.append(np.abs(b_bp_orig[b_rep] - p_bp_orig[p_rep]))
    if not b_parts:
        return pd.DataFrame(columns=["b_row", "p_row", "bp_delta"])
    return pd.DataFrame(
        {
            "b_row": np.concatenate(b_parts).astype(np.int64),
            "p_row": np.concatenate(p_parts).astype(np.int64),
            "bp_delta": np.concatenate(d_parts).astype(np.int64),
        }
    )


def _int_key_dict(raw: dict[Any, Any]) -> dict[int, int]:
    return {int(k): int(v) for k, v in raw.items()}


def best_matches(pairs: pd.DataFrame, base: Run, perf: Run) -> tuple[dict[int, int], dict[int, int]]:
    """Closest perf match per baseline row and closest baseline match per perf row.

    Every baseline row that has at least one perf candidate within the window is
    matched to its nearest one (ties broken by lower perf rank).
    """
    if pairs.empty:
        return {}, {}
    p = pairs.copy()
    p["p_rank"] = perf.rank.to_numpy()[p["p_row"].to_numpy()]
    p["b_rank"] = base.rank.to_numpy()[p["b_row"].to_numpy()]
    best_b = _int_key_dict(
        p.sort_values(["bp_delta", "p_rank"], kind="mergesort")
        .drop_duplicates("b_row")
        .set_index("b_row")["p_row"]
        .to_dict()
    )
    best_p = _int_key_dict(
        p.sort_values(["bp_delta", "b_rank"], kind="mergesort")
        .drop_duplicates("p_row")
        .set_index("p_row")["b_row"]
        .to_dict()
    )
    return best_b, best_p


def build_summary(base: Run, perf: Run, pairs: pd.DataFrame, top_rank: int) -> pd.DataFrame:
    """Vectorized locus-by-locus comparison table.

    Row layout: one row per baseline candidate (in file order), followed by one
    row per perf candidate that no baseline candidate matched (in file order).
    """
    best_b, best_p = best_matches(pairs, base, perf)
    n_b = len(base)
    n_p = len(perf)
    b_to_p = np.full(n_b, -1, dtype=np.int64)
    for k, v in best_b.items():
        b_to_p[k] = v
    matched_perf = set(best_b.values())
    perf_only = np.array([j for j in range(n_p) if j not in matched_perf], dtype=np.int64)

    # Per-output-row source indices (-1 = side absent).
    out_b_idx = np.concatenate(
        [np.arange(n_b, dtype=np.int64), np.array([best_p.get(int(j), -1) for j in perf_only], dtype=np.int64)]
    )
    out_p_idx = np.concatenate([b_to_p, perf_only])
    has_b = out_b_idx >= 0
    has_p = out_p_idx >= 0
    b_safe = np.where(has_b, out_b_idx, 0)
    p_safe = np.where(has_p, out_p_idx, 0)

    b_chrom = base.chrom.to_numpy()
    p_chrom = perf.chrom.to_numpy()
    b_bp = base.bp.to_numpy().astype("float64")
    p_bp = perf.bp.to_numpy().astype("float64")
    b_rank = base.rank.to_numpy().astype("float64")
    p_rank = perf.rank.to_numpy().astype("float64")
    b_fam = base.family.to_numpy()
    p_fam = perf.family.to_numpy()

    b_bp_row = np.where(has_b, b_bp[b_safe], np.nan)
    p_bp_row = np.where(has_p, p_bp[p_safe], np.nan)
    b_rank_row = np.where(has_b, b_rank[b_safe], np.nan)
    p_rank_row = np.where(has_p, p_rank[p_safe], np.nan)
    top_b_row = np.where(has_b, b_rank[b_safe] <= top_rank, False)
    top_p_row = np.where(has_p, p_rank[p_safe] <= top_rank, False)
    both = has_b & has_p

    # All-by-all match counts within the window (per locus, per side).
    n_bm = np.zeros(n_b, dtype=np.int64)
    n_pm = np.zeros(n_p, dtype=np.int64)
    if not pairs.empty:
        np.add.at(n_bm, pairs["b_row"].to_numpy(), 1)
        np.add.at(n_pm, pairs["p_row"].to_numpy(), 1)

    data: dict[str, object] = {
        "match_status": np.select([both, has_b], ["shared", "baseline_only"], default="perf_only"),
        "chrom": np.where(has_b, b_chrom[b_safe], p_chrom[p_safe]),
        f"bp_{base.label}": b_bp_row,
        f"bp_{perf.label}": p_bp_row,
        "bp_delta_bp": np.where(both, np.abs(b_bp_row - p_bp_row), np.nan),
        f"rank_{base.label}": b_rank_row,
        f"rank_{perf.label}": p_rank_row,
        "rank_shift": np.where(both, p_rank_row - b_rank_row, np.nan),
        f"top{top_rank}_baseline": top_b_row,
        f"top{top_rank}_perf": top_p_row,
        "top_tier_change": np.select(
            [top_b_row & top_p_row, top_b_row, top_p_row],
            ["stayed_in_top", "exited_top", "entered_top"],
            default="outside_top",
        ),
        "n_perf_matches_within_window": np.concatenate([n_bm, np.zeros(perf_only.size, dtype=np.int64)]),
        "n_baseline_matches_within_window": np.concatenate(
            [np.where(b_to_p >= 0, n_pm[np.where(b_to_p >= 0, b_to_p, 0)], 0), n_pm[perf_only]]
        ),
    }
    for name in METRIC_FALLBACKS:
        b_m = base.metrics[name].to_numpy().astype("float64")
        p_m = perf.metrics[name].to_numpy().astype("float64")
        b_m_row = np.where(has_b, b_m[b_safe], np.nan)
        p_m_row = np.where(has_p, p_m[p_safe], np.nan)
        data[f"{name}_{base.label}"] = b_m_row
        data[f"{name}_{perf.label}"] = p_m_row
        data[f"delta_{name}"] = np.where(both, p_m_row - b_m_row, np.nan)

    fam_b_row = np.where(has_b, b_fam[b_safe], "")
    fam_p_row = np.where(has_p, p_fam[p_safe], "")
    data[f"family_{base.label}"] = fam_b_row
    data[f"family_{perf.label}"] = fam_p_row
    data["family_changed"] = both & (fam_b_row != fam_p_row)

    gold_b = np.where(has_b, base.is_gold.to_numpy()[b_safe], False)
    gold_p = np.where(has_p, perf.is_gold.to_numpy()[p_safe], False)
    hc_b = np.where(has_b, base.is_high_conf.to_numpy()[b_safe], False)
    hc_p = np.where(has_p, perf.is_high_conf.to_numpy()[p_safe], False)
    g1k_b = np.where(has_b, base.is_1000g.to_numpy()[b_safe], False)
    g1k_p = np.where(has_p, perf.is_1000g.to_numpy()[p_safe], False)
    data[f"gold_{base.label}"] = gold_b
    data[f"gold_{perf.label}"] = gold_p
    data[f"high_conf_{base.label}"] = hc_b
    data[f"high_conf_{perf.label}"] = hc_p
    data[f"1000g_{base.label}"] = g1k_b
    data[f"1000g_{perf.label}"] = g1k_p
    data["1000g_overlap_change"] = np.select(
        [g1k_b & g1k_p, g1k_p, g1k_b],
        ["concordant", "gained_in_perf", "lost_in_perf"],
        default="neither",
    )
    gold_tok = np.select([~gold_b & gold_p, gold_b & ~gold_p], ["gained_gold", "lost_gold"], default="")
    hc_tok = np.select([~hc_b & hc_p, hc_b & ~hc_p], ["gained_high_conf", "lost_high_conf"], default="")
    data["status_change"] = np.select(
        [
            (gold_tok != "") & (hc_tok != ""),
            (gold_tok != "") & (hc_tok == ""),
            (gold_tok == "") & (hc_tok != ""),
        ],
        [gold_tok + ";" + hc_tok, gold_tok, hc_tok],
        default="unchanged",
    )

    columns = [
        "match_status",
        "chrom",
        f"bp_{base.label}",
        f"bp_{perf.label}",
        "bp_delta_bp",
        f"rank_{base.label}",
        f"rank_{perf.label}",
        "rank_shift",
        f"top{top_rank}_baseline",
        f"top{top_rank}_perf",
        "top_tier_change",
        "n_perf_matches_within_window",
        "n_baseline_matches_within_window",
    ]
    for name in METRIC_FALLBACKS:
        columns += [f"{name}_{base.label}", f"{name}_{perf.label}", f"delta_{name}"]
    columns += [
        f"family_{base.label}",
        f"family_{perf.label}",
        "family_changed",
        f"gold_{base.label}",
        f"gold_{perf.label}",
        f"high_conf_{base.label}",
        f"high_conf_{perf.label}",
        f"1000g_{base.label}",
        f"1000g_{perf.label}",
        "1000g_overlap_change",
        "status_change",
    ]
    out = pd.DataFrame(data, columns=columns)
    return out.sort_values(f"rank_{base.label}", kind="mergesort").reset_index(drop=True)


def venn_counts(base: Run, perf: Run, pairs: pd.DataFrame) -> dict[str, int]:
    """1000G overlap set sizes using the same spatial matching as the join."""
    best_b, _ = best_matches(pairs, base, perf)
    g_b = base.is_1000g.to_numpy()
    g_p = perf.is_1000g.to_numpy()
    b_idx = [i for i in range(len(base)) if bool(g_b[i])]
    concordant_b = 0
    matched_perf: set[int] = set()
    for i in b_idx:
        j = best_b.get(i)
        if j is not None and bool(g_p[j]):
            concordant_b += 1
            matched_perf.add(j)
    n_perf = int(g_p.sum())
    return {
        "baseline": len(b_idx),
        "perf": n_perf,
        "shared": concordant_b,
        "baseline_only": len(b_idx) - concordant_b,
        "perf_only": n_perf - len(matched_perf),
    }


def _lens_area(r1: float, r2: float, d: float) -> float:
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a1 = r1**2 * math.acos((d**2 + r1**2 - r2**2) / (2.0 * d * r1))
    a2 = r2**2 * math.acos((d**2 + r2**2 - r1**2) / (2.0 * d * r2))
    a3 = 0.5 * math.sqrt(max((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2), 0.0))
    return a1 + a2 - a3


def _solve_venn_geometry(area_a: float, area_b: float, area_ab: float) -> tuple[float, float, float]:
    """Radii and center distance for two area-accurate overlapping circles."""
    r1 = math.sqrt(area_a / math.pi) if area_a > 0 else 0.0
    r2 = math.sqrt(area_b / math.pi) if area_b > 0 else 0.0
    if area_ab <= 0.0 or r1 == 0.0 or r2 == 0.0:
        return r1, r2, r1 + r2
    if area_ab >= min(area_a, area_b) - 1e-9:
        return r1, r2, 0.0  # concentric (subset)
    lo, hi = 0.0, r1 + r2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _lens_area(r1, r2, mid) > area_ab:
            lo = mid
        else:
            hi = mid
    return r1, r2, 0.5 * (lo + hi)


def plot_venn(
    out_png: Path,
    label_a: str,
    label_b: str,
    n_a: int,
    n_b: int,
    n_both: int,
    window_bp: int,
) -> None:
    n_both = max(0, min(int(n_both), int(n_a), int(n_b)))
    only_a = int(n_a) - n_both
    only_b = int(n_b) - n_both
    both = n_both
    r1, r2, d = _solve_venn_geometry(only_a + both, only_b + both, both)
    scale = max(r1, r2, 1.0)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.add_patch(Circle((0.0, 0.0), r1, fill=True, facecolor="#4c8bf5", alpha=0.35, edgecolor="#1f4e96", lw=1.5))
    ax.add_patch(Circle((d, 0.0), r2, fill=True, facecolor="#f5a74c", alpha=0.35, edgecolor="#96561f", lw=1.5))

    def _place(x: float, y: float, text: str, color: str = "#222222") -> None:
        ax.text(x, y, text, ha="center", va="center", fontsize=11, color=color)

    if d > 1e-6:
        _place(-0.5 * d - 0.45 * r1, 0.0, f"only\n{only_a}")
        _place(0.5 * d, 0.0, f"shared\n{both}", color="#111111")
        _place(d + 0.5 * d + 0.45 * r2, 0.0, f"only\n{only_b}")
        _place(-0.5 * d - 0.45 * r1, -r1 - 0.35 * scale, label_a, color="#1f4e96")
        _place(d + 0.5 * d + 0.45 * r2, -r2 - 0.35 * scale, label_b, color="#96561f")
    else:
        outer_only = only_a if r1 >= r2 else only_b
        _place(0.0, 0.75 * max(r1, r2), f"only\n{outer_only}")
        _place(0.0, 0.0, f"shared\n{both}")
        ax.text(
            0.0,
            -max(r1, r2) - 0.35 * scale,
            f"{label_a} / {label_b} (concentric: one set is a subset of the other)",
            ha="center",
            va="top",
            fontsize=9,
        )
    ax.set_xlim(-r1 - 1.6 * scale, d + r2 + 1.6 * scale)
    ax.set_ylim(-max(r1, r2) - 0.9 * scale, max(r1, r2) + 0.7 * scale)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"1000G MEI overlap: {label_a} vs {label_b}\n(total {n_a} / {n_b}, window {window_bp} bp)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_rank_decay(out_png: Path, runs: list[Run], window: int = ROLLING_WINDOW) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    plotted = False
    for run in runs:
        if len(run) == 0:
            continue
        order = np.argsort(run.rank.to_numpy(), kind="stable")
        ranks = np.arange(1, len(run) + 1)
        g = run.is_1000g.to_numpy()[order].astype(float)
        density = pd.Series(g).rolling(window, min_periods=1).mean().to_numpy()
        ax.plot(ranks, density, lw=1.6, label=f"{run.label} (n={len(run)})")
        plotted = True
    ax.set_xlabel("review rank (1 = top candidate)")
    ax.set_ylabel(f"rolling 1000G overlap density (window={window})")
    ax.set_title("1000G overlap density vs review rank")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    if plotted:
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def print_summary(
    summary: pd.DataFrame, venn_counts: dict[str, int], outdir: Path, window_bp: int, top_rank: int
) -> None:
    shared = summary[summary["match_status"] == "shared"]
    rank_shift = shared["rank_shift"].dropna()
    lines: list[tuple[str, str]] = [
        ("baseline candidates", str(len(summary[summary["match_status"].isin(["shared", "baseline_only"])]))),
        ("perf candidates", str(len(summary[summary["match_status"].isin(["shared", "perf_only"])]))),
        ("matched loci (within window)", str(len(shared))),
        ("baseline-only loci", str((summary["match_status"] == "baseline_only").sum())),
        ("perf-only loci", str((summary["match_status"] == "perf_only").sum())),
        ("1000G overlap baseline", str(venn_counts["baseline"])),
        ("1000G overlap perf", str(venn_counts["perf"])),
        ("concordant 1000G overlaps", str(venn_counts["shared"])),
        ("1000G gained in perf", str(venn_counts["perf_only"])),
        ("1000G lost in perf", str(venn_counts["baseline_only"])),
    ]
    if len(rank_shift):
        lines += [
            ("rank shift min", str(int(rank_shift.min()))),
            ("rank shift max", str(int(rank_shift.max()))),
            ("rank shift median", str(float(rank_shift.median()))),
            ("loci with rank shift", str(int((rank_shift != 0).sum()))),
        ]
    entered = int(summary["top_tier_change"].eq("entered_top").sum())
    exited = int(summary["top_tier_change"].eq("exited_top").sum())
    lines.append((f"top-{top_rank} tier entered/exited", f"{entered}/{exited}"))
    for name in METRIC_FALLBACKS:
        delta = shared[f"delta_{name}"].dropna()
        lines.append((f"delta_{name} != 0", str(int((delta != 0).sum())) + f" (of {len(delta)})"))
    if "family_changed" in summary.columns:
        lines.append(("family changed (shared)", str(int(shared["family_changed"].sum()) if len(shared) else 0)))
    for status in ("gold", "high_conf"):
        gained = int(
            summary["status_change"].fillna("").str.contains(f"gained_{status}", regex=False).sum()
        )
        lost = int(summary["status_change"].fillna("").str.contains(f"lost_{status}", regex=False).sum())
        lines.append((f"{status} status gained/lost", f"{gained}/{lost}"))
    width = max(len(k) for k, _ in lines) + 2
    print()
    print("== variant run comparison summary ==")
    for key, value in lines:
        print(f"  {key:<{width}} {value}")
    print(f"  {'output dir':<{width}} {outdir}")
    print(f"  {'window (bp)':<{width}} {window_bp}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare baseline vs refactored/vectorized annotate runs across all candidate loci.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--baseline-tsv", type=Path, required=True, help="Baseline candidate_loci.mei.gold_review.tsv")
    parser.add_argument("--perf-tsv", type=Path, required=True, help="Refactored candidate_loci.mei.gold_review.tsv")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory for summary CSVs and plot PNGs")
    parser.add_argument(
        "--window-bp",
        type=int,
        default=50,
        help="Breakpoint coordinate matching tolerance in base pairs (default: 50)",
    )
    parser.add_argument(
        "--top-rank",
        type=int,
        default=10,
        help="Top-rank tier size used to flag loci entering/leaving the top ranks (default: 10)",
    )
    args = parser.parse_args()
    if args.window_bp < 0:
        parser.error("--window-bp must be >= 0")
    if args.top_rank < 1:
        parser.error("--top-rank must be >= 1")
    args.outdir.mkdir(parents=True, exist_ok=True)

    base = load_run(args.baseline_tsv, "baseline")
    perf = load_run(args.perf_tsv, "perf")

    pairs = spatial_pairs(base, perf, args.window_bp)
    summary = build_summary(base, perf, pairs, args.top_rank)
    summary_path = args.outdir / "variant_discordance_summary.csv"
    summary.to_csv(summary_path, index=False)

    venn = venn_counts(base, perf, pairs)
    venn_path = args.outdir / "1000g_overlap_venn.png"
    plot_venn(venn_path, "baseline", "perf", venn["baseline"], venn["perf"], venn["shared"], args.window_bp)

    decay_path = args.outdir / "rank_decay_1000g_density.png"
    plot_rank_decay(decay_path, [base, perf])

    print_summary(summary, venn, args.outdir, args.window_bp, args.top_rank)
    print(f"wrote {summary_path}")
    print(f"wrote {venn_path}")
    print(f"wrote {decay_path}")


if __name__ == "__main__":
    main()
