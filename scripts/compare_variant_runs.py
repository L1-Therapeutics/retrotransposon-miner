#!/usr/bin/env python3
"""Compare baseline vs refactored/vectorized annotate runs across ALL candidate loci.

Reads two ``candidate_loci.mei.gold_review.tsv`` files (one per run) and performs
an all-by-all spatial outer join on ``chrom`` + ``consensus_insertion_breakpoint_pos``
within ``--window-bp`` base pairs. Every candidate row in each TSV is covered (the
full ranked candidate list, not just top-ranked subsets).

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

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

CHROM_COL = "chrom"
BP_COL = "consensus_insertion_breakpoint_pos"
FAMILY_COL = "consensus_mei_family"
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
    tier_col = _first_col(df, ("consensus_breakpoint_confidence_tier", "breakpoint_confidence_tier"))
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
    return hit


def load_run(path: Path, label: str) -> Run:
    if not path.exists():
        raise SystemExit(f"{label} TSV not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in (CHROM_COL, BP_COL) if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing required column(s): {', '.join(missing)}")

    chrom = df[CHROM_COL].fillna("").astype(str).str.strip()
    bp = pd.to_numeric(df[BP_COL], errors="coerce").fillna(0).astype("int64")
    gold = _gold_mask(df)
    metrics: dict[str, pd.Series] = {}
    for name in METRIC_FALLBACKS:
        series = _resolve_metric(df, name)
        metrics[name] = series if series is not None else pd.Series(np.nan, index=df.index)
    fam_col = _first_col(df, (FAMILY_COL, "mei_family"))
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
    """All (baseline_row, perf_row) pairs on the same chrom within ``window_bp``."""
    b_chrom = base.chrom.to_numpy()
    p_chrom = perf.chrom.to_numpy()
    b_bp = base.bp.to_numpy()
    p_bp = perf.bp.to_numpy()
    rows: list[tuple[int, int, int]] = []
    for chrom_name in np.unique(b_chrom):
        p_idx_all = np.where(p_chrom == chrom_name)[0]
        if p_idx_all.size == 0:
            continue
        order = np.argsort(p_bp[p_idx_all], kind="stable")
        p_bp_sorted = p_bp[p_idx_all][order]
        p_idx_sorted = p_idx_all[order]
        b_idx = np.where(b_chrom == chrom_name)[0]
        b_bp_c = b_bp[b_idx]
        lo = np.searchsorted(p_bp_sorted, b_bp_c - window_bp, side="left")
        hi = np.searchsorted(p_bp_sorted, b_bp_c + window_bp, side="right")
        for k in range(b_idx.size):
            for j in range(int(lo[k]), int(hi[k])):
                j_pos = int(p_idx_sorted[j])
                rows.append((int(b_idx[k]), j_pos, abs(int(b_bp_c[k]) - int(p_bp[j_pos]))))
    return pd.DataFrame(rows, columns=["b_row", "p_row", "bp_delta"])


def best_matches(pairs: pd.DataFrame, base: Run, perf: Run) -> tuple[dict[int, int], dict[int, int]]:
    """Closest perf match per baseline row and closest baseline match per perf row."""
    if pairs.empty:
        return {}, {}
    p = pairs.copy()
    p["p_rank"] = perf.rank.to_numpy()[p["p_row"].to_numpy()]
    p["b_rank"] = base.rank.to_numpy()[p["b_row"].to_numpy()]
    best_b = (
        p.sort_values(["bp_delta", "p_rank"], kind="mergesort")
        .drop_duplicates("b_row")
        .set_index("b_row")["p_row"]
        .to_dict()
    )
    best_p = (
        p.sort_values(["bp_delta", "b_rank"], kind="mergesort")
        .drop_duplicates("p_row")
        .set_index("p_row")["b_row"]
        .to_dict()
    )
    return best_b, best_p


def _status_change(label: str, base_flag: bool, perf_flag: bool) -> list[str]:
    if base_flag == perf_flag:
        return []
    return [f"{'gained' if perf_flag else 'lost'}_{label}"]


def build_summary(base: Run, perf: Run, pairs: pd.DataFrame) -> pd.DataFrame:
    best_b, best_p = best_matches(pairs, base, perf)
    matched_perf_rows = set(best_b.values())
    n_b_matches = pairs.groupby("b_row").size().to_dict() if not pairs.empty else {}
    n_p_matches = pairs.groupby("p_row").size().to_dict() if not pairs.empty else {}

    metric_names = list(METRIC_FALLBACKS)
    columns = [
        "match_status",
        "chrom",
        f"bp_{base.label}",
        f"bp_{perf.label}",
        "bp_delta_bp",
        f"rank_{base.label}",
        f"rank_{perf.label}",
        "rank_shift",
        "n_perf_matches_within_window",
        "n_baseline_matches_within_window",
    ]
    for name in metric_names:
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

    rows: list[dict[str, object]] = []
    b_arr = {
        "chrom": base.chrom.to_numpy(),
        "bp": base.bp.to_numpy(),
        "rank": base.rank.to_numpy(),
        "gold": base.is_gold.to_numpy(),
        "high_conf": base.is_high_conf.to_numpy(),
        "g1k": base.is_1000g.to_numpy(),
        "family": base.family.to_numpy(),
        **{f"m_{m}": base.metrics[m].to_numpy() for m in metric_names},
    }
    p_arr = {
        "chrom": perf.chrom.to_numpy(),
        "bp": perf.bp.to_numpy(),
        "rank": perf.rank.to_numpy(),
        "gold": perf.is_gold.to_numpy(),
        "high_conf": perf.is_high_conf.to_numpy(),
        "g1k": perf.is_1000g.to_numpy(),
        "family": perf.family.to_numpy(),
        **{f"m_{m}": perf.metrics[m].to_numpy() for m in metric_names},
    }

    def _row(
        *,
        status: str,
        b: int | None,
        p: int | None,
        delta_bp: int | None,
        n_b: int,
        n_p: int,
    ) -> dict[str, object]:
        b_rec = b_arr if b is not None else None
        p_rec = p_arr if p is not None else None
        rec: dict[str, object] = {
            "match_status": status,
            "chrom": (b_rec or p_rec)["chrom"][b if b is not None else p],
            f"bp_{base.label}": b_arr["bp"][b] if b is not None else np.nan,
            f"bp_{perf.label}": p_arr["bp"][p] if p is not None else np.nan,
            "bp_delta_bp": delta_bp,
            f"rank_{base.label}": b_arr["rank"][b] if b is not None else np.nan,
            f"rank_{perf.label}": p_arr["rank"][p] if p is not None else np.nan,
            "rank_shift": (p_arr["rank"][p] - b_arr["rank"][b]) if (b is not None and p is not None) else np.nan,
            "n_perf_matches_within_window": n_b,
            "n_baseline_matches_within_window": n_p,
        }
        for name in metric_names:
            b_val = float(b_arr[f"m_{name}"][b]) if b is not None else np.nan
            p_val = float(p_arr[f"m_{name}"][p]) if p is not None else np.nan
            rec[f"{name}_{base.label}"] = b_val
            rec[f"{name}_{perf.label}"] = p_val
            rec[f"delta_{name}"] = (p_val - b_val) if (b is not None and p is not None) else np.nan
        fam_b = b_arr["family"][b] if b is not None else ""
        fam_p = p_arr["family"][p] if p is not None else ""
        rec[f"family_{base.label}"] = fam_b
        rec[f"family_{perf.label}"] = fam_p
        rec["family_changed"] = (b is not None and p is not None) and bool(fam_b != fam_p)
        gold_b = bool(b_arr["gold"][b]) if b is not None else False
        gold_p = bool(p_arr["gold"][p]) if p is not None else False
        hc_b = bool(b_arr["high_conf"][b]) if b is not None else False
        hc_p = bool(p_arr["high_conf"][p]) if p is not None else False
        g_b = bool(b_arr["g1k"][b]) if b is not None else False
        g_p = bool(p_arr["g1k"][p]) if p is not None else False
        rec[f"gold_{base.label}"] = gold_b
        rec[f"gold_{perf.label}"] = gold_p
        rec[f"high_conf_{base.label}"] = hc_b
        rec[f"high_conf_{perf.label}"] = hc_p
        rec[f"1000g_{base.label}"] = g_b
        rec[f"1000g_{perf.label}"] = g_p
        if g_b and g_p:
            rec["1000g_overlap_change"] = "concordant"
        elif g_p:
            rec["1000g_overlap_change"] = "gained_in_perf"
        elif g_b:
            rec["1000g_overlap_change"] = "lost_in_perf"
        else:
            rec["1000g_overlap_change"] = "neither"
        changes = _status_change("gold", gold_b, gold_p) + _status_change("high_conf", hc_b, hc_p)
        rec["status_change"] = ";".join(changes) if changes else "unchanged"
        return rec

    for i in range(len(base)):
        j = best_b.get(i)
        if j is not None:
            rows.append(
                _row(
                    status="shared",
                    b=i,
                    p=j,
                    delta_bp=abs(int(base.bp.iloc[i]) - int(perf.bp.iloc[j])),
                    n_b=int(n_b_matches.get(i, 1)),
                    n_p=int(n_p_matches.get(j, 1)),
                )
            )
        else:
            rows.append(_row(status="baseline_only", b=i, p=None, delta_bp=None, n_b=0, n_p=0))
    for j in range(len(perf)):
        if j in matched_perf_rows:
            continue
        i = best_p.get(j)
        rows.append(
            _row(
                status="perf_only",
                b=i,
                p=j,
                delta_bp=(abs(int(base.bp.iloc[i]) - int(perf.bp.iloc[j])) if i is not None else None),
                n_b=0,
                n_p=int(n_p_matches.get(j, 0)),
            )
        )
    out = pd.DataFrame(rows, columns=columns)
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


def print_summary(summary: pd.DataFrame, venn_counts: dict[str, int], outdir: Path, window_bp: int) -> None:
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
    args = parser.parse_args()
    if args.window_bp < 0:
        parser.error("--window-bp must be >= 0")
    args.outdir.mkdir(parents=True, exist_ok=True)

    base = load_run(args.baseline_tsv, "baseline")
    perf = load_run(args.perf_tsv, "perf")

    pairs = spatial_pairs(base, perf, args.window_bp)
    summary = build_summary(base, perf, pairs)
    summary_path = args.outdir / "variant_discordance_summary.csv"
    summary.to_csv(summary_path, index=False)

    venn = venn_counts(base, perf, pairs)
    venn_path = args.outdir / "1000g_overlap_venn.png"
    plot_venn(venn_path, "baseline", "perf", venn["baseline"], venn["perf"], venn["shared"], args.window_bp)

    decay_path = args.outdir / "rank_decay_1000g_density.png"
    plot_rank_decay(decay_path, [base, perf])

    print_summary(summary, venn, args.outdir, args.window_bp)
    print(f"wrote {summary_path}")
    print(f"wrote {venn_path}")
    print(f"wrote {decay_path}")


if __name__ == "__main__":
    main()
