from pathlib import Path
import sys

import pandas as pd

#: Canonical score column resolved against both the post-refactor schema and the
#: legacy pre-refactor schema (first column present wins).
SCORE_COL_FALLBACKS = (
    "mei_score",
    "insertion_model_score",
    "read_support_heuristic_score",
    "coherence_score",
    "mei_score_enrichment_ratio",
)


def _first_col(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in df.columns:
            return name
    return None


def _resolve_score(df: pd.DataFrame) -> tuple[str, pd.Series] | None:
    col = _first_col(df, SCORE_COL_FALLBACKS)
    if col is None:
        return None
    return col, pd.to_numeric(df[col], errors="coerce")


def _tier_counts(df: pd.DataFrame) -> dict[str, int]:
    col = _first_col(df, ("analysis_stage_tier", "stage_tier", "tier"))
    if col is None:
        return {}
    return (
        df[col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .value_counts()
        .to_dict()
    )


def _high_conf_count(df: pd.DataFrame) -> int:
    tier_col = _first_col(
        df,
        (
            "consensus_breakpoint_confidence_tier",
            "breakpoint_confidence_tier",
            "insertion_breakpoint_confidence_tier",
        ),
    )
    if tier_col is None:
        return 0
    high = df[tier_col].fillna("").astype(str).str.strip().str.lower().eq("high")
    conf_col = _first_col(df, ("assembly_confidence_score",))
    if conf_col is not None:
        high = high | (pd.to_numeric(df[conf_col], errors="coerce").fillna(0.0) >= 0.5)
    return int(high.sum())


def _bool_from_str(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "t"}


def _concordance_stats(out_md: Path) -> dict[str, object]:
    """Best-effort enrichment from a sibling ``variant_discordance_summary.csv``."""
    summary_csv = out_md.parent / "variant_discordance_summary.csv"
    if not summary_csv.exists():
        return {}
    summary = pd.read_csv(summary_csv, low_memory=False)
    status = summary.get("match_status")
    stats: dict[str, object] = {
        "matched_loci": int((status == "shared").sum()) if status is not None else None,
        "baseline_only": int((status == "baseline_only").sum()) if status is not None else None,
        "perf_only": int((status == "perf_only").sum()) if status is not None else None,
    }
    for col in ("gold_baseline", "gold_perf", "high_conf_baseline", "high_conf_perf"):
        if col in summary.columns:
            stats[col] = int(summary[col].astype(bool).sum())
    if "status_change" in summary.columns:
        sc = summary["status_change"].fillna("").astype(str)
        stats["gained_gold"] = int(sc.str.contains("gained_gold", regex=False).sum())
        stats["lost_gold"] = int(sc.str.contains("lost_gold", regex=False).sum())
        stats["gained_high_conf"] = int(sc.str.contains("gained_high_conf", regex=False).sum())
        stats["lost_high_conf"] = int(sc.str.contains("lost_high_conf", regex=False).sum())
    if "family_changed" in summary.columns:
        stats["family_changed"] = int(summary["family_changed"].astype(bool).sum())
    if "rank_shift" in summary.columns:
        rs = pd.to_numeric(summary["rank_shift"], errors="coerce").dropna()
        if len(rs):
            stats["rank_shift_min"] = int(rs.min())
            stats["rank_shift_max"] = int(rs.max())
            stats["rank_shift_median"] = float(rs.median())
    return stats


def _marker(match: bool) -> str:
    return "OK" if match else "MISMATCH"


def build_report(baseline_path: str, perf_path: str, out_md: str) -> None:
    b_df = pd.read_csv(baseline_path, sep="\t", low_memory=False)
    p_df = pd.read_csv(perf_path, sep="\t", low_memory=False)

    b_score = _resolve_score(b_df)
    p_score = _resolve_score(p_df)
    score_col = b_score[0] if b_score else (p_score[0] if p_score else None)
    b_mean = b_score[1].mean() if b_score else float("nan")
    p_mean = p_score[1].mean() if p_score else float("nan")

    n_match = len(b_df) == len(p_df)
    score_known = b_score is not None and p_score is not None
    score_match = score_known and abs(b_mean - p_mean) < 1e-4

    b_tiers = _tier_counts(b_df)
    p_tiers = _tier_counts(p_df)
    b_hc = _high_conf_count(b_df)
    p_hc = _high_conf_count(p_df)

    lines = [
        "# Pipeline Benchmark Performance Summary",
        "",
        "## Candidate / Annotation Parity",
        "",
        "| Metric | Baseline (`main`) | Optimized (`perf`) | Status |",
        "| :--- | :--- | :--- | :--- |",
        f"| **Total Candidates** | {len(b_df):,} | {len(p_df):,} | {_marker(n_match)} |",
    ]
    if score_col is not None:
        lines.append(
            f"| **Mean {score_col}** | {b_mean:.3f} | {p_mean:.3f} | {_marker(score_match)} |"
        )
    for tier in sorted(set(b_tiers) | set(p_tiers)):
        lines.append(
            f"| **{tier.title()} tier** | {b_tiers.get(tier, 0):,} | {p_tiers.get(tier, 0):,} | {_marker(b_tiers.get(tier, 0) == p_tiers.get(tier, 0))} |"
        )
    lines.append(
        f"| **High-confidence loci** | {b_hc:,} | {p_hc:,} | {_marker(b_hc == p_hc)} |"
    )

    stats = _concordance_stats(Path(out_md))
    if stats and stats.get("matched_loci") is not None:
        matched = stats["matched_loci"]
        bo = stats.get("baseline_only", 0)
        po = stats.get("perf_only", 0)
        concordance = 100.0 * matched / (matched + bo + po) if (matched + bo + po) else 0.0
        lines += [
            "",
            "## Variant Concordance Audit",
            "",
            "| Metric | Value |",
            "| :--- | :--- |",
            f"| Matched loci (50 bp window) | {matched:,} |",
            f"| Baseline-only loci | {bo:,} |",
            f"| Perf-only loci | {po:,} |",
            f"| Per-locus concordance | {concordance:.2f}% |",
        ]
        for key, label in (
            ("gained_gold", "Gold gained in perf"),
            ("lost_gold", "Gold lost in perf"),
            ("gained_high_conf", "High-confidence gained in perf"),
            ("lost_high_conf", "High-confidence lost in perf"),
            ("family_changed", "Family changed (matched)"),
            ("rank_shift_min", "Rank shift min"),
            ("rank_shift_max", "Rank shift max"),
            ("rank_shift_median", "Rank shift median"),
        ):
            if key in stats:
                lines.append(f"| {label} | {stats[key]:,} |")

    lines.append("")
    Path(out_md).write_text("\n".join(lines))


if __name__ == "__main__":
    if len(sys.argv) > 3:
        build_report(sys.argv[1], sys.argv[2], sys.argv[3])