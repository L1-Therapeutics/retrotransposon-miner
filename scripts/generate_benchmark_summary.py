#!/usr/bin/env python3
"""Generate markdown + terminal benchmark reports for pipeline comparison runs.

Compares a *baseline* pipeline run against a *performance* run.  Both runs are
expected to share the same output layout produced by
`scripts/run_candidate_discovery_and_annotation.sh`:

  <dir>/pipeline.log                       # stage wall-clock timings
  <dir>/split_evidence.summary.tsv        # extract-split-evidence summary
  <dir>/candidate_loci.mei.gold_review.tsv  # consolidated gold-review loci

Metrics reported:

  * Stage speedups  - stage-by-stage wall-clock durations (T_baseline / T_perf)
  * Candidate concordance - join gold-review loci within a 200 bp window and
    report the Jaccard index, MATCHED count, and branch-only counts.
  * Truth-set density - overlap rates against g1k_melt / lr_svan hits.

Usage (from repo root, inside the rtm-miner env):

  python3 scripts/generate_benchmark_summary.py \\
    --baseline-dir results/1000g_hg00100_baseline \\
    --perf-dir     results/1000g_hg00100_perf \\
    --outdir       results/benchmark_summary
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.table import Table

MATCH_WINDOW_BP = 200

GOLD_REVIEW_FILENAME = "candidate_loci.mei.gold_review.tsv"
SPLIT_SUMMARY_FILENAME = "split_evidence.summary.tsv"
PIPELINE_LOG_FILENAME = "pipeline.log"

# Pipeline stages recorded in pipeline.log as:
#   [candidate-pipeline] stage=<name> done region=<region> elapsed=<N>s
KNOWN_STAGES = (
    "extract-split-evidence",
    "build-candidate-loci",
    "annotate-mei-support",
)


def _find_gold_review(directory: Path) -> Path | None:
    path = directory / GOLD_REVIEW_FILENAME
    if path.is_file():
        return path
    # Fall back to per-chromosome logs under a nested logs/ layout.
    logs = directory / "logs"
    if logs.is_dir():
        found = sorted(logs.glob(f"*/{GOLD_REVIEW_FILENAME}"))
        if found:
            return found[0]
    return None


# --------------------------------------------------------------------------- #
# Stage timing parsing
# --------------------------------------------------------------------------- #
def _parse_pipeline_log(path: Path) -> dict[str, float]:
    """Return cumulative per-stage wall-clock seconds from a pipeline.log."""
    stage_elapsed: dict[str, float] = defaultdict(float)
    if not path.is_file():
        return dict(stage_elapsed)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            _accumulate_stage_line(line, stage_elapsed)
    return dict(stage_elapsed)


def _accumulate_stage_line(line: str, acc: dict[str, float]) -> None:
    if "stage=" not in line or " done " not in line or "elapsed=" not in line:
        return
    stage = None
    for token in line.split():
        if token.startswith("stage="):
            stage = token.split("=", 1)[1].strip()
        elif token.startswith("elapsed="):
            value = token.split("=", 1)[1].rstrip("s")
            try:
                seconds = float(value)
            except ValueError:
                return
            if stage and stage in KNOWN_STAGES:
                acc[stage] += seconds
            return


def collect_stage_times(directory: Path) -> dict[str, float]:
    """Aggregate per-stage wall times from all pipeline.log files in a dir."""
    agg: dict[str, float] = defaultdict(float)

    log = directory / PIPELINE_LOG_FILENAME
    if log.is_file():
        for stage, seconds in _parse_pipeline_log(log).items():
            agg[stage] += seconds

    logs = directory / "logs"
    if logs.is_dir():
        for log in sorted(logs.glob("*.log")):
            for stage, seconds in _parse_pipeline_log(log).items():
                agg[stage] += seconds

    return dict(agg)


# --------------------------------------------------------------------------- #
# Gold-review loading & normalization
# --------------------------------------------------------------------------- #
def _normalize_gold_review(df: pd.DataFrame) -> pd.DataFrame:
    """Rename location/hit columns to a canonical benchmark schema."""
    out = df.copy()
    if "start_base" not in out.columns and "window_start" in out.columns:
        out["start_base"] = pd.to_numeric(out["window_start"], errors="coerce")
    if "start_base" not in out.columns and "pos" in out.columns:
        out["start_base"] = pd.to_numeric(out["pos"], errors="coerce")
    if "start_base" not in out.columns:
        out["start_base"] = float("nan")

    for src, dst in (("g1k_melt_id", "g1k_melt_hit"), ("lr_svan_id", "lr_svan_hit")):
        if dst not in out.columns:
            if src in out.columns:
                out[dst] = out[src].astype(str).str.strip().ne("") & out[src].notna()
            else:
                out[dst] = False

    if "chrom" not in out.columns:
        out["chrom"] = ""
    return out


def _read_gold_review(directory: Path) -> pd.DataFrame:
    path = _find_gold_review(directory)
    if path is None:
        return _normalize_gold_review(pd.DataFrame())
    df = pd.read_csv(path, sep="\t", low_memory=False)
    return _normalize_gold_review(df)


def _is_truth_hit(row: pd.Series, col: str) -> bool:
    val = row.get(col)
    if val is None or pd.isna(val):
        return False
    if isinstance(val, bool):
        return bool(val)
    if isinstance(val, (int, float)):
        return bool(val != 0)
    return str(val).strip().lower() not in {"", "0", "false", "nan", "none"}


# --------------------------------------------------------------------------- #
# Candidate concordance
# --------------------------------------------------------------------------- #
def concordance_stats(baseline: pd.DataFrame, perf: pd.DataFrame) -> dict[str, object]:
    """Join baseline/perf candidates within a 200 bp same-chrom window."""
    b = baseline[["chrom", "start_base"]].dropna(subset=["start_base"]).copy()
    p = perf[["chrom", "start_base"]].dropna(subset=["start_base"]).copy()
    b["start_base"] = b["start_base"].astype(int)
    p["start_base"] = p["start_base"].astype(int)

    baseline_only: list[tuple[str, int]] = []
    perf_only: list[tuple[str, int]] = []
    matched: list[tuple[str, int, int]] = []

    for chrom in sorted(set(b["chrom"]).union(p["chrom"])):
        b_starts = sorted(b.loc[b["chrom"] == chrom, "start_base"].tolist())
        p_starts = sorted(p.loc[p["chrom"] == chrom, "start_base"].tolist())

        used_b: set[int] = set()
        used_p: set[int] = set()
        pairs: list[tuple[int, int]] = []

        b_index = 0
        for ps in p_starts:
            while b_index < len(b_starts) and b_starts[b_index] < ps - MATCH_WINDOW_BP:
                b_index += 1
            j = b_index
            while j < len(b_starts) and abs(b_starts[j] - ps) <= MATCH_WINDOW_BP:
                if j not in used_b and ps not in used_p:
                    used_b.add(j)
                    used_p.add(ps)
                    pairs.append((b_starts[j], ps))
                    break
                j += 1

        for i, bs in enumerate(b_starts):
            if i not in used_b:
                baseline_only.append((chrom, bs))
        for ps in p_starts:
            if ps not in used_p:
                perf_only.append((chrom, ps))
        matched.extend((chrom, bs, ps) for bs, ps in pairs)

    matched = sorted(set(matched))
    baseline_only = sorted(set(baseline_only))
    perf_only = sorted(set(perf_only))

    matched_by_chrom: dict[str, int] = defaultdict(int)
    for chrom, _bs, _ps in matched:
        matched_by_chrom[chrom] += 1

    return {
        "baseline_candidates": len(b),
        "perf_candidates": len(p),
        "matched": len(matched),
        "baseline_only": len(baseline_only),
        "perf_only": len(perf_only),
        "baseline_only_count": len(baseline_only),
        "perf_only_count": len(perf_only),
        "matched_by_chrom": dict(matched_by_chrom),
    }


def _jaccard(stats: dict[str, object]) -> float:
    union = stats["matched"] + stats["baseline_only"] + stats["perf_only"]
    if union <= 0:
        return 0.0
    return float(stats["matched"]) / float(union)


# --------------------------------------------------------------------------- #
# Truth-set density
# --------------------------------------------------------------------------- #
def truth_density(df: pd.DataFrame) -> dict[str, object]:
    if df.empty:
        return {
            "total": 0,
            "g1k_melt_hit": 0.0,
            "lr_svan_hit": 0.0,
            "g1k_melt_count": 0,
            "lr_svan_count": 0,
        }
    total = float(len(df))
    g1k = sum(1 for _, row in df.iterrows() if _is_truth_hit(row, "g1k_melt_hit"))
    svan = sum(1 for _, row in df.iterrows() if _is_truth_hit(row, "lr_svan_hit"))
    return {
        "total": len(df),
        "g1k_melt_count": g1k,
        "lr_svan_count": svan,
        "g1k_melt_hit": 100.0 * g1k / total,
        "lr_svan_hit": 100.0 * svan / total,
    }


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def _fmt_speedup(ratio: float) -> str:
    return f"{ratio:.2f}x"


def _build_stage_table_stub(baseline: dict[str, float], perf: dict[str, float]) -> list[dict[str, object]]:
    stages = [s for s in KNOWN_STAGES if s in baseline or s in perf]
    rows: list[dict[str, object]] = []
    for stage in stages:
        tb = baseline.get(stage, 0.0)
        tp = perf.get(stage, 0.0)
        speedup = (tb / tp) if tp > 0 else float("inf")
        rows.append(
            {
                "stage": stage,
                "baseline_s": tb,
                "perf_s": tp,
                "speedup": speedup,
            }
        )
    return rows


def render_terminal_report(
    *,
    baseline_dir_name: str,
    perf_dir_name: str,
    stages: list[dict[str, object]],
    concordance: dict[str, object],
    baseline_truth: dict[str, object],
    perf_truth: dict[str, object],
) -> None:
    console = Console()

    console.rule(f"{baseline_dir_name}  vs  {perf_dir_name}")
    console.print()

    stage_table = Table(title="Stage Speedups", show_lines=True)
    stage_table.add_column("Stage")
    stage_table.add_column("Baseline (s)")
    stage_table.add_column("Perf (s)")
    stage_table.add_column("Speedup")
    for row in stages:
        stage_table.add_row(
            str(row["stage"]),
            f"{row['baseline_s']:.1f}",
            f"{row['perf_s']:.1f}",
            _fmt_speedup(float(row["speedup"])),
        )
    console.print(stage_table)
    console.print()

    conc_table = Table(title=f"Candidate Concordance (within {MATCH_WINDOW_BP} bp)", show_lines=True)
    conc_table.add_column("Metric")
    conc_table.add_column("Count")
    conc_rows = [
        ("Baseline candidates", concordance["baseline_candidates"]),
        ("Perf candidates", concordance["perf_candidates"]),
        ("MATCHED", concordance["matched"]),
        (f"Baseline-only ({baseline_dir_name})", concordance["baseline_only"]),
        (f"Perf-only ({perf_dir_name})", concordance["perf_only"]),
        ("", ""),
        ("Jaccard index", f"{_jaccard(concordance):.3f}"),
    ]
    for label, value in conc_rows:
        conc_table.add_row(str(label), str(value))
    console.print(conc_table)
    console.print()

    truth_table = Table(title="Truth-Set Density", show_lines=True)
    truth_table.add_column("Run")
    truth_table.add_column("Total")
    truth_table.add_column("g1k_melt hit")
    truth_table.add_column("lr_svan hit")
    truth_table.add_column("g1k_melt count")
    truth_table.add_column("lr_svan count")
    truth_table.add_row(
        baseline_dir_name,
        str(baseline_truth["total"]),
        f"{baseline_truth['g1k_melt_hit']:.1f}%",
        f"{baseline_truth['lr_svan_hit']:.1f}%",
        str(baseline_truth["g1k_melt_count"]),
        str(baseline_truth["lr_svan_count"]),
    )
    truth_table.add_row(
        perf_dir_name,
        str(perf_truth["total"]),
        f"{perf_truth['g1k_melt_hit']:.1f}%",
        f"{perf_truth['lr_svan_hit']:.1f}%",
        str(perf_truth["g1k_melt_count"]),
        str(perf_truth["lr_svan_count"]),
    )
    console.print(truth_table)


def write_markdown_report(
    path: Path,
    *,
    baseline_dir_name: str,
    perf_dir_name: str,
    stages: list[dict[str, object]],
    concordance: dict[str, object],
    baseline_truth: dict[str, object],
    perf_truth: dict[str, object],
) -> None:
    lines: list[str] = []
    lines.append(f"# Benchmark Summary: `{baseline_dir_name}` vs `{perf_dir_name}`")
    lines.append("")

    lines.append("## Stage Speedups")
    lines.append("")
    lines.append("| Stage | Baseline (s) | Perf (s) | Speedup |")
    lines.append("|---|---:|---:|---:|")
    for row in stages:
        lines.append(
            f"| {row['stage']} | {row['baseline_s']:.1f} | {row['perf_s']:.1f} | "
            f"{_fmt_speedup(float(row['speedup']))} |"
        )
    lines.append("")

    lines.append(f"## Candidate Concordance (within {MATCH_WINDOW_BP} bp)")
    lines.append("")
    lines.append(f"- Baseline candidates: **{concordance['baseline_candidates']}**")
    lines.append(f"- Perf candidates: **{concordance['perf_candidates']}**")
    lines.append(f"- MATCHED: **{concordance['matched']}**")
    lines.append(f"- Baseline-only (`{baseline_dir_name}`): **{concordance['baseline_only']}**")
    lines.append(f"- Perf-only (`{perf_dir_name}`): **{concordance['perf_only']}**")
    lines.append(f"- **Jaccard index: {_jaccard(concordance):.3f}**")
    lines.append("")
    lines.append("### Matched per chromosome")
    lines.append("")
    lines.append("| Chromosome | MATCHED |")
    lines.append("|---|---:|")
    for chrom, count in sorted(concordance["matched_by_chrom"].items()):
        lines.append(f"| {chrom} | {count} |")
    lines.append("")

    lines.append("## Truth-Set Density")
    lines.append("")
    lines.append("| Run | Total | g1k_melt hit | lr_svan hit | g1k_melt count | lr_svan count |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    lines.append(
        f"| `{baseline_dir_name}` | {baseline_truth['total']} | "
        f"{baseline_truth['g1k_melt_hit']:.1f}% | {baseline_truth['lr_svan_hit']:.1f}% | "
        f"{baseline_truth['g1k_melt_count']} | {baseline_truth['lr_svan_count']} |"
    )
    lines.append(
        f"| `{perf_dir_name}` | {perf_truth['total']} | "
        f"{perf_truth['g1k_melt_hit']:.1f}% | {perf_truth['lr_svan_hit']:.1f}% | "
        f"{perf_truth['g1k_melt_count']} | {perf_truth['lr_svan_count']} |"
    )
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True, help="Directory with baseline run results.")
    parser.add_argument("--perf-dir", type=Path, required=True, help="Directory with performance run results.")
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory to store generated reports and figures.",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()

    baseline_stage = collect_stage_times(args.baseline_dir)
    perf_stage = collect_stage_times(args.perf_dir)

    baseline_gold = _read_gold_review(args.baseline_dir)
    perf_gold = _read_gold_review(args.perf_dir)

    concordance = concordance_stats(baseline_gold, perf_gold)
    baseline_truth = truth_density(baseline_gold)
    perf_truth = truth_density(perf_gold)

    stages = _build_stage_table_stub(baseline_stage, perf_stage)

    render_terminal_report(
        baseline_dir_name=args.baseline_dir.name,
        perf_dir_name=args.perf_dir.name,
        stages=stages,
        concordance=concordance,
        baseline_truth=baseline_truth,
        perf_truth=perf_truth,
    )

    report_path = args.outdir / "BENCHMARK_REPORT.md"
    write_markdown_report(
        report_path,
        baseline_dir_name=args.baseline_dir.name,
        perf_dir_name=args.perf_dir.name,
        stages=stages,
        concordance=concordance,
        baseline_truth=baseline_truth,
        perf_truth=perf_truth,
    )
    print(f"\nWrote {report_path}")


if __name__ == "__main__":
    main()
