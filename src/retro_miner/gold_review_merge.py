"""Concatenate per-chromosome gold calls and re-rank them as one genome.

Last step of ``--chr all``, after every chromosome has finished annotation and
while staged alignments are still on disk. Only rows with
``analysis_stage_tier`` gold are kept. Silver and bronze candidates stay in
the per-chromosome tables. The sort is ``_prioritize_mei_candidates``, the
same order a single chromosome already writes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd

from retro_miner.mei_support import _prioritize_mei_candidates
from retro_miner.vcf_export import export_vcf_from_tsv, lookup_reference_build

# Columns the review sort treats as booleans. CSV re-reads them as strings,
# and bool("False") is True, which would scramble stage and catalog flags.
_BOOL_COLUMNS = (
    "gold_stage_pass",
    "silver_stage_pass",
    "known_mei_polymorphism",
    "consensus_tsd_detected",
    "tsd_detected",
    "consensus_poly_at_supported",
    "tsd_poly_at_filter_applied",
    "poly_at_supported",
    "two_sided_support",
)


def _parse_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "t", "yes"})


def _gold_tier_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep calls labeled gold. Silver and bronze are not genome-wide output."""
    if "analysis_stage_tier" not in frame.columns:
        raise ValueError("gold review table is missing analysis_stage_tier")
    tier = frame["analysis_stage_tier"].fillna("").astype(str).str.strip().str.lower()
    return frame.loc[tier.eq("gold")].copy()


def merge_gold_review_tables(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Stack gold-tier rows and sort them like one chromosome."""
    if not frames:
        return pd.DataFrame()
    merged = _gold_tier_rows(pd.concat(list(frames), ignore_index=True))
    if merged.empty:
        return merged
    for col in _BOOL_COLUMNS:
        if col in merged.columns:
            merged[col] = _parse_bool(merged[col])
    ranked = _prioritize_mei_candidates(merged, stage_first=True)
    for col in ranked.columns:
        if ranked[col].dtype == bool:
            ranked[col] = ranked[col].astype(int)
    return ranked


GOLD_REVIEW_NAME = "candidate_loci.mei.gold_review.tsv"
GOLD_REVIEW_VCF_NAME = "candidate_loci.mei.gold_review.vcf"


def annotation_done_marker(chrom: str) -> str:
    """Log line ``run_annotate_mei_support`` prints when a chromosome finishes."""
    return f"stage=annotate-mei-support done region={chrom}"


def incomplete_chromosomes(base_outdir: Path, chroms: Sequence[str]) -> list[str]:
    """Chromosomes whose gold-review table or finished annotation log is missing."""
    missing: list[str] = []
    for chrom in chroms:
        gold = base_outdir / chrom / GOLD_REVIEW_NAME
        log = base_outdir / "logs" / f"{chrom}.log"
        if not gold.is_file() or gold.stat().st_size == 0 or not log.is_file():
            missing.append(chrom)
            continue
        text = log.read_text(errors="replace")
        if annotation_done_marker(chrom) not in text:
            missing.append(chrom)
    return missing


def write_merged_gold_review(inputs: Sequence[Path], output: Path) -> int:
    """Read per-chromosome TSVs, re-rank, and write one genome-wide table."""
    frames: list[pd.DataFrame] = []
    for path in inputs:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        frame = frame.loc[:, ~frame.columns.duplicated()]
        for col in frame.columns:
            frame[col] = frame[col].where(frame[col] != "", other=pd.NA)
        frames.append(frame)
    ranked = merge_gold_review_tables(frames)
    if ranked.empty:
        return 0
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(output, sep="\t", index=False)
    return int(len(ranked))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Re-rank per-chromosome gold reviews as one genome.")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--base-outdir",
        type=Path,
        help="Sample directory. Positional arguments are chromosome names; every one must be finished.",
    )
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args(argv)
    if args.base_outdir is not None:
        chroms = [str(item) for item in args.inputs]
        missing = incomplete_chromosomes(args.base_outdir, chroms)
        if missing:
            joined = ", ".join(missing)
            raise SystemExit(
                f"refusing to aggregate gold review; chromosomes not complete: {joined}"
            )
        paths = [args.base_outdir / chrom / GOLD_REVIEW_NAME for chrom in chroms]
    else:
        paths = [Path(item) for item in args.inputs]
        absent = [str(path) for path in paths if not path.is_file()]
        if absent:
            raise SystemExit(f"gold review table not found: {absent[0]}")
    n_rows = write_merged_gold_review(paths, args.output)
    print(f"[gold-review-merge] rows={n_rows} output={args.output}")
    if n_rows:
        vcf_path = args.output.with_suffix(".vcf")
        build, _fasta = lookup_reference_build(args.output)
        if build:
            print(f"[gold-review-merge] reference_build={build}")
        else:
            print(
                "[gold-review-merge] reference_build not found in pipeline_params.env; "
                "VCF header will omit ##reference"
            )
        n_vcf = export_vcf_from_tsv(args.output, vcf_path, reference_build=build or None)
        print(f"[gold-review-merge] vcf_rows={n_vcf} output={vcf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
