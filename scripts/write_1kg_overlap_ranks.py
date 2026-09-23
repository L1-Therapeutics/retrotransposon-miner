#!/usr/bin/env python3
"""Rank gold-review calls against HG00100 1000 Genomes MELT and ONT sites.

Writes a rank-sorted overlap table. Unrecovered catalog sites are plotted by
default: an IGV snapshot and a read-architecture figure centered on the catalog
coordinate, using a candidate window when one already covers that base.

Example (HG00100 chr17 on the VM):

  PYTHONPATH=src python scripts/write_1kg_overlap_ranks.py \\
    --calls ~/retrotransposon-workdir/results/hg00100_germline_chr17_span/candidate_loci.mei.gold_review.tsv \\
    --chrom chr17 \\
    --sample HG00100 \\
    --melt-vcf $DATA/polymorphism/hg38/1kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz \\
    --ont-svim $DATA/polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.bcf \\
    --ont-svan $DATA/polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.noGt.SVAN_1.3.bcf \\
    --reference-fasta $DATA/reference/hg38/Homo_sapiens_assembly38.fasta \\
    --disease-bam ~/retrotransposon-workdir/data/bam_stage/HG00100.final.cram \\
    --outdir ~/retrotransposon-workdir/results/hg00100_germline_chr17_span

MELT and ONT sites inside the caller junk mask are dropped by default
(segdup, low mappability, gaps, ENCODE blacklist). The bed is
$RTM_PUBLIC_DATA_DIR/annotation/hg38/junk/junk_exclusion_merged.bed.
Pass --no-junk-filter to keep those catalog sites.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from retro_miner.callset_eval import (
    Variant,
    filter_variants_to_regions,
    load_melt_sample_meis,
    load_ont_sample_meis,
    match_variants,
    merge_unique_events,
    overlap_metrics,
    resolve_exclude_beds,
)

PAD_BP = 300
SYNTHETIC_FLANK_BP = 1500


def _gold_frame(path: Path, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    gold = df.loc[df["analysis_stage_tier"].astype(str).str.lower().eq("gold")].copy()
    if chrom:
        gold = gold.loc[gold["chrom"].astype(str).eq(chrom)]
    gold = gold.reset_index(drop=True)
    gold["rank"] = gold.index + 1
    gold["bp"] = pd.to_numeric(gold["consensus_insertion_breakpoint_pos"], errors="coerce")
    return gold


def _gold_calls(gold: pd.DataFrame, chrom: str) -> list[Variant]:
    calls: list[Variant] = []
    for rec in gold.to_dict(orient="records"):
        if pd.isna(rec["bp"]) or int(rec["bp"]) <= 0:
            continue
        pos = int(rec["bp"])
        calls.append(
            Variant(
                chrom=str(rec.get("chrom") or chrom),
                pos=pos,
                end=pos,
                variant_id=f"rtm_{pos}_{int(rec['rank'])}",
                family=str(rec.get("consensus_mei_family") or ""),
                source="rtm",
            )
        )
    return calls


def _truth_label(sources: dict) -> str:
    has_melt = any("melt" in source for source in sources)
    has_ont = any(token in source for source in sources for token in ("ont", "svim", "long_read"))
    if has_melt and has_ont:
        return "MELT+ONT"
    if has_melt:
        return "MELT"
    if has_ont:
        return "ONT"
    return ",".join(sorted(sources))


def overlap_table(gold: pd.DataFrame, events: list[Variant], calls: list[Variant]) -> pd.DataFrame:
    hits = match_variants(events, calls, pad_bp=PAD_BP, require_family=False)
    rows = []
    for event in events:
        sources = event.extra.get("source_ids", {})
        matched = hits.get(event.key) or []
        hit = min(matched, key=lambda variant: abs(variant.pos - event.pos)) if matched else None
        rank = "-"
        rtm_pos = "-"
        if hit is not None:
            _prefix, rtm_pos, rank = hit.variant_id.split("_")
        rows.append(
            {
                "rank": rank,
                "pos": event.pos,
                "rtm_pos": rtm_pos,
                "family": event.family,
                "truth": _truth_label(sources),
                "ids": ", ".join(variant_id for ids in sources.values() for variant_id in ids),
                "recovered": "yes" if hit is not None else "no",
                "chrom": event.chrom,
            }
        )
    table = pd.DataFrame(rows)
    table["_rank"] = pd.to_numeric(table["rank"], errors="coerce")
    table = table.sort_values(["_rank", "pos"], na_position="last").drop(columns="_rank")
    return table.reset_index(drop=True)


def _plot_stem(chrom: str, pos: int, family: str, variant_id: str) -> str:
    safe_id = re.sub(r"[^A-Za-z0-9._+-]+", "_", variant_id).strip("_") or "na"
    fam = re.sub(r"[^A-Za-z0-9._+-]+", "_", family) or "NA"
    return f"{chrom}_{pos}_{fam}_{safe_id}"


def _containing_or_nearest(review: pd.DataFrame, chrom: str, pos: int) -> tuple[pd.Series, int, bool]:
    sub = review.loc[review["chrom"].astype(str).eq(chrom)].copy()
    if sub.empty:
        raise ValueError(f"No candidate rows on {chrom}")
    ws = pd.to_numeric(sub["window_start"], errors="coerce")
    we = pd.to_numeric(sub["window_end"], errors="coerce")
    bp = pd.to_numeric(sub["consensus_insertion_breakpoint_pos"], errors="coerce")
    contains = (ws <= pos) & (we >= pos)
    if bool(contains.any()):
        cand = sub.loc[contains].copy()
        cand["_dist"] = (pd.to_numeric(cand["consensus_insertion_breakpoint_pos"], errors="coerce") - pos).abs()
        cand["_span"] = (
            pd.to_numeric(cand["window_end"], errors="coerce")
            - pd.to_numeric(cand["window_start"], errors="coerce")
        )
        row = cand.sort_values(["_dist", "_span"], na_position="last").iloc[0]
        dist = int(row["_dist"]) if pd.notna(row["_dist"]) else -1
        return row.drop(labels=["_dist", "_span"]), dist, True
    sub["_dist"] = (bp - pos).abs()
    row = sub.sort_values("_dist", na_position="last").iloc[0]
    dist = int(row["_dist"]) if pd.notna(row["_dist"]) else -1
    return row.drop(labels=["_dist"]), dist, False


def _view_row(source: pd.Series, *, chrom: str, pos: int, contains: bool) -> pd.Series:
    row = source.copy()
    if contains:
        return row
    start = max(1, pos - SYNTHETIC_FLANK_BP)
    end = pos + SYNTHETIC_FLANK_BP
    row["chrom"] = chrom
    row["window_start"] = start
    row["window_end"] = end
    row["discovery_window_start"] = start
    row["discovery_window_end"] = end
    row["consensus_insertion_breakpoint_pos"] = pos
    row["insertion_breakpoint_pos"] = pos
    for col in row.index:
        if str(col).endswith("_supporting_reads"):
            row[col] = ""
    return row


def _reuse_chromosome_bam(snapshot_dir: Path, calls_path: Path) -> None:
    """Point IGV at the chromosome BAM from the gold-review run when it exists."""
    existing = calls_path.parent / "candidate_loci.mei.gold_review.igv" / "igv_disease.bam"
    bai = Path(str(existing) + ".bai")
    if not existing.exists() or not bai.exists():
        return
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    dest = snapshot_dir / "igv_disease.bam"
    dest_bai = Path(str(dest) + ".bai")
    if not dest.exists():
        dest.symlink_to(existing.resolve())
    if not dest_bai.exists():
        dest_bai.symlink_to(bai.resolve())


def plot_unrecovered(
    review: pd.DataFrame,
    missed: pd.DataFrame,
    *,
    calls_path: Path,
    outdir: Path,
    reference_fasta: Path,
    disease_bam: Path,
    control_bam: Path | None,
) -> Path:
    from retro_miner.igv_plots import generate_gold_review_igv_plots
    from retro_miner.read_architecture import ReadArchitectureCache, plot_locus_architecture

    igv_dir = outdir / "missed_1kg" / "igv"
    arch_dir = outdir / "missed_1kg" / "read_architecture"
    igv_dir.mkdir(parents=True, exist_ok=True)
    arch_dir.mkdir(parents=True, exist_ok=True)

    view_rows: list[pd.Series] = []
    manifest_rows: list[dict[str, object]] = []
    for rec in missed.to_dict(orient="records"):
        chrom = str(rec["chrom"])
        pos = int(rec["pos"])
        ids = str(rec["ids"])
        primary_id = ids.split(",")[0].strip() or "na"
        source, dist, contains = _containing_or_nearest(review, chrom, pos)
        view = _view_row(source, chrom=chrom, pos=pos, contains=contains)
        view_rows.append(view)
        stem = _plot_stem(chrom, pos, str(rec["family"]), primary_id)
        manifest_rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "family": rec["family"],
                "truth": rec["truth"],
                "ids": ids,
                "candidate_tier": source.get("analysis_stage_tier", ""),
                "candidate_bp": source.get("consensus_insertion_breakpoint_pos", ""),
                "candidate_distance_bp": dist,
                "window_contains_catalog": int(contains),
                "window_start": int(view["window_start"]),
                "window_end": int(view["window_end"]),
                "igv_png": str(igv_dir / f"{stem}.igv.png"),
                "readarch_png": str(arch_dir / f"{stem}.readarch.png"),
                "stem": stem,
            }
        )
        print(
            f"miss {chrom}:{pos} {rec['family']} {primary_id} "
            f"tier={source.get('analysis_stage_tier', '')} "
            f"bp={source.get('consensus_insertion_breakpoint_pos', '')} "
            f"dist={dist} contains={int(contains)}",
            flush=True,
        )

    views = pd.DataFrame(view_rows).reset_index(drop=True)
    _reuse_chromosome_bam(igv_dir, calls_path)
    generate_gold_review_igv_plots(
        views,
        reference_fasta=reference_fasta,
        disease_bam=disease_bam,
        control_bam=control_bam or disease_bam,
        snapshot_dir=igv_dir,
        top_n=0,
        gold_only=False,
    )
    rank_pngs = sorted(igv_dir.glob("rank*.png"))
    for manifest, rank_png in zip(manifest_rows, rank_pngs):
        dest = Path(str(manifest["igv_png"]))
        dest.write_bytes(rank_png.read_bytes())

    cache = ReadArchitectureCache.from_paths(calls_path, load_split_evidence=True)
    for manifest, view in zip(manifest_rows, view_rows):
        plot_locus_architecture(
            chrom=str(view["chrom"]),
            pos=int(pd.to_numeric(view["consensus_insertion_breakpoint_pos"], errors="coerce")),
            sample="auto",
            out_png=Path(str(manifest["readarch_png"])),
            cache=cache,
            row=view,
        )

    manifest_path = outdir / "missed_1kg" / "missed_1kg_manifest.tsv"
    pd.DataFrame(manifest_rows).drop(columns=["stem"]).to_csv(manifest_path, sep="\t", index=False)
    print(f"wrote {manifest_path}", flush=True)
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True, help="candidate_loci.mei.gold_review.tsv")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--sample", default="HG00100")
    parser.add_argument("--melt-vcf", type=Path, required=True)
    parser.add_argument("--ont-svim", type=Path, required=True)
    parser.add_argument("--ont-svan", type=Path, required=True)
    parser.add_argument(
        "--exclude-bed",
        type=Path,
        action="append",
        default=None,
        help=(
            "Drop MELT and ONT sites inside this BED before overlap (repeatable). "
            "Default: junk_exclusion_merged.bed for hg38."
        ),
    )
    parser.add_argument(
        "--no-junk-filter",
        action="store_true",
        help="Do not apply the default caller junk mask to MELT and ONT.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, default=None)
    parser.add_argument("--disease-bam", type=Path, default=None)
    parser.add_argument("--control-bam", type=Path, default=None)
    parser.add_argument(
        "--no-missed-plots",
        action="store_true",
        help="Write the overlap table only",
    )
    args = parser.parse_args()

    review = pd.read_csv(args.calls, sep="\t", low_memory=False)
    gold = _gold_frame(args.calls, args.chrom)
    calls = _gold_calls(gold, args.chrom)
    melt = load_melt_sample_meis(args.melt_vcf, args.sample, args.chrom)
    ont = load_ont_sample_meis(args.ont_svim, args.ont_svan, args.sample, args.chrom)
    try:
        exclude_beds = resolve_exclude_beds(
            args.exclude_bed,
            apply_default_junk=not args.no_junk_filter,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
    if exclude_beds:
        raw_melt, raw_ont = len(melt), len(ont)
        melt = filter_variants_to_regions(melt, exclude_beds=exclude_beds)
        ont = filter_variants_to_regions(ont, exclude_beds=exclude_beds)
        print(
            "exclude "
            + ",".join(str(path) for path in exclude_beds)
            + f" melt {raw_melt}->{len(melt)} ont {raw_ont}->{len(ont)}",
            flush=True,
        )
    events = merge_unique_events([melt, ont], pad_bp=PAD_BP, require_family=False)
    melt_m = overlap_metrics(melt, calls, pad_bp=PAD_BP, require_family=False)
    ont_m = overlap_metrics(ont, calls, pad_bp=PAD_BP, require_family=False)
    union_m = overlap_metrics(events, calls, pad_bp=PAD_BP, require_family=False)
    print(
        f"gold {len(gold)} melt {melt_m['truth_recovered']}/{melt_m['n_truth']} "
        f"ont {ont_m['truth_recovered']}/{ont_m['n_truth']} "
        f"union {union_m['truth_recovered']}/{union_m['n_truth']}",
        flush=True,
    )

    table = overlap_table(gold, events, calls)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / f"{args.sample.lower()}_1kg_overlap_ranks.tsv"
    table.drop(columns=["chrom"]).to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv} rows={len(table)}", flush=True)

    missed = table.loc[table["recovered"].eq("no")].copy()
    print(f"unrecovered {len(missed)}", flush=True)
    if args.no_missed_plots or missed.empty:
        return
    if args.reference_fasta is None or args.disease_bam is None:
        raise SystemExit("Missed-site plots need --reference-fasta and --disease-bam")
    plot_unrecovered(
        review,
        missed,
        calls_path=args.calls,
        outdir=args.outdir,
        reference_fasta=args.reference_fasta,
        disease_bam=args.disease_bam,
        control_bam=args.control_bam,
    )


if __name__ == "__main__":
    sys.exit(main())
