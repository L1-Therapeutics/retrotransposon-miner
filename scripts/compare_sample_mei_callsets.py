#!/usr/bin/env python3
"""Compare one-sample RTM MEI calls to 1000 Genomes MELT and ONT genotypes.

Truth records keep ID, MEI class (ALU / SVA / LINE1), strand, and length.
Inserted sequences are not loaded.

Example (HG00100 chr22 on the VM):

  PYTHONPATH=src python scripts/compare_sample_mei_callsets.py \\
    --calls ~/retrotransposon-workdir/results/hg00100_germline_chr22/candidate_loci.mei.gold_review.tsv \\
    --sample HG00100 \\
    --chrom chr22 \\
    --melt-vcf $DATA/polymorphism/hg38/1kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz \\
    --ont-svim $DATA/polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.bcf \\
    --ont-svan $DATA/polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.noGt.SVAN_1.3.bcf \\
    --outdir ~/retrotransposon-workdir/results/hg00100_germline_chr22/callset_eval

Junk sites are dropped by default. The bed is
$RTM_PUBLIC_DATA_DIR/annotation/hg38/junk/junk_exclusion_merged.bed
(segdup, low mappability, gaps, ENCODE blacklist), the same mask the caller
uses before silver. Pass --no-junk-filter to score the unfiltered catalogs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from retro_miner.callset_eval import (
    catalog_overlap,
    filter_variants_to_regions,
    label_rtm_calls,
    resolve_exclude_beds,
    load_melt_sample_meis,
    load_ont_sample_meis,
    load_rtm_calls,
    merge_unique_events,
    overlap_metrics,
    suggest_heuristic_cutoff,
    sweep_cutoffs,
)


def _write_variants(path: Path, variants) -> None:
    columns = [
        "chrom",
        "pos",
        "end",
        "id",
        "family",
        "strand",
        "mei_length",
        "source",
        "genotype",
        "svtype",
        "source_ids",
        "member_count",
    ]
    rows = [
        {
            "chrom": v.chrom,
            "pos": v.pos,
            "end": v.end,
            "id": v.variant_id,
            "family": v.family,
            "strand": v.strand,
            "mei_length": v.svlen,
            "source": v.source,
            "genotype": v.genotype,
            "svtype": v.svtype,
            "source_ids": json.dumps(v.extra.get("source_ids", {}), sort_keys=True),
            "member_count": v.extra.get("member_count", 1),
        }
        for v in variants
    ]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--calls", type=Path, required=True, help="RTM gold_review or candidate_loci.mei TSV")
    p.add_argument("--sample", default="HG00100", help="Genotype sample ID (1000G uses HG00100, not HG001)")
    p.add_argument(
        "--chrom",
        default="all",
        help="Chromosome such as chr22, or 'all' for a whole-genome evaluation",
    )
    p.add_argument("--melt-vcf", type=Path, required=True, help="Phase 3 integrated SV genotypes VCF")
    p.add_argument("--ont-svim", type=Path, required=True, help="SVIM-asm multi-sample BCF/VCF")
    p.add_argument("--ont-svan", type=Path, required=True, help="SVAN noGt annotation BCF for those IDs")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--pad-bp", type=int, default=300, help="Overlap pad around the RTM breakpoint")
    p.add_argument("--require-family", action="store_true", help="Require ALU/LINE1/SVA to match")
    p.add_argument("--min-ont-recall", type=float, default=0.80)
    p.add_argument(
        "--exclude-bed",
        type=Path,
        action="append",
        default=None,
        help=(
            "Drop RTM, MELT, and ONT sites inside this BED (repeatable). "
            "Default: junk_exclusion_merged.bed for hg38."
        ),
    )
    p.add_argument(
        "--no-junk-filter",
        action="store_true",
        help="Do not apply the default caller junk mask.",
    )
    p.add_argument(
        "--callable-bed",
        "--highconf-bed",
        dest="callable_bed",
        type=Path,
        action="append",
        default=[],
        help="Keep RTM, MELT, and ONT only inside this callable BED (repeatable)",
    )
    args = p.parse_args()
    query_chrom = None if args.chrom.lower() in {"all", "*"} else args.chrom

    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f"[callset-eval] loading RTM calls from {args.calls}", flush=True)
    calls, _ = load_rtm_calls(args.calls, query_chrom)
    raw_rtm_count = len(calls)
    print(f"[callset-eval] rtm_calls_raw={raw_rtm_count}", flush=True)
    print(f"[callset-eval] loading MELT genotypes {args.melt_vcf}", flush=True)
    melt = load_melt_sample_meis(args.melt_vcf, args.sample, query_chrom)
    raw_melt_count = len(melt)
    print(f"[callset-eval] melt_insertions_raw={raw_melt_count}", flush=True)
    print(f"[callset-eval] loading ONT SVIM+SVAN for {args.sample}", flush=True)
    ont = load_ont_sample_meis(args.ont_svim, args.ont_svan, args.sample, query_chrom)
    raw_ont_count = len(ont)
    print(f"[callset-eval] ont_insertions_raw={raw_ont_count}", flush=True)

    try:
        exclude_beds = resolve_exclude_beds(
            args.exclude_bed,
            apply_default_junk=not args.no_junk_filter,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
    if args.callable_bed or exclude_beds:
        region_args = {
            "include_beds": args.callable_bed,
            "exclude_beds": exclude_beds,
        }
        filtered = filter_variants_to_regions(calls + melt + ont, **region_args)
        calls = [variant for variant in filtered if variant.source == "rtm"]
        melt = [variant for variant in filtered if variant.source == "melt"]
        ont = [variant for variant in filtered if variant.source == "ont"]
        print(
            "[callset-eval] shared callable filter "
            f"rtm={raw_rtm_count}->{len(calls)} "
            f"melt={raw_melt_count}->{len(melt)} "
            f"ont={raw_ont_count}->{len(ont)}",
            flush=True,
        )

    catalogs = catalog_overlap(
        melt, ont, pad_bp=args.pad_bp, require_family=args.require_family
    )
    truth_union = merge_unique_events(
        [melt, ont], pad_bp=args.pad_bp, require_family=args.require_family
    )
    melt_m = overlap_metrics(melt, calls, pad_bp=args.pad_bp, require_family=args.require_family)
    ont_m = overlap_metrics(ont, calls, pad_bp=args.pad_bp, require_family=args.require_family)
    union_m = overlap_metrics(
        truth_union, calls, pad_bp=args.pad_bp, require_family=args.require_family
    )
    labeled = label_rtm_calls(calls, melt, ont, pad_bp=args.pad_bp, require_family=args.require_family)
    sweep = sweep_cutoffs(labeled, melt, ont, pad_bp=args.pad_bp, require_family=args.require_family)
    suggestion = suggest_heuristic_cutoff(sweep, min_ont_recall=args.min_ont_recall)

    _write_variants(args.outdir / "truth_melt.tsv", melt)
    _write_variants(args.outdir / "truth_ont.tsv", ont)
    _write_variants(args.outdir / "truth_unique.tsv", truth_union)
    _write_variants(args.outdir / "truth_shared.tsv", catalogs["shared"])
    _write_variants(args.outdir / "truth_melt_only.tsv", catalogs["left_only"])
    _write_variants(args.outdir / "truth_ont_only.tsv", catalogs["right_only"])
    pd.DataFrame(labeled).to_csv(args.outdir / "rtm_calls.labeled.tsv", sep="\t", index=False)
    pd.DataFrame(sweep).to_csv(args.outdir / "cutoff_sweep.tsv", sep="\t", index=False)

    summary = {
        "sample": args.sample,
        "chrom": args.chrom,
        "pad_bp": args.pad_bp,
        "require_family": args.require_family,
        "callable_beds": [str(path) for path in args.callable_bed],
        "exclude_beds": [str(path) for path in exclude_beds],
        "raw_counts": {
            "rtm": raw_rtm_count,
            "melt_insertions": raw_melt_count,
            "ont_insertions": raw_ont_count,
        },
        "n_rtm_calls": len(calls),
        "melt": melt_m,
        "ont": ont_m,
        "unique_truth": union_m,
        "catalog_overlap": {
            "shared": len(catalogs["shared"]),
            "melt_only": len(catalogs["left_only"]),
            "ont_only": len(catalogs["right_only"]),
            "unique_events": len(truth_union),
        },
        "novel_vs_union": int(sum(1 for r in labeled if r["novel"])),
        "suggested_cutoff": suggestion,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"sample={args.sample} chrom={args.chrom} pad_bp={args.pad_bp}")
    print(f"rtm_calls={len(calls)}")
    print(
        f"MELT insertions={melt_m['n_truth']} recovered={melt_m['truth_recovered']} "
        f"recall={melt_m['recall']:.3f}"
    )
    print(
        f"ONT insertions={ont_m['n_truth']} recovered={ont_m['truth_recovered']} "
        f"recall={ont_m['recall']:.3f}"
    )
    print(
        f"catalog_overlap shared={len(catalogs['shared'])} "
        f"melt_only={len(catalogs['left_only'])} "
        f"ont_only={len(catalogs['right_only'])} "
        f"unique_events={len(truth_union)}"
    )
    print(f"novel_rtm_vs_either={summary['novel_vs_union']}")
    if suggestion:
        print(
            "suggested_cutoff "
            f"{suggestion['cutoff']}={suggestion['value']} "
            f"ont_recall={suggestion['ont_recall']:.3f} "
            f"melt_recall={suggestion['melt_recall']:.3f} "
            f"n_calls={suggestion['n_calls']} "
            f"novel={suggestion['novel_calls']} "
            f"meets_min_ont_recall={suggestion['meets_min_ont_recall']}"
        )
    print(f"wrote {args.outdir}")


if __name__ == "__main__":
    main()
