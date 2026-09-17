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
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from retro_miner.callset_eval import (
    label_rtm_calls,
    load_melt_sample_meis,
    load_ont_sample_meis,
    load_rtm_calls,
    overlap_metrics,
    suggest_heuristic_cutoff,
    sweep_cutoffs,
)


def _write_variants(path: Path, variants) -> None:
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
        }
        for v in variants
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--calls", type=Path, required=True, help="RTM gold_review or candidate_loci.mei TSV")
    p.add_argument("--sample", default="HG00100", help="Genotype sample ID (1000G uses HG00100, not HG001)")
    p.add_argument("--chrom", default="chr22")
    p.add_argument("--melt-vcf", type=Path, required=True, help="Phase 3 integrated SV genotypes VCF")
    p.add_argument("--ont-svim", type=Path, required=True, help="SVIM-asm multi-sample BCF/VCF")
    p.add_argument("--ont-svan", type=Path, required=True, help="SVAN noGt annotation BCF for those IDs")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--pad-bp", type=int, default=200, help="Overlap pad around the RTM breakpoint")
    p.add_argument("--require-family", action="store_true", help="Require ALU/LINE1/SVA to match")
    p.add_argument("--min-ont-recall", type=float, default=0.80)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f"[callset-eval] loading RTM calls from {args.calls}", flush=True)
    calls, _ = load_rtm_calls(args.calls, args.chrom)
    print(f"[callset-eval] rtm_calls={len(calls)}", flush=True)
    print(f"[callset-eval] loading MELT genotypes {args.melt_vcf}", flush=True)
    melt = load_melt_sample_meis(args.melt_vcf, args.sample, args.chrom)
    print(f"[callset-eval] melt_carriers={len(melt)}", flush=True)
    print(f"[callset-eval] loading ONT SVIM+SVAN for {args.sample}", flush=True)
    ont = load_ont_sample_meis(args.ont_svim, args.ont_svan, args.sample, args.chrom)
    print(f"[callset-eval] ont_mei_carriers={len(ont)}", flush=True)

    melt_m = overlap_metrics(melt, calls, pad_bp=args.pad_bp, require_family=args.require_family)
    ont_m = overlap_metrics(ont, calls, pad_bp=args.pad_bp, require_family=args.require_family)
    labeled = label_rtm_calls(calls, melt, ont, pad_bp=args.pad_bp, require_family=args.require_family)
    sweep = sweep_cutoffs(labeled, melt, ont, pad_bp=args.pad_bp, require_family=args.require_family)
    suggestion = suggest_heuristic_cutoff(sweep, min_ont_recall=args.min_ont_recall)

    _write_variants(args.outdir / "truth_melt.tsv", melt)
    _write_variants(args.outdir / "truth_ont.tsv", ont)
    pd.DataFrame(labeled).to_csv(args.outdir / "rtm_calls.labeled.tsv", sep="\t", index=False)
    pd.DataFrame(sweep).to_csv(args.outdir / "cutoff_sweep.tsv", sep="\t", index=False)

    summary = {
        "sample": args.sample,
        "chrom": args.chrom,
        "pad_bp": args.pad_bp,
        "require_family": args.require_family,
        "n_rtm_calls": len(calls),
        "melt": melt_m,
        "ont": ont_m,
        "novel_vs_union": int(sum(1 for r in labeled if r["novel"])),
        "suggested_cutoff": suggestion,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"sample={args.sample} chrom={args.chrom} pad_bp={args.pad_bp}")
    print(f"rtm_calls={len(calls)}")
    print(
        f"MELT carriers={melt_m['n_truth']} recovered={melt_m['truth_recovered']} "
        f"recall={melt_m['recall']:.3f}"
    )
    print(
        f"ONT MEI carriers={ont_m['n_truth']} recovered={ont_m['truth_recovered']} "
        f"recall={ont_m['recall']:.3f}"
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
