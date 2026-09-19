#!/usr/bin/env python3
"""Cut a small BAM snippet plus extract tables for a real-locus pytest fixture."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import pysam


def _slice_by_pos(df: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    pos = pd.to_numeric(df["pos"], errors="coerce").fillna(0).astype(int)
    soft = (
        pd.to_numeric(df["soft_clip_pos"], errors="coerce").fillna(0).astype(int)
        if "soft_clip_pos" in df.columns
        else pd.Series(0, index=df.index)
    )
    keep = df["chrom"].astype(str).eq(chrom) & (pos.between(start, end) | soft.between(start, end))
    return df.loc[keep].copy()


def _write_sorted_bam(
    *,
    cram: Path,
    reference: Path,
    chrom: str,
    start: int,
    end: int,
    dest: Path,
) -> int:
    raw_path = dest.with_suffix(".unsorted.bam")
    src = pysam.AlignmentFile(str(cram), "rc" if str(cram).endswith(".cram") else "rb", reference_filename=str(reference))
    raw = pysam.AlignmentFile(str(raw_path), "wb", header=src.header)
    written = 0
    seen: set[tuple[str, int, int, int]] = set()

    def _write(read: pysam.AlignedSegment) -> None:
        nonlocal written
        key = (read.query_name or "", int(read.flag), int(read.reference_id), int(read.reference_start))
        if key in seen:
            return
        seen.add(key)
        raw.write(read)
        written += 1

    region_reads: list[pysam.AlignedSegment] = []
    for read in src.fetch(chrom, start - 1, end):
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail:
            continue
        region_reads.append(read)
        _write(read)

    for read in region_reads:
        if not read.is_paired or read.mate_is_unmapped:
            continue
        mate_chrom = read.next_reference_name
        mate_pos = int(read.next_reference_start)
        if mate_chrom is None or mate_pos < 0:
            continue
        if mate_chrom == chrom and start - 1 <= mate_pos < end:
            continue
        try:
            for mate in src.fetch(mate_chrom, max(0, mate_pos - 5), mate_pos + 5):
                if mate.query_name == read.query_name:
                    _write(mate)
                    break
        except ValueError:
            continue

    src.close()
    raw.close()
    pysam.sort("-o", str(dest), str(raw_path))
    pysam.index(str(dest))
    raw_path.unlink(missing_ok=True)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--cram", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--breakpoint", required=True, type=int)
    parser.add_argument("--window-start", required=True, type=int)
    parser.add_argument("--window-end", required=True, type=int)
    parser.add_argument("--catalog-id", required=True)
    parser.add_argument("--sample", default="HG00100")
    parser.add_argument("--pad-bp", type=int, default=800)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    fetch_start = int(args.window_start) - int(args.pad_bp)
    fetch_end = int(args.window_end) + int(args.pad_bp)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    tables: dict[str, pd.DataFrame] = {}
    for name in (
        "split_evidence.disease.parquet",
        "split_evidence.control.parquet",
        "discordant_evidence.disease.parquet",
        "discordant_evidence.control.parquet",
    ):
        src = args.results_dir / name
        if not src.exists():
            raise SystemExit(f"missing {src}")
        sl = _slice_by_pos(pd.read_parquet(src), args.chrom, fetch_start, fetch_end)
        sl.to_parquet(out / name, index=False)
        tables[name] = sl

    detail_src = args.results_dir / "supporting_reads_detail.mei.tsv"
    if detail_src.exists():
        detail = pd.read_csv(detail_src, sep="\t", low_memory=False)
        dpos = pd.to_numeric(detail["genomic_pos"], errors="coerce").fillna(0).astype(int)
        dsoft = (
            pd.to_numeric(detail["soft_clip_pos"], errors="coerce").fillna(0).astype(int)
            if "soft_clip_pos" in detail.columns
            else pd.Series(0, index=detail.index)
        )
        keep = detail["chrom"].astype(str).eq(args.chrom) & (
            dpos.between(fetch_start, fetch_end) | dsoft.between(fetch_start, fetch_end)
        )
        dsub = detail.loc[keep].copy()
        dsub.to_csv(out / "supporting_reads_detail.mei.tsv", sep="\t", index=False)
        dpe_names = set(
            dsub.loc[dsub["evidence_type"].astype(str) == "DPE", "read_name"].dropna().astype(str)
        )
        for name in ("discordant_evidence.disease.parquet", "discordant_evidence.control.parquet"):
            sl = tables[name].copy()
            sl["mei_hit"] = sl["read_name"].astype(str).isin(dpe_names)
            sl["mate_mei_hit"] = sl["mei_hit"]
            sl.to_parquet(out / name, index=False)
            tables[name] = sl

    pd.DataFrame(
        [
            {
                "chrom": args.chrom,
                "window_start": int(args.window_start),
                "window_end": int(args.window_end),
                "tsd_left_breakpoint": 0,
                "tsd_right_breakpoint": 0,
                "tsd_len_estimate": 0,
                "tsd_detected": False,
                "tsd_evidence_source": "",
            }
        ]
    ).to_csv(out / "candidate.tsv", sep="\t", index=False)

    bam_name = f"{args.sample}.{args.chrom}_{fetch_start}_{fetch_end}.bam"
    n_reads = _write_sorted_bam(
        cram=args.cram,
        reference=args.reference_fasta,
        chrom=args.chrom,
        start=fetch_start,
        end=fetch_end,
        dest=out / bam_name,
    )
    manifest = {
        "catalog_id": args.catalog_id,
        "sample": args.sample,
        "chrom": args.chrom,
        "expected_breakpoint": int(args.breakpoint),
        "discovery_window_start": int(args.window_start),
        "discovery_window_end": int(args.window_end),
        "fetch_start": fetch_start,
        "fetch_end": fetch_end,
        "bam": bam_name,
        "n_bam_reads": n_reads,
        "source_results": str(args.results_dir),
        "source_cram": str(args.cram),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
