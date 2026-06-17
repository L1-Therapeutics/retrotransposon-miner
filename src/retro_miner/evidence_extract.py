from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam


@dataclass
class ExtractionSummary:
    sample: str
    total_reads_scanned: int
    passing_reads: int
    split_evidence_rows: int
    discordant_evidence_rows: int = 0
    insert_size_threshold: int = 0
    weak_only_discordant_filtered_rows: int = 0


def _normalize_regions(regions: list[str] | str) -> list[str]:
    if isinstance(regions, str):
        return [regions]
    clean = [r.strip() for r in regions if r and r.strip()]
    if not clean:
        raise ValueError("No valid regions provided.")
    return clean


def _iter_reads_for_regions(bam: pysam.AlignmentFile, regions: list[str]):
    for region in regions:
        for read in bam.fetch(region=region):
            yield read


def _collect_soft_clips(read: pysam.AlignedSegment, min_clip_len: int) -> list[tuple[str, int]]:
    if read.cigartuples is None:
        return []

    clips: list[tuple[str, int]] = []
    first_op, first_len = read.cigartuples[0]
    if first_op == 4 and first_len >= min_clip_len:
        clips.append(("L", first_len))

    last_op, last_len = read.cigartuples[-1]
    if last_op == 4 and last_len >= min_clip_len:
        clips.append(("R", last_len))

    return clips


def _poly_at_stats(seq: str) -> tuple[int, float, str]:
    s = (seq or "").upper()
    if not s:
        return (0, 0.0, "")
    best = 0
    best_base = ""
    cur = 0
    prev = ""
    at_bases = 0
    for ch in s:
        if ch in {"A", "T"}:
            at_bases += 1
            if ch == prev:
                cur += 1
            else:
                cur = 1
                prev = ch
            if cur > best:
                best = cur
                best_base = ch
        else:
            cur = 0
            prev = ""
    return (best, float(at_bases) / float(len(s)), best_base)


def _poly_at_breakpoint_proximal_stats(
    read_seq: str,
    window_bases: int,
) -> tuple[int, float, str, str]:
    seq = (read_seq or "").upper()
    if not seq:
        return (0, 0.0, "", "")
    win = max(1, int(window_bases))
    left = seq[:win]
    right = seq[-win:]
    l_run, l_frac, l_base = _poly_at_stats(left)
    r_run, r_frac, r_base = _poly_at_stats(right)
    if (l_run, l_frac) >= (r_run, r_frac):
        return (int(l_run), float(l_frac), l_base, "L")
    return (int(r_run), float(r_frac), r_base, "R")


def extract_split_evidence(
    bam_path: Path,
    sample_name: str,
    outdir: Path,
    regions: list[str] | str,
    min_mapq: int = 20,
    min_clip_len: int = 20,
    poly_tail_rescue_min_clip_len: int = 8,
    poly_tail_rescue_min_run: int = 8,
    poly_tail_rescue_min_frac: float = 0.8,
) -> ExtractionSummary:
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    total_reads_scanned = 0
    passing_reads = 0
    region_list = _normalize_regions(regions)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        # Fetch explicit regions to support targeted chromosome subsets.
        for read in _iter_reads_for_regions(bam, region_list):
            total_reads_scanned += 1

            if read.is_unmapped:
                continue
            if read.is_qcfail or read.is_duplicate or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            passing_reads += 1
            collect_clip_len = min(int(min_clip_len), max(1, int(poly_tail_rescue_min_clip_len)))
            clips = _collect_soft_clips(read, min_clip_len=collect_clip_len)
            if not clips:
                continue

            has_sa = read.has_tag("SA")
            sa_raw = read.get_tag("SA") if has_sa else ""
            nm = int(read.get_tag("NM")) if read.has_tag("NM") else -1
            chrom = bam.get_reference_name(read.reference_id)
            for clip_side, clip_len in clips:
                # Breakpoint coordinate should depend on clipping side:
                # - Left clip: mapped segment starts at breakpoint (reference_start + 1)
                # - Right clip: mapped segment ends at breakpoint (reference_end)
                if clip_side == "L":
                    pos_1based = read.reference_start + 1
                else:
                    pos_1based = read.reference_end
                query_seq = read.query_sequence or ""
                clip_seq = ""
                if query_seq:
                    if clip_side == "L":
                        clip_seq = query_seq[:clip_len]
                    else:
                        clip_seq = query_seq[-clip_len:]
                poly_run, poly_frac, poly_base = _poly_at_stats(clip_seq)
                poly_tail_rescued = (
                    clip_len < min_clip_len
                    and clip_len >= max(1, int(poly_tail_rescue_min_clip_len))
                    and poly_run >= max(1, int(poly_tail_rescue_min_run))
                    and poly_frac >= float(poly_tail_rescue_min_frac)
                )
                if clip_len < min_clip_len and not poly_tail_rescued:
                    continue
                rows.append(
                    {
                        "sample": sample_name,
                        "chrom": chrom,
                        "pos": pos_1based,
                        "clip_side": clip_side,
                        "clip_len": int(clip_len),
                        "mapq": int(read.mapping_quality),
                        "is_reverse": bool(read.is_reverse),
                        "read_name": read.query_name,
                        "has_sa": bool(has_sa),
                        "sa_raw": sa_raw,
                        "clip_seq": clip_seq,
                        "nm": nm,
                        "clip_poly_at_run": int(poly_run),
                        "clip_poly_at_fraction": float(poly_frac),
                        "clip_poly_base": poly_base,
                        "poly_tail_rescued": bool(poly_tail_rescued),
                    }
                )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chrom", "pos", "read_name", "clip_side"], kind="mergesort")
    else:
        df = pd.DataFrame(
            columns=[
                "sample",
                "chrom",
                "pos",
                "clip_side",
                "clip_len",
                "mapq",
                "is_reverse",
                "read_name",
                "has_sa",
                "sa_raw",
                "clip_seq",
                "nm",
                "clip_poly_at_run",
                "clip_poly_at_fraction",
                "clip_poly_base",
                "poly_tail_rescued",
            ]
        )

    tsv_path = outdir / f"split_evidence.{sample_name}.tsv"
    parquet_path = outdir / f"split_evidence.{sample_name}.parquet"
    df.to_csv(tsv_path, sep="\t", index=False)
    df.to_parquet(parquet_path, index=False)

    return ExtractionSummary(
        sample=sample_name,
        total_reads_scanned=total_reads_scanned,
        passing_reads=passing_reads,
        split_evidence_rows=len(df),
    )


def _estimate_insert_size_threshold(
    bam_path: Path,
    regions: list[str] | str,
    min_mapq: int,
    quantile: float,
    fallback_threshold: int,
) -> int:
    insert_sizes: list[int] = []
    region_list = _normalize_regions(regions)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in _iter_reads_for_regions(bam, region_list):
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_qcfail or read.is_duplicate or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_id != read.next_reference_id:
                continue
            abs_tlen = abs(read.template_length)
            if abs_tlen > 0:
                insert_sizes.append(abs_tlen)

    if not insert_sizes:
        return fallback_threshold

    threshold = int(np.quantile(np.asarray(insert_sizes), quantile))
    return max(threshold, fallback_threshold)


def extract_discordant_evidence(
    bam_path: Path,
    sample_name: str,
    outdir: Path,
    regions: list[str] | str,
    min_mapq: int = 20,
    insert_quantile: float = 0.995,
    min_abs_tlen: int = 1000,
    poly_tail_rescue_window_bases: int = 25,
    poly_tail_rescue_min_run: int = 10,
    poly_tail_rescue_min_frac: float = 0.8,
    poly_tail_rescue_min_abs_tlen: int = 500,
    require_strong_discordant_reason: bool = True,
) -> ExtractionSummary:
    outdir.mkdir(parents=True, exist_ok=True)
    insert_threshold = _estimate_insert_size_threshold(
        bam_path=bam_path,
        regions=regions,
        min_mapq=min_mapq,
        quantile=insert_quantile,
        fallback_threshold=min_abs_tlen,
    )

    rows: list[dict[str, Any]] = []
    weak_only_filtered_rows = 0
    total_reads_scanned = 0
    passing_reads = 0
    region_list = _normalize_regions(regions)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in _iter_reads_for_regions(bam, region_list):
            total_reads_scanned += 1

            if not read.is_paired:
                continue
            if read.is_unmapped:
                continue
            if read.is_qcfail or read.is_duplicate or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            passing_reads += 1
            reasons: list[str] = []
            mate_chrom = "*"
            mate_pos_1based = 0
            if read.next_reference_id >= 0:
                mate_chrom = bam.get_reference_name(read.next_reference_id)
                mate_pos_1based = read.next_reference_start + 1

            abs_tlen = abs(read.template_length)
            if read.mate_is_unmapped:
                reasons.append("mate_unmapped")
            elif read.reference_id != read.next_reference_id:
                reasons.append("interchrom")
            else:
                if abs_tlen >= insert_threshold:
                    reasons.append("large_insert")

            # Orientation is weak as a stand-alone MEI signal but useful context.
            if read.is_reverse == read.mate_is_reverse:
                reasons.append("same_strand")
            if not read.is_proper_pair:
                reasons.append("improper_pair")

            read_seq = read.query_sequence or ""
            poly_run, poly_frac, poly_base, poly_side = _poly_at_breakpoint_proximal_stats(
                read_seq,
                window_bases=poly_tail_rescue_window_bases,
            )
            has_structural_context = (
                read.mate_is_unmapped
                or (read.reference_id != read.next_reference_id)
                or (abs_tlen >= max(1, int(poly_tail_rescue_min_abs_tlen)))
                or (read.is_reverse == read.mate_is_reverse)
                or (not read.is_proper_pair)
            )
            poly_tail_anchor_rescued = (
                poly_run >= max(1, int(poly_tail_rescue_min_run))
                and poly_frac >= float(poly_tail_rescue_min_frac)
                and has_structural_context
            )
            if poly_tail_anchor_rescued:
                reasons.append("poly_tail_anchor_rescue")

            if not reasons:
                continue
            strong_reasons = {"mate_unmapped", "interchrom", "large_insert", "poly_tail_anchor_rescue"}
            has_strong_reason = any(r in strong_reasons for r in reasons)
            if require_strong_discordant_reason and not has_strong_reason:
                weak_only_filtered_rows += 1
                continue

            chrom = bam.get_reference_name(read.reference_id)
            pos_1based = read.reference_start + 1
            rows.append(
                {
                    "sample": sample_name,
                    "chrom": chrom,
                    "pos": pos_1based,
                    "mate_chrom": mate_chrom,
                    "mate_pos": mate_pos_1based,
                    "mapq": int(read.mapping_quality),
                    "template_len": int(read.template_length),
                    "is_reverse": bool(read.is_reverse),
                    "mate_is_reverse": bool(read.mate_is_reverse),
                    "is_proper_pair": bool(read.is_proper_pair),
                    "is_read1": bool(read.is_read1),
                    "read_name": read.query_name,
                    "discordant_reasons": ",".join(sorted(set(reasons))),
                    "nm": int(read.get_tag("NM")) if read.has_tag("NM") else -1,
                    "read_seq": read_seq,
                    "anchor_poly_at_run": int(poly_run),
                    "anchor_poly_at_fraction": float(poly_frac),
                    "anchor_poly_base": poly_base,
                    "anchor_poly_side": poly_side,
                    "poly_tail_anchor_rescued": bool(poly_tail_anchor_rescued),
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chrom", "pos", "read_name"], kind="mergesort")
    else:
        df = pd.DataFrame(
            columns=[
                "sample",
                "chrom",
                "pos",
                "mate_chrom",
                "mate_pos",
                "mapq",
                "template_len",
                "is_reverse",
                "mate_is_reverse",
                "is_proper_pair",
                "is_read1",
                "read_name",
                "discordant_reasons",
                "nm",
                "read_seq",
                "anchor_poly_at_run",
                "anchor_poly_at_fraction",
                "anchor_poly_base",
                "anchor_poly_side",
                "poly_tail_anchor_rescued",
            ]
        )

    tsv_path = outdir / f"discordant_evidence.{sample_name}.tsv"
    parquet_path = outdir / f"discordant_evidence.{sample_name}.parquet"
    df.to_csv(tsv_path, sep="\t", index=False)
    df.to_parquet(parquet_path, index=False)

    return ExtractionSummary(
        sample=sample_name,
        total_reads_scanned=total_reads_scanned,
        passing_reads=passing_reads,
        split_evidence_rows=0,
        discordant_evidence_rows=len(df),
        insert_size_threshold=insert_threshold,
        weak_only_discordant_filtered_rows=weak_only_filtered_rows,
    )
