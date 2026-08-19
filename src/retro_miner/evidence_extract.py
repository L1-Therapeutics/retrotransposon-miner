from __future__ import annotations

from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import click
import numpy as np
import pandas as pd
import pysam

from ._utils import _longest_poly_at_span, _poly_at_stats


@dataclass
class ExtractionSummary:
    sample: str
    total_reads_scanned: int
    passing_reads: int
    split_evidence_rows: int
    discordant_evidence_rows: int = 0
    insert_size_threshold: int = 0
    weak_only_discordant_filtered_rows: int = 0
    mate_seq_fetched_rows: int = 0
    mate_seq_missing_interchrom_rows: int = 0


def _normalize_regions(regions: list[str] | str) -> list[str]:
    if isinstance(regions, str):
        return [regions]
    clean = [r.strip() for r in regions if r and r.strip()]
    if not clean:
        raise ValueError("No valid regions provided.")
    return clean


def _parse_region_to_bounds(region: str) -> tuple[str, int | None, int | None]:
    """Parse a samtools-style region string into (contig, start0, end0).

    Returned *start0* and *end0* are **0-based half-open** for
    ``pysam.AlignmentFile.fetch(contig, start0, end0)`` — identical to what
    pysam's own region parser produces but computed once per region rather
    than on every BAM fetch call.  A bare chromosome name (no ``:``)
    returns ``(chrom, None, None)`` so the full contig is fetched.

    Samtools region syntax uses 1-based inclusive coordinates, e.g.
    ``chr1:1000-2000`` covers positions 1000–2000 (1-based).  The
    0-based half-open equivalent passed to HTSlib is [999, 2000).
    """
    r = (region or "").strip()
    if not r:
        raise ValueError("Empty region string.")
    if ":" not in r:
        # Bare chromosome name — fetch the full contig.
        return (r, None, None)
    chrom, coords = r.split(":", 1)
    coords = coords.replace(",", "")  # tolerate thousand-separator commas
    if "-" in coords:
        s_str, e_str = coords.split("-", 1)
        # 1-based inclusive → 0-based half-open
        start0 = max(0, int(s_str) - 1)
        end0 = int(e_str)
    else:
        # Single-position region (e.g. "chr1:5000")
        start0 = max(0, int(coords) - 1)
        end0 = start0 + 1
    return (chrom, start0, end0)


def _iter_reads_for_regions(bam: pysam.AlignmentFile, regions: list[str]):
    """Yield reads from a BAM over a list of samtools-style region strings.

    Region strings are parsed to integer (contig, start0, end0) bounds **once**
    before entering the HTSlib fetch loop.  Using integer bounds avoids the
    per-call string-parse overhead inside ``pysam.AlignmentFile.fetch`` and
    matches the recommended HTSlib calling convention for repeated fetches.
    """
    for region in regions:
        contig, start0, end0 = _parse_region_to_bounds(region)
        # pysam.fetch(contig) with no bounds fetches the full contig;
        # pysam.fetch(contig, start0, end0) uses 0-based half-open HTSlib coords.
        for read in bam.fetch(contig, start0, end0):
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


def _soft_clip_query_seq(query_seq: str, side: str, clip_len: int) -> str:
    """Extract leftmost (L) or rightmost (R) soft-clip bases from a query sequence."""
    seq = query_seq or ""
    side_u = (side or "").upper()[:1]
    n = int(clip_len or 0)
    if n <= 0 or len(seq) < n:
        return ""
    if side_u == "L":
        return seq[:n]
    if side_u == "R":
        return seq[-n:]
    return ""


def _longest_soft_clip_from_read(
    read: pysam.AlignedSegment,
) -> tuple[str, int, int, str]:
    """Return ``(side, clip_len, clip_pos_1based, clip_seq)`` for the longest soft clip.

    ``clip_pos_1based`` is the genomic junction tip (left clip → alignment start;
    right clip → alignment end). Empty when the read has no soft clip.
    """
    if read.cigartuples is None:
        return ("", 0, 0, "")
    query_seq = read.query_sequence or ""
    pos_1based = int(read.reference_start) + 1 if read.reference_start is not None else 0
    ref_end_1based = int(read.reference_end) if read.reference_end is not None else pos_1based
    first_op, first_len = read.cigartuples[0]
    last_op, last_len = read.cigartuples[-1]
    candidates: list[tuple[str, int, int]] = []
    if first_op == 4 and first_len > 0:
        candidates.append(("L", int(first_len), pos_1based))
    if last_op == 4 and last_len > 0:
        candidates.append(("R", int(last_len), ref_end_1based))
    if not candidates:
        return ("", 0, 0, "")
    side, clip_len, clip_pos = max(candidates, key=lambda x: x[1])
    return (side, clip_len, clip_pos, _soft_clip_query_seq(query_seq, side, clip_len))


def _clip_to_poly_at_region(seq: str, *, min_dom_frac: float = 0.90) -> str:
    """Return the longest mostly-A or mostly-T substring (empty if none)."""
    _length, _frac, _base, span = _longest_poly_at_span(seq, min_frac=float(min_dom_frac))
    return span


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
    short_mei_rescue_min_clip_len: int = 12,
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
            mate_chrom = "*"
            mate_pos_1based = 0
            if read.is_paired and not read.mate_is_unmapped and read.next_reference_id >= 0:
                mate_chrom = bam.get_reference_name(read.next_reference_id)
                mate_pos_1based = read.next_reference_start + 1
            collect_clip_len = min(
                int(min_clip_len),
                max(1, int(poly_tail_rescue_min_clip_len)),
                max(1, int(short_mei_rescue_min_clip_len)),
            )
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
                span_len, poly_frac, poly_base, poly_region = _longest_poly_at_span(clip_seq, min_frac=0.90, min_len=8)
                if span_len <= 0:
                    poly_run, poly_frac, poly_base = _poly_at_stats(clip_seq)
                else:
                    poly_run = span_len
                    # keep poly_frac/base from span
                poly_tail_rescued = (
                    clip_len < min_clip_len
                    and clip_len >= max(1, int(poly_tail_rescue_min_clip_len))
                    and poly_run >= max(1, int(poly_tail_rescue_min_run))
                    and poly_frac >= float(poly_tail_rescue_min_frac)
                )
                short_mei_candidate = (
                    clip_len < int(min_clip_len)
                    and clip_len >= max(1, int(short_mei_rescue_min_clip_len))
                )
                if clip_len < min_clip_len and not poly_tail_rescued and not short_mei_candidate:
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
                        "mate_chrom": mate_chrom,
                        "mate_pos": mate_pos_1based,
                        "has_sa": bool(has_sa),
                        "sa_raw": sa_raw,
                        "clip_seq": clip_seq,
                        "nm": nm,
                        "short_mei_candidate": bool(short_mei_candidate and not poly_tail_rescued),
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
                "mate_chrom",
                "mate_pos",
                "has_sa",
                "sa_raw",
                "clip_seq",
                "nm",
                "short_mei_candidate",
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


def _fetch_mate_sequence_from_bam(
    bam: pysam.AlignmentFile,
    read_name: str,
    mate_chrom: str,
    mate_pos_1based: int,
    *,
    fetch_window_bp: int = 500,
) -> tuple[str, int, int, str, int, str]:
    """Best-effort mate sequence fetch for discordant-pair MEI remapping.

    Returns
    ``(mate_seq, mate_ref_start, mate_ref_end, mate_soft_clip_side,
    mate_soft_clip_len, mate_soft_clip_seq)``.
    """
    empty = ("", 0, 0, "", 0, "")
    if mate_chrom in {"", "*"} or mate_pos_1based <= 0:
        return empty
    if mate_chrom not in bam.references:
        return empty

    start0 = max(0, int(mate_pos_1based) - 1)
    end0 = start0 + max(1, int(fetch_window_bp))
    for mate_read in bam.fetch(mate_chrom, start0, end0):
        if mate_read.query_name != read_name:
            continue
        if mate_read.is_secondary or mate_read.is_supplementary:
            continue
        if mate_read.is_unmapped:
            continue
        seq = mate_read.query_sequence or ""
        if not seq:
            continue
        clip_side, clip_len, _clip_pos, clip_seq = _longest_soft_clip_from_read(mate_read)
        return (
            seq,
            int(mate_read.reference_start) + 1 if mate_read.reference_start is not None else 0,
            int(mate_read.reference_end) if mate_read.reference_end is not None else 0,
            clip_side,
            int(clip_len),
            clip_seq,
        )
    return empty


def _validate_mate_fetch_bam(
    scan_bam_path: Path,
    mate_bam_path: Path | None,
    regions: list[str],
) -> None:
    """Warn when region-scanned BAM cannot resolve interchrom mate sequences."""
    if mate_bam_path is not None and mate_bam_path != scan_bam_path:
        click.echo(
            f"[extract-discordant] using mate-resolution BAM {mate_bam_path} "
            f"(scan BAM {scan_bam_path})"
        )
        return

    scan_chroms = {region.split(":", 1)[0] for region in regions}
    try:
        with pysam.AlignmentFile(str(scan_bam_path), "rb") as bam:
            stats = bam.get_index_statistics()
            if not stats:
                return
            other_with_reads = 0
            for i, ref in enumerate(bam.references):
                if ref in scan_chroms:
                    continue
                if i < len(stats) and (int(stats[i].mapped) + int(stats[i].unmapped)) > 0:
                    other_with_reads += 1
            if other_with_reads == 0:
                click.echo(
                    "[extract-discordant] warning: scan BAM appears chromosome-subset "
                    f"({scan_bam_path}); interchrom mate sequences will be missing unless "
                    "--disease-mate-bam/--control-mate-bam points to a full-genome BAM."
                )
    except (OSError, ValueError, RuntimeError):
        return


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
    mate_bam_path: Path | None = None,
    mate_fetch_window_bp: int = 500,
    fetch_mate_seq: bool = True,
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
    mate_seq_fetched_rows = 0
    mate_seq_missing_interchrom_rows = 0
    region_list = _normalize_regions(regions)
    mate_bam_resolved = mate_bam_path if mate_bam_path is not None else bam_path
    if fetch_mate_seq:
        _validate_mate_fetch_bam(bam_path, mate_bam_path, region_list)

    # Mate BAM is only opened when fetch_mate_seq is enabled. Annotate can re-fetch
    # mates for candidate loci later; skipping here is the main extract speedup.
    mate_bam_ctx = (
        pysam.AlignmentFile(str(mate_bam_resolved), "rb") if fetch_mate_seq else nullcontext(None)
    )
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, mate_bam_ctx as mate_bam:
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
            soft_clip_side, soft_clip_len, soft_clip_pos, soft_clip_seq = _longest_soft_clip_from_read(read)
            if fetch_mate_seq and mate_bam is not None:
                (
                    mate_seq,
                    mate_ref_start,
                    mate_ref_end,
                    mate_soft_clip_side,
                    mate_soft_clip_len,
                    mate_soft_clip_seq,
                ) = _fetch_mate_sequence_from_bam(
                    mate_bam,
                    read.query_name,
                    mate_chrom,
                    mate_pos_1based,
                    fetch_window_bp=mate_fetch_window_bp,
                )
                if mate_seq:
                    mate_seq_fetched_rows += 1
                elif "interchrom" in reasons:
                    mate_seq_missing_interchrom_rows += 1
            else:
                mate_seq = ""
                mate_ref_start = 0
                mate_ref_end = 0
                mate_soft_clip_side = ""
                mate_soft_clip_len = 0
                mate_soft_clip_seq = ""
            rows.append(
                {
                    "sample": sample_name,
                    "chrom": chrom,
                    "pos": int(read.reference_start) + 1 if read.reference_start is not None else 0,
                    "ref_end": int(read.reference_end) if read.reference_end is not None else (
                        int(read.reference_start) + 1 if read.reference_start is not None else 0
                    ),
                    "soft_clip_side": soft_clip_side,
                    "soft_clip_len": int(soft_clip_len),
                    "soft_clip_pos": int(soft_clip_pos),
                    "soft_clip_seq": soft_clip_seq,
                    "mate_chrom": mate_chrom,
                    "mate_pos": mate_pos_1based,
                    "mate_seq": mate_seq,
                    "mate_ref_start": int(mate_ref_start),
                    "mate_ref_end": int(mate_ref_end),
                    "mate_soft_clip_side": mate_soft_clip_side,
                    "mate_soft_clip_len": int(mate_soft_clip_len),
                    "mate_soft_clip_seq": mate_soft_clip_seq,
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
                "ref_end",
                "soft_clip_side",
                "soft_clip_len",
                "soft_clip_pos",
                "soft_clip_seq",
                "mate_chrom",
                "mate_pos",
                "mate_seq",
                "mate_ref_start",
                "mate_ref_end",
                "mate_soft_clip_side",
                "mate_soft_clip_len",
                "mate_soft_clip_seq",
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

    if not df.empty and mate_seq_missing_interchrom_rows > 0:
        interchrom_total = int(df["discordant_reasons"].fillna("").astype(str).str.contains("interchrom").sum())
        missing_frac = mate_seq_missing_interchrom_rows / max(interchrom_total, 1)
        click.echo(
            f"[extract-discordant] sample={sample_name} interchrom_rows={interchrom_total} "
            f"mate_seq_missing={mate_seq_missing_interchrom_rows} ({missing_frac:.1%}); "
            "MEI_MAPPED discordant support may be undercounted."
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
        mate_seq_fetched_rows=mate_seq_fetched_rows,
        mate_seq_missing_interchrom_rows=mate_seq_missing_interchrom_rows,
    )
