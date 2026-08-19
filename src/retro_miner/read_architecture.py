#!/usr/bin/env python3
"""Schematic read-pair architecture plot for a gold-review MEI locus.

Composite axis (forward / ``+`` orientation):

    left flank  |  MEI 5′→3′  |  polyA zone  |  right flank

Reverse (``-``) swaps MEI direction and puts the polyA zone on the 3′ side
(left of the MEI body). The MEI segment is consensus 5′–3′ only; polyA/T
rescues and polyA split clips are drawn in the dedicated polyA zone whose
width is the longest observed polyA/T support at the locus.
"""

from __future__ import annotations

import argparse
import hashlib
import random
import re
import textwrap
import time
from dataclasses import dataclass, field
from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

# L1 palette (only these four colors are used in the figure).
COLOR_LIGHT_ORANGE = "#feb24c"
COLOR_DARK_ORANGE = "#fd8d3c"
COLOR_BLACK = "black"
COLOR_WHITE = "white"

# Illumina PE ~150 bp: a single polyA mate/clip is at most one read.
# Anchor-clip + polyA-mate on the same DPE can cover non-overlapping stretches,
# so allow up to ~2 reads when summing those lengths for the plot zone/bar.
_MAX_POLYA_SINGLE_BP = 151
_MAX_POLYA_COMBINED_BP = 280

# Slim columns pulled from candidate_loci.mei.tsv when gold_review lacks them.
_MEI_ENRICH_COLS = (
    "chrom",
    "window_start",
    "window_end",
    "insertion_breakpoint_pos",
    "insertion_breakpoint_interval_start",
    "insertion_breakpoint_interval_end",
    "disease_discordant_mei_left_target_pos_median",
    "disease_discordant_mei_right_target_pos_median",
    "disease_discordant_mei_insertion_span_estimate",
    "control_discordant_mei_left_target_pos_median",
    "control_discordant_mei_right_target_pos_median",
    "control_discordant_mei_insertion_span_estimate",
    "g1k_melt_insertion_length",
    "consensus_insertion_mei_5p_coord",
    "consensus_insertion_mei_3p_coord",
    "consensus_insertion_mei_span",
    "consensus_insertion_mei_5p_coord_full",
    "consensus_insertion_mei_3p_coord_full",
    "consensus_insertion_mei_span_full",
)

_DETAIL_COLS = (
    "sample",
    "evidence_type",
    "read_name",
    "chrom",
    "window_start",
    "window_end",
    "anchor_side",
    "genomic_pos",
    "mate_chrom",
    "mate_genomic_pos",
    "mei_start",
    "mei_end",
    "mate_mei_start",
    "mate_mei_end",
    "mei_hit",
    "mate_mei_hit",
    "short_mei_seed_rescued",
    "vntr_rescue",
    "polya_rescue",
    "mei_hit_source",
    "mate_seq_len",
    "polya_base",
    "clip_poly_base",
    "clip_len",
    "soft_clip_side",
    "soft_clip_len",
    "soft_clip_pos",
    "poly_tail_rescued",
    "clip_poly_at_run",
    "anchor_poly_at_run",
    "anchor_poly_base",
    "anchor_poly_side",
    "poly_tail_anchor_rescued",
)


@dataclass(frozen=True)
class LocusLayout:
    chrom: str
    window_start: int
    window_end: int
    breakpoint: int
    breakpoint_left: int
    breakpoint_right: int
    mei_5p: int
    mei_3p: int
    flank_bp: int
    mei_span_bp: int
    span_source: str = ""
    insert_size_estimates: tuple[int, ...] = ()
    # "-" means the insertion is reverse-oriented: left flank abuts MEI 3', right abuts 5'.
    orientation: str = "+"
    # Dedicated polyA/T tail zone beyond MEI 3' (max observed polyA length).
    polya_zone_bp: int = 0

    @property
    def left_flank_bp(self) -> int:
        return self.flank_bp

    @property
    def right_flank_bp(self) -> int:
        return self.flank_bp

    @property
    def reverse_oriented(self) -> bool:
        return str(self.orientation).strip() in {"-", "−", "rev", "reverse", "-1"}

    @property
    def total_width(self) -> int:
        return int(self.flank_bp + self.polya_zone_bp + self.mei_span_bp + self.flank_bp)

    @property
    def polya_region_start_x(self) -> float:
        """Start of the polyA zone (abuts MEI 3')."""
        if self.polya_zone_bp <= 0:
            # Degenerate: collapse onto the 3' MEI edge.
            return self.mei_3p_junction_x
        if self.reverse_oriented:
            return float(self.flank_bp)
        return float(self.flank_bp + self.mei_span_bp)

    @property
    def polya_region_end_x(self) -> float:
        return self.polya_region_start_x + float(max(0, self.polya_zone_bp))

    @property
    def mei_region_start_x(self) -> float:
        if self.reverse_oriented and self.polya_zone_bp > 0:
            return float(self.flank_bp + self.polya_zone_bp)
        return float(self.flank_bp)

    @property
    def mei_region_end_x(self) -> float:
        return self.mei_region_start_x + float(self.mei_span_bp)

    @property
    def mei_3p_junction_x(self) -> float:
        """X coordinate of the MEI 3′ edge (where the polyA zone abuts)."""
        if self.reverse_oriented:
            return self.mei_region_start_x
        return self.mei_region_end_x

    @property
    def insertion_left_x(self) -> float:
        """Left edge of the insertion block (polyA or MEI) abutting left flank."""
        if self.reverse_oriented and self.polya_zone_bp > 0:
            return self.polya_region_start_x
        return self.mei_region_start_x

    @property
    def insertion_right_x(self) -> float:
        """Right edge of the insertion block abutting right flank."""
        if (not self.reverse_oriented) and self.polya_zone_bp > 0:
            return self.polya_region_end_x
        return self.mei_region_end_x

    def genomic_to_x(self, pos: int) -> float:
        """Map genomic coordinate onto the composite axis (clamped to visible flanks)."""
        pos = int(pos)
        left_end = self.insertion_left_x
        right_start = self.insertion_right_x
        if pos <= self.breakpoint:
            x = left_end - max(0, self.breakpoint - pos)
            return max(0.0, min(x, left_end))
        x = right_start + max(0, pos - self.breakpoint)
        return max(right_start, min(x, float(self.total_width)))

    def mei_coord_to_x(self, coord: int) -> float:
        """Map MEI consensus coordinate onto the MEI body segment (not polyA zone).

        Forward (+): 5' at left junction, 3' at right junction.
        Reverse (-): 3' at left junction, 5' at right junction.
        """
        coord = int(coord)
        lo = int(self.mei_5p)
        hi = int(self.mei_3p)
        if hi <= lo:
            return self.mei_region_start_x
        coord = max(lo, min(coord, hi))
        if self.reverse_oriented:
            rel = hi - coord
        else:
            rel = coord - lo
        rel = max(0, min(rel, self.mei_span_bp))
        return self.mei_region_start_x + rel

    def polya_span_x(self, width_bp: int) -> tuple[float, float]:
        """Place a polyA/T bar of ``width_bp`` at the distal (flank) end of the polyA zone.

        The zone width is the locus max polyA. A shorter bar of length L in a zone
        of width W sits at offset ``W - L`` from MEI 3′ (abutting the flank), not at
        ``0..L``. Example: max=180, read=60 → span 120–180 from MEI 3′.
        """
        width = float(max(1, min(int(width_bp), max(1, self.polya_zone_bp or int(width_bp)))))
        if self.polya_zone_bp <= 0:
            # Fallback: tiny stub just past 3'.
            if self.reverse_oriented:
                return max(0.0, self.mei_3p_junction_x - width), self.mei_3p_junction_x
            return self.mei_3p_junction_x, self.mei_3p_junction_x + width
        if self.reverse_oriented:
            # left flank | polyA | MEI — abut left flank, grow toward MEI.
            x0 = self.polya_region_start_x
            x1 = min(self.polya_region_end_x, x0 + width)
        else:
            # MEI | polyA | right flank — abut right flank, grow toward MEI.
            x1 = self.polya_region_end_x
            x0 = max(self.polya_region_start_x, x1 - width)
        if x1 <= x0:
            x1 = x0 + 1.0
        return x0, x1


def _row_int(row: pd.Series, col: str, default: int = 0) -> int:
    if col not in row.index:
        return default
    val = row[col]
    if pd.isna(val):
        return default
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return default


def _target_family(target: str) -> str:
    t = (target or "").upper()
    if "SVA" in t:
        return "SVA"
    if "ALU" in t:
        return "ALU"
    if "LINE1" in t or "LINE/L1" in t or re.search(r"(?:^|[^A-Z0-9])L1(?:[^A-Z0-9]|$)", t):
        return "LINE1"
    if "HERV" in t or "ERV" in t:
        return "ERV"
    return "OTHER"


def _row_mei_family(row: pd.Series | None) -> str:
    if row is None:
        return ""
    for col in (
        "consensus_mei_family",
        "mei_family",
        "disease_discordant_mei_family",
        "control_discordant_mei_family",
    ):
        fam = _target_family(str(row.get(col, "") or ""))
        if fam in {"ALU", "LINE1", "SVA"}:
            return fam
    return ""


def _mei_coords_from_detail(
    detail: pd.DataFrame,
    *,
    family: str = "",
) -> tuple[int, int] | None:
    """MEI body extent from consensus-mapped hits only (exclude polyA rescues).

    When ``family`` is ALU/LINE1/SVA, ignore off-family hits (e.g. stray L1
    mates at an Alu locus) so the plot axis is not inflated to ~6 kb.
    """
    coords: list[int] = []
    work = detail
    # PolyA-rescue mates use synthetic coords past MEI 3′; exclude from body span.
    if "polya_rescue" in detail.columns:
        work = detail.loc[~detail["polya_rescue"].fillna(False).astype(bool)]
    if "mei_hit_source" in work.columns:
        src = work["mei_hit_source"].fillna("").astype(str)
        work = work.loc[~src.eq("polya_rescue")]
    fam = (family or "").upper()
    if fam in {"ALU", "LINE1", "SVA"}:
        keep = pd.Series(True, index=work.index)
        for tcol in ("mei_target", "mate_mei_target", "target"):
            if tcol not in work.columns:
                continue
            tfam = work[tcol].fillna("").astype(str).map(_target_family)
            # Keep rows whose target is empty/OTHER (unknown) or matches locus family.
            keep = keep & (tfam.eq("") | tfam.eq("OTHER") | tfam.eq(fam))
        # If filtering would drop everything, fall back to unfiltered (better than empty).
        filtered = work.loc[keep]
        if not filtered.empty:
            work = filtered
    for start_col, end_col, hit_col in (
        ("mei_start", "mei_end", "mei_hit"),
        ("mate_mei_start", "mate_mei_end", "mate_mei_hit"),
    ):
        if start_col not in work.columns or end_col not in work.columns:
            continue
        hits = work
        if hit_col in work.columns:
            hits = work.loc[work[hit_col].astype(bool)]
        # Per-hit family gate using the column that owns these coords.
        tcol = "mei_target" if start_col.startswith("mei_") and not start_col.startswith("mate_") else "mate_mei_target"
        if fam in {"ALU", "LINE1", "SVA"} and tcol in hits.columns:
            tfam = hits[tcol].fillna("").astype(str).map(_target_family)
            hits = hits.loc[tfam.eq(fam)]
        starts = pd.to_numeric(hits[start_col], errors="coerce").fillna(0).astype(int)
        ends = pd.to_numeric(hits[end_col], errors="coerce").fillna(0).astype(int)
        for start, end in zip(starts.tolist(), ends.tolist()):
            if start > 0:
                coords.append(start)
            if end > 0:
                coords.append(end)
    if len(coords) < 2:
        return None
    return min(coords), max(coords)


def _mate_polya_width_bp(row: pd.Series) -> int:
    """Observed polyA/T width for a mate — never raw mate_seq_len alone.

    Prefer purity-checked ``mate_seq``, then imputed polyA rescue spans.
    ``mate_seq_len`` is only a fallback when the row is explicitly a polyA
    rescue (annotate already required a mostly-A/T mate).
    """
    mate_seq = str(row.get("mate_seq", "") or "")
    if len(mate_seq) >= 12:
        s = mate_seq.upper()
        n_a = s.count("A")
        n_t = s.count("T")
        dom = max(n_a, n_t)
        if dom / float(len(s)) >= 0.90:
            return int(min(_MAX_POLYA_SINGLE_BP, len(s)))
    for start_col, end_col in (("mate_mei_start", "mate_mei_end"), ("mei_start", "mei_end")):
        if start_col not in row.index or end_col not in row.index:
            continue
        try:
            start = int(float(row.get(start_col) or 0))
            end = int(float(row.get(end_col) or 0))
        except (TypeError, ValueError):
            continue
        if start > 0 and end >= start:
            # Only trust synthetic/imputed spans on polyA-rescue rows.
            is_rescue = bool(row.get("polya_rescue", False)) or str(
                row.get("mei_hit_source", "") or ""
            ).strip().lower() in {"polya_rescue", "poly_tail_rescue"}
            if is_rescue:
                return int(max(0, min(_MAX_POLYA_SINGLE_BP, end - start + 1)))
    is_rescue = bool(row.get("polya_rescue", False)) or str(
        row.get("mei_hit_source", "") or ""
    ).strip().lower() in {"polya_rescue", "poly_tail_rescue"}
    if is_rescue and "mate_seq_len" in row.index:
        try:
            n = int(float(row.get("mate_seq_len") or 0))
        except (TypeError, ValueError):
            n = 0
        if n >= 12:
            return int(min(_MAX_POLYA_SINGLE_BP, n))
    return 0


def _max_polya_zone_bp(
    detail: pd.DataFrame | None,
    *,
    consensus_poly_at_min_bp: int = 0,
) -> int:
    """Width of the dedicated polyA zone = longest observed polyA/T support."""
    widths: list[int] = []
    gold = int(consensus_poly_at_min_bp or 0)
    if gold >= 8:
        widths.append(gold)
    if detail is None or detail.empty:
        return int(max(0, min(_MAX_POLYA_COMBINED_BP, max(widths)))) if widths else 0

    def _add_span_widths(frame: pd.DataFrame) -> None:
        for start_col, end_col in (("mate_mei_start", "mate_mei_end"), ("mei_start", "mei_end")):
            if start_col not in frame.columns or end_col not in frame.columns:
                continue
            s = pd.to_numeric(frame[start_col], errors="coerce").fillna(0).astype(int)
            e = pd.to_numeric(frame[end_col], errors="coerce").fillna(0).astype(int)
            widths.extend(
                int(ee - ss + 1) for ss, ee in zip(s.tolist(), e.tolist()) if ss > 0 and ee >= ss
            )

    if "polya_rescue" in detail.columns:
        resc = detail.loc[detail["polya_rescue"].fillna(False).astype(bool)]
        if not resc.empty:
            # Use imputed polyA spans / purity-checked mate width — not raw mate_seq_len
            # on non-rescue MEI mates (that inflated zones to ~clip+151).
            widths.extend(int(v) for v in resc.apply(_mate_polya_width_bp, axis=1).tolist() if int(v) >= 12)
            _add_span_widths(resc)
    if "mei_hit_source" in detail.columns:
        src = detail["mei_hit_source"].fillna("").astype(str)
        poly_src = detail.loc[src.eq("polya_rescue")]
        if not poly_src.empty:
            widths.extend(
                int(v) for v in poly_src.apply(_mate_polya_width_bp, axis=1).tolist() if int(v) >= 12
            )
            _add_span_widths(poly_src)
    if "polya_base" in detail.columns:
        poly_base = detail.loc[detail["polya_base"].fillna("").astype(str).isin({"A", "T", "a", "t"})]
        if not poly_base.empty:
            if "clip_len" in poly_base.columns:
                cl = pd.to_numeric(poly_base["clip_len"], errors="coerce").fillna(0).astype(int)
                widths.extend(int(v) for v in cl.tolist() if int(v) >= 12)
            widths.extend(
                int(v) for v in poly_base.apply(_mate_polya_width_bp, axis=1).tolist() if int(v) >= 12
            )
    if "poly_tail_rescued" in detail.columns:
        rescued = detail.loc[detail["poly_tail_rescued"].fillna(False).astype(bool)]
        if not rescued.empty and "clip_len" in rescued.columns:
            cl = pd.to_numeric(rescued["clip_len"], errors="coerce").fillna(0).astype(int)
            widths.extend(int(v) for v in cl.tolist() if int(v) >= 12)
    # DPE anchors clipped into polyA at the flank junction.
    if "poly_tail_anchor_rescued" in detail.columns:
        anc = detail.loc[detail["poly_tail_anchor_rescued"].fillna(False).astype(bool)]
    else:
        anc = detail.iloc[0:0].copy()
    if anc.empty and "anchor_poly_at_run" in detail.columns:
        runs = pd.to_numeric(detail["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int)
        anc = detail.loc[runs.ge(8)]
    if not anc.empty:
        if "anchor_poly_at_run" in anc.columns:
            runs = pd.to_numeric(anc["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int)
            widths.extend(int(v) for v in runs.tolist() if int(v) >= 8)
        if "soft_clip_len" in anc.columns:
            cl = pd.to_numeric(anc["soft_clip_len"], errors="coerce").fillna(0).astype(int)
            widths.extend(int(v) for v in cl.tolist() if int(v) >= 8)
        # Same DPE: junction polyA clip + polyA mate → sum only when mate is
        # actually polyA (rescue / purity), never raw mate_seq_len of MEI mates.
        clip_lens = (
            pd.to_numeric(anc["soft_clip_len"], errors="coerce").fillna(0).astype(int)
            if "soft_clip_len" in anc.columns
            else pd.Series(0, index=anc.index)
        )
        run_lens = (
            pd.to_numeric(anc["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int)
            if "anchor_poly_at_run" in anc.columns
            else pd.Series(0, index=anc.index)
        )
        anchor_w = clip_lens.combine(run_lens, max)
        mate_w = anc.apply(_mate_polya_width_bp, axis=1)
        for aw, mw in zip(anchor_w.tolist(), mate_w.tolist()):
            if int(aw) >= 8 and int(mw) >= 12:
                widths.append(int(aw) + int(mw))
    if not widths:
        return 0
    return int(max(12, min(_MAX_POLYA_COMBINED_BP, max(widths))))


def _discordant_mei_axis(row: pd.Series, sample: str) -> tuple[int, int, int] | None:
    for prefix in (sample, "disease", "control"):
        left = _row_int(row, f"{prefix}_discordant_mei_left_target_pos_median")
        right = _row_int(row, f"{prefix}_discordant_mei_right_target_pos_median")
        if left > 0 and right > 0:
            mei_5p = min(left, right)
            mei_3p = max(left, right)
            return mei_5p, mei_3p, mei_3p - mei_5p + 1
        span = _row_int(row, f"{prefix}_discordant_mei_insertion_span_estimate")
        if span > 0:
            return 1, span, span
    return None


def _breakpoint_interval(row: pd.Series) -> tuple[int, int, int]:
    bp = _row_int(row, "consensus_insertion_breakpoint_pos")
    if bp <= 0:
        bp = _row_int(row, "insertion_breakpoint_pos")
    if bp <= 0:
        ws = _row_int(row, "window_start")
        we = _row_int(row, "window_end")
        bp = (ws + we) // 2
    bp_l = _row_int(row, "insertion_breakpoint_interval_start")
    bp_r = _row_int(row, "insertion_breakpoint_interval_end")
    if bp_l <= 0 or bp_r < bp_l:
        bp_l = _row_int(row, "consensus_breakpoint_interval_start")
        bp_r = _row_int(row, "consensus_breakpoint_interval_end")
    if bp_l <= 0 or bp_r < bp_l:
        bp_l = bp_r = bp
    return bp_l, bp_r, bp


def _plausible_local_mate_pos(mate_pos: int, *, window_start: int, window_end: int, slack: int = 600) -> bool:
    return window_start - slack <= mate_pos <= window_end + slack


def _genomic_insert_size_estimates(
    detail: pd.DataFrame | None,
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    bp_l: int,
    bp_r: int,
    breakpoint: int,
) -> list[int]:
    """Estimate insertion sizes from discordant pairs with locally plausible mates."""
    if detail is None or detail.empty:
        return []

    sizes: list[int] = []
    dpe = detail.loc[detail["evidence_type"].astype(str) == "DPE"]
    for rec in dpe.itertuples(index=False):
        anchor = int(rec.genomic_pos)
        mate_chrom = str(getattr(rec, "mate_chrom", "") or "")
        mate_pos = int(getattr(rec, "mate_genomic_pos", 0) or 0)
        if mate_chrom != chrom or mate_pos <= 0:
            continue
        if not _plausible_local_mate_pos(mate_pos, window_start=window_start, window_end=window_end):
            continue
        if mate_pos == anchor:
            continue
        if anchor <= bp_l and mate_pos >= bp_r:
            sizes.append(mate_pos - anchor)
        elif anchor >= bp_r and mate_pos <= bp_l:
            sizes.append(anchor - mate_pos)
        elif anchor >= bp_r and mate_pos > anchor:
            sizes.append(mate_pos - breakpoint)
        elif anchor <= bp_l and mate_pos < anchor:
            sizes.append(breakpoint - mate_pos)
    return sizes


def _insertion_span_from_evidence(
    row: pd.Series,
    *,
    sample: str,
    detail: pd.DataFrame | None,
    bp_l: int,
    bp_r: int,
    breakpoint: int,
) -> tuple[int, int, int, str]:
    """Return (mei_5p, mei_3p, span, source_label) for the plotted MEI axis.

    Axis width is the extent of same-family MEI consensus mappings only.
    Off-family minimap hits (e.g. L1 mates at an Alu locus) are ignored.
    Genomic insert-size / TLEN estimates are never used to set the axis.
    """
    del bp_l, bp_r, breakpoint  # retained for call-site compatibility
    family = _row_mei_family(row)

    def _consensus_span() -> tuple[int, int, int, str] | None:
        full_5p = _row_int(row, "consensus_insertion_mei_5p_coord_full")
        full_3p = _row_int(row, "consensus_insertion_mei_3p_coord_full")
        full_span = _row_int(row, "consensus_insertion_mei_span_full")
        if full_5p > 0 and full_3p > full_5p:
            return full_5p, full_3p, full_3p - full_5p + 1, "consensus_coords_full"
        if full_span > 0 and full_5p > 0:
            return full_5p, full_5p + full_span - 1, full_span, "consensus_span_full"
        if full_span > 0 and full_3p > 0:
            return max(1, full_3p - full_span + 1), full_3p, full_span, "consensus_span_full"
        mei_5p = _row_int(row, "consensus_insertion_mei_5p_coord")
        mei_3p = _row_int(row, "consensus_insertion_mei_3p_coord")
        span = _row_int(row, "consensus_insertion_mei_span")
        if mei_5p > 0 and mei_3p > mei_5p:
            return mei_5p, mei_3p, mei_3p - mei_5p + 1, "consensus_coords"
        if span > 0:
            if mei_5p > 0 and mei_3p <= 0:
                mei_3p = mei_5p + span - 1
            elif mei_3p > 0 and mei_5p <= 0:
                mei_5p = max(1, mei_3p - span + 1)
            else:
                mei_5p, mei_3p = 1, span
            return mei_5p, mei_3p, span, "consensus_span"
        return None

    # Family-plausible upper bounds for read-derived axes (guards residual noise).
    max_plausible = {"ALU": 400, "SVA": 2000, "LINE1": 7000}.get(family, 7000)

    detail_extent = (
        _mei_coords_from_detail(detail, family=family)
        if detail is not None and not detail.empty
        else None
    )
    if detail_extent is not None:
        detail_5p, detail_3p = detail_extent
        detail_span = detail_3p - detail_5p + 1
        if detail_span >= 1 and detail_span <= max_plausible:
            return detail_5p, detail_3p, detail_span, "read_mei_coords"

    consensus = _consensus_span()
    if consensus is not None:
        return consensus

    disc = _discordant_mei_axis(row, sample)
    if disc is not None and disc[2] <= max_plausible:
        return disc[0], disc[1], disc[2], "discordant_mei_targets"

    melt_len = _row_int(row, "g1k_melt_insertion_length")
    if melt_len > 0:
        return 1, melt_len, melt_len, "g1k_melt_length"

    # Last-resort axis width when no SR/DPE/consensus evidence exists.
    # Same default for ALU/SVA/LINE1 — do not invent family-specific spans.
    return 1, 300, 300, "mei_default"


def _read_tsv(path: Path, *, usecols: tuple[str, ...] | None = None) -> pd.DataFrame:
    """Load a TSV, optionally restricting to existing columns in ``usecols``."""
    if usecols is None:
        return pd.read_csv(path, sep="\t", low_memory=False)
    wanted = set(usecols)

    def _keep(name: str) -> bool:
        return name in wanted

    return pd.read_csv(path, sep="\t", usecols=_keep, low_memory=False)


def _load_slim_mei_table(mei_tsv: Path) -> pd.DataFrame:
    return _read_tsv(mei_tsv, usecols=_MEI_ENRICH_COLS)


def _enrich_row_from_mei(
    row: pd.Series,
    mei_df: pd.DataFrame | None,
    *,
    mei_index: dict[tuple[str, int, int], pd.Series] | None = None,
) -> pd.Series:
    """Fill missing plot fields from a slim candidate_loci.mei table."""
    chrom = str(row["chrom"])
    ws = int(row["window_start"])
    we = int(row["window_end"])
    src: pd.Series | None = None
    if mei_index is not None:
        src = mei_index.get((chrom, ws, we))
    elif mei_df is not None and not mei_df.empty:
        mei_hit = mei_df.loc[
            (mei_df["chrom"].astype(str) == chrom)
            & (pd.to_numeric(mei_df["window_start"], errors="coerce") == ws)
            & (pd.to_numeric(mei_df["window_end"], errors="coerce") == we)
        ]
        if not mei_hit.empty:
            src = mei_hit.iloc[0]
    if src is None:
        return row
    out = row.copy()
    for col in src.index:
        if col in {"chrom", "window_start", "window_end"}:
            continue
        if col not in out.index or pd.isna(out.get(col)) or out.get(col) in {0, -1, "", "nan"}:
            out[col] = src[col]
    return out


def _build_mei_index(mei_df: pd.DataFrame | None) -> dict[tuple[str, int, int], pd.Series]:
    if mei_df is None or mei_df.empty:
        return {}
    index: dict[tuple[str, int, int], pd.Series] = {}
    chrom_s = mei_df["chrom"].astype(str)
    ws = pd.to_numeric(mei_df["window_start"], errors="coerce").fillna(-1).astype(int)
    we = pd.to_numeric(mei_df["window_end"], errors="coerce").fillna(-1).astype(int)
    for i in range(len(mei_df)):
        key = (str(chrom_s.iloc[i]), int(ws.iloc[i]), int(we.iloc[i]))
        if key not in index:
            index[key] = mei_df.iloc[i]
    return index


def _select_locus_row(gold_df: pd.DataFrame, chrom: str, pos: int) -> pd.Series:
    subset = gold_df.loc[gold_df["chrom"].astype(str) == chrom]
    if subset.empty:
        raise ValueError(f"No rows for {chrom} in gold review table")
    hit = subset.loc[
        (pd.to_numeric(subset["window_start"], errors="coerce") <= pos)
        & (pd.to_numeric(subset["window_end"], errors="coerce") >= pos)
    ]
    if hit.empty:
        hit = subset.assign(
            dist=(pd.to_numeric(subset["window_start"], errors="coerce") - pos).abs()
            + (pd.to_numeric(subset["window_end"], errors="coerce") - pos).abs()
        ).sort_values("dist").head(1)
    return hit.iloc[0].copy()


def _load_locus_row(
    gold_tsv: Path,
    chrom: str,
    pos: int,
    *,
    gold_df: pd.DataFrame | None = None,
    mei_df: pd.DataFrame | None = None,
) -> pd.Series:
    df = gold_df if gold_df is not None else _read_tsv(gold_tsv)
    row = _select_locus_row(df, chrom, pos)
    if mei_df is None:
        mei_tsv = gold_tsv.parent / "candidate_loci.mei.tsv"
        if mei_tsv.exists():
            mei_df = _load_slim_mei_table(mei_tsv)
    return _enrich_row_from_mei(row, mei_df)


def _flank_side(genomic_pos: int, *, bp_l: int, bp_r: int) -> str | None:
    if genomic_pos <= bp_l:
        return "L"
    if genomic_pos >= bp_r:
        return "R"
    return None


def _layout_from_row(
    row: pd.Series,
    *,
    sample: str,
    detail: pd.DataFrame | None,
    flank_bp: int,
) -> LocusLayout:
    chrom = str(row["chrom"])
    ws = int(row["window_start"])
    we = int(row["window_end"])
    bp_l, bp_r, bp = _breakpoint_interval(row)
    insert_sizes = tuple(
        sorted(
            _genomic_insert_size_estimates(
                detail,
                chrom=chrom,
                window_start=ws,
                window_end=we,
                bp_l=bp_l,
                bp_r=bp_r,
                breakpoint=bp,
            )
        )
    )
    mei_5p, mei_3p, mei_span, span_source = _insertion_span_from_evidence(
        row,
        sample=sample,
        detail=detail,
        bp_l=bp_l,
        bp_r=bp_r,
        breakpoint=bp,
    )
    orientation = _infer_orientation(
        detail,
        table_orientation=str(row.get("consensus_insertion_orientation", "") or ""),
        breakpoint=bp,
        breakpoint_left=bp_l,
        breakpoint_right=bp_r,
        window_start=ws,
        window_end=we,
    )
    gold_poly = 0
    try:
        gold_poly = int(float(row.get("consensus_poly_at_min_bp", 0) or 0))
    except (TypeError, ValueError):
        gold_poly = 0
    polya_zone_bp = _max_polya_zone_bp(detail, consensus_poly_at_min_bp=gold_poly)
    return LocusLayout(
        chrom=chrom,
        window_start=ws,
        window_end=we,
        breakpoint=bp,
        breakpoint_left=bp_l,
        breakpoint_right=bp_r,
        mei_5p=mei_5p,
        mei_3p=mei_3p,
        flank_bp=flank_bp,
        mei_span_bp=mei_span,
        span_source=span_source,
        insert_size_estimates=insert_sizes,
        orientation=orientation,
        polya_zone_bp=polya_zone_bp,
    )


def _infer_orientation(
    detail: pd.DataFrame | None,
    *,
    table_orientation: str,
    breakpoint: int | None = None,
    breakpoint_left: int | None = None,
    breakpoint_right: int | None = None,
    window_start: int | None = None,
    window_end: int | None = None,
) -> str:
    """Infer insertion orientation from SR/DPE MEI coordinates.

    Family-agnostic (Alu/SVA/LINE-1 use the same geometry):

    Forward (+): left genomic junction abuts MEI 5' (low coords), right abuts 3'.
    Reverse (-): left abuts MEI 3' (high), right abuts 5'.

    Evidence mapping onto junctions:
      - DPE: genomic anchor flank is the junction the mate should abut.
      - SR: BAM L-softclip maps on the right flank → clip abuts the *right*
        junction; R-softclip maps on the left flank → *left* junction.

    Prefer exact-locus rows (neighbor windows often belong to other insertions
    and can invert DPE medians). Prefer SR when both junctions are observed.
    """
    table = (table_orientation or "+").strip() or "+"
    if detail is None or detail.empty:
        return table

    work = detail
    if (
        window_start is not None
        and window_end is not None
        and "window_start" in detail.columns
        and "window_end" in detail.columns
    ):
        exact = detail.loc[
            (pd.to_numeric(detail["window_start"], errors="coerce") == int(window_start))
            & (pd.to_numeric(detail["window_end"], errors="coerce") == int(window_end))
        ]
        if not exact.empty:
            work = exact

    bp_l = int(breakpoint_left) if breakpoint_left is not None else None
    bp_r = int(breakpoint_right) if breakpoint_right is not None else None
    bp = int(breakpoint) if breakpoint is not None else None

    def _junction_side(rec, *, evidence: str) -> str | None:
        if evidence == "DPE":
            genomic_pos = int(getattr(rec, "genomic_pos", 0) or 0)
            if bp_l is not None and bp_r is not None and genomic_pos > 0:
                side = _flank_side(genomic_pos, bp_l=bp_l, bp_r=bp_r)
                if side is not None:
                    return side
            if bp is not None and genomic_pos > 0:
                return "L" if genomic_pos <= bp else "R"
            recorded = str(getattr(rec, "anchor_side", "") or "").upper()
            return recorded if recorded in {"L", "R"} else None

        # SR: clip_side is BAM soft-clip side, not the insertion junction.
        clip = str(
            getattr(rec, "clip_side", "") or getattr(rec, "anchor_side", "") or ""
        ).upper()
        if clip == "L":
            return "R"
        if clip == "R":
            return "L"
        return None

    def _collect(evidence: str) -> tuple[list[float], list[float]]:
        left_mids: list[float] = []
        right_mids: list[float] = []
        for rec in work.itertuples(index=False):
            if bool(getattr(rec, "polya_rescue", False)):
                continue
            if str(getattr(rec, "mei_hit_source", "") or "") == "polya_rescue":
                continue
            et = str(rec.evidence_type)
            if evidence == "DPE":
                if et != "DPE" or not bool(getattr(rec, "mate_mei_hit", False)):
                    continue
                start = int(getattr(rec, "mate_mei_start", 0) or 0)
                end = int(getattr(rec, "mate_mei_end", 0) or 0)
            else:
                if et != "SR" or not bool(getattr(rec, "mei_hit", False)):
                    continue
                start = int(getattr(rec, "mei_start", 0) or 0)
                end = int(getattr(rec, "mei_end", 0) or 0)
            if start <= 0 or end <= 0:
                continue
            side = _junction_side(rec, evidence=evidence)
            if side not in {"L", "R"}:
                continue
            mid = (start + end) / 2.0
            if side == "L":
                left_mids.append(mid)
            else:
                right_mids.append(mid)
        return left_mids, right_mids

    def _call(left_mids: list[float], right_mids: list[float]) -> str | None:
        if not left_mids or not right_mids:
            return None
        left_med = sorted(left_mids)[len(left_mids) // 2]
        right_med = sorted(right_mids)[len(right_mids) // 2]
        if abs(left_med - right_med) < 50:
            return None
        # Left junction nearer the high (3') end ⇒ reverse insertion.
        return "-" if left_med > right_med else "+"

    # Split reads abut the junctions directly — prefer when both sides exist.
    sr_call = _call(*_collect("SR"))
    if sr_call is not None:
        return sr_call
    dpe_call = _call(*_collect("DPE"))
    if dpe_call is not None:
        return dpe_call
    return table


def _parse_support_counts(support_str: object) -> dict[str, int]:
    """Parse ``SR_L=..,SR_R=..,DPE_L=..,DPE_R=..,MEI_MAPPED=..,polyA_MAPPED=..,VNTR_MAPPED=..``."""
    text = "" if support_str is None or (isinstance(support_str, float) and pd.isna(support_str)) else str(support_str)
    out = {
        "SR_L": 0,
        "SR_R": 0,
        "DPE_L": 0,
        "DPE_R": 0,
        "MEI_MAPPED": 0,
        "polyA_MAPPED": 0,
        "VNTR_MAPPED": 0,
    }
    for key in out:
        m = re.search(rf"{re.escape(key)}=([0-9]+)", text)
        if m:
            out[key] = int(m.group(1))
    return out


def _support_score(support_str: object) -> tuple[int, int, int, int]:
    """Rank sample support: MEI_MAPPED, polyA_MAPPED, VNTR_MAPPED, then SR+DPE."""
    counts = _parse_support_counts(support_str)
    flank = counts["SR_L"] + counts["SR_R"] + counts["DPE_L"] + counts["DPE_R"]
    return counts["MEI_MAPPED"], counts["polyA_MAPPED"], counts["VNTR_MAPPED"], flank


def _choose_sample(row: pd.Series, sample: str) -> str:
    """Pick disease/control with the most support when sample is auto/empty."""
    requested = (sample or "auto").strip().lower()
    if requested in {"disease", "control"}:
        return requested
    disease_score = _support_score(row.get("disease_supporting_reads", ""))
    control_score = _support_score(row.get("control_supporting_reads", ""))
    if control_score > disease_score:
        return "control"
    if disease_score > control_score:
        return "disease"
    # Tie: prefer disease when both present, else whichever has any flank support.
    if any(disease_score):
        return "disease"
    if any(control_score):
        return "control"
    return "disease"


def _sample_status_label(row: pd.Series) -> str:
    """Locus call status for plot labels: shared / disease_only / control_only / …"""
    raw = row.get("sample_status_label", "")
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return ""
    text = str(raw).strip()
    return text if text else ""


def _normalize_mei_coords(
    start: int,
    end: int,
    layout: LocusLayout,
) -> tuple[int, int] | None:
    """Clamp MEI consensus coords into the plotted insertion footprint.

    Does not relocate hits to the opposite end — if the axis is too short,
    extend it via `_insertion_span_from_evidence` instead.
    """
    start = int(start)
    end = int(end)
    if start <= 0 and end <= 0:
        return None
    if end < start:
        start, end = end, start

    lo = int(layout.mei_5p)
    hi = int(layout.mei_3p)
    if hi < lo:
        return None
    start_c = max(lo, min(start, hi))
    end_c = max(lo, min(end, hi))
    if end_c < start_c:
        start_c, end_c = end_c, start_c
    return start_c, end_c


def _mei_span_x(
    start: int,
    end: int,
    layout: LocusLayout,
) -> tuple[float, float] | None:
    normalized = _normalize_mei_coords(start, end, layout)
    if normalized is None:
        return None
    start_n, end_n = normalized
    x0 = layout.mei_coord_to_x(start_n)
    x1 = layout.mei_coord_to_x(end_n)
    if layout.reverse_oriented:
        # Higher consensus coords map leftward; keep (x0, x1) left-to-right.
        return min(x0, x1), max(x0, x1)
    return min(x0, x1), max(x0, x1)


def _flag_true(rec, *names: str) -> bool:
    """True if any named attribute is a non-null truthy value (NaN → False)."""
    for name in names:
        val = getattr(rec, name, False)
        if val is None:
            continue
        try:
            if pd.isna(val):
                continue
        except (TypeError, ValueError):
            pass
        if bool(val):
            return True
    return False


def _rescue_kind(rec) -> str:
    """Return 'polya', 'vntr', or '' for second-pass MEI rescues."""
    source = str(getattr(rec, "mei_hit_source", "") or "").strip().lower()
    if _flag_true(rec, "polya_rescue") or source == "polya_rescue":
        return "polya"
    if _flag_true(rec, "vntr_rescue") or source == "vntr_rescue":
        return "vntr"
    return ""


def _polya_base_from_rec(rec) -> str:
    """Dominant A/T base for a polyA-rescue mate (default A)."""
    base = str(getattr(rec, "polya_base", "") or "").upper()
    if base in {"A", "T"}:
        return base
    mate_seq = str(getattr(rec, "mate_seq", "") or "").upper()
    if mate_seq:
        n_a = mate_seq.count("A")
        n_t = mate_seq.count("T")
        if n_t > n_a:
            return "T"
        if n_a > 0:
            return "A"
    return "A"


def _polya_tail_label(width_bp: int, base: str = "A") -> str:
    """Plot label as polyA length (e.g. ``66A``), never strand-dependent ``T``."""
    del base  # sequenced mate may be polyT; display as polyA length
    return f"{int(max(1, width_bp))}A"


def _draw_polya_label(
    ax,
    *,
    x0: float,
    x1: float,
    yi: float,
    label: str,
) -> None:
    """Compact length+base label centered on the polyA bar (small-plot friendly)."""
    span = float(x1) - float(x0)
    if span <= 0 or not label:
        return
    # ~7" figure / ~700 bp axis ⇒ short tails are tiny; keep type readable but small.
    if span >= 50:
        fs = 5.0
    elif span >= 25:
        fs = 4.5
    else:
        fs = 4.0
    ax.text(
        (float(x0) + float(x1)) / 2.0,
        float(yi),
        label,
        ha="center",
        va="center",
        fontsize=fs,
        color=COLOR_BLACK,
        zorder=5,
        clip_on=True,
        fontweight="bold",
    )


def _polya_tail_width_from_rec(rec, *, default: int = 0) -> int:
    """Width of a polyA-rescue bar from purity-checked mate / imputed span."""
    if isinstance(rec, pd.Series):
        w = _mate_polya_width_bp(rec)
        if w >= 12:
            return int(w)
    else:
        # namedtuple / row-like from itertuples
        row = {
            "mate_seq": getattr(rec, "mate_seq", ""),
            "mate_seq_len": getattr(rec, "mate_seq_len", 0),
            "mate_mei_start": getattr(rec, "mate_mei_start", 0),
            "mate_mei_end": getattr(rec, "mate_mei_end", 0),
            "mei_start": getattr(rec, "mei_start", 0),
            "mei_end": getattr(rec, "mei_end", 0),
            "polya_rescue": getattr(rec, "polya_rescue", False),
            "mei_hit_source": getattr(rec, "mei_hit_source", ""),
        }
        w = _mate_polya_width_bp(pd.Series(row))
        if w >= 12:
            return int(w)
        try:
            pw = int(float(getattr(rec, "polya_width", 0) or 0))
        except (TypeError, ValueError):
            pw = 0
        if pw >= 12:
            return int(min(_MAX_POLYA_SINGLE_BP, pw))
    return int(default)


def _combined_dpe_polya_width(anchor_width: int, mate_width: int) -> int:
    """Minimum polyA length when junction clip and mate are both polyA.

    Treat as non-overlapping stretches of the same tail (clip at the flank
    junction + mate further into the A-tail), so length ≥ clip + mate.
    """
    a = max(0, int(anchor_width))
    m = max(0, int(mate_width))
    if a >= 8 and m >= 12:
        return int(min(_MAX_POLYA_COMBINED_BP, a + m))
    return int(min(_MAX_POLYA_COMBINED_BP, max(a, m)))


def _polya_tail_span_x(layout: LocusLayout, width_bp: int) -> tuple[float, float]:
    """Place a polyA/T bar inside the dedicated polyA zone (distal / flank-aligned)."""
    return layout.polya_span_x(width_bp)


def _dpe_has_anchor_polya_clip(rec) -> bool:
    """True when the genomic DPE anchor itself is clipped / polyA-tailed at the junction."""
    if bool(getattr(rec, "poly_tail_anchor_rescued", False)):
        return True
    try:
        run = int(float(getattr(rec, "anchor_poly_at_run", 0) or 0))
    except (TypeError, ValueError):
        run = 0
    if run >= 8:
        return True
    try:
        soft = int(float(getattr(rec, "soft_clip_len", 0) or 0))
    except (TypeError, ValueError):
        soft = 0
    base = str(getattr(rec, "anchor_poly_base", "") or "").upper()
    return soft >= 8 and base in {"A", "T"}


def _dpe_anchor_polya_width(rec, layout: LocusLayout) -> int:
    """Width of the orange polyA overhang drawn from the flank junction into the polyA zone."""
    try:
        soft = int(float(getattr(rec, "soft_clip_len", 0) or 0))
    except (TypeError, ValueError):
        soft = 0
    try:
        run = int(float(getattr(rec, "anchor_poly_at_run", 0) or 0))
    except (TypeError, ValueError):
        run = 0
    width = max(soft, run, 0)
    if width < 8:
        return 0
    if layout.polya_zone_bp > 0:
        width = min(width, int(layout.polya_zone_bp))
    return int(min(_MAX_POLYA_SINGLE_BP, max(8, width)))


def _dpe_anchor_polya_base(rec) -> str:
    base = str(getattr(rec, "anchor_poly_base", "") or "").upper()
    if base in {"A", "T"}:
        return "A"  # display as polyA length
    return _polya_base_from_rec(rec) or "A"


def _dpe_anchor_abuts_polya_zone(side: str, layout: LocusLayout) -> bool:
    """True when the genomic flank abuts the dedicated polyA zone (3′ side)."""
    side_u = str(side or "").upper()
    if layout.reverse_oriented:
        return side_u == "L"
    return side_u == "R"


def _polya_span_from_flank_x(layout: LocusLayout, width_bp: int) -> tuple[float, float]:
    """Place a polyA bar at the flank|polyA boundary (same as ``polya_span_x``)."""
    return layout.polya_span_x(int(width_bp))


def _sr_is_polya_clip(rec, layout: LocusLayout) -> bool:
    """True when a split clip should be drawn in the polyA zone, not the MEI body.

    Requires an explicit polyA signal (rescue flag or long A/T run). A lone
    ``clip_poly_base`` of A/T is not enough — most softclips report a dominant
    base even when they are ordinary MEI remaps.
    """
    if _rescue_kind(rec) == "polya":
        return True
    if _flag_true(rec, "polya_rescue", "poly_tail_rescued"):
        return True
    try:
        run = getattr(rec, "clip_poly_at_run", 0)
        if run is not None and not (isinstance(run, float) and pd.isna(run)):
            if int(run) >= 8:
                return True
    except (TypeError, ValueError):
        pass
    try:
        mei_start = int(getattr(rec, "mei_start", 0) or 0)
    except (TypeError, ValueError):
        mei_start = 0
    # Synthetic polyA coords are tacked past consensus 3′.
    if layout.polya_zone_bp > 0 and mei_start > int(layout.mei_3p):
        return True
    return False


def _brk_clp_stub_span_x(layout: LocusLayout, clip_side: str, clip_len: int) -> tuple[float, float]:
    """Junction-anchored orange stub for BRK_CLP clips without MEI remap coords.

    Clip into the full insertion block (MEI + polyA zone). Clipping to the MEI
    body alone collapses right-junction stubs into a ~1 bp sliver whenever the
    polyA zone occupies the 3' edge.
    """
    width = max(1, int(clip_len))
    left = float(layout.insertion_left_x)
    right = float(layout.insertion_right_x)
    if str(clip_side).upper() == "L":
        x0 = left
        x1 = x0 + width
    else:
        x1 = right
        x0 = x1 - width
    return _clip_span(x0, x1, left, right)


def _imputed_rescue_mei_span(rec, layout: LocusLayout) -> tuple[float, float] | None:
    """Place polyA/VNTR-rescued mates that lack real consensus coords."""
    kind = _rescue_kind(rec)
    if not kind:
        return None
    if kind == "polya":
        return _polya_tail_span_x(layout, _polya_tail_width_from_rec(rec))
    lo = int(layout.mei_5p)
    hi = int(layout.mei_3p)
    if hi < lo:
        return None
    span = max(1, hi - lo + 1)
    width = min(40, span)
    start = lo + span // 3
    end = lo + (2 * span) // 3
    if end < start:
        end = min(hi, start + max(1, width - 1))
    return _mei_span_x(start, end, layout)


def _resolve_anchor_side(rec, layout: LocusLayout) -> str | None:
    """Flank side from genomic position; fall back to detail anchor_side at the BP."""
    anchor_pos = int(rec.genomic_pos)
    side = _flank_side(anchor_pos, bp_l=layout.breakpoint_left, bp_r=layout.breakpoint_right)
    if side is not None:
        return side
    recorded = str(getattr(rec, "anchor_side", "") or "").upper()
    if recorded in {"L", "R"}:
        return recorded
    if anchor_pos <= layout.breakpoint:
        return "L"
    return "R"


def _clip_len_from_rec(rec) -> int:
    """Soft-clip length from evidence columns, else MEI alignment span."""
    for col in ("clip_len", "soft_clip_len"):
        val = getattr(rec, col, None)
        if val is not None and not (isinstance(val, float) and pd.isna(val)):
            try:
                n = int(float(val))
            except (TypeError, ValueError):
                n = 0
            if n > 0:
                return n
    mei_start = int(getattr(rec, "mei_start", 0) or 0)
    mei_end = int(getattr(rec, "mei_end", 0) or 0)
    if mei_start > 0 and mei_end >= mei_start:
        return mei_end - mei_start + 1
    return 0


def _pair_from_sr_row(rec, layout: LocusLayout, *, read_width_bp: int = 150) -> dict[str, object] | None:
    """Split read using real soft-clip length from evidence.

    BAM geometry:
      L softclip at junction J, clip_len C, schematic read width R:
        orange (MEI/polyA) = C bp into insertion from the left junction
        black (ref)  = (R-C) bp on the RIGHT flank starting at J
      R softclip at junction J:
        orange (MEI/polyA) = C bp into insertion from the right junction
        black (ref)  = (R-C) bp on the LEFT flank ending at J

    ``side`` is the genomic flank holding the mapped (black) portion.
    PolyA/T clips are placed in the dedicated polyA zone beyond MEI 3′.
    """
    mei_hit = _flag_true(rec, "mei_hit")
    polya_clip = _sr_is_polya_clip(rec, layout)
    if not (mei_hit or polya_clip):
        return None
    junction = int(rec.genomic_pos)
    clip_side = str(getattr(rec, "anchor_side", "") or getattr(rec, "clip_side", "") or "").upper()
    if clip_side not in {"L", "R"}:
        clip_side = _resolve_anchor_side(rec, layout) or ""
    if clip_side not in {"L", "R"}:
        return None

    clip_len = _clip_len_from_rec(rec)
    if clip_len <= 0:
        return None
    read_w = max(int(read_width_bp), clip_len + 1)
    mapped_len = max(1, read_w - clip_len)

    if clip_side == "L":
        plot_side = "R"
        ref_start = junction
        ref_end = junction + mapped_len - 1
    else:
        plot_side = "L"
        ref_end = junction
        ref_start = junction - mapped_len + 1

    polya_width = 0
    polya_base = ""
    polya_label = ""
    if polya_clip:
        polya_width = int(max(clip_len, _polya_tail_width_from_rec(rec, default=clip_len)))
        if layout.polya_zone_bp > 0:
            polya_width = min(polya_width, int(layout.polya_zone_bp))
        polya_base = _polya_base_from_rec(rec)
        if not polya_base:
            polya_base = str(getattr(rec, "clip_poly_base", "") or "A").upper() or "A"
            if polya_base not in {"A", "T"}:
                polya_base = "A"
        polya_label = _polya_tail_label(polya_width, polya_base)
        mei = layout.polya_span_x(polya_width)
        remote_kind = "sr_polya"
    else:
        mei = _mei_span_x(int(getattr(rec, "mei_start", 0) or 0), int(getattr(rec, "mei_end", 0) or 0), layout)
        if mei is None:
            return None
        remote_kind = "sr_mei"

    pair: dict[str, object] = {
        "anchor_pos": junction,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "clip_len": clip_len,
        "mapped_len": mapped_len,
        "clip_side": clip_side,
        "anchor_x": layout.genomic_to_x(junction),
        "remote_x": (mei[0] + mei[1]) / 2.0,
        "remote_mei_x0": mei[0],
        "remote_mei_x1": mei[1],
        "mei_start": int(getattr(rec, "mei_start", 0) or 0),
        "mei_end": int(getattr(rec, "mei_end", 0) or 0),
        "side": plot_side,
        "interchrom": False,
        "mei_mapped": bool(mei_hit),
        "read_name": str(rec.read_name),
        "evidence_type": "SR",
        "remote_kind": remote_kind,
        "anchor_on_reference": True,
        "is_split": True,
        "rescue_kind": "polya" if polya_clip else "",
        "polya_width": polya_width,
        "polya_base": polya_base,
        "polya_label": polya_label,
        "short_mei_seed_rescued": _flag_true(rec, "short_mei_seed_rescued"),
    }
    if bool(getattr(rec, "mate_mei_hit", False)):
        mate_mei = _mei_span_x(
            int(getattr(rec, "mate_mei_start", 0) or 0),
            int(getattr(rec, "mate_mei_end", 0) or 0),
            layout,
        )
        if mate_mei is not None:
            pair["mate_mei_x0"] = mate_mei[0]
            pair["mate_mei_x1"] = mate_mei[1]
            pair["mate_in_mei"] = True
    return pair


def _pair_from_dpe_row(rec, layout: LocusLayout) -> dict[str, object] | None:
    """Discordant pair: reference anchor + MEI-mapped mate.

    If the schematic 150 bp anchor would cross the breakpoint, the overhang is
    drawn inside the MEI (clipped sequence placed in the insertion).
    """
    rescue_kind = _rescue_kind(rec)
    if not (
        _flag_true(rec, "mate_mei_hit", "mei_hit", "polya_rescue", "vntr_rescue")
        or bool(rescue_kind)
    ):
        return None

    anchor_pos = int(rec.genomic_pos)
    anchor_side = _resolve_anchor_side(rec, layout)
    if anchor_side is None:
        return None
    if layout.breakpoint_left < anchor_pos < layout.breakpoint_right:
        return None
    # Require a real visible genomic anchor. BP-parked discordant rows that are
    # really soft-clipped splits (pos≈BP, almost no flank span) produce fake
    # MEI-only drawings when a schematic overhang is added.
    if anchor_side == "L" and layout.breakpoint - anchor_pos < 30:
        return None
    if anchor_side == "R" and anchor_pos - layout.breakpoint < 0:
        return None

    mate_chrom = str(getattr(rec, "mate_chrom", "") or "")
    remote_kind = ""
    remote_mei_x0 = 0.0
    remote_mei_x1 = 0.0
    mei_start = 0
    mei_end = 0
    polya_width = 0
    polya_base = ""
    polya_label = ""

    # PolyA rescues are not real consensus hits — always tack a polyA/T bar onto
    # the MEI 3′ end instead of trusting clamped/zero remap coords.
    if rescue_kind == "polya":
        polya_width = _polya_tail_width_from_rec(rec)
        polya_base = _polya_base_from_rec(rec)
        polya_label = _polya_tail_label(polya_width, polya_base)
        remote_mei_x0, remote_mei_x1 = _polya_tail_span_x(layout, polya_width)
        remote_kind = "mate_mei"
        # Keep synthetic coords as a tail past consensus 3′ for detail tables.
        mei_start = int(layout.mei_3p) + 1
        mei_end = int(layout.mei_3p) + int(polya_width)
    else:
        if bool(getattr(rec, "mate_mei_hit", False)):
            mei_start = int(getattr(rec, "mate_mei_start", 0) or 0)
            mei_end = int(getattr(rec, "mate_mei_end", 0) or 0)
            mate_mei = _mei_span_x(mei_start, mei_end, layout)
            if mate_mei is not None:
                remote_kind = "mate_mei"
                remote_mei_x0, remote_mei_x1 = mate_mei

        if remote_kind == "" and bool(getattr(rec, "mei_hit", False)):
            mei_start = int(getattr(rec, "mei_start", 0) or 0)
            mei_end = int(getattr(rec, "mei_end", 0) or 0)
            anchor_mei = _mei_span_x(mei_start, mei_end, layout)
            if anchor_mei is not None:
                remote_kind = "anchor_mei"
                remote_mei_x0, remote_mei_x1 = anchor_mei

        if remote_kind == "" and rescue_kind:
            imputed = _imputed_rescue_mei_span(rec, layout)
            if imputed is not None:
                remote_kind = "mate_mei"
                remote_mei_x0, remote_mei_x1 = imputed
                lo, hi = int(layout.mei_5p), int(layout.mei_3p)
                span = max(1, hi - lo + 1)
                width = min(40, span)
                mei_start = lo + span // 3
                mei_end = lo + (2 * span) // 3
                if mei_end < mei_start:
                    mei_end = mei_start + max(1, width - 1)

    if remote_kind == "":
        return None

    anchor_polya_clip = False
    anchor_polya_width = 0
    anchor_polya_base = ""
    anchor_polya_label = ""
    anchor_polya_x0 = 0.0
    anchor_polya_x1 = 0.0
    if _dpe_has_anchor_polya_clip(rec):
        # Ensure a polyA zone exists so the overhang has somewhere to land.
        # (Zone width is usually set from locus-wide max; fall back to this clip.)
        aw = _dpe_anchor_polya_width(rec, layout)
        if aw <= 0 and layout.polya_zone_bp <= 0:
            aw = max(12, int(getattr(rec, "anchor_poly_at_run", 0) or 12))
        if aw > 0:
            # Junction clip + polyA mate: minimum tail length is the sum
            # (non-overlapping stretches of the same A-tail).
            if rescue_kind == "polya" and polya_width >= 12:
                combined = _combined_dpe_polya_width(aw, polya_width)
                if layout.polya_zone_bp > 0:
                    combined = min(combined, int(layout.polya_zone_bp))
                polya_width = combined
                polya_label = _polya_tail_label(combined, polya_base or _dpe_anchor_polya_base(rec))
                remote_mei_x0, remote_mei_x1 = _polya_tail_span_x(layout, combined)
                mei_start = int(layout.mei_3p) + 1
                mei_end = int(layout.mei_3p) + int(combined)
                aw = combined
            elif layout.polya_zone_bp > 0:
                aw = min(aw, int(layout.polya_zone_bp))
            if _dpe_anchor_abuts_polya_zone(anchor_side, layout):
                anchor_polya_x0, anchor_polya_x1 = _polya_span_from_flank_x(layout, aw)
            else:
                anchor_polya_x0, anchor_polya_x1 = _polya_tail_span_x(layout, aw)
            anchor_polya_clip = True
            anchor_polya_width = aw
            anchor_polya_base = _dpe_anchor_polya_base(rec)
            anchor_polya_label = _polya_tail_label(aw, anchor_polya_base)

    return {
        "anchor_pos": anchor_pos,
        "anchor_x": layout.genomic_to_x(anchor_pos),
        "remote_x": (remote_mei_x0 + remote_mei_x1) / 2.0,
        "remote_mei_x0": remote_mei_x0,
        "remote_mei_x1": remote_mei_x1,
        "mei_start": mei_start,
        "mei_end": mei_end,
        "side": anchor_side,
        "interchrom": mate_chrom not in {"", "*", layout.chrom},
        "mei_mapped": True,
        "read_name": str(rec.read_name),
        "evidence_type": "DPE",
        "remote_kind": remote_kind,
        "anchor_on_reference": True,
        "rescue_kind": rescue_kind,
        "polya_width": polya_width,
        "polya_base": polya_base,
        "polya_label": polya_label,
        "anchor_polya_clip": anchor_polya_clip,
        "anchor_polya_width": anchor_polya_width,
        "anchor_polya_base": anchor_polya_base,
        "anchor_polya_label": anchor_polya_label,
        "anchor_polya_x0": anchor_polya_x0,
        "anchor_polya_x1": anchor_polya_x1,
    }


def _default_supporting_reads_detail(gold_tsv: Path) -> Path | None:
    # Prefer TSV: parquet needs pyarrow/fastparquet, which may be absent.
    for name in ("supporting_reads_detail.mei.tsv", "supporting_reads_detail.mei.parquet"):
        candidate = gold_tsv.parent / name
        if candidate.exists():
            return candidate
    return None


def _load_detail_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        try:
            return pd.read_parquet(path)
        except ImportError:
            tsv = path.with_suffix(".tsv")
            if tsv.exists():
                return _read_tsv(tsv, usecols=_DETAIL_COLS)
            raise
    return _read_tsv(path, usecols=_DETAIL_COLS)


def _enrich_detail_with_clip_lens(
    detail: pd.DataFrame,
    gold_tsv: Path | None,
    sample: str,
    *,
    split_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Attach clip_len from split_evidence.{sample}.tsv when missing from detail."""
    if detail.empty:
        return detail
    out = detail.copy()
    if "clip_len" in out.columns and pd.to_numeric(out["clip_len"], errors="coerce").fillna(0).gt(0).any():
        return out

    if split_df is None and gold_tsv is not None:
        split_path = gold_tsv.parent / f"split_evidence.{sample}.tsv"
        if split_path.exists():
            split_df = pd.read_csv(
                split_path,
                sep="\t",
                usecols=lambda c: c in {"read_name", "clip_len", "clip_side", "pos"},
            )

    if split_df is None or split_df.empty or "clip_len" not in split_df.columns:
        out["clip_len"] = 0
        mei_span = (
            pd.to_numeric(out.get("mei_end", 0), errors="coerce").fillna(0).astype(int)
            - pd.to_numeric(out.get("mei_start", 0), errors="coerce").fillna(0).astype(int)
            + 1
        )
        sr_mask = out["evidence_type"].astype(str).eq("SR")
        out.loc[sr_mask, "clip_len"] = mei_span.loc[sr_mask].clip(lower=0)
        return out

    split_use = (
        split_df.sort_values("clip_len", ascending=False)
        .drop_duplicates("read_name", keep="first")
    )
    merged = out.merge(split_use[["read_name", "clip_len"]], on="read_name", how="left", suffixes=("", "_ev"))
    if "clip_len_ev" in merged.columns:
        merged["clip_len"] = pd.to_numeric(merged["clip_len_ev"], errors="coerce")
        merged = merged.drop(columns=["clip_len_ev"])
    mei_span = (
        pd.to_numeric(merged.get("mei_end", 0), errors="coerce").fillna(0).astype(int)
        - pd.to_numeric(merged.get("mei_start", 0), errors="coerce").fillna(0).astype(int)
        + 1
    )
    need = merged["evidence_type"].astype(str).eq("SR") & (
        pd.to_numeric(merged["clip_len"], errors="coerce").fillna(0) <= 0
    )
    merged.loc[need, "clip_len"] = mei_span.loc[need].clip(lower=0)
    merged["clip_len"] = pd.to_numeric(merged["clip_len"], errors="coerce").fillna(0).astype(int)
    return merged


def _dedupe_detail_rows(detail: pd.DataFrame) -> pd.DataFrame:
    """Keep one row per read_name + evidence_type (prefer MEI-mapped mates)."""
    if detail.empty:
        return detail
    out = detail.copy()

    def _flag(col: str) -> pd.Series:
        s = out.get(col, False)
        if not isinstance(s, pd.Series):
            return pd.Series(False, index=out.index)
        return s.fillna(False).infer_objects(copy=False).astype(bool)

    out["_mei_rank"] = (
        _flag("mate_mei_hit").astype(int) * 4
        + _flag("mei_hit").astype(int) * 3
        + (_flag("polya_rescue") | _flag("poly_tail_rescued")).astype(int) * 2
        + _flag("short_mei_seed_rescued").astype(int)
    )
    out = out.sort_values(["_mei_rank"], ascending=False)
    out = out.drop_duplicates(subset=["read_name", "evidence_type"], keep="first")
    return out.drop(columns=["_mei_rank"]).reset_index(drop=True)


def _detail_support_mask(detail: pd.DataFrame) -> pd.Series:
    """Rows eligible for read-architecture plots (MEI or polyA)."""
    if detail.empty:
        return pd.Series(dtype=bool)

    def _flag(col: str) -> pd.Series:
        s = detail.get(col, False)
        if not isinstance(s, pd.Series):
            return pd.Series(False, index=detail.index)
        return s.fillna(False).infer_objects(copy=False).astype(bool)

    poly_run = pd.Series(False, index=detail.index)
    if "clip_poly_at_run" in detail.columns:
        poly_run = pd.to_numeric(detail["clip_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
    return (
        _flag("mei_hit")
        | _flag("mate_mei_hit")
        | _flag("polya_rescue")
        | _flag("poly_tail_rescued")
        | _flag("vntr_rescue")
        | poly_run
    )


def _filter_detail_for_locus(
    detail: pd.DataFrame,
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    sample: str,
    breakpoint: int | None = None,
    neighbor_bp: int = 600,
) -> pd.DataFrame:
    """Select MEI / polyA supporting reads for one locus."""
    if detail.empty:
        return detail.iloc[0:0].copy()

    chrom_s = detail["chrom"].astype(str)
    sample_s = detail["sample"].astype(str)
    same_sample = (chrom_s == chrom) & (sample_s == sample)
    exact = detail.loc[
        same_sample
        & (pd.to_numeric(detail["window_start"], errors="coerce") == window_start)
        & (pd.to_numeric(detail["window_end"], errors="coerce") == window_end)
    ].copy()

    bp = int(breakpoint) if breakpoint is not None else (window_start + window_end) // 2
    gpos = pd.to_numeric(detail["genomic_pos"], errors="coerce")
    support = _detail_support_mask(detail)
    nearby = detail.loc[
        same_sample
        & (gpos >= bp - int(neighbor_bp))
        & (gpos <= bp + int(neighbor_bp))
        & support
    ].copy()
    return _dedupe_detail_rows(pd.concat([exact, nearby], ignore_index=True))


def _build_read_table_for_locus(
    *,
    chrom: str,
    window_start: int,
    window_end: int,
    sample: str,
    supporting_reads_detail: Path | None = None,
    breakpoint: int | None = None,
    neighbor_bp: int = 600,
    gold_tsv: Path | None = None,
    detail_df: pd.DataFrame | None = None,
    split_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Load MEI-supporting reads for the locus, including nearby windows."""
    if detail_df is None:
        if supporting_reads_detail is None:
            raise FileNotFoundError("supporting_reads_detail path or detail_df required")
        detail_df = _load_detail_table(supporting_reads_detail)
    combined = _filter_detail_for_locus(
        detail_df,
        chrom=chrom,
        window_start=window_start,
        window_end=window_end,
        sample=sample,
        breakpoint=breakpoint,
        neighbor_bp=neighbor_bp,
    )
    if gold_tsv is not None or split_df is not None:
        combined = _enrich_detail_with_clip_lens(combined, gold_tsv, sample, split_df=split_df)
    return combined


@dataclass
class ReadArchitectureCache:
    """In-memory tables shared across many locus plots."""

    gold_df: pd.DataFrame
    detail_df: pd.DataFrame
    mei_df: pd.DataFrame | None = None
    split_by_sample: dict[str, pd.DataFrame] = field(default_factory=dict)
    gold_tsv: Path | None = None
    _window_index: dict[tuple[str, str, int, int], object] = field(default_factory=dict, repr=False)
    _hit_gpos: dict[tuple[str, str], object] = field(default_factory=dict, repr=False)
    _hit_rows: dict[tuple[str, str], object] = field(default_factory=dict, repr=False)
    _mei_index: dict[tuple[str, int, int], pd.Series] = field(default_factory=dict, repr=False)
    _indexed: bool = field(default=False, repr=False)

    def ensure_indexes(self) -> None:
        """Build exact-window and MEI-hit indexes once for batch plotting."""
        if self._indexed:
            return
        self._mei_index = _build_mei_index(self.mei_df)
        if self.detail_df.empty:
            self._indexed = True
            return
        detail = self.detail_df
        sample_s = detail["sample"].astype(str).to_numpy()
        chrom_s = detail["chrom"].astype(str).to_numpy()
        ws = pd.to_numeric(detail["window_start"], errors="coerce").fillna(-1).astype(int).to_numpy()
        we = pd.to_numeric(detail["window_end"], errors="coerce").fillna(-1).astype(int).to_numpy()
        # groupby.indices returns ndarray positions — much cheaper than materializing frames.
        keys = pd.MultiIndex.from_arrays([sample_s, chrom_s, ws, we])
        self._window_index = {
            (str(k[0]), str(k[1]), int(k[2]), int(k[3])): idxs
            for k, idxs in detail.groupby(keys, sort=False).indices.items()
        }

        hit_mask = _detail_support_mask(detail).to_numpy()
        gpos = pd.to_numeric(detail["genomic_pos"], errors="coerce").fillna(-1).astype(int).to_numpy()
        hit_gpos: dict[tuple[str, str], object] = {}
        hit_rows: dict[tuple[str, str], object] = {}
        if hit_mask.any():
            hit_idx = np.flatnonzero(hit_mask)
            hit_sample = sample_s[hit_idx]
            hit_chrom = chrom_s[hit_idx]
            hit_pos = gpos[hit_idx]
            # Sort within each (sample, chrom) for searchsorted neighbor queries.
            order = np.lexsort((hit_pos, hit_chrom, hit_sample))
            hit_idx = hit_idx[order]
            hit_sample = hit_sample[order]
            hit_chrom = hit_chrom[order]
            hit_pos = hit_pos[order]
            # Split contiguous runs of (sample, chrom).
            if len(hit_idx):
                change = np.ones(len(hit_idx), dtype=bool)
                change[1:] = (hit_sample[1:] != hit_sample[:-1]) | (hit_chrom[1:] != hit_chrom[:-1])
                starts = np.flatnonzero(change)
                ends = np.append(starts[1:], len(hit_idx))
                for start, end in zip(starts.tolist(), ends.tolist()):
                    key = (str(hit_sample[start]), str(hit_chrom[start]))
                    hit_gpos[key] = hit_pos[start:end]
                    hit_rows[key] = hit_idx[start:end]
        self._hit_gpos = hit_gpos
        self._hit_rows = hit_rows
        self._indexed = True

    def detail_for_locus(
        self,
        *,
        chrom: str,
        window_start: int,
        window_end: int,
        sample: str,
        breakpoint: int | None = None,
        neighbor_bp: int = 600,
    ) -> pd.DataFrame:
        self.ensure_indexes()
        parts: list[pd.DataFrame] = []
        exact_idx = self._window_index.get((sample, chrom, int(window_start), int(window_end)))
        if exact_idx is not None and len(exact_idx):
            parts.append(self.detail_df.iloc[exact_idx])

        bp = int(breakpoint) if breakpoint is not None else (window_start + window_end) // 2
        key = (sample, chrom)
        gpos = self._hit_gpos.get(key)
        rows = self._hit_rows.get(key)
        if gpos is not None and rows is not None and len(gpos):
            lo = bp - int(neighbor_bp)
            hi = bp + int(neighbor_bp)
            left = int(np.searchsorted(gpos, lo, side="left"))
            right = int(np.searchsorted(gpos, hi, side="right"))
            if right > left:
                parts.append(self.detail_df.iloc[rows[left:right]])
        if not parts:
            return self.detail_df.iloc[0:0].copy()
        return _dedupe_detail_rows(pd.concat(parts, ignore_index=True))

    @classmethod
    def from_paths(
        cls,
        gold_tsv: Path,
        *,
        supporting_reads_detail: Path | None = None,
        mei_tsv: Path | None = None,
        load_split_evidence: bool = True,
    ) -> "ReadArchitectureCache":
        gold_df = _read_tsv(gold_tsv)
        detail_path = supporting_reads_detail or _default_supporting_reads_detail(gold_tsv)
        if detail_path is None or not detail_path.exists():
            raise FileNotFoundError(
                f"Missing supporting_reads_detail.mei.tsv beside {gold_tsv}. "
                "Re-run annotate-mei-support after candidate generation."
            )
        detail_df = _load_detail_table(detail_path)
        mei_path = mei_tsv if mei_tsv is not None else gold_tsv.parent / "candidate_loci.mei.tsv"
        mei_df = _load_slim_mei_table(mei_path) if mei_path.exists() else None
        split_by_sample: dict[str, pd.DataFrame] = {}
        if load_split_evidence:
            for sample in ("disease", "control"):
                split_path = gold_tsv.parent / f"split_evidence.{sample}.tsv"
                if split_path.exists():
                    split_by_sample[sample] = pd.read_csv(
                        split_path,
                        sep="\t",
                        usecols=lambda c: c in {"read_name", "clip_len", "clip_side", "pos"},
                    )
        cache = cls(
            gold_df=gold_df,
            detail_df=detail_df,
            mei_df=mei_df,
            split_by_sample=split_by_sample,
            gold_tsv=gold_tsv,
        )
        cache.ensure_indexes()
        return cache

    @classmethod
    def from_frames(
        cls,
        gold_df: pd.DataFrame,
        detail_df: pd.DataFrame,
        *,
        mei_df: pd.DataFrame | None = None,
        split_by_sample: dict[str, pd.DataFrame] | None = None,
        gold_tsv: Path | None = None,
    ) -> "ReadArchitectureCache":
        cache = cls(
            gold_df=gold_df,
            detail_df=detail_df,
            mei_df=mei_df,
            split_by_sample=dict(split_by_sample or {}),
            gold_tsv=gold_tsv,
        )
        cache.ensure_indexes()
        return cache


def _auto_flank_bp(detail: pd.DataFrame, layout_breakpoint: int, *, min_flank: int, max_flank: int = 500) -> int:
    """Widen flanks so mapped anchors are not all clamped to the junction."""
    if detail.empty:
        return int(min_flank)
    positions = pd.to_numeric(detail["genomic_pos"], errors="coerce").dropna().astype(int)
    if positions.empty:
        return int(min_flank)
    left_span = int((layout_breakpoint - positions[positions <= layout_breakpoint]).max()) if (positions <= layout_breakpoint).any() else 0
    right_span = int((positions[positions >= layout_breakpoint] - layout_breakpoint).max()) if (positions >= layout_breakpoint).any() else 0
    # Include room for a full 150 bp read beyond the farthest start.
    needed = max(left_span, right_span) + 150
    return int(max(min_flank, min(max_flank, needed)))


def _pair_segments(
    detail: pd.DataFrame,
    layout: LocusLayout,
    *,
    max_pairs: int,
    rng: random.Random | None = None,
    read_width_bp: int = 150,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    pairs: list[dict[str, object]] = []
    stats = {
        "detail_rows": len(detail),
        "sr_skipped": 0,
        "dpe_skipped": 0,
        "sr_plotted": 0,
        "dpe_plotted": 0,
        "polya_rescue_plotted": 0,
        "vntr_rescue_plotted": 0,
    }
    for rec in detail.itertuples(index=False):
        if str(rec.evidence_type) == "SR":
            pair = _pair_from_sr_row(rec, layout, read_width_bp=read_width_bp)
            if pair is None:
                stats["sr_skipped"] += 1
            else:
                stats["sr_plotted"] += 1
                pairs.append(pair)
                rk = str(pair.get("rescue_kind", "") or "")
                if rk == "polya":
                    stats["polya_rescue_plotted"] += 1
                elif rk == "vntr":
                    stats["vntr_rescue_plotted"] += 1
        else:
            pair = _pair_from_dpe_row(rec, layout)
            if pair is None:
                stats["dpe_skipped"] += 1
            else:
                stats["dpe_plotted"] += 1
                rk = str(pair.get("rescue_kind", "") or "")
                if rk == "polya":
                    stats["polya_rescue_plotted"] += 1
                elif rk == "vntr":
                    stats["vntr_rescue_plotted"] += 1
                pairs.append(pair)

    stats["pairs_before_cap"] = len(pairs)
    cap = max(1, int(max_pairs))
    sampler = rng if rng is not None else random.Random(0)
    if len(pairs) > cap:
        sr_pairs = [p for p in pairs if p.get("evidence_type") == "SR"]
        dpe_pairs = [p for p in pairs if p.get("evidence_type") != "SR"]

        def _sample_stratified(items: list[dict[str, object]], n: int) -> list[dict[str, object]]:
            """Keep a mix of MEI / polyA SR (and both clip sides) under the cap."""
            if n <= 0 or not items:
                return []
            if len(items) <= n:
                return list(items)
            buckets: dict[str, list[dict[str, object]]] = {
                "sr_mei": [],
                "sr_polya": [],
                "other": [],
            }
            for p in items:
                kind = str(p.get("remote_kind", "") or "other")
                buckets[kind if kind in buckets else "other"].append(p)
            # Prefer covering each evidence class and both L/R clip sides.
            chosen: list[dict[str, object]] = []
            chosen_ids: set[int] = set()

            def _take(pool: list[dict[str, object]], k: int) -> None:
                if k <= 0 or not pool:
                    return
                avail = [p for p in pool if id(p) not in chosen_ids]
                if not avail:
                    return
                pick = avail if len(avail) <= k else sampler.sample(avail, k)
                for p in pick:
                    chosen.append(p)
                    chosen_ids.add(id(p))

            # Guarantee at least one L and one R MEI SR when present.
            for side in ("L", "R"):
                side_mei = [
                    p
                    for p in buckets["sr_mei"]
                    if str(p.get("clip_side", "") or "").upper() == side
                ]
                _take(side_mei, 1)
            # Round-robin a few from each SR class, then fill randomly.
            per_class = max(1, (n - len(chosen)) // 2)
            for key in ("sr_mei", "sr_polya"):
                _take(buckets[key], per_class)
            leftover = [p for p in items if id(p) not in chosen_ids]
            _take(leftover, n - len(chosen))
            return chosen[:n]

        # Prefer keeping all split reads; reserve a small DPE slice only when SR
        # alone would consume the whole cap (common after polyA expansion).
        dpe_reserve = min(len(dpe_pairs), max(0, min(30, cap // 5))) if dpe_pairs else 0
        sr_budget = min(len(sr_pairs), max(0, cap - dpe_reserve))
        if len(sr_pairs) > sr_budget:
            keep_sr = _sample_stratified(sr_pairs, sr_budget)
        else:
            keep_sr = list(sr_pairs)
        remaining = max(0, cap - len(keep_sr))
        dpe_left = [p for p in dpe_pairs if p["side"] == "L"]
        dpe_right = [p for p in dpe_pairs if p["side"] == "R"]
        n_left = min(len(dpe_left), remaining // 2)
        n_right = min(len(dpe_right), remaining - n_left)
        # If one side is short, give leftover slots to the other.
        leftover = remaining - n_left - n_right
        if leftover > 0 and len(dpe_left) > n_left:
            extra = min(leftover, len(dpe_left) - n_left)
            n_left += extra
            leftover -= extra
        if leftover > 0 and len(dpe_right) > n_right:
            n_right += min(leftover, len(dpe_right) - n_right)
        chosen = list(keep_sr)
        if n_left:
            chosen.extend(sampler.sample(dpe_left, n_left))
        if n_right:
            chosen.extend(sampler.sample(dpe_right, n_right))
        # Unused DPE slots go back to SR.
        if len(chosen) < cap and len(sr_pairs) > len(keep_sr):
            keep_ids = {id(p) for p in keep_sr}
            extra_sr = _sample_stratified(
                [p for p in sr_pairs if id(p) not in keep_ids],
                cap - len(chosen),
            )
            chosen.extend(extra_sr)
        pairs = chosen
    # Stable visual order: left then right, SR before DPE, then genomic position.
    pairs.sort(
        key=lambda p: (
            0 if p["side"] == "L" else 1,
            0 if p.get("evidence_type") == "SR" else 1,
            int(p.get("anchor_pos", 0)),
            str(p.get("read_name", "")),
        )
    )
    stats["pairs_shown"] = len(pairs)
    stats["pairs_capped"] = max(0, stats["pairs_before_cap"] - len(pairs))
    return pairs, stats


def _anchor_fully_on_reference(genomic_pos: int, layout: LocusLayout) -> bool:
    """True when the anchor maps to reference flank without edge clamping."""
    pos = int(genomic_pos)
    if pos <= layout.breakpoint:
        return layout.breakpoint - pos <= layout.flank_bp
    return pos - layout.breakpoint <= layout.flank_bp


def _split_junction_x(layout: LocusLayout, side: str) -> float:
    if side == "L":
        return layout.insertion_left_x
    return layout.insertion_right_x


def _clip_span(x0: float, x1: float, x_min: float, x_max: float) -> tuple[float, float]:
    x0 = max(float(x_min), min(float(x0), float(x_max)))
    x1 = max(float(x_min), min(float(x1), float(x_max)))
    if x1 < x0:
        x0, x1 = x1, x0
    if x1 - x0 < 1.0:
        mid = (x0 + x1) / 2.0
        x0 = max(float(x_min), mid - 0.5)
        x1 = min(float(x_max), mid + 0.5)
    return x0, x1


def _genomic_read_bar(
    genomic_pos: int,
    *,
    side: str,
    width_bp: float,
    layout: LocusLayout,
) -> tuple[float, float]:
    """Reference bar at BAM-mapped genomic start; clipped at the breakpoint.

    Uses exclusive end (start + width) so the drawn x-span equals width_bp.
    """
    start = int(genomic_pos)
    width = int(width_bp)
    end_excl = start + width
    left_edge = layout.insertion_left_x
    right_edge = layout.insertion_right_x
    if side == "L":
        end_excl = min(end_excl, layout.breakpoint)
        if end_excl <= start:
            end_excl = start + 1
        x0 = layout.genomic_to_x(start)
        # Map exclusive end: last included base is end_excl-1, but keep width in x.
        x1 = x0 + (end_excl - start)
        x1 = min(x1, left_edge)
        return _clip_span(x0, x1, 0.0, left_edge)
    start = max(start, layout.breakpoint)
    end_excl = start + width
    x0 = layout.genomic_to_x(start)
    x1 = x0 + width
    return _clip_span(x0, x1, right_edge, float(layout.total_width))


def _split_segments(
    *,
    clip_side: str,
    clip_len: int,
    mapped_len: int,
    mei_x0: float,
    mei_x1: float,
    layout: LocusLayout,
    polya_clip: bool = False,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Black|orange split-read segments from real clip/mapped lengths + MEI coords.

    BAM: L softclip → mapped extends right of junction; R softclip → mapped ends at junction.
    Orange width is the real clip_len, centered on the mapped MEI hit (so a 52 bp clip
    at consensus 1327-1378 is drawn there, not as a generic stub). PolyA clips are
    placed in the dedicated polyA zone abutting MEI 3′.
    """
    clip_len = max(1, int(clip_len))
    mapped_len = max(1, int(mapped_len))
    left_edge = layout.insertion_left_x
    right_edge = layout.insertion_right_x
    if clip_side == "L":
        # Mapped genomic portion starts at BP and extends into the right flank.
        j_ref = right_edge
        ref_span = _clip_span(j_ref, j_ref + mapped_len, right_edge, float(layout.total_width))
    else:
        # Mapped genomic portion ends at BP on the left flank.
        j_ref = left_edge
        ref_span = _clip_span(j_ref - mapped_len, j_ref, 0.0, left_edge)

    if polya_clip and layout.polya_zone_bp > 0:
        return ref_span, layout.polya_span_x(clip_len)

    # Orange: clip_len bp centered on the mapped MEI / BRK stub interval.
    # Bound to the full insertion (MEI + polyA): BRK stubs often sit on the
    # polyA-side junction and collapse if clipped to the MEI body alone.
    mx0 = min(float(mei_x0), float(mei_x1))
    mx1 = max(float(mei_x0), float(mei_x1))
    center = (mx0 + mx1) / 2.0 if mx1 > mx0 else mx0
    half = clip_len / 2.0
    mei_span = _clip_span(center - half, center + half, left_edge, right_edge)
    # If centering clipped the width, fall back to the true mapped interval.
    if mei_span[1] - mei_span[0] < clip_len * 0.9 and mx1 > mx0:
        mei_span = _clip_span(mx0, mx1, left_edge, right_edge)
    return ref_span, mei_span

def _mei_read_bar(
    *,
    mei_x0: float,
    mei_x1: float,
    width_bp: float,
    layout: LocusLayout,
    anchor_side: str | None = None,
) -> tuple[float, float]:
    """Fixed-width MEI bar at the mapped consensus interval.

    Prefer growing the schematic ``width_bp`` bar *into* the insertion from the
    true mapped hit, rather than centering (centering a short 5' hit to 150 bp
    spills back to the junction and looks like insert-size 0).
    """
    mapped_x0 = min(float(mei_x0), float(mei_x1))
    mapped_x1 = max(float(mei_x0), float(mei_x1))
    if mapped_x1 <= mapped_x0:
        mapped_x1 = mapped_x0 + 1.0
    width = float(width_bp)
    side = (anchor_side or "").upper()
    # Left-flank mates usually hit near MEI 5' (left on forward axis): grow rightward.
    # Right-flank mates usually hit near MEI 3' (right): grow leftward.
    if side == "L":
        return _clip_span(
            mapped_x0,
            mapped_x0 + width,
            layout.mei_region_start_x,
            layout.mei_region_end_x,
        )
    if side == "R":
        return _clip_span(
            mapped_x1 - width,
            mapped_x1,
            layout.mei_region_start_x,
            layout.mei_region_end_x,
        )
    center = (mapped_x0 + mapped_x1) / 2.0
    return _clip_span(
        center - width / 2.0,
        center + width / 2.0,
        layout.mei_region_start_x,
        layout.mei_region_end_x,
    )


def _dpe_anchor_and_overhang(
    *,
    anchor_pos: int,
    side: str,
    width_bp: float,
    layout: LocusLayout,
    mate_mei_x0: float,
    mate_mei_x1: float,
) -> tuple[tuple[float, float], tuple[float, float] | None]:
    """Schematic width_bp DPE anchor at the true BAM start, clipped at the BP.

    Does **not** invent an MEI overhang to pad to 150 bp. Illumina DPE anchors
    that extend past the breakpoint are simply truncated on the reference; the
    mate's real MEI mapping is drawn separately. Invented overhangs were creating
    fake MEI–MEI pairs with implausible insert sizes.
    """
    del mate_mei_x0, mate_mei_x1
    width = int(width_bp)
    start = int(anchor_pos)
    left_edge = layout.insertion_left_x
    right_edge = layout.insertion_right_x
    if side == "L":
        end_excl = min(start + width, layout.breakpoint)
        if end_excl <= start:
            return (left_edge, left_edge), None
        x0 = layout.genomic_to_x(start)
        return _clip_span(x0, x0 + (end_excl - start), 0.0, left_edge), None
    start = max(start, layout.breakpoint)
    x0 = layout.genomic_to_x(start)
    return _clip_span(x0, x0 + width, right_edge, float(layout.total_width)), None


def _connector_endpoints(span_a: tuple[float, float], span_b: tuple[float, float]) -> tuple[float, float]:
    """Nearest edges between two bars (for the mate-pair connector)."""
    a0, a1 = span_a
    b0, b1 = span_b
    if a1 <= b0:
        return a1, b0
    if b1 <= a0:
        return a0, b1
    return (a0 + a1) / 2.0, (b0 + b1) / 2.0


def _draw_bar_span(
    ax,
    *,
    x0: float,
    x1: float,
    yi: float,
    color: str,
    height: float = 0.45,
) -> None:
    if x1 <= x0:
        return
    # Match SR/DPE bar thickness. Skip contrasting edges: white strokes on black
    # bars change perceived height under antialiasing on short/long spans.
    rect = Rectangle(
        (x0, yi - height / 2.0),
        x1 - x0,
        height,
        facecolor=color,
        edgecolor="none",
        linewidth=0.0,
        alpha=1.0,
        zorder=3,
    )
    ax.add_patch(rect)


def _draw_connector(ax, *, x0: float, x1: float, yi: float) -> None:
    if abs(x1 - x0) < 0.5:
        return
    ax.plot([x0, x1], [yi, yi], color=COLOR_BLACK, lw=1.0, alpha=0.8, zorder=2)


def _draw_read_pair(
    ax,
    *,
    yi: float,
    pair: dict[str, object],
    layout: LocusLayout,
    read_width_bp: float,
) -> None:
    remote_kind = str(pair.get("remote_kind", ""))
    is_split = bool(pair.get("is_split", False)) or remote_kind.startswith("sr_")

    if is_split and remote_kind in {"sr_mei", "sr_polya"}:
        clip_side = str(pair.get("clip_side", "L"))
        clip_len = int(pair.get("clip_len", 0) or 0)
        mapped_len = int(pair.get("mapped_len", 0) or 0)
        if mapped_len <= 0:
            mapped_len = max(1, int(read_width_bp) - max(1, clip_len))
        polya_clip = remote_kind == "sr_polya" or str(pair.get("rescue_kind", "")) == "polya"
        ref_span, mei_span = _split_segments(
            clip_side=clip_side,
            clip_len=clip_len,
            mapped_len=mapped_len,
            mei_x0=float(pair.get("remote_mei_x0", 0)),
            mei_x1=float(pair.get("remote_mei_x1", 0)),
            layout=layout,
            polya_clip=polya_clip,
        )
        mei_color = COLOR_LIGHT_ORANGE if polya_clip else COLOR_DARK_ORANGE
        _draw_bar_span(ax, x0=ref_span[0], x1=ref_span[1], yi=yi, color=COLOR_BLACK)
        _draw_bar_span(ax, x0=mei_span[0], x1=mei_span[1], yi=yi, color=mei_color)
        if polya_clip:
            label = str(pair.get("polya_label", "") or "")
            if not label:
                label = _polya_tail_label(
                    int(pair.get("polya_width", 0) or max(1, int(mei_span[1] - mei_span[0]))),
                    str(pair.get("polya_base", "A") or "A"),
                )
            _draw_polya_label(ax, x0=mei_span[0], x1=mei_span[1], yi=yi, label=label)
        # Tick at the genomic junction the mapped portion abuts.
        j = layout.insertion_right_x if clip_side == "L" else layout.insertion_left_x
        ax.plot([j, j], [yi - 0.32, yi + 0.32], color=COLOR_BLACK, lw=2.2, zorder=4)
        c0, c1 = _connector_endpoints(mei_span, ref_span)
        _draw_connector(ax, x0=c0, x1=c1, yi=yi)
        return

    side = str(pair["side"])
    mate_mei = None
    if remote_kind in {"mate_mei", "anchor_mei"}:
        mate_mei = (
            float(pair.get("remote_mei_x0", 0)),
            float(pair.get("remote_mei_x1", 0)),
        )
    ref_span, overhang = _dpe_anchor_and_overhang(
        anchor_pos=int(pair.get("anchor_pos", 0)),
        side=side,
        width_bp=read_width_bp,
        layout=layout,
        mate_mei_x0=float(pair.get("remote_mei_x0", 0)),
        mate_mei_x1=float(pair.get("remote_mei_x1", 0)),
    )
    _draw_bar_span(ax, x0=ref_span[0], x1=ref_span[1], yi=yi, color=COLOR_BLACK)

    # Anchor softclip into polyA: orange from flank|polyA boundary into the polyA zone
    # (SR-style), even when the mate is also a polyA rescue.
    anchor_polya = bool(pair.get("anchor_polya_clip"))
    anchor_polya_span: tuple[float, float] | None = None
    mate_is_polya = str(pair.get("rescue_kind", "")) == "polya"
    if anchor_polya and layout.polya_zone_bp > 0:
        aw = int(pair.get("anchor_polya_width", 0) or 0)
        if aw <= 0:
            aw = max(1, int(round(float(pair.get("anchor_polya_x1", 0)) - float(pair.get("anchor_polya_x0", 0)))))
        # Junction polyA clip + polyA mate → sum (minimum non-overlapping tail).
        if mate_is_polya:
            mate_w = int(pair.get("polya_width", 0) or 0)
            # polya_width may already be the combined sum from pair-building;
            # if anchor_polya_width was also set to that sum, don't double-add.
            if mate_w > aw:
                aw = mate_w
            elif aw == mate_w:
                pass
            else:
                # If pair already combined, polya_width == anchor_polya_width == sum.
                aw = max(aw, mate_w)
        aw = min(max(1, aw), int(layout.polya_zone_bp))
        if _dpe_anchor_abuts_polya_zone(side, layout):
            anchor_polya_span = _polya_span_from_flank_x(layout, aw)
            j = layout.insertion_left_x if layout.reverse_oriented else layout.insertion_right_x
        else:
            anchor_polya_span = layout.polya_span_x(aw)
            j = layout.insertion_left_x if side == "L" else layout.insertion_right_x
        _draw_bar_span(
            ax,
            x0=anchor_polya_span[0],
            x1=anchor_polya_span[1],
            yi=yi,
            color=COLOR_LIGHT_ORANGE,
        )
        label = str(pair.get("polya_label", "") or "") or str(pair.get("anchor_polya_label", "") or "")
        if not label:
            label = _polya_tail_label(aw, str(pair.get("anchor_polya_base", "A") or "A"))
        _draw_polya_label(ax, x0=anchor_polya_span[0], x1=anchor_polya_span[1], yi=yi, label=label)
        ax.plot([j, j], [yi - 0.32, yi + 0.32], color=COLOR_BLACK, lw=2.2, zorder=4)

    mei_color = COLOR_LIGHT_ORANGE if str(pair.get("rescue_kind", "")) else COLOR_DARK_ORANGE
    mei_span = None
    if mate_mei is not None:
        if mate_is_polya:
            if anchor_polya_span is not None:
                # Already drew the polyA zone bar from the anchor clip (merged width).
                mei_span = anchor_polya_span
            else:
                # Mate-only polyA rescue (no anchor polyA clip).
                mei_span = (
                    min(float(mate_mei[0]), float(mate_mei[1])),
                    max(float(mate_mei[0]), float(mate_mei[1])),
                )
                _draw_bar_span(ax, x0=mei_span[0], x1=mei_span[1], yi=yi, color=mei_color)
                label = str(pair.get("polya_label", "") or "")
                if not label:
                    label = _polya_tail_label(
                        int(pair.get("polya_width", 0) or max(1, int(mei_span[1] - mei_span[0]))),
                        str(pair.get("polya_base", "A") or "A"),
                    )
                _draw_polya_label(ax, x0=mei_span[0], x1=mei_span[1], yi=yi, label=label)
        else:
            mei_span = _mei_read_bar(
                mei_x0=mate_mei[0],
                mei_x1=mate_mei[1],
                width_bp=read_width_bp,
                layout=layout,
                anchor_side=side,
            )
            _draw_bar_span(ax, x0=mei_span[0], x1=mei_span[1], yi=yi, color=mei_color)

    if overhang is not None:
        _draw_bar_span(ax, x0=overhang[0], x1=overhang[1], yi=yi, color=COLOR_DARK_ORANGE)
        # Connect ref -> overhang (contiguous at junction) then overhang/mate.
        c0, c1 = _connector_endpoints(ref_span, overhang)
        _draw_connector(ax, x0=c0, x1=c1, yi=yi)
        if mei_span is not None:
            c0, c1 = _connector_endpoints(overhang, mei_span)
            _draw_connector(ax, x0=c0, x1=c1, yi=yi)
    elif anchor_polya_span is not None:
        c0, c1 = _connector_endpoints(ref_span, anchor_polya_span)
        _draw_connector(ax, x0=c0, x1=c1, yi=yi)
        # If mate is a real MEI body hit (not polyA), also connect polyA -> mate.
        if mei_span is not None and str(pair.get("rescue_kind", "")) != "polya":
            c0, c1 = _connector_endpoints(anchor_polya_span, mei_span)
            _draw_connector(ax, x0=c0, x1=c1, yi=yi)
    elif mei_span is not None:
        c0, c1 = _connector_endpoints(ref_span, mei_span)
        _draw_connector(ax, x0=c0, x1=c1, yi=yi)

    if pair.get("mate_in_mei"):
        mate_mei_span = _mei_read_bar(
            mei_x0=float(pair.get("mate_mei_x0", 0)),
            mei_x1=float(pair.get("mate_mei_x1", 0)),
            width_bp=read_width_bp,
            layout=layout,
        )
        _draw_bar_span(ax, x0=mate_mei_span[0], x1=mate_mei_span[1], yi=yi, color=COLOR_DARK_ORANGE)


def _mei_axis_ticks(layout: LocusLayout) -> tuple[list[float], list[str]]:
    """Tick positions/labels along the MEI segment (orientation-aware)."""
    left_coord = layout.mei_3p if layout.reverse_oriented else layout.mei_5p
    right_coord = layout.mei_5p if layout.reverse_oriented else layout.mei_3p
    ticks = [layout.mei_region_start_x, layout.mei_region_end_x]
    labels = [str(left_coord), str(right_coord)]
    if layout.mei_span_bp >= 300:
        mid_coord = (layout.mei_5p + layout.mei_3p) // 2
        mid_x = layout.mei_coord_to_x(mid_coord)
        ticks.insert(1, mid_x)
        labels.insert(1, str(mid_coord))
    return ticks, labels


def _safe_plot_stem(
    chrom: str,
    window_start: int,
    window_end: int,
    *,
    sample: str,
    rank: int | None = None,
) -> str:
    """Return a filesystem-safe read-architecture plot filename stem.

    Special characters in *chrom* and *sample* are replaced with underscores.
    When sanitisation changes *chrom*, a 6-hex-character SHA-256 digest of
    the original string is embedded between the sanitised name and the
    coordinates to prevent distinct chromosome names that sanitise identically
    (e.g. ``chr1/a`` and ``chr1_a``) from producing the same filename and
    silently overwriting each other.  Canonical names such as ``chr22`` are
    unaffected.
    """
    chrom_str = str(chrom)
    safe_chrom = re.sub(r"[^A-Za-z0-9_.-]+", "_", chrom_str)
    if safe_chrom != chrom_str:
        digest = hashlib.sha256(chrom_str.encode("utf-8")).hexdigest()[:6]
        safe_chrom = f"{safe_chrom}_{digest}"
    safe_sample = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(sample))
    body = f"read_arch_{safe_sample}_{safe_chrom}_{int(window_start)}_{int(window_end)}"
    if rank is None:
        return body
    return f"rank{int(rank):03d}_{body}"


def _select_gold_rows(
    gold_df: pd.DataFrame,
    *,
    gold_only: bool = True,
    top_n: int = 0,
) -> pd.DataFrame:
    if gold_df.empty:
        return gold_df.iloc[0:0].copy()
    subset = gold_df
    if gold_only and "analysis_stage_tier" in subset.columns:
        subset = subset.loc[subset["analysis_stage_tier"].fillna("").astype(str).str.lower() == "gold"]
    if int(top_n) > 0:
        return subset.head(int(top_n)).copy()
    return subset.copy()


def plot_locus_architecture(
    *,
    gold_tsv: Path | None = None,
    chrom: str,
    pos: int,
    sample: str = "auto",
    out_png: Path,
    supporting_reads_detail: Path | None = None,
    max_pairs: int = 50,
    flank_bp: int = 200,
    read_width_bp: int = 150,
    seed: int = 0,
    cache: ReadArchitectureCache | None = None,
    row: pd.Series | None = None,
) -> tuple[Path, pd.DataFrame]:
    if cache is not None:
        gold_df = cache.gold_df
        detail_df = cache.detail_df
        mei_df = cache.mei_df
        gold_path = cache.gold_tsv or gold_tsv
    else:
        if gold_tsv is None:
            raise ValueError("gold_tsv is required when cache is not provided")
        gold_df = None
        detail_df = None
        mei_df = None
        gold_path = gold_tsv

    if row is None:
        if gold_df is not None:
            row = _select_locus_row(gold_df, chrom, pos)
            row = _enrich_row_from_mei(
                row,
                mei_df,
                mei_index=cache._mei_index if cache is not None else None,
            )
        else:
            assert gold_path is not None
            row = _load_locus_row(gold_path, chrom, pos, mei_df=mei_df)
    elif mei_df is not None or (cache is not None and cache._mei_index):
        row = _enrich_row_from_mei(
            row,
            mei_df,
            mei_index=cache._mei_index if cache is not None else None,
        )

    sample = _choose_sample(row, sample)
    detail_path = supporting_reads_detail
    if detail_path is None and gold_path is not None:
        detail_path = _default_supporting_reads_detail(gold_path)
    if detail_df is None and (detail_path is None or not detail_path.exists()):
        raise FileNotFoundError(
            f"Missing supporting_reads_detail.mei.tsv beside {gold_path}. "
            "Re-run annotate-mei-support after candidate generation."
        )

    bp_l, bp_r, bp = _breakpoint_interval(row)
    discovery_start = _row_int(row, "discovery_window_start")
    discovery_end = _row_int(row, "discovery_window_end")
    if discovery_start <= 0 or discovery_end <= 0:
        discovery_start = int(row["window_start"])
        discovery_end = int(row["window_end"])
    split_df = None
    if cache is not None:
        split_df = cache.split_by_sample.get(sample)
    if cache is not None:
        detail = cache.detail_for_locus(
            chrom=str(row["chrom"]),
            window_start=discovery_start,
            window_end=discovery_end,
            sample=sample,
            breakpoint=bp,
        )
        if gold_path is not None or split_df is not None:
            detail = _enrich_detail_with_clip_lens(detail, gold_path, sample, split_df=split_df)
    else:
        detail = _build_read_table_for_locus(
            chrom=str(row["chrom"]),
            window_start=discovery_start,
            window_end=discovery_end,
            sample=sample,
            supporting_reads_detail=detail_path,
            breakpoint=bp,
            gold_tsv=gold_path,
            detail_df=detail_df,
            split_df=split_df,
        )
    flank_bp = _auto_flank_bp(detail, bp, min_flank=flank_bp)
    layout = _layout_from_row(row, sample=sample, detail=detail, flank_bp=flank_bp)
    pairs, pair_stats = _pair_segments(
        detail,
        layout,
        max_pairs=max_pairs,
        rng=random.Random(seed),
        read_width_bp=read_width_bp,
    )

    insert_hint = ""
    if layout.insert_size_estimates:
        insert_hint = f"  genomic_insert_sizes={list(layout.insert_size_estimates)}"

    # Plot reads from one BAM (disease or control); label by call status.
    status_label = _sample_status_label(row) or sample
    support_str = str(row.get(f"{sample}_supporting_reads", "") or "")
    ori_label = "reverse (−)" if layout.reverse_oriented else "forward (+)"

    title_lines = [
        (
            f"SAMPLE: {status_label}  |  "
            f"{layout.chrom}:{layout.window_start}-{layout.window_end}  "
            f"BP={layout.breakpoint_left}-{layout.breakpoint_right}  "
            f"{row.get('consensus_mei_family', '')}/{row.get('consensus_mei_subfamily', '')}  "
            f"tier={row.get('analysis_stage_tier', '')}"
        ),
        (
            f"{status_label} support: {support_str}  "
            f"mei={layout.mei_5p}-{layout.mei_3p} ({layout.mei_span_bp} bp, {layout.span_source}, "
            f"ori={layout.orientation} {ori_label})"
            f"{insert_hint}"
        ),
        (
            f"pairs_shown={pair_stats['pairs_shown']}/{pair_stats['pairs_before_cap']} "
            f"(SR={pair_stats['sr_plotted']}, DPE={pair_stats['dpe_plotted']}; "
            f"polyA_rescue={pair_stats.get('polya_rescue_plotted', 0)}, "
            f"VNTR_rescue={pair_stats.get('vntr_rescue_plotted', 0)}; "
            f"skipped SR={pair_stats['sr_skipped']}, DPE={pair_stats['dpe_skipped']}; "
            f"detail_rows={pair_stats['detail_rows']})"
        ),
    ]
    # Half-width figure (~700px @ 100dpi): wrap conservatively — bold 8.4pt is
    # wider than monospace char estimates, so 95 chars still clips on the right.
    title_fs = 8.4
    wrap_width = 86
    wrapped_parts: list[str] = []
    for line in title_lines:
        wrapped_parts.extend(
            textwrap.wrap(line, width=wrap_width, break_long_words=True, break_on_hyphens=False) or [""]
        )
    wrapped_title = "\n".join(wrapped_parts)
    n_title_rows = max(1, len(wrapped_parts))

    # Reserve header/footer in inches so short (few-pair) figures still clear the
    # wrapped title + flank labels. Fractional top margins shrink too much at fig_h~3.
    title_line_in = (title_fs * 1.35) / 72.0
    flank_label_in = 0.50 if layout.polya_zone_bp > 0 else 0.42
    header_pad_in = 0.12
    header_in = n_title_rows * title_line_in + flank_label_in + header_pad_in
    bottom_in = 0.48
    axes_in = max(1.35, len(pairs) * 0.09 + 0.55)
    fig_h = header_in + axes_in + bottom_in
    fig, ax = plt.subplots(figsize=(7, fig_h))
    fig.patch.set_facecolor(COLOR_WHITE)
    ax.set_facecolor(COLOR_WHITE)
    ax.set_xlim(-read_width_bp / 2, layout.total_width + read_width_bp / 2)
    # Tight axes: only reads (+ small pad). Labels live above the axes box.
    ax.set_ylim(-1, max(1, len(pairs)) + 0.8)

    # Left flank | (polyA?) MEI (polyA?) | right flank
    ax.axvspan(0, layout.insertion_left_x, color=COLOR_WHITE, alpha=1.0, zorder=0)
    if layout.polya_zone_bp > 0:
        ax.axvspan(
            layout.polya_region_start_x,
            layout.polya_region_end_x,
            color=COLOR_LIGHT_ORANGE,
            alpha=0.28,
            zorder=0,
        )
    ax.axvspan(
        layout.mei_region_start_x,
        layout.mei_region_end_x,
        color=COLOR_LIGHT_ORANGE,
        alpha=0.55,
        zorder=0,
    )
    ax.axvspan(layout.insertion_right_x, layout.total_width, color=COLOR_WHITE, alpha=1.0, zorder=0)
    ax.axvline(layout.mei_region_start_x, color=COLOR_DARK_ORANGE, ls="--", lw=1.2, alpha=0.9)
    ax.axvline(layout.mei_region_end_x, color=COLOR_DARK_ORANGE, ls="--", lw=1.2, alpha=0.9)
    if layout.polya_zone_bp > 0:
        outer = layout.polya_region_start_x if layout.reverse_oriented else layout.polya_region_end_x
        ax.axvline(outer, color=COLOR_DARK_ORANGE, ls=":", lw=1.0, alpha=0.7)
    ax.axvline(layout.genomic_to_x(layout.breakpoint), color=COLOR_BLACK, ls="-", lw=2, alpha=1.0)

    for yi, pair in enumerate(reversed(pairs)):
        _draw_read_pair(ax, yi=yi, pair=pair, layout=layout, read_width_bp=read_width_bp)

    mei_ticks, mei_labels = _mei_axis_ticks(layout)
    ax.set_xticks(mei_ticks)
    ax.set_xticklabels(mei_labels, color=COLOR_BLACK)
    ax.tick_params(axis="x", colors=COLOR_BLACK)
    ax.set_yticks([])
    polya_note = (
        f"; polyA zone {layout.polya_zone_bp} bp beyond 3′" if layout.polya_zone_bp > 0 else ""
    )
    ax.set_xlabel(
        f"MEI insertion 5′–3′ ({layout.mei_5p}–{layout.mei_3p}, {layout.mei_span_bp} bp, {ori_label})"
        f"{polya_note}; "
        f"flanks are {layout.flank_bp} bp each side of breakpoint",
        color=COLOR_BLACK,
        fontsize=8,
    )

    top = 1.0 - header_in / fig_h
    bottom = bottom_in / fig_h
    fig.subplots_adjust(left=0.04, right=0.99, bottom=bottom, top=top)

    # Title in figure coordinates (above axes), not ax.set_title (which collides when top is tight).
    fig.text(
        0.02,
        0.995,
        wrapped_title,
        ha="left",
        va="top",
        fontsize=title_fs,
        fontweight="bold",
        color=COLOR_BLACK,
        linespacing=1.15,
    )

    # Region labels sit just above the axes frame (axes-fraction y > 1), not in the data area.
    label_trans = blended_transform_factory(ax.transData, ax.transAxes)
    left_lab = layout.mei_3p if layout.reverse_oriented else layout.mei_5p
    right_lab = layout.mei_5p if layout.reverse_oriented else layout.mei_3p
    left_end = "3′" if layout.reverse_oriented else "5′"
    right_end = "5′" if layout.reverse_oriented else "3′"
    region_labels: list[tuple[float, str]] = [
        (layout.insertion_left_x / 2.0, f"{layout.chrom}\nleft flank"),
        (
            layout.mei_region_start_x + layout.mei_span_bp / 2.0,
            f"MEI ({ori_label})\n{left_end} {left_lab}…{right_lab} {right_end}",
        ),
    ]
    if layout.polya_zone_bp > 0:
        region_labels.append(
            (
                (layout.polya_region_start_x + layout.polya_region_end_x) / 2.0,
                f"polyA\n{layout.polya_zone_bp} bp",
            )
        )
    region_labels.append(
        (
            layout.insertion_right_x + layout.flank_bp / 2.0,
            f"{layout.chrom}\nright flank",
        )
    )
    for x_data, label_text in region_labels:
        ax.text(
            x_data,
            1.02,
            label_text,
            transform=label_trans,
            ha="center",
            va="bottom",
            fontsize=8,
            color=COLOR_BLACK,
            clip_on=False,
        )

    # Call-status badge inside the plot (corner): shared / disease_only / control_only.
    ax.text(
        0.01,
        0.98,
        status_label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=COLOR_BLACK,
        bbox=dict(boxstyle="round,pad=0.25", facecolor=COLOR_LIGHT_ORANGE, edgecolor=COLOR_BLACK, linewidth=1.0),
        zorder=10,
    )

    for spine in ax.spines.values():
        spine.set_color(COLOR_BLACK)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(out_png, dpi=100)
        return out_png, detail
    finally:
        plt.close(fig)


def generate_gold_read_architecture_plots(
    gold_review: pd.DataFrame,
    *,
    supporting_reads_detail: pd.DataFrame | Path,
    out_dir: Path,
    mei_table: pd.DataFrame | Path | None = None,
    gold_tsv: Path | None = None,
    gold_only: bool = True,
    top_n: int = 0,
    sample: str = "auto",
    max_pairs: int = 50,
    flank_bp: int = 200,
    read_width_bp: int = 150,
    seed: int = 0,
    progress_every: int = 50,
    cache: ReadArchitectureCache | None = None,
) -> pd.DataFrame:
    """Plot read architecture for gold-tier loci (batch; tables loaded once)."""
    variants = _select_gold_rows(gold_review, gold_only=gold_only, top_n=top_n)
    if variants.empty:
        click.echo("[read-arch] no variants selected for plots; skipping")
        return pd.DataFrame()

    if cache is None:
        if isinstance(supporting_reads_detail, Path):
            detail_df = _load_detail_table(supporting_reads_detail)
        else:
            detail_df = supporting_reads_detail

        if isinstance(mei_table, Path):
            mei_df = _load_slim_mei_table(mei_table) if mei_table.exists() else None
        elif mei_table is not None and not mei_table.empty:
            keep = [c for c in _MEI_ENRICH_COLS if c in mei_table.columns]
            mei_df = mei_table.loc[:, keep].copy() if keep else None
        else:
            mei_df = None

        split_by_sample: dict[str, pd.DataFrame] = {}
        if gold_tsv is not None:
            for samp in ("disease", "control"):
                split_path = gold_tsv.parent / f"split_evidence.{samp}.tsv"
                if split_path.exists():
                    split_by_sample[samp] = pd.read_csv(
                        split_path,
                        sep="\t",
                        usecols=lambda c: c in {"read_name", "clip_len", "clip_side", "pos"},
                    )

        cache = ReadArchitectureCache.from_frames(
            gold_df=gold_review,
            detail_df=detail_df,
            mei_df=mei_df,
            split_by_sample=split_by_sample,
            gold_tsv=gold_tsv,
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    index_rows: list[dict[str, object]] = []
    t0 = time.monotonic()
    n = len(variants)
    click.echo(f"[read-arch] generating {n} plots in {out_dir}")
    for rank, (_, row) in enumerate(variants.iterrows(), start=1):
        chrom = str(row.get("chrom", "") or "")
        ws = int(row.get("window_start", 0) or 0)
        we = int(row.get("window_end", 0) or 0)
        bp = _row_int(row, "consensus_insertion_breakpoint_pos")
        if bp <= 0:
            bp = (ws + we) // 2 if ws > 0 and we >= ws else 0
        if not chrom or bp <= 0:
            continue
        chosen = _choose_sample(row, sample)
        stem = _safe_plot_stem(chrom, ws, we, sample=chosen, rank=rank)
        out_png = out_dir / f"{stem}.png"
        try:
            plot_locus_architecture(
                chrom=chrom,
                pos=bp,
                sample=chosen,
                out_png=out_png,
                max_pairs=max_pairs,
                flank_bp=flank_bp,
                read_width_bp=read_width_bp,
                seed=seed + rank,
                cache=cache,
                row=row,
                gold_tsv=gold_tsv,
            )
            status = "ok"
            err = ""
        except Exception as exc:  # noqa: BLE001 - keep batch going
            status = "error"
            err = str(exc)
            click.echo(f"[read-arch] failed {chrom}:{ws}-{we}: {exc}")
        index_rows.append(
            {
                "plot_rank": rank,
                "analysis_stage_tier": row.get("analysis_stage_tier", ""),
                "chrom": chrom,
                "window_start": ws,
                "window_end": we,
                "discovery_window_start": _row_int(row, "discovery_window_start"),
                "discovery_window_end": _row_int(row, "discovery_window_end"),
                "consensus_insertion_breakpoint_pos": bp,
                "sample": chosen,
                "sample_status_label": _sample_status_label(row),
                "consensus_mei_family": row.get("consensus_mei_family", ""),
                "consensus_mei_subfamily": row.get("consensus_mei_subfamily", ""),
                "png": str(out_png) if status == "ok" else "",
                "status": status,
                "error": err,
            }
        )
        if progress_every > 0 and (rank % progress_every == 0 or rank == n):
            elapsed = time.monotonic() - t0
            rate = rank / elapsed if elapsed > 0 else 0.0
            click.echo(
                f"[read-arch] {rank}/{n} plots "
                f"({rate:.1f}/s, elapsed={elapsed:.1f}s)"
            )

    index = pd.DataFrame(index_rows)
    index_path = out_dir / "read_architecture_index.tsv"
    index.to_csv(index_path, sep="\t", index=False)
    ok = int((index["status"] == "ok").sum()) if not index.empty else 0
    click.echo(
        f"[read-arch] wrote {ok}/{len(index)} plots to {out_dir} "
        f"(index={index_path}, elapsed={time.monotonic() - t0:.1f}s)"
    )
    return index


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gold-review-tsv", type=Path, required=True)
    parser.add_argument("--chrom", default=None, help="Required unless --all-gold.")
    parser.add_argument(
        "--pos",
        type=int,
        default=None,
        help="1-based position inside locus window (required unless --all-gold).",
    )
    parser.add_argument(
        "--all-gold",
        action="store_true",
        help="Plot all gold-tier loci (loads tables once).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=0,
        help="With --all-gold, limit to first N prioritized gold loci (<=0 = all).",
    )
    parser.add_argument(
        "--supporting-reads-detail",
        type=Path,
        default=None,
        help="Per-read MEI-remap detail TSV/parquet from annotate-mei-support (default: beside gold review TSV).",
    )
    parser.add_argument(
        "--sample",
        default="auto",
        help="disease, control, or auto (choose the sample with more MEI/flank support).",
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=None,
        help="Output PNG for a single locus (required unless --all-gold).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for --all-gold plots (default: <gold-dir>/read_architecture).",
    )
    parser.add_argument("--out-detail-tsv", type=Path, default=None)
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=50,
        help="Maximum supporting read pairs to plot (stratified sample if more).",
    )
    parser.add_argument("--seed", type=int, default=0, help="RNG seed for random pair sampling.")
    parser.add_argument("--flank-bp", type=int, default=200, help="Symmetric flank width on each side of breakpoint.")
    parser.add_argument("--read-width-bp", type=int, default=150, help="Rendered read width in bp.")
    args = parser.parse_args()

    t0 = time.time()
    cache = ReadArchitectureCache.from_paths(
        args.gold_review_tsv,
        supporting_reads_detail=args.supporting_reads_detail,
    )

    if args.all_gold:
        out_dir = args.out_dir or (args.gold_review_tsv.parent / "read_architecture")
        index = generate_gold_read_architecture_plots(
            cache.gold_df,
            supporting_reads_detail=cache.detail_df,
            out_dir=out_dir,
            mei_table=cache.mei_df,
            gold_tsv=args.gold_review_tsv,
            gold_only=True,
            top_n=args.top_n,
            sample=args.sample,
            max_pairs=args.max_pairs,
            flank_bp=args.flank_bp,
            read_width_bp=args.read_width_bp,
            seed=args.seed,
            cache=cache,
        )
        click.echo(str(out_dir))
        click.echo(f"plots={len(index)} elapsed_sec={time.time() - t0:.2f}")
        return

    if args.chrom is None or args.pos is None or args.out_png is None:
        parser.error("--chrom, --pos, and --out-png are required unless --all-gold")

    out, detail = plot_locus_architecture(
        gold_tsv=args.gold_review_tsv,
        chrom=args.chrom,
        pos=args.pos,
        sample=args.sample,
        out_png=args.out_png,
        supporting_reads_detail=args.supporting_reads_detail,
        max_pairs=args.max_pairs,
        flank_bp=args.flank_bp,
        read_width_bp=args.read_width_bp,
        seed=args.seed,
        cache=cache,
    )
    if args.out_detail_tsv is not None:
        args.out_detail_tsv.parent.mkdir(parents=True, exist_ok=True)
        detail.to_csv(args.out_detail_tsv, sep="\t", index=False)
    click.echo(str(out))
    click.echo(f"detail_rows={len(detail)} elapsed_sec={time.time() - t0:.2f}")


if __name__ == "__main__":
    main()
