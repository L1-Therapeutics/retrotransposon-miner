from __future__ import annotations

import gzip
import hashlib
import json
import math
import random
import re
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import click
import pysam
from intervaltree import IntervalTree

from retro_miner.igv_plots import generate_gold_review_igv_plots
from retro_miner.read_architecture import generate_gold_read_architecture_plots
from retro_miner.local_assembly import annotate_silver_with_local_assembly
from retro_miner.evidence_extract import _longest_soft_clip_from_read, _soft_clip_query_seq


@dataclass
class ClipAlignmentSummary:
    sample: str
    clip_count: int
    paf_hits: int


_MIN_MEI_ANCHOR_BP = 25
_MIN_POLYA_RUN_FOR_END_IMPUTE = 12
_MIN_MEI_ANCHOR_BP_RELAXED = 15
_MIN_REPORTABLE_MEI_SPAN_BP = 20


@dataclass(frozen=True)
class FragmentToFullMap:
    """Linear projection of a Dfam panel fragment onto a full-length consensus."""

    fragment_name: str
    fragment_length: int
    full_name: str
    full_length: int
    fragment_aln_start: int  # 1-based inclusive
    fragment_aln_end: int
    full_aln_start: int
    full_aln_end: int
    strand: str
    family: str = ""


def _fragment_name_keys(name: str) -> list[str]:
    text = str(name or "").strip()
    if not text:
        return []
    keys = [text]
    base = text.split("|", 1)[0]
    if base and base not in keys:
        keys.append(base)
    no_hash = base.split("#", 1)[0]
    if no_hash and no_hash not in keys:
        keys.append(no_hash)
    return keys


def _resolve_fragment_to_full_map_tsv(
    mei_fasta: Path,
    mei_full_fasta: Path | None = None,
) -> Path | None:
    """Locate ``mei_fragment_to_full_coords.tsv`` beside the MEI / full FASTA."""
    candidates: list[Path] = []
    if mei_full_fasta is not None:
        full = Path(mei_full_fasta)
        candidates.append(full.with_name("mei_fragment_to_full_coords.tsv"))
        candidates.append(full.parent / "mei_fragment_to_full_coords.tsv")
    mei_fasta = Path(mei_fasta)
    candidates.extend(
        [
            mei_fasta.parent.parent / "full_consensus" / "mei_fragment_to_full_coords.tsv",
            mei_fasta.parent / "full_consensus" / "mei_fragment_to_full_coords.tsv",
        ]
    )
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def _load_fragment_to_full_map(path: Path | None) -> dict[str, FragmentToFullMap]:
    """Load fragment→full projection table keyed by fragment name variants."""
    if path is None or not path.exists():
        return {}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return {}
    required = {
        "fragment_name",
        "full_name",
        "fragment_aln_start",
        "fragment_aln_end",
        "full_aln_start",
        "full_aln_end",
        "strand",
    }
    if df.empty or not required.issubset(set(df.columns)):
        return {}
    out: dict[str, FragmentToFullMap] = {}
    for rec in df.itertuples(index=False):
        frag = str(getattr(rec, "fragment_name", "") or "").strip()
        full = str(getattr(rec, "full_name", "") or "").strip()
        if not frag or not full or full.lower() == "nan":
            continue
        try:
            entry = FragmentToFullMap(
                fragment_name=frag,
                fragment_length=int(float(getattr(rec, "fragment_length", 0) or 0)),
                full_name=full,
                full_length=int(float(getattr(rec, "full_length", 0) or 0)),
                fragment_aln_start=int(float(getattr(rec, "fragment_aln_start", 0) or 0)),
                fragment_aln_end=int(float(getattr(rec, "fragment_aln_end", 0) or 0)),
                full_aln_start=int(float(getattr(rec, "full_aln_start", 0) or 0)),
                full_aln_end=int(float(getattr(rec, "full_aln_end", 0) or 0)),
                strand=str(getattr(rec, "strand", "+") or "+"),
                family=str(getattr(rec, "family", "") or ""),
            )
        except (TypeError, ValueError):
            continue
        if entry.fragment_aln_end < entry.fragment_aln_start or entry.full_aln_end < entry.full_aln_start:
            continue
        for key in _fragment_name_keys(frag):
            out.setdefault(key, entry)
    return out


def _project_panel_coords_to_full(
    start: int,
    end: int,
    fragment_name: str,
    frag_map: dict[str, FragmentToFullMap],
) -> tuple[int, int, str] | None:
    """Project 1-based inclusive panel coords onto the full-length consensus.

    Returns ``(full_start, full_end, full_name)`` or None if unmapped.
    """
    if not frag_map:
        return None
    entry = None
    for key in _fragment_name_keys(fragment_name):
        entry = frag_map.get(key)
        if entry is not None:
            break
    if entry is None:
        return None
    a = int(start)
    b = int(end)
    if a <= 0 and b <= 0:
        return None
    if b < a:
        a, b = b, a
    lo = max(a, entry.fragment_aln_start)
    hi = min(b, entry.fragment_aln_end)
    if hi < lo:
        return None

    def _one(pos: int) -> int:
        offset = pos - entry.fragment_aln_start
        if entry.strand == "-":
            return entry.full_aln_end - offset
        return entry.full_aln_start + offset

    full_a = _one(lo)
    full_b = _one(hi)
    return min(full_a, full_b), max(full_a, full_b), entry.full_name


def _row_int_field(row: pd.Series, col: str) -> int:
    if col not in row.index:
        return 0
    val = pd.to_numeric(row.get(col), errors="coerce")
    if pd.isna(val):
        return 0
    return int(val)


def _collect_panel_fragment_intervals(row: pd.Series) -> list[tuple[int, int, str]]:
    """Gather (start, end, fragment_name) intervals that live on panel axes.

    L1 5′-end and 3′-end hits use different short references; they must be
    projected separately before any min/max on the full-length axis.
    """
    intervals: list[tuple[int, int, str]] = []

    def _add(start: int, end: int, frag: str) -> None:
        frag = str(frag or "").strip()
        if not frag or start <= 0 or end <= 0:
            return
        lo, hi = (start, end) if end >= start else (end, start)
        intervals.append((int(lo), int(hi), frag))

    for sample in ("disease", "control"):
        for side in ("L", "R"):
            start = _row_int_field(row, f"{sample}_{side}_mei_start")
            end = _row_int_field(row, f"{sample}_{side}_mei_end")
            if start <= 0 or end <= 0:
                start = _row_int_field(row, f"{sample}_{side}_mei_start_x")
                end = _row_int_field(row, f"{sample}_{side}_mei_end_x")
            if start <= 0 or end <= 0:
                start = _row_int_field(row, f"{sample}_{side}_detail_mei_start")
                end = _row_int_field(row, f"{sample}_{side}_detail_mei_end")
            frag = str(row.get(f"{sample}_{side}_mei_subfamily", "") or "")
            _add(start, end, frag)
        for side in ("left", "right"):
            start = _row_int_field(row, f"{sample}_discordant_mei_{side}_target_start_min")
            end = _row_int_field(row, f"{sample}_discordant_mei_{side}_target_end_max")
            frag = str(row.get(f"{sample}_discordant_mei_{side}_subfamily", "") or "")
            _add(start, end, frag)
    return intervals


def _row_consensus_mei_family(row: pd.Series) -> str:
    """Preferred MEI family for full-axis projection (never cross-family)."""
    for col in (
        "consensus_mei_family",
        "mei_family",
        "asm_mei_family",
        "consensus_mei_subfamily",
        "mei_subfamily",
    ):
        if col not in row.index:
            continue
        fam = _normalize_mei_family_token(str(row.get(col, "") or ""))
        if fam:
            return fam
    return ""


def _full_axis_union_from_panel_fragments(
    row: pd.Series,
    frag_map: dict[str, FragmentToFullMap],
) -> tuple[int, int] | None:
    """Project panel-fragment intervals onto one full-consensus axis, then union.

    Never mix families (Alu tip + L1 body) or different ``full_name`` axes.
    When the call has a consensus family, only same-family projections count;
    among those, keep the ``full_name`` with the most intervals (then longest span).
    """
    if not frag_map:
        return None
    preferred_fam = _row_consensus_mei_family(row)
    by_full: dict[str, list[int]] = {}
    for start, end, frag in _collect_panel_fragment_intervals(row):
        frag_fam = _normalize_mei_family_token(frag)
        if preferred_fam and frag_fam and frag_fam != preferred_fam:
            continue
        projected = _project_panel_coords_to_full(start, end, frag, frag_map)
        if projected is None:
            continue
        full_start, full_end, full_name = projected
        full_fam = _normalize_mei_family_token(full_name)
        if preferred_fam and full_fam and full_fam != preferred_fam:
            continue
        if not preferred_fam and frag_fam and full_fam and frag_fam != full_fam:
            continue
        key = str(full_name or "").strip()
        if not key:
            continue
        by_full.setdefault(key, []).extend([int(full_start), int(full_end)])
    if not by_full:
        return None

    preferred_full = ""
    for col in ("consensus_mei_subfamily", "mei_subfamily"):
        if col not in row.index:
            continue
        sub = str(row.get(col, "") or "").strip()
        if not sub:
            continue
        base = sub.split("#", 1)[0]
        base = re.sub(r"_(3end|5end|orf2)$", "", base, flags=re.IGNORECASE)
        base = re.sub(r"_short_?$", "", base, flags=re.IGNORECASE)
        if base:
            preferred_full = f"{base}_full"
            break

    # Prefer consensus-subfamily full axis, then most intervals, then longest span.
    def _axis_score(name: str) -> tuple[int, int, int]:
        coords = by_full[name]
        n_intervals = len(coords) // 2
        span = int(max(coords)) - int(min(coords))
        t_base = name.split("#", 1)[0]
        prefer = 1 if preferred_full and (
            t_base == preferred_full or name.startswith(preferred_full + "#") or name.startswith(preferred_full)
        ) else 0
        return prefer, n_intervals, span

    best_full = max(by_full, key=_axis_score)
    coords = by_full[best_full]
    if len(coords) < 2:
        return None
    return int(min(coords)), int(max(coords))


def _reference_contig_aliases(chrom: str) -> list[str]:
    c = str(chrom or "").strip()
    if not c:
        return []
    out: list[str] = [c]
    low = c.lower()
    if low.startswith("chr"):
        bare = c[3:]
        if bare:
            out.append(bare)
            if bare.lower() == "m":
                out.extend(["MT", "mt", "M", "m"])
    else:
        out.append(f"chr{c}")
        if low in {"m", "mt"}:
            out.extend(["chrM", "chrm", "MT", "mt", "M"])
    if low == "chrm":
        out.extend(["MT", "mt", "M", "m"])
    elif low == "mt":
        out.extend(["chrM", "chrm", "M", "m"])
    dedup: list[str] = []
    seen: set[str] = set()
    for name in out:
        if name in seen:
            continue
        seen.add(name)
        dedup.append(name)
    return dedup


def _make_reference_fetcher(ref: pysam.FastaFile):
    ref_names = set(ref.references)
    resolved_cache: dict[str, str] = {}

    def fetch(chrom: str, start0: int, end0: int) -> str:
        c = str(chrom or "").strip()
        if not c:
            return ""
        resolved = resolved_cache.get(c, "")
        if not resolved:
            for cand in _reference_contig_aliases(c):
                if cand in ref_names:
                    resolved = cand
                    break
            resolved_cache[c] = resolved
        if not resolved:
            return ""
        try:
            return ref.fetch(resolved, int(start0), int(end0)).upper()
        except Exception:
            return ""

    return fetch


def _load_table(base_dir: Path, stem: str, sample: str) -> pd.DataFrame:
    parquet_path = base_dir / f"{stem}.{sample}.parquet"
    tsv_path = base_dir / f"{stem}.{sample}.tsv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    raise FileNotFoundError(f"Missing {stem} for sample={sample}")


def _family_from_target(target: str) -> str:
    t = target.upper()
    if "SVA" in t:
        return "SVA"
    if "ALU" in t:
        return "ALU"
    if "LINE1" in t or "L1" in t:
        return "LINE1"
    if "HERV" in t or "ERV" in t:
        return "ERV"
    return "OTHER"


# Same-chr DPE mates closer than this often land in adjacent reference MEIs and
# must not drive insertion family/subfamily or MEI_MAPPED. Interchromosomal mates
# are always OK. No upper bound: far same-chr mates can still be valid homologous
# MEI placements; deletion-vs-insertion is handled by other call logic.
_DPE_MEI_IDENTITY_MIN_SAME_CHR_BP = 1000

# Short soft-clips may count as SR only when consistent with a strict MEI SR seed
# (preferred) or a same-side DPE MEI seed (fallback).
_STRICT_MEI_CLIP_MIN_BP = 20
_SHORT_MEI_RESCUE_MIN_CLIP_BP = 12
_SHORT_MEI_RESCUE_MIN_ALN_BP = 10
_SHORT_MEI_RESCUE_MAX_COORD_GAP_BP = 5
_SHORT_MEI_DPE_MIN_SUPPORT = 2
_SHORT_MEI_DPE_MAX_AXIS_GAP_BP = 500
_SHORT_MEI_DPE_PROXIMAL_SLACK_BP = 80

# Minimum soft-clip length for DPE→MEI consensus remap. Shorter clips are too
# ambiguous among Alu/SVA/L1; full reference-matching read bodies must never be
# remapped (∼43% of the genome is MEI-derived).
_DPE_MEI_REMAP_MIN_CLIP_BP = 20
# Unclipped mates may still be fully inside the insertion / homologous MEI; allow
# full-mate remap only when no soft clip meets the threshold above.
_DPE_MEI_REMAP_MIN_FULL_MATE_BP = 30


def _discordant_anchor_mei_query_seq(row: pd.Series | object) -> str:
    """Soft-clipped insert-like bases from a DPE anchor; never the ref-matched body."""
    clip_seq = str(getattr(row, "soft_clip_seq", "") or "")
    if len(clip_seq) >= _DPE_MEI_REMAP_MIN_CLIP_BP:
        return clip_seq
    # Derive from stored side/len + full read_seq when soft_clip_seq is absent
    # (older evidence tables).
    derived = _soft_clip_query_seq(
        str(getattr(row, "read_seq", "") or ""),
        str(getattr(row, "soft_clip_side", "") or ""),
        int(getattr(row, "soft_clip_len", 0) or 0),
    )
    return derived if len(derived) >= _DPE_MEI_REMAP_MIN_CLIP_BP else ""


def _discordant_mate_mei_query_seq(row: pd.Series | object) -> str:
    """Mate query for MEI remap: soft clip when present, else full mate sequence.

    Clipped mates at the opposite junction must not contribute their
    reference-matched bases. Unclipped mates may lie entirely in the insertion
    (or a distant homologous MEI), so the full mate sequence is allowed.
    """
    mate_clip = str(getattr(row, "mate_soft_clip_seq", "") or "")
    derived = _soft_clip_query_seq(
        str(getattr(row, "mate_seq", "") or ""),
        str(getattr(row, "mate_soft_clip_side", "") or ""),
        int(getattr(row, "mate_soft_clip_len", 0) or 0),
    )
    best_clip = mate_clip if len(mate_clip) >= len(derived) else derived
    clip_len = max(int(getattr(row, "mate_soft_clip_len", 0) or 0), len(best_clip))
    if clip_len > 0:
        # Any soft clip ⇒ never remap the ref-matched body; clip must be long enough.
        return best_clip if len(best_clip) >= _DPE_MEI_REMAP_MIN_CLIP_BP else ""
    mate_seq = str(getattr(row, "mate_seq", "") or "")
    return mate_seq if len(mate_seq) >= _DPE_MEI_REMAP_MIN_FULL_MATE_BP else ""


def _discordant_mate_ok_for_mei_identity(df: pd.DataFrame) -> pd.Series:
    """True when the mate is interchromosomal or >1 kb away on the same chrom.

    Same-chr mates within 1 kb are treated as landing in a nearby reference MEI
    nest and are excluded from insertion identity and MEI_MAPPED support.
    """
    chrom = df["chrom"].fillna("").astype(str) if "chrom" in df.columns else pd.Series("", index=df.index)
    mate_chrom = (
        df["mate_chrom"].fillna("").astype(str)
        if "mate_chrom" in df.columns
        else pd.Series("", index=df.index)
    )
    pos = pd.to_numeric(df["pos"] if "pos" in df.columns else 0, errors="coerce").fillna(0).astype(int)
    mate_pos = pd.to_numeric(
        df["mate_pos"] if "mate_pos" in df.columns else 0,
        errors="coerce",
    ).fillna(0).astype(int)
    valid_mate = (mate_chrom != "") & (mate_chrom != "*") & mate_pos.gt(0)
    interchrom = valid_mate & (mate_chrom != chrom)
    same_chr_far = (
        valid_mate
        & (mate_chrom == chrom)
        & (mate_pos - pos).abs().gt(_DPE_MEI_IDENTITY_MIN_SAME_CHR_BP)
    )
    return interchrom | same_chr_far


def _discordant_rows_for_mei_mapped_support(mei_df: pd.DataFrame) -> pd.DataFrame:
    """DPE rows that count toward MEI_MAPPED.

    Requires a consensus MEI hit (already filtered by callers via ``_mei_rows_only``
    or equivalent) and a mate that is interchromosomal or >1 kb away — i.e. not a
    same-chr mate into a nearby reference MEI.
    """
    if mei_df is None or mei_df.empty:
        return mei_df.iloc[0:0].copy() if mei_df is not None else pd.DataFrame()
    return mei_df.loc[_discordant_mate_ok_for_mei_identity(mei_df)].copy()


# NOTE: short Alu fragments can remap into SVA's Alu-like 5' domain (~80% ID over
# ~300 bp), so SVA labels may appear at true Alu insertions. Family-first + pooled
# disease/control voting usually resolves this; a stricter SVA-domain gate is TBD.
def _discordant_rows_for_mei_identity(mei_df: pd.DataFrame) -> pd.DataFrame:
    """DPE rows eligible to vote on insertion family/subfamily.

    Requires a mate that is interchromosomal or >1 kb away (nearby same-chr mates
    into adjacent reference MEIs are excluded). Anchor-only consensus hits still
    vote; when a mate consensus hit exists, prefer mate labels over anchor labels.
    """
    if mei_df is None or mei_df.empty:
        return mei_df.iloc[0:0].copy() if mei_df is not None else pd.DataFrame()
    out = mei_df.loc[_discordant_mate_ok_for_mei_identity(mei_df)].copy()
    if out.empty:
        return out
    if "mate_mei_target" in out.columns:
        mate_t = out["mate_mei_target"].fillna("").astype(str)
        use = mate_t.str.len().gt(0)
        out.loc[use, "target"] = mate_t.loc[use]
    if "mate_mei_family" in out.columns:
        mate_f = out["mate_mei_family"].fillna("").astype(str)
        use = mate_f.str.len().gt(0)
        out.loc[use, "family"] = mate_f.loc[use]
    elif "target" in out.columns:
        out["family"] = out["target"].fillna("").astype(str).map(_family_from_target)
    return out


def _format_vote_map(counts: dict[str, int]) -> str:
    """Serialize ``{label: weight}`` as ``LABEL:N,LABEL:N`` (stable key order)."""
    parts = []
    for key in sorted(counts):
        weight = int(counts[key])
        label = str(key or "").strip()
        if label and weight > 0:
            parts.append(f"{label}:{weight}")
    return ",".join(parts)


def _parse_vote_map(text: object) -> dict[str, int]:
    """Parse ``LABEL:N,LABEL:N`` vote maps; invalid entries are ignored."""
    raw = str(text or "").strip()
    if not raw or raw.lower() in {"nan", "none"}:
        return {}
    out: dict[str, int] = {}
    for part in raw.split(","):
        part = part.strip()
        if not part or ":" not in part:
            continue
        label, weight_s = part.rsplit(":", 1)
        label = label.strip()
        try:
            weight = int(float(weight_s))
        except (TypeError, ValueError):
            continue
        if label and weight > 0:
            out[label] = out.get(label, 0) + weight
    return out


def _identity_vote_maps(
    identity_df: pd.DataFrame,
    *,
    sample_prefix: str,
) -> pd.DataFrame:
    """Per-locus unique-read vote maps for family and subfamily (identity rows)."""
    cols = [
        "chrom",
        "window_start",
        "window_end",
        f"{sample_prefix}_discordant_mei_family_votes",
        f"{sample_prefix}_discordant_mei_subfamily_votes",
    ]
    if identity_df is None or identity_df.empty:
        return pd.DataFrame(columns=cols)
    work = identity_df.copy()
    if "family" not in work.columns and "target" in work.columns:
        work["family"] = work["target"].fillna("").astype(str).map(_family_from_target)
    if "family" not in work.columns or "target" not in work.columns or "read_name" not in work.columns:
        return pd.DataFrame(columns=cols)
    work["family"] = work["family"].fillna("").astype(str).map(_normalize_mei_family_token)
    work["target"] = work["target"].fillna("").astype(str)
    work["read_name"] = work["read_name"].fillna("").astype(str)
    work = work.loc[
        work["family"].ne("")
        & work["target"].ne("")
        & work["read_name"].ne("")
    ].copy()
    if work.empty:
        return pd.DataFrame(columns=cols)

    rows: list[dict[str, object]] = []
    for (chrom, ws, we), grp in work.groupby(["chrom", "window_start", "window_end"], sort=False):
        fam_counts = grp.groupby("family", sort=False)["read_name"].nunique().astype(int).to_dict()
        sub_counts = grp.groupby("target", sort=False)["read_name"].nunique().astype(int).to_dict()
        rows.append(
            {
                "chrom": chrom,
                "window_start": ws,
                "window_end": we,
                f"{sample_prefix}_discordant_mei_family_votes": _format_vote_map(fam_counts),
                f"{sample_prefix}_discordant_mei_subfamily_votes": _format_vote_map(sub_counts),
            }
        )
    return pd.DataFrame(rows, columns=cols)


def _top_family_then_subfamily(
    df: pd.DataFrame,
    *,
    group_cols: list[str],
    family_col: str,
    subfamily_col: str,
    score_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pick family by pooled score, then subfamily only within that family.

    Returns ``(family_top, subfamily_top)`` with one row per group.
    """
    empty_fam = pd.DataFrame(columns=group_cols + [family_col, score_col])
    empty_sub = pd.DataFrame(columns=group_cols + [subfamily_col, score_col])
    if df is None or df.empty:
        return empty_fam, empty_sub
    work = df.copy()
    work[family_col] = work[family_col].fillna("").astype(str)
    work[subfamily_col] = work[subfamily_col].fillna("").astype(str)
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce").fillna(0.0)
    work = work.loc[
        work[family_col].ne("")
        & work[family_col].ne("OTHER")
        & work[subfamily_col].ne("")
        & work[score_col].gt(0)
    ].copy()
    if work.empty:
        return empty_fam, empty_sub

    ascending = [True] * len(group_cols) + [False]
    family_top = (
        work.groupby(group_cols + [family_col], as_index=False)[score_col]
        .sum()
        .sort_values(group_cols + [score_col], ascending=ascending)
        .drop_duplicates(group_cols, keep="first")
    )
    within = work.merge(family_top[group_cols + [family_col]], on=group_cols + [family_col], how="inner")
    subfamily_top = (
        within.groupby(group_cols + [subfamily_col], as_index=False)[score_col]
        .sum()
        .sort_values(group_cols + [score_col], ascending=ascending)
        .drop_duplicates(group_cols, keep="first")
    )
    return family_top, subfamily_top


def _resolve_full_consensus_fasta(
    mei_fasta: Path,
    out_dir: Path,
    mei_full_fasta: Path | None = None,
) -> Path | None:
    """Resolve full-consensus FASTA for coordinate-normalized remapping.

    Preference order:
      1. Explicit ``mei_full_fasta``
      2. Packaged full-length panel (true L1HS ~6 kb, AluY, SVA_F)
      3. Auto-build from longest/canonical panel records (avoid short 5'/3' stubs)
    """
    if mei_full_fasta is not None:
        full = Path(mei_full_fasta)
        return full if full.exists() else None

    mei_fasta = Path(mei_fasta)
    packaged = [
        # Prefer the full Dfam-matched panel (polyA-trimmed at prep) over the
        # tiny 3-sequence RepeatBrowser stub.
        mei_fasta.parent.parent / "full_consensus" / "mei_full_canonical.panel.fa",
        mei_fasta.parent.parent / "full_consensus" / "mei_full_canonical.ucsc_repeatbrowser.panel.fa",
        mei_fasta.parent.parent / "full_consensus" / "mei_full_canonical.ucsc_repeatbrowser.fa",
        mei_fasta.parent / "full_consensus" / "mei_full_canonical.panel.fa",
        mei_fasta.parent / "full_consensus" / "mei_full_canonical.ucsc_repeatbrowser.fa",
    ]
    for cand in packaged:
        if cand.exists() and cand.stat().st_size > 0:
            fai = cand.with_suffix(cand.suffix + ".fai")
            if not fai.exists():
                try:
                    subprocess.run(["samtools", "faidx", str(cand)], check=False, capture_output=True, text=True)
                except Exception:
                    pass
            return cand

    if not mei_fasta.exists():
        return None
    fai = mei_fasta.with_suffix(mei_fasta.suffix + ".fai")
    if not fai.exists():
        return None
    try:
        fai_df = pd.read_csv(
            fai,
            sep="\t",
            header=None,
            names=["name", "length", "offset", "line_bases", "line_width"],
            usecols=[0, 1],
        )
    except Exception:
        return None
    if fai_df.empty:
        return None
    fai_df["name"] = fai_df["name"].astype(str)
    fai_df["family"] = fai_df["name"].map(_family_from_target)
    fai_df["length"] = pd.to_numeric(fai_df["length"], errors="coerce").fillna(0).astype(int)

    def _pick_target(family: str, preferred: list[str]) -> str:
        fam = fai_df.loc[fai_df["family"] == family].copy()
        if fam.empty:
            return ""
        names = set(fam["name"].tolist())
        for cand_name in preferred:
            if cand_name in names:
                return cand_name
        fam = fam.sort_values("length", ascending=False)
        for _, row in fam.iterrows():
            name = str(row["name"])
            # Skip short fragment consensus records for LINE1 when possible.
            if family == "LINE1" and int(row["length"]) < 5000 and any(
                tag in name for tag in ("_5end", "_3end", "_orf2")
            ):
                continue
            return name
        return str(fam.iloc[0]["name"])

    picks = {
        "ALU": _pick_target("ALU", ["AluY#SINE/Alu", "AluYa5#SINE/Alu", "AluYb8#SINE/Alu"]),
        "SVA": _pick_target("SVA", ["SVA_F#Retroposon/SVA", "SVA_D#Retroposon/SVA"]),
        "LINE1": _pick_target(
            "LINE1",
            [
                "L1HS_full#LINE/L1",
                "L1HS#LINE/L1",
            ],
        ),
    }
    if not any(picks.values()):
        return None
    full_fa = out_dir / "mei_full_consensus.auto.fasta"
    with pysam.FastaFile(str(mei_fasta)) as ref, full_fa.open("w", encoding="utf-8") as out:
        for fam, target in picks.items():
            if not target:
                continue
            try:
                seq = ref.fetch(target)
            except Exception:
                continue
            if not seq:
                continue
            out.write(f">{fam}_FULL|source={target}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
    if not full_fa.exists() or full_fa.stat().st_size == 0:
        return None
    try:
        subprocess.run(["samtools", "faidx", str(full_fa)], check=False, capture_output=True, text=True)
    except Exception:
        pass
    return full_fa



_MEI_HIT_COLUMNS = [
    "qname",
    "target",
    "target_start",
    "target_end",
    "target_len",
    "target_strand",
    "alnlen",
    "mapq",
    "pid",
    "qcov",
    "mei_score",
    "family",
]

# Short MEI tips (~20–30 bp) need smaller seeds than bwa mem defaults (k=19, T=30).
# Empirically (scripts/benchmark_mei_aligners.py on Alu/LINE1/SVA consensus samples):
# bwa mem -k10 -T10 recovers ~100% of 20 bp+ tips; default minimap2 -x sr recovers 0% ≤30 bp.
# Primary alignments only (no -a): we keep the best hit per query anyway, and -a
# inflates SAM ~8× with secondaries that are discarded after scoring.
# Thread count is per invocation via --bwa-threads (default 1). The pipeline wrapper
# raises it for single-chrom runs and keeps 1 under multi-chrom concurrency.
def _bwa_mem_mei_args(*, bwa_threads: int = 1) -> tuple[str, ...]:
    threads = max(1, int(bwa_threads))
    return ("mem", "-t", str(threads), "-k", "10", "-T", "10")
# After polyA/T trim, length-aware hit gate (qcov vs trimmed query length):
# short tips must cover most of the query; longer clips often tip-align only.
_MEI_ALIGN_MIN_QCOV_SHORT = 0.80
_MEI_ALIGN_MIN_PID = 0.90
_MEI_ALIGN_MIN_TRIMMED_BP = 12
# Short path: qcov≥0.80 on trimmed len≥12 already implies ≳10 bp aligned.
_MEI_ALIGN_MIN_ALN_BP_LONG = 20
_MEI_ALIGN_SHORT_QLEN_MAX = 30


def _trim_poly_at_from_clip(seq: str, side: str = "") -> str:
    """Remove the longest A- or T-homopolymer run (≥8 bp) from a soft-clip.

    Used before MEI consensus remap so polyA+tip clips are scored on the tip.
    Prefers the non-poly segment consistent with clip side (L→3' of clip /
    right residual; R→5' of clip / left residual).
    """
    s = (seq or "").upper().replace("U", "T")
    if not s:
        return ""
    best_start = -1
    best_end = -1
    best_len = 0
    cur_start = -1
    cur_base = ""
    for i, ch in enumerate(s):
        if ch not in {"A", "T"}:
            if cur_start >= 0:
                run_len = i - cur_start
                if run_len > best_len:
                    best_start, best_end, best_len = cur_start, i, run_len
                cur_start = -1
                cur_base = ""
            continue
        if cur_start < 0:
            cur_start = i
            cur_base = ch
        elif ch != cur_base:
            run_len = i - cur_start
            if run_len > best_len:
                best_start, best_end, best_len = cur_start, i, run_len
            cur_start = i
            cur_base = ch
    if cur_start >= 0:
        run_len = len(s) - cur_start
        if run_len > best_len:
            best_start, best_end, best_len = cur_start, len(s), run_len
    if best_len < 8:
        return s
    left = s[:best_start]
    right = s[best_end:]
    side_u = (side or "").upper()
    preferred = right if side_u == "L" else left if side_u == "R" else (
        right if len(right) >= len(left) else left
    )
    backup = left if preferred is right else right
    if len(preferred) >= _MEI_ALIGN_MIN_TRIMMED_BP:
        return preferred
    if len(backup) >= _MEI_ALIGN_MIN_TRIMMED_BP:
        return backup
    joined = (left + right).strip()
    # Pure / near-pure polyA clips: do not fall back to the untrimmed sequence
    # (that remaps onto consensus A-tails and inflates MEI_MAPPED).
    return joined if len(joined) >= _MEI_ALIGN_MIN_TRIMMED_BP else ""


def _mei_align_quality_ok(
    df: pd.DataFrame,
    *,
    query_len: pd.Series | None = None,
) -> pd.Series:
    """Length-aware MEI-hit gate on trimmed-query alignment metrics.

    Short trimmed queries (≤30 bp): require high qcov + pid.
    Longer queries: tip alignments are allowed — require pid + min alnlen,
    not full-query coverage (long clips often have only an MEI tip).
    """
    if df is None or df.empty:
        return pd.Series(dtype=bool)
    has = (
        df["target"].notna() & df["target"].fillna("").astype(str).ne("")
        if "target" in df.columns
        else pd.Series(False, index=df.index)
    )
    qcov = (
        pd.to_numeric(df["qcov"], errors="coerce").fillna(0.0)
        if "qcov" in df.columns
        else pd.Series(0.0, index=df.index)
    )
    pid = (
        pd.to_numeric(df["pid"], errors="coerce").fillna(0.0)
        if "pid" in df.columns
        else pd.Series(0.0, index=df.index)
    )
    aln = (
        pd.to_numeric(df["alnlen"], errors="coerce").fillna(0.0)
        if "alnlen" in df.columns
        else pd.Series(0.0, index=df.index)
    )
    if query_len is None:
        # Recover approximate query length from q_aln / qcov when available.
        qlen = pd.Series(0.0, index=df.index)
        recoverable = qcov.gt(0)
        qlen.loc[recoverable] = (aln.loc[recoverable] / qcov.loc[recoverable]).astype(float)
    else:
        qlen = pd.to_numeric(query_len, errors="coerce").reindex(df.index).fillna(0.0)
    short = qlen.le(float(_MEI_ALIGN_SHORT_QLEN_MAX))
    short_ok = (
        short
        & qcov.ge(float(_MEI_ALIGN_MIN_QCOV_SHORT))
        & pid.ge(float(_MEI_ALIGN_MIN_PID))
    )
    long_ok = (
        (~short)
        & pid.ge(float(_MEI_ALIGN_MIN_PID))
        & aln.ge(float(_MEI_ALIGN_MIN_ALN_BP_LONG))
    )
    return has & (short_ok | long_ok)


def _apply_mei_align_quality_gate(
    df: pd.DataFrame,
    *,
    query_len: pd.Series | None = None,
) -> pd.DataFrame:
    """Clear MEI hit fields that fail the length-aware qcov/pid gate."""
    if df is None or df.empty or "target" not in df.columns:
        return df
    out = df.copy()
    ok = _mei_align_quality_ok(out, query_len=query_len)
    bad = out["target"].fillna("").astype(str).ne("") & ~ok
    if not bool(bad.any()):
        return out
    clear_empty = ("target", "family", "target_strand")
    clear_zero_int = ("target_start", "target_end", "target_len", "alnlen", "mapq")
    clear_zero_float = ("pid", "qcov", "mei_score")
    for col in clear_empty:
        if col in out.columns:
            out.loc[bad, col] = ""
    for col in clear_zero_int:
        if col in out.columns:
            out.loc[bad, col] = 0
    for col in clear_zero_float:
        if col in out.columns:
            out.loc[bad, col] = 0.0
    return out


def _empty_mei_hits() -> pd.DataFrame:
    return pd.DataFrame(columns=_MEI_HIT_COLUMNS)


def _mei_score_from_metrics(*, pid: float, qcov: float, mapq: int) -> float:
    raw_score = (0.45 * float(pid)) + (0.35 * float(qcov)) + (0.2 * (float(mapq) / 60.0))
    return float(max(0.0, min(1.0, raw_score)))


def _pick_best_mei_hits(rows: list[dict]) -> pd.DataFrame:
    if not rows:
        return _empty_mei_hits()
    hits = pd.DataFrame(rows)
    hits = hits.sort_values(
        ["qname", "mei_score", "alnlen", "mapq", "pid", "qcov"],
        ascending=[True, False, False, False, False, False],
    )
    return hits.drop_duplicates(subset=["qname"], keep="first")


def _best_hits_from_paf(paf_path: Path) -> pd.DataFrame:
    rows = []
    with paf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7])
            tend = int(parts[8])
            nmatch = int(parts[9])
            alnlen = int(parts[10])
            mapq = int(parts[11])
            qcov = (qend - qstart) / qlen if qlen > 0 else 0.0
            pid = (nmatch / alnlen) if alnlen > 0 else 0.0
            rows.append(
                {
                    "qname": qname,
                    "target": tname,
                    "target_start": tstart + 1,
                    "target_end": tend,
                    "target_len": tlen,
                    "target_strand": strand,
                    "alnlen": alnlen,
                    "mapq": mapq,
                    "pid": pid,
                    "qcov": qcov,
                    "mei_score": _mei_score_from_metrics(pid=pid, qcov=qcov, mapq=mapq),
                    "family": _family_from_target(tname),
                }
            )
    return _pick_best_mei_hits(rows)


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _cigar_alignment_spans(cigar: str) -> tuple[int, int, int]:
    """Return (query_aligned_bp, ref_span_bp, alnlen) from a CIGAR string."""
    q_aln = 0
    r_span = 0
    alnlen = 0
    for n_s, op in _CIGAR_RE.findall(cigar or ""):
        n = int(n_s)
        if op in {"M", "=", "X"}:
            q_aln += n
            r_span += n
            alnlen += n
        elif op == "I":
            q_aln += n
            alnlen += n
        elif op in {"D", "N"}:
            r_span += n
            alnlen += n
    return q_aln, r_span, alnlen


def _cigar_query_len(cigar: str) -> int:
    """Full query length implied by CIGAR (M/I/=/X/S/H)."""
    qlen = 0
    for n_s, op in _CIGAR_RE.findall(cigar or ""):
        if op in {"M", "I", "=", "X", "S", "H"}:
            qlen += int(n_s)
    return qlen


def _best_hits_from_sam(sam_text: str, *, target_lengths: dict[str, int] | None = None) -> pd.DataFrame:
    """Best *primary* SAM hit per query (same columns as PAF path).

    Skips secondary (0x100) and supplementary (0x800). Those often hard-clip SEQ
    to the aligned tip only, which falsely makes qcov look like 1.0 and beats the
    real soft-clipped primary when ranking by mei_score.
    """
    tlen_map = dict(target_lengths or {})
    rows: list[dict] = []
    for line in (sam_text or "").splitlines():
        if not line or line.startswith("@"):
            if line.startswith("@SQ"):
                # @SQ SN:name LN:len
                sn = ln = None
                for field in line.split("\t")[1:]:
                    if field.startswith("SN:"):
                        sn = field[3:]
                    elif field.startswith("LN:"):
                        ln = int(field[3:])
                if sn and ln:
                    tlen_map[sn] = ln
            continue
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        qname = parts[0]
        flag = int(parts[1])
        if flag & 0x4:  # unmapped
            continue
        if flag & 0x100 or flag & 0x800:  # secondary / supplementary
            continue
        tname = parts[2]
        if not tname or tname == "*":
            continue
        tstart = int(parts[3])  # 1-based
        mapq = int(parts[4])
        cigar = parts[5]
        seq = parts[9]
        tags: dict[str, str] = {}
        for field in parts[11:]:
            try:
                tag, _typ, val = field.split(":", 2)
            except ValueError:
                continue
            tags[tag] = val
        q_aln, r_span, alnlen = _cigar_alignment_spans(cigar)
        if alnlen <= 0:
            continue
        nm = int(tags["NM"]) if "NM" in tags and str(tags["NM"]).lstrip("-").isdigit() else 0
        nmatch = max(0, alnlen - nm)
        pid = (nmatch / alnlen) if alnlen > 0 else 0.0
        # Prefer CIGAR-implied full query length so soft/hard clips count in qcov.
        qlen = _cigar_query_len(cigar)
        if qlen <= 0:
            qlen = len(seq) if seq and seq != "*" else q_aln
        qcov = (q_aln / qlen) if qlen > 0 else 0.0
        strand = "-" if (flag & 0x10) else "+"
        tend = tstart + max(r_span - 1, 0)
        rows.append(
            {
                "qname": qname,
                "target": tname,
                "target_start": tstart,
                "target_end": tend,
                "target_len": int(tlen_map.get(tname, 0)),
                "target_strand": strand,
                "alnlen": int(alnlen),
                "mapq": mapq,
                "pid": float(pid),
                "qcov": float(qcov),
                "mei_score": _mei_score_from_metrics(pid=pid, qcov=qcov, mapq=mapq),
                "family": _family_from_target(tname),
            }
        )
    return _pick_best_mei_hits(rows)


# Terminal polyA on Alu/SVA/L1 consensus is not MEI body — remap against a
# trimmed copy so polyA soft-clips cannot "MEI-map" onto the consensus tail.
_MEI_CONSENSUS_POLYA_MIN_TRIM = 8


def trim_mei_consensus_terminal_polya(seq: str, *, min_run: int = _MEI_CONSENSUS_POLYA_MIN_TRIM) -> str:
    """Strip terminal polyA (and leading polyT) from a consensus sequence."""
    s = (seq or "").upper().replace("U", "T")
    if not s:
        return ""
    end = len(s)
    while end > 0 and s[end - 1] == "A":
        end -= 1
    if (len(s) - end) >= int(min_run):
        s = s[:end]
    start = 0
    while start < len(s) and s[start] == "T":
        start += 1
    if start >= int(min_run):
        s = s[start:]
    return s


def _write_polya_trimmed_fasta(src: Path, dst: Path, *, min_run: int = _MEI_CONSENSUS_POLYA_MIN_TRIM) -> int:
    """Write ``dst`` with terminal polyA/T stripped from each record. Returns n trimmed."""
    src = Path(src)
    dst = Path(dst)
    n_trimmed = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open("r", encoding="utf-8") as hin, dst.open("w", encoding="utf-8") as hout:
        name: str | None = None
        seq_parts: list[str] = []

        def _flush() -> None:
            nonlocal n_trimmed, name, seq_parts
            if name is None:
                return
            raw = "".join(seq_parts)
            trimmed = trim_mei_consensus_terminal_polya(raw, min_run=min_run)
            if len(trimmed) != len(raw.upper().replace("U", "T")):
                n_trimmed += 1
            if not trimmed:
                name = None
                seq_parts = []
                return
            # ``name`` is the full FASTA header line including the leading '>'.
            hout.write(name if name.endswith("\n") else f"{name}\n")
            for i in range(0, len(trimmed), 60):
                hout.write(trimmed[i : i + 60] + "\n")
            name = None
            seq_parts = []

        for line in hin:
            if line.startswith(">"):
                _flush()
                name = line.rstrip("\n")
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        _flush()
    return n_trimmed


def _ensure_polya_trimmed_mei_fasta(mei_fasta: Path) -> Path:
    """Return a polyA-trimmed MEI FASTA (sidecar ``*.nopolya.fa``), refreshing when stale."""
    src = Path(mei_fasta)
    if not src.exists():
        raise FileNotFoundError(f"MEI FASTA not found: {src}")
    dst = src.with_name(f"{src.stem}.nopolya{src.suffix}")
    src_mtime = src.stat().st_mtime
    needs = (not dst.exists()) or (dst.stat().st_mtime < src_mtime) or (dst.stat().st_size <= 0)
    if needs:
        _write_polya_trimmed_fasta(src, dst)
        # Drop stale bwa index so _ensure_bwa_index rebuilds.
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            idx = Path(f"{dst}{suffix}")
            if idx.exists():
                idx.unlink()
    return dst


def _ensure_bwa_index(mei_fasta: Path) -> None:
    """Build classic bwa index next to ``mei_fasta`` when missing."""
    if Path(f"{mei_fasta}.bwt").exists() and Path(f"{mei_fasta}.sa").exists():
        return
    subprocess.run(
        ["bwa", "index", str(mei_fasta)],
        check=True,
        capture_output=True,
        text=True,
    )


def _align_queries_to_mei_bwa(
    query_fa: Path,
    mei_fasta: Path,
    *,
    bwa_threads: int = 1,
) -> pd.DataFrame:
    """Remap query FASTA to MEI consensus with ``bwa mem -k10 -T10``."""
    remap_fa = _ensure_polya_trimmed_mei_fasta(mei_fasta)
    _ensure_bwa_index(remap_fa)
    cmd = ["bwa", *_bwa_mem_mei_args(bwa_threads=bwa_threads), str(remap_fa), str(query_fa)]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return _best_hits_from_sam(proc.stdout or "")


def _canonicalize_alignment_hit_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize MEI-alignment hit columns after merges with possible name collisions."""
    out = df.copy()
    expected_cols = [
        "target",
        "family",
        "target_strand",
        "target_start",
        "target_end",
        "target_len",
        "alnlen",
        "mapq",
        "pid",
        "qcov",
        "mei_score",
    ]
    for col in expected_cols:
        preferred_sources = [f"{col}_y", f"{col}_mei", col]
        src = next((name for name in preferred_sources if name in out.columns), None)
        if src is None:
            out[col] = pd.NA
            continue
        if src != col:
            out[col] = out[src]
    return out


def _align_clips_with_minimap2(
    split_df: pd.DataFrame,
    mei_fasta: Path,
    sample: str,
    *,
    bwa_threads: int = 1,
) -> tuple[pd.DataFrame, ClipAlignmentSummary]:
    if "clip_seq" not in split_df.columns:
        raise ValueError(
            "split evidence table is missing 'clip_seq'. Re-run extract-split-evidence with the current code."
        )

    # Keep every input clip (including pure polyA) so polyA_MAPPED counting still
    # sees them; only non-empty trimmed queries are remapped to MEI.
    out = split_df.copy()
    out = out.loc[out["clip_seq"].fillna("").astype(str).str.len() > 0].copy()
    out["clip_id"] = [f"{sample}_{i}" for i in range(len(out))]
    out["clip_seq_coord"] = [
        _trim_poly_at_from_clip(str(seq), str(side))
        for seq, side in zip(out["clip_seq"].fillna(""), out["clip_side"].fillna(""))
    ]
    out["mei_query_seq"] = out["clip_seq_coord"].fillna("").astype(str)
    for col, default in (
        ("mei_hit", False),
        ("family", ""),
        ("target", ""),
        ("target_strand", ""),
        ("target_start", 0),
        ("target_end", 0),
        ("target_len", 0),
        ("alnlen", 0),
        ("mapq", 0),
        ("pid", 0.0),
        ("qcov", 0.0),
        ("mei_score", 0.0),
    ):
        if col not in out.columns:
            out[col] = default

    remap = out.loc[out["mei_query_seq"].str.len() >= _MEI_ALIGN_MIN_TRIMMED_BP].copy()
    if not remap.empty:
        with tempfile.TemporaryDirectory(prefix=f"rtm_mei_{sample}_") as tmpdir:
            tmp = Path(tmpdir)
            query_fa = tmp / "clips.fa"

            with query_fa.open("w", encoding="utf-8") as handle:
                for row in remap.itertuples(index=False):
                    handle.write(f">{row.clip_id}\n{row.mei_query_seq}\n")

            # bwa mem -k10 -T10: recovers ~20 bp Alu/L1/SVA tips that minimap2 -x sr misses.
            best_hits = _align_queries_to_mei_bwa(query_fa, mei_fasta, bwa_threads=bwa_threads)
            hit_rows = remap.merge(best_hits, left_on="clip_id", right_on="qname", how="left")
            hit_rows = _canonicalize_alignment_hit_columns(hit_rows)
            hit_rows = _apply_mei_align_quality_gate(
                hit_rows,
                query_len=hit_rows["mei_query_seq"].fillna("").astype(str).str.len(),
            )
            hit_rows["mei_hit"] = hit_rows["target"].notna() & hit_rows["target"].fillna("").astype(str).ne("")
            hit_cols = [
                "clip_id",
                "mei_hit",
                "family",
                "target",
                "target_strand",
                "target_start",
                "target_end",
                "target_len",
                "alnlen",
                "mapq",
                "pid",
                "qcov",
                "mei_score",
            ]
            hit_cols = [c for c in hit_cols if c in hit_rows.columns]
            hit_slim = hit_rows.loc[:, hit_cols].copy()
            out = out.drop(columns=[c for c in hit_cols if c != "clip_id" and c in out.columns], errors="ignore")
            out = out.merge(hit_slim, on="clip_id", how="left")

    out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool)
    out["family"] = out["family"].fillna("")
    out["target"] = out["target"].fillna("")
    out["target_strand"] = out["target_strand"].fillna("")
    out["target_start"] = out["target_start"].fillna(0).astype(int)
    out["target_end"] = out["target_end"].fillna(0).astype(int)
    out["target_len"] = out["target_len"].fillna(0).astype(int)
    out["alnlen"] = out["alnlen"].fillna(0).astype(int)
    out["mapq"] = out["mapq"].fillna(0).astype(int)
    out["pid"] = out["pid"].fillna(0.0).astype(float)
    out["qcov"] = out["qcov"].fillna(0.0).astype(float)
    out["mei_score"] = out["mei_score"].fillna(0.0).astype(float)

    # Coord columns: same trimmed+gated alignment (legacy consumers).
    out["target_coord"] = out["target"]
    out["family_coord"] = out["family"]
    out["target_strand_coord"] = out["target_strand"]
    out["target_start_coord"] = out["target_start"]
    out["target_end_coord"] = out["target_end"]
    out["target_len_coord"] = out["target_len"]
    out["alnlen_coord"] = out["alnlen"]
    out["mapq_coord"] = out["mapq"]
    out["pid_coord"] = out["pid"]
    out["qcov_coord"] = out["qcov"]
    out["mei_score_coord"] = out["mei_score"]
    out["mei_hit_coord"] = out["mei_hit"]

    summary = ClipAlignmentSummary(
        sample=sample,
        clip_count=len(out),
        paf_hits=int(out["mei_hit"].sum()),
    )
    return out, summary


def _align_discordant_reads_with_minimap2(
    discordant_df: pd.DataFrame,
    mei_fasta: Path,
    sample: str,
    *,
    bwa_threads: int = 1,
) -> tuple[pd.DataFrame, ClipAlignmentSummary]:
    """Remap DPE *soft-clipped* bases (not ref-matched bodies) to MEI consensus."""
    if discordant_df is None or discordant_df.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return pd.DataFrame() if discordant_df is None else discordant_df.copy(), summary

    reads = discordant_df.copy()
    reason = (
        reads["discordant_reasons"].fillna("").astype(str)
        if "discordant_reasons" in reads.columns
        else pd.Series("", index=reads.index)
    )
    # poly_tail_anchor_rescue is a strong discordant class; mates/anchors from those
    # pairs can still carry MEI sequence and must be eligible for consensus remap.
    reason_ok = (
        reason.str.contains("interchrom", regex=False)
        | reason.str.contains("mate_unmapped", regex=False)
        | reason.str.contains("large_insert", regex=False)
        | reason.str.contains("poly_tail_anchor_rescue", regex=False)
    )
    reads = reads.loc[reason_ok].copy()
    if reads.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return discordant_df.copy(), summary

    # Ensure soft_clip_seq exists (derive from side/len + read_seq when needed).
    if "soft_clip_seq" not in reads.columns:
        reads["soft_clip_seq"] = ""
    missing_clip = reads["soft_clip_seq"].fillna("").astype(str).str.len().lt(_DPE_MEI_REMAP_MIN_CLIP_BP)
    if missing_clip.any():
        derived = [
            _discordant_anchor_mei_query_seq(row)
            for row in reads.loc[missing_clip].itertuples(index=False)
        ]
        reads.loc[missing_clip, "soft_clip_seq"] = derived

    raw_query = reads["soft_clip_seq"].fillna("").astype(str)
    clip_side = (
        reads["soft_clip_side"].fillna("").astype(str)
        if "soft_clip_side" in reads.columns
        else pd.Series("", index=reads.index)
    )
    reads["mei_query_seq"] = [
        _trim_poly_at_from_clip(seq, side) for seq, side in zip(raw_query, clip_side)
    ]
    reads = reads.loc[
        reads["mei_query_seq"].str.len() >= max(_DPE_MEI_REMAP_MIN_CLIP_BP, _MEI_ALIGN_MIN_TRIMMED_BP)
    ].copy()
    reads["discordant_id"] = [f"{sample}_disc_{i}" for i in range(len(reads))]
    if reads.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        # Preserve original columns; no MEI hits from full-body remap.
        out = discordant_df.copy()
        out["mei_hit"] = False
        out["family"] = ""
        out["target"] = ""
        out["target_strand"] = ""
        out["target_start"] = 0
        out["target_end"] = 0
        out["target_len"] = 0
        out["alnlen"] = 0
        out["mei_score"] = 0.0
        return out, summary

    with tempfile.TemporaryDirectory(prefix=f"rtm_mei_disc_{sample}_") as tmpdir:
        tmp = Path(tmpdir)
        query_fa = tmp / "discordant_clips.fa"

        with query_fa.open("w", encoding="utf-8") as handle:
            for row in reads.itertuples(index=False):
                handle.write(f">{row.discordant_id}\n{row.mei_query_seq}\n")

        best_hits = _align_queries_to_mei_bwa(query_fa, mei_fasta, bwa_threads=bwa_threads)
        hit_rows = reads.merge(best_hits, left_on="discordant_id", right_on="qname", how="left")
        hit_rows = _canonicalize_alignment_hit_columns(hit_rows)
        hit_rows = _apply_mei_align_quality_gate(
            hit_rows,
            query_len=hit_rows["mei_query_seq"].fillna("").astype(str).str.len(),
        )

    # Left-join hits back onto the full discordant table so unclipped anchors stay.
    out = discordant_df.copy()
    merge_keys = [c for c in ["read_name", "chrom", "pos", "mate_chrom", "mate_pos"] if c in out.columns and c in hit_rows.columns]
    hit_cols = [
        "mei_hit",
        "target",
        "family",
        "target_strand",
        "target_start",
        "target_end",
        "target_len",
        "alnlen",
        "mapq",
        "pid",
        "qcov",
        "mei_score",
    ]
    if not merge_keys:
        # Fallback: index-unsafe; attach by read_name only.
        merge_keys = [c for c in ["read_name"] if c in out.columns and c in hit_rows.columns]
    for col in hit_cols:
        if col not in hit_rows.columns:
            if col == "mei_hit":
                hit_rows[col] = False
            elif col in {"target", "family", "target_strand"}:
                hit_rows[col] = ""
            elif col in {"mei_score", "pid", "qcov"}:
                hit_rows[col] = 0.0
            else:
                hit_rows[col] = 0
    hit_rows["mei_hit"] = hit_rows["target"].notna() & hit_rows["target"].fillna("").astype(str).ne("")
    keep = merge_keys + hit_cols
    keep = [c for c in keep if c in hit_rows.columns]
    hits_slim = hit_rows.loc[:, keep].drop_duplicates(subset=merge_keys, keep="first")
    out = out.merge(hits_slim, on=merge_keys, how="left", suffixes=("", "_mei"))
    out = _canonicalize_alignment_hit_columns(out)
    out["mei_hit"] = out["target"].notna() & out["target"].fillna("").astype(str).ne("")
    out["family"] = out["family"].fillna("")
    out["target"] = out["target"].fillna("")
    out["target_strand"] = out["target_strand"].fillna("")
    out["target_start"] = out["target_start"].fillna(0).astype(int)
    out["target_end"] = out["target_end"].fillna(0).astype(int)
    out["target_len"] = out["target_len"].fillna(0).astype(int)
    out["alnlen"] = out["alnlen"].fillna(0).astype(int)
    out["mei_score"] = out["mei_score"].fillna(0.0).astype(float)

    summary = ClipAlignmentSummary(
        sample=sample,
        clip_count=len(reads),
        paf_hits=int(out["mei_hit"].sum()),
    )
    return out, summary


def _merge_fetched_mate_sequences(
    out: pd.DataFrame,
    fetched: dict[str, tuple[str, int, int, str, int, str]],
) -> pd.DataFrame:
    """Apply BAM-fetched mate fields in O(rows) via map — not per-qname masks."""
    if not fetched:
        return out

    names = out["read_name"].astype(str)
    seq_map = {q: t[0] for q, t in fetched.items()}
    start_map = {q: int(t[1]) for q, t in fetched.items()}
    end_map = {q: int(t[2]) for q, t in fetched.items()}
    side_map = {q: t[3] for q, t in fetched.items()}
    len_map = {q: int(t[4]) for q, t in fetched.items()}
    clip_map = {q: t[5] for q, t in fetched.items()}

    fetched_seq = names.map(seq_map)
    hit = fetched_seq.notna()
    if not bool(hit.any()):
        return out

    # Keep existing mate_seq when already present; always refresh soft-clip fields.
    empty_seq = out["mate_seq"].astype(str).str.len().eq(0)
    fill_seq = hit & empty_seq
    out.loc[fill_seq, "mate_seq"] = fetched_seq.loc[fill_seq].astype(str)
    out.loc[hit, "mate_ref_start"] = names.loc[hit].map(start_map).astype(int)
    out.loc[hit, "mate_ref_end"] = names.loc[hit].map(end_map).astype(int)
    out.loc[hit, "mate_soft_clip_side"] = names.loc[hit].map(side_map).astype(str)
    out.loc[hit, "mate_soft_clip_len"] = names.loc[hit].map(len_map).astype(int)
    out.loc[hit, "mate_soft_clip_seq"] = names.loc[hit].map(clip_map).astype(str)
    return out


def _fetch_discordant_mate_sequences(
    discordant_df: pd.DataFrame,
    bam_path: Path | None,
) -> pd.DataFrame:
    """Attach mate_seq (+ soft-clip fields) to discordant rows."""
    if discordant_df.empty:
        return discordant_df.copy()

    out = discordant_df.copy()
    for col, default in (
        ("mate_seq", ""),
        ("mate_ref_start", 0),
        ("mate_ref_end", 0),
        ("mate_soft_clip_side", ""),
        ("mate_soft_clip_len", 0),
        ("mate_soft_clip_seq", ""),
    ):
        if col not in out.columns:
            out[col] = default

    out["mate_seq"] = out["mate_seq"].fillna("").astype(str)
    out["mate_ref_start"] = pd.to_numeric(out["mate_ref_start"], errors="coerce").fillna(0).astype(int)
    out["mate_ref_end"] = pd.to_numeric(out["mate_ref_end"], errors="coerce").fillna(0).astype(int)
    out["mate_soft_clip_side"] = out["mate_soft_clip_side"].fillna("").astype(str)
    out["mate_soft_clip_len"] = pd.to_numeric(out["mate_soft_clip_len"], errors="coerce").fillna(0).astype(int)
    out["mate_soft_clip_seq"] = out["mate_soft_clip_seq"].fillna("").astype(str)

    if bam_path is None or not bam_path.exists():
        return out

    # Refetch when mate_seq missing OR soft-clip fields missing (older extracts).
    need_fetch = out.loc[
        (
            (out["mate_seq"].str.len() == 0)
            | (
                out["mate_soft_clip_seq"].str.len().eq(0)
                & out["mate_soft_clip_len"].eq(0)
                & out["mate_seq"].str.len().gt(0)
            )
        )
        & (out["mate_chrom"].fillna("").astype(str) != "*")
        & (pd.to_numeric(out["mate_pos"], errors="coerce").fillna(0).astype(int) > 0)
    ].copy()
    if need_fetch.empty:
        return out

    fetched: dict[str, tuple[str, int, int, str, int, str]] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in need_fetch.itertuples(index=False):
            qname = str(rec.read_name)
            if qname in fetched:
                continue
            mate_chrom = str(rec.mate_chrom)
            mate_pos = int(rec.mate_pos)
            start0 = max(0, mate_pos - 1)
            end0 = start0 + 500
            for mate_read in bam.fetch(mate_chrom, start0, end0):
                if mate_read.query_name != qname:
                    continue
                if mate_read.is_secondary or mate_read.is_supplementary:
                    continue
                if mate_read.is_unmapped:
                    continue
                seq = mate_read.query_sequence or ""
                if not seq:
                    continue
                clip_side, clip_len, _clip_pos, clip_seq = _longest_soft_clip_from_read(mate_read)
                fetched[qname] = (
                    seq,
                    int(mate_read.reference_start) + 1 if mate_read.reference_start is not None else 0,
                    int(mate_read.reference_end) if mate_read.reference_end is not None else 0,
                    clip_side,
                    int(clip_len),
                    clip_seq,
                )
                break

    return _merge_fetched_mate_sequences(out, fetched)


def _align_discordant_mates_with_minimap2(
    discordant_df: pd.DataFrame,
    mei_fasta: Path,
    sample: str,
    bam_path: Path | None = None,
    *,
    bwa_threads: int = 1,
) -> tuple[pd.DataFrame, ClipAlignmentSummary]:
    """Map discordant mates to MEI consensus (clip-only when mate is soft-clipped)."""
    fetch_t0 = time.monotonic()
    enriched = _fetch_discordant_mate_sequences(discordant_df, bam_path)
    click.echo(
        f"[mei-annotate] sample={sample} mate_fetch rows={len(discordant_df)} "
        f"elapsed={time.monotonic() - fetch_t0:.1f}s"
    )
    if enriched.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return enriched, summary

    mates = enriched.copy()
    mates["mei_query_seq"] = [
        _discordant_mate_mei_query_seq(row) for row in mates.itertuples(index=False)
    ]
    mates = mates.loc[mates["mei_query_seq"].fillna("").astype(str).str.len() >= _DPE_MEI_REMAP_MIN_CLIP_BP].copy()
    if mates.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        out = enriched.copy()
        out["mate_mei_hit"] = False
        return out, summary

    mates["mate_query_id"] = [f"{sample}_mate_{i}" for i in range(len(mates))]
    mate_query = mates.copy()
    # Shared discordant aligner historically expected read_seq; pass clip/full query.
    # PolyA trim uses mate soft-clip side when the query is a clip; full-mate queries
    # leave side empty so trim keeps the longer non-poly residual.
    mate_query["read_seq"] = mate_query["mei_query_seq"].fillna("").astype(str)
    mate_query["soft_clip_seq"] = mate_query["mei_query_seq"].fillna("").astype(str)
    mate_has_clip = (
        pd.to_numeric(mate_query.get("mate_soft_clip_len", 0), errors="coerce").fillna(0).astype(int).gt(0)
        | mate_query.get("mate_soft_clip_seq", pd.Series("", index=mate_query.index))
        .fillna("")
        .astype(str)
        .str.len()
        .gt(0)
    )
    mate_side = (
        mate_query["mate_soft_clip_side"].fillna("").astype(str)
        if "mate_soft_clip_side" in mate_query.columns
        else pd.Series("", index=mate_query.index)
    )
    mate_query["soft_clip_side"] = mate_side.where(mate_has_clip, "")
    mate_query["soft_clip_len"] = mate_query["mei_query_seq"].str.len().astype(int)
    # Force a strong reason so the shared discordant aligner does not drop mates
    # whose only original reason was poly_tail_anchor_rescue (or empty).
    mate_query["discordant_reasons"] = "interchrom"
    aln_t0 = time.monotonic()
    hits, summary = _align_discordant_reads_with_minimap2(
        mate_query, mei_fasta, sample=sample, bwa_threads=bwa_threads
    )
    click.echo(
        f"[mei-annotate] sample={sample} mate_bwa_mem queries={len(mate_query)} "
        f"mei_hits={summary.paf_hits} elapsed={time.monotonic() - aln_t0:.1f}s"
    )
    keep_cols = [
        "read_name",
        "mate_query_id",
        "mei_hit",
        "target",
        "family",
        "target_strand",
        "target_start",
        "target_end",
        "target_len",
        "alnlen",
        "mei_score",
    ]
    keep_cols = [c for c in keep_cols if c in hits.columns]
    hits = hits.loc[:, keep_cols].rename(
        columns={
            "mei_hit": "mate_mei_hit",
            "target": "mate_mei_target",
            "family": "mate_mei_family",
            "target_strand": "mate_mei_strand",
            "target_start": "mate_mei_start",
            "target_end": "mate_mei_end",
            "target_len": "mate_mei_target_len",
            "alnlen": "mate_mei_alnlen",
            "mei_score": "mate_mei_score",
        }
    )
    # hits from shared aligner may be the full mate_query table; merge by read_name.
    hit_merge_cols = [c for c in hits.columns if c != "mate_query_id"]
    out = enriched.merge(hits.loc[:, hit_merge_cols], on="read_name", how="left")
    for col in [
        "mate_mei_hit",
        "mate_mei_start",
        "mate_mei_end",
        "mate_mei_target_len",
        "mate_mei_alnlen",
        "mate_mei_score",
    ]:
        if col not in out.columns:
            if col == "mate_mei_hit":
                out[col] = False
            elif col in {"mate_mei_start", "mate_mei_end", "mate_mei_target_len", "mate_mei_alnlen"}:
                out[col] = 0
            else:
                out[col] = 0.0
    out["mate_mei_hit"] = out["mate_mei_hit"].fillna(False).astype(bool)
    out["mate_mei_start"] = pd.to_numeric(out["mate_mei_start"], errors="coerce").fillna(0).astype(int)
    out["mate_mei_end"] = pd.to_numeric(out["mate_mei_end"], errors="coerce").fillna(0).astype(int)
    out["mate_mei_score"] = pd.to_numeric(out["mate_mei_score"], errors="coerce").fillna(0.0).astype(float)
    return out, summary


def _resolve_supporting_reads_detail_path(reuse_dir: Path) -> Path:
    """Locate ``supporting_reads_detail.mei.{parquet,tsv}`` beside a prior annotate outdir."""
    root = Path(reuse_dir)
    for name in ("supporting_reads_detail.mei.parquet", "supporting_reads_detail.mei.tsv"):
        path = root / name
        if path.exists():
            return path
    # Also accept a gold_review sibling layout where detail sits next to the TSV.
    for name in ("supporting_reads_detail.mei.parquet", "supporting_reads_detail.mei.tsv"):
        path = root.parent / name if root.is_file() else root / name
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No supporting_reads_detail.mei.{{parquet,tsv}} under {root}. "
        "Pass --reuse-mei-annotate-dir to a prior annotate output directory."
    )


def _load_supporting_reads_detail_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def _hydrate_sample_mei_hits_from_detail(
    *,
    sample: str,
    split_df: pd.DataFrame,
    discordant_df: pd.DataFrame,
    detail: pd.DataFrame,
) -> dict[str, object]:
    """Rebuild per-read MEI hit frames from a prior supporting_reads_detail table.

    Skips MEI consensus remap and mate-BAM fetch. Non-hit evidence rows are retained so
    residual complex fractions stay correct.
    """
    t0 = time.monotonic()
    key_cols = ["chrom", "window_start", "window_end", "read_name"]
    sample_detail = (
        detail.loc[detail["sample"].fillna("").astype(str) == str(sample)].copy()
        if (not detail.empty and "sample" in detail.columns)
        else detail.iloc[0:0].copy()
    )

    def _ensure_int(series: pd.Series) -> pd.Series:
        return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)

    # --- split ---
    split_hits = split_df.copy() if split_df is not None else pd.DataFrame()
    if not split_hits.empty:
        for col, default in (
            ("mei_hit", False),
            ("mei_hit_coord", False),
            ("target", ""),
            ("family", ""),
            ("target_strand", ""),
            ("target_start", 0),
            ("target_end", 0),
            ("target_len", 0),
            ("mei_score", 0.0),
            ("alnlen", 0),
            ("mapq", 0),
            ("pid", 0.0),
            ("qcov", 0.0),
        ):
            if col not in split_hits.columns:
                split_hits[col] = default
        sr = sample_detail.loc[sample_detail.get("evidence_type", pd.Series(dtype=object)).astype(str) == "SR"].copy()
        if not sr.empty and all(c in split_hits.columns for c in key_cols):
            for c in key_cols:
                sr[c] = sr[c]
                if c in ("window_start", "window_end"):
                    sr[c] = _ensure_int(sr[c])
                    split_hits[c] = _ensure_int(split_hits[c])
                else:
                    sr[c] = sr[c].fillna("").astype(str)
                    split_hits[c] = split_hits[c].fillna("").astype(str)
            hit_cols = {
                "mei_target": "target",
                "mei_start": "target_start",
                "mei_end": "target_end",
                "mei_strand": "target_strand",
            }
            keep = key_cols + [c for c in hit_cols if c in sr.columns]
            m = sr.loc[:, keep].drop_duplicates(key_cols, keep="first").rename(columns=hit_cols)
            m["target"] = m.get("target", pd.Series("", index=m.index)).fillna("").astype(str)
            m["target_start"] = _ensure_int(m.get("target_start", 0))
            m["target_end"] = _ensure_int(m.get("target_end", 0))
            strand = m.get("target_strand", pd.Series("", index=m.index)).fillna("").astype(str)
            m["target_strand"] = strand.where(strand.isin(["+", "-"]), "")
            m["family"] = m["target"].map(_normalize_mei_family_token)
            m["mei_hit"] = m["target"].str.len() > 0
            m["mei_score"] = m["mei_hit"].astype(float)
            drop_overlap = [c for c in m.columns if c not in key_cols and c in split_hits.columns]
            split_hits = split_hits.drop(columns=drop_overlap, errors="ignore").merge(
                m, on=key_cols, how="left"
            )
            split_hits["target"] = split_hits["target"].fillna("").astype(str)
            split_hits["mei_hit"] = split_hits["target"].str.len().gt(0)
            split_hits["mei_score"] = pd.to_numeric(split_hits["mei_score"], errors="coerce").fillna(0.0)
            split_hits["family"] = split_hits["family"].fillna("").astype(str)
            split_hits["target_start"] = _ensure_int(split_hits["target_start"])
            split_hits["target_end"] = _ensure_int(split_hits["target_end"])
            strand_hit = split_hits.get(
                "target_strand", pd.Series("", index=split_hits.index)
            ).fillna("").astype(str)
            split_hits["target_strand"] = strand_hit.where(strand_hit.isin(["+", "-"]), "")

    split_paf = int(split_hits["mei_hit"].fillna(False).astype(bool).sum()) if (
        not split_hits.empty and "mei_hit" in split_hits.columns
    ) else 0
    split_summary = ClipAlignmentSummary(
        sample=sample,
        clip_count=int(len(split_df) if split_df is not None else 0),
        paf_hits=split_paf,
    )

    # --- discordant ---
    disc_hits = discordant_df.copy() if discordant_df is not None else pd.DataFrame()
    if not disc_hits.empty:
        for col, default in (
            ("mei_hit", False),
            ("mate_mei_hit", False),
            ("vntr_rescue", False),
            ("polya_rescue", False),
            ("mei_hit_source", ""),
            ("target", ""),
            ("family", ""),
            ("target_strand", ""),
            ("target_start", 0),
            ("target_end", 0),
            ("target_len", 0),
            ("mate_mei_target", ""),
            ("mate_mei_strand", ""),
            ("mate_mei_start", 0),
            ("mate_mei_end", 0),
            ("mei_score", 0.0),
            ("alnlen", 0),
            ("mapq", 0),
            ("pid", 0.0),
            ("qcov", 0.0),
        ):
            if col not in disc_hits.columns:
                disc_hits[col] = default
        dpe = sample_detail.loc[sample_detail.get("evidence_type", pd.Series(dtype=object)).astype(str) == "DPE"].copy()
        if not dpe.empty and all(c in disc_hits.columns for c in key_cols):
            for c in key_cols:
                if c in ("window_start", "window_end"):
                    dpe[c] = _ensure_int(dpe[c])
                    disc_hits[c] = _ensure_int(disc_hits[c])
                else:
                    dpe[c] = dpe[c].fillna("").astype(str)
                    disc_hits[c] = disc_hits[c].fillna("").astype(str)
            rename = {
                "mei_target": "target",
                "mei_start": "target_start",
                "mei_end": "target_end",
                "mei_strand": "target_strand",
            }
            keep = key_cols + [
                c
                for c in (
                    "mei_target",
                    "mei_start",
                    "mei_end",
                    "mei_strand",
                    "mate_mei_target",
                    "mate_mei_start",
                    "mate_mei_end",
                    "mate_mei_strand",
                    "mei_hit",
                    "mate_mei_hit",
                    "vntr_rescue",
                    "polya_rescue",
                    "mei_hit_source",
                )
                if c in dpe.columns
            ]
            m = dpe.loc[:, keep].drop_duplicates(key_cols, keep="first").rename(columns=rename)
            if "target" in m.columns:
                m["target"] = m["target"].fillna("").astype(str)
                m["family"] = m["target"].map(_normalize_mei_family_token)
            for c in ("target_start", "target_end", "mate_mei_start", "mate_mei_end"):
                if c in m.columns:
                    m[c] = _ensure_int(m[c])
            for c in ("target_strand", "mate_mei_strand"):
                if c in m.columns:
                    strand = m[c].fillna("").astype(str)
                    m[c] = strand.where(strand.isin(["+", "-"]), "")
            for c in ("mei_hit", "mate_mei_hit", "vntr_rescue", "polya_rescue"):
                if c in m.columns:
                    m[c] = m[c].fillna(False).infer_objects(copy=False).astype(bool)
            if "mei_hit_source" in m.columns:
                m["mei_hit_source"] = m["mei_hit_source"].fillna("").astype(str)
            m["mei_score"] = (
                m.get("mei_hit", pd.Series(False, index=m.index)).fillna(False).astype(bool)
                | m.get("mate_mei_hit", pd.Series(False, index=m.index)).fillna(False).astype(bool)
            ).astype(float)
            drop_overlap = [c for c in m.columns if c not in key_cols and c in disc_hits.columns]
            disc_hits = disc_hits.drop(columns=drop_overlap, errors="ignore").merge(
                m, on=key_cols, how="left"
            )
            for c in ("mei_hit", "mate_mei_hit", "vntr_rescue", "polya_rescue"):
                if c in disc_hits.columns:
                    disc_hits[c] = disc_hits[c].fillna(False).infer_objects(copy=False).astype(bool)
            disc_hits["target"] = disc_hits["target"].fillna("").astype(str)
            disc_hits["mate_mei_target"] = disc_hits.get(
                "mate_mei_target", pd.Series("", index=disc_hits.index)
            ).fillna("").astype(str)
            disc_hits["mei_hit_source"] = disc_hits.get(
                "mei_hit_source", pd.Series("", index=disc_hits.index)
            ).fillna("").astype(str)
            disc_hits["family"] = disc_hits.get("family", pd.Series("", index=disc_hits.index)).fillna("").astype(str)
            disc_hits["mei_score"] = pd.to_numeric(disc_hits["mei_score"], errors="coerce").fillna(0.0)
            for c in ("target_strand", "mate_mei_strand"):
                strand = disc_hits.get(c, pd.Series("", index=disc_hits.index)).fillna("").astype(str)
                disc_hits[c] = strand.where(strand.isin(["+", "-"]), "")

    disc_paf = 0
    mate_paf = 0
    if not disc_hits.empty:
        if "mei_hit" in disc_hits.columns:
            disc_paf = int(disc_hits["mei_hit"].fillna(False).astype(bool).sum())
        if "mate_mei_hit" in disc_hits.columns:
            mate_paf = int(disc_hits["mate_mei_hit"].fillna(False).astype(bool).sum())
    disc_summary = ClipAlignmentSummary(
        sample=sample,
        clip_count=int(len(discordant_df) if discordant_df is not None else 0),
        paf_hits=disc_paf,
    )
    disc_mate_summary = ClipAlignmentSummary(
        sample=f"{sample}_mate",
        clip_count=int(len(discordant_df) if discordant_df is not None else 0),
        paf_hits=mate_paf,
    )
    click.echo(
        f"[mei-annotate] sample={sample} reused detail hits "
        f"split_mei={split_paf} disc_anchor_mei={disc_paf} disc_mate_mei={mate_paf} "
        f"elapsed={time.monotonic() - t0:.1f}s"
    )
    return {
        "sample": sample,
        "split_hits": split_hits,
        "split_summary": split_summary,
        "disc_hits": disc_hits,
        "disc_summary": disc_summary,
        "disc_mate_summary": disc_mate_summary,
    }


def _remap_one_sample_mei_evidence(
    *,
    sample: str,
    split_df: pd.DataFrame,
    discordant_df: pd.DataFrame,
    mei_fasta: Path,
    bam_path: Path | None,
    mate_bam_path: Path | None,
    bwa_threads: int = 1,
) -> dict[str, object]:
    """Split + discordant (anchor/mate) MEI remaps for one sample."""
    t0 = time.monotonic()
    click.echo(f"[mei-annotate] sample={sample} remap start bwa_threads={max(1, int(bwa_threads))}")
    split_t0 = time.monotonic()
    split_hits, split_summary = _align_clips_with_minimap2(
        split_df, mei_fasta, sample=sample, bwa_threads=bwa_threads
    )
    click.echo(
        f"[mei-annotate] sample={sample} split_bwa_mem rows={len(split_df)} "
        f"mei_hits={split_summary.paf_hits} "
        f"elapsed={time.monotonic() - split_t0:.1f}s"
    )
    disc_t0 = time.monotonic()
    disc_anchor_hits, disc_summary = _align_discordant_reads_with_minimap2(
        discordant_df, mei_fasta, sample=sample, bwa_threads=bwa_threads
    )
    click.echo(
        f"[mei-annotate] sample={sample} disc_anchor_bwa_mem rows={len(discordant_df)} "
        f"mei_hits={disc_summary.paf_hits} elapsed={time.monotonic() - disc_t0:.1f}s"
    )
    disc_mate_hits, disc_mate_summary = _align_discordant_mates_with_minimap2(
        discordant_df,
        mei_fasta,
        sample=f"{sample}_mate",
        bam_path=mate_bam_path or bam_path,
        bwa_threads=bwa_threads,
    )
    post_t0 = time.monotonic()
    disc_hits = _attach_mei_hits_to_discordant_rows(discordant_df, disc_anchor_hits, disc_mate_hits)
    disc_hits = _rescue_vntr_like_discordant_mei_hits(disc_hits)
    disc_hits = _rescue_polya_like_discordant_mei_hits(disc_hits)
    split_hits = _enrich_split_hits_with_mate_positions(split_hits, bam_path)
    # PolyA/T soft-clips → polyA_MAPPED only (never MEI_MAPPED/SR/side coords).
    split_hits = _demote_polya_split_mei_hits(split_hits)
    # CCCTCT / SVA-VNTR soft-clips → VNTR_MAPPED (never MEI-SR), before short rescue.
    split_hits = _annotate_vntr_like_split_clips(split_hits, discordant_df=disc_hits)
    # Short-clip rescue after DPE remap so DPE can seed when SR≥20 is absent.
    split_hits = _annotate_short_mei_seed_rescue(split_hits, discordant_df=disc_hits)
    n_short_rescue = (
        int(split_hits["short_mei_seed_rescued"].fillna(False).astype(bool).sum())
        if "short_mei_seed_rescued" in split_hits.columns
        else 0
    )
    n_vntr_split = (
        int(split_hits["vntr_rescue"].fillna(False).astype(bool).sum())
        if "vntr_rescue" in split_hits.columns
        else 0
    )
    click.echo(
        f"[mei-annotate] sample={sample} attach/rescue/enrich "
        f"short_mei_rescued={n_short_rescue} vntr_split={n_vntr_split} "
        f"elapsed={time.monotonic() - post_t0:.1f}s"
    )
    click.echo(
        f"[mei-annotate] sample={sample} remap done elapsed={time.monotonic() - t0:.1f}s "
        f"split_mei={getattr(split_summary, 'paf_hits', 0)} "
        f"disc_anchor_mei={getattr(disc_summary, 'paf_hits', 0)} "
        f"disc_mate_mei={getattr(disc_mate_summary, 'paf_hits', 0)}"
    )
    return {
        "sample": sample,
        "split_hits": split_hits,
        "split_summary": split_summary,
        "disc_hits": disc_hits,
        "disc_summary": disc_summary,
        "disc_mate_summary": disc_mate_summary,
    }


def _enrich_discordant_anchor_hits_with_mate_mei(
    anchor_hits: pd.DataFrame,
    mate_hits: pd.DataFrame,
) -> pd.DataFrame:
    """Promote mate MEI mappings into discordant anchor hit rows for metrics/MEI_MAPPED."""
    if anchor_hits.empty:
        return anchor_hits.copy()
    out = anchor_hits.copy()
    if mate_hits.empty or "mate_mei_hit" not in mate_hits.columns:
        return out

    merge_keys = [
        c for c in ["read_name", "chrom", "window_start", "window_end"] if c in out.columns and c in mate_hits.columns
    ]
    mate_cols = [
        c
        for c in [
            "mate_mei_hit",
            "mate_mei_start",
            "mate_mei_end",
            "mate_mei_target",
            "mate_mei_family",
            "mate_mei_strand",
            "mate_mei_score",
        ]
        if c in mate_hits.columns
    ]
    if not merge_keys or not mate_cols:
        return out

    merged = out.merge(
        mate_hits.loc[:, merge_keys + mate_cols].drop_duplicates(merge_keys),
        on=merge_keys,
        how="left",
    )
    anchor_hit = merged["mei_hit"].fillna(False).astype(bool) if "mei_hit" in merged.columns else pd.Series(False, index=merged.index)
    mate_hit = merged["mate_mei_hit"].fillna(False).astype(bool)
    use_mate = (~anchor_hit) & mate_hit
    if use_mate.any():
        if "target_start" in merged.columns:
            merged.loc[use_mate, "target_start"] = merged.loc[use_mate, "mate_mei_start"]
        if "target_end" in merged.columns:
            merged.loc[use_mate, "target_end"] = merged.loc[use_mate, "mate_mei_end"]
        if "target" in merged.columns and "mate_mei_target" in merged.columns:
            merged.loc[use_mate, "target"] = merged.loc[use_mate, "mate_mei_target"]
        if "family" in merged.columns and "mate_mei_family" in merged.columns:
            merged.loc[use_mate, "family"] = merged.loc[use_mate, "mate_mei_family"]
        if "target_strand" in merged.columns and "mate_mei_strand" in merged.columns:
            merged.loc[use_mate, "target_strand"] = merged.loc[use_mate, "mate_mei_strand"]
        if "mei_score" in merged.columns and "mate_mei_score" in merged.columns:
            merged.loc[use_mate, "mei_score"] = merged.loc[use_mate, "mate_mei_score"]
        if "mei_hit" in merged.columns:
            merged.loc[use_mate, "mei_hit"] = True
    if "mei_hit" in merged.columns:
        merged["mei_hit"] = merged["mei_hit"].fillna(False).astype(bool) | mate_hit
    return merged


def _discordant_row_mei_mapped(df: pd.DataFrame) -> pd.Series:
    """True when anchor and/or mate remapped to MEI consensus (not polyA/VNTR rescue)."""
    if df.empty:
        return pd.Series(dtype=bool)
    mapped = (
        df["mei_hit"].fillna(False).astype(bool)
        if "mei_hit" in df.columns
        else pd.Series(False, index=df.index)
    )
    if "mate_mei_hit" in df.columns:
        mapped = mapped | df["mate_mei_hit"].fillna(False).astype(bool)
    # Second-pass rescues are tracked separately (polyA_MAPPED / VNTR_MAPPED).
    if "polya_rescue" in df.columns:
        mapped = mapped & ~df["polya_rescue"].fillna(False).astype(bool)
    if "vntr_rescue" in df.columns:
        mapped = mapped & ~df["vntr_rescue"].fillna(False).astype(bool)
    if "mei_hit_source" in df.columns:
        src = df["mei_hit_source"].fillna("").astype(str)
        mapped = mapped & ~src.isin({"polya_rescue", "vntr_rescue"})
    return mapped


def _discordant_row_rescue_mapped(df: pd.DataFrame) -> pd.Series:
    """True for polyA/VNTR second-pass rescues (excluded from residual complex)."""
    if df.empty:
        return pd.Series(dtype=bool)
    flagged = pd.Series(False, index=df.index)
    if "polya_rescue" in df.columns:
        flagged = flagged | df["polya_rescue"].fillna(False).astype(bool)
    if "vntr_rescue" in df.columns:
        flagged = flagged | df["vntr_rescue"].fillna(False).astype(bool)
    if "mei_hit_source" in df.columns:
        src = df["mei_hit_source"].fillna("").astype(str)
        flagged = flagged | src.isin({"polya_rescue", "vntr_rescue"})
    return flagged


def _discordant_row_mei_related(df: pd.DataFrame) -> pd.Series:
    """Consensus MEI hit or polyA/VNTR rescue (for residual-complex exclusion)."""
    if df.empty:
        return pd.Series(dtype=bool)
    return _discordant_row_mei_mapped(df) | _discordant_row_rescue_mapped(df)


def _attach_mei_hits_to_discordant_rows(
    discordant_df: pd.DataFrame,
    anchor_hits: pd.DataFrame,
    mate_hits: pd.DataFrame,
) -> pd.DataFrame:
    """
    Left-join MEI hit status onto all locus-assigned discordant rows.

    Unlike enrich-on-anchor-subset, this keeps poly-tail-only (and other) pairs so
    mate MEI hits contribute to MEI_MAPPED counts and residual complex fractions.
    """
    if discordant_df.empty:
        return discordant_df.copy()

    out = discordant_df.copy()
    merge_keys = [
        c
        for c in ["read_name", "chrom", "window_start", "window_end"]
        if c in out.columns
    ]
    if not merge_keys:
        out["mei_hit"] = False
        out["mate_mei_hit"] = False
        return out

    out["mei_hit"] = False
    out["mate_mei_hit"] = False
    for col, default in (
        ("target", ""),
        ("family", ""),
        ("target_strand", ""),
        ("target_start", 0),
        ("target_end", 0),
        ("mei_score", 0.0),
    ):
        if col not in out.columns:
            out[col] = default

    if not anchor_hits.empty:
        a_keys = [c for c in merge_keys if c in anchor_hits.columns]
        a_cols = [
            c
            for c in [
                "mei_hit",
                "target",
                "family",
                "target_strand",
                "target_start",
                "target_end",
                "mei_score",
            ]
            if c in anchor_hits.columns
        ]
        if a_keys and a_cols:
            a = anchor_hits.loc[:, a_keys + a_cols].drop_duplicates(a_keys)
            out = out.drop(columns=[c for c in a_cols if c in out.columns], errors="ignore")
            out = out.merge(a, on=a_keys, how="left")
            out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool)

    if not mate_hits.empty and "mate_mei_hit" in mate_hits.columns:
        m_keys = [c for c in merge_keys if c in mate_hits.columns]
        m_cols = [
            c
            for c in [
                "mate_mei_hit",
                "mate_mei_start",
                "mate_mei_end",
                "mate_mei_target",
                "mate_mei_family",
                "mate_mei_strand",
                "mate_mei_score",
            ]
            if c in mate_hits.columns
        ]
        if m_keys and m_cols:
            m = mate_hits.loc[:, m_keys + m_cols].drop_duplicates(m_keys)
            out = out.drop(columns=[c for c in m_cols if c in out.columns], errors="ignore")
            out = out.merge(m, on=m_keys, how="left")
            out["mate_mei_hit"] = out["mate_mei_hit"].fillna(False).astype(bool)
            use_mate = (~out["mei_hit"].fillna(False).astype(bool)) & out["mate_mei_hit"]
            if use_mate.any():
                if "mate_mei_start" in out.columns:
                    out.loc[use_mate, "target_start"] = out.loc[use_mate, "mate_mei_start"]
                if "mate_mei_end" in out.columns:
                    out.loc[use_mate, "target_end"] = out.loc[use_mate, "mate_mei_end"]
                if "mate_mei_target" in out.columns:
                    out.loc[use_mate, "target"] = out.loc[use_mate, "mate_mei_target"]
                if "mate_mei_family" in out.columns:
                    out.loc[use_mate, "family"] = out.loc[use_mate, "mate_mei_family"]
                if "mate_mei_strand" in out.columns:
                    out.loc[use_mate, "target_strand"] = out.loc[use_mate, "mate_mei_strand"]
                if "mate_mei_score" in out.columns:
                    out.loc[use_mate, "mei_score"] = out.loc[use_mate, "mate_mei_score"]
                out.loc[use_mate, "mei_hit"] = True
            out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool) | out["mate_mei_hit"].fillna(False).astype(bool)
    else:
        out["mate_mei_hit"] = out.get("mate_mei_hit", pd.Series(False, index=out.index))
        out["mate_mei_hit"] = out["mate_mei_hit"].fillna(False).astype(bool)

    out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool)
    out["target"] = out["target"].fillna("")
    out["family"] = out["family"].fillna("")
    out["target_strand"] = out["target_strand"].fillna("")
    out["target_start"] = pd.to_numeric(out.get("target_start", 0), errors="coerce").fillna(0).astype(int)
    out["target_end"] = pd.to_numeric(out.get("target_end", 0), errors="coerce").fillna(0).astype(int)
    out["mei_score"] = pd.to_numeric(out.get("mei_score", 0.0), errors="coerce").fillna(0.0).astype(float)
    if "mei_hit_source" not in out.columns:
        out["mei_hit_source"] = ""
    out.loc[out["mei_hit"] & (out["mei_hit_source"].fillna("").astype(str).str.len() == 0), "mei_hit_source"] = "consensus"
    return out


_VNTR_HEXAMERS = ("CCCTCT", "CTCTCC", "CTCCCT", "TCCCTC", "CCCTCTC", "GGGAGA")
_VNTR_MIN_SEQ_LEN = 40
# Soft-clips can be shorter than discordant mates; still require clear hexamer content.
_VNTR_MIN_CLIP_SEQ_LEN = 20
_VNTR_MIN_HEXAMER_HITS = 3
_VNTR_MIN_HEXAMER_COV = 0.30
_VNTR_RESCUE_MIN_SUPPORT_READS = 3
_VNTR_RESCUE_MIN_FAMILY_PURITY = 0.60
_VNTR_RESCUE_FAMILIES = frozenset({"SVA"})


def _sva_vntr_like_score(seq: str, *, min_seq_len: int | None = None) -> float:
    """
    Score SVA VNTR / hexamer-repeat content in [0, 1].

    SVA central VNTR is (CCCTCT)n-like. Short-read remap to Dfam consensus often
    fails on allele-specific VNTR even when BAM placed the mate on a genomic SVA.
    """
    s = (seq or "").upper().replace("U", "T")
    min_len = int(_VNTR_MIN_SEQ_LEN if min_seq_len is None else min_seq_len)
    if len(s) < min_len:
        return 0.0
    hits = 0
    covered = [False] * len(s)
    for motif in _VNTR_HEXAMERS:
        start = 0
        mlen = len(motif)
        while True:
            i = s.find(motif, start)
            if i < 0:
                break
            hits += 1
            for j in range(i, min(len(s), i + mlen)):
                covered[j] = True
            start = i + 1
    cov = sum(1 for x in covered if x) / float(len(s))
    # Require both enough motif hits and substantial coverage.
    if hits < _VNTR_MIN_HEXAMER_HITS and cov < _VNTR_MIN_HEXAMER_COV:
        return 0.0
    hit_term = min(1.0, hits / 6.0)
    return float(max(0.0, min(1.0, 0.55 * cov + 0.45 * hit_term)))


def _is_sva_vntr_like(
    seq: str, *, min_score: float = 0.35, min_seq_len: int | None = None
) -> bool:
    return _sva_vntr_like_score(seq, min_seq_len=min_seq_len) >= float(min_score)


def _rescue_vntr_like_discordant_mei_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Second-pass VNTR_MAPPED rescue for VNTR-like discordant sequences.

    If a locus already has consistent SVA consensus support, flag non-hitting
    mates/anchors whose sequence looks like SVA VNTR. These count toward
    VNTR_MAPPED (not MEI_MAPPED) and are excluded from residual complex.
    """
    if df.empty:
        return df.copy()
    out = df.copy()
    if "mei_hit" not in out.columns:
        out["mei_hit"] = False
    out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool)
    if "mei_hit_source" not in out.columns:
        out["mei_hit_source"] = ""
    out["mei_hit_source"] = out["mei_hit_source"].fillna("").astype(str)
    if "family" not in out.columns:
        out["family"] = ""
    if "target" not in out.columns:
        out["target"] = ""
    if "mei_score" not in out.columns:
        out["mei_score"] = 0.0
    if "vntr_like_score" not in out.columns:
        out["vntr_like_score"] = 0.0
    if "vntr_rescue" not in out.columns:
        out["vntr_rescue"] = False

    # Prefer mate_seq (insert side); fall back to anchor read_seq.
    mate_seq = out["mate_seq"].fillna("").astype(str) if "mate_seq" in out.columns else pd.Series("", index=out.index)
    read_seq = out["read_seq"].fillna("").astype(str) if "read_seq" in out.columns else pd.Series("", index=out.index)
    prefer_mate = mate_seq.str.len() >= _VNTR_MIN_SEQ_LEN
    seqs = mate_seq.where(prefer_mate, read_seq)
    out["vntr_like_score"] = seqs.map(_sva_vntr_like_score).astype(float)
    # Consensus-only gate: do not treat prior rescues as mapped support.
    consensus_mapped = _discordant_row_mei_mapped(out)
    vntr_cand = (~consensus_mapped) & (~out["vntr_rescue"].fillna(False).astype(bool)) & (out["vntr_like_score"] >= 0.35)
    if not vntr_cand.any():
        return out

    # Locus-level SVA support from consensus hits already present.
    mapped = consensus_mapped
    fam = out["family"].fillna("").astype(str).map(_family_from_target)
    # Also derive family from target name when family column empty.
    empty_fam = fam.eq("") | fam.eq("OTHER")
    if empty_fam.any() and "target" in out.columns:
        fam = fam.where(~empty_fam, out["target"].fillna("").astype(str).map(_family_from_target))
    out["_fam_norm"] = fam
    support = (
        out.loc[mapped]
        .groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            support_reads=("read_name", "nunique"),
            top_family=(
                "_fam_norm",
                lambda s: s.value_counts().index[0] if len(s) else "",
            ),
        )
    )
    if support.empty:
        out = out.drop(columns=["_fam_norm"], errors="ignore")
        return out

    # Family purity among mapped reads at locus.
    mapped_rows = out.loc[mapped, ["chrom", "window_start", "window_end", "_fam_norm", "read_name"]].copy()
    purity_parts = []
    for (chrom, ws, we), grp in mapped_rows.groupby(["chrom", "window_start", "window_end"]):
        vc = grp["_fam_norm"].value_counts()
        top = str(vc.index[0]) if len(vc) else ""
        pur = float(vc.iloc[0] / max(len(grp), 1)) if len(vc) else 0.0
        purity_parts.append(
            {"chrom": chrom, "window_start": ws, "window_end": we, "top_family": top, "family_purity": pur}
        )
    purity = pd.DataFrame(purity_parts)
    support = support.drop(columns=["top_family"], errors="ignore").merge(
        purity, on=["chrom", "window_start", "window_end"], how="left"
    )
    eligible = support.loc[
        (support["support_reads"] >= _VNTR_RESCUE_MIN_SUPPORT_READS)
        & (support["family_purity"].fillna(0.0) >= _VNTR_RESCUE_MIN_FAMILY_PURITY)
        & (support["top_family"].fillna("").astype(str).isin(_VNTR_RESCUE_FAMILIES))
    ].copy()
    if eligible.empty:
        out = out.drop(columns=["_fam_norm"], errors="ignore")
        return out

    out = out.merge(
        eligible[["chrom", "window_start", "window_end", "top_family", "support_reads"]],
        on=["chrom", "window_start", "window_end"],
        how="left",
        suffixes=("", "_elig"),
    )
    # Recompute mask after merge so boolean alignment follows the merged index.
    vntr_cand = (
        (~_discordant_row_mei_mapped(out))
        & (~out["vntr_rescue"].fillna(False).astype(bool))
        & (out["vntr_like_score"].fillna(0.0) >= 0.35)
    )
    rescue = vntr_cand & out["top_family"].fillna("").astype(str).isin(_VNTR_RESCUE_FAMILIES)
    n_rescue = int(rescue.sum())
    if n_rescue:
        prefer_mate = (
            out["mate_seq"].fillna("").astype(str).str.len() >= _VNTR_MIN_SEQ_LEN
            if "mate_seq" in out.columns
            else pd.Series(False, index=out.index)
        )
        # Snapshot pre-rescue mapped rows for coordinate imputation.
        mapped_before = _discordant_row_mei_mapped(out)
        # Do NOT set mei_hit — VNTR rescues are counted as VNTR_MAPPED separately.
        out.loc[rescue, "vntr_rescue"] = True
        out.loc[rescue, "mei_hit_source"] = "vntr_rescue"
        out.loc[rescue, "family"] = out.loc[rescue, "top_family"]
        # Keep a stable synthetic target label for downstream family parsing.
        out.loc[rescue, "target"] = out.loc[rescue, "top_family"].map(
            lambda f: f"{f}_VNTR_rescue#Retroposon/{f}" if f == "SVA" else f"{f}_VNTR_rescue"
        )
        # Soft score from VNTR strength; not a real alignment score.
        out.loc[rescue, "mei_score"] = out.loc[rescue, "vntr_like_score"].clip(lower=0.35, upper=0.85)
        if "mate_mei_family" in out.columns:
            used_mate = rescue & prefer_mate
            out.loc[used_mate, "mate_mei_family"] = out.loc[used_mate, "top_family"]
            if "mate_mei_target" in out.columns:
                out.loc[used_mate, "mate_mei_target"] = out.loc[used_mate, "target"]
            if "mate_mei_score" in out.columns:
                out.loc[used_mate, "mate_mei_score"] = out.loc[used_mate, "mei_score"]
        _impute_rescue_mei_coords(out, rescue, kind="vntr", mapped=mapped_before)
        click.echo(
            f"[mei-annotate] VNTR-like rescue: flagged {n_rescue} discordant rows as VNTR_MAPPED "
            f"across {int(eligible.drop_duplicates(['chrom','window_start','window_end']).shape[0])} SVA-supported loci"
        )
    out = out.drop(columns=["_fam_norm", "top_family", "support_reads"], errors="ignore")
    return out


def _demote_polya_split_mei_hits(split_df: pd.DataFrame) -> pd.DataFrame:
    """Clear MEI hits on polyA/T junction clips (polyA_MAPPED only).

    Soft-clips with a polyA/T run (≥8) or ``poly_tail_rescued`` must not
    contribute to MEI_MAPPED, SR, side MEI coords, or family votes — even if a
    residual tip (or untrimmed fallback) remapped to consensus.
    """
    if split_df is None or split_df.empty:
        return split_df
    out = split_df.copy()
    polya = _split_polya_member_mask(out)
    if not bool(polya.any()):
        return out
    mei_hit = (
        out["mei_hit"].fillna(False).astype(bool)
        if "mei_hit" in out.columns
        else pd.Series(False, index=out.index)
    )
    mei_hit_coord = (
        out["mei_hit_coord"].fillna(False).astype(bool)
        if "mei_hit_coord" in out.columns
        else pd.Series(False, index=out.index)
    )
    demote = polya & (mei_hit | mei_hit_coord)
    n = int(demote.sum())
    if not n:
        return out
    out.loc[demote, "mei_hit"] = False
    if "mei_hit_coord" in out.columns:
        out.loc[demote, "mei_hit_coord"] = False
    if "short_mei_seed_rescued" in out.columns:
        out.loc[demote, "short_mei_seed_rescued"] = False
    for col in (
        "target",
        "family",
        "target_strand",
        "target_coord",
        "family_coord",
        "target_strand_coord",
    ):
        if col in out.columns:
            out.loc[demote, col] = ""
    for col in (
        "target_start",
        "target_end",
        "target_len",
        "alnlen",
        "mapq",
        "target_start_coord",
        "target_end_coord",
        "target_len_coord",
        "alnlen_coord",
        "mapq_coord",
    ):
        if col in out.columns:
            out.loc[demote, col] = 0
    for col in ("pid", "qcov", "mei_score", "pid_coord", "qcov_coord", "mei_score_coord"):
        if col in out.columns:
            out.loc[demote, col] = 0.0
    click.echo(
        f"[mei-annotate] polyA split demote: cleared MEI hits on {n} polyA/T soft-clips "
        "(counted as polyA_MAPPED only)"
    )
    return out


def _annotate_vntr_like_split_clips(
    split_df: pd.DataFrame,
    discordant_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Move CCCTCT/VNTR-like soft-clips out of MEI-SR into VNTR_MAPPED.

    - VNTR-like clips that minimap assigned to SVA are demoted (not MEI_MAPPED/SR).
    - VNTR-like non-hits at SVA-supported loci are rescued as ``vntr_rescue``.
    """
    if split_df is None or split_df.empty:
        return split_df
    out = split_df.copy()
    if "vntr_rescue" not in out.columns:
        out["vntr_rescue"] = False
    else:
        out["vntr_rescue"] = out["vntr_rescue"].fillna(False).astype(bool)
    if "vntr_like_score" not in out.columns:
        out["vntr_like_score"] = 0.0
    if "mei_hit" not in out.columns:
        out["mei_hit"] = False
    if "mei_hit_source" not in out.columns:
        out["mei_hit_source"] = ""

    seq = (
        out["clip_seq"].fillna("").astype(str)
        if "clip_seq" in out.columns
        else pd.Series("", index=out.index, dtype=str)
    )
    out["vntr_like_score"] = seq.map(
        lambda s: _sva_vntr_like_score(s, min_seq_len=_VNTR_MIN_CLIP_SEQ_LEN)
    ).astype(float)
    vntr_like = out["vntr_like_score"].fillna(0.0) >= 0.35
    if not bool(vntr_like.any()):
        return out

    fam = (
        out["family"].fillna("").astype(str).map(_family_from_target)
        if "family" in out.columns
        else pd.Series("", index=out.index, dtype=str)
    )
    if "target" in out.columns:
        empty_fam = fam.eq("") | fam.eq("OTHER")
        fam = fam.where(~empty_fam, out["target"].fillna("").astype(str).map(_family_from_target))
    mei_hit = out["mei_hit"].fillna(False).astype(bool)
    if "mei_hit_coord" in out.columns:
        mei_hit = mei_hit | out["mei_hit_coord"].fillna(False).astype(bool)

    # Always demote VNTR-like clips that aligned to SVA (hexamer smash-hits).
    demote = vntr_like & mei_hit & fam.isin(_VNTR_RESCUE_FAMILIES)
    # Rescue non-MEI VNTR-like clips at loci with SVA consensus support.
    sva_loci: set[tuple[str, int, int]] = set()
    key_cols = ["chrom", "window_start", "window_end"]
    if all(c in out.columns for c in key_cols):
        split_sva = out.loc[mei_hit & ~demote & fam.isin(_VNTR_RESCUE_FAMILIES), key_cols]
        for row in split_sva.itertuples(index=False):
            sva_loci.add((str(row.chrom), int(row.window_start), int(row.window_end)))
        if discordant_df is not None and not discordant_df.empty:
            d = discordant_df
            d_mapped = _discordant_row_mei_mapped(d) if not d.empty else pd.Series(dtype=bool)
            if d_mapped.any() and all(c in d.columns for c in key_cols):
                d_fam = (
                    d["family"].fillna("").astype(str).map(_family_from_target)
                    if "family" in d.columns
                    else pd.Series("", index=d.index, dtype=str)
                )
                if "mate_mei_family" in d.columns:
                    mf = d["mate_mei_family"].fillna("").astype(str).map(_family_from_target)
                    d_fam = d_fam.where(d_fam.isin(_VNTR_RESCUE_FAMILIES), mf)
                if "target" in d.columns:
                    empty = d_fam.eq("") | d_fam.eq("OTHER")
                    d_fam = d_fam.where(
                        ~empty, d["target"].fillna("").astype(str).map(_family_from_target)
                    )
                for row in d.loc[d_mapped & d_fam.isin(_VNTR_RESCUE_FAMILIES), key_cols].itertuples(
                    index=False
                ):
                    sva_loci.add((str(row.chrom), int(row.window_start), int(row.window_end)))
    at_sva_locus = pd.Series(False, index=out.index)
    if sva_loci and all(c in out.columns for c in key_cols):
        at_sva_locus = pd.Series(
            [
                (str(c), int(ws), int(we)) in sva_loci
                for c, ws, we in zip(out["chrom"], out["window_start"], out["window_end"])
            ],
            index=out.index,
        )
    rescue = vntr_like & (~mei_hit | demote) & (demote | at_sva_locus)
    rescue = rescue & ~out["vntr_rescue"]
    n = int(rescue.sum())
    if not n:
        return out

    out.loc[rescue, "vntr_rescue"] = True
    out.loc[rescue, "mei_hit_source"] = "vntr_rescue"
    # Remove from MEI-SR / MEI_MAPPED; keep as VNTR_MAPPED only.
    out.loc[rescue, "mei_hit"] = False
    if "mei_hit_coord" in out.columns:
        out.loc[rescue, "mei_hit_coord"] = False
    if "short_mei_seed_rescued" in out.columns:
        out.loc[rescue, "short_mei_seed_rescued"] = False
    out.loc[rescue, "family"] = "SVA"
    if "target" in out.columns:
        out.loc[rescue, "target"] = "SVA_VNTR_rescue#Retroposon/SVA"
    out.loc[rescue, "mei_score"] = out.loc[rescue, "vntr_like_score"].clip(lower=0.35, upper=0.85)
    click.echo(
        f"[mei-annotate] VNTR-like split rescue: flagged {n} soft-clips as VNTR_MAPPED "
        f"(demoted_mei_hits={int(demote.sum())})"
    )
    return out


_POLYA_MIN_SEQ_LEN = 25
_POLYA_MIN_RUN = 12
# Empirically (chr22 disease discordant mates with A/T homopolymer run ≥40):
#   whole-mate dominant-base purity median ≈ 0.91, p10 ≈ 0.79.
# Longest-substring search at 0.90 recovers those polyA/T bodies (median span
# ≈50 bp, full-read pure polyT → 151) with ~10% mismatch tolerance and low FP.
_POLYA_MIN_FRAC = 0.90
_POLYA_RESCUE_MIN_SUPPORT_READS = 3
_POLYA_RESCUE_MIN_FAMILY_PURITY = 0.60
# SVA also carries a 3' polyA tail (after the Alu-like region); include it here
# in addition to the separate central-VNTR rescue pass.
_POLYA_RESCUE_FAMILIES = frozenset({"ALU", "LINE1", "SVA"})
# VNTR-only schematic stub width for synthetic consensus coords (NOT a polyA length).
_VNTR_RESCUE_IMPUTE_WIDTH = 40
# One Illumina read ≈ 150 bp of observable polyA. A DPE with junction polyA
# clip + polyA mate can cover non-overlapping stretches; absolute ceiling is
# ~pair insert (300) minus a short uniquely mapped flank anchor (~20).
_MAX_POLYA_SINGLE_BP = 151
_MAX_POLYA_PAIR_BP = 280


def _longest_poly_at_span(
    seq: str,
    *,
    min_frac: float = _POLYA_MIN_FRAC,
    min_len: int = _POLYA_MIN_SEQ_LEN,
) -> tuple[int, float, str, str]:
    """Longest substring that is mostly polyA **or** mostly polyT.

    For each base in {A,T}, two-pointer search for the longest window with
    ``count(base) / window_len ≥ min_frac``. Returns
    ``(length, purity, base, span_seq)``. Length 0 if none found.

    If the span covers essentially the whole read (≥140 bp or within 2 bp of
    read length), that is a full-read polyA/T (few mismatches allowed).
    """
    s = "".join(ch for ch in (seq or "").upper() if ch in {"A", "C", "G", "T"})
    n = len(s)
    if n < int(min_len):
        return (0, 0.0, "", "")
    best_len = 0
    best_frac = 0.0
    best_base = ""
    best_ij = (0, 0)
    thr = float(min_frac)
    for base in ("A", "T"):
        left = 0
        n_base = 0
        for right in range(n):
            if s[right] == base:
                n_base += 1
            while left <= right and (n_base / float(right - left + 1)) < thr:
                if s[left] == base:
                    n_base -= 1
                left += 1
            cur_len = right - left + 1
            if cur_len >= int(min_len) and n_base > 0:
                frac = float(n_base) / float(cur_len)
                if cur_len > best_len or (cur_len == best_len and frac > best_frac):
                    best_len = cur_len
                    best_frac = frac
                    best_base = base
                    best_ij = (left, right + 1)
    if best_len <= 0:
        return (0, 0.0, "", "")
    span = s[best_ij[0] : best_ij[1]]
    # Whole-read polyA/T with a few end mismatches → report full read length.
    if best_len >= 140 or best_len >= n - 2:
        best_len = n
        best_frac = float(span.count(best_base)) / float(len(span)) if span else best_frac
        span = s
    return (int(best_len), float(best_frac), best_base, span)


def _polya_tail_width_bp(seq: str) -> int:
    """Observed polyA/T length = longest mostly-A or mostly-T span (else 0)."""
    length, _frac, _base, _span = _longest_poly_at_span(seq)
    if length < _POLYA_MIN_SEQ_LEN:
        return 0
    return int(min(length, _MAX_POLYA_SINGLE_BP))


def _homopolymer_at_run(seq: str) -> int:
    """Longest consecutive A-only or T-only run (legacy clip metric)."""
    s = (seq or "").upper()
    best = 0
    cur = 0
    prev = ""
    for ch in s:
        if ch not in {"A", "T"}:
            cur = 0
            prev = ""
            continue
        if ch == prev:
            cur += 1
        else:
            cur = 1
            prev = ch
        if cur > best:
            best = cur
    return int(best)


def _observed_poly_at_len_bp(seq: str) -> int:
    """Best observed polyA/T length in a clip/read (read-length lower bound).

    Takes the max of: 90%-purity span (≥25), looser span (≥8), and consecutive
    A/T homopolymer run. Used for consensus_poly_at_min_bp inputs.
    """
    s = seq or ""
    span25 = _polya_tail_width_bp(s)
    span8, _frac, _base, _span = _longest_poly_at_span(s, min_len=8)
    run = _homopolymer_at_run(s)
    return int(min(_MAX_POLYA_SINGLE_BP, max(int(span25), int(span8), int(run), 0)))


def _dpe_polya_observed_len_bp(
    *,
    mate_len: int,
    anchor_clip_len: int = 0,
    anchor_poly_run: int = 0,
    poly_tail_anchor_rescued: bool = False,
) -> int:
    """Observed polyA length for one DPE row (gold + plots).

    Mate-only → mate length (≤1 read). Junction polyA clip + polyA mate on the
    same pair → sum as a minimum non-overlapping tail (≤ pair insert − flank).
    """
    mate = max(0, int(mate_len))
    anchor = max(int(anchor_clip_len), int(anchor_poly_run), 0)
    anchor_ok = bool(poly_tail_anchor_rescued) or anchor >= 8
    mate_ok = mate >= _POLYA_MIN_SEQ_LEN
    if anchor_ok and mate_ok and anchor >= 8:
        return int(min(_MAX_POLYA_PAIR_BP, anchor + mate))
    if mate_ok:
        return int(min(_MAX_POLYA_SINGLE_BP, mate))
    if anchor_ok and anchor >= 8:
        return int(min(_MAX_POLYA_SINGLE_BP, anchor))
    return 0


def _impute_rescue_mei_coords(
    out: pd.DataFrame,
    rescue: pd.Series,
    *,
    kind: str,
    mapped: pd.Series | None = None,
) -> None:
    """Assign synthetic consensus coords so rescued mates appear in architecture plots.

    polyA → tack a polyA/T tail onto the locus 3' end (just past max mapped end).
    VNTR → central third of the mapped MEI footprint.
    """
    if not bool(rescue.any()):
        return
    mapped_mask = mapped if mapped is not None else out["mei_hit"].fillna(False).astype(bool)
    # Prefer real alignment ends; fall back to mate-side coords.
    end_cols = [c for c in ("target_end", "mate_mei_end") if c in out.columns]
    start_cols = [c for c in ("target_start", "mate_mei_start") if c in out.columns]
    if not end_cols:
        return

    mate_seq = (
        out["mate_seq"].fillna("").astype(str)
        if "mate_seq" in out.columns
        else pd.Series("", index=out.index)
    )
    read_seq = (
        out["read_seq"].fillna("").astype(str)
        if "read_seq" in out.columns
        else pd.Series("", index=out.index)
    )

    for (chrom, ws, we), idx in out.loc[rescue].groupby(
        ["chrom", "window_start", "window_end"]
    ).groups.items():
        locus = (
            (out["chrom"] == chrom)
            & (out["window_start"].astype(int) == int(ws))
            & (out["window_end"].astype(int) == int(we))
        )
        mapped_locus = out.loc[locus & mapped_mask]
        ends = []
        starts = []
        for col in end_cols:
            vals = pd.to_numeric(mapped_locus[col], errors="coerce").fillna(0).astype(int)
            ends.extend(int(v) for v in vals.tolist() if int(v) > 0)
        for col in start_cols:
            vals = pd.to_numeric(mapped_locus[col], errors="coerce").fillna(0).astype(int)
            starts.extend(int(v) for v in vals.tolist() if int(v) > 0)
        if not ends:
            continue
        hi = max(ends)
        lo = min(starts) if starts else max(1, hi - 200)
        if hi < lo:
            lo, hi = hi, lo
        span = max(1, hi - lo + 1)
        rescue_idx = out.index.intersection(idx)
        if kind == "polya":
            # Tack polyA onto the MEI 3' end: coords just past the consensus body.
            for i in rescue_idx:
                seq = str(mate_seq.loc[i] or "")
                if len(seq) < _POLYA_MIN_SEQ_LEN:
                    seq = str(read_seq.loc[i] or "")
                width = _polya_tail_width_bp(seq)
                if width <= 0:
                    # No inventing lengths — skip synthetic coords when unmeasurable.
                    continue
                start_i = int(hi) + 1
                end_i = int(hi) + int(width)
                if "target_start" in out.columns:
                    out.at[i, "target_start"] = start_i
                if "target_end" in out.columns:
                    out.at[i, "target_end"] = end_i
                if "mate_mei_start" in out.columns:
                    out.at[i, "mate_mei_start"] = start_i
                if "mate_mei_end" in out.columns:
                    out.at[i, "mate_mei_end"] = end_i
            continue
        # VNTR: central third of the mapped footprint (schematic stub only).
        width = min(_VNTR_RESCUE_IMPUTE_WIDTH, span)
        start_i = lo + span // 3
        end_i = lo + (2 * span) // 3
        if end_i < start_i:
            start_i, end_i = start_i, start_i + max(1, width - 1)
        if "target_start" in out.columns:
            out.loc[rescue_idx, "target_start"] = int(start_i)
        if "target_end" in out.columns:
            out.loc[rescue_idx, "target_end"] = int(end_i)
        if "mate_mei_start" in out.columns:
            out.loc[rescue_idx, "mate_mei_start"] = int(start_i)
        if "mate_mei_end" in out.columns:
            out.loc[rescue_idx, "mate_mei_end"] = int(end_i)


def _polya_like_stats(seq: str) -> tuple[float, int, str]:
    """Return (score, span_len, base) from the longest mostly-A/T span."""
    length, frac, base, span = _longest_poly_at_span(seq)
    if length < _POLYA_MIN_SEQ_LEN or not base:
        return (0.0, 0, "")
    run, _f, _b = _poly_at_stats(span)
    run_term = min(1.0, float(run) / 25.0)
    score = float(max(0.0, min(1.0, 0.55 * float(frac) + 0.45 * run_term)))
    if frac >= _POLYA_MIN_FRAC and length >= _POLYA_MIN_SEQ_LEN:
        score = max(score, 0.40)
    return (score, int(length), base)


def _rescue_polya_like_discordant_mei_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Second-pass polyA_MAPPED rescue for polyA/T-like discordant mates.

    Existing polyA logic rescues short split clips / anchors into evidence, but
    mates that are almost pure polyA/T often fail consensus remap. If the locus
    already has consistent ALU/LINE1/SVA support, flag polyA/T-like mates from
    either flank.

    Both flanks are allowed: for Alu-sized inserts a mate can land in the 3'
    polyA tail from either side, and the sequenced mate may be polyA or polyT
    depending on strand. Orientation (when known) is stored for labeling only.
    These count toward polyA_MAPPED (not MEI_MAPPED).
    """
    if df.empty:
        return df.copy()
    out = df.copy()
    if "mei_hit" not in out.columns:
        out["mei_hit"] = False
    out["mei_hit"] = out["mei_hit"].fillna(False).astype(bool)
    if "mei_hit_source" not in out.columns:
        out["mei_hit_source"] = ""
    out["mei_hit_source"] = out["mei_hit_source"].fillna("").astype(str)
    for col, default in (
        ("family", ""),
        ("target", ""),
        ("target_strand", ""),
        ("mei_score", 0.0),
        ("polya_like_score", 0.0),
        ("polya_rescue", False),
    ):
        if col not in out.columns:
            out[col] = default

    mate_seq = out["mate_seq"].fillna("").astype(str) if "mate_seq" in out.columns else pd.Series("", index=out.index)
    read_seq = out["read_seq"].fillna("").astype(str) if "read_seq" in out.columns else pd.Series("", index=out.index)
    prefer_mate = mate_seq.str.len() >= _POLYA_MIN_SEQ_LEN
    seqs = mate_seq.where(prefer_mate, read_seq)
    stats = seqs.map(_polya_like_stats)
    out["polya_like_score"] = stats.map(lambda x: float(x[0])).astype(float)
    out["_polya_base"] = stats.map(lambda x: str(x[2] or ""))
    consensus_mapped = _discordant_row_mei_mapped(out)
    already_rescue = out["polya_rescue"].fillna(False).astype(bool) | out.get(
        "vntr_rescue", pd.Series(False, index=out.index)
    ).fillna(False).astype(bool)
    poly_cand = (~consensus_mapped) & (~already_rescue) & (out["polya_like_score"] >= 0.40) & prefer_mate
    if not poly_cand.any():
        out = out.drop(columns=["_polya_base"], errors="ignore")
        return out

    mapped = consensus_mapped
    fam = out["family"].fillna("").astype(str).map(_family_from_target)
    empty_fam = fam.eq("") | fam.eq("OTHER")
    if empty_fam.any() and "target" in out.columns:
        fam = fam.where(~empty_fam, out["target"].fillna("").astype(str).map(_family_from_target))
    out["_fam_norm"] = fam

    # Locus family + orientation from consensus MEI hits.
    locus_rows = []
    mapped_rows = out.loc[mapped].copy()
    if mapped_rows.empty:
        out = out.drop(columns=["_fam_norm", "_polya_base"], errors="ignore")
        return out
    for (chrom, ws, we), grp in mapped_rows.groupby(["chrom", "window_start", "window_end"]):
        vc = grp["_fam_norm"].value_counts()
        top = str(vc.index[0]) if len(vc) else ""
        pur = float(vc.iloc[0] / max(len(grp), 1)) if len(vc) else 0.0
        support_reads = int(grp["read_name"].nunique())
        strand_vc = (
            grp["target_strand"].fillna("").astype(str).replace({"": pd.NA}).dropna().value_counts()
            if "target_strand" in grp.columns
            else pd.Series(dtype=int)
        )
        if "mate_mei_strand" in grp.columns and strand_vc.empty:
            strand_vc = grp["mate_mei_strand"].fillna("").astype(str).replace({"": pd.NA}).dropna().value_counts()
        orient = ""
        if len(strand_vc):
            top_strand = str(strand_vc.index[0])
            if top_strand in {"+", "-"} and float(strand_vc.iloc[0] / max(int(strand_vc.sum()), 1)) >= 0.55:
                orient = top_strand
        locus_rows.append(
            {
                "chrom": chrom,
                "window_start": int(ws),
                "window_end": int(we),
                "top_family": top,
                "family_purity": pur,
                "support_reads": support_reads,
                "locus_orientation": orient,
            }
        )
    support = pd.DataFrame(locus_rows)
    # Orientation is optional (used for labeling only); both flanks are eligible.
    eligible = support.loc[
        (support["support_reads"] >= _POLYA_RESCUE_MIN_SUPPORT_READS)
        & (support["family_purity"].fillna(0.0) >= _POLYA_RESCUE_MIN_FAMILY_PURITY)
        & (support["top_family"].fillna("").astype(str).isin(_POLYA_RESCUE_FAMILIES))
    ].copy()
    if eligible.empty:
        out = out.drop(columns=["_fam_norm", "_polya_base"], errors="ignore")
        return out

    out = out.merge(
        eligible[
            ["chrom", "window_start", "window_end", "top_family", "support_reads", "locus_orientation"]
        ],
        on=["chrom", "window_start", "window_end"],
        how="left",
        suffixes=("", "_elig"),
    )
    prefer_mate = (
        out["mate_seq"].fillna("").astype(str).str.len() >= _POLYA_MIN_SEQ_LEN
        if "mate_seq" in out.columns
        else pd.Series(False, index=out.index)
    )

    already_rescue = out["polya_rescue"].fillna(False).astype(bool)
    if "vntr_rescue" in out.columns:
        already_rescue = already_rescue | out["vntr_rescue"].fillna(False).astype(bool)
    poly_cand = (
        (~_discordant_row_mei_mapped(out))
        & (~already_rescue)
        & (out["polya_like_score"].fillna(0.0) >= 0.40)
        & prefer_mate
    )
    rescue = poly_cand & out["top_family"].fillna("").astype(str).isin(_POLYA_RESCUE_FAMILIES)
    n_rescue = int(rescue.sum())
    if n_rescue:
        mapped_before = _discordant_row_mei_mapped(out)
        # Do NOT set mei_hit — polyA rescues are counted as polyA_MAPPED separately.
        out.loc[rescue, "polya_rescue"] = True
        out.loc[rescue, "mei_hit_source"] = "polya_rescue"
        out.loc[rescue, "family"] = out.loc[rescue, "top_family"]
        out.loc[rescue, "target"] = out.loc[rescue, "top_family"].map(
            lambda f: (
                f"{f}_polyA_rescue#SINE/Alu"
                if f == "ALU"
                else f"{f}_polyA_rescue#LINE/L1"
                if f == "LINE1"
                else f"{f}_polyA_rescue#Retroposon/SVA"
                if f == "SVA"
                else f"{f}_polyA_rescue"
            )
        )
        # Label with locus orientation when known; leave blank otherwise.
        orient = out["locus_orientation"].fillna("").astype(str)
        out.loc[rescue & orient.isin({"+", "-"}), "target_strand"] = orient.loc[
            rescue & orient.isin({"+", "-"})
        ]
        out.loc[rescue, "mei_score"] = out.loc[rescue, "polya_like_score"].clip(lower=0.40, upper=0.85)
        if "mate_mei_family" in out.columns:
            out.loc[rescue, "mate_mei_family"] = out.loc[rescue, "top_family"]
            if "mate_mei_target" in out.columns:
                out.loc[rescue, "mate_mei_target"] = out.loc[rescue, "target"]
            if "mate_mei_strand" in out.columns:
                out.loc[rescue & orient.isin({"+", "-"}), "mate_mei_strand"] = orient.loc[
                    rescue & orient.isin({"+", "-"})
                ]
            if "mate_mei_score" in out.columns:
                out.loc[rescue, "mate_mei_score"] = out.loc[rescue, "mei_score"]
        _impute_rescue_mei_coords(out, rescue, kind="polya", mapped=mapped_before)
        # Side breakdown for logging (both flanks allowed).
        mid = (out["window_start"].astype(int) + out["window_end"].astype(int)) // 2
        sides = (
            out.loc[rescue, "pos"]
            .astype(int)
            .le(mid.loc[rescue])
            .map({True: "L", False: "R"})
            .value_counts()
            .to_dict()
        )
        click.echo(
            f"[mei-annotate] polyA-like rescue: flagged {n_rescue} discordant rows as polyA_MAPPED "
            f"across {int(eligible.drop_duplicates(['chrom','window_start','window_end']).shape[0])} "
            f"ALU/LINE1/SVA-supported loci (both flanks; sides={sides})"
        )
    out = out.drop(
        columns=["_fam_norm", "_polya_base", "top_family", "support_reads", "locus_orientation"],
        errors="ignore",
    )
    return out


def _aggregate_discordant_residual_complex_metrics(df: pd.DataFrame, sample_prefix: str) -> pd.DataFrame:
    """
    Complex SV reason fractions among discordants that did NOT map to MEI.

    Classic MEI insertions often have interchrom/large-insert mates that remap to
    MEI consensus; those must not drive complex_mei_event. Only the non-MEI
    residual can support a companion SV signature.
    """
    empty_cols = [
        "chrom",
        "window_start",
        "window_end",
        f"{sample_prefix}_discordant_mei_mapped_unique_reads",
        f"{sample_prefix}_discordant_residual_unique_reads",
        f"{sample_prefix}_discordant_mei_mapped_fraction",
        f"{sample_prefix}_discordant_residual_interchrom_fraction",
        f"{sample_prefix}_discordant_residual_large_insert_fraction",
        f"{sample_prefix}_discordant_residual_mate_unmapped_fraction",
        f"{sample_prefix}_discordant_residual_same_strand_fraction",
        f"{sample_prefix}_discordant_residual_improper_pair_fraction",
        f"{sample_prefix}_discordant_residual_complex_any_fraction",
    ]
    if df.empty:
        return pd.DataFrame(columns=empty_cols)

    tmp = df.copy()
    # Consensus MEI hits + polyA/VNTR rescues are MEI-related evidence; only the
    # remaining residual can support a companion complex SV signature.
    related = _discordant_row_mei_related(tmp)
    tmp["_mei_mapped"] = _discordant_row_mei_mapped(tmp).astype(bool)
    tmp["_mei_related"] = related.astype(bool)
    reason_col = tmp["discordant_reasons"].fillna("").astype(str) if "discordant_reasons" in tmp.columns else pd.Series("", index=tmp.index)
    tmp["reason_interchrom"] = reason_col.str.contains("interchrom", regex=False).astype(int)
    tmp["reason_mate_unmapped"] = reason_col.str.contains("mate_unmapped", regex=False).astype(int)
    tmp["reason_large_insert"] = reason_col.str.contains("large_insert", regex=False).astype(int)
    tmp["reason_same_strand"] = reason_col.str.contains("same_strand", regex=False).astype(int)
    tmp["reason_improper_pair"] = reason_col.str.contains("improper_pair", regex=False).astype(int)
    tmp["reason_complex_any"] = (
        (tmp["reason_interchrom"] == 1)
        | (tmp["reason_mate_unmapped"] == 1)
        | (tmp["reason_large_insert"] == 1)
        | (tmp["reason_same_strand"] == 1)
        | (tmp["reason_improper_pair"] == 1)
    ).astype(int)

    mei_u = (
        tmp.loc[tmp["_mei_mapped"]]
        .groupby(["chrom", "window_start", "window_end"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": f"{sample_prefix}_discordant_mei_mapped_unique_reads"})
    )
    res_u = (
        tmp.loc[~tmp["_mei_related"]]
        .groupby(["chrom", "window_start", "window_end"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": f"{sample_prefix}_discordant_residual_unique_reads"})
    )
    all_u = (
        tmp.groupby(["chrom", "window_start", "window_end"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": "_total_unique"})
    )
    residual = tmp.loc[~tmp["_mei_related"]].copy()
    if residual.empty:
        out = all_u.rename(columns={"_total_unique": "_total_unique"})
        out = out.merge(mei_u, on=["chrom", "window_start", "window_end"], how="left")
        out = out.merge(res_u, on=["chrom", "window_start", "window_end"], how="left")
        out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"] = (
            out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"].fillna(0).astype(int)
        )
        out[f"{sample_prefix}_discordant_residual_unique_reads"] = 0
        out[f"{sample_prefix}_discordant_mei_mapped_fraction"] = (
            out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"].astype(float)
            / out["_total_unique"].astype(float).clip(lower=1.0)
        )
        for col in [
            f"{sample_prefix}_discordant_residual_interchrom_fraction",
            f"{sample_prefix}_discordant_residual_large_insert_fraction",
            f"{sample_prefix}_discordant_residual_mate_unmapped_fraction",
            f"{sample_prefix}_discordant_residual_same_strand_fraction",
            f"{sample_prefix}_discordant_residual_improper_pair_fraction",
            f"{sample_prefix}_discordant_residual_complex_any_fraction",
        ]:
            out[col] = 0.0
        return out.drop(columns=["_total_unique"])

    res_frac = (
        residual.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_discordant_residual_interchrom_fraction": ("reason_interchrom", "mean"),
                f"{sample_prefix}_discordant_residual_large_insert_fraction": ("reason_large_insert", "mean"),
                f"{sample_prefix}_discordant_residual_mate_unmapped_fraction": ("reason_mate_unmapped", "mean"),
                f"{sample_prefix}_discordant_residual_same_strand_fraction": ("reason_same_strand", "mean"),
                f"{sample_prefix}_discordant_residual_improper_pair_fraction": ("reason_improper_pair", "mean"),
                f"{sample_prefix}_discordant_residual_complex_any_fraction": ("reason_complex_any", "mean"),
            }
        )
    )
    out = all_u.merge(mei_u, on=["chrom", "window_start", "window_end"], how="left")
    out = out.merge(res_u, on=["chrom", "window_start", "window_end"], how="left")
    out = out.merge(res_frac, on=["chrom", "window_start", "window_end"], how="left")
    out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"] = (
        out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"].fillna(0).astype(int)
    )
    out[f"{sample_prefix}_discordant_residual_unique_reads"] = (
        out[f"{sample_prefix}_discordant_residual_unique_reads"].fillna(0).astype(int)
    )
    out[f"{sample_prefix}_discordant_mei_mapped_fraction"] = (
        out[f"{sample_prefix}_discordant_mei_mapped_unique_reads"].astype(float)
        / out["_total_unique"].astype(float).clip(lower=1.0)
    )
    for col in [
        f"{sample_prefix}_discordant_residual_interchrom_fraction",
        f"{sample_prefix}_discordant_residual_large_insert_fraction",
        f"{sample_prefix}_discordant_residual_mate_unmapped_fraction",
        f"{sample_prefix}_discordant_residual_same_strand_fraction",
        f"{sample_prefix}_discordant_residual_improper_pair_fraction",
        f"{sample_prefix}_discordant_residual_complex_any_fraction",
    ]:
        out[col] = out[col].fillna(0.0).astype(float)
    return out.drop(columns=["_total_unique"])


def _enrich_split_hits_with_mate_positions(
    split_hits: pd.DataFrame,
    bam_path: Path | None,
) -> pd.DataFrame:
    """Attach mate_chrom/mate_pos to split hits from evidence columns or BAM fallback."""
    if split_hits.empty:
        return split_hits.copy()

    out = split_hits.copy()
    if "mate_chrom" not in out.columns:
        out["mate_chrom"] = ""
    if "mate_pos" not in out.columns:
        out["mate_pos"] = 0
    out["mate_chrom"] = out["mate_chrom"].fillna("").astype(str)
    out["mate_pos"] = pd.to_numeric(out["mate_pos"], errors="coerce").fillna(0).astype(int)

    if bam_path is None or not bam_path.exists():
        return out

    need_fetch = out.loc[(out["mate_chrom"].isin({"", "*"})) | (out["mate_pos"] <= 0)].copy()
    if need_fetch.empty:
        return out

    fetched: dict[str, tuple[str, int]] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in need_fetch.itertuples(index=False):
            qname = str(rec.read_name)
            if qname in fetched:
                continue
            chrom = str(rec.chrom)
            pos = int(rec.pos)
            if chrom not in bam.references:
                continue
            start0 = max(0, pos - 1)
            end0 = start0 + 500
            for read in bam.fetch(chrom, start0, end0):
                if read.query_name != qname:
                    continue
                if read.is_secondary or read.is_supplementary:
                    continue
                if not read.is_paired or read.mate_is_unmapped:
                    continue
                mate_chrom = "*"
                mate_pos = 0
                if read.next_reference_id >= 0:
                    mate_chrom = bam.get_reference_name(read.next_reference_id)
                    mate_pos = read.next_reference_start + 1
                fetched[qname] = (mate_chrom, mate_pos)
                break

    for qname, (mate_chrom, mate_pos) in fetched.items():
        mask = out["read_name"].astype(str) == qname
        out.loc[mask, "mate_chrom"] = mate_chrom
        out.loc[mask, "mate_pos"] = mate_pos
    return out


def _split_clip_mei_fields(split_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize per-split MEI alignment fields used by short-clip rescue."""
    mei_hit = (
        split_df["mei_hit"].fillna(False).astype(bool)
        if "mei_hit" in split_df.columns
        else pd.Series(False, index=split_df.index)
    )
    mei_hit_coord = (
        split_df["mei_hit_coord"].fillna(False).astype(bool)
        if "mei_hit_coord" in split_df.columns
        else pd.Series(False, index=split_df.index)
    )
    has_mei = mei_hit | mei_hit_coord
    clip_len = pd.to_numeric(split_df.get("clip_len", 0), errors="coerce").fillna(0).astype(int)
    if "target_start_coord" in split_df.columns:
        t_start = pd.to_numeric(split_df["target_start_coord"], errors="coerce").fillna(0).astype(int)
        t_end = pd.to_numeric(split_df["target_end_coord"], errors="coerce").fillna(0).astype(int)
        t_strand = split_df.get("target_strand_coord", pd.Series("", index=split_df.index)).fillna("").astype(str)
        t_family = split_df.get("family_coord", pd.Series("", index=split_df.index)).fillna("").astype(str)
        t_aln = pd.to_numeric(split_df.get("alnlen_coord", 0), errors="coerce").fillna(0).astype(int)
        use_coord = (t_start > 0) | (t_end > 0)
        if "target_start" in split_df.columns:
            t_start = t_start.where(
                use_coord, pd.to_numeric(split_df["target_start"], errors="coerce").fillna(0).astype(int)
            )
            t_end = t_end.where(
                use_coord, pd.to_numeric(split_df["target_end"], errors="coerce").fillna(0).astype(int)
            )
            raw_strand = split_df.get("target_strand", pd.Series("", index=split_df.index)).fillna("").astype(str)
            raw_family = split_df.get("family", pd.Series("", index=split_df.index)).fillna("").astype(str)
            raw_aln = pd.to_numeric(split_df.get("alnlen", 0), errors="coerce").fillna(0).astype(int)
            t_strand = t_strand.where(use_coord & t_strand.ne(""), raw_strand)
            t_family = t_family.where(use_coord & t_family.ne(""), raw_family)
            t_aln = t_aln.where(use_coord & t_aln.gt(0), raw_aln)
    else:
        t_start = pd.to_numeric(split_df.get("target_start", 0), errors="coerce").fillna(0).astype(int)
        t_end = pd.to_numeric(split_df.get("target_end", 0), errors="coerce").fillna(0).astype(int)
        t_strand = split_df.get("target_strand", pd.Series("", index=split_df.index)).fillna("").astype(str)
        t_family = split_df.get("family", pd.Series("", index=split_df.index)).fillna("").astype(str)
        t_aln = pd.to_numeric(split_df.get("alnlen", 0), errors="coerce").fillna(0).astype(int)
    side = (
        split_df["clip_side"].fillna("").astype(str).str.upper().str[:1]
        if "clip_side" in split_df.columns
        else pd.Series("", index=split_df.index)
    )
    work = pd.DataFrame(
        {
            "clip_len": clip_len,
            "has_mei": has_mei,
            "t_strand": t_strand,
            "t_family": t_family,
            "t_aln": t_aln,
            "side": side,
            "t_lo": pd.concat([t_start, t_end], axis=1).min(axis=1),
            "t_hi": pd.concat([t_start, t_end], axis=1).max(axis=1),
        },
        index=split_df.index,
    )
    return work


def _dpe_mei_seed_table(discordant_df: pd.DataFrame, key_cols: list[str]) -> pd.DataFrame:
    """Build per-locus/side DPE MEI seeds for short-clip rescue."""
    empty_cols = key_cols + [
        "side",
        "seed_lo",
        "seed_hi",
        "seed_mid",
        "seed_strand",
        "seed_family",
        "seed_n",
        "seed_source",
    ]
    if discordant_df is None or discordant_df.empty:
        return pd.DataFrame(columns=empty_cols)
    if not all(c in discordant_df.columns for c in key_cols + ["pos"]):
        return pd.DataFrame(columns=empty_cols)

    mapped = _discordant_row_mei_mapped(discordant_df)
    if not mapped.any():
        return pd.DataFrame(columns=empty_cols)
    d = discordant_df.loc[mapped].copy()
    mate_hit = (
        d["mate_mei_hit"].fillna(False).astype(bool)
        if "mate_mei_hit" in d.columns
        else pd.Series(False, index=d.index)
    )
    # Prefer mate MEI coords (classic DPE→MEI); fall back to anchor remap.
    t_start = (
        pd.to_numeric(d["target_start"], errors="coerce").fillna(0).astype(int)
        if "target_start" in d.columns
        else pd.Series(0, index=d.index, dtype=int)
    )
    t_end = (
        pd.to_numeric(d["target_end"], errors="coerce").fillna(0).astype(int)
        if "target_end" in d.columns
        else pd.Series(0, index=d.index, dtype=int)
    )
    t_strand = (
        d["target_strand"].fillna("").astype(str)
        if "target_strand" in d.columns
        else pd.Series("", index=d.index, dtype=str)
    )
    t_family = (
        d["family"].fillna("").astype(str)
        if "family" in d.columns
        else pd.Series("", index=d.index, dtype=str)
    )
    if "mate_mei_start" in d.columns:
        m_start = pd.to_numeric(d["mate_mei_start"], errors="coerce").fillna(0).astype(int)
        m_end = (
            pd.to_numeric(d["mate_mei_end"], errors="coerce").fillna(0).astype(int)
            if "mate_mei_end" in d.columns
            else pd.Series(0, index=d.index, dtype=int)
        )
        m_strand = (
            d["mate_mei_strand"].fillna("").astype(str)
            if "mate_mei_strand" in d.columns
            else pd.Series("", index=d.index, dtype=str)
        )
        m_family = (
            d["mate_mei_family"].fillna("").astype(str)
            if "mate_mei_family" in d.columns
            else pd.Series("", index=d.index, dtype=str)
        )
        t_start = m_start.where(mate_hit & ((m_start > 0) | (m_end > 0)), t_start)
        t_end = m_end.where(mate_hit & ((m_start > 0) | (m_end > 0)), t_end)
        t_strand = m_strand.where(mate_hit & m_strand.ne(""), t_strand)
        t_family = m_family.where(mate_hit & m_family.ne(""), t_family)
    d["_lo"] = pd.concat([t_start, t_end], axis=1).min(axis=1)
    d["_hi"] = pd.concat([t_start, t_end], axis=1).max(axis=1)
    d = d.loc[d["_hi"].gt(0)].copy()
    if d.empty:
        return pd.DataFrame(columns=empty_cols)

    mid = (
        pd.to_numeric(d["window_start"], errors="coerce").fillna(0).astype(int)
        + pd.to_numeric(d["window_end"], errors="coerce").fillna(0).astype(int)
    ) // 2
    pos = pd.to_numeric(d["pos"], errors="coerce").fillna(mid).astype(int)
    d["_side"] = pd.Series(["L"] * len(d), index=d.index).where(pos <= mid, "R")
    d["_strand"] = t_strand.loc[d.index]
    d["_family"] = t_family.loc[d.index]

    seeds = (
        d.groupby(key_cols + ["_side"], as_index=False)
        .agg(
            seed_lo=("_lo", "min"),
            seed_hi=("_hi", "max"),
            seed_strand=("_strand", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
            seed_family=("_family", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
            seed_n=("pos", "count"),
        )
        .rename(columns={"_side": "side"})
    )
    seeds = seeds.loc[seeds["seed_n"] >= int(_SHORT_MEI_DPE_MIN_SUPPORT)].copy()
    if seeds.empty:
        return pd.DataFrame(columns=empty_cols)
    seeds["seed_mid"] = (seeds["seed_lo"] + seeds["seed_hi"]) / 2.0
    seeds["seed_source"] = "dpe"
    return seeds.loc[:, empty_cols]


def _annotate_short_mei_seed_rescue(
    split_df: pd.DataFrame,
    discordant_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Mark short MEI-mapped clips consistent with SR≥20 or DPE MEI seeds.

    Priority per locus/side:
      1. Strict SR MEI clips (≥20bp): short clip must overlap seed MEI coords.
      2. Else DPE MEI (≥``_SHORT_MEI_DPE_MIN_SUPPORT`` mates): same family/strand,
         near DPE span along the MEI axis, and junction-proximal to the DPE
         cluster (clips need not overlap DPE mate coords).
    """
    if split_df is None or split_df.empty:
        return split_df
    out = split_df.copy()
    out["short_mei_seed_rescued"] = False
    if "clip_len" not in out.columns:
        return out
    key_cols = ["chrom", "window_start", "window_end"]
    if not all(c in out.columns for c in key_cols):
        return out

    work = _split_clip_mei_fields(out)
    for c in key_cols:
        work[c] = out[c]
    work["t_mid"] = (work["t_lo"] + work["t_hi"]) / 2.0

    vntr_blocked = (
        out["vntr_rescue"].fillna(False).astype(bool)
        if "vntr_rescue" in out.columns
        else pd.Series(False, index=out.index)
    )
    short_mask = (
        work["has_mei"]
        & ~vntr_blocked
        & work["clip_len"].ge(int(_SHORT_MEI_RESCUE_MIN_CLIP_BP))
        & work["clip_len"].lt(int(_STRICT_MEI_CLIP_MIN_BP))
        & work["t_aln"].ge(int(_SHORT_MEI_RESCUE_MIN_ALN_BP))
        & work["side"].isin(["L", "R"])
        & work["t_hi"].gt(0)
    )
    if not short_mask.any():
        return out

    rescued = pd.Series(False, index=out.index)

    # --- 1) Strict SR seeds (overlap) ---
    sr_seed_mask = (
        work["has_mei"]
        & ~vntr_blocked
        & work["clip_len"].ge(int(_STRICT_MEI_CLIP_MIN_BP))
        & work["t_aln"].ge(int(_SHORT_MEI_RESCUE_MIN_ALN_BP))
        & work["side"].isin(["L", "R"])
        & work["t_hi"].gt(0)
    )
    sr_sides: set[tuple[str, int, int, str]] = set()
    if sr_seed_mask.any():
        sr_seeds = (
            work.loc[sr_seed_mask]
            .groupby(key_cols + ["side"], as_index=False)
            .agg(
                seed_lo=("t_lo", "min"),
                seed_hi=("t_hi", "max"),
                seed_strand=("t_strand", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
                seed_family=("t_family", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
            )
        )
        for row in sr_seeds.itertuples(index=False):
            sr_sides.add((str(row.chrom), int(row.window_start), int(row.window_end), str(row.side)))
        cand = work.loc[short_mask].reset_index().merge(sr_seeds, on=key_cols + ["side"], how="inner")
        if not cand.empty:
            gap = int(_SHORT_MEI_RESCUE_MAX_COORD_GAP_BP)
            strand_ok = cand["t_strand"].eq(cand["seed_strand"]) | cand["seed_strand"].eq("") | cand["t_strand"].eq("")
            fam_ok = (
                cand["t_family"].eq(cand["seed_family"])
                | cand["seed_family"].eq("")
                | cand["t_family"].eq("")
            )
            overlap = (cand["t_lo"] <= cand["seed_hi"] + gap) & (cand["t_hi"] >= cand["seed_lo"] - gap)
            hit_idx = cand.loc[strand_ok & fam_ok & overlap, "index"].astype(int)
            rescued.loc[hit_idx] = True

    # --- 2) DPE seeds on sides without SR≥20 ---
    still_short = short_mask & ~rescued
    if still_short.any() and discordant_df is not None and not discordant_df.empty:
        dpe_seeds = _dpe_mei_seed_table(discordant_df, key_cols)
        if not dpe_seeds.empty:
            if sr_sides:
                keep_rows = []
                for row in dpe_seeds.itertuples(index=False):
                    key = (str(row.chrom), int(row.window_start), int(row.window_end), str(row.side))
                    if key not in sr_sides:
                        keep_rows.append(True)
                    else:
                        keep_rows.append(False)
                dpe_seeds = dpe_seeds.loc[keep_rows].copy()
            if not dpe_seeds.empty:
                cand = work.loc[still_short].reset_index().merge(
                    dpe_seeds, on=key_cols + ["side"], how="inner"
                )
                if not cand.empty:
                    axis_gap = int(_SHORT_MEI_DPE_MAX_AXIS_GAP_BP)
                    slack = int(_SHORT_MEI_DPE_PROXIMAL_SLACK_BP)
                    strand_ok = (
                        cand["t_strand"].eq(cand["seed_strand"])
                        | cand["seed_strand"].eq("")
                        | cand["t_strand"].eq("")
                    )
                    fam_ok = (
                        cand["t_family"].eq(cand["seed_family"])
                        | cand["seed_family"].eq("")
                        | cand["t_family"].eq("")
                    )
                    # Distance between intervals along MEI (0 if overlap).
                    sep = pd.concat(
                        [
                            cand["t_lo"] - cand["seed_hi"],
                            cand["seed_lo"] - cand["t_hi"],
                        ],
                        axis=1,
                    ).max(axis=1).clip(lower=0)
                    near = sep <= axis_gap
                    # Junction-proximal: clip toward the flank entry end vs DPE cluster.
                    # +/L or -/R → enter at MEI 5' (lower coords); else enter at 3'.
                    entry_5p = (
                        (cand["side"].eq("L") & cand["seed_strand"].isin(["+", ""]))
                        | (cand["side"].eq("R") & cand["seed_strand"].eq("-"))
                    )
                    proximal = pd.Series(False, index=cand.index)
                    if entry_5p.any():
                        proximal.loc[entry_5p] = cand.loc[entry_5p, "t_mid"] <= (
                            cand.loc[entry_5p, "seed_mid"] + slack
                        )
                    if (~entry_5p).any():
                        proximal.loc[~entry_5p] = cand.loc[~entry_5p, "t_mid"] >= (
                            cand.loc[~entry_5p, "seed_mid"] - slack
                        )
                    hit_idx = cand.loc[strand_ok & fam_ok & near & proximal, "index"].astype(int)
                    rescued.loc[hit_idx] = True

    out.loc[rescued, "short_mei_seed_rescued"] = True
    return out


def _split_mei_support_eligible_mask(split_df: pd.DataFrame) -> pd.Series:
    """True for MEI-mapped splits eligible for SR / MEI_MAPPED counting."""
    if split_df is None or split_df.empty:
        return pd.Series(dtype=bool)
    mei_hit = (
        split_df["mei_hit"].fillna(False).astype(bool)
        if "mei_hit" in split_df.columns
        else pd.Series(False, index=split_df.index)
    )
    mei_hit_coord = (
        split_df["mei_hit_coord"].fillna(False).astype(bool)
        if "mei_hit_coord" in split_df.columns
        else pd.Series(False, index=split_df.index)
    )
    has_mei = mei_hit | mei_hit_coord
    if "vntr_rescue" in split_df.columns:
        has_mei = has_mei & ~split_df["vntr_rescue"].fillna(False).astype(bool)
    if "mei_hit_source" in split_df.columns:
        src = split_df["mei_hit_source"].fillna("").astype(str)
        has_mei = has_mei & ~src.eq("vntr_rescue")
    # PolyA/T junction clips are polyA_MAPPED only — never also MEI_MAPPED/SR.
    has_mei = has_mei & ~_split_polya_member_mask(split_df)
    if "clip_len" not in split_df.columns:
        return has_mei
    clip_len = pd.to_numeric(split_df["clip_len"], errors="coerce").fillna(int(_STRICT_MEI_CLIP_MIN_BP)).astype(int)
    rescued = (
        split_df["short_mei_seed_rescued"].fillna(False).astype(bool)
        if "short_mei_seed_rescued" in split_df.columns
        else pd.Series(False, index=split_df.index)
    )
    return has_mei & (clip_len.ge(int(_STRICT_MEI_CLIP_MIN_BP)) | rescued)


def _split_polya_member_mask(split_df: pd.DataFrame) -> pd.Series:
    """True for split rows that contribute to polyA_MAPPED."""
    if split_df is None or split_df.empty:
        return pd.Series(dtype=bool)
    mask = pd.Series(False, index=split_df.index)
    if "poly_tail_rescued" in split_df.columns:
        mask = mask | split_df["poly_tail_rescued"].fillna(False).infer_objects(copy=False).astype(bool)
    if "clip_poly_at_run" in split_df.columns:
        mask = mask | (
            pd.to_numeric(split_df["clip_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
        )
    return mask


def _build_supporting_reads_detail_table(
    *,
    split_hits: pd.DataFrame,
    discordant_hits: pd.DataFrame,
    discordant_mate_hits: pd.DataFrame,
    sample: str,
) -> pd.DataFrame:
    """Flatten per-read coordinates for architecture plots and review.

    Includes MEI-eligible splits (strict or short-seed-rescued) and polyA-
    supporting splits so read-architecture plots match the support-string evidence.
    """
    rows: list[dict[str, object]] = []

    def _as_int(val: object, default: int = 0) -> int:
        try:
            if val is None or (isinstance(val, float) and pd.isna(val)):
                return default
            return int(val)
        except (TypeError, ValueError):
            return default

    def _as_bool(val: object, default: bool = False) -> bool:
        if val is None:
            return default
        try:
            if pd.isna(val):
                return default
        except (TypeError, ValueError):
            pass
        return bool(val)

    if not split_hits.empty:
        mei_eligible = _split_mei_support_eligible_mask(split_hits)
        polya_mask = _split_polya_member_mask(split_hits)
        for i, rec in enumerate(split_hits.itertuples(index=False)):
            idx = split_hits.index[i]
            mei_hit = bool(mei_eligible.loc[idx]) if idx in mei_eligible.index else False
            polya_split = bool(polya_mask.loc[idx]) if idx in polya_mask.index else False
            if not (mei_hit or polya_split):
                continue
            start_col = (
                "target_start_coord"
                if _as_int(getattr(rec, "target_start_coord", 0), 0) > 0
                else "target_start"
            )
            end_col = (
                "target_end_coord"
                if _as_int(getattr(rec, "target_end_coord", 0), 0) > 0
                else "target_end"
            )
            strand = str(getattr(rec, "target_strand", "") or "")
            strand_coord = str(getattr(rec, "target_strand_coord", "") or "")
            if strand not in {"+", "-"} and strand_coord in {"+", "-"}:
                strand = strand_coord
            clip_len = _as_int(getattr(rec, "clip_len", 0), 0)
            if clip_len <= 0:
                clip_len = _as_int(getattr(rec, "soft_clip_len", 0), 0)
            poly_run = _as_int(getattr(rec, "clip_poly_at_run", 0), 0)
            poly_rescued = _as_bool(getattr(rec, "poly_tail_rescued", False)) or polya_split
            # Only persist A/T base for true polyA contributors (not every softclip).
            poly_base = ""
            if polya_split:
                poly_base = str(getattr(rec, "clip_poly_base", "") or "").upper()
                if poly_base in {"NAN", ""} or poly_base not in {"A", "T"}:
                    poly_base = ""
            hit_source = ""
            if mei_hit:
                hit_source = (
                    "short_mei_rescue"
                    if _as_bool(getattr(rec, "short_mei_seed_rescued", False))
                    else "mei"
                )
            elif polya_split:
                hit_source = "polya_clip"
            rows.append(
                {
                    "sample": sample,
                    "evidence_type": "SR",
                    "read_name": str(rec.read_name),
                    "chrom": str(rec.chrom),
                    "window_start": int(rec.window_start),
                    "window_end": int(rec.window_end),
                    "anchor_side": str(getattr(rec, "clip_side", "") or ""),
                    "genomic_pos": int(rec.pos),
                    "mate_chrom": str(getattr(rec, "mate_chrom", "") or ""),
                    "mate_genomic_pos": _as_int(getattr(rec, "mate_pos", 0), 0),
                    "mei_start": _as_int(getattr(rec, start_col, 0), 0) if mei_hit else 0,
                    "mei_end": _as_int(getattr(rec, end_col, 0), 0) if mei_hit else 0,
                    "mate_mei_start": 0,
                    "mate_mei_end": 0,
                    "mei_target": str(getattr(rec, "target", "") or "") if mei_hit else "",
                    "mate_mei_target": "",
                    "mei_strand": strand if mei_hit and strand in {"+", "-"} else "",
                    "mate_mei_strand": "",
                    "mei_hit": bool(mei_hit),
                    "mate_mei_hit": False,
                    "short_mei_seed_rescued": _as_bool(getattr(rec, "short_mei_seed_rescued", False)),
                    "polya_rescue": bool(poly_rescued),
                    "mei_hit_source": hit_source,
                    "clip_len": int(clip_len),
                    "soft_clip_len": int(clip_len),
                    "soft_clip_side": str(getattr(rec, "clip_side", "") or ""),
                    "soft_clip_pos": int(rec.pos),
                    "poly_tail_rescued": bool(poly_rescued),
                    "clip_poly_at_run": int(poly_run),
                    "polya_base": poly_base,
                    "clip_poly_base": poly_base,
                }
            )

    disc = discordant_hits if not discordant_hits.empty else discordant_mate_hits
    if not disc.empty:
        mid_cache: dict[tuple[str, int, int], int] = {}
        for rec in disc.itertuples(index=False):
            mei_hit = bool(getattr(rec, "mei_hit", False))
            mate_mei_hit = bool(getattr(rec, "mate_mei_hit", False))
            vntr_rescue = bool(getattr(rec, "vntr_rescue", False))
            polya_rescue = bool(getattr(rec, "polya_rescue", False))
            # Keep detail focused on MEI-supporting discordants (consensus or rescue).
            if not (mei_hit or mate_mei_hit or vntr_rescue or polya_rescue):
                continue
            key = (str(rec.chrom), int(rec.window_start), int(rec.window_end))
            if key not in mid_cache:
                mid_cache[key] = (int(rec.window_start) + int(rec.window_end)) // 2
            anchor_side = "L" if int(rec.pos) <= mid_cache[key] else "R"
            soft_clip_pos = int(getattr(rec, "soft_clip_pos", 0) or 0)
            ref_end = int(getattr(rec, "ref_end", 0) or 0)
            if ref_end <= 0:
                read_seq = str(getattr(rec, "read_seq", "") or "")
                ref_end = int(rec.pos) + max(0, len(read_seq) - 1) if read_seq else int(rec.pos)
            if soft_clip_pos > 0:
                genomic_pos = soft_clip_pos
            elif anchor_side == "L":
                genomic_pos = ref_end
            else:
                genomic_pos = int(rec.pos)
            hit_source = str(getattr(rec, "mei_hit_source", "") or "")
            if not hit_source:
                if polya_rescue:
                    hit_source = "polya_rescue"
                elif vntr_rescue:
                    hit_source = "vntr_rescue"
            mate_seq = str(getattr(rec, "mate_seq", "") or "")
            polya_base = ""
            if polya_rescue and mate_seq:
                n_a = mate_seq.upper().count("A")
                n_t = mate_seq.upper().count("T")
                polya_base = "T" if n_t > n_a else ("A" if n_a else "")
            rows.append(
                {
                    "sample": sample,
                    "evidence_type": "DPE",
                    "read_name": str(rec.read_name),
                    "chrom": str(rec.chrom),
                    "window_start": int(rec.window_start),
                    "window_end": int(rec.window_end),
                    "anchor_side": anchor_side,
                    "genomic_pos": int(genomic_pos),
                    "mate_chrom": str(getattr(rec, "mate_chrom", "")),
                    "mate_genomic_pos": int(getattr(rec, "mate_pos", 0) or 0),
                    "mei_start": int(getattr(rec, "target_start", 0) or 0),
                    "mei_end": int(getattr(rec, "target_end", 0) or 0),
                    "mate_mei_start": int(getattr(rec, "mate_mei_start", 0) or 0),
                    "mate_mei_end": int(getattr(rec, "mate_mei_end", 0) or 0),
                    "mei_target": str(getattr(rec, "target", "") or ""),
                    "mate_mei_target": str(getattr(rec, "mate_mei_target", "") or ""),
                    "mei_strand": (
                        str(getattr(rec, "target_strand", "") or "")
                        if str(getattr(rec, "target_strand", "") or "") in {"+", "-"}
                        else ""
                    ),
                    "mate_mei_strand": (
                        str(getattr(rec, "mate_mei_strand", "") or "")
                        if str(getattr(rec, "mate_mei_strand", "") or "") in {"+", "-"}
                        else ""
                    ),
                    "mei_hit": bool(mei_hit or vntr_rescue or polya_rescue),
                    "mate_mei_hit": bool(mate_mei_hit or vntr_rescue or polya_rescue),
                    "vntr_rescue": bool(vntr_rescue),
                    "polya_rescue": bool(polya_rescue),
                    "mei_hit_source": hit_source,
                    "mate_seq_len": int(len(mate_seq)),
                    "polya_base": polya_base,
                    "soft_clip_side": str(getattr(rec, "soft_clip_side", "") or ""),
                    "soft_clip_len": int(getattr(rec, "soft_clip_len", 0) or 0),
                    "soft_clip_pos": int(getattr(rec, "soft_clip_pos", 0) or 0),
                    "anchor_poly_at_run": int(getattr(rec, "anchor_poly_at_run", 0) or 0),
                    "anchor_poly_base": str(getattr(rec, "anchor_poly_base", "") or ""),
                    "anchor_poly_side": str(getattr(rec, "anchor_poly_side", "") or ""),
                    "poly_tail_anchor_rescued": bool(getattr(rec, "poly_tail_anchor_rescued", False)),
                }
            )

    _detail_empty_cols = [
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
        "mei_target",
        "mate_mei_target",
        "mei_strand",
        "mate_mei_strand",
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
    ]
    if not rows:
        return pd.DataFrame(columns=_detail_empty_cols)
    out = pd.DataFrame(rows)
    for col, default in (
        ("mei_strand", ""),
        ("mate_mei_strand", ""),
        ("short_mei_seed_rescued", False),
        ("vntr_rescue", False),
        ("polya_rescue", False),
        ("poly_tail_rescued", False),
        ("clip_poly_at_run", 0),
        ("clip_len", 0),
        ("soft_clip_len", 0),
    ):
        if col not in out.columns:
            out[col] = default
    out["mei_strand"] = out["mei_strand"].fillna("").astype(str)
    out["mate_mei_strand"] = out["mate_mei_strand"].fillna("").astype(str)
    out["short_mei_seed_rescued"] = out["short_mei_seed_rescued"].fillna(False).infer_objects(copy=False).astype(bool)
    out["polya_rescue"] = out["polya_rescue"].fillna(False).astype(bool)
    out["poly_tail_rescued"] = out["poly_tail_rescued"].fillna(False).infer_objects(copy=False).astype(bool)
    return out



def _project_detail_coords_to_full(
    detail: pd.DataFrame,
    frag_map: dict[str, FragmentToFullMap],
) -> pd.DataFrame:
    """Project panel MEI coords in supporting-read detail onto the full-length axis.

    Uses the one-time fragment→full map from public-data download. PolyA/VNTR
    rescue rows are left unchanged (synthetic / non-consensus).
    """
    if detail is None or detail.empty or not frag_map:
        return detail
    out = detail.copy()
    if "mei_target" not in out.columns:
        out["mei_target"] = ""
    if "mate_mei_target" not in out.columns:
        out["mate_mei_target"] = ""

    rescue = pd.Series(False, index=out.index)
    if "polya_rescue" in out.columns:
        rescue = rescue | out["polya_rescue"].fillna(False).astype(bool)
    if "vntr_rescue" in out.columns:
        rescue = rescue | out["vntr_rescue"].fillna(False).astype(bool)
    if "mei_hit_source" in out.columns:
        src = out["mei_hit_source"].fillna("").astype(str)
        rescue = rescue | src.isin({"polya_rescue", "vntr_rescue"})

    for start_col, end_col, target_col in (
        ("mei_start", "mei_end", "mei_target"),
        ("mate_mei_start", "mate_mei_end", "mate_mei_target"),
    ):
        if start_col not in out.columns or end_col not in out.columns:
            continue
        starts = pd.to_numeric(out[start_col], errors="coerce").fillna(0).astype(int)
        ends = pd.to_numeric(out[end_col], errors="coerce").fillna(0).astype(int)
        targets = out[target_col].fillna("").astype(str)
        new_starts = starts.copy()
        new_ends = ends.copy()
        for idx in out.index:
            if bool(rescue.loc[idx]):
                continue
            if int(starts.loc[idx]) <= 0 and int(ends.loc[idx]) <= 0:
                continue
            projected = _project_panel_coords_to_full(
                int(starts.loc[idx]),
                int(ends.loc[idx]),
                str(targets.loc[idx]),
                frag_map,
            )
            if projected is None:
                continue
            new_starts.loc[idx] = int(projected[0])
            new_ends.loc[idx] = int(projected[1])
            if target_col in out.columns and projected[2]:
                out.loc[idx, target_col] = projected[2]
        out[start_col] = new_starts
        out[end_col] = new_ends
    return out


def _overlay_full_consensus_coords_onto_detail(
    detail: pd.DataFrame,
    full_detail: pd.DataFrame,
) -> pd.DataFrame:
    """Deprecated path: kept for callers; prefer ``_project_detail_coords_to_full``."""
    if detail is None or detail.empty or full_detail is None or full_detail.empty:
        return detail
    keys = [
        c
        for c in ("sample", "evidence_type", "read_name", "chrom", "window_start", "window_end")
        if c in detail.columns and c in full_detail.columns
    ]
    if len(keys) < 6:
        return detail
    keep_cols = keys + [
        c
        for c in (
            "mei_start",
            "mei_end",
            "mate_mei_start",
            "mate_mei_end",
            "mei_hit",
            "mate_mei_hit",
        )
        if c in full_detail.columns
    ]
    full = full_detail.loc[:, keep_cols].copy()
    full = full.rename(
        columns={
            "mei_start": "mei_start_full",
            "mei_end": "mei_end_full",
            "mate_mei_start": "mate_mei_start_full",
            "mate_mei_end": "mate_mei_end_full",
            "mei_hit": "mei_hit_full",
            "mate_mei_hit": "mate_mei_hit_full",
        }
    )
    # One full row per read key (prefer rows with any MEI hit).
    if "mei_hit_full" in full.columns or "mate_mei_hit_full" in full.columns:
        hit = pd.Series(False, index=full.index)
        if "mei_hit_full" in full.columns:
            hit = hit | full["mei_hit_full"].fillna(False).astype(bool)
        if "mate_mei_hit_full" in full.columns:
            hit = hit | full["mate_mei_hit_full"].fillna(False).astype(bool)
        full["_hit_rank"] = hit.astype(int)
        full = full.sort_values("_hit_rank", ascending=False).drop_duplicates(keys, keep="first")
        full = full.drop(columns=["_hit_rank"])
    else:
        full = full.drop_duplicates(keys, keep="first")

    out = detail.merge(full, on=keys, how="left")
    rescue = pd.Series(False, index=out.index)
    if "polya_rescue" in out.columns:
        rescue = rescue | out["polya_rescue"].fillna(False).astype(bool)
    if "vntr_rescue" in out.columns:
        rescue = rescue | out["vntr_rescue"].fillna(False).astype(bool)
    if "mei_hit_source" in out.columns:
        src = out["mei_hit_source"].fillna("").astype(str)
        rescue = rescue | src.isin({"polya_rescue", "vntr_rescue"})

    for src_col, dst_col in (
        ("mei_start_full", "mei_start"),
        ("mei_end_full", "mei_end"),
        ("mate_mei_start_full", "mate_mei_start"),
        ("mate_mei_end_full", "mate_mei_end"),
    ):
        if src_col not in out.columns or dst_col not in out.columns:
            continue
        src = pd.to_numeric(out[src_col], errors="coerce").fillna(0).astype(int)
        use = (~rescue) & src.gt(0)
        if use.any():
            out.loc[use, dst_col] = src.loc[use]
    drop_tmp = [c for c in out.columns if c.endswith("_full") and c not in detail.columns]
    if drop_tmp:
        out = out.drop(columns=drop_tmp)
    return out


def _robust_coord_extent(lo_values: pd.Series, hi_values: pd.Series) -> tuple[float, float]:
    """Return outlier-resistant min/max MEI coords for one locus/sample group.

    Uses Tukey fences (k=3) on the pooled start/end endpoints when enough
    points exist; otherwise falls back to raw min/max. This keeps true full-
    length SVA/LINE1 footprints while dropping rare off-target mates that can
    inflate an Alu-sized insertion to >1 kb.
    """
    pts = pd.concat(
        [
            pd.to_numeric(lo_values, errors="coerce"),
            pd.to_numeric(hi_values, errors="coerce"),
        ],
        ignore_index=True,
    )
    pts = pts[pts.gt(0)].astype(float)
    if pts.empty:
        return float("nan"), float("nan")
    if len(pts) < 8:
        return float(pts.min()), float(pts.max())
    q1 = float(pts.quantile(0.25))
    q3 = float(pts.quantile(0.75))
    iqr = max(q3 - q1, 1.0)
    lo_fence = q1 - 3.0 * iqr
    hi_fence = q3 + 3.0 * iqr
    kept = pts[(pts >= lo_fence) & (pts <= hi_fence)]
    if kept.empty:
        return float(pts.min()), float(pts.max())
    return float(kept.min()), float(kept.max())


def _candidate_mei_target_lengths(candidates: pd.DataFrame) -> pd.DataFrame:
    """Per-locus consensus MEI target length used to keep on-target mappings only.

    Prefer assembly consensus target length when present. Side-level
    ``*_mei_target_len`` values can come from off-family hits (e.g. LINE1 3294
    on an Alu locus) and must not override the assembled element.
    """
    key_cols = ["chrom", "window_start", "window_end"]
    if candidates is None or candidates.empty or not set(key_cols).issubset(candidates.columns):
        return pd.DataFrame(columns=key_cols + ["mei_target_length"])
    fallback_cols = [
        "disease_L_mei_target_len",
        "disease_R_mei_target_len",
        "control_L_mei_target_len",
        "control_R_mei_target_len",
        "disease_full_L_mei_target_len",
        "disease_full_R_mei_target_len",
        "control_full_L_mei_target_len",
        "control_full_R_mei_target_len",
    ]
    work = candidates.loc[:, key_cols].copy()
    asm = (
        pd.to_numeric(candidates["asm_mei_target_length"], errors="coerce")
        if "asm_mei_target_length" in candidates.columns
        else pd.Series(float("nan"), index=candidates.index)
    )
    present = [c for c in fallback_cols if c in candidates.columns]
    if present:
        fallback = pd.concat(
            [pd.to_numeric(candidates[c], errors="coerce") for c in present],
            axis=1,
        )
        fallback_len = fallback.where(fallback.gt(0)).max(axis=1, skipna=True)
    else:
        fallback_len = pd.Series(float("nan"), index=candidates.index)
    work["mei_target_length"] = asm.where(asm.gt(0), fallback_len)
    return (
        work.groupby(key_cols, as_index=False)["mei_target_length"]
        .max()
        .loc[:, key_cols + ["mei_target_length"]]
    )


def _keep_on_target_mei_interval(
    start: pd.Series,
    end: pd.Series,
    target_length: pd.Series,
    *,
    slack: int = 50,
) -> tuple[pd.Series, pd.Series]:
    """Zero intervals that fall outside the consensus target element.

    When ``target_length`` is known, only keep mappings with both endpoints in
    ``[1, target_length + slack]``. Slack covers minor end overhangs (e.g. SVA
    1378 vs consensus 1375) without admitting wrong-element mates (e.g. Alu
    mate at 1330-1462 against a ~312 bp target).
    """
    start_n = pd.to_numeric(start, errors="coerce").fillna(0).astype(int)
    end_n = pd.to_numeric(end, errors="coerce").fillna(0).astype(int)
    tlen = pd.to_numeric(target_length, errors="coerce")
    has_tlen = tlen.gt(0).fillna(False)
    max_pos = (tlen + float(slack)).where(has_tlen, float("inf"))
    on_target = start_n.gt(0) & end_n.ge(start_n) & (~has_tlen | (start_n.le(max_pos) & end_n.le(max_pos)))
    return start_n.where(on_target, 0), end_n.where(on_target, 0)


def _on_target_extent_ok(
    lo: pd.Series,
    hi: pd.Series,
    target_length: pd.Series,
    *,
    slack: int = 50,
) -> pd.Series:
    """True when both extent endpoints map within the consensus target element."""
    lo_n = pd.to_numeric(lo, errors="coerce")
    hi_n = pd.to_numeric(hi, errors="coerce")
    tlen = pd.to_numeric(target_length, errors="coerce")
    has_tlen = tlen.gt(0).fillna(False)
    max_pos = (tlen + float(slack)).where(has_tlen, float("inf"))
    return lo_n.gt(0) & hi_n.ge(lo_n) & (~has_tlen | (lo_n.le(max_pos) & hi_n.le(max_pos)))


def _aggregate_detail_mei_extents(
    detail: pd.DataFrame,
    target_lengths: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Per-locus min/max MEI consensus coords from supporting-read detail rows.

    Matches the plot-path footprint: robust min/max over SR ``mei_start/end``
    and DPE ``mate_mei_start/end`` (and any direct ``mei_*`` hits) for disease
    and control samples separately, plus a combined locus extent.

    Only intervals that map within the consensus target element length are
    included when that length is known.

    Also emits per-side SR extents (``{sample}_{L|R}_detail_mei_start/end``) so
    gold/annotation rebuilds can restore zeroed L/R aggregated coords, and
    per-sample mapped-read counts for sample selection.
    """
    key_cols = ["chrom", "window_start", "window_end"]
    empty_cols = key_cols + [
        "disease_detail_mei_start_min",
        "disease_detail_mei_end_max",
        "control_detail_mei_start_min",
        "control_detail_mei_end_max",
        "detail_mei_start_min",
        "detail_mei_end_max",
        "disease_detail_mei_mapped_reads",
        "control_detail_mei_mapped_reads",
        "disease_L_detail_mei_start",
        "disease_L_detail_mei_end",
        "disease_R_detail_mei_start",
        "disease_R_detail_mei_end",
        "control_L_detail_mei_start",
        "control_L_detail_mei_end",
        "control_R_detail_mei_start",
        "control_R_detail_mei_end",
    ]
    if detail is None or detail.empty:
        return pd.DataFrame(columns=empty_cols)
    required = set(key_cols + ["sample", "mei_start", "mei_end", "mate_mei_start", "mate_mei_end"])
    if not required.issubset(set(detail.columns)):
        return pd.DataFrame(columns=empty_cols)

    work = detail.loc[
        :, list(required | {"mei_hit", "mate_mei_hit", "evidence_type", "anchor_side", "read_name"})
    ].copy()
    work["sample"] = work["sample"].fillna("").astype(str).str.lower()
    work = work.loc[work["sample"].isin(["disease", "control"])].copy()
    if work.empty:
        return pd.DataFrame(columns=empty_cols)

    if target_lengths is not None and not target_lengths.empty and "mei_target_length" in target_lengths.columns:
        work = work.merge(
            target_lengths.loc[:, key_cols + ["mei_target_length"]],
            on=key_cols,
            how="left",
        )
    else:
        work["mei_target_length"] = float("nan")

    mei_hit = (
        work["mei_hit"].fillna(False).astype(bool)
        if "mei_hit" in work.columns
        else pd.Series(True, index=work.index)
    )
    mate_hit = (
        work["mate_mei_hit"].fillna(False).astype(bool)
        if "mate_mei_hit" in work.columns
        else pd.Series(True, index=work.index)
    )
    mei_start = pd.to_numeric(work["mei_start"], errors="coerce").fillna(0).astype(int).where(mei_hit, 0)
    mei_end = pd.to_numeric(work["mei_end"], errors="coerce").fillna(0).astype(int).where(mei_hit, 0)
    mate_start = (
        pd.to_numeric(work["mate_mei_start"], errors="coerce").fillna(0).astype(int).where(mate_hit, 0)
    )
    mate_end = pd.to_numeric(work["mate_mei_end"], errors="coerce").fillna(0).astype(int).where(mate_hit, 0)
    mei_start, mei_end = _keep_on_target_mei_interval(mei_start, mei_end, work["mei_target_length"])
    mate_start, mate_end = _keep_on_target_mei_interval(mate_start, mate_end, work["mei_target_length"])
    work["extent_lo"] = pd.concat(
        [mei_start.where(mei_start.gt(0)), mate_start.where(mate_start.gt(0))],
        axis=1,
    ).min(axis=1, skipna=True)
    work["extent_hi"] = pd.concat(
        [mei_end.where(mei_end.gt(0)), mate_end.where(mate_end.gt(0))],
        axis=1,
    ).max(axis=1, skipna=True)
    work = work.loc[work["extent_lo"].gt(0) & work["extent_hi"].ge(work["extent_lo"])].copy()
    if work.empty:
        return pd.DataFrame(columns=empty_cols)

    rows: list[dict[str, object]] = []
    for (chrom, ws, we, sample), grp in work.groupby(key_cols + ["sample"], sort=False):
        lo, hi = _robust_coord_extent(grp["extent_lo"], grp["extent_hi"])
        n_reads = (
            grp["read_name"].fillna("").astype(str).nunique()
            if "read_name" in grp.columns
            else int(len(grp))
        )
        rows.append(
            {
                "chrom": chrom,
                "window_start": int(ws),
                "window_end": int(we),
                "sample": sample,
                "extent_lo": lo,
                "extent_hi": hi,
                "mapped_reads": int(n_reads),
            }
        )
    per_sample = pd.DataFrame(rows)
    disease = per_sample.loc[per_sample["sample"].eq("disease"), key_cols + ["extent_lo", "extent_hi", "mapped_reads"]].rename(
        columns={
            "extent_lo": "disease_detail_mei_start_min",
            "extent_hi": "disease_detail_mei_end_max",
            "mapped_reads": "disease_detail_mei_mapped_reads",
        }
    )
    control = per_sample.loc[per_sample["sample"].eq("control"), key_cols + ["extent_lo", "extent_hi", "mapped_reads"]].rename(
        columns={
            "extent_lo": "control_detail_mei_start_min",
            "extent_hi": "control_detail_mei_end_max",
            "mapped_reads": "control_detail_mei_mapped_reads",
        }
    )
    combined_rows: list[dict[str, object]] = []
    for (chrom, ws, we), grp in work.groupby(key_cols, sort=False):
        lo, hi = _robust_coord_extent(grp["extent_lo"], grp["extent_hi"])
        combined_rows.append(
            {
                "chrom": chrom,
                "window_start": int(ws),
                "window_end": int(we),
                "detail_mei_start_min": lo,
                "detail_mei_end_max": hi,
            }
        )
    combined = pd.DataFrame(combined_rows)
    out = combined.merge(disease, on=key_cols, how="left").merge(control, on=key_cols, how="left")

    # Per-side SR extents (used to restore zeroed disease/control_L/R_mei_start/end).
    if "evidence_type" in work.columns and "anchor_side" in work.columns:
        sr = work.loc[
            work["evidence_type"].astype(str).str.upper().eq("SR")
            & work["anchor_side"].astype(str).str.upper().isin(["L", "R"])
        ].copy()
        if not sr.empty:
            sr["anchor_side"] = sr["anchor_side"].astype(str).str.upper().str[:1]
            side_rows: list[dict[str, object]] = []
            for (chrom, ws, we, sample, side), grp in sr.groupby(
                key_cols + ["sample", "anchor_side"], sort=False
            ):
                lo, hi = _robust_coord_extent(grp["extent_lo"], grp["extent_hi"])
                side_rows.append(
                    {
                        "chrom": chrom,
                        "window_start": int(ws),
                        "window_end": int(we),
                        "sample": sample,
                        "anchor_side": side,
                        "extent_lo": lo,
                        "extent_hi": hi,
                    }
                )
            side_agg = pd.DataFrame(side_rows)
            for sample in ("disease", "control"):
                for side in ("L", "R"):
                    part = side_agg.loc[
                        side_agg["sample"].eq(sample) & side_agg["anchor_side"].eq(side),
                        key_cols + ["extent_lo", "extent_hi"],
                    ].rename(
                        columns={
                            "extent_lo": f"{sample}_{side}_detail_mei_start",
                            "extent_hi": f"{sample}_{side}_detail_mei_end",
                        }
                    )
                    out = out.merge(part, on=key_cols, how="left")

    for col in empty_cols:
        if col not in out.columns:
            out[col] = pd.NA
    return out.loc[:, empty_cols]


def _merge_detail_mei_extents(candidates: pd.DataFrame, detail: pd.DataFrame | None) -> pd.DataFrame:
    """Attach detail-derived MEI extents onto candidate rows when available.

    When aggregated ``{sample}_{L|R}_mei_start/end`` are zero/missing but detail
    SR extents exist, restore those L/R fields from detail (same min/max the
    plots use).
    """
    if detail is None or detail.empty or candidates.empty:
        return candidates
    extents = _aggregate_detail_mei_extents(
        detail,
        target_lengths=_candidate_mei_target_lengths(candidates),
    )
    if extents.empty:
        return candidates
    out = candidates.copy()
    drop_cols = [c for c in extents.columns if c not in {"chrom", "window_start", "window_end"} and c in out.columns]
    if drop_cols:
        out = out.drop(columns=drop_cols)
    out = out.merge(extents, on=["chrom", "window_start", "window_end"], how="left")

    for sample in ("disease", "control"):
        for side in ("L", "R"):
            start_col = f"{sample}_{side}_mei_start"
            end_col = f"{sample}_{side}_mei_end"
            d_start = f"{sample}_{side}_detail_mei_start"
            d_end = f"{sample}_{side}_detail_mei_end"
            if d_start not in out.columns or d_end not in out.columns:
                continue
            if start_col not in out.columns:
                out[start_col] = 0
            if end_col not in out.columns:
                out[end_col] = 0
            cur_start = pd.to_numeric(out[start_col], errors="coerce").fillna(0)
            cur_end = pd.to_numeric(out[end_col], errors="coerce").fillna(0)
            new_start = pd.to_numeric(out[d_start], errors="coerce")
            new_end = pd.to_numeric(out[d_end], errors="coerce")
            replace = cur_start.le(0) & new_start.gt(0) & new_end.ge(new_start)
            if replace.any():
                out.loc[replace, start_col] = new_start.loc[replace]
                out.loc[replace, end_col] = new_end.loc[replace]
    return out


def _assign_rows_to_candidate_loci(split_df: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    if split_df.empty or candidates.empty:
        return pd.DataFrame(columns=list(split_df.columns) + ["window_start", "window_end"])

    trees: dict[str, IntervalTree] = {}
    loci = candidates.loc[:, ["chrom", "window_start", "window_end"]].drop_duplicates()
    for row in loci.itertuples(index=False):
        chrom = str(row.chrom)
        tree = trees.setdefault(chrom, IntervalTree())
        tree.addi(int(row.window_start), int(row.window_end) + 1, (int(row.window_start), int(row.window_end)))

    assigned_rows: list[dict[str, object]] = []
    for row in split_df.itertuples(index=False):
        chrom = str(row.chrom)
        pos = int(row.pos)
        tree = trees.get(chrom)
        if tree is None:
            continue
        overlaps = list(tree.at(pos))
        if not overlaps:
            continue
        best = min(overlaps, key=lambda iv: (iv.end - iv.begin, abs(((iv.begin + iv.end) // 2) - pos)))
        locus_start, locus_end = best.data
        as_dict = row._asdict()
        as_dict["window_start"] = locus_start
        as_dict["window_end"] = locus_end
        assigned_rows.append(as_dict)
    return pd.DataFrame(assigned_rows)


def _poly_at_artifact_tsd_mask(tsd_seq: pd.Series) -> pd.Series:
    """True for sequences that are polyA/polyT tails, not real TSDs."""
    tsd_seq_s = tsd_seq.fillna("").astype(str).str.upper()
    tsd_len_s = tsd_seq_s.str.len().astype(int)
    a_fraction = (tsd_seq_s.str.count("A") / tsd_len_s.replace(0, pd.NA)).fillna(0.0).infer_objects(copy=False).astype(float)
    t_fraction = (tsd_seq_s.str.count("T") / tsd_len_s.replace(0, pd.NA)).fillna(0.0).infer_objects(copy=False).astype(float)
    dominant_poly_fraction = pd.concat([a_fraction, t_fraction], axis=1).max(axis=1)
    longest_at_run = tsd_seq_s.str.findall(r"[AT]+").map(
        lambda parts: max((len(p) for p in parts), default=0)
    ).astype(int)
    longest_a_run = tsd_seq_s.str.findall(r"A+").map(
        lambda parts: max((len(p) for p in parts), default=0)
    ).astype(int)
    longest_t_run = tsd_seq_s.str.findall(r"T+").map(
        lambda parts: max((len(p) for p in parts), default=0)
    ).astype(int)
    poly_at_only = tsd_len_s.ge(4) & tsd_seq_s.str.fullmatch(r"[AT]+", na=False)
    near_poly_at = tsd_len_s.ge(8) & dominant_poly_fraction.ge(0.85) & longest_at_run.ge(6)
    # Long A/T homopolymer (≥12) is a poly-tail fragment even with a short GC tip.
    long_homopolymer = (longest_a_run.ge(12) | longest_t_run.ge(12)) & tsd_len_s.ge(12)
    return (poly_at_only | near_poly_at | long_homopolymer).fillna(False).astype(bool)


def _clear_poly_at_artifact_tsd_fields(
    out: pd.DataFrame,
    *,
    seq_col: str = "tsd_seq",
    len_col: str = "tsd_len_estimate",
    detected_col: str = "tsd_detected",
    source_col: str = "tsd_evidence_source",
    filter_flag_col: str = "tsd_poly_at_filter_applied",
) -> pd.Series:
    """Blank polyA/T artifact TSDs (all evidence sources). Returns filter mask."""
    if out is None or out.empty or seq_col not in out.columns:
        empty = pd.Series(False, index=(out.index if out is not None else None))
        if out is not None and filter_flag_col not in out.columns:
            out[filter_flag_col] = False
        return empty
    mask = _poly_at_artifact_tsd_mask(out[seq_col])
    prev = (
        out[filter_flag_col].fillna(False).astype(bool)
        if filter_flag_col in out.columns
        else pd.Series(False, index=out.index)
    )
    out[filter_flag_col] = (prev | mask).astype(bool)
    if not bool(mask.any()):
        if detected_col in out.columns and len_col in out.columns:
            out[detected_col] = pd.to_numeric(out[len_col], errors="coerce").fillna(0).astype(int) >= 4
        return mask
    if "tsd_left_breakpoint" in out.columns:
        out.loc[mask, "tsd_left_breakpoint"] = 0
    if "tsd_right_breakpoint" in out.columns:
        out.loc[mask, "tsd_right_breakpoint"] = 0
    if len_col in out.columns:
        out.loc[mask, len_col] = 0
    out.loc[mask, seq_col] = ""
    if detected_col in out.columns:
        if len_col in out.columns:
            out[detected_col] = pd.to_numeric(out[len_col], errors="coerce").fillna(0).astype(int) >= 4
        else:
            out.loc[mask, detected_col] = False
    if source_col in out.columns:
        src = out.loc[mask, source_col].fillna("").astype(str)
        already = src.str.contains("filtered_poly_at", regex=False)
        out.loc[mask & ~already, source_col] = src.loc[mask & ~already] + "_filtered_poly_at_only"
    return mask


def _poly_at_stats(seq: str) -> tuple[int, float, str]:
    """PolyA/T stats: longest dominant-base run, dominant-base fraction, base.

    Purity is ``max(n_A, n_T) / length`` — mostly A **or** mostly T — not
    combined A+T. Mixed AT sequence scores ~0.5 and fails typical thresholds.
    """
    s = (seq or "").upper()
    if not s:
        return (0, 0.0, "")
    n_a = s.count("A")
    n_t = s.count("T")
    if n_a <= 0 and n_t <= 0:
        return (0, 0.0, "")
    if n_a >= n_t:
        base = "A"
        n_dom = n_a
    else:
        base = "T"
        n_dom = n_t
    frac = float(n_dom) / float(len(s))
    best = 0
    cur = 0
    for ch in s:
        if ch == base:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return (int(best), float(frac), base)


def _clip_to_poly_at_region(seq: str, *, min_dom_frac: float | None = None) -> str:
    """Return the longest mostly-A or mostly-T substring (empty if none)."""
    thr = float(_POLYA_MIN_FRAC if min_dom_frac is None else min_dom_frac)
    _length, _frac, _base, span = _longest_poly_at_span(seq, min_frac=thr)
    return span


def _normalized_clip_seq(seq: str, *, max_len: int = 80) -> str:
    s = "".join(ch for ch in (seq or "").upper() if ch in {"A", "C", "G", "T"})
    if not s:
        return ""
    return s[: int(max_len)]


def _clip_shannon_entropy(seq: str) -> float:
    s = _normalized_clip_seq(seq)
    if not s:
        return 0.0
    n = float(len(s))
    counts = {b: s.count(b) for b in ("A", "C", "G", "T")}
    ent = 0.0
    for v in counts.values():
        if v <= 0:
            continue
        p = float(v) / n
        ent -= p * math.log2(p)
    return float(ent)


def _breakpoint_proximal_clip_seq(seq: str, side: str, *, max_len: int = 40) -> str:
    s = _normalized_clip_seq(seq, max_len=200)
    if not s:
        return ""
    side_u = (side or "").upper()
    if side_u == "L":
        # L-clips are left-anchored; breakpoint-proximal bases are near clip end.
        return s[-int(max_len) :]
    # R-clips are right-anchored; breakpoint-proximal bases are near clip start.
    return s[: int(max_len)]


def _is_informative_split_clip(seq: str, *, min_len: int = 20, min_non_at_fraction: float = 0.15, min_entropy: float = 1.20) -> bool:
    s = _normalized_clip_seq(seq)
    if len(s) < int(min_len):
        return False
    non_at = sum(1 for ch in s if ch in {"C", "G"})
    non_at_fraction = float(non_at) / float(len(s))
    if non_at_fraction < float(min_non_at_fraction):
        return False
    return _clip_shannon_entropy(s) >= float(min_entropy)


def _pair_clip_similarity(a: str, b: str) -> float:
    """Breakpoint-proximal basewise identity.

    Compare only the overlapping breakpoint-proximal span. This reflects
    "do reads agree at the same reference-relative offset?" rather than
    generic global sequence similarity.
    """
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    if n <= 0:
        return 0.0
    aa = a[:n]
    bb = b[:n]
    matches = sum(1 for x, y in zip(aa, bb) if x == y)
    return float(matches) / float(n)


def _clip_overlap_consistency_stats(group: pd.DataFrame) -> tuple[int, float, float]:
    if group.empty:
        return (0, 0.0, 0.0)
    work = group.copy()
    # Keep one representative clip per read (highest MEI score first) so one
    # read with multiple alignments does not dominate overlap statistics.
    if "read_name" in work.columns:
        work["mei_score_effective"] = pd.to_numeric(work.get("mei_score_effective", 0.0), errors="coerce").fillna(0.0)
        work = work.sort_values("mei_score_effective", ascending=False).drop_duplicates("read_name", keep="first")
    work = work.head(40).copy()
    side_hint = ""
    if "clip_side" in work.columns and not work.empty:
        side_hint = str(work["clip_side"].iloc[0] or "")
    informative_rows = []
    for row in work.itertuples(index=False):
        seq = str(getattr(row, "clip_seq", "") or "")
        side = str(getattr(row, "clip_side", side_hint) or side_hint)
        if not _is_informative_split_clip(seq):
            continue
        prox = _breakpoint_proximal_clip_seq(seq, side)
        if len(prox) < 12:
            continue
        pos_val = int(getattr(row, "pos", 0) or 0) if hasattr(row, "pos") else 0
        informative_rows.append((pos_val, prox))
    informative_reads = int(len(informative_rows))
    if informative_reads < 2:
        return (informative_reads, 0.0, 0.0)
    pos_counts: dict[int, int] = defaultdict(int)
    for p, _ in informative_rows:
        pos_counts[int(p)] += 1
    mode_support = max(pos_counts.values()) if pos_counts else 0
    mode_fraction = float(mode_support) / float(informative_reads) if informative_reads > 0 else 0.0
    sims: list[float] = []
    ge20 = 0
    total_pairs = 0
    for i in range(informative_reads):
        _, a = informative_rows[i]
        for j in range(i + 1, informative_reads):
            _, b = informative_rows[j]
            raw_identity = _pair_clip_similarity(a, b)
            # Random DNA baseline is 0.25; scale to [0,1] above chance.
            sim = max(0.0, min(1.0, (raw_identity - 0.25) / 0.75))
            sims.append(float(sim))
            total_pairs += 1
            if sim >= 0.20:
                ge20 += 1
    if not sims or total_pairs <= 0:
        return (informative_reads, 0.0, 0.0)
    sims_sorted = sorted(sims)
    mid = len(sims_sorted) // 2
    if len(sims_sorted) % 2 == 1:
        median_sim = float(sims_sorted[mid])
    else:
        median_sim = float((sims_sorted[mid - 1] + sims_sorted[mid]) / 2.0)
    frac_ge20 = float(ge20) / float(total_pairs)
    # Penalize diffuse breakpoint-position clouds: real events usually have
    # consistent split positions, noisy loci often do not.
    median_sim = float(median_sim) * float(mode_fraction)
    return (informative_reads, float(median_sim), float(frac_ge20))


def _collect_indel_breakpoint_evidence(
    bam_path: Path,
    candidates: pd.DataFrame,
    *,
    sample: str,
    min_mapq: int = 20,
    min_indel_bp: int = 12,
    query_context_bases: int = 12,
) -> pd.DataFrame:
    """Collect breakpoint-proximal CIGAR indel evidence assigned to candidate loci."""
    key_cols = ["chrom", "window_start", "window_end"]
    if candidates.empty:
        return pd.DataFrame(
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
                "evidence_type",
                "indel_type",
                "indel_len",
            ]
        )

    loci = candidates.loc[:, key_cols].drop_duplicates()
    trees: dict[str, IntervalTree] = {}
    span_by_chrom: dict[str, tuple[int, int]] = {}
    for row in loci.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.window_start)
        end = int(row.window_end)
        tree = trees.setdefault(chrom, IntervalTree())
        tree.addi(start, end + 1, (start, end))
        if chrom not in span_by_chrom:
            span_by_chrom[chrom] = (start, end)
        else:
            lo, hi = span_by_chrom[chrom]
            span_by_chrom[chrom] = (min(lo, start), max(hi, end))

    rows: list[dict[str, object]] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chrom, tree in trees.items():
            lo, hi = span_by_chrom[chrom]
            fetch_start0 = max(0, int(lo) - 1)
            fetch_end0 = max(fetch_start0 + 1, int(hi))
            for read in bam.fetch(chrom, fetch_start0, fetch_end0):
                if read.is_unmapped:
                    continue
                if read.is_qcfail or read.is_duplicate or read.is_secondary:
                    continue
                if read.mapping_quality < int(min_mapq):
                    continue
                if not read.cigartuples:
                    continue
                query_seq = str(read.query_sequence or "")
                if not query_seq:
                    continue

                ref_pos = int(read.reference_start) + 1  # 1-based
                query_pos = 0
                for op, length in read.cigartuples:
                    l = int(length)
                    if op in {0, 7, 8}:  # M/=/X
                        ref_pos += l
                        query_pos += l
                        continue
                    if op == 1:  # insertion relative to reference
                        if l >= int(min_indel_bp):
                            pos = max(1, int(ref_pos))
                            overlaps = list(tree.at(pos))
                            if overlaps:
                                best = min(
                                    overlaps,
                                    key=lambda iv: (iv.end - iv.begin, abs(((iv.begin + iv.end) // 2) - pos)),
                                )
                                window_start, window_end = best.data
                                ins_seq = query_seq[query_pos : query_pos + l]
                                if len(ins_seq) < 8:
                                    q0 = max(0, query_pos - int(query_context_bases))
                                    q1 = min(len(query_seq), query_pos + l + int(query_context_bases))
                                    ins_seq = query_seq[q0:q1]
                                poly_run, poly_frac, poly_base = _poly_at_stats(ins_seq)
                                rows.append(
                                    {
                                        "sample": sample,
                                        "chrom": chrom,
                                        "window_start": int(window_start),
                                        "window_end": int(window_end),
                                        "pos": int(pos),
                                        "clip_side": "",
                                        "clip_len": int(l),
                                        "mapq": int(read.mapping_quality),
                                        "is_reverse": bool(read.is_reverse),
                                        "read_name": str(read.query_name or ""),
                                        "has_sa": bool(read.has_tag("SA")),
                                        "sa_raw": str(read.get_tag("SA")) if read.has_tag("SA") else "",
                                        "clip_seq": str(ins_seq),
                                        "nm": int(read.get_tag("NM")) if read.has_tag("NM") else -1,
                                        "clip_poly_at_run": int(poly_run),
                                        "clip_poly_at_fraction": float(poly_frac),
                                        "clip_poly_base": str(poly_base),
                                        "poly_tail_rescued": bool(poly_run >= 8 and poly_frac >= 0.8),
                                        "evidence_type": "indel",
                                        "indel_type": "I",
                                        "indel_len": int(l),
                                    }
                                )
                        query_pos += l
                        continue
                    if op == 2:  # deletion relative to reference
                        if l >= int(min_indel_bp):
                            pos = max(1, int(ref_pos + (l // 2)))
                            overlaps = list(tree.at(pos))
                            if overlaps:
                                best = min(
                                    overlaps,
                                    key=lambda iv: (iv.end - iv.begin, abs(((iv.begin + iv.end) // 2) - pos)),
                                )
                                window_start, window_end = best.data
                                q0 = max(0, query_pos - int(query_context_bases))
                                q1 = min(len(query_seq), query_pos + int(query_context_bases))
                                del_seq = query_seq[q0:q1]
                                poly_run, poly_frac, poly_base = _poly_at_stats(del_seq)
                                rows.append(
                                    {
                                        "sample": sample,
                                        "chrom": chrom,
                                        "window_start": int(window_start),
                                        "window_end": int(window_end),
                                        "pos": int(pos),
                                        "clip_side": "",
                                        "clip_len": int(l),
                                        "mapq": int(read.mapping_quality),
                                        "is_reverse": bool(read.is_reverse),
                                        "read_name": str(read.query_name or ""),
                                        "has_sa": bool(read.has_tag("SA")),
                                        "sa_raw": str(read.get_tag("SA")) if read.has_tag("SA") else "",
                                        "clip_seq": str(del_seq),
                                        "nm": int(read.get_tag("NM")) if read.has_tag("NM") else -1,
                                        "clip_poly_at_run": int(poly_run),
                                        "clip_poly_at_fraction": float(poly_frac),
                                        "clip_poly_base": str(poly_base),
                                        "poly_tail_rescued": bool(poly_run >= 8 and poly_frac >= 0.8),
                                        "evidence_type": "indel",
                                        "indel_type": "D",
                                        "indel_len": int(l),
                                    }
                                )
                        ref_pos += l
                        continue
                    if op == 3:  # N
                        ref_pos += l
                        continue
                    if op == 4:  # S
                        query_pos += l
                        continue
                    if op == 5:  # H
                        continue
                    if op == 6:  # P
                        continue

    if not rows:
        return pd.DataFrame(
            columns=[
                "sample",
                "chrom",
                "window_start",
                "window_end",
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
                "evidence_type",
                "indel_type",
                "indel_len",
            ]
        )
    out = pd.DataFrame(rows)
    out = out.loc[out["read_name"].fillna("").astype(str).str.len() > 0].copy()
    if out.empty:
        return out
    return out.sort_values(["chrom", "window_start", "window_end", "pos", "read_name"], kind="mergesort").reset_index(drop=True)


def _build_locus_read_name_map(df: pd.DataFrame) -> dict[tuple[str, int, int], set[str]]:
    if df.empty or "read_name" not in df.columns:
        return {}
    cols = {"chrom", "window_start", "window_end", "read_name"}
    if not cols.issubset(set(df.columns)):
        return {}
    out: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    subset = df.loc[:, ["chrom", "window_start", "window_end", "read_name"]].dropna(subset=["read_name"])
    for row in subset.itertuples(index=False):
        read_name = str(row.read_name).strip()
        if not read_name:
            continue
        key = (str(row.chrom), int(row.window_start), int(row.window_end))
        out[key].add(read_name)
    return dict(out)


def _add_candidate_support_info_fields(
    candidates: pd.DataFrame,
    *,
    split_disease: pd.DataFrame,
    split_control: pd.DataFrame,
    discordant_disease: pd.DataFrame,
    discordant_control: pd.DataFrame,
    split_disease_mei: pd.DataFrame,
    split_control_mei: pd.DataFrame,
    discordant_disease_mei: pd.DataFrame,
    discordant_control_mei: pd.DataFrame,
    indel_disease: pd.DataFrame | None = None,
    indel_control: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Populate pre-assembly support strings using MEI-gated raw side counts.

    For each locus:
    - if either sample has >=1 MEI-supporting read (split or discordant),
      count assigned SR/DPE reads on each side for both disease and control.
    - otherwise set support counts to zero.
    - SR counts only MEI-mapped split clips: strict (≥20bp) or short clips
      rescued when consistent with a MEI seed. PolyA and CIGAR indels do
      not count toward SR (polyA is reported separately as polyA_MAPPED).
    - PolyA/T junction clips are never MEI_MAPPED (even if a residual tip
      remapped); they count only toward polyA_MAPPED.
    - support strings also include MEI_MAPPED (consensus remap), then
      polyA_MAPPED (mate polyA rescue ∪ junction-clip / anchor polyA) and
      VNTR_MAPPED (discordant VNTR rescue ∪ VNTR-like soft-clips demoted
      out of MEI-SR), plus polyA_side=L|R when polyA evidence has a clear
      majority flank.
    """

    out = candidates.copy()
    key_cols = ["chrom", "window_start", "window_end"]
    if out.empty:
        out["disease_supporting_reads"] = ""
        out["control_supporting_reads"] = ""
        return out

    if "insertion_breakpoint_pos" in out.columns:
        bp_tbl = out.loc[:, key_cols + ["insertion_breakpoint_pos"]].copy()
    else:
        bp_tbl = out.loc[:, key_cols].copy()
        bp_tbl["insertion_breakpoint_pos"] = 0
    bp_tbl["insertion_breakpoint_pos"] = pd.to_numeric(bp_tbl["insertion_breakpoint_pos"], errors="coerce").fillna(0).astype(int)
    midpoint = (bp_tbl["window_start"].astype(int) + bp_tbl["window_end"].astype(int)) // 2
    bp_tbl["insertion_breakpoint_pos"] = bp_tbl["insertion_breakpoint_pos"].where(bp_tbl["insertion_breakpoint_pos"] > 0, midpoint)

    def _counts_from_split(
        df: pd.DataFrame,
        prefix: str,
        informative_reads: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"])
        cols = key_cols + ["read_name"]
        if "clip_side" in df.columns:
            cols.append("clip_side")
        if "pos" in df.columns:
            cols.append("pos")
        work = df.loc[:, [c for c in cols if c in df.columns]].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"])
        if informative_reads is not None:
            # Empty informative set means zero SR (do not fall back to all clips).
            if informative_reads.empty:
                return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"])
            work = work.merge(informative_reads, on=key_cols + ["read_name"], how="inner")
            if work.empty:
                return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"])
        work = work.merge(bp_tbl, on=key_cols, how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"])

        if "clip_side" in work.columns:
            side = work["clip_side"].fillna("").astype(str).str.upper().str[:1]
            if "pos" in work.columns:
                pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
                fallback = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
                side = side.where(side.isin(["L", "R"]), fallback)
            else:
                side = side.where(side.isin(["L", "R"]), "L")
        elif "pos" in work.columns:
            pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
            side = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
        else:
            side = pd.Series(["L"] * len(work), index=work.index)
        work["raw_side"] = side

        agg = (
            work.groupby(key_cols + ["raw_side"], as_index=False)["read_name"]
            .nunique()
            .pivot_table(index=key_cols, columns="raw_side", values="read_name", fill_value=0)
            .reset_index()
        )
        agg.columns = [str(c) for c in agg.columns]
        if "L" not in agg.columns:
            agg["L"] = 0
        if "R" not in agg.columns:
            agg["R"] = 0
        agg[f"{prefix}_sr_l_raw"] = pd.to_numeric(agg["L"], errors="coerce").fillna(0).astype(int)
        agg[f"{prefix}_sr_r_raw"] = pd.to_numeric(agg["R"], errors="coerce").fillna(0).astype(int)
        return agg[key_cols + [f"{prefix}_sr_l_raw", f"{prefix}_sr_r_raw"]]

    def _counts_from_discordant(
        df: pd.DataFrame,
        prefix: str,
        informative_reads: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns or "pos" not in df.columns:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"])
        work = df.loc[:, key_cols + ["read_name", "pos"]].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"])
        if informative_reads is not None:
            if informative_reads.empty:
                return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"])
            work = work.merge(informative_reads, on=key_cols + ["read_name"], how="inner")
            if work.empty:
                return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"])
        work = work.merge(bp_tbl, on=key_cols, how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"])
        pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
        work["raw_side"] = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
        agg = (
            work.groupby(key_cols + ["raw_side"], as_index=False)["read_name"]
            .nunique()
            .pivot_table(index=key_cols, columns="raw_side", values="read_name", fill_value=0)
            .reset_index()
        )
        agg.columns = [str(c) for c in agg.columns]
        if "L" not in agg.columns:
            agg["L"] = 0
        if "R" not in agg.columns:
            agg["R"] = 0
        agg[f"{prefix}_dpe_l_raw"] = pd.to_numeric(agg["L"], errors="coerce").fillna(0).astype(int)
        agg[f"{prefix}_dpe_r_raw"] = pd.to_numeric(agg["R"], errors="coerce").fillna(0).astype(int)
        return agg[key_cols + [f"{prefix}_dpe_l_raw", f"{prefix}_dpe_r_raw"]]

    def _polya_side_table(
        disc_df: pd.DataFrame,
        split_df: pd.DataFrame,
        prefix: str,
    ) -> pd.DataFrame:
        """Majority polyA flank (L/R) among polyA_MAPPED contributors."""
        out_col = f"{prefix}_polya_side"
        parts: list[pd.DataFrame] = []

        if not split_df.empty and "read_name" in split_df.columns:
            work = split_df.copy()
            keep = pd.Series(False, index=work.index)
            if "poly_tail_rescued" in work.columns:
                keep = keep | work["poly_tail_rescued"].fillna(False).astype(bool)
            if "clip_poly_at_run" in work.columns:
                keep = keep | (
                    pd.to_numeric(work["clip_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
                )
            if bool(keep.any()):
                tmp = work.loc[keep].copy()
                if "clip_side" in tmp.columns:
                    side = tmp["clip_side"].fillna("").astype(str).str.upper().str[:1]
                else:
                    side = pd.Series("", index=tmp.index, dtype="object")
                tmp = tmp.loc[:, [c for c in key_cols + ["read_name"] if c in tmp.columns]].copy()
                tmp["poly_side"] = side.values
                tmp = tmp.loc[tmp["poly_side"].isin(["L", "R"])]
                if not tmp.empty:
                    parts.append(tmp)

        if not disc_df.empty and "read_name" in disc_df.columns:
            work = disc_df.copy()
            keep = pd.Series(False, index=work.index)
            if "polya_rescue" in work.columns:
                keep = keep | work["polya_rescue"].fillna(False).astype(bool)
            if "poly_tail_anchor_rescued" in work.columns:
                keep = keep | work["poly_tail_anchor_rescued"].fillna(False).astype(bool)
            if "anchor_poly_at_run" in work.columns:
                keep = keep | (
                    pd.to_numeric(work["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
                )
            if bool(keep.any()):
                tmp = work.loc[keep].copy()
                side = pd.Series("", index=tmp.index, dtype="object")
                if "anchor_poly_side" in tmp.columns:
                    side = tmp["anchor_poly_side"].fillna("").astype(str).str.upper().str[:1]
                if "soft_clip_side" in tmp.columns:
                    soft = tmp["soft_clip_side"].fillna("").astype(str).str.upper().str[:1]
                    side = side.where(side.isin(["L", "R"]), soft)
                need_bp = ~side.isin(["L", "R"])
                if bool(need_bp.any()) and "pos" in tmp.columns:
                    bp_map = (
                        bp_tbl.drop_duplicates(key_cols)
                        .set_index(key_cols)["insertion_breakpoint_pos"]
                    )
                    keys = list(
                        zip(
                            tmp.loc[need_bp, "chrom"].tolist(),
                            tmp.loc[need_bp, "window_start"].tolist(),
                            tmp.loc[need_bp, "window_end"].tolist(),
                        )
                    )
                    bp = pd.Series(
                        [bp_map.get(k, float("nan")) for k in keys],
                        index=tmp.loc[need_bp].index,
                        dtype="float64",
                    )
                    pos = pd.to_numeric(tmp.loc[need_bp, "pos"], errors="coerce")
                    fallback = pd.Series("L", index=pos.index).where(pos <= bp, "R")
                    side.loc[need_bp] = fallback
                out_tmp = tmp.loc[:, [c for c in key_cols + ["read_name"] if c in tmp.columns]].copy()
                out_tmp["poly_side"] = side.values
                out_tmp = out_tmp.loc[out_tmp["poly_side"].isin(["L", "R"])]
                if not out_tmp.empty:
                    parts.append(out_tmp)

        if not parts:
            return pd.DataFrame(columns=key_cols + [out_col])
        all_sides = pd.concat(parts, ignore_index=True)
        all_sides["read_name"] = all_sides["read_name"].fillna("").astype(str)
        all_sides = all_sides.loc[all_sides["read_name"].str.len() > 0]
        if all_sides.empty:
            return pd.DataFrame(columns=key_cols + [out_col])
        # Unique reads per side, then majority.
        per_side = (
            all_sides.drop_duplicates(key_cols + ["read_name", "poly_side"])
            .groupby(key_cols + ["poly_side"], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "n"})
        )
        pivot = (
            per_side.pivot_table(index=key_cols, columns="poly_side", values="n", fill_value=0)
            .reset_index()
        )
        pivot.columns = [str(c) for c in pivot.columns]
        if "L" not in pivot.columns:
            pivot["L"] = 0
        if "R" not in pivot.columns:
            pivot["R"] = 0
        l_n = pd.to_numeric(pivot["L"], errors="coerce").fillna(0).astype(int)
        r_n = pd.to_numeric(pivot["R"], errors="coerce").fillna(0).astype(int)
        pivot[out_col] = ""
        pivot.loc[l_n > r_n, out_col] = "L"
        pivot.loc[r_n > l_n, out_col] = "R"
        return pivot[key_cols + [out_col]]

    def _mei_gate(split_mei_df: pd.DataFrame, disc_mei_df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        parts: list[pd.DataFrame] = []
        for src in (split_mei_df, disc_mei_df):
            if src.empty or "read_name" not in src.columns:
                continue
            work = src.loc[:, [c for c in key_cols + ["read_name"] if c in src.columns]].copy()
            work["read_name"] = work["read_name"].fillna("").astype(str)
            work = work.loc[work["read_name"].str.len() > 0].copy()
            if not work.empty:
                parts.append(work)
        if not parts:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_has_mei_support"])
        all_mei = pd.concat(parts, ignore_index=True)
        gate = (
            all_mei.groupby(key_cols, as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": f"{prefix}_mei_unique_reads"})
        )
        gate[f"{prefix}_has_mei_support"] = gate[f"{prefix}_mei_unique_reads"].astype(int) >= 1
        return gate[key_cols + [f"{prefix}_has_mei_support"]]

    def _mei_mapped_counts(split_mei_df: pd.DataFrame, disc_mei_df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        parts: list[pd.DataFrame] = []
        for src in (split_mei_df, disc_mei_df):
            if src.empty or "read_name" not in src.columns:
                continue
            work = src.loc[:, [c for c in key_cols + ["read_name"] if c in src.columns]].copy()
            work["read_name"] = work["read_name"].fillna("").astype(str)
            work = work.loc[work["read_name"].str.len() > 0]
            if not work.empty:
                parts.append(work)
        if not parts:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_mei_mapped"])
        all_mei = pd.concat(parts, ignore_index=True)
        out_mei = (
            all_mei.groupby(key_cols, as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": f"{prefix}_mei_mapped"})
        )
        out_mei[f"{prefix}_mei_mapped"] = pd.to_numeric(out_mei[f"{prefix}_mei_mapped"], errors="coerce").fillna(0).astype(int)
        return out_mei

    def _rescue_mapped_counts(disc_df: pd.DataFrame, prefix: str, flag_col: str, out_col: str) -> pd.DataFrame:
        if disc_df.empty or "read_name" not in disc_df.columns or flag_col not in disc_df.columns:
            return pd.DataFrame(columns=key_cols + [out_col])
        work = disc_df.loc[disc_df[flag_col].fillna(False).astype(bool)].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + [out_col])
        work = work.loc[:, [c for c in key_cols + ["read_name"] if c in work.columns]]
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0]
        if work.empty:
            return pd.DataFrame(columns=key_cols + [out_col])
        out_r = (
            work.groupby(key_cols, as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": out_col})
        )
        out_r[out_col] = pd.to_numeric(out_r[out_col], errors="coerce").fillna(0).astype(int)
        return out_r

    def _polya_mapped_counts(
        disc_df: pd.DataFrame,
        split_df: pd.DataFrame,
        prefix: str,
    ) -> pd.DataFrame:
        """polyA_MAPPED = unique mate-rescue ∪ junction-clip polyA reads."""
        out_col = f"{prefix}_polya_mapped"
        parts: list[pd.DataFrame] = []

        def _take(df: pd.DataFrame, mask: pd.Series) -> None:
            if df.empty or "read_name" not in df.columns or not bool(mask.any()):
                return
            work = df.loc[mask, [c for c in key_cols + ["read_name"] if c in df.columns]].copy()
            work["read_name"] = work["read_name"].fillna("").astype(str)
            work = work.loc[work["read_name"].str.len() > 0]
            if not work.empty:
                parts.append(work)

        if not disc_df.empty:
            if "polya_rescue" in disc_df.columns:
                _take(disc_df, disc_df["polya_rescue"].fillna(False).astype(bool))
            anchor_mask = pd.Series(False, index=disc_df.index)
            if "poly_tail_anchor_rescued" in disc_df.columns:
                anchor_mask = anchor_mask | disc_df["poly_tail_anchor_rescued"].fillna(False).astype(bool)
            if "anchor_poly_at_run" in disc_df.columns:
                anchor_mask = anchor_mask | (
                    pd.to_numeric(disc_df["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
                )
            _take(disc_df, anchor_mask)

        if not split_df.empty:
            split_mask = pd.Series(False, index=split_df.index)
            if "poly_tail_rescued" in split_df.columns:
                split_mask = split_mask | split_df["poly_tail_rescued"].fillna(False).astype(bool)
            if "clip_poly_at_run" in split_df.columns:
                split_mask = split_mask | (
                    pd.to_numeric(split_df["clip_poly_at_run"], errors="coerce").fillna(0).astype(int) >= 8
                )
            _take(split_df, split_mask)

        if not parts:
            return pd.DataFrame(columns=key_cols + [out_col])
        all_reads = pd.concat(parts, ignore_index=True)
        out_r = (
            all_reads.groupby(key_cols, as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": out_col})
        )
        out_r[out_col] = pd.to_numeric(out_r[out_col], errors="coerce").fillna(0).astype(int)
        return out_r

    def _rescue_polya_max_len(
        disc_df: pd.DataFrame,
        split_df: pd.DataFrame,
        prefix: str,
    ) -> pd.DataFrame:
        """Longest observed polyA/T among polyA_MAPPED contributors.

        Includes mate polyA rescues and junction-clip / anchor polyA reads.
        When a DPE has both a junction polyA clip and a polyA mate, length is
        clip + mate (minimum non-overlapping tail), capped at
        ``_MAX_POLYA_PAIR_BP``.
        """
        out_col = f"{prefix}_polya_rescue_max_len"
        lens: list[pd.DataFrame] = []

        if not disc_df.empty and "read_name" in disc_df.columns:
            work = disc_df.copy()
            mate_lens = (
                work["mate_seq"].fillna("").astype(str).map(_polya_tail_width_bp)
                if "mate_seq" in work.columns
                else pd.Series(0, index=work.index, dtype="int64")
            )
            soft = (
                pd.to_numeric(work["soft_clip_len"], errors="coerce").fillna(0).astype(int)
                if "soft_clip_len" in work.columns
                else pd.Series(0, index=work.index, dtype="int64")
            )
            run = (
                pd.to_numeric(work["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int)
                if "anchor_poly_at_run" in work.columns
                else pd.Series(0, index=work.index, dtype="int64")
            )
            rescued = (
                work["poly_tail_anchor_rescued"].fillna(False).astype(bool)
                if "poly_tail_anchor_rescued" in work.columns
                else pd.Series(False, index=work.index)
            )
            read_lens = (
                work["read_seq"].fillna("").astype(str).map(_observed_poly_at_len_bp)
                if "read_seq" in work.columns
                else pd.Series(0, index=work.index, dtype="int64")
            )
            # Prefer soft-clip length when present; otherwise use observed poly
            # length on the anchor read (covers extracts without soft_clip_len).
            anchor_clip = pd.concat([soft, run, read_lens], axis=1).max(axis=1).astype(int)
            dpe_len = [
                _dpe_polya_observed_len_bp(
                    mate_len=int(m),
                    anchor_clip_len=int(a),
                    anchor_poly_run=int(r),
                    poly_tail_anchor_rescued=bool(p),
                )
                for m, a, r, p in zip(
                    mate_lens.tolist(),
                    anchor_clip.tolist(),
                    run.tolist(),
                    rescued.tolist(),
                )
            ]
            polya_rescue = (
                work["polya_rescue"].fillna(False).astype(bool)
                if "polya_rescue" in work.columns
                else pd.Series(False, index=work.index)
            )
            keep_mask = polya_rescue | rescued | run.ge(8) | pd.Series(dpe_len, index=work.index).gt(0)
            if bool(keep_mask.any()):
                tmp = work.loc[keep_mask, [c for c in key_cols if c in work.columns]].copy()
                tmp[out_col] = [
                    max(int(d), int(rl))
                    for d, rl in zip(
                        pd.Series(dpe_len, index=work.index).loc[keep_mask].tolist(),
                        read_lens.loc[keep_mask].tolist(),
                    )
                ]
                tmp = tmp.loc[pd.to_numeric(tmp[out_col], errors="coerce").fillna(0).astype(int) > 0]
                if not tmp.empty:
                    lens.append(tmp)

        if not split_df.empty and "read_name" in split_df.columns:
            work = split_df.copy()
            if "clip_seq" in work.columns:
                clip_lens = work["clip_seq"].fillna("").astype(str).map(_observed_poly_at_len_bp)
            elif "clip_poly_at_run" in work.columns:
                clip_lens = pd.to_numeric(work["clip_poly_at_run"], errors="coerce").fillna(0).astype(int)
            else:
                clip_lens = pd.Series(0, index=work.index, dtype="int64")
            rescued = (
                work["poly_tail_rescued"].fillna(False).astype(bool)
                if "poly_tail_rescued" in work.columns
                else pd.Series(False, index=work.index)
            )
            run = (
                pd.to_numeric(work["clip_poly_at_run"], errors="coerce").fillna(0).astype(int)
                if "clip_poly_at_run" in work.columns
                else pd.Series(0, index=work.index, dtype="int64")
            )
            keep_mask = rescued | run.ge(8) | clip_lens.ge(8)
            if bool(keep_mask.any()):
                tmp = work.loc[keep_mask, [c for c in key_cols if c in work.columns]].copy()
                tmp[out_col] = clip_lens.loc[keep_mask].astype(int).clip(upper=_MAX_POLYA_SINGLE_BP)
                tmp = tmp.loc[pd.to_numeric(tmp[out_col], errors="coerce").fillna(0).astype(int) > 0]
                if not tmp.empty:
                    lens.append(tmp)

        if not lens:
            return pd.DataFrame(columns=key_cols + [out_col])
        all_lens = pd.concat(lens, ignore_index=True)
        return all_lens.groupby(key_cols, as_index=False)[out_col].max()

    def _polyA_split_read_table(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        cols = [c for c in key_cols + ["read_name", "clip_poly_at_run", "poly_tail_rescued"] if c in df.columns]
        work = df.loc[:, cols].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        if "clip_poly_at_run" in work.columns:
            poly_run = pd.to_numeric(work["clip_poly_at_run"], errors="coerce").fillna(0).astype(int)
        else:
            poly_run = pd.Series(0, index=work.index, dtype="int64")
        if "poly_tail_rescued" in work.columns:
            rescued = work["poly_tail_rescued"].fillna(False).astype(bool)
        else:
            rescued = pd.Series(False, index=work.index)
        keep = (poly_run >= 8) | rescued
        return work.loc[keep, key_cols + ["read_name"]].drop_duplicates()

    def _polyA_discordant_read_table(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        cols = [c for c in key_cols + ["read_name", "anchor_poly_at_run", "poly_tail_anchor_rescued"] if c in df.columns]
        work = df.loc[:, cols].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        if "anchor_poly_at_run" in work.columns:
            poly_run = pd.to_numeric(work["anchor_poly_at_run"], errors="coerce").fillna(0).astype(int)
        else:
            poly_run = pd.Series(0, index=work.index, dtype="int64")
        if "poly_tail_anchor_rescued" in work.columns:
            rescued = work["poly_tail_anchor_rescued"].fillna(False).astype(bool)
        else:
            rescued = pd.Series(False, index=work.index)
        keep = (poly_run >= 8) | rescued
        return work.loc[keep, key_cols + ["read_name"]].drop_duplicates()

    def _mei_read_table(split_mei_df: pd.DataFrame, disc_mei_df: pd.DataFrame) -> pd.DataFrame:
        parts: list[pd.DataFrame] = []
        for src in (split_mei_df, disc_mei_df):
            if src.empty or "read_name" not in src.columns:
                continue
            work = src.loc[:, [c for c in key_cols + ["read_name"] if c in src.columns]].copy()
            work["read_name"] = work["read_name"].fillna("").astype(str)
            work = work.loc[work["read_name"].str.len() > 0, key_cols + ["read_name"]].drop_duplicates()
            if not work.empty:
                parts.append(work)
        if not parts:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        return pd.concat(parts, ignore_index=True).drop_duplicates()

    def _unique_read_count_table(df: pd.DataFrame, out_col: str) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns:
            return pd.DataFrame(columns=key_cols + [out_col])
        work = df.loc[:, [c for c in key_cols + ["read_name"] if c in df.columns]].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0]
        if work.empty:
            return pd.DataFrame(columns=key_cols + [out_col])
        return work.groupby(key_cols, as_index=False)["read_name"].nunique().rename(columns={"read_name": out_col})

    disease_gate = _mei_gate(split_disease_mei, discordant_disease_mei, "disease")
    control_gate = _mei_gate(split_control_mei, discordant_control_mei, "control")
    locus_gate = (
        bp_tbl.loc[:, key_cols]
        .drop_duplicates()
        .merge(disease_gate, on=key_cols, how="left")
        .merge(control_gate, on=key_cols, how="left")
    )
    disease_gate_series = (
        locus_gate["disease_has_mei_support"].fillna(False).astype(bool)
        if "disease_has_mei_support" in locus_gate.columns
        else pd.Series(False, index=locus_gate.index)
    )
    control_gate_series = (
        locus_gate["control_has_mei_support"].fillna(False).infer_objects(copy=False).astype(bool)
        if "control_has_mei_support" in locus_gate.columns
        else pd.Series(False, index=locus_gate.index)
    )
    locus_has_mei_support = disease_gate_series | control_gate_series
    locus_gate["locus_has_mei_support"] = locus_has_mei_support

    sample_frames = {
        "disease": (
            split_disease,
            discordant_disease,
            split_disease_mei,
            discordant_disease_mei,
            indel_disease if indel_disease is not None else pd.DataFrame(),
        ),
        "control": (
            split_control,
            discordant_control,
            split_control_mei,
            discordant_control_mei,
            indel_control if indel_control is not None else pd.DataFrame(),
        ),
    }

    for prefix in ("disease", "control"):
        split_df, disc_df, split_mei_df, disc_mei_df, indel_df = sample_frames[prefix]
        other_prefix = "control" if prefix == "disease" else "disease"
        split_support_df = (
            pd.concat([split_df, indel_df], ignore_index=True)
            if (indel_df is not None and not indel_df.empty)
            else split_df
        )

        poly_split_reads = _polyA_split_read_table(split_support_df)
        poly_disc_reads = _polyA_discordant_read_table(disc_df)
        mei_split_reads = _mei_read_table(split_mei_df, split_mei_df.iloc[0:0].copy())
        mei_disc_reads = _mei_read_table(disc_mei_df.iloc[0:0].copy(), disc_mei_df)
        indel_reads = (
            split_support_df.loc[
                split_support_df.get("evidence_type", "").fillna("").astype(str) == "indel",
                key_cols + ["read_name"],
            ].drop_duplicates()
            if ("evidence_type" in split_support_df.columns and "read_name" in split_support_df.columns)
            else pd.DataFrame(columns=key_cols + ["read_name"])
        )
        strict_reads = pd.concat(
            [
                poly_split_reads,
                poly_disc_reads,
                mei_split_reads,
                mei_disc_reads,
                indel_reads,
            ],
            ignore_index=True,
        ).drop_duplicates()

        # SR = MEI-mapped split clips only (strict ≥20bp or short-seed-rescued).
        # PolyA stays in polyA_MAPPED; CIGAR indels are not SR.
        sr_raw = _counts_from_split(split_df, prefix, informative_reads=mei_split_reads)
        dpe_raw = _counts_from_discordant(disc_df, prefix)
        sr_strict = _counts_from_split(split_df, prefix, informative_reads=mei_split_reads)
        dpe_strict = _counts_from_discordant(disc_df, prefix, informative_reads=strict_reads)
        mei_mapped = _mei_mapped_counts(split_mei_df, disc_mei_df, prefix)
        polya_mapped = _polya_mapped_counts(disc_df, split_support_df, prefix)
        polya_side = _polya_side_table(disc_df, split_support_df, prefix)
        # VNTR_MAPPED = discordant VNTR rescue ∪ demoted/rescued VNTR-like soft-clips.
        vntr_parts = []
        for src in (disc_df, split_df):
            if src is None or src.empty or "vntr_rescue" not in src.columns or "read_name" not in src.columns:
                continue
            if not all(c in src.columns for c in key_cols):
                continue
            part = src.loc[src["vntr_rescue"].fillna(False).astype(bool), key_cols + ["read_name"]].copy()
            if not part.empty:
                vntr_parts.append(part)
        if not vntr_parts:
            vntr_mapped = pd.DataFrame(columns=key_cols + [f"{prefix}_vntr_mapped"])
        else:
            vntr_mapped = (
                pd.concat(vntr_parts, ignore_index=True)
                .drop_duplicates()
                .groupby(key_cols, as_index=False)["read_name"]
                .nunique()
                .rename(columns={"read_name": f"{prefix}_vntr_mapped"})
            )
            vntr_mapped[f"{prefix}_vntr_mapped"] = (
                pd.to_numeric(vntr_mapped[f"{prefix}_vntr_mapped"], errors="coerce")
                .fillna(0)
                .astype(int)
            )
        polya_max_len = _rescue_polya_max_len(disc_df, split_support_df, prefix)

        merged = (
            bp_tbl.loc[:, key_cols]
            .drop_duplicates()
            .merge(sr_raw, on=key_cols, how="left")
            .merge(dpe_raw, on=key_cols, how="left")
            .merge(
                sr_strict.rename(
                    columns={
                        f"{prefix}_sr_l_raw": f"{prefix}_sr_l_strict",
                        f"{prefix}_sr_r_raw": f"{prefix}_sr_r_strict",
                    }
                ),
                on=key_cols,
                how="left",
            )
            .merge(
                dpe_strict.rename(
                    columns={
                        f"{prefix}_dpe_l_raw": f"{prefix}_dpe_l_strict",
                        f"{prefix}_dpe_r_raw": f"{prefix}_dpe_r_strict",
                    }
                ),
                on=key_cols,
                how="left",
            )
            .merge(mei_mapped, on=key_cols, how="left")
            .merge(polya_mapped, on=key_cols, how="left")
            .merge(polya_side, on=key_cols, how="left")
            .merge(vntr_mapped, on=key_cols, how="left")
            .merge(polya_max_len, on=key_cols, how="left")
            .merge(locus_gate.loc[:, key_cols + ["locus_has_mei_support"]], on=key_cols, how="left")
        )

        sample_has_mei_col = f"{prefix}_has_mei_support"
        other_has_mei_col = f"{other_prefix}_has_mei_support"
        status_tbl = locus_gate.loc[:, key_cols].copy()
        status_tbl["sample_has_mei_support"] = (
            locus_gate[sample_has_mei_col].fillna(False).infer_objects(copy=False).astype(bool)
            if sample_has_mei_col in locus_gate.columns
            else pd.Series(False, index=locus_gate.index)
        )
        status_tbl["other_has_mei_support"] = (
            locus_gate[other_has_mei_col].fillna(False).infer_objects(copy=False).astype(bool)
            if other_has_mei_col in locus_gate.columns
            else pd.Series(False, index=locus_gate.index)
        )
        merged = merged.merge(status_tbl, on=key_cols, how="left")

        has_locus_mei_support = merged["locus_has_mei_support"].fillna(False).astype(bool)
        weakness = (
            bp_tbl.loc[:, key_cols]
            .drop_duplicates()
            .merge(_unique_read_count_table(poly_split_reads, "poly_split_reads"), on=key_cols, how="left")
            .merge(_unique_read_count_table(mei_split_reads, "mei_split_reads"), on=key_cols, how="left")
            .merge(_unique_read_count_table(mei_disc_reads, "mei_disc_reads"), on=key_cols, how="left")
        )
        merged = merged.merge(weakness, on=key_cols, how="left")
        poly_split_ct = pd.to_numeric(merged.get("poly_split_reads", 0), errors="coerce").fillna(0).astype(int)
        mei_split_ct = pd.to_numeric(merged.get("mei_split_reads", 0), errors="coerce").fillna(0).astype(int)
        mei_disc_ct = pd.to_numeric(merged.get("mei_disc_reads", 0), errors="coerce").fillna(0).astype(int)
        weak_mei_only_discordant = (mei_split_ct <= 0) & (poly_split_ct <= 0) & (mei_disc_ct <= 2)
        sr_l_raw = pd.to_numeric(merged.get(f"{prefix}_sr_l_raw", 0), errors="coerce").fillna(0).astype(int)
        sr_r_raw = pd.to_numeric(merged.get(f"{prefix}_sr_r_raw", 0), errors="coerce").fillna(0).astype(int)
        dpe_l_raw = pd.to_numeric(merged.get(f"{prefix}_dpe_l_raw", 0), errors="coerce").fillna(0).astype(int)
        dpe_r_raw = pd.to_numeric(merged.get(f"{prefix}_dpe_r_raw", 0), errors="coerce").fillna(0).astype(int)
        sr_total_raw = sr_l_raw + sr_r_raw
        dpe_total_raw = dpe_l_raw + dpe_r_raw
        weak_window_only_dpe = (sr_total_raw <= 0) & (dpe_total_raw <= 2) & (poly_split_ct <= 0)
        strict_mode = merged["other_has_mei_support"].fillna(False).astype(bool) & (
            (~merged["sample_has_mei_support"].fillna(False).astype(bool))
            | weak_mei_only_discordant
            | weak_window_only_dpe
        )

        sr_l_strict = pd.to_numeric(merged.get(f"{prefix}_sr_l_strict", 0), errors="coerce").fillna(0).astype(int)
        sr_r_strict = pd.to_numeric(merged.get(f"{prefix}_sr_r_strict", 0), errors="coerce").fillna(0).astype(int)
        dpe_l_strict = pd.to_numeric(merged.get(f"{prefix}_dpe_l_strict", 0), errors="coerce").fillna(0).astype(int)
        dpe_r_strict = pd.to_numeric(merged.get(f"{prefix}_dpe_r_strict", 0), errors="coerce").fillna(0).astype(int)
        mei_mapped_total = pd.to_numeric(merged.get(f"{prefix}_mei_mapped", 0), errors="coerce").fillna(0).astype(int)
        polya_mapped_total = pd.to_numeric(merged.get(f"{prefix}_polya_mapped", 0), errors="coerce").fillna(0).astype(int)
        vntr_mapped_total = pd.to_numeric(merged.get(f"{prefix}_vntr_mapped", 0), errors="coerce").fillna(0).astype(int)

        # SR is always informative-only; strict_mode still gates DPE raw vs strict.
        sr_l = sr_l_raw.where(has_locus_mei_support, 0)
        sr_r = sr_r_raw.where(has_locus_mei_support, 0)
        dpe_l = dpe_l_raw.where(~strict_mode, dpe_l_strict).where(has_locus_mei_support, 0)
        dpe_r = dpe_r_raw.where(~strict_mode, dpe_r_strict).where(has_locus_mei_support, 0)

        # VNTR_MAPPED is SVA-only; omit the token for ALU/LINE1.
        fam_series = pd.Series("", index=merged.index, dtype="object")
        for fam_col in (
            "consensus_mei_family",
            "mei_family",
            "disease_discordant_mei_family",
            "control_discordant_mei_family",
        ):
            if fam_col in out.columns:
                fam_tbl = out.loc[:, key_cols + [fam_col]].drop_duplicates(key_cols, keep="first")
                merged = merged.merge(fam_tbl, on=key_cols, how="left", suffixes=("", "_fam"))
                use_col = fam_col if fam_col in merged.columns else f"{fam_col}_fam"
                if use_col in merged.columns:
                    fam_series = merged[use_col].fillna("").astype(str)
                    break
        is_sva = fam_series.str.upper().eq("SVA")
        polya_side_vals = merged.get(f"{prefix}_polya_side", pd.Series("", index=merged.index)).fillna("").astype(str)
        polya_side_vals = polya_side_vals.where(polya_side_vals.isin(["L", "R"]), "")

        merged[f"{prefix}_supporting_reads"] = [
            (
                f"SR_L={sl},SR_R={srx},DPE_L={dl},DPE_R={dr},"
                f"MEI_MAPPED={mm},polyA_MAPPED={pa}"
                + (f",VNTR_MAPPED={vn}" if sva else "")
                + (f",polyA_side={ps}" if ps else "")
            )
            for sl, srx, dl, dr, mm, pa, vn, sva, ps in zip(
                sr_l.tolist(),
                sr_r.tolist(),
                dpe_l.tolist(),
                dpe_r.tolist(),
                mei_mapped_total.tolist(),
                polya_mapped_total.tolist(),
                vntr_mapped_total.tolist(),
                is_sva.tolist(),
                polya_side_vals.tolist(),
            )
        ]
        if f"{prefix}_supporting_reads" in out.columns:
            out = out.drop(columns=[f"{prefix}_supporting_reads"])
        merge_cols = key_cols + [f"{prefix}_supporting_reads"]
        max_len_col = f"{prefix}_polya_rescue_max_len"
        if max_len_col in merged.columns:
            merge_cols.append(max_len_col)
        out = out.merge(merged[merge_cols], on=key_cols, how="left")
        out[f"{prefix}_supporting_reads"] = (
            out[f"{prefix}_supporting_reads"]
            .fillna(
                "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0"
            )
            .astype(str)
        )
        if max_len_col in out.columns:
            out[max_len_col] = pd.to_numeric(out[max_len_col], errors="coerce").fillna(0).astype(int)
    return out


def _aggregate_side_metrics(
    df: pd.DataFrame,
    sample_prefix: str,
    side: str,
    preferred_subfamily_by_locus: dict[tuple[str, int, int], str] | None = None,
) -> pd.DataFrame:
    if df is None or df.empty or "clip_side" not in df.columns:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_{side}_mei_supported_reads",
                f"{sample_prefix}_{side}_mei_score_sum",
                f"{sample_prefix}_{side}_mei_family",
                f"{sample_prefix}_{side}_mei_subfamily",
                f"{sample_prefix}_{side}_mei_strand",
                f"{sample_prefix}_{side}_mei_start",
                f"{sample_prefix}_{side}_mei_end",
                f"{sample_prefix}_{side}_mei_anchor_bp_max",
                f"{sample_prefix}_{side}_mei_target_len",
                f"{sample_prefix}_{side}_mei_subfamily_purity",
                f"{sample_prefix}_{side}_mei_breakpoint_mode",
                f"{sample_prefix}_{side}_mei_breakpoint_mode_fraction",
                f"{sample_prefix}_{side}_mei_breakpoint_unique_positions",
                f"{sample_prefix}_{side}_poly_at_reads",
                f"{sample_prefix}_{side}_poly_at_fraction",
                f"{sample_prefix}_{side}_poly_at_max_run",
                f"{sample_prefix}_{side}_clip_overlap_informative_reads",
                f"{sample_prefix}_{side}_clip_overlap_jaccard_median",
                f"{sample_prefix}_{side}_clip_overlap_pair_ge20_fraction",
            ]
        )
    mei_hit_raw = df["mei_hit"].fillna(False).astype(bool) if "mei_hit" in df.columns else pd.Series(False, index=df.index)
    mei_hit_coord = (
        df["mei_hit_coord"].fillna(False).astype(bool) if "mei_hit_coord" in df.columns else pd.Series(False, index=df.index)
    )
    side_df = df.loc[(df["clip_side"] == side) & (mei_hit_raw | mei_hit_coord)].copy()
    if side_df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_{side}_mei_supported_reads",
                f"{sample_prefix}_{side}_mei_score_sum",
                f"{sample_prefix}_{side}_mei_family",
                f"{sample_prefix}_{side}_mei_subfamily",
                f"{sample_prefix}_{side}_mei_strand",
                f"{sample_prefix}_{side}_mei_start",
                f"{sample_prefix}_{side}_mei_end",
                f"{sample_prefix}_{side}_mei_anchor_bp_max",
                f"{sample_prefix}_{side}_mei_target_len",
                f"{sample_prefix}_{side}_mei_subfamily_purity",
                f"{sample_prefix}_{side}_mei_breakpoint_mode",
                f"{sample_prefix}_{side}_mei_breakpoint_mode_fraction",
                f"{sample_prefix}_{side}_mei_breakpoint_unique_positions",
                f"{sample_prefix}_{side}_poly_at_reads",
                f"{sample_prefix}_{side}_poly_at_fraction",
                f"{sample_prefix}_{side}_poly_at_max_run",
                f"{sample_prefix}_{side}_clip_overlap_informative_reads",
                f"{sample_prefix}_{side}_clip_overlap_jaccard_median",
                f"{sample_prefix}_{side}_clip_overlap_pair_ge20_fraction",
            ]
        )

    def poly_at_max_run(seq: str) -> int:
        return _observed_poly_at_len_bp(seq)

    side_df["poly_at_max_run"] = side_df["clip_seq"].fillna("").astype(str).map(poly_at_max_run)
    side_df["poly_at_flag"] = (side_df["poly_at_max_run"] >= 8).astype(int)
    side_df["mei_score_effective"] = pd.to_numeric(side_df.get("mei_score", 0.0), errors="coerce").fillna(0.0)
    if "mei_score_coord" in side_df.columns:
        score_coord = pd.to_numeric(side_df["mei_score_coord"], errors="coerce").fillna(0.0)
        side_df.loc[side_df["mei_score_effective"] <= 0.0, "mei_score_effective"] = score_coord.loc[
            side_df["mei_score_effective"] <= 0.0
        ]
    side_df["family_effective"] = side_df.get("family", "").fillna("").astype(str)
    if "family_coord" in side_df.columns:
        coord_family = side_df["family_coord"].fillna("").astype(str)
        side_df.loc[side_df["family_effective"] == "", "family_effective"] = coord_family.loc[side_df["family_effective"] == ""]
    side_df["target_effective"] = side_df.get("target", "").fillna("").astype(str)
    if "target_coord" in side_df.columns:
        coord_target = side_df["target_coord"].fillna("").astype(str)
        side_df.loc[side_df["target_effective"] == "", "target_effective"] = coord_target.loc[side_df["target_effective"] == ""]
    side_df["target_strand_effective"] = side_df.get("target_strand", "").fillna("").astype(str)
    if "target_strand_coord" in side_df.columns:
        coord_strand = side_df["target_strand_coord"].fillna("").astype(str)
        side_df.loc[side_df["target_strand_effective"] == "", "target_strand_effective"] = coord_strand.loc[
            side_df["target_strand_effective"] == ""
        ]
    side_df["alnlen_effective"] = pd.to_numeric(side_df.get("alnlen", 0), errors="coerce").fillna(0).astype(int)
    if "alnlen_coord" in side_df.columns:
        coord_alnlen = pd.to_numeric(side_df["alnlen_coord"], errors="coerce").fillna(0).astype(int)
        side_df["alnlen_effective"] = pd.concat(
            [side_df["alnlen_effective"], coord_alnlen], axis=1
        ).max(axis=1).astype(int)
    side_df["target_len_effective"] = pd.to_numeric(side_df.get("target_len", 0), errors="coerce").fillna(0).astype(int)
    if "target_len_coord" in side_df.columns:
        coord_tlen = pd.to_numeric(side_df["target_len_coord"], errors="coerce").fillna(0).astype(int)
        side_df["target_len_effective"] = pd.concat(
            [side_df["target_len_effective"], coord_tlen], axis=1
        ).max(axis=1).astype(int)
    # Prefer per-row coord-remap fields only when they carry a real hit; otherwise
    # fall back to the primary target_* columns. Blindly selecting *_coord when the
    # column exists (even if all zeros) previously zeroed L/R mei_start/end despite
    # valid split-read MEI mappings.
    def _coalesce_coord_numeric(primary: str, coord: str) -> pd.Series:
        base = pd.to_numeric(side_df.get(primary, 0), errors="coerce").fillna(0)
        if coord in side_df.columns:
            alt = pd.to_numeric(side_df[coord], errors="coerce").fillna(0)
            return alt.where(alt.gt(0), base)
        return base

    def _coalesce_coord_string(primary: str, coord: str) -> pd.Series:
        base = side_df.get(primary, pd.Series("", index=side_df.index)).fillna("").astype(str)
        if coord in side_df.columns:
            alt = side_df[coord].fillna("").astype(str)
            return alt.where(alt.str.len().gt(0), base)
        return base

    side_df["coord_target"] = _coalesce_coord_string("target", "target_coord")
    side_df["coord_target_start"] = _coalesce_coord_numeric("target_start", "target_start_coord").astype(int)
    side_df["coord_target_end"] = _coalesce_coord_numeric("target_end", "target_end_coord").astype(int)
    side_df["coord_target_len"] = _coalesce_coord_numeric("target_len", "target_len_coord").astype(int)
    side_df["coord_mapq"] = _coalesce_coord_numeric("mapq", "mapq_coord")
    side_df["coord_qcov"] = _coalesce_coord_numeric("qcov", "qcov_coord")
    side_df["coord_pid"] = _coalesce_coord_numeric("pid", "pid_coord")
    side_df["coord_alnlen"] = _coalesce_coord_numeric("alnlen", "alnlen_coord").astype(int)
    if "mei_hit_coord" in side_df.columns:
        coord_hit_mask = (
            side_df["mei_hit"].fillna(False).astype(bool)
            | side_df["mei_hit_coord"].fillna(False).astype(bool)
        )
    else:
        coord_hit_mask = side_df["mei_hit"].fillna(False).astype(bool)

    family_top, subfamily_top = _top_family_then_subfamily(
        side_df,
        group_cols=["chrom", "window_start", "window_end"],
        family_col="family_effective",
        subfamily_col="target_effective",
        score_col="mei_score_effective",
    )
    family_top = family_top.rename(columns={"family_effective": f"{sample_prefix}_{side}_mei_family"})
    subfamily_top = subfamily_top.rename(columns={"target_effective": f"{sample_prefix}_{side}_mei_subfamily"})
    # Coordinate-estimation subset: avoid low-confidence/polyA-only hits that can
    # collapse inferred spans to tail-length artifacts. Prefer non-poly clips, but
    # do not hard-drop poly-containing clips that still have a real MEI alignment
    # (common at the 3' junction where the softclip includes the polyA tail).
    coord_mapq = side_df["coord_mapq"]
    coord_qcov = side_df["coord_qcov"]
    coord_pid = side_df["coord_pid"]
    coord_alnlen = side_df["coord_alnlen"]
    has_mei_coords = side_df["coord_target_start"].gt(0) & side_df["coord_target_end"].ge(
        side_df["coord_target_start"]
    )
    # MEI clips can be highly repetitive with low MAPQ despite very high identity/coverage.
    # Permit these as coordinate candidates if alignment quality itself is strong.
    strong_repetitive_clip = (coord_qcov >= 0.90) & (coord_pid >= 0.90) & (coord_alnlen >= 40)
    quality_ok_strict = (
        coord_hit_mask
        & has_mei_coords
        & ((coord_mapq >= 20) | strong_repetitive_clip)
        & (coord_qcov >= 0.60)
        & (coord_pid >= 0.85)
        & (coord_alnlen >= _MIN_MEI_ANCHOR_BP)
    )
    quality_ok_relaxed = (
        coord_hit_mask
        & has_mei_coords
        & ((coord_mapq >= 10) | strong_repetitive_clip)
        & (coord_qcov >= 0.35)
        & (coord_pid >= 0.75)
        & (coord_alnlen >= _MIN_MEI_ANCHOR_BP_RELAXED)
    )
    non_poly = side_df["poly_at_flag"] == 0
    coord_df = side_df.loc[quality_ok_strict & non_poly].copy()
    if coord_df.empty:
        coord_df = side_df.loc[quality_ok_relaxed & non_poly].copy()
    if coord_df.empty:
        # Keep poly-tailed junction clips with otherwise strong MEI mappings.
        coord_df = side_df.loc[quality_ok_strict].copy()
    if coord_df.empty:
        coord_df = side_df.loc[quality_ok_relaxed].copy()
    if not coord_df.empty:
        # Columns already coalesced above.
        pass
    if not coord_df.empty:
        coord_df = coord_df.merge(
            subfamily_top[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_{side}_mei_subfamily",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        coord_df["coord_side_target"] = coord_df[f"{sample_prefix}_{side}_mei_subfamily"].fillna("").astype(str)
        if preferred_subfamily_by_locus:
            pref_tbl = pd.DataFrame(
                [
                    {
                        "chrom": chrom,
                        "window_start": int(ws),
                        "window_end": int(we),
                        "coord_preferred_target": str(target or ""),
                    }
                    for (chrom, ws, we), target in preferred_subfamily_by_locus.items()
                    if str(target or "") != ""
                ]
            )
            if not pref_tbl.empty:
                coord_df = coord_df.merge(
                    pref_tbl,
                    on=["chrom", "window_start", "window_end"],
                    how="left",
                )
            else:
                coord_df["coord_preferred_target"] = ""
        else:
            coord_df["coord_preferred_target"] = ""
        coord_df["coord_preferred_target"] = coord_df["coord_preferred_target"].fillna("").astype(str)
        pref_match = (
            coord_df["coord_preferred_target"].str.len().gt(0)
            & coord_df["coord_target"].eq(coord_df["coord_preferred_target"])
        )
        side_match = coord_df["coord_target"].eq(coord_df["coord_side_target"])
        keep_mask = pref_match | side_match
        filtered = coord_df.loc[keep_mask].copy()
        # If preferred/side-subfamily filtering removes every hit (common for
        # multi-subfamily SVA/L1 loci), keep the unfiltered coordinate set so
        # 5'/3' footprint estimation is not left empty.
        if filtered.empty:
            coord_df = coord_df.copy()
        else:
            coord_df = filtered
            coord_df["coord_target_rank"] = pd.Series([2] * len(coord_df), index=coord_df.index)
            coord_df.loc[side_match.loc[coord_df.index], "coord_target_rank"] = 1
            coord_df.loc[pref_match.loc[coord_df.index], "coord_target_rank"] = 0
            best_rank = coord_df.groupby(["chrom", "window_start", "window_end"])["coord_target_rank"].transform("min")
            coord_df = coord_df.loc[coord_df["coord_target_rank"] == best_rank].copy()
    coord_agg = (
        coord_df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_{side}_mei_start": ("coord_target_start", "min"),
                f"{sample_prefix}_{side}_mei_end": ("coord_target_end", "max"),
            }
        )
        if not coord_df.empty
        else pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_{side}_mei_start",
                f"{sample_prefix}_{side}_mei_end",
            ]
        )
    )
    strand_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "target_strand_effective"], as_index=False)["mei_score_effective"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score_effective"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target_strand_effective": f"{sample_prefix}_{side}_mei_strand"})
    )
    subfamily_totals = (
        side_df.groupby(["chrom", "window_start", "window_end", "target_effective"], as_index=False)["mei_score_effective"]
        .sum()
        .rename(columns={"mei_score_effective": "subfamily_score_sum"})
    )
    subfamily_sum = (
        side_df.groupby(["chrom", "window_start", "window_end"], as_index=False)["mei_score_effective"]
        .sum()
        .rename(columns={"mei_score_effective": "all_subfamily_score_sum"})
    )
    purity = (
        subfamily_top.rename(columns={"mei_score_effective": "top_subfamily_score_sum"})
        .merge(
            subfamily_sum,
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
    )
    purity[f"{sample_prefix}_{side}_mei_subfamily_purity"] = (
        purity["top_subfamily_score_sum"] / purity["all_subfamily_score_sum"]
    ).fillna(0.0)

    pos_counts = (
        side_df.groupby(["chrom", "window_start", "window_end", "pos"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": "support_reads"})
    )
    pos_counts = pos_counts.sort_values(
        ["chrom", "window_start", "window_end", "support_reads", "pos"],
        ascending=[True, True, True, False, True],
    )
    pos_mode = pos_counts.drop_duplicates(["chrom", "window_start", "window_end"], keep="first").rename(
        columns={"pos": f"{sample_prefix}_{side}_mei_breakpoint_mode", "support_reads": "mode_support_reads"}
    )
    pos_unique = (
        pos_counts.groupby(["chrom", "window_start", "window_end"], as_index=False)["pos"]
        .nunique()
        .rename(columns={"pos": f"{sample_prefix}_{side}_mei_breakpoint_unique_positions"})
    )
    overlap_rows: list[dict[str, object]] = []
    for (chrom, ws, we), grp in side_df.groupby(["chrom", "window_start", "window_end"], sort=False):
        informative_reads, jaccard_median, pair_ge20_fraction = _clip_overlap_consistency_stats(grp)
        overlap_rows.append(
            {
                "chrom": str(chrom),
                "window_start": int(ws),
                "window_end": int(we),
                f"{sample_prefix}_{side}_clip_overlap_informative_reads": int(informative_reads),
                f"{sample_prefix}_{side}_clip_overlap_jaccard_median": float(jaccard_median),
                f"{sample_prefix}_{side}_clip_overlap_pair_ge20_fraction": float(pair_ge20_fraction),
            }
        )
    overlap_tbl = pd.DataFrame(overlap_rows)
    agg = (
        side_df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_{side}_mei_supported_reads": ("read_name", "nunique"),
                f"{sample_prefix}_{side}_mei_score_sum": ("mei_score_effective", "sum"),
                f"{sample_prefix}_{side}_mei_anchor_bp_max": ("alnlen_effective", "max"),
                f"{sample_prefix}_{side}_mei_target_len": ("target_len_effective", "max"),
                f"{sample_prefix}_{side}_poly_at_reads": ("poly_at_flag", "sum"),
                f"{sample_prefix}_{side}_poly_at_fraction": ("poly_at_flag", "mean"),
                f"{sample_prefix}_{side}_poly_at_max_run": ("poly_at_max_run", "max"),
            }
        )
        .merge(family_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_family"]], on=["chrom", "window_start", "window_end"], how="left")
        .merge(subfamily_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_subfamily"]], on=["chrom", "window_start", "window_end"], how="left")
        .merge(strand_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_strand"]], on=["chrom", "window_start", "window_end"], how="left")
        .merge(
            purity[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_{side}_mei_subfamily_purity",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            pos_mode[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_{side}_mei_breakpoint_mode",
                    "mode_support_reads",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            pos_unique,
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            coord_agg,
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            overlap_tbl,
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
    )
    agg[f"{sample_prefix}_{side}_mei_breakpoint_mode_fraction"] = (
        agg["mode_support_reads"] / agg[f"{sample_prefix}_{side}_mei_supported_reads"]
    ).fillna(0.0)
    agg[f"{sample_prefix}_{side}_mei_start"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_mei_start"],
        errors="coerce",
    ).fillna(0).astype(int)
    agg[f"{sample_prefix}_{side}_mei_end"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_mei_end"],
        errors="coerce",
    ).fillna(0).astype(int)
    agg[f"{sample_prefix}_{side}_mei_anchor_bp_max"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_mei_anchor_bp_max"],
        errors="coerce",
    ).fillna(0).astype(int)
    agg[f"{sample_prefix}_{side}_mei_target_len"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_mei_target_len"],
        errors="coerce",
    ).fillna(0).astype(int)
    agg[f"{sample_prefix}_{side}_clip_overlap_informative_reads"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_clip_overlap_informative_reads"],
        errors="coerce",
    ).fillna(0).astype(int)
    agg[f"{sample_prefix}_{side}_clip_overlap_jaccard_median"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_clip_overlap_jaccard_median"],
        errors="coerce",
    ).fillna(0.0).astype(float)
    agg[f"{sample_prefix}_{side}_clip_overlap_pair_ge20_fraction"] = pd.to_numeric(
        agg[f"{sample_prefix}_{side}_clip_overlap_pair_ge20_fraction"],
        errors="coerce",
    ).fillna(0.0).astype(float)
    agg = agg.drop(columns=["mode_support_reads"])
    return agg


def _aggregate_discordant_mei_metrics(df: pd.DataFrame, sample_prefix: str) -> pd.DataFrame:
    empty_cols = [
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_discordant_mei_supported_reads",
                f"{sample_prefix}_discordant_mei_score_sum",
                f"{sample_prefix}_discordant_mei_family",
                f"{sample_prefix}_discordant_mei_subfamily",
                f"{sample_prefix}_discordant_mei_family_votes",
                f"{sample_prefix}_discordant_mei_subfamily_votes",
                f"{sample_prefix}_discordant_mei_strand",
                f"{sample_prefix}_discordant_mei_family_purity",
                f"{sample_prefix}_discordant_mei_strand_purity",
                f"{sample_prefix}_discordant_mei_left_supported_reads",
                f"{sample_prefix}_discordant_mei_right_supported_reads",
                f"{sample_prefix}_discordant_mei_left_target_pos_median",
                f"{sample_prefix}_discordant_mei_right_target_pos_median",
                f"{sample_prefix}_discordant_mei_insertion_span_estimate",
                f"{sample_prefix}_discordant_mei_orientation_order_consistent",
                f"{sample_prefix}_discordant_mei_geometry_consistent",
                f"{sample_prefix}_discordant_mei_left_subfamily",
                f"{sample_prefix}_discordant_mei_right_subfamily",
                f"{sample_prefix}_discordant_mei_side_subfamily_consistent",
                f"{sample_prefix}_discordant_mei_left_anchor_bin_mode_fraction",
                f"{sample_prefix}_discordant_mei_right_anchor_bin_mode_fraction",
                f"{sample_prefix}_discordant_mei_left_target_bin_mode_fraction",
                f"{sample_prefix}_discordant_mei_right_target_bin_mode_fraction",
                f"{sample_prefix}_discordant_mei_left_side_coherence",
                f"{sample_prefix}_discordant_mei_right_side_coherence",
                f"{sample_prefix}_discordant_mei_side_coherence_min",
                f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
                f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
                f"{sample_prefix}_discordant_mei_anchor_target_spearman_abs_min",
                f"{sample_prefix}_discordant_mei_left_local_jump_violation",
                f"{sample_prefix}_discordant_mei_right_local_jump_violation",
                f"{sample_prefix}_discordant_mei_any_local_jump_violation",
                f"{sample_prefix}_discordant_mei_insert_sd_proxy",
                f"{sample_prefix}_discordant_mei_max_pair_swing",
                f"{sample_prefix}_discordant_mei_self_consistent",
    ]
    if df is None or df.empty or "mei_hit" not in df.columns:
        return pd.DataFrame(columns=empty_cols)
    mei_df = df.loc[df["mei_hit"]].copy()
    if mei_df.empty:
        return pd.DataFrame(columns=empty_cols)

    mei_df["locus_midpoint"] = (mei_df["window_start"].astype(int) + mei_df["window_end"].astype(int)) // 2
    mei_df["anchor_side"] = mei_df.apply(
        lambda r: "L" if int(r["pos"]) <= int(r["locus_midpoint"]) else "R",
        axis=1,
    )
    mei_df["anchor_bin_10bp"] = (mei_df["pos"].astype(int) // 10).astype(int)

    # Family/subfamily identity: only mates that are interchromosomal or >1 kb away.
    # Prefer mate consensus labels when present. Support counts / geometry below
    # still use the full mei_df.
    identity_df = _discordant_rows_for_mei_identity(mei_df)
    family_top, subfamily_top = _top_family_then_subfamily(
        identity_df,
        group_cols=["chrom", "window_start", "window_end"],
        family_col="family",
        subfamily_col="target",
        score_col="mei_score",
    )
    family_top = family_top.rename(columns={"family": f"{sample_prefix}_discordant_mei_family"})
    subfamily_top = subfamily_top.rename(columns={"target": f"{sample_prefix}_discordant_mei_subfamily"})
    identity_votes = _identity_vote_maps(identity_df, sample_prefix=sample_prefix)
    strand_top = (
        mei_df.groupby(["chrom", "window_start", "window_end", "target_strand"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target_strand": f"{sample_prefix}_discordant_mei_strand"})
    )
    # Purity for the identity-eligible subset (same rows that chose family/subfamily).
    purity_src = identity_df if not identity_df.empty else mei_df.iloc[0:0].copy()
    family_sum = (
        purity_src.groupby(["chrom", "window_start", "window_end"], as_index=False)["mei_score"]
        .sum()
        .rename(columns={"mei_score": "all_family_score_sum"})
    )
    family_purity = family_top.rename(columns={"mei_score": "top_family_score_sum"}).merge(
        family_sum,
        on=["chrom", "window_start", "window_end"],
        how="left",
    )
    family_purity[f"{sample_prefix}_discordant_mei_family_purity"] = (
        family_purity["top_family_score_sum"] / family_purity["all_family_score_sum"]
    ).fillna(0.0)

    strand_sum = (
        mei_df.groupby(["chrom", "window_start", "window_end"], as_index=False)["mei_score"]
        .sum()
        .rename(columns={"mei_score": "all_strand_score_sum"})
    )
    strand_purity = strand_top.rename(columns={"mei_score": "top_strand_score_sum"}).merge(
        strand_sum,
        on=["chrom", "window_start", "window_end"],
        how="left",
    )
    strand_purity[f"{sample_prefix}_discordant_mei_strand_purity"] = (
        strand_purity["top_strand_score_sum"] / strand_purity["all_strand_score_sum"]
    ).fillna(0.0)

    side_counts = (
        mei_df.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": "side_unique_reads"})
    )
    side_pivot = (
        side_counts.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values="side_unique_reads",
            fill_value=0,
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    if "L" not in side_pivot.columns:
        side_pivot["L"] = 0
    if "R" not in side_pivot.columns:
        side_pivot["R"] = 0
    side_pivot = side_pivot.rename(
        columns={
            "L": f"{sample_prefix}_discordant_mei_left_supported_reads",
            "R": f"{sample_prefix}_discordant_mei_right_supported_reads",
        }
    )

    mei_df["target_mid"] = ((mei_df["target_start"].astype(int) + mei_df["target_end"].astype(int)) // 2).astype(int)
    mei_df["target_bin_25bp"] = (mei_df["target_mid"].astype(int) // 25).astype(int)
    side_target_mid = (
        mei_df.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)["target_mid"]
        .median()
        .rename(columns={"target_mid": "target_mid_median"})
    )
    side_target_extent = (
        mei_df.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)
        .agg(
            target_start_min=("target_start", "min"),
            target_end_max=("target_end", "max"),
        )
    )
    side_mid_pivot = (
        side_target_mid.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values="target_mid_median",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    if "L" not in side_mid_pivot.columns:
        side_mid_pivot["L"] = 0
    if "R" not in side_mid_pivot.columns:
        side_mid_pivot["R"] = 0
    side_mid_pivot = side_mid_pivot.rename(
        columns={
            "L": f"{sample_prefix}_discordant_mei_left_target_pos_median",
            "R": f"{sample_prefix}_discordant_mei_right_target_pos_median",
        }
    )
    side_extent_start = (
        side_target_extent.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values="target_start_min",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    side_extent_end = (
        side_target_extent.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values="target_end_max",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for side_col in ("L", "R"):
        if side_col not in side_extent_start.columns:
            side_extent_start[side_col] = 0
        if side_col not in side_extent_end.columns:
            side_extent_end[side_col] = 0
    side_extent_start = side_extent_start.rename(
        columns={
            "L": f"{sample_prefix}_discordant_mei_left_target_start_min",
            "R": f"{sample_prefix}_discordant_mei_right_target_start_min",
        }
    )
    side_extent_end = side_extent_end.rename(
        columns={
            "L": f"{sample_prefix}_discordant_mei_left_target_end_max",
            "R": f"{sample_prefix}_discordant_mei_right_target_end_max",
        }
    )

    side_subfamily_src = identity_df if not identity_df.empty else mei_df.iloc[0:0].copy()
    if not side_subfamily_src.empty:
        # Constrain per-side subfamily to the locus winning family when available.
        if not family_top.empty:
            side_subfamily_src = side_subfamily_src.merge(
                family_top[["chrom", "window_start", "window_end", f"{sample_prefix}_discordant_mei_family"]],
                on=["chrom", "window_start", "window_end"],
                how="inner",
            )
            side_subfamily_src = side_subfamily_src.loc[
                side_subfamily_src["family"].astype(str)
                == side_subfamily_src[f"{sample_prefix}_discordant_mei_family"].astype(str)
            ]
        side_subfamily_top = (
            side_subfamily_src.groupby(
                ["chrom", "window_start", "window_end", "anchor_side", "target"],
                as_index=False,
            )["mei_score"]
            .sum()
            .sort_values(
                ["chrom", "window_start", "window_end", "anchor_side", "mei_score"],
                ascending=[True, True, True, True, False],
            )
            .drop_duplicates(["chrom", "window_start", "window_end", "anchor_side"], keep="first")
            .rename(columns={"target": "side_top_subfamily"})
        )
    else:
        side_subfamily_top = pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "anchor_side", "side_top_subfamily"]
        )
    side_subfamily_pivot = (
        side_subfamily_top.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values="side_top_subfamily",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    if "L" not in side_subfamily_pivot.columns:
        side_subfamily_pivot["L"] = ""
    if "R" not in side_subfamily_pivot.columns:
        side_subfamily_pivot["R"] = ""
    side_subfamily_pivot = side_subfamily_pivot.rename(
        columns={
            "L": f"{sample_prefix}_discordant_mei_left_subfamily",
            "R": f"{sample_prefix}_discordant_mei_right_subfamily",
        }
    )

    def _side_bin_mode_fraction(
        data: pd.DataFrame,
        key_col: str,
        out_col: str,
    ) -> pd.DataFrame:
        counts = (
            data.groupby(["chrom", "window_start", "window_end", "anchor_side", key_col], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "n_reads"})
        )
        top = (
            counts.sort_values(
                ["chrom", "window_start", "window_end", "anchor_side", "n_reads"],
                ascending=[True, True, True, True, False],
            )
            .drop_duplicates(["chrom", "window_start", "window_end", "anchor_side"], keep="first")
            .rename(columns={"n_reads": "n_reads_mode"})
        )
        totals = (
            data.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "n_reads_total"})
        )
        merged = top.merge(
            totals,
            on=["chrom", "window_start", "window_end", "anchor_side"],
            how="inner",
        )
        merged[out_col] = (merged["n_reads_mode"] / merged["n_reads_total"]).fillna(0.0).astype(float)
        return merged[["chrom", "window_start", "window_end", "anchor_side", out_col]]

    side_anchor_mode_frac = _side_bin_mode_fraction(
        data=mei_df,
        key_col="anchor_bin_10bp",
        out_col="side_anchor_bin_mode_fraction",
    )
    side_target_mode_frac = _side_bin_mode_fraction(
        data=mei_df,
        key_col="target_bin_25bp",
        out_col="side_target_bin_mode_fraction",
    )
    side_mode_frac = side_anchor_mode_frac.merge(
        side_target_mode_frac,
        on=["chrom", "window_start", "window_end", "anchor_side"],
        how="outer",
    ).fillna(0.0)
    side_mode_frac["side_coherence"] = side_mode_frac[
        ["side_anchor_bin_mode_fraction", "side_target_bin_mode_fraction"]
    ].min(axis=1)
    # Side-wise monotonicity between genomic anchor position and MEI target position.
    # This captures insert-size-driven spread (including inverse ordering) better than
    # strict local bin concentration alone.
    side_spearman = (
        mei_df.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)
        .apply(
            lambda g: pd.Series(
                {
                    "side_anchor_target_spearman_abs": abs(
                        float(
                            g.loc[:, ["pos", "target_mid"]]
                            .corr(method="spearman")
                            .iloc[0, 1]
                        )
                    )
                    if len(g) >= 3
                    else 1.0
                }
            ),
            include_groups=False,
        )
        .reset_index(drop=True)
    )
    side_mode_frac = side_mode_frac.merge(
        side_spearman,
        on=["chrom", "window_start", "window_end", "anchor_side"],
        how="left",
    )
    side_mode_frac["side_anchor_target_spearman_abs"] = (
        side_mode_frac["side_anchor_target_spearman_abs"].fillna(0.0).astype(float)
    )
    side_mode_pivot = (
        side_mode_frac.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values=[
                "side_anchor_bin_mode_fraction",
                "side_target_bin_mode_fraction",
                "side_coherence",
                "side_anchor_target_spearman_abs",
            ],
            aggfunc="first",
        )
        .reset_index()
    )
    side_mode_pivot.columns = [
        (
            col
            if isinstance(col, str)
            else col[0]
            if len(col) > 1 and col[1] in {"", None}
            else f"{col[0]}_{col[1]}"
        )
        for col in side_mode_pivot.columns
    ]
    side_mode_pivot = side_mode_pivot.rename(
        columns={
            "side_anchor_bin_mode_fraction_L": f"{sample_prefix}_discordant_mei_left_anchor_bin_mode_fraction",
            "side_anchor_bin_mode_fraction_R": f"{sample_prefix}_discordant_mei_right_anchor_bin_mode_fraction",
            "side_target_bin_mode_fraction_L": f"{sample_prefix}_discordant_mei_left_target_bin_mode_fraction",
            "side_target_bin_mode_fraction_R": f"{sample_prefix}_discordant_mei_right_target_bin_mode_fraction",
            "side_coherence_L": f"{sample_prefix}_discordant_mei_left_side_coherence",
            "side_coherence_R": f"{sample_prefix}_discordant_mei_right_side_coherence",
            "side_anchor_target_spearman_abs_L": f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
            "side_anchor_target_spearman_abs_R": f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
        }
    )
    for col in [
        f"{sample_prefix}_discordant_mei_left_anchor_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_right_anchor_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_left_target_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_right_target_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_left_side_coherence",
        f"{sample_prefix}_discordant_mei_right_side_coherence",
        f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
        f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
    ]:
        if col not in side_mode_pivot.columns:
            side_mode_pivot[col] = 0.0

    # Proxy for expected fragment-size variation. Prefer sample-observed spread;
    # clamp to a sensible lower bound so very small estimates do not over-penalize.
    insert_sd_proxy = float(mei_df["template_len"].abs().astype(float).std(ddof=0))
    if not math.isfinite(insert_sd_proxy):
        insert_sd_proxy = 100.0
    insert_sd_proxy = max(50.0, insert_sd_proxy)
    swing_sigma_cutoff = 3.0 * insert_sd_proxy

    def _side_local_jump_violation(data: pd.DataFrame) -> pd.DataFrame:
        # Side-internal mapping incoherence test based on relative-position swing.
        #
        # For adjacent reads sorted by genomic anchor:
        #   d_anchor = pos_i - pos_{i-1}
        #   d_target = target_i - target_{i-1}
        #
        # Expected consistent behavior can look direct (d_target ~= d_anchor)
        # or inverse (d_target ~= -d_anchor), depending on orientation/mapping frame.
        # We therefore use:
        #   swing = min(|d_target - d_anchor|, |d_target + d_anchor|)
        #
        # Flag violation if any adjacent pair exceeds 3*insert_sd_proxy.
        rows: list[dict[str, object]] = []
        key_cols = ["chrom", "window_start", "window_end", "anchor_side"]
        for key, g in data.groupby(key_cols, sort=False):
            gg = g.loc[:, ["pos", "target_mid"]].copy()
            gg["pos"] = gg["pos"].astype(int)
            gg["target_mid"] = gg["target_mid"].astype(int)
            gg = gg.sort_values("pos", kind="mergesort").drop_duplicates()
            violated = False
            max_pair_swing = 0.0
            if len(gg) >= 2:
                pos_vals = gg["pos"].tolist()
                tgt_vals = gg["target_mid"].tolist()
                for i in range(1, len(pos_vals)):
                    d_anchor_signed = float(int(pos_vals[i]) - int(pos_vals[i - 1]))
                    d_target_signed = float(int(tgt_vals[i]) - int(tgt_vals[i - 1]))
                    swing_direct = abs(d_target_signed - d_anchor_signed)
                    swing_inverse = abs(d_target_signed + d_anchor_signed)
                    pair_swing = min(swing_direct, swing_inverse)
                    if pair_swing > max_pair_swing:
                        max_pair_swing = float(pair_swing)
                    if pair_swing > swing_sigma_cutoff:
                        violated = True
                        break
            rows.append(
                {
                    "chrom": key[0],
                    "window_start": key[1],
                    "window_end": key[2],
                    "anchor_side": key[3],
                    "side_local_jump_violation": bool(violated),
                    "side_max_pair_swing": float(max_pair_swing),
                }
            )
        if not rows:
            return pd.DataFrame(columns=key_cols + ["side_local_jump_violation", "side_max_pair_swing"])
        return pd.DataFrame(rows)

    side_jump_violation = _side_local_jump_violation(mei_df)
    side_jump_pivot = (
        side_jump_violation.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values=["side_local_jump_violation", "side_max_pair_swing"],
            aggfunc="first",
        )
        .reset_index()
    )
    side_jump_pivot.columns = [
        (
            col
            if isinstance(col, str)
            else col[0]
            if len(col) > 1 and col[1] in {"", None}
            else f"{col[0]}_{col[1]}"
        )
        for col in side_jump_pivot.columns
    ]
    side_jump_pivot = side_jump_pivot.rename(
        columns={
            "side_local_jump_violation_L": f"{sample_prefix}_discordant_mei_left_local_jump_violation",
            "side_local_jump_violation_R": f"{sample_prefix}_discordant_mei_right_local_jump_violation",
            "side_max_pair_swing_L": f"{sample_prefix}_discordant_mei_left_max_pair_swing",
            "side_max_pair_swing_R": f"{sample_prefix}_discordant_mei_right_max_pair_swing",
        }
    )
    for col, default in [
        (f"{sample_prefix}_discordant_mei_left_local_jump_violation", False),
        (f"{sample_prefix}_discordant_mei_right_local_jump_violation", False),
        (f"{sample_prefix}_discordant_mei_left_max_pair_swing", 0.0),
        (f"{sample_prefix}_discordant_mei_right_max_pair_swing", 0.0),
    ]:
        if col not in side_jump_pivot.columns:
            side_jump_pivot[col] = default

    agg = (
        mei_df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_discordant_mei_supported_reads": ("read_name", "nunique"),
                f"{sample_prefix}_discordant_mei_score_sum": ("mei_score", "sum"),
            }
        )
        .merge(
            family_top[["chrom", "window_start", "window_end", f"{sample_prefix}_discordant_mei_family"]],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            subfamily_top[["chrom", "window_start", "window_end", f"{sample_prefix}_discordant_mei_subfamily"]],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            identity_votes,
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            strand_top[["chrom", "window_start", "window_end", f"{sample_prefix}_discordant_mei_strand"]],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            family_purity[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_family_purity",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            strand_purity[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_strand_purity",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_pivot[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_supported_reads",
                    f"{sample_prefix}_discordant_mei_right_supported_reads",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_mid_pivot[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_target_pos_median",
                    f"{sample_prefix}_discordant_mei_right_target_pos_median",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_extent_start[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_target_start_min",
                    f"{sample_prefix}_discordant_mei_right_target_start_min",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_extent_end[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_target_end_max",
                    f"{sample_prefix}_discordant_mei_right_target_end_max",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_subfamily_pivot[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_subfamily",
                    f"{sample_prefix}_discordant_mei_right_subfamily",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_mode_pivot[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_anchor_bin_mode_fraction",
                    f"{sample_prefix}_discordant_mei_right_anchor_bin_mode_fraction",
                    f"{sample_prefix}_discordant_mei_left_target_bin_mode_fraction",
                    f"{sample_prefix}_discordant_mei_right_target_bin_mode_fraction",
                    f"{sample_prefix}_discordant_mei_left_side_coherence",
                    f"{sample_prefix}_discordant_mei_right_side_coherence",
                    f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
                    f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
        .merge(
            side_jump_pivot[
                [
                    "chrom",
                    "window_start",
                    "window_end",
                    f"{sample_prefix}_discordant_mei_left_local_jump_violation",
                    f"{sample_prefix}_discordant_mei_right_local_jump_violation",
                    f"{sample_prefix}_discordant_mei_left_max_pair_swing",
                    f"{sample_prefix}_discordant_mei_right_max_pair_swing",
                ]
            ],
            on=["chrom", "window_start", "window_end"],
            how="left",
        )
    )
    agg[f"{sample_prefix}_discordant_mei_left_supported_reads"] = (
        agg[f"{sample_prefix}_discordant_mei_left_supported_reads"].fillna(0).astype(int)
    )
    agg[f"{sample_prefix}_discordant_mei_right_supported_reads"] = (
        agg[f"{sample_prefix}_discordant_mei_right_supported_reads"].fillna(0).astype(int)
    )
    agg[f"{sample_prefix}_discordant_mei_family_purity"] = (
        agg[f"{sample_prefix}_discordant_mei_family_purity"].fillna(0.0).astype(float)
    )
    agg[f"{sample_prefix}_discordant_mei_strand_purity"] = (
        agg[f"{sample_prefix}_discordant_mei_strand_purity"].fillna(0.0).astype(float)
    )

    agg[f"{sample_prefix}_discordant_mei_left_target_pos_median"] = (
        agg[f"{sample_prefix}_discordant_mei_left_target_pos_median"].fillna(0).astype(float)
    )
    agg[f"{sample_prefix}_discordant_mei_right_target_pos_median"] = (
        agg[f"{sample_prefix}_discordant_mei_right_target_pos_median"].fillna(0).astype(float)
    )
    for extent_col in (
        f"{sample_prefix}_discordant_mei_left_target_start_min",
        f"{sample_prefix}_discordant_mei_right_target_start_min",
        f"{sample_prefix}_discordant_mei_left_target_end_max",
        f"{sample_prefix}_discordant_mei_right_target_end_max",
    ):
        agg[extent_col] = pd.to_numeric(agg.get(extent_col, 0), errors="coerce").fillna(0).astype(int)
    agg[f"{sample_prefix}_discordant_mei_left_subfamily"] = (
        agg[f"{sample_prefix}_discordant_mei_left_subfamily"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_right_subfamily"] = (
        agg[f"{sample_prefix}_discordant_mei_right_subfamily"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_family"] = (
        agg[f"{sample_prefix}_discordant_mei_family"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_subfamily"] = (
        agg[f"{sample_prefix}_discordant_mei_subfamily"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_family_votes"] = (
        agg.get(f"{sample_prefix}_discordant_mei_family_votes", pd.Series("", index=agg.index))
        .fillna("")
        .astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_subfamily_votes"] = (
        agg.get(f"{sample_prefix}_discordant_mei_subfamily_votes", pd.Series("", index=agg.index))
        .fillna("")
        .astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_strand"] = (
        agg[f"{sample_prefix}_discordant_mei_strand"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_side_subfamily_consistent"] = (
        (agg[f"{sample_prefix}_discordant_mei_left_subfamily"] != "")
        & (agg[f"{sample_prefix}_discordant_mei_right_subfamily"] != "")
        & (agg[f"{sample_prefix}_discordant_mei_left_subfamily"] == agg[f"{sample_prefix}_discordant_mei_right_subfamily"])
    )
    for col in [
        f"{sample_prefix}_discordant_mei_left_anchor_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_right_anchor_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_left_target_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_right_target_bin_mode_fraction",
        f"{sample_prefix}_discordant_mei_left_side_coherence",
        f"{sample_prefix}_discordant_mei_right_side_coherence",
        f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
        f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
    ]:
        agg[col] = agg[col].fillna(0.0).astype(float)
    agg[f"{sample_prefix}_discordant_mei_side_coherence_min"] = agg[
        [
            f"{sample_prefix}_discordant_mei_left_side_coherence",
            f"{sample_prefix}_discordant_mei_right_side_coherence",
        ]
    ].min(axis=1)
    agg[f"{sample_prefix}_discordant_mei_anchor_target_spearman_abs_min"] = agg[
        [
            f"{sample_prefix}_discordant_mei_left_anchor_target_spearman_abs",
            f"{sample_prefix}_discordant_mei_right_anchor_target_spearman_abs",
        ]
    ].min(axis=1)
    agg[f"{sample_prefix}_discordant_mei_left_local_jump_violation"] = (
        agg[f"{sample_prefix}_discordant_mei_left_local_jump_violation"].fillna(False).astype(bool)
    )
    agg[f"{sample_prefix}_discordant_mei_right_local_jump_violation"] = (
        agg[f"{sample_prefix}_discordant_mei_right_local_jump_violation"].fillna(False).astype(bool)
    )
    agg[f"{sample_prefix}_discordant_mei_left_max_pair_swing"] = (
        agg[f"{sample_prefix}_discordant_mei_left_max_pair_swing"].fillna(0.0).astype(float)
    )
    agg[f"{sample_prefix}_discordant_mei_right_max_pair_swing"] = (
        agg[f"{sample_prefix}_discordant_mei_right_max_pair_swing"].fillna(0.0).astype(float)
    )
    agg[f"{sample_prefix}_discordant_mei_insert_sd_proxy"] = float(insert_sd_proxy)
    agg[f"{sample_prefix}_discordant_mei_max_pair_swing"] = agg[
        [
            f"{sample_prefix}_discordant_mei_left_max_pair_swing",
            f"{sample_prefix}_discordant_mei_right_max_pair_swing",
        ]
    ].max(axis=1)
    agg[f"{sample_prefix}_discordant_mei_any_local_jump_violation"] = (
        agg[f"{sample_prefix}_discordant_mei_left_local_jump_violation"]
        | agg[f"{sample_prefix}_discordant_mei_right_local_jump_violation"]
    )
    left_mid = agg[f"{sample_prefix}_discordant_mei_left_target_pos_median"]
    right_mid = agg[f"{sample_prefix}_discordant_mei_right_target_pos_median"]
    span = (right_mid - left_mid).abs() + 1.0
    agg[f"{sample_prefix}_discordant_mei_insertion_span_estimate"] = span
    dominant_strand = agg[f"{sample_prefix}_discordant_mei_strand"].fillna("").astype(str)
    order_consistent = (
        ((dominant_strand == "+") & (right_mid >= left_mid))
        | ((dominant_strand == "-") & (left_mid >= right_mid))
    )
    agg[f"{sample_prefix}_discordant_mei_orientation_order_consistent"] = order_consistent
    # Geometry-consistent DPE insertion footprint:
    # - bilateral support
    # - expected left/right order by orientation
    # - plausible insertion span on consensus (exclude tiny/noise and huge artifacts)
    agg[f"{sample_prefix}_discordant_mei_geometry_consistent"] = (
        (agg[f"{sample_prefix}_discordant_mei_left_supported_reads"] >= 1)
        & (agg[f"{sample_prefix}_discordant_mei_right_supported_reads"] >= 1)
        & agg[f"{sample_prefix}_discordant_mei_orientation_order_consistent"]
        & (agg[f"{sample_prefix}_discordant_mei_insertion_span_estimate"] >= 30.0)
        & (agg[f"{sample_prefix}_discordant_mei_insertion_span_estimate"] <= 8000.0)
    )
    side_reads_min = agg[
        [
            f"{sample_prefix}_discordant_mei_left_supported_reads",
            f"{sample_prefix}_discordant_mei_right_supported_reads",
        ]
    ].min(axis=1)
    # Position coherence on each side:
    # - concentration in local bins (good for tight clusters), OR
    # - strong monotonic anchor<->target relationship (good for broader insert-size spread,
    #   including inverse ordering).
    # For low-support sides (<3 reads), avoid over-penalizing by treating coherence as pass.
    agg[f"{sample_prefix}_discordant_mei_self_consistent"] = (
        (~agg[f"{sample_prefix}_discordant_mei_any_local_jump_violation"])
        & (
            (side_reads_min < 3)
            | (agg[f"{sample_prefix}_discordant_mei_side_coherence_min"] >= 0.5)
            | (agg[f"{sample_prefix}_discordant_mei_anchor_target_spearman_abs_min"] >= 0.6)
        )
    )
    return agg


def _aggregate_discordant_anchor_side_metrics(df: pd.DataFrame, sample_prefix: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_discordant_anchor_left_unique_reads",
                f"{sample_prefix}_discordant_anchor_right_unique_reads",
                f"{sample_prefix}_discordant_anchor_left_complex_reason_max_fraction",
                f"{sample_prefix}_discordant_anchor_right_complex_reason_max_fraction",
                f"{sample_prefix}_discordant_anchor_left_mapq_mean",
                f"{sample_prefix}_discordant_anchor_right_mapq_mean",
            ]
        )

    tmp = df.copy()
    tmp["locus_midpoint"] = (tmp["window_start"].astype(int) + tmp["window_end"].astype(int)) // 2
    tmp["anchor_side"] = tmp.apply(
        lambda r: "L" if int(r["pos"]) <= int(r["locus_midpoint"]) else "R",
        axis=1,
    )
    # Complex sidepair evidence should ignore MEI-mapped / polyA/VNTR-rescued
    # discordants (expected for simple insertions).
    mei_related = _discordant_row_mei_related(tmp)
    if len(mei_related) == len(tmp) and bool(mei_related.any()):
        complex_src = tmp.loc[~mei_related.astype(bool)].copy()
    else:
        complex_src = tmp
    reason_col = complex_src["discordant_reasons"].fillna("").astype(str)
    complex_src["reason_interchrom"] = reason_col.str.contains("interchrom", regex=False).astype(int)
    complex_src["reason_mate_unmapped"] = reason_col.str.contains("mate_unmapped", regex=False).astype(int)
    complex_src["reason_large_insert"] = reason_col.str.contains("large_insert", regex=False).astype(int)
    complex_src["reason_same_strand"] = reason_col.str.contains("same_strand", regex=False).astype(int)
    complex_src["reason_improper_pair"] = reason_col.str.contains("improper_pair", regex=False).astype(int)
    complex_src["reason_complex_any"] = (
        (complex_src["reason_interchrom"] == 1)
        | (complex_src["reason_mate_unmapped"] == 1)
        | (complex_src["reason_large_insert"] == 1)
        | (complex_src["reason_same_strand"] == 1)
        | (complex_src["reason_improper_pair"] == 1)
    ).astype(int)
    complex_src["mapq"] = complex_src["mapq"].astype(float)

    if complex_src.empty:
        side = pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                "anchor_side",
                "side_unique_reads",
                "side_mapq_mean",
                "side_interchrom_fraction",
                "side_mate_unmapped_fraction",
                "side_large_insert_fraction",
                "side_same_strand_fraction",
                "side_improper_pair_fraction",
                "side_complex_any_fraction",
                "side_complex_reason_max_fraction",
            ]
        )
    else:
        side = (
            complex_src.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)
            .agg(
                side_unique_reads=("read_name", "nunique"),
                side_mapq_mean=("mapq", "mean"),
                side_interchrom_fraction=("reason_interchrom", "mean"),
                side_mate_unmapped_fraction=("reason_mate_unmapped", "mean"),
                side_large_insert_fraction=("reason_large_insert", "mean"),
                side_same_strand_fraction=("reason_same_strand", "mean"),
                side_improper_pair_fraction=("reason_improper_pair", "mean"),
                side_complex_any_fraction=("reason_complex_any", "mean"),
            )
            .sort_values(["chrom", "window_start", "window_end", "anchor_side"], kind="mergesort")
        )
        side["side_complex_reason_max_fraction"] = side[
            [
                "side_interchrom_fraction",
                "side_mate_unmapped_fraction",
                "side_large_insert_fraction",
                "side_same_strand_fraction",
                "side_improper_pair_fraction",
                "side_complex_any_fraction",
            ]
        ].max(axis=1)

    if side.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_discordant_anchor_left_unique_reads",
                f"{sample_prefix}_discordant_anchor_right_unique_reads",
                f"{sample_prefix}_discordant_anchor_left_complex_reason_max_fraction",
                f"{sample_prefix}_discordant_anchor_right_complex_reason_max_fraction",
                f"{sample_prefix}_discordant_anchor_left_mapq_mean",
                f"{sample_prefix}_discordant_anchor_right_mapq_mean",
            ]
        )

    pivot = (
        side.pivot_table(
            index=["chrom", "window_start", "window_end"],
            columns="anchor_side",
            values=["side_unique_reads", "side_complex_reason_max_fraction", "side_mapq_mean"],
            aggfunc="first",
        )
        .reset_index()
    )
    pivot.columns = [
        (
            col
            if isinstance(col, str)
            else col[0]
            if len(col) > 1 and col[1] in {"", None}
            else f"{col[0]}_{col[1]}"
        )
        for col in pivot.columns
    ]
    pivot = pivot.rename(
        columns={
            "side_unique_reads_L": f"{sample_prefix}_discordant_anchor_left_unique_reads",
            "side_unique_reads_R": f"{sample_prefix}_discordant_anchor_right_unique_reads",
            "side_complex_reason_max_fraction_L": f"{sample_prefix}_discordant_anchor_left_complex_reason_max_fraction",
            "side_complex_reason_max_fraction_R": f"{sample_prefix}_discordant_anchor_right_complex_reason_max_fraction",
            "side_mapq_mean_L": f"{sample_prefix}_discordant_anchor_left_mapq_mean",
            "side_mapq_mean_R": f"{sample_prefix}_discordant_anchor_right_mapq_mean",
        }
    )

    defaults: list[tuple[str, float | int]] = [
        (f"{sample_prefix}_discordant_anchor_left_unique_reads", 0),
        (f"{sample_prefix}_discordant_anchor_right_unique_reads", 0),
        (f"{sample_prefix}_discordant_anchor_left_complex_reason_max_fraction", 0.0),
        (f"{sample_prefix}_discordant_anchor_right_complex_reason_max_fraction", 0.0),
        (f"{sample_prefix}_discordant_anchor_left_mapq_mean", 0.0),
        (f"{sample_prefix}_discordant_anchor_right_mapq_mean", 0.0),
    ]
    for col, default in defaults:
        if col not in pivot.columns:
            pivot[col] = default
    pivot[f"{sample_prefix}_discordant_anchor_left_unique_reads"] = (
        pivot[f"{sample_prefix}_discordant_anchor_left_unique_reads"].fillna(0).astype(int)
    )
    pivot[f"{sample_prefix}_discordant_anchor_right_unique_reads"] = (
        pivot[f"{sample_prefix}_discordant_anchor_right_unique_reads"].fillna(0).astype(int)
    )
    for col in [
        f"{sample_prefix}_discordant_anchor_left_complex_reason_max_fraction",
        f"{sample_prefix}_discordant_anchor_right_complex_reason_max_fraction",
        f"{sample_prefix}_discordant_anchor_left_mapq_mean",
        f"{sample_prefix}_discordant_anchor_right_mapq_mean",
    ]:
        pivot[col] = pivot[col].fillna(0.0).astype(float)
    return pivot


def _infer_disease_insertion_metrics(
    candidates: pd.DataFrame,
    reference_fasta: Path | None = None,
    split_disease: pd.DataFrame | None = None,
    split_control: pd.DataFrame | None = None,
) -> pd.DataFrame:
    out = candidates.copy()
    for col in [
        "disease_L_mei_start",
        "disease_R_mei_start",
        "disease_L_mei_end",
        "disease_R_mei_end",
        "disease_L_mei_breakpoint_mode",
        "disease_R_mei_breakpoint_mode",
        "control_L_mei_breakpoint_mode",
        "control_R_mei_breakpoint_mode",
        "disease_L_mei_supported_reads",
        "disease_R_mei_supported_reads",
        "control_L_mei_supported_reads",
        "control_R_mei_supported_reads",
    ]:
        if col not in out.columns:
            out[col] = 0
        out[col] = out[col].fillna(0).astype(int)

    disease_metrics = out.apply(
        lambda r: _sample_insertion_span_and_orientation(r, "disease"),
        axis=1,
        result_type="expand",
    )
    disease_metrics.columns = [
        "disease_insertion_mei_start",
        "disease_insertion_mei_end",
        "disease_insertion_mei_span",
        "disease_insertion_orientation",
    ]
    for col in disease_metrics.columns:
        out[col] = disease_metrics[col]

    control_metrics = out.apply(
        lambda r: _sample_insertion_span_and_orientation(r, "control"),
        axis=1,
        result_type="expand",
    )
    control_metrics.columns = [
        "control_insertion_mei_start",
        "control_insertion_mei_end",
        "control_insertion_mei_span",
        "control_insertion_orientation",
    ]
    for col in control_metrics.columns:
        out[col] = control_metrics[col]

    def _pick_tsd_pair(row: pd.Series) -> tuple[int, int, int, str]:
        # Prefer strict bilateral pairs from either disease or control.
        # If strict length is invalid, allow a small ±2 bp breakpoint rescue.
        candidates: list[tuple[int, int, int, str]] = []
        t_l = int(row.get("disease_L_mei_breakpoint_mode", 0))
        t_r = int(row.get("disease_R_mei_breakpoint_mode", 0))
        t_support = int(row.get("disease_L_mei_supported_reads", 0)) + int(row.get("disease_R_mei_supported_reads", 0))
        if t_l > 0 and t_r > 0:
            candidates.append((t_l, t_r, t_support, "tsd_disease"))
        n_l = int(row.get("control_L_mei_breakpoint_mode", 0))
        n_r = int(row.get("control_R_mei_breakpoint_mode", 0))
        n_support = int(row.get("control_L_mei_supported_reads", 0)) + int(row.get("control_R_mei_supported_reads", 0))
        if n_l > 0 and n_r > 0:
            candidates.append((n_l, n_r, n_support, "tsd_control"))
        if not candidates:
            return (0, 0, 0, "")

        # Try strict first (no coordinate adjustment).
        strict_ok: list[tuple[int, int, int, str]] = []
        for l, r, support, source in candidates:
            tsd_len = int(r - l + 1)
            if 2 <= tsd_len <= 30:
                strict_ok.append((support, l, r, source))
        if strict_ok:
            strict_ok.sort(key=lambda x: (x[0], x[2] - x[1]), reverse=True)
            _support, best_l, best_r, source = strict_ok[0]
            return (best_l, best_r, int(best_r - best_l + 1), source)

        # Rescue with ±2 bp shift when strict pairing misses by a few bases.
        rescue: list[tuple[int, int, int, str, int, int]] = []
        for l, r, support, source in candidates:
            sample_priority = 0 if source == "tsd_disease" else 1
            for dl in (-2, -1, 0, 1, 2):
                for dr in (-2, -1, 0, 1, 2):
                    ll = int(l + dl)
                    rr = int(r + dr)
                    if ll <= 0 or rr <= 0 or rr < ll:
                        continue
                    tsd_len = int(rr - ll + 1)
                    if 2 <= tsd_len <= 30:
                        shift_penalty = abs(dl) + abs(dr)
                        rescue.append((shift_penalty, -support, sample_priority, source, ll, rr))
        if not rescue:
            return (0, 0, 0, "")
        rescue.sort()
        _shift_penalty, _neg_support, _sample_priority, source, best_l, best_r = rescue[0]
        return (best_l, best_r, int(best_r - best_l + 1), source)

    tsd_pairs = out.apply(_pick_tsd_pair, axis=1, result_type="expand")
    tsd_pairs.columns = ["tsd_left_breakpoint", "tsd_right_breakpoint", "tsd_len_estimate", "tsd_evidence_source"]
    out["tsd_left_breakpoint"] = tsd_pairs["tsd_left_breakpoint"].astype(int)
    out["tsd_right_breakpoint"] = tsd_pairs["tsd_right_breakpoint"].astype(int)
    out["tsd_len_estimate"] = tsd_pairs["tsd_len_estimate"].astype(int)
    out["tsd_evidence_source"] = tsd_pairs["tsd_evidence_source"].fillna("").astype(str)
    # Strict TSD evidence threshold: 4 bp or longer.
    out["tsd_detected"] = out["tsd_len_estimate"] >= 4

    def _build_poly_breakpoint_mode_table(
        split_df: pd.DataFrame | None,
        *,
        prefix: str,
        min_poly_run: int = 8,
        min_poly_frac: float = 0.80,
    ) -> pd.DataFrame:
        key_cols = ["chrom", "window_start", "window_end"]
        out_cols = key_cols + [
            f"{prefix}_L_poly_breakpoint_mode",
            f"{prefix}_R_poly_breakpoint_mode",
            f"{prefix}_L_poly_breakpoint_support",
            f"{prefix}_R_poly_breakpoint_support",
        ]
        if split_df is None or split_df.empty:
            return pd.DataFrame(columns=out_cols)
        required = set(key_cols + ["clip_side", "pos", "read_name"])
        if not required.issubset(set(split_df.columns)):
            return pd.DataFrame(columns=out_cols)

        work = split_df.loc[:, key_cols + ["clip_side", "pos", "read_name"]].copy()
        work["clip_side"] = work["clip_side"].fillna("").astype(str).str.upper().str[:1]
        work = work.loc[work["clip_side"].isin(["L", "R"])].copy()
        if work.empty:
            return pd.DataFrame(columns=out_cols)
        work["pos"] = pd.to_numeric(work["pos"], errors="coerce").fillna(0).astype(int)
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[(work["pos"] > 0) & (work["read_name"].str.len() > 0)].copy()
        if work.empty:
            return pd.DataFrame(columns=out_cols)

        poly_flag = pd.Series(False, index=work.index)
        if "poly_tail_rescued" in split_df.columns:
            poly_tail_rescued = split_df.loc[work.index, "poly_tail_rescued"].fillna(False).astype(bool)
            poly_flag = poly_flag | poly_tail_rescued
        if "clip_poly_at_run" in split_df.columns:
            poly_run = pd.to_numeric(split_df.loc[work.index, "clip_poly_at_run"], errors="coerce").fillna(0).astype(int)
            poly_flag = poly_flag | poly_run.ge(int(min_poly_run))
        if "clip_poly_at_fraction" in split_df.columns:
            poly_frac = pd.to_numeric(
                split_df.loc[work.index, "clip_poly_at_fraction"],
                errors="coerce",
            ).fillna(0.0)
            poly_flag = poly_flag | poly_frac.ge(float(min_poly_frac))
        work = work.loc[poly_flag].copy()
        if work.empty:
            return pd.DataFrame(columns=out_cols)

        pos_support = (
            work.groupby(key_cols + ["clip_side", "pos"], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "support_reads"})
            .sort_values(
                key_cols + ["clip_side", "support_reads", "pos"],
                ascending=[True, True, True, True, False, True],
            )
        )
        modes = pos_support.drop_duplicates(key_cols + ["clip_side"], keep="first").copy()
        modes = modes.rename(columns={"pos": "poly_breakpoint_mode"})

        mode_pivot = (
            modes.pivot_table(
                index=key_cols,
                columns="clip_side",
                values="poly_breakpoint_mode",
                aggfunc="first",
            )
            .reset_index()
        )
        support_pivot = (
            modes.pivot_table(
                index=key_cols,
                columns="clip_side",
                values="support_reads",
                aggfunc="first",
            )
            .reset_index()
        )
        mode_pivot.columns = [str(c) for c in mode_pivot.columns]
        support_pivot.columns = [str(c) for c in support_pivot.columns]
        for side in ("L", "R"):
            if side not in mode_pivot.columns:
                mode_pivot[side] = 0
            if side not in support_pivot.columns:
                support_pivot[side] = 0
        merged = mode_pivot.merge(support_pivot, on=key_cols, suffixes=("_mode", "_support"), how="outer")
        for side in ("L", "R"):
            merged[f"{prefix}_{side}_poly_breakpoint_mode"] = (
                pd.to_numeric(merged.get(f"{side}_mode", 0), errors="coerce").fillna(0).astype(int)
            )
            merged[f"{prefix}_{side}_poly_breakpoint_support"] = (
                pd.to_numeric(merged.get(f"{side}_support", 0), errors="coerce").fillna(0).astype(int)
            )
        return merged[out_cols]

    key_cols = ["chrom", "window_start", "window_end"]
    disease_poly_modes = _build_poly_breakpoint_mode_table(split_disease, prefix="disease")
    control_poly_modes = _build_poly_breakpoint_mode_table(split_control, prefix="control")
    for tbl in (disease_poly_modes, control_poly_modes):
        if tbl.empty:
            continue
        out = out.merge(tbl, on=key_cols, how="left")
    for col in [
        "disease_L_poly_breakpoint_mode",
        "disease_R_poly_breakpoint_mode",
        "disease_L_poly_breakpoint_support",
        "disease_R_poly_breakpoint_support",
        "control_L_poly_breakpoint_mode",
        "control_R_poly_breakpoint_mode",
        "control_L_poly_breakpoint_support",
        "control_R_poly_breakpoint_support",
    ]:
        if col not in out.columns:
            out[col] = 0
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)

    def _build_split_breakpoint_mode_table(
        split_df: pd.DataFrame | None,
        *,
        prefix: str,
    ) -> pd.DataFrame:
        key_cols = ["chrom", "window_start", "window_end"]
        out_cols = key_cols + [
            f"{prefix}_L_split_breakpoint_mode",
            f"{prefix}_R_split_breakpoint_mode",
            f"{prefix}_L_split_breakpoint_support",
            f"{prefix}_R_split_breakpoint_support",
            f"{prefix}_L_split_breakpoint_total_reads",
            f"{prefix}_R_split_breakpoint_total_reads",
        ]
        if split_df is None or split_df.empty:
            return pd.DataFrame(columns=out_cols)
        required = set(key_cols + ["clip_side", "pos", "read_name"])
        if not required.issubset(set(split_df.columns)):
            return pd.DataFrame(columns=out_cols)

        work = split_df.loc[:, key_cols + ["clip_side", "pos", "read_name"]].copy()
        work["clip_side"] = work["clip_side"].fillna("").astype(str).str.upper().str[:1]
        work = work.loc[work["clip_side"].isin(["L", "R"])].copy()
        if work.empty:
            return pd.DataFrame(columns=out_cols)
        work["pos"] = pd.to_numeric(work["pos"], errors="coerce").fillna(0).astype(int)
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[(work["pos"] > 0) & (work["read_name"].str.len() > 0)].copy()
        if work.empty:
            return pd.DataFrame(columns=out_cols)

        pos_support = (
            work.groupby(key_cols + ["clip_side", "pos"], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "support_reads"})
            .sort_values(
                key_cols + ["clip_side", "support_reads", "pos"],
                ascending=[True, True, True, True, False, True],
            )
        )
        modes = pos_support.drop_duplicates(key_cols + ["clip_side"], keep="first").copy()
        modes = modes.rename(columns={"pos": "split_breakpoint_mode"})

        side_totals = (
            work.groupby(key_cols + ["clip_side"], as_index=False)["read_name"]
            .nunique()
            .rename(columns={"read_name": "total_reads"})
        )

        mode_pivot = (
            modes.pivot_table(
                index=key_cols,
                columns="clip_side",
                values="split_breakpoint_mode",
                aggfunc="first",
            )
            .reset_index()
        )
        support_pivot = (
            modes.pivot_table(
                index=key_cols,
                columns="clip_side",
                values="support_reads",
                aggfunc="first",
            )
            .reset_index()
        )
        total_pivot = (
            side_totals.pivot_table(
                index=key_cols,
                columns="clip_side",
                values="total_reads",
                aggfunc="first",
            )
            .reset_index()
        )
        mode_pivot.columns = [str(c) for c in mode_pivot.columns]
        support_pivot.columns = [str(c) for c in support_pivot.columns]
        total_pivot.columns = [str(c) for c in total_pivot.columns]
        for side in ("L", "R"):
            if side not in mode_pivot.columns:
                mode_pivot[side] = 0
            if side not in support_pivot.columns:
                support_pivot[side] = 0
            if side not in total_pivot.columns:
                total_pivot[side] = 0

        merged = mode_pivot.merge(support_pivot, on=key_cols, suffixes=("_mode", "_support"), how="outer")
        merged = merged.merge(total_pivot, on=key_cols, suffixes=("", "_total"), how="outer")
        for side in ("L", "R"):
            merged[f"{prefix}_{side}_split_breakpoint_mode"] = (
                pd.to_numeric(merged.get(f"{side}_mode", 0), errors="coerce").fillna(0).astype(int)
            )
            merged[f"{prefix}_{side}_split_breakpoint_support"] = (
                pd.to_numeric(merged.get(f"{side}_support", 0), errors="coerce").fillna(0).astype(int)
            )
            merged[f"{prefix}_{side}_split_breakpoint_total_reads"] = (
                pd.to_numeric(merged.get(side, 0), errors="coerce").fillna(0).astype(int)
            )
        return merged[out_cols]

    disease_split_modes = _build_split_breakpoint_mode_table(split_disease, prefix="disease")
    control_split_modes = _build_split_breakpoint_mode_table(split_control, prefix="control")
    for tbl in (disease_split_modes, control_split_modes):
        if tbl.empty:
            continue
        out = out.merge(tbl, on=key_cols, how="left")
    for col in [
        "disease_L_split_breakpoint_mode",
        "disease_R_split_breakpoint_mode",
        "disease_L_split_breakpoint_support",
        "disease_R_split_breakpoint_support",
        "disease_L_split_breakpoint_total_reads",
        "disease_R_split_breakpoint_total_reads",
        "control_L_split_breakpoint_mode",
        "control_R_split_breakpoint_mode",
        "control_L_split_breakpoint_support",
        "control_R_split_breakpoint_support",
        "control_L_split_breakpoint_total_reads",
        "control_R_split_breakpoint_total_reads",
    ]:
        if col not in out.columns:
            out[col] = 0
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)

    def _split_gap_tsd_rescue(
        row: pd.Series,
        *,
        min_len: int = 4,
        max_len: int = 40,
        min_side_support: int = 2,
        min_total_support: int = 6,
    ) -> tuple[int, int, int, str]:
        def _to_int_safe(value: object, default: int = 0) -> int:
            try:
                if pd.isna(value):
                    return int(default)
                return int(float(value))
            except (TypeError, ValueError):
                return int(default)

        if _to_int_safe(row.get("tsd_len_estimate", 0), 0) >= int(min_len):
            return (
                _to_int_safe(row.get("tsd_left_breakpoint", 0), 0),
                _to_int_safe(row.get("tsd_right_breakpoint", 0), 0),
                _to_int_safe(row.get("tsd_len_estimate", 0), 0),
                str(row.get("tsd_evidence_source", "") or ""),
            )

        any_mei_signal = max(
            _to_int_safe(row.get("disease_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("control_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("disease_discordant_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("control_discordant_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("disease_split_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("control_split_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("disease_full_mei_supported_reads", 0), 0),
            _to_int_safe(row.get("control_full_mei_supported_reads", 0), 0),
        )
        if any_mei_signal < 1:
            return (0, 0, 0, "")

        candidates: list[tuple[int, int, int, int, int, int, str]] = []
        for sample in ("disease", "control"):
            sample_pri = 0 if sample == "disease" else 1
            l_mode = _to_int_safe(row.get(f"{sample}_L_split_breakpoint_mode", 0), 0)
            r_mode = _to_int_safe(row.get(f"{sample}_R_split_breakpoint_mode", 0), 0)
            l_support = _to_int_safe(row.get(f"{sample}_L_split_breakpoint_support", 0), 0)
            r_support = _to_int_safe(row.get(f"{sample}_R_split_breakpoint_support", 0), 0)
            if l_mode <= 0 or r_mode <= 0:
                continue
            if l_support < int(min_side_support) or r_support < int(min_side_support):
                continue
            total_support = int(l_support + r_support)
            if total_support < int(min_total_support):
                continue
            if r_mode < l_mode:
                continue
            tsd_len = int(r_mode - l_mode + 1)
            if tsd_len < int(min_len) or tsd_len > int(max_len):
                continue
            l_poly_support = _to_int_safe(row.get(f"{sample}_L_poly_breakpoint_support", 0), 0)
            r_poly_support = _to_int_safe(row.get(f"{sample}_R_poly_breakpoint_support", 0), 0)
            poly_support = int(l_poly_support + r_poly_support)
            candidates.append(
                (
                    sample_pri,
                    -total_support,
                    -poly_support,
                    -tsd_len,
                    l_mode,
                    r_mode,
                    f"tsd_{sample}_split_clip_gap_rescue",
                )
            )
        if not candidates:
            return (0, 0, 0, "")
        candidates.sort()
        _sample_pri, _neg_total, _neg_poly, _neg_len, ll, rr, src = candidates[0]
        return (int(ll), int(rr), int(rr - ll + 1), str(src))

    gap_seed = out.apply(
        lambda r: _split_gap_tsd_rescue(r),
        axis=1,
        result_type="expand",
    )
    gap_seed.columns = [
        "gap_tsd_left_breakpoint",
        "gap_tsd_right_breakpoint",
        "gap_tsd_len_estimate",
        "gap_tsd_evidence_source",
    ]
    gap_len = pd.to_numeric(gap_seed["gap_tsd_len_estimate"], errors="coerce").fillna(0).astype(int)
    gap_mask = (out["tsd_len_estimate"].fillna(0).astype(int) < 4) & (gap_len >= 4)
    out.loc[gap_mask, "tsd_left_breakpoint"] = (
        pd.to_numeric(gap_seed.loc[gap_mask, "gap_tsd_left_breakpoint"], errors="coerce")
        .fillna(0)
        .astype(int)
    )
    out.loc[gap_mask, "tsd_right_breakpoint"] = (
        pd.to_numeric(gap_seed.loc[gap_mask, "gap_tsd_right_breakpoint"], errors="coerce")
        .fillna(0)
        .astype(int)
    )
    out.loc[gap_mask, "tsd_len_estimate"] = gap_len.loc[gap_mask].astype(int)
    out.loc[gap_mask, "tsd_evidence_source"] = (
        gap_seed.loc[gap_mask, "gap_tsd_evidence_source"].fillna("").astype(str)
    )
    out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

    def _one_sided_polyA_tsd_bridge_rescue(
        row: pd.Series,
        *,
        min_len: int = 4,
        max_len: int = 40,
    ) -> tuple[int, int, int, str]:
        def _to_int_safe(value: object, default: int = 0) -> int:
            try:
                if pd.isna(value):
                    return int(default)
                return int(float(value))
            except (TypeError, ValueError):
                return int(default)

        if _to_int_safe(row.get("tsd_len_estimate", 0), 0) >= int(min_len):
            return (
                _to_int_safe(row.get("tsd_left_breakpoint", 0), 0),
                _to_int_safe(row.get("tsd_right_breakpoint", 0), 0),
                _to_int_safe(row.get("tsd_len_estimate", 0), 0),
                str(row.get("tsd_evidence_source", "") or ""),
            )

        candidates: list[tuple[int, int, int, int, int, int, str]] = []
        for sample in ("disease", "control"):
            sample_pri = 0 if sample == "disease" else 1
            l_bp = _to_int_safe(row.get(f"{sample}_L_mei_breakpoint_mode", 0), 0)
            r_bp = _to_int_safe(row.get(f"{sample}_R_mei_breakpoint_mode", 0), 0)
            l_support = _to_int_safe(row.get(f"{sample}_L_mei_supported_reads", 0), 0)
            r_support = _to_int_safe(row.get(f"{sample}_R_mei_supported_reads", 0), 0)
            l_poly_bp = _to_int_safe(row.get(f"{sample}_L_poly_breakpoint_mode", 0), 0)
            r_poly_bp = _to_int_safe(row.get(f"{sample}_R_poly_breakpoint_mode", 0), 0)
            l_poly_support = _to_int_safe(row.get(f"{sample}_L_poly_breakpoint_support", 0), 0)
            r_poly_support = _to_int_safe(row.get(f"{sample}_R_poly_breakpoint_support", 0), 0)
            l_poly = _to_int_safe(row.get(f"{sample}_L_poly_at_reads", 0), 0)
            r_poly = _to_int_safe(row.get(f"{sample}_R_poly_at_reads", 0), 0)
            split_poly = _to_int_safe(row.get(f"split_{sample}_poly_tail_rescued_unique_reads", 0), 0)
            disc_poly = _to_int_safe(row.get(f"discordant_{sample}_poly_tail_rescued_unique_reads", 0), 0)
            poly_run = _to_int_safe(row.get(f"{sample}_poly_at_max_run", 0), 0)
            any_poly = max(l_poly, r_poly, split_poly, disc_poly, 1 if poly_run >= 8 else 0)
            if any_poly < 1:
                continue

            if l_bp > 0 and r_support <= 1 and l_support >= 1:
                if r_poly_bp > 0 and r_poly_bp >= l_bp and r_poly_support >= 1:
                    ll = int(l_bp)
                    rr = int(r_poly_bp)
                    tsd_len = int(rr - ll + 1)
                    if int(min_len) <= tsd_len <= int(max_len):
                        support = int(l_support + r_poly_support)
                        candidates.append(
                            (
                                sample_pri,
                                0,
                                -support,
                                -tsd_len,
                                ll,
                                rr,
                                f"tsd_{sample}_L_poly_mode_bridge_rescue",
                            )
                        )

            if r_bp > 0 and l_support <= 1 and r_support >= 1:
                if l_poly_bp > 0 and l_poly_bp <= r_bp and l_poly_support >= 1:
                    ll = int(l_poly_bp)
                    rr = int(r_bp)
                    tsd_len = int(rr - ll + 1)
                    if int(min_len) <= tsd_len <= int(max_len):
                        support = int(r_support + l_poly_support)
                        candidates.append(
                            (
                                sample_pri,
                                0,
                                -support,
                                -tsd_len,
                                ll,
                                rr,
                                f"tsd_{sample}_R_poly_mode_bridge_rescue",
                            )
                        )

        if not candidates:
            return (0, 0, 0, "")
        candidates.sort()
        _sample_pri, _fallback_rank, _neg_support, _neg_len, ll, rr, src = candidates[0]
        if ll <= 0 or rr < ll:
            return (0, 0, 0, "")
        return (int(ll), int(rr), int(rr - ll + 1), str(src))

    one_sided_seed = out.apply(
        lambda r: _one_sided_polyA_tsd_bridge_rescue(r),
        axis=1,
        result_type="expand",
    )
    one_sided_seed.columns = [
        "seed_tsd_left_breakpoint",
        "seed_tsd_right_breakpoint",
        "seed_tsd_len_estimate",
        "seed_tsd_evidence_source",
    ]
    seed_len = pd.to_numeric(one_sided_seed["seed_tsd_len_estimate"], errors="coerce").fillna(0).astype(int)
    seed_mask = (out["tsd_len_estimate"].fillna(0).astype(int) < 4) & (seed_len >= 4)
    out.loc[seed_mask, "tsd_left_breakpoint"] = (
        pd.to_numeric(one_sided_seed.loc[seed_mask, "seed_tsd_left_breakpoint"], errors="coerce")
        .fillna(0)
        .astype(int)
    )
    out.loc[seed_mask, "tsd_right_breakpoint"] = (
        pd.to_numeric(one_sided_seed.loc[seed_mask, "seed_tsd_right_breakpoint"], errors="coerce")
        .fillna(0)
        .astype(int)
    )
    out.loc[seed_mask, "tsd_len_estimate"] = seed_len.loc[seed_mask].astype(int)
    out.loc[seed_mask, "tsd_evidence_source"] = (
        one_sided_seed.loc[seed_mask, "seed_tsd_evidence_source"].fillna("").astype(str)
    )
    out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

    def _rescue_tsd_pair_with_reference(
        row: pd.Series,
        fetch_ref,
        *,
        shift_bp: int = 12,
        min_len: int = 4,
        max_len: int = 40,
    ) -> tuple[int, int, int, str]:
        if int(row.get("tsd_len_estimate", 0)) >= int(min_len):
            l = int(row.get("tsd_left_breakpoint", 0))
            r = int(row.get("tsd_right_breakpoint", 0))
            src = str(row.get("tsd_evidence_source", "") or "")
            return (l, r, int(max(0, row.get("tsd_len_estimate", 0))), src)

        chrom = str(row.get("chrom", "") or "").strip()
        if not chrom:
            return (0, 0, 0, "")

        seed_pairs: list[tuple[int, int, int, str]] = []
        t_l = int(row.get("disease_L_mei_breakpoint_mode", 0))
        t_r = int(row.get("disease_R_mei_breakpoint_mode", 0))
        t_support = int(row.get("disease_L_mei_supported_reads", 0)) + int(row.get("disease_R_mei_supported_reads", 0))
        if t_l > 0 and t_r > 0:
            seed_pairs.append((t_l, t_r, t_support, "tsd_disease"))
        n_l = int(row.get("control_L_mei_breakpoint_mode", 0))
        n_r = int(row.get("control_R_mei_breakpoint_mode", 0))
        n_support = int(row.get("control_L_mei_supported_reads", 0)) + int(row.get("control_R_mei_supported_reads", 0))
        if n_l > 0 and n_r > 0:
            seed_pairs.append((n_l, n_r, n_support, "tsd_control"))
        if not seed_pairs:
            return (0, 0, 0, "")

        seed_midpoints = [int((l + r) // 2) for l, r, _support, _source in seed_pairs if l > 0 and r > 0]
        bp_seed = int(seed_midpoints[0]) if seed_midpoints else 0

        best_key: tuple[int, int, int, int, int] | None = None
        best_value: tuple[int, int, int, str] | None = None
        for l0, r0, support, source in seed_pairs:
            src_priority = 0 if source == "tsd_disease" else 1
            for dl in range(-int(shift_bp), int(shift_bp) + 1):
                for dr in range(-int(shift_bp), int(shift_bp) + 1):
                    ll = int(l0 + dl)
                    rr = int(r0 + dr)
                    if ll <= 0 or rr <= 0 or rr < ll:
                        continue
                    tsd_len = int(rr - ll + 1)
                    if tsd_len < int(min_len) or tsd_len > int(max_len):
                        continue
                    seq = fetch_ref(chrom, ll - 1, rr)
                    if len(seq) != tsd_len or not seq:
                        continue
                    if "N" in seq:
                        continue
                    shift_penalty = abs(dl) + abs(dr)
                    mid = int((ll + rr) // 2)
                    mid_penalty = abs(mid - bp_seed) if bp_seed > 0 else 0
                    key = (shift_penalty, src_priority, mid_penalty, -int(support), -tsd_len)
                    if best_key is None or key < best_key:
                        best_key = key
                        best_value = (ll, rr, tsd_len, f"{source}_seq_rescue")

        if best_value is None:
            return (0, 0, 0, "")
        return best_value

    if reference_fasta is not None:
        def _build_clip_lookup(split_df: pd.DataFrame | None, *, max_probe_bp: int = 50) -> dict[tuple[str, int, int], dict[str, list[str]]]:
            lookup: dict[tuple[str, int, int], dict[str, list[str]]] = {}
            if split_df is None or split_df.empty:
                return lookup
            required = {"chrom", "window_start", "window_end", "clip_side", "clip_seq"}
            if not required.issubset(set(split_df.columns)):
                return lookup
            work = split_df.loc[:, ["chrom", "window_start", "window_end", "clip_side", "clip_seq"]].copy()
            work["clip_seq"] = work["clip_seq"].fillna("").astype(str)
            work = work.loc[work["clip_seq"].str.len() >= 4].copy()
            if work.empty:
                return lookup
            for row in work.itertuples(index=False):
                chrom = str(getattr(row, "chrom", "") or "")
                if not chrom:
                    continue
                try:
                    ws = int(getattr(row, "window_start"))
                    we = int(getattr(row, "window_end"))
                except (TypeError, ValueError):
                    continue
                side = str(getattr(row, "clip_side", "") or "").strip().upper()
                if side not in {"L", "R"}:
                    continue
                seq = str(getattr(row, "clip_seq", "") or "").upper()
                if len(seq) < 4:
                    continue
                if side == "L":
                    prox = seq[-int(max_probe_bp) :]
                else:
                    prox = seq[: int(max_probe_bp)]
                key = (chrom, ws, we)
                rec = lookup.setdefault(key, {"L": [], "R": []})
                rec[side].append(prox)
            return lookup

        with pysam.FastaFile(str(reference_fasta)) as ref:
            fetch_ref = _make_reference_fetcher(ref)
            rescued = out.apply(
                lambda r: _rescue_tsd_pair_with_reference(r, fetch_ref),
                axis=1,
                result_type="expand",
            )
        rescued.columns = [
            "resc_tsd_left_breakpoint",
            "resc_tsd_right_breakpoint",
            "resc_tsd_len_estimate",
            "resc_tsd_evidence_source",
        ]
        resc_len = pd.to_numeric(rescued["resc_tsd_len_estimate"], errors="coerce").fillna(0).astype(int)
        replace_mask = (out["tsd_len_estimate"].fillna(0).astype(int) < 4) & (resc_len >= 4)
        out.loc[replace_mask, "tsd_left_breakpoint"] = (
            pd.to_numeric(rescued.loc[replace_mask, "resc_tsd_left_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[replace_mask, "tsd_right_breakpoint"] = (
            pd.to_numeric(rescued.loc[replace_mask, "resc_tsd_right_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[replace_mask, "tsd_len_estimate"] = resc_len.loc[replace_mask].astype(int)
        out.loc[replace_mask, "tsd_evidence_source"] = (
            rescued.loc[replace_mask, "resc_tsd_evidence_source"].fillna("").astype(str)
        )
        out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

        disease_clip_lookup = _build_clip_lookup(split_disease)
        control_clip_lookup = _build_clip_lookup(split_control)

        def _best_match_stats(
            window_seq: str,
            motif: str,
            center_idx: int,
            *,
            max_mismatch: int = 1,
        ) -> tuple[int, int] | None:
            if not motif or not window_seq:
                return None
            m = len(motif)
            if m <= 0 or m > len(window_seq):
                return None
            best: tuple[int, int] | None = None
            for i in range(0, len(window_seq) - m + 1):
                mm = 0
                w = window_seq[i : i + m]
                for a, b in zip(w, motif):
                    if a != b:
                        mm += 1
                        if mm > int(max_mismatch):
                            break
                if mm > int(max_mismatch):
                    continue
                dist = abs(int(i) - int(center_idx))
                rank = (mm, dist)
                if best is None or rank < best:
                    best = rank
            if best is None:
                return None
            return (int(best[1]), int(best[0]))

        def _clip_exact_tsd_rescue(
            row: pd.Series,
            fetch_ref,
            *,
            min_len: int = 4,
            max_len: int = 30,
            flank_bp: int = 80,
            max_mismatch: int = 1,
            min_support_clips: int = 2,
        ) -> tuple[int, int, int, str]:
            chrom = str(row.get("chrom", "") or "").strip()
            if not chrom:
                return (0, 0, 0, "")
            try:
                ws = int(row.get("window_start", -1))
                we = int(row.get("window_end", -1))
            except (TypeError, ValueError):
                return (0, 0, 0, "")
            key = (chrom, ws, we)
            if key not in disease_clip_lookup and key not in control_clip_lookup:
                return (0, 0, 0, "")

            best_key: tuple[int, int, int, int, int] | None = None
            best_value: tuple[int, int, int, str] | None = None
            event_ori = str(row.get("insertion_orientation", "") or "").strip()
            if event_ori not in {"+", "-"}:
                event_ori = str(_choose_consolidated_insertion_orientation(row) or "").strip()
            for sample, sample_pri, lookup in (
                ("disease", 0, disease_clip_lookup),
                ("control", 1, control_clip_lookup),
            ):
                side_map = lookup.get(key, {})
                for side, side_pri in (("L", 0), ("R", 1)):
                    seqs = side_map.get(side, [])
                    if not seqs:
                        continue
                    try:
                        bp = int(row.get(f"{sample}_{side}_mei_breakpoint_mode", 0))
                    except (TypeError, ValueError):
                        bp = 0
                    if bp <= 0:
                        continue
                    start0 = max(0, bp - int(flank_bp) - 1)
                    end0 = max(start0 + 1, bp + int(flank_bp))
                    window_seq = fetch_ref(chrom, start0, end0)
                    if not window_seq:
                        continue
                    center_idx = int(bp - 1 - start0)
                    cand_votes: dict[tuple[int, int, str], int] = {}
                    cand_rank: dict[tuple[int, int, str], tuple[int, int, int, int, int]] = {}
                    for prox in seqs:
                        if len(prox) < int(min_len):
                            continue
                        oriented_variants: list[str]
                        if event_ori == "+":
                            oriented_variants = [prox]
                        elif event_ori == "-":
                            oriented_variants = [_revcomp(prox)]
                        else:
                            oriented_variants = [prox, _revcomp(prox)]
                        clip_best_key: tuple[int, int, str] | None = None
                        clip_best_rank: tuple[int, int, int, int, int] | None = None
                        for prox_oriented in oriented_variants:
                            if len(prox_oriented) < int(min_len):
                                continue
                            max_l = min(int(max_len), len(prox_oriented))
                            for L in range(max_l, int(min_len) - 1, -1):
                                motif = prox_oriented[-L:] if side == "L" else prox_oriented[:L]
                                if len(motif) != L or "N" in motif:
                                    continue
                                match = _best_match_stats(
                                    window_seq,
                                    motif,
                                    center_idx,
                                    max_mismatch=int(max_mismatch),
                                )
                                if match is None:
                                    continue
                                dist, mm = match
                                if side == "L":
                                    ll = int(bp)
                                    rr = int(bp + L - 1)
                                else:
                                    ll = int(bp - L + 1)
                                    rr = int(bp)
                                if ll <= 0 or rr < ll:
                                    continue
                                if mm == 0:
                                    src = f"tsd_{sample}_{side}_clip_exact_rescue"
                                else:
                                    src = f"tsd_{sample}_{side}_clip_near_exact_rescue"
                                ckey = (ll, rr, src)
                                crank = (-L, mm, dist, sample_pri, side_pri)
                                if clip_best_rank is None or crank < clip_best_rank:
                                    clip_best_rank = crank
                                    clip_best_key = ckey
                                break
                        if clip_best_key is None:
                            continue
                        cand_votes[clip_best_key] = int(cand_votes.get(clip_best_key, 0)) + 1
                        prev_rank = cand_rank.get(clip_best_key)
                        if prev_rank is None or (clip_best_rank is not None and clip_best_rank < prev_rank):
                            cand_rank[clip_best_key] = clip_best_rank if clip_best_rank is not None else (0, 0, 0, 0, 0)

                    for ckey, votes in cand_votes.items():
                        if int(votes) < int(min_support_clips):
                            continue
                        ll, rr, src = ckey
                        L = int(rr - ll + 1)
                        base_rank = cand_rank.get(ckey, (-L, 0, 0, sample_pri, side_pri))
                        rank = (base_rank[0], -int(votes), base_rank[1], base_rank[2], sample_pri)
                        if best_key is None or rank < best_key:
                            best_key = rank
                            best_value = (ll, rr, int(L), src)
            if best_value is None:
                return (0, 0, 0, "")
            return best_value

        clip_exact = out.apply(
            lambda r: _clip_exact_tsd_rescue(r, fetch_ref),
            axis=1,
            result_type="expand",
        )
        clip_exact.columns = [
            "clip_tsd_left_breakpoint",
            "clip_tsd_right_breakpoint",
            "clip_tsd_len_estimate",
            "clip_tsd_evidence_source",
        ]
        clip_len = pd.to_numeric(clip_exact["clip_tsd_len_estimate"], errors="coerce").fillna(0).astype(int)
        # Keep diagnostic columns so we can audit clip-rescue behavior on all loci.
        out["clip_tsd_len_estimate"] = clip_len.astype(int)
        out["clip_tsd_evidence_source"] = clip_exact["clip_tsd_evidence_source"].fillna("").astype(str)
        cur_len = out["tsd_len_estimate"].fillna(0).astype(int)
        cur_src = out["tsd_evidence_source"].fillna("").astype(str)
        cur_is_primary_pair = cur_src.isin(["tsd_disease", "tsd_control"])
        clip_is_primary_side = out["clip_tsd_evidence_source"].str.contains("_L_clip_|_R_clip_", regex=True)
        clip_mask = (
            (clip_len >= 4)
            & (
                (cur_len < 4)
                | ((clip_len > cur_len) & (~cur_is_primary_pair))
                | cur_src.str.contains("one_sided_polyA_rescue", regex=False)
                | (
                    (clip_len == cur_len)
                    & clip_is_primary_side
                    & cur_src.str.contains("seq_rescue", regex=False)
                )
            )
        )
        out["clip_tsd_applied"] = clip_mask.astype(bool)
        out.loc[clip_mask, "tsd_left_breakpoint"] = (
            pd.to_numeric(clip_exact.loc[clip_mask, "clip_tsd_left_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[clip_mask, "tsd_right_breakpoint"] = (
            pd.to_numeric(clip_exact.loc[clip_mask, "clip_tsd_right_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[clip_mask, "tsd_len_estimate"] = clip_len.loc[clip_mask].astype(int)
        out.loc[clip_mask, "tsd_evidence_source"] = (
            clip_exact.loc[clip_mask, "clip_tsd_evidence_source"].fillna("").astype(str)
        )
        out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

        def _one_sided_polyA_tsd_rescue(
            row: pd.Series,
            fetch_ref,
            *,
            min_len: int = 4,
        ) -> tuple[int, int, int, str]:
            def _to_int_safe(value: object, default: int = 0) -> int:
                try:
                    if pd.isna(value):
                        return int(default)
                    return int(float(value))
                except (TypeError, ValueError):
                    return int(default)

            if _to_int_safe(row.get("tsd_len_estimate", 0), 0) >= int(min_len):
                return (
                    _to_int_safe(row.get("tsd_left_breakpoint", 0), 0),
                    _to_int_safe(row.get("tsd_right_breakpoint", 0), 0),
                    _to_int_safe(row.get("tsd_len_estimate", 0), 0),
                    str(row.get("tsd_evidence_source", "") or ""),
                )

            chrom = str(row.get("chrom", "") or "").strip()
            if not chrom:
                return (0, 0, 0, "")

            candidates: list[tuple[int, int, int, int, str]] = []
            for sample in ("disease", "control"):
                l_bp = _to_int_safe(row.get(f"{sample}_L_mei_breakpoint_mode", 0), 0)
                r_bp = _to_int_safe(row.get(f"{sample}_R_mei_breakpoint_mode", 0), 0)
                l_support = _to_int_safe(row.get(f"{sample}_L_mei_supported_reads", 0), 0)
                r_support = _to_int_safe(row.get(f"{sample}_R_mei_supported_reads", 0), 0)
                sample_pri = 0 if sample == "disease" else 1

                # Disable fixed-length one-sided rescue; keep only interval-resolved
                # rescues (split-gap/bridge/clip-exact) to avoid noisy 4 bp artifacts.
                _ = sample_pri  # keep deterministic sample ordering if logic is re-enabled

            if not candidates:
                return (0, 0, 0, "")

            candidates.sort()
            for _sample_pri, _neg_support, ll, rr, src in candidates:
                if ll <= 0 or rr < ll:
                    continue
                seq = fetch_ref(chrom, ll - 1, rr)
                if len(seq) != int(min_len) or not seq or "N" in seq:
                    continue
                return (ll, rr, int(min_len), src)
            return (0, 0, 0, "")

        one_sided = out.apply(
            lambda r: _one_sided_polyA_tsd_rescue(r, fetch_ref),
            axis=1,
            result_type="expand",
        )
        one_sided.columns = [
            "one_tsd_left_breakpoint",
            "one_tsd_right_breakpoint",
            "one_tsd_len_estimate",
            "one_tsd_evidence_source",
        ]
        one_len = pd.to_numeric(one_sided["one_tsd_len_estimate"], errors="coerce").fillna(0).astype(int)
        one_mask = (out["tsd_len_estimate"].fillna(0).astype(int) < 4) & (one_len >= 4)
        out.loc[one_mask, "tsd_left_breakpoint"] = (
            pd.to_numeric(one_sided.loc[one_mask, "one_tsd_left_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[one_mask, "tsd_right_breakpoint"] = (
            pd.to_numeric(one_sided.loc[one_mask, "one_tsd_right_breakpoint"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        out.loc[one_mask, "tsd_len_estimate"] = one_len.loc[one_mask].astype(int)
        out.loc[one_mask, "tsd_evidence_source"] = (
            one_sided.loc[one_mask, "one_tsd_evidence_source"].fillna("").astype(str)
        )
        out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

    def _breakpoint_pos_and_source(row: pd.Series) -> tuple[int, str]:
        l = int(row.get("tsd_left_breakpoint", 0))
        r = int(row.get("tsd_right_breakpoint", 0))
        if l > 0 and r > 0:
            source = str(row.get("tsd_evidence_source", "") or "").strip() or "tsd_unknown"
            return int((l + r) // 2), source
        # Prefer MEI-mapped split modes, then raw split-clip modes. Soft-clipped
        # anchors (including one-sided) are junction-resolving.
        for prefix, label in (
            ("disease", "disease"),
            ("control", "control"),
        ):
            l = int(row.get(f"{prefix}_L_mei_breakpoint_mode", 0))
            r = int(row.get(f"{prefix}_R_mei_breakpoint_mode", 0))
            if l > 0 and r > 0:
                return int((l + r) // 2), f"{label}_split"
            if l > 0:
                return l, f"{label}_single"
            if r > 0:
                return r, f"{label}_single"
        for prefix, label in (
            ("disease", "disease"),
            ("control", "control"),
        ):
            l = int(row.get(f"{prefix}_L_split_breakpoint_mode", 0))
            r = int(row.get(f"{prefix}_R_split_breakpoint_mode", 0))
            l_sup = int(row.get(f"{prefix}_L_split_breakpoint_support", 0))
            r_sup = int(row.get(f"{prefix}_R_split_breakpoint_support", 0))
            if l > 0 and r > 0:
                return int((l + r) // 2), f"{label}_split_clip"
            # One-sided soft-clip / split-clip is enough for a point estimate.
            if l > 0 and l_sup >= 1:
                return l, f"{label}_single_clip"
            if r > 0 and r_sup >= 1:
                return r, f"{label}_single_clip"
        return 0, ""

    bp_fields = out.apply(_breakpoint_pos_and_source, axis=1, result_type="expand")
    bp_fields.columns = ["insertion_breakpoint_pos", "breakpoint_evidence_source"]
    out["insertion_breakpoint_pos"] = bp_fields["insertion_breakpoint_pos"].astype(int)
    out["breakpoint_evidence_source"] = bp_fields["breakpoint_evidence_source"].fillna("").astype(str)
    out["tsd_seq"] = ""
    out["breakpoint_context_11bp"] = ""
    out["breakpoint_l1_en_hexamer"] = ""
    out["breakpoint_l1_en_pattern"] = ""
    out["breakpoint_context_11bp_oriented"] = ""
    out["breakpoint_l1_en_hexamer_oriented"] = ""
    out["breakpoint_l1_en_pattern_yy_rrrr"] = ""
    out["breakpoint_l1_en_orientation_source"] = "unknown"
    out["breakpoint_l1_en_motif_like"] = False
    out["breakpoint_l1_en_best_motif"] = ""
    out["breakpoint_l1_en_motif_type"] = ""
    out["breakpoint_l1_en_mismatches"] = 99
    out["breakpoint_l1_en_mismatch_tolerance"] = 0
    out["breakpoint_l1_en_best_match_seq"] = ""
    out["breakpoint_l1_en_best_match_offset"] = 0
    out["breakpoint_l1_en_best_match_strand"] = "unknown"
    out["breakpoint_l1_en_best_match_anchor_6mer"] = ""
    out["breakpoint_l1_en_best_match_pattern_yy_rrrr"] = ""
    out["breakpoint_yyrrrr_logodds"] = float("nan")
    out["breakpoint_yyrrrr_logodds_shift1_max"] = float("nan")
    out["breakpoint_yyrrrr_best_offset"] = -1
    out["breakpoint_yyrrrr_logodds_shift1_mt_adj"] = float("nan")
    if reference_fasta is not None:
        with pysam.FastaFile(str(reference_fasta)) as ref:
            fetch_ref = _make_reference_fetcher(ref)
            seqs = []
            contexts_11bp: list[str] = []
            l1_hexamers: list[str] = []
            l1_patterns: list[str] = []
            contexts_11bp_oriented: list[str] = []
            l1_hexamers_oriented: list[str] = []
            l1_patterns_yy_rrrr: list[str] = []
            l1_orientation_source: list[str] = []
            l1_like: list[bool] = []
            l1_best_motif: list[str] = []
            l1_motif_type: list[str] = []
            l1_mismatches: list[int] = []
            l1_tolerance: list[int] = []
            l1_best_match_seq: list[str] = []
            l1_best_match_offset: list[int] = []
            l1_best_match_strand: list[str] = []
            l1_best_match_anchor_6mer: list[str] = []
            l1_best_match_pattern: list[str] = []
            yyrrrr_scores: list[float] = []
            yyrrrr_shift1_scores: list[float] = []
            yyrrrr_best_offsets: list[int] = []
            yyrrrr_shift1_mt_adj_scores: list[float] = []
            for row in out.itertuples(index=False):
                if int(row.tsd_len_estimate) <= 0:
                    seqs.append("")
                else:
                    chrom = str(row.chrom)
                    start0 = int(getattr(row, "tsd_left_breakpoint", 0)) - 1
                    end0 = int(getattr(row, "tsd_right_breakpoint", 0))
                    seqs.append(fetch_ref(chrom, start0, end0))

                bp = int(getattr(row, "insertion_breakpoint_pos", 0))
                if bp <= 0:
                    contexts_11bp.append("")
                    l1_hexamers.append("")
                    l1_patterns.append("")
                    contexts_11bp_oriented.append("")
                    l1_hexamers_oriented.append("")
                    l1_patterns_yy_rrrr.append("")
                    l1_orientation_source.append("unknown")
                    l1_like.append(False)
                    l1_best_motif.append("")
                    l1_motif_type.append("")
                    l1_mismatches.append(99)
                    l1_tolerance.append(0)
                    l1_best_match_seq.append("")
                    l1_best_match_offset.append(0)
                    l1_best_match_strand.append("unknown")
                    l1_best_match_anchor_6mer.append("")
                    l1_best_match_pattern.append("")
                    yyrrrr_scores.append(float("nan"))
                    yyrrrr_shift1_scores.append(float("nan"))
                    yyrrrr_best_offsets.append(-1)
                    yyrrrr_shift1_mt_adj_scores.append(float("nan"))
                    continue
                chrom = str(row.chrom)
                # 11 bp centered on breakpoint base (5 upstream + anchor + 5 downstream).
                start0_11 = max(0, bp - 6)
                end0_11 = max(start0_11 + 1, bp + 5)
                # 6 bp motif window near cleavage preference (4 upstream + 2 downstream).
                start0_6 = max(0, bp - 5)
                end0_6 = max(start0_6 + 1, bp + 1)
                ctx11 = fetch_ref(chrom, start0_11, end0_11)
                hex6 = fetch_ref(chrom, start0_6, end0_6)
                patt = f"{hex6[:4]}/{hex6[4:6]}" if len(hex6) == 6 else ""
                oriented_hex6, oriented_ctx11, orientation_source = _orient_to_insertion_strand(
                    hexamer=hex6,
                    context11bp=ctx11,
                    orientation=str(
                        _choose_consolidated_insertion_orientation(pd.Series(row._asdict()))
                    ),
                )
                patt_yy_rrrr = f"{oriented_hex6[:2]}/{oriented_hex6[2:6]}" if len(oriented_hex6) == 6 else ""
                allow_reverse_scan = orientation_source == "unknown"
                (
                    motif_like,
                    motif,
                    motif_type,
                    motif_mm,
                    motif_tol,
                    best_seq,
                    best_off,
                    best_strand,
                    best_anchor6,
                    best_pattern,
                ) = _match_l1_endonuclease_motif(
                    context11bp_oriented=oriented_ctx11,
                    allow_reverse_scan=allow_reverse_scan,
                )
                yyrrrr_score, yyrrrr_shift1_score, yyrrrr_best_off = _yyrrrr_logodds_with_shift_tolerance(
                    oriented_ctx11=oriented_ctx11
                )
                yyrrrr_shift1_mt_adj = _yyrrrr_shift1_logodds_mt_adjusted(yyrrrr_shift1_score)
                contexts_11bp.append(ctx11)
                l1_hexamers.append(hex6)
                l1_patterns.append(patt)
                contexts_11bp_oriented.append(oriented_ctx11)
                l1_hexamers_oriented.append(oriented_hex6)
                l1_patterns_yy_rrrr.append(patt_yy_rrrr)
                l1_orientation_source.append(orientation_source)
                l1_like.append(bool(motif_like))
                l1_best_motif.append(motif)
                l1_motif_type.append(motif_type)
                l1_mismatches.append(int(motif_mm))
                l1_tolerance.append(int(motif_tol))
                l1_best_match_seq.append(best_seq)
                l1_best_match_offset.append(int(best_off))
                l1_best_match_strand.append(best_strand)
                l1_best_match_anchor_6mer.append(best_anchor6)
                l1_best_match_pattern.append(best_pattern)
                yyrrrr_scores.append(float(yyrrrr_score))
                yyrrrr_shift1_scores.append(float(yyrrrr_shift1_score))
                yyrrrr_best_offsets.append(int(yyrrrr_best_off))
                yyrrrr_shift1_mt_adj_scores.append(float(yyrrrr_shift1_mt_adj))
            out["tsd_seq"] = seqs
            out["breakpoint_context_11bp"] = contexts_11bp
            out["breakpoint_l1_en_hexamer"] = l1_hexamers
            out["breakpoint_l1_en_pattern"] = l1_patterns
            out["breakpoint_context_11bp_oriented"] = contexts_11bp_oriented
            out["breakpoint_l1_en_hexamer_oriented"] = l1_hexamers_oriented
            out["breakpoint_l1_en_pattern_yy_rrrr"] = l1_patterns_yy_rrrr
            out["breakpoint_l1_en_orientation_source"] = l1_orientation_source
            out["breakpoint_l1_en_motif_like"] = l1_like
            out["breakpoint_l1_en_best_motif"] = l1_best_motif
            out["breakpoint_l1_en_motif_type"] = l1_motif_type
            out["breakpoint_l1_en_mismatches"] = l1_mismatches
            out["breakpoint_l1_en_mismatch_tolerance"] = l1_tolerance
            out["breakpoint_l1_en_best_match_seq"] = l1_best_match_seq
            out["breakpoint_l1_en_best_match_offset"] = l1_best_match_offset
            out["breakpoint_l1_en_best_match_strand"] = l1_best_match_strand
            out["breakpoint_l1_en_best_match_anchor_6mer"] = l1_best_match_anchor_6mer
            out["breakpoint_l1_en_best_match_pattern_yy_rrrr"] = l1_best_match_pattern
            out["breakpoint_yyrrrr_logodds"] = yyrrrr_scores
            out["breakpoint_yyrrrr_logodds_shift1_max"] = yyrrrr_shift1_scores
            out["breakpoint_yyrrrr_best_offset"] = yyrrrr_best_offsets
            out["breakpoint_yyrrrr_logodds_shift1_mt_adj"] = yyrrrr_shift1_mt_adj_scores

    # Guardrail: pure/near-pure polyA/polyT "TSDs" are poly-tail / reference-tail
    # artifacts, not target-site duplications. Apply to ALL evidence sources
    # (tsd_disease/tsd_control as well as *_rescue) — previously only rescue
    # rows were filtered, so primary pairs could publish AAAA…/TTTT… as TSD.
    poly_at_filter_mask = _clear_poly_at_artifact_tsd_fields(out)
    if bool(poly_at_filter_mask.any()):
        bp_fields = out.apply(_breakpoint_pos_and_source, axis=1, result_type="expand")
        bp_fields.columns = ["insertion_breakpoint_pos", "breakpoint_evidence_source"]
        out["insertion_breakpoint_pos"] = bp_fields["insertion_breakpoint_pos"].astype(int)
        out["breakpoint_evidence_source"] = bp_fields["breakpoint_evidence_source"].fillna("").astype(str)
    else:
        out["tsd_detected"] = out["tsd_len_estimate"].fillna(0).astype(int) >= 4

    # Weighted coherence metrics for ranking (annotation-only, no hard filtering).
    out["disease_breakpoint_mode_fraction_weighted"] = (
        out.get("disease_L_mei_breakpoint_mode_fraction", 0.0) * out.get("disease_L_mei_supported_reads", 0)
        + out.get("disease_R_mei_breakpoint_mode_fraction", 0.0) * out.get("disease_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "disease_mei_supported_reads", 0).replace(0, 1)
    out["control_breakpoint_mode_fraction_weighted"] = (
        out.get("control_L_mei_breakpoint_mode_fraction", 0.0) * out.get("control_L_mei_supported_reads", 0)
        + out.get("control_R_mei_breakpoint_mode_fraction", 0.0) * out.get("control_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "control_mei_supported_reads", 0).replace(0, 1)
    out["disease_subfamily_purity_weighted"] = (
        out.get("disease_L_mei_subfamily_purity", 0.0) * out.get("disease_L_mei_supported_reads", 0)
        + out.get("disease_R_mei_subfamily_purity", 0.0) * out.get("disease_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "disease_mei_supported_reads", 0).replace(0, 1)
    out["control_subfamily_purity_weighted"] = (
        out.get("control_L_mei_subfamily_purity", 0.0) * out.get("control_L_mei_supported_reads", 0)
        + out.get("control_R_mei_subfamily_purity", 0.0) * out.get("control_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "control_mei_supported_reads", 0).replace(0, 1)
    mapq_scaled = (_df_col_series(out, "split_disease_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0)
    out["coherence_score"] = (
        0.4 * out["disease_breakpoint_mode_fraction_weighted"].fillna(0.0)
        + 0.4 * out["disease_subfamily_purity_weighted"].fillna(0.0)
        + 0.2 * mapq_scaled.fillna(0.0)
    )
    out["control_background_score"] = (
        _df_col_series(out, "control_mei_supported_reads", 0).astype(float)
        + _df_col_series(out, "control_total_rows", 0).astype(float)
    )

    out["disease_poly_at_reads"] = _df_col_series(out, "disease_L_poly_at_reads", 0).fillna(0).astype(int) + _df_col_series(
        out, "disease_R_poly_at_reads", 0
    ).fillna(0).astype(int)
    out["control_poly_at_reads"] = _df_col_series(out, "control_L_poly_at_reads", 0).fillna(0).astype(int) + _df_col_series(
        out, "control_R_poly_at_reads", 0
    ).fillna(0).astype(int)
    out["disease_poly_at_max_run"] = (
        _df_col_series(out, "disease_L_poly_at_max_run", 0).fillna(0).astype(int).combine(
            _df_col_series(out, "disease_R_poly_at_max_run", 0).fillna(0).astype(int), max
        )
    )
    out["control_poly_at_max_run"] = (
        _df_col_series(out, "control_L_poly_at_max_run", 0).fillna(0).astype(int).combine(
            _df_col_series(out, "control_R_poly_at_max_run", 0).fillna(0).astype(int), max
        )
    )
    out["disease_poly_at_fraction_weighted"] = (
        _df_col_series(out, "disease_L_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * _df_col_series(out, "disease_L_mei_supported_reads", 0)
        + _df_col_series(out, "disease_R_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * _df_col_series(out, "disease_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "disease_mei_supported_reads", 0).replace(0, 1)
    out["control_poly_at_fraction_weighted"] = (
        _df_col_series(out, "control_L_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * _df_col_series(out, "control_L_mei_supported_reads", 0)
        + _df_col_series(out, "control_R_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * _df_col_series(out, "control_R_mei_supported_reads", 0)
    ) / _df_col_series(out, "control_mei_supported_reads", 0).replace(0, 1)

    out["insertion_orientation"] = out.apply(_choose_consolidated_insertion_orientation, axis=1)
    out["insertion_mei_span"] = out.apply(_choose_consolidated_insertion_mei_span, axis=1).astype(int)
    return out


def _apply_assembly_refinement_overrides(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    s = lambda col, default: _df_col_series(out, col, default)
    asm_source = s("asm_breakpoint_source", "").fillna("").astype(str)
    asm_has_mei = asm_source.isin(["disease", "control"])

    asm_bp = pd.to_numeric(s("asm_consensus_breakpoint_pos", float("nan")), errors="coerce")
    # Only override when assembly actually resolved a positive breakpoint.
    use_asm_bp = asm_has_mei & asm_bp.notna() & asm_bp.gt(0)
    out["insertion_breakpoint_pos"] = asm_bp.where(use_asm_bp, s("insertion_breakpoint_pos", 0)).fillna(0).astype(int)
    out.loc[use_asm_bp, "breakpoint_evidence_source"] = asm_source.loc[use_asm_bp]

    asm_tsd_seq = s("asm_tsd_seq", "").fillna("").astype(str)
    asm_tsd_len = pd.to_numeric(s("asm_tsd_len", float("nan")), errors="coerce")
    asm_tsd_detected = asm_has_mei & asm_tsd_len.notna() & asm_tsd_len.ge(4)
    out["tsd_seq"] = asm_tsd_seq.where(asm_tsd_detected & (asm_tsd_seq.str.len() > 0), s("tsd_seq", "").fillna("").astype(str))
    out["tsd_len_estimate"] = asm_tsd_len.where(asm_tsd_detected, s("tsd_len_estimate", 0)).fillna(0).astype(int)
    out["tsd_detected"] = out["tsd_len_estimate"].astype(float) >= 4.0
    # Assembly can also emit polyA/T tails as TSD; reject those too.
    _clear_poly_at_artifact_tsd_fields(out)

    asm_poly = pd.to_numeric(s("asm_polyA_max_run", float("nan")), errors="coerce")
    base_poly = pd.to_numeric(s("poly_at_max_run", 0), errors="coerce")
    picked_poly = asm_poly.where(asm_has_mei & asm_poly.notna(), base_poly)
    out["poly_at_max_run"] = (
        pd.concat([picked_poly, base_poly], axis=1).max(axis=1).fillna(0).astype(int)
    )

    # Do not override insertion orientation or MEI 5'/3' coords from local
    # assembly: contigs often fail to span the element and can invert/truncate
    # the footprint. Keep SR/DPE-derived orientation and coordinates.

    return out


def _recompute_breakpoint_sequence_metrics(candidates: pd.DataFrame, reference_fasta: Path | None) -> pd.DataFrame:
    out = candidates.copy()
    if reference_fasta is None:
        return out
    with pysam.FastaFile(str(reference_fasta)) as ref:
        fetch_ref = _make_reference_fetcher(ref)
        contexts_11bp: list[str] = []
        l1_hexamers: list[str] = []
        l1_patterns: list[str] = []
        contexts_11bp_oriented: list[str] = []
        l1_hexamers_oriented: list[str] = []
        l1_patterns_yy_rrrr: list[str] = []
        l1_orientation_source: list[str] = []
        l1_like: list[bool] = []
        l1_best_motif: list[str] = []
        l1_motif_type: list[str] = []
        l1_mismatches: list[int] = []
        l1_tolerance: list[int] = []
        l1_best_match_seq: list[str] = []
        l1_best_match_offset: list[int] = []
        l1_best_match_strand: list[str] = []
        l1_best_match_anchor_6mer: list[str] = []
        l1_best_match_pattern: list[str] = []
        yyrrrr_scores: list[float] = []
        yyrrrr_shift1_scores: list[float] = []
        yyrrrr_best_offsets: list[int] = []
        yyrrrr_shift1_mt_adj_scores: list[float] = []
        for row in out.itertuples(index=False):
            bp = int(getattr(row, "insertion_breakpoint_pos", 0) or 0)
            if bp <= 0:
                contexts_11bp.append("")
                l1_hexamers.append("")
                l1_patterns.append("")
                contexts_11bp_oriented.append("")
                l1_hexamers_oriented.append("")
                l1_patterns_yy_rrrr.append("")
                l1_orientation_source.append("unknown")
                l1_like.append(False)
                l1_best_motif.append("")
                l1_motif_type.append("")
                l1_mismatches.append(99)
                l1_tolerance.append(0)
                l1_best_match_seq.append("")
                l1_best_match_offset.append(0)
                l1_best_match_strand.append("unknown")
                l1_best_match_anchor_6mer.append("")
                l1_best_match_pattern.append("")
                yyrrrr_scores.append(float("nan"))
                yyrrrr_shift1_scores.append(float("nan"))
                yyrrrr_best_offsets.append(-1)
                yyrrrr_shift1_mt_adj_scores.append(float("nan"))
                continue
            chrom = str(getattr(row, "chrom", ""))
            start0_11 = max(0, bp - 6)
            end0_11 = max(start0_11 + 1, bp + 5)
            start0_6 = max(0, bp - 5)
            end0_6 = max(start0_6 + 1, bp + 1)
            ctx11 = fetch_ref(chrom, start0_11, end0_11)
            hex6 = fetch_ref(chrom, start0_6, end0_6)
            patt = f"{hex6[:4]}/{hex6[4:6]}" if len(hex6) == 6 else ""
            oriented_hex6, oriented_ctx11, orientation_source = _orient_to_insertion_strand(
                hexamer=hex6,
                context11bp=ctx11,
                orientation=str(getattr(row, "insertion_orientation", "")),
            )
            patt_yy_rrrr = f"{oriented_hex6[:2]}/{oriented_hex6[2:6]}" if len(oriented_hex6) == 6 else ""
            allow_reverse_scan = orientation_source == "unknown"
            motif_like, motif, motif_type, motif_mm, motif_tol, best_seq, best_off, best_strand, best_anchor6, best_pattern = (
                _match_l1_endonuclease_motif(
                    context11bp_oriented=oriented_ctx11,
                    allow_reverse_scan=allow_reverse_scan,
                )
            )
            yyrrrr_score, yyrrrr_shift1_score, yyrrrr_best_off = _yyrrrr_logodds_with_shift_tolerance(oriented_ctx11=oriented_ctx11)
            yyrrrr_shift1_mt_adj = _yyrrrr_shift1_logodds_mt_adjusted(yyrrrr_shift1_score)
            contexts_11bp.append(ctx11)
            l1_hexamers.append(hex6)
            l1_patterns.append(patt)
            contexts_11bp_oriented.append(oriented_ctx11)
            l1_hexamers_oriented.append(oriented_hex6)
            l1_patterns_yy_rrrr.append(patt_yy_rrrr)
            l1_orientation_source.append(orientation_source)
            l1_like.append(bool(motif_like))
            l1_best_motif.append(motif)
            l1_motif_type.append(motif_type)
            l1_mismatches.append(int(motif_mm))
            l1_tolerance.append(int(motif_tol))
            l1_best_match_seq.append(best_seq)
            l1_best_match_offset.append(int(best_off))
            l1_best_match_strand.append(best_strand)
            l1_best_match_anchor_6mer.append(best_anchor6)
            l1_best_match_pattern.append(best_pattern)
            yyrrrr_scores.append(float(yyrrrr_score))
            yyrrrr_shift1_scores.append(float(yyrrrr_shift1_score))
            yyrrrr_best_offsets.append(int(yyrrrr_best_off))
            yyrrrr_shift1_mt_adj_scores.append(float(yyrrrr_shift1_mt_adj))
    out["breakpoint_context_11bp"] = contexts_11bp
    out["breakpoint_l1_en_hexamer"] = l1_hexamers
    out["breakpoint_l1_en_pattern"] = l1_patterns
    out["breakpoint_context_11bp_oriented"] = contexts_11bp_oriented
    out["breakpoint_l1_en_hexamer_oriented"] = l1_hexamers_oriented
    out["breakpoint_l1_en_pattern_yy_rrrr"] = l1_patterns_yy_rrrr
    out["breakpoint_l1_en_orientation_source"] = l1_orientation_source
    out["breakpoint_l1_en_motif_like"] = l1_like
    out["breakpoint_l1_en_best_motif"] = l1_best_motif
    out["breakpoint_l1_en_motif_type"] = l1_motif_type
    out["breakpoint_l1_en_mismatches"] = l1_mismatches
    out["breakpoint_l1_en_mismatch_tolerance"] = l1_tolerance
    out["breakpoint_l1_en_best_match_seq"] = l1_best_match_seq
    out["breakpoint_l1_en_best_match_offset"] = l1_best_match_offset
    out["breakpoint_l1_en_best_match_strand"] = l1_best_match_strand
    out["breakpoint_l1_en_best_match_anchor_6mer"] = l1_best_match_anchor_6mer
    out["breakpoint_l1_en_best_match_pattern_yy_rrrr"] = l1_best_match_pattern
    out["breakpoint_yyrrrr_logodds"] = yyrrrr_scores
    out["breakpoint_yyrrrr_logodds_shift1_max"] = yyrrrr_shift1_scores
    out["breakpoint_yyrrrr_best_offset"] = yyrrrr_best_offsets
    out["breakpoint_yyrrrr_logodds_shift1_mt_adj"] = yyrrrr_shift1_mt_adj_scores
    return out


def _add_post_assembly_support_info_fields(
    candidates: pd.DataFrame,
    *,
    split_disease: pd.DataFrame,
    split_control: pd.DataFrame,
    discordant_disease: pd.DataFrame,
    discordant_control: pd.DataFrame,
) -> pd.DataFrame:
    out = candidates.copy()
    key_cols = ["chrom", "window_start", "window_end"]
    if out.empty:
        out["disease_supporting_reads_post_assembly"] = ""
        out["control_supporting_reads_post_assembly"] = ""
        return out

    bp_tbl = out.loc[:, key_cols + ["insertion_breakpoint_pos"]].copy()
    bp_tbl["insertion_breakpoint_pos"] = pd.to_numeric(bp_tbl["insertion_breakpoint_pos"], errors="coerce").fillna(0).astype(int)
    midpoint = (bp_tbl["window_start"].astype(int) + bp_tbl["window_end"].astype(int)) // 2
    bp_tbl["insertion_breakpoint_pos"] = bp_tbl["insertion_breakpoint_pos"].where(bp_tbl["insertion_breakpoint_pos"] > 0, midpoint)

    def _assembly_recruited_name_table(prefix: str) -> pd.DataFrame:
        names_col = f"asm_{prefix}_recruited_evidence_read_names"
        primary_col = f"asm_{prefix}_primary_contig_id"
        if names_col not in out.columns or primary_col not in out.columns:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        names = out.loc[:, key_cols + [names_col, primary_col]].copy()
        names[names_col] = names[names_col].fillna("").astype(str)
        has_mei_contig = names[primary_col].fillna("").astype(str).str.len() > 0
        names = names.loc[has_mei_contig & (names[names_col].str.len() > 0), key_cols + [names_col]].copy()
        if names.empty:
            return pd.DataFrame(columns=key_cols + ["read_name"])
        names["read_name"] = names[names_col].str.split(",")
        names = names.explode("read_name")
        names["read_name"] = names["read_name"].fillna("").astype(str).str.strip()
        names = names.loc[names["read_name"].str.len() > 0, key_cols + ["read_name"]].drop_duplicates()
        return names

    def _counts_from_split(df: pd.DataFrame, prefix: str, recruited: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns or recruited.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_post_asm", f"{prefix}_sr_r_post_asm"])
        cols = key_cols + ["read_name"]
        if "clip_side" in df.columns:
            cols.append("clip_side")
        if "pos" in df.columns:
            cols.append("pos")
        work = df.loc[:, [c for c in cols if c in df.columns]].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_post_asm", f"{prefix}_sr_r_post_asm"])
        work = work.merge(recruited, on=key_cols + ["read_name"], how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_post_asm", f"{prefix}_sr_r_post_asm"])
        work = work.merge(bp_tbl, on=key_cols, how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_sr_l_post_asm", f"{prefix}_sr_r_post_asm"])

        if "clip_side" in work.columns:
            side = work["clip_side"].fillna("").astype(str).str.upper().str[:1]
            if "pos" in work.columns:
                pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
                fallback = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
                side = side.where(side.isin(["L", "R"]), fallback)
            else:
                side = side.where(side.isin(["L", "R"]), "L")
        elif "pos" in work.columns:
            pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
            side = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
        else:
            side = pd.Series(["L"] * len(work), index=work.index)
        work["post_side"] = side

        agg = (
            work.groupby(key_cols + ["post_side"], as_index=False)["read_name"]
            .nunique()
            .pivot_table(index=key_cols, columns="post_side", values="read_name", fill_value=0)
            .reset_index()
        )
        agg.columns = [str(c) for c in agg.columns]
        if "L" not in agg.columns:
            agg["L"] = 0
        if "R" not in agg.columns:
            agg["R"] = 0
        agg[f"{prefix}_sr_l_post_asm"] = pd.to_numeric(agg["L"], errors="coerce").fillna(0).astype(int)
        agg[f"{prefix}_sr_r_post_asm"] = pd.to_numeric(agg["R"], errors="coerce").fillna(0).astype(int)
        return agg[key_cols + [f"{prefix}_sr_l_post_asm", f"{prefix}_sr_r_post_asm"]]

    def _counts_from_discordant(df: pd.DataFrame, prefix: str, recruited: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "read_name" not in df.columns or "pos" not in df.columns or recruited.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_post_asm", f"{prefix}_dpe_r_post_asm"])
        work = df.loc[:, key_cols + ["read_name", "pos"]].copy()
        work["read_name"] = work["read_name"].fillna("").astype(str)
        work = work.loc[work["read_name"].str.len() > 0].copy()
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_post_asm", f"{prefix}_dpe_r_post_asm"])
        work = work.merge(recruited, on=key_cols + ["read_name"], how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_post_asm", f"{prefix}_dpe_r_post_asm"])
        work = work.merge(bp_tbl, on=key_cols, how="inner")
        if work.empty:
            return pd.DataFrame(columns=key_cols + [f"{prefix}_dpe_l_post_asm", f"{prefix}_dpe_r_post_asm"])
        pos = pd.to_numeric(work["pos"], errors="coerce").fillna(work["insertion_breakpoint_pos"]).astype(int)
        work["post_side"] = pd.Series(["L"] * len(work), index=work.index).where(pos <= work["insertion_breakpoint_pos"], "R")
        agg = (
            work.groupby(key_cols + ["post_side"], as_index=False)["read_name"]
            .nunique()
            .pivot_table(index=key_cols, columns="post_side", values="read_name", fill_value=0)
            .reset_index()
        )
        agg.columns = [str(c) for c in agg.columns]
        if "L" not in agg.columns:
            agg["L"] = 0
        if "R" not in agg.columns:
            agg["R"] = 0
        agg[f"{prefix}_dpe_l_post_asm"] = pd.to_numeric(agg["L"], errors="coerce").fillna(0).astype(int)
        agg[f"{prefix}_dpe_r_post_asm"] = pd.to_numeric(agg["R"], errors="coerce").fillna(0).astype(int)
        return agg[key_cols + [f"{prefix}_dpe_l_post_asm", f"{prefix}_dpe_r_post_asm"]]

    for prefix, split_df, disc_df in (
        ("disease", split_disease, discordant_disease),
        ("control", split_control, discordant_control),
    ):
        recruited = _assembly_recruited_name_table(prefix)
        sr = _counts_from_split(split_df, prefix, recruited)
        dpe = _counts_from_discordant(disc_df, prefix, recruited)
        merged = bp_tbl.loc[:, key_cols].drop_duplicates().merge(sr, on=key_cols, how="left").merge(dpe, on=key_cols, how="left")
        sr_l = pd.to_numeric(merged.get(f"{prefix}_sr_l_post_asm", 0), errors="coerce").fillna(0).astype(int)
        sr_r = pd.to_numeric(merged.get(f"{prefix}_sr_r_post_asm", 0), errors="coerce").fillna(0).astype(int)
        dpe_l = pd.to_numeric(merged.get(f"{prefix}_dpe_l_post_asm", 0), errors="coerce").fillna(0).astype(int)
        dpe_r = pd.to_numeric(merged.get(f"{prefix}_dpe_r_post_asm", 0), errors="coerce").fillna(0).astype(int)
        merged[f"{prefix}_supporting_reads_post_assembly"] = [
            f"SR_L={sl},SR_R={srx},DPE_L={dl},DPE_R={dr}"
            for sl, srx, dl, dr in zip(sr_l.tolist(), sr_r.tolist(), dpe_l.tolist(), dpe_r.tolist())
        ]
        if f"{prefix}_supporting_reads_post_assembly" in out.columns:
            out = out.drop(columns=[f"{prefix}_supporting_reads_post_assembly"])
        out = out.merge(merged[key_cols + [f"{prefix}_supporting_reads_post_assembly"]], on=key_cols, how="left")
        out[f"{prefix}_supporting_reads_post_assembly"] = (
            out[f"{prefix}_supporting_reads_post_assembly"].fillna("SR_L=0,SR_R=0,DPE_L=0,DPE_R=0").astype(str)
        )
    return out


def _row_int(row: pd.Series, key: str, default: int = 0) -> int:
    val = row.get(key, default)
    if pd.isna(val):
        return default
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


def _row_bool(row: pd.Series, key: str, default: bool = False) -> bool:
    val = row.get(key, default)
    if pd.isna(val):
        return default
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        return bool(val)
    txt = str(val).strip().lower()
    if txt in {"true", "t", "1", "yes", "y"}:
        return True
    if txt in {"false", "f", "0", "no", "n"}:
        return False
    return default


def _sample_insertion_span_and_orientation(row: pd.Series, prefix: str) -> tuple[int, int, int, str]:
    l_start = _row_int(row, f"{prefix}_L_mei_start")
    r_start = _row_int(row, f"{prefix}_R_mei_start")
    l_end = _row_int(row, f"{prefix}_L_mei_end")
    r_end = _row_int(row, f"{prefix}_R_mei_end")
    l_support = _row_int(row, f"{prefix}_L_mei_supported_reads")
    r_support = _row_int(row, f"{prefix}_R_mei_supported_reads")
    l_anchor_bp = _row_int(row, f"{prefix}_L_mei_anchor_bp_max")
    r_anchor_bp = _row_int(row, f"{prefix}_R_mei_anchor_bp_max")
    l_target_len = _row_int(row, f"{prefix}_L_mei_target_len")
    r_target_len = _row_int(row, f"{prefix}_R_mei_target_len")
    l_poly_reads = _row_int(row, f"{prefix}_L_poly_at_reads")
    r_poly_reads = _row_int(row, f"{prefix}_R_poly_at_reads")
    l_poly_run = _row_int(row, f"{prefix}_L_poly_at_max_run")
    r_poly_run = _row_int(row, f"{prefix}_R_poly_at_max_run")
    raw_start = 0
    raw_end = 0
    strands = [
        s
        for s in [row.get(f"{prefix}_L_mei_strand", ""), row.get(f"{prefix}_R_mei_strand", "")]
        if s in {"+", "-"}
    ]
    if not strands:
        orient = ""
    elif len(set(strands)) == 1:
        orient = strands[0]
    else:
        orient = "mixed"

    l_strong = l_support >= 1 and l_anchor_bp >= _MIN_MEI_ANCHOR_BP and l_start > 0 and l_end >= l_start
    r_strong = r_support >= 1 and r_anchor_bp >= _MIN_MEI_ANCHOR_BP and r_start > 0 and r_end >= r_start
    l_relaxed = (
        l_support >= 1
        and l_anchor_bp >= _MIN_MEI_ANCHOR_BP_RELAXED
        and l_start > 0
        and l_end >= l_start
    )
    r_relaxed = (
        r_support >= 1
        and r_anchor_bp >= _MIN_MEI_ANCHOR_BP_RELAXED
        and r_start > 0
        and r_end >= r_start
    )
    l_poly_strong = l_poly_reads >= 1 and l_poly_run >= _MIN_POLYA_RUN_FOR_END_IMPUTE
    r_poly_strong = r_poly_reads >= 1 and r_poly_run >= _MIN_POLYA_RUN_FOR_END_IMPUTE

    if l_strong and r_strong:
        raw_start = min(l_start, r_start)
        raw_end = max(l_end, r_end)
    elif l_relaxed and r_relaxed:
        raw_start = min(l_start, r_start)
        raw_end = max(l_end, r_end)
    elif l_strong and r_poly_strong:
        tlen = max(l_target_len, r_target_len)
        if tlen > 0:
            raw_start = min(l_start, l_end)
            raw_end = max(tlen, max(l_start, l_end))
    elif l_relaxed and r_poly_strong:
        tlen = max(l_target_len, r_target_len)
        if tlen > 0:
            raw_start = min(l_start, l_end)
            raw_end = max(tlen, max(l_start, l_end))
    elif r_strong and l_poly_strong:
        tlen = max(r_target_len, l_target_len)
        if tlen > 0:
            raw_start = min(r_start, r_end)
            raw_end = max(tlen, max(r_start, r_end))
    elif r_relaxed and l_poly_strong:
        tlen = max(r_target_len, l_target_len)
        if tlen > 0:
            raw_start = min(r_start, r_end)
            raw_end = max(tlen, max(r_start, r_end))

    if raw_start <= 0 or raw_end < raw_start:
        d_left = _row_int(row, f"{prefix}_discordant_mei_left_supported_reads")
        d_right = _row_int(row, f"{prefix}_discordant_mei_right_supported_reads")
        d_two_sided = _row_bool(row, f"{prefix}_discordant_mei_two_sided_support", (d_left >= 1 and d_right >= 1))
        d_geom = _row_bool(row, f"{prefix}_discordant_mei_geometry_consistent", False)
        l_target = _row_int(row, f"{prefix}_discordant_mei_left_target_pos_median")
        r_target = _row_int(row, f"{prefix}_discordant_mei_right_target_pos_median")
        if d_two_sided and d_geom and l_target > 0 and r_target > 0:
            raw_start = min(l_target, r_target)
            raw_end = max(l_target, r_target)

    if orient not in {"+", "-"}:
        discordant_strand = str(row.get(f"{prefix}_discordant_mei_strand", "") or "").strip()
        if discordant_strand in {"+", "-"}:
            orient = discordant_strand

    if raw_start <= 0 or raw_end < raw_start:
        return 0, 0, 0, orient

    # Consensus coordinates are on the MEI reference axis; insertion strand does
    # not change which coordinate corresponds to element 3' vs 5'.
    # Under the project's 3'->5' convention, start is the higher coordinate.
    start = raw_end
    end = raw_start
    span = abs(end - start) + 1
    return start, end, span, orient


def _sample_has_bilateral_split_support(row: pd.Series, prefix: str) -> bool:
    left = _row_int(row, f"{prefix}_L_mei_supported_reads")
    right = _row_int(row, f"{prefix}_R_mei_supported_reads")
    return left >= 1 and right >= 1


def _sample_has_bilateral_discordant_support(row: pd.Series, prefix: str) -> bool:
    left = _row_int(row, f"{prefix}_discordant_mei_left_supported_reads")
    right = _row_int(row, f"{prefix}_discordant_mei_right_supported_reads")
    return left >= 1 and right >= 1


def _choose_consolidated_insertion_orientation(row: pd.Series) -> str:
    disease_orient = str(row.get("disease_insertion_orientation", "") or "").strip()
    control_orient = str(row.get("control_insertion_orientation", "") or "").strip()
    disease_bilateral = _sample_has_bilateral_split_support(row, "disease") or _sample_has_bilateral_discordant_support(
        row, "disease"
    )
    control_bilateral = _sample_has_bilateral_split_support(row, "control") or _sample_has_bilateral_discordant_support(
        row, "control"
    )
    if disease_bilateral and disease_orient in {"+", "-"}:
        return disease_orient
    if control_bilateral and control_orient in {"+", "-"}:
        return control_orient
    if disease_orient in {"+", "-"}:
        return disease_orient
    if control_orient in {"+", "-"}:
        return control_orient
    return _choose_event_orientation(row)


def _choose_consolidated_insertion_mei_span(row: pd.Series) -> int:
    disease_span = _row_int(row, "disease_insertion_mei_span")
    control_span = _row_int(row, "control_insertion_mei_span")
    disease_bilateral = _sample_has_bilateral_split_support(row, "disease") or _sample_has_bilateral_discordant_support(
        row, "disease"
    )
    control_bilateral = _sample_has_bilateral_split_support(row, "control") or _sample_has_bilateral_discordant_support(
        row, "control"
    )
    if disease_bilateral and disease_span > 0:
        return disease_span
    if control_bilateral and control_span > 0:
        return control_span
    if disease_span > 0 and control_span > 0:
        disease_reads = _row_int(row, "disease_mei_supported_reads")
        control_reads = _row_int(row, "control_mei_supported_reads")
        return disease_span if disease_reads >= control_reads else control_span
    return max(disease_span, control_span)


def _broaden_poly_at_fields(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()

    def _col_int(col: str) -> pd.Series:
        if col in out.columns:
            return out[col].fillna(0).astype(int)
        return pd.Series(0, index=out.index, dtype=int)

    def _max_int_series(*cols: str) -> pd.Series:
        parts = [_col_int(col) for col in cols]
        if not parts:
            return pd.Series(0, index=out.index, dtype=int)
        return pd.concat(parts, axis=1).max(axis=1).astype(int)

    for prefix in ("disease", "control"):
        out[f"{prefix}_poly_at_max_run"] = _max_int_series(
            f"{prefix}_poly_at_max_run",
            f"split_{prefix}_poly_tail_at_run_max",
            f"discordant_{prefix}_poly_tail_at_run_max",
            f"{prefix}_polya_rescue_max_len",
        )
        mei_poly_reads = _col_int(f"{prefix}_poly_at_reads")
        split_poly_reads = _col_int(f"split_{prefix}_poly_tail_rescued_unique_reads")
        discordant_poly_reads = _col_int(f"discordant_{prefix}_poly_tail_rescued_unique_reads")
        out[f"{prefix}_poly_at_reads"] = (
            pd.concat([mei_poly_reads, split_poly_reads], axis=1).max(axis=1).astype(int)
            + discordant_poly_reads
        )

    out["poly_at_max_run"] = _max_int_series("disease_poly_at_max_run", "control_poly_at_max_run")
    out["poly_at_reads"] = _col_int("disease_poly_at_reads") + _col_int("control_poly_at_reads")
    return out


def _add_consolidated_event_fields(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    s = lambda col, default: _df_col_series(out, col, default)
    for prefix in ("disease", "control"):
        out[f"{prefix}_left_supported_reads"] = (
            s(f"{prefix}_L_mei_supported_reads", 0).fillna(0).astype(int)
            + s(f"{prefix}_discordant_mei_left_supported_reads", 0).fillna(0).astype(int)
        )
        out[f"{prefix}_right_supported_reads"] = (
            s(f"{prefix}_R_mei_supported_reads", 0).fillna(0).astype(int)
            + s(f"{prefix}_discordant_mei_right_supported_reads", 0).fillna(0).astype(int)
        )
    out["mei_subfamily"] = out.apply(_choose_event_subfamily, axis=1)
    out["mei_family"] = out.apply(_choose_event_family, axis=1)
    # COMPLEX_INS is not a retrotransposon MEI call — do not carry nest-driven family.
    complex_ins = (
        _df_col_series(out, "insertion_event_class", "")
        .fillna("")
        .astype(str)
        .eq("COMPLEX_INS")
    )
    if complex_ins.any():
        out.loc[complex_ins, "mei_family"] = ""
        out.loc[complex_ins, "mei_subfamily"] = ""
    return out


def _agreement_flag(a: str, b: str) -> int:
    a = (a or "").strip()
    b = (b or "").strip()
    if not a and not b:
        return 0
    if not a or not b:
        # One-sided support can still be valid for low-support/subclonal events.
        return 1
    return 1 if a == b else 0


_COMPLEX_ANCHOR_MIN_UNIQUE_READS = 2
_COMPLEX_ANCHOR_MIN_FRACTION = 0.60
_COMPLEX_SPLIT_MIN_READS = 2
_COMPLEX_SPLIT_MIN_PURITY = 0.70
_COMPLEX_SPLIT_MIN_MODE_FRAC = 0.50
_COMPLEX_LOCUS_STRONG_MIN_FRACTION = 0.60
_COMPLEX_LOCUS_WEAK_MIN_FRACTION = 0.50
_COMPLEX_RESIDUAL_MIN_UNIQUE_READS = 2
# COMPLEX_INS: MEI_MAPPED must be weak relative to residual discordants.
_COMPLEX_INS_MAX_MEI_MAPPED = 2
_COMPLEX_INS_MEI_FRAC_OF_RESIDUAL = 0.25
_COMPLEX_INS_MIN_SPLIT_MODE_SUPPORT = 2


def _support_string_token(series: pd.Series, label: str) -> pd.Series:
    """Parse ``LABEL=N`` integer tokens from supporting-reads summary strings."""
    text = series.fillna("").astype(str)
    vals = pd.to_numeric(text.str.extract(rf"{re.escape(label)}=([0-9]+)", expand=False), errors="coerce")
    return vals.fillna(0).astype(int)


def _discordant_anchor_side_is_complex(
    unique_reads: pd.Series,
    complex_frac: pd.Series,
    mei_supported_on_side: pd.Series,
) -> pd.Series:
    return (
        (unique_reads.fillna(0).astype(float) >= _COMPLEX_ANCHOR_MIN_UNIQUE_READS)
        & (complex_frac.fillna(0.0).astype(float) >= _COMPLEX_ANCHOR_MIN_FRACTION)
        & (mei_supported_on_side.fillna(0).astype(float) <= 1)
    )


def _split_side_mei_for_complex(reads: pd.Series, purity: pd.Series, mode_frac: pd.Series) -> pd.Series:
    return (
        (reads.fillna(0).astype(float) >= _COMPLEX_SPLIT_MIN_READS)
        & (purity.fillna(0.0).astype(float) >= _COMPLEX_SPLIT_MIN_PURITY)
        & (mode_frac.fillna(0.0).astype(float) >= _COMPLEX_SPLIT_MIN_MODE_FRAC)
    )


def _df_col_float(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col in df.columns:
        return df[col].astype(float).fillna(default)
    return pd.Series(default, index=df.index, dtype=float)


def _df_col_series(df: pd.DataFrame, col: str, default: object) -> pd.Series:
    if col in df.columns:
        return df[col]
    return pd.Series([default] * len(df), index=df.index)


def _ensure_candidate_schema_defaults(candidates: pd.DataFrame) -> pd.DataFrame:
    """Guarantee optional evidence columns exist with safe defaults."""
    out = candidates.copy()
    defaults: dict[str, object] = {
        "disease_L_mei_supported_reads": 0,
        "disease_R_mei_supported_reads": 0,
        "control_L_mei_supported_reads": 0,
        "control_R_mei_supported_reads": 0,
        "disease_discordant_mei_left_supported_reads": 0,
        "disease_discordant_mei_right_supported_reads": 0,
        "control_discordant_mei_left_supported_reads": 0,
        "control_discordant_mei_right_supported_reads": 0,
        "disease_discordant_mei_supported_reads": 0,
        "control_discordant_mei_supported_reads": 0,
        "disease_discordant_mei_score_sum": 0.0,
        "control_discordant_mei_score_sum": 0.0,
        "disease_mei_supported_reads": 0,
        "control_mei_supported_reads": 0,
        "disease_total_rows": 0,
        "control_total_rows": 0,
        "disease_left_supported_reads": 0,
        "disease_right_supported_reads": 0,
        "control_left_supported_reads": 0,
        "control_right_supported_reads": 0,
        "disease_two_sided_support": False,
        "control_two_sided_support": False,
        "disease_family_agreement": 0,
        "control_family_agreement": 0,
        "disease_strand_agreement": 0,
        "control_strand_agreement": 0,
        "silver_stage_pass": False,
        "junk_flag_count": 999,
        "disease_poly_at_reads": 0,
        "control_poly_at_reads": 0,
        "poly_at_reads": 0,
        "poly_at_max_run": 0,
        "tsd_detected": False,
        "insertion_breakpoint_pos": 0,
        "asm_status": "",
        "known_mei_polymorphism": False,
        "known_mei_polymorphism_source": "",
        "known_mei_polymorphism_family": "",
        "known_mei_polymorphism_subfamily": "",
        "known_mei_polymorphism_id": "",
    }
    for col, default in defaults.items():
        if col not in out.columns:
            out[col] = default
    return out


def _complex_locus_fraction_cols(*, strong: bool) -> list[str]:
    """Prefer non-MEI residual fractions; fall back to all-discordant legacy columns."""
    residual_strong = [
        "disease_discordant_residual_large_insert_fraction",
        "disease_discordant_residual_interchrom_fraction",
        "disease_discordant_residual_mate_unmapped_fraction",
        "control_discordant_residual_large_insert_fraction",
        "control_discordant_residual_interchrom_fraction",
        "control_discordant_residual_mate_unmapped_fraction",
    ]
    residual_weak = [
        "disease_discordant_residual_same_strand_fraction",
        "disease_discordant_residual_improper_pair_fraction",
        "control_discordant_residual_same_strand_fraction",
        "control_discordant_residual_improper_pair_fraction",
    ]
    legacy_strong = [
        "discordant_disease_large_insert_fraction",
        "discordant_disease_interchrom_fraction",
        "discordant_disease_mate_unmapped_fraction",
        "discordant_control_large_insert_fraction",
        "discordant_control_interchrom_fraction",
        "discordant_control_mate_unmapped_fraction",
    ]
    legacy_weak = [
        "discordant_disease_same_strand_fraction",
        "discordant_disease_improper_pair_fraction",
        "discordant_control_same_strand_fraction",
        "discordant_control_improper_pair_fraction",
    ]
    if strong:
        return residual_strong + legacy_strong
    return residual_strong + residual_weak + legacy_strong + legacy_weak


def _complex_locus_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = _complex_locus_fraction_cols(strong=False)
    has_residual = any(
        c.startswith(("disease_discordant_residual_", "control_discordant_residual_")) and c in df.columns
        for c in fraction_cols
    )
    if has_residual:
        fraction_cols = [
            c for c in fraction_cols if c.startswith(("disease_discordant_residual_", "control_discordant_residual_"))
        ]
    else:
        fraction_cols = [c for c in fraction_cols if c.startswith("discordant_")]
    parts = [_df_col_float(df, col) for col in fraction_cols if col in df.columns]
    if not parts:
        return pd.Series(0.0, index=df.index)
    return pd.concat(parts, axis=1).max(axis=1)


def _complex_locus_strong_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = _complex_locus_fraction_cols(strong=True)
    has_residual = any(
        c.startswith(("disease_discordant_residual_", "control_discordant_residual_")) and c in df.columns
        for c in fraction_cols
    )
    if has_residual:
        fraction_cols = [
            c for c in fraction_cols if c.startswith(("disease_discordant_residual_", "control_discordant_residual_"))
        ]
    else:
        fraction_cols = [c for c in fraction_cols if c.startswith("discordant_")]
    parts = [_df_col_float(df, col) for col in fraction_cols if col in df.columns]
    if not parts:
        return pd.Series(0.0, index=df.index)
    return pd.concat(parts, axis=1).max(axis=1)


def _classic_polya_mei_sidepair(df: pd.DataFrame) -> pd.Series:
    """
    Classic simple MEI geometry: polyA/T on one flank + MEI support on the other.

    These should not be labeled complex even if residual discordants look noisy.
    """
    def _side_poly(prefix: str, side: str) -> pd.Series:
        reads = _df_col_float(df, f"{prefix}_{side}_poly_at_reads")
        run = _df_col_float(df, f"{prefix}_{side}_poly_at_max_run")
        return (reads >= 1) | (run >= 8)

    def _side_mei(prefix: str, side: str) -> pd.Series:
        split = _df_col_float(df, f"{prefix}_{side}_mei_supported_reads")
        disc = _df_col_float(
            df, f"{prefix}_discordant_mei_{'left' if side == 'L' else 'right'}_supported_reads"
        )
        return (split >= 1) | (disc >= 1)

    disease = (
        (_side_poly("disease", "L") & _side_mei("disease", "R"))
        | (_side_poly("disease", "R") & _side_mei("disease", "L"))
    )
    control = (
        (_side_poly("control", "L") & _side_mei("control", "R"))
        | (_side_poly("control", "R") & _side_mei("control", "L"))
    )
    return disease | control


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return (seq or "").translate(tr)[::-1]


def _hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


# Motif examples from published analyses; these are supportive mechanism hints,
# not strict pass/fail requirements (non-classical insertions may diverge).
_L1_EN_PAPER_MOTIFS: dict[str, str] = {
    "TTAAAA": "l1_en_canonical",
    "TTTAAA": "l1_en_canonical",
    "TTTTAA": "l1_en_canonical",
    "AAACTT": "l1_en_alternative",
    "CTGGG": "l1_en_alternative",
    "CCATT": "nested_novel_like",
}

# Motif-specific mismatch tolerance:
# - canonical 6bp motifs: allow up to 1 mismatch
# - alternative 6bp motif (AAACTT): allow up to 1 mismatch
# - shorter/novel-like 5bp motifs: allow up to 2 mismatches
_L1_EN_MOTIF_ALLOWED_MISMATCHES: dict[str, int] = {
    "TTAAAA": 1,
    "TTTAAA": 1,
    "TTTTAA": 1,
    "AAACTT": 1,
    "CTGGG": 2,
    "CCATT": 2,
}


def _yyrrrr_logodds(seq6: str) -> float:
    s = (seq6 or "").upper()
    if len(s) != 6:
        return 0.0
    favored = [
        {"C", "T"},
        {"C", "T"},
        {"A", "G"},
        {"A", "G"},
        {"A", "G"},
        {"A", "G"},
    ]
    score = 0.0
    for i, base in enumerate(s):
        p = 0.45 if base in favored[i] else 0.05
        score += float(math.log2(p / 0.25))
    return score


def _yyrrrr_logodds_with_shift_tolerance(oriented_ctx11: str) -> tuple[float, float, int]:
    ctx = (oriented_ctx11 or "").upper()
    if len(ctx) < 8:
        return (0.0, 0.0, 0)
    candidates: list[tuple[int, str]] = []
    for offset, start in [(-1, 0), (0, 1), (1, 2)]:
        end = start + 6
        if end <= len(ctx):
            candidates.append((offset, ctx[start:end]))
    if not candidates:
        return (0.0, 0.0, 0)
    scores = [(offset, _yyrrrr_logodds(seq)) for offset, seq in candidates]
    strict_score = next((sc for off, sc in scores if off == 0), scores[0][1])
    best_off, best_score = max(scores, key=lambda x: x[1])
    return (strict_score, best_score, int(best_off))


def _yyrrrr_shift1_logodds_mt_adjusted(best_score: float) -> float:
    # Multiple-testing adjustment for evaluating three offsets (-1, 0, +1).
    return float(best_score) - float(math.log2(3.0))


def _orient_to_insertion_strand(hexamer: str, context11bp: str, orientation: str) -> tuple[str, str, str]:
    ori = (orientation or "").strip()
    h = (hexamer or "").upper()
    c = (context11bp or "").upper()
    if ori == "+":
        return (h, c, "+")
    if ori == "-":
        return (_revcomp(h), _revcomp(c), "-")
    # Unknown/mixed orientation: keep reference orientation.
    return (h, c, "unknown")


def _match_l1_endonuclease_motif(
    context11bp_oriented: str,
    allow_reverse_scan: bool = True,
) -> tuple[bool, str, str, int, int, str, int, str, str, str]:
    q11 = (context11bp_oriented or "").upper()
    if len(q11) < 8:
        return (False, "", "", 99, 0, "", 0, "unknown", "", "")

    # Use the observed breakpoint-anchor 6-mer (offset 0) as the source of truth
    # for "best motif", so it always reflects what is closest to observed sequence.
    anchor6 = q11[1:7]
    if len(anchor6) != 6:
        return (False, "", "", 99, 0, "", 0, "unknown", "", "")

    best_motif = ""
    best_type = ""
    best_mm = 99
    best_seq = ""
    best_offset = 0
    best_strand = "forward"
    best_anchor6 = anchor6
    best_pattern = ""

    for motif, mtype in _L1_EN_PAPER_MOTIFS.items():
        mlen = len(motif)
        windows = [anchor6] if mlen == 6 else [anchor6[:5], anchor6[1:6]]
        for w_idx, win in enumerate(windows):
            if len(win) != mlen:
                continue
            mm = _hamming(win, motif)
            if mm < best_mm:
                best_mm = mm
                best_motif = motif
                best_type = mtype
                best_seq = win
                best_pattern = f"{win[:2]}/{win[2:]}" if len(win) >= 2 else win
                # For 5bp windows, index 0 is left-shifted slice and index 1 is right-shifted slice.
                best_offset = 0 if mlen == 6 else (0 if w_idx == 0 else 1)

    allowed_mm = _L1_EN_MOTIF_ALLOWED_MISMATCHES.get(best_motif, 0)
    motif_like = bool(best_motif) and best_mm <= allowed_mm
    return (
        motif_like,
        best_motif,
        best_type,
        best_mm,
        allowed_mm,
        best_seq,
        best_offset,
        best_strand,
        best_anchor6,
        best_pattern,
    )


def _compute_insertion_model_scores(candidates: pd.DataFrame) -> pd.DataFrame:
    out = _ensure_candidate_schema_defaults(candidates)
    s = lambda col, default: _df_col_series(out, col, default)

    for col in [
        "disease_L_mei_family",
        "disease_R_mei_family",
        "disease_L_mei_subfamily",
        "disease_R_mei_subfamily",
        "disease_L_mei_strand",
        "disease_R_mei_strand",
    ]:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str)

    out["disease_family_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["disease_L_mei_family"], out["disease_R_mei_family"])
    ]
    out["disease_subfamily_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["disease_L_mei_subfamily"], out["disease_R_mei_subfamily"])
    ]
    out["disease_strand_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["disease_L_mei_strand"], out["disease_R_mei_strand"])
    ]
    out["control_family_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            s("control_L_mei_family", "").fillna("").astype(str),
            s("control_R_mei_family", "").fillna("").astype(str),
        )
    ]
    out["control_subfamily_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            s("control_L_mei_subfamily", "").fillna("").astype(str),
            s("control_R_mei_subfamily", "").fillna("").astype(str),
        )
    ]
    out["control_strand_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            s("control_L_mei_strand", "").fillna("").astype(str),
            s("control_R_mei_strand", "").fillna("").astype(str),
        )
    ]

    disease_mei_reads = s("disease_mei_supported_reads", 0).astype(float)
    control_mei_reads = s("control_mei_supported_reads", 0).astype(float)
    total_rows = s("disease_total_rows", 0).astype(float).replace(0, 1.0)
    mei_enrichment = s("mei_score_enrichment_ratio", 0.0).astype(float)
    mei_enrichment_scaled = (mei_enrichment / (mei_enrichment + 1.0)).clip(lower=0.0, upper=1.0)
    mei_read_fraction = (disease_mei_reads / total_rows).clip(lower=0.0, upper=1.0)

    # Event-centric confidence score: do not bias to disease-only support.
    event_subfamily_purity = pd.concat(
        [
            s("disease_subfamily_purity_weighted", 0.0).astype(float).fillna(0.0),
            s("control_subfamily_purity_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    event_breakpoint_consistency = pd.concat(
        [
            s("disease_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
            s("control_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    event_family_agreement = pd.concat(
        [out["disease_family_agreement"].astype(float), out["control_family_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    event_subfamily_agreement = pd.concat(
        [out["disease_subfamily_agreement"].astype(float), out["control_subfamily_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    event_strand_agreement = pd.concat(
        [out["disease_strand_agreement"].astype(float), out["control_strand_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    control_mei_fraction = (
        control_mei_reads / s("control_total_rows", 0).astype(float).replace(0, 1.0)
    ).clip(lower=0.0, upper=1.0)
    event_mei_fraction = pd.concat([mei_read_fraction.fillna(0.0), control_mei_fraction.fillna(0.0)], axis=1).max(axis=1)
    mapq_event = pd.concat(
        [
            (s("split_disease_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0),
            (s("split_control_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0),
        ],
        axis=1,
    ).max(axis=1)

    tsd_boost = s("tsd_detected", False).fillna(False).astype(bool).astype(float)
    polyA_event = pd.concat(
        [
            s("disease_poly_at_fraction_weighted", 0.0).astype(float).fillna(0.0),
            s("control_poly_at_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1).clip(lower=0.0, upper=1.0)
    motif_boost = s("breakpoint_l1_en_motif_like", False).fillna(False).astype(bool).astype(float)
    motif_logodds = s("breakpoint_yyrrrr_logodds_shift1_mt_adj", 0.0).astype(float).fillna(0.0)
    motif_logodds_scaled = (motif_logodds / 6.0).clip(lower=0.0, upper=1.0)
    # Split-read overlap consistency across clipped sequences (k-mer Jaccard):
    # true loci tend to have multiple clips with mutually consistent overlap,
    # while slippage/noise loci often do not.
    overlap_med_cols = [
        "disease_L_clip_overlap_jaccard_median",
        "disease_R_clip_overlap_jaccard_median",
        "control_L_clip_overlap_jaccard_median",
        "control_R_clip_overlap_jaccard_median",
    ]
    overlap_n_cols = [
        "disease_L_clip_overlap_informative_reads",
        "disease_R_clip_overlap_informative_reads",
        "control_L_clip_overlap_informative_reads",
        "control_R_clip_overlap_informative_reads",
    ]
    overlap_med_tbl = pd.concat([s(c, 0.0).astype(float).fillna(0.0) for c in overlap_med_cols], axis=1)
    overlap_n_tbl = pd.concat([s(c, 0).astype(float).fillna(0.0) for c in overlap_n_cols], axis=1)
    # Keep side-wise alignment between median and informative-read tables.
    # The raw column names differ (jaccard_median vs informative_reads), so
    # rename the mask columns to the median-table columns before where().
    overlap_valid = overlap_n_tbl.copy()
    overlap_valid.columns = overlap_med_tbl.columns
    overlap_valid = overlap_valid >= 2.0
    overlap_med_masked = overlap_med_tbl.where(overlap_valid, 0.0)
    event_clip_overlap_consistency = overlap_med_masked.max(axis=1).fillna(0.0).clip(lower=0.0, upper=1.0)
    event_clip_overlap_informative_reads_max = overlap_n_tbl.max(axis=1).fillna(0.0)
    # One-good-side anchor: allow true asymmetric events where one side is
    # consistently aligned but the opposite (often polyA/T-rich) side is noisy.
    overlap_side_strong = ((overlap_n_tbl >= 3.0) & (overlap_med_masked >= 0.55)).any(axis=1)
    overlap_side_very_strong = ((overlap_n_tbl >= 2.0) & (overlap_med_masked >= 0.70)).any(axis=1)
    event_one_side_overlap_anchor = overlap_side_strong | overlap_side_very_strong
    out["event_clip_overlap_consistency"] = event_clip_overlap_consistency
    out["event_clip_overlap_informative_reads_max"] = event_clip_overlap_informative_reads_max.astype(int)
    out["event_one_side_overlap_anchor"] = event_one_side_overlap_anchor.astype(bool)

    base_score = (
        0.20 * event_subfamily_purity
        + 0.16 * event_breakpoint_consistency
        + 0.15 * mei_enrichment_scaled.fillna(0.0)
        + 0.10 * event_mei_fraction
        + 0.12 * event_family_agreement
        + 0.04 * event_subfamily_agreement
        + 0.06 * event_strand_agreement
        + 0.07 * tsd_boost
        + 0.05 * polyA_event
        + 0.03 * motif_boost
        + 0.02 * motif_logodds_scaled
        + 0.08 * event_clip_overlap_consistency
    )
    base_score = (base_score + 0.05 * mapq_event).clip(lower=0.0, upper=1.0)

    # Track complex SV-like companion signatures without suppressing MEI detection.
    # Prefer non-MEI residual discordant fractions when available (MEI-mapped
    # interchrom/large-insert mates are expected for simple insertions).
    complex_companion_fraction = _complex_locus_companion_fraction(out)
    complex_strong_companion_fraction = _complex_locus_strong_companion_fraction(out)
    has_residual = (
        "disease_discordant_residual_interchrom_fraction" in out.columns
        or "control_discordant_residual_interchrom_fraction" in out.columns
    )
    if has_residual:
        large_insert_fraction = _df_col_float(out, "disease_discordant_residual_large_insert_fraction")
        interchrom_fraction = _df_col_float(out, "disease_discordant_residual_interchrom_fraction")
        mate_unmapped_fraction = _df_col_float(out, "disease_discordant_residual_mate_unmapped_fraction")
        same_strand_fraction = _df_col_float(out, "disease_discordant_residual_same_strand_fraction")
        improper_pair_fraction = _df_col_float(out, "disease_discordant_residual_improper_pair_fraction")
        residual_unique = pd.concat(
            [
                _df_col_float(out, "disease_discordant_residual_unique_reads"),
                _df_col_float(out, "control_discordant_residual_unique_reads"),
            ],
            axis=1,
        ).max(axis=1)
    else:
        large_insert_fraction = _df_col_float(out, "discordant_disease_large_insert_fraction")
        interchrom_fraction = _df_col_float(out, "discordant_disease_interchrom_fraction")
        mate_unmapped_fraction = _df_col_float(out, "discordant_disease_mate_unmapped_fraction")
        same_strand_fraction = _df_col_float(out, "discordant_disease_same_strand_fraction")
        improper_pair_fraction = _df_col_float(out, "discordant_disease_improper_pair_fraction")
        residual_unique = pd.Series(float(_COMPLEX_RESIDUAL_MIN_UNIQUE_READS), index=out.index)

    residual_enough = residual_unique >= float(_COMPLEX_RESIDUAL_MIN_UNIQUE_READS)
    mei_mapped_fraction = pd.concat(
        [
            _df_col_float(out, "disease_discordant_mei_mapped_fraction"),
            _df_col_float(out, "control_discordant_mei_mapped_fraction"),
        ],
        axis=1,
    ).max(axis=1)
    # Fallback when residual metrics were patched from older annotations:
    # use supported-MEI DPE counts over all discordant unique reads.
    if "discordant_disease_unique_reads" in out.columns:
        stored_frac = (
            _df_col_float(out, "disease_discordant_mei_supported_reads")
            / _df_col_float(out, "discordant_disease_unique_reads").clip(lower=1.0)
        )
        mei_mapped_fraction = mei_mapped_fraction.combine(stored_frac, max)
    if "discordant_control_unique_reads" in out.columns:
        stored_frac_n = (
            _df_col_float(out, "control_discordant_mei_supported_reads")
            / _df_col_float(out, "discordant_control_unique_reads").clip(lower=1.0)
        )
        mei_mapped_fraction = mei_mapped_fraction.combine(stored_frac_n, max)
    # If a majority of DPEs remap to MEI, treat as simple insertion geometry.
    # Residual interchrom/large-insert mates are then insufficient for complex.
    mei_majority = mei_mapped_fraction >= 0.50
    out["discordant_mei_majority"] = mei_majority
    out["complex_sv_large_insert_flag"] = (
        residual_enough & ~mei_majority & (large_insert_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION)
    )
    out["complex_sv_interchrom_flag"] = (
        residual_enough & ~mei_majority & (interchrom_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION)
    )
    out["complex_sv_mate_unmapped_flag"] = (
        residual_enough & ~mei_majority & (mate_unmapped_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION)
    )
    out["complex_sv_companion_signal"] = (
        residual_enough
        & ~mei_majority
        & (complex_strong_companion_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION)
    )
    out["complex_sv_signal_score"] = complex_companion_fraction
    out["mei_with_complex_sv_signature"] = out["complex_sv_companion_signal"] & (disease_mei_reads >= 2)
    out["complex_sv_signature_label"] = "none"
    out.loc[out["complex_sv_large_insert_flag"], "complex_sv_signature_label"] = "large_insert"
    out.loc[out["complex_sv_interchrom_flag"], "complex_sv_signature_label"] = "interchrom"
    out.loc[
        out["complex_sv_large_insert_flag"] & out["complex_sv_interchrom_flag"],
        "complex_sv_signature_label",
    ] = "large_insert+interchrom"
    out.loc[
        out["complex_sv_mate_unmapped_flag"] & (out["complex_sv_signature_label"] == "none"),
        "complex_sv_signature_label",
    ] = "mate_unmapped"
    out.loc[
        out["complex_sv_mate_unmapped_flag"]
        & (out["complex_sv_signature_label"] != "none")
        & ~out["complex_sv_signature_label"].astype(str).str.contains("mate_unmapped", regex=False),
        "complex_sv_signature_label",
    ] = out["complex_sv_signature_label"] + "+mate_unmapped"
    out.loc[
        residual_enough
        & ~mei_majority
        & (same_strand_fraction >= _COMPLEX_LOCUS_WEAK_MIN_FRACTION)
        & (out["complex_sv_signature_label"] == "none"),
        "complex_sv_signature_label",
    ] = "same_strand"
    out.loc[
        residual_enough
        & ~mei_majority
        & (improper_pair_fraction >= _COMPLEX_LOCUS_WEAK_MIN_FRACTION)
        & (out["complex_sv_signature_label"] == "none"),
        "complex_sv_signature_label",
    ] = "improper_pair"

    score = base_score.clip(lower=0.0, upper=1.0)
    out["insertion_model_score"] = score
    left_reads = s("disease_L_mei_supported_reads", 0).astype(float)
    right_reads = s("disease_R_mei_supported_reads", 0).astype(float)
    discordant_mei_reads = s("disease_discordant_mei_supported_reads", 0).astype(float)
    left_mode_frac = s("disease_L_mei_breakpoint_mode_fraction", 0.0).astype(float).fillna(0.0)
    right_mode_frac = s("disease_R_mei_breakpoint_mode_fraction", 0.0).astype(float).fillna(0.0)
    left_purity = s("disease_L_mei_subfamily_purity", 0.0).astype(float).fillna(0.0)
    right_purity = s("disease_R_mei_subfamily_purity", 0.0).astype(float).fillna(0.0)
    out["disease_two_sided_support"] = (left_reads >= 1) & (right_reads >= 1)
    out["disease_two_sided_strong_support"] = (left_reads >= 2) & (right_reads >= 2)
    out["disease_one_sided_split_support"] = ((left_reads >= 2) & (right_reads < 2)) | (
        (right_reads >= 2) & (left_reads < 2)
    )
    out["disease_discordant_mei_strong_support"] = discordant_mei_reads >= 3
    dpe_left = s("disease_discordant_mei_left_supported_reads", 0).astype(float)
    dpe_right = s("disease_discordant_mei_right_supported_reads", 0).astype(float)
    dpe_family_purity = s("disease_discordant_mei_family_purity", 0.0).astype(float).fillna(0.0)
    dpe_strand_purity = s("disease_discordant_mei_strand_purity", 0.0).astype(float).fillna(0.0)
    dpe_geometry_consistent = (
        s("disease_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    )
    dpe_self_consistent = (
        s("disease_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    out["disease_discordant_mei_two_sided_support"] = (dpe_left >= 1) & (dpe_right >= 1)
    out["disease_discordant_mei_consistent_support"] = (
        out["disease_discordant_mei_two_sided_support"]
        & (dpe_family_purity >= 0.60)
        & (dpe_strand_purity >= 0.60)
        & dpe_geometry_consistent
        & dpe_self_consistent
    )
    control_left_reads = s("control_L_mei_supported_reads", 0).astype(float)
    control_right_reads = s("control_R_mei_supported_reads", 0).astype(float)
    out["control_two_sided_support"] = (control_left_reads >= 1) & (control_right_reads >= 1)
    control_dpe_left = s("control_discordant_mei_left_supported_reads", 0).astype(float)
    control_dpe_right = s("control_discordant_mei_right_supported_reads", 0).astype(float)
    control_dpe_family_purity = s("control_discordant_mei_family_purity", 0.0).astype(float).fillna(0.0)
    control_dpe_strand_purity = s("control_discordant_mei_strand_purity", 0.0).astype(float).fillna(0.0)
    control_dpe_geometry_consistent = (
        s("control_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    )
    control_dpe_self_consistent = (
        s("control_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    out["control_discordant_mei_two_sided_support"] = (control_dpe_left >= 1) & (control_dpe_right >= 1)
    out["control_discordant_mei_consistent_support"] = (
        out["control_discordant_mei_two_sided_support"]
        & (control_dpe_family_purity >= 0.60)
        & (control_dpe_strand_purity >= 0.60)
        & control_dpe_geometry_consistent
        & control_dpe_self_consistent
    )
    disease_left_mei_consistent = _split_side_mei_for_complex(
        left_reads,
        s("disease_L_mei_subfamily_purity", 0.0).astype(float),
        s("disease_L_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    disease_right_mei_consistent = _split_side_mei_for_complex(
        right_reads,
        s("disease_R_mei_subfamily_purity", 0.0).astype(float),
        s("disease_R_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    disease_left_anchor_complex = _discordant_anchor_side_is_complex(
        s("disease_discordant_anchor_left_unique_reads", 0).astype(float),
        s("disease_discordant_anchor_left_complex_reason_max_fraction", 0.0).astype(float),
        s("disease_discordant_mei_left_supported_reads", 0).astype(float),
    )
    disease_right_anchor_complex = _discordant_anchor_side_is_complex(
        s("disease_discordant_anchor_right_unique_reads", 0).astype(float),
        s("disease_discordant_anchor_right_complex_reason_max_fraction", 0.0).astype(float),
        s("disease_discordant_mei_right_supported_reads", 0).astype(float),
    )
    out["disease_discordant_anchor_left_complex_side"] = disease_left_anchor_complex
    out["disease_discordant_anchor_right_complex_side"] = disease_right_anchor_complex
    out["disease_mei_with_complex_sidepair"] = (
        (disease_left_mei_consistent & disease_right_anchor_complex)
        | (disease_right_mei_consistent & disease_left_anchor_complex)
    )

    control_left_mei_consistent = _split_side_mei_for_complex(
        control_left_reads,
        s("control_L_mei_subfamily_purity", 0.0).astype(float),
        s("control_L_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    control_right_mei_consistent = _split_side_mei_for_complex(
        control_right_reads,
        s("control_R_mei_subfamily_purity", 0.0).astype(float),
        s("control_R_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    control_left_anchor_complex = _discordant_anchor_side_is_complex(
        s("control_discordant_anchor_left_unique_reads", 0).astype(float),
        s("control_discordant_anchor_left_complex_reason_max_fraction", 0.0).astype(float),
        s("control_discordant_mei_left_supported_reads", 0).astype(float),
    )
    control_right_anchor_complex = _discordant_anchor_side_is_complex(
        s("control_discordant_anchor_right_unique_reads", 0).astype(float),
        s("control_discordant_anchor_right_complex_reason_max_fraction", 0.0).astype(float),
        s("control_discordant_mei_right_supported_reads", 0).astype(float),
    )
    out["control_discordant_anchor_left_complex_side"] = control_left_anchor_complex
    out["control_discordant_anchor_right_complex_side"] = control_right_anchor_complex
    out["control_mei_with_complex_sidepair"] = (
        (control_left_mei_consistent & control_right_anchor_complex)
        | (control_right_mei_consistent & control_left_anchor_complex)
    )
    out["disease_two_sided_like_support"] = out["disease_two_sided_strong_support"] | (
        out["disease_one_sided_split_support"] & out["disease_discordant_mei_strong_support"]
    ) | out["disease_discordant_mei_consistent_support"]
    out["disease_side_breakpoint_consistency"] = left_mode_frac.combine(right_mode_frac, min)
    out["disease_side_subfamily_purity"] = left_purity.combine(right_purity, min)
    out["disease_two_sided_family_consistent"] = out["disease_two_sided_support"] & (out["disease_family_agreement"] == 1)
    out["disease_two_sided_subfamily_consistent"] = out["disease_two_sided_support"] & (
        out["disease_subfamily_agreement"] == 1
    )
    out["event_two_sided_like_support"] = (
        out["disease_two_sided_like_support"]
        | out["control_two_sided_support"]
        | out["control_discordant_mei_consistent_support"]
    )
    out["event_family_consistent"] = (out["disease_family_agreement"] == 1) | (out["control_family_agreement"] == 1)
    out["event_strand_consistent"] = (out["disease_strand_agreement"] == 1) | (out["control_strand_agreement"] == 1)
    out["event_side_breakpoint_consistency"] = pd.concat(
        [
            out["disease_side_breakpoint_consistency"].astype(float).fillna(0.0),
            s("control_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    out["event_polyA_or_tsd_or_motif"] = (
        (tsd_boost >= 1.0) | (polyA_event >= 0.20) | (motif_boost >= 1.0) | (motif_logodds_scaled >= 0.25)
    )
    tsd_poly_filtered = s("tsd_poly_at_filter_applied", False).fillna(False).astype(bool)
    low_complexity_context = (polyA_event >= 0.20) | tsd_poly_filtered
    one_side_overlap_rescue = (
        event_one_side_overlap_anchor
        & (
            (tsd_boost >= 1.0)
            | out["disease_discordant_mei_consistent_support"]
            | out["control_discordant_mei_consistent_support"]
        )
    )
    overlap_consistent_for_context = (
        (~low_complexity_context)
        | (event_clip_overlap_consistency >= 0.20)
        | one_side_overlap_rescue
    )
    out["event_quality_clean"] = (
        (s("junk_flag_count", 0).fillna(0).astype(int) == 0)
        & (mapq_event >= 0.30)
    )

    # Structural confidence gates (sample-status agnostic, no explicit minimum read-count gate).
    high_conf_pass = (
        (out["insertion_model_score"] >= 0.60)
        & out["event_two_sided_like_support"]
        & out["event_family_consistent"]
        & out["event_strand_consistent"]
        & (out["event_side_breakpoint_consistency"] >= 0.55)
        & out["event_quality_clean"]
        & out["event_polyA_or_tsd_or_motif"]
        & (s("coherence_score", 0.0).astype(float) >= 0.55)
        & overlap_consistent_for_context
    )
    provisional_one_sided = (
        (~high_conf_pass)
        & (out["insertion_model_score"] >= 0.55)
        & (
            out["event_two_sided_like_support"]
            | ((disease_mei_reads + control_mei_reads) >= 1)
            | ((discordant_mei_reads + s("control_discordant_mei_supported_reads", 0).astype(float)) >= 1)
        )
        & out["event_family_consistent"]
        & (s("coherence_score", 0.0).astype(float) >= 0.50)
        & (
            (~low_complexity_context)
            | (event_clip_overlap_consistency >= 0.12)
            | one_side_overlap_rescue
        )
    )
    complex_sidepair_event = (
        out["disease_mei_with_complex_sidepair"] | out["control_mei_with_complex_sidepair"]
    )
    # Classic simple MEI: polyA/T on one flank + MEI support on the other.
    # Do not label these complex even if residual discordants look SV-like.
    out["classic_polya_mei_sidepair"] = _classic_polya_mei_sidepair(out)
    complex_sidepair_event = complex_sidepair_event & ~out["classic_polya_mei_sidepair"]
    out["mei_with_complex_sv_signature"] = (
        out["mei_with_complex_sv_signature"] & ~out["classic_polya_mei_sidepair"]
    )
    complex_sidepair_pass = (
        (~high_conf_pass)
        & (~provisional_one_sided)
        & complex_sidepair_event
        & out["event_family_consistent"]
        & out["event_strand_consistent"]
        & out["event_quality_clean"]
        & (out["insertion_model_score"] >= 0.50)
        & (s("coherence_score", 0.0).astype(float) >= 0.45)
    )
    out["complex_mei_event"] = complex_sidepair_event | out["mei_with_complex_sv_signature"]
    out["passes_insertion_model"] = high_conf_pass
    out["passes_insertion_model_provisional"] = provisional_one_sided
    out["passes_insertion_model_complex"] = complex_sidepair_pass
    out["insertion_call_tier"] = "none"
    out.loc[provisional_one_sided, "insertion_call_tier"] = "provisional_one_sided"
    out.loc[complex_sidepair_event, "insertion_call_tier"] = "mei_with_complex"
    out.loc[high_conf_pass, "insertion_call_tier"] = "high_conf_two_sided"

    # Careful public label: COMPLEX_INS when breakpoint pileup is real but the
    # majority of discordants are non-MEI residual (interchrom / large-insert).
    # Distinct from MEI_WITH_COMPLEX (real MEI + complex companion).
    disease_support = s("disease_supporting_reads", "").fillna("").astype(str)
    control_support = s("control_supporting_reads", "").fillna("").astype(str)
    d_sr_l = _support_string_token(disease_support, "SR_L")
    d_sr_r = _support_string_token(disease_support, "SR_R")
    d_dpe_l = _support_string_token(disease_support, "DPE_L")
    d_dpe_r = _support_string_token(disease_support, "DPE_R")
    c_sr_l = _support_string_token(control_support, "SR_L")
    c_sr_r = _support_string_token(control_support, "SR_R")
    c_dpe_l = _support_string_token(control_support, "DPE_L")
    c_dpe_r = _support_string_token(control_support, "DPE_R")
    bilateral_pileup = (
        ((d_sr_l + d_dpe_l) >= 1) & ((d_sr_r + d_dpe_r) >= 1)
    ) | (
        ((c_sr_l + c_dpe_l) >= 1) & ((c_sr_r + c_dpe_r) >= 1)
    )
    split_mode_support = pd.concat(
        [
            _df_col_float(out, "disease_L_split_breakpoint_support"),
            _df_col_float(out, "disease_R_split_breakpoint_support"),
            _df_col_float(out, "control_L_split_breakpoint_support"),
            _df_col_float(out, "control_R_split_breakpoint_support"),
        ],
        axis=1,
    ).max(axis=1)
    split_mode_pileup = split_mode_support >= float(_COMPLEX_INS_MIN_SPLIT_MODE_SUPPORT)
    breakpoint_pileup = bilateral_pileup | split_mode_pileup

    mei_mapped_max = pd.concat(
        [
            _support_string_token(disease_support, "MEI_MAPPED"),
            _support_string_token(control_support, "MEI_MAPPED"),
        ],
        axis=1,
    ).max(axis=1)
    mei_of_residual = mei_mapped_max.astype(float) / (
        mei_mapped_max.astype(float) + residual_unique.astype(float)
    ).clip(lower=1.0)
    weak_mei_for_complex_ins = (mei_mapped_max <= int(_COMPLEX_INS_MAX_MEI_MAPPED)) | (
        mei_of_residual < float(_COMPLEX_INS_MEI_FRAC_OF_RESIDUAL)
    )
    strong_residual_ins = (
        out["complex_sv_interchrom_flag"].fillna(False).astype(bool)
        | out["complex_sv_large_insert_flag"].fillna(False).astype(bool)
    )
    complex_ins = (
        residual_enough
        & (~mei_majority)
        & strong_residual_ins
        & breakpoint_pileup
        & (~out["classic_polya_mei_sidepair"].fillna(False).astype(bool))
        & weak_mei_for_complex_ins
    )
    mei_with_complex_class = (
        (~complex_ins)
        & (
            out["mei_with_complex_sv_signature"].fillna(False).astype(bool)
            | complex_sidepair_event
        )
        & (
            mei_majority
            | (mei_mapped_max >= 3)
        )
    )
    simple_mei = (
        (~complex_ins)
        & (~mei_with_complex_class)
        & (
            out["classic_polya_mei_sidepair"].fillna(False).astype(bool)
            | mei_majority
            | high_conf_pass
        )
    )
    out["insertion_event_class"] = "NONE"
    out.loc[simple_mei, "insertion_event_class"] = "SIMPLE_MEI"
    out.loc[mei_with_complex_class, "insertion_event_class"] = "MEI_WITH_COMPLEX"
    out.loc[complex_ins, "insertion_event_class"] = "COMPLEX_INS"
    out.loc[complex_ins, "insertion_call_tier"] = "complex_ins"

    # Sample presence: genotype on MEI_MAPPED (not SR/DPE). Soft-clip / discordant
    # pileup without MEI mate mapping is not treated as sample presence.
    #
    # Extreme imbalance (≥90% of total MEI_MAPPED in one sample) is also labeled
    # disease_only / control_only. Tumor/normal pairs often show a few bleed-
    # through MEI-mapped reads from contamination, adjacent tissue, or noise;
    # true shared germline events are usually much closer to balanced. Mosaic
    # intermediate cases stay shared when neither side reaches 90%.
    _mei_status_enrichment_frac = 0.90
    disease_mei = _mei_mapped_from_support_string(s("disease_supporting_reads", ""))
    control_mei = _mei_mapped_from_support_string(s("control_supporting_reads", ""))
    mei_total = (disease_mei + control_mei).astype(float)
    disease_present = disease_mei.ge(1)
    control_present = control_mei.ge(1)
    disease_frac = disease_mei.astype(float) / mei_total.clip(lower=1.0)
    control_frac = control_mei.astype(float) / mei_total.clip(lower=1.0)
    disease_enriched = disease_present & control_present & disease_frac.ge(_mei_status_enrichment_frac)
    control_enriched = disease_present & control_present & control_frac.ge(_mei_status_enrichment_frac)

    out["sample_status_label"] = "low_support"
    out.loc[disease_present & control_present, "sample_status_label"] = "shared"
    out.loc[disease_present & (~control_present), "sample_status_label"] = "disease_only"
    out.loc[(~disease_present) & control_present, "sample_status_label"] = "control_only"
    out.loc[disease_enriched, "sample_status_label"] = "disease_only"
    out.loc[control_enriched, "sample_status_label"] = "control_only"

    # Explicit convenience flag for downstream filtering.
    out["likely_false_positive_control_only"] = out["sample_status_label"] == "control_only"
    return out


def _consistent_family_mask(df: pd.DataFrame) -> pd.Series:
    family_cols = [
        "disease_L_mei_family",
        "disease_R_mei_family",
        "control_L_mei_family",
        "control_R_mei_family",
    ]
    missing = [c for c in family_cols if c not in df.columns]
    if missing:
        return pd.Series(False, index=df.index)

    def _is_consistent(row: pd.Series) -> bool:
        fams = [str(row[c]).strip() for c in family_cols]
        fams = [f for f in fams if f]
        if not fams:
            return False
        return len(set(fams)) == 1

    return df.apply(_is_consistent, axis=1)


def _depth_stats_for_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
) -> tuple[float, float]:
    if end_1based < start_1based:
        return (0.0, 0.0)
    start0 = max(0, int(start_1based) - 1)
    end0 = max(start0 + 1, int(end_1based))
    try:
        cov = bam.count_coverage(chrom, start0, end0, quality_threshold=0, read_callback="all")
    except ValueError:
        return (0.0, 0.0)
    span = end0 - start0
    if span <= 0:
        return (0.0, 0.0)
    depths = [a + c + g + t for a, c, g, t in zip(cov[0], cov[1], cov[2], cov[3])]
    if not depths:
        return (0.0, 0.0)
    total_depth = float(sum(depths))
    mean_depth = float(total_depth) / float(span)
    peak_depth = float(max(depths))
    return (mean_depth, peak_depth)


def _mean_depth_for_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
) -> float:
    mean_depth, _ = _depth_stats_for_interval(
        bam=bam,
        chrom=chrom,
        start_1based=start_1based,
        end_1based=end_1based,
    )
    return mean_depth


def _has_long_soft_clip(read: pysam.AlignedSegment, min_softclip: int = 20) -> bool:
    cigar = read.cigartuples
    if not cigar:
        return False
    first_op, first_len = cigar[0]
    if first_op == 4 and int(first_len) >= int(min_softclip):
        return True
    last_op, last_len = cigar[-1]
    return last_op == 4 and int(last_len) >= int(min_softclip)


def _is_non_sv_context_read(
    read: pysam.AlignedSegment,
    min_softclip: int = 20,
    discordant_abs_tlen_threshold: int = 1000,
) -> bool:
    if read.is_unmapped:
        return False
    if read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_supplementary:
        return False
    if _has_long_soft_clip(read, min_softclip=min_softclip):
        return False
    if read.has_tag("SA"):
        return False

    if read.is_paired:
        if read.mate_is_unmapped:
            return False
        if read.reference_id != read.next_reference_id:
            return False
        if abs(int(read.template_length)) >= int(discordant_abs_tlen_threshold):
            return False
        if not read.is_proper_pair:
            return False
    return True


def _context_quality_metrics_for_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
) -> dict[str, float]:
    start0 = max(0, int(start_1based) - 1)
    end0 = max(start0 + 1, int(end_1based))
    mapqs: list[int] = []
    nm_per_100bp: list[float] = []

    try:
        iterator = bam.fetch(chrom, start0, end0)
    except ValueError:
        return {
            "local_bam_mean_depth": 0.0,
            "local_bam_peak_depth": 0.0,
            "context_non_sv_reads": 0.0,
            "context_mapq_mean": 0.0,
            "context_mapq_lt20_fraction": 0.0,
            "context_nm_per_100bp_mean": 0.0,
            "context_nm_per_100bp_p90": 0.0,
        }

    for read in iterator:
        if not _is_non_sv_context_read(read):
            continue
        mapq = int(read.mapping_quality)
        mapqs.append(mapq)
        if read.has_tag("NM"):
            nm = int(read.get_tag("NM"))
            aligned_len = int(read.query_alignment_length or 0)
            if aligned_len > 0:
                nm_per_100bp.append((100.0 * float(nm)) / float(aligned_len))

    mapq_mean = float(sum(mapqs) / len(mapqs)) if mapqs else 0.0
    low_mapq_frac = float(sum(1 for q in mapqs if q < 20) / len(mapqs)) if mapqs else 0.0
    nm_mean = float(sum(nm_per_100bp) / len(nm_per_100bp)) if nm_per_100bp else 0.0
    nm_p90 = float(pd.Series(nm_per_100bp).quantile(0.9)) if nm_per_100bp else 0.0
    mean_depth, peak_depth = _depth_stats_for_interval(
        bam=bam, chrom=chrom, start_1based=start_1based, end_1based=end_1based
    )
    return {
        "local_bam_mean_depth": mean_depth,
        "local_bam_peak_depth": peak_depth,
        "context_non_sv_reads": float(len(mapqs)),
        "context_mapq_mean": mapq_mean,
        "context_mapq_lt20_fraction": low_mapq_frac,
        "context_nm_per_100bp_mean": nm_mean,
        "context_nm_per_100bp_p90": nm_p90,
    }


def _load_bed_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    def _open_textmaybe_gz(p: Path):
        if str(p).endswith(".gz"):
            return gzip.open(p, "rt", encoding="utf-8")
        return p.open("r", encoding="utf-8")

    def _parse_interval_parts(parts: list[str]) -> tuple[str, int, int] | None:
        if len(parts) >= 3 and parts[0].startswith("chr"):
            chrom = parts[0]
            start_idx = 1
            end_idx = 2
        elif len(parts) >= 4 and parts[1].startswith("chr"):
            chrom = parts[1]
            start_idx = 2
            end_idx = 3
        else:
            return None
        try:
            start0 = int(parts[start_idx])
            end0 = int(parts[end_idx])
        except ValueError:
            return None
        if end0 <= start0:
            return None
        return (chrom, start0 + 1, end0)

    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with _open_textmaybe_gz(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            parsed = _parse_interval_parts(parts)
            if parsed is None:
                continue
            chrom, start1, end1 = parsed
            intervals[chrom].append((int(start1), int(end1)))
    for chrom in list(intervals):
        intervals[chrom] = sorted(intervals[chrom], key=lambda x: (x[0], x[1]))
    return intervals


def _load_low_mappability_intervals(path: Path, threshold: float) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[0].startswith("chr"):
                chrom = parts[0]
                start_idx = 1
                end_idx = 2
                score_idx = 3
            elif len(parts) >= 5 and parts[1].startswith("chr"):
                chrom = parts[1]
                start_idx = 2
                end_idx = 3
                score_idx = 4
            else:
                continue
            try:
                start0 = int(parts[start_idx])
                end0 = int(parts[end_idx])
                score = float(parts[score_idx])
            except ValueError:
                continue
            if end0 <= start0 or score >= float(threshold):
                continue
            intervals[chrom].append((start0 + 1, end0))
    for chrom in list(intervals):
        intervals[chrom] = sorted(intervals[chrom], key=lambda x: (x[0], x[1]))
    return intervals


def _build_junk_interval_trees(
    segdup_bed: Path | None,
    low_mappability_bedgraph: Path | None,
    low_mappability_threshold: float,
    gap_bed: Path | None,
    encode_blacklist_bed: Path | None,
) -> dict[str, IntervalTree]:
    trees: dict[str, IntervalTree] = {}

    def _add_intervals(intervals: dict[str, list[tuple[int, int]]]) -> None:
        for chrom, rows in intervals.items():
            tree = trees.setdefault(str(chrom), IntervalTree())
            for start1, end1 in rows:
                tree.addi(int(start1), int(end1) + 1, 1)

    if segdup_bed is not None and segdup_bed.exists():
        _add_intervals(_load_bed_intervals(segdup_bed))
    if low_mappability_bedgraph is not None and low_mappability_bedgraph.exists():
        low_map_name = str(low_mappability_bedgraph).lower()
        if low_map_name.endswith(".bed") or low_map_name.endswith(".bed.gz"):
            _add_intervals(_load_bed_intervals(low_mappability_bedgraph))
        else:
            _add_intervals(_load_low_mappability_intervals(low_mappability_bedgraph, threshold=low_mappability_threshold))
    if gap_bed is not None and gap_bed.exists():
        _add_intervals(_load_bed_intervals(gap_bed))
    if encode_blacklist_bed is not None and encode_blacklist_bed.exists():
        _add_intervals(_load_bed_intervals(encode_blacklist_bed))
    return trees


def _write_interval_trees_to_bed(interval_trees: dict[str, IntervalTree], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for chrom in sorted(interval_trees):
            for iv in sorted(interval_trees[chrom], key=lambda x: (x.begin, x.end)):
                start0 = max(0, int(iv.begin) - 1)
                end0 = max(start0 + 1, int(iv.end) - 1)
                handle.write(f"{chrom}\t{start0}\t{end0}\n")


def _write_intervals_dict_to_bed(intervals: dict[str, list[tuple[int, int]]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for chrom in sorted(intervals):
            for start1, end1 in intervals[chrom]:
                start0 = max(0, int(start1) - 1)
                end0 = max(start0 + 1, int(end1))
                handle.write(f"{chrom}\t{start0}\t{end0}\n")


def _sample_random_windows_with_bedtools(
    target_chroms: list[str],
    reference_lengths: dict[str, int],
    sampled_span: int,
    n_windows: int,
    scope: str,
    random_seed: int,
    excluded_trees: dict[str, IntervalTree],
    highconf_bed: Path | None,
    junk_trees: dict[str, IntervalTree] | None = None,
    junk_exclusion_bed: Path | None = None,
) -> pd.DataFrame:
    if not target_chroms or sampled_span <= 0 or n_windows <= 0:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    rng = random.Random(int(random_seed))
    click.echo(
        f"[mei-annotate] empirical stage: bedtools shuffle start "
        f"scope={scope} n={n_windows} chroms={len(target_chroms)} span={sampled_span}"
    )
    with tempfile.TemporaryDirectory(prefix="rtm_empirical_shuffle_") as tmpdir:
        tmp = Path(tmpdir)
        genome_path = tmp / "genome.txt"
        seed_windows_path = tmp / "seed_windows.bed"
        excl_path = tmp / "exclude.bed"
        incl_path = tmp / "include.bed"
        shuffled_path = tmp / "shuffled.bed"

        with genome_path.open("w", encoding="utf-8") as gh:
            for chrom in target_chroms:
                gh.write(f"{chrom}\t{int(reference_lengths[chrom])}\n")

        # Seed intervals define count and lengths; shuffle randomizes positions.
        seeds: list[tuple[str, int, int]] = []
        if scope == "chromosome":
            for chrom in target_chroms:
                for _ in range(int(n_windows)):
                    seeds.append((chrom, 0, sampled_span))
        else:
            for _ in range(int(n_windows)):
                chrom = rng.choice(target_chroms)
                seeds.append((chrom, 0, sampled_span))
        with seed_windows_path.open("w", encoding="utf-8") as sh:
            for chrom, s0, e0 in seeds:
                sh.write(f"{chrom}\t{s0}\t{e0}\n")

        merged_excl: dict[str, IntervalTree] = {}
        for chrom, tree in excluded_trees.items():
            merged_excl.setdefault(chrom, IntervalTree()).update(tree)
        if junk_exclusion_bed is None and junk_trees is not None:
            for chrom, tree in junk_trees.items():
                merged_excl.setdefault(chrom, IntervalTree()).update(tree)
        for chrom in list(merged_excl):
            merged_excl[chrom].merge_overlaps()
        _write_interval_trees_to_bed(merged_excl, excl_path)
        if junk_exclusion_bed is not None and Path(junk_exclusion_bed).exists():
            with _open_textmaybe_gz(Path(junk_exclusion_bed)) as jh, excl_path.open("a", encoding="utf-8") as oh:
                for line in jh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    oh.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

        cmd = [
            "bedtools",
            "shuffle",
            "-i",
            str(seed_windows_path),
            "-g",
            str(genome_path),
            "-seed",
            str(int(random_seed)),
            "-chrom",
            "-excl",
            str(excl_path),
        ]
        if highconf_bed is not None:
            allowed = _load_bed_intervals(highconf_bed)
            if not allowed:
                return pd.DataFrame(columns=["chrom", "window_start", "window_end"])
            _write_intervals_dict_to_bed(allowed, incl_path)
            cmd.extend(["-incl", str(incl_path)])

        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        shuffled_path.write_text(proc.stdout, encoding="utf-8")
        rows: list[dict[str, int | str]] = []
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = str(parts[0])
            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            rows.append({"chrom": chrom, "window_start": start0 + 1, "window_end": end0})
        click.echo(
            f"[mei-annotate] empirical stage: bedtools shuffle done windows={len(rows)}"
        )
        return pd.DataFrame(rows)


def _sample_random_windows(
    candidates: pd.DataFrame,
    bam: pysam.AlignmentFile,
    n_windows: int,
    scope: str,
    random_seed: int,
    highconf_bed: Path | None,
    junk_trees: dict[str, IntervalTree] | None = None,
    junk_exclusion_bed: Path | None = None,
) -> pd.DataFrame:
    if n_windows <= 0 or candidates.empty:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    rng = random.Random(int(random_seed))
    spans = (candidates["window_end"].astype(int) - candidates["window_start"].astype(int) + 1).clip(lower=50)
    sampled_span = int(spans.median()) if len(spans) else 200

    excluded_trees: dict[str, IntervalTree] = {}
    for row in candidates.loc[:, ["chrom", "window_start", "window_end"]].itertuples(index=False):
        chrom = str(row.chrom)
        tree = excluded_trees.setdefault(chrom, IntervalTree())
        tree.addi(int(row.window_start), int(row.window_end) + 1, 1)

    reference_lengths = {str(chrom): int(length) for chrom, length in zip(bam.references, bam.lengths)}
    target_chroms = [str(c) for c in candidates["chrom"].astype(str).unique().tolist() if str(c) in reference_lengths]
    if not target_chroms:
        target_chroms = [str(c) for c in bam.references if str(c) in reference_lengths]

    # Prefer bedtools shuffle for interval randomization speed/reliability.
    try:
        click.echo("[mei-annotate] empirical stage: trying bedtools-based random sampling")
        sampled = _sample_random_windows_with_bedtools(
            target_chroms=target_chroms,
            reference_lengths=reference_lengths,
            sampled_span=sampled_span,
            n_windows=n_windows,
            scope=scope,
            random_seed=random_seed,
            excluded_trees=excluded_trees,
            highconf_bed=highconf_bed,
            junk_trees=junk_trees,
            junk_exclusion_bed=junk_exclusion_bed,
        )
        if not sampled.empty:
            return sampled
    except Exception:
        # Fall back to pure-Python sampling if bedtools shuffle is unavailable/fails.
        click.echo("[mei-annotate] empirical stage: bedtools sampling unavailable; using python fallback")
        pass

    allowed_intervals = _load_bed_intervals(highconf_bed) if highconf_bed is not None else {}
    if highconf_bed is not None:
        target_chroms = [c for c in target_chroms if c in allowed_intervals]
        if not target_chroms:
            return pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    targets: list[str] = []
    if scope == "chromosome":
        # Interpret n_windows as per-chromosome count when scope is chromosome.
        for chrom in target_chroms:
            for _ in range(int(n_windows)):
                targets.append(chrom)
        rng.shuffle(targets)
    else:
        for _ in range(int(n_windows)):
            targets.append(rng.choice(target_chroms))

    windows: list[dict[str, int | str]] = []
    target_total = len(targets)
    max_attempts = max(1000, target_total * 50)
    attempts = 0
    while len(windows) < target_total and attempts < max_attempts:
        attempts += 1
        chrom = targets[len(windows)] if len(windows) < len(targets) else rng.choice(target_chroms)
        chrom_len = int(reference_lengths.get(chrom, 0))
        if chrom_len < sampled_span:
            continue

        if highconf_bed is not None:
            intervals = [iv for iv in allowed_intervals.get(chrom, []) if (iv[1] - iv[0] + 1) >= sampled_span]
            if not intervals:
                continue
            iv_start, iv_end = rng.choice(intervals)
            max_start = iv_end - sampled_span + 1
            if max_start < iv_start:
                continue
            start = rng.randint(iv_start, max_start)
        else:
            start = rng.randint(1, chrom_len - sampled_span + 1)
        end = start + sampled_span - 1

        tree = excluded_trees.get(chrom)
        if tree is not None and tree.overlaps(start, end + 1):
            continue
        if junk_trees is not None:
            junk_tree = junk_trees.get(chrom)
            if junk_tree is not None and junk_tree.overlaps(start, end + 1):
                continue
        windows.append({"chrom": chrom, "window_start": int(start), "window_end": int(end)})

    click.echo(
        f"[mei-annotate] empirical stage: python fallback sampling done windows={len(windows)} "
        f"attempts={attempts}"
    )
    return pd.DataFrame(windows)


def _empirical_tail_prob(values: pd.Series, value: float, tail: str) -> float:
    arr = values.dropna().astype(float)
    n = int(len(arr))
    if n <= 0:
        return 1.0
    if tail == "high":
        k = int((arr >= float(value)).sum())
    else:
        k = int((arr <= float(value)).sum())
    return float(k + 1) / float(n + 1)


def _empirical_percentile(values: pd.Series, value: float) -> float:
    arr = values.dropna().astype(float)
    n = int(len(arr))
    if n <= 0:
        return 0.0
    k = int((arr <= float(value)).sum())
    return float(k) / float(n)


def _apply_empirical_context_scores(
    loci_metrics: pd.DataFrame,
    random_metrics: pd.DataFrame,
    sample_prefix: str,
    scope: str,
    progress_every: int = 0,
) -> pd.DataFrame:
    out = loci_metrics.copy()
    metric_specs: list[tuple[str, str]] = [
        ("local_bam_mean_depth", "high"),
        ("context_mapq_mean", "low"),
        ("context_mapq_lt20_fraction", "high"),
        ("context_nm_per_100bp_mean", "high"),
        ("context_nm_per_100bp_p90", "high"),
    ]

    out[f"{sample_prefix}_empirical_random_n"] = 0
    for metric, tail in metric_specs:
        out[f"{sample_prefix}_empirical_{metric}_percentile"] = 0.0
        out[f"{sample_prefix}_empirical_{metric}_p_{tail}"] = 1.0

    if random_metrics.empty:
        return out

    out[f"{sample_prefix}_empirical_random_n"] = int(len(random_metrics))
    global_lookup = {metric: random_metrics[metric] for metric, _ in metric_specs}
    by_chrom_lookup: dict[str, dict[str, pd.Series]] = {}
    if scope == "chromosome":
        for chrom, cdf in random_metrics.groupby("chrom", sort=False):
            by_chrom_lookup[str(chrom)] = {metric: cdf[metric] for metric, _ in metric_specs}

    total = int(len(out))
    for i, (idx, row) in enumerate(out.iterrows(), start=1):
        chrom = str(row.get("chrom", ""))
        for metric, tail in metric_specs:
            series = global_lookup[metric]
            if scope == "chromosome":
                chrom_series = by_chrom_lookup.get(chrom, {}).get(metric)
                if chrom_series is not None and len(chrom_series) >= 50:
                    series = chrom_series
            value = float(row.get(metric, 0.0) or 0.0)
            out.at[idx, f"{sample_prefix}_empirical_{metric}_percentile"] = _empirical_percentile(series, value)
            out.at[idx, f"{sample_prefix}_empirical_{metric}_p_{tail}"] = _empirical_tail_prob(series, value, tail=tail)
        if progress_every > 0 and (i % progress_every == 0 or i == total):
            click.echo(f"[mei-annotate] empirical scoring {sample_prefix}: {i}/{total} loci")
    return out


def _file_stamp(path: Path | None) -> dict[str, object]:
    if path is None:
        return {"path": "", "exists": False}
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {"path": str(p), "exists": True, "size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}


def _empirical_cache_key(
    loci: pd.DataFrame,
    disease_bam_path: Path,
    control_bam_path: Path,
    empirical_random_windows: int,
    empirical_random_scope: str,
    empirical_random_seed: int,
    empirical_highconf_bed: Path | None,
    empirical_exclude_merged_bed: Path | None,
    empirical_exclude_segdup_bed: Path | None,
    empirical_exclude_mappability_bedgraph: Path | None,
    empirical_exclude_mappability_threshold: float,
    empirical_exclude_gap_bed: Path | None,
    empirical_exclude_blacklist_bed: Path | None,
) -> str:
    loci_view = loci.loc[:, ["chrom", "window_start", "window_end"]].copy()
    chrom_counts = (
        loci_view.groupby("chrom", sort=True).size().to_dict() if not loci_view.empty else {}
    )
    spans = (
        (loci_view["window_end"].astype(int) - loci_view["window_start"].astype(int) + 1).tolist()
        if not loci_view.empty
        else []
    )
    payload = {
        "version": "empirical_cache_v1",
        "loci_count": int(len(loci_view)),
        "chrom_counts": {str(k): int(v) for k, v in chrom_counts.items()},
        "span_median": float(pd.Series(spans).median()) if spans else 0.0,
        "disease_bam": _file_stamp(disease_bam_path),
        "control_bam": _file_stamp(control_bam_path),
        "random_windows": int(empirical_random_windows),
        "random_scope": str(empirical_random_scope),
        "random_seed": int(empirical_random_seed),
        "highconf": _file_stamp(empirical_highconf_bed),
        "merged_exclusion": _file_stamp(empirical_exclude_merged_bed),
        "segdup": _file_stamp(empirical_exclude_segdup_bed),
        "mappability": _file_stamp(empirical_exclude_mappability_bedgraph),
        "mappability_threshold": float(empirical_exclude_mappability_threshold),
        "gap": _file_stamp(empirical_exclude_gap_bed),
        "blacklist": _file_stamp(empirical_exclude_blacklist_bed),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]


def _annotate_bam_depth_for_consistent_loci(
    candidates: pd.DataFrame,
    disease_bam_path: Path,
    control_bam_path: Path,
    empirical_random_windows: int = 1000,
    empirical_random_scope: str = "chromosome",
    empirical_random_seed: int = 13,
    empirical_highconf_bed: Path | None = None,
    empirical_exclude_merged_bed: Path | None = None,
    empirical_exclude_segdup_bed: Path | None = None,
    empirical_exclude_mappability_bedgraph: Path | None = None,
    empirical_exclude_mappability_threshold: float = 0.5,
    empirical_exclude_gap_bed: Path | None = None,
    empirical_exclude_blacklist_bed: Path | None = None,
    empirical_cache_dir: Path | None = None,
) -> pd.DataFrame:
    stage_start = time.monotonic()
    out = candidates.copy()
    s = lambda col, default: _df_col_series(out, col, default)
    out["depth_filter_family_consistent"] = False
    out["depth_filter_two_sided_consistent"] = False
    out["depth_filter_pass"] = False
    out["disease_local_bam_mean_depth"] = 0.0
    out["control_local_bam_mean_depth"] = 0.0
    out["disease_local_bam_peak_depth"] = 0.0
    out["control_local_bam_peak_depth"] = 0.0
    out["disease_context_non_sv_reads"] = 0
    out["control_context_non_sv_reads"] = 0
    out["disease_context_mapq_mean"] = 0.0
    out["control_context_mapq_mean"] = 0.0
    out["disease_context_mapq_lt20_fraction"] = 0.0
    out["control_context_mapq_lt20_fraction"] = 0.0
    out["disease_context_nm_per_100bp_mean"] = 0.0
    out["control_context_nm_per_100bp_mean"] = 0.0
    out["disease_context_nm_per_100bp_p90"] = 0.0
    out["control_context_nm_per_100bp_p90"] = 0.0
    out["disease_mei_support_per_100x_bam_depth"] = 0.0
    out["control_mei_support_per_100x_bam_depth"] = 0.0
    out["mei_support_per_100x_bam_depth_delta"] = 0.0
    out["mei_support_per_100x_bam_depth_ratio"] = 1.0

    if out.empty:
        return out
    t0 = time.monotonic()

    family_consistent = _consistent_family_mask(out)
    out["depth_filter_family_consistent"] = family_consistent
    disease_two_sided = s("disease_two_sided_support", False).fillna(False).astype(bool)
    control_two_sided = s("control_two_sided_support", False).fillna(False).astype(bool)
    disease_family_consistent = s("disease_family_agreement", 0).fillna(0).astype(int) == 1
    control_family_consistent = s("control_family_agreement", 0).fillna(0).astype(int) == 1
    disease_orientation_consistent = s("disease_strand_agreement", 0).fillna(0).astype(int) == 1
    control_orientation_consistent = s("control_strand_agreement", 0).fillna(0).astype(int) == 1
    two_sided_consistent = (
        (disease_two_sided & disease_family_consistent & disease_orientation_consistent)
        | (control_two_sided & control_family_consistent & control_orientation_consistent)
    )
    out["depth_filter_two_sided_consistent"] = two_sided_consistent
    silver_mask = s("silver_stage_pass", False).fillna(False).astype(bool)
    if silver_mask.any():
        depth_mask = silver_mask
    else:
        depth_mask = s("junk_flag_count", 999).fillna(999).astype(int) == 0
    out["depth_filter_pass"] = depth_mask
    idxs = out.index[depth_mask].tolist()
    if not idxs:
        click.echo("[mei-annotate] empirical stage skipped: no loci passed empirical prefilter")
        return out

    loci_for_empirical = out.loc[depth_mask, ["chrom", "window_start", "window_end"]].copy()
    cache_key = _empirical_cache_key(
        loci=loci_for_empirical,
        disease_bam_path=disease_bam_path,
        control_bam_path=control_bam_path,
        empirical_random_windows=empirical_random_windows,
        empirical_random_scope=empirical_random_scope,
        empirical_random_seed=empirical_random_seed,
        empirical_highconf_bed=empirical_highconf_bed,
        empirical_exclude_merged_bed=empirical_exclude_merged_bed,
        empirical_exclude_segdup_bed=empirical_exclude_segdup_bed,
        empirical_exclude_mappability_bedgraph=empirical_exclude_mappability_bedgraph,
        empirical_exclude_mappability_threshold=empirical_exclude_mappability_threshold,
        empirical_exclude_gap_bed=empirical_exclude_gap_bed,
        empirical_exclude_blacklist_bed=empirical_exclude_blacklist_bed,
    )
    random_disease_df = pd.DataFrame()
    random_control_df = pd.DataFrame()
    cache_hit = False
    if empirical_cache_dir is not None:
        empirical_cache_dir.mkdir(parents=True, exist_ok=True)
        disease_cache_path = empirical_cache_dir / f"{cache_key}.disease.parquet"
        control_cache_path = empirical_cache_dir / f"{cache_key}.control.parquet"
        if disease_cache_path.exists() and control_cache_path.exists():
            try:
                random_disease_df = pd.read_parquet(disease_cache_path)
                random_control_df = pd.read_parquet(control_cache_path)
                cache_hit = True
                click.echo(
                    f"[mei-annotate] empirical cache hit key={cache_key} "
                    f"rows={len(random_disease_df)}"
                )
            except Exception:
                cache_hit = False
        else:
            click.echo(f"[mei-annotate] empirical cache miss key={cache_key}")

    prep_t0 = time.monotonic()
    click.echo("[mei-annotate] empirical stage: preparing junk exclusion masks")
    merged_exclusion_ready = (
        empirical_exclude_merged_bed is not None and Path(empirical_exclude_merged_bed).exists()
    )
    junk_trees = {}
    if not cache_hit and not merged_exclusion_ready:
        junk_trees = _build_junk_interval_trees(
            segdup_bed=empirical_exclude_segdup_bed,
            low_mappability_bedgraph=empirical_exclude_mappability_bedgraph,
            low_mappability_threshold=empirical_exclude_mappability_threshold,
            gap_bed=empirical_exclude_gap_bed,
            encode_blacklist_bed=empirical_exclude_blacklist_bed,
        )
        junk_interval_count = sum(len(tree) for tree in junk_trees.values())
        click.echo(
            f"[mei-annotate] empirical stage: junk masks ready chroms={len(junk_trees)} intervals={junk_interval_count}"
        )
    elif not cache_hit and merged_exclusion_ready:
        click.echo(
            f"[mei-annotate] empirical stage: using merged exclusion bed {empirical_exclude_merged_bed}"
        )
    click.echo(
        f"[mei-annotate] empirical stage: exclusion mask prep elapsed={time.monotonic() - prep_t0:.1f}s"
    )

    with pysam.AlignmentFile(str(disease_bam_path), "rb") as disease_bam, pysam.AlignmentFile(
        str(control_bam_path), "rb"
    ) as control_bam:
        total_loci = int(len(idxs))
        loci_progress_every = 100
        click.echo(f"[mei-annotate] empirical stage: computing context metrics for {total_loci} loci")
        for i, idx in enumerate(idxs, start=1):
            row = out.loc[idx]
            chrom = str(row["chrom"])
            start = int(row["window_start"])
            end = int(row["window_end"])
            t_metrics = _context_quality_metrics_for_interval(disease_bam, chrom=chrom, start_1based=start, end_1based=end)
            n_metrics = _context_quality_metrics_for_interval(control_bam, chrom=chrom, start_1based=start, end_1based=end)
            out.at[idx, "disease_local_bam_mean_depth"] = float(t_metrics["local_bam_mean_depth"])
            out.at[idx, "control_local_bam_mean_depth"] = float(n_metrics["local_bam_mean_depth"])
            out.at[idx, "disease_local_bam_peak_depth"] = float(t_metrics.get("local_bam_peak_depth", 0.0))
            out.at[idx, "control_local_bam_peak_depth"] = float(n_metrics.get("local_bam_peak_depth", 0.0))
            out.at[idx, "disease_context_non_sv_reads"] = int(t_metrics["context_non_sv_reads"])
            out.at[idx, "control_context_non_sv_reads"] = int(n_metrics["context_non_sv_reads"])
            out.at[idx, "disease_context_mapq_mean"] = float(t_metrics["context_mapq_mean"])
            out.at[idx, "control_context_mapq_mean"] = float(n_metrics["context_mapq_mean"])
            out.at[idx, "disease_context_mapq_lt20_fraction"] = float(t_metrics["context_mapq_lt20_fraction"])
            out.at[idx, "control_context_mapq_lt20_fraction"] = float(n_metrics["context_mapq_lt20_fraction"])
            out.at[idx, "disease_context_nm_per_100bp_mean"] = float(t_metrics["context_nm_per_100bp_mean"])
            out.at[idx, "control_context_nm_per_100bp_mean"] = float(n_metrics["context_nm_per_100bp_mean"])
            out.at[idx, "disease_context_nm_per_100bp_p90"] = float(t_metrics["context_nm_per_100bp_p90"])
            out.at[idx, "control_context_nm_per_100bp_p90"] = float(n_metrics["context_nm_per_100bp_p90"])
            if i % loci_progress_every == 0 or i == total_loci:
                elapsed = time.monotonic() - t0
                click.echo(
                    f"[mei-annotate] empirical stage: locus metrics {i}/{total_loci} "
                    f"(elapsed={elapsed:.1f}s)"
                )

        if not cache_hit:
            click.echo("[mei-annotate] empirical stage: building random-window background metrics")
            random_windows = _sample_random_windows(
                candidates=out.loc[depth_mask].copy() if depth_mask.any() else out.copy(),
                bam=disease_bam,
                n_windows=int(empirical_random_windows),
                scope=str(empirical_random_scope),
                random_seed=int(empirical_random_seed),
                highconf_bed=empirical_highconf_bed,
                junk_trees=junk_trees,
                junk_exclusion_bed=empirical_exclude_merged_bed if merged_exclusion_ready else None,
            )
            click.echo(
                f"[mei-annotate] empirical stage: sampled {len(random_windows)} random windows "
                f"(scope={empirical_random_scope}, n={empirical_random_windows})"
            )
            random_disease_rows: list[dict[str, float | int | str]] = []
            random_control_rows: list[dict[str, float | int | str]] = []
            random_progress_every = 200
            total_random = int(len(random_windows))
            for i, rw in enumerate(random_windows.itertuples(index=False), start=1):
                chrom = str(rw.chrom)
                start = int(rw.window_start)
                end = int(rw.window_end)
                t_metrics = _context_quality_metrics_for_interval(
                    disease_bam, chrom=chrom, start_1based=start, end_1based=end
                )
                n_metrics = _context_quality_metrics_for_interval(
                    control_bam, chrom=chrom, start_1based=start, end_1based=end
                )
                random_disease_rows.append(
                    {
                        "chrom": chrom,
                        "local_bam_mean_depth": float(t_metrics["local_bam_mean_depth"]),
                        "context_mapq_mean": float(t_metrics["context_mapq_mean"]),
                        "context_mapq_lt20_fraction": float(t_metrics["context_mapq_lt20_fraction"]),
                        "context_nm_per_100bp_mean": float(t_metrics["context_nm_per_100bp_mean"]),
                        "context_nm_per_100bp_p90": float(t_metrics["context_nm_per_100bp_p90"]),
                    }
                )
                random_control_rows.append(
                    {
                        "chrom": chrom,
                        "local_bam_mean_depth": float(n_metrics["local_bam_mean_depth"]),
                        "context_mapq_mean": float(n_metrics["context_mapq_mean"]),
                        "context_mapq_lt20_fraction": float(n_metrics["context_mapq_lt20_fraction"]),
                        "context_nm_per_100bp_mean": float(n_metrics["context_nm_per_100bp_mean"]),
                        "context_nm_per_100bp_p90": float(n_metrics["context_nm_per_100bp_p90"]),
                    }
                )
                if i % random_progress_every == 0 or i == total_random:
                    elapsed = time.monotonic() - t0
                    click.echo(
                        f"[mei-annotate] empirical stage: random-window metrics {i}/{total_random} "
                        f"(elapsed={elapsed:.1f}s)"
                    )
            random_disease_df = pd.DataFrame(random_disease_rows)
            random_control_df = pd.DataFrame(random_control_rows)
            if empirical_cache_dir is not None:
                disease_cache_path = empirical_cache_dir / f"{cache_key}.disease.parquet"
                control_cache_path = empirical_cache_dir / f"{cache_key}.control.parquet"
                random_disease_df.to_parquet(disease_cache_path, index=False)
                random_control_df.to_parquet(control_cache_path, index=False)
                click.echo(
                    f"[mei-annotate] empirical cache write key={cache_key} "
                    f"rows={len(random_disease_df)}"
                )
        else:
            click.echo("[mei-annotate] empirical stage: using cached random-window metrics")

    t_mei = _df_col_series(out, "disease_mei_supported_reads", 0).astype(float)
    n_mei = _df_col_series(out, "control_mei_supported_reads", 0).astype(float)
    t_depth = out["disease_local_bam_mean_depth"].astype(float)
    n_depth = out["control_local_bam_mean_depth"].astype(float)
    out["disease_mei_support_per_100x_bam_depth"] = (t_mei * 100.0) / t_depth.replace(0, 1.0)
    out["control_mei_support_per_100x_bam_depth"] = (n_mei * 100.0) / n_depth.replace(0, 1.0)
    out["mei_support_per_100x_bam_depth_delta"] = (
        out["disease_mei_support_per_100x_bam_depth"] - out["control_mei_support_per_100x_bam_depth"]
    )
    out["mei_support_per_100x_bam_depth_ratio"] = (
        (out["disease_mei_support_per_100x_bam_depth"] + 1e-3)
        / (out["control_mei_support_per_100x_bam_depth"] + 1e-3)
    )

    # Empirical scoring should be applied only to the evaluated subset (depth_mask),
    # while leaving default neutral values for non-evaluated rows.
    metric_specs: list[tuple[str, str]] = [
        ("local_bam_mean_depth", "high"),
        ("context_mapq_mean", "low"),
        ("context_mapq_lt20_fraction", "high"),
        ("context_nm_per_100bp_mean", "high"),
        ("context_nm_per_100bp_p90", "high"),
    ]
    out["disease_empirical_random_n"] = 0
    out["control_empirical_random_n"] = 0
    for metric, tail in metric_specs:
        out[f"disease_empirical_{metric}_percentile"] = 0.0
        out[f"disease_empirical_{metric}_p_{tail}"] = 1.0
        out[f"control_empirical_{metric}_percentile"] = 0.0
        out[f"control_empirical_{metric}_p_{tail}"] = 1.0

    score_idx = out.index[depth_mask]
    if len(score_idx) > 0:
        disease_for_scoring = out.loc[score_idx, ["chrom"]].copy()
        disease_for_scoring["local_bam_mean_depth"] = out.loc[score_idx, "disease_local_bam_mean_depth"].astype(float)
        disease_for_scoring["context_mapq_mean"] = out.loc[score_idx, "disease_context_mapq_mean"].astype(float)
        disease_for_scoring["context_mapq_lt20_fraction"] = out.loc[score_idx, "disease_context_mapq_lt20_fraction"].astype(
            float
        )
        disease_for_scoring["context_nm_per_100bp_mean"] = out.loc[score_idx, "disease_context_nm_per_100bp_mean"].astype(
            float
        )
        disease_for_scoring["context_nm_per_100bp_p90"] = out.loc[score_idx, "disease_context_nm_per_100bp_p90"].astype(float)
        control_for_scoring = out.loc[score_idx, ["chrom"]].copy()
        control_for_scoring["local_bam_mean_depth"] = out.loc[score_idx, "control_local_bam_mean_depth"].astype(float)
        control_for_scoring["context_mapq_mean"] = out.loc[score_idx, "control_context_mapq_mean"].astype(float)
        control_for_scoring["context_mapq_lt20_fraction"] = out.loc[score_idx, "control_context_mapq_lt20_fraction"].astype(
            float
        )
        control_for_scoring["context_nm_per_100bp_mean"] = out.loc[score_idx, "control_context_nm_per_100bp_mean"].astype(
            float
        )
        control_for_scoring["context_nm_per_100bp_p90"] = out.loc[score_idx, "control_context_nm_per_100bp_p90"].astype(
            float
        )

        disease_scored = _apply_empirical_context_scores(
            loci_metrics=disease_for_scoring,
            random_metrics=random_disease_df,
            sample_prefix="disease",
            scope=str(empirical_random_scope),
            progress_every=200,
        )
        control_scored = _apply_empirical_context_scores(
            loci_metrics=control_for_scoring,
            random_metrics=random_control_df,
            sample_prefix="control",
            scope=str(empirical_random_scope),
            progress_every=200,
        )
        for col in disease_scored.columns:
            if col.startswith("disease_empirical_"):
                out.loc[score_idx, col] = disease_scored[col].values
        for col in control_scored.columns:
            if col.startswith("control_empirical_"):
                out.loc[score_idx, col] = control_scored[col].values
    click.echo("[mei-annotate] empirical stage: applying empirical p-value scoring complete")
    elapsed_total = time.monotonic() - t0
    click.echo(f"[mei-annotate] empirical stage complete (elapsed={elapsed_total:.1f}s)")
    click.echo(f"[mei-annotate] empirical stage walltime={time.monotonic() - stage_start:.1f}s")
    return out


def _add_local_depth_normalized_support(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    disease_total = _df_col_series(out, "disease_total_rows", 0).astype(float)
    control_total = _df_col_series(out, "control_total_rows", 0).astype(float)
    disease_mei = _df_col_series(out, "disease_mei_supported_reads", 0).astype(float)
    control_mei = _df_col_series(out, "control_mei_supported_reads", 0).astype(float)

    # Local informative depth proxy from candidate-building stage.
    out["disease_local_informative_rows"] = disease_total.fillna(0.0).astype(int)
    out["control_local_informative_rows"] = control_total.fillna(0.0).astype(int)

    out["disease_mei_support_local_frac"] = (disease_mei / disease_total.replace(0, 1)).fillna(0.0)
    out["control_mei_support_local_frac"] = (control_mei / control_total.replace(0, 1)).fillna(0.0)
    out["disease_mei_support_per_100_local_rows"] = out["disease_mei_support_local_frac"] * 100.0
    out["control_mei_support_per_100_local_rows"] = out["control_mei_support_local_frac"] * 100.0
    out["mei_local_support_frac_delta"] = (
        out["disease_mei_support_local_frac"] - out["control_mei_support_local_frac"]
    )
    out["mei_local_support_frac_ratio"] = (
        (out["disease_mei_support_local_frac"] + 1e-4) / (out["control_mei_support_local_frac"] + 1e-4)
    )
    return out


def _normal_ci_bounds_from_soft_counts(
    p: pd.Series,
    n_eff: pd.Series,
    z: float = 1.96,
) -> tuple[pd.Series, pd.Series]:
    # Heuristic uncertainty bounds for weighted-support VAF.
    n_pos = n_eff.astype(float).where(n_eff.astype(float) > 0.0)
    se = ((p * (1.0 - p)) / n_pos).pow(0.5)
    low = (p - (z * se)).clip(lower=0.0, upper=1.0)
    high = (p + (z * se)).clip(lower=0.0, upper=1.0)
    return low.where(n_pos.notna()), high.where(n_pos.notna())


def _add_heuristic_assembly_like_vaf_fields(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    s = lambda col, default: _df_col_series(out, col, default)

    sr_t = s("disease_split_mei_supported_reads", 0).fillna(0).astype(float)
    sr_n = s("control_split_mei_supported_reads", 0).fillna(0).astype(float)
    dpe_t = s("disease_discordant_mei_supported_reads", 0).fillna(0).astype(float)
    dpe_n = s("control_discordant_mei_supported_reads", 0).fillna(0).astype(float)

    # TODO(v2): replace this heuristic weighted model with RF-based TE genotyping/AF
    # inference (xTea-style feature model using SR/DRP/reference-support evidence).
    out["asm_disease_sr_alt_reads"] = sr_t
    out["asm_control_sr_alt_reads"] = sr_n
    out["asm_disease_dpe_alt_reads"] = dpe_t
    out["asm_control_dpe_alt_reads"] = dpe_n
    out["asm_disease_alt_soft_reads"] = sr_t + (0.5 * dpe_t)
    out["asm_control_alt_soft_reads"] = sr_n + (0.5 * dpe_n)
    out["asm_vaf_method"] = "heuristic_sr_plus_half_dpe_over_alt_plus_ref"

    if "disease_context_non_sv_reads" in out.columns and "control_context_non_sv_reads" in out.columns:
        out["asm_disease_ref_support_reads"] = out["disease_context_non_sv_reads"].fillna(0).astype(float)
        out["asm_control_ref_support_reads"] = out["control_context_non_sv_reads"].fillna(0).astype(float)
        out["asm_reference_support_source"] = "context_non_sv_reads"
    else:
        out["asm_disease_ref_support_reads"] = float("nan")
        out["asm_control_ref_support_reads"] = float("nan")
        out["asm_reference_support_source"] = "unavailable"

    disease_total = out["asm_disease_alt_soft_reads"] + out["asm_disease_ref_support_reads"]
    control_total = out["asm_control_alt_soft_reads"] + out["asm_control_ref_support_reads"]
    out["asm_disease_callable_reads"] = disease_total
    out["asm_control_callable_reads"] = control_total
    out["asm_disease_vaf"] = out["asm_disease_alt_soft_reads"] / disease_total.where(disease_total > 0.0)
    out["asm_control_vaf"] = out["asm_control_alt_soft_reads"] / control_total.where(control_total > 0.0)
    out["asm_vaf_delta"] = out["asm_disease_vaf"] - out["asm_control_vaf"]

    d_low, d_high = _normal_ci_bounds_from_soft_counts(
        out["asm_disease_vaf"].fillna(0.0),
        out["asm_disease_callable_reads"].fillna(0.0),
    )
    n_low, n_high = _normal_ci_bounds_from_soft_counts(
        out["asm_control_vaf"].fillna(0.0),
        out["asm_control_callable_reads"].fillna(0.0),
    )
    out["asm_disease_vaf_ci_low"] = d_low
    out["asm_disease_vaf_ci_high"] = d_high
    out["asm_control_vaf_ci_low"] = n_low
    out["asm_control_vaf_ci_high"] = n_high

    disease_width = (out["asm_disease_vaf_ci_high"] - out["asm_disease_vaf_ci_low"]).astype(float)
    control_width = (out["asm_control_vaf_ci_high"] - out["asm_control_vaf_ci_low"]).astype(float)
    out["assembly_confidence_score"] = (
        1.0 - ((disease_width.fillna(1.0) + control_width.fillna(1.0)) / 2.0)
    ).clip(lower=0.0, upper=1.0)

    silver_mask = s("silver_stage_pass", False).fillna(False).astype(bool)
    existing_status = out.get("asm_status", pd.Series([""] * len(out), index=out.index)).fillna("").astype(str)
    out["asm_status"] = existing_status
    out.loc[silver_mask & (out["asm_status"] == ""), "asm_status"] = "heuristic_estimated"
    no_ref = out["asm_reference_support_source"] == "unavailable"
    no_evidence = silver_mask & (
        (out["asm_disease_callable_reads"].fillna(0.0) <= 0.0)
        & (out["asm_control_callable_reads"].fillna(0.0) <= 0.0)
    )
    out.loc[silver_mask & no_ref & (out["asm_status"] == "heuristic_estimated"), "asm_status"] = (
        "heuristic_no_reference_support"
    )
    out.loc[no_evidence & (out["asm_status"].str.startswith("heuristic")), "asm_status"] = "heuristic_no_callable_reads"
    return out


def _assign_bronze_silver_stages(candidates: pd.DataFrame) -> pd.DataFrame:
    out = _ensure_candidate_schema_defaults(candidates)
    s = lambda col, default: _df_col_series(out, col, default)
    out["bronze_stage_pass"] = True

    junk_clean = s("junk_flag_count", 999).fillna(999).astype(int) == 0
    t_left_split = s("disease_L_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    t_right_split = s("disease_R_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    t_left_disc = s("disease_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1
    t_right_disc = s("disease_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1
    n_left_split = s("control_L_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    n_right_split = s("control_R_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    n_left_disc = s("control_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1
    n_right_disc = s("control_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1

    disease_bilateral_any = (t_left_split | t_left_disc) & (t_right_split | t_right_disc)
    control_bilateral_any = (n_left_split | n_left_disc) & (n_right_split | n_right_disc)
    out["silver_bilateral_support_any"] = disease_bilateral_any | control_bilateral_any
    t_left_poly = s("disease_L_poly_at_reads", 0).fillna(0).astype(float) >= 1
    t_right_poly = s("disease_R_poly_at_reads", 0).fillna(0).astype(float) >= 1
    n_left_poly = s("control_L_poly_at_reads", 0).fillna(0).astype(float) >= 1
    n_right_poly = s("control_R_poly_at_reads", 0).fillna(0).astype(float) >= 1

    disease_split_consistent = (
        s("disease_two_sided_support", False).fillna(False).astype(bool)
        & (s("disease_family_agreement", 0).fillna(0).astype(int) == 1)
        & (s("disease_strand_agreement", 0).fillna(0).astype(int) == 1)
    )
    control_split_consistent = (
        s("control_two_sided_support", False).fillna(False).astype(bool)
        & (s("control_family_agreement", 0).fillna(0).astype(int) == 1)
        & (s("control_strand_agreement", 0).fillna(0).astype(int) == 1)
    )

    disease_disc_consistent = (
        s("disease_discordant_mei_two_sided_support", False).fillna(False).astype(bool)
        & (s("disease_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & (s("disease_discordant_mei_strand_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & s("disease_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
        & s("disease_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    control_disc_consistent = (
        (s("control_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1)
        & (s("control_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1)
        & (s("control_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & (s("control_discordant_mei_strand_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & s("control_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
        & s("control_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )

    event_family_consistent = s("event_family_consistent", False).fillna(False).astype(bool)
    event_strand_consistent = s("event_strand_consistent", False).fillna(False).astype(bool)

    # PolyA-rescue bilateral support:
    # one side has MEI anchor support and the opposite side has polyA-clipped support.
    # If orientation is known, enforce expected tail side:
    # + insertion => right-side polyA; - insertion => left-side polyA.
    disease_ori = s("disease_insertion_orientation", "").fillna("").astype(str)
    control_ori = s("control_discordant_mei_strand", "").fillna("").astype(str)
    t_poly_mei_any = (t_left_poly & (t_right_split | t_right_disc)) | (t_right_poly & (t_left_split | t_left_disc))
    n_poly_mei_any = (n_left_poly & (n_right_split | n_right_disc)) | (n_right_poly & (n_left_split | n_left_disc))
    t_poly_oriented = (
        ((disease_ori == "+") & t_right_poly & (t_left_split | t_left_disc))
        | ((disease_ori == "-") & t_left_poly & (t_right_split | t_right_disc))
    )
    n_poly_oriented = (
        ((control_ori == "+") & n_right_poly & (n_left_split | n_left_disc))
        | ((control_ori == "-") & n_left_poly & (n_right_split | n_right_disc))
    )
    poly_sidepair_support = (t_poly_mei_any & ((disease_ori == "") | t_poly_oriented)) | (
        n_poly_mei_any & ((control_ori == "") | n_poly_oriented)
    )
    out["silver_polyA_sidepair_support"] = poly_sidepair_support

    t_left_anchor_complex = s("disease_discordant_anchor_left_complex_side", False).fillna(False).astype(bool)
    t_right_anchor_complex = s("disease_discordant_anchor_right_complex_side", False).fillna(False).astype(bool)
    n_left_anchor_complex = s("control_discordant_anchor_left_complex_side", False).fillna(False).astype(bool)
    n_right_anchor_complex = s("control_discordant_anchor_right_complex_side", False).fillna(False).astype(bool)
    t_left_structural = t_left_split | t_left_disc | t_left_poly | t_left_anchor_complex
    t_right_structural = t_right_split | t_right_disc | t_right_poly | t_right_anchor_complex
    n_left_structural = n_left_split | n_left_disc | n_left_poly | n_left_anchor_complex
    n_right_structural = n_right_split | n_right_disc | n_right_poly | n_right_anchor_complex
    disease_bilateral_structural = t_left_structural & t_right_structural
    control_bilateral_structural = n_left_structural & n_right_structural
    out["silver_bilateral_structural_support"] = disease_bilateral_structural | control_bilateral_structural

    disease_complex_sidepair = (
        (t_left_split | t_left_disc) & (t_right_anchor_complex | t_right_poly)
    ) | ((t_right_split | t_right_disc) & (t_left_anchor_complex | t_left_poly))
    control_complex_sidepair = (
        (n_left_split | n_left_disc) & (n_right_anchor_complex | n_right_poly)
    ) | ((n_right_split | n_right_disc) & (n_left_anchor_complex | n_left_poly))
    out["silver_complex_sidepair_support"] = disease_complex_sidepair | control_complex_sidepair
    out["silver_complex_structural_consistent"] = (
        out["silver_bilateral_structural_support"]
        & out["silver_complex_sidepair_support"]
        & (
            s("disease_mei_with_complex_sidepair", False).fillna(False).astype(bool)
            | s("control_mei_with_complex_sidepair", False).fillna(False).astype(bool)
            | s("mei_with_complex_sv_signature", False).fillna(False).astype(bool)
        )
    )

    silver_consistency = (
        disease_split_consistent | control_split_consistent | disease_disc_consistent | control_disc_consistent
    )
    out["silver_consistency_pass"] = (
        silver_consistency
        | (event_family_consistent & event_strand_consistent)
        | (poly_sidepair_support & event_family_consistent)
        | out["silver_complex_structural_consistent"]
    )
    out["silver_discordant_two_sided_consistent"] = disease_disc_consistent | control_disc_consistent

    disease_l_bp = s("disease_L_mei_breakpoint_mode", 0).fillna(0).astype(int)
    disease_r_bp = s("disease_R_mei_breakpoint_mode", 0).fillna(0).astype(int)
    control_l_bp = s("control_L_mei_breakpoint_mode", 0).fillna(0).astype(int)
    control_r_bp = s("control_R_mei_breakpoint_mode", 0).fillna(0).astype(int)
    out["silver_split_breakpoint_resolved"] = (
        (t_left_split & (disease_l_bp > 0))
        | (t_right_split & (disease_r_bp > 0))
        | (n_left_split & (control_l_bp > 0))
        | (n_right_split & (control_r_bp > 0))
    )
    out["silver_breakpoint_interval_resolved"] = (
        pd.to_numeric(s("insertion_breakpoint_interval_width_bp", 0), errors="coerce").fillna(0).astype(int) > 0
    )
    out["silver_insertion_span_resolved"] = s("insertion_mei_span", 0).fillna(0).astype(int) > 0
    out["silver_breakpoint_or_span_resolved"] = (
        out["silver_split_breakpoint_resolved"]
        | out["silver_breakpoint_interval_resolved"]
        | out["silver_insertion_span_resolved"]
    )
    disease_family_consistent = (
        (s("disease_family_agreement", 0).fillna(0).astype(int) == 1)
        | (s("disease_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
    )
    control_family_consistent = (
        (s("control_family_agreement", 0).fillna(0).astype(int) == 1)
        | (s("control_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
    )
    # Asymmetric insertion-like support:
    # one side has MEI anchor support, opposite side has split/discordant support
    # (or orientation-consistent polyA/T), with per-sample family/self-consistency.
    # This path intentionally does not require resolved breakpoint/span because
    # real MEI events can present as one-sided anchor + opposite polyA/T/DPE.
    disease_left_anchor = t_left_split | t_left_disc
    disease_right_anchor = t_right_split | t_right_disc
    disease_left_opposite = t_left_split | t_left_disc | (((disease_ori == "-") | (disease_ori == "")) & t_left_poly)
    disease_right_opposite = t_right_split | t_right_disc | (((disease_ori == "+") | (disease_ori == "")) & t_right_poly)
    disease_self_consistent = s("disease_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    disease_asymmetric = (
        ((disease_left_anchor & disease_right_opposite) | (disease_right_anchor & disease_left_opposite))
        & disease_family_consistent
        & disease_self_consistent
    )
    control_left_anchor = n_left_split | n_left_disc
    control_right_anchor = n_right_split | n_right_disc
    control_left_opposite = n_left_split | n_left_disc | (((control_ori == "-") | (control_ori == "")) & n_left_poly)
    control_right_opposite = n_right_split | n_right_disc | (((control_ori == "+") | (control_ori == "")) & n_right_poly)
    control_self_consistent = s("control_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    control_asymmetric = (
        ((control_left_anchor & control_right_opposite) | (control_right_anchor & control_left_opposite))
        & control_family_consistent
        & control_self_consistent
    )
    out["silver_asymmetric_insertion_like"] = disease_asymmetric | control_asymmetric
    out["silver_consistency_pass"] = out["silver_consistency_pass"] | out["silver_asymmetric_insertion_like"]

    disease_any_mei = (
        s("disease_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    ) | (s("disease_full_mei_supported_reads", 0).fillna(0).astype(float) >= 1)
    control_any_mei = (
        s("control_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    ) | (s("control_full_mei_supported_reads", 0).fillna(0).astype(float) >= 1)
    out["silver_any_mei_support"] = disease_any_mei | control_any_mei

    out["silver_stage_pass"] = junk_clean & out["silver_any_mei_support"] & (
        (
            out["silver_bilateral_support_any"]
            | poly_sidepair_support
            | out["silver_complex_structural_consistent"]
            | out["silver_breakpoint_or_span_resolved"]
            | out["silver_asymmetric_insertion_like"]
        )
        & (out["silver_consistency_pass"] | out["silver_breakpoint_or_span_resolved"])
    )
    silver_fail = ~out["silver_stage_pass"]
    has_support_or_resolution = (
        out["silver_bilateral_support_any"]
        | poly_sidepair_support
        | out["silver_complex_structural_consistent"]
        | out["silver_breakpoint_or_span_resolved"]
        | out["silver_asymmetric_insertion_like"]
    )
    has_consistency_or_resolution = out["silver_consistency_pass"] | out["silver_breakpoint_or_span_resolved"]
    silver_fail_reason = pd.Series("", index=out.index, dtype="object")

    def _append_reason(mask: pd.Series, reason: str) -> None:
        target = mask.fillna(False).astype(bool)
        if not target.any():
            return
        empty = silver_fail_reason.eq("")
        silver_fail_reason.loc[target & empty] = reason
        silver_fail_reason.loc[target & ~empty] = silver_fail_reason.loc[target & ~empty] + ";" + reason

    _append_reason(silver_fail & (~junk_clean), "junk_region_flagged")
    _append_reason(silver_fail & (~out["silver_any_mei_support"]), "no_mei_signal")
    _append_reason(silver_fail & (~has_support_or_resolution), "no_bilateral_or_breakpoint_span_support")
    _append_reason(silver_fail & (~has_consistency_or_resolution), "no_consistency_signal")
    out["silver_stage_fail_reason"] = ""
    out.loc[silver_fail, "silver_stage_fail_reason"] = silver_fail_reason.loc[silver_fail]
    out["stage_fail_reason"] = out["silver_stage_fail_reason"].fillna("").astype(str)

    out["analysis_stage_tier"] = "bronze"
    out.loc[out["silver_stage_pass"], "analysis_stage_tier"] = "silver"
    click.echo(
        "[mei-annotate] stage counts "
        f"bronze={len(out)} silver={int(out['silver_stage_pass'].sum())}"
    )
    return out


def _mei_mapped_from_support_string(series: pd.Series) -> pd.Series:
    """Parse ``MEI_MAPPED=N`` from supporting-reads summary strings."""
    text = series.fillna("").astype(str)
    vals = pd.to_numeric(text.str.extract(r"MEI_MAPPED=([0-9]+)", expand=False), errors="coerce")
    return vals.fillna(0).astype(int)


def _assign_gold_stage(
    candidates: pd.DataFrame,
    empirical_p_threshold: float = 0.001,
    empirical_stage: bool = False,
    min_mei_mapped: int = 3,
) -> pd.DataFrame:
    out = candidates.copy()
    out["gold_empirical_p_threshold"] = float(empirical_p_threshold)
    out["gold_empirical_eval_available"] = False
    out["gold_empirical_outlier"] = False
    out["gold_stage_pass"] = False
    out["gold_stage_fail_reason"] = ""
    out["gold_min_mei_mapped"] = int(min_mei_mapped)

    p_cols = [
        "disease_empirical_local_bam_mean_depth_p_high",
        "disease_empirical_context_mapq_mean_p_low",
        "disease_empirical_context_mapq_lt20_fraction_p_high",
        "disease_empirical_context_nm_per_100bp_mean_p_high",
        "disease_empirical_context_nm_per_100bp_p90_p_high",
        "control_empirical_local_bam_mean_depth_p_high",
        "control_empirical_context_mapq_mean_p_low",
        "control_empirical_context_mapq_lt20_fraction_p_high",
        "control_empirical_context_nm_per_100bp_mean_p_high",
        "control_empirical_context_nm_per_100bp_p90_p_high",
    ]
    available_cols = [c for c in p_cols if c in out.columns] if empirical_stage else []
    silver = _df_col_series(out, "silver_stage_pass", False).fillna(False).astype(bool)
    if available_cols:
        out["gold_empirical_eval_available"] = True
        pvals = out.loc[:, available_cols].fillna(1.0).astype(float)
        out["gold_empirical_outlier"] = (pvals < float(empirical_p_threshold)).any(axis=1)
        out["gold_stage_pass"] = silver & (~out["gold_empirical_outlier"])
        out.loc[silver & out["gold_empirical_outlier"], "gold_stage_fail_reason"] = "empirical_outlier"
    else:
        out["gold_stage_pass"] = silver
        if empirical_stage:
            out.loc[silver, "gold_stage_fail_reason"] = "empirical_not_available"

    # Require enough MEI-mapped support in at least one sample. Low MEI_MAPPED
    # silver calls are dominated by disease/control-only noise and flood review plots.
    disease_mei = pd.to_numeric(_df_col_series(out, "disease_mei_mapped", float("nan")), errors="coerce")
    control_mei = pd.to_numeric(_df_col_series(out, "control_mei_mapped", float("nan")), errors="coerce")
    if disease_mei.isna().all():
        disease_mei = _mei_mapped_from_support_string(_df_col_series(out, "disease_supporting_reads", ""))
    else:
        disease_mei = disease_mei.fillna(0)
    if control_mei.isna().all():
        control_mei = _mei_mapped_from_support_string(_df_col_series(out, "control_supporting_reads", ""))
    else:
        control_mei = control_mei.fillna(0)
    disease_mei = disease_mei.astype(int)
    control_mei = control_mei.astype(int)
    out["disease_mei_mapped"] = disease_mei
    out["control_mei_mapped"] = control_mei
    mei_mapped_ok = (disease_mei >= int(min_mei_mapped)) | (control_mei >= int(min_mei_mapped))
    low_mei = silver & out["gold_stage_pass"] & (~mei_mapped_ok)
    if low_mei.any():
        out.loc[low_mei, "gold_stage_pass"] = False
        prev = _df_col_series(out, "gold_stage_fail_reason", "").fillna("").astype(str)
        fail_tag = f"mei_mapped_lt_{int(min_mei_mapped)}"
        need_append = low_mei & prev.ne("")
        need_set = low_mei & prev.eq("")
        out.loc[need_set, "gold_stage_fail_reason"] = fail_tag
        out.loc[need_append, "gold_stage_fail_reason"] = prev.loc[need_append] + ";" + fail_tag

    # COMPLEX_INS is a non-MEI structural insertion class — keep out of gold MEI review.
    complex_ins = (
        _df_col_series(out, "insertion_event_class", "").fillna("").astype(str).eq("COMPLEX_INS")
    )
    complex_ins_gold = silver & out["gold_stage_pass"] & complex_ins
    if complex_ins_gold.any():
        out.loc[complex_ins_gold, "gold_stage_pass"] = False
        prev = _df_col_series(out, "gold_stage_fail_reason", "").fillna("").astype(str)
        fail_tag = "complex_ins_non_mei"
        need_append = complex_ins_gold & prev.ne("")
        need_set = complex_ins_gold & prev.eq("")
        out.loc[need_set, "gold_stage_fail_reason"] = fail_tag
        out.loc[need_append, "gold_stage_fail_reason"] = prev.loc[need_append] + ";" + fail_tag

    # Non-empirical depth-outlier guard for obvious pileup artifacts.
    # Use a run-adaptive 3-sigma threshold and keep known overlaps exempt.
    known_poly = _df_col_series(out, "known_mei_polymorphism", False).fillna(False).astype(bool)
    d_depth_peak = pd.to_numeric(_df_col_series(out, "disease_local_bam_peak_depth", 0.0), errors="coerce").fillna(0.0)
    c_depth_peak = pd.to_numeric(_df_col_series(out, "control_local_bam_peak_depth", 0.0), errors="coerce").fillna(0.0)
    max_depth = pd.concat([d_depth_peak, c_depth_peak], axis=1).max(axis=1)
    # Fallback for older outputs that may not carry peak-depth fields.
    if float(max_depth.max()) <= 0.0:
        d_depth_mean = pd.to_numeric(_df_col_series(out, "disease_local_bam_mean_depth", 0.0), errors="coerce").fillna(0.0)
        c_depth_mean = pd.to_numeric(_df_col_series(out, "control_local_bam_mean_depth", 0.0), errors="coerce").fillna(0.0)
        max_depth = pd.concat([d_depth_mean, c_depth_mean], axis=1).max(axis=1)
    depth_ref = max_depth.loc[silver & max_depth.gt(0.0)]
    if depth_ref.empty:
        depth_ref = max_depth.loc[max_depth.gt(0.0)]
    depth_mean = float(depth_ref.mean()) if not depth_ref.empty else 0.0
    depth_sigma = float(depth_ref.std(ddof=0)) if not depth_ref.empty else 0.0
    if depth_sigma > 1e-6:
        depth_z = (max_depth - depth_mean) / depth_sigma
        depth_outlier = depth_z >= 3.0
    else:
        depth_outlier = pd.Series(False, index=out.index)
    # Depth-only artifact gate (user requested): reject extreme peak-depth
    # outliers among silver loci, except known polymorphism overlaps.
    depth_pileup_artifact = (
        silver
        & (~known_poly)
        & depth_outlier
    )
    if depth_pileup_artifact.any():
        out.loc[depth_pileup_artifact, "gold_stage_pass"] = False
        prev = _df_col_series(out, "gold_stage_fail_reason", "").fillna("").astype(str)
        fail_tag = "depth_pileup_artifact"
        need_append = depth_pileup_artifact & prev.ne("")
        need_set = depth_pileup_artifact & prev.eq("")
        out.loc[need_set, "gold_stage_fail_reason"] = fail_tag
        out.loc[need_append, "gold_stage_fail_reason"] = prev.loc[need_append] + ";" + fail_tag

    stage_fail_reason = _df_col_series(out, "stage_fail_reason", "").fillna("").astype(str)
    gold_fail_reason = _df_col_series(out, "gold_stage_fail_reason", "").fillna("").astype(str)
    silver_failed_gold = silver & (~out["gold_stage_pass"]) & gold_fail_reason.ne("")
    if silver_failed_gold.any():
        idx = silver_failed_gold[silver_failed_gold].index
        prefixed_gold_reason = "gold:" + gold_fail_reason.loc[idx]
        empty = stage_fail_reason.loc[idx].eq("")
        stage_fail_reason.loc[idx[empty]] = prefixed_gold_reason.loc[idx[empty]]
        stage_fail_reason.loc[idx[~empty]] = stage_fail_reason.loc[idx[~empty]] + ";" + prefixed_gold_reason.loc[idx[~empty]]
    stage_fail_reason.loc[out["gold_stage_pass"]] = ""
    out["stage_fail_reason"] = stage_fail_reason

    out.loc[out["gold_stage_pass"], "analysis_stage_tier"] = "gold"
    click.echo(
        "[mei-annotate] stage counts "
        f"silver={int(silver.sum())} gold={int(out['gold_stage_pass'].sum())} "
        f"(min_mei_mapped>={int(min_mei_mapped)})"
    )
    return out


def _two_sided_support_mask(df: pd.DataFrame) -> pd.Series:
    required_cols = [
        "disease_left_supported_reads",
        "disease_right_supported_reads",
        "control_left_supported_reads",
        "control_right_supported_reads",
    ]
    if (not any(col in df.columns for col in required_cols)) and ("two_sided_support" in df.columns):
        return _df_col_series(df, "two_sided_support", False).fillna(False).astype(bool)
    disease_left = _df_col_series(df, "disease_left_supported_reads", 0).fillna(0).astype(int)
    disease_right = _df_col_series(df, "disease_right_supported_reads", 0).fillna(0).astype(int)
    control_left = _df_col_series(df, "control_left_supported_reads", 0).fillna(0).astype(int)
    control_right = _df_col_series(df, "control_right_supported_reads", 0).fillna(0).astype(int)
    bilateral = ((disease_left >= 1) & (disease_right >= 1)) | ((control_left >= 1) & (control_right >= 1))
    if "silver_bilateral_support_any" in df.columns:
        bilateral = bilateral | df["silver_bilateral_support_any"].fillna(False).astype(bool)
    return bilateral


def _poly_at_supported_mask(df: pd.DataFrame) -> pd.Series:
    return (_df_col_series(df, "poly_at_reads", 0).fillna(0).astype(int) > 0) | (
        _df_col_series(df, "poly_at_max_run", 0).fillna(0).astype(int) > 0
    )


def _prioritize_mei_candidates(candidates: pd.DataFrame, *, stage_first: bool = True) -> pd.DataFrame:
    """Rank loci by evidence strength for manual review."""
    out = candidates.copy()
    out["two_sided_support"] = _two_sided_support_mask(out)
    out["poly_at_supported"] = _poly_at_supported_mask(out)

    def _extract_support_total(series: pd.Series, label: str) -> pd.Series:
        text = series.fillna("").astype(str)
        vals = pd.to_numeric(text.str.extract(rf"{label}=([0-9]+)", expand=False), errors="coerce").fillna(0)
        return vals.astype(int)

    if "consensus_tsd_detected" in out.columns:
        tsd_signal = _df_col_series(out, "consensus_tsd_detected", False).fillna(False).astype(bool)
    elif "tsd_detected" in out.columns:
        tsd_signal = _df_col_series(out, "tsd_detected", False).fillna(False).astype(bool)
    else:
        tsd_signal = _df_col_series(out, "consensus_tsd_len_estimate", 0).fillna(0).astype(float) >= 4.0
    out["_prio_tsd"] = tsd_signal
    out["_prio_poly_at"] = (
        _df_col_series(out, "consensus_poly_at_supported", False).fillna(False).astype(bool)
        | out["poly_at_supported"].astype(bool)
    )
    disease_split_total = (
        _df_col_series(out, "disease_split_mei_supported_reads", float("nan")).astype(float)
    )
    control_split_total = (
        _df_col_series(out, "control_split_mei_supported_reads", float("nan")).astype(float)
    )
    if disease_split_total.isna().all():
        disease_sr = _extract_support_total(_df_col_series(out, "disease_supporting_reads", ""), "SR_L") + _extract_support_total(
            _df_col_series(out, "disease_supporting_reads", ""), "SR_R"
        )
    else:
        disease_sr = disease_split_total.fillna(0).astype(int)
    if control_split_total.isna().all():
        control_sr = _extract_support_total(_df_col_series(out, "control_supporting_reads", ""), "SR_L") + _extract_support_total(
            _df_col_series(out, "control_supporting_reads", ""), "SR_R"
        )
    else:
        control_sr = control_split_total.fillna(0).astype(int)
    out["_prio_split_reads_max"] = pd.concat([disease_sr, control_sr], axis=1).max(axis=1).astype(int)

    disease_disc_total = _df_col_series(out, "disease_discordant_mei_supported_reads", float("nan")).astype(float)
    control_disc_total = _df_col_series(out, "control_discordant_mei_supported_reads", float("nan")).astype(float)
    if disease_disc_total.isna().all():
        disease_dpe = _extract_support_total(_df_col_series(out, "disease_supporting_reads", ""), "DPE_L") + _extract_support_total(
            _df_col_series(out, "disease_supporting_reads", ""), "DPE_R"
        )
    else:
        disease_dpe = disease_disc_total.fillna(0).astype(int)
    if control_disc_total.isna().all():
        control_dpe = _extract_support_total(_df_col_series(out, "control_supporting_reads", ""), "DPE_L") + _extract_support_total(
            _df_col_series(out, "control_supporting_reads", ""), "DPE_R"
        )
    else:
        control_dpe = control_disc_total.fillna(0).astype(int)
    out["_prio_discordant_reads_max"] = pd.concat([disease_dpe, control_dpe], axis=1).max(axis=1).astype(int)
    disease_mei_mapped = _extract_support_total(_df_col_series(out, "disease_supporting_reads", ""), "MEI_MAPPED")
    control_mei_mapped = _extract_support_total(_df_col_series(out, "control_supporting_reads", ""), "MEI_MAPPED")
    out["_prio_mei_mapped_max"] = pd.concat([disease_mei_mapped, control_mei_mapped], axis=1).max(axis=1).astype(int)
    disease_polya_mapped = _extract_support_total(_df_col_series(out, "disease_supporting_reads", ""), "polyA_MAPPED")
    control_polya_mapped = _extract_support_total(_df_col_series(out, "control_supporting_reads", ""), "polyA_MAPPED")
    out["_prio_polya_mapped_max"] = pd.concat([disease_polya_mapped, control_polya_mapped], axis=1).max(axis=1).astype(int)
    disease_vntr_mapped = _extract_support_total(_df_col_series(out, "disease_supporting_reads", ""), "VNTR_MAPPED")
    control_vntr_mapped = _extract_support_total(_df_col_series(out, "control_supporting_reads", ""), "VNTR_MAPPED")
    vntr_raw = pd.concat([disease_vntr_mapped, control_vntr_mapped], axis=1).max(axis=1).astype(int)
    # VNTR rescue is biologically meaningful for SVA only; ignore elsewhere.
    fam = _df_col_series(out, "consensus_mei_family", "").fillna("").astype(str).str.upper()
    out["_prio_vntr_mapped_max"] = vntr_raw.where(fam.eq("SVA"), 0).astype(int)
    # Prefer true MEI-mapped support over raw SR/DPE pileups (which can be complex SVs).
    # MEI_MAPPED dominates; polyA/VNTR rescues next; SR/DPE are secondary tie-breakers.
    out["read_support_heuristic_score"] = (
        1.00 * (out["_prio_mei_mapped_max"].astype(float).map(math.log1p))
        + 0.50 * (out["_prio_polya_mapped_max"].astype(float).map(math.log1p))
        + 0.50 * (out["_prio_vntr_mapped_max"].astype(float).map(math.log1p))
        + 0.20 * (out["_prio_split_reads_max"].astype(float).map(math.log1p))
        + 0.10 * (out["_prio_discordant_reads_max"].astype(float).map(math.log1p))
    )
    # Down-rank low-complexity pileup artifacts:
    # - strong A/T-rich/polyA signatures,
    # - weak cross-read coherence,
    # - very high support pileups that can arise in noisy simple-repeat contexts.
    # Keep known polymorphisms exempt from this penalty.
    tsd_seq_s = _df_col_series(out, "tsd_seq", "").fillna("").astype(str).str.upper()
    tsd_len_s = tsd_seq_s.str.len().astype(int)
    tsd_at_fraction = (
        (tsd_seq_s.str.count("A") + tsd_seq_s.str.count("T")) / tsd_len_s.replace(0, pd.NA)
    ).fillna(0.0).astype(float)
    tsd_longest_at_run = tsd_seq_s.str.findall(r"[AT]+").map(lambda runs: max((len(run) for run in runs), default=0)).astype(int)
    tsd_poly_filtered = _df_col_series(out, "tsd_poly_at_filter_applied", False).fillna(False).astype(bool)
    poly_run = _df_col_series(out, "poly_at_max_run", 0).fillna(0).astype(int)
    poly_reads = _df_col_series(out, "poly_at_reads", 0).fillna(0).astype(int)
    coherence = _df_col_series(out, "coherence_score", 0.0).fillna(0.0).astype(float)
    support_pileup = (out["_prio_split_reads_max"] + out["_prio_discordant_reads_max"]).astype(int)
    complex_signal = _df_col_series(out, "complex_sv_signal_score", 0.0).fillna(0.0).astype(float)
    known_poly = _df_col_series(out, "known_mei_polymorphism", False).fillna(False).astype(bool)
    at_rich_tsd_like = tsd_len_s.ge(8) & tsd_at_fraction.ge(0.85) & tsd_longest_at_run.ge(6)
    low_complexity_signature = tsd_poly_filtered | at_rich_tsd_like
    out["low_complexity_noisy_artifact_flag"] = (
        (~known_poly)
        & (coherence < 0.45)
        & (poly_run >= 30)
        & (poly_reads >= 8)
        & (support_pileup >= 35)
        & (low_complexity_signature | (complex_signal >= 0.80))
    )
    if "consensus_insertion_breakpoint_pos" in out.columns:
        bp_consensus = _df_col_series(out, "consensus_insertion_breakpoint_pos", 0)
    else:
        bp_consensus = _df_col_series(out, "insertion_breakpoint_pos", 0)
    out["_prio_breakpoint_consensus_available"] = bp_consensus.fillna(0).astype(int) > 0
    out["_prio_insertion_model_score"] = pd.to_numeric(
        _df_col_series(out, "insertion_model_score", 0.0), errors="coerce"
    ).fillna(0.0)
    out["_prio_clip_overlap_consistency"] = pd.to_numeric(
        _df_col_series(out, "event_clip_overlap_consistency", 0.0), errors="coerce"
    ).fillna(0.0)

    sort_cols: list[str] = []
    ascending: list[bool] = []
    if stage_first:
        if "gold_stage_pass" in out.columns and "silver_stage_pass" in out.columns:
            for col in ("gold_stage_pass", "silver_stage_pass"):
                sort_cols.append(col)
                ascending.append(False)
        elif "analysis_stage_tier" in out.columns:
            out["_prio_stage"] = (
                out["analysis_stage_tier"]
                .fillna("")
                .astype(str)
                .map({"gold": 0, "silver": 1, "bronze": 2})
                .fillna(3)
                .astype(int)
            )
            sort_cols.append("_prio_stage")
            ascending.append(True)
    sort_cols.extend(
        [
            "low_complexity_noisy_artifact_flag",
            "_prio_mei_mapped_max",
            "_prio_polya_mapped_max",
            "_prio_vntr_mapped_max",
            "read_support_heuristic_score",
            "_prio_split_reads_max",
            "_prio_discordant_reads_max",
            "_prio_tsd",
            "_prio_poly_at",
            "_prio_breakpoint_consensus_available",
            "_prio_clip_overlap_consistency",
            "_prio_insertion_model_score",
        ]
    )
    # Artifact flag ascending (False first), then MEI / polyA / VNTR(SVA) support.
    ascending.extend([True] + [False] * 11)
    sorted_out = out.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    return sorted_out.drop(
        columns=[c for c in sorted_out.columns if c.startswith("_prio_")],
        errors="ignore",
    ).reset_index(drop=True)


_YYRRRR_MT_ADJ_REPORT_MIN = 0.0  # report motif fields only when MT-adjusted log-odds are positive.


def _yyrrrr_mt_adj_value(row: pd.Series) -> float:
    val = row.get("breakpoint_yyrrrr_logodds_shift1_mt_adj", float("nan"))
    if pd.isna(val):
        return float("nan")
    try:
        return float(val)
    except (TypeError, ValueError):
        return float("nan")


def _yyrrrr_mt_adj_reportable(row: pd.Series) -> bool:
    val = _yyrrrr_mt_adj_value(row)
    return not pd.isna(val) and val > _YYRRRR_MT_ADJ_REPORT_MIN


def _apply_breakpoint_motif_report_gating(df: pd.DataFrame) -> pd.DataFrame:
    """Mask motif report fields unless breakpoint log-odds pass the report threshold."""
    out = df.copy()
    observed_hex = _df_col_series(out, "breakpoint_l1_en_hexamer_oriented", "").fillna("").astype(str)
    observed_pattern = _df_col_series(out, "breakpoint_l1_en_pattern_yy_rrrr", "").fillna("").astype(str)
    out["breakpoint_l1_en_observed_motif"] = observed_hex
    out["breakpoint_l1_en_observed_motif_pattern"] = observed_pattern
    mt_adj = _df_col_series(out, "breakpoint_yyrrrr_logodds_shift1_mt_adj", float("nan")).astype(float)
    reportable = mt_adj.notna() & (mt_adj > _YYRRRR_MT_ADJ_REPORT_MIN)
    # Keep raw observed breakpoint motif fields visible; only gate derived interpretation fields.
    for col in (
        "breakpoint_l1_en_best_motif",
        "breakpoint_l1_en_motif_type",
    ):
        if col not in out.columns:
            continue
        out.loc[~reportable, col] = ""
    return out


def _consensus_retrotransposition_class(row: pd.Series) -> str:
    if not _yyrrrr_mt_adj_reportable(row):
        return ""
    if bool(row.get("breakpoint_l1_en_motif_like", False)):
        motif_type = str(row.get("breakpoint_l1_en_motif_type", "") or "").strip()
        if motif_type == "l1_en_canonical":
            return "classical"
        if motif_type in {"l1_en_alternative", "nested_novel_like"}:
            return "non_classical"
    return "classical"


def _consensus_sequence_signature(row: pd.Series, *, retro_class: str = "") -> str:
    if not retro_class:
        return ""
    for col in (
        "breakpoint_l1_en_observed_motif_pattern",
        "breakpoint_l1_en_observed_motif",
        "breakpoint_l1_en_pattern_yy_rrrr",
        "breakpoint_l1_en_best_match_pattern_yy_rrrr",
        "breakpoint_l1_en_best_motif",
    ):
        pattern = str(row.get(col, "") or "").strip()
        if pattern:
            return pattern
    return ""


def _annotate_consensus_retrotransposition_fields(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    classes = out.apply(_consensus_retrotransposition_class, axis=1)
    out["consensus_retrotransposition_class"] = classes
    out["consensus_sequence_signature"] = [
        _consensus_sequence_signature(row, retro_class=retro_class)
        for (_, row), retro_class in zip(out.iterrows(), classes)
    ]
    return out


def _stable_tsv_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Return a TSV-export copy with stable column dtypes for pandas re-import."""
    out = df.copy()
    for col in out.columns:
        series = out[col]
        if pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series):
            out[col] = series.where(series.notna(), "").astype(str)
        elif pd.api.types.is_bool_dtype(series):
            out[col] = series.astype(int)
    return out


def _round_sig_value(value: float, sig: int) -> float:
    if pd.isna(value):
        return value
    if value == 0:
        return 0.0
    return round(float(value), int(sig - math.floor(math.log10(abs(float(value)))) - 1))


def _round_sig_series(series: pd.Series, sig: int) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    rounded = numeric.apply(lambda v: _round_sig_value(v, sig))
    return rounded.where(numeric.notna(), series)


def _estimate_discordant_ref_end(df: pd.DataFrame) -> pd.Series:
    """1-based inclusive mapped end; prefer stored ref_end, else pos+read_len-1."""
    pos = pd.to_numeric(df.get("pos", 0), errors="coerce").fillna(0).astype(int)
    if "ref_end" in df.columns:
        ref_end = pd.to_numeric(df["ref_end"], errors="coerce").fillna(0).astype(int)
    else:
        ref_end = pd.Series(0, index=df.index, dtype=int)
    if "read_seq" in df.columns:
        read_len = df["read_seq"].fillna("").astype(str).str.len().astype(int)
    else:
        read_len = pd.Series(0, index=df.index, dtype=int)
    estimated = (pos + read_len - 1).where(read_len > 0, pos)
    return ref_end.where(ref_end > 0, estimated).astype(int)


def _discordant_junction_tips(df: pd.DataFrame) -> pd.DataFrame:
    """Per-read junction-proximal genomic tip for discordant anchors."""
    if df is None or df.empty:
        return pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "anchor_side", "junction_tip", "soft_clipped"]
        )
    required = {"chrom", "window_start", "window_end", "pos", "read_name"}
    if not required.issubset(df.columns):
        return pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "anchor_side", "junction_tip", "soft_clipped"]
        )

    work = df.loc[
        :,
        [
            c
            for c in list(required)
            + ["ref_end", "soft_clip_pos", "soft_clip_side", "soft_clip_len", "read_seq", "mei_hit", "mate_mei_hit"]
            if c in df.columns
        ],
    ].copy()
    # Discordant MEI support is often on the mate (interchrom / large-insert), not the anchor.
    # Soft-clipped anchors are junction-resolving even before MEI remapping succeeds.
    mei_hit = work["mei_hit"].fillna(False).astype(bool) if "mei_hit" in work.columns else pd.Series(False, index=work.index)
    mate_mei = (
        work["mate_mei_hit"].fillna(False).astype(bool)
        if "mate_mei_hit" in work.columns
        else pd.Series(False, index=work.index)
    )
    soft_pos = (
        pd.to_numeric(work["soft_clip_pos"], errors="coerce").fillna(0).astype(int)
        if "soft_clip_pos" in work.columns
        else pd.Series(0, index=work.index, dtype=int)
    )
    if "mei_hit" in work.columns or "mate_mei_hit" in work.columns or "soft_clip_pos" in work.columns:
        work = work.loc[mei_hit | mate_mei | soft_pos.gt(0)].copy()
    if work.empty:
        return pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "anchor_side", "junction_tip", "soft_clipped"]
        )

    work["pos"] = pd.to_numeric(work["pos"], errors="coerce").fillna(0).astype(int)
    work = work.loc[work["pos"] > 0].copy()
    if work.empty:
        return pd.DataFrame(
            columns=["chrom", "window_start", "window_end", "read_name", "anchor_side", "junction_tip", "soft_clipped"]
        )

    work["ref_end"] = _estimate_discordant_ref_end(work)
    if "soft_clip_pos" in work.columns:
        work["soft_clip_pos"] = pd.to_numeric(work["soft_clip_pos"], errors="coerce").fillna(0).astype(int)
    else:
        work["soft_clip_pos"] = 0
    if "soft_clip_side" in work.columns:
        work["soft_clip_side"] = work["soft_clip_side"].fillna("").astype(str).str.upper().str[:1]
    else:
        work["soft_clip_side"] = ""

    work["anchor_side"] = ""
    work["junction_tip"] = 0
    work["soft_clipped"] = work["soft_clip_pos"].gt(0) & work["soft_clip_side"].isin(["L", "R"])

    # Soft-clip geometry: right clip → left flank; left clip → right flank.
    clipped = work["soft_clipped"]
    if clipped.any():
        work.loc[clipped, "anchor_side"] = work.loc[clipped, "soft_clip_side"].map({"R": "L", "L": "R"})
        work.loc[clipped, "junction_tip"] = work.loc[clipped, "soft_clip_pos"]

    # Assign unclipped sides per locus so multi-locus batches do not share one BP estimate.
    unclipped = ~work["soft_clipped"]
    key_cols = ["chrom", "window_start", "window_end"]
    for (_, ws, we), grp in work.loc[unclipped].groupby(key_cols, sort=False):
        if grp.empty:
            continue
        if clipped.any():
            clip_tips = work.loc[
                clipped
                & work["chrom"].eq(grp["chrom"].iloc[0])
                & work["window_start"].eq(ws)
                & work["window_end"].eq(we),
                "junction_tip",
            ]
            if not clip_tips.empty:
                bp_est = float(clip_tips.median())
            else:
                bp_est = None
        else:
            bp_est = None
        if bp_est is None:
            positions = sorted({int(p) for p in grp["pos"].tolist()})
            if len(positions) >= 2:
                gaps = [(positions[i + 1] - positions[i], i) for i in range(len(positions) - 1)]
                _gap, idx = max(gaps)
                bp_est = (positions[idx] + positions[idx + 1]) / 2.0
            else:
                bp_est = float((int(ws) + int(we)) / 2.0)
        left_idx = grp.index[grp["pos"].le(bp_est)]
        right_idx = grp.index[~grp["pos"].le(bp_est)]
        work.loc[left_idx, "anchor_side"] = "L"
        work.loc[right_idx, "anchor_side"] = "R"
        work.loc[left_idx, "junction_tip"] = work.loc[left_idx, "ref_end"]
        work.loc[right_idx, "junction_tip"] = work.loc[right_idx, "pos"]

    work["read_name"] = work["read_name"].fillna("").astype(str)
    work = work.loc[work["read_name"].str.len() > 0].copy()
    work = work.loc[work["anchor_side"].isin(["L", "R"]) & work["junction_tip"].gt(0)].copy()
    return work.loc[:, ["chrom", "window_start", "window_end", "read_name", "anchor_side", "junction_tip", "soft_clipped"]]


def _estimate_discordant_gap_intervals(
    discordant_frames: list[pd.DataFrame],
    *,
    min_reads_per_side: int = 2,
    max_gap_bp: int = 250,
    min_soft_clip_support: int = 1,
) -> pd.DataFrame:
    """Estimate insertion intervals from discordant L/R tips and soft-clip modes.

    Soft-clipped discordant anchors are junction-resolving even when support is
    one-sided: the clip coordinate is the insertion breakpoint.
    """
    key_cols = ["chrom", "window_start", "window_end"]
    empty = pd.DataFrame(
        columns=key_cols
        + [
            "dpe_gap_left",
            "dpe_gap_right",
            "dpe_gap_pos",
            "dpe_gap_n_left",
            "dpe_gap_n_right",
            "dpe_gap_n_soft_clip",
            "dpe_soft_clip_mode",
            "dpe_soft_clip_support",
            "dpe_side_window_start",
            "dpe_side_window_end",
            "dpe_gap_two_sided",
            "dpe_soft_clip_resolved",
        ]
    )
    parts = [_discordant_junction_tips(frame) for frame in discordant_frames if frame is not None and not frame.empty]
    if not parts:
        return empty
    tips = pd.concat(parts, ignore_index=True)
    if tips.empty:
        return empty

    rows: list[dict[str, object]] = []
    for (chrom, ws, we), grp in tips.groupby(key_cols, sort=False):
        left = grp.loc[grp["anchor_side"].eq("L"), "junction_tip"].astype(int)
        right = grp.loc[grp["anchor_side"].eq("R"), "junction_tip"].astype(int)
        n_l = int(grp.loc[grp["anchor_side"].eq("L"), "read_name"].nunique())
        n_r = int(grp.loc[grp["anchor_side"].eq("R"), "read_name"].nunique())
        clipped = grp.loc[grp["soft_clipped"]].copy()
        n_clip = int(clipped["read_name"].nunique()) if not clipped.empty else 0
        soft_mode = 0
        soft_support = 0
        if not clipped.empty:
            # Mode of soft-clip junction positions (unique reads per position).
            clip_counts = (
                clipped.groupby("junction_tip", as_index=False)["read_name"]
                .nunique()
                .sort_values(["read_name", "junction_tip"], ascending=[False, True])
            )
            if not clip_counts.empty:
                soft_mode = int(clip_counts.iloc[0]["junction_tip"])
                soft_support = int(clip_counts.iloc[0]["read_name"])
        soft_resolved = soft_mode > 0 and soft_support >= int(min_soft_clip_support)

        side_start = int(grp["junction_tip"].min())
        side_end = int(grp["junction_tip"].max())
        if soft_resolved:
            # Soft-clip mode is a point estimate; shrink side window to that base.
            side_start = soft_mode
            side_end = soft_mode

        two_sided = n_l >= int(min_reads_per_side) and n_r >= int(min_reads_per_side) and (not left.empty) and (not right.empty)
        dpe_left = 0
        dpe_right = 0
        dpe_pos = 0
        if soft_resolved:
            dpe_left = soft_mode
            dpe_right = soft_mode
            dpe_pos = soft_mode
        elif two_sided:
            # Junction-facing edges: rightmost left-anchor tip, leftmost right-anchor tip.
            left_edge = int(left.max())
            right_edge = int(right.min())
            lo = min(left_edge, right_edge)
            hi = max(left_edge, right_edge)
            if (hi - lo) <= int(max_gap_bp):
                dpe_left = lo
                dpe_right = hi
                dpe_pos = int((lo + hi) // 2)
                side_start = lo
                side_end = hi
            else:
                two_sided = False
        rows.append(
            {
                "chrom": str(chrom),
                "window_start": int(ws),
                "window_end": int(we),
                "dpe_gap_left": int(dpe_left),
                "dpe_gap_right": int(dpe_right),
                "dpe_gap_pos": int(dpe_pos),
                "dpe_gap_n_left": int(n_l),
                "dpe_gap_n_right": int(n_r),
                "dpe_gap_n_soft_clip": int(n_clip),
                "dpe_soft_clip_mode": int(soft_mode),
                "dpe_soft_clip_support": int(soft_support),
                "dpe_side_window_start": int(side_start),
                "dpe_side_window_end": int(side_end),
                "dpe_gap_two_sided": bool(two_sided and dpe_pos > 0 and not soft_resolved),
                "dpe_soft_clip_resolved": bool(soft_resolved),
            }
        )
    if not rows:
        return empty
    return pd.DataFrame(rows)


def _apply_discordant_gap_breakpoint_fallback(
    candidates: pd.DataFrame,
    *,
    discordant_disease: pd.DataFrame | None = None,
    discordant_control: pd.DataFrame | None = None,
    min_reads_per_side: int = 2,
    max_gap_bp: int = 250,
    min_soft_clip_support: int = 1,
) -> pd.DataFrame:
    """Fill unresolved breakpoints from discordant soft-clips or L/R gaps.

    Priority for unresolved loci:
      1. soft-clipped DPE mode (one- or two-sided) → discordant_soft_clip
      2. two-sided DPE gap midpoint → discordant_gap_midpoint
    Never overwrites SR/TSD/assembly-resolved breakpoints.
    """
    out = candidates.copy()
    key_cols = ["chrom", "window_start", "window_end"]
    for col, default in [
        ("dpe_gap_left", 0),
        ("dpe_gap_right", 0),
        ("dpe_gap_pos", 0),
        ("dpe_gap_n_left", 0),
        ("dpe_gap_n_right", 0),
        ("dpe_gap_n_soft_clip", 0),
        ("dpe_soft_clip_mode", 0),
        ("dpe_soft_clip_support", 0),
        ("dpe_side_window_start", 0),
        ("dpe_side_window_end", 0),
        ("dpe_gap_two_sided", False),
        ("dpe_soft_clip_resolved", False),
    ]:
        if col not in out.columns:
            out[col] = default

    gap = _estimate_discordant_gap_intervals(
        [
            discordant_disease if discordant_disease is not None else pd.DataFrame(),
            discordant_control if discordant_control is not None else pd.DataFrame(),
        ],
        min_reads_per_side=min_reads_per_side,
        max_gap_bp=max_gap_bp,
        min_soft_clip_support=min_soft_clip_support,
    )
    if gap.empty:
        return out

    drop_cols = [c for c in gap.columns if c not in key_cols and c in out.columns]
    if drop_cols:
        out = out.drop(columns=drop_cols)
    out = out.merge(gap, on=key_cols, how="left")
    for col, default in [
        ("dpe_gap_left", 0),
        ("dpe_gap_right", 0),
        ("dpe_gap_pos", 0),
        ("dpe_gap_n_left", 0),
        ("dpe_gap_n_right", 0),
        ("dpe_gap_n_soft_clip", 0),
        ("dpe_soft_clip_mode", 0),
        ("dpe_soft_clip_support", 0),
        ("dpe_side_window_start", 0),
        ("dpe_side_window_end", 0),
        ("dpe_gap_two_sided", False),
        ("dpe_soft_clip_resolved", False),
    ]:
        if col not in out.columns:
            out[col] = default
        if col in {"dpe_gap_two_sided", "dpe_soft_clip_resolved"}:
            out[col] = out[col].fillna(False).astype(bool)
        else:
            out[col] = pd.to_numeric(out[col], errors="coerce").fillna(default).astype(type(default))

    bp = pd.to_numeric(out.get("insertion_breakpoint_pos", 0), errors="coerce").fillna(0).astype(int)
    unresolved = bp.le(0)
    if "breakpoint_evidence_source" not in out.columns:
        out["breakpoint_evidence_source"] = ""

    use_clip = unresolved & out["dpe_soft_clip_resolved"] & out["dpe_soft_clip_mode"].gt(0)
    out.loc[use_clip, "insertion_breakpoint_pos"] = out.loc[use_clip, "dpe_soft_clip_mode"].astype(int)
    out.loc[use_clip, "breakpoint_evidence_source"] = "discordant_soft_clip"

    still_unresolved = pd.to_numeric(out["insertion_breakpoint_pos"], errors="coerce").fillna(0).le(0)
    use_gap = still_unresolved & out["dpe_gap_two_sided"] & out["dpe_gap_pos"].gt(0)
    out.loc[use_gap, "insertion_breakpoint_pos"] = out.loc[use_gap, "dpe_gap_pos"].astype(int)
    out.loc[use_gap, "breakpoint_evidence_source"] = "discordant_gap_midpoint"
    return out


def _tighten_windows_to_breakpoint_interval(
    df: pd.DataFrame,
    *,
    breakpoint_pos_col: str,
    interval_start_col: str,
    interval_end_col: str,
) -> pd.DataFrame:
    """Set window_start/end to the resolved breakpoint interval (no pad).

    - TSD / two-sided split / DPE gap: window = interval
    - Point-resolved BP: window = [pos, pos]
    - Unresolved one-sided DPE: keep pos unresolved, shrink window to that side's tips

    Preserves the original discovery span in discovery_window_{start,end} so
    supporting-read detail tables (keyed on discovery windows) still join.
    """
    out = df.copy()
    if out.empty or "window_start" not in out.columns or "window_end" not in out.columns:
        return out

    if "discovery_window_start" not in out.columns:
        out["discovery_window_start"] = pd.to_numeric(out["window_start"], errors="coerce").fillna(0).astype(int)
    if "discovery_window_end" not in out.columns:
        out["discovery_window_end"] = pd.to_numeric(out["window_end"], errors="coerce").fillna(0).astype(int)

    ws = pd.to_numeric(out["discovery_window_start"], errors="coerce")
    we = pd.to_numeric(out["discovery_window_end"], errors="coerce")
    bp = pd.to_numeric(out.get(breakpoint_pos_col, 0), errors="coerce")
    lo = pd.to_numeric(out.get(interval_start_col, float("nan")), errors="coerce")
    hi = pd.to_numeric(out.get(interval_end_col, float("nan")), errors="coerce")

    # Prefer explicit interval when valid.
    valid_interval = lo.notna() & hi.notna() & lo.gt(0) & hi.gt(0) & hi.ge(lo)
    new_ws = ws.copy()
    new_we = we.copy()
    new_ws = new_ws.where(~valid_interval, lo)
    new_we = new_we.where(~valid_interval, hi)

    # Point estimate with no interval → 1 bp window.
    point = (~valid_interval) & bp.notna() & bp.gt(0)
    new_ws = new_ws.where(~point, bp)
    new_we = new_we.where(~point, bp)

    # One-sided / unresolved DPE: shrink to observed junction-tip span on that side.
    side_lo = pd.to_numeric(out.get("dpe_side_window_start", float("nan")), errors="coerce")
    side_hi = pd.to_numeric(out.get("dpe_side_window_end", float("nan")), errors="coerce")
    unresolved = (~valid_interval) & (bp.isna() | bp.le(0))
    use_side = unresolved & side_lo.notna() & side_hi.notna() & side_lo.gt(0) & side_hi.ge(side_lo)
    new_ws = new_ws.where(~use_side, side_lo)
    new_we = new_we.where(~use_side, side_hi)

    # Keep windows ordered and within original discovery span when possible.
    ordered_ws = pd.concat([new_ws, new_we], axis=1).min(axis=1)
    ordered_we = pd.concat([new_ws, new_we], axis=1).max(axis=1)
    # Do not expand beyond the discovery window.
    ordered_ws = pd.concat([ordered_ws, ws], axis=1).max(axis=1)
    ordered_we = pd.concat([ordered_we, we], axis=1).min(axis=1)
    # If clipping inverted the interval, fall back to discovery window.
    bad = ordered_ws.isna() | ordered_we.isna() | ordered_we.lt(ordered_ws)
    ordered_ws = ordered_ws.where(~bad, ws)
    ordered_we = ordered_we.where(~bad, we)

    out["window_start"] = ordered_ws.fillna(ws).astype(int)
    out["window_end"] = ordered_we.fillna(we).astype(int)
    return out


def _derive_breakpoint_interval_fields(
    df: pd.DataFrame,
    *,
    breakpoint_pos_col: str,
    output_prefix: str,
) -> pd.DataFrame:
    out = df.copy()

    def s(col: str, default: object) -> pd.Series:
        if col in out.columns:
            return out[col]
        return pd.Series([default] * len(out), index=out.index)

    tsd_left = pd.to_numeric(s("tsd_left_breakpoint", float("nan")), errors="coerce")
    tsd_right = pd.to_numeric(s("tsd_right_breakpoint", float("nan")), errors="coerce")
    dpe_left = pd.to_numeric(s("dpe_gap_left", float("nan")), errors="coerce")
    dpe_right = pd.to_numeric(s("dpe_gap_right", float("nan")), errors="coerce")
    split_candidates = pd.concat(
        [
            pd.to_numeric(s("disease_L_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("disease_R_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_L_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_R_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("disease_L_split_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("disease_R_split_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_L_split_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_R_split_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("dpe_soft_clip_mode", float("nan")), errors="coerce"),
        ],
        axis=1,
    )
    split_candidates = split_candidates.where(split_candidates.gt(0))
    tsd_ok = tsd_left.gt(0) & tsd_right.gt(0)
    dpe_ok = dpe_left.gt(0) & dpe_right.gt(0)
    split_lo = split_candidates.min(axis=1, skipna=True)
    split_hi = split_candidates.max(axis=1, skipna=True)
    split_ok = split_lo.notna() & split_hi.notna()

    # Prefer TSD span, then split-read modes, then discordant gap. Do not mix
    # weaker evidence into a tighter higher-confidence interval.
    bp_lo = pd.Series(float("nan"), index=out.index, dtype=float)
    bp_hi = pd.Series(float("nan"), index=out.index, dtype=float)
    bp_lo = bp_lo.where(~tsd_ok, tsd_left)
    bp_hi = bp_hi.where(~tsd_ok, tsd_right)
    use_split = (~tsd_ok) & split_ok
    bp_lo = bp_lo.where(~use_split, split_lo)
    bp_hi = bp_hi.where(~use_split, split_hi)
    use_dpe = (~tsd_ok) & (~split_ok) & dpe_ok
    bp_lo = bp_lo.where(~use_dpe, dpe_left)
    bp_hi = bp_hi.where(~use_dpe, dpe_right)

    bp_pos = pd.to_numeric(s(breakpoint_pos_col, float("nan")), errors="coerce")
    use_pos_single = bp_lo.isna() & bp_hi.isna() & bp_pos.gt(0)
    bp_lo = bp_lo.where(~use_pos_single, bp_pos)
    bp_hi = bp_hi.where(~use_pos_single, bp_pos)

    ws = pd.to_numeric(s("window_start", float("nan")), errors="coerce")
    we = pd.to_numeric(s("window_end", float("nan")), errors="coerce")
    bp_lo = bp_lo.clip(lower=ws, upper=we)
    bp_hi = bp_hi.clip(lower=ws, upper=we)
    valid_interval = bp_lo.notna() & bp_hi.notna() & bp_hi.ge(bp_lo)
    bp_width = (bp_hi - bp_lo + 1.0).where(valid_interval, float("nan"))

    out[f"{output_prefix}breakpoint_interval_start"] = bp_lo.where(valid_interval, -1).fillna(-1).astype(int)
    out[f"{output_prefix}breakpoint_interval_end"] = bp_hi.where(valid_interval, -1).fillna(-1).astype(int)
    out[f"{output_prefix}breakpoint_interval_width_bp"] = bp_width.where(valid_interval, -1).fillna(-1).astype(int)
    out[f"{output_prefix}breakpoint_confidence_tier"] = "unresolved"
    out.loc[valid_interval & bp_width.le(3.0), f"{output_prefix}breakpoint_confidence_tier"] = "high"
    out.loc[
        valid_interval & bp_width.gt(3.0) & bp_width.le(15.0),
        f"{output_prefix}breakpoint_confidence_tier",
    ] = "medium"
    out.loc[valid_interval & bp_width.gt(15.0), f"{output_prefix}breakpoint_confidence_tier"] = "low"

    # Use interval midpoint as point estimate when interval exists.
    interval_midpoint = ((bp_lo + bp_hi) / 2.0).round()
    use_interval_midpoint = valid_interval & interval_midpoint.notna()
    out.loc[use_interval_midpoint, breakpoint_pos_col] = interval_midpoint.loc[use_interval_midpoint].astype(int)
    return out


def _build_gold_review_table(candidates: pd.DataFrame, empirical_stage: bool = False, fragment_to_full_map: dict[str, FragmentToFullMap] | None = None) -> pd.DataFrame:
    out = _ensure_candidate_schema_defaults(candidates)
    def _series_or_default(col: str, default: object) -> pd.Series:
        if col in out.columns:
            return out[col]
        return pd.Series([default] * len(out), index=out.index)

    out["tsd_or_polyA_supported"] = (
        _series_or_default("tsd_detected", False).fillna(False).astype(bool)
        | (_series_or_default("disease_poly_at_reads", 0).fillna(0).astype(float) >= 1)
        | (_series_or_default("control_poly_at_reads", 0).fillna(0).astype(float) >= 1)
        | (_series_or_default("poly_at_reads", 0).fillna(0).astype(float) >= 1)
    )
    out = _add_known_mei_polymorphism_consensus(out)
    out = _annotate_consensus_retrotransposition_fields(out)
    out["two_sided_support"] = _two_sided_support_mask(out)
    out["poly_at_supported"] = _poly_at_supported_mask(out)
    out["disease_vaf"] = _series_or_default("asm_disease_vaf", float("nan")).astype(float)
    out["control_vaf"] = _series_or_default("asm_control_vaf", float("nan")).astype(float)
    out["vaf_delta"] = _series_or_default("asm_vaf_delta", float("nan")).astype(float)
    out["assembly_status"] = _series_or_default("asm_status", "not_run").fillna("not_run").astype(str)
    out["assembly_confidence_score"] = _series_or_default("assembly_confidence_score", 0.0).fillna(0.0).astype(float)
    asm_source = _series_or_default("asm_breakpoint_source", "").fillna("").astype(str)
    asm_has_mei = asm_source.isin(["disease", "control"])
    bp_pos = _series_or_default("insertion_breakpoint_pos", 0).fillna(0).astype(int)
    out["insertion_breakpoint_pos"] = bp_pos.where(bp_pos > 0, -1)

    # Assembly-preferred consensus fields for breakpoint/TSD/polyA only.
    # MEI 5'/3' coords and orientation are set from SR/DPE below (not assembly).
    asm_bp = pd.to_numeric(_series_or_default("asm_consensus_breakpoint_pos", float("nan")), errors="coerce")
    use_asm_bp = asm_has_mei & asm_bp.notna() & asm_bp.gt(0)
    out["consensus_insertion_breakpoint_pos"] = asm_bp.where(use_asm_bp, out["insertion_breakpoint_pos"]).astype(int)
    out["consensus_breakpoint_source"] = ""
    out.loc[use_asm_bp, "consensus_breakpoint_source"] = asm_source.loc[use_asm_bp]
    out.loc[out["consensus_breakpoint_source"] == "", "consensus_breakpoint_source"] = _series_or_default(
        "breakpoint_evidence_source", ""
    ).fillna("").astype(str)
    out = _derive_breakpoint_interval_fields(
        out,
        breakpoint_pos_col="consensus_insertion_breakpoint_pos",
        output_prefix="consensus_",
    )
    empty_source = out["consensus_breakpoint_source"].fillna("").astype(str).str.len().eq(0)
    out.loc[
        (out["consensus_breakpoint_interval_width_bp"].astype(int) > 0) & (~asm_has_mei) & empty_source,
        "consensus_breakpoint_source",
    ] = "interval_midpoint"
    out = _tighten_windows_to_breakpoint_interval(
        out,
        breakpoint_pos_col="consensus_insertion_breakpoint_pos",
        interval_start_col="consensus_breakpoint_interval_start",
        interval_end_col="consensus_breakpoint_interval_end",
    )

    asm_tsd_seq = _series_or_default("asm_tsd_seq", "").fillna("").astype(str)
    asm_tsd_len = pd.to_numeric(_series_or_default("asm_tsd_len", float("nan")), errors="coerce")
    asm_tsd_detected = asm_has_mei & asm_tsd_len.notna() & asm_tsd_len.ge(4)
    out["consensus_tsd_seq"] = asm_tsd_seq.where(
        asm_tsd_detected & (asm_tsd_seq.str.len() > 0),
        _series_or_default("tsd_seq", "").fillna("").astype(str),
    )
    base_tsd_len = pd.to_numeric(_series_or_default("tsd_len_estimate", float("nan")), errors="coerce")
    out["consensus_tsd_len_estimate"] = asm_tsd_len.where(asm_tsd_detected, base_tsd_len)
    # Final consensus guard: never publish polyA/T tails as TSD.
    poly_consensus = _poly_at_artifact_tsd_mask(out["consensus_tsd_seq"])
    if bool(poly_consensus.any()):
        out.loc[poly_consensus, "consensus_tsd_seq"] = ""
        out.loc[poly_consensus, "consensus_tsd_len_estimate"] = 0
    # Keep TSD sequence/length internally consistent in review output. Some upstream
    # rows can carry sequence but a zero/missing length estimate.
    consensus_tsd_seq_len = out["consensus_tsd_seq"].fillna("").astype(str).str.len().astype(float)
    need_len_from_seq = (
        consensus_tsd_seq_len.gt(0)
        & (
            out["consensus_tsd_len_estimate"].isna()
            | pd.to_numeric(out["consensus_tsd_len_estimate"], errors="coerce").fillna(0).le(0)
        )
    )
    out.loc[need_len_from_seq, "consensus_tsd_len_estimate"] = consensus_tsd_seq_len.loc[need_len_from_seq]
    out["consensus_tsd_detected"] = out["consensus_tsd_len_estimate"].fillna(0).astype(float) >= 4.0
    out["junction_overlap_sequence"] = (
        _series_or_default("asm_junction_overlap_sequence", "").fillna("").astype(str)
    )
    # Backward-compatible fallback for prior cache artifacts that may only carry
    # the strict microhomology key.
    mh_fallback = _series_or_default("asm_microhomology_sequence", "").fillna("").astype(str)
    out.loc[out["junction_overlap_sequence"] == "", "junction_overlap_sequence"] = mh_fallback.loc[
        out["junction_overlap_sequence"] == ""
    ]

    asm_poly = pd.to_numeric(_series_or_default("asm_polyA_max_run", float("nan")), errors="coerce")
    base_poly = pd.to_numeric(out.get("poly_at_max_run", 0), errors="coerce")
    picked_poly = asm_poly.where(asm_has_mei & asm_poly.notna(), base_poly)
    # Observed polyA/T length is a lower bound (read-length limited), not a true max.
    out["consensus_poly_at_min_bp"] = (
        pd.concat([picked_poly, base_poly], axis=1).max(axis=1).fillna(0).round().astype(int)
    )
    out["consensus_poly_at_supported"] = out["consensus_poly_at_min_bp"].fillna(0).astype(float) >= 8.0

    # MEI 5'/3' coords and orientation come from split reads first, then DPE.
    # Local assembly is intentionally excluded: contigs rarely span the full MEI
    # and previously truncated/inverted consensus footprints (e.g. SVA stubs).
    asm_span = pd.to_numeric(_series_or_default("asm_insertion_length", float("nan")), errors="coerce")
    base_span = pd.to_numeric(_series_or_default("insertion_mei_span", float("nan")), errors="coerce")
    asm_mei_start = pd.to_numeric(_series_or_default("asm_insertion_mei_start", float("nan")), errors="coerce")
    asm_mei_end = pd.to_numeric(_series_or_default("asm_insertion_mei_end", float("nan")), errors="coerce")
    disease_start = pd.to_numeric(_series_or_default("disease_insertion_mei_start", float("nan")), errors="coerce")
    control_start = pd.to_numeric(_series_or_default("control_insertion_mei_start", float("nan")), errors="coerce")
    disease_end = pd.to_numeric(_series_or_default("disease_insertion_mei_end", float("nan")), errors="coerce")
    control_end = pd.to_numeric(_series_or_default("control_insertion_mei_end", float("nan")), errors="coerce")

    d_l_sr_start = pd.to_numeric(_series_or_default("disease_L_mei_start", float("nan")), errors="coerce")
    d_l_sr_end = pd.to_numeric(_series_or_default("disease_L_mei_end", float("nan")), errors="coerce")
    d_r_sr_start = pd.to_numeric(_series_or_default("disease_R_mei_start", float("nan")), errors="coerce")
    d_r_sr_end = pd.to_numeric(_series_or_default("disease_R_mei_end", float("nan")), errors="coerce")
    n_l_sr_start = pd.to_numeric(_series_or_default("control_L_mei_start", float("nan")), errors="coerce")
    n_l_sr_end = pd.to_numeric(_series_or_default("control_L_mei_end", float("nan")), errors="coerce")
    n_r_sr_start = pd.to_numeric(_series_or_default("control_R_mei_start", float("nan")), errors="coerce")
    n_r_sr_end = pd.to_numeric(_series_or_default("control_R_mei_end", float("nan")), errors="coerce")
    d_sr_bilateral = (
        d_l_sr_start.gt(0) & d_l_sr_end.ge(d_l_sr_start) & d_r_sr_start.gt(0) & d_r_sr_end.ge(d_r_sr_start)
    )
    n_sr_bilateral = (
        n_l_sr_start.gt(0) & n_l_sr_end.ge(n_l_sr_start) & n_r_sr_start.gt(0) & n_r_sr_end.ge(n_r_sr_start)
    )
    d_sr_lo = pd.concat([d_l_sr_start, d_l_sr_end, d_r_sr_start, d_r_sr_end], axis=1).min(axis=1, skipna=True)
    d_sr_hi = pd.concat([d_l_sr_start, d_l_sr_end, d_r_sr_start, d_r_sr_end], axis=1).max(axis=1, skipna=True)
    n_sr_lo = pd.concat([n_l_sr_start, n_l_sr_end, n_r_sr_start, n_r_sr_end], axis=1).min(axis=1, skipna=True)
    n_sr_hi = pd.concat([n_l_sr_start, n_l_sr_end, n_r_sr_start, n_r_sr_end], axis=1).max(axis=1, skipna=True)

    # Prefer supporting-read detail min/max (same footprint used by read-architecture
    # plots). Aggregated L/R start fields and DPE medians often under-call span.
    d_detail_lo = pd.to_numeric(_series_or_default("disease_detail_mei_start_min", float("nan")), errors="coerce")
    d_detail_hi = pd.to_numeric(_series_or_default("disease_detail_mei_end_max", float("nan")), errors="coerce")
    n_detail_lo = pd.to_numeric(_series_or_default("control_detail_mei_start_min", float("nan")), errors="coerce")
    n_detail_hi = pd.to_numeric(_series_or_default("control_detail_mei_end_max", float("nan")), errors="coerce")
    combined_detail_lo = pd.to_numeric(_series_or_default("detail_mei_start_min", float("nan")), errors="coerce")
    combined_detail_hi = pd.to_numeric(_series_or_default("detail_mei_end_max", float("nan")), errors="coerce")
    d_detail_ok = d_detail_lo.gt(0) & d_detail_hi.ge(d_detail_lo)
    n_detail_ok = n_detail_lo.gt(0) & n_detail_hi.ge(n_detail_lo)
    combined_detail_ok = combined_detail_lo.gt(0) & combined_detail_hi.ge(combined_detail_lo)
    d_detail_span = (d_detail_hi - d_detail_lo + 1.0).where(d_detail_ok, 0.0)
    n_detail_span = (n_detail_hi - n_detail_lo + 1.0).where(n_detail_ok, 0.0)
    d_detail_reads = pd.to_numeric(
        _series_or_default("disease_detail_mei_mapped_reads", 0), errors="coerce"
    ).fillna(0.0)
    n_detail_reads = pd.to_numeric(
        _series_or_default("control_detail_mei_mapped_reads", 0), errors="coerce"
    ).fillna(0.0)

    # Consensus target length: only keep footprint sources that map onto this
    # element (drops off-family DPE mates / medians that inflate Alu spans).
    target_length = pd.to_numeric(_series_or_default("asm_mei_target_length", float("nan")), errors="coerce")
    if "mei_target_length" in out.columns:
        target_length = target_length.where(
            target_length.gt(0),
            pd.to_numeric(out["mei_target_length"], errors="coerce"),
        )
    fallback_tlen = pd.concat(
        [
            pd.to_numeric(_series_or_default("disease_L_mei_target_len", float("nan")), errors="coerce"),
            pd.to_numeric(_series_or_default("disease_R_mei_target_len", float("nan")), errors="coerce"),
            pd.to_numeric(_series_or_default("control_L_mei_target_len", float("nan")), errors="coerce"),
            pd.to_numeric(_series_or_default("control_R_mei_target_len", float("nan")), errors="coerce"),
        ],
        axis=1,
    ).where(lambda x: x.gt(0)).max(axis=1, skipna=True)
    # Prefer asm; only fall back when asm is missing. Do not take max(asm, side)
    # because side lengths can be off-family (LINE1 3294 on an Alu call).
    target_length = target_length.where(target_length.gt(0), fallback_tlen)
    d_detail_ok = d_detail_ok & _on_target_extent_ok(d_detail_lo, d_detail_hi, target_length)
    n_detail_ok = n_detail_ok & _on_target_extent_ok(n_detail_lo, n_detail_hi, target_length)
    combined_detail_ok = combined_detail_ok & _on_target_extent_ok(
        combined_detail_lo, combined_detail_hi, target_length
    )
    d_sr_bilateral = d_sr_bilateral & _on_target_extent_ok(d_sr_lo, d_sr_hi, target_length)
    n_sr_bilateral = n_sr_bilateral & _on_target_extent_ok(n_sr_lo, n_sr_hi, target_length)

    d_left_t_early = pd.to_numeric(
        _series_or_default("disease_discordant_mei_left_target_pos_median", float("nan")), errors="coerce"
    )
    d_right_t_early = pd.to_numeric(
        _series_or_default("disease_discordant_mei_right_target_pos_median", float("nan")), errors="coerce"
    )
    n_left_t_early = pd.to_numeric(
        _series_or_default("control_discordant_mei_left_target_pos_median", float("nan")), errors="coerce"
    )
    n_right_t_early = pd.to_numeric(
        _series_or_default("control_discordant_mei_right_target_pos_median", float("nan")), errors="coerce"
    )
    d_left_start_min = pd.to_numeric(
        _series_or_default("disease_discordant_mei_left_target_start_min", float("nan")), errors="coerce"
    )
    d_right_start_min = pd.to_numeric(
        _series_or_default("disease_discordant_mei_right_target_start_min", float("nan")), errors="coerce"
    )
    d_left_end_max = pd.to_numeric(
        _series_or_default("disease_discordant_mei_left_target_end_max", float("nan")), errors="coerce"
    )
    d_right_end_max = pd.to_numeric(
        _series_or_default("disease_discordant_mei_right_target_end_max", float("nan")), errors="coerce"
    )
    n_left_start_min = pd.to_numeric(
        _series_or_default("control_discordant_mei_left_target_start_min", float("nan")), errors="coerce"
    )
    n_right_start_min = pd.to_numeric(
        _series_or_default("control_discordant_mei_right_target_start_min", float("nan")), errors="coerce"
    )
    n_left_end_max = pd.to_numeric(
        _series_or_default("control_discordant_mei_left_target_end_max", float("nan")), errors="coerce"
    )
    n_right_end_max = pd.to_numeric(
        _series_or_default("control_discordant_mei_right_target_end_max", float("nan")), errors="coerce"
    )
    d_left_reads_early = pd.to_numeric(
        _series_or_default("disease_discordant_mei_left_supported_reads", 0), errors="coerce"
    ).fillna(0.0)
    d_right_reads_early = pd.to_numeric(
        _series_or_default("disease_discordant_mei_right_supported_reads", 0), errors="coerce"
    ).fillna(0.0)
    n_left_reads_early = pd.to_numeric(
        _series_or_default("control_discordant_mei_left_supported_reads", 0), errors="coerce"
    ).fillna(0.0)
    n_right_reads_early = pd.to_numeric(
        _series_or_default("control_discordant_mei_right_supported_reads", 0), errors="coerce"
    ).fillna(0.0)
    d_dpe_extent_lo = pd.concat([d_left_start_min, d_right_start_min], axis=1).min(axis=1, skipna=True)
    d_dpe_extent_hi = pd.concat([d_left_end_max, d_right_end_max], axis=1).max(axis=1, skipna=True)
    n_dpe_extent_lo = pd.concat([n_left_start_min, n_right_start_min], axis=1).min(axis=1, skipna=True)
    n_dpe_extent_hi = pd.concat([n_left_end_max, n_right_end_max], axis=1).max(axis=1, skipna=True)
    d_dpe_bilateral = (
        d_left_reads_early.ge(1.0)
        & d_right_reads_early.ge(1.0)
        & d_left_t_early.gt(0)
        & d_right_t_early.gt(0)
    )
    n_dpe_bilateral = (
        n_left_reads_early.ge(1.0)
        & n_right_reads_early.ge(1.0)
        & n_left_t_early.gt(0)
        & n_right_t_early.gt(0)
    )
    d_dpe_extent_ok = (
        d_dpe_bilateral
        & d_dpe_extent_lo.gt(0)
        & d_dpe_extent_hi.ge(d_dpe_extent_lo)
        & _on_target_extent_ok(d_dpe_extent_lo, d_dpe_extent_hi, target_length)
    )
    n_dpe_extent_ok = (
        n_dpe_bilateral
        & n_dpe_extent_lo.gt(0)
        & n_dpe_extent_hi.ge(n_dpe_extent_lo)
        & _on_target_extent_ok(n_dpe_extent_lo, n_dpe_extent_hi, target_length)
    )
    # Median-based DPE footprint is only a last resort (orientation signal, not span).
    d_dpe_lo = pd.concat([d_left_t_early, d_right_t_early], axis=1).min(axis=1, skipna=True)
    d_dpe_hi = pd.concat([d_left_t_early, d_right_t_early], axis=1).max(axis=1, skipna=True)
    n_dpe_lo = pd.concat([n_left_t_early, n_right_t_early], axis=1).min(axis=1, skipna=True)
    n_dpe_hi = pd.concat([n_left_t_early, n_right_t_early], axis=1).max(axis=1, skipna=True)
    d_dpe_bilateral = d_dpe_bilateral & _on_target_extent_ok(d_dpe_lo, d_dpe_hi, target_length)
    n_dpe_bilateral = n_dpe_bilateral & _on_target_extent_ok(n_dpe_lo, n_dpe_hi, target_length)

    # Drop per-sample insertion coords that were previously copied from assembly.
    asm_polluted_disease = (
        asm_has_mei
        & (asm_source == "disease")
        & disease_start.gt(0)
        & asm_mei_start.gt(0)
        & disease_start.eq(asm_mei_start)
        & disease_end.eq(asm_mei_end)
    )
    asm_polluted_control = (
        asm_has_mei
        & (asm_source == "control")
        & control_start.gt(0)
        & asm_mei_start.gt(0)
        & control_start.eq(asm_mei_start)
        & control_end.eq(asm_mei_end)
    )
    disease_pair_lo = pd.concat([disease_start, disease_end], axis=1).min(axis=1, skipna=True)
    disease_pair_hi = pd.concat([disease_start, disease_end], axis=1).max(axis=1, skipna=True)
    control_pair_lo = pd.concat([control_start, control_end], axis=1).min(axis=1, skipna=True)
    control_pair_hi = pd.concat([control_start, control_end], axis=1).max(axis=1, skipna=True)
    disease_pair_valid = (
        disease_start.gt(0)
        & disease_end.gt(0)
        & ~asm_polluted_disease
        & _on_target_extent_ok(disease_pair_lo, disease_pair_hi, target_length)
    )
    control_pair_valid = (
        control_start.gt(0)
        & control_end.gt(0)
        & ~asm_polluted_control
        & _on_target_extent_ok(control_pair_lo, control_pair_hi, target_length)
    )

    raw_start = pd.Series([float("nan")] * len(out), index=out.index)
    raw_end = pd.Series([float("nan")] * len(out), index=out.index)

    # 0) Supporting-read detail footprint (min/max of all mapped SR+DPE MEI coords).
    # Prefer the sample with more mapped MEI-supporting reads; only use span as a
    # tie-breaker. Choosing the larger span alone let a single off-target control
    # DPE mate inflate Alu footprints (e.g. 197-1462).
    choose_n_detail = n_detail_ok & (
        ~d_detail_ok
        | n_detail_reads.gt(d_detail_reads)
        | (n_detail_reads.eq(d_detail_reads) & n_detail_span.gt(d_detail_span))
    )
    choose_d_detail = d_detail_ok & ~choose_n_detail
    raw_start = raw_start.where(~choose_d_detail, d_detail_lo)
    raw_end = raw_end.where(~choose_d_detail, d_detail_hi)
    raw_start = raw_start.where(~choose_n_detail, n_detail_lo)
    raw_end = raw_end.where(~choose_n_detail, n_detail_hi)
    unresolved = raw_start.isna() | raw_end.isna()
    choose_combined_detail = unresolved & combined_detail_ok
    raw_start = raw_start.where(~choose_combined_detail, combined_detail_lo)
    raw_end = raw_end.where(~choose_combined_detail, combined_detail_hi)

    # 1) Bilateral split-read footprint = full mapped extent (min start, max end).
    unresolved = raw_start.isna() | raw_end.isna()
    d_sr_support = (
        pd.to_numeric(_series_or_default("disease_L_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("disease_R_mei_supported_reads", 0), errors="coerce").fillna(0.0)
    )
    n_sr_support = (
        pd.to_numeric(_series_or_default("control_L_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("control_R_mei_supported_reads", 0), errors="coerce").fillna(0.0)
    )
    choose_n_sr = unresolved & n_sr_bilateral & (~d_sr_bilateral | n_sr_support.gt(d_sr_support))
    choose_d_sr = unresolved & d_sr_bilateral & ~choose_n_sr
    raw_start = raw_start.where(~choose_d_sr, d_sr_lo)
    raw_end = raw_end.where(~choose_d_sr, d_sr_hi)
    raw_start = raw_start.where(~choose_n_sr, n_sr_lo)
    raw_end = raw_end.where(~choose_n_sr, n_sr_hi)

    # 2) Bilateral DPE full mapped extent (not medians) when SR footprint missing.
    unresolved = raw_start.isna() | raw_end.isna()
    d_dpe_support = d_left_reads_early + d_right_reads_early
    n_dpe_support = n_left_reads_early + n_right_reads_early
    choose_n_dpe_ext = unresolved & n_dpe_extent_ok & (~d_dpe_extent_ok | n_dpe_support.gt(d_dpe_support))
    choose_d_dpe_ext = unresolved & d_dpe_extent_ok & ~choose_n_dpe_ext
    raw_start = raw_start.where(~choose_d_dpe_ext, d_dpe_extent_lo)
    raw_end = raw_end.where(~choose_d_dpe_ext, d_dpe_extent_hi)
    raw_start = raw_start.where(~choose_n_dpe_ext, n_dpe_extent_lo)
    raw_end = raw_end.where(~choose_n_dpe_ext, n_dpe_extent_hi)

    # 3) Previously fused per-sample insertion coords (SR/DPE), excluding asm copies.
    unresolved = raw_start.isna() | raw_end.isna()
    disease_pick = unresolved & disease_pair_valid
    raw_start = raw_start.where(~disease_pick, disease_start)
    raw_end = raw_end.where(~disease_pick, disease_end)
    unresolved = raw_start.isna() | raw_end.isna()
    control_pick = unresolved & control_pair_valid
    raw_start = raw_start.where(~control_pick, control_start)
    raw_end = raw_end.where(~control_pick, control_end)

    # 4) Last resort: DPE side medians (orientation signal only; underestimates span).
    unresolved = raw_start.isna() | raw_end.isna()
    choose_n_dpe = unresolved & n_dpe_bilateral & (~d_dpe_bilateral | n_dpe_support.gt(d_dpe_support))
    choose_d_dpe = unresolved & d_dpe_bilateral & ~choose_n_dpe
    raw_start = raw_start.where(~choose_d_dpe, d_dpe_lo)
    raw_end = raw_end.where(~choose_d_dpe, d_dpe_hi)
    raw_start = raw_start.where(~choose_n_dpe, n_dpe_lo)
    raw_end = raw_end.where(~choose_n_dpe, n_dpe_hi)

    # Orientation: SR/DPE consolidated call only (never assembly).
    # Footprint stays strictly min-max of mapped SR/DPE coords (no target_len expansion).
    consolidated_orient = out.apply(_choose_consolidated_insertion_orientation, axis=1)
    read_orient = _series_or_default("insertion_orientation", "").fillna("").astype(str)
    # If insertion_orientation was previously overwritten by assembly, prefer
    # disease/control insertion orientations via the consolidated chooser.
    asm_orient = _series_or_default("asm_insertion_orientation", "").fillna("").astype(str)
    read_orient_clean = read_orient.where(
        ~(asm_has_mei & asm_orient.isin(["+", "-"]) & read_orient.eq(asm_orient)),
        "",
    )
    out["consensus_insertion_orientation"] = consolidated_orient.where(
        consolidated_orient.isin(["+", "-"]),
        read_orient_clean.where(read_orient_clean.isin(["+", "-"]), ""),
    )
    raw_start_num = pd.to_numeric(raw_start, errors="coerce")
    raw_end_num = pd.to_numeric(raw_end, errors="coerce")
    valid_coords = raw_start_num.gt(0) & raw_end_num.gt(0)
    out["consensus_insertion_mei_3p_coord"] = raw_start_num.where(
        raw_start_num >= raw_end_num,
        raw_end_num,
    ).where(valid_coords, -1)
    out["consensus_insertion_mei_5p_coord"] = raw_start_num.where(
        raw_start_num <= raw_end_num,
        raw_end_num,
    ).where(valid_coords, -1)
    out["consensus_insertion_mei_start"] = out["consensus_insertion_mei_3p_coord"]
    out["consensus_insertion_mei_end"] = out["consensus_insertion_mei_5p_coord"]
    span_from_coords = (
        out["consensus_insertion_mei_3p_coord"].astype(float) - out["consensus_insertion_mei_5p_coord"].astype(float) + 1.0
    )
    out["consensus_insertion_mei_span"] = span_from_coords.where(
        out["consensus_insertion_mei_3p_coord"].astype(float).gt(0)
        & out["consensus_insertion_mei_5p_coord"].astype(float).gt(0),
        base_span.where(base_span.notna() & base_span.gt(0), float("nan")),
    )
    # Non-assembly fallback: impute MEI-axis coords only when discordant evidence
    # provides both-side MEI target mapping (no placeholder 1..span fallback).
    consensus_3p = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
    consensus_5p = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
    consensus_span = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce")
    bp_resolved = pd.to_numeric(out["consensus_insertion_breakpoint_pos"], errors="coerce").fillna(-1).gt(0)
    d_left_t = pd.to_numeric(_series_or_default("disease_discordant_mei_left_target_pos_median", float("nan")), errors="coerce")
    d_right_t = pd.to_numeric(_series_or_default("disease_discordant_mei_right_target_pos_median", float("nan")), errors="coerce")
    n_left_t = pd.to_numeric(_series_or_default("control_discordant_mei_left_target_pos_median", float("nan")), errors="coerce")
    n_right_t = pd.to_numeric(_series_or_default("control_discordant_mei_right_target_pos_median", float("nan")), errors="coerce")
    d_left_reads = pd.to_numeric(_series_or_default("disease_discordant_mei_left_supported_reads", 0), errors="coerce").fillna(0.0)
    d_right_reads = pd.to_numeric(_series_or_default("disease_discordant_mei_right_supported_reads", 0), errors="coerce").fillna(0.0)
    n_left_reads = pd.to_numeric(_series_or_default("control_discordant_mei_left_supported_reads", 0), errors="coerce").fillna(0.0)
    n_right_reads = pd.to_numeric(_series_or_default("control_discordant_mei_right_supported_reads", 0), errors="coerce").fillna(0.0)
    d_two_sided = _series_or_default("disease_discordant_mei_two_sided_support", (d_left_reads.ge(1.0) & d_right_reads.ge(1.0)))
    n_two_sided = _series_or_default("control_discordant_mei_two_sided_support", (n_left_reads.ge(1.0) & n_right_reads.ge(1.0)))
    d_two_sided = d_two_sided.fillna(False).astype(bool)
    n_two_sided = n_two_sided.fillna(False).astype(bool)
    d_pair_valid = d_two_sided & d_left_t.gt(0) & d_right_t.gt(0)
    n_pair_valid = n_two_sided & n_left_t.gt(0) & n_right_t.gt(0)

    # Mixed split + DPE fallback (opposite-side pairing):
    # allow one split-side MEI coordinate on one side paired with opposite-side
    # DPE target coordinate, with sample-level consistency checks.
    d_l_split = pd.to_numeric(_series_or_default("disease_L_mei_start", float("nan")), errors="coerce")
    d_r_split = pd.to_numeric(_series_or_default("disease_R_mei_start", float("nan")), errors="coerce")
    n_l_split = pd.to_numeric(_series_or_default("control_L_mei_start", float("nan")), errors="coerce")
    n_r_split = pd.to_numeric(_series_or_default("control_R_mei_start", float("nan")), errors="coerce")
    d_family_ok = pd.to_numeric(_series_or_default("disease_discordant_mei_family_purity", 0.0), errors="coerce").fillna(0.0).ge(0.60)
    n_family_ok = pd.to_numeric(_series_or_default("control_discordant_mei_family_purity", 0.0), errors="coerce").fillna(0.0).ge(0.60)
    d_geom_ok = _series_or_default("disease_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    n_geom_ok = _series_or_default("control_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    d_self_ok = _series_or_default("disease_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    n_self_ok = _series_or_default("control_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    d_mixed_ok = d_family_ok & d_geom_ok & d_self_ok
    n_mixed_ok = n_family_ok & n_geom_ok & n_self_ok
    d_mix_lsplit_rdisc = d_mixed_ok & d_l_split.gt(0) & d_right_t.gt(0)
    d_mix_rsplit_ldisc = d_mixed_ok & d_r_split.gt(0) & d_left_t.gt(0)
    n_mix_lsplit_rdisc = n_mixed_ok & n_l_split.gt(0) & n_right_t.gt(0)
    n_mix_rsplit_ldisc = n_mixed_ok & n_r_split.gt(0) & n_left_t.gt(0)

    d_total = d_left_reads + d_right_reads
    n_total = n_left_reads + n_right_reads
    choose_control = n_pair_valid & (~d_pair_valid | n_total.gt(d_total))
    choose_disease = d_pair_valid & ~choose_control

    pair_5p = pd.Series([float("nan")] * len(out), index=out.index)
    pair_3p = pd.Series([float("nan")] * len(out), index=out.index)
    pair_5p = pair_5p.where(~choose_disease, pd.concat([d_left_t, d_right_t], axis=1).min(axis=1, skipna=True))
    pair_3p = pair_3p.where(~choose_disease, pd.concat([d_left_t, d_right_t], axis=1).max(axis=1, skipna=True))
    pair_5p = pair_5p.where(~choose_control, pd.concat([n_left_t, n_right_t], axis=1).min(axis=1, skipna=True))
    pair_3p = pair_3p.where(~choose_control, pd.concat([n_left_t, n_right_t], axis=1).max(axis=1, skipna=True))

    # Fill unresolved rows with mixed split+DPE opposite-side pairs.
    unresolved_pair = pair_5p.isna() | pair_3p.isna()
    d_mix = unresolved_pair & (d_mix_lsplit_rdisc | d_mix_rsplit_ldisc)
    n_mix = unresolved_pair & (~d_mix) & (n_mix_lsplit_rdisc | n_mix_rsplit_ldisc)
    if d_mix.any():
        d_mix_5p = pd.concat(
            [
                d_l_split.where(d_mix_lsplit_rdisc),
                d_left_t.where(d_mix_rsplit_ldisc),
            ],
            axis=1,
        ).min(axis=1, skipna=True)
        d_mix_3p = pd.concat(
            [
                d_right_t.where(d_mix_lsplit_rdisc),
                d_r_split.where(d_mix_rsplit_ldisc),
            ],
            axis=1,
        ).max(axis=1, skipna=True)
        pair_5p = pair_5p.where(~d_mix, d_mix_5p)
        pair_3p = pair_3p.where(~d_mix, d_mix_3p)
    if n_mix.any():
        n_mix_5p = pd.concat(
            [
                n_l_split.where(n_mix_lsplit_rdisc),
                n_left_t.where(n_mix_rsplit_ldisc),
            ],
            axis=1,
        ).min(axis=1, skipna=True)
        n_mix_3p = pd.concat(
            [
                n_right_t.where(n_mix_lsplit_rdisc),
                n_r_split.where(n_mix_rsplit_ldisc),
            ],
            axis=1,
        ).max(axis=1, skipna=True)
        pair_5p = pair_5p.where(~n_mix, n_mix_5p)
        pair_3p = pair_3p.where(~n_mix, n_mix_3p)
    pair_span = (pair_3p - pair_5p + 1.0).where(pair_5p.notna() & pair_3p.notna(), float("nan"))

    missing_any = consensus_span.fillna(0).le(0) | consensus_3p.fillna(-1).le(0) | consensus_5p.fillna(-1).le(0)
    can_impute = bp_resolved & missing_any & pair_span.fillna(0).gt(0)
    if can_impute.any():
        out.loc[can_impute, "consensus_insertion_mei_5p_coord"] = pair_5p.loc[can_impute]
        out.loc[can_impute, "consensus_insertion_mei_3p_coord"] = pair_3p.loc[can_impute]
        out.loc[can_impute, "consensus_insertion_mei_span"] = pair_span.loc[can_impute]

    # Secondary fallback (regression guard):
    # if breakpoint is resolved but strict paired-coordinate imputation is still
    # unavailable, use best available span hint and complete coords from any
    # existing side; as last resort use 1..span to avoid unresolved 0/-1 fields.
    consensus_3p = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
    consensus_5p = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
    consensus_span = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce")
    still_missing = bp_resolved & (
        consensus_span.fillna(0).le(0) | consensus_3p.fillna(-1).le(0) | consensus_5p.fillna(-1).le(0)
    )
    if still_missing.any():
        target_len_hint = pd.to_numeric(
            _series_or_default("asm_mei_target_length", float("nan")), errors="coerce"
        )
        if target_len_hint.fillna(0).le(0).any():
            side_hint = pd.concat(
                [
                    pd.to_numeric(_series_or_default("disease_L_mei_target_len", float("nan")), errors="coerce"),
                    pd.to_numeric(_series_or_default("disease_R_mei_target_len", float("nan")), errors="coerce"),
                    pd.to_numeric(_series_or_default("control_L_mei_target_len", float("nan")), errors="coerce"),
                    pd.to_numeric(_series_or_default("control_R_mei_target_len", float("nan")), errors="coerce"),
                    pd.to_numeric(base_span, errors="coerce"),
                ],
                axis=1,
            ).where(lambda x: x.gt(0)).max(axis=1, skipna=True)
            target_len_hint = target_len_hint.where(target_len_hint.gt(0), side_hint)
        relaxed_span = pd.concat(
            [
                consensus_span.where(consensus_span.gt(0)),
                target_len_hint.where(target_len_hint.gt(0)),
            ],
            axis=1,
        ).max(axis=1, skipna=True)
        can_relax = still_missing & relaxed_span.fillna(0).gt(0)
        if can_relax.any():
            out.loc[can_relax, "consensus_insertion_mei_span"] = relaxed_span.loc[can_relax]
            span_now = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce")
            c3_now = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
            c5_now = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
            only_3p = can_relax & c3_now.fillna(-1).gt(0) & c5_now.fillna(-1).le(0)
            only_5p = can_relax & c5_now.fillna(-1).gt(0) & c3_now.fillna(-1).le(0)
            out.loc[only_3p, "consensus_insertion_mei_5p_coord"] = (
                c3_now.loc[only_3p] - span_now.loc[only_3p] + 1.0
            ).clip(lower=1.0)
            out.loc[only_5p, "consensus_insertion_mei_3p_coord"] = (
                c5_now.loc[only_5p] + span_now.loc[only_5p] - 1.0
            ).clip(lower=1.0)

        # If still unresolved, prefer explicit mapped MEI coordinates first
        # (split/disc medians), then polyA + target-length projection.
        c3_now = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
        c5_now = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
        unresolved_both = still_missing & c3_now.fillna(-1).le(0) & c5_now.fillna(-1).le(0)
        mapped_tbl = pd.concat(
            [
                d_l_split,
                d_r_split,
                n_l_split,
                n_r_split,
                d_left_t,
                d_right_t,
                n_left_t,
                n_right_t,
            ],
            axis=1,
        )
        mapped_pos = mapped_tbl.where(mapped_tbl.gt(0))
        mapped_min = mapped_pos.min(axis=1, skipna=True)
        mapped_max = mapped_pos.max(axis=1, skipna=True)
        mapped_count = mapped_pos.count(axis=1)
        can_map_impute = unresolved_both & bp_resolved & mapped_count.ge(2) & mapped_min.notna() & mapped_max.notna()
        if can_map_impute.any():
            out.loc[can_map_impute, "consensus_insertion_mei_5p_coord"] = mapped_min.loc[can_map_impute]
            out.loc[can_map_impute, "consensus_insertion_mei_3p_coord"] = mapped_max.loc[can_map_impute]
            out.loc[can_map_impute, "consensus_insertion_mei_span"] = (
                mapped_max.loc[can_map_impute] - mapped_min.loc[can_map_impute] + 1.0
            )

        c3_now = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
        c5_now = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
        unresolved_both = still_missing & c3_now.fillna(-1).le(0) & c5_now.fillna(-1).le(0)
        poly_reads_any = pd.concat(
            [
                pd.to_numeric(_series_or_default("disease_poly_at_reads", 0), errors="coerce").fillna(0.0),
                pd.to_numeric(_series_or_default("control_poly_at_reads", 0), errors="coerce").fillna(0.0),
                pd.to_numeric(_series_or_default("poly_at_reads", 0), errors="coerce").fillna(0.0),
            ],
            axis=1,
        ).max(axis=1)
        poly_run_any = pd.to_numeric(_series_or_default("consensus_poly_at_min_bp", 0), errors="coerce").fillna(0.0)
        poly_support = poly_reads_any.ge(1.0) | poly_run_any.ge(8.0)
        anchor_5p = mapped_min
        can_poly_project = (
            unresolved_both
            & bp_resolved
            & poly_support
            & target_len_hint.fillna(0).gt(0)
            & anchor_5p.notna()
            & anchor_5p.gt(0)
        )
        if can_poly_project.any():
            out.loc[can_poly_project, "consensus_insertion_mei_5p_coord"] = anchor_5p.loc[can_poly_project]
            out.loc[can_poly_project, "consensus_insertion_mei_3p_coord"] = pd.concat(
                [target_len_hint.loc[can_poly_project], anchor_5p.loc[can_poly_project]],
                axis=1,
            ).max(axis=1)
            out.loc[can_poly_project, "consensus_insertion_mei_span"] = (
                pd.to_numeric(out.loc[can_poly_project, "consensus_insertion_mei_3p_coord"], errors="coerce")
                - pd.to_numeric(out.loc[can_poly_project, "consensus_insertion_mei_5p_coord"], errors="coerce")
                + 1.0
            )

        # One-sided coordinate retention:
        # if polyA evidence indicates the tail side, keep a 3p coord even when
        # 5p stays unresolved; span will be clamped to 0 later.
        c3_now = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
        c5_now = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
        one_sided_poly_3p = (
            poly_support
            & mapped_max.notna()
            & mapped_max.gt(0)
            & c3_now.fillna(-1).le(0)
            & c5_now.fillna(-1).le(0)
        )
        if one_sided_poly_3p.any():
            out.loc[one_sided_poly_3p, "consensus_insertion_mei_3p_coord"] = mapped_max.loc[one_sided_poly_3p]

    c3_final = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce")
    c5_final = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce")
    coords_valid_final = c3_final.gt(0) & c5_final.gt(0)
    out["consensus_insertion_mei_span"] = (c3_final - c5_final + 1.0).where(coords_valid_final, 0.0)

    out["consensus_insertion_mei_3p_coord"] = pd.to_numeric(
        out["consensus_insertion_mei_3p_coord"],
        errors="coerce",
    ).fillna(-1).astype(int)
    out["consensus_insertion_mei_5p_coord"] = pd.to_numeric(
        out["consensus_insertion_mei_5p_coord"],
        errors="coerce",
    ).fillna(-1).astype(int)
    out["consensus_insertion_mei_start"] = out["consensus_insertion_mei_3p_coord"].astype(int)
    out["consensus_insertion_mei_end"] = out["consensus_insertion_mei_5p_coord"].astype(int)

    disease_full_start = pd.to_numeric(_series_or_default("disease_full_insertion_mei_start", float("nan")), errors="coerce")
    disease_full_end = pd.to_numeric(_series_or_default("disease_full_insertion_mei_end", float("nan")), errors="coerce")
    control_full_start = pd.to_numeric(_series_or_default("control_full_insertion_mei_start", float("nan")), errors="coerce")
    control_full_end = pd.to_numeric(_series_or_default("control_full_insertion_mei_end", float("nan")), errors="coerce")
    disease_full_valid = disease_full_start.gt(0) & disease_full_end.gt(0)
    control_full_valid = control_full_start.gt(0) & control_full_end.gt(0)
    disease_full_support = (
        pd.to_numeric(_series_or_default("disease_full_L_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("disease_full_R_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("disease_full_discordant_mei_supported_reads", 0), errors="coerce").fillna(0.0)
    )
    control_full_support = (
        pd.to_numeric(_series_or_default("control_full_L_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("control_full_R_mei_supported_reads", 0), errors="coerce").fillna(0.0)
        + pd.to_numeric(_series_or_default("control_full_discordant_mei_supported_reads", 0), errors="coerce").fillna(0.0)
    )
    choose_control_full = control_full_valid & (~disease_full_valid | control_full_support.gt(disease_full_support))
    choose_disease_full = disease_full_valid & ~choose_control_full
    raw_full_start = pd.Series([float("nan")] * len(out), index=out.index)
    raw_full_end = pd.Series([float("nan")] * len(out), index=out.index)
    raw_full_start = raw_full_start.where(~choose_disease_full, disease_full_start)
    raw_full_end = raw_full_end.where(~choose_disease_full, disease_full_end)
    raw_full_start = raw_full_start.where(~choose_control_full, control_full_start)
    raw_full_end = raw_full_end.where(~choose_control_full, control_full_end)
    full_3p = raw_full_start.where(raw_full_start >= raw_full_end, raw_full_end)
    full_5p = raw_full_start.where(raw_full_start <= raw_full_end, raw_full_end)
    full_valid = full_3p.gt(0) & full_5p.gt(0)
    out["consensus_insertion_mei_3p_coord_full"] = full_3p.where(full_valid, -1).fillna(-1).astype(int)
    out["consensus_insertion_mei_5p_coord_full"] = full_5p.where(full_valid, -1).fillna(-1).astype(int)
    out["consensus_insertion_mei_span_full"] = (
        pd.to_numeric(out["consensus_insertion_mei_3p_coord_full"], errors="coerce")
        - pd.to_numeric(out["consensus_insertion_mei_5p_coord_full"], errors="coerce")
        + 1.0
    ).where(full_valid, 0.0).fillna(0.0).astype(int)

    asm_subfamily = _series_or_default("asm_mei_subfamily", "").fillna("").astype(str)
    out["consensus_mei_subfamily"] = asm_subfamily.where(
        asm_subfamily.str.len() > 0,
        _series_or_default("mei_subfamily", "").fillna("").astype(str),
    )
    asm_family = _series_or_default("asm_mei_family", "").fillna("").astype(str)
    out["consensus_mei_family"] = asm_family.where(
        asm_family.str.len() > 0,
        _series_or_default("mei_family", "").fillna("").astype(str),
    )
    # COMPLEX_INS: blank MEI identity (nest/context stays in nested_* fields only).
    complex_ins = (
        _series_or_default("insertion_event_class", "").fillna("").astype(str).eq("COMPLEX_INS")
    )
    if complex_ins.any():
        out.loc[complex_ins, "consensus_mei_family"] = ""
        out.loc[complex_ins, "consensus_mei_subfamily"] = ""
    # Project each panel-fragment interval (e.g. L1HS_5end and L1HS_3end) onto
    # one full-length axis, then take the union. Never min/max raw panel coords
    # across fragment refs, and never mix Alu/L1 (or different full_name axes).
    frag_map = fragment_to_full_map or {}
    if frag_map:
        proj_5 = pd.Series(-1, index=out.index, dtype=int)
        proj_3 = pd.Series(-1, index=out.index, dtype=int)
        for idx in out.index:
            # Prefer already-valid per-read full remaps when present.
            cur_5 = int(pd.to_numeric(out.at[idx, "consensus_insertion_mei_5p_coord_full"], errors="coerce") or -1)
            cur_3 = int(pd.to_numeric(out.at[idx, "consensus_insertion_mei_3p_coord_full"], errors="coerce") or -1)
            union = _full_axis_union_from_panel_fragments(out.loc[idx], frag_map)
            if union is not None:
                u5, u3 = int(union[0]), int(union[1])
                if cur_5 > 0 and cur_3 >= cur_5:
                    # Expand existing full coords with the fragment-union footprint.
                    proj_5.at[idx] = min(cur_5, u5)
                    proj_3.at[idx] = max(cur_3, u3)
                else:
                    proj_5.at[idx] = u5
                    proj_3.at[idx] = u3
            elif cur_5 > 0 and cur_3 >= cur_5:
                proj_5.at[idx] = cur_5
                proj_3.at[idx] = cur_3
        out["consensus_insertion_mei_5p_coord_full"] = proj_5.astype(int)
        out["consensus_insertion_mei_3p_coord_full"] = proj_3.astype(int)
        full_ok = proj_5.gt(0) & proj_3.gt(0) & proj_3.ge(proj_5)
        out["consensus_insertion_mei_span_full"] = (proj_3 - proj_5 + 1).where(full_ok, 0).astype(int)
    # Keep *_full strictly on the full-length consensus axis. Do NOT copy panel /
    # partial-fragment coords (e.g. L1HS_3end 1–300) into *_full — that made
    # truncated L1s look like they started at genomic L1 position 1.
    full_5p_num = pd.to_numeric(out["consensus_insertion_mei_5p_coord_full"], errors="coerce").fillna(-1.0)
    full_3p_num = pd.to_numeric(out["consensus_insertion_mei_3p_coord_full"], errors="coerce").fillna(-1.0)
    full_pair_valid = full_5p_num.gt(0) & full_3p_num.gt(0)
    full_span_from_pair = (full_3p_num - full_5p_num + 1.0).where(full_pair_valid, 0.0)
    out["consensus_insertion_mei_span_full"] = pd.concat(
        [
            pd.to_numeric(out["consensus_insertion_mei_span_full"], errors="coerce").fillna(0.0),
            full_span_from_pair,
        ],
        axis=1,
    ).max(axis=1)

    # Prefer full-length coordinates for the primary consensus fields used by
    # gold tables and read-architecture plots.
    full_5p_num = pd.to_numeric(out["consensus_insertion_mei_5p_coord_full"], errors="coerce").fillna(-1.0)
    full_3p_num = pd.to_numeric(out["consensus_insertion_mei_3p_coord_full"], errors="coerce").fillna(-1.0)
    full_span_num = pd.to_numeric(out["consensus_insertion_mei_span_full"], errors="coerce").fillna(0.0)
    prefer_full = full_5p_num.gt(0) & full_3p_num.gt(0) & full_span_num.gt(0)
    if prefer_full.any():
        out.loc[prefer_full, "consensus_insertion_mei_5p_coord"] = full_5p_num.loc[prefer_full].astype(int)
        out.loc[prefer_full, "consensus_insertion_mei_3p_coord"] = full_3p_num.loc[prefer_full].astype(int)
        out.loc[prefer_full, "consensus_insertion_mei_span"] = full_span_num.loc[prefer_full].round().astype(int)
        out.loc[prefer_full, "consensus_insertion_mei_start"] = out.loc[
            prefer_full, "consensus_insertion_mei_3p_coord"
        ]
        out.loc[prefer_full, "consensus_insertion_mei_end"] = out.loc[
            prefer_full, "consensus_insertion_mei_5p_coord"
        ]

    # Enforce minimum reportable MEI span for both base and full coordinate
    # fields; shorter spans are treated as unresolved.
    base_span_num = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce").fillna(0.0)
    short_base = base_span_num.gt(0) & base_span_num.lt(_MIN_REPORTABLE_MEI_SPAN_BP)
    if short_base.any():
        out.loc[short_base, "consensus_insertion_mei_span"] = 0
        out.loc[short_base, "consensus_insertion_mei_5p_coord"] = -1
        out.loc[short_base, "consensus_insertion_mei_3p_coord"] = -1
    full_span_num = pd.to_numeric(out["consensus_insertion_mei_span_full"], errors="coerce").fillna(0.0)
    short_full = full_span_num.gt(0) & full_span_num.lt(_MIN_REPORTABLE_MEI_SPAN_BP)
    if short_full.any():
        out.loc[short_full, "consensus_insertion_mei_span_full"] = 0
        out.loc[short_full, "consensus_insertion_mei_5p_coord_full"] = -1
        out.loc[short_full, "consensus_insertion_mei_3p_coord_full"] = -1

    out["assembly_best_contig_id"] = _series_or_default("asm_consensus_primary_contig_id", "").fillna("").astype(str)
    out.loc[out["assembly_best_contig_id"] == "", "assembly_best_contig_id"] = _series_or_default(
        "asm_disease_primary_contig_id", ""
    ).fillna("").astype(str)
    out.loc[out["assembly_best_contig_id"] == "", "assembly_best_contig_id"] = _series_or_default(
        "asm_control_primary_contig_id", ""
    ).fillna("").astype(str)

    def _support_info_field(prefix: str) -> pd.Series:
        sr_l = pd.to_numeric(_series_or_default(f"{prefix}_L_mei_supported_reads", 0), errors="coerce").fillna(0).astype(int)
        sr_r = pd.to_numeric(_series_or_default(f"{prefix}_R_mei_supported_reads", 0), errors="coerce").fillna(0).astype(int)
        dpe_l = pd.to_numeric(
            _series_or_default(f"{prefix}_discordant_mei_left_supported_reads", 0), errors="coerce"
        ).fillna(0).astype(int)
        dpe_r = pd.to_numeric(
            _series_or_default(f"{prefix}_discordant_mei_right_supported_reads", 0), errors="coerce"
        ).fillna(0).astype(int)
        mei_mapped = sr_l + sr_r + dpe_l + dpe_r
        fam = _series_or_default("consensus_mei_family", "").fillna("").astype(str)
        if fam.eq("").all():
            fam = _series_or_default("mei_family", "").fillna("").astype(str)
        is_sva = fam.str.upper().eq("SVA")
        return pd.Series(
            [
                (
                    f"SR_L={sl},SR_R={sr},DPE_L={dl},DPE_R={dr},"
                    f"MEI_MAPPED={mm},polyA_MAPPED=0"
                    + (",VNTR_MAPPED=0" if sva else "")
                )
                for sl, sr, dl, dr, mm, sva in zip(
                    sr_l.tolist(),
                    sr_r.tolist(),
                    dpe_l.tolist(),
                    dpe_r.tolist(),
                    mei_mapped.tolist(),
                    is_sva.tolist(),
                )
            ],
            index=out.index,
        )

    existing_disease_support = _series_or_default("disease_supporting_reads", "").fillna("").astype(str)
    existing_control_support = _series_or_default("control_supporting_reads", "").fillna("").astype(str)
    out["disease_supporting_reads"] = existing_disease_support.where(
        existing_disease_support.str.len() > 0,
        _support_info_field("disease"),
    )
    out["control_supporting_reads"] = existing_control_support.where(
        existing_control_support.str.len() > 0,
        _support_info_field("control"),
    )
    def _with_mei_mapped(series: pd.Series, prefix: str) -> pd.Series:
        text = series.fillna("").astype(str)
        mapped_from_cols = pd.to_numeric(_series_or_default(f"{prefix}_mei_supported_reads", float("nan")), errors="coerce")
        mapped_from_token = pd.to_numeric(
            text.str.extract(r"MEI_MAPPED=([0-9]+)", expand=False),
            errors="coerce",
        )
        mapped = mapped_from_cols.where(mapped_from_cols.notna(), mapped_from_token).fillna(0).astype(int)
        polya = pd.to_numeric(
            text.str.extract(r"polyA_MAPPED=([0-9]+)", expand=False),
            errors="coerce",
        ).fillna(0).astype(int)
        vntr = pd.to_numeric(
            text.str.extract(r"VNTR_MAPPED=([0-9]+)", expand=False),
            errors="coerce",
        ).fillna(0).astype(int)
        fam = _series_or_default("consensus_mei_family", "").fillna("").astype(str)
        if fam.eq("").all():
            fam = _series_or_default("mei_family", "").fillna("").astype(str)
        is_sva = fam.str.upper().eq("SVA")
        base = text.str.replace(r",?MEI_MAPPED=[0-9]+", "", regex=True)
        base = base.str.replace(r",?polyA_MAPPED=[0-9]+", "", regex=True)
        base = base.str.replace(r",?VNTR_MAPPED=[0-9]+", "", regex=True)
        # polyA_side is rebuilt from the existing token when present.
        polya_side = text.str.extract(r"polyA_side=([LR])", expand=False).fillna("")
        base = base.str.replace(r",?polyA_side=[LR]", "", regex=True).str.strip(",")
        base = base.where(base.str.len() > 0, "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0")
        # Drop legacy BRK_CLP tokens if present in older support strings.
        base = base.str.replace(r",?BRK_CLP_[LR]=[0-9]+", "", regex=True).str.strip(",")
        out_s = base + ",MEI_MAPPED=" + mapped.astype(str) + ",polyA_MAPPED=" + polya.astype(str)
        out_s = out_s.where(~is_sva, out_s + ",VNTR_MAPPED=" + vntr.astype(str))
        out_s = out_s.where(polya_side.eq(""), out_s + ",polyA_side=" + polya_side)
        return out_s

    out["disease_supporting_reads"] = _with_mei_mapped(out["disease_supporting_reads"], "disease")
    out["control_supporting_reads"] = _with_mei_mapped(out["control_supporting_reads"], "control")
    out["nested_in_same_MEI"] = _series_or_default("nested_same_class_orientation", "").fillna("").astype(str)
    span_numeric = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce")
    out["consensus_insertion_mei_span"] = span_numeric.round().where(span_numeric.notna(), span_numeric)

    # Show compact breakpoint motif interpretation fields when motif signal is
    # significant by MT-adjusted YYRRRR log-odds.
    mt_adj = pd.to_numeric(_series_or_default("breakpoint_yyrrrr_logodds_shift1_mt_adj", float("nan")), errors="coerce")
    motif_reportable = mt_adj.notna() & (mt_adj > _YYRRRR_MT_ADJ_REPORT_MIN)
    # Keep observed breakpoint pattern visible for all resolved breakpoints.
    out["breakpoint_l1_en_observed_motif_pattern"] = (
        _series_or_default("breakpoint_l1_en_observed_motif_pattern", "").fillna("").astype(str)
    )
    # Gate derived motif interpretation fields to high-confidence motif calls.
    for col in (
        "breakpoint_l1_en_best_match_pattern_yy_rrrr",
        "breakpoint_l1_en_motif_type",
        "consensus_retrotransposition_class",
    ):
        out[col] = _series_or_default(col, "").fillna("").astype(str)
        out.loc[~motif_reportable, col] = ""

    empirical_cols = [
        "disease_empirical_local_bam_mean_depth_p_high",
        "disease_empirical_context_mapq_mean_p_low",
        "disease_empirical_context_mapq_lt20_fraction_p_high",
        "disease_empirical_context_nm_per_100bp_mean_p_high",
        "disease_empirical_context_nm_per_100bp_p90_p_high",
        "control_empirical_local_bam_mean_depth_p_high",
        "control_empirical_context_mapq_mean_p_low",
        "control_empirical_context_mapq_lt20_fraction_p_high",
        "control_empirical_context_nm_per_100bp_mean_p_high",
        "control_empirical_context_nm_per_100bp_p90_p_high",
        "gold_empirical_outlier",
    ]
    priority_cols = [
        "chrom",
        "consensus_insertion_breakpoint_pos",
        "window_start",
        "window_end",
        "discovery_window_start",
        "discovery_window_end",
        "control_supporting_reads",
        "disease_supporting_reads",
        "sample_status_label",
        "consensus_tsd_seq",
        "consensus_poly_at_min_bp",
        "consensus_mei_family",
        "consensus_mei_subfamily",
        "known_mei_polymorphism_family",
        "known_mei_polymorphism_id",
        "known_mei_polymorphism_source",
        "insertion_event_class",
        "complex_mei_event",
        "classic_polya_mei_sidepair",
        "complex_sv_signature_label",
        "discordant_mei_majority",
        "consensus_insertion_orientation",
        "nested_in_same_MEI",
        "consensus_insertion_mei_span_full",
        "consensus_insertion_mei_5p_coord_full",
        "consensus_insertion_mei_3p_coord_full",
    ]
    full_cols = [
        "analysis_stage_tier",
        "stage_fail_reason",
        "sample_status_label",
        "insertion_call_tier",
        "insertion_event_class",
        "complex_mei_event",
        "asm_complex_class",
        "chrom",
        "window_start",
        "window_end",
        "discovery_window_start",
        "discovery_window_end",
        "consensus_insertion_breakpoint_pos",
        "consensus_breakpoint_source",
        "consensus_breakpoint_interval_start",
        "consensus_breakpoint_interval_end",
        "consensus_breakpoint_interval_width_bp",
        "consensus_breakpoint_confidence_tier",
        "consensus_insertion_orientation",
        "consensus_insertion_mei_span_full",
        "consensus_insertion_mei_5p_coord_full",
        "consensus_insertion_mei_3p_coord_full",
        "asm_mei_target_length",
        "asm_insertion_length_observed",
        "asm_insertion_length_imputed",
        "asm_insertion_length_confidence_tier",
        "nested_in_same_MEI",
        "consensus_poly_at_min_bp",
        "mei_subfamily",
        "consensus_mei_family",
        "consensus_mei_subfamily",
        "consensus_tsd_seq",
        "consensus_tsd_len_estimate",
        "consensus_breakpoint_interval_width_bp",
        "junction_overlap_sequence",
        "known_mei_polymorphism",
        "known_mei_polymorphism_source",
        "known_mei_polymorphism_family",
        "known_mei_polymorphism_subfamily",
        "known_mei_polymorphism_id",
        "breakpoint_l1_en_observed_motif_pattern",
        "breakpoint_l1_en_best_match_pattern_yy_rrrr",
        "breakpoint_l1_en_motif_type",
        "consensus_retrotransposition_class",
        "disease_supporting_reads",
        "control_supporting_reads",
        "disease_supporting_reads_post_assembly",
        "control_supporting_reads_post_assembly",
        "disease_family_agreement",
        "disease_strand_agreement",
        "control_family_agreement",
        "control_strand_agreement",
        "two_sided_support",
        "assembly_best_contig_id",
        "asm_insertion_mei_start",
        "asm_insertion_mei_end",
        "asm_non_mei_partner_chrom",
        "asm_non_mei_partner_pos",
        "asm_non_mei_partner_type",
        "asm_breakpoint_side_status",
        "asm_complexity_source",
        "asm_top_contigs",
        "asm_mei_alignment_preset",
        "asm_left_support_contig_id",
        "asm_right_support_contig_id",
        "asm_left_support_mei_start",
        "asm_left_support_mei_end",
        "asm_right_support_mei_start",
        "asm_right_support_mei_end",
        "asm_left_support_mei_aln_len",
        "asm_right_support_mei_aln_len",
        "asm_coord_model",
        "disease_poly_at_reads",
        "disease_poly_at_max_run",
        "disease_poly_at_fraction_weighted",
        "control_poly_at_reads",
        "control_poly_at_max_run",
        "control_poly_at_fraction_weighted",
        "poly_at_reads",
        "poly_at_supported",
        "tsd_or_polyA_supported",
        "gold_stage_fail_reason",
        "insertion_model_score",
        "coherence_score",
        "mei_score_enrichment_ratio",
        "read_support_heuristic_score",
        "consensus_insertion_mei_span",
        "consensus_insertion_mei_5p_coord",
        "consensus_insertion_mei_3p_coord",
    ]
    selected_cols = list(priority_cols) + [c for c in full_cols if c not in set(priority_cols)]
    if empirical_stage:
        selected_cols = selected_cols[:-4] + empirical_cols + selected_cols[-4:]
    for col in selected_cols:
        if col not in out.columns:
            out[col] = ""
    review = out.loc[:, selected_cols].copy()
    if "consensus_insertion_mei_span" in review.columns:
        span = pd.to_numeric(review["consensus_insertion_mei_span"], errors="coerce")
        review["consensus_insertion_mei_span"] = span.round().fillna(-1).astype(int)
    if "consensus_insertion_mei_span_full" in review.columns:
        span_full = pd.to_numeric(review["consensus_insertion_mei_span_full"], errors="coerce")
        review["consensus_insertion_mei_span_full"] = span_full.round().fillna(-1).astype(int)

    sig4_cols = [
        "disease_vaf",
        "control_vaf",
        "vaf_delta",
    ]
    if empirical_stage:
        sig4_cols.extend(empirical_cols[:-1])
    sig3_cols = [
        "assembly_confidence_score",
        "disease_poly_at_fraction_weighted",
        "control_poly_at_fraction_weighted",
        "breakpoint_yyrrrr_logodds_shift1_mt_adj",
        "insertion_model_score",
        "coherence_score",
        "mei_score_enrichment_ratio",
        "read_support_heuristic_score",
    ]
    for col in sig4_cols:
        if col in review.columns:
            review[col] = _round_sig_series(review[col], sig=4)
    for col in sig3_cols:
        if col in review.columns:
            review[col] = _round_sig_series(review[col], sig=3)

    prioritized = _prioritize_mei_candidates(review, stage_first=True)
    return prioritized.loc[:, selected_cols]


def _infer_mei_family_from_fields(hit_id: str, family_hint: str, extra_hint: str) -> str:
    txt = " ".join([hit_id or "", family_hint or "", extra_hint or ""]).strip().upper()
    if not txt:
        return ""
    if "ALU" in txt:
        return "ALU"
    if "SVA" in txt:
        return "SVA"
    if "LINE1" in txt or "L1" in txt:
        return "LINE1"
    if "HERV" in txt or "ERV" in txt:
        return "ERV"
    return ""


def _extract_float_from_info(value: object, default: float = -1.0) -> float:
    if value is None:
        return default
    if isinstance(value, tuple):
        vals = [v for v in value if v is not None]
        if not vals:
            return default
        try:
            return float(max(vals))
        except Exception:
            return default
    try:
        return float(value)
    except Exception:
        return default


def _extract_int_from_info(value: object, default: int = -1) -> int:
    if value is None:
        return default
    if isinstance(value, tuple):
        vals = [v for v in value if v is not None]
        if not vals:
            return default
        try:
            return int(max(vals))
        except Exception:
            return default
    try:
        return int(value)
    except Exception:
        return default


def _is_mei_like_variant(vid: str, alt_txt: str, svtype: str, meinfo: str) -> bool:
    txt = " ".join([vid or "", alt_txt or "", svtype or "", meinfo or ""]).upper()
    markers = ("ALU", "SVA", "LINE", "L1", "MEI", "INS:ME")
    return any(m in txt for m in markers)


def _first_info_str(info: object) -> str:
    if info is None:
        return ""
    if isinstance(info, tuple):
        vals = [str(v) for v in info if v is not None]
        return vals[0] if vals else ""
    return str(info)


def _safe_info_get(info_map: object, key: str, default: object = None) -> object:
    # pysam raises ValueError for INFO keys absent from header definitions.
    if not hasattr(info_map, "get"):
        return default
    try:
        return info_map.get(key, default)
    except (KeyError, ValueError):
        return default


def _infer_subfamily_from_alt_meinfo(alt_txt: str, meinfo: str) -> str:
    # Prefer explicit MEINFO first token (e.g. "SVA,48,1315,-"), then ALT tags.
    me_first = (meinfo.split(",")[0] if meinfo else "").strip()
    if me_first:
        return me_first
    alt = (alt_txt or "").upper()
    for token in alt.replace("<", "").replace(">", "").split(":"):
        t = token.strip()
        if t and t not in {"INS", "ME"}:
            return t
    return ""


def _infer_tsd_from_info(info_map: object) -> str:
    # MELT/dbVar exports vary; TSD may appear as dedicated INFO or embedded in DESC.
    if hasattr(info_map, "get"):
        tsd = _first_info_str(_safe_info_get(info_map, "TSD", ""))
        if tsd:
            return tsd
        desc = _first_info_str(_safe_info_get(info_map, "DESC", ""))
        m = re.search(r"TSD(?:=|%3D)([A-Za-z]+)", desc)
        if m:
            return m.group(1).upper()
    return ""


def _build_g1k_mei_bed_from_vcf(vcf_path: Path, out_bed_path: Path) -> int:
    kept = 0
    prev_verbosity = pysam.set_verbosity(0)
    try:
        with pysam.VariantFile(str(vcf_path)) as vf, out_bed_path.open("w", encoding="utf-8") as oh:
            for rec in vf:
                chrom = str(rec.contig or "")
                if not chrom:
                    continue
                pos1 = int(rec.pos)
                ref = str(rec.ref or "")
                rid = str(rec.id or "")
                alts = [str(a) for a in (rec.alts or ())]
                alt_txt = ",".join(alts)
                info = rec.info
                svtype = _first_info_str(_safe_info_get(info, "SVTYPE", "")).upper()
                meinfo = _first_info_str(_safe_info_get(info, "MEINFO", ""))
                # Restrict to insertion MEI records only.
                is_insertion = svtype == "INS" or "INS:ME" in alt_txt.upper()
                if not is_insertion:
                    continue
                if not _is_mei_like_variant(vid=rid, alt_txt=alt_txt, svtype=svtype, meinfo=meinfo):
                    continue

                end1 = _extract_int_from_info(_safe_info_get(info, "END", None), default=pos1 + max(1, len(ref)) - 1)
                end1 = max(end1, pos1)
                start0 = max(0, pos1 - 1)
                end0 = max(start0 + 1, end1)

                melt_ins_type = "INS"
                melt_ins_subfamily = _infer_subfamily_from_alt_meinfo(alt_txt=alt_txt, meinfo=meinfo)
                melt_ins_len = abs(_extract_int_from_info(_safe_info_get(info, "SVLEN", None), default=-1))
                melt_tsd = _infer_tsd_from_info(info_map=info)
                melt_region_id = _first_info_str(_safe_info_get(info, "REGIONID", ""))
                rec_id = rid if rid and rid != "." else f"{chrom}:{pos1}:INS"
                # Keep a strict minimal schema to avoid downstream column drift.
                oh.write(
                    f"{chrom}\t{start0}\t{end0}\t{rec_id}\t{melt_ins_type}\t"
                    f"{melt_ins_subfamily}\t{melt_ins_len}\t{melt_tsd}\t{melt_region_id}\n"
                )
                kept += 1
    finally:
        pysam.set_verbosity(prev_verbosity)
    return kept


def _run_bedtools_checked(
    cmd: list[str],
    *,
    label: str,
) -> subprocess.CompletedProcess[str]:
    """Run bedtools and raise with stderr/stdout if it fails."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if int(proc.returncode) == 0:
        return proc
    err = (proc.stderr or "").strip()
    out = (proc.stdout or "").strip()
    click.echo(
        f"[mei-annotate] {label} failed exit={proc.returncode} cmd={' '.join(cmd)}"
    )
    if err:
        click.echo(f"[mei-annotate] {label} stderr:\n{err}")
    if out:
        click.echo(f"[mei-annotate] {label} stdout (truncated):\n{out[:2000]}")
    raise RuntimeError(
        f"{label} failed (exit={proc.returncode}): {' '.join(cmd)}"
        + (f"\nstderr:\n{err}" if err else "")
        + (f"\nstdout (truncated):\n{out[:2000]}" if out else "")
    )


def _normalize_bed_chrom_style(input_bed: Path, output_bed: Path, target_has_chr_prefix: bool) -> None:
    with input_bed.open("r", encoding="utf-8") as ih, output_bed.open("w", encoding="utf-8") as oh:
        for line in ih:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if target_has_chr_prefix:
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
            else:
                if chrom.startswith("chr"):
                    chrom = chrom[3:]
            parts[0] = chrom
            oh.write("\t".join(parts) + "\n")


def _g1k_query_interval_for_row(
    row: pd.Series,
    split_padding_bp: int,
    dpe_padding_min_bp: int,
    dpe_padding_max_bp: int,
    dpe_padding_tlen_factor: float,
) -> tuple[int, int]:
    window_start = int(row.get("window_start", 1))
    window_end = int(row.get("window_end", window_start))
    midpoint = (window_start + window_end) // 2
    breakpoint_pos = int(row.get("insertion_breakpoint_pos", 0))
    if breakpoint_pos <= 0:
        breakpoint_pos = midpoint

    left_split = int(row.get("disease_L_mei_supported_reads", 0))
    right_split = int(row.get("disease_R_mei_supported_reads", 0))
    split_total = int(row.get("disease_split_mei_supported_reads", 0))
    split_resolved = (split_total >= 2) or ((left_split >= 1) and (right_split >= 1))

    disease_dpe = int(row.get("disease_discordant_mei_supported_reads", 0))
    control_dpe = int(row.get("control_discordant_mei_supported_reads", 0))
    dpe_present = (disease_dpe + control_dpe) > 0
    dpe_tlen_mean = max(
        float(row.get("discordant_disease_abs_tlen_mean", 0.0) or 0.0),
        float(row.get("discordant_control_abs_tlen_mean", 0.0) or 0.0),
    )

    if split_resolved:
        pad = max(1, int(split_padding_bp))
        center = int(breakpoint_pos)
    elif dpe_present:
        dynamic_pad = max(int(dpe_padding_min_bp), int(round(dpe_tlen_mean * float(dpe_padding_tlen_factor))))
        pad = max(1, min(int(dpe_padding_max_bp), dynamic_pad))
        center = int(breakpoint_pos if breakpoint_pos > 0 else midpoint)
    else:
        pad = max(1, int(split_padding_bp) * 2)
        center = int(midpoint)

    start_1based = max(1, center - pad)
    end_1based = max(start_1based, center + pad)
    return start_1based, end_1based


def _annotate_g1k_mei_overlap(
    candidates: pd.DataFrame,
    g1k_mei_vcf: Path | None,
    split_padding_bp: int,
    dpe_padding_min_bp: int,
    dpe_padding_max_bp: int,
    dpe_padding_tlen_factor: float,
) -> pd.DataFrame:
    if g1k_mei_vcf is None:
        return candidates.copy()

    out = candidates.copy().reset_index(drop=True)
    out["g1k_melt_id"] = ""
    out["g1k_melt_insertion_type"] = ""
    out["g1k_melt_insertion_subfamily"] = ""
    out["g1k_melt_insertion_length"] = -1
    out["g1k_melt_tsd"] = ""
    out["g1k_melt_region_id"] = ""
    if out.empty:
        return out

    out["row_id"] = out.index.astype(int)
    row_by_id = {int(row.row_id): pd.Series(row._asdict()) for row in out.itertuples(index=False)}
    best_hits: dict[int, dict[str, object]] = {}

    with tempfile.TemporaryDirectory(prefix="rtm_g1k_mei_") as tmpdir:
        tmp = Path(tmpdir)
        source_bed = tmp / "g1k_mei_from_vcf.bed"
        kept = _build_g1k_mei_bed_from_vcf(g1k_mei_vcf, source_bed)
        click.echo(f"[mei-annotate] parsed g1k MEI VCF records kept={kept} path={g1k_mei_vcf}")

        query_bed = tmp / "candidate_g1k_query.bed"
        with query_bed.open("w", encoding="utf-8") as handle:
            for row in out.itertuples(index=False):
                start_1based, end_1based = _g1k_query_interval_for_row(
                    pd.Series(row._asdict()),
                    split_padding_bp=split_padding_bp,
                    dpe_padding_min_bp=dpe_padding_min_bp,
                    dpe_padding_max_bp=dpe_padding_max_bp,
                    dpe_padding_tlen_factor=dpe_padding_tlen_factor,
                )
                start0 = max(0, int(start_1based) - 1)
                end0 = max(start0 + 1, int(end_1based))
                handle.write(f"{row.chrom}\t{start0}\t{end0}\t{row.row_id}\n")

        query_has_chr_prefix = out["chrom"].astype(str).str.startswith("chr").any()
        source_bed_norm = tmp / "g1k_mei.chromnorm.bed"
        _normalize_bed_chrom_style(
            input_bed=source_bed,
            output_bed=source_bed_norm,
            target_has_chr_prefix=bool(query_has_chr_prefix),
        )

        intersect_cmd = ["bedtools", "intersect", "-a", str(query_bed), "-b", str(source_bed_norm), "-wa", "-wb"]
        proc = _run_bedtools_checked(intersect_cmd, label="g1k-mei bedtools intersect")
        for line in proc.stdout.splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            try:
                row_id = int(parts[3])
            except ValueError:
                continue
            b_cols = parts[4:]
            hit_id = b_cols[3] if len(b_cols) >= 4 and b_cols[3] not in {"", "."} else ""
            if not hit_id and len(b_cols) >= 3:
                hit_id = f"{b_cols[0]}:{b_cols[1]}-{b_cols[2]}"
            ins_type_s = b_cols[4] if len(b_cols) >= 5 else ""
            subfamily_s = b_cols[5] if len(b_cols) >= 6 else ""
            ins_len_i = -1
            if len(b_cols) >= 7:
                try:
                    ins_len_i = int(float(b_cols[6]))
                except ValueError:
                    ins_len_i = -1
            tsd_s = b_cols[7] if len(b_cols) >= 8 else ""
            region_s = b_cols[8] if len(b_cols) >= 9 else ""
            try:
                a_start0 = int(parts[1])
                a_end0 = int(parts[2])
                b_start0 = int(b_cols[1])
                b_end0 = int(b_cols[2])
                overlap_bp = max(0, min(a_end0, b_end0) - max(a_start0, b_start0))
            except (ValueError, IndexError):
                overlap_bp = 0
            row = row_by_id.get(row_id)
            if row is None:
                continue
            event_family = _choose_event_family(row)
            g1k_family = _normalize_mei_family_token(f"{hit_id} {ins_type_s} {subfamily_s}")
            if not event_family or g1k_family != event_family:
                continue
            current = best_hits.get(row_id)
            if (current is None) or (int(current.get("overlap_bp", -1)) < overlap_bp):
                best_hits[row_id] = {
                    "overlap_bp": overlap_bp,
                    "id": hit_id,
                    "ins_type": ins_type_s,
                    "subfamily": subfamily_s,
                    "ins_len": ins_len_i,
                    "tsd": tsd_s,
                    "region_id": region_s,
                }

    if best_hits:
        out["g1k_melt_id"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("id", "")))
        out["g1k_melt_insertion_type"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("ins_type", "")))
        out["g1k_melt_insertion_subfamily"] = out["row_id"].map(
            lambda i: str(best_hits.get(i, {}).get("subfamily", ""))
        )
        out["g1k_melt_insertion_length"] = (
            out["row_id"].map(lambda i: int(best_hits.get(i, {}).get("ins_len", -1))).fillna(-1).astype(int)
        )
        out["g1k_melt_tsd"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("tsd", "")))
        out["g1k_melt_region_id"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("region_id", "")))
    return out.drop(columns=["row_id"])


def _build_lr_mei_bed_from_vcf(vcf_path: Path, out_bed_path: Path) -> int:
    kept = 0
    prev_verbosity = pysam.set_verbosity(0)
    try:
        with pysam.VariantFile(str(vcf_path)) as vf, out_bed_path.open("w", encoding="utf-8") as oh:
            for rec in vf:
                chrom = str(rec.contig or "")
                if not chrom:
                    continue
                info = rec.info
                rid = str(rec.id or "")
                fam_n = _first_info_str(_safe_info_get(info, "FAM_N", "")).strip()
                if not fam_n:
                    continue
                norm_family = _normalize_mei_family_token(fam_n)
                if norm_family not in {"ALU", "SVA", "LINE1"}:
                    continue

                # Strict SVAN mode:
                # keep only insertion-like records by ID signature (e.g. chrXX-YYYYYY-INS)
                # and ignore DTYPE_N fallback (which can include deletion-like classes).
                itype_n = _first_info_str(_safe_info_get(info, "ITYPE_N", "")).strip()
                itype_u = itype_n.upper()
                rid_u = rid.upper()
                is_insertion_like = ("-INS" in rid_u) or rid_u.endswith("INS") or (":INS" in rid_u)
                if not is_insertion_like:
                    continue

                pos1 = int(rec.pos)
                ref = str(rec.ref or "")
                end1 = _extract_int_from_info(_safe_info_get(info, "END", None), default=pos1 + max(1, len(ref)) - 1)
                end1 = max(end1, pos1)
                start0 = max(0, pos1 - 1)
                end0 = max(start0 + 1, end1)

                rec_id = rid if rid and rid != "." else f"{chrom}:{pos1}:SVAN_MEI"
                event_type = itype_n if itype_n else "INS"
                subfamily = fam_n
                ins_len = abs(_extract_int_from_info(_safe_info_get(info, "INS_LEN", None), default=-1))
                tsd_len = abs(_extract_int_from_info(_safe_info_get(info, "TSD_LEN", None), default=-1))
                polya_len = abs(_extract_int_from_info(_safe_info_get(info, "POLYA_LEN", None), default=-1))
                conformation = _first_info_str(_safe_info_get(info, "CONFORMATION", "")).strip()
                not_canonical = 1 if bool(_safe_info_get(info, "NOT_CANONICAL", False)) else 0
                oh.write(
                    f"{chrom}\t{start0}\t{end0}\t{rec_id}\t{event_type}\t{subfamily}\t{ins_len}\t"
                    f"{tsd_len}\t{polya_len}\t{conformation}\t{not_canonical}\n"
                )
                kept += 1
    finally:
        pysam.set_verbosity(prev_verbosity)
    return kept


def _annotate_lr_mei_overlap(
    candidates: pd.DataFrame,
    lr_mei_vcf: Path | None,
    split_padding_bp: int,
    dpe_padding_min_bp: int,
    dpe_padding_max_bp: int,
    dpe_padding_tlen_factor: float,
) -> pd.DataFrame:
    if lr_mei_vcf is None:
        return candidates.copy()

    out = candidates.copy().reset_index(drop=True)
    out["lr_svan_id"] = ""
    out["lr_svan_event_type"] = ""
    out["lr_svan_subfamily"] = ""
    out["lr_svan_insertion_length"] = -1
    out["lr_svan_tsd_len"] = -1
    out["lr_svan_polya_len"] = -1
    out["lr_svan_conformation"] = ""
    out["lr_svan_not_canonical"] = False
    if out.empty:
        return out

    out["row_id"] = out.index.astype(int)
    row_by_id = {int(row.row_id): pd.Series(row._asdict()) for row in out.itertuples(index=False)}
    best_hits: dict[int, dict[str, object]] = {}

    with tempfile.TemporaryDirectory(prefix="rtm_lr_mei_") as tmpdir:
        tmp = Path(tmpdir)
        source_bed = tmp / "lr_mei_from_vcf.bed"
        kept = _build_lr_mei_bed_from_vcf(lr_mei_vcf, source_bed)
        click.echo(f"[mei-annotate] parsed long-read SVAN MEI VCF records kept={kept} path={lr_mei_vcf}")

        query_bed = tmp / "candidate_lr_query.bed"
        with query_bed.open("w", encoding="utf-8") as handle:
            for row in out.itertuples(index=False):
                start_1based, end_1based = _g1k_query_interval_for_row(
                    pd.Series(row._asdict()),
                    split_padding_bp=split_padding_bp,
                    dpe_padding_min_bp=dpe_padding_min_bp,
                    dpe_padding_max_bp=dpe_padding_max_bp,
                    dpe_padding_tlen_factor=dpe_padding_tlen_factor,
                )
                start0 = max(0, int(start_1based) - 1)
                end0 = max(start0 + 1, int(end_1based))
                handle.write(f"{row.chrom}\t{start0}\t{end0}\t{row.row_id}\n")

        query_has_chr_prefix = out["chrom"].astype(str).str.startswith("chr").any()
        source_bed_norm = tmp / "lr_mei.chromnorm.bed"
        _normalize_bed_chrom_style(
            input_bed=source_bed,
            output_bed=source_bed_norm,
            target_has_chr_prefix=bool(query_has_chr_prefix),
        )

        intersect_cmd = ["bedtools", "intersect", "-a", str(query_bed), "-b", str(source_bed_norm), "-wa", "-wb"]
        proc = _run_bedtools_checked(intersect_cmd, label="lr-mei bedtools intersect")
        for line in proc.stdout.splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            try:
                row_id = int(parts[3])
            except ValueError:
                continue
            b_cols = parts[4:]
            hit_id = b_cols[3] if len(b_cols) >= 4 and b_cols[3] not in {"", "."} else ""
            if not hit_id and len(b_cols) >= 3:
                hit_id = f"{b_cols[0]}:{b_cols[1]}-{b_cols[2]}"
            event_type = b_cols[4] if len(b_cols) >= 5 else ""
            subfamily = b_cols[5] if len(b_cols) >= 6 else ""
            ins_len = -1
            tsd_len = -1
            polya_len = -1
            not_canonical = False
            if len(b_cols) >= 7:
                try:
                    ins_len = int(float(b_cols[6]))
                except ValueError:
                    ins_len = -1
            if len(b_cols) >= 8:
                try:
                    tsd_len = int(float(b_cols[7]))
                except ValueError:
                    tsd_len = -1
            if len(b_cols) >= 9:
                try:
                    polya_len = int(float(b_cols[8]))
                except ValueError:
                    polya_len = -1
            conformation = b_cols[9] if len(b_cols) >= 10 else ""
            if len(b_cols) >= 11:
                try:
                    not_canonical = int(float(b_cols[10])) > 0
                except ValueError:
                    not_canonical = False
            try:
                a_start0 = int(parts[1])
                a_end0 = int(parts[2])
                b_start0 = int(b_cols[1])
                b_end0 = int(b_cols[2])
                overlap_bp = max(0, min(a_end0, b_end0) - max(a_start0, b_start0))
            except (ValueError, IndexError):
                overlap_bp = 0
            row = row_by_id.get(row_id)
            if row is None:
                continue
            event_family = _choose_event_family(row)
            lr_family = _normalize_mei_family_token(f"{subfamily} {hit_id}")
            if not event_family or lr_family != event_family:
                continue
            current = best_hits.get(row_id)
            if (current is None) or (int(current.get("overlap_bp", -1)) < overlap_bp):
                best_hits[row_id] = {
                    "overlap_bp": overlap_bp,
                    "id": hit_id,
                    "event_type": event_type,
                    "subfamily": subfamily,
                    "ins_len": ins_len,
                    "tsd_len": tsd_len,
                    "polya_len": polya_len,
                    "conformation": conformation,
                    "not_canonical": bool(not_canonical),
                }

    if best_hits:
        out["lr_svan_id"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("id", "")))
        out["lr_svan_event_type"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("event_type", "")))
        out["lr_svan_subfamily"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("subfamily", "")))
        out["lr_svan_insertion_length"] = (
            out["row_id"].map(lambda i: int(best_hits.get(i, {}).get("ins_len", -1))).fillna(-1).astype(int)
        )
        out["lr_svan_tsd_len"] = (
            out["row_id"].map(lambda i: int(best_hits.get(i, {}).get("tsd_len", -1))).fillna(-1).astype(int)
        )
        out["lr_svan_polya_len"] = (
            out["row_id"].map(lambda i: int(best_hits.get(i, {}).get("polya_len", -1))).fillna(-1).astype(int)
        )
        out["lr_svan_conformation"] = out["row_id"].map(lambda i: str(best_hits.get(i, {}).get("conformation", "")))
        out["lr_svan_not_canonical"] = (
            out["row_id"].map(lambda i: bool(best_hits.get(i, {}).get("not_canonical", False))).fillna(False).astype(bool)
        )
    return out.drop(columns=["row_id"])


def _add_known_mei_polymorphism_consensus(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    g1k_id = _df_col_series(out, "g1k_melt_id", "").fillna("").astype(str).str.strip()
    lr_id = _df_col_series(out, "lr_svan_id", "").fillna("").astype(str).str.strip()
    has_g1k = g1k_id != ""
    has_lr = lr_id != ""

    out["known_mei_polymorphism"] = has_g1k | has_lr
    out["known_mei_polymorphism_source"] = ""
    out.loc[has_g1k & ~has_lr, "known_mei_polymorphism_source"] = "melt_1kg"
    out.loc[~has_g1k & has_lr, "known_mei_polymorphism_source"] = "long_read_1kg_ont_vienna"
    out.loc[has_g1k & has_lr, "known_mei_polymorphism_source"] = "melt_1kg,long_read_1kg_ont_vienna"

    g1k_family = _df_col_series(out, "g1k_melt_insertion_subfamily", "").fillna("").astype(str).apply(
        lambda x: _infer_mei_family_from_fields(hit_id=x, family_hint=x, extra_hint="")
    )
    lr_family = _df_col_series(out, "lr_svan_subfamily", "").fillna("").astype(str).apply(
        lambda x: _infer_mei_family_from_fields(hit_id=x, family_hint=x, extra_hint="")
    )
    out["known_mei_polymorphism_family"] = ""
    out.loc[has_g1k & ~has_lr, "known_mei_polymorphism_family"] = g1k_family.loc[has_g1k & ~has_lr]
    out.loc[~has_g1k & has_lr, "known_mei_polymorphism_family"] = lr_family.loc[~has_g1k & has_lr]
    both = has_g1k & has_lr
    same_family = both & (g1k_family == lr_family) & (g1k_family != "")
    out.loc[same_family, "known_mei_polymorphism_family"] = g1k_family.loc[same_family]
    out.loc[both & ~same_family, "known_mei_polymorphism_family"] = "MIXED"

    g1k_subfamily = _df_col_series(out, "g1k_melt_insertion_subfamily", "").fillna("").astype(str).str.strip()
    lr_subfamily = _df_col_series(out, "lr_svan_subfamily", "").fillna("").astype(str).str.strip()
    out["known_mei_polymorphism_subfamily"] = ""
    out.loc[has_g1k & ~has_lr, "known_mei_polymorphism_subfamily"] = g1k_subfamily.loc[has_g1k & ~has_lr]
    out.loc[~has_g1k & has_lr, "known_mei_polymorphism_subfamily"] = lr_subfamily.loc[~has_g1k & has_lr]
    same_subfamily = both & (g1k_subfamily == lr_subfamily) & (g1k_subfamily != "")
    out.loc[same_subfamily, "known_mei_polymorphism_subfamily"] = g1k_subfamily.loc[same_subfamily]
    out.loc[both & ~same_subfamily, "known_mei_polymorphism_subfamily"] = "MULTI_SOURCE"

    out["known_mei_polymorphism_id"] = ""
    out.loc[has_g1k & ~has_lr, "known_mei_polymorphism_id"] = g1k_id.loc[has_g1k & ~has_lr]
    out.loc[~has_g1k & has_lr, "known_mei_polymorphism_id"] = lr_id.loc[~has_g1k & has_lr]
    out.loc[both, "known_mei_polymorphism_id"] = (
        "g1k:" + g1k_id.loc[both] + "|lr:" + lr_id.loc[both]
    )
    return out


def _normalize_mei_family_token(token: str) -> str:
    t = (token or "").upper()
    if "ALU" in t:
        return "ALU"
    if "SVA" in t:
        return "SVA"
    if "LINE1" in t or "L1" in t:
        return "LINE1"
    return ""


# Event-level identity: pool disease+control family/subfamily evidence.
# Discordant vote-map columns carry all families (ALU:N,SVA:N), not winner-take-all.
# Split L/R columns still contribute their per-side subfamily × support weights.
_EVENT_DISCORDANT_VOTE_SOURCES: list[tuple[str, str]] = [
    ("disease_discordant_mei_family_votes", "disease_discordant_mei_subfamily_votes"),
    ("control_discordant_mei_family_votes", "control_discordant_mei_subfamily_votes"),
]
_EVENT_SUBFAMILY_WEIGHT_SOURCES: list[tuple[str, str]] = [
    ("disease_L_mei_subfamily", "disease_L_mei_supported_reads"),
    ("disease_R_mei_subfamily", "disease_R_mei_supported_reads"),
    ("control_L_mei_subfamily", "control_L_mei_supported_reads"),
    ("control_R_mei_subfamily", "control_R_mei_supported_reads"),
]
# Fallback only when vote maps are absent (older rows / unit tests).
_EVENT_DISCORDANT_FALLBACK_SOURCES: list[tuple[str, str]] = [
    ("disease_discordant_mei_subfamily", "disease_discordant_mei_supported_reads"),
    ("control_discordant_mei_subfamily", "control_discordant_mei_supported_reads"),
]


def _collect_event_family_and_subfamily_weights(
    row: pd.Series,
) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    """Pool disease+control family and subfamily weights.

    Returns ``(family_weights, subfamily_weights_by_family)`` where family weights
    are ``ALU/SVA/LINE1 -> count`` summed across samples.
    """
    fam_weights: dict[str, int] = {}
    sub_weights: dict[str, dict[str, int]] = {}

    def _add_family(fam: str, weight: int) -> None:
        if not fam or weight <= 0:
            return
        fam_weights[fam] = fam_weights.get(fam, 0) + int(weight)

    def _add_subfamily(label: str, weight: int) -> None:
        fam = _normalize_mei_family_token(label)
        if not fam or weight <= 0:
            return
        _add_family(fam, weight)
        sub_weights.setdefault(fam, {})
        sub_weights[fam][label] = sub_weights[fam].get(label, 0) + int(weight)

    has_vote_maps = False
    for fam_col, sub_col in _EVENT_DISCORDANT_VOTE_SOURCES:
        fam_votes = _parse_vote_map(row.get(fam_col, ""))
        sub_votes = _parse_vote_map(row.get(sub_col, ""))
        if fam_votes or sub_votes:
            has_vote_maps = True
        if fam_votes:
            for fam, weight in fam_votes.items():
                _add_family(_normalize_mei_family_token(fam) or fam, weight)
            for label, weight in sub_votes.items():
                fam = _normalize_mei_family_token(label)
                if not fam or weight <= 0:
                    continue
                sub_weights.setdefault(fam, {})
                sub_weights[fam][label] = sub_weights[fam].get(label, 0) + int(weight)
        else:
            # Subfamily map only: derive family weights from labels.
            for label, weight in sub_votes.items():
                _add_subfamily(label, weight)

    # Split L/R: each side contributes its winning subfamily × side support.
    for subfamily_col, weight_col in _EVENT_SUBFAMILY_WEIGHT_SOURCES:
        label = str(row.get(subfamily_col, "") or "").strip()
        weight = _row_int(row, weight_col)
        if label and weight > 0:
            _add_subfamily(label, weight)

    # Older rows / tests without vote maps: fall back to discordant winner×support.
    if not has_vote_maps:
        for subfamily_col, weight_col in _EVENT_DISCORDANT_FALLBACK_SOURCES:
            label = str(row.get(subfamily_col, "") or "").strip()
            weight = _row_int(row, weight_col)
            if label and weight > 0:
                _add_subfamily(label, weight)

    for col in (
        "g1k_melt_insertion_subfamily",
        "lr_svan_subfamily",
        "known_mei_polymorphism_subfamily",
    ):
        label = str(row.get(col, "") or "").strip()
        if label and label not in {"MULTI_SOURCE"}:
            _add_subfamily(label, 1)

    return fam_weights, sub_weights


def _collect_event_subfamily_weights(row: pd.Series) -> list[tuple[str, int]]:
    """Flatten pooled subfamily weights for callers that expect a label list."""
    _fam_weights, sub_weights = _collect_event_family_and_subfamily_weights(row)
    out: list[tuple[str, int]] = []
    for fam_map in sub_weights.values():
        for label, weight in fam_map.items():
            out.append((label, int(weight)))
    return out


def _choose_event_family_and_subfamily(row: pd.Series) -> tuple[str, str]:
    """Family-first event identity from pooled disease+control evidence.

    1. Sum family weights across samples: disease ALU + control ALU, etc.
    2. Pick the family with the largest pooled weight.
    3. Pick the subfamily only among labels of that family.
    """
    fam_weights, sub_weights = _collect_event_family_and_subfamily_weights(row)
    if fam_weights:
        best_family = max(fam_weights.items(), key=lambda item: item[1])[0]
        fam_subs = sub_weights.get(best_family) or {}
        if fam_subs:
            best_subfamily = max(fam_subs.items(), key=lambda item: item[1])[0]
        else:
            best_subfamily = ""
        return best_family, best_subfamily

    # No resolvable weights: fall back to explicit family columns only.
    for col in (
        "disease_discordant_mei_family",
        "disease_L_mei_family",
        "disease_R_mei_family",
        "control_discordant_mei_family",
        "control_L_mei_family",
        "control_R_mei_family",
        "known_mei_polymorphism_family",
        "g1k_melt_id",
        "lr_svan_id",
    ):
        fam = _normalize_mei_family_token(str(row.get(col, "") or ""))
        if fam:
            return fam, ""
    return "", ""


def _choose_event_subfamily(row: pd.Series) -> str:
    return _choose_event_family_and_subfamily(row)[1]


def _choose_event_family(row: pd.Series) -> str:
    return _choose_event_family_and_subfamily(row)[0]


def _choose_event_orientation(row: pd.Series) -> str:
    candidates = [
        str(row.get("consensus_insertion_orientation", "")),
        str(row.get("insertion_orientation", "")),
        str(row.get("asm_insertion_orientation", "")),
        str(row.get("disease_insertion_orientation", "")),
        str(row.get("control_insertion_orientation", "")),
        str(row.get("disease_discordant_mei_strand", "")),
        str(row.get("disease_L_mei_strand", "")),
        str(row.get("disease_R_mei_strand", "")),
        str(row.get("control_discordant_mei_strand", "")),
        str(row.get("control_L_mei_strand", "")),
        str(row.get("control_R_mei_strand", "")),
    ]
    for c in candidates:
        cc = (c or "").strip()
        if cc in {"+", "-"}:
            return cc
    return ""


def _bed_field(value: object, *, default: str = ".") -> str:
    """Format one BED column for bedtools.

    bedtools type-checking rejects trailing empty fields (``...\\t\\n``). Missing
    optional values must be the BED placeholder ``.``, and embedded tabs/newlines
    must never appear inside a field.
    """
    if value is None:
        text = ""
    elif isinstance(value, float) and pd.isna(value):
        text = ""
    else:
        text = str(value)
        if text.lower() == "nan":
            text = ""
    text = text.replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()
    return text if text else default


def _write_bed_row(handle, fields: list[object]) -> None:
    """Write a BED row with sanitized, non-empty fields."""
    handle.write("\t".join(_bed_field(v) for v in fields) + "\n")


def _write_rmsk_mei_bed(
    rmsk_table_path: Path,
    out_bed: Path,
    *,
    chroms: set[str] | None = None,
) -> int:
    """Write MEI-normalizable rmsk intervals as BED for bedtools intersect.

    Columns: chrom start0 end0 repName length strand repClass repFamily normFamily
    """
    n = 0
    opener = gzip.open if str(rmsk_table_path).endswith(".gz") else open
    with opener(rmsk_table_path, "rt", encoding="utf-8") as hin, out_bed.open("w", encoding="utf-8") as hout:
        for line in hin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                # UCSC rmsk table with leading bin.
                if len(parts) >= 13 and parts[5].startswith("chr"):
                    chrom = parts[5]
                    start0 = int(parts[6])
                    end0 = int(parts[7])
                    strand = parts[9]
                    rep_name = parts[10]
                    rep_class = parts[11]
                    rep_family = parts[12]
                # BED-like fallback: chrom start end ... strand repName repClass repFamily
                elif len(parts) >= 8 and parts[0].startswith("chr"):
                    chrom = parts[0]
                    start0 = int(parts[1])
                    end0 = int(parts[2])
                    strand = parts[5]
                    rep_name = parts[6]
                    rep_class = parts[7]
                    rep_family = parts[8] if len(parts) > 8 else ""
                else:
                    continue
            except (ValueError, IndexError):
                continue
            if chroms is not None and chrom not in chroms:
                continue
            if end0 <= start0:
                continue
            fam = _normalize_mei_family_token(f"{rep_name} {rep_class} {rep_family}")
            if not fam:
                continue
            length = int(end0) - int(start0)
            _write_bed_row(
                hout,
                [chrom, start0, end0, rep_name, length, strand, rep_class, rep_family, fam],
            )
            n += 1
    return n


def _annotate_nested_retrotransposon(candidates: pd.DataFrame, rmsk_table_path: Path) -> pd.DataFrame:
    """Flag candidates whose breakpoint sits in a same-family rmsk MEI (nested).

    Uses ``bedtools intersect`` against an MEI-filtered rmsk BED (chrom-restricted
    to the candidate set). Requires ``bedtools`` on PATH.
    """
    out = candidates.copy().reset_index(drop=True)
    out["nested_repeat_overlap"] = False
    out["nested_repeat_name"] = ""
    out["nested_repeat_class"] = ""
    out["nested_repeat_family"] = ""
    out["nested_repeat_strand"] = ""
    out["nested_mei_family"] = ""
    out["nested_insertion_orientation"] = ""
    out["nested_same_class"] = False
    out["nested_same_orientation"] = False
    out["nested_same_class_orientation"] = "unnested"
    if out.empty:
        return out

    chroms = set(out["chrom"].fillna("").astype(str))
    chroms.discard("")

    bedtools_bin = shutil.which("bedtools")
    if bedtools_bin is None:
        raise RuntimeError(
            "nested-rmsk annotation requires bedtools on PATH "
            "(install bedtools or activate the rtm-miner env)."
        )

    with tempfile.TemporaryDirectory(prefix="rtm_nested_rmsk_") as tmp_dir:
        tmp = Path(tmp_dir)
        rmsk_bed = tmp / "rmsk.mei.bed"
        cand_bed = tmp / "candidate.points.bed"

        bed_t0 = time.monotonic()
        n_rmsk = _write_rmsk_mei_bed(rmsk_table_path, rmsk_bed, chroms=chroms)
        click.echo(
            f"[mei-annotate] nested-rmsk wrote MEI BED rows={n_rmsk} "
            f"chroms={len(chroms)} elapsed={time.monotonic() - bed_t0:.1f}s"
        )
        if n_rmsk == 0:
            raise ValueError(
                "RepeatMasker annotation lacks repName/repClass/repFamily fields required "
                "for nested_same_class_orientation. Use a UCSC rmsk table (e.g. rmsk.txt.gz), "
                "not a stripped 3-4 column BED."
            )

        cand_t0 = time.monotonic()
        with cand_bed.open("w", encoding="utf-8") as hout:
            for i, row in enumerate(out.itertuples(index=False)):
                as_row = pd.Series(row._asdict())
                chrom = str(getattr(row, "chrom"))
                pos_1based = int(getattr(row, "insertion_breakpoint_pos", 0) or 0)
                if pos_1based <= 0:
                    pos_1based = int(
                        (int(getattr(row, "window_start", 1)) + int(getattr(row, "window_end", 1))) // 2
                    )
                pos0 = max(0, pos_1based - 1)
                event_family = _choose_event_family(as_row)
                event_orientation = _choose_event_orientation(as_row)
                # name=row index; keep event family/orient for intersect-side filter.
                # Missing family/orient must be "." — empty trailing BED fields make
                # bedtools exit 1 ("extra TAB at the end of your line").
                _write_bed_row(
                    hout,
                    [chrom, pos0, pos0 + 1, i, 0, ".", event_family, event_orientation],
                )
        click.echo(
            f"[mei-annotate] nested-rmsk wrote candidate point BED "
            f"rows={len(out)} elapsed={time.monotonic() - cand_t0:.1f}s"
        )

        inter_t0 = time.monotonic()
        proc = _run_bedtools_checked(
            [bedtools_bin, "intersect", "-wa", "-wb", "-a", str(cand_bed), "-b", str(rmsk_bed)],
            label="nested-rmsk bedtools intersect",
        )
        click.echo(
            f"[mei-annotate] nested-rmsk bedtools intersect "
            f"elapsed={time.monotonic() - inter_t0:.1f}s"
        )

        # a(8) + b(9): pick best same-family hit (prefer same-orient, then longer).
        post_t0 = time.monotonic()
        best: dict[int, tuple[tuple[int, int], dict[str, object]]] = {}
        for line in proc.stdout.splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 17:
                continue
            try:
                idx = int(parts[3])
            except ValueError:
                continue
            event_fam = parts[6]
            event_orient = parts[7]
            if event_fam in {"", "."}:
                event_fam = ""
            if event_orient in {"", "."}:
                event_orient = ""
            try:
                length = int(parts[12])
            except ValueError:
                length = max(0, int(parts[10]) - int(parts[9]))
            rep_name, strand = parts[11], parts[13]
            rep_class, rep_family, norm_fam = parts[14], parts[15], parts[16]
            if not event_fam or norm_fam != event_fam:
                continue
            same_orient = int(
                event_orient in {"+", "-"} and strand in {"+", "-"} and strand == event_orient
            )
            score = (same_orient, length)
            rec: dict[str, object] = {
                "nested_repeat_overlap": True,
                "nested_repeat_name": rep_name,
                "nested_repeat_class": rep_class,
                "nested_repeat_family": rep_family,
                "nested_repeat_strand": strand,
                "nested_mei_family": event_fam,
                "nested_insertion_orientation": event_orient,
                "nested_same_class": True,
                "nested_same_orientation": bool(same_orient),
                "nested_same_class_orientation": "nested" if same_orient else "unnested",
            }
            prev = best.get(idx)
            if prev is None or score > prev[0]:
                best[idx] = (score, rec)

        if best:
            for idx, (_score, rec) in best.items():
                for key, val in rec.items():
                    out.at[idx, key] = val
        click.echo(
            f"[mei-annotate] nested-rmsk post-select hits={len(best)}/{len(out)} "
            f"elapsed={time.monotonic() - post_t0:.1f}s"
        )
    return out



def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    reference_fasta: Path | None = None,
    disease_bam_path: Path | None = None,
    control_bam_path: Path | None = None,
    disease_mate_bam_path: Path | None = None,
    control_mate_bam_path: Path | None = None,
    rmsk_table_path: Path | None = None,
    g1k_mei_vcf: Path | None = None,
    lr_mei_vcf: Path | None = None,
    g1k_split_padding_bp: int = 200,
    g1k_dpe_padding_min_bp: int = 200,
    g1k_dpe_padding_max_bp: int = 200,
    g1k_dpe_padding_tlen_factor: float = 0.0,
    empirical_stage: bool = False,
    empirical_random_windows: int = 1000,
    empirical_random_scope: str = "chromosome",
    empirical_random_seed: int = 13,
    empirical_highconf_bed: Path | None = None,
    empirical_exclude_merged_bed: Path | None = None,
    empirical_exclude_segdup_bed: Path | None = None,
    empirical_exclude_mappability_bedgraph: Path | None = None,
    empirical_exclude_mappability_threshold: float = 0.5,
    empirical_exclude_gap_bed: Path | None = None,
    empirical_exclude_blacklist_bed: Path | None = None,
    empirical_cache_dir: Path | None = None,
    progress_every: int = 20000,
    igv_plots: bool = True,
    igv_top_n: int = 0,
    igv_snapshot_dir: Path | None = None,
    igv_launcher: Path | None = None,
    igv_gold_only: bool = True,
    igv_panel_height_min: int = 250,
    igv_panel_height_max: int = 8000,
    igv_timeout_sec: int | None = None,
    read_architecture_plots: bool = True,
    read_architecture_top_n: int = 0,
    read_architecture_dir: Path | None = None,
    local_assembly: bool = False,
    assembly_cache_dir: Path | None = None,
    assembly_interval_pad_bp: int = 250,
    assembly_retry_pad_bp: int = 600,
    assembly_max_reads_per_sample: int = 600,
    assembly_spades_threads: int = 1,
    assembly_spades_memory_gb: int = 8,
    assembly_minimap2_threads: int = 1,
    assembly_locus_workers: int = 0,
    assembly_reuse_cache_only: bool = False,
    mei_full_fasta: Path | None = None,
    reuse_mei_annotate_dir: Path | None = None,
    bwa_threads: int = 1,
) -> Path:
    total_t0 = time.monotonic()
    reuse_dir = Path(reuse_mei_annotate_dir) if reuse_mei_annotate_dir is not None else None
    bwa_threads = max(1, int(bwa_threads))
    load_t0 = time.monotonic()
    candidate = pd.read_csv(candidate_loci_path, sep="\t")
    split_disease_raw = _load_table(evidence_dir, "split_evidence", "disease")
    split_control_raw = _load_table(evidence_dir, "split_evidence", "control")
    discordant_disease_raw = _load_table(evidence_dir, "discordant_evidence", "disease")
    discordant_control_raw = _load_table(evidence_dir, "discordant_evidence", "control")
    click.echo(
        f"[mei-annotate] loaded evidence "
        f"candidates={len(candidate)} "
        f"split_d={len(split_disease_raw)} split_c={len(split_control_raw)} "
        f"disc_d={len(discordant_disease_raw)} disc_c={len(discordant_control_raw)} "
        f"elapsed={time.monotonic() - load_t0:.1f}s"
    )

    assign_t0 = time.monotonic()
    split_disease = _assign_rows_to_candidate_loci(split_disease_raw, candidate)
    split_control = _assign_rows_to_candidate_loci(split_control_raw, candidate)
    discordant_disease = _assign_rows_to_candidate_loci(discordant_disease_raw, candidate)
    discordant_control = _assign_rows_to_candidate_loci(discordant_control_raw, candidate)
    click.echo(
        f"[mei-annotate] assigned evidence to loci "
        f"split_d={len(split_disease)} split_c={len(split_control)} "
        f"disc_d={len(discordant_disease)} disc_c={len(discordant_control)} "
        f"elapsed={time.monotonic() - assign_t0:.1f}s"
    )

    if reuse_dir is not None:
        # Re-label path: hydrate MEI hits from a prior annotate detail table.
        # Skips indel BAM scan, mate fetch, and MEI consensus remaps.
        detail_path = _resolve_supporting_reads_detail_path(reuse_dir)
        click.echo(
            f"[mei-annotate] reuse-mei-annotate-dir={reuse_dir} "
            f"detail={detail_path.name} (skip indel+remap)"
        )
        reuse_detail = _load_supporting_reads_detail_table(detail_path)
        indel_disease = pd.DataFrame()
        indel_control = pd.DataFrame()
        remap_t0 = time.monotonic()
        remap_by_sample = {
            "disease": _hydrate_sample_mei_hits_from_detail(
                sample="disease",
                split_df=split_disease,
                discordant_df=discordant_disease,
                detail=reuse_detail,
            ),
            "control": _hydrate_sample_mei_hits_from_detail(
                sample="control",
                split_df=split_control,
                discordant_df=discordant_control,
                detail=reuse_detail,
            ),
        }
        click.echo(
            f"[mei-annotate] hydrated MEI hits from detail wall elapsed={time.monotonic() - remap_t0:.1f}s"
        )
        # Re-apply sequence-based rescues (mate_seq is on discordant evidence; prior
        # detail often omitted polyA-only DPE rows, so hydrate alone cannot restore them).
        for sample in ("disease", "control"):
            disc = remap_by_sample[sample]["disc_hits"]  # type: ignore[index]
            if isinstance(disc, pd.DataFrame) and not disc.empty:
                disc = _rescue_vntr_like_discordant_mei_hits(disc)
                disc = _rescue_polya_like_discordant_mei_hits(disc)
                remap_by_sample[sample]["disc_hits"] = disc  # type: ignore[index]
    else:
        # Indel collection + MEI remaps: disease∥control (I/O and bwa mem release the GIL).
        click.echo("[mei-annotate] running disease∥control indel + MEI remaps")
        indel_t0 = time.monotonic()
        indel_jobs: dict[str, object] = {}
        with ThreadPoolExecutor(max_workers=2) as pool:
            indel_futs = {}
            if disease_bam_path is not None:
                indel_futs[
                    pool.submit(
                        _collect_indel_breakpoint_evidence,
                        disease_bam_path,
                        candidate,
                        sample="disease",
                    )
                ] = "disease"
            if control_bam_path is not None:
                indel_futs[
                    pool.submit(
                        _collect_indel_breakpoint_evidence,
                        control_bam_path,
                        candidate,
                        sample="control",
                    )
                ] = "control"
            for fut in as_completed(indel_futs):
                indel_jobs[indel_futs[fut]] = fut.result()
        indel_disease = indel_jobs.get("disease", pd.DataFrame())
        indel_control = indel_jobs.get("control", pd.DataFrame())
        click.echo(
            f"[mei-annotate] indel collection disease={len(indel_disease)} "
            f"control={len(indel_control)} elapsed={time.monotonic() - indel_t0:.1f}s"
        )

        remap_t0 = time.monotonic()
        # disease∥control remaps can run together; split bwa threads across them
        # when the caller asked for more than one thread.
        per_sample_bwa_threads = max(1, bwa_threads // 2) if bwa_threads > 1 else 1
        click.echo(
            f"[mei-annotate] bwa_threads total={bwa_threads} "
            f"per_sample={per_sample_bwa_threads} (disease∥control)"
        )
        remap_by_sample = {}
        with ThreadPoolExecutor(max_workers=2) as pool:
            remap_futs = {
                pool.submit(
                    _remap_one_sample_mei_evidence,
                    sample="disease",
                    split_df=split_disease,
                    discordant_df=discordant_disease,
                    mei_fasta=mei_fasta,
                    bam_path=disease_bam_path,
                    mate_bam_path=disease_mate_bam_path,
                    bwa_threads=per_sample_bwa_threads,
                ): "disease",
                pool.submit(
                    _remap_one_sample_mei_evidence,
                    sample="control",
                    split_df=split_control,
                    discordant_df=discordant_control,
                    mei_fasta=mei_fasta,
                    bam_path=control_bam_path,
                    mate_bam_path=control_mate_bam_path,
                    bwa_threads=per_sample_bwa_threads,
                ): "control",
            }
            for fut in as_completed(remap_futs):
                sample = remap_futs[fut]
                remap_by_sample[sample] = fut.result()
        click.echo(
            f"[mei-annotate] disease∥control MEI remaps wall elapsed={time.monotonic() - remap_t0:.1f}s"
        )

    disease_hits = remap_by_sample["disease"]["split_hits"]  # type: ignore[assignment]
    control_hits = remap_by_sample["control"]["split_hits"]  # type: ignore[assignment]
    disease_summary = remap_by_sample["disease"]["split_summary"]  # type: ignore[assignment]
    control_summary = remap_by_sample["control"]["split_summary"]  # type: ignore[assignment]
    disease_disc_hits = remap_by_sample["disease"]["disc_hits"]  # type: ignore[assignment]
    control_disc_hits = remap_by_sample["control"]["disc_hits"]  # type: ignore[assignment]
    disease_disc_summary = remap_by_sample["disease"]["disc_summary"]  # type: ignore[assignment]
    control_disc_summary = remap_by_sample["control"]["disc_summary"]  # type: ignore[assignment]
    disease_disc_mate_summary = remap_by_sample["disease"]["disc_mate_summary"]  # type: ignore[assignment]
    control_disc_mate_summary = remap_by_sample["control"]["disc_mate_summary"]  # type: ignore[assignment]
    detail_t0 = time.monotonic()
    # Always rebuild from hit frames so short-MEI rescue / polyA split evidence is included
    # even when MEI hits were hydrated from a prior detail table (reuse mode).
    supporting_reads_detail = pd.concat(
        [
            _build_supporting_reads_detail_table(
                split_hits=disease_hits,
                discordant_hits=disease_disc_hits,
                discordant_mate_hits=pd.DataFrame(),
                sample="disease",
            ),
            _build_supporting_reads_detail_table(
                split_hits=control_hits,
                discordant_hits=control_disc_hits,
                discordant_mate_hits=pd.DataFrame(),
                sample="control",
            ),
        ],
        ignore_index=True,
    )
    if not supporting_reads_detail.empty:
        detail_tsv = out_path.with_name("supporting_reads_detail.mei.tsv")
        detail_parquet = out_path.with_name("supporting_reads_detail.mei.parquet")
        supporting_reads_detail.to_csv(detail_tsv, sep="\t", index=False)
        supporting_reads_detail.to_parquet(detail_parquet, index=False)
        candidate = _merge_detail_mei_extents(candidate, supporting_reads_detail)
        click.echo(
            f"[mei-annotate] wrote supporting read detail table to {detail_tsv} "
            f"rows={len(supporting_reads_detail)} elapsed={time.monotonic() - detail_t0:.1f}s"
        )
    else:
        click.echo(
            f"[mei-annotate] supporting read detail empty elapsed={time.monotonic() - detail_t0:.1f}s"
        )
    frag_map_tsv = _resolve_fragment_to_full_map_tsv(mei_fasta=mei_fasta, mei_full_fasta=mei_full_fasta)
    frag_map = _load_fragment_to_full_map(frag_map_tsv)
    # Prefer one-time fragment→full projection over per-read full-consensus remaps.
    disease_hits_full = pd.DataFrame()
    control_hits_full = pd.DataFrame()
    disease_disc_hits_full = pd.DataFrame()
    control_disc_hits_full = pd.DataFrame()
    disease_summary_full = ClipAlignmentSummary(sample="disease_full", clip_count=0, paf_hits=0)
    control_summary_full = ClipAlignmentSummary(sample="control_full", clip_count=0, paf_hits=0)
    disease_disc_summary_full = ClipAlignmentSummary(sample="disease_full", clip_count=0, paf_hits=0)
    control_disc_summary_full = ClipAlignmentSummary(sample="control_full", clip_count=0, paf_hits=0)
    if frag_map:
        # Project detail for plots; then refresh candidate detail extents on the
        # full-length axis (panel extents were merged above before projection).
        if not supporting_reads_detail.empty:
            detail_full = _project_detail_coords_to_full(supporting_reads_detail, frag_map)
            detail_tsv = out_path.with_name("supporting_reads_detail.mei.tsv")
            detail_parquet = out_path.with_name("supporting_reads_detail.mei.parquet")
            detail_full.to_csv(detail_tsv, sep="\t", index=False)
            detail_full.to_parquet(detail_parquet, index=False)
            supporting_reads_detail = detail_full
            # Re-aggregate without panel target-length caps so L1 5′+3′ unions
            # survive onto disease/control_detail_mei_* used by gold.
            candidate = _merge_detail_mei_extents(candidate, supporting_reads_detail)
        n_frags = len({e.fragment_name for e in frag_map.values()})
        click.echo(
            f"[mei-annotate] projected panel MEI coords onto full-length consensus via "
            f"{frag_map_tsv.name if frag_map_tsv is not None else 'fragment map'} "
            f"(entries={n_frags})"
        )
    else:
        full_consensus_fasta = _resolve_full_consensus_fasta(
            mei_fasta=mei_fasta,
            out_dir=out_path.parent,
            mei_full_fasta=mei_full_fasta,
        )
        if reuse_dir is not None:
            click.echo(
                "[mei-annotate] reuse mode: skipping full-consensus remaps "
                "(no fragment map; keeping prior detail coords)"
            )
        elif full_consensus_fasta is not None:
            click.echo(
                "[mei-annotate] warning: mei_fragment_to_full_coords.tsv missing; "
                "falling back to per-read full-consensus remaps. Re-run "
                "scripts/download_public_data.py to build the fragment map."
            )
            disease_hits_full, disease_summary_full = _align_clips_with_minimap2(
                split_disease,
                full_consensus_fasta,
                sample="disease_full",
                bwa_threads=bwa_threads,
            )
            control_hits_full, control_summary_full = _align_clips_with_minimap2(
                split_control,
                full_consensus_fasta,
                sample="control_full",
                bwa_threads=bwa_threads,
            )
            disease_disc_anchor_hits_full, disease_disc_summary_full = _align_discordant_reads_with_minimap2(
                discordant_disease,
                full_consensus_fasta,
                sample="disease_full",
                bwa_threads=bwa_threads,
            )
            control_disc_anchor_hits_full, control_disc_summary_full = _align_discordant_reads_with_minimap2(
                discordant_control,
                full_consensus_fasta,
                sample="control_full",
                bwa_threads=bwa_threads,
            )
            disease_disc_mate_hits_full, _disease_disc_mate_summary_full = _align_discordant_mates_with_minimap2(
                discordant_disease,
                full_consensus_fasta,
                sample="disease_full_mate",
                bam_path=disease_mate_bam_path or disease_bam_path,
                bwa_threads=bwa_threads,
            )
            control_disc_mate_hits_full, _control_disc_mate_summary_full = _align_discordant_mates_with_minimap2(
                discordant_control,
                full_consensus_fasta,
                sample="control_full_mate",
                bam_path=control_mate_bam_path or control_bam_path,
                bwa_threads=bwa_threads,
            )
            disease_disc_hits_full = _attach_mei_hits_to_discordant_rows(
                discordant_disease, disease_disc_anchor_hits_full, disease_disc_mate_hits_full
            )
            control_disc_hits_full = _attach_mei_hits_to_discordant_rows(
                discordant_control, control_disc_anchor_hits_full, control_disc_mate_hits_full
            )
            if not supporting_reads_detail.empty:
                full_detail = pd.concat(
                    [
                        _build_supporting_reads_detail_table(
                            split_hits=disease_hits_full,
                            discordant_hits=disease_disc_hits_full,
                            discordant_mate_hits=pd.DataFrame(),
                            sample="disease",
                        ),
                        _build_supporting_reads_detail_table(
                            split_hits=control_hits_full,
                            discordant_hits=control_disc_hits_full,
                            discordant_mate_hits=pd.DataFrame(),
                            sample="control",
                        ),
                    ],
                    ignore_index=True,
                )
                supporting_reads_detail = _overlay_full_consensus_coords_onto_detail(
                    supporting_reads_detail, full_detail
                )
                detail_tsv = out_path.with_name("supporting_reads_detail.mei.tsv")
                detail_parquet = out_path.with_name("supporting_reads_detail.mei.parquet")
                supporting_reads_detail.to_csv(detail_tsv, sep="\t", index=False)
                supporting_reads_detail.to_parquet(detail_parquet, index=False)
                candidate = _merge_detail_mei_extents(candidate, supporting_reads_detail)
                click.echo(
                    f"[mei-annotate] overlaid full-consensus coords onto supporting read detail "
                    f"({full_consensus_fasta.name})"
                )

    click.echo(
        f"[mei-annotate] disease clips={disease_summary.clip_count} hits={disease_summary.paf_hits}; "
        f"control clips={control_summary.clip_count} hits={control_summary.paf_hits}; "
        f"disease discordant reads={disease_disc_summary.clip_count} hits={disease_disc_summary.paf_hits}; "
        f"control discordant reads={control_disc_summary.clip_count} hits={control_disc_summary.paf_hits}; "
        f"disease discordant mates={disease_disc_mate_summary.clip_count} mei_hits={disease_disc_mate_summary.paf_hits}; "
        f"control discordant mates={control_disc_mate_summary.clip_count} mei_hits={control_disc_mate_summary.paf_hits}; "
        f"disease full clips={disease_summary_full.clip_count} hits={disease_summary_full.paf_hits}; "
        f"control full clips={control_summary_full.clip_count} hits={control_summary_full.paf_hits}; "
        f"disease full discordant reads={disease_disc_summary_full.clip_count} hits={disease_disc_summary_full.paf_hits}; "
        f"control full discordant reads={control_disc_summary_full.clip_count} hits={control_disc_summary_full.paf_hits}; "
        f"disease indel reads={len(indel_disease)} control indel reads={len(indel_control)}"
    )

    def _mei_rows_only(df: pd.DataFrame, *, is_split: bool = False) -> pd.DataFrame:
        if ("mei_hit" not in df.columns) and ("mei_hit_coord" not in df.columns) and (
            "mate_mei_hit" not in df.columns
        ):
            return df.iloc[0:0].copy()
        work = df
        if is_split:
            # Short-clip rescue should already be annotated during remap; re-run
            # without DPE only if the flag is missing (e.g. unit tests).
            if "short_mei_seed_rescued" not in work.columns:
                work = _annotate_short_mei_seed_rescue(work)
            keep = _split_mei_support_eligible_mask(work)
        else:
            keep = _discordant_row_mei_mapped(work)
        # Exclude second-pass rescues from MEI_MAPPED row sets.
        if "polya_rescue" in work.columns:
            keep = keep & ~work["polya_rescue"].fillna(False).astype(bool)
        if "vntr_rescue" in work.columns:
            keep = keep & ~work["vntr_rescue"].fillna(False).astype(bool)
        if "mei_hit_source" in work.columns:
            src = work["mei_hit_source"].fillna("").astype(str)
            keep = keep & ~src.isin({"polya_rescue", "vntr_rescue"})
        return work.loc[keep].copy()

    split_disease_mei = _mei_rows_only(disease_hits, is_split=True)
    split_control_mei = _mei_rows_only(control_hits, is_split=True)
    # MEI_MAPPED for DPE: exclude same-chr mates within 1 kb (nearby ref MEIs).
    discordant_disease_mei = _discordant_rows_for_mei_mapped_support(
        _mei_rows_only(disease_disc_hits, is_split=False)
    )
    discordant_control_mei = _discordant_rows_for_mei_mapped_support(
        _mei_rows_only(control_disc_hits, is_split=False)
    )

    metrics_t0 = time.monotonic()
    disc_t = _aggregate_discordant_mei_metrics(disease_disc_hits, sample_prefix="disease")
    disc_n = _aggregate_discordant_mei_metrics(control_disc_hits, sample_prefix="control")
    disc_t_full = _aggregate_discordant_mei_metrics(disease_disc_hits_full, sample_prefix="disease_full")
    disc_n_full = _aggregate_discordant_mei_metrics(control_disc_hits_full, sample_prefix="control_full")
    disc_anchor_t = _aggregate_discordant_anchor_side_metrics(disease_disc_hits, sample_prefix="disease")
    disc_anchor_n = _aggregate_discordant_anchor_side_metrics(control_disc_hits, sample_prefix="control")
    disc_residual_t = _aggregate_discordant_residual_complex_metrics(disease_disc_hits, sample_prefix="disease")
    disc_residual_n = _aggregate_discordant_residual_complex_metrics(control_disc_hits, sample_prefix="control")
    click.echo(
        f"[mei-annotate] aggregated discordant MEI metrics "
        f"elapsed={time.monotonic() - metrics_t0:.1f}s"
    )

    def _preferred_map_from_discordant(df: pd.DataFrame, target_col: str) -> dict[tuple[str, int, int], str]:
        if df.empty or (target_col not in df.columns):
            return {}
        subset = df.loc[:, ["chrom", "window_start", "window_end", target_col]].copy()
        subset[target_col] = subset[target_col].fillna("").astype(str)
        subset = subset.loc[subset[target_col].str.len() > 0]
        out_map: dict[tuple[str, int, int], str] = {}
        for row in subset.itertuples(index=False):
            out_map[(str(row.chrom), int(row.window_start), int(row.window_end))] = str(getattr(row, target_col))
        return out_map

    disease_pref_target = _preferred_map_from_discordant(disc_t, "disease_discordant_mei_subfamily")
    control_pref_target = _preferred_map_from_discordant(disc_n, "control_discordant_mei_subfamily")
    disease_pref_target_full = _preferred_map_from_discordant(disc_t_full, "disease_full_discordant_mei_subfamily")
    control_pref_target_full = _preferred_map_from_discordant(disc_n_full, "control_full_discordant_mei_subfamily")

    anno_parts = []
    side_t0 = time.monotonic()
    for sample_prefix, df, pref_map in (
        ("disease", disease_hits, disease_pref_target),
        ("control", control_hits, control_pref_target),
        ("disease_full", disease_hits_full, disease_pref_target_full),
        ("control_full", control_hits_full, control_pref_target_full),
    ):
        for side in ("L", "R"):
            anno_parts.append(
                _aggregate_side_metrics(
                    df,
                    sample_prefix=sample_prefix,
                    side=side,
                    preferred_subfamily_by_locus=pref_map,
                )
            )
    click.echo(
        f"[mei-annotate] aggregated side metrics parts={len(anno_parts)} "
        f"elapsed={time.monotonic() - side_t0:.1f}s"
    )

    merge_t0 = time.monotonic()
    for idx, part in enumerate(anno_parts):
        if part.empty:
            continue
        candidate = candidate.merge(part, on=["chrom", "window_start", "window_end"], how="left")
        if (idx + 1) % 2 == 0:
            click.echo(f"[mei-annotate] merged side metrics {idx + 1}/{len(anno_parts)}")

    if not disc_t.empty:
        candidate = candidate.merge(disc_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_n.empty:
        candidate = candidate.merge(disc_n, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_t_full.empty:
        candidate = candidate.merge(disc_t_full, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_n_full.empty:
        candidate = candidate.merge(disc_n_full, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_anchor_t.empty:
        candidate = candidate.merge(disc_anchor_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_anchor_n.empty:
        candidate = candidate.merge(disc_anchor_n, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_residual_t.empty:
        candidate = candidate.merge(disc_residual_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_residual_n.empty:
        candidate = candidate.merge(disc_residual_n, on=["chrom", "window_start", "window_end"], how="left")
    click.echo(
        f"[mei-annotate] merged discordant MEI support metrics "
        f"elapsed={time.monotonic() - merge_t0:.1f}s"
    )

    support_t0 = time.monotonic()
    candidate = _add_candidate_support_info_fields(
        candidate,
        split_disease=split_disease,
        split_control=split_control,
        # Pass hit-annotated frames so polyA_MAPPED / VNTR_MAPPED flags are visible.
        discordant_disease=disease_disc_hits,
        discordant_control=control_disc_hits,
        split_disease_mei=split_disease_mei,
        split_control_mei=split_control_mei,
        discordant_disease_mei=discordant_disease_mei,
        discordant_control_mei=discordant_control_mei,
        indel_disease=indel_disease,
        indel_control=indel_control,
    )
    click.echo(
        f"[mei-annotate] added candidate support info fields "
        f"elapsed={time.monotonic() - support_t0:.1f}s"
    )

    for col in candidate.columns:
        if re.search(r"_mei_supported_reads$|_mei_start$|_mei_end$", col):
            candidate[col] = candidate[col].fillna(0).astype(int)
        if re.search(r"_mei_breakpoint_mode$|_mei_breakpoint_unique_positions$", col):
            candidate[col] = candidate[col].fillna(0).astype(int)
        if col.endswith("_mei_score_sum"):
            candidate[col] = candidate[col].fillna(0.0).astype(float)
        if col.endswith("_mei_subfamily_purity") or col.endswith("_mei_breakpoint_mode_fraction"):
            candidate[col] = candidate[col].fillna(0.0).astype(float)
        if re.search(r"_mei_family$|_mei_subfamily$|_mei_strand$", col):
            candidate[col] = candidate[col].fillna("")

    # De-fragment frame before adding many derived columns to avoid pandas PerformanceWarning.
    candidate = _ensure_candidate_schema_defaults(candidate.copy())

    candidate["disease_split_mei_score_sum"] = candidate.get("disease_L_mei_score_sum", 0.0) + candidate.get(
        "disease_R_mei_score_sum", 0.0
    )
    candidate["control_split_mei_score_sum"] = candidate.get("control_L_mei_score_sum", 0.0) + candidate.get(
        "control_R_mei_score_sum", 0.0
    )
    candidate["disease_split_mei_supported_reads"] = candidate.get("disease_L_mei_supported_reads", 0) + candidate.get(
        "disease_R_mei_supported_reads", 0
    )
    candidate["control_split_mei_supported_reads"] = candidate.get("control_L_mei_supported_reads", 0) + candidate.get(
        "control_R_mei_supported_reads", 0
    )
    candidate["disease_discordant_mei_supported_reads"] = (
        candidate.get("disease_discordant_mei_supported_reads", pd.Series(0, index=candidate.index)).fillna(0).astype(int)
    )
    candidate["control_discordant_mei_supported_reads"] = (
        candidate.get("control_discordant_mei_supported_reads", pd.Series(0, index=candidate.index)).fillna(0).astype(int)
    )
    candidate["disease_discordant_mei_score_sum"] = (
        candidate.get("disease_discordant_mei_score_sum", pd.Series(0.0, index=candidate.index)).fillna(0.0).astype(float)
    )
    candidate["control_discordant_mei_score_sum"] = (
        candidate.get("control_discordant_mei_score_sum", pd.Series(0.0, index=candidate.index)).fillna(0.0).astype(float)
    )
    candidate["disease_full_split_mei_supported_reads"] = candidate.get("disease_full_L_mei_supported_reads", 0) + candidate.get(
        "disease_full_R_mei_supported_reads", 0
    )
    candidate["control_full_split_mei_supported_reads"] = candidate.get("control_full_L_mei_supported_reads", 0) + candidate.get(
        "control_full_R_mei_supported_reads", 0
    )
    candidate["disease_full_discordant_mei_supported_reads"] = (
        candidate.get("disease_full_discordant_mei_supported_reads", pd.Series(0, index=candidate.index)).fillna(0).astype(int)
    )
    candidate["control_full_discordant_mei_supported_reads"] = (
        candidate.get("control_full_discordant_mei_supported_reads", pd.Series(0, index=candidate.index)).fillna(0).astype(int)
    )
    candidate["disease_full_mei_supported_reads"] = (
        candidate["disease_full_split_mei_supported_reads"] + candidate["disease_full_discordant_mei_supported_reads"]
    )
    candidate["control_full_mei_supported_reads"] = (
        candidate["control_full_split_mei_supported_reads"] + candidate["control_full_discordant_mei_supported_reads"]
    )

    candidate["disease_mei_score_sum"] = candidate["disease_split_mei_score_sum"] + candidate["disease_discordant_mei_score_sum"]
    candidate["control_mei_score_sum"] = candidate["control_split_mei_score_sum"] + candidate["control_discordant_mei_score_sum"]
    candidate["disease_mei_supported_reads"] = (
        candidate["disease_split_mei_supported_reads"] + candidate["disease_discordant_mei_supported_reads"]
    )
    candidate["control_mei_supported_reads"] = (
        candidate["control_split_mei_supported_reads"] + candidate["control_discordant_mei_supported_reads"]
    )
    candidate["mei_score_enrichment_ratio"] = (candidate["disease_mei_score_sum"] + 0.1) / (
        candidate["control_mei_score_sum"] + 0.1
    )

    candidate = _add_local_depth_normalized_support(candidate)
    candidate = _infer_disease_insertion_metrics(
        candidate,
        reference_fasta=reference_fasta,
        split_disease=split_disease,
        split_control=split_control,
    )
    candidate = _apply_discordant_gap_breakpoint_fallback(
        candidate,
        discordant_disease=discordant_disease_mei,
        discordant_control=discordant_control_mei,
    )
    candidate = _derive_breakpoint_interval_fields(
        candidate,
        breakpoint_pos_col="insertion_breakpoint_pos",
        output_prefix="insertion_",
    )
    for full_prefix in ("disease_full", "control_full"):
        full_metrics = candidate.apply(
            lambda r: _sample_insertion_span_and_orientation(r, full_prefix),
            axis=1,
            result_type="expand",
        )
        full_metrics.columns = [
            f"{full_prefix}_insertion_mei_start",
            f"{full_prefix}_insertion_mei_end",
            f"{full_prefix}_insertion_mei_span",
            f"{full_prefix}_insertion_orientation",
        ]
        for col in full_metrics.columns:
            candidate[col] = full_metrics[col]
    candidate = _compute_insertion_model_scores(candidate)
    stage_t0 = time.monotonic()
    candidate = _assign_bronze_silver_stages(candidate)
    click.echo(
        f"[mei-annotate] bronze/silver staging elapsed={time.monotonic() - stage_t0:.1f}s"
    )
    if local_assembly and disease_bam_path is not None and control_bam_path is not None:
        asm_t0 = time.monotonic()
        asm_dir = assembly_cache_dir if assembly_cache_dir is not None else out_path.parent / "assembly_cache"
        disease_preferred_read_names_by_locus = _build_locus_read_name_map(
            pd.concat([split_disease, discordant_disease], ignore_index=True)
        )
        control_preferred_read_names_by_locus = _build_locus_read_name_map(
            pd.concat([split_control, discordant_control], ignore_index=True)
        )
        asm_df = annotate_silver_with_local_assembly(
            candidate,
            disease_bam_path=disease_bam_path,
            control_bam_path=control_bam_path,
            assembly_cache_dir=asm_dir,
            mei_fasta=mei_fasta,
            reference_fasta=reference_fasta,
            interval_pad_bp=assembly_interval_pad_bp,
            retry_pad_bp=assembly_retry_pad_bp,
            max_reads_per_sample=assembly_max_reads_per_sample,
            spades_threads=assembly_spades_threads,
            spades_memory_gb=assembly_spades_memory_gb,
            minimap2_threads=assembly_minimap2_threads,
            locus_workers=assembly_locus_workers,
            reuse_cache_only=assembly_reuse_cache_only,
            disease_preferred_read_names_by_locus=disease_preferred_read_names_by_locus,
            control_preferred_read_names_by_locus=control_preferred_read_names_by_locus,
        )
        if not asm_df.empty:
            candidate = candidate.merge(asm_df, on=["chrom", "window_start", "window_end"], how="left")
            candidate = _apply_assembly_refinement_overrides(candidate)
            candidate = _recompute_breakpoint_sequence_metrics(candidate, reference_fasta=reference_fasta)
            candidate = _derive_breakpoint_interval_fields(
                candidate,
                breakpoint_pos_col="insertion_breakpoint_pos",
                output_prefix="insertion_",
            )
        click.echo(
            f"[mei-annotate] local assembly complete loci={len(asm_df)} "
            f"cache={asm_dir} elapsed={time.monotonic() - asm_t0:.1f}s"
        )
    candidate = _add_post_assembly_support_info_fields(
        candidate,
        split_disease=split_disease,
        split_control=split_control,
        discordant_disease=discordant_disease,
        discordant_control=discordant_control,
    )
    if g1k_mei_vcf is not None:
        g1k_t0 = time.monotonic()
        candidate = _annotate_g1k_mei_overlap(
            candidate,
            g1k_mei_vcf=g1k_mei_vcf,
            split_padding_bp=g1k_split_padding_bp,
            dpe_padding_min_bp=g1k_dpe_padding_min_bp,
            dpe_padding_max_bp=g1k_dpe_padding_max_bp,
            dpe_padding_tlen_factor=g1k_dpe_padding_tlen_factor,
        )
        click.echo(
            f"[mei-annotate] added 1000G/MELT polymorphism overlap fields "
            f"(elapsed={time.monotonic() - g1k_t0:.1f}s)"
        )
    if lr_mei_vcf is not None:
        lr_t0 = time.monotonic()
        candidate = _annotate_lr_mei_overlap(
            candidate,
            lr_mei_vcf=lr_mei_vcf,
            split_padding_bp=g1k_split_padding_bp,
            dpe_padding_min_bp=g1k_dpe_padding_min_bp,
            dpe_padding_max_bp=g1k_dpe_padding_max_bp,
            dpe_padding_tlen_factor=g1k_dpe_padding_tlen_factor,
        )
        click.echo(
            f"[mei-annotate] added long-read SVAN polymorphism overlap fields "
            f"(elapsed={time.monotonic() - lr_t0:.1f}s)"
        )
    candidate = _add_known_mei_polymorphism_consensus(candidate)
    candidate = _add_consolidated_event_fields(candidate)
    candidate = _broaden_poly_at_fields(candidate)
    if rmsk_table_path is not None:
        rmsk_t0 = time.monotonic()
        candidate = _annotate_nested_retrotransposon(candidate, rmsk_table_path=rmsk_table_path)
        click.echo(
            f"[mei-annotate] added nested-retrotransposon overlap annotation "
            f"(elapsed={time.monotonic() - rmsk_t0:.1f}s)"
        )
    if empirical_stage and disease_bam_path is not None and control_bam_path is not None:
        emp_t0 = time.monotonic()
        candidate = _annotate_bam_depth_for_consistent_loci(
            candidate,
            disease_bam_path=disease_bam_path,
            control_bam_path=control_bam_path,
            empirical_random_windows=empirical_random_windows,
            empirical_random_scope=empirical_random_scope,
            empirical_random_seed=empirical_random_seed,
            empirical_highconf_bed=empirical_highconf_bed,
            empirical_exclude_merged_bed=empirical_exclude_merged_bed,
            empirical_exclude_segdup_bed=empirical_exclude_segdup_bed,
            empirical_exclude_mappability_bedgraph=empirical_exclude_mappability_bedgraph,
            empirical_exclude_mappability_threshold=empirical_exclude_mappability_threshold,
            empirical_exclude_gap_bed=empirical_exclude_gap_bed,
            empirical_exclude_blacklist_bed=empirical_exclude_blacklist_bed,
            empirical_cache_dir=empirical_cache_dir if empirical_cache_dir is not None else out_path.parent / "empirical_cache",
        )
        click.echo(
            f"[mei-annotate] added BAM-depth controlization for family-consistent, junk-clean loci "
            f"(elapsed={time.monotonic() - emp_t0:.1f}s)"
        )
    elif not empirical_stage:
        click.echo("[mei-annotate] empirical stage disabled (--no-empirical-stage)")
    candidate = _add_heuristic_assembly_like_vaf_fields(candidate)
    candidate = _assign_gold_stage(candidate, empirical_stage=empirical_stage)

    candidate = _apply_breakpoint_motif_report_gating(candidate)
    candidate = _prioritize_mei_candidates(candidate, stage_first=True)

    # Final window shrink to resolved breakpoint interval (no pad). Keep discovery
    # windows as merge keys until here so all upstream joins stay stable.
    if "insertion_breakpoint_interval_start" not in candidate.columns:
        candidate = _derive_breakpoint_interval_fields(
            candidate,
            breakpoint_pos_col="insertion_breakpoint_pos",
            output_prefix="insertion_",
        )
    candidate = _tighten_windows_to_breakpoint_interval(
        candidate,
        breakpoint_pos_col="insertion_breakpoint_pos",
        interval_start_col="insertion_breakpoint_interval_start",
        interval_end_col="insertion_breakpoint_interval_end",
    )

    candidate_tsv = _stable_tsv_export_frame(candidate)
    write_t0 = time.monotonic()
    candidate_tsv.to_csv(out_path, sep="\t", index=False)
    candidate.to_parquet(out_path.with_suffix(".parquet"), index=False)
    gold_review = _build_gold_review_table(candidate, empirical_stage=empirical_stage, fragment_to_full_map=frag_map)
    gold_review_path = out_path.with_name(out_path.stem + ".gold_review.tsv")
    gold_review_tsv = _stable_tsv_export_frame(gold_review)
    gold_review_tsv.to_csv(gold_review_path, sep="\t", index=False)
    click.echo(
        f"[mei-annotate] wrote gold review table to {gold_review_path} "
        f"(candidate_rows={len(candidate)} gold_rows={len(gold_review)} "
        f"elapsed={time.monotonic() - write_t0:.1f}s)"
    )
    if (
        igv_plots
        and disease_bam_path is not None
        and control_bam_path is not None
        and reference_fasta is not None
    ):
        igv_dir = igv_snapshot_dir if igv_snapshot_dir is not None else out_path.with_name(out_path.stem + ".gold_review.igv")
        try:
            generate_gold_review_igv_plots(
                gold_review,
                reference_fasta=reference_fasta,
                disease_bam=disease_bam_path,
                control_bam=control_bam_path,
                snapshot_dir=igv_dir,
                top_n=igv_top_n,
                gold_only=igv_gold_only,
                launcher=igv_launcher,
                panel_height_min=igv_panel_height_min,
                panel_height_max=igv_panel_height_max,
                timeout_sec=igv_timeout_sec,
                assembly_cache_dir=asm_dir if local_assembly else None,
            )
        except FileNotFoundError as exc:
            click.echo(f"[mei-annotate] IGV snapshot generation skipped: {exc}")
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            click.echo(
                f"[mei-annotate] IGV snapshot generation failed: {exc}; "
                f"batch script remains at {igv_dir / 'igv_batch.txt'}"
            )
    elif igv_plots and (disease_bam_path is None or control_bam_path is None or reference_fasta is None):
        click.echo(
            "[mei-annotate] IGV snapshot generation skipped: require --reference-fasta, "
            "--disease-bam-depth, and --control-bam-depth"
        )
    if read_architecture_plots:
        detail_path = out_path.with_name("supporting_reads_detail.mei.tsv")
        detail_parquet = out_path.with_name("supporting_reads_detail.mei.parquet")
        detail_source: pd.DataFrame | Path | None = None
        if not supporting_reads_detail.empty:
            detail_source = supporting_reads_detail
        elif detail_parquet.exists():
            detail_source = detail_parquet
        elif detail_path.exists():
            detail_source = detail_path
        if detail_source is None:
            click.echo(
                "[mei-annotate] read-architecture plots skipped: missing supporting_reads_detail.mei.tsv"
            )
        else:
            arch_dir = (
                read_architecture_dir
                if read_architecture_dir is not None
                else out_path.with_name(out_path.stem + ".read_architecture")
            )
            try:
                generate_gold_read_architecture_plots(
                    gold_review,
                    supporting_reads_detail=detail_source,
                    out_dir=arch_dir,
                    mei_table=candidate,
                    gold_tsv=gold_review_path,
                    gold_only=True,
                    top_n=read_architecture_top_n,
                )
            except Exception as exc:  # noqa: BLE001 - do not fail annotate on plot errors
                click.echo(f"[mei-annotate] read-architecture plot generation failed: {exc}")
    click.echo(f"[mei-annotate] wrote {len(candidate)} rows to {out_path}")
    click.echo(f"[mei-annotate] total annotate walltime={time.monotonic() - total_t0:.1f}s")
    return out_path
