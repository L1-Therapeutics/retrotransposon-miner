from __future__ import annotations

import gzip
import hashlib
import json
import math
import random
import re
import subprocess
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam
from intervaltree import IntervalTree

from retro_miner.igv_plots import generate_gold_review_igv_plots
from retro_miner.local_assembly import annotate_silver_with_local_assembly


@dataclass
class ClipAlignmentSummary:
    sample: str
    clip_count: int
    paf_hits: int


_MIN_MEI_ANCHOR_BP = 25
_MIN_POLYA_RUN_FOR_END_IMPUTE = 12
_MIN_MEI_ANCHOR_BP_RELAXED = 15
_MIN_REPORTABLE_MEI_SPAN_BP = 20


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


def _resolve_full_consensus_fasta(
    mei_fasta: Path,
    out_dir: Path,
    mei_full_fasta: Path | None = None,
) -> Path | None:
    """Resolve full-consensus FASTA for coordinate-normalized remapping.

    If `mei_full_fasta` is provided, it is used directly. Otherwise, build a
    small 3-record FASTA (ALU/SVA/LINE1) from the current panel by selecting
    canonical representatives.
    """
    if mei_full_fasta is not None:
        full = Path(mei_full_fasta)
        return full if full.exists() else None
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
        for cand in preferred:
            if cand in names:
                return cand
        fam = fam.sort_values("length", ascending=False)
        return str(fam.iloc[0]["name"])

    picks = {
        "ALU": _pick_target("ALU", ["AluY#SINE/Alu", "AluYa5#SINE/Alu", "AluYb8#SINE/Alu"]),
        "SVA": _pick_target("SVA", ["SVA_F#Retroposon/SVA", "SVA_D#Retroposon/SVA"]),
        "LINE1": _pick_target(
            "LINE1",
            [
                "L1HS#LINE/L1",
                "L1HS_full#LINE/L1",
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
            raw_score = (0.45 * pid) + (0.35 * qcov) + (0.2 * (mapq / 60.0))
            score = max(0.0, min(1.0, raw_score))
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
                    "mei_score": score,
                    "family": _family_from_target(tname),
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
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
        )

    paf = pd.DataFrame(rows)
    paf = paf.sort_values(
        ["qname", "mei_score", "alnlen", "mapq", "pid", "qcov"],
        ascending=[True, False, False, False, False, False],
    )
    return paf.drop_duplicates(subset=["qname"], keep="first")


def _canonicalize_alignment_hit_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize minimap2 hit columns after merges with possible name collisions."""
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
) -> tuple[pd.DataFrame, ClipAlignmentSummary]:
    if "clip_seq" not in split_df.columns:
        raise ValueError(
            "split evidence table is missing 'clip_seq'. Re-run extract-split-evidence with the current code."
        )

    clips = split_df.copy()
    clips = clips.loc[clips["clip_seq"].fillna("").astype(str).str.len() > 0].copy()
    clips["clip_id"] = [f"{sample}_{i}" for i in range(len(clips))]
    def _trim_clip_for_coord_alignment(seq: str, side: str) -> str:
        s = (seq or "").upper()
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
        preferred = right if side_u == "L" else left if side_u == "R" else (right if len(right) >= len(left) else left)
        backup = left if preferred is right else right
        if len(preferred) >= 12:
            return preferred
        if len(backup) >= 12:
            return backup
        joined = (left + right).strip()
        return joined if len(joined) >= 12 else s
    clips["clip_seq_coord"] = [
        _trim_clip_for_coord_alignment(str(seq), str(side))
        for seq, side in zip(clips["clip_seq"].fillna(""), clips["clip_side"].fillna(""))
    ]
    if clips.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return clips, summary

    with tempfile.TemporaryDirectory(prefix=f"rtm_mei_{sample}_") as tmpdir:
        tmp = Path(tmpdir)
        query_fa = tmp / "clips.fa"
        paf_path = tmp / "clips_vs_mei.paf"
        query_coord_fa = tmp / "clips_coord.fa"
        paf_coord_path = tmp / "clips_coord_vs_mei.paf"

        with query_fa.open("w", encoding="utf-8") as handle:
            for row in clips.itertuples(index=False):
                handle.write(f">{row.clip_id}\n{row.clip_seq}\n")

        cmd = [
            "minimap2",
            "-x",
            "sr",
            "--secondary=yes",
            "-c",
            str(mei_fasta),
            str(query_fa),
        ]
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        paf_path.write_text(proc.stdout, encoding="utf-8")

        best_hits = _best_hits_from_paf(paf_path)
        out = clips.merge(best_hits, left_on="clip_id", right_on="qname", how="left")
        out = _canonicalize_alignment_hit_columns(out)
        out["mei_hit"] = out["target"].notna()
        out["family"] = out["family"].fillna("")
        out["target"] = out["target"].fillna("")
        out["target_strand"] = out["target_strand"].fillna("")
        out["target_start"] = out["target_start"].fillna(0).astype(int)
        out["target_end"] = out["target_end"].fillna(0).astype(int)
        out["target_len"] = out["target_len"].fillna(0).astype(int)
        out["alnlen"] = out["alnlen"].fillna(0).astype(int)
        out["mei_score"] = out["mei_score"].fillna(0.0).astype(float)

        # Coordinate-only alignment path from poly-trimmed clips. This keeps raw
        # clip alignment for support scoring but improves MEI-axis coordinate recovery.
        with query_coord_fa.open("w", encoding="utf-8") as handle:
            for row in clips.itertuples(index=False):
                coord_seq = str(getattr(row, "clip_seq_coord", "") or "")
                if len(coord_seq) < 8:
                    continue
                handle.write(f">{row.clip_id}\n{coord_seq}\n")
        coord_hits = pd.DataFrame(columns=["qname", "target", "family", "target_strand", "target_start", "target_end", "target_len", "alnlen", "mapq", "pid", "qcov", "mei_score"])
        if query_coord_fa.exists() and query_coord_fa.stat().st_size > 0:
            proc_coord = subprocess.run(cmd[:-1] + [str(query_coord_fa)], check=True, capture_output=True, text=True)
            paf_coord_path.write_text(proc_coord.stdout, encoding="utf-8")
            coord_hits = _best_hits_from_paf(paf_coord_path)
        coord_hits = coord_hits.rename(
            columns={
                "target": "target_coord",
                "family": "family_coord",
                "target_strand": "target_strand_coord",
                "target_start": "target_start_coord",
                "target_end": "target_end_coord",
                "target_len": "target_len_coord",
                "alnlen": "alnlen_coord",
                "mapq": "mapq_coord",
                "pid": "pid_coord",
                "qcov": "qcov_coord",
                "mei_score": "mei_score_coord",
            }
        )
        out = out.merge(coord_hits, left_on="clip_id", right_on="qname", how="left", suffixes=("", "_coord_merge"))
        out["mei_hit_coord"] = out["target_coord"].notna()
        out["target_coord"] = out["target_coord"].fillna("")
        out["target_strand_coord"] = out["target_strand_coord"].fillna("")
        out["target_start_coord"] = out["target_start_coord"].fillna(0).astype(int)
        out["target_end_coord"] = out["target_end_coord"].fillna(0).astype(int)
        out["target_len_coord"] = out["target_len_coord"].fillna(0).astype(int)
        out["alnlen_coord"] = out["alnlen_coord"].fillna(0).astype(int)
        out["mapq_coord"] = out["mapq_coord"].fillna(0).astype(int)
        out["pid_coord"] = out["pid_coord"].fillna(0.0).astype(float)
        out["qcov_coord"] = out["qcov_coord"].fillna(0.0).astype(float)
        out["family_coord"] = out["family_coord"].fillna("")
        out["mei_score_coord"] = out["mei_score_coord"].fillna(0.0).astype(float)

    summary = ClipAlignmentSummary(
        sample=sample,
        clip_count=len(clips),
        paf_hits=int(out["mei_hit"].sum()),
    )
    return out, summary


def _align_discordant_reads_with_minimap2(
    discordant_df: pd.DataFrame,
    mei_fasta: Path,
    sample: str,
) -> tuple[pd.DataFrame, ClipAlignmentSummary]:
    if "read_seq" not in discordant_df.columns:
        raise ValueError(
            "discordant evidence table is missing 'read_seq'. Re-run extract-split-evidence with the current code."
        )

    reads = discordant_df.copy()
    reason = reads["discordant_reasons"].fillna("").astype(str)
    reads = reads.loc[
        (reads["read_seq"].fillna("").astype(str).str.len() >= 30)
        & (
            reason.str.contains("interchrom", regex=False)
            | reason.str.contains("mate_unmapped", regex=False)
            | reason.str.contains("large_insert", regex=False)
        )
    ].copy()
    reads["discordant_id"] = [f"{sample}_disc_{i}" for i in range(len(reads))]
    if reads.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return reads, summary

    with tempfile.TemporaryDirectory(prefix=f"rtm_mei_disc_{sample}_") as tmpdir:
        tmp = Path(tmpdir)
        query_fa = tmp / "discordant_reads.fa"
        paf_path = tmp / "discordant_vs_mei.paf"

        with query_fa.open("w", encoding="utf-8") as handle:
            for row in reads.itertuples(index=False):
                handle.write(f">{row.discordant_id}\n{row.read_seq}\n")

        cmd = [
            "minimap2",
            "-x",
            "sr",
            "--secondary=yes",
            "-c",
            str(mei_fasta),
            str(query_fa),
        ]
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        paf_path.write_text(proc.stdout, encoding="utf-8")
        best_hits = _best_hits_from_paf(paf_path)
        out = reads.merge(best_hits, left_on="discordant_id", right_on="qname", how="left")
        out = _canonicalize_alignment_hit_columns(out)
        out["mei_hit"] = out["target"].notna()
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
    return (int(best), float(at_bases) / float(len(s)), best_base)


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
    - support strings also include MEI_MAPPED, the count of split+discordant
      reads that map to MEI for that sample/locus.
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
        if informative_reads is not None and not informative_reads.empty:
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
        if informative_reads is not None and not informative_reads.empty:
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
        locus_gate["control_has_mei_support"].fillna(False).astype(bool)
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

        sr_raw = _counts_from_split(split_support_df, prefix)
        dpe_raw = _counts_from_discordant(disc_df, prefix)
        sr_strict = _counts_from_split(split_support_df, prefix, informative_reads=strict_reads)
        dpe_strict = _counts_from_discordant(disc_df, prefix, informative_reads=strict_reads)
        mei_mapped = _mei_mapped_counts(split_mei_df, disc_mei_df, prefix)

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
            .merge(locus_gate.loc[:, key_cols + ["locus_has_mei_support"]], on=key_cols, how="left")
        )

        sample_has_mei_col = f"{prefix}_has_mei_support"
        other_has_mei_col = f"{other_prefix}_has_mei_support"
        status_tbl = locus_gate.loc[:, key_cols].copy()
        status_tbl["sample_has_mei_support"] = (
            locus_gate[sample_has_mei_col].fillna(False).astype(bool)
            if sample_has_mei_col in locus_gate.columns
            else pd.Series(False, index=locus_gate.index)
        )
        status_tbl["other_has_mei_support"] = (
            locus_gate[other_has_mei_col].fillna(False).astype(bool)
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

        sr_l = sr_l_raw.where(~strict_mode, sr_l_strict).where(has_locus_mei_support, 0)
        sr_r = sr_r_raw.where(~strict_mode, sr_r_strict).where(has_locus_mei_support, 0)
        dpe_l = dpe_l_raw.where(~strict_mode, dpe_l_strict).where(has_locus_mei_support, 0)
        dpe_r = dpe_r_raw.where(~strict_mode, dpe_r_strict).where(has_locus_mei_support, 0)

        merged[f"{prefix}_supporting_reads"] = [
            f"SR_L={sl},SR_R={srx},DPE_L={dl},DPE_R={dr},MEI_MAPPED={mm}"
            for sl, srx, dl, dr, mm in zip(
                sr_l.tolist(),
                sr_r.tolist(),
                dpe_l.tolist(),
                dpe_r.tolist(),
                mei_mapped_total.tolist(),
            )
        ]
        if f"{prefix}_supporting_reads" in out.columns:
            out = out.drop(columns=[f"{prefix}_supporting_reads"])
        out = out.merge(merged[key_cols + [f"{prefix}_supporting_reads"]], on=key_cols, how="left")
        out[f"{prefix}_supporting_reads"] = (
            out[f"{prefix}_supporting_reads"].fillna("SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0").astype(str)
        )
    return out


def _aggregate_side_metrics(
    df: pd.DataFrame,
    sample_prefix: str,
    side: str,
    preferred_subfamily_by_locus: dict[tuple[str, int, int], str] | None = None,
) -> pd.DataFrame:
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
        return best

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
    coord_target_col = "target_coord" if "target_coord" in side_df.columns else "target"
    coord_mapq_col = "mapq_coord" if "mapq_coord" in side_df.columns else "mapq"
    coord_qcov_col = "qcov_coord" if "qcov_coord" in side_df.columns else "qcov"
    coord_pid_col = "pid_coord" if "pid_coord" in side_df.columns else "pid"
    coord_alnlen_col = "alnlen_coord" if "alnlen_coord" in side_df.columns else "alnlen"
    coord_start_col = "target_start_coord" if "target_start_coord" in side_df.columns else "target_start"
    coord_end_col = "target_end_coord" if "target_end_coord" in side_df.columns else "target_end"
    coord_tlen_col = "target_len_coord" if "target_len_coord" in side_df.columns else "target_len"
    if "mei_hit_coord" in side_df.columns:
        coord_hit_mask = side_df["mei_hit_coord"].fillna(False).astype(bool)
    else:
        coord_hit_mask = side_df["mei_hit"].fillna(False).astype(bool)

    family_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "family_effective"], as_index=False)["mei_score_effective"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score_effective"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"family_effective": f"{sample_prefix}_{side}_mei_family"})
    )
    subfamily_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "target_effective"], as_index=False)["mei_score_effective"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score_effective"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target_effective": f"{sample_prefix}_{side}_mei_subfamily"})
    )
    # Coordinate-estimation subset: avoid low-confidence/polyA-only hits that can
    # collapse inferred spans to tail-length artifacts.
    coord_mapq = pd.to_numeric(side_df[coord_mapq_col], errors="coerce").fillna(0)
    coord_qcov = pd.to_numeric(side_df[coord_qcov_col], errors="coerce").fillna(0.0)
    coord_pid = pd.to_numeric(side_df[coord_pid_col], errors="coerce").fillna(0.0)
    coord_alnlen = pd.to_numeric(side_df[coord_alnlen_col], errors="coerce").fillna(0)
    # MEI clips can be highly repetitive with low MAPQ despite very high identity/coverage.
    # Permit these as coordinate candidates if alignment quality itself is strong.
    strong_repetitive_clip = (coord_qcov >= 0.90) & (coord_pid >= 0.90) & (coord_alnlen >= 40)
    coord_df = side_df.loc[
        coord_hit_mask
        & ((coord_mapq >= 20) | strong_repetitive_clip)
        & (coord_qcov >= 0.60)
        & (coord_pid >= 0.85)
        & (coord_alnlen >= _MIN_MEI_ANCHOR_BP)
        & (side_df["poly_at_flag"] == 0)
    ].copy()
    if coord_df.empty:
        # Fallback for split-supported loci where strict filters are too harsh:
        # still require uniquely mappable-ish, non-polyA anchors.
        coord_df = side_df.loc[
            coord_hit_mask
            & ((coord_mapq >= 10) | strong_repetitive_clip)
            & (coord_qcov >= 0.35)
            & (coord_pid >= 0.75)
            & (coord_alnlen >= _MIN_MEI_ANCHOR_BP_RELAXED)
            & (side_df["poly_at_flag"] == 0)
        ].copy()
    if not coord_df.empty:
        coord_df["coord_target"] = coord_df[coord_target_col].fillna("").astype(str)
        coord_df["coord_target_start"] = pd.to_numeric(coord_df[coord_start_col], errors="coerce").fillna(0).astype(int)
        coord_df["coord_target_end"] = pd.to_numeric(coord_df[coord_end_col], errors="coerce").fillna(0).astype(int)
        coord_df["coord_target_len"] = pd.to_numeric(coord_df[coord_tlen_col], errors="coerce").fillna(0).astype(int)
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
        coord_df["coord_target"] = coord_df["coord_target"].fillna("").astype(str)
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
        coord_df = coord_df.loc[keep_mask].copy()
        if not coord_df.empty:
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
    mei_df = df.loc[df["mei_hit"]].copy()
    if mei_df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{sample_prefix}_discordant_mei_supported_reads",
                f"{sample_prefix}_discordant_mei_score_sum",
                f"{sample_prefix}_discordant_mei_family",
                f"{sample_prefix}_discordant_mei_subfamily",
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
        )

    mei_df["locus_midpoint"] = (mei_df["window_start"].astype(int) + mei_df["window_end"].astype(int)) // 2
    mei_df["anchor_side"] = mei_df.apply(
        lambda r: "L" if int(r["pos"]) <= int(r["locus_midpoint"]) else "R",
        axis=1,
    )
    mei_df["anchor_bin_10bp"] = (mei_df["pos"].astype(int) // 10).astype(int)

    family_top = (
        mei_df.groupby(["chrom", "window_start", "window_end", "family"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"family": f"{sample_prefix}_discordant_mei_family"})
    )
    subfamily_top = (
        mei_df.groupby(["chrom", "window_start", "window_end", "target"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target": f"{sample_prefix}_discordant_mei_subfamily"})
    )
    strand_top = (
        mei_df.groupby(["chrom", "window_start", "window_end", "target_strand"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target_strand": f"{sample_prefix}_discordant_mei_strand"})
    )
    family_sum = (
        mei_df.groupby(["chrom", "window_start", "window_end"], as_index=False)["mei_score"]
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

    side_subfamily_top = (
        mei_df.groupby(["chrom", "window_start", "window_end", "anchor_side", "target"], as_index=False)["mei_score"]
        .sum()
        .sort_values(
            ["chrom", "window_start", "window_end", "anchor_side", "mei_score"],
            ascending=[True, True, True, True, False],
        )
        .drop_duplicates(["chrom", "window_start", "window_end", "anchor_side"], keep="first")
        .rename(columns={"target": "side_top_subfamily"})
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
            )
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
    agg[f"{sample_prefix}_discordant_mei_left_subfamily"] = (
        agg[f"{sample_prefix}_discordant_mei_left_subfamily"].fillna("").astype(str)
    )
    agg[f"{sample_prefix}_discordant_mei_right_subfamily"] = (
        agg[f"{sample_prefix}_discordant_mei_right_subfamily"].fillna("").astype(str)
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
    reason_col = tmp["discordant_reasons"].fillna("").astype(str)
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
    tmp["mapq"] = tmp["mapq"].astype(float)

    side = (
        tmp.groupby(["chrom", "window_start", "window_end", "anchor_side"], as_index=False)
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
        l = int(row.get("disease_L_mei_breakpoint_mode", 0))
        r = int(row.get("disease_R_mei_breakpoint_mode", 0))
        if l > 0 and r > 0:
            return int((l + r) // 2), "disease_split"
        l = int(row.get("control_L_mei_breakpoint_mode", 0))
        r = int(row.get("control_R_mei_breakpoint_mode", 0))
        if l > 0 and r > 0:
            return int((l + r) // 2), "control_split"
        if l > 0:
            return l, "control_single"
        if r > 0:
            return r, "control_single"
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

    # Guardrail: pure polyA/polyT-only TSDs from rescue heuristics are often
    # reference-tail artifacts (especially in low-MEI-mapping contexts). Treat
    # these as non-TSD and fall back to non-TSD breakpoint evidence.
    tsd_src = out.get("tsd_evidence_source", "").fillna("").astype(str)
    tsd_seq_s = out.get("tsd_seq", "").fillna("").astype(str).str.upper()
    rescue_src = tsd_src.str.contains("_rescue", regex=False)
    tsd_len_s = tsd_seq_s.str.len().astype(int)
    a_fraction = (tsd_seq_s.str.count("A") / tsd_len_s.replace(0, pd.NA)).fillna(0.0).astype(float)
    t_fraction = (tsd_seq_s.str.count("T") / tsd_len_s.replace(0, pd.NA)).fillna(0.0).astype(float)
    dominant_poly_fraction = pd.concat([a_fraction, t_fraction], axis=1).max(axis=1)
    longest_at_run = tsd_seq_s.str.findall(r"[AT]+").map(
        lambda parts: max((len(p) for p in parts), default=0)
    ).astype(int)
    # Filter both strict polyA/T-only and near-pure A/T tails that contain only
    # a tiny number of non-A/T bases (common alignment jitter around homopolymers).
    poly_at_only_tsd = tsd_len_s.ge(4) & tsd_seq_s.str.fullmatch(r"[AT]+", na=False)
    # Keep rescue-derived TSDs conservative in low-complexity sequence:
    # a long pure/poly-A/T-like signature with >=6bp dominant run is often a
    # repeat-context artifact rather than a true short duplication footprint.
    near_poly_at_tsd = tsd_len_s.ge(8) & dominant_poly_fraction.ge(0.90) & longest_at_run.ge(6)
    poly_at_filter_mask = rescue_src & (poly_at_only_tsd | near_poly_at_tsd)
    out["tsd_poly_at_filter_applied"] = poly_at_filter_mask.astype(bool)
    if poly_at_filter_mask.any():
        out.loc[poly_at_filter_mask, "tsd_left_breakpoint"] = 0
        out.loc[poly_at_filter_mask, "tsd_right_breakpoint"] = 0
        out.loc[poly_at_filter_mask, "tsd_len_estimate"] = 0
        out.loc[poly_at_filter_mask, "tsd_seq"] = ""
        out.loc[poly_at_filter_mask, "tsd_detected"] = False
        out.loc[poly_at_filter_mask, "tsd_evidence_source"] = (
            out.loc[poly_at_filter_mask, "tsd_evidence_source"].astype(str) + "_filtered_poly_at_only"
        )
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
    out["insertion_breakpoint_pos"] = asm_bp.where(asm_has_mei & asm_bp.notna(), s("insertion_breakpoint_pos", 0)).fillna(0).astype(int)
    out.loc[asm_has_mei, "breakpoint_evidence_source"] = asm_source.loc[asm_has_mei]

    asm_tsd_seq = s("asm_tsd_seq", "").fillna("").astype(str)
    asm_tsd_len = pd.to_numeric(s("asm_tsd_len", float("nan")), errors="coerce")
    asm_tsd_detected = asm_has_mei & asm_tsd_len.notna() & asm_tsd_len.ge(4)
    out["tsd_seq"] = asm_tsd_seq.where(asm_tsd_detected & (asm_tsd_seq.str.len() > 0), s("tsd_seq", "").fillna("").astype(str))
    out["tsd_len_estimate"] = asm_tsd_len.where(asm_tsd_detected, s("tsd_len_estimate", 0)).fillna(0).astype(int)
    out["tsd_detected"] = out["tsd_len_estimate"].astype(float) >= 4.0

    asm_poly = pd.to_numeric(s("asm_polyA_max_run", float("nan")), errors="coerce")
    base_poly = pd.to_numeric(s("poly_at_max_run", 0), errors="coerce")
    picked_poly = asm_poly.where(asm_has_mei & asm_poly.notna(), base_poly)
    out["poly_at_max_run"] = (
        pd.concat([picked_poly, base_poly], axis=1).max(axis=1).fillna(0).astype(int)
    )

    asm_orient = s("asm_insertion_orientation", "").fillna("").astype(str)
    out["insertion_orientation"] = asm_orient.where(asm_has_mei & asm_orient.isin(["+", "-"]), s("insertion_orientation", "").fillna("").astype(str))

    asm_span = pd.to_numeric(s("asm_insertion_length", float("nan")), errors="coerce")
    out["insertion_mei_span"] = asm_span.where(asm_has_mei & asm_span.notna(), s("insertion_mei_span", 0)).fillna(0).astype(int)

    asm_start = pd.to_numeric(s("asm_insertion_mei_start", float("nan")), errors="coerce")
    asm_end = pd.to_numeric(s("asm_insertion_mei_end", float("nan")), errors="coerce")
    out["disease_insertion_mei_start"] = asm_start.where(
        asm_has_mei & (asm_source == "disease") & asm_start.gt(0),
        s("disease_insertion_mei_start", 0),
    ).fillna(0).astype(int)
    out["disease_insertion_mei_end"] = asm_end.where(
        asm_has_mei & (asm_source == "disease") & asm_end.gt(0),
        s("disease_insertion_mei_end", 0),
    ).fillna(0).astype(int)
    out["control_insertion_mei_start"] = asm_start.where(
        asm_has_mei & (asm_source == "control") & asm_start.gt(0),
        s("control_insertion_mei_start", 0),
    ).fillna(0).astype(int)
    out["control_insertion_mei_end"] = asm_end.where(
        asm_has_mei & (asm_source == "control") & asm_end.gt(0),
        s("control_insertion_mei_end", 0),
    ).fillna(0).astype(int)

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
    empty_family = out["mei_family"].fillna("").astype(str) == ""
    out.loc[empty_family, "mei_family"] = (
        out.loc[empty_family, "mei_subfamily"].fillna("").astype(str).map(_normalize_mei_family_token)
    )
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


def _complex_locus_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = [
        "discordant_disease_large_insert_fraction",
        "discordant_disease_interchrom_fraction",
        "discordant_disease_mate_unmapped_fraction",
        "discordant_disease_same_strand_fraction",
        "discordant_disease_improper_pair_fraction",
        "discordant_control_large_insert_fraction",
        "discordant_control_interchrom_fraction",
        "discordant_control_mate_unmapped_fraction",
        "discordant_control_same_strand_fraction",
        "discordant_control_improper_pair_fraction",
    ]
    parts = [_df_col_float(df, col) for col in fraction_cols if col in df.columns]
    if not parts:
        return pd.Series(0.0, index=df.index)
    return pd.concat(parts, axis=1).max(axis=1)


def _complex_locus_strong_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = [
        "discordant_disease_large_insert_fraction",
        "discordant_disease_interchrom_fraction",
        "discordant_disease_mate_unmapped_fraction",
        "discordant_control_large_insert_fraction",
        "discordant_control_interchrom_fraction",
        "discordant_control_mate_unmapped_fraction",
    ]
    parts = [_df_col_float(df, col) for col in fraction_cols if col in df.columns]
    if not parts:
        return pd.Series(0.0, index=df.index)
    return pd.concat(parts, axis=1).max(axis=1)


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
    complex_companion_fraction = _complex_locus_companion_fraction(out)
    complex_strong_companion_fraction = _complex_locus_strong_companion_fraction(out)
    large_insert_fraction = _df_col_float(out, "discordant_disease_large_insert_fraction")
    interchrom_fraction = _df_col_float(out, "discordant_disease_interchrom_fraction")
    mate_unmapped_fraction = _df_col_float(out, "discordant_disease_mate_unmapped_fraction")

    out["complex_sv_large_insert_flag"] = large_insert_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_interchrom_flag"] = interchrom_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_mate_unmapped_flag"] = mate_unmapped_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_companion_signal"] = complex_strong_companion_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
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
        out["complex_sv_mate_unmapped_flag"] & (out["complex_sv_signature_label"] != "none"),
        "complex_sv_signature_label",
    ] = out["complex_sv_signature_label"] + "+mate_unmapped"
    same_strand_fraction = _df_col_float(out, "discordant_disease_same_strand_fraction")
    improper_pair_fraction = _df_col_float(out, "discordant_disease_improper_pair_fraction")
    out.loc[
        (same_strand_fraction >= _COMPLEX_LOCUS_WEAK_MIN_FRACTION) & (out["complex_sv_signature_label"] == "none"),
        "complex_sv_signature_label",
    ] = "same_strand"
    out.loc[
        (improper_pair_fraction >= _COMPLEX_LOCUS_WEAK_MIN_FRACTION) & (out["complex_sv_signature_label"] == "none"),
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

    # Sample presence: genotype on per-locus SR/DPE support strings so a disease-
    # discovered locus with control-side SR/DPE evidence is marked shared.
    def _support_total_from_string_col(col_name: str) -> pd.Series:
        text = s(col_name, "").fillna("").astype(str)
        total = pd.Series(0, index=out.index, dtype="int64")
        for label in ("SR_L", "SR_R", "DPE_L", "DPE_R"):
            vals = pd.to_numeric(text.str.extract(rf"{label}=([0-9]+)", expand=False), errors="coerce").fillna(0).astype(int)
            total = total + vals
        return total.astype(float)

    disease_window_support_reads = _support_total_from_string_col("disease_supporting_reads")
    control_window_support_reads = _support_total_from_string_col("control_supporting_reads")
    disease_present = disease_window_support_reads.ge(1.0)
    control_present = control_window_support_reads.ge(1.0)

    out["sample_status_label"] = "low_support"
    out.loc[disease_present & control_present, "sample_status_label"] = "shared"
    out.loc[disease_present & (~control_present), "sample_status_label"] = "disease_only"
    out.loc[(~disease_present) & control_present, "sample_status_label"] = "control_only"

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
    print(
        f"[mei-annotate] empirical stage: bedtools shuffle start "
        f"scope={scope} n={n_windows} chroms={len(target_chroms)} span={sampled_span}",
        flush=True,
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
        print(
            f"[mei-annotate] empirical stage: bedtools shuffle done windows={len(rows)}",
            flush=True,
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
        print("[mei-annotate] empirical stage: trying bedtools-based random sampling", flush=True)
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
        print("[mei-annotate] empirical stage: bedtools sampling unavailable; using python fallback", flush=True)
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

    print(
        f"[mei-annotate] empirical stage: python fallback sampling done windows={len(windows)} "
        f"attempts={attempts}",
        flush=True,
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
            print(f"[mei-annotate] empirical scoring {sample_prefix}: {i}/{total} loci", flush=True)
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
        print("[mei-annotate] empirical stage skipped: no loci passed empirical prefilter", flush=True)
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
                print(
                    f"[mei-annotate] empirical cache hit key={cache_key} "
                    f"rows={len(random_disease_df)}",
                    flush=True,
                )
            except Exception:
                cache_hit = False
        else:
            print(f"[mei-annotate] empirical cache miss key={cache_key}", flush=True)

    prep_t0 = time.monotonic()
    print("[mei-annotate] empirical stage: preparing junk exclusion masks", flush=True)
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
        print(
            f"[mei-annotate] empirical stage: junk masks ready chroms={len(junk_trees)} intervals={junk_interval_count}",
            flush=True,
        )
    elif not cache_hit and merged_exclusion_ready:
        print(
            f"[mei-annotate] empirical stage: using merged exclusion bed {empirical_exclude_merged_bed}",
            flush=True,
        )
    print(
        f"[mei-annotate] empirical stage: exclusion mask prep elapsed={time.monotonic() - prep_t0:.1f}s",
        flush=True,
    )

    with pysam.AlignmentFile(str(disease_bam_path), "rb") as disease_bam, pysam.AlignmentFile(
        str(control_bam_path), "rb"
    ) as control_bam:
        total_loci = int(len(idxs))
        loci_progress_every = 100
        print(f"[mei-annotate] empirical stage: computing context metrics for {total_loci} loci", flush=True)
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
                print(
                    f"[mei-annotate] empirical stage: locus metrics {i}/{total_loci} "
                    f"(elapsed={elapsed:.1f}s)",
                    flush=True,
                )

        if not cache_hit:
            print("[mei-annotate] empirical stage: building random-window background metrics", flush=True)
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
            print(
                f"[mei-annotate] empirical stage: sampled {len(random_windows)} random windows "
                f"(scope={empirical_random_scope}, n={empirical_random_windows})",
                flush=True,
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
                    print(
                        f"[mei-annotate] empirical stage: random-window metrics {i}/{total_random} "
                        f"(elapsed={elapsed:.1f}s)",
                        flush=True,
                    )
            random_disease_df = pd.DataFrame(random_disease_rows)
            random_control_df = pd.DataFrame(random_control_rows)
            if empirical_cache_dir is not None:
                disease_cache_path = empirical_cache_dir / f"{cache_key}.disease.parquet"
                control_cache_path = empirical_cache_dir / f"{cache_key}.control.parquet"
                random_disease_df.to_parquet(disease_cache_path, index=False)
                random_control_df.to_parquet(control_cache_path, index=False)
                print(
                    f"[mei-annotate] empirical cache write key={cache_key} "
                    f"rows={len(random_disease_df)}",
                    flush=True,
                )
        else:
            print("[mei-annotate] empirical stage: using cached random-window metrics", flush=True)

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
    print("[mei-annotate] empirical stage: applying empirical p-value scoring complete", flush=True)
    elapsed_total = time.monotonic() - t0
    print(f"[mei-annotate] empirical stage complete (elapsed={elapsed_total:.1f}s)", flush=True)
    print(f"[mei-annotate] empirical stage walltime={time.monotonic() - stage_start:.1f}s", flush=True)
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
    print(
        "[mei-annotate] stage counts "
        f"bronze={len(out)} silver={int(out['silver_stage_pass'].sum())}",
        flush=True,
    )
    return out


def _assign_gold_stage(
    candidates: pd.DataFrame,
    empirical_p_threshold: float = 0.001,
    empirical_stage: bool = False,
) -> pd.DataFrame:
    out = candidates.copy()
    out["gold_empirical_p_threshold"] = float(empirical_p_threshold)
    out["gold_empirical_eval_available"] = False
    out["gold_empirical_outlier"] = False
    out["gold_stage_pass"] = False
    out["gold_stage_fail_reason"] = ""

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
    print(
        "[mei-annotate] stage counts "
        f"silver={int(silver.sum())} gold={int(out['gold_stage_pass'].sum())}",
        flush=True,
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
    support_total = (out["_prio_split_reads_max"] + out["_prio_discordant_reads_max"]).astype(float)
    mei_mapped_frac = (out["_prio_mei_mapped_max"].astype(float) / support_total.where(support_total > 0, 1.0)).clip(
        lower=0.0, upper=1.0
    )
    out["read_support_heuristic_score"] = (
        0.50 * (out["_prio_split_reads_max"].astype(float).map(math.log1p))
        + 0.35 * (out["_prio_discordant_reads_max"].astype(float).map(math.log1p))
        + 0.15 * mei_mapped_frac
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
            "read_support_heuristic_score",
            "_prio_mei_mapped_max",
            "_prio_split_reads_max",
            "_prio_discordant_reads_max",
            "_prio_tsd",
            "_prio_poly_at",
            "_prio_breakpoint_consensus_available",
            "_prio_clip_overlap_consistency",
            "_prio_insertion_model_score",
        ]
    )
    # Artifact flag sorted ascending (False first), then strongest support.
    ascending.extend([True] + [False] * 9)
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
    bp_candidates = pd.concat(
        [
            tsd_left,
            tsd_right,
            pd.to_numeric(s("disease_L_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("disease_R_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_L_mei_breakpoint_mode", float("nan")), errors="coerce"),
            pd.to_numeric(s("control_R_mei_breakpoint_mode", float("nan")), errors="coerce"),
        ],
        axis=1,
    )
    bp_candidates = bp_candidates.where(bp_candidates.gt(0))
    bp_lo = bp_candidates.min(axis=1, skipna=True)
    bp_hi = bp_candidates.max(axis=1, skipna=True)

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


def _build_gold_review_table(candidates: pd.DataFrame, empirical_stage: bool = False) -> pd.DataFrame:
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

    # Assembly-preferred consensus fields (non-destructive): use assembly-derived
    # values when present, else fall back to current evidence-derived fields.
    asm_bp = pd.to_numeric(_series_or_default("asm_consensus_breakpoint_pos", float("nan")), errors="coerce")
    out["consensus_insertion_breakpoint_pos"] = asm_bp.where(asm_has_mei & asm_bp.notna(), out["insertion_breakpoint_pos"]).astype(int)
    out["consensus_breakpoint_source"] = asm_source.copy()
    out.loc[out["consensus_breakpoint_source"] == "", "consensus_breakpoint_source"] = _series_or_default(
        "breakpoint_evidence_source", ""
    ).fillna("").astype(str)
    out = _derive_breakpoint_interval_fields(
        out,
        breakpoint_pos_col="consensus_insertion_breakpoint_pos",
        output_prefix="consensus_",
    )
    out.loc[
        (out["consensus_breakpoint_interval_width_bp"].astype(int) > 0) & (~asm_has_mei),
        "consensus_breakpoint_source",
    ] = "interval_midpoint"

    asm_tsd_seq = _series_or_default("asm_tsd_seq", "").fillna("").astype(str)
    asm_tsd_len = pd.to_numeric(_series_or_default("asm_tsd_len", float("nan")), errors="coerce")
    asm_tsd_detected = asm_has_mei & asm_tsd_len.notna() & asm_tsd_len.ge(4)
    out["consensus_tsd_seq"] = asm_tsd_seq.where(
        asm_tsd_detected & (asm_tsd_seq.str.len() > 0),
        _series_or_default("tsd_seq", "").fillna("").astype(str),
    )
    base_tsd_len = pd.to_numeric(_series_or_default("tsd_len_estimate", float("nan")), errors="coerce")
    out["consensus_tsd_len_estimate"] = asm_tsd_len.where(asm_tsd_detected, base_tsd_len)
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
    out["consensus_poly_at_max_run"] = (
        pd.concat([picked_poly, base_poly], axis=1).max(axis=1).fillna(0).round().astype(int)
    )
    out["consensus_poly_at_supported"] = out["consensus_poly_at_max_run"].fillna(0).astype(float) >= 8.0

    asm_span = pd.to_numeric(_series_or_default("asm_insertion_length", float("nan")), errors="coerce")
    base_span = pd.to_numeric(_series_or_default("insertion_mei_span", float("nan")), errors="coerce")
    asm_mei_start = pd.to_numeric(_series_or_default("asm_insertion_mei_start", float("nan")), errors="coerce")
    asm_mei_end = pd.to_numeric(_series_or_default("asm_insertion_mei_end", float("nan")), errors="coerce")
    disease_start = pd.to_numeric(_series_or_default("disease_insertion_mei_start", float("nan")), errors="coerce")
    control_start = pd.to_numeric(_series_or_default("control_insertion_mei_start", float("nan")), errors="coerce")
    disease_end = pd.to_numeric(_series_or_default("disease_insertion_mei_end", float("nan")), errors="coerce")
    control_end = pd.to_numeric(_series_or_default("control_insertion_mei_end", float("nan")), errors="coerce")
    asm_pair_valid = asm_has_mei & asm_mei_start.gt(0) & asm_mei_end.gt(0)
    disease_pair_valid = disease_start.gt(0) & disease_end.gt(0)
    control_pair_valid = control_start.gt(0) & control_end.gt(0)
    raw_start = pd.Series([float("nan")] * len(out), index=out.index)
    raw_end = pd.Series([float("nan")] * len(out), index=out.index)
    raw_start = raw_start.where(~asm_pair_valid, asm_mei_start)
    raw_end = raw_end.where(~asm_pair_valid, asm_mei_end)
    disease_pick = (~asm_pair_valid) & disease_pair_valid
    raw_start = raw_start.where(~disease_pick, disease_start)
    raw_end = raw_end.where(~disease_pick, disease_end)
    control_pick = (~asm_pair_valid) & (~disease_pair_valid) & control_pair_valid
    raw_start = raw_start.where(~control_pick, control_start)
    raw_end = raw_end.where(~control_pick, control_end)

    asm_orient = _series_or_default("asm_insertion_orientation", "").fillna("").astype(str)
    out["consensus_insertion_orientation"] = asm_orient.where(
        asm_orient.isin(["+", "-"]),
        _series_or_default("insertion_orientation", "").fillna("").astype(str),
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
        asm_span.where(asm_has_mei & asm_span.notna(), base_span),
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
        target_len_hint = pd.concat(
            [
                pd.to_numeric(_series_or_default("asm_mei_target_length", float("nan")), errors="coerce"),
                pd.to_numeric(_series_or_default("disease_L_mei_target_len", float("nan")), errors="coerce"),
                pd.to_numeric(_series_or_default("disease_R_mei_target_len", float("nan")), errors="coerce"),
                pd.to_numeric(_series_or_default("control_L_mei_target_len", float("nan")), errors="coerce"),
                pd.to_numeric(_series_or_default("control_R_mei_target_len", float("nan")), errors="coerce"),
                pd.to_numeric(base_span, errors="coerce"),
            ],
            axis=1,
        ).where(lambda x: x.gt(0)).max(axis=1, skipna=True)
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
        poly_run_any = pd.to_numeric(_series_or_default("consensus_poly_at_max_run", 0), errors="coerce").fillna(0.0)
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
    # Consolidate legacy consensus coordinates into *_full when full remap is
    # unresolved. For ALU, also carry over one-sided coordinate fields.
    base_span_num = pd.to_numeric(out["consensus_insertion_mei_span"], errors="coerce").fillna(0.0)
    base_5p_num = pd.to_numeric(out["consensus_insertion_mei_5p_coord"], errors="coerce").fillna(-1.0)
    base_3p_num = pd.to_numeric(out["consensus_insertion_mei_3p_coord"], errors="coerce").fillna(-1.0)
    full_span_num = pd.to_numeric(out["consensus_insertion_mei_span_full"], errors="coerce").fillna(0.0)
    full_5p_num = pd.to_numeric(out["consensus_insertion_mei_5p_coord_full"], errors="coerce").fillna(-1.0)
    full_3p_num = pd.to_numeric(out["consensus_insertion_mei_3p_coord_full"], errors="coerce").fillna(-1.0)
    full_span_missing = full_span_num.le(0) & base_span_num.gt(0)
    if full_span_missing.any():
        out.loc[full_span_missing, "consensus_insertion_mei_span_full"] = (
            base_span_num.loc[full_span_missing].round().astype(int)
        )
    fam_upper = out["consensus_mei_family"].fillna("").astype(str).str.upper()
    alu_mask = fam_upper.eq("ALU")
    fill_full_5p = alu_mask & full_5p_num.le(0) & base_5p_num.gt(0)
    fill_full_3p = alu_mask & full_3p_num.le(0) & base_3p_num.gt(0)
    if fill_full_5p.any():
        out.loc[fill_full_5p, "consensus_insertion_mei_5p_coord_full"] = base_5p_num.loc[fill_full_5p]
    if fill_full_3p.any():
        out.loc[fill_full_3p, "consensus_insertion_mei_3p_coord_full"] = base_3p_num.loc[fill_full_3p]
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
        return pd.Series(
            [
                f"SR_L={sl},SR_R={sr},DPE_L={dl},DPE_R={dr},MEI_MAPPED={mm}"
                for sl, sr, dl, dr, mm in zip(
                    sr_l.tolist(),
                    sr_r.tolist(),
                    dpe_l.tolist(),
                    dpe_r.tolist(),
                    mei_mapped.tolist(),
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
        base = text.str.replace(r",?MEI_MAPPED=[0-9]+", "", regex=True).str.strip(",")
        base = base.where(base.str.len() > 0, "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0")
        return base + ",MEI_MAPPED=" + mapped.astype(str)

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
        "control_supporting_reads",
        "disease_supporting_reads",
        "sample_status_label",
        "consensus_tsd_seq",
        "consensus_poly_at_max_run",
        "consensus_mei_family",
        "consensus_mei_subfamily",
        "known_mei_polymorphism_family",
        "known_mei_polymorphism_id",
        "known_mei_polymorphism_source",
        "complex_mei_event",
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
        "complex_mei_event",
        "asm_complex_class",
        "chrom",
        "window_start",
        "window_end",
        "consensus_insertion_breakpoint_pos",
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
        "consensus_poly_at_max_run",
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
        print(f"[mei-annotate] parsed g1k MEI VCF records kept={kept} path={g1k_mei_vcf}")

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
        proc = subprocess.run(intersect_cmd, check=True, capture_output=True, text=True)
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
        print(f"[mei-annotate] parsed long-read SVAN MEI VCF records kept={kept} path={lr_mei_vcf}")

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
        proc = subprocess.run(intersect_cmd, check=True, capture_output=True, text=True)
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


def _choose_event_subfamily(row: pd.Series) -> str:
    weighted: list[tuple[str, int]] = []
    for subfamily_col, weight_col in [
        ("disease_discordant_mei_subfamily", "disease_discordant_mei_supported_reads"),
        ("disease_discordant_mei_left_subfamily", "disease_discordant_mei_left_supported_reads"),
        ("disease_discordant_mei_right_subfamily", "disease_discordant_mei_right_supported_reads"),
        ("disease_L_mei_subfamily", "disease_L_mei_supported_reads"),
        ("disease_R_mei_subfamily", "disease_R_mei_supported_reads"),
        ("control_discordant_mei_subfamily", "control_discordant_mei_supported_reads"),
        ("control_discordant_mei_left_subfamily", "control_discordant_mei_left_supported_reads"),
        ("control_discordant_mei_right_subfamily", "control_discordant_mei_right_supported_reads"),
        ("control_L_mei_subfamily", "control_L_mei_supported_reads"),
        ("control_R_mei_subfamily", "control_R_mei_supported_reads"),
    ]:
        label = str(row.get(subfamily_col, "") or "").strip()
        weight = _row_int(row, weight_col)
        if label and weight > 0:
            weighted.append((label, weight))
    g1k_subfamily = str(row.get("g1k_melt_insertion_subfamily", "") or "").strip()
    if g1k_subfamily:
        weighted.append((g1k_subfamily, 1))
    lr_subfamily = str(row.get("lr_svan_subfamily", "") or "").strip()
    if lr_subfamily:
        weighted.append((lr_subfamily, 1))
    known_subfamily = str(row.get("known_mei_polymorphism_subfamily", "") or "").strip()
    if known_subfamily and known_subfamily not in {"MULTI_SOURCE"}:
        weighted.append((known_subfamily, 1))
    if not weighted:
        return ""
    weighted.sort(key=lambda item: item[1], reverse=True)
    return weighted[0][0]


def _choose_event_family(row: pd.Series) -> str:
    candidates = [
        str(row.get("disease_discordant_mei_family", "")),
        str(row.get("disease_L_mei_family", "")),
        str(row.get("disease_R_mei_family", "")),
        str(row.get("control_discordant_mei_family", "")),
        str(row.get("control_L_mei_family", "")),
        str(row.get("control_R_mei_family", "")),
        str(row.get("disease_discordant_mei_subfamily", "")),
        str(row.get("disease_L_mei_subfamily", "")),
        str(row.get("disease_R_mei_subfamily", "")),
        str(row.get("control_discordant_mei_subfamily", "")),
        str(row.get("control_L_mei_subfamily", "")),
        str(row.get("control_R_mei_subfamily", "")),
        str(row.get("g1k_melt_insertion_subfamily", "")),
        str(row.get("g1k_melt_id", "")),
        str(row.get("lr_svan_subfamily", "")),
        str(row.get("lr_svan_id", "")),
        str(row.get("known_mei_polymorphism_family", "")),
        str(row.get("known_mei_polymorphism_subfamily", "")),
    ]
    for c in candidates:
        fam = _normalize_mei_family_token(c)
        if fam:
            return fam
    return ""


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


def _load_rmsk_interval_trees(rmsk_table_path: Path) -> tuple[dict[str, IntervalTree], int]:
    trees: dict[str, IntervalTree] = {}
    intervals_with_family = 0
    opener = gzip.open if str(rmsk_table_path).endswith(".gz") else open
    with opener(rmsk_table_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = ""
            start0 = -1
            end0 = -1
            strand = ""
            rep_name = ""
            rep_class = ""
            rep_family = ""
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
            if end0 <= start0:
                continue
            tree = trees.setdefault(chrom, IntervalTree())
            tree.addi(
                int(start0),
                int(end0),
                {
                    "rep_name": rep_name,
                    "rep_class": rep_class,
                    "rep_family": rep_family,
                    "strand": strand,
                },
            )
            if _rmsk_interval_family_norm(
                {
                    "rep_name": rep_name,
                    "rep_class": rep_class,
                    "rep_family": rep_family,
                }
            ):
                intervals_with_family += 1
    return trees, intervals_with_family


def _rmsk_interval_family_norm(data: dict[str, object]) -> str:
    return _normalize_mei_family_token(
        f"{data.get('rep_name', '')} {data.get('rep_class', '')} {data.get('rep_family', '')}"
    )


def _annotate_nested_retrotransposon(candidates: pd.DataFrame, rmsk_table_path: Path) -> pd.DataFrame:
    out = candidates.copy()
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

    trees, intervals_with_family = _load_rmsk_interval_trees(rmsk_table_path)
    if not trees or intervals_with_family == 0:
        raise ValueError(
            "RepeatMasker annotation lacks repName/repClass/repFamily fields required "
            "for nested_same_class_orientation. Use a UCSC rmsk table (e.g. rmsk.txt.gz), "
            "not a stripped 3-4 column BED."
        )
    selected: list[dict[str, object]] = []
    for row in out.itertuples(index=False):
        as_row = pd.Series(row._asdict())
        chrom = str(getattr(row, "chrom"))
        pos_1based = int(getattr(row, "insertion_breakpoint_pos", 0))
        if pos_1based <= 0:
            pos_1based = int((int(getattr(row, "window_start", 1)) + int(getattr(row, "window_end", 1))) // 2)
        pos0 = max(0, pos_1based - 1)
        event_family = _choose_event_family(as_row)
        event_orientation = _choose_event_orientation(as_row)
        tree = trees.get(chrom)
        overlaps = list(tree.at(pos0)) if tree is not None else []
        if not overlaps or not event_family:
            selected.append({})
            continue

        same_family_overlaps = [
            iv for iv in overlaps if _rmsk_interval_family_norm(iv.data) == event_family
        ]
        if not same_family_overlaps:
            selected.append({})
            continue

        def score(iv) -> tuple[int, int]:
            d = iv.data
            rep_strand = str(d.get("strand", "")).strip()
            same_orient = int(
                bool(event_orientation)
                and rep_strand in {"+", "-"}
                and rep_strand == event_orientation
            )
            return (same_orient, int(iv.end - iv.begin))

        best = max(same_family_overlaps, key=score)
        d = best.data
        rep_strand = str(d.get("strand", "")).strip()
        same_orient = bool(event_orientation) and rep_strand in {"+", "-"} and (rep_strand == event_orientation)
        selected.append(
            {
                "nested_repeat_overlap": True,
                "nested_repeat_name": str(d.get("rep_name", "")),
                "nested_repeat_class": str(d.get("rep_class", "")),
                "nested_repeat_family": str(d.get("rep_family", "")),
                "nested_repeat_strand": rep_strand,
                "nested_mei_family": event_family,
                "nested_insertion_orientation": event_orientation,
                "nested_same_class": True,
                "nested_same_orientation": same_orient,
                "nested_same_class_orientation": "nested" if same_orient else "unnested",
            }
        )

    sel_df = pd.DataFrame(selected)
    if not sel_df.empty:
        for col in [
            "nested_repeat_overlap",
            "nested_repeat_name",
            "nested_repeat_class",
            "nested_repeat_family",
            "nested_repeat_strand",
            "nested_mei_family",
            "nested_insertion_orientation",
            "nested_same_class",
            "nested_same_orientation",
            "nested_same_class_orientation",
        ]:
            if col in sel_df.columns:
                out[col] = sel_df[col].where(sel_df[col].notna(), out[col])
    return out


def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    reference_fasta: Path | None = None,
    disease_bam_path: Path | None = None,
    control_bam_path: Path | None = None,
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
) -> Path:
    total_t0 = time.monotonic()
    candidate = pd.read_csv(candidate_loci_path, sep="\t")

    split_disease_raw = _load_table(evidence_dir, "split_evidence", "disease")
    split_control_raw = _load_table(evidence_dir, "split_evidence", "control")
    discordant_disease_raw = _load_table(evidence_dir, "discordant_evidence", "disease")
    discordant_control_raw = _load_table(evidence_dir, "discordant_evidence", "control")
    split_disease = _assign_rows_to_candidate_loci(split_disease_raw, candidate)
    split_control = _assign_rows_to_candidate_loci(split_control_raw, candidate)
    discordant_disease = _assign_rows_to_candidate_loci(discordant_disease_raw, candidate)
    discordant_control = _assign_rows_to_candidate_loci(discordant_control_raw, candidate)
    if disease_bam_path is not None:
        indel_disease = _collect_indel_breakpoint_evidence(
            disease_bam_path,
            candidate,
            sample="disease",
        )
    else:
        indel_disease = pd.DataFrame()
    if control_bam_path is not None:
        indel_control = _collect_indel_breakpoint_evidence(
            control_bam_path,
            candidate,
            sample="control",
        )
    else:
        indel_control = pd.DataFrame()
    disease_hits, disease_summary = _align_clips_with_minimap2(split_disease, mei_fasta, sample="disease")
    control_hits, control_summary = _align_clips_with_minimap2(split_control, mei_fasta, sample="control")
    disease_disc_hits, disease_disc_summary = _align_discordant_reads_with_minimap2(
        discordant_disease, mei_fasta, sample="disease"
    )
    control_disc_hits, control_disc_summary = _align_discordant_reads_with_minimap2(
        discordant_control, mei_fasta, sample="control"
    )
    full_consensus_fasta = _resolve_full_consensus_fasta(
        mei_fasta=mei_fasta,
        out_dir=out_path.parent,
        mei_full_fasta=mei_full_fasta,
    )
    if full_consensus_fasta is not None:
        fai_full = full_consensus_fasta.with_suffix(full_consensus_fasta.suffix + ".fai")
        if fai_full.exists():
            try:
                fai_tbl = pd.read_csv(fai_full, sep="\t", header=None, names=["name", "length", "offset", "line_bases", "line_width"])
                line1_rows = fai_tbl.loc[fai_tbl["name"].astype(str).str.contains("LINE1", na=False)]
                if not line1_rows.empty:
                    max_l1_len = pd.to_numeric(line1_rows["length"], errors="coerce").fillna(0).max()
                    if max_l1_len < 5000:
                        print(
                            "[mei-annotate] warning: full-consensus LINE1 target is shorter than 5kb "
                            f"(max={int(max_l1_len)}). Provide --mei-full-fasta with full-length L1HS for "
                            "6kb-relative coordinates."
                        )
            except Exception:
                pass
    if full_consensus_fasta is not None:
        disease_hits_full, disease_summary_full = _align_clips_with_minimap2(
            split_disease,
            full_consensus_fasta,
            sample="disease_full",
        )
        control_hits_full, control_summary_full = _align_clips_with_minimap2(
            split_control,
            full_consensus_fasta,
            sample="control_full",
        )
        disease_disc_hits_full, disease_disc_summary_full = _align_discordant_reads_with_minimap2(
            discordant_disease,
            full_consensus_fasta,
            sample="disease_full",
        )
        control_disc_hits_full, control_disc_summary_full = _align_discordant_reads_with_minimap2(
            discordant_control,
            full_consensus_fasta,
            sample="control_full",
        )
    else:
        disease_hits_full = pd.DataFrame()
        control_hits_full = pd.DataFrame()
        disease_disc_hits_full = pd.DataFrame()
        control_disc_hits_full = pd.DataFrame()
        disease_summary_full = ClipAlignmentSummary(sample="disease_full", clip_count=0, paf_hits=0)
        control_summary_full = ClipAlignmentSummary(sample="control_full", clip_count=0, paf_hits=0)
        disease_disc_summary_full = ClipAlignmentSummary(sample="disease_full", clip_count=0, paf_hits=0)
        control_disc_summary_full = ClipAlignmentSummary(sample="control_full", clip_count=0, paf_hits=0)

    print(
        f"[mei-annotate] disease clips={disease_summary.clip_count} hits={disease_summary.paf_hits}; "
        f"control clips={control_summary.clip_count} hits={control_summary.paf_hits}; "
        f"disease discordant reads={disease_disc_summary.clip_count} hits={disease_disc_summary.paf_hits}; "
        f"control discordant reads={control_disc_summary.clip_count} hits={control_disc_summary.paf_hits}; "
        f"disease full clips={disease_summary_full.clip_count} hits={disease_summary_full.paf_hits}; "
        f"control full clips={control_summary_full.clip_count} hits={control_summary_full.paf_hits}; "
        f"disease full discordant reads={disease_disc_summary_full.clip_count} hits={disease_disc_summary_full.paf_hits}; "
        f"control full discordant reads={control_disc_summary_full.clip_count} hits={control_disc_summary_full.paf_hits}; "
        f"disease indel reads={len(indel_disease)} control indel reads={len(indel_control)}"
    )

    def _mei_rows_only(df: pd.DataFrame) -> pd.DataFrame:
        if ("mei_hit" not in df.columns) and ("mei_hit_coord" not in df.columns):
            return df.iloc[0:0].copy()
        raw_hit = df["mei_hit"].fillna(False).astype(bool) if "mei_hit" in df.columns else pd.Series(False, index=df.index)
        coord_hit = (
            df["mei_hit_coord"].fillna(False).astype(bool)
            if "mei_hit_coord" in df.columns
            else pd.Series(False, index=df.index)
        )
        return df.loc[raw_hit | coord_hit].copy()

    split_disease_mei = _mei_rows_only(disease_hits)
    split_control_mei = _mei_rows_only(control_hits)
    discordant_disease_mei = _mei_rows_only(disease_disc_hits)
    discordant_control_mei = _mei_rows_only(control_disc_hits)

    disc_t = _aggregate_discordant_mei_metrics(disease_disc_hits, sample_prefix="disease")
    disc_n = _aggregate_discordant_mei_metrics(control_disc_hits, sample_prefix="control")
    disc_t_full = _aggregate_discordant_mei_metrics(disease_disc_hits_full, sample_prefix="disease_full")
    disc_n_full = _aggregate_discordant_mei_metrics(control_disc_hits_full, sample_prefix="control_full")
    disc_anchor_t = _aggregate_discordant_anchor_side_metrics(discordant_disease, sample_prefix="disease")
    disc_anchor_n = _aggregate_discordant_anchor_side_metrics(discordant_control, sample_prefix="control")

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

    for idx, part in enumerate(anno_parts):
        if part.empty:
            continue
        candidate = candidate.merge(part, on=["chrom", "window_start", "window_end"], how="left")
        if (idx + 1) % 2 == 0:
            print(f"[mei-annotate] merged side metrics {idx + 1}/{len(anno_parts)}")

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
    print("[mei-annotate] merged discordant MEI support metrics")

    candidate = _add_candidate_support_info_fields(
        candidate,
        split_disease=split_disease,
        split_control=split_control,
        discordant_disease=discordant_disease,
        discordant_control=discordant_control,
        split_disease_mei=split_disease_mei,
        split_control_mei=split_control_mei,
        discordant_disease_mei=discordant_disease_mei,
        discordant_control_mei=discordant_control_mei,
        indel_disease=indel_disease,
        indel_control=indel_control,
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
    candidate = _assign_bronze_silver_stages(candidate)
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
        print(
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
        print(
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
        print(
            f"[mei-annotate] added long-read SVAN polymorphism overlap fields "
            f"(elapsed={time.monotonic() - lr_t0:.1f}s)"
        )
    candidate = _add_known_mei_polymorphism_consensus(candidate)
    candidate = _add_consolidated_event_fields(candidate)
    candidate = _broaden_poly_at_fields(candidate)
    if rmsk_table_path is not None:
        rmsk_t0 = time.monotonic()
        candidate = _annotate_nested_retrotransposon(candidate, rmsk_table_path=rmsk_table_path)
        print(
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
        print(
            f"[mei-annotate] added BAM-depth controlization for family-consistent, junk-clean loci "
            f"(elapsed={time.monotonic() - emp_t0:.1f}s)"
        )
    elif not empirical_stage:
        print("[mei-annotate] empirical stage disabled (--no-empirical-stage)", flush=True)
    candidate = _add_heuristic_assembly_like_vaf_fields(candidate)
    candidate = _assign_gold_stage(candidate, empirical_stage=empirical_stage)

    candidate = _apply_breakpoint_motif_report_gating(candidate)
    candidate = _prioritize_mei_candidates(candidate, stage_first=True)

    candidate_tsv = _stable_tsv_export_frame(candidate)
    candidate_tsv.to_csv(out_path, sep="\t", index=False)
    candidate.to_parquet(out_path.with_suffix(".parquet"), index=False)
    gold_review = _build_gold_review_table(candidate, empirical_stage=empirical_stage)
    gold_review_path = out_path.with_name(out_path.stem + ".gold_review.tsv")
    gold_review_tsv = _stable_tsv_export_frame(gold_review)
    gold_review_tsv.to_csv(gold_review_path, sep="\t", index=False)
    print(f"[mei-annotate] wrote gold review table to {gold_review_path}")
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
            print(f"[mei-annotate] IGV snapshot generation skipped: {exc}", flush=True)
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            print(
                f"[mei-annotate] IGV snapshot generation failed: {exc}; "
                f"batch script remains at {igv_dir / 'igv_batch.txt'}",
                flush=True,
            )
    elif igv_plots and (disease_bam_path is None or control_bam_path is None or reference_fasta is None):
        print(
            "[mei-annotate] IGV snapshot generation skipped: require --reference-fasta, "
            "--disease-bam-depth, and --control-bam-depth",
            flush=True,
        )
    print(f"[mei-annotate] wrote {len(candidate)} rows to {out_path}")
    print(f"[mei-annotate] total annotate walltime={time.monotonic() - total_t0:.1f}s")
    return out_path
