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


@dataclass
class ClipAlignmentSummary:
    sample: str
    clip_count: int
    paf_hits: int


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
                    "target_strand": strand,
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
                "target_strand",
                "mapq",
                "pid",
                "qcov",
                "mei_score",
                "family",
            ]
        )

    paf = pd.DataFrame(rows)
    paf = paf.sort_values(["qname", "mei_score", "mapq", "pid", "qcov"], ascending=[True, False, False, False, False])
    return paf.drop_duplicates(subset=["qname"], keep="first")


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
    if clips.empty:
        summary = ClipAlignmentSummary(sample=sample, clip_count=0, paf_hits=0)
        return clips, summary

    with tempfile.TemporaryDirectory(prefix=f"rtm_mei_{sample}_") as tmpdir:
        tmp = Path(tmpdir)
        query_fa = tmp / "clips.fa"
        paf_path = tmp / "clips_vs_mei.paf"

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
        out["mei_hit"] = out["target"].notna()
        out["family"] = out["family"].fillna("")
        out["target"] = out["target"].fillna("")
        out["target_strand"] = out["target_strand"].fillna("")
        out["target_start"] = out["target_start"].fillna(0).astype(int)
        out["target_end"] = out["target_end"].fillna(0).astype(int)
        out["mei_score"] = out["mei_score"].fillna(0.0).astype(float)

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
        out["mei_hit"] = out["target"].notna()
        out["family"] = out["family"].fillna("")
        out["target"] = out["target"].fillna("")
        out["target_strand"] = out["target_strand"].fillna("")
        out["target_start"] = out["target_start"].fillna(0).astype(int)
        out["target_end"] = out["target_end"].fillna(0).astype(int)
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


def _aggregate_side_metrics(df: pd.DataFrame, sample_prefix: str, side: str) -> pd.DataFrame:
    side_df = df.loc[(df["clip_side"] == side) & (df["mei_hit"])].copy()
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
                f"{sample_prefix}_{side}_mei_subfamily_purity",
                f"{sample_prefix}_{side}_mei_breakpoint_mode",
                f"{sample_prefix}_{side}_mei_breakpoint_mode_fraction",
                f"{sample_prefix}_{side}_mei_breakpoint_unique_positions",
                f"{sample_prefix}_{side}_poly_at_reads",
                f"{sample_prefix}_{side}_poly_at_fraction",
                f"{sample_prefix}_{side}_poly_at_max_run",
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

    family_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "family"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"family": f"{sample_prefix}_{side}_mei_family"})
    )
    subfamily_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "target"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target": f"{sample_prefix}_{side}_mei_subfamily"})
    )
    strand_top = (
        side_df.groupby(["chrom", "window_start", "window_end", "target_strand"], as_index=False)["mei_score"]
        .sum()
        .sort_values(["chrom", "window_start", "window_end", "mei_score"], ascending=[True, True, True, False])
        .drop_duplicates(["chrom", "window_start", "window_end"], keep="first")
        .rename(columns={"target_strand": f"{sample_prefix}_{side}_mei_strand"})
    )
    subfamily_totals = (
        side_df.groupby(["chrom", "window_start", "window_end", "target"], as_index=False)["mei_score"]
        .sum()
        .rename(columns={"mei_score": "subfamily_score_sum"})
    )
    subfamily_sum = (
        side_df.groupby(["chrom", "window_start", "window_end"], as_index=False)["mei_score"]
        .sum()
        .rename(columns={"mei_score": "all_subfamily_score_sum"})
    )
    purity = (
        subfamily_top.rename(columns={"mei_score": "top_subfamily_score_sum"})
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
    agg = (
        side_df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_{side}_mei_supported_reads": ("read_name", "nunique"),
                f"{sample_prefix}_{side}_mei_score_sum": ("mei_score", "sum"),
                f"{sample_prefix}_{side}_mei_start": ("target_start", "min"),
                f"{sample_prefix}_{side}_mei_end": ("target_end", "max"),
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
    )
    agg[f"{sample_prefix}_{side}_mei_breakpoint_mode_fraction"] = (
        agg["mode_support_reads"] / agg[f"{sample_prefix}_{side}_mei_supported_reads"]
    ).fillna(0.0)
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


def _infer_tumor_insertion_metrics(candidates: pd.DataFrame, reference_fasta: Path | None = None) -> pd.DataFrame:
    out = candidates.copy()
    for col in [
        "tumor_L_mei_start",
        "tumor_R_mei_start",
        "tumor_L_mei_end",
        "tumor_R_mei_end",
        "tumor_L_mei_breakpoint_mode",
        "tumor_R_mei_breakpoint_mode",
        "normal_L_mei_breakpoint_mode",
        "normal_R_mei_breakpoint_mode",
        "tumor_L_mei_supported_reads",
        "tumor_R_mei_supported_reads",
        "normal_L_mei_supported_reads",
        "normal_R_mei_supported_reads",
    ]:
        if col not in out.columns:
            out[col] = 0
        out[col] = out[col].fillna(0).astype(int)

    out["tumor_insertion_mei_start"] = out.apply(
        lambda r: min([x for x in [r["tumor_L_mei_start"], r["tumor_R_mei_start"]] if x > 0], default=0),
        axis=1,
    )
    out["tumor_insertion_mei_end"] = out.apply(
        lambda r: max([x for x in [r["tumor_L_mei_end"], r["tumor_R_mei_end"]] if x > 0], default=0),
        axis=1,
    )
    out["tumor_insertion_mei_span"] = (
        out["tumor_insertion_mei_end"] - out["tumor_insertion_mei_start"] + 1
    ).where((out["tumor_insertion_mei_start"] > 0) & (out["tumor_insertion_mei_end"] >= out["tumor_insertion_mei_start"]), 0)

    def orient(row: pd.Series) -> str:
        strands = [s for s in [row.get("tumor_L_mei_strand", ""), row.get("tumor_R_mei_strand", "")] if s in {"+", "-"}]
        if not strands:
            return ""
        if len(set(strands)) == 1:
            return strands[0]
        return "mixed"

    out["tumor_insertion_orientation"] = out.apply(orient, axis=1)

    def _pick_tsd_pair(row: pd.Series) -> tuple[int, int, int]:
        # Prefer strict bilateral pairs from either tumor or normal.
        # If strict length is invalid, allow a small ±2 bp breakpoint rescue.
        candidates: list[tuple[int, int, int, int]] = []
        t_l = int(row.get("tumor_L_mei_breakpoint_mode", 0))
        t_r = int(row.get("tumor_R_mei_breakpoint_mode", 0))
        t_support = int(row.get("tumor_L_mei_supported_reads", 0)) + int(row.get("tumor_R_mei_supported_reads", 0))
        if t_l > 0 and t_r > 0:
            candidates.append((t_l, t_r, t_support, 0))
        n_l = int(row.get("normal_L_mei_breakpoint_mode", 0))
        n_r = int(row.get("normal_R_mei_breakpoint_mode", 0))
        n_support = int(row.get("normal_L_mei_supported_reads", 0)) + int(row.get("normal_R_mei_supported_reads", 0))
        if n_l > 0 and n_r > 0:
            candidates.append((n_l, n_r, n_support, 1))
        if not candidates:
            return (0, 0, 0)

        # Try strict first (no coordinate adjustment).
        strict_ok: list[tuple[int, int, int]] = []
        for l, r, support, _sample_priority in candidates:
            tsd_len = int(r - l + 1)
            if 2 <= tsd_len <= 30:
                strict_ok.append((support, l, r))
        if strict_ok:
            strict_ok.sort(key=lambda x: (x[0], x[2] - x[1]), reverse=True)
            _support, best_l, best_r = strict_ok[0]
            return (best_l, best_r, int(best_r - best_l + 1))

        # Rescue with ±2 bp shift when strict pairing misses by a few bases.
        rescue: list[tuple[int, int, int, int, int]] = []
        for l, r, support, sample_priority in candidates:
            for dl in (-2, -1, 0, 1, 2):
                for dr in (-2, -1, 0, 1, 2):
                    ll = int(l + dl)
                    rr = int(r + dr)
                    if ll <= 0 or rr <= 0 or rr < ll:
                        continue
                    tsd_len = int(rr - ll + 1)
                    if 2 <= tsd_len <= 30:
                        shift_penalty = abs(dl) + abs(dr)
                        rescue.append((shift_penalty, -support, sample_priority, ll, rr))
        if not rescue:
            return (0, 0, 0)
        rescue.sort()
        _shift_penalty, _neg_support, _sample_priority, best_l, best_r = rescue[0]
        return (best_l, best_r, int(best_r - best_l + 1))

    tsd_pairs = out.apply(_pick_tsd_pair, axis=1, result_type="expand")
    tsd_pairs.columns = ["tsd_left_breakpoint", "tsd_right_breakpoint", "tsd_len_estimate"]
    out["tsd_left_breakpoint"] = tsd_pairs["tsd_left_breakpoint"].astype(int)
    out["tsd_right_breakpoint"] = tsd_pairs["tsd_right_breakpoint"].astype(int)
    out["tsd_len_estimate"] = tsd_pairs["tsd_len_estimate"].astype(int)
    # Strict TSD evidence threshold: 4 bp or longer.
    out["tsd_detected"] = out["tsd_len_estimate"] >= 4

    def _breakpoint_pos(row: pd.Series) -> int:
        l = int(row.get("tsd_left_breakpoint", 0))
        r = int(row.get("tsd_right_breakpoint", 0))
        if l > 0 and r > 0:
            return int((l + r) // 2)
        l = int(row.get("tumor_L_mei_breakpoint_mode", 0))
        r = int(row.get("tumor_R_mei_breakpoint_mode", 0))
        if l > 0 and r > 0:
            return int((l + r) // 2)
        l = int(row.get("normal_L_mei_breakpoint_mode", 0))
        r = int(row.get("normal_R_mei_breakpoint_mode", 0))
        if l > 0 and r > 0:
            return int((l + r) // 2)
        if l > 0:
            return l
        if r > 0:
            return r
        return 0

    out["tumor_insertion_breakpoint_pos"] = out.apply(_breakpoint_pos, axis=1).astype(int)
    out["tsd_seq"] = ""
    out["tumor_breakpoint_context_11bp"] = ""
    out["tumor_breakpoint_l1_en_hexamer"] = ""
    out["tumor_breakpoint_l1_en_pattern"] = ""
    out["tumor_breakpoint_context_11bp_oriented"] = ""
    out["tumor_breakpoint_l1_en_hexamer_oriented"] = ""
    out["tumor_breakpoint_l1_en_pattern_yy_rrrr"] = ""
    out["tumor_breakpoint_l1_en_orientation_source"] = "unknown"
    out["tumor_breakpoint_l1_en_motif_like"] = False
    out["tumor_breakpoint_l1_en_best_motif"] = ""
    out["tumor_breakpoint_l1_en_motif_type"] = ""
    out["tumor_breakpoint_l1_en_mismatches"] = 99
    out["tumor_breakpoint_l1_en_mismatch_tolerance"] = 0
    out["tumor_breakpoint_l1_en_best_match_seq"] = ""
    out["tumor_breakpoint_l1_en_best_match_offset"] = 0
    out["tumor_breakpoint_l1_en_best_match_strand"] = "unknown"
    out["tumor_breakpoint_l1_en_best_match_anchor_6mer"] = ""
    out["tumor_breakpoint_l1_en_best_match_pattern_yy_rrrr"] = ""
    out["tumor_breakpoint_yyrrrr_logodds"] = 0.0
    out["tumor_breakpoint_yyrrrr_logodds_shift1_max"] = 0.0
    out["tumor_breakpoint_yyrrrr_best_offset"] = 0
    out["tumor_breakpoint_yyrrrr_logodds_shift1_mt_adj"] = 0.0
    if reference_fasta is not None:
        with pysam.FastaFile(str(reference_fasta)) as ref:
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
                    try:
                        seqs.append(ref.fetch(chrom, start0, end0).upper())
                    except Exception:
                        seqs.append("")

                bp = int(getattr(row, "tumor_insertion_breakpoint_pos", 0))
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
                    yyrrrr_scores.append(0.0)
                    yyrrrr_shift1_scores.append(0.0)
                    yyrrrr_best_offsets.append(0)
                    yyrrrr_shift1_mt_adj_scores.append(0.0)
                    continue
                chrom = str(row.chrom)
                # 11 bp centered on breakpoint base (5 upstream + anchor + 5 downstream).
                start0_11 = max(0, bp - 6)
                end0_11 = max(start0_11 + 1, bp + 5)
                # 6 bp motif window near cleavage preference (4 upstream + 2 downstream).
                start0_6 = max(0, bp - 5)
                end0_6 = max(start0_6 + 1, bp + 1)
                try:
                    ctx11 = ref.fetch(chrom, start0_11, end0_11).upper()
                except Exception:
                    ctx11 = ""
                try:
                    hex6 = ref.fetch(chrom, start0_6, end0_6).upper()
                except Exception:
                    hex6 = ""
                patt = f"{hex6[:4]}/{hex6[4:6]}" if len(hex6) == 6 else ""
                oriented_hex6, oriented_ctx11, orientation_source = _orient_to_insertion_strand(
                    hexamer=hex6,
                    context11bp=ctx11,
                    orientation=str(getattr(row, "tumor_insertion_orientation", "")),
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
            out["tumor_breakpoint_context_11bp"] = contexts_11bp
            out["tumor_breakpoint_l1_en_hexamer"] = l1_hexamers
            out["tumor_breakpoint_l1_en_pattern"] = l1_patterns
            out["tumor_breakpoint_context_11bp_oriented"] = contexts_11bp_oriented
            out["tumor_breakpoint_l1_en_hexamer_oriented"] = l1_hexamers_oriented
            out["tumor_breakpoint_l1_en_pattern_yy_rrrr"] = l1_patterns_yy_rrrr
            out["tumor_breakpoint_l1_en_orientation_source"] = l1_orientation_source
            out["tumor_breakpoint_l1_en_motif_like"] = l1_like
            out["tumor_breakpoint_l1_en_best_motif"] = l1_best_motif
            out["tumor_breakpoint_l1_en_motif_type"] = l1_motif_type
            out["tumor_breakpoint_l1_en_mismatches"] = l1_mismatches
            out["tumor_breakpoint_l1_en_mismatch_tolerance"] = l1_tolerance
            out["tumor_breakpoint_l1_en_best_match_seq"] = l1_best_match_seq
            out["tumor_breakpoint_l1_en_best_match_offset"] = l1_best_match_offset
            out["tumor_breakpoint_l1_en_best_match_strand"] = l1_best_match_strand
            out["tumor_breakpoint_l1_en_best_match_anchor_6mer"] = l1_best_match_anchor_6mer
            out["tumor_breakpoint_l1_en_best_match_pattern_yy_rrrr"] = l1_best_match_pattern
            out["tumor_breakpoint_yyrrrr_logodds"] = yyrrrr_scores
            out["tumor_breakpoint_yyrrrr_logodds_shift1_max"] = yyrrrr_shift1_scores
            out["tumor_breakpoint_yyrrrr_best_offset"] = yyrrrr_best_offsets
            out["tumor_breakpoint_yyrrrr_logodds_shift1_mt_adj"] = yyrrrr_shift1_mt_adj_scores

    # Weighted coherence metrics for ranking (annotation-only, no hard filtering).
    out["tumor_breakpoint_mode_fraction_weighted"] = (
        out.get("tumor_L_mei_breakpoint_mode_fraction", 0.0) * out.get("tumor_L_mei_supported_reads", 0)
        + out.get("tumor_R_mei_breakpoint_mode_fraction", 0.0) * out.get("tumor_R_mei_supported_reads", 0)
    ) / out.get("tumor_mei_supported_reads", 0).replace(0, 1)
    out["normal_breakpoint_mode_fraction_weighted"] = (
        out.get("normal_L_mei_breakpoint_mode_fraction", 0.0) * out.get("normal_L_mei_supported_reads", 0)
        + out.get("normal_R_mei_breakpoint_mode_fraction", 0.0) * out.get("normal_R_mei_supported_reads", 0)
    ) / out.get("normal_mei_supported_reads", 0).replace(0, 1)
    out["tumor_subfamily_purity_weighted"] = (
        out.get("tumor_L_mei_subfamily_purity", 0.0) * out.get("tumor_L_mei_supported_reads", 0)
        + out.get("tumor_R_mei_subfamily_purity", 0.0) * out.get("tumor_R_mei_supported_reads", 0)
    ) / out.get("tumor_mei_supported_reads", 0).replace(0, 1)
    out["normal_subfamily_purity_weighted"] = (
        out.get("normal_L_mei_subfamily_purity", 0.0) * out.get("normal_L_mei_supported_reads", 0)
        + out.get("normal_R_mei_subfamily_purity", 0.0) * out.get("normal_R_mei_supported_reads", 0)
    ) / out.get("normal_mei_supported_reads", 0).replace(0, 1)
    mapq_scaled = (out.get("split_tumor_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0)
    out["coherence_score"] = (
        0.4 * out["tumor_breakpoint_mode_fraction_weighted"].fillna(0.0)
        + 0.4 * out["tumor_subfamily_purity_weighted"].fillna(0.0)
        + 0.2 * mapq_scaled.fillna(0.0)
    )
    out["normal_background_score"] = (
        out.get("normal_mei_supported_reads", 0).astype(float)
        + out.get("normal_total_rows", 0).astype(float)
    )

    out["tumor_poly_at_reads"] = out.get("tumor_L_poly_at_reads", 0).fillna(0).astype(int) + out.get(
        "tumor_R_poly_at_reads", 0
    ).fillna(0).astype(int)
    out["normal_poly_at_reads"] = out.get("normal_L_poly_at_reads", 0).fillna(0).astype(int) + out.get(
        "normal_R_poly_at_reads", 0
    ).fillna(0).astype(int)
    out["tumor_poly_at_max_run"] = (
        out.get("tumor_L_poly_at_max_run", 0).fillna(0).astype(int).combine(
            out.get("tumor_R_poly_at_max_run", 0).fillna(0).astype(int), max
        )
    )
    out["normal_poly_at_max_run"] = (
        out.get("normal_L_poly_at_max_run", 0).fillna(0).astype(int).combine(
            out.get("normal_R_poly_at_max_run", 0).fillna(0).astype(int), max
        )
    )
    out["tumor_poly_at_fraction_weighted"] = (
        out.get("tumor_L_poly_at_fraction", 0.0).fillna(0.0).astype(float) * out.get("tumor_L_mei_supported_reads", 0)
        + out.get("tumor_R_poly_at_fraction", 0.0).fillna(0.0).astype(float) * out.get("tumor_R_mei_supported_reads", 0)
    ) / out.get("tumor_mei_supported_reads", 0).replace(0, 1)
    out["normal_poly_at_fraction_weighted"] = (
        out.get("normal_L_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * out.get("normal_L_mei_supported_reads", 0)
        + out.get("normal_R_poly_at_fraction", 0.0).fillna(0.0).astype(float)
        * out.get("normal_R_mei_supported_reads", 0)
    ) / out.get("normal_mei_supported_reads", 0).replace(0, 1)

    normal_metrics = out.apply(
        lambda r: _sample_insertion_span_and_orientation(r, "normal"),
        axis=1,
        result_type="expand",
    )
    normal_metrics.columns = [
        "normal_insertion_mei_start",
        "normal_insertion_mei_end",
        "normal_insertion_mei_span",
        "normal_insertion_orientation",
    ]
    for col in normal_metrics.columns:
        out[col] = normal_metrics[col]
    out["insertion_orientation"] = out.apply(_choose_consolidated_insertion_orientation, axis=1)
    out["insertion_mei_span"] = out.apply(_choose_consolidated_insertion_mei_span, axis=1).astype(int)
    return out


def _row_int(row: pd.Series, key: str, default: int = 0) -> int:
    val = row.get(key, default)
    if pd.isna(val):
        return default
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


def _sample_insertion_span_and_orientation(row: pd.Series, prefix: str) -> tuple[int, int, int, str]:
    l_start = _row_int(row, f"{prefix}_L_mei_start")
    r_start = _row_int(row, f"{prefix}_R_mei_start")
    l_end = _row_int(row, f"{prefix}_L_mei_end")
    r_end = _row_int(row, f"{prefix}_R_mei_end")
    starts = [x for x in [l_start, r_start] if x > 0]
    ends = [x for x in [l_end, r_end] if x > 0]
    start = min(starts) if starts else 0
    end = max(ends) if ends else 0
    span = (end - start + 1) if start > 0 and end >= start else 0
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

    if span <= 0:
        l_target = _row_int(row, f"{prefix}_discordant_mei_left_target_pos_median")
        r_target = _row_int(row, f"{prefix}_discordant_mei_right_target_pos_median")
        if l_target > 0 and r_target > 0:
            start = min(l_target, r_target)
            end = max(l_target, r_target)
            span = end - start + 1

    if orient not in {"+", "-"}:
        discordant_strand = str(row.get(f"{prefix}_discordant_mei_strand", "") or "").strip()
        if discordant_strand in {"+", "-"}:
            orient = discordant_strand
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
    tumor_orient = str(row.get("tumor_insertion_orientation", "") or "").strip()
    normal_orient = str(row.get("normal_insertion_orientation", "") or "").strip()
    tumor_bilateral = _sample_has_bilateral_split_support(row, "tumor") or _sample_has_bilateral_discordant_support(
        row, "tumor"
    )
    normal_bilateral = _sample_has_bilateral_split_support(row, "normal") or _sample_has_bilateral_discordant_support(
        row, "normal"
    )
    if tumor_bilateral and tumor_orient in {"+", "-"}:
        return tumor_orient
    if normal_bilateral and normal_orient in {"+", "-"}:
        return normal_orient
    if tumor_orient in {"+", "-"}:
        return tumor_orient
    if normal_orient in {"+", "-"}:
        return normal_orient
    return _choose_event_orientation(row)


def _choose_consolidated_insertion_mei_span(row: pd.Series) -> int:
    tumor_span = _row_int(row, "tumor_insertion_mei_span")
    normal_span = _row_int(row, "normal_insertion_mei_span")
    tumor_bilateral = _sample_has_bilateral_split_support(row, "tumor") or _sample_has_bilateral_discordant_support(
        row, "tumor"
    )
    normal_bilateral = _sample_has_bilateral_split_support(row, "normal") or _sample_has_bilateral_discordant_support(
        row, "normal"
    )
    if tumor_bilateral and tumor_span > 0:
        return tumor_span
    if normal_bilateral and normal_span > 0:
        return normal_span
    if tumor_span > 0 and normal_span > 0:
        tumor_reads = _row_int(row, "tumor_mei_supported_reads")
        normal_reads = _row_int(row, "normal_mei_supported_reads")
        return tumor_span if tumor_reads >= normal_reads else normal_span
    return max(tumor_span, normal_span)


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

    for prefix in ("tumor", "normal"):
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

    out["poly_at_max_run"] = _max_int_series("tumor_poly_at_max_run", "normal_poly_at_max_run")
    out["poly_at_reads"] = _col_int("tumor_poly_at_reads") + _col_int("normal_poly_at_reads")
    return out


def _add_consolidated_event_fields(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    for prefix in ("tumor", "normal"):
        out[f"{prefix}_left_supported_reads"] = (
            out.get(f"{prefix}_L_mei_supported_reads", 0).fillna(0).astype(int)
            + out.get(f"{prefix}_discordant_mei_left_supported_reads", 0).fillna(0).astype(int)
        )
        out[f"{prefix}_right_supported_reads"] = (
            out.get(f"{prefix}_R_mei_supported_reads", 0).fillna(0).astype(int)
            + out.get(f"{prefix}_discordant_mei_right_supported_reads", 0).fillna(0).astype(int)
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


def _complex_locus_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = [
        "discordant_tumor_large_insert_fraction",
        "discordant_tumor_interchrom_fraction",
        "discordant_tumor_mate_unmapped_fraction",
        "discordant_tumor_same_strand_fraction",
        "discordant_tumor_improper_pair_fraction",
        "discordant_normal_large_insert_fraction",
        "discordant_normal_interchrom_fraction",
        "discordant_normal_mate_unmapped_fraction",
        "discordant_normal_same_strand_fraction",
        "discordant_normal_improper_pair_fraction",
    ]
    parts = [_df_col_float(df, col) for col in fraction_cols if col in df.columns]
    if not parts:
        return pd.Series(0.0, index=df.index)
    return pd.concat(parts, axis=1).max(axis=1)


def _complex_locus_strong_companion_fraction(df: pd.DataFrame) -> pd.Series:
    fraction_cols = [
        "discordant_tumor_large_insert_fraction",
        "discordant_tumor_interchrom_fraction",
        "discordant_tumor_mate_unmapped_fraction",
        "discordant_normal_large_insert_fraction",
        "discordant_normal_interchrom_fraction",
        "discordant_normal_mate_unmapped_fraction",
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

    best_motif = ""
    best_type = ""
    best_mm = 99
    best_seq = ""
    best_offset = 0
    best_strand = "forward"
    best_anchor6 = ""
    best_pattern = ""
    strand_queries = [("forward", q11)]
    if allow_reverse_scan:
        strand_queries.append(("reverse", _revcomp(q11)))

    for strand_label, q in strand_queries:
        # Candidate 6bp windows anchored around breakpoint with offsets -1/0/+1.
        for offset, start in [(-1, 0), (0, 1), (1, 2)]:
            end = start + 6
            if end > len(q):
                continue
            seq6 = q[start:end]
            for motif, mtype in _L1_EN_PAPER_MOTIFS.items():
                mlen = len(motif)
                windows = [seq6] if mlen == 6 else [seq6[:5], seq6[1:6]]
                for win in windows:
                    if len(win) != mlen:
                        continue
                    mm = _hamming(win, motif)
                    if mm < best_mm:
                        best_mm = mm
                        best_motif = motif
                        best_type = mtype
                        best_seq = win
                        best_offset = int(offset)
                        best_strand = strand_label
                        best_anchor6 = seq6
                        best_pattern = f"{seq6[:2]}/{seq6[2:6]}" if len(seq6) == 6 else ""

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
    out = candidates.copy()

    for col in [
        "tumor_L_mei_family",
        "tumor_R_mei_family",
        "tumor_L_mei_subfamily",
        "tumor_R_mei_subfamily",
        "tumor_L_mei_strand",
        "tumor_R_mei_strand",
    ]:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str)

    out["tumor_family_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["tumor_L_mei_family"], out["tumor_R_mei_family"])
    ]
    out["tumor_subfamily_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["tumor_L_mei_subfamily"], out["tumor_R_mei_subfamily"])
    ]
    out["tumor_strand_agreement"] = [
        _agreement_flag(a, b) for a, b in zip(out["tumor_L_mei_strand"], out["tumor_R_mei_strand"])
    ]
    out["normal_family_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            out.get("normal_L_mei_family", "").fillna("").astype(str),
            out.get("normal_R_mei_family", "").fillna("").astype(str),
        )
    ]
    out["normal_subfamily_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            out.get("normal_L_mei_subfamily", "").fillna("").astype(str),
            out.get("normal_R_mei_subfamily", "").fillna("").astype(str),
        )
    ]
    out["normal_strand_agreement"] = [
        _agreement_flag(a, b)
        for a, b in zip(
            out.get("normal_L_mei_strand", "").fillna("").astype(str),
            out.get("normal_R_mei_strand", "").fillna("").astype(str),
        )
    ]

    tumor_mei_reads = out.get("tumor_mei_supported_reads", 0).astype(float)
    normal_mei_reads = out.get("normal_mei_supported_reads", 0).astype(float)
    total_rows = out.get("tumor_total_rows", 0).astype(float).replace(0, 1.0)
    mei_enrichment = out.get("mei_score_enrichment_ratio", 0.0).astype(float)
    mei_enrichment_scaled = (mei_enrichment / (mei_enrichment + 1.0)).clip(lower=0.0, upper=1.0)
    mei_read_fraction = (tumor_mei_reads / total_rows).clip(lower=0.0, upper=1.0)

    # Event-centric confidence score: do not bias to tumor-only support.
    event_subfamily_purity = pd.concat(
        [
            out.get("tumor_subfamily_purity_weighted", 0.0).astype(float).fillna(0.0),
            out.get("normal_subfamily_purity_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    event_breakpoint_consistency = pd.concat(
        [
            out.get("tumor_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
            out.get("normal_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    event_family_agreement = pd.concat(
        [out["tumor_family_agreement"].astype(float), out["normal_family_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    event_subfamily_agreement = pd.concat(
        [out["tumor_subfamily_agreement"].astype(float), out["normal_subfamily_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    event_strand_agreement = pd.concat(
        [out["tumor_strand_agreement"].astype(float), out["normal_strand_agreement"].astype(float)],
        axis=1,
    ).max(axis=1)
    normal_mei_fraction = (
        normal_mei_reads / out.get("normal_total_rows", 0).astype(float).replace(0, 1.0)
    ).clip(lower=0.0, upper=1.0)
    event_mei_fraction = pd.concat([mei_read_fraction.fillna(0.0), normal_mei_fraction.fillna(0.0)], axis=1).max(axis=1)
    mapq_event = pd.concat(
        [
            (out.get("split_tumor_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0),
            (out.get("split_normal_mapq_mean", 0.0).astype(float) / 60.0).clip(lower=0.0, upper=1.0),
        ],
        axis=1,
    ).max(axis=1)

    tsd_boost = out.get("tsd_detected", False).fillna(False).astype(bool).astype(float)
    polyA_event = pd.concat(
        [
            out.get("tumor_poly_at_fraction_weighted", 0.0).astype(float).fillna(0.0),
            out.get("normal_poly_at_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1).clip(lower=0.0, upper=1.0)
    motif_boost = out.get("tumor_breakpoint_l1_en_motif_like", False).fillna(False).astype(bool).astype(float)
    motif_logodds = out.get("tumor_breakpoint_yyrrrr_logodds_shift1_mt_adj", 0.0).astype(float).fillna(0.0)
    motif_logodds_scaled = (motif_logodds / 6.0).clip(lower=0.0, upper=1.0)

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
    )
    base_score = (base_score + 0.05 * mapq_event).clip(lower=0.0, upper=1.0)

    # Track complex SV-like companion signatures without suppressing MEI detection.
    complex_companion_fraction = _complex_locus_companion_fraction(out)
    complex_strong_companion_fraction = _complex_locus_strong_companion_fraction(out)
    large_insert_fraction = _df_col_float(out, "discordant_tumor_large_insert_fraction")
    interchrom_fraction = _df_col_float(out, "discordant_tumor_interchrom_fraction")
    mate_unmapped_fraction = _df_col_float(out, "discordant_tumor_mate_unmapped_fraction")

    out["complex_sv_large_insert_flag"] = large_insert_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_interchrom_flag"] = interchrom_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_mate_unmapped_flag"] = mate_unmapped_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_companion_signal"] = complex_strong_companion_fraction >= _COMPLEX_LOCUS_STRONG_MIN_FRACTION
    out["complex_sv_signal_score"] = complex_companion_fraction
    out["mei_with_complex_sv_signature"] = out["complex_sv_companion_signal"] & (tumor_mei_reads >= 2)
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
    same_strand_fraction = _df_col_float(out, "discordant_tumor_same_strand_fraction")
    improper_pair_fraction = _df_col_float(out, "discordant_tumor_improper_pair_fraction")
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
    left_reads = out.get("tumor_L_mei_supported_reads", 0).astype(float)
    right_reads = out.get("tumor_R_mei_supported_reads", 0).astype(float)
    discordant_mei_reads = out.get("tumor_discordant_mei_supported_reads", 0).astype(float)
    left_mode_frac = out.get("tumor_L_mei_breakpoint_mode_fraction", 0.0).astype(float).fillna(0.0)
    right_mode_frac = out.get("tumor_R_mei_breakpoint_mode_fraction", 0.0).astype(float).fillna(0.0)
    left_purity = out.get("tumor_L_mei_subfamily_purity", 0.0).astype(float).fillna(0.0)
    right_purity = out.get("tumor_R_mei_subfamily_purity", 0.0).astype(float).fillna(0.0)
    out["tumor_two_sided_support"] = (left_reads >= 1) & (right_reads >= 1)
    out["tumor_two_sided_strong_support"] = (left_reads >= 2) & (right_reads >= 2)
    out["tumor_one_sided_split_support"] = ((left_reads >= 2) & (right_reads < 2)) | (
        (right_reads >= 2) & (left_reads < 2)
    )
    out["tumor_discordant_mei_strong_support"] = discordant_mei_reads >= 3
    dpe_left = out.get("tumor_discordant_mei_left_supported_reads", 0).astype(float)
    dpe_right = out.get("tumor_discordant_mei_right_supported_reads", 0).astype(float)
    dpe_family_purity = out.get("tumor_discordant_mei_family_purity", 0.0).astype(float).fillna(0.0)
    dpe_strand_purity = out.get("tumor_discordant_mei_strand_purity", 0.0).astype(float).fillna(0.0)
    dpe_geometry_consistent = (
        out.get("tumor_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    )
    dpe_self_consistent = (
        out.get("tumor_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    out["tumor_discordant_mei_two_sided_support"] = (dpe_left >= 1) & (dpe_right >= 1)
    out["tumor_discordant_mei_consistent_support"] = (
        out["tumor_discordant_mei_two_sided_support"]
        & (dpe_family_purity >= 0.60)
        & (dpe_strand_purity >= 0.60)
        & dpe_geometry_consistent
        & dpe_self_consistent
    )
    normal_left_reads = out.get("normal_L_mei_supported_reads", 0).astype(float)
    normal_right_reads = out.get("normal_R_mei_supported_reads", 0).astype(float)
    out["normal_two_sided_support"] = (normal_left_reads >= 1) & (normal_right_reads >= 1)
    normal_dpe_left = out.get("normal_discordant_mei_left_supported_reads", 0).astype(float)
    normal_dpe_right = out.get("normal_discordant_mei_right_supported_reads", 0).astype(float)
    normal_dpe_family_purity = out.get("normal_discordant_mei_family_purity", 0.0).astype(float).fillna(0.0)
    normal_dpe_strand_purity = out.get("normal_discordant_mei_strand_purity", 0.0).astype(float).fillna(0.0)
    normal_dpe_geometry_consistent = (
        out.get("normal_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
    )
    normal_dpe_self_consistent = (
        out.get("normal_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    out["normal_discordant_mei_consistent_support"] = (
        (normal_dpe_left >= 1)
        & (normal_dpe_right >= 1)
        & (normal_dpe_family_purity >= 0.60)
        & (normal_dpe_strand_purity >= 0.60)
        & normal_dpe_geometry_consistent
        & normal_dpe_self_consistent
    )
    tumor_left_mei_consistent = _split_side_mei_for_complex(
        left_reads,
        out.get("tumor_L_mei_subfamily_purity", 0.0).astype(float),
        out.get("tumor_L_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    tumor_right_mei_consistent = _split_side_mei_for_complex(
        right_reads,
        out.get("tumor_R_mei_subfamily_purity", 0.0).astype(float),
        out.get("tumor_R_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    tumor_left_anchor_complex = _discordant_anchor_side_is_complex(
        out.get("tumor_discordant_anchor_left_unique_reads", 0).astype(float),
        out.get("tumor_discordant_anchor_left_complex_reason_max_fraction", 0.0).astype(float),
        out.get("tumor_discordant_mei_left_supported_reads", 0).astype(float),
    )
    tumor_right_anchor_complex = _discordant_anchor_side_is_complex(
        out.get("tumor_discordant_anchor_right_unique_reads", 0).astype(float),
        out.get("tumor_discordant_anchor_right_complex_reason_max_fraction", 0.0).astype(float),
        out.get("tumor_discordant_mei_right_supported_reads", 0).astype(float),
    )
    out["tumor_discordant_anchor_left_complex_side"] = tumor_left_anchor_complex
    out["tumor_discordant_anchor_right_complex_side"] = tumor_right_anchor_complex
    out["tumor_mei_with_complex_sidepair"] = (
        (tumor_left_mei_consistent & tumor_right_anchor_complex)
        | (tumor_right_mei_consistent & tumor_left_anchor_complex)
    )

    normal_left_mei_consistent = _split_side_mei_for_complex(
        normal_left_reads,
        out.get("normal_L_mei_subfamily_purity", 0.0).astype(float),
        out.get("normal_L_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    normal_right_mei_consistent = _split_side_mei_for_complex(
        normal_right_reads,
        out.get("normal_R_mei_subfamily_purity", 0.0).astype(float),
        out.get("normal_R_mei_breakpoint_mode_fraction", 0.0).astype(float),
    )
    normal_left_anchor_complex = _discordant_anchor_side_is_complex(
        out.get("normal_discordant_anchor_left_unique_reads", 0).astype(float),
        out.get("normal_discordant_anchor_left_complex_reason_max_fraction", 0.0).astype(float),
        out.get("normal_discordant_mei_left_supported_reads", 0).astype(float),
    )
    normal_right_anchor_complex = _discordant_anchor_side_is_complex(
        out.get("normal_discordant_anchor_right_unique_reads", 0).astype(float),
        out.get("normal_discordant_anchor_right_complex_reason_max_fraction", 0.0).astype(float),
        out.get("normal_discordant_mei_right_supported_reads", 0).astype(float),
    )
    out["normal_discordant_anchor_left_complex_side"] = normal_left_anchor_complex
    out["normal_discordant_anchor_right_complex_side"] = normal_right_anchor_complex
    out["normal_mei_with_complex_sidepair"] = (
        (normal_left_mei_consistent & normal_right_anchor_complex)
        | (normal_right_mei_consistent & normal_left_anchor_complex)
    )
    out["tumor_two_sided_like_support"] = out["tumor_two_sided_strong_support"] | (
        out["tumor_one_sided_split_support"] & out["tumor_discordant_mei_strong_support"]
    ) | out["tumor_discordant_mei_consistent_support"]
    out["tumor_side_breakpoint_consistency"] = left_mode_frac.combine(right_mode_frac, min)
    out["tumor_side_subfamily_purity"] = left_purity.combine(right_purity, min)
    out["tumor_two_sided_family_consistent"] = out["tumor_two_sided_support"] & (out["tumor_family_agreement"] == 1)
    out["tumor_two_sided_subfamily_consistent"] = out["tumor_two_sided_support"] & (
        out["tumor_subfamily_agreement"] == 1
    )
    out["event_two_sided_like_support"] = (
        out["tumor_two_sided_like_support"]
        | out["normal_two_sided_support"]
        | out["normal_discordant_mei_consistent_support"]
    )
    out["event_family_consistent"] = (out["tumor_family_agreement"] == 1) | (out["normal_family_agreement"] == 1)
    out["event_strand_consistent"] = (out["tumor_strand_agreement"] == 1) | (out["normal_strand_agreement"] == 1)
    out["event_side_breakpoint_consistency"] = pd.concat(
        [
            out["tumor_side_breakpoint_consistency"].astype(float).fillna(0.0),
            out.get("normal_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0),
        ],
        axis=1,
    ).max(axis=1)
    out["event_polyA_or_tsd_or_motif"] = (
        (tsd_boost >= 1.0) | (polyA_event >= 0.20) | (motif_boost >= 1.0) | (motif_logodds_scaled >= 0.25)
    )
    out["event_quality_clean"] = (
        (out.get("junk_flag_count", 0).fillna(0).astype(int) == 0)
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
        & (out.get("coherence_score", 0.0).astype(float) >= 0.55)
    )
    provisional_one_sided = (
        (~high_conf_pass)
        & (out["insertion_model_score"] >= 0.55)
        & (
            out["event_two_sided_like_support"]
            | ((tumor_mei_reads + normal_mei_reads) >= 1)
            | ((discordant_mei_reads + out.get("normal_discordant_mei_supported_reads", 0).astype(float)) >= 1)
        )
        & out["event_family_consistent"]
        & (out.get("coherence_score", 0.0).astype(float) >= 0.50)
    )
    complex_sidepair_event = (
        out["tumor_mei_with_complex_sidepair"] | out["normal_mei_with_complex_sidepair"]
    )
    complex_sidepair_pass = (
        (~high_conf_pass)
        & (~provisional_one_sided)
        & complex_sidepair_event
        & out["event_family_consistent"]
        & out["event_strand_consistent"]
        & out["event_quality_clean"]
        & (out["insertion_model_score"] >= 0.50)
        & (out.get("coherence_score", 0.0).astype(float) >= 0.45)
    )
    out["complex_mei_event"] = complex_sidepair_event | out["mei_with_complex_sv_signature"]
    out["passes_insertion_model"] = high_conf_pass
    out["passes_insertion_model_provisional"] = provisional_one_sided
    out["passes_insertion_model_complex"] = complex_sidepair_pass
    out["insertion_call_tier"] = "none"
    out.loc[provisional_one_sided, "insertion_call_tier"] = "provisional_one_sided"
    out.loc[complex_sidepair_event, "insertion_call_tier"] = "mei_with_complex"
    out.loc[high_conf_pass, "insertion_call_tier"] = "high_conf_two_sided"

    # Sample presence: shared when both tumor and normal have MEI support (>=1 read each).
    out["sample_status_label"] = "low_support"
    out.loc[(tumor_mei_reads >= 1) & (normal_mei_reads >= 1), "sample_status_label"] = "shared"
    out.loc[(tumor_mei_reads >= 1) & (normal_mei_reads == 0), "sample_status_label"] = "somatic_only"
    out.loc[(tumor_mei_reads == 0) & (normal_mei_reads >= 1), "sample_status_label"] = "germline_only"

    # Explicit convenience flag for downstream filtering.
    out["likely_false_positive_germline_only"] = out["sample_status_label"] == "germline_only"
    return out


def _consistent_family_mask(df: pd.DataFrame) -> pd.Series:
    family_cols = [
        "tumor_L_mei_family",
        "tumor_R_mei_family",
        "normal_L_mei_family",
        "normal_R_mei_family",
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


def _mean_depth_for_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
) -> float:
    if end_1based < start_1based:
        return 0.0
    start0 = max(0, int(start_1based) - 1)
    end0 = max(start0 + 1, int(end_1based))
    try:
        cov = bam.count_coverage(chrom, start0, end0, quality_threshold=0, read_callback="all")
    except ValueError:
        return 0.0
    span = end0 - start0
    if span <= 0:
        return 0.0
    total_depth = float(sum(cov[0]) + sum(cov[1]) + sum(cov[2]) + sum(cov[3]))
    return float(total_depth) / float(span)


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
    return {
        "local_bam_mean_depth": _mean_depth_for_interval(bam, chrom=chrom, start_1based=start_1based, end_1based=end_1based),
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
    tumor_bam_path: Path,
    normal_bam_path: Path,
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
        "tumor_bam": _file_stamp(tumor_bam_path),
        "normal_bam": _file_stamp(normal_bam_path),
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
    tumor_bam_path: Path,
    normal_bam_path: Path,
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
    out["depth_filter_family_consistent"] = False
    out["depth_filter_two_sided_consistent"] = False
    out["depth_filter_pass"] = False
    out["tumor_local_bam_mean_depth"] = 0.0
    out["normal_local_bam_mean_depth"] = 0.0
    out["tumor_context_non_sv_reads"] = 0
    out["normal_context_non_sv_reads"] = 0
    out["tumor_context_mapq_mean"] = 0.0
    out["normal_context_mapq_mean"] = 0.0
    out["tumor_context_mapq_lt20_fraction"] = 0.0
    out["normal_context_mapq_lt20_fraction"] = 0.0
    out["tumor_context_nm_per_100bp_mean"] = 0.0
    out["normal_context_nm_per_100bp_mean"] = 0.0
    out["tumor_context_nm_per_100bp_p90"] = 0.0
    out["normal_context_nm_per_100bp_p90"] = 0.0
    out["tumor_mei_support_per_100x_bam_depth"] = 0.0
    out["normal_mei_support_per_100x_bam_depth"] = 0.0
    out["mei_support_per_100x_bam_depth_delta"] = 0.0
    out["mei_support_per_100x_bam_depth_ratio"] = 1.0

    if out.empty:
        return out
    t0 = time.monotonic()

    family_consistent = _consistent_family_mask(out)
    out["depth_filter_family_consistent"] = family_consistent
    tumor_two_sided = out.get("tumor_two_sided_support", False).fillna(False).astype(bool)
    normal_two_sided = out.get("normal_two_sided_support", False).fillna(False).astype(bool)
    tumor_family_consistent = out.get("tumor_family_agreement", 0).fillna(0).astype(int) == 1
    normal_family_consistent = out.get("normal_family_agreement", 0).fillna(0).astype(int) == 1
    tumor_orientation_consistent = out.get("tumor_strand_agreement", 0).fillna(0).astype(int) == 1
    normal_orientation_consistent = out.get("normal_strand_agreement", 0).fillna(0).astype(int) == 1
    two_sided_consistent = (
        (tumor_two_sided & tumor_family_consistent & tumor_orientation_consistent)
        | (normal_two_sided & normal_family_consistent & normal_orientation_consistent)
    )
    out["depth_filter_two_sided_consistent"] = two_sided_consistent
    silver_mask = out.get("silver_stage_pass", False).fillna(False).astype(bool)
    if silver_mask.any():
        depth_mask = silver_mask
    else:
        depth_mask = out.get("junk_flag_count", 999).fillna(999).astype(int) == 0
    out["depth_filter_pass"] = depth_mask
    idxs = out.index[depth_mask].tolist()
    if not idxs:
        print("[mei-annotate] empirical stage skipped: no loci passed empirical prefilter", flush=True)
        return out

    loci_for_empirical = out.loc[depth_mask, ["chrom", "window_start", "window_end"]].copy()
    cache_key = _empirical_cache_key(
        loci=loci_for_empirical,
        tumor_bam_path=tumor_bam_path,
        normal_bam_path=normal_bam_path,
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
    random_tumor_df = pd.DataFrame()
    random_normal_df = pd.DataFrame()
    cache_hit = False
    if empirical_cache_dir is not None:
        empirical_cache_dir.mkdir(parents=True, exist_ok=True)
        tumor_cache_path = empirical_cache_dir / f"{cache_key}.tumor.parquet"
        normal_cache_path = empirical_cache_dir / f"{cache_key}.normal.parquet"
        if tumor_cache_path.exists() and normal_cache_path.exists():
            try:
                random_tumor_df = pd.read_parquet(tumor_cache_path)
                random_normal_df = pd.read_parquet(normal_cache_path)
                cache_hit = True
                print(
                    f"[mei-annotate] empirical cache hit key={cache_key} "
                    f"rows={len(random_tumor_df)}",
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

    with pysam.AlignmentFile(str(tumor_bam_path), "rb") as tumor_bam, pysam.AlignmentFile(
        str(normal_bam_path), "rb"
    ) as normal_bam:
        total_loci = int(len(idxs))
        loci_progress_every = 100
        print(f"[mei-annotate] empirical stage: computing context metrics for {total_loci} loci", flush=True)
        for i, idx in enumerate(idxs, start=1):
            row = out.loc[idx]
            chrom = str(row["chrom"])
            start = int(row["window_start"])
            end = int(row["window_end"])
            t_metrics = _context_quality_metrics_for_interval(tumor_bam, chrom=chrom, start_1based=start, end_1based=end)
            n_metrics = _context_quality_metrics_for_interval(normal_bam, chrom=chrom, start_1based=start, end_1based=end)
            out.at[idx, "tumor_local_bam_mean_depth"] = float(t_metrics["local_bam_mean_depth"])
            out.at[idx, "normal_local_bam_mean_depth"] = float(n_metrics["local_bam_mean_depth"])
            out.at[idx, "tumor_context_non_sv_reads"] = int(t_metrics["context_non_sv_reads"])
            out.at[idx, "normal_context_non_sv_reads"] = int(n_metrics["context_non_sv_reads"])
            out.at[idx, "tumor_context_mapq_mean"] = float(t_metrics["context_mapq_mean"])
            out.at[idx, "normal_context_mapq_mean"] = float(n_metrics["context_mapq_mean"])
            out.at[idx, "tumor_context_mapq_lt20_fraction"] = float(t_metrics["context_mapq_lt20_fraction"])
            out.at[idx, "normal_context_mapq_lt20_fraction"] = float(n_metrics["context_mapq_lt20_fraction"])
            out.at[idx, "tumor_context_nm_per_100bp_mean"] = float(t_metrics["context_nm_per_100bp_mean"])
            out.at[idx, "normal_context_nm_per_100bp_mean"] = float(n_metrics["context_nm_per_100bp_mean"])
            out.at[idx, "tumor_context_nm_per_100bp_p90"] = float(t_metrics["context_nm_per_100bp_p90"])
            out.at[idx, "normal_context_nm_per_100bp_p90"] = float(n_metrics["context_nm_per_100bp_p90"])
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
                bam=tumor_bam,
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
            random_tumor_rows: list[dict[str, float | int | str]] = []
            random_normal_rows: list[dict[str, float | int | str]] = []
            random_progress_every = 200
            total_random = int(len(random_windows))
            for i, rw in enumerate(random_windows.itertuples(index=False), start=1):
                chrom = str(rw.chrom)
                start = int(rw.window_start)
                end = int(rw.window_end)
                t_metrics = _context_quality_metrics_for_interval(
                    tumor_bam, chrom=chrom, start_1based=start, end_1based=end
                )
                n_metrics = _context_quality_metrics_for_interval(
                    normal_bam, chrom=chrom, start_1based=start, end_1based=end
                )
                random_tumor_rows.append(
                    {
                        "chrom": chrom,
                        "local_bam_mean_depth": float(t_metrics["local_bam_mean_depth"]),
                        "context_mapq_mean": float(t_metrics["context_mapq_mean"]),
                        "context_mapq_lt20_fraction": float(t_metrics["context_mapq_lt20_fraction"]),
                        "context_nm_per_100bp_mean": float(t_metrics["context_nm_per_100bp_mean"]),
                        "context_nm_per_100bp_p90": float(t_metrics["context_nm_per_100bp_p90"]),
                    }
                )
                random_normal_rows.append(
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
            random_tumor_df = pd.DataFrame(random_tumor_rows)
            random_normal_df = pd.DataFrame(random_normal_rows)
            if empirical_cache_dir is not None:
                tumor_cache_path = empirical_cache_dir / f"{cache_key}.tumor.parquet"
                normal_cache_path = empirical_cache_dir / f"{cache_key}.normal.parquet"
                random_tumor_df.to_parquet(tumor_cache_path, index=False)
                random_normal_df.to_parquet(normal_cache_path, index=False)
                print(
                    f"[mei-annotate] empirical cache write key={cache_key} "
                    f"rows={len(random_tumor_df)}",
                    flush=True,
                )
        else:
            print("[mei-annotate] empirical stage: using cached random-window metrics", flush=True)

    t_mei = out.get("tumor_mei_supported_reads", 0).astype(float)
    n_mei = out.get("normal_mei_supported_reads", 0).astype(float)
    t_depth = out["tumor_local_bam_mean_depth"].astype(float)
    n_depth = out["normal_local_bam_mean_depth"].astype(float)
    out["tumor_mei_support_per_100x_bam_depth"] = (t_mei * 100.0) / t_depth.replace(0, 1.0)
    out["normal_mei_support_per_100x_bam_depth"] = (n_mei * 100.0) / n_depth.replace(0, 1.0)
    out["mei_support_per_100x_bam_depth_delta"] = (
        out["tumor_mei_support_per_100x_bam_depth"] - out["normal_mei_support_per_100x_bam_depth"]
    )
    out["mei_support_per_100x_bam_depth_ratio"] = (
        (out["tumor_mei_support_per_100x_bam_depth"] + 1e-3)
        / (out["normal_mei_support_per_100x_bam_depth"] + 1e-3)
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
    out["tumor_empirical_random_n"] = 0
    out["normal_empirical_random_n"] = 0
    for metric, tail in metric_specs:
        out[f"tumor_empirical_{metric}_percentile"] = 0.0
        out[f"tumor_empirical_{metric}_p_{tail}"] = 1.0
        out[f"normal_empirical_{metric}_percentile"] = 0.0
        out[f"normal_empirical_{metric}_p_{tail}"] = 1.0

    score_idx = out.index[depth_mask]
    if len(score_idx) > 0:
        tumor_for_scoring = out.loc[score_idx, ["chrom"]].copy()
        tumor_for_scoring["local_bam_mean_depth"] = out.loc[score_idx, "tumor_local_bam_mean_depth"].astype(float)
        tumor_for_scoring["context_mapq_mean"] = out.loc[score_idx, "tumor_context_mapq_mean"].astype(float)
        tumor_for_scoring["context_mapq_lt20_fraction"] = out.loc[score_idx, "tumor_context_mapq_lt20_fraction"].astype(
            float
        )
        tumor_for_scoring["context_nm_per_100bp_mean"] = out.loc[score_idx, "tumor_context_nm_per_100bp_mean"].astype(
            float
        )
        tumor_for_scoring["context_nm_per_100bp_p90"] = out.loc[score_idx, "tumor_context_nm_per_100bp_p90"].astype(float)
        normal_for_scoring = out.loc[score_idx, ["chrom"]].copy()
        normal_for_scoring["local_bam_mean_depth"] = out.loc[score_idx, "normal_local_bam_mean_depth"].astype(float)
        normal_for_scoring["context_mapq_mean"] = out.loc[score_idx, "normal_context_mapq_mean"].astype(float)
        normal_for_scoring["context_mapq_lt20_fraction"] = out.loc[score_idx, "normal_context_mapq_lt20_fraction"].astype(
            float
        )
        normal_for_scoring["context_nm_per_100bp_mean"] = out.loc[score_idx, "normal_context_nm_per_100bp_mean"].astype(
            float
        )
        normal_for_scoring["context_nm_per_100bp_p90"] = out.loc[score_idx, "normal_context_nm_per_100bp_p90"].astype(
            float
        )

        tumor_scored = _apply_empirical_context_scores(
            loci_metrics=tumor_for_scoring,
            random_metrics=random_tumor_df,
            sample_prefix="tumor",
            scope=str(empirical_random_scope),
            progress_every=200,
        )
        normal_scored = _apply_empirical_context_scores(
            loci_metrics=normal_for_scoring,
            random_metrics=random_normal_df,
            sample_prefix="normal",
            scope=str(empirical_random_scope),
            progress_every=200,
        )
        for col in tumor_scored.columns:
            if col.startswith("tumor_empirical_"):
                out.loc[score_idx, col] = tumor_scored[col].values
        for col in normal_scored.columns:
            if col.startswith("normal_empirical_"):
                out.loc[score_idx, col] = normal_scored[col].values
    print("[mei-annotate] empirical stage: applying empirical p-value scoring complete", flush=True)
    elapsed_total = time.monotonic() - t0
    print(f"[mei-annotate] empirical stage complete (elapsed={elapsed_total:.1f}s)", flush=True)
    print(f"[mei-annotate] empirical stage walltime={time.monotonic() - stage_start:.1f}s", flush=True)
    return out


def _add_local_depth_normalized_support(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    tumor_total = out.get("tumor_total_rows", 0).astype(float)
    normal_total = out.get("normal_total_rows", 0).astype(float)
    tumor_mei = out.get("tumor_mei_supported_reads", 0).astype(float)
    normal_mei = out.get("normal_mei_supported_reads", 0).astype(float)

    # Local informative depth proxy from candidate-building stage.
    out["tumor_local_informative_rows"] = tumor_total.fillna(0.0).astype(int)
    out["normal_local_informative_rows"] = normal_total.fillna(0.0).astype(int)

    out["tumor_mei_support_local_frac"] = (tumor_mei / tumor_total.replace(0, 1)).fillna(0.0)
    out["normal_mei_support_local_frac"] = (normal_mei / normal_total.replace(0, 1)).fillna(0.0)
    out["tumor_mei_support_per_100_local_rows"] = out["tumor_mei_support_local_frac"] * 100.0
    out["normal_mei_support_per_100_local_rows"] = out["normal_mei_support_local_frac"] * 100.0
    out["mei_local_support_frac_delta"] = (
        out["tumor_mei_support_local_frac"] - out["normal_mei_support_local_frac"]
    )
    out["mei_local_support_frac_ratio"] = (
        (out["tumor_mei_support_local_frac"] + 1e-4) / (out["normal_mei_support_local_frac"] + 1e-4)
    )
    return out


def _assign_bronze_silver_stages(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    out["bronze_stage_pass"] = True

    junk_clean = out.get("junk_flag_count", 999).fillna(999).astype(int) == 0
    t_left_split = out.get("tumor_L_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    t_right_split = out.get("tumor_R_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    t_left_disc = out.get("tumor_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1
    t_right_disc = out.get("tumor_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1
    n_left_split = out.get("normal_L_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    n_right_split = out.get("normal_R_mei_supported_reads", 0).fillna(0).astype(float) >= 1
    n_left_disc = out.get("normal_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1
    n_right_disc = out.get("normal_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1

    tumor_bilateral_any = (t_left_split | t_left_disc) & (t_right_split | t_right_disc)
    normal_bilateral_any = (n_left_split | n_left_disc) & (n_right_split | n_right_disc)
    out["silver_bilateral_support_any"] = tumor_bilateral_any | normal_bilateral_any
    t_left_poly = out.get("tumor_L_poly_at_reads", 0).fillna(0).astype(float) >= 1
    t_right_poly = out.get("tumor_R_poly_at_reads", 0).fillna(0).astype(float) >= 1
    n_left_poly = out.get("normal_L_poly_at_reads", 0).fillna(0).astype(float) >= 1
    n_right_poly = out.get("normal_R_poly_at_reads", 0).fillna(0).astype(float) >= 1

    tumor_split_consistent = (
        out.get("tumor_two_sided_support", False).fillna(False).astype(bool)
        & (out.get("tumor_family_agreement", 0).fillna(0).astype(int) == 1)
        & (out.get("tumor_strand_agreement", 0).fillna(0).astype(int) == 1)
    )
    normal_split_consistent = (
        out.get("normal_two_sided_support", False).fillna(False).astype(bool)
        & (out.get("normal_family_agreement", 0).fillna(0).astype(int) == 1)
        & (out.get("normal_strand_agreement", 0).fillna(0).astype(int) == 1)
    )

    tumor_disc_consistent = (
        out.get("tumor_discordant_mei_two_sided_support", False).fillna(False).astype(bool)
        & (out.get("tumor_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & (out.get("tumor_discordant_mei_strand_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & out.get("tumor_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
        & out.get("tumor_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )
    normal_disc_consistent = (
        (out.get("normal_discordant_mei_left_supported_reads", 0).fillna(0).astype(float) >= 1)
        & (out.get("normal_discordant_mei_right_supported_reads", 0).fillna(0).astype(float) >= 1)
        & (out.get("normal_discordant_mei_family_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & (out.get("normal_discordant_mei_strand_purity", 0.0).fillna(0.0).astype(float) >= 0.95)
        & out.get("normal_discordant_mei_geometry_consistent", False).fillna(False).astype(bool)
        & out.get("normal_discordant_mei_self_consistent", True).fillna(True).astype(bool)
    )

    event_family_consistent = out.get("event_family_consistent", False).fillna(False).astype(bool)
    event_strand_consistent = out.get("event_strand_consistent", False).fillna(False).astype(bool)

    # PolyA-rescue bilateral support:
    # one side has MEI anchor support and the opposite side has polyA-clipped support.
    # If orientation is known, enforce expected tail side:
    # + insertion => right-side polyA; - insertion => left-side polyA.
    tumor_ori = out.get("tumor_insertion_orientation", "").fillna("").astype(str)
    normal_ori = out.get("normal_discordant_mei_strand", "").fillna("").astype(str)
    t_poly_mei_any = (t_left_poly & (t_right_split | t_right_disc)) | (t_right_poly & (t_left_split | t_left_disc))
    n_poly_mei_any = (n_left_poly & (n_right_split | n_right_disc)) | (n_right_poly & (n_left_split | n_left_disc))
    t_poly_oriented = (
        ((tumor_ori == "+") & t_right_poly & (t_left_split | t_left_disc))
        | ((tumor_ori == "-") & t_left_poly & (t_right_split | t_right_disc))
    )
    n_poly_oriented = (
        ((normal_ori == "+") & n_right_poly & (n_left_split | n_left_disc))
        | ((normal_ori == "-") & n_left_poly & (n_right_split | n_right_disc))
    )
    poly_sidepair_support = (t_poly_mei_any & ((tumor_ori == "") | t_poly_oriented)) | (
        n_poly_mei_any & ((normal_ori == "") | n_poly_oriented)
    )
    out["silver_polyA_sidepair_support"] = poly_sidepair_support

    t_left_anchor_complex = out.get("tumor_discordant_anchor_left_complex_side", False).fillna(False).astype(bool)
    t_right_anchor_complex = out.get("tumor_discordant_anchor_right_complex_side", False).fillna(False).astype(bool)
    n_left_anchor_complex = out.get("normal_discordant_anchor_left_complex_side", False).fillna(False).astype(bool)
    n_right_anchor_complex = out.get("normal_discordant_anchor_right_complex_side", False).fillna(False).astype(bool)
    t_left_structural = t_left_split | t_left_disc | t_left_poly | t_left_anchor_complex
    t_right_structural = t_right_split | t_right_disc | t_right_poly | t_right_anchor_complex
    n_left_structural = n_left_split | n_left_disc | n_left_poly | n_left_anchor_complex
    n_right_structural = n_right_split | n_right_disc | n_right_poly | n_right_anchor_complex
    tumor_bilateral_structural = t_left_structural & t_right_structural
    normal_bilateral_structural = n_left_structural & n_right_structural
    out["silver_bilateral_structural_support"] = tumor_bilateral_structural | normal_bilateral_structural

    tumor_complex_sidepair = (
        (t_left_split | t_left_disc) & (t_right_anchor_complex | t_right_poly)
    ) | ((t_right_split | t_right_disc) & (t_left_anchor_complex | t_left_poly))
    normal_complex_sidepair = (
        (n_left_split | n_left_disc) & (n_right_anchor_complex | n_right_poly)
    ) | ((n_right_split | n_right_disc) & (n_left_anchor_complex | n_left_poly))
    out["silver_complex_sidepair_support"] = tumor_complex_sidepair | normal_complex_sidepair
    out["silver_complex_structural_consistent"] = (
        out["silver_bilateral_structural_support"]
        & out["silver_complex_sidepair_support"]
        & (
            out.get("tumor_mei_with_complex_sidepair", False).fillna(False).astype(bool)
            | out.get("normal_mei_with_complex_sidepair", False).fillna(False).astype(bool)
            | out.get("mei_with_complex_sv_signature", False).fillna(False).astype(bool)
        )
    )

    silver_consistency = (
        tumor_split_consistent | normal_split_consistent | tumor_disc_consistent | normal_disc_consistent
    )
    out["silver_consistency_pass"] = (
        silver_consistency
        | (event_family_consistent & event_strand_consistent)
        | (poly_sidepair_support & event_family_consistent)
        | out["silver_complex_structural_consistent"]
    )
    out["silver_discordant_two_sided_consistent"] = tumor_disc_consistent | normal_disc_consistent

    tumor_l_bp = out.get("tumor_L_mei_breakpoint_mode", 0).fillna(0).astype(int)
    tumor_r_bp = out.get("tumor_R_mei_breakpoint_mode", 0).fillna(0).astype(int)
    normal_l_bp = out.get("normal_L_mei_breakpoint_mode", 0).fillna(0).astype(int)
    normal_r_bp = out.get("normal_R_mei_breakpoint_mode", 0).fillna(0).astype(int)
    out["silver_split_breakpoint_resolved"] = (
        (t_left_split & (tumor_l_bp > 0))
        | (t_right_split & (tumor_r_bp > 0))
        | (n_left_split & (normal_l_bp > 0))
        | (n_right_split & (normal_r_bp > 0))
    )
    out["silver_insertion_span_resolved"] = out.get("insertion_mei_span", 0).fillna(0).astype(int) > 0
    out["silver_breakpoint_or_span_resolved"] = (
        out["silver_split_breakpoint_resolved"] | out["silver_insertion_span_resolved"]
    )

    out["silver_stage_pass"] = junk_clean & (
        (
            out["silver_bilateral_support_any"]
            | poly_sidepair_support
            | out["silver_complex_structural_consistent"]
            | out["silver_breakpoint_or_span_resolved"]
        )
        & (out["silver_consistency_pass"] | out["silver_breakpoint_or_span_resolved"])
    )
    out["analysis_stage_tier"] = "bronze"
    out.loc[out["silver_stage_pass"], "analysis_stage_tier"] = "silver"
    print(
        "[mei-annotate] stage counts "
        f"bronze={len(out)} silver={int(out['silver_stage_pass'].sum())}",
        flush=True,
    )
    return out


def _assign_gold_stage(candidates: pd.DataFrame, empirical_p_threshold: float = 0.001) -> pd.DataFrame:
    out = candidates.copy()
    out["gold_empirical_p_threshold"] = float(empirical_p_threshold)
    out["gold_empirical_eval_available"] = False
    out["gold_empirical_outlier"] = False
    out["gold_stage_pass"] = False
    out["gold_stage_fail_reason"] = ""

    p_cols = [
        "tumor_empirical_local_bam_mean_depth_p_high",
        "tumor_empirical_context_mapq_mean_p_low",
        "tumor_empirical_context_mapq_lt20_fraction_p_high",
        "tumor_empirical_context_nm_per_100bp_mean_p_high",
        "tumor_empirical_context_nm_per_100bp_p90_p_high",
        "normal_empirical_local_bam_mean_depth_p_high",
        "normal_empirical_context_mapq_mean_p_low",
        "normal_empirical_context_mapq_lt20_fraction_p_high",
        "normal_empirical_context_nm_per_100bp_mean_p_high",
        "normal_empirical_context_nm_per_100bp_p90_p_high",
    ]
    available_cols = [c for c in p_cols if c in out.columns]
    silver = out.get("silver_stage_pass", False).fillna(False).astype(bool)
    if available_cols:
        out["gold_empirical_eval_available"] = True
        pvals = out.loc[:, available_cols].fillna(1.0).astype(float)
        out["gold_empirical_outlier"] = (pvals < float(empirical_p_threshold)).any(axis=1)
        out["gold_stage_pass"] = silver & (~out["gold_empirical_outlier"])
        out.loc[silver & out["gold_empirical_outlier"], "gold_stage_fail_reason"] = "empirical_outlier"
    else:
        out["gold_stage_pass"] = silver
        out.loc[silver, "gold_stage_fail_reason"] = "empirical_not_available"

    out.loc[out["gold_stage_pass"], "analysis_stage_tier"] = "gold"
    print(
        "[mei-annotate] stage counts "
        f"silver={int(silver.sum())} gold={int(out['gold_stage_pass'].sum())}",
        flush=True,
    )
    return out


def _two_sided_support_mask(df: pd.DataFrame) -> pd.Series:
    tumor_left = df.get("tumor_left_supported_reads", 0).fillna(0).astype(int)
    tumor_right = df.get("tumor_right_supported_reads", 0).fillna(0).astype(int)
    normal_left = df.get("normal_left_supported_reads", 0).fillna(0).astype(int)
    normal_right = df.get("normal_right_supported_reads", 0).fillna(0).astype(int)
    bilateral = ((tumor_left >= 1) & (tumor_right >= 1)) | ((normal_left >= 1) & (normal_right >= 1))
    if "silver_bilateral_support_any" in df.columns:
        bilateral = bilateral | df["silver_bilateral_support_any"].fillna(False).astype(bool)
    return bilateral


def _prioritize_mei_candidates(candidates: pd.DataFrame, *, stage_first: bool = True) -> pd.DataFrame:
    """Rank loci by evidence strength for manual review."""
    out = candidates.copy()
    out["two_sided_support"] = _two_sided_support_mask(out)

    out["_prio_tsd"] = out.get("tsd_detected", False).fillna(False).astype(bool)
    out["_prio_two_sided"] = out["two_sided_support"].astype(bool)
    out["_prio_poly_at_max_run"] = out.get("poly_at_max_run", 0).fillna(0).astype(int)
    out["_prio_breakpoint_resolved"] = out.get("tumor_insertion_breakpoint_pos", 0).fillna(0).astype(int) > 0
    out["_prio_g1k_region"] = out.get("g1k_melt_region_id", "").fillna("").astype(str).str.strip() != ""
    out["_prio_insertion_mei_span"] = out.get("insertion_mei_span", 0).fillna(0).astype(int)
    out["_prio_split_reads"] = (
        out.get("tumor_split_mei_supported_reads", 0).fillna(0).astype(int)
        + out.get("normal_split_mei_supported_reads", 0).fillna(0).astype(int)
    )
    out["_prio_discordant_reads"] = (
        out.get("tumor_discordant_mei_supported_reads", 0).fillna(0).astype(int)
        + out.get("normal_discordant_mei_supported_reads", 0).fillna(0).astype(int)
    )

    sort_cols: list[str] = []
    ascending: list[bool] = []
    if stage_first:
        for col in ("gold_stage_pass", "silver_stage_pass"):
            if col in out.columns:
                sort_cols.append(col)
                ascending.append(False)
    sort_cols.extend(
        [
            "_prio_tsd",
            "_prio_two_sided",
            "_prio_poly_at_max_run",
            "_prio_breakpoint_resolved",
            "_prio_g1k_region",
            "_prio_insertion_mei_span",
            "_prio_split_reads",
            "_prio_discordant_reads",
        ]
    )
    ascending.extend([False] * 8)
    sorted_out = out.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    return sorted_out.drop(
        columns=[c for c in sorted_out.columns if c.startswith("_prio_")],
        errors="ignore",
    ).reset_index(drop=True)


def _build_gold_review_table(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    def _series_or_default(col: str, default: object) -> pd.Series:
        if col in out.columns:
            return out[col]
        return pd.Series([default] * len(out), index=out.index)

    out["breakpoint_resolved"] = out.get("tumor_insertion_breakpoint_pos", 0).fillna(0).astype(int) > 0
    out["tsd_or_polyA_supported"] = (
        out.get("tsd_detected", False).fillna(False).astype(bool)
        | (out.get("tumor_poly_at_reads", 0).fillna(0).astype(float) >= 1)
        | (out.get("normal_poly_at_reads", 0).fillna(0).astype(float) >= 1)
        | (out.get("poly_at_reads", 0).fillna(0).astype(float) >= 1)
    )
    out["consensus_sequence_label"] = _series_or_default("mei_subfamily", "").fillna("").astype(str)
    out["g1k_overlap"] = _series_or_default("g1k_melt_id", "").fillna("").astype(str) != ""
    out["g1k_mei_family"] = _series_or_default("g1k_melt_insertion_subfamily", "").fillna("").astype(str).apply(
        lambda x: _infer_mei_family_from_fields(hit_id=x, family_hint=x, extra_hint="")
    )
    out["two_sided_support"] = _two_sided_support_mask(out)

    selected_cols = [
        "analysis_stage_tier",
        "gold_stage_pass",
        "silver_stage_pass",
        "sample_status_label",
        "insertion_call_tier",
        "complex_mei_event",
        "complex_sv_signature_label",
        "mei_with_complex_sv_signature",
        "chrom",
        "window_start",
        "window_end",
        "tumor_split_mei_supported_reads",
        "normal_split_mei_supported_reads",
        "tumor_discordant_mei_supported_reads",
        "normal_discordant_mei_supported_reads",
        "tumor_left_supported_reads",
        "tumor_right_supported_reads",
        "normal_left_supported_reads",
        "normal_right_supported_reads",
        "tumor_family_agreement",
        "tumor_strand_agreement",
        "normal_family_agreement",
        "normal_strand_agreement",
        "tumor_insertion_breakpoint_pos",
        "breakpoint_resolved",
        "two_sided_support",
        "insertion_orientation",
        "insertion_mei_span",
        "tsd_detected",
        "tsd_len_estimate",
        "tsd_seq",
        "tumor_poly_at_reads",
        "tumor_poly_at_max_run",
        "tumor_poly_at_fraction_weighted",
        "normal_poly_at_reads",
        "normal_poly_at_max_run",
        "normal_poly_at_fraction_weighted",
        "poly_at_reads",
        "poly_at_max_run",
        "tsd_or_polyA_supported",
        "mei_family",
        "mei_subfamily",
        "consensus_sequence_label",
        "nested_same_class_orientation",
        "g1k_overlap",
        "g1k_melt_id",
        "g1k_mei_family",
        "g1k_melt_insertion_subfamily",
        "g1k_melt_insertion_length",
        "g1k_melt_tsd",
        "g1k_melt_region_id",
        "tumor_empirical_local_bam_mean_depth_p_high",
        "tumor_empirical_context_mapq_mean_p_low",
        "tumor_empirical_context_mapq_lt20_fraction_p_high",
        "tumor_empirical_context_nm_per_100bp_mean_p_high",
        "tumor_empirical_context_nm_per_100bp_p90_p_high",
        "normal_empirical_local_bam_mean_depth_p_high",
        "normal_empirical_context_mapq_mean_p_low",
        "normal_empirical_context_mapq_lt20_fraction_p_high",
        "normal_empirical_context_nm_per_100bp_mean_p_high",
        "normal_empirical_context_nm_per_100bp_p90_p_high",
        "gold_empirical_outlier",
        "gold_stage_fail_reason",
        "insertion_model_score",
        "coherence_score",
        "mei_score_enrichment_ratio",
    ]
    for col in selected_cols:
        if col not in out.columns:
            out[col] = ""
    review = out.loc[:, selected_cols].copy()
    return _prioritize_mei_candidates(review, stage_first=True)


def _infer_mei_family_from_fields(hit_id: str, family_hint: str, extra_hint: str) -> str:
    txt = " ".join([hit_id or "", family_hint or "", extra_hint or ""]).upper()
    if "ALU" in txt:
        return "ALU"
    if "SVA" in txt:
        return "SVA"
    if "LINE1" in txt or "L1" in txt:
        return "LINE1"
    if "HERV" in txt or "ERV" in txt:
        return "ERV"
    return "OTHER"


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
    breakpoint_pos = int(row.get("tumor_insertion_breakpoint_pos", 0))
    if breakpoint_pos <= 0:
        breakpoint_pos = midpoint

    left_split = int(row.get("tumor_L_mei_supported_reads", 0))
    right_split = int(row.get("tumor_R_mei_supported_reads", 0))
    split_total = int(row.get("tumor_split_mei_supported_reads", 0))
    split_resolved = (split_total >= 2) or ((left_split >= 1) and (right_split >= 1))

    tumor_dpe = int(row.get("tumor_discordant_mei_supported_reads", 0))
    normal_dpe = int(row.get("normal_discordant_mei_supported_reads", 0))
    dpe_present = (tumor_dpe + normal_dpe) > 0
    dpe_tlen_mean = max(
        float(row.get("discordant_tumor_abs_tlen_mean", 0.0) or 0.0),
        float(row.get("discordant_normal_abs_tlen_mean", 0.0) or 0.0),
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
        ("tumor_discordant_mei_subfamily", "tumor_discordant_mei_supported_reads"),
        ("tumor_discordant_mei_left_subfamily", "tumor_discordant_mei_left_supported_reads"),
        ("tumor_discordant_mei_right_subfamily", "tumor_discordant_mei_right_supported_reads"),
        ("tumor_L_mei_subfamily", "tumor_L_mei_supported_reads"),
        ("tumor_R_mei_subfamily", "tumor_R_mei_supported_reads"),
        ("normal_discordant_mei_subfamily", "normal_discordant_mei_supported_reads"),
        ("normal_discordant_mei_left_subfamily", "normal_discordant_mei_left_supported_reads"),
        ("normal_discordant_mei_right_subfamily", "normal_discordant_mei_right_supported_reads"),
        ("normal_L_mei_subfamily", "normal_L_mei_supported_reads"),
        ("normal_R_mei_subfamily", "normal_R_mei_supported_reads"),
    ]:
        label = str(row.get(subfamily_col, "") or "").strip()
        weight = _row_int(row, weight_col)
        if label and weight > 0:
            weighted.append((label, weight))
    g1k_subfamily = str(row.get("g1k_melt_insertion_subfamily", "") or "").strip()
    if g1k_subfamily:
        weighted.append((g1k_subfamily, 1))
    if not weighted:
        return ""
    weighted.sort(key=lambda item: item[1], reverse=True)
    return weighted[0][0]


def _choose_event_family(row: pd.Series) -> str:
    candidates = [
        str(row.get("tumor_discordant_mei_family", "")),
        str(row.get("tumor_L_mei_family", "")),
        str(row.get("tumor_R_mei_family", "")),
        str(row.get("normal_discordant_mei_family", "")),
        str(row.get("normal_L_mei_family", "")),
        str(row.get("normal_R_mei_family", "")),
        str(row.get("tumor_discordant_mei_subfamily", "")),
        str(row.get("tumor_L_mei_subfamily", "")),
        str(row.get("tumor_R_mei_subfamily", "")),
        str(row.get("normal_discordant_mei_subfamily", "")),
        str(row.get("normal_L_mei_subfamily", "")),
        str(row.get("normal_R_mei_subfamily", "")),
        str(row.get("g1k_melt_insertion_subfamily", "")),
        str(row.get("g1k_melt_id", "")),
    ]
    for c in candidates:
        fam = _normalize_mei_family_token(c)
        if fam:
            return fam
    return ""


def _choose_event_orientation(row: pd.Series) -> str:
    candidates = [
        str(row.get("insertion_orientation", "")),
        str(row.get("tumor_insertion_orientation", "")),
        str(row.get("normal_insertion_orientation", "")),
        str(row.get("tumor_discordant_mei_strand", "")),
        str(row.get("tumor_L_mei_strand", "")),
        str(row.get("tumor_R_mei_strand", "")),
        str(row.get("normal_discordant_mei_strand", "")),
        str(row.get("normal_L_mei_strand", "")),
        str(row.get("normal_R_mei_strand", "")),
    ]
    for c in candidates:
        cc = (c or "").strip()
        if cc in {"+", "-"}:
            return cc
    return ""


def _load_rmsk_interval_trees(rmsk_table_path: Path) -> dict[str, IntervalTree]:
    trees: dict[str, IntervalTree] = {}
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
    return trees


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

    trees = _load_rmsk_interval_trees(rmsk_table_path)
    selected: list[dict[str, object]] = []
    for row in out.itertuples(index=False):
        as_row = pd.Series(row._asdict())
        chrom = str(getattr(row, "chrom"))
        pos_1based = int(getattr(row, "tumor_insertion_breakpoint_pos", 0))
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
    tumor_bam_path: Path | None = None,
    normal_bam_path: Path | None = None,
    rmsk_table_path: Path | None = None,
    g1k_mei_vcf: Path | None = None,
    g1k_split_padding_bp: int = 200,
    g1k_dpe_padding_min_bp: int = 200,
    g1k_dpe_padding_max_bp: int = 200,
    g1k_dpe_padding_tlen_factor: float = 0.0,
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
) -> Path:
    total_t0 = time.monotonic()
    candidate = pd.read_csv(candidate_loci_path, sep="\t")

    split_tumor_raw = _load_table(evidence_dir, "split_evidence", "tumor")
    split_normal_raw = _load_table(evidence_dir, "split_evidence", "normal")
    discordant_tumor_raw = _load_table(evidence_dir, "discordant_evidence", "tumor")
    discordant_normal_raw = _load_table(evidence_dir, "discordant_evidence", "normal")
    split_tumor = _assign_rows_to_candidate_loci(split_tumor_raw, candidate)
    split_normal = _assign_rows_to_candidate_loci(split_normal_raw, candidate)
    discordant_tumor = _assign_rows_to_candidate_loci(discordant_tumor_raw, candidate)
    discordant_normal = _assign_rows_to_candidate_loci(discordant_normal_raw, candidate)

    tumor_hits, tumor_summary = _align_clips_with_minimap2(split_tumor, mei_fasta, sample="tumor")
    normal_hits, normal_summary = _align_clips_with_minimap2(split_normal, mei_fasta, sample="normal")
    tumor_disc_hits, tumor_disc_summary = _align_discordant_reads_with_minimap2(
        discordant_tumor, mei_fasta, sample="tumor"
    )
    normal_disc_hits, normal_disc_summary = _align_discordant_reads_with_minimap2(
        discordant_normal, mei_fasta, sample="normal"
    )

    print(
        f"[mei-annotate] tumor clips={tumor_summary.clip_count} hits={tumor_summary.paf_hits}; "
        f"normal clips={normal_summary.clip_count} hits={normal_summary.paf_hits}; "
        f"tumor discordant reads={tumor_disc_summary.clip_count} hits={tumor_disc_summary.paf_hits}; "
        f"normal discordant reads={normal_disc_summary.clip_count} hits={normal_disc_summary.paf_hits}"
    )

    anno_parts = []
    for sample_prefix, df in (("tumor", tumor_hits), ("normal", normal_hits)):
        for side in ("L", "R"):
            anno_parts.append(_aggregate_side_metrics(df, sample_prefix=sample_prefix, side=side))

    for idx, part in enumerate(anno_parts):
        if part.empty:
            continue
        candidate = candidate.merge(part, on=["chrom", "window_start", "window_end"], how="left")
        if (idx + 1) % 2 == 0:
            print(f"[mei-annotate] merged side metrics {idx + 1}/{len(anno_parts)}")

    disc_t = _aggregate_discordant_mei_metrics(tumor_disc_hits, sample_prefix="tumor")
    disc_n = _aggregate_discordant_mei_metrics(normal_disc_hits, sample_prefix="normal")
    disc_anchor_t = _aggregate_discordant_anchor_side_metrics(discordant_tumor, sample_prefix="tumor")
    disc_anchor_n = _aggregate_discordant_anchor_side_metrics(discordant_normal, sample_prefix="normal")
    if not disc_t.empty:
        candidate = candidate.merge(disc_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_n.empty:
        candidate = candidate.merge(disc_n, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_anchor_t.empty:
        candidate = candidate.merge(disc_anchor_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_anchor_n.empty:
        candidate = candidate.merge(disc_anchor_n, on=["chrom", "window_start", "window_end"], how="left")
    print("[mei-annotate] merged discordant MEI support metrics")

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
    candidate = candidate.copy()

    candidate["tumor_split_mei_score_sum"] = candidate.get("tumor_L_mei_score_sum", 0.0) + candidate.get(
        "tumor_R_mei_score_sum", 0.0
    )
    candidate["normal_split_mei_score_sum"] = candidate.get("normal_L_mei_score_sum", 0.0) + candidate.get(
        "normal_R_mei_score_sum", 0.0
    )
    candidate["tumor_split_mei_supported_reads"] = candidate.get("tumor_L_mei_supported_reads", 0) + candidate.get(
        "tumor_R_mei_supported_reads", 0
    )
    candidate["normal_split_mei_supported_reads"] = candidate.get("normal_L_mei_supported_reads", 0) + candidate.get(
        "normal_R_mei_supported_reads", 0
    )
    candidate["tumor_discordant_mei_supported_reads"] = candidate.get("tumor_discordant_mei_supported_reads", 0).fillna(0).astype(int)
    candidate["normal_discordant_mei_supported_reads"] = candidate.get("normal_discordant_mei_supported_reads", 0).fillna(0).astype(int)
    candidate["tumor_discordant_mei_score_sum"] = candidate.get("tumor_discordant_mei_score_sum", 0.0).fillna(0.0).astype(float)
    candidate["normal_discordant_mei_score_sum"] = candidate.get("normal_discordant_mei_score_sum", 0.0).fillna(0.0).astype(float)

    candidate["tumor_mei_score_sum"] = candidate["tumor_split_mei_score_sum"] + candidate["tumor_discordant_mei_score_sum"]
    candidate["normal_mei_score_sum"] = candidate["normal_split_mei_score_sum"] + candidate["normal_discordant_mei_score_sum"]
    candidate["tumor_mei_supported_reads"] = (
        candidate["tumor_split_mei_supported_reads"] + candidate["tumor_discordant_mei_supported_reads"]
    )
    candidate["normal_mei_supported_reads"] = (
        candidate["normal_split_mei_supported_reads"] + candidate["normal_discordant_mei_supported_reads"]
    )
    candidate["mei_score_enrichment_ratio"] = (candidate["tumor_mei_score_sum"] + 0.1) / (
        candidate["normal_mei_score_sum"] + 0.1
    )

    candidate = _add_local_depth_normalized_support(candidate)
    candidate = _infer_tumor_insertion_metrics(candidate, reference_fasta=reference_fasta)
    candidate = _compute_insertion_model_scores(candidate)
    candidate = _assign_bronze_silver_stages(candidate)
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
    candidate = _add_consolidated_event_fields(candidate)
    candidate = _broaden_poly_at_fields(candidate)
    if rmsk_table_path is not None:
        rmsk_t0 = time.monotonic()
        candidate = _annotate_nested_retrotransposon(candidate, rmsk_table_path=rmsk_table_path)
        print(
            f"[mei-annotate] added nested-retrotransposon overlap annotation "
            f"(elapsed={time.monotonic() - rmsk_t0:.1f}s)"
        )
    if tumor_bam_path is not None and normal_bam_path is not None:
        emp_t0 = time.monotonic()
        candidate = _annotate_bam_depth_for_consistent_loci(
            candidate,
            tumor_bam_path=tumor_bam_path,
            normal_bam_path=normal_bam_path,
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
            f"[mei-annotate] added BAM-depth normalization for family-consistent, junk-clean loci "
            f"(elapsed={time.monotonic() - emp_t0:.1f}s)"
        )
    candidate = _assign_gold_stage(candidate)

    candidate = _prioritize_mei_candidates(candidate, stage_first=True)

    candidate.to_csv(out_path, sep="\t", index=False)
    candidate.to_parquet(out_path.with_suffix(".parquet"), index=False)
    gold_review = _build_gold_review_table(candidate)
    gold_review_path = out_path.with_name(out_path.stem + ".gold_review.tsv")
    gold_review.to_csv(gold_review_path, sep="\t", index=False)
    print(f"[mei-annotate] wrote gold review table to {gold_review_path}")
    print(f"[mei-annotate] wrote {len(candidate)} rows to {out_path}")
    print(f"[mei-annotate] total annotate walltime={time.monotonic() - total_t0:.1f}s")
    return out_path
