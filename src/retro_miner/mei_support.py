from __future__ import annotations

import math
import re
import subprocess
import tempfile
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
            ]
        )

    mei_df["locus_midpoint"] = (mei_df["window_start"].astype(int) + mei_df["window_end"].astype(int)) // 2
    mei_df["anchor_side"] = mei_df.apply(
        lambda r: "L" if int(r["pos"]) <= int(r["locus_midpoint"]) else "R",
        axis=1,
    )

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
    return agg


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
    large_insert_fraction = out.get("discordant_tumor_large_insert_fraction", 0.0).astype(float).fillna(0.0)
    interchrom_fraction = out.get("discordant_tumor_interchrom_fraction", 0.0).astype(float).fillna(0.0)
    mate_unmapped_fraction = out.get("discordant_tumor_mate_unmapped_fraction", 0.0).astype(float).fillna(0.0)

    out["complex_sv_large_insert_flag"] = large_insert_fraction >= 0.60
    out["complex_sv_interchrom_flag"] = interchrom_fraction >= 0.60
    out["complex_sv_mate_unmapped_flag"] = mate_unmapped_fraction >= 0.60
    out["complex_sv_companion_signal"] = (
        out["complex_sv_large_insert_flag"]
        | out["complex_sv_interchrom_flag"]
        | out["complex_sv_mate_unmapped_flag"]
    )
    out["complex_sv_signal_score"] = (
        pd.concat([large_insert_fraction, interchrom_fraction, mate_unmapped_fraction], axis=1).max(axis=1)
    )
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
    out["tumor_discordant_mei_two_sided_support"] = (dpe_left >= 1) & (dpe_right >= 1)
    out["tumor_discordant_mei_consistent_support"] = (
        out["tumor_discordant_mei_two_sided_support"] & (dpe_family_purity >= 0.60) & (dpe_strand_purity >= 0.60)
    )
    normal_left_reads = out.get("normal_L_mei_supported_reads", 0).astype(float)
    normal_right_reads = out.get("normal_R_mei_supported_reads", 0).astype(float)
    out["normal_two_sided_support"] = (normal_left_reads >= 1) & (normal_right_reads >= 1)
    normal_dpe_left = out.get("normal_discordant_mei_left_supported_reads", 0).astype(float)
    normal_dpe_right = out.get("normal_discordant_mei_right_supported_reads", 0).astype(float)
    normal_dpe_family_purity = out.get("normal_discordant_mei_family_purity", 0.0).astype(float).fillna(0.0)
    normal_dpe_strand_purity = out.get("normal_discordant_mei_strand_purity", 0.0).astype(float).fillna(0.0)
    out["normal_discordant_mei_consistent_support"] = (
        (normal_dpe_left >= 1)
        & (normal_dpe_right >= 1)
        & (normal_dpe_family_purity >= 0.60)
        & (normal_dpe_strand_purity >= 0.60)
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

    out["passes_insertion_model"] = high_conf_pass
    out["passes_insertion_model_provisional"] = provisional_one_sided
    out["insertion_call_tier"] = "none"
    out.loc[provisional_one_sided, "insertion_call_tier"] = "provisional_one_sided"
    out.loc[high_conf_pass, "insertion_call_tier"] = "high_conf_two_sided"

    # Separate sample-status label from structural confidence.
    # We allow +/-1 read slack to avoid overcalling "shared" or "germline-only"
    # from minor depth/mapping fluctuations.
    out["sample_status_label"] = "low_support"
    out.loc[(tumor_mei_reads >= 2) & (normal_mei_reads <= 1), "sample_status_label"] = "somatic_only"
    out.loc[(tumor_mei_reads >= 2) & (normal_mei_reads >= 2), "sample_status_label"] = "shared"
    out.loc[(tumor_mei_reads <= 1) & (normal_mei_reads >= 2), "sample_status_label"] = "germline_only"

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


def _annotate_bam_depth_for_consistent_loci(
    candidates: pd.DataFrame,
    tumor_bam_path: Path,
    normal_bam_path: Path,
) -> pd.DataFrame:
    out = candidates.copy()
    out["depth_filter_family_consistent"] = False
    out["depth_filter_pass"] = False
    out["tumor_local_bam_mean_depth"] = 0.0
    out["normal_local_bam_mean_depth"] = 0.0
    out["tumor_mei_support_per_100x_bam_depth"] = 0.0
    out["normal_mei_support_per_100x_bam_depth"] = 0.0
    out["mei_support_per_100x_bam_depth_delta"] = 0.0
    out["mei_support_per_100x_bam_depth_ratio"] = 1.0

    if out.empty:
        return out

    family_consistent = _consistent_family_mask(out)
    out["depth_filter_family_consistent"] = family_consistent
    depth_mask = (out.get("junk_flag_count", 999).fillna(999).astype(int) == 0) & family_consistent
    out["depth_filter_pass"] = depth_mask
    idxs = out.index[depth_mask].tolist()
    if not idxs:
        return out

    with pysam.AlignmentFile(str(tumor_bam_path), "rb") as tumor_bam, pysam.AlignmentFile(
        str(normal_bam_path), "rb"
    ) as normal_bam:
        for idx in idxs:
            row = out.loc[idx]
            chrom = str(row["chrom"])
            start = int(row["window_start"])
            end = int(row["window_end"])
            t_depth = _mean_depth_for_interval(tumor_bam, chrom=chrom, start_1based=start, end_1based=end)
            n_depth = _mean_depth_for_interval(normal_bam, chrom=chrom, start_1based=start, end_1based=end)
            out.at[idx, "tumor_local_bam_mean_depth"] = t_depth
            out.at[idx, "normal_local_bam_mean_depth"] = n_depth

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
    g1k_mei_bed: Path | None,
    g1k_mei_vcf: Path | None,
    split_padding_bp: int,
    dpe_padding_min_bp: int,
    dpe_padding_max_bp: int,
    dpe_padding_tlen_factor: float,
) -> pd.DataFrame:
    if g1k_mei_bed is not None and g1k_mei_vcf is not None:
        raise ValueError("Provide only one of g1k_mei_bed or g1k_mei_vcf.")
    if g1k_mei_bed is None and g1k_mei_vcf is None:
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
    best_hits: dict[int, dict[str, object]] = {}

    with tempfile.TemporaryDirectory(prefix="rtm_g1k_mei_") as tmpdir:
        tmp = Path(tmpdir)
        source_bed = g1k_mei_bed
        if g1k_mei_vcf is not None:
            source_bed = tmp / "g1k_mei_from_vcf.bed"
            kept = _build_g1k_mei_bed_from_vcf(g1k_mei_vcf, source_bed)
            print(f"[mei-annotate] parsed g1k MEI VCF records kept={kept} path={g1k_mei_vcf}")
        assert source_bed is not None

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


def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    reference_fasta: Path | None = None,
    tumor_bam_path: Path | None = None,
    normal_bam_path: Path | None = None,
    g1k_mei_bed: Path | None = None,
    g1k_mei_vcf: Path | None = None,
    g1k_split_padding_bp: int = 200,
    g1k_dpe_padding_min_bp: int = 200,
    g1k_dpe_padding_max_bp: int = 200,
    g1k_dpe_padding_tlen_factor: float = 0.0,
    progress_every: int = 20000,
) -> Path:
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
    if not disc_t.empty:
        candidate = candidate.merge(disc_t, on=["chrom", "window_start", "window_end"], how="left")
    if not disc_n.empty:
        candidate = candidate.merge(disc_n, on=["chrom", "window_start", "window_end"], how="left")
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
    if g1k_mei_bed is not None or g1k_mei_vcf is not None:
        candidate = _annotate_g1k_mei_overlap(
            candidate,
            g1k_mei_bed=g1k_mei_bed,
            g1k_mei_vcf=g1k_mei_vcf,
            split_padding_bp=g1k_split_padding_bp,
            dpe_padding_min_bp=g1k_dpe_padding_min_bp,
            dpe_padding_max_bp=g1k_dpe_padding_max_bp,
            dpe_padding_tlen_factor=g1k_dpe_padding_tlen_factor,
        )
        print("[mei-annotate] added 1000G/MELT polymorphism overlap fields")
    if tumor_bam_path is not None and normal_bam_path is not None:
        candidate = _annotate_bam_depth_for_consistent_loci(
            candidate,
            tumor_bam_path=tumor_bam_path,
            normal_bam_path=normal_bam_path,
        )
        print("[mei-annotate] added BAM-depth normalization for family-consistent, junk-clean loci")

    candidate = candidate.sort_values(
        [
            "passes_insertion_model",
            "insertion_model_score",
            "mei_score_enrichment_ratio",
            "coherence_score",
            "tumor_mei_supported_reads",
            "enrichment_ratio",
        ],
        ascending=[False, False, False, False, False, False],
        kind="mergesort",
    )

    candidate.to_csv(out_path, sep="\t", index=False)
    candidate.to_parquet(out_path.with_suffix(".parquet"), index=False)
    print(f"[mei-annotate] wrote {len(candidate)} rows to {out_path}")
    return out_path
