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
    left_mode = out["tumor_L_mei_breakpoint_mode"].astype(int)
    right_mode = out["tumor_R_mei_breakpoint_mode"].astype(int)
    tsd_len = (right_mode - left_mode + 1).where((left_mode > 0) & (right_mode > 0), 0)
    out["tsd_len_estimate"] = tsd_len.where((tsd_len >= 2) & (tsd_len <= 30), 0).astype(int)
    out["tsd_detected"] = out["tsd_len_estimate"] > 0

    def _breakpoint_pos(row: pd.Series) -> int:
        l = int(row.get("tumor_L_mei_breakpoint_mode", 0))
        r = int(row.get("tumor_R_mei_breakpoint_mode", 0))
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
                    start0 = int(row.tumor_L_mei_breakpoint_mode) - 1
                    end0 = int(row.tumor_R_mei_breakpoint_mode)
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

    tumor_mei_reads = out.get("tumor_mei_supported_reads", 0).astype(float)
    normal_mei_reads = out.get("normal_mei_supported_reads", 0).astype(float)
    total_rows = out.get("tumor_total_rows", 0).astype(float).replace(0, 1.0)
    mei_enrichment = out.get("mei_score_enrichment_ratio", 0.0).astype(float)
    mei_enrichment_scaled = (mei_enrichment / (mei_enrichment + 1.0)).clip(lower=0.0, upper=1.0)
    mei_read_fraction = (tumor_mei_reads / total_rows).clip(lower=0.0, upper=1.0)

    base_score = (
        0.24 * out.get("tumor_subfamily_purity_weighted", 0.0).astype(float).fillna(0.0)
        + 0.18 * out.get("tumor_breakpoint_mode_fraction_weighted", 0.0).astype(float).fillna(0.0)
        + 0.20 * mei_enrichment_scaled.fillna(0.0)
        + 0.12 * mei_read_fraction.fillna(0.0)
        + 0.15 * out["tumor_family_agreement"].astype(float)
        # Subfamily agreement is informative but noisy with short-read clips:
        # keep as a small positive boost only (no explicit penalty on mismatch).
        + 0.03 * out["tumor_subfamily_agreement"].astype(float)
        + 0.08 * out["tumor_strand_agreement"].astype(float)
    )

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
    out["tumor_two_sided_like_support"] = out["tumor_two_sided_strong_support"] | (
        out["tumor_one_sided_split_support"] & out["tumor_discordant_mei_strong_support"]
    ) | out["tumor_discordant_mei_consistent_support"]
    out["tumor_side_breakpoint_consistency"] = left_mode_frac.combine(right_mode_frac, min)
    out["tumor_side_subfamily_purity"] = left_purity.combine(right_purity, min)
    out["tumor_two_sided_family_consistent"] = out["tumor_two_sided_support"] & (out["tumor_family_agreement"] == 1)
    out["tumor_two_sided_subfamily_consistent"] = out["tumor_two_sided_support"] & (
        out["tumor_subfamily_agreement"] == 1
    )

    # Structural confidence gates (sample-status agnostic).
    high_conf_pass = (
        (out["insertion_model_score"] >= 0.60)
        & (tumor_mei_reads >= 4)
        & out["tumor_two_sided_like_support"]
        & out["tumor_two_sided_family_consistent"]
        & (out["tumor_side_breakpoint_consistency"] >= 0.60)
        & (out.get("coherence_score", 0.0).astype(float) >= 0.55)
    )
    provisional_one_sided = (
        (~high_conf_pass)
        & (out["insertion_model_score"] >= 0.55)
        & (tumor_mei_reads >= 2)
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


def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    reference_fasta: Path | None = None,
    tumor_bam_path: Path | None = None,
    normal_bam_path: Path | None = None,
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
