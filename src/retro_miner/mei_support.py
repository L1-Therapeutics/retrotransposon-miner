from __future__ import annotations

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
            ]
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
    out["tsd_seq"] = ""
    if reference_fasta is not None:
        with pysam.FastaFile(str(reference_fasta)) as ref:
            seqs = []
            for row in out.itertuples(index=False):
                if int(row.tsd_len_estimate) <= 0:
                    seqs.append("")
                    continue
                chrom = str(row.chrom)
                start0 = int(row.tumor_L_mei_breakpoint_mode) - 1
                end0 = int(row.tumor_R_mei_breakpoint_mode)
                try:
                    seqs.append(ref.fetch(chrom, start0, end0).upper())
                except Exception:
                    seqs.append("")
            out["tsd_seq"] = seqs

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
        + 0.10 * out["tumor_family_agreement"].astype(float)
        + 0.08 * out["tumor_subfamily_agreement"].astype(float)
        + 0.08 * out["tumor_strand_agreement"].astype(float)
    )

    # Penalize patterns more consistent with deletion/complex SV hotspots than MEI insertion.
    large_insert_fraction = out.get("discordant_tumor_large_insert_fraction", 0.0).astype(float).fillna(0.0)
    interchrom_fraction = out.get("discordant_tumor_interchrom_fraction", 0.0).astype(float).fillna(0.0)
    deletion_complex_penalty = 0.0
    deletion_complex_penalty += (
        ((large_insert_fraction > 0.6) & (tumor_mei_reads < 4)).astype(float) * 0.18
    )
    deletion_complex_penalty += (
        ((interchrom_fraction > 0.6) & (tumor_mei_reads < 4)).astype(float) * 0.18
    )
    deletion_complex_penalty += ((normal_mei_reads >= tumor_mei_reads) & (tumor_mei_reads > 0)).astype(float) * 0.15

    score = (base_score - deletion_complex_penalty).clip(lower=0.0, upper=1.0)
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
    out["tumor_two_sided_like_support"] = out["tumor_two_sided_strong_support"] | (
        out["tumor_one_sided_split_support"] & out["tumor_discordant_mei_strong_support"]
    )
    out["tumor_side_breakpoint_consistency"] = left_mode_frac.combine(right_mode_frac, min)
    out["tumor_side_subfamily_purity"] = left_purity.combine(right_purity, min)
    out["tumor_two_sided_family_consistent"] = out["tumor_two_sided_support"] & (out["tumor_family_agreement"] == 1)
    out["tumor_two_sided_subfamily_consistent"] = out["tumor_two_sided_support"] & (
        out["tumor_subfamily_agreement"] == 1
    )

    high_conf_pass = (
        (out["insertion_model_score"] >= 0.60)
        & (tumor_mei_reads >= 4)
        & out["tumor_two_sided_like_support"]
        & out["tumor_two_sided_family_consistent"]
        & out["tumor_two_sided_subfamily_consistent"]
        & (out["tumor_side_breakpoint_consistency"] >= 0.60)
        & (out["tumor_side_subfamily_purity"] >= 0.60)
        & (mei_enrichment > 1.10)
        & (out.get("coherence_score", 0.0).astype(float) >= 0.55)
    )
    provisional_one_sided = (
        (~high_conf_pass)
        & (out["insertion_model_score"] >= 0.55)
        & (tumor_mei_reads >= 2)
        & (mei_enrichment > 1.10)
        & (out.get("coherence_score", 0.0).astype(float) >= 0.50)
    )

    out["passes_insertion_model"] = high_conf_pass
    out["passes_insertion_model_provisional"] = provisional_one_sided
    out["insertion_call_tier"] = "none"
    out.loc[provisional_one_sided, "insertion_call_tier"] = "provisional_one_sided"
    out.loc[high_conf_pass, "insertion_call_tier"] = "high_conf_two_sided"
    return out


def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    reference_fasta: Path | None = None,
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

    candidate = _infer_tumor_insertion_metrics(candidate, reference_fasta=reference_fasta)
    candidate = _compute_insertion_model_scores(candidate)

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
