from __future__ import annotations

import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


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
            ]
        )

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
    agg = (
        side_df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{sample_prefix}_{side}_mei_supported_reads": ("read_name", "nunique"),
                f"{sample_prefix}_{side}_mei_score_sum": ("mei_score", "sum"),
                f"{sample_prefix}_{side}_mei_start": ("target_start", "min"),
                f"{sample_prefix}_{side}_mei_end": ("target_end", "max"),
            }
        )
        .merge(family_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_family"]], on=["chrom", "window_start", "window_end"], how="left")
        .merge(subfamily_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_subfamily"]], on=["chrom", "window_start", "window_end"], how="left")
        .merge(strand_top[["chrom", "window_start", "window_end", f"{sample_prefix}_{side}_mei_strand"]], on=["chrom", "window_start", "window_end"], how="left")
    )
    return agg


def _infer_tumor_insertion_metrics(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    for col in ["tumor_L_mei_start", "tumor_R_mei_start", "tumor_L_mei_end", "tumor_R_mei_end"]:
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
    return out


def annotate_candidate_loci_with_mei(
    evidence_dir: Path,
    candidate_loci_path: Path,
    mei_fasta: Path,
    out_path: Path,
    progress_every: int = 20000,
) -> Path:
    candidate = pd.read_csv(candidate_loci_path, sep="\t")
    window_size = int(candidate["window_end"].iloc[0] - candidate["window_start"].iloc[0] + 1) if not candidate.empty else 200

    split_tumor = _load_table(evidence_dir, "split_evidence", "tumor")
    split_normal = _load_table(evidence_dir, "split_evidence", "normal")
    for df in (split_tumor, split_normal):
        df["window_start"] = ((df["pos"].astype(int) - 1) // window_size) * window_size + 1
        df["window_end"] = df["window_start"] + window_size - 1

    tumor_hits, tumor_summary = _align_clips_with_minimap2(split_tumor, mei_fasta, sample="tumor")
    normal_hits, normal_summary = _align_clips_with_minimap2(split_normal, mei_fasta, sample="normal")

    print(
        f"[mei-annotate] tumor clips={tumor_summary.clip_count} hits={tumor_summary.paf_hits}; "
        f"normal clips={normal_summary.clip_count} hits={normal_summary.paf_hits}"
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

    for col in candidate.columns:
        if re.search(r"_mei_supported_reads$|_mei_start$|_mei_end$", col):
            candidate[col] = candidate[col].fillna(0).astype(int)
        if col.endswith("_mei_score_sum"):
            candidate[col] = candidate[col].fillna(0.0).astype(float)
        if re.search(r"_mei_family$|_mei_subfamily$|_mei_strand$", col):
            candidate[col] = candidate[col].fillna("")

    candidate["tumor_mei_score_sum"] = candidate.get("tumor_L_mei_score_sum", 0.0) + candidate.get(
        "tumor_R_mei_score_sum", 0.0
    )
    candidate["normal_mei_score_sum"] = candidate.get("normal_L_mei_score_sum", 0.0) + candidate.get(
        "normal_R_mei_score_sum", 0.0
    )
    candidate["tumor_mei_supported_reads"] = candidate.get("tumor_L_mei_supported_reads", 0) + candidate.get(
        "tumor_R_mei_supported_reads", 0
    )
    candidate["normal_mei_supported_reads"] = candidate.get("normal_L_mei_supported_reads", 0) + candidate.get(
        "normal_R_mei_supported_reads", 0
    )
    candidate["mei_score_enrichment_ratio"] = (candidate["tumor_mei_score_sum"] + 0.1) / (
        candidate["normal_mei_score_sum"] + 0.1
    )

    candidate = _infer_tumor_insertion_metrics(candidate)

    candidate = candidate.sort_values(
        ["mei_score_enrichment_ratio", "tumor_mei_supported_reads", "enrichment_ratio"],
        ascending=[False, False, False],
        kind="mergesort",
    )

    candidate.to_csv(out_path, sep="\t", index=False)
    candidate.to_parquet(out_path.with_suffix(".parquet"), index=False)
    print(f"[mei-annotate] wrote {len(candidate)} rows to {out_path}")
    return out_path
