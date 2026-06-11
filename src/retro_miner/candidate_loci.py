from __future__ import annotations

import gzip
from collections.abc import Iterable
from pathlib import Path
import subprocess
import tempfile

import pandas as pd


def _load_evidence_table(base_dir: Path, stem: str, sample: str) -> pd.DataFrame:
    parquet_path = base_dir / f"{stem}.{sample}.parquet"
    tsv_path = base_dir / f"{stem}.{sample}.tsv"
    if parquet_path.exists():
        return pd.read_parquet(parquet_path)
    if tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    raise FileNotFoundError(
        f"Missing evidence table for stem={stem} sample={sample}; "
        f"looked for {parquet_path} and {tsv_path}"
    )


def _read_passing_counts(summary_path: Path) -> dict[str, int]:
    summary = pd.read_csv(summary_path, sep="\t")
    counts: dict[str, int] = {}
    for _, row in summary.iterrows():
        counts[str(row["sample"])] = int(row["passing_reads"])
    return counts


def _windowize(df: pd.DataFrame, window_size: int) -> pd.DataFrame:
    out = df.copy()
    out["window_start"] = ((out["pos"].astype(int) - 1) // window_size) * window_size + 1
    out["window_end"] = out["window_start"] + window_size - 1
    return out


def _aggregate_evidence(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{prefix}_rows",
                f"{prefix}_unique_reads",
            ]
        )

    grouped = (
        df.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{prefix}_rows": ("read_name", "size"),
                f"{prefix}_unique_reads": ("read_name", "nunique"),
            }
        )
        .sort_values(["chrom", "window_start"], kind="mergesort")
    )
    return grouped


def _safe_cpm(count_series: pd.Series, denominator: int) -> pd.Series:
    if denominator <= 0:
        return pd.Series(0.0, index=count_series.index)
    return count_series.astype(float) * 1_000_000.0 / float(denominator)


def _open_textmaybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_interval_parts(parts: list[str]) -> tuple[str, int, int] | None:
    if len(parts) >= 3 and parts[0].startswith("chr"):
        chrom = parts[0]
        start_idx = 1
        end_idx = 2
    elif len(parts) >= 4 and parts[1].startswith("chr"):
        # UCSC table format with leading `bin`.
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
    return (chrom, start0, end0)


def _normalize_track_to_bed(
    input_path: Path,
    output_path: Path,
    mappability_threshold: float | None = None,
) -> None:
    with _open_textmaybe_gz(input_path) as in_handle, output_path.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            parsed = _parse_interval_parts(parts)
            if parsed is None:
                continue
            chrom, start0, end0 = parsed

            if mappability_threshold is not None:
                if len(parts) < 4:
                    continue
                try:
                    score = float(parts[3] if parts[0].startswith("chr") else parts[4])
                except ValueError:
                    continue
                if score >= mappability_threshold:
                    continue

            out_handle.write(f"{chrom}\t{start0}\t{end0}\n")


def _write_candidate_windows_bed(loci: pd.DataFrame, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for row in loci.itertuples(index=False):
            # BED uses 0-based start, half-open end.
            start0 = int(row.window_start) - 1
            end0 = int(row.window_end)
            handle.write(f"{row.chrom}\t{start0}\t{end0}\t{row.row_id}\n")


def _get_overlapping_row_ids(a_bed: Path, b_bed: Path) -> set[int]:
    cmd = ["bedtools", "intersect", "-a", str(a_bed), "-b", str(b_bed), "-wa", "-u"]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    row_ids: set[int] = set()
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            row_ids.add(int(parts[3]))
        except ValueError:
            continue
    return row_ids


def _get_overlapping_row_ids_with_fraction(a_bed: Path, b_bed: Path, min_fraction: float) -> set[int]:
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(a_bed),
        "-b",
        str(b_bed),
        "-wa",
        "-u",
        "-f",
        str(min_fraction),
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    row_ids: set[int] = set()
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            row_ids.add(int(parts[3]))
        except ValueError:
            continue
    return row_ids


def _annotate_junk_flags(
    loci: pd.DataFrame,
    discordant_tumor: pd.DataFrame,
    discordant_normal: pd.DataFrame,
    segdup_bed: Path | None,
    segdup_min_fraction: float,
    low_mappability_bedgraph: Path | None,
    low_mappability_threshold: float,
    low_mappability_min_fraction: float,
    giab_highconf_bed: Path | None,
    gap_bed: Path | None,
    gap_min_fraction: float,
    encode_blacklist_bed: Path | None,
    encode_blacklist_min_fraction: float,
) -> pd.DataFrame:
    out = loci.copy()
    out = out.reset_index(drop=True)
    out["row_id"] = out.index.astype(int)
    out["flag_segdup"] = False
    out["flag_low_mappability"] = False
    out["flag_outside_giab_highconf"] = False
    out["flag_gap_region"] = False
    out["flag_encode_blacklist"] = False
    out["flag_mate_in_segdup"] = False
    out["flag_mate_low_mappability"] = False
    out["flag_mate_in_gap"] = False
    out["flag_mate_in_blacklist"] = False
    out["junk_flag_count"] = 0
    out["mate_junk_flag_count"] = 0

    with tempfile.TemporaryDirectory(prefix="rtm_candidate_loci_") as tmpdir:
        tmp = Path(tmpdir)
        candidate_bed = tmp / "candidate_windows.bed"
        _write_candidate_windows_bed(out, candidate_bed)

        key_cols = ["chrom", "window_start", "window_end"]
        candidate_key = out.loc[:, key_cols + ["row_id"]]
        discordant_all = pd.concat([discordant_tumor, discordant_normal], ignore_index=True)
        if discordant_all.empty:
            mate_with_id = pd.DataFrame(columns=key_cols + ["row_id", "mate_chrom", "mate_pos"])
        else:
            mate_with_id = discordant_all.merge(candidate_key, on=key_cols, how="inner")
            mate_with_id = mate_with_id.loc[
                (mate_with_id["mate_chrom"].astype(str) != "*") & (mate_with_id["mate_pos"].astype(int) > 0),
                ["row_id", "mate_chrom", "mate_pos"],
            ].drop_duplicates()

        mate_bed = tmp / "candidate_mate_positions.bed"
        with mate_bed.open("w", encoding="utf-8") as handle:
            for row in mate_with_id.itertuples(index=False):
                start0 = int(row.mate_pos) - 1
                end0 = int(row.mate_pos)
                handle.write(f"{row.mate_chrom}\t{start0}\t{end0}\t{row.row_id}\n")

        if segdup_bed is not None and segdup_bed.exists():
            segdup_norm = tmp / "segdup.norm.bed"
            _normalize_track_to_bed(segdup_bed, segdup_norm)
            segdup_hits = _get_overlapping_row_ids_with_fraction(candidate_bed, segdup_norm, segdup_min_fraction)
            out["flag_segdup"] = out["row_id"].isin(segdup_hits)
            if not mate_with_id.empty:
                segdup_mate_hits = _get_overlapping_row_ids_with_fraction(mate_bed, segdup_norm, segdup_min_fraction)
                out["flag_mate_in_segdup"] = out["row_id"].isin(segdup_mate_hits)

        if low_mappability_bedgraph is not None and low_mappability_bedgraph.exists():
            mappability_norm = tmp / "mappability_low.norm.bed"
            _normalize_track_to_bed(
                low_mappability_bedgraph,
                mappability_norm,
                mappability_threshold=low_mappability_threshold,
            )
            low_map_hits = _get_overlapping_row_ids_with_fraction(
                candidate_bed, mappability_norm, low_mappability_min_fraction
            )
            out["flag_low_mappability"] = out["row_id"].isin(low_map_hits)
            if not mate_with_id.empty:
                low_map_mate_hits = _get_overlapping_row_ids_with_fraction(
                    mate_bed, mappability_norm, low_mappability_min_fraction
                )
                out["flag_mate_low_mappability"] = out["row_id"].isin(low_map_mate_hits)

        if giab_highconf_bed is not None and giab_highconf_bed.exists():
            giab_norm = tmp / "giab_highconf.norm.bed"
            _normalize_track_to_bed(giab_highconf_bed, giab_norm)
            giab_hits = _get_overlapping_row_ids(candidate_bed, giab_norm)
            out["flag_outside_giab_highconf"] = ~out["row_id"].isin(giab_hits)

        if gap_bed is not None and gap_bed.exists():
            gap_norm = tmp / "gap.norm.bed"
            _normalize_track_to_bed(gap_bed, gap_norm)
            gap_hits = _get_overlapping_row_ids_with_fraction(candidate_bed, gap_norm, gap_min_fraction)
            out["flag_gap_region"] = out["row_id"].isin(gap_hits)
            if not mate_with_id.empty:
                gap_mate_hits = _get_overlapping_row_ids_with_fraction(mate_bed, gap_norm, gap_min_fraction)
                out["flag_mate_in_gap"] = out["row_id"].isin(gap_mate_hits)

        if encode_blacklist_bed is not None and encode_blacklist_bed.exists():
            blacklist_norm = tmp / "blacklist.norm.bed"
            _normalize_track_to_bed(encode_blacklist_bed, blacklist_norm)
            blacklist_hits = _get_overlapping_row_ids_with_fraction(
                candidate_bed, blacklist_norm, encode_blacklist_min_fraction
            )
            out["flag_encode_blacklist"] = out["row_id"].isin(blacklist_hits)
            if not mate_with_id.empty:
                blacklist_mate_hits = _get_overlapping_row_ids_with_fraction(
                    mate_bed, blacklist_norm, encode_blacklist_min_fraction
                )
                out["flag_mate_in_blacklist"] = out["row_id"].isin(blacklist_mate_hits)

    flag_cols: Iterable[str] = (
        "flag_segdup",
        "flag_low_mappability",
        "flag_outside_giab_highconf",
        "flag_gap_region",
        "flag_encode_blacklist",
    )
    out["junk_flag_count"] = out.loc[:, list(flag_cols)].sum(axis=1)
    mate_flag_cols: Iterable[str] = (
        "flag_mate_in_segdup",
        "flag_mate_low_mappability",
        "flag_mate_in_gap",
        "flag_mate_in_blacklist",
    )
    out["mate_junk_flag_count"] = out.loc[:, list(mate_flag_cols)].sum(axis=1)
    out = out.drop(columns=["row_id"])
    return out


def build_candidate_loci(
    evidence_dir: Path,
    outdir: Path,
    window_size: int = 200,
    pseudocount: float = 1.0,
    segdup_bed: Path | None = None,
    segdup_min_fraction: float = 0.1,
    low_mappability_bedgraph: Path | None = None,
    low_mappability_threshold: float = 0.5,
    low_mappability_min_fraction: float = 0.5,
    giab_highconf_bed: Path | None = None,
    gap_bed: Path | None = None,
    gap_min_fraction: float = 0.1,
    encode_blacklist_bed: Path | None = None,
    encode_blacklist_min_fraction: float = 0.1,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    passing_counts = _read_passing_counts(evidence_dir / "split_evidence.summary.tsv")

    split_tumor = _windowize(_load_evidence_table(evidence_dir, "split_evidence", "tumor"), window_size)
    split_normal = _windowize(_load_evidence_table(evidence_dir, "split_evidence", "normal"), window_size)
    discordant_tumor = _windowize(_load_evidence_table(evidence_dir, "discordant_evidence", "tumor"), window_size)
    discordant_normal = _windowize(_load_evidence_table(evidence_dir, "discordant_evidence", "normal"), window_size)

    split_t = _aggregate_evidence(split_tumor, "split_tumor")
    split_n = _aggregate_evidence(split_normal, "split_normal")
    disc_t = _aggregate_evidence(discordant_tumor, "discordant_tumor")
    disc_n = _aggregate_evidence(discordant_normal, "discordant_normal")

    merged = split_t.merge(split_n, on=["chrom", "window_start", "window_end"], how="outer")
    merged = merged.merge(disc_t, on=["chrom", "window_start", "window_end"], how="outer")
    merged = merged.merge(disc_n, on=["chrom", "window_start", "window_end"], how="outer")
    merged = merged.fillna(0)

    int_cols = [c for c in merged.columns if c.endswith("_rows") or c.endswith("_unique_reads")]
    for col in int_cols:
        merged[col] = merged[col].astype(int)

    merged["tumor_total_rows"] = merged["split_tumor_rows"] + merged["discordant_tumor_rows"]
    merged["normal_total_rows"] = merged["split_normal_rows"] + merged["discordant_normal_rows"]

    tumor_den = passing_counts.get("tumor", 0)
    normal_den = passing_counts.get("normal", 0)
    merged["tumor_total_cpm"] = _safe_cpm(merged["tumor_total_rows"], tumor_den)
    merged["normal_total_cpm"] = _safe_cpm(merged["normal_total_rows"], normal_den)

    merged["enrichment_ratio"] = (merged["tumor_total_cpm"] + pseudocount) / (
        merged["normal_total_cpm"] + pseudocount
    )
    merged["delta_cpm"] = merged["tumor_total_cpm"] - merged["normal_total_cpm"]

    merged = _annotate_junk_flags(
        loci=merged,
        discordant_tumor=discordant_tumor,
        discordant_normal=discordant_normal,
        segdup_bed=segdup_bed,
        segdup_min_fraction=segdup_min_fraction,
        low_mappability_bedgraph=low_mappability_bedgraph,
        low_mappability_threshold=low_mappability_threshold,
        low_mappability_min_fraction=low_mappability_min_fraction,
        giab_highconf_bed=giab_highconf_bed,
        gap_bed=gap_bed,
        gap_min_fraction=gap_min_fraction,
        encode_blacklist_bed=encode_blacklist_bed,
        encode_blacklist_min_fraction=encode_blacklist_min_fraction,
    )

    merged = merged.sort_values(["enrichment_ratio", "tumor_total_rows"], ascending=[False, False], kind="mergesort")

    out_tsv = outdir / "candidate_loci.tsv"
    out_parquet = outdir / "candidate_loci.parquet"
    merged.to_csv(out_tsv, sep="\t", index=False)
    merged.to_parquet(out_parquet, index=False)
    return out_tsv
