from __future__ import annotations

import gzip
import time
from collections.abc import Iterable
from pathlib import Path
import subprocess
import tempfile

import pandas as pd
from intervaltree import IntervalTree


_RUN_T0: float | None = None


def _progress(msg: str) -> None:
    if _RUN_T0 is None:
        print(f"[candidate-loci] {msg}", flush=True)
    else:
        print(f"[candidate-loci] +{(time.monotonic() - _RUN_T0):.1f}s {msg}", flush=True)


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


def _cluster_sorted_positions(positions: list[int], max_gap_bp: int) -> list[list[int]]:
    if not positions:
        return []
    if max_gap_bp < 0:
        max_gap_bp = 0
    clusters: list[list[int]] = [[int(positions[0])]]
    for pos in positions[1:]:
        p = int(pos)
        if p - clusters[-1][-1] <= max_gap_bp:
            clusters[-1].append(p)
        else:
            clusters.append([p])
    return clusters


def _split_cluster_positions(
    positions: list[int],
    valley_gap_bp: int,
    max_locus_span_bp: int,
) -> list[list[int]]:
    if not positions:
        return []
    if valley_gap_bp < 0:
        valley_gap_bp = 0
    if max_locus_span_bp <= 0:
        max_locus_span_bp = 1

    sorted_pos = sorted(int(p) for p in positions)
    out: list[list[int]] = [[sorted_pos[0]]]
    current_start = sorted_pos[0]
    for pos in sorted_pos[1:]:
        prev = out[-1][-1]
        span_if_added = pos - current_start + 1
        if (pos - prev > valley_gap_bp) or (span_if_added > max_locus_span_bp):
            out.append([pos])
            current_start = pos
        else:
            out[-1].append(pos)
    return out


def _distance_to_closed_interval(pos: int, start: int, end: int) -> int:
    if pos < start:
        return start - pos
    if pos > end:
        return pos - end
    return 0


def _build_loci_from_evidence(
    split_tumor: pd.DataFrame,
    split_normal: pd.DataFrame,
    discordant_tumor: pd.DataFrame,
    discordant_normal: pd.DataFrame,
    split_cluster_bp: int,
    discordant_cluster_bp: int,
    max_locus_span_bp: int,
) -> pd.DataFrame:
    _progress(
        "clustering evidence "
        f"(split_cluster_bp={split_cluster_bp}, discordant_cluster_bp={discordant_cluster_bp}, "
        f"max_locus_span_bp={max_locus_span_bp})"
    )
    split_all = pd.concat([split_tumor, split_normal], ignore_index=True)
    discordant_all = pd.concat([discordant_tumor, discordant_normal], ignore_index=True)
    if split_all.empty and discordant_all.empty:
        _progress("no split/discordant evidence found; returning empty loci")
        return pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    loci_by_chrom: dict[str, list[dict[str, object]]] = {}

    if not split_all.empty:
        split_pos = split_all.loc[:, ["chrom", "pos"]].copy()
        split_pos["pos"] = split_pos["pos"].astype(int)
        for chrom, chrom_df in split_pos.groupby("chrom", sort=False):
            positions = sorted(chrom_df["pos"].tolist())
            for cluster in _cluster_sorted_positions(positions, max_gap_bp=int(split_cluster_bp)):
                if not cluster:
                    continue
                loci_by_chrom.setdefault(str(chrom), []).append(
                    {
                        "positions": list(cluster),
                        "min_pos": int(cluster[0]),
                        "max_pos": int(cluster[-1]),
                        "seeded_by_split": True,
                    }
                )

    unassigned_discordant: dict[str, list[int]] = {}
    if not discordant_all.empty:
        disc_pos = discordant_all.loc[:, ["chrom", "pos"]].copy()
        disc_pos["pos"] = disc_pos["pos"].astype(int)
        for chrom, chrom_df in disc_pos.groupby("chrom", sort=False):
            chrom_key = str(chrom)
            seeds = loci_by_chrom.get(chrom_key, [])
            seed_tree: IntervalTree | None = None
            if seeds:
                seed_tree = IntervalTree()
                radius = int(discordant_cluster_bp)
                for idx, seed in enumerate(seeds):
                    seed_start = int(seed["min_pos"])
                    seed_end = int(seed["max_pos"])
                    seed_tree.addi(max(1, seed_start - radius), seed_end + radius + 1, idx)
            for pos in sorted(chrom_df["pos"].tolist()):
                if not seeds:
                    unassigned_discordant.setdefault(chrom_key, []).append(int(pos))
                    continue

                candidate_seed_idxs = list(seed_tree.at(int(pos))) if seed_tree is not None else []
                if not candidate_seed_idxs:
                    unassigned_discordant.setdefault(chrom_key, []).append(int(pos))
                    continue

                best_idx = min(
                    (int(iv.data) for iv in candidate_seed_idxs),
                    key=lambda idx: _distance_to_closed_interval(
                        int(pos),
                        int(seeds[idx]["min_pos"]),
                        int(seeds[idx]["max_pos"]),
                    ),
                )
                best_dist = _distance_to_closed_interval(
                    int(pos),
                    int(seeds[best_idx]["min_pos"]),
                    int(seeds[best_idx]["max_pos"]),
                )
                if best_dist <= int(discordant_cluster_bp):
                    chosen = seeds[best_idx]
                    chosen["positions"].append(int(pos))
                else:
                    unassigned_discordant.setdefault(chrom_key, []).append(int(pos))

    for chrom, positions in unassigned_discordant.items():
        for cluster in _cluster_sorted_positions(sorted(positions), max_gap_bp=int(discordant_cluster_bp)):
            if not cluster:
                continue
            loci_by_chrom.setdefault(chrom, []).append(
                {
                    "positions": list(cluster),
                    "min_pos": int(cluster[0]),
                    "max_pos": int(cluster[-1]),
                    "seeded_by_split": False,
                }
            )

    rows: list[dict[str, int | str]] = []
    valley_gap_bp = max(int(split_cluster_bp), min(int(discordant_cluster_bp), int(max_locus_span_bp)))
    for chrom in sorted(loci_by_chrom):
        for locus in loci_by_chrom[chrom]:
            segments = _split_cluster_positions(
                positions=[int(p) for p in locus["positions"]],
                valley_gap_bp=valley_gap_bp,
                max_locus_span_bp=int(max_locus_span_bp),
            )
            flank = int(split_cluster_bp) if bool(locus["seeded_by_split"]) else int(discordant_cluster_bp)
            for seg in segments:
                if not seg:
                    continue
                seg_start = int(min(seg))
                seg_end = int(max(seg))
                rows.append(
                    {
                        "chrom": chrom,
                        "window_start": max(1, seg_start - flank),
                        "window_end": seg_end + flank,
                    }
                )

    if not rows:
        _progress("no loci generated after clustering")
        return pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    out = pd.DataFrame(rows).drop_duplicates()
    out = out.sort_values(["chrom", "window_start", "window_end"], kind="mergesort").reset_index(drop=True)
    _progress(f"built {len(out)} candidate loci from {len(split_all)} split and {len(discordant_all)} discordant rows")
    return out


def _assign_rows_to_loci(df: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    if df.empty or loci.empty:
        return pd.DataFrame(columns=list(df.columns) + ["window_start", "window_end"])

    trees: dict[str, IntervalTree] = {}
    for row in loci.itertuples(index=False):
        chrom = str(row.chrom)
        tree = trees.setdefault(chrom, IntervalTree())
        tree.addi(int(row.window_start), int(row.window_end) + 1, (int(row.window_start), int(row.window_end)))

    assigned_rows: list[dict[str, object]] = []
    for row in df.itertuples(index=False):
        chrom = str(row.chrom)
        pos = int(row.pos)
        tree = trees.get(chrom)
        if tree is None:
            continue
        overlaps = list(tree.at(pos))
        if not overlaps:
            continue
        # If a position lands in multiple merged loci, choose the tightest interval.
        best = min(overlaps, key=lambda iv: (iv.end - iv.begin, abs(((iv.begin + iv.end) // 2) - pos)))
        locus_start, locus_end = best.data
        as_dict = row._asdict()
        as_dict["window_start"] = locus_start
        as_dict["window_end"] = locus_end
        assigned_rows.append(as_dict)
    return pd.DataFrame(assigned_rows)


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


def _aggregate_split_metrics(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{prefix}_rows",
                f"{prefix}_unique_reads",
                f"{prefix}_mapq_mean",
                f"{prefix}_mapq_min",
                f"{prefix}_clip_len_mean",
                f"{prefix}_clip_len_max",
                f"{prefix}_sa_fraction",
                f"{prefix}_poly_tail_rescued_rows",
                f"{prefix}_poly_tail_rescued_unique_reads",
                f"{prefix}_poly_tail_at_fraction_mean",
                f"{prefix}_poly_tail_at_run_max",
            ]
        )
    tmp = df.copy()
    if "nm" not in tmp.columns:
        tmp["nm"] = -1
    if "poly_tail_rescued" not in tmp.columns:
        tmp["poly_tail_rescued"] = False
    if "clip_poly_at_fraction" not in tmp.columns:
        tmp["clip_poly_at_fraction"] = 0.0
    if "clip_poly_at_run" not in tmp.columns:
        tmp["clip_poly_at_run"] = 0
    tmp["has_sa_int"] = tmp["has_sa"].astype(bool).astype(int)
    tmp["poly_tail_rescued_int"] = tmp["poly_tail_rescued"].astype(bool).astype(int)
    grouped = (
        tmp.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{prefix}_rows": ("read_name", "size"),
                f"{prefix}_unique_reads": ("read_name", "nunique"),
                f"{prefix}_mapq_mean": ("mapq", "mean"),
                f"{prefix}_mapq_min": ("mapq", "min"),
                f"{prefix}_clip_len_mean": ("clip_len", "mean"),
                f"{prefix}_clip_len_max": ("clip_len", "max"),
                f"{prefix}_sa_fraction": ("has_sa_int", "mean"),
                f"{prefix}_nm_mean": ("nm", "mean"),
                f"{prefix}_nm_max": ("nm", "max"),
                f"{prefix}_poly_tail_rescued_rows": ("poly_tail_rescued_int", "sum"),
                f"{prefix}_poly_tail_at_fraction_mean": ("clip_poly_at_fraction", "mean"),
                f"{prefix}_poly_tail_at_run_max": ("clip_poly_at_run", "max"),
            }
        )
        .sort_values(["chrom", "window_start"], kind="mergesort")
    )
    rescued_unique = (
        tmp.loc[tmp["poly_tail_rescued_int"] == 1]
        .groupby(["chrom", "window_start", "window_end"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": f"{prefix}_poly_tail_rescued_unique_reads"})
    )
    grouped = grouped.merge(rescued_unique, on=["chrom", "window_start", "window_end"], how="left")
    grouped[f"{prefix}_poly_tail_rescued_unique_reads"] = (
        grouped[f"{prefix}_poly_tail_rescued_unique_reads"].fillna(0).astype(int)
    )
    return grouped


def _aggregate_discordant_metrics(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                f"{prefix}_rows",
                f"{prefix}_unique_reads",
                f"{prefix}_mapq_mean",
                f"{prefix}_mapq_min",
                f"{prefix}_abs_tlen_mean",
                f"{prefix}_mate_unmapped_fraction",
                f"{prefix}_interchrom_fraction",
                f"{prefix}_large_insert_fraction",
                f"{prefix}_same_strand_fraction",
                f"{prefix}_improper_pair_fraction",
                f"{prefix}_poly_tail_rescued_rows",
                f"{prefix}_poly_tail_rescued_unique_reads",
                f"{prefix}_poly_tail_at_fraction_mean",
                f"{prefix}_poly_tail_at_run_max",
            ]
        )
    tmp = df.copy()
    if "nm" not in tmp.columns:
        tmp["nm"] = -1
    if "poly_tail_anchor_rescued" not in tmp.columns:
        tmp["poly_tail_anchor_rescued"] = False
    if "anchor_poly_at_fraction" not in tmp.columns:
        tmp["anchor_poly_at_fraction"] = 0.0
    if "anchor_poly_at_run" not in tmp.columns:
        tmp["anchor_poly_at_run"] = 0
    reason_col = tmp["discordant_reasons"].fillna("").astype(str)
    tmp["reason_mate_unmapped"] = reason_col.str.contains("mate_unmapped", regex=False).astype(int)
    tmp["reason_interchrom"] = reason_col.str.contains("interchrom", regex=False).astype(int)
    tmp["reason_large_insert"] = reason_col.str.contains("large_insert", regex=False).astype(int)
    tmp["reason_same_strand"] = reason_col.str.contains("same_strand", regex=False).astype(int)
    tmp["reason_improper_pair"] = reason_col.str.contains("improper_pair", regex=False).astype(int)
    tmp["poly_tail_anchor_rescued_int"] = tmp["poly_tail_anchor_rescued"].astype(bool).astype(int)
    tmp["abs_tlen"] = tmp["template_len"].abs()
    grouped = (
        tmp.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            **{
                f"{prefix}_rows": ("read_name", "size"),
                f"{prefix}_unique_reads": ("read_name", "nunique"),
                f"{prefix}_mapq_mean": ("mapq", "mean"),
                f"{prefix}_mapq_min": ("mapq", "min"),
                f"{prefix}_abs_tlen_mean": ("abs_tlen", "mean"),
                f"{prefix}_mate_unmapped_fraction": ("reason_mate_unmapped", "mean"),
                f"{prefix}_interchrom_fraction": ("reason_interchrom", "mean"),
                f"{prefix}_large_insert_fraction": ("reason_large_insert", "mean"),
                f"{prefix}_same_strand_fraction": ("reason_same_strand", "mean"),
                f"{prefix}_improper_pair_fraction": ("reason_improper_pair", "mean"),
                f"{prefix}_nm_mean": ("nm", "mean"),
                f"{prefix}_nm_max": ("nm", "max"),
                f"{prefix}_poly_tail_rescued_rows": ("poly_tail_anchor_rescued_int", "sum"),
                f"{prefix}_poly_tail_at_fraction_mean": ("anchor_poly_at_fraction", "mean"),
                f"{prefix}_poly_tail_at_run_max": ("anchor_poly_at_run", "max"),
            }
        )
        .sort_values(["chrom", "window_start"], kind="mergesort")
    )
    rescued_unique = (
        tmp.loc[tmp["poly_tail_anchor_rescued_int"] == 1]
        .groupby(["chrom", "window_start", "window_end"], as_index=False)["read_name"]
        .nunique()
        .rename(columns={"read_name": f"{prefix}_poly_tail_rescued_unique_reads"})
    )
    grouped = grouped.merge(rescued_unique, on=["chrom", "window_start", "window_end"], how="left")
    grouped[f"{prefix}_poly_tail_rescued_unique_reads"] = (
        grouped[f"{prefix}_poly_tail_rescued_unique_reads"].fillna(0).astype(int)
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
    _progress("starting junk-region annotation")
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
            _progress("annotating segdup overlaps")
            segdup_norm = tmp / "segdup.norm.bed"
            _normalize_track_to_bed(segdup_bed, segdup_norm)
            segdup_hits = _get_overlapping_row_ids_with_fraction(candidate_bed, segdup_norm, segdup_min_fraction)
            out["flag_segdup"] = out["row_id"].isin(segdup_hits)
            _progress(f"segdup hits={len(segdup_hits)}")
            if not mate_with_id.empty:
                segdup_mate_hits = _get_overlapping_row_ids_with_fraction(mate_bed, segdup_norm, segdup_min_fraction)
                out["flag_mate_in_segdup"] = out["row_id"].isin(segdup_mate_hits)
                _progress(f"segdup mate hits={len(segdup_mate_hits)}")

        if low_mappability_bedgraph is not None and low_mappability_bedgraph.exists():
            _progress("annotating low-mappability overlaps")
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
            _progress(f"low-mappability hits={len(low_map_hits)}")
            if not mate_with_id.empty:
                low_map_mate_hits = _get_overlapping_row_ids_with_fraction(
                    mate_bed, mappability_norm, low_mappability_min_fraction
                )
                out["flag_mate_low_mappability"] = out["row_id"].isin(low_map_mate_hits)
                _progress(f"low-mappability mate hits={len(low_map_mate_hits)}")

        if giab_highconf_bed is not None and giab_highconf_bed.exists():
            _progress("annotating GIAB high-confidence coverage")
            giab_norm = tmp / "giab_highconf.norm.bed"
            _normalize_track_to_bed(giab_highconf_bed, giab_norm)
            giab_hits = _get_overlapping_row_ids(candidate_bed, giab_norm)
            out["flag_outside_giab_highconf"] = ~out["row_id"].isin(giab_hits)
            _progress(f"inside GIAB high-confidence hits={len(giab_hits)}")

        if gap_bed is not None and gap_bed.exists():
            _progress("annotating gap overlaps")
            gap_norm = tmp / "gap.norm.bed"
            _normalize_track_to_bed(gap_bed, gap_norm)
            gap_hits = _get_overlapping_row_ids_with_fraction(candidate_bed, gap_norm, gap_min_fraction)
            out["flag_gap_region"] = out["row_id"].isin(gap_hits)
            _progress(f"gap hits={len(gap_hits)}")
            if not mate_with_id.empty:
                gap_mate_hits = _get_overlapping_row_ids_with_fraction(mate_bed, gap_norm, gap_min_fraction)
                out["flag_mate_in_gap"] = out["row_id"].isin(gap_mate_hits)
                _progress(f"gap mate hits={len(gap_mate_hits)}")

        if encode_blacklist_bed is not None and encode_blacklist_bed.exists():
            _progress("annotating ENCODE blacklist overlaps")
            blacklist_norm = tmp / "blacklist.norm.bed"
            _normalize_track_to_bed(encode_blacklist_bed, blacklist_norm)
            blacklist_hits = _get_overlapping_row_ids_with_fraction(
                candidate_bed, blacklist_norm, encode_blacklist_min_fraction
            )
            out["flag_encode_blacklist"] = out["row_id"].isin(blacklist_hits)
            _progress(f"blacklist hits={len(blacklist_hits)}")
            if not mate_with_id.empty:
                blacklist_mate_hits = _get_overlapping_row_ids_with_fraction(
                    mate_bed, blacklist_norm, encode_blacklist_min_fraction
                )
                out["flag_mate_in_blacklist"] = out["row_id"].isin(blacklist_mate_hits)
                _progress(f"blacklist mate hits={len(blacklist_mate_hits)}")

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
    _progress("finished junk-region annotation")
    return out


def build_candidate_loci(
    evidence_dir: Path,
    outdir: Path,
    window_size: int = 200,
    split_cluster_bp: int = 30,
    discordant_cluster_bp: int = 400,
    max_locus_span_bp: int = 2000,
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
    global _RUN_T0
    _RUN_T0 = time.monotonic()
    try:
        _progress(f"start build_candidate_loci evidence_dir={evidence_dir} outdir={outdir}")
        outdir.mkdir(parents=True, exist_ok=True)
        _progress("reading split evidence summary")
        passing_counts = _read_passing_counts(evidence_dir / "split_evidence.summary.tsv")

        _progress("loading split/discordant evidence tables")
        split_tumor_raw = _load_evidence_table(evidence_dir, "split_evidence", "tumor")
        split_normal_raw = _load_evidence_table(evidence_dir, "split_evidence", "normal")
        discordant_tumor_raw = _load_evidence_table(evidence_dir, "discordant_evidence", "tumor")
        discordant_normal_raw = _load_evidence_table(evidence_dir, "discordant_evidence", "normal")
        _progress(
            "loaded evidence rows "
            f"split_tumor={len(split_tumor_raw)}, split_normal={len(split_normal_raw)}, "
            f"discordant_tumor={len(discordant_tumor_raw)}, discordant_normal={len(discordant_normal_raw)}"
        )

        loci = _build_loci_from_evidence(
        split_tumor=split_tumor_raw,
        split_normal=split_normal_raw,
        discordant_tumor=discordant_tumor_raw,
        discordant_normal=discordant_normal_raw,
        split_cluster_bp=split_cluster_bp,
        discordant_cluster_bp=discordant_cluster_bp,
        max_locus_span_bp=max_locus_span_bp,
    )
        _progress(f"locus assignment target count={len(loci)}")

        _progress("assigning evidence rows to loci")
        split_tumor = _assign_rows_to_loci(split_tumor_raw, loci)
        split_normal = _assign_rows_to_loci(split_normal_raw, loci)
        discordant_tumor = _assign_rows_to_loci(discordant_tumor_raw, loci)
        discordant_normal = _assign_rows_to_loci(discordant_normal_raw, loci)
        _progress(
            "assigned rows "
            f"split_tumor={len(split_tumor)}, split_normal={len(split_normal)}, "
            f"discordant_tumor={len(discordant_tumor)}, discordant_normal={len(discordant_normal)}"
        )

        _progress("aggregating split/discordant metrics per locus")
        split_t = _aggregate_split_metrics(split_tumor, "split_tumor")
        split_n = _aggregate_split_metrics(split_normal, "split_normal")
        disc_t = _aggregate_discordant_metrics(discordant_tumor, "discordant_tumor")
        disc_n = _aggregate_discordant_metrics(discordant_normal, "discordant_normal")

        _progress("merging per-sample locus metrics")
        merged = split_t.merge(split_n, on=["chrom", "window_start", "window_end"], how="outer")
        merged = merged.merge(disc_t, on=["chrom", "window_start", "window_end"], how="outer")
        merged = merged.merge(disc_n, on=["chrom", "window_start", "window_end"], how="outer")
        merged = merged.fillna(0)

        int_cols = [
        c
        for c in merged.columns
        if c.endswith("_rows")
        or c.endswith("_unique_reads")
        or c.endswith("_mapq_min")
        or c.endswith("_clip_len_max")
        or c.endswith("_run_max")
        or c.endswith("_nm_max")
    ]
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

        _progress("computing junk-region flags")
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

        _progress("sorting final candidate loci")
        merged = merged.sort_values(["enrichment_ratio", "tumor_total_rows"], ascending=[False, False], kind="mergesort")

        out_tsv = outdir / "candidate_loci.tsv"
        out_parquet = outdir / "candidate_loci.parquet"
        _progress(f"writing outputs tsv={out_tsv} parquet={out_parquet}")
        merged.to_csv(out_tsv, sep="\t", index=False)
        merged.to_parquet(out_parquet, index=False)
        _progress(f"done build_candidate_loci rows={len(merged)}")
        return out_tsv
    finally:
        _RUN_T0 = None
