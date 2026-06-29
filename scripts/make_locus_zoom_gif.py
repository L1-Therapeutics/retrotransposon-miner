#!/usr/bin/env python3
"""Generate jump-cut GIFs across random clean 1 kb windows, ending at a MEI locus."""

from __future__ import annotations

import argparse
import os
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from intervaltree import IntervalTree
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from retro_miner.igv_plots import (  # noqa: E402
    _estimate_panel_height,
    _resolve_bam_index,
    resolve_igv_launcher,
    run_igv_batch,
)

CHR22_LENGTH = 50_818_468
DEFAULT_PUBLIC = Path.home() / "retrotransposon-workdir" / "data" / "public"
DEFAULT_RESULTS = Path.home() / "retrotransposon-workdir" / "results"

# All frames use fixed 1 kb windows (IGV shows reads/coverage well at this scale).
DEFAULT_WINDOW_BP = 1_000
FIXED_PANEL_HEIGHT = 380
MAX_CONTENT_HEIGHT = 460
CANVAS_DEFAULT = 1080
DEFAULT_TOTAL_MS = 24_000
FINAL_CAPTION_LINES = ("Genomic scar from", "retrotransposon insertion")


@dataclass(frozen=True)
class LocusSpec:
    chrom: str
    center: int
    final_start: int
    final_end: int
    caption_plain: str
    family: str
    subfamily: str
    tsd: str
    poly_a_run: int
    orientation: str
    control_support: str
    known_source: str

    @property
    def coord_label(self) -> str:
        return f"{self.chrom}:{self.final_start:,}-{self.final_end:,} (GRCh38)"


@dataclass(frozen=True)
class ViewWindow:
    start: int
    end: int
    phase: str  # scan | zoom | final

    @property
    def width(self) -> int:
        return self.end - self.start


def _load_locus_from_tsv(tsv_path: Path, chrom: str, pos: int) -> LocusSpec:
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    subset = df.loc[df["chrom"].astype(str) == chrom].copy()
    if subset.empty:
        raise ValueError(f"No rows for {chrom} in {tsv_path}")

    hit = subset.loc[(subset["window_start"] <= pos) & (subset["window_end"] >= pos)]
    if hit.empty:
        hit = subset.iloc[(subset["consensus_insertion_breakpoint_pos"] - pos).abs().argsort()[:1]]
    row = hit.iloc[0]

    final_start = int(row.get("window_start", pos) or pos)
    final_end = int(row.get("window_end", pos) or pos)
    center = int(row.get("consensus_insertion_breakpoint_pos", 0) or 0)
    if center <= 0:
        center = (final_start + final_end) // 2

    return LocusSpec(
        chrom=chrom,
        center=center,
        final_start=final_start,
        final_end=final_end,
        caption_plain="Genome scar from a viral-like retrotransposon (Alu)",
        family=str(row.get("consensus_mei_family", "") or "MEI"),
        subfamily=str(row.get("consensus_mei_subfamily", "") or ""),
        tsd=str(row.get("consensus_tsd_seq", "") or "").strip(),
        poly_a_run=int(row.get("consensus_poly_at_max_run", 0) or 0),
        orientation=str(row.get("consensus_insertion_orientation", "") or ""),
        control_support=str(row.get("control_supporting_reads", "") or ""),
        known_source=str(row.get("known_mei_polymorphism_source", "") or ""),
    )


@dataclass(frozen=True)
class AnnotationPaths:
    junk_bed: Path
    segdup_bed: Path
    low_mappability_bed: Path
    gap_bed: Path


def _default_annotation_paths(reference_build: str = "hg38") -> AnnotationPaths:
    base = DEFAULT_PUBLIC / "annotation" / reference_build
    return AnnotationPaths(
        junk_bed=base / "junk" / "junk_exclusion_merged.bed",
        segdup_bed=base / "segdup" / "genomicSuperDups.bed",
        low_mappability_bed=base
        / "mappability"
        / "k100.Umap.MultiTrackMappability.low_lt0.5.bed",
        gap_bed=base / "masks" / "gap.bed",
    )


def _load_chr_interval_tree(path: Path, chrom: str) -> IntervalTree:
    tree: IntervalTree = IntervalTree()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[0] != chrom:
                continue
            start, end = int(parts[1]), int(parts[2])
            if end > start:
                tree[start:end] = True
    return tree


def _window_is_clean(start: int, end: int, exclusion_trees: list[IntervalTree]) -> bool:
    return all(not tree.overlap(start, end) for tree in exclusion_trees)


def _find_clean_windows(
    *,
    chrom: str,
    chrom_length: int,
    window_bp: int,
    exclusion_trees: list[IntervalTree],
    cache_path: Path | None = None,
) -> list[tuple[int, int, int]]:
    if cache_path is not None and cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t")
        return [(int(r.start), int(r.end), int(r.center)) for r in df.itertuples(index=False)]

    step = max(250, window_bp // 4)
    clean: list[tuple[int, int, int]] = []
    for start in range(1, chrom_length - window_bp, step):
        end = start + window_bp
        if _window_is_clean(start, end, exclusion_trees):
            clean.append((start, end, (start + end) // 2))

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(clean, columns=["start", "end", "center"]).to_csv(
            cache_path, sep="\t", index=False
        )
    return clean


def _sample_random_clean_windows(
    clean_windows: list[tuple[int, int, int]],
    *,
    n: int,
    locus_center: int,
    locus_buffer: int,
    rng: random.Random,
) -> list[tuple[int, int, int]]:
    """Sample random clean 1 kb windows away from the target locus."""
    candidates = [w for w in clean_windows if abs(w[2] - locus_center) >= locus_buffer]
    if len(candidates) < n:
        raise RuntimeError(
            f"Only {len(candidates)} clean windows available (need {n}); "
            "relax filters or reduce --n-random."
        )
    return rng.sample(candidates, n)


def _locus_window_1kb(locus: LocusSpec, chrom_length: int) -> tuple[int, int, int]:
    half = DEFAULT_WINDOW_BP // 2
    start = max(1, locus.center - half)
    end = min(chrom_length, start + DEFAULT_WINDOW_BP)
    start = max(1, end - DEFAULT_WINDOW_BP)
    return start, end, locus.center


def _build_view_windows(
    locus: LocusSpec,
    *,
    annotation: AnnotationPaths,
    chrom_length: int,
    window_bp: int,
    n_random: int,
    random_seed: int,
    work_dir: Path,
) -> list[ViewWindow]:
    """Random clean 1 kb jump cuts, hard cut to 1 kb variant window (no zoom)."""
    for path in (
        annotation.junk_bed,
        annotation.segdup_bed,
        annotation.low_mappability_bed,
        annotation.gap_bed,
    ):
        if not path.exists():
            raise FileNotFoundError(f"Required annotation track not found: {path}")

    exclusion_trees = [
        _load_chr_interval_tree(annotation.junk_bed, locus.chrom),
        _load_chr_interval_tree(annotation.segdup_bed, locus.chrom),
        _load_chr_interval_tree(annotation.low_mappability_bed, locus.chrom),
        _load_chr_interval_tree(annotation.gap_bed, locus.chrom),
    ]
    cache_path = work_dir / f"clean_windows_{locus.chrom}_{window_bp}.tsv"
    clean_windows = _find_clean_windows(
        chrom=locus.chrom,
        chrom_length=chrom_length,
        window_bp=window_bp,
        exclusion_trees=exclusion_trees,
        cache_path=cache_path,
    )
    if not clean_windows:
        raise RuntimeError(f"No clean {window_bp} bp windows found on {locus.chrom}.")

    rng = random.Random(random_seed)
    random_regions = _sample_random_clean_windows(
        clean_windows,
        n=n_random,
        locus_center=locus.center,
        locus_buffer=max(window_bp * 2, 5_000),
        rng=rng,
    )
    variant_start, variant_end, _ = _locus_window_1kb(locus, chrom_length)

    path_tsv = work_dir / "random_clean_path.tsv"
    rows = [{"start": s, "end": e, "center": c, "phase": "random"} for s, e, c in random_regions]
    rows.append({"start": variant_start, "end": variant_end, "center": locus.center, "phase": "final"})
    pd.DataFrame(rows).to_csv(path_tsv, sep="\t", index=False)
    print(
        f"[zoom-gif] random clean path ({window_bp} bp, n={n_random}+1, seed={random_seed}): "
        f"{path_tsv}",
        flush=True,
    )
    for start, end, center in random_regions:
        print(f"  random {locus.chrom}:{start:,}-{end:,} (center {center:,})", flush=True)
    print(
        f"  final  {locus.chrom}:{variant_start:,}-{variant_end:,} (center {locus.center:,})",
        flush=True,
    )

    windows: list[ViewWindow] = [
        ViewWindow(start, end, "scan") for start, end, _ in random_regions
    ]
    windows.append(ViewWindow(variant_start, variant_end, "final"))
    return windows


def build_single_bam_igv_batch(
    *,
    reference_fasta: Path,
    bam_path: Path,
    snapshot_dir: Path,
    windows: list[tuple[str, int, int]],
    panel_heights: list[int],
    panel_height: int = FIXED_PANEL_HEIGHT,
) -> Path:
    bam_index = _resolve_bam_index(bam_path)
    if bam_index is None:
        raise FileNotFoundError(f"Missing BAM index for {bam_path}")

    lines: list[str] = [
        "new",
        f"genome {reference_fasta.resolve()}",
        f"snapshotDirectory {snapshot_dir.resolve()}",
        "preference SAM.SHOW_SOFT_CLIPPED true",
        f"load {bam_path.resolve()} index={bam_index.resolve()}",
    ]

    for idx, (chrom, start, end) in enumerate(windows, start=1):
        height = panel_heights[idx - 1] if idx - 1 < len(panel_heights) else panel_height
        stem = f"frame_{idx:02d}_{chrom}_{start}_{end}"
        lines.extend(
            [
                f"goto {chrom}:{start}-{end}",
                "expand",
                "sort position",
                f"maxPanelHeight {height}",
                f"snapshot {stem}.png",
            ]
        )

    lines.append("exit")
    batch_path = snapshot_dir / "igv_zoom_batch.txt"
    batch_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return batch_path


def _uniform_panel_height(
    bam_path: Path,
    chrom: str,
    view_windows: list[ViewWindow],
) -> int:
    """One shared maxPanelHeight for every frame (matches the deepest read pile)."""
    heights = [
        _estimate_panel_height(
            bam_path,
            bam_path,
            chrom,
            w.start,
            w.end,
            min_height=FIXED_PANEL_HEIGHT,
            max_height=8000,
        )
        for w in view_windows
    ]
    return max(heights)


def render_igv_frames(
    *,
    reference_fasta: Path,
    bam_path: Path,
    work_dir: Path,
    view_windows: list[ViewWindow],
    chrom: str,
    igv_launcher: Path | None,
    timeout_sec: int,
) -> list[Path]:
    igv_dir = work_dir / "igv_frames"
    igv_dir.mkdir(parents=True, exist_ok=True)
    igv_windows = [(chrom, w.start, w.end) for w in view_windows]
    panel_height = _uniform_panel_height(bam_path, chrom, view_windows)
    panel_heights = [panel_height] * len(view_windows)
    print(
        f"[zoom-gif] uniform maxPanelHeight={panel_height} for all {len(view_windows)} frames",
        flush=True,
    )
    batch_path = build_single_bam_igv_batch(
        reference_fasta=reference_fasta,
        bam_path=bam_path,
        snapshot_dir=igv_dir,
        windows=igv_windows,
        panel_heights=panel_heights,
    )
    print(f"[zoom-gif] running IGV batch ({len(view_windows)} snapshots)...", flush=True)
    run_igv_batch(batch_path, launcher=igv_launcher, timeout_sec=timeout_sec)

    paths: list[Path] = []
    for idx, window in enumerate(view_windows, start=1):
        stem = f"frame_{idx:02d}_{chrom}_{window.start}_{window.end}"
        png = igv_dir / f"{stem}.png"
        if not png.exists():
            raise RuntimeError(f"IGV snapshot missing: {png}")
        paths.append(png)
    return paths


def render_schematic_frame(
    *,
    gold_review_tsv: Path,
    locus: LocusSpec,
    out_path: Path,
    size: int,
) -> None:
    df = pd.read_csv(gold_review_tsv, sep="\t", low_memory=False)
    c22 = df.loc[df["chrom"].astype(str) == locus.chrom].copy()
    if "analysis_stage_tier" in c22.columns:
        c22 = c22.loc[c22["analysis_stage_tier"].astype(str).str.lower() == "gold"]

    pos = c22["consensus_insertion_breakpoint_pos"].fillna(-1).astype(int)
    xs = pos.where(pos > 0, (c22["window_start"] + c22["window_end"]) // 2).astype(float)

    dpi = 100
    fig_size = size / dpi
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=dpi)
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    ax.scatter(xs / 1e6, np.ones(len(xs)), s=12, c="#5b6b8a", alpha=0.55, linewidths=0, zorder=2)
    ax.scatter(
        [locus.center / 1e6],
        [1.0],
        s=160,
        c="#ffb020",
        edgecolors="#ffffff",
        linewidths=1.2,
        zorder=5,
    )

    ax.set_xlim(0, CHR22_LENGTH / 1e6)
    ax.set_ylim(0.84, 1.16)
    ax.set_yticks([])
    ax.set_xlabel("Chromosome 22 position (Mb)", color="#d7dce8", fontsize=12)
    ax.tick_params(colors="#aab3c5", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#3a4256")

    ax.set_title(
        locus.caption_plain,
        color="#f3f5fa",
        fontsize=16,
        pad=16,
        loc="center",
        fontweight="bold",
        wrap=True,
    )
    ax.text(
        locus.center / 1e6,
        0.875,
        f"target · {locus.center / 1e6:.2f} Mb",
        ha="center",
        va="top",
        color="#ffb020",
        fontsize=11,
    )
    ax.text(
        0.5,
        0.06,
        f"{len(c22):,} MEI candidate loci on chr22 (gold tier)",
        transform=ax.transAxes,
        ha="center",
        color="#9aa6bc",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(out_path, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def _detect_igv_content_band(arr: np.ndarray, *, tall_panel: bool = False) -> tuple[int, int, int, int]:
    """Return (left, top, right, bottom) for ruler+coverage+reads band."""
    h, w, _ = arr.shape
    gray = arr.mean(axis=2)

    col_active = (gray < 248).any(axis=0)
    active_cols = np.where(col_active)[0]
    left = int(active_cols[0]) if active_cols.size else 0

    align_left = left
    if tall_panel:
        for x in range(left, min(left + 220, w)):
            if (gray[:, x] > 248).mean() > 0.82:
                align_left = x + 1
        align_left = max(align_left, 100 if left == 0 else left)

    frac = (gray[:, align_left:] < 240).mean(axis=1)
    active_rows = np.where(frac > 0.02)[0]
    if active_rows.size == 0:
        cap = h if tall_panel else min(h, MAX_CONTENT_HEIGHT)
        return align_left, 0, w, cap

    top = max(0, int(active_rows[0]) - 2)
    if tall_panel:
        bottom = top
        for y in range(h - 1, top, -1):
            if frac[y] > 0.12:
                bottom = y
                break
        bottom = min(h, bottom + 20)
    else:
        bottom = min(h, top + MAX_CONTENT_HEIGHT)
    return align_left, top, w, bottom


def _load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates: list[Path] = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        fonts_dir = Path(conda_prefix) / "fonts"
        if bold:
            candidates.extend(
                [
                    fonts_dir / "Ubuntu-B.ttf",
                    fonts_dir / "DejaVuSans-Bold.ttf",
                ]
            )
        candidates.append(fonts_dir / "DejaVuSans.ttf")

    candidates.extend(
        [
            Path("/usr/share/fonts/google-noto-vf/NotoSans[wght].ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
            if bold
            else Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf")
            if bold
            else Path("/usr/share/fonts/dejavu/DejaVuSans.ttf"),
        ]
    )
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def _fit_caption_font(
    draw: ImageDraw.ImageDraw,
    lines: tuple[str, ...],
    canvas_size: int,
) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    max_w = canvas_size * 0.90
    for font_size in range(180, 36, -2):
        font = _load_font(font_size, bold=True)
        if max(draw.textlength(line, font=font) for line in lines) <= max_w:
            return font
    return _load_font(48, bold=True)


def _overlay_scar_caption(base: Image.Image, lines: tuple[str, ...] = FINAL_CAPTION_LINES) -> Image.Image:
    """Two-line caption anchored slightly left of centre."""
    img = base.copy()
    draw = ImageDraw.Draw(img)
    size = img.width
    font = _fit_caption_font(draw, lines, size)
    font_size = max(12, getattr(font, "size", 48) // 4)
    font = _load_font(font_size, bold=True)
    line_gap = max(4, int(font_size * 0.12))
    block_h = len(lines) * font_size + line_gap * (len(lines) - 1)

    anchor_x = int(size * 0.38)
    anchor_y = int(size * 0.45)
    y = anchor_y - block_h // 2

    for line in lines:
        tw = draw.textlength(line, font=font)
        x = max(8, int(anchor_x - tw // 2))
        draw.text(
            (x, y),
            line,
            font=font,
            fill=(0, 0, 0),
        )
        y += font_size + line_gap
    return img


def _fit_to_square(img: Image.Image, canvas_size: int) -> Image.Image:
    scale = canvas_size / max(img.width, img.height)
    new_w = max(1, int(img.width * scale))
    new_h = max(1, int(img.height * scale))
    resized = img.resize((new_w, new_h), Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (canvas_size, canvas_size), (15, 17, 23))
    offset = ((canvas_size - new_w) // 2, (canvas_size - new_h) // 2)
    canvas.paste(resized, offset)
    return canvas


def _format_igv_frame(
    path: Path,
    *,
    canvas_size: int,
) -> Image.Image:
    img = Image.open(path).convert("RGB")
    arr = np.array(img)
    left, top, right, bottom = _detect_igv_content_band(arr, tall_panel=True)
    cropped = img.crop((left, top, right, bottom))
    return _fit_to_square(cropped, canvas_size)


def _assemble_gif(
    frames: list[Image.Image],
    durations_ms: list[int],
    out_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=durations_ms,
        loop=0,
        optimize=True,
    )


def _durations_with_caption(
    frame_count: int,
    *,
    total_ms: int,
    final_hold_ms: int,
    caption_hold_ms: int,
) -> list[int]:
    """Last two frames: clean locus hold, then caption hold."""
    if frame_count <= 0:
        return []
    if frame_count == 1:
        return [total_ms]
    if frame_count == 2:
        clean = max(400, total_ms - caption_hold_ms)
        return [clean, caption_hold_ms]
    body_ms = max(0, total_ms - final_hold_ms - caption_hold_ms)
    per = max(400, body_ms // (frame_count - 2))
    durs = [per] * (frame_count - 2) + [final_hold_ms, caption_hold_ms]
    durs[-1] += total_ms - sum(durs)
    return durs


def build_gifs(
    *,
    locus: LocusSpec,
    reference_fasta: Path,
    bam_path: Path,
    gold_review_tsv: Path,
    annotation: AnnotationPaths,
    chrom_length: int,
    out_dir: Path,
    work_dir: Path,
    igv_launcher: Path | None,
    igv_timeout_sec: int,
    canvas_size: int,
    window_bp: int,
    n_random: int,
    random_seed: int,
    total_ms: int,
    final_hold_ms: int,
    caption_hold_ms: int,
) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    view_windows = _build_view_windows(
        locus,
        annotation=annotation,
        chrom_length=chrom_length,
        window_bp=window_bp,
        n_random=n_random,
        random_seed=random_seed,
        work_dir=work_dir,
    )
    igv_paths = render_igv_frames(
        reference_fasta=reference_fasta,
        bam_path=bam_path,
        work_dir=work_dir,
        view_windows=view_windows,
        chrom=locus.chrom,
        igv_launcher=igv_launcher,
        timeout_sec=igv_timeout_sec,
    )

    motion = [
        _format_igv_frame(path, canvas_size=canvas_size)
        for path in igv_paths
    ]

    final_clean = motion[-1]
    final_caption = _overlay_scar_caption(final_clean)

    outreach_frames = [*motion, final_caption]
    technical_frames = outreach_frames

    stem = (
        f"grch38_hg00100_known_longread_alu_{locus.chrom}_{locus.final_start}_{locus.final_end}"
    )
    outreach_path = out_dir / f"{stem}_outreach.gif"
    technical_path = out_dir / f"{stem}_technical.gif"

    _assemble_gif(
        outreach_frames,
        _durations_with_caption(
            len(outreach_frames),
            total_ms=total_ms,
            final_hold_ms=final_hold_ms,
            caption_hold_ms=caption_hold_ms,
        ),
        outreach_path,
    )
    _assemble_gif(
        technical_frames,
        _durations_with_caption(
            len(technical_frames),
            total_ms=total_ms,
            final_hold_ms=final_hold_ms,
            caption_hold_ms=caption_hold_ms,
        ),
        technical_path,
    )

    static_path = out_dir / f"{stem}.png"
    final_caption.save(static_path)

    print(f"[zoom-gif] wrote {outreach_path}", flush=True)
    print(f"[zoom-gif] wrote {technical_path}", flush=True)
    print(f"[zoom-gif] wrote {static_path}", flush=True)
    return outreach_path, technical_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", default="chr22")
    parser.add_argument("--center-pos", type=int, default=17_567_655)
    parser.add_argument("--final-start", type=int, default=17_567_258)
    parser.add_argument("--final-end", type=int, default=17_567_931)
    parser.add_argument(
        "--gold-review-tsv",
        type=Path,
        default=DEFAULT_RESULTS
        / "mei_step1_hg38_chr22_hg00100_control_only"
        / "candidate_loci.mei.gold_review.tsv",
    )
    parser.add_argument(
        "--bam",
        type=Path,
        default=DEFAULT_PUBLIC / "test_data/1kg_hg00100/chr22/hg00100.shortread.chr22.hg38.bam",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=DEFAULT_PUBLIC / "reference/hg38/Homo_sapiens_assembly38.fasta",
    )
    parser.add_argument("--out-dir", type=Path, default=ROOT / "docs" / "examples")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "workdir-results" / "locus_zoom_gif")
    parser.add_argument("--reference-build", default="hg38", choices=["hg38", "hg19", "hs1"])
    parser.add_argument("--annotation-dir", type=Path, default=None)
    parser.add_argument("--chrom-length", type=int, default=CHR22_LENGTH)
    parser.add_argument(
        "--window-bp",
        type=int,
        default=DEFAULT_WINDOW_BP,
        help="Fixed window size for every frame (default 1 kb)",
    )
    parser.add_argument(
        "--n-random",
        type=int,
        default=7,
        help="Number of random clean 1 kb windows before the variant frame",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        help="Seed for reproducible random window selection",
    )
    parser.add_argument("--igv-launcher", type=Path, default=None)
    parser.add_argument("--igv-timeout-sec", type=int, default=900)
    parser.add_argument(
        "--canvas-size",
        type=int,
        default=CANVAS_DEFAULT,
        help="Square output size (1080 works for LinkedIn/iPhone/README)",
    )
    parser.add_argument(
        "--total-ms",
        type=int,
        default=DEFAULT_TOTAL_MS,
        help="Target GIF duration in milliseconds (default 24 s)",
    )
    parser.add_argument(
        "--final-hold-ms",
        type=int,
        default=1800,
        help="Clean locus frame before caption (default 1.8 s)",
    )
    parser.add_argument(
        "--caption-hold-ms",
        type=int,
        default=7500,
        help="Pause on final frame with caption overlay",
    )
    args = parser.parse_args()

    if not args.gold_review_tsv.exists():
        raise SystemExit(f"gold review TSV not found: {args.gold_review_tsv}")
    if not args.bam.exists():
        raise SystemExit(f"BAM not found: {args.bam}")
    if not args.reference_fasta.exists():
        raise SystemExit(f"reference FASTA not found: {args.reference_fasta}")

    annotation = (
        _default_annotation_paths(args.reference_build)
        if args.annotation_dir is None
        else AnnotationPaths(
            junk_bed=args.annotation_dir / "junk" / "junk_exclusion_merged.bed",
            segdup_bed=args.annotation_dir / "segdup" / "genomicSuperDups.bed",
            low_mappability_bed=args.annotation_dir
            / "mappability"
            / "k100.Umap.MultiTrackMappability.low_lt0.5.bed",
            gap_bed=args.annotation_dir / "masks" / "gap.bed",
        )
    )

    locus = _load_locus_from_tsv(args.gold_review_tsv, args.chrom, args.center_pos)
    locus = LocusSpec(
        chrom=locus.chrom,
        center=locus.center,
        final_start=args.final_start,
        final_end=args.final_end,
        caption_plain=locus.caption_plain,
        family=locus.family,
        subfamily=locus.subfamily,
        tsd=locus.tsd,
        poly_a_run=locus.poly_a_run,
        orientation=locus.orientation,
        control_support=locus.control_support,
        known_source=locus.known_source,
    )
    launcher = args.igv_launcher or resolve_igv_launcher(None)

    build_gifs(
        locus=locus,
        reference_fasta=args.reference_fasta,
        bam_path=args.bam,
        gold_review_tsv=args.gold_review_tsv,
        annotation=annotation,
        chrom_length=args.chrom_length,
        out_dir=args.out_dir,
        work_dir=args.work_dir,
        igv_launcher=launcher,
        igv_timeout_sec=args.igv_timeout_sec,
        canvas_size=args.canvas_size,
        window_bp=args.window_bp,
        n_random=args.n_random,
        random_seed=args.random_seed,
        total_ms=args.total_ms,
        final_hold_ms=args.final_hold_ms,
        caption_hold_ms=args.caption_hold_ms,
    )


if __name__ == "__main__":
    main()
