from __future__ import annotations

import atexit
import json
import os
import platform
import re
import shutil
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pysam

_XVFB_PROC: subprocess.Popen[bytes] | None = None


def _headless_display_help() -> str:
    return (
        "IGV requires a graphical display. On headless Linux run: "
        "bash scripts/install_headless_igv_deps.sh "
        "(or: sudo dnf install -y xorg-x11-server-Xvfb xorg-x11-xauth)"
    )


def _find_xvfb_binary() -> str | None:
    for candidate in (
        shutil.which("Xvfb"),
        "/usr/bin/Xvfb",
        "/usr/local/bin/Xvfb",
    ):
        if candidate and Path(candidate).exists():
            return str(candidate)
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_xvfb = Path(conda_prefix) / "bin" / "Xvfb"
        if conda_xvfb.exists():
            return str(conda_xvfb)
    return None


def _needs_virtual_display() -> bool:
    if platform.system() != "Linux":
        return False
    display = os.environ.get("DISPLAY", "").strip()
    return not display


def _start_xvfb(display: str = ":99") -> subprocess.Popen[bytes]:
    global _XVFB_PROC
    xvfb = _find_xvfb_binary()
    if xvfb is None:
        raise RuntimeError(_headless_display_help())

    proc = subprocess.Popen(
        [xvfb, display, "-screen", "0", "1920x1080x24", "-nolisten", "tcp"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    time.sleep(0.5)
    if proc.poll() is not None:
        raise RuntimeError(f"Failed to start Xvfb for virtual display {display}")
    _XVFB_PROC = proc

    def _stop_xvfb() -> None:
        global _XVFB_PROC
        if _XVFB_PROC is not None and _XVFB_PROC.poll() is None:
            _XVFB_PROC.terminate()
            try:
                _XVFB_PROC.wait(timeout=5)
            except subprocess.TimeoutExpired:
                _XVFB_PROC.kill()
        _XVFB_PROC = None

    atexit.register(_stop_xvfb)
    return proc


@contextmanager
def _headless_display_env():
    env = os.environ.copy()
    if not _needs_virtual_display():
        yield env
        return

    if shutil.which("xvfb-run"):
        yield env
        return

    display = ":99"
    _start_xvfb(display)
    env["DISPLAY"] = display
    try:
        yield env
    finally:
        global _XVFB_PROC
        if _XVFB_PROC is not None and _XVFB_PROC.poll() is None:
            _XVFB_PROC.terminate()
            try:
                _XVFB_PROC.wait(timeout=5)
            except subprocess.TimeoutExpired:
                _XVFB_PROC.kill()
            _XVFB_PROC = None


def _resolve_bam_index(bam_path: Path) -> Path | None:
    for candidate in (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
        Path(f"{bam_path}.csi"),
        bam_path.with_suffix(".csi"),
    ):
        if candidate.exists():
            return candidate
    return None


def resolve_igv_launcher(launcher: Path | None = None) -> Path:
    if launcher is not None:
        if not launcher.exists():
            raise FileNotFoundError(f"IGV launcher not found: {launcher}")
        return launcher

    for name in ("igv", "igv.sh"):
        found = shutil.which(name)
        if found:
            return Path(found)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        for name in ("igv", "igv.sh"):
            candidate = Path(conda_prefix) / "bin" / name
            if candidate.exists():
                return candidate

    igv_home = os.environ.get("IGV_HOME")
    if igv_home:
        candidate = Path(igv_home) / "igv.sh"
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "IGV launcher not found. Install with: conda install -c bioconda igv "
        "(or set IGV_HOME / pass --igv-launcher)."
    )


def _count_reads_in_window(bam_path: Path, chrom: str, start: int, end: int) -> int:
    if end <= start:
        return 0
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            if chrom not in bam.references and not chrom.startswith("chr"):
                alt = f"chr{chrom}"
                if alt in bam.references:
                    chrom = alt
            return int(bam.count(chrom, max(0, start - 1), end))
    except (ValueError, OSError):
        return 0


def _estimate_panel_height(
    disease_bam: Path,
    control_bam: Path,
    chrom: str,
    start: int,
    end: int,
    *,
    pixels_per_read: int = 15,
    min_height: int = 250,
    max_height: int = 8000,
) -> int:
    disease_reads = _count_reads_in_window(disease_bam, chrom, start, end)
    control_reads = _count_reads_in_window(control_bam, chrom, start, end)
    stacked = disease_reads + control_reads + 2
    return max(min_height, min(max_height, stacked * pixels_per_read + 80))


def _safe_snapshot_stem(rank: int, chrom: str, start: int, end: int, *, contig_id: str = "") -> str:
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    contig_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(contig_id or "")).strip("_")
    if contig_safe:
        contig_safe = contig_safe[:32]
        return f"rank{rank:03d}_{chrom_safe}_{start}_{end}_{contig_safe}"
    return f"rank{rank:03d}_{chrom_safe}_{start}_{end}"


def _safe_locus_id(chrom: str, start: int, end: int) -> str:
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    return f"{chrom_safe}_{int(start)}_{int(end)}"


def _window_locus_id(chrom: str, window_start: int, window_end: int) -> str:
    s = int(window_start)
    e = int(window_end)
    if e < s:
        s, e = e, s
    return _safe_locus_id(chrom, max(1, s), max(1, e))


def _read_json_dict(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return data if isinstance(data, dict) else None


def _iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    out: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if name:
        out.append((name, "".join(seq_parts)))
    return out


def _build_assembly_contig_track(
    variants: pd.DataFrame,
    *,
    assembly_cache_dir: Path,
    reference_fasta: Path,
    snapshot_dir: Path,
) -> Path | None:
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        print("[igv-plots] minimap2/samtools unavailable; skipping contig alignment track", flush=True)
        return None

    contig_entries: list[tuple[str, str]] = []
    for rank, row in enumerate(variants.itertuples(index=False), start=1):
        chrom = str(getattr(row, "chrom", "") or "")
        window_start = int(getattr(row, "window_start", 0) or 0)
        window_end = int(getattr(row, "window_end", 0) or 0)
        contig_id = str(getattr(row, "assembly_best_contig_id", "") or "")
        if not chrom or window_start <= 0 or window_end <= 0 or not contig_id:
            continue

        stable_locus_dir = assembly_cache_dir / _window_locus_id(chrom, window_start, window_end)
        legacy_locus_dir = assembly_cache_dir / _safe_locus_id(chrom, window_start, window_end)
        locus_dir = stable_locus_dir if stable_locus_dir.exists() else legacy_locus_dir
        manifest = _read_json_dict(locus_dir / "assembly_manifest.json")
        if manifest is None:
            continue
        pad_bp = int(manifest.get("interval", {}).get("pad_bp", 250)) if isinstance(manifest.get("interval", {}), dict) else 250

        candidate_fastas = [
            locus_dir / f"disease.spades.pad{pad_bp}" / "contigs.fasta",
            locus_dir / f"control.spades.pad{pad_bp}" / "contigs.fasta",
        ]
        seq = ""
        source = ""
        for fa in candidate_fastas:
            recs = dict(_iter_fasta_records(fa))
            if contig_id in recs:
                seq = recs[contig_id]
                source = "disease" if "disease." in fa.name or "disease.spades" in str(fa.parent) else "control"
                break
        if not seq:
            continue
        header = f"rank{rank:03d}|{chrom}:{window_start}-{window_end}|{source}|{contig_id}"
        contig_entries.append((header, seq))

    if not contig_entries:
        print("[igv-plots] no assembly contigs resolved for IGV track", flush=True)
        return None

    query_fa = snapshot_dir / "assembly_selected_contigs.fa"
    with query_fa.open("w", encoding="utf-8") as oh:
        for name, seq in contig_entries:
            oh.write(f">{name}\n{seq}\n")

    bam_path = snapshot_dir / "assembly_selected_contigs.bam"
    cmd = f'minimap2 -a -x asm5 "{reference_fasta.resolve()}" "{query_fa.resolve()}" | samtools sort -o "{bam_path.resolve()}" -'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not bam_path.exists():
        detail = ((proc.stderr or "") + "\n" + (proc.stdout or "")).strip()[-2000:]
        print(f"[igv-plots] failed contig track alignment; skipping ({detail})", flush=True)
        return None
    idx_proc = subprocess.run(["samtools", "index", str(bam_path)], capture_output=True, text=True, check=False)
    if idx_proc.returncode != 0 or _resolve_bam_index(bam_path) is None:
        detail = ((idx_proc.stderr or "") + "\n" + (idx_proc.stdout or "")).strip()[-1000:]
        print(f"[igv-plots] failed indexing contig track; skipping ({detail})", flush=True)
        return None
    return bam_path


def _select_variants_for_plots(
    gold_review: pd.DataFrame,
    *,
    top_n: int,
    gold_only: bool,
) -> pd.DataFrame:
    if gold_review.empty:
        return gold_review.iloc[0:0].copy()
    subset = gold_review
    if gold_only and "analysis_stage_tier" in subset.columns:
        subset = subset.loc[subset["analysis_stage_tier"].fillna("").astype(str).str.lower() == "gold"]
    return subset.head(max(0, int(top_n))).copy()


def _write_contig_annotation_bed(variants: pd.DataFrame, snapshot_dir: Path) -> Path | None:
    rows: list[str] = []
    for row in variants.itertuples(index=False):
        chrom = str(getattr(row, "chrom", "") or "")
        start = int(getattr(row, "window_start", 0) or 0)
        end = int(getattr(row, "window_end", 0) or 0)
        contig = str(getattr(row, "assembly_best_contig_id", "") or "")
        complex_class = str(getattr(row, "asm_complex_class", "") or "")
        if not chrom or start <= 0 or end <= start or not contig:
            continue
        label = contig
        if complex_class:
            label = f"{contig}|{complex_class}"
        safe_label = re.sub(r"[\t\n\r]+", "_", label)
        rows.append(f"{chrom}\t{start - 1}\t{end}\t{safe_label}")
    if not rows:
        return None
    bed_path = snapshot_dir / "assembly_best_contigs.bed"
    bed_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return bed_path


def build_igv_batch_script(
    variants: pd.DataFrame,
    *,
    reference_fasta: Path,
    disease_bam: Path,
    control_bam: Path,
    snapshot_dir: Path,
    contig_annotation_bed: Path | None = None,
    contig_alignment_bam: Path | None = None,
    panel_height_min: int = 250,
    panel_height_max: int = 8000,
) -> str:
    disease_index = _resolve_bam_index(disease_bam)
    control_index = _resolve_bam_index(control_bam)
    if disease_index is None:
        raise FileNotFoundError(f"Missing BAM index for disease BAM: {disease_bam}")
    if control_index is None:
        raise FileNotFoundError(f"Missing BAM index for control BAM: {control_bam}")

    lines: list[str] = [
        "new",
        f"genome {reference_fasta.resolve()}",
        f"snapshotDirectory {snapshot_dir.resolve()}",
        "preference SAM.SHOW_SOFT_CLIPPED true",
        f"load {disease_bam.resolve()} index={disease_index.resolve()}",
        f"load {control_bam.resolve()} index={control_index.resolve()}",
    ]
    if contig_alignment_bam is not None:
        contig_idx = _resolve_bam_index(contig_alignment_bam)
        if contig_idx is not None:
            lines.append(f"load {contig_alignment_bam.resolve()} index={contig_idx.resolve()}")
    if contig_annotation_bed is not None and contig_annotation_bed.exists():
        lines.append(f"load {contig_annotation_bed.resolve()}")

    for rank, row in enumerate(variants.itertuples(index=False), start=1):
        chrom = str(getattr(row, "chrom", "") or "")
        start = int(getattr(row, "window_start", 0) or 0)
        end = int(getattr(row, "window_end", 0) or 0)
        if not chrom or start <= 0 or end <= start:
            continue
        panel_height = _estimate_panel_height(
            disease_bam,
            control_bam,
            chrom,
            start,
            end,
            min_height=panel_height_min,
            max_height=panel_height_max,
        )
        best_contig = str(getattr(row, "assembly_best_contig_id", "") or "")
        snapshot_name = _safe_snapshot_stem(rank, chrom, start, end, contig_id=best_contig)
        lines.extend(
            [
                f"goto {chrom}:{start}-{end}",
                "expand",
                "sort position",
                f"maxPanelHeight {panel_height}",
                f"snapshot {snapshot_name}.png",
            ]
        )

    lines.append("exit")
    return "\n".join(lines) + "\n"


def _wrap_headless_command(launcher: Path, batch_script: Path) -> list[str]:
    base = [str(launcher), "-b", str(batch_script)]
    if _needs_virtual_display() and shutil.which("xvfb-run"):
        return ["xvfb-run", "--auto-servernum", "--server-num=1", *base]
    return base


def _verify_snapshot_pngs(index_rows: list[dict[str, object]]) -> int:
    paths = [Path(str(row["snapshot_png"])) for row in index_rows if row.get("snapshot_png")]
    created = sum(1 for path in paths if path.exists() and path.stat().st_size > 0)
    if created == 0 and paths:
        raise RuntimeError(
            f"IGV produced 0/{len(paths)} snapshot PNGs. {_headless_display_help()}"
        )
    return created


def run_igv_batch(
    batch_script_path: Path,
    *,
    launcher: Path | None = None,
    timeout_sec: int | None = None,
) -> None:
    igv = resolve_igv_launcher(launcher)
    cmd = _wrap_headless_command(igv, batch_script_path)
    if _needs_virtual_display() and not shutil.which("xvfb-run"):
        if _find_xvfb_binary() is None:
            raise RuntimeError(_headless_display_help())

    with _headless_display_env() as env:
        result = subprocess.run(
            cmd,
            check=False,
            timeout=timeout_sec,
            env=env,
            capture_output=True,
            text=True,
        )
    combined = "\n".join(part for part in (result.stdout, result.stderr) if part)
    if result.returncode != 0 or "HeadlessException" in combined:
        detail = combined.strip()[-2000:] if combined.strip() else f"exit code {result.returncode}"
        raise RuntimeError(f"IGV batch run failed: {detail}")


def generate_gold_review_igv_plots(
    gold_review: pd.DataFrame,
    *,
    reference_fasta: Path,
    disease_bam: Path,
    control_bam: Path,
    snapshot_dir: Path,
    top_n: int = 100,
    gold_only: bool = True,
    launcher: Path | None = None,
    panel_height_min: int = 250,
    panel_height_max: int = 8000,
    timeout_sec: int | None = None,
    assembly_cache_dir: Path | None = None,
) -> Path | None:
    variants = _select_variants_for_plots(gold_review, top_n=top_n, gold_only=gold_only)
    if variants.empty:
        print("[igv-plots] no variants selected for snapshots; skipping", flush=True)
        return None

    snapshot_dir.mkdir(parents=True, exist_ok=True)
    contig_annotation_bed = _write_contig_annotation_bed(variants, snapshot_dir)
    contig_alignment_bam: Path | None = None
    if assembly_cache_dir is not None and assembly_cache_dir.exists():
        contig_alignment_bam = _build_assembly_contig_track(
            variants,
            assembly_cache_dir=assembly_cache_dir,
            reference_fasta=reference_fasta,
            snapshot_dir=snapshot_dir,
        )
    batch_script_path = snapshot_dir / "igv_batch.txt"
    batch_text = build_igv_batch_script(
        variants,
        reference_fasta=reference_fasta,
        disease_bam=disease_bam,
        control_bam=control_bam,
        snapshot_dir=snapshot_dir,
        contig_annotation_bed=contig_annotation_bed,
        contig_alignment_bam=contig_alignment_bam,
        panel_height_min=panel_height_min,
        panel_height_max=panel_height_max,
    )
    batch_script_path.write_text(batch_text, encoding="utf-8")

    index_rows: list[dict[str, object]] = []
    for rank, row in enumerate(variants.itertuples(index=False), start=1):
        chrom = str(getattr(row, "chrom", "") or "")
        start = int(getattr(row, "window_start", 0) or 0)
        end = int(getattr(row, "window_end", 0) or 0)
        if not chrom or start <= 0 or end <= start:
            continue
        best_contig = str(getattr(row, "assembly_best_contig_id", "") or "")
        snapshot_name = _safe_snapshot_stem(rank, chrom, start, end, contig_id=best_contig)
        index_rows.append(
            {
                "plot_rank": rank,
                "analysis_stage_tier": getattr(row, "analysis_stage_tier", ""),
                "chrom": chrom,
                "window_start": start,
                "window_end": end,
                "insertion_breakpoint_pos": getattr(row, "insertion_breakpoint_pos", -1),
                "mei_family": getattr(row, "mei_family", ""),
                "mei_subfamily": getattr(row, "mei_subfamily", ""),
                "insertion_model_score": getattr(row, "insertion_model_score", ""),
                "assembly_best_contig_id": best_contig,
                "asm_complex_class": getattr(row, "asm_complex_class", ""),
                "asm_breakpoint_side_status": getattr(row, "asm_breakpoint_side_status", ""),
                "snapshot_png": str(snapshot_dir / f"{snapshot_name}.png"),
            }
        )

    index_path = snapshot_dir / "igv_snapshot_index.tsv"
    pd.DataFrame(index_rows).to_csv(index_path, sep="\t", index=False)

    t0 = time.monotonic()
    print(
        f"[igv-plots] generating {len(index_rows)} snapshots in {snapshot_dir} "
        f"(top_n={top_n}, gold_only={gold_only})",
        flush=True,
    )
    run_igv_batch(batch_script_path, launcher=launcher, timeout_sec=timeout_sec)
    created = _verify_snapshot_pngs(index_rows)
    print(
        f"[igv-plots] wrote {created}/{len(index_rows)} snapshot PNGs; "
        f"index at {index_path} elapsed={time.monotonic() - t0:.1f}s",
        flush=True,
    )
    return index_path
