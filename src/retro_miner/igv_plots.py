from __future__ import annotations

import atexit
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

import click
import pandas as pd
import pysam

from ._utils import _iter_fasta_records, safe_locus_id as _safe_locus_id

_XVFB_PROC: subprocess.Popen[bytes] | None = None


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


@contextmanager
def _igv_singleton_lock(
    *,
    lock_path: Path | None = None,
    timeout_sec: float = 600.0,
    poll_sec: float = 1.0,
):
    """Serialize IGV batch runs to avoid java.net.BindException port collisions.

    The default lock path includes the current user's UID so distinct OS users
    cannot interfere with each other's IGV singleton coordination.
    """
    lock_file = lock_path or Path(tempfile.gettempdir()) / f"retro_miner_igv_{os.getuid()}.lock"
    start = time.monotonic()
    owner_fd: int | None = None
    while owner_fd is None:
        try:
            owner_fd = os.open(str(lock_file), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(owner_fd, f"{os.getpid()}\t{int(time.time())}\n".encode("utf-8"))
            break
        except FileExistsError:
            stale = False
            try:
                raw = lock_file.read_text(encoding="utf-8").strip()
                owner_pid = int(raw.split("\t")[0]) if raw else 0
                stale = not _pid_is_alive(owner_pid)
            except (OSError, ValueError, IndexError):
                stale = True
            if stale:
                try:
                    lock_file.unlink()
                except FileNotFoundError:
                    pass
                continue
            if time.monotonic() - start >= float(timeout_sec):
                raise RuntimeError(
                    f"Timed out waiting for IGV singleton lock: {lock_file}. "
                    "Another IGV batch run may still be active."
                )
            time.sleep(max(0.1, float(poll_sec)))
    try:
        yield
    finally:
        if owner_fd is not None:
            try:
                os.close(owner_fd)
            except OSError:
                pass
        try:
            if lock_file.exists():
                raw = lock_file.read_text(encoding="utf-8").strip()
                owner_pid = int(raw.split("\t")[0]) if raw else 0
                if owner_pid == os.getpid():
                    lock_file.unlink()
        except (OSError, ValueError, IndexError):
            pass


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
    """Return a filesystem-safe IGV snapshot stem.

    Special characters in *chrom* are replaced with underscores. When
    sanitisation changes the value, a 6-hex-character SHA-256 digest of the
    original string is embedded to prevent distinct chromosome names that
    sanitise identically from producing the same snapshot filename.
    """
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    if chrom_safe != str(chrom):
        digest = hashlib.sha256(str(chrom).encode()).hexdigest()[:6]
        chrom_safe = f"{chrom_safe}_{digest}"
    contig_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(contig_id or "")).strip("_")
    if contig_safe:
        contig_safe = contig_safe[:32]
        return f"rank{rank:03d}_{chrom_safe}_{start}_{end}_{contig_safe}"
    return f"rank{rank:03d}_{chrom_safe}_{start}_{end}"


_IGV_BATCH_FORBIDDEN: frozenset[str] = frozenset("\n\r\t")


def _validate_igv_chrom(chrom: str) -> None:
    """Raise ValueError if *chrom* contains characters that would corrupt an IGV batch command.

    Newlines, carriage returns, and tabs are forbidden because a ``goto`` line such as
    ``goto chr22\\nBAD:100-200`` would be split into two separate batch commands, causing
    IGV to receive an unintended command and to navigate to the wrong locus.
    """
    for ch in chrom:
        if ch in _IGV_BATCH_FORBIDDEN:
            raise ValueError(
                f"IGV batch chromosome contains a forbidden character (0x{ord(ch):02x}): {chrom!r}"
            )


def _quote_igv_path(p: Path) -> str:
    """Return *p* wrapped in double quotes for use in an IGV batch command."""
    s = str(p)
    if '"' in s:
        raise ValueError(
            f"IGV batch path contains a forbidden double-quote character: {s!r}"
        )
    return f'"{s}"'


def _window_locus_id(chrom: str, window_start: int, window_end: int) -> str:
    s = int(window_start)
    e = int(window_end)
    if e < s:
        s, e = e, s
    return _safe_locus_id(chrom, max(1, s), max(1, e))


def _row_discovery_window(row: object) -> tuple[str, int, int]:
    """Prefer discovery span for IGV/assembly-cache keys; fall back to reported window."""
    chrom = str(getattr(row, "chrom", "") or "")
    disc_start = int(getattr(row, "discovery_window_start", 0) or 0)
    disc_end = int(getattr(row, "discovery_window_end", 0) or 0)
    if disc_start > 0 and disc_end >= disc_start:
        return chrom, disc_start, disc_end
    start = int(getattr(row, "window_start", 0) or 0)
    end = int(getattr(row, "window_end", 0) or 0)
    return chrom, start, end


def _read_json_dict(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def _build_assembly_contig_track(
    variants: pd.DataFrame,
    *,
    assembly_cache_dir: Path,
    reference_fasta: Path,
    snapshot_dir: Path,
) -> Path | None:
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        click.echo("[igv-plots] minimap2/samtools unavailable; skipping contig alignment track")
        return None

    contig_entries: list[tuple[str, str]] = []
    for rank, row in enumerate(variants.itertuples(index=False), start=1):
        chrom, window_start, window_end = _row_discovery_window(row)
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
        click.echo("[igv-plots] no assembly contigs resolved for IGV track")
        return None

    query_fa = snapshot_dir / "assembly_selected_contigs.fa"
    with query_fa.open("w", encoding="utf-8") as oh:
        for name, seq in contig_entries:
            oh.write(f">{name}\n{seq}\n")

    bam_path = snapshot_dir / "assembly_selected_contigs.bam"
    try:
        minimap2_proc = subprocess.Popen(
            [
                "minimap2",
                "-a",
                "-x",
                "asm5",
                str(reference_fasta.resolve()),
                str(query_fa.resolve()),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        samtools_proc = subprocess.Popen(
            [
                "samtools",
                "sort",
                "-o",
                str(bam_path.resolve()),
                "-",
            ],
            stdin=minimap2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as exc:
        click.echo(f"[igv-plots] failed contig track alignment; skipping ({exc})")
        return None
    assert minimap2_proc.stdout is not None  # always set when stdout=PIPE; satisfies type checkers
    minimap2_proc.stdout.close()  # release our reference so minimap2 receives SIGPIPE if samtools exits early
    _, samtools_stderr = samtools_proc.communicate()
    minimap2_stderr = minimap2_proc.stderr.read() if minimap2_proc.stderr is not None else b""
    minimap2_rc = minimap2_proc.wait()
    samtools_rc = samtools_proc.returncode
    if minimap2_rc != 0 or samtools_rc != 0 or not bam_path.exists():
        detail = (
            samtools_stderr.decode("utf-8", errors="replace")
            + "\n"
            + minimap2_stderr.decode("utf-8", errors="replace")
        ).strip()[-2000:]
        click.echo(f"[igv-plots] failed contig track alignment; skipping ({detail})")
        if bam_path.exists():
            try:
                bam_path.unlink()
            except OSError:
                pass
        return None
    idx_proc = subprocess.run(["samtools", "index", str(bam_path)], capture_output=True, text=True, check=False)
    if idx_proc.returncode != 0 or _resolve_bam_index(bam_path) is None:
        detail = ((idx_proc.stderr or "") + "\n" + (idx_proc.stdout or "")).strip()[-1000:]
        click.echo(f"[igv-plots] failed indexing contig track; skipping ({detail})")
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
    limit = int(top_n)
    if limit > 0:
        return subset.head(limit).copy()
    return subset.copy()


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
        f"genome {_quote_igv_path(reference_fasta.resolve())}",
        f"snapshotDirectory {_quote_igv_path(snapshot_dir.resolve())}",
        "preference SAM.SHOW_SOFT_CLIPPED true",
        f"load {_quote_igv_path(disease_bam.resolve())} index={_quote_igv_path(disease_index.resolve())}",
        f"load {_quote_igv_path(control_bam.resolve())} index={_quote_igv_path(control_index.resolve())}",
    ]
    if contig_alignment_bam is not None:
        contig_idx = _resolve_bam_index(contig_alignment_bam)
        if contig_idx is not None:
            lines.append(
                f"load {_quote_igv_path(contig_alignment_bam.resolve())} index={_quote_igv_path(contig_idx.resolve())}"
            )
    if contig_annotation_bed is not None and contig_annotation_bed.exists():
        lines.append(f"load {_quote_igv_path(contig_annotation_bed.resolve())}")

    for rank, row in enumerate(variants.itertuples(index=False), start=1):
        chrom, start, end = _row_discovery_window(row)
        # Tight 1bp reported windows are hard to review; pad tiny spans for IGV view.
        view_start, view_end = start, end
        if view_end - view_start < 50:
            mid = (view_start + view_end) // 2
            view_start = max(1, mid - 100)
            view_end = mid + 100
        if not chrom or start <= 0 or end < start:
            continue
        try:
            _validate_igv_chrom(chrom)
        except ValueError as exc:
            click.echo(f"[igv-plots] skipping row {rank}: invalid chromosome \u2014 {exc}")
            continue
        panel_height = _estimate_panel_height(
            disease_bam,
            control_bam,
            chrom,
            view_start,
            view_end,
            min_height=panel_height_min,
            max_height=panel_height_max,
        )
        best_contig = str(getattr(row, "assembly_best_contig_id", "") or "")
        snapshot_name = _safe_snapshot_stem(rank, chrom, start, end, contig_id=best_contig)
        lines.extend(
            [
                f"goto {chrom}:{view_start}-{view_end}",
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
    created = 0
    for path in paths:
        try:
            if path.stat().st_size > 0:
                created += 1
        except OSError:
            continue
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
    bind_retry_attempts: int = 3,
    bind_retry_sleep_sec: float = 2.0,
) -> None:
    igv = resolve_igv_launcher(launcher)
    cmd = _wrap_headless_command(igv, batch_script_path)
    if _needs_virtual_display() and not shutil.which("xvfb-run"):
        if _find_xvfb_binary() is None:
            raise RuntimeError(_headless_display_help())

    max_attempts = max(1, int(bind_retry_attempts))
    bind_sleep = max(0.1, float(bind_retry_sleep_sec))
    last_detail = ""
    with _igv_singleton_lock():
        for attempt in range(1, max_attempts + 1):
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
            bind_error = ("java.net.BindException" in combined) or ("Address already in use" in combined)
            headless_error = "HeadlessException" in combined
            failed = (result.returncode != 0) or headless_error or bind_error
            if not failed:
                return
            detail = combined.strip()[-2000:] if combined.strip() else f"exit code {result.returncode}"
            last_detail = detail
            if bind_error and attempt < max_attempts:
                wait_sec = bind_sleep * float(attempt)
                click.echo(
                    f"[igv-plots] transient BindException on attempt {attempt}/{max_attempts}; "
                    f"retrying in {wait_sec:.1f}s"
                )
                time.sleep(wait_sec)
                continue
            break
    raise RuntimeError(f"IGV batch run failed: {last_detail}")


def generate_gold_review_igv_plots(
    gold_review: pd.DataFrame,
    *,
    reference_fasta: Path,
    disease_bam: Path,
    control_bam: Path,
    snapshot_dir: Path,
    top_n: int = 0,
    gold_only: bool = True,
    launcher: Path | None = None,
    panel_height_min: int = 250,
    panel_height_max: int = 8000,
    timeout_sec: int | None = None,
    assembly_cache_dir: Path | None = None,
) -> Path | None:
    variants = _select_variants_for_plots(gold_review, top_n=top_n, gold_only=gold_only)
    if variants.empty:
        click.echo("[igv-plots] no variants selected for snapshots; skipping")
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
        chrom, start, end = _row_discovery_window(row)
        if not chrom or start <= 0 or end < start:
            continue
        best_contig = str(getattr(row, "assembly_best_contig_id", "") or "")
        snapshot_name = _safe_snapshot_stem(rank, chrom, start, end, contig_id=best_contig)
        index_rows.append(
            {
                "plot_rank": rank,
                "analysis_stage_tier": getattr(row, "analysis_stage_tier", ""),
                "chrom": chrom,
                "window_start": int(getattr(row, "window_start", start) or start),
                "window_end": int(getattr(row, "window_end", end) or end),
                "discovery_window_start": start,
                "discovery_window_end": end,
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
    click.echo(
        f"[igv-plots] generating {len(index_rows)} snapshots in {snapshot_dir} "
        f"(top_n={'all' if top_n <= 0 else top_n}, gold_only={gold_only})"
    )
    run_igv_batch(batch_script_path, launcher=launcher, timeout_sec=timeout_sec)
    created = _verify_snapshot_pngs(index_rows)
    click.echo(
        f"[igv-plots] wrote {created}/{len(index_rows)} snapshot PNGs; "
        f"index at {index_path} elapsed={time.monotonic() - t0:.1f}s"
    )
    return index_path
