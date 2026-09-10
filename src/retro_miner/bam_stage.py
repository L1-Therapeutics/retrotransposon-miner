"""Stage remote WGS BAMs to local disk for multi-chrom / full-genome runs.

Single-chromosome runs may keep streaming ``s3://`` or ``http(s)://`` BAMs.
When the pipeline covers multiple chromosomes, ``--chr all``, or
``--chr_concurrency > 1``, remote alignments are copied next to their
BAI/CSI/CRAI so extract/annotate hit local files.

S3 copies use :mod:`retro_miner.s3_transfer` (64 concurrent 64 MiB parts).
No live object-store calls happen unless ``apply_bam_stage`` / ``--apply``.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

from retro_miner.s3_transfer import download_s3_uri

HEADROOM_RATIO = 0.10
MIN_HEADROOM_BYTES = 5 * 1024**3
REMOTE_SCHEMES = ("s3://", "http://", "https://")


def is_remote_alignment_uri(path: str | None) -> bool:
    text = (path or "").strip()
    return text.startswith(REMOTE_SCHEMES)


def default_bam_stage_dir(workdir: str | Path | None = None) -> Path:
    override = (os.environ.get("RTM_BAM_STAGE_DIR") or "").strip()
    if override:
        return Path(override).expanduser()
    root = Path(workdir or os.environ.get("RTM_WORKDIR") or Path.home() / "retrotransposon-workdir")
    return root.expanduser() / "data" / "bam_stage"


def bam_stage_enabled(raw: str | None = None) -> bool:
    text = (raw if raw is not None else os.environ.get("RTM_BAM_STAGE", "1")).strip().lower()
    return text not in {"0", "false", "no", "off"}


def should_stage_remote_bams(
    *,
    remote_bam_present: bool,
    chromosome_count: int,
    chr_all: bool = False,
    chr_concurrency: int = 1,
    enabled: bool = True,
) -> bool:
    """Return True when remote BAMs should be copied before extract/annotate."""
    if not enabled or not remote_bam_present:
        return False
    if chr_all or int(chromosome_count) > 1 or int(chr_concurrency) > 1:
        return True
    return False


def alignment_basename(uri: str) -> str:
    text = (uri or "").strip()
    if text.startswith("s3://"):
        path = text[len("s3://") :].split("/", 1)[-1] if "/" in text[len("s3://") :] else text
        name = Path(path).name
    else:
        name = Path(urllib.parse.urlparse(text).path).name
    if not name:
        raise ValueError(f"Could not derive alignment basename from {uri!r}")
    return name


def _split_alignment_name(name: str) -> tuple[str, str]:
    for suffix in (".bam", ".cram"):
        if name.endswith(suffix):
            return name[: -len(suffix)], suffix
    stem, ext = os.path.splitext(name)
    return stem, ext


def _uri_hash(uri: str) -> str:
    return hashlib.sha1(uri.encode("utf-8")).hexdigest()[:8]


def staged_alignment_dest(
    stage_dir: str | Path,
    remote_uri: str,
    *,
    colliding_basenames: set[str] | None = None,
) -> Path:
    """Local destination path for one remote BAM/CRAM."""
    name = alignment_basename(remote_uri)
    if colliding_basenames and name in colliding_basenames:
        stem, suffix = _split_alignment_name(name)
        name = f"{stem}.{_uri_hash(remote_uri)}{suffix}"
    return Path(stage_dir) / name


def plan_staged_dests(stage_dir: str | Path, uris: Sequence[str]) -> dict[str, Path]:
    """Map unique remote URIs to local dests; hash-suffix only on basename clashes."""
    unique = [u.strip() for u in uris if (u or "").strip()]
    unique = list(dict.fromkeys(unique))
    names = [alignment_basename(u) for u in unique]
    counts: dict[str, int] = {}
    for name in names:
        counts[name] = counts.get(name, 0) + 1
    colliding = {name for name, n in counts.items() if n > 1}
    return {u: staged_alignment_dest(stage_dir, u, colliding_basenames=colliding) for u in unique}


def index_sidecar_uris(alignment_uri: str) -> list[str]:
    """Likely remote index URLs/keys next to a BAM/CRAM."""
    uri = alignment_uri.rstrip("/")
    cands: list[str] = []
    if uri.endswith(".bam"):
        cands.extend([uri + ".bai", uri[:-4] + ".bai", uri + ".csi", uri[:-4] + ".csi"])
    elif uri.endswith(".cram"):
        cands.extend([uri + ".crai", uri[:-5] + ".crai"])
    else:
        cands.extend([uri + ".bai", uri + ".csi", uri + ".crai"])
    seen: set[str] = set()
    out: list[str] = []
    for cand in cands:
        if cand not in seen:
            seen.add(cand)
            out.append(cand)
    return out


def index_sidecar_dests(local_alignment: Path) -> list[Path]:
    name = local_alignment.name
    parent = local_alignment.parent
    if name.endswith(".bam"):
        return [
            Path(str(local_alignment) + ".bai"),
            parent / f"{name[:-4]}.bai",
            Path(str(local_alignment) + ".csi"),
            parent / f"{name[:-4]}.csi",
        ]
    if name.endswith(".cram"):
        return [
            Path(str(local_alignment) + ".crai"),
            parent / f"{name[:-5]}.crai",
        ]
    return [Path(str(local_alignment) + ".bai")]


def local_copy_is_complete(dest: str | Path, expected_size: int | None) -> bool:
    path = Path(dest)
    if expected_size is None or expected_size < 0:
        return False
    try:
        return path.is_file() and path.stat().st_size == int(expected_size)
    except OSError:
        return False


def required_free_bytes(
    pending_copy_sizes: Iterable[int],
    *,
    headroom_ratio: float = HEADROOM_RATIO,
    min_headroom_bytes: int = MIN_HEADROOM_BYTES,
) -> int:
    total = int(sum(int(x) for x in pending_copy_sizes if int(x) > 0))
    if total <= 0:
        return 0
    headroom = max(int(total * headroom_ratio), int(min_headroom_bytes))
    return total + headroom


def format_bytes(n: int) -> str:
    value = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{n} B"


def _local_index_present(local_alignment: Path) -> bool:
    return any(p.is_file() and p.stat().st_size > 0 for p in index_sidecar_dests(local_alignment))


@dataclass
class BamStagePlan:
    should_stage: bool
    stage_dir: Path
    dest_by_uri: dict[str, Path] = field(default_factory=dict)
    rewritten: dict[str, str] = field(default_factory=dict)


def plan_bam_stage(
    *,
    disease_bam: str,
    control_bam: str,
    disease_mate_bam: str = "",
    control_mate_bam: str = "",
    stage_dir: str | Path | None = None,
    chromosome_count: int,
    chr_all: bool = False,
    chr_concurrency: int = 1,
    enabled: bool = True,
) -> BamStagePlan:
    labels = {
        "DISEASE_BAM": disease_bam,
        "CONTROL_BAM": control_bam,
        "DISEASE_MATE_BAM": disease_mate_bam,
        "CONTROL_MATE_BAM": control_mate_bam,
    }
    remotes = [p for p in labels.values() if is_remote_alignment_uri(p)]
    dest_dir = Path(stage_dir) if stage_dir is not None else default_bam_stage_dir()
    do_stage = should_stage_remote_bams(
        remote_bam_present=bool(remotes),
        chromosome_count=chromosome_count,
        chr_all=chr_all,
        chr_concurrency=chr_concurrency,
        enabled=enabled,
    )
    dest_by_uri = plan_staged_dests(dest_dir, remotes) if do_stage else {}
    rewritten = {}
    for key, value in labels.items():
        if do_stage and is_remote_alignment_uri(value):
            rewritten[key] = str(dest_by_uri[value])
        else:
            rewritten[key] = value
    return BamStagePlan(
        should_stage=do_stage,
        stage_dir=dest_dir,
        dest_by_uri=dest_by_uri,
        rewritten=rewritten,
    )


def _run_cmd(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{err}")


def _default_head_size(uri: str) -> int | None:
    if uri.startswith("s3://"):
        if shutil.which("aws") is None:
            raise RuntimeError("aws CLI is required to stage s3:// BAMs")
        rest = uri[len("s3://") :]
        bucket, _, key = rest.partition("/")
        proc = subprocess.run(
            [
                "aws",
                "s3api",
                "head-object",
                "--bucket",
                bucket,
                "--key",
                key,
                "--query",
                "ContentLength",
                "--output",
                "text",
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if proc.returncode != 0:
            return None
        try:
            return int((proc.stdout or "").strip())
        except ValueError:
            return None
        request = urllib.request.Request(
        uri,
        method="HEAD",
        headers={"User-Agent": "retrotransposon-miner/bam-stage"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as resp:
            cl = resp.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:  # noqa: BLE001
        return None


def _default_copy(src: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.startswith("s3://"):
        download_s3_uri(src, dest)
        return
    if shutil.which("curl") is None:
        raise RuntimeError("curl is required to stage http(s):// BAMs")
    _run_cmd(["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-o", str(dest), src])


def _default_disk_free(path: Path) -> int:
    path.mkdir(parents=True, exist_ok=True)
    return int(shutil.disk_usage(path).free)


def apply_bam_stage(
    plan: BamStagePlan,
    *,
    head_size: Callable[[str], int | None] | None = None,
    copy_object: Callable[[str, Path], None] | None = None,
    disk_free: Callable[[Path], int] | None = None,
    log: Callable[[str], None] | None = None,
) -> dict[str, str]:
    """Copy planned remotes if dest size mismatches. Returns rewritten paths."""
    emit = log or (lambda msg: print(msg, file=sys.stderr))
    if not plan.should_stage:
        emit("[bam-stage] skip (single-chrom stream or staging disabled)")
        return dict(plan.rewritten)

    head = head_size or _default_head_size
    copy_fn = copy_object or _default_copy
    free_fn = disk_free or _default_disk_free
    plan.stage_dir.mkdir(parents=True, exist_ok=True)

    pending: list[tuple[str, Path, int]] = []
    for uri, dest in plan.dest_by_uri.items():
        size = head(uri)
        if size is None:
            raise RuntimeError(f"Could not determine size of remote BAM: {uri}")
        if local_copy_is_complete(dest, size) and _local_index_present(dest):
            emit(f"[bam-stage] reuse {dest} ({format_bytes(size)})")
            continue
        pending.append((uri, dest, size))

    need = required_free_bytes(size for _uri, _dest, size in pending)
    free = free_fn(plan.stage_dir)
    if need > free:
        raise RuntimeError(
            "not enough free space to stage remote BAMs.\n"
            f"  dest: {plan.stage_dir}\n"
            f"  need: {format_bytes(need)} "
            f"({format_bytes(sum(s for _u, _d, s in pending))} BAMs + headroom)\n"
            f"  free: {format_bytes(free)}\n"
            "Use --bam-stage-dir / RTM_BAM_STAGE_DIR on a larger volume, "
            "or disable with --no-bam-stage / RTM_BAM_STAGE=0."
        )

    def _copy_with_index(uri: str, dest: Path, size: int) -> None:
        if not local_copy_is_complete(dest, size):
            emit(f"[bam-stage] copy {uri} -> {dest} ({format_bytes(size)})")
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                dest.unlink()
            copy_fn(uri, dest)
            if not local_copy_is_complete(dest, size):
                got = dest.stat().st_size if dest.exists() else 0
                raise RuntimeError(
                    f"staged BAM size mismatch for {dest}: expected {size}, got {got}"
                )
        if _local_index_present(dest):
            return
        last_err: str | None = None
        for idx_uri in index_sidecar_uris(uri):
            idx_name = alignment_basename(idx_uri)
            idx_dest = dest.parent / idx_name
            try:
                emit(f"[bam-stage] copy index {idx_uri} -> {idx_dest}")
                copy_fn(idx_uri, idx_dest)
                if idx_dest.is_file() and idx_dest.stat().st_size > 0:
                    return
            except Exception as exc:  # noqa: BLE001
                last_err = str(exc)
                continue
        raise RuntimeError(
            f"Could not stage BAI/CSI/CRAI next to {dest}"
            + (f" ({last_err})" if last_err else "")
        )

    if pending:
        workers = min(4, len(pending))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(_copy_with_index, uri, dest, size) for uri, dest, size in pending]
            for fut in as_completed(futs):
                fut.result()
    emit(f"[bam-stage] ready under {plan.stage_dir}")
    return dict(plan.rewritten)


def write_env_file(path: str | Path, rewritten: Mapping[str, str]) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for key in ("DISEASE_BAM", "CONTROL_BAM", "DISEASE_MATE_BAM", "CONTROL_MATE_BAM"):
        val = rewritten.get(key, "")
        lines.append(f"{key}={_shell_single_quote(val)}")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _shell_single_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Stage remote BAMs for multi-chrom RTM runs.")
    p.add_argument("--disease-bam", required=True)
    p.add_argument("--control-bam", required=True)
    p.add_argument("--disease-mate-bam", default="")
    p.add_argument("--control-mate-bam", default="")
    p.add_argument("--stage-dir", default="")
    p.add_argument("--chr-count", type=int, required=True)
    p.add_argument("--chr-concurrency", type=int, default=1)
    p.add_argument("--chr-all", action="store_true")
    p.add_argument("--disabled", action="store_true", help="Do not stage (RTM_BAM_STAGE=0).")
    p.add_argument("--out-env", required=True, help="Write DISEASE_BAM=... shell assignments.")
    p.add_argument("--apply", action="store_true", help="Copy remotes when the plan says to stage.")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    # Wrapper already merged RTM_BAM_STAGE / --no-bam-stage into --disabled.
    enabled = not args.disabled
    stage_dir = Path(args.stage_dir) if args.stage_dir else default_bam_stage_dir()
    plan = plan_bam_stage(
        disease_bam=args.disease_bam,
        control_bam=args.control_bam,
        disease_mate_bam=args.disease_mate_bam,
        control_mate_bam=args.control_mate_bam,
        stage_dir=stage_dir,
        chromosome_count=args.chr_count,
        chr_all=args.chr_all,
        chr_concurrency=args.chr_concurrency,
        enabled=enabled,
    )
    rewritten = apply_bam_stage(plan) if args.apply else dict(plan.rewritten)
    write_env_file(args.out_env, rewritten)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
