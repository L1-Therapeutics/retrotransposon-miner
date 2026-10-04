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
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence

from retro_miner.s3_transfer import download_s3_uri, split_s3_uri

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


def _uri_hash(uri: str) -> str:
    # Never persist signed URLs or access tokens in cache metadata.
    return hashlib.sha256(uri.strip().encode("utf-8")).hexdigest()


def staged_alignment_dest(
    stage_dir: str | Path,
    remote_uri: str,
    *,
    colliding_basenames: set[str] | None = None,
) -> Path:
    """Source-isolated destination; basename collisions cannot cross runs.

    ``colliding_basenames`` is retained for API compatibility. Every source is
    namespaced, regardless of the other inputs in the current plan.
    """
    name = alignment_basename(remote_uri)
    if name in {".", "..", "manifest.json"}:
        raise ValueError("Alignment basename must be a file name distinct from manifest.json")
    return Path(stage_dir) / "objects" / _uri_hash(remote_uri) / name


def plan_staged_dests(stage_dir: str | Path, uris: Sequence[str]) -> dict[str, Path]:
    """Map unique remote URIs to stable source-isolated destinations.

    Apply resolves these provisional paths to immutable remote generations.
    Planning does not access remote storage or create local files.
    """
    unique = list(dict.fromkeys(u.strip() for u in uris if (u or "").strip()))
    return {u: staged_alignment_dest(stage_dir, u) for u in unique}


def index_sidecar_uris(alignment_uri: str) -> list[str]:
    """Likely remote index URLs/keys next to a BAM/CRAM."""
    uri = alignment_uri.strip()
    # HTTP query strings (including signed access tokens) belong after the
    # modified path. S3 keys, in contrast, can contain literal '?' and '#'.
    parsed = urllib.parse.urlsplit(uri) if uri.startswith(("http://", "https://")) else None
    path = parsed.path if parsed is not None else uri
    cands: list[str] = []
    if path.endswith(".bam"):
        cands.extend([path + ".bai", path[:-4] + ".bai", path + ".csi", path[:-4] + ".csi"])
    elif path.endswith(".cram"):
        cands.extend([path + ".crai", path[:-5] + ".crai"])
    else:
        cands.extend([path + ".bai", path + ".csi", path + ".crai"])
    if parsed is not None:
        cands = [urllib.parse.urlunsplit(parsed._replace(path=cand)) for cand in cands]
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
    return [Path(str(local_alignment) + suffix) for suffix in (".bai", ".csi", ".crai")]


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
            rewritten[key] = str(dest_by_uri[value.strip()])
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


@dataclass(frozen=True)
class RemoteMetadata:
    size: int
    etag: str = ""
    last_modified: str = ""
    version_id: str = ""

    @property
    def reusable(self) -> bool:
        # Weak HTTP ETags are not byte-identity validators. Last-Modified is
        # useful as an additional change signal, but not sufficient alone.
        return bool(self.version_id or (self.etag and not self.etag.startswith("W/")))


def _default_head_metadata(uri: str) -> RemoteMetadata | None:
    if uri.startswith("s3://"):
        if shutil.which("aws") is None:
            raise RuntimeError("aws CLI is required to stage s3:// BAMs")
        bucket, key = split_s3_uri(uri)
        proc = subprocess.run(
            ["aws", "s3api", "head-object", "--bucket", bucket, "--key", key, "--output", "json"],
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60,
        )
        if proc.returncode != 0:
            return None
        try:
            obj = json.loads(proc.stdout)
            version = obj.get("VersionId") or ""
            return RemoteMetadata(
                size=int(obj["ContentLength"]), etag=str(obj.get("ETag") or ""),
                last_modified=str(obj.get("LastModified") or ""),
                version_id="" if version == "null" else str(version),
            )
        except (ValueError, KeyError, TypeError, AttributeError):
            return None
    if not uri.startswith(("http://", "https://")):
        return None
    request = urllib.request.Request(
        uri, method="HEAD", headers={"User-Agent": "retrotransposon-miner/bam-stage", "Accept-Encoding": "identity"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as resp:
            size = resp.headers.get("Content-Length")
            return RemoteMetadata(
                size=int(size), etag=resp.headers.get("ETag", ""),
                last_modified=resp.headers.get("Last-Modified", ""),
            ) if size is not None else None
    except (OSError, ValueError):
        return None


def _default_head_size(uri: str) -> int | None:
    metadata = _default_head_metadata(uri)
    return metadata.size if metadata is not None else None


@contextmanager
def _stage_lock(stage_dir: Path, timeout: float) -> Iterator[None]:
    """Serialize staging transactions; the kernel releases locks on exit/crash."""
    if not math.isfinite(timeout) or timeout < 0:
        raise ValueError("lock_timeout must be finite and >= 0")
    try:
        import fcntl
    except ImportError as exc:
        raise RuntimeError("BAM staging locks require a POSIX filesystem (Linux/macOS)") from exc
    stage_dir.mkdir(parents=True, exist_ok=True)
    # Keep this inode: unlinking a lock file lets different processes lock
    # different inodes and defeats mutual exclusion.
    with (stage_dir / ".rtm-stage.lock").open("a+b") as handle:
        deadline = time.monotonic() + timeout
        while True:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise RuntimeError(f"Timed out waiting for BAM staging lock under {stage_dir}") from None
                time.sleep(min(0.05, max(0.0, deadline - time.monotonic())))
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _file_fingerprint(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "inode": stat.st_ino}


def _cache_valid(
    dest: Path, metadata: RemoteMetadata, manifest: dict[str, object],
    head: Callable[[str], RemoteMetadata | None], uri: str,
) -> bool:
    if not metadata.reusable:
        return False
    try:
        saved = json.loads((dest.parent / "manifest.json").read_text(encoding="utf-8"))
        if not isinstance(saved, dict) or any(saved.get(k) != v for k, v in manifest.items()):
            return False
        index_name = saved.get("index_name")
        if not isinstance(index_name, str) or index_name not in {p.name for p in index_sidecar_dests(dest)}:
            return False
        index = dest.parent / index_name
        if not (
            local_copy_is_complete(dest, metadata.size) and index.is_file() and index.stat().st_size > 0
            and saved.get("alignment_file") == _file_fingerprint(dest)
            and saved.get("index_file") == _file_fingerprint(index)
        ):
            return False
        candidates = dict(zip((p.name for p in index_sidecar_dests(dest)), index_sidecar_uris(uri)))
        index_metadata = head(candidates[index_name])
        return bool(index_metadata is not None and index_metadata.reusable
                    and index_metadata.size == index.stat().st_size
                    and saved.get("index_remote") == asdict(index_metadata))
    except (OSError, ValueError, TypeError):
        return False


def _default_copy(src: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.startswith("s3://"):
        download_s3_uri(src, dest)
        return
    if shutil.which("curl") is None:
        raise RuntimeError("curl is required to stage http(s):// BAMs")
    try:
        _run_cmd(["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-o", str(dest), src])
    except RuntimeError:
        # curl stderr and command lines may contain signed URL credentials.
        raise RuntimeError(f"HTTP alignment/index transfer failed for {dest}") from None


def _default_disk_free(path: Path) -> int:
    path.mkdir(parents=True, exist_ok=True)
    return int(shutil.disk_usage(path).free)


def apply_bam_stage(
    plan: BamStagePlan,
    *,
    head_size: Callable[[str], int | None] | None = None,
    head_metadata: Callable[[str], RemoteMetadata | None] | None = None,
    copy_object: Callable[[str, Path], None] | None = None,
    disk_free: Callable[[Path], int] | None = None,
    log: Callable[[str], None] | None = None,
    lock_timeout: float = 3600,
) -> dict[str, str]:
    """Resolve immutable remote generations and stage under a bounded POSIX lock.

    A cache hit requires remote validators plus a matching completion manifest
    and local file fingerprints. Size-only providers are supported but cannot
    reuse cached data. The lock covers metadata, disk budgeting and publication.
    """
    emit = log or (lambda msg: print(msg, file=sys.stderr))
    if not plan.should_stage:
        emit("[bam-stage] skip (single-chrom stream or staging disabled)")
        return dict(plan.rewritten)
    if head_metadata is not None and head_size is not None:
        raise ValueError("Provide head_metadata or head_size, not both")

    def head(uri: str) -> RemoteMetadata | None:
        if head_metadata is not None:
            return head_metadata(uri)
        if head_size is not None:
            size = head_size(uri)
            return RemoteMetadata(size) if size is not None else None
        return _default_head_metadata(uri)

    with _stage_lock(plan.stage_dir, lock_timeout):
        return _apply_locked(plan, head, copy_object or _default_copy, disk_free or _default_disk_free, emit)


def _apply_locked(
    plan: BamStagePlan,
    head: Callable[[str], RemoteMetadata | None],
    copy_fn: Callable[[str, Path], None],
    free_fn: Callable[[Path], int],
    emit: Callable[[str], None],
) -> dict[str, str]:
    pending: list[tuple[str, Path, RemoteMetadata, dict[str, object]]] = []
    rewritten = dict(plan.rewritten)
    for uri in plan.dest_by_uri:
        metadata = head(uri)
        if metadata is None or metadata.size <= 0:
            raise RuntimeError("Could not determine a positive size of remote BAM")
        source = staged_alignment_dest(plan.stage_dir, uri)
        generation = hashlib.sha256(json.dumps(asdict(metadata), sort_keys=True).encode()).hexdigest()
        if not metadata.reusable:
            emit("[bam-stage] remote has no strong identity validator; cache reuse disabled")
            # Without a byte-identity validator even same-size objects may
            # differ. Give every transfer a new generation and never reuse.
            generation = uuid.uuid4().hex
        dest = source.parent / generation / source.name
        manifest: dict[str, object] = {"schema": 1, "source_hash": _uri_hash(uri), "remote": asdict(metadata)}
        valid = _cache_valid(dest, metadata, manifest, head, uri)
        if not valid and metadata.reusable and source.parent.exists():
            for candidate in sorted(source.parent.glob(f"{generation}-*")):
                if candidate.is_dir() and _cache_valid(candidate / source.name, metadata, manifest, head, uri):
                    dest = candidate / source.name
                    valid = True
                    break
        if not valid and dest.parent.exists():
            # Never mutate a published generation: a previous run may be
            # reading it. Invalid/incomplete generations get a fresh namespace.
            dest = source.parent / f"{generation}-{uuid.uuid4().hex}" / source.name
        for label, original in plan.rewritten.items():
            if original == str(plan.dest_by_uri[uri]):
                rewritten[label] = str(dest)
        if valid:
            emit(f"[bam-stage] reuse {dest} ({format_bytes(metadata.size)})")
        else:
            pending.append((uri, dest, metadata, manifest))

    need = required_free_bytes(metadata.size for _uri, _dest, metadata, _manifest in pending)
    free = free_fn(plan.stage_dir)
    if need > free:
        raise RuntimeError(
            "not enough free space to stage remote BAMs.\n"
            f"  dest: {plan.stage_dir}\n  need: {format_bytes(need)} (BAMs + headroom)\n"
            f"  free: {format_bytes(free)}\n"
            "Use --bam-stage-dir / RTM_BAM_STAGE_DIR on a larger volume, "
            "or disable with --no-bam-stage / RTM_BAM_STAGE=0."
        )

    def copy_pair(uri: str, dest: Path, metadata: RemoteMetadata, manifest: dict[str, object]) -> None:
        dest.parent.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".rtm-stage-", dir=dest.parent.parent) as tmp:
            staged_bam = Path(tmp) / dest.name
            emit(f"[bam-stage] copy alignment -> {dest} ({format_bytes(metadata.size)})")
            copy_fn(uri, staged_bam)
            if not local_copy_is_complete(staged_bam, metadata.size):
                raise RuntimeError(f"staged BAM size mismatch for {dest}")
            for idx_uri, idx_dest in zip(index_sidecar_uris(uri), index_sidecar_dests(dest)):
                staged_index = Path(tmp) / idx_dest.name
                try:
                    index_metadata = head(idx_uri)
                    copy_fn(idx_uri, staged_index)
                    if staged_index.is_file() and staged_index.stat().st_size > 0:
                        if index_metadata is not None and not local_copy_is_complete(staged_index, index_metadata.size):
                            raise ValueError("Index size does not match remote metadata")
                        if head(idx_uri) != index_metadata:
                            raise RuntimeError("Remote index changed during staging")
                        break
                except Exception:  # noqa: BLE001 - try the next naming convention; do not expose signed URLs
                    pass
                staged_index.unlink(missing_ok=True)
            else:
                raise RuntimeError(f"Could not stage BAI/CSI/CRAI next to {dest}")
            # Detect changes during the transfer, including same-size updates.
            if head(uri) != metadata:
                raise RuntimeError(f"Remote BAM changed during staging for {dest}; retry the run")
            saved = dict(manifest)
            saved.update(index_name=staged_index.name, alignment_file=_file_fingerprint(staged_bam),
                         index_file=_file_fingerprint(staged_index),
                         index_remote=asdict(index_metadata) if index_metadata is not None else None)
            (Path(tmp) / "manifest.json").write_text(json.dumps(saved, sort_keys=True) + "\n", encoding="utf-8")
            # Publish the complete directory (alignment, index and manifest)
            # in one rename. No reader can observe half of the pair.
            Path(tmp).rename(dest.parent)

    if pending:
        with ThreadPoolExecutor(max_workers=min(4, len(pending))) as pool:
            futures = [pool.submit(copy_pair, *item) for item in pending]
            for future in as_completed(futures):
                future.result()
    emit(f"[bam-stage] ready under {plan.stage_dir}")
    return rewritten


def write_env_file(path: str | Path, rewritten: Mapping[str, str]) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for key in ("DISEASE_BAM", "CONTROL_BAM", "DISEASE_MATE_BAM", "CONTROL_MATE_BAM"):
        val = rewritten.get(key, "")
        lines.append(f"{key}={_shell_single_quote(val)}")
    with tempfile.TemporaryDirectory(prefix=".rtm-env-", dir=dest.parent) as tmp:
        staged = Path(tmp) / "env"
        staged.write_text("\n".join(lines) + "\n", encoding="utf-8")
        staged.replace(dest)


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
    p.add_argument("--lock-timeout", type=float, default=3600, help="Seconds to wait for the staging lock (default: 3600).")
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
    if plan.should_stage and not args.apply:
        raise ValueError("--apply is required to publish staged paths; a plan has no verified local files")
    rewritten = apply_bam_stage(plan, lock_timeout=args.lock_timeout) if args.apply else dict(plan.rewritten)
    write_env_file(args.out_env, rewritten)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
