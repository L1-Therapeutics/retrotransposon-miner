"""High-concurrency S3 multipart copies (process-local AWS config).

AWS CLI defaults (~10 × 8 MiB) cap an r6i.4xlarge near ~75 MiB/s. This helper
uses 64 concurrent 64 MiB parts instead. Override concurrency with
``RTM_S3_MAX_CONCURRENCY``. ``max_bandwidth`` is never set.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

DEFAULT_S3_MAX_CONCURRENCY = 64
DEFAULT_S3_MULTIPART_CHUNKSIZE = 64 * 1024 * 1024
DEFAULT_S3_MULTIPART_THRESHOLD = 64 * 1024 * 1024
DEFAULT_S3_MAX_QUEUE_SIZE = 10_000
_MIB = 1024 * 1024


def split_s3_uri(uri: str) -> tuple[str, str]:
    text = (uri or "").strip()
    if not text.startswith("s3://"):
        raise ValueError(f"Not an s3 URI: {uri}")
    rest = text[len("s3://") :].lstrip("/")
    bucket, _, key = rest.partition("/")
    if not bucket or not key:
        raise ValueError(f"s3 URI must include bucket and key: {uri}")
    return bucket, key


def s3_max_concurrency(raw: str | None = None) -> int:
    text = raw if raw is not None else os.environ.get("RTM_S3_MAX_CONCURRENCY", "")
    if text is None or not str(text).strip():
        return DEFAULT_S3_MAX_CONCURRENCY
    try:
        value = int(str(text).strip())
    except ValueError as exc:
        raise ValueError(f"RTM_S3_MAX_CONCURRENCY must be a positive integer, got {text!r}") from exc
    if value < 1:
        raise ValueError(f"RTM_S3_MAX_CONCURRENCY must be >= 1, got {value}")
    return value


def s3_transfer_settings(*, max_concurrency: int | None = None) -> dict[str, int]:
    """Return TransferConfig / AWS CLI knobs. Does not include max_bandwidth."""
    return {
        "max_concurrency": int(max_concurrency) if max_concurrency is not None else s3_max_concurrency(),
        "multipart_chunksize": DEFAULT_S3_MULTIPART_CHUNKSIZE,
        "multipart_threshold": DEFAULT_S3_MULTIPART_THRESHOLD,
        "max_queue_size": DEFAULT_S3_MAX_QUEUE_SIZE,
    }


def _aws_region() -> str:
    for key in ("AWS_DEFAULT_REGION", "AWS_REGION"):
        val = (os.environ.get(key) or "").strip()
        if val:
            return val
    cfg_path = (os.environ.get("AWS_CONFIG_FILE") or "").strip()
    candidates = [Path(cfg_path)] if cfg_path else []
    candidates.append(Path.home() / ".aws" / "config")
    for path in candidates:
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        in_default = False
        for raw_line in text.splitlines():
            line = raw_line.strip()
            if line.startswith("[") and line.endswith("]"):
                in_default = line.lower() in {"[default]", "[profile default]"}
                continue
            if in_default and line.lower().startswith("region"):
                _, _, value = line.partition("=")
                value = value.strip()
                if value:
                    return value
    return ""


def write_process_local_aws_config(
    path: str | Path,
    *,
    settings: Mapping[str, int] | None = None,
    region: str | None = None,
) -> Path:
    """Write a temp AWS CLI config. Does not modify ``~/.aws/config``."""
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    knobs = s3_transfer_settings() if settings is None else dict(settings)
    chunk_mib = max(1, int(knobs["multipart_chunksize"]) // _MIB)
    thresh_mib = max(1, int(knobs["multipart_threshold"]) // _MIB)
    resolved_region = (region if region is not None else _aws_region()).strip()
    lines = ["[default]"]
    if resolved_region:
        lines.append(f"region = {resolved_region}")
    lines.extend(
        [
            "s3 =",
            f"    max_concurrent_requests = {int(knobs['max_concurrency'])}",
            f"    multipart_chunksize = {chunk_mib}MB",
            f"    multipart_threshold = {thresh_mib}MB",
            f"    max_queue_size = {int(knobs['max_queue_size'])}",
            "",
        ]
    )
    profile = (os.environ.get("AWS_PROFILE") or "").strip()
    if profile and profile.lower() != "default":
        lines.append(f"[profile {profile}]")
        if resolved_region:
            lines.append(f"region = {resolved_region}")
        lines.extend(
            [
                "s3 =",
                f"    max_concurrent_requests = {int(knobs['max_concurrency'])}",
                f"    multipart_chunksize = {chunk_mib}MB",
                f"    multipart_threshold = {thresh_mib}MB",
                f"    max_queue_size = {int(knobs['max_queue_size'])}",
                "",
            ]
        )
    dest.write_text("\n".join(lines), encoding="utf-8")
    return dest


@contextmanager
def process_local_aws_config(
    *,
    settings: Mapping[str, int] | None = None,
) -> Iterator[dict[str, str]]:
    """Yield an env mapping with ``AWS_CONFIG_FILE`` pointing at a temp config."""
    with tempfile.TemporaryDirectory(prefix="rtm_aws_cfg_") as tmp:
        cfg = write_process_local_aws_config(Path(tmp) / "config", settings=settings)
        env = os.environ.copy()
        env["AWS_CONFIG_FILE"] = str(cfg)
        region = _aws_region()
        if region and not (env.get("AWS_DEFAULT_REGION") or "").strip():
            env["AWS_DEFAULT_REGION"] = region
        yield env


def boto3_available() -> bool:
    try:
        import boto3  # noqa: F401
    except ImportError:
        return False
    return True


def boto3_transfer_config(settings: Mapping[str, int] | None = None) -> Any:
    from boto3.s3.transfer import TransferConfig

    knobs = s3_transfer_settings() if settings is None else dict(settings)
    return TransferConfig(
        multipart_threshold=int(knobs["multipart_threshold"]),
        max_concurrency=int(knobs["max_concurrency"]),
        multipart_chunksize=int(knobs["multipart_chunksize"]),
        max_io_queue=int(knobs["max_queue_size"]),
        use_threads=True,
    )


def _boto3_client():
    import boto3

    kwargs: dict[str, str] = {}
    region = _aws_region()
    if region:
        kwargs["region_name"] = region
    return boto3.client("s3", **kwargs)


def _run_aws_cli(args: list[str], *, env: Mapping[str, str] | None = None) -> None:
    if shutil.which("aws") is None:
        raise RuntimeError(
            "aws CLI is required for S3 copies when boto3 is unavailable. "
            "On EC2, attach an instance profile instead of copying local AWS keys."
        )
    cmd = ["aws", *args]
    proc = subprocess.run(
        cmd,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=dict(env) if env else None,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{err}")


def download_s3_uri(s3_uri: str, dest: str | Path) -> None:
    """Copy an S3 object to a local path using boto3 TransferConfig, else aws CLI."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if boto3_available():
        bucket, key = split_s3_uri(s3_uri)
        _boto3_client().download_file(bucket, key, str(dest_path), Config=boto3_transfer_config())
        return
    with process_local_aws_config() as env:
        _run_aws_cli(["s3", "cp", s3_uri, str(dest_path)], env=env)


def copy_s3_uri(src_uri: str, dst: str | Path) -> None:
    """Copy S3→S3 or S3→local with the high-concurrency transfer settings."""
    dst_text = str(dst)
    if dst_text.startswith("s3://"):
        if boto3_available():
            src_bucket, src_key = split_s3_uri(src_uri)
            dst_bucket, dst_key = split_s3_uri(dst_text)
            _boto3_client().copy(
                {"Bucket": src_bucket, "Key": src_key},
                dst_bucket,
                dst_key,
                Config=boto3_transfer_config(),
            )
            return
        with process_local_aws_config() as env:
            _run_aws_cli(["s3", "cp", src_uri, dst_text], env=env)
        return
    download_s3_uri(src_uri, dst)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="High-concurrency S3 object copy.")
    p.add_argument("src", help="s3://bucket/key")
    p.add_argument("dest", help="Local path or s3://bucket/key")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    copy_s3_uri(args.src, args.dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
