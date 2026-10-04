"""High-concurrency S3 multipart copies (process-local AWS config).

AWS CLI defaults (~10 × 8 MiB) cap an m7i.8xlarge near ~75 MiB/s. This helper
uses 64 concurrent 64 MiB parts instead. Override concurrency with
``RTM_S3_MAX_CONCURRENCY``. ``max_bandwidth`` is never set.
"""

from __future__ import annotations

import argparse
import configparser
import os
import shutil
import subprocess
import sys
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
    concurrency = s3_max_concurrency(str(max_concurrency) if max_concurrency is not None else None)
    return {
        "max_concurrency": concurrency,
        "multipart_chunksize": DEFAULT_S3_MULTIPART_CHUNKSIZE,
        "multipart_threshold": DEFAULT_S3_MULTIPART_THRESHOLD,
        "max_queue_size": DEFAULT_S3_MAX_QUEUE_SIZE,
    }


def _aws_profile_section() -> str:
    profile = (os.environ.get("AWS_DEFAULT_PROFILE") or os.environ.get("AWS_PROFILE") or "default").strip()
    return "default" if not profile or profile == "default" else f"profile {profile}"


def _read_aws_config() -> configparser.ConfigParser:
    # Disable interpolation: credential_process commands and other settings
    # may contain literal percent signs. Do not fall back to a different
    # config when the caller explicitly selected AWS_CONFIG_FILE.
    config = configparser.ConfigParser(interpolation=None)
    raw_path = (os.environ.get("AWS_CONFIG_FILE") or "").strip()
    path = Path(raw_path).expanduser() if raw_path else Path.home() / ".aws" / "config"
    if path.is_file():
        with path.open(encoding="utf-8") as handle:
            config.read_file(handle)
    return config


def _aws_region() -> str:
    for key in ("AWS_DEFAULT_REGION", "AWS_REGION"):
        val = (os.environ.get(key) or "").strip()
        if val:
            return val
    return _read_aws_config().get(_aws_profile_section(), "region", fallback="").strip()


def write_process_local_aws_config(
    path: str | Path,
    *,
    settings: Mapping[str, int] | None = None,
    region: str | None = None,
) -> Path:
    """Copy AWS config and override transfer knobs without changing authentication.

    Role, credential-process and SSO profiles are retained. The source config
    (including a caller-selected ``AWS_CONFIG_FILE``) is never rewritten.
    """
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    knobs = s3_transfer_settings() if settings is None else dict(settings)
    chunk_mib = max(1, int(knobs["multipart_chunksize"]) // _MIB)
    thresh_mib = max(1, int(knobs["multipart_threshold"]) // _MIB)
    resolved_region = (region if region is not None else _aws_region()).strip()
    config = _read_aws_config()
    selected = _aws_profile_section()
    transfer_options = {
        "max_concurrent_requests": str(s3_max_concurrency(str(knobs["max_concurrency"]))),
        "multipart_chunksize": f"{chunk_mib}MB",
        "multipart_threshold": f"{thresh_mib}MB",
        "max_queue_size": str(int(knobs["max_queue_size"])),
    }
    for section in dict.fromkeys(("default", selected)):
        if not config.has_section(section):
            config.add_section(section)
        if section == selected and resolved_region:
            config.set(section, "region", resolved_region)
        # AWS stores nested S3 options as an indented multiline value. Keep
        # caller options such as addressing_style and endpoint behavior.
        options: dict[str, str] = {}
        for line in config.get(section, "s3", fallback="").splitlines():
            key, sep, value = line.strip().partition("=")
            if sep:
                options[key.strip()] = value.strip()
        options.update(transfer_options)
        config.set(section, "s3", "\n" + "\n".join(f"{key} = {value}" for key, value in options.items()))
    with dest.open("w", encoding="utf-8") as handle:
        config.write(handle)
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
    """Download with boto3 or the AWS CLI, atomically publishing on success.

    Failed transfers leave an existing destination untouched. Temporary files
    stay on the destination filesystem so replacement does not copy the object
    a second time. Empty objects are valid for this generic transfer helper.
    """
    bucket, key = split_s3_uri(s3_uri)
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".rtm-s3-", dir=dest_path.parent) as tmp:
        staged = Path(tmp) / dest_path.name
        if boto3_available():
            _boto3_client().download_file(bucket, key, str(staged), Config=boto3_transfer_config())
        else:
            with process_local_aws_config() as env:
                _run_aws_cli(["s3", "cp", s3_uri, str(staged)], env=env)
        if not staged.is_file():
            raise RuntimeError(f"S3 download output is missing for {dest_path}")
        staged.replace(dest_path)


def copy_s3_uri(src_uri: str, dst: str | Path) -> None:
    """Copy S3→S3 or S3→local with the high-concurrency transfer settings."""
    src_bucket, src_key = split_s3_uri(src_uri)
    dst_text = str(dst)
    if dst_text.startswith("s3://"):
        dst_bucket, dst_key = split_s3_uri(dst_text)
        if boto3_available():
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
