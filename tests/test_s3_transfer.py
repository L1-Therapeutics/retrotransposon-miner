"""Unit tests for high-concurrency S3 transfer settings (no live S3)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from retro_miner.s3_transfer import (
    DEFAULT_S3_MAX_CONCURRENCY,
    DEFAULT_S3_MAX_QUEUE_SIZE,
    DEFAULT_S3_MULTIPART_CHUNKSIZE,
    DEFAULT_S3_MULTIPART_THRESHOLD,
    boto3_available,
    boto3_transfer_config,
    copy_s3_uri,
    download_s3_uri,
    process_local_aws_config,
    s3_max_concurrency,
    s3_transfer_settings,
    split_s3_uri,
    write_process_local_aws_config,
)


def test_default_s3_transfer_settings() -> None:
    knobs = s3_transfer_settings()
    assert knobs["max_concurrency"] == DEFAULT_S3_MAX_CONCURRENCY == 64
    assert knobs["multipart_chunksize"] == DEFAULT_S3_MULTIPART_CHUNKSIZE == 64 * 1024 * 1024
    assert knobs["multipart_threshold"] == DEFAULT_S3_MULTIPART_THRESHOLD == 64 * 1024 * 1024
    assert knobs["max_queue_size"] == DEFAULT_S3_MAX_QUEUE_SIZE == 10_000
    assert "max_bandwidth" not in knobs


def test_s3_max_concurrency_env_override(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("RTM_S3_MAX_CONCURRENCY", raising=False)
    assert s3_max_concurrency() == 64
    monkeypatch.setenv("RTM_S3_MAX_CONCURRENCY", "16")
    assert s3_max_concurrency() == 16
    assert s3_transfer_settings()["max_concurrency"] == 16
    assert s3_max_concurrency("8") == 8
    with pytest.raises(ValueError, match=">= 1"):
        s3_max_concurrency("0")
    with pytest.raises(ValueError, match="positive"):
        s3_max_concurrency("nope")


def test_split_s3_uri() -> None:
    assert split_s3_uri("s3://l1tx-data/public/foo.bam") == ("l1tx-data", "public/foo.bam")
    with pytest.raises(ValueError):
        split_s3_uri("https://example.com/foo.bam")
    with pytest.raises(ValueError):
        split_s3_uri("s3://bucket-only")


def test_write_process_local_aws_config_does_not_touch_home(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    home = tmp_path / "home"
    aws_dir = home / ".aws"
    aws_dir.mkdir(parents=True)
    global_cfg = aws_dir / "config"
    global_cfg.write_text("[default]\nregion = us-east-1\n", encoding="utf-8")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("AWS_CONFIG_FILE", raising=False)
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-west-2")
    dest = tmp_path / "cfg" / "config"
    write_process_local_aws_config(dest)
    text = dest.read_text(encoding="utf-8")
    assert "max_concurrent_requests = 64" in text
    assert "multipart_chunksize = 64MB" in text
    assert "multipart_threshold = 64MB" in text
    assert "max_queue_size = 10000" in text
    assert "max_bandwidth" not in text
    assert "region = us-west-2" in text
    assert global_cfg.read_text(encoding="utf-8") == "[default]\nregion = us-east-1\n"


def test_process_local_aws_config_leaves_process_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("AWS_CONFIG_FILE", raising=False)
    with process_local_aws_config() as env:
        cfg = Path(env["AWS_CONFIG_FILE"])
        assert cfg.is_file()
        assert "max_concurrent_requests = 64" in cfg.read_text(encoding="utf-8")
        assert "AWS_CONFIG_FILE" not in os.environ
    assert not cfg.exists()


def test_boto3_transfer_config_values() -> None:
    if not boto3_available():
        pytest.skip("boto3 is not installed")
    cfg = boto3_transfer_config()
    assert cfg.max_concurrency == 64
    assert cfg.multipart_chunksize == 64 * 1024 * 1024
    assert cfg.multipart_threshold == 64 * 1024 * 1024
    assert cfg.max_io_queue == 10_000
    assert getattr(cfg, "max_bandwidth", None) is None


def test_download_s3_uri_uses_boto3_when_available(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dest = tmp_path / "out.bam"
    calls: list[tuple[str, str, str]] = []

    class _Cfg:
        max_concurrency = 64
        max_bandwidth = None

    class _Client:
        def download_file(self, bucket, key, filename, Config=None):
            calls.append((bucket, key, filename))
            Path(filename).write_bytes(b"ok")
            assert Config is not None
            assert Config.max_concurrency == 64
            assert getattr(Config, "max_bandwidth", None) is None

    monkeypatch.setattr("retro_miner.s3_transfer.boto3_available", lambda: True)
    monkeypatch.setattr("retro_miner.s3_transfer.boto3_transfer_config", lambda settings=None: _Cfg())
    monkeypatch.setattr("retro_miner.s3_transfer._boto3_client", lambda: _Client())
    download_s3_uri("s3://l1tx-data/public/foo.bam", dest)
    assert calls == [("l1tx-data", "public/foo.bam", str(dest))]
    assert dest.read_bytes() == b"ok"


def test_copy_s3_uri_server_side_uses_boto3(monkeypatch: pytest.MonkeyPatch) -> None:
    copies: list[tuple[dict[str, str], str, str]] = []

    class _Cfg:
        max_concurrency = 16
        max_bandwidth = None

    class _Client:
        def copy(self, src, bucket, key, Config=None):
            copies.append((src, bucket, key))
            assert Config.max_concurrency == 16

    monkeypatch.setattr("retro_miner.s3_transfer.boto3_available", lambda: True)
    monkeypatch.setattr("retro_miner.s3_transfer.boto3_transfer_config", lambda settings=None: _Cfg())
    monkeypatch.setattr("retro_miner.s3_transfer._boto3_client", lambda: _Client())
    monkeypatch.setenv("RTM_S3_MAX_CONCURRENCY", "16")
    copy_s3_uri("s3://src-bkt/a.bam", "s3://dst-bkt/b.bam")
    assert copies == [({"Bucket": "src-bkt", "Key": "a.bam"}, "dst-bkt", "b.bam")]


def test_copy_s3_uri_local_routes_to_download(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[str] = []

    def _fake_download(uri: str, dest) -> None:
        seen.append(uri)
        Path(dest).write_bytes(b"x")

    monkeypatch.setattr("retro_miner.s3_transfer.download_s3_uri", _fake_download)
    dest = tmp_path / "local.bam"
    copy_s3_uri("s3://bucket/key.bam", dest)
    assert seen == ["s3://bucket/key.bam"]
