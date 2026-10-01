"""Regression tests for remote BAM Content-Length probing and AWS region parsing.

Both areas were previously exercised only indirectly (via injected ``head_size``
lambdas and the happy-path ``region =`` config), which is how a
mis-indented ``request =`` assignment and a prefix-based region match survived.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

from retro_miner import bam_stage, s3_transfer


class _Resp:
    def __init__(self, headers):
        self.headers = headers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _Proc:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_head_size_https_issues_head_request(monkeypatch):
    """An https:// BAM must actually be probed; previously `request` was
    unbound inside the s3 branch, so urlopen was never called and the
    UnboundLocalError was swallowed, always returning None."""
    seen = {}

    def fake_urlopen(request, timeout=None):
        seen["url"] = request.full_url
        seen["method"] = request.get_method()
        return _Resp({"Content-Length": "4096"})

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    assert bam_stage._default_head_size("https://example.com/x.bam") == 4096
    assert seen == {"url": "https://example.com/x.bam", "method": "HEAD"}


def test_head_size_http_works_too(monkeypatch):
    monkeypatch.setattr(
        urllib.request,
        "urlopen",
        lambda request, timeout=None: _Resp({"Content-Length": "17"}),
    )
    assert bam_stage._default_head_size("http://example.com/x.bam") == 17


def test_head_size_non_http_uri_does_not_call_urlopen(monkeypatch):
    """A non-http(s) URI must short-circuit rather than build a Request."""
    calls = []

    def fake_urlopen(*a, **k):
        calls.append(a)
        raise AssertionError("urlopen must not be called for a non-http(s) URI")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    assert bam_stage._default_head_size("ftp://example.com/x.bam") is None
    assert calls == []


def test_head_size_missing_content_length_returns_none(monkeypatch):
    monkeypatch.setattr(
        urllib.request, "urlopen", lambda request, timeout=None: _Resp({})
    )
    assert bam_stage._default_head_size("https://example.com/x.bam") is None


def test_head_size_s3_uses_aws_cli(monkeypatch):
    """The s3:// branch must keep short-circuiting via the aws CLI."""

    def fail_urlopen(*a, **k):
        raise AssertionError("s3:// must not fall through to an HTTP HEAD")

    monkeypatch.setattr(urllib.request, "urlopen", fail_urlopen)
    monkeypatch.setattr(bam_stage.shutil, "which", lambda _n: "/usr/bin/aws")
    monkeypatch.setattr(
        bam_stage.subprocess, "run", lambda *a, **k: _Proc(stdout="2048\n")
    )
    assert bam_stage._default_head_size("s3://bucket/key.bam") == 2048


def _isolate_aws_env(monkeypatch, tmp_path, *, config_text: str) -> Path:
    """Point _aws_region at a throwaway config and an empty home, so the
    developer's real ~/.aws/config cannot satisfy the lookup."""
    home = tmp_path / "home"
    home.mkdir()
    cfg = tmp_path / "config"
    cfg.write_text(config_text, encoding="utf-8")
    for var in ("AWS_DEFAULT_REGION", "AWS_REGION"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("AWS_CONFIG_FILE", str(cfg))
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    return cfg


def test_aws_region_plain_key(tmp_path, monkeypatch):
    _isolate_aws_env(monkeypatch, tmp_path, config_text="[default]\nregion = us-east-1\n")
    assert s3_transfer._aws_region() == "us-east-1"


def test_aws_region_ignores_plural_regions_key(tmp_path, monkeypatch):
    """`regions` is a real AWS config key (space-separated list). A prefix
    match would wrongly return 'us-east-1 us-west-2' as the region."""
    _isolate_aws_env(
        monkeypatch, tmp_path, config_text="[default]\nregions = us-east-1 us-west-2\n"
    )
    assert s3_transfer._aws_region() == ""


def test_aws_region_ignores_region_prefixed_other_keys(tmp_path, monkeypatch):
    _isolate_aws_env(
        monkeypatch,
        tmp_path,
        config_text="[default]\nregion_id = not-a-region\noutput = json\n",
    )
    assert s3_transfer._aws_region() == ""


def test_aws_region_falls_through_to_real_region_after_similar_key(tmp_path, monkeypatch):
    _isolate_aws_env(
        monkeypatch,
        tmp_path,
        config_text="[default]\nregions = us-east-1 us-west-2\nregion = eu-west-1\n",
    )
    assert s3_transfer._aws_region() == "eu-west-1"


def test_aws_region_env_var_wins(tmp_path, monkeypatch):
    monkeypatch.setenv("AWS_DEFAULT_REGION", "ap-south-1")
    assert s3_transfer._aws_region() == "ap-south-1"
