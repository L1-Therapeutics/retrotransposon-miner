"""Tests for scripts/generate_benchmark_report.py."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "generate_benchmark_report.py"


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
    )


def _write_json(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_valid_json_generates_markdown_table(tmp_path: Path) -> None:
    json_in = tmp_path / "bench.json"
    md_out = tmp_path / "report.md"
    _write_json(
        json_in,
        {
            "read_count": 1000,
            "runtime_seconds": 12.345,
            "peak_memory_mb": 42.5,
        },
    )

    result = _run("--json-input", str(json_in), "--markdown-out", str(md_out))
    assert result.returncode == 0, result.stderr
    assert md_out.exists()

    markdown = md_out.read_text(encoding="utf-8")
    assert "## Benchmark Results" in markdown
    assert "PR #32 One-Pass Sweep" in markdown
    assert "Legacy 2-Pass Sort" in markdown
    assert "Peak Memory (MB)" in markdown
    assert "Runtime (sec)" in markdown
    assert "1,000" in markdown  # formatted read_count
    assert "42.50 MB" in markdown
    assert "12.345s" in markdown


def test_malformed_json_exits_nonzero(tmp_path: Path) -> None:
    json_in = tmp_path / "bad.json"
    json_in.write_text("not valid json", encoding="utf-8")
    md_out = tmp_path / "report.md"

    result = _run("--json-input", str(json_in), "--markdown-out", str(md_out))
    assert result.returncode != 0
    assert "malformed JSON" in result.stderr


def test_missing_keys_exits_nonzero(tmp_path: Path) -> None:
    json_in = tmp_path / "partial.json"
    _write_json(json_in, {"read_count": 10})
    md_out = tmp_path / "report.md"

    result = _run("--json-input", str(json_in), "--markdown-out", str(md_out))
    assert result.returncode != 0
    assert "missing required keys" in result.stderr


def test_missing_json_input_exits_nonzero(tmp_path: Path) -> None:
    md_out = tmp_path / "report.md"
    result = _run("--json-input", str(tmp_path / "nonexistent.json"), "--markdown-out", str(md_out))
    assert result.returncode != 0
    assert "JSON input not found" in result.stderr
