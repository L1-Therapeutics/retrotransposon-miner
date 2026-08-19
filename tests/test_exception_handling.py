"""Tests for exception handling in retro_miner helper functions."""

from __future__ import annotations

import pytest

from retro_miner.mei_support import _extract_float_from_info
from retro_miner.local_assembly import _parse_existing_manifest


def test_extract_float_from_info_valid():
    _extract_float_from_info(42.0)


def test_extract_float_from_info_int():
    result = _extract_float_from_info(42)
    assert result == 42.0


def test_extract_float_from_info_none():
    result = _extract_float_from_info(None)
    assert result == -1.0


def test_extract_float_from_info_tuple_empty():
    result = _extract_float_from_info(())
    assert result == -1.0


def test_extract_float_from_info_tuple_with_none():
    result = _extract_float_from_info((None,))
    assert result == -1.0


def test_extract_float_from_info_bad_string():
    result = _extract_float_from_info("not_a_number")
    assert result == -1.0


def test_extract_float_from_info_bad_tuple():
    result = _extract_float_from_info(("abc", "def"))
    assert result == -1.0


def test_extract_float_from_info_mixed_tuple():
    result = _extract_float_from_info((42.0, "abc"))
    assert result == -1.0


def test_parse_existing_manifest_valid():
    import json
    import tempfile
    from pathlib import Path

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump({"key": "value"}, f)
        f.flush()
        result = _parse_existing_manifest(Path(f.name))
    assert result == {"key": "value"}


def test_parse_existing_manifest_invalid_json():
    import tempfile
    from pathlib import Path

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        f.write("not valid json{{{")
        f.flush()
        result = _parse_existing_manifest(Path(f.name))
    assert result is None


def test_parse_existing_manifest_not_dict():
    import tempfile
    from pathlib import Path

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        f.write('["list", not, "dict"]')
        f.flush()
        result = _parse_existing_manifest(Path(f.name))
    assert result is None


def test_parse_existing_manifest_missing_file():
    from pathlib import Path

    result = _parse_existing_manifest(Path("/nonexistent/manifest.json"))
    assert result is None