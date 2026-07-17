"""README example table cells must never contain raw '|' (markdown column shift)."""

from __future__ import annotations

import pytest

from retro_miner.readme_example_table import (
    assert_markdown_table_shape,
    markdown_table_cell,
)


def test_markdown_table_cell_replaces_pipes_and_newlines():
    assert markdown_table_cell("g1k:a|lr:b") == "g1k:a;lr:b"
    assert markdown_table_cell("line1\nline2") == "line1 line2"
    assert markdown_table_cell(None) == ""
    assert markdown_table_cell(1386.0) == "1386"


def test_assert_markdown_table_shape_catches_pipe_split_rows():
    good = (
        "| a | b | c |\n"
        "| --- | --- | --- |\n"
        "| x | g1k:a;lr:b | z |\n"
    )
    assert_markdown_table_shape(good)

    bad = (
        "| a | b | c |\n"
        "| --- | --- | --- |\n"
        "| x | g1k:a|lr:b | z |\n"
    )
    with pytest.raises(ValueError, match="column shift"):
        assert_markdown_table_shape(bad)


def test_current_readme_example_table_has_stable_columns():
    from pathlib import Path

    from retro_miner.readme_example_table import (
        EXAMPLE_SECTION_END,
        EXAMPLE_SECTION_START,
    )

    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text()
    start = text.index(EXAMPLE_SECTION_START)
    end = text.index(EXAMPLE_SECTION_END)
    assert_markdown_table_shape(text[start:end])
