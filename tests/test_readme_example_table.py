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


def test_readme_example_is_the_chr22_classifier_vcf():
    from pathlib import Path

    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text()
    start = text.index("## Example Variant Calls")
    end = text.index("\n## Examples\n")
    section = text[start:end]
    assert "#CHROM\tPOS\tID\t" in section
    assert section.count("\nchr22\t") == 21
    assert "Gold-tier calls from the SEQC2" not in section
