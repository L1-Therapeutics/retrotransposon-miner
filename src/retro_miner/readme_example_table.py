"""Helpers for the README gold-example markdown table.

Root cause of past column shifts: fields like ``known_mei_polymorphism_id`` can
contain literal ``|`` (e.g. ``g1k:…|lr:…``), which splits markdown table cells.
Always format cells with ``markdown_table_cell`` and validate with
``assert_markdown_table_shape``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

EXAMPLE_SECTION_START = "Gold-tier calls from the SEQC2"
EXAMPLE_SECTION_END = "\n## Examples\n"

TABLE_COLS = [
    "chrom",
    "consensus_insertion_breakpoint_pos",
    "window_start",
    "window_end",
    "control_supporting_reads",
    "disease_supporting_reads",
    "sample_status_label",
    "consensus_tsd_seq",
    "consensus_poly_at_min_bp",
    "consensus_mei_family",
    "consensus_mei_subfamily",
    "known_mei_polymorphism_id",
    "known_mei_polymorphism_source",
    "consensus_insertion_orientation",
    "nested_in_same_MEI",
    "consensus_insertion_mei_span_full",
    "consensus_insertion_mei_5p_coord_full",
    "consensus_insertion_mei_3p_coord_full",
]

FILL_COLS = [
    "consensus_insertion_mei_span_full",
    "consensus_insertion_mei_5p_coord_full",
    "consensus_insertion_mei_3p_coord_full",
    "consensus_insertion_orientation",
    "nested_in_same_MEI",
    "consensus_mei_subfamily",
]


def markdown_table_cell(value: object) -> str:
    """Format one markdown table cell; never emit an unescaped ``|``."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float) and float(value).is_integer():
        text = str(int(value))
    else:
        text = str(value)
    text = text.replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    text = text.replace("|", ";")
    return text.strip()


def assert_markdown_table_shape(table_md: str) -> None:
    """Raise if any data row's cell count differs from the header."""
    rows = [ln for ln in table_md.splitlines() if ln.startswith("|")]
    if len(rows) < 2:
        raise ValueError("markdown table needs a header and separator")

    def n_cells(line: str) -> int:
        body = line.strip()
        if not (body.startswith("|") and body.endswith("|")):
            raise ValueError(f"malformed table row: {line[:80]!r}")
        return len([c for c in body[1:-1].split("|")])

    expected = n_cells(rows[0])
    for line in rows[2:]:
        got = n_cells(line)
        if got != expected:
            raise ValueError(
                f"markdown table column shift: expected {expected} cells, "
                f"got {got} in row starting {line[:100]!r}. "
                "Cell values must not contain unescaped '|' "
                "(use markdown_table_cell / scripts/update_readme_example_table.py)."
            )


def _clean_subfamily(value: object) -> str:
    text = "" if value is None or (isinstance(value, float) and pd.isna(value)) else str(value)
    return text.replace("_full#", "#")


def _prefer(row: pd.Series, col: str) -> object:
    cur = row.get(col)
    alt = row.get(f"{col}_fill")
    if col.endswith("_full") or "coord" in col:
        try:
            cv = float(cur) if pd.notna(cur) else 0.0
        except (TypeError, ValueError):
            cv = 0.0
        if cv > 0:
            return cur
        return alt if pd.notna(alt) else cur
    if pd.isna(cur) or str(cur).strip() in {"", "nan", "."}:
        return alt if pd.notna(alt) else cur
    return cur


def build_example_table_markdown(
    *,
    gold_review: Path,
    rank_index: Path | None = None,
    fill_gold_review: Path | None = None,
    top_n: int = 25,
    caption: str | None = None,
) -> str:
    """Build the README example table markdown (caption + table)."""
    gold = pd.read_csv(gold_review, sep="\t", low_memory=False)
    g = gold.loc[gold["analysis_stage_tier"].astype(str).str.lower() == "gold"].copy()
    if g.empty:
        raise ValueError(f"no gold rows in {gold_review}")

    keys = ["chrom", "window_start", "window_end"]
    for c in ("window_start", "window_end"):
        g[c] = pd.to_numeric(g[c], errors="coerce").fillna(0).astype(int)
    g["chrom"] = g["chrom"].astype(str)

    if rank_index is not None and Path(rank_index).exists():
        idx = pd.read_csv(rank_index, sep="\t")
        idx["chrom"] = idx["chrom"].astype(str)
        for c in ("window_start", "window_end"):
            idx[c] = pd.to_numeric(idx[c], errors="coerce").fillna(0).astype(int)
        selected = idx.head(int(top_n)).merge(g, on=keys, how="left", suffixes=("_idx", ""))
    else:
        selected = g.head(int(top_n)).copy()

    if fill_gold_review is not None and Path(fill_gold_review).exists():
        fill = pd.read_csv(fill_gold_review, sep="\t", low_memory=False)
        fill["chrom"] = fill["chrom"].astype(str)
        for c in ("window_start", "window_end"):
            fill[c] = pd.to_numeric(fill[c], errors="coerce").fillna(0).astype(int)
        keep = keys + [c for c in FILL_COLS if c in fill.columns]
        f = fill.loc[:, keep].drop_duplicates(keys, keep="first")
        selected = selected.merge(f, on=keys, how="left", suffixes=("", "_fill"))
        for c in FILL_COLS:
            if c in selected.columns:
                selected[c] = [_prefer(r, c) for _, r in selected.iterrows()]

    for full_c, panel_c in (
        ("consensus_insertion_mei_span_full", "consensus_insertion_mei_span"),
        ("consensus_insertion_mei_5p_coord_full", "consensus_insertion_mei_5p_coord"),
        ("consensus_insertion_mei_3p_coord_full", "consensus_insertion_mei_3p_coord"),
    ):
        if full_c in selected.columns and panel_c in selected.columns:
            full = pd.to_numeric(selected[full_c], errors="coerce").fillna(0)
            panel = pd.to_numeric(selected[panel_c], errors="coerce")
            selected.loc[full <= 0, full_c] = panel.loc[full <= 0]

    if "consensus_mei_subfamily" in selected.columns:
        selected["consensus_mei_subfamily"] = selected["consensus_mei_subfamily"].map(_clean_subfamily)

    n_gold = int(len(g))
    if caption is None:
        caption = (
            "Gold-tier calls from the SEQC2 tumor/normal chr22 annotate run after "
            f"`COMPLEX_INS` demotion (top {int(top_n)} of n={n_gold} by review rank)."
        )

    lines = [
        caption,
        "",
        "| " + " | ".join(TABLE_COLS) + " |",
        "| " + " | ".join(["---"] * len(TABLE_COLS)) + " |",
    ]
    for _, row in selected.iterrows():
        cells = [markdown_table_cell(row.get(c)) for c in TABLE_COLS]
        lines.append("| " + " | ".join(cells) + " |")
    table_md = "\n".join(lines) + "\n"
    assert_markdown_table_shape(table_md)
    return table_md


def replace_readme_example_table(readme_text: str, table_md: str) -> str:
    start = readme_text.index(EXAMPLE_SECTION_START)
    end = readme_text.index(EXAMPLE_SECTION_END)
    return readme_text[:start] + table_md.rstrip() + "\n\n" + readme_text[end + 1 :]
