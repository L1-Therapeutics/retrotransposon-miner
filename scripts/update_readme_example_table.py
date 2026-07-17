#!/usr/bin/env python3
"""Refresh the README gold-example markdown table from an annotate run.

Always use this script (or ``retro_miner.readme_example_table``) so cell values
with literal ``|`` cannot shift markdown columns.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from retro_miner.readme_example_table import (  # noqa: E402
    EXAMPLE_SECTION_END,
    EXAMPLE_SECTION_START,
    assert_markdown_table_shape,
    build_example_table_markdown,
    replace_readme_example_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gold-review-tsv",
        type=Path,
        required=True,
        help="candidate_loci.mei.gold_review.tsv from the annotate run",
    )
    parser.add_argument(
        "--rank-index-tsv",
        type=Path,
        default=None,
        help="Optional read_architecture_index.tsv (defines review rank order)",
    )
    parser.add_argument(
        "--fill-gold-review-tsv",
        type=Path,
        default=None,
        help="Optional prior gold_review for span/orientation when re-label left blanks",
    )
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument(
        "--readme",
        type=Path,
        default=_REPO_ROOT / "README.md",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print table to stdout instead of rewriting README.md",
    )
    args = parser.parse_args()

    table_md = build_example_table_markdown(
        gold_review=args.gold_review_tsv,
        rank_index=args.rank_index_tsv,
        fill_gold_review=args.fill_gold_review_tsv,
        top_n=args.top_n,
    )
    if args.dry_run:
        print(table_md, end="")
        return

    readme = Path(args.readme)
    updated = replace_readme_example_table(readme.read_text(), table_md)
    start = updated.index(EXAMPLE_SECTION_START)
    end = updated.index(EXAMPLE_SECTION_END)
    assert_markdown_table_shape(updated[start:end])
    readme.write_text(updated)
    print(f"updated {readme} top_n={args.top_n}")


if __name__ == "__main__":
    main()
