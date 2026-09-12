#!/usr/bin/env python3
"""Callset Parity Comparator for retrotransposon-miner.

Validates 100% concordance between baseline and PR #32 candidate outputs
(VCF or TSV call files).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Repo-root sys.path resolution for standalone execution
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def parse_callset(filepath: Path) -> set[tuple[str, ...]]:
    """Parse a VCF/TSV callset file and return a set of call key tuples.

    Each call key is ``(CHROM, POS, REF, ALT, MEI_FAMILY)``.

    Args:
        filepath: Path to the VCF or TSV callset file.

    Returns:
        A set of call key tuples.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
    """
    calls: set[tuple[str, ...]] = set()
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with filepath.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                calls.add((parts[0], parts[1], parts[3], parts[4]))
    return calls


def compare_parity(baseline_path: Path, candidate_path: Path) -> dict[str, Any]:
    """Compare two callsets and return a parity summary dictionary.

    Args:
        baseline_path: Path to the baseline callset.
        candidate_path: Path to the candidate (PR #32) callset.

    Returns:
        A dict with keys ``concordance_pct``, ``total_baseline``,
        ``total_candidate``, ``shared``, ``missing_in_candidate``,
        ``extra_in_candidate``, and ``exact_match``.
    """
    base_calls = parse_callset(baseline_path)
    cand_calls = parse_callset(candidate_path)

    shared = base_calls.intersection(cand_calls)
    missing = base_calls - cand_calls
    extra = cand_calls - base_calls

    total_unique = len(base_calls.union(cand_calls))
    concordance = (len(shared) / total_unique * 100.0) if total_unique > 0 else 100.0

    return {
        "concordance_pct": round(concordance, 4),
        "total_baseline": len(base_calls),
        "total_candidate": len(cand_calls),
        "shared": len(shared),
        "missing_in_candidate": len(missing),
        "extra_in_candidate": len(extra),
        "exact_match": len(missing) == 0 and len(extra) == 0,
    }


def main() -> int:
    """Entry point: parse args, run parity check, print results, exit with status."""
    parser = argparse.ArgumentParser(description="Compare VCF/TSV callset parity for RTM.")
    parser.add_argument("--baseline", type=Path, required=True, help="Baseline callset path")
    parser.add_argument("--candidate", type=Path, required=True, help="PR #32 candidate callset path")
    parser.add_argument("--json", action="store_true", help="Output raw JSON summary")
    args = parser.parse_args()

    results = compare_parity(args.baseline, args.candidate)

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("=== RTM Callset Parity Analysis ===")
        print(f"Concordance Rate : {results['concordance_pct']}%")
        print(f"Baseline Calls   : {results['total_baseline']}")
        print(f"Candidate Calls  : {results['total_candidate']}")
        print(f"Shared Calls     : {results['shared']}")
        print(f"Divergence       : Missing={results['missing_in_candidate']}, Extra={results['extra_in_candidate']}")
        print(f"Exact Parity     : {'PASSED' if results['exact_match'] else 'FAILED'}")

    return 0 if results["exact_match"] else 1


if __name__ == "__main__":
    sys.exit(main())
