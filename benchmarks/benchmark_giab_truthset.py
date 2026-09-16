#!/usr/bin/env python3
"""Genome in a Bottle (GIAB) HG002 Structural Variant Benchmark Harness.

Evaluates retrotransposon-miner callset precision, recall, and F1-score against
truth sets within defined positional window tolerances (default ±50 bp).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_vcf_or_tsv_locus(file_path: Path) -> list[tuple[str, int, str]]:
    """Parse candidate callset loci from VCF or TSV file."""
    loci = []
    content = file_path.read_text().splitlines()
    for line in content:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            chrom = parts[0]
            try:
                pos = int(parts[1])
                mei_type = parts[4] if len(parts) > 4 else "MEI"
                loci.append((chrom, pos, mei_type))
            except ValueError:
                continue
    return loci


def evaluate_benchmark_parity(
    truth_loci: list[tuple[str, int, str]],
    candidate_loci: list[tuple[str, int, str]],
    window_bp: int = 50,
) -> dict[str, float]:
    """Calculate Precision, Recall, and F1-score between candidate and truth set."""
    true_positives = 0
    matched_candidates = set()

    for t_chrom, t_pos, _ in truth_loci:
        for idx, (c_chrom, c_pos, _) in enumerate(candidate_loci):
            if idx in matched_candidates:
                continue
            if t_chrom == c_chrom and abs(t_pos - c_pos) <= window_bp:
                true_positives += 1
                matched_candidates.add(idx)
                break

    fn = len(truth_loci) - true_positives
    fp = len(candidate_loci) - true_positives

    precision = true_positives / len(candidate_loci) if candidate_loci else 0.0
    recall = true_positives / len(truth_loci) if truth_loci else 0.0
    f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "true_positives": true_positives,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1_score": round(f1, 4),
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark MEI calls against GIAB truth set.")
    parser.add_argument("--truth", type=Path, required=True, help="Path to GIAB truth VCF/TSV")
    parser.add_argument("--candidate", type=Path, required=True, help="Path to candidate VCF/TSV")
    parser.add_argument("--window", type=int, default=50, help="Matching tolerance window in bp")
    args = parser.parse_args()

    if not args.truth.exists() or not args.candidate.exists():
        print("[-] Error: Input benchmark files do not exist.", file=sys.stderr)
        sys.exit(1)

    truth_loci = parse_vcf_or_tsv_locus(args.truth)
    cand_loci = parse_vcf_or_tsv_locus(args.candidate)

    metrics = evaluate_benchmark_parity(truth_loci, cand_loci, window_bp=args.window)
    print(f"[+] GIAB Benchmark Results (Window ±{args.window} bp):")
    print(f"  |-- Precision : {metrics['precision']:.4f}")
    print(f"  |-- Recall    : {metrics['recall']:.4f}")
    print(f"  |-- F1-Score  : {metrics['f1_score']:.4f}")
    print(f"  |-- TP: {metrics['true_positives']} | FP: {metrics['false_positives']} | FN: {metrics['false_negatives']}")


if __name__ == "__main__":
    main()
