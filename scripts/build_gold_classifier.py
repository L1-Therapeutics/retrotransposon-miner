#!/usr/bin/env python3
"""Build a gold-call training set and fit the gradient-boosted score.

Each --sample argument is NAME=DIR, where DIR contains per-chromosome
candidate_loci.mei.gold_review.tsv tables. Those tables are concatenated and
re-sorted with the single-chromosome review priority. Positives are 1000
Genomes overlaps in the high-rate prefix of that genome-wide list. Known
overlaps below the prefix are excluded. Negatives are a family-matched sample
of the other gold calls below it.

Example:
  python scripts/build_gold_classifier.py \
    --sample HG00100=/path/hg00100_germline_wgs \
    --sample NA18525=/path/NA18525_germline_wgs \
    --exclude-chrom HG00100=chrY \
    --out-dir /path/gold_classifier
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from retro_miner.gold_classifier import (  # noqa: E402
    CATEGORICAL_FEATURES,
    NUMERIC_FEATURES,
    SampleGold,
    build_labeled_gold,
    dump_json,
    encode_feature_matrix,
    fit_gold_classifier,
    leave_one_sample_out_metrics,
    role_summary,
    training_table_from_labels,
)


def _parse_sample(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError(f"expected NAME=DIR, got {text}")
    name, raw = text.split("=", 1)
    return name, Path(raw)


def _parse_exclude(text: str) -> tuple[str, str]:
    if "=" not in text:
        raise argparse.ArgumentTypeError(f"expected SAMPLE=CHROM, got {text}")
    sample, chrom = text.split("=", 1)
    return sample, chrom


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", required=True, type=_parse_sample, dest="samples")
    parser.add_argument(
        "--exclude-chrom",
        action="append",
        default=[],
        type=_parse_exclude,
        help="Drop a chromosome for one sample, e.g. HG00100=chrY",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument(
        "--holdout-sample",
        action="append",
        default=[],
        help="Score these samples and leave them out of the saved fit",
    )
    args = parser.parse_args()
    excluded: dict[str, set[str]] = {}
    for sample, chrom in args.exclude_chrom:
        excluded.setdefault(sample, set()).add(chrom)
    specs = [
        SampleGold(sample=name, directory=path, exclude_chroms=frozenset(excluded.get(name, ())))
        for name, path in args.samples
    ]
    labeled = build_labeled_gold(specs)
    if labeled.empty:
        raise SystemExit("no gold rows found")
    training = training_table_from_labels(labeled, seed=args.seed)
    holdout = set(args.holdout_sample)
    fit_rows = training.loc[~training["sample"].isin(holdout)].copy()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = role_summary(labeled)
    summary.to_csv(args.out_dir / "role_summary.tsv", sep="\t", index=False)
    meta_cols = ["sample", "chrom", "gold_rank", "mei_family", "train_role", "train_label", "source_table"]
    keep_meta = [c for c in meta_cols if c in training.columns]
    feature_cols = [c for c in list(NUMERIC_FEATURES) + list(CATEGORICAL_FEATURES) if c in training.columns]
    training[keep_meta + feature_cols].to_csv(args.out_dir / "training_calls.tsv", sep="\t", index=False)
    outlier_cols = [c for c in meta_cols if c in labeled.columns]
    outliers = labeled.loc[labeled["train_role"].eq("excluded_outlier_known"), outlier_cols]
    outliers.to_csv(args.out_dir / "excluded_outlier_known.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(
        f"training_rows={len(training)} fit_rows={len(fit_rows)} "
        f"positives={int((fit_rows['train_label'] == 1).sum())} "
        f"negatives={int((fit_rows['train_label'] == 0).sum())}"
    )
    logo = leave_one_sample_out_metrics(fit_rows, seed=args.seed)
    dump_json(args.out_dir / "leave_one_sample_out.json", logo)
    for row in logo:
        print(
            f"holdout {row['held_out_sample']}: n={row['n']} "
            f"roc_auc={row['roc_auc']} average_precision={row['average_precision']}"
        )
    matrix, levels = encode_feature_matrix(fit_rows)
    model = fit_gold_classifier(matrix, fit_rows["train_label"], seed=args.seed)
    import joblib

    joblib.dump(
        {
            "model": model,
            "categorical_levels": levels,
            "numeric_features": list(NUMERIC_FEATURES),
            "categorical_features": list(CATEGORICAL_FEATURES),
            "seed": args.seed,
            "fit_samples": sorted(fit_rows["sample"].astype(str).unique()),
            "holdout_samples": sorted(holdout),
        },
        args.out_dir / "gold_classifier.joblib",
    )
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
