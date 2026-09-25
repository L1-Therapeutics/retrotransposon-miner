"""Copy a run's gold table and plot directories to the public-data bucket.

The destination is ``s3://<bucket>/results/<outdir name>/``. The bucket is the
one that already holds public reference data (``RTM_S3_BUCKET``), not the
``public/`` prefix inside it. Per-chromosome candidate tables stay local.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from retro_miner.s3_transfer import _run_aws_cli, process_local_aws_config

GENOME_GOLD_NAME = "candidate_loci.mei.gold_review.tsv"
GENOME_VCF_NAME = "candidate_loci.mei.gold_review.vcf"
IGV_DIR_NAME = "candidate_loci.mei.gold_review.igv"
READ_ARCH_DIR_NAME = "candidate_loci.mei.read_architecture"


def bucket_uri(raw: str) -> str:
    """Return ``s3://bucket`` with no key."""
    text = (raw or "").strip()
    if not text:
        raise ValueError("S3 bucket is empty")
    if not text.startswith("s3://"):
        text = "s3://" + text
    name = text[len("s3://") :].split("/", 1)[0].strip()
    if not name:
        raise ValueError(f"S3 bucket is empty: {raw}")
    return f"s3://{name}"


def resolve_bucket(explicit: str | None = None) -> str:
    """Bucket from the flag, else ``RTM_S3_BUCKET``, else ``S3_BUCKET``."""
    raw = (explicit or "").strip()
    if not raw:
        raw = (os.environ.get("RTM_S3_BUCKET") or os.environ.get("S3_BUCKET") or "").strip()
    if not raw:
        raise ValueError("S3 bucket is empty; set RTM_S3_BUCKET or pass --bucket")
    return bucket_uri(raw)


def results_destination(bucket: str, sample_name: str) -> str:
    """``s3://bucket/results/<sample>``."""
    name = sample_name.strip().strip("/")
    if not name or "/" in name:
        raise ValueError(f"sample name must be a single path segment: {sample_name!r}")
    return f"{bucket_uri(bucket)}/results/{name}"


def sync_argv(outdir: Path, dest_prefix: str) -> list[str]:
    """Sync the genome gold table plus IGV and read-architecture plots.

    The include for the gold table is the file at the run root. A per-chromosome
    ``chr1/candidate_loci.mei.gold_review.tsv`` does not match that name, so the
    full silver/bronze tables are not uploaded.
    """
    return [
        "s3",
        "sync",
        str(outdir),
        dest_prefix.rstrip("/"),
        "--exclude",
        "*",
        "--include",
        GENOME_GOLD_NAME,
        "--include",
        GENOME_VCF_NAME,
        "--include",
        f"*/{IGV_DIR_NAME}/*",
        "--include",
        f"*/{READ_ARCH_DIR_NAME}/*",
        "--include",
        f"{IGV_DIR_NAME}/*",
        "--include",
        f"{READ_ARCH_DIR_NAME}/*",
        "--only-show-errors",
    ]


def publish_outdir(
    outdir: Path,
    bucket: str,
    *,
    run_aws=None,
) -> str:
    """Upload gold calls and plots. Returns the destination prefix."""
    dest = results_destination(bucket, outdir.name)
    argv = sync_argv(outdir, dest)
    if run_aws is None:
        with process_local_aws_config() as env:
            _run_aws_cli(argv, env=env)
    else:
        run_aws(argv)
    return dest


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Upload the genome gold table, IGV plots, and read-architecture plots."
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--bucket",
        default="",
        help="s3://bucket that holds public data. Default: RTM_S3_BUCKET or S3_BUCKET.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the destination and the aws sync arguments without uploading.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    outdir = args.outdir
    if not outdir.is_dir():
        raise SystemExit(f"ERROR: outdir is not a directory: {outdir}")
    bucket = resolve_bucket(args.bucket)
    dest = results_destination(bucket, outdir.name)
    if args.dry_run:
        print(dest)
        print(" ".join(sync_argv(outdir, dest)))
        return 0
    published = publish_outdir(outdir, bucket)
    print(f"[publish-results] uploaded gold calls and plots to {published}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
