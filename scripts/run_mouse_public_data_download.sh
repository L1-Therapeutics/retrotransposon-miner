#!/usr/bin/env bash
# Run on the mouse EC2 instance after this repo is copied there.
# Reference tracks are prepared before the multi-terabyte FASTQ streams start.
set -euo pipefail

export PATH="${HOME}/.local/bin:${PATH}"
set -a
# shellcheck disable=SC1090
. "${HOME}/.config/rtm/s3.env"
set +a

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"
# The micromamba shell hook expands MAMBA_ROOT_PREFIX before it is set, which
# aborts this script under `set -u`. `micromamba run` uses the env directly.
run() {
  micromamba run -n rtm-miner "$@"
}

OUT="${RTM_PUBLIC_DATA_DIR:-${HOME}/retrotransposon-workdir/data/public}"
mkdir -p "${OUT}"

run python scripts/download_public_data.py \
  --references mm39 \
  --s3-cache-prefix "${RTM_S3_CACHE}" \
  --outdir "${OUT}" \
  --download-workers 4

# ENA HTTPS drops the FASTQ streams. The same reads are the public SRA object
# in us-east-1; one CRT copy fills the disk, then fasterq-dump feeds bwa mem.
# This run is the test batch only: Hazen clone B2, its matched thymus, and
# two Sugo sister clones from mouse A. The rest of the catalog waits.
run python scripts/align_mouse_sra.py \
  --config resources/public_datasets.yaml \
  --reference "${OUT}/reference/mm39/mm39.fa" \
  --workdir "${HOME}/retrotransposon-workdir/mouse-align" \
  --s3-prefix "${RTM_S3_BUCKET}/mouse/align/mm39" \
  --threads 16 \
  --accession SRR3184315 \
  --accession SRR3184314 \
  --accession DRR697908 \
  --accession DRR697909
