#!/usr/bin/env bash
# Run on the mouse EC2 instance after scripts/ec2_mouse_public_data.sh sync.
# Reference tracks are prepared before the multi-terabyte FASTQ streams start.
set -euo pipefail

export PATH="${HOME}/.local/bin:${PATH}"
set -a
# shellcheck disable=SC1090
. "${HOME}/.config/rtm/s3.env"
set +a

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rtm-miner

OUT="${RTM_PUBLIC_DATA_DIR:-${HOME}/retrotransposon-workdir/data/public}"
mkdir -p "${OUT}"

python scripts/download_public_data.py \
  --references mm39 \
  --s3-cache-prefix "${RTM_S3_CACHE}" \
  --outdir "${OUT}" \
  --download-workers 4

python scripts/download_public_data.py \
  --references mm39 \
  --include-mouse-test-fastqs \
  --skip-postprocess \
  --s3-cache-prefix "${RTM_S3_CACHE}" \
  --outdir "${OUT}" \
  --download-workers 2
