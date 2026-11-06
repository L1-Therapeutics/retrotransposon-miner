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
# in us-east-1; one CRT copy fills the disk, then fasterq-dump feeds bwa-mem3.
# The mem3 index is already built. Do not run bwa-mem3 index on a 64 GB host.
# This run is the test batch only: Hazen clone B2, its matched thymus, and
# two Sugo sister clones from mouse A. The rest of the catalog waits.
MEM3="${OUT}/reference/mm39-mem3"
mkdir -p "${MEM3}"
aws s3 sync "${RTM_S3_BUCKET}/mouse/public/reference/mm39-mem3/" "${MEM3}/" \
  --exclude "INDEX_COMPLETE"
for suffix in .amb .ann .bwt.2bit.64 .pac; do
  test -s "${MEM3}/mm39${suffix}"
done
if ! run bash -c 'command -v bwa-mem3 >/dev/null'; then
  micromamba install -y -n rtm-miner -c bioconda -c conda-forge 'bwa-mem3=0.16.0'
fi
run python scripts/align_mouse_sra.py \
  --config resources/public_datasets.yaml \
  --reference "${MEM3}/mm39" \
  --workdir "${HOME}/retrotransposon-workdir/mouse-align" \
  --s3-prefix "${RTM_S3_BUCKET}/mouse/align/mm39" \
  --threads 16 \
  --accession SRR3184315 \
  --accession SRR3184314 \
  --accession DRR697908 \
  --accession DRR697909
