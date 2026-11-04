#!/usr/bin/env bash
# Pull the shared bwa-mem3 mm39 index from S3 and align one accession.
# Does not build the index. A 64 GB host cannot.
set -euo pipefail

export PATH="${HOME}/.local/bin:${PATH}"
BOOT="s3://l1tx-data/mouse/align/bootstrap"
REF_SRC="s3://l1tx-data/mouse/public/reference/mm39-mem3"
REF_DIR="${HOME}/retrotransposon-workdir/data/public/reference/mm39-mem3"
REPO="${HOME}/retrotransposon-miner"
ACC="$(tr -d '[:space:]' < "${HOME}/accession.txt")"

mkdir -p "${HOME}/.local" "${REPO}/scripts" "${REPO}/src/retro_miner" "${REPO}/resources" "${HOME}/logs" "${REF_DIR}"

for _ in 1 2 3 4 5 6 7 8 9 10; do
  if aws s3 ls "${BOOT}/rtm-miner-env.tar" >/dev/null 2>&1; then
    break
  fi
  sleep 15
done
aws s3 cp "${BOOT}/rtm-miner-env.tar" - | tar -C "${HOME}/.local" -xf -
aws s3 cp "${BOOT}/align_mouse_sra.py" "${REPO}/scripts/align_mouse_sra.py"
aws s3 cp "${BOOT}/mouse_public_data.py" "${REPO}/src/retro_miner/mouse_public_data.py"
aws s3 cp "${BOOT}/retro_miner_init.py" "${REPO}/src/retro_miner/__init__.py"
aws s3 cp "${BOOT}/public_datasets.yaml" "${REPO}/resources/public_datasets.yaml"
chmod +x "${HOME}/.local/bin/micromamba"

# Existing bwa-mem3 index. Sync does not rebuild it.
aws s3 sync "${REF_SRC}/" "${REF_DIR}/" --exclude "INDEX_COMPLETE"
for suffix in .amb .ann .bwt.2bit.64 .pac; do
  test -s "${REF_DIR}/mm39${suffix}"
done
if ! micromamba run -n rtm-miner bash -c 'command -v bwa-mem3 >/dev/null'; then
  micromamba install -y -n rtm-miner -c bioconda -c conda-forge 'bwa-mem3=0.16.0'
fi

cd "${REPO}"
setsid nohup micromamba run -n rtm-miner python scripts/align_mouse_sra.py \
  --config resources/public_datasets.yaml \
  --reference "${REF_DIR}/mm39" \
  --workdir "${HOME}/retrotransposon-workdir/mouse-align" \
  --s3-prefix s3://l1tx-data/mouse/align/mm39 \
  --threads 16 \
  --accession "${ACC}" \
  > "${HOME}/logs/mouse-align.log" 2>&1 < /dev/null &
echo $! > "${HOME}/logs/mouse-align.pid"
echo "started ${ACC} pid $(cat "${HOME}/logs/mouse-align.pid")"
setsid nohup bash "${HOME}/retrotransposon-miner/scripts/shutdown_when_aligned.sh" "${ACC}" \
  > "${HOME}/logs/mouse-align-shutdown.nohup" 2>&1 < /dev/null &
