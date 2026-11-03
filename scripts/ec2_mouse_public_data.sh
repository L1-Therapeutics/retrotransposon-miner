#!/usr/bin/env bash
# Dedicated mouse public-data instance.
#
# Uses a separate name, SSH alias, state file, and S3 prefix from the human
# box launched by scripts/ec2_jlab.sh. Default action brings the instance up,
# copies this working tree, installs the miner env, then starts two downloads:
# mm39 reference/repeats/polymorphisms first, then Hazen 2016 and Sugo 2025
# FASTQs streamed to S3.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export INSTANCE_NAME="${INSTANCE_NAME:-retrotransposon-miner-mouse}"
export HOST_ALIAS="${HOST_ALIAS:-retro-ec2-mouse}"
export JLAB_ALIAS="${JLAB_ALIAS:-jlab-mouse}"
export JLAB_LOCAL_PORT="${JLAB_LOCAL_PORT:-8891}"
export INSTANCE_STATE_FILE="${INSTANCE_STATE_FILE:-${REPO_ROOT}/.ec2-mouse.env}"
export ROOT_VOLUME_GB="${ROOT_VOLUME_GB:-500}"
export S3_BUCKET="${S3_BUCKET:-s3://l1tx-data}"
export S3_CACHE_PREFIX="${S3_CACHE_PREFIX:-s3://l1tx-data/mouse/public}"
export INSTANCE_TYPE="${INSTANCE_TYPE:-m7i.4xlarge}"
export SPOT="${SPOT:-0}"
export SSH_USER="${SSH_USER:-ec2-user}"

REMOTE_REPO="${REMOTE_REPO:-retrotransposon-miner}"
LOG_REL="logs/mouse-public-data.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }

ec2() {
  bash "${REPO_ROOT}/scripts/ec2_jlab.sh" "$@"
}

sync_tree() {
  log "Copying working tree to ${HOST_ALIAS}:~/${REMOTE_REPO} (no .git, data, or results)"
  tar -C "${REPO_ROOT}" \
    --exclude .git \
    --exclude data \
    --exclude results \
    --exclude workdir-results \
    --exclude .ec2-instance.env \
    --exclude .ec2-mouse.env \
    --exclude '__pycache__' \
    --exclude '.pytest_cache' \
    --exclude '*.pem' \
    -czf - . \
    | ssh "${HOST_ALIAS}" "mkdir -p ~/${REMOTE_REPO} && tar -xzf - -C ~/${REMOTE_REPO}"
}

install_remote_env() {
  log "Installing rtm-miner and UCSC tools on ${HOST_ALIAS}"
  ssh "${HOST_ALIAS}" "bash -lc 'set -euo pipefail; export PATH=\"\$HOME/.local/bin:\$PATH\"; cd ~/${REMOTE_REPO}; bash scripts/bootstrap_env.sh; bash scripts/install_ucsc_tools.sh'"
}

start_download() {
  log "Starting mm39 public-data download, then Hazen/Sugo FASTQ streams"
  ssh "${HOST_ALIAS}" "bash -lc 'set -euo pipefail
export PATH=\"\$HOME/.local/bin:\$PATH\"
mkdir -p \"\$HOME/logs\" \"\$HOME/retrotransposon-workdir/data/public\"
if [[ -f \"\$HOME/logs/mouse-public-data.pid\" ]] && kill -0 \"\$(cat \"\$HOME/logs/mouse-public-data.pid\")\" 2>/dev/null; then
  echo \"Download already running pid \$(cat \"\$HOME/logs/mouse-public-data.pid\")\"
  exit 0
fi
chmod 755 \"\$HOME/${REMOTE_REPO}/scripts/run_mouse_public_data_download.sh\"
nohup \"\$HOME/${REMOTE_REPO}/scripts/run_mouse_public_data_download.sh\" > \"\$HOME/${LOG_REL}\" 2>&1 &
echo \$! > \"\$HOME/logs/mouse-public-data.pid\"
echo \"pid \$(cat \"\$HOME/logs/mouse-public-data.pid\") log ~/${LOG_REL}\"
'"
}

show_log() {
  ssh "${HOST_ALIAS}" "tail -n ${1:-80} ~/${LOG_REL}"
}

usage() {
  cat <<EOF
Usage: $0 [bootstrap|sync|install|download|log|status|help]

  bootstrap   Launch or reuse the mouse instance, sync this tree, install the
              env, and start the mm39 download plus FASTQ streams (default).
  sync        Copy the working tree to the instance.
  install     bootstrap_env.sh and install_ucsc_tools.sh on the instance.
  download    Start the remote download if it is not already running.
  log         Tail the remote download log.
  status      EC2 status for this mouse instance.

Other arguments are passed to scripts/ec2_jlab.sh with the mouse instance
env (stop-instance, up, ssh-config, ...).

Env overrides: INSTANCE_TYPE, ROOT_VOLUME_GB, SPOT, S3_BUCKET, S3_CACHE_PREFIX,
REGION, KEY_PATH. State file: ${INSTANCE_STATE_FILE}
SSH: ssh ${HOST_ALIAS}
S3:  ${S3_CACHE_PREFIX}
EOF
}

cmd="${1:-bootstrap}"
shift || true

case "${cmd}" in
  help|-h|--help) usage ;;
  bootstrap)
    ec2 bootstrap
    sync_tree
    install_remote_env
    start_download
    log "Mouse public-data download is running on ${HOST_ALIAS}."
    log "S3 prefix: ${S3_CACHE_PREFIX}"
    log "Follow it with: $0 log"
    ;;
  sync) sync_tree ;;
  install) install_remote_env ;;
  download) start_download ;;
  log) show_log "${1:-80}" ;;
  status)
    ec2 status
    ssh "${HOST_ALIAS}" "if [[ -f \$HOME/logs/mouse-public-data.pid ]]; then echo -n 'download pid '; cat \$HOME/logs/mouse-public-data.pid; echo; tail -n 20 ~/${LOG_REL} || true; else echo 'download not started'; fi" || true
    ;;
  *)
    ec2 "${cmd}" "$@"
    ;;
esac
