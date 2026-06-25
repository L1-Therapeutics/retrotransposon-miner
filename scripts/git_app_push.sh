#!/usr/bin/env bash
set -euo pipefail

# Mint a short-lived GitHub App installation token and push.
#
# Usage:
#   scripts/git_app_push.sh
#   scripts/git_app_push.sh origin main
#   scripts/git_app_push.sh origin HEAD:main
#
# Required env (or defaults below):
#   GITHUB_APP_ID
#   GITHUB_APP_INSTALLATION_ID
#   GITHUB_APP_KEY_PATH

GITHUB_APP_ID="${GITHUB_APP_ID:-4147126}"
GITHUB_APP_INSTALLATION_ID="${GITHUB_APP_INSTALLATION_ID:-142633660}"
GITHUB_APP_KEY_PATH="${GITHUB_APP_KEY_PATH:-$HOME/.config/gh-app/retrotransposon-miner-vm-pusher.pem}"
GITHUB_OWNER="${GITHUB_OWNER:-L1-Therapeutics}"
GITHUB_REPO="${GITHUB_REPO:-retrotransposon-miner}"

if [[ ! -f "${GITHUB_APP_KEY_PATH}" ]]; then
  echo "ERROR: missing GitHub App private key: ${GITHUB_APP_KEY_PATH}" >&2
  exit 1
fi

b64url() {
  openssl base64 -A | tr '+/' '-_' | tr -d '='
}

now_epoch="$(date +%s)"
iat="$((now_epoch - 60))"
exp="$((now_epoch + 540))"

header="$(printf '{"alg":"RS256","typ":"JWT"}' | b64url)"
payload="$(printf '{"iat":%s,"exp":%s,"iss":"%s"}' "${iat}" "${exp}" "${GITHUB_APP_ID}" | b64url)"
unsigned="${header}.${payload}"
sig="$(printf %s "${unsigned}" | openssl dgst -binary -sha256 -sign "${GITHUB_APP_KEY_PATH}" | b64url)"
jwt="${unsigned}.${sig}"

token="$(
  curl -fsSL -X POST \
    -H "Authorization: Bearer ${jwt}" \
    -H "Accept: application/vnd.github+json" \
    "https://api.github.com/app/installations/${GITHUB_APP_INSTALLATION_ID}/access_tokens" \
    | /home/ec2-user/.local/share/mamba/envs/rtm-miner/bin/python -c 'import json,sys; print(json.load(sys.stdin)["token"])'
)"

remote="${1:-origin}"
shift || true

if [[ "${remote}" != "origin" ]]; then
  echo "ERROR: this helper currently supports pushing to 'origin' only" >&2
  exit 1
fi

push_url="https://x-access-token:${token}@github.com/${GITHUB_OWNER}/${GITHUB_REPO}.git"
git push "${push_url}" "$@"
