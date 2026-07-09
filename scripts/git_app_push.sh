#!/usr/bin/env bash
set -euo pipefail

# Mint a short-lived GitHub App installation token, push the current HEAD, and
# open/reuse a pull request into the default base branch.
#
# Usage:
#   scripts/git_app_push.sh
#   scripts/git_app_push.sh --branch feat/my-change
#   scripts/git_app_push.sh --branch feat/my-change --title "My change" --body "Details"
#   scripts/git_app_push.sh --no-pr
#   scripts/git_app_push.sh origin HEAD:refs/heads/feat/manual   # legacy raw git-push mode
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
GITHUB_BASE_BRANCH="${GITHUB_BASE_BRANCH:-main}"
PYTHON_BIN="${PYTHON_BIN:-/home/ec2-user/.local/share/mamba/envs/rtm-miner/bin/python}"
if [[ ! -x "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3)"
fi

if [[ ! -f "${GITHUB_APP_KEY_PATH}" ]]; then
  echo "ERROR: missing GitHub App private key: ${GITHUB_APP_KEY_PATH}" >&2
  exit 1
fi

usage() {
  cat <<'EOF'
Usage:
  scripts/git_app_push.sh [--branch NAME] [--title TEXT] [--body TEXT] [--base BRANCH] [--no-pr]
  scripts/git_app_push.sh origin <git-push-refspec...>

Default mode pushes HEAD to a writable feature branch and opens/reuses a PR into
the base branch (default: main). Protected branches are never pushed directly.
EOF
}

b64url() {
  openssl base64 -A | tr '+/' '-_' | tr -d '='
}

mint_token() {
  local now_epoch iat exp header payload unsigned sig jwt
  now_epoch="$(date +%s)"
  iat="$((now_epoch - 60))"
  exp="$((now_epoch + 540))"
  header="$(printf '{"alg":"RS256","typ":"JWT"}' | b64url)"
  payload="$(printf '{"iat":%s,"exp":%s,"iss":"%s"}' "${iat}" "${exp}" "${GITHUB_APP_ID}" | b64url)"
  unsigned="${header}.${payload}"
  sig="$(printf %s "${unsigned}" | openssl dgst -binary -sha256 -sign "${GITHUB_APP_KEY_PATH}" | b64url)"
  jwt="${unsigned}.${sig}"
  curl -fsSL -X POST \
    -H "Authorization: Bearer ${jwt}" \
    -H "Accept: application/vnd.github+json" \
    "https://api.github.com/app/installations/${GITHUB_APP_INSTALLATION_ID}/access_tokens" \
    | "${PYTHON_BIN}" -c 'import json,sys; print(json.load(sys.stdin)["token"])'
}

slugify() {
  printf '%s' "$1" \
    | tr '[:upper:]' '[:lower:]' \
    | sed -E 's/[^a-z0-9]+/-/g; s/^-+//; s/-+$//; s/-+/-/g' \
    | cut -c1-48
}

default_branch_name() {
  local local_branch subject slug stamp
  local_branch="$(git rev-parse --abbrev-ref HEAD)"
  if [[ "${local_branch}" != "HEAD" && "${local_branch}" != "main" && "${local_branch}" != "master" ]]; then
    # Prefer a feat/ prefix when the local branch looks protected/legacy.
    if [[ "${local_branch}" == feat/* || "${local_branch}" == fix/* || "${local_branch}" == chore/* ]]; then
      # Avoid pushing to known-protected fix/* branches by rewriting to feat/.
      if [[ "${local_branch}" == fix/* ]]; then
        printf 'feat/%s\n' "${local_branch#fix/}"
      else
        printf '%s\n' "${local_branch}"
      fi
      return 0
    fi
  fi
  subject="$(git log -1 --pretty=%s | slugify)"
  if [[ -z "${subject}" ]]; then
    subject="update"
  fi
  stamp="$(date +%Y%m%d-%H%M%S)"
  printf 'feat/%s-%s\n' "${subject}" "${stamp}"
}

default_pr_title() {
  git log -1 --pretty=%s
}

default_pr_body() {
  local base="$1"
  local range
  if git rev-parse --verify "origin/${base}" >/dev/null 2>&1; then
    range="origin/${base}..HEAD"
  else
    range="HEAD"
  fi
  {
    echo "## Summary"
    git log --reverse --pretty='- %s' "${range}"
    echo
    echo "## Test plan"
    echo "- [ ] CI passes"
    echo "- [ ] Spot-check changed behavior / docs"
  }
}

api_json() {
  local token="$1" method="$2" url="$3" body="${4-}"
  if [[ -n "${body}" ]]; then
    curl -fsSL -X "${method}" \
      -H "Authorization: Bearer ${token}" \
      -H "Accept: application/vnd.github+json" \
      -H "Content-Type: application/json" \
      --data "${body}" \
      "${url}"
  else
    curl -fsSL -X "${method}" \
      -H "Authorization: Bearer ${token}" \
      -H "Accept: application/vnd.github+json" \
      "${url}"
  fi
}

ensure_pr() {
  local token="$1" head_branch="$2" base_branch="$3" title="$4" body="$5"
  local existing pr_json
  existing="$(
    api_json "${token}" GET \
      "https://api.github.com/repos/${GITHUB_OWNER}/${GITHUB_REPO}/pulls?state=open&head=${GITHUB_OWNER}:${head_branch}&base=${base_branch}"
  )"
  pr_json="$(
    TITLE="${title}" BODY="${body}" HEAD_BRANCH="${head_branch}" BASE_BRANCH="${base_branch}" \
    EXISTING="${existing}" OWNER="${GITHUB_OWNER}" REPO="${GITHUB_REPO}" TOKEN="${token}" \
    "${PYTHON_BIN}" - <<'PY'
import json, os, urllib.request

existing = json.loads(os.environ["EXISTING"])
if existing:
    print(existing[0]["html_url"])
    raise SystemExit(0)

payload = {
    "title": os.environ["TITLE"],
    "body": os.environ["BODY"],
    "head": os.environ["HEAD_BRANCH"],
    "base": os.environ["BASE_BRANCH"],
}
req = urllib.request.Request(
    f"https://api.github.com/repos/{os.environ['OWNER']}/{os.environ['REPO']}/pulls",
    data=json.dumps(payload).encode(),
    method="POST",
    headers={
        "Authorization": f"Bearer {os.environ['TOKEN']}",
        "Accept": "application/vnd.github+json",
        "Content-Type": "application/json",
        "User-Agent": "rtm-git-app-push",
    },
)
with urllib.request.urlopen(req) as resp:
    print(json.load(resp)["html_url"])
PY
  )"
  printf '%s\n' "${pr_json}"
}

BRANCH=""
TITLE=""
BODY=""
BASE_BRANCH="${GITHUB_BASE_BRANCH}"
CREATE_PR=1
LEGACY_MODE=0

if [[ "${1-}" == "origin" ]]; then
  LEGACY_MODE=1
fi

if [[ "${LEGACY_MODE}" -eq 0 ]]; then
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --branch)
        BRANCH="${2-}"
        shift 2
        ;;
      --title)
        TITLE="${2-}"
        shift 2
        ;;
      --body)
        BODY="${2-}"
        shift 2
        ;;
      --base)
        BASE_BRANCH="${2-}"
        shift 2
        ;;
      --no-pr)
        CREATE_PR=0
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        echo "ERROR: unknown argument: $1" >&2
        usage >&2
        exit 1
        ;;
    esac
  done
fi

token="$(mint_token)"
push_url="https://x-access-token:${token}@github.com/${GITHUB_OWNER}/${GITHUB_REPO}.git"

if [[ "${LEGACY_MODE}" -eq 1 ]]; then
  remote="${1}"
  shift
  if [[ "${remote}" != "origin" ]]; then
    echo "ERROR: this helper currently supports pushing to 'origin' only" >&2
    exit 1
  fi
  git push "${push_url}" "$@"
  exit 0
fi

if [[ -z "${BRANCH}" ]]; then
  BRANCH="$(default_branch_name)"
fi
if [[ "${BRANCH}" == "main" || "${BRANCH}" == "master" ]]; then
  echo "ERROR: refusing to push directly to protected branch '${BRANCH}'" >&2
  exit 1
fi
if [[ -z "${TITLE}" ]]; then
  TITLE="$(default_pr_title)"
fi
if [[ -z "${BODY}" ]]; then
  BODY="$(default_pr_body "${BASE_BRANCH}")"
fi

push_head_to_branch() {
  local target_branch="$1"
  # Do not use `git push -u` with a tokenized URL: that stores the token in
  # branch.<name>.remote. Push, then set upstream to the normal origin remote.
  git push "${push_url}" "HEAD:refs/heads/${target_branch}"
}

echo "Pushing HEAD -> origin/${BRANCH}"
if ! push_head_to_branch "${BRANCH}" 2>/tmp/rtm_git_app_push.err; then
  if grep -q "Changes must be made through a pull request\|protected branch\|GH013" /tmp/rtm_git_app_push.err; then
    fallback="feat/$(slugify "$(git log -1 --pretty=%s)")-$(date +%Y%m%d-%H%M%S)"
    echo "Branch origin/${BRANCH} is protected; falling back to origin/${fallback}" >&2
    cat /tmp/rtm_git_app_push.err >&2 || true
    BRANCH="${fallback}"
    push_head_to_branch "${BRANCH}"
  else
    cat /tmp/rtm_git_app_push.err >&2 || true
    exit 1
  fi
fi
rm -f /tmp/rtm_git_app_push.err

# Ensure the remote-tracking ref exists, then track the normal origin remote.
git fetch origin "refs/heads/${BRANCH}:refs/remotes/origin/${BRANCH}" >/dev/null 2>&1 || true
git branch --set-upstream-to="origin/${BRANCH}" >/dev/null 2>&1 || true

if [[ "${CREATE_PR}" -eq 1 ]]; then
  echo "Ensuring pull request ${BRANCH} -> ${BASE_BRANCH}"
  if ! pr_url="$(ensure_pr "${token}" "${BRANCH}" "${BASE_BRANCH}" "${TITLE}" "${BODY}")"; then
    echo "ERROR: push succeeded but PR create/reuse failed for ${BRANCH} -> ${BASE_BRANCH}" >&2
    echo "Open manually: https://github.com/${GITHUB_OWNER}/${GITHUB_REPO}/compare/${BASE_BRANCH}...${BRANCH}" >&2
    exit 1
  fi
  echo "PR: ${pr_url}"
fi
