#!/usr/bin/env bash
set -euo pipefail

# Bootstrap environment on Linux even if conda/mamba are missing.
# - Uses mamba/conda if available.
# - Falls back to installing micromamba into ~/.local/bin.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/environment.yml"
ENV_NAME="${1:-rtm-miner}"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "ERROR: missing ${ENV_FILE}" >&2
  exit 1
fi

echo "Bootstrapping environment '${ENV_NAME}' from ${ENV_FILE}"

ensure_system_git() {
  if command -v git >/dev/null 2>&1; then
    return 0
  fi

  if [[ "$(uname -s)" != "Linux" ]]; then
    echo "WARN: git is missing and automatic install is only configured for Linux." >&2
    return 0
  fi

  echo "git not found; attempting system install..."
  if command -v dnf >/dev/null 2>&1; then
    sudo dnf install -y git
  elif command -v yum >/dev/null 2>&1; then
    sudo yum install -y git
  elif command -v apt-get >/dev/null 2>&1; then
    sudo apt-get update
    sudo apt-get install -y git
  else
    echo "WARN: no supported package manager found for git install." >&2
    return 0
  fi

  if command -v git >/dev/null 2>&1; then
    echo "git installed successfully"
  else
    echo "WARN: git installation attempted but git is still missing." >&2
  fi
}

ensure_system_git

create_or_update_with_mamba() {
  if mamba env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    mamba env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
  else
    mamba env create -n "${ENV_NAME}" -f "${ENV_FILE}"
  fi
}

create_or_update_with_conda() {
  if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
  else
    conda env create -n "${ENV_NAME}" -f "${ENV_FILE}"
  fi
}

create_or_update_with_micromamba() {
  local mm_cmd="$1"
  if "${mm_cmd}" env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    "${mm_cmd}" env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
  else
    "${mm_cmd}" create -y -n "${ENV_NAME}" -f "${ENV_FILE}"
  fi
}

if command -v mamba >/dev/null 2>&1; then
  echo "Found mamba"
  create_or_update_with_mamba
elif command -v conda >/dev/null 2>&1; then
  echo "Found conda"
  create_or_update_with_conda
else
  echo "mamba/conda not found. Installing micromamba..."
  mkdir -p "${HOME}/.local/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "${HOME}/.local/bin" --strip-components=1 bin/micromamba

  if command -v micromamba >/dev/null 2>&1; then
    MM_CMD="micromamba"
  elif [[ -x "${HOME}/.local/bin/micromamba" ]]; then
    MM_CMD="${HOME}/.local/bin/micromamba"
  else
    echo "ERROR: micromamba installation failed." >&2
    exit 1
  fi

  create_or_update_with_micromamba "${MM_CMD}"
fi

if [[ "$(uname -s)" == "Linux" ]]; then
  echo
  echo "Checking headless IGV dependencies (system Xvfb)..."
  bash "${ROOT_DIR}/scripts/install_headless_igv_deps.sh" || \
    echo "WARN: could not install headless IGV deps; run: bash scripts/install_headless_igv_deps.sh" >&2
fi

echo
echo "Environment ready."
echo "Activate it with one of:"
echo "  micromamba activate ${ENV_NAME}"
echo "  mamba activate ${ENV_NAME}"
echo "  conda activate ${ENV_NAME}"
echo
echo "Then run:"
echo "  bash scripts/validate_environment.sh"
