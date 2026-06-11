#!/usr/bin/env bash
set -euo pipefail

# Installs selected UCSC command-line binaries.
# Preferred path:
#   1) Install from bioconda into an existing env (best shared-lib compatibility)
# Fallback path:
#   2) Download UCSC standalone Linux binaries into ~/.local/ucsc

INSTALL_DIR="${1:-$HOME/.local/ucsc}"
BASE_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64"
ENV_NAME="${ENV_NAME:-rtm-miner}"

install_with_pkg_manager() {
  local installer="$1"
  echo "Attempting UCSC install via ${installer} into env ${ENV_NAME}..."
  "${installer}" install -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
    ucsc-liftover ucsc-twobittofa ucsc-fatotwobit ucsc-bigbedtobed ucsc-bigwigtobedgraph
}

if command -v micromamba >/dev/null 2>&1; then
  if install_with_pkg_manager micromamba; then
    echo "Installed UCSC tools via micromamba."
    echo "Run with:"
    echo "  micromamba activate ${ENV_NAME}"
    echo "  liftOver"
    exit 0
  fi
elif command -v mamba >/dev/null 2>&1; then
  if install_with_pkg_manager mamba; then
    echo "Installed UCSC tools via mamba."
    exit 0
  fi
elif command -v conda >/dev/null 2>&1; then
  if install_with_pkg_manager conda; then
    echo "Installed UCSC tools via conda."
    exit 0
  fi
fi

echo "Package-manager UCSC install not available/failed; falling back to standalone binaries."
echo "Note: some Linux images emit a libcurl warning for standalone liftOver."

mkdir -p "${INSTALL_DIR}"

for tool in liftOver twoBitToFa faToTwoBit bigBedToBed bigWigToBedGraph; do
  echo "Installing ${tool} -> ${INSTALL_DIR}"
  curl -fsSL "${BASE_URL}/${tool}" -o "${INSTALL_DIR}/${tool}"
  chmod +x "${INSTALL_DIR}/${tool}"
done

echo "Done."
echo "Add to PATH:"
echo "  export PATH=\"${INSTALL_DIR}:\$PATH\""
