#!/usr/bin/env bash
set -euo pipefail

# Configure external data/output directories outside the repo checkout.
# Usage (recommended):
#   source scripts/use_external_workdir.sh
#   source scripts/use_external_workdir.sh /custom/workdir/path

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "This script should be sourced so exports persist in your shell." >&2
  echo "Run: source scripts/use_external_workdir.sh [optional-workdir]" >&2
  exit 1
fi

RTM_WORKDIR_INPUT="${1:-${RTM_WORKDIR:-${HOME}/retrotransposon-workdir}}"

export RTM_WORKDIR="${RTM_WORKDIR_INPUT}"
export RTM_PUBLIC_DATA_DIR="${RTM_PUBLIC_DATA_DIR:-${RTM_WORKDIR}/data/public}"
export RTM_RESULTS_DIR="${RTM_RESULTS_DIR:-${RTM_WORKDIR}/results}"

mkdir -p "${RTM_PUBLIC_DATA_DIR}" "${RTM_RESULTS_DIR}" "${RTM_WORKDIR}/logs"

echo "Configured external runtime directories:"
echo "  RTM_WORKDIR=${RTM_WORKDIR}"
echo "  RTM_PUBLIC_DATA_DIR=${RTM_PUBLIC_DATA_DIR}"
echo "  RTM_RESULTS_DIR=${RTM_RESULTS_DIR}"
