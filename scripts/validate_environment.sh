#!/usr/bin/env bash
set -euo pipefail

echo "Validating retrotransposon-miner environment..."

required_bins=(
  git
  jupyter
  samtools
  bedtools
  minimap2
  bwa
  bwa-mem2
  bcftools
  liftOver
  bigBedToBed
  bigWigToBedGraph
)

for bin in "${required_bins[@]}"; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    echo "ERROR: missing binary: ${bin}" >&2
    exit 1
  fi
done

optional_bins=(
  igv
  spades.py
  java
  Xvfb
  xvfb-run
)

for bin in "${optional_bins[@]}"; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    echo "WARN: optional binary not found: ${bin}" >&2
  fi
done

if [[ "$(uname -s)" == "Linux" ]] && [[ -z "${DISPLAY:-}" ]]; then
  if ! command -v xvfb-run >/dev/null 2>&1 && ! command -v Xvfb >/dev/null 2>&1; then
    echo "WARN: headless Linux without Xvfb; IGV snapshots will fail." >&2
    echo "      Run: bash scripts/install_headless_igv_deps.sh" >&2
  fi
fi

report_version() {
  local label="$1"
  shift
  local out
  out="$("$@" 2>&1)" || true
  if [[ -n "${out}" ]]; then
    printf '%s: %s\n' "${label}" "${out}" | awk 'NR==1'
  fi
}

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  echo "ERROR: ${PYTHON_BIN} not found in PATH." >&2
  exit 1
fi

"${PYTHON_BIN}" - <<'PY'
from importlib.util import find_spec
import sys

mods = [
    "matplotlib",
    "plotly",
    "ipykernel",
    "pysam",
    "pandas",
    "pyarrow",
    "h5py",
    "click",
    "Bio",
]

missing = [m for m in mods if find_spec(m) is None]
if missing:
    print(f"ERROR: missing python modules: {', '.join(missing)}", file=sys.stderr)
    raise SystemExit(1)

print("Python modules OK")
PY

echo "All required tools detected."
report_version "samtools" samtools --version
report_version "bedtools" bedtools --version
report_version "minimap2" minimap2 --version
report_version "bwa" bwa
report_version "bwa-mem2" bwa-mem2 version
report_version "bcftools" bcftools --version
if command -v liftOver >/dev/null 2>&1; then
  report_version "liftOver" liftOver
fi

echo "Environment validation complete."
