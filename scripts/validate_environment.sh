#!/usr/bin/env bash
set -euo pipefail

echo "Validating retrotransposon-miner environment..."

required_bins=(
  git
  jupyter
  samtools
  bedtools
  minimap2
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
samtools --version | awk 'NR==1 {print "samtools:", $0}'
bedtools --version | awk '{print "bedtools:", $0}'
minimap2 --version | awk '{print "minimap2:", $0}'
bwa-mem2 version 2>&1 | awk 'NR==1 {print "bwa-mem2:", $0}'
bcftools --version | awk 'NR==1 {print "bcftools:", $0}'
if command -v liftOver >/dev/null 2>&1; then
  liftOver 2>&1 | awk 'NR==1 {print "liftOver:", $0}'
fi

echo "Environment validation complete."
