#!/usr/bin/env bash
set -euo pipefail

# Smoke test:
# - Accepts disease/control BAM paths or HTTP/FTP URLs
# - Optional local download mode for remote BAMs (recommended)
# - Extracts one chromosome (default chr22) into small test BAMs
# - Indexes outputs
#
# Usage:
#   bash scripts/smoke_test_chr.sh \
#     --disease-bam <path-or-url> \
#     --control-bam <path-or-url> \
#     --download-local \
#     --chrom chr22 \
#     --outdir smoke-test

DISEASE_BAM=""
CONTROL_BAM=""
CHROM="chr22"
OUTDIR="smoke-test"
THREADS="${THREADS:-4}"
DOWNLOAD_LOCAL="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --disease-bam)
      DISEASE_BAM="$2"
      shift 2
      ;;
    --control-bam)
      CONTROL_BAM="$2"
      shift 2
      ;;
    --chrom)
      CHROM="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --download-local)
      DOWNLOAD_LOCAL="1"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${DISEASE_BAM}" || -z "${CONTROL_BAM}" ]]; then
  echo "ERROR: --disease-bam and --control-bam are required." >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools is required but not found in PATH." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

is_url() {
  local x="$1"
  [[ "${x}" =~ ^https?:// ]] || [[ "${x}" =~ ^ftp:// ]]
}

download_file() {
  local url="$1"
  local dest="$2"

  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --retry-delay 3 -o "${dest}" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${dest}" "${url}"
  else
    echo "ERROR: need curl or wget to download ${url}" >&2
    exit 1
  fi
}

ensure_local_bam_and_index() {
  local input_bam="$1"
  local label="$2"
  local dl_dir="${OUTDIR}/downloads"
  local out_bam=""

  mkdir -p "${dl_dir}"

  if is_url "${input_bam}"; then
    out_bam="${dl_dir}/${label}.bam"
    local out_bai="${out_bam}.bai"
    local bam_url="${input_bam}"
    local bai_url_a="${input_bam%.bam}.bai"
    local bai_url_b="${input_bam}.bai"

    echo "Downloading ${label} BAM..."
    download_file "${bam_url}" "${out_bam}"

    echo "Downloading ${label} BAM index..."
    if ! download_file "${bai_url_a}" "${out_bai}"; then
      download_file "${bai_url_b}" "${out_bai}"
    fi
    echo "${out_bam}"
    return
  fi

  if [[ ! -f "${input_bam}" ]]; then
    echo "ERROR: local BAM not found: ${input_bam}" >&2
    exit 1
  fi

  out_bam="${input_bam}"
  if [[ ! -f "${out_bam}.bai" && ! -f "${out_bam%.bam}.bai" ]]; then
    echo "Index not found for ${out_bam}; creating one..."
    samtools index -@ "${THREADS}" "${out_bam}"
  fi
  echo "${out_bam}"
}

slice_bam() {
  local in_bam="$1"
  local out_bam="$2"

  # This uses indexed random access if index is available.
  # For remote URLs, ensure corresponding .bai/.csi exists and server supports range requests.
  samtools view -@ "${THREADS}" -b "${in_bam}" "${CHROM}" -o "${out_bam}"
  samtools index -@ "${THREADS}" "${out_bam}"
}

if [[ "${DOWNLOAD_LOCAL}" == "1" ]]; then
  echo "Local download mode enabled."
  DISEASE_BAM="$(ensure_local_bam_and_index "${DISEASE_BAM}" "disease")"
  CONTROL_BAM="$(ensure_local_bam_and_index "${CONTROL_BAM}" "control")"
fi

echo "Extracting ${CHROM} from disease BAM..."
slice_bam "${DISEASE_BAM}" "${OUTDIR}/disease.${CHROM}.bam"

echo "Extracting ${CHROM} from control BAM..."
slice_bam "${CONTROL_BAM}" "${OUTDIR}/control.${CHROM}.bam"

echo "Done."
echo "disease subset:  ${OUTDIR}/disease.${CHROM}.bam"
echo "control subset: ${OUTDIR}/control.${CHROM}.bam"
echo "Suggested quick checks:"
echo "  samtools idxstats ${OUTDIR}/disease.${CHROM}.bam | head"
echo "  samtools idxstats ${OUTDIR}/control.${CHROM}.bam | head"
