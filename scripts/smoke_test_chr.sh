#!/usr/bin/env bash
set -euo pipefail

# Smoke test:
# - Accepts tumor/normal BAM paths or HTTP/FTP URLs
# - Optional local download mode for remote BAMs (recommended)
# - Extracts one chromosome (default chr22) into small test BAMs
# - Indexes outputs
#
# Usage:
#   bash scripts/smoke_test_chr.sh \
#     --tumor-bam <path-or-url> \
#     --normal-bam <path-or-url> \
#     --download-local \
#     --chrom chr22 \
#     --outdir smoke-test

TUMOR_BAM=""
NORMAL_BAM=""
CHROM="chr22"
OUTDIR="smoke-test"
THREADS="${THREADS:-4}"
DOWNLOAD_LOCAL="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tumor-bam)
      TUMOR_BAM="$2"
      shift 2
      ;;
    --normal-bam)
      NORMAL_BAM="$2"
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

if [[ -z "${TUMOR_BAM}" || -z "${NORMAL_BAM}" ]]; then
  echo "ERROR: --tumor-bam and --normal-bam are required." >&2
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
  TUMOR_BAM="$(ensure_local_bam_and_index "${TUMOR_BAM}" "tumor")"
  NORMAL_BAM="$(ensure_local_bam_and_index "${NORMAL_BAM}" "normal")"
fi

echo "Extracting ${CHROM} from tumor BAM..."
slice_bam "${TUMOR_BAM}" "${OUTDIR}/tumor.${CHROM}.bam"

echo "Extracting ${CHROM} from normal BAM..."
slice_bam "${NORMAL_BAM}" "${OUTDIR}/normal.${CHROM}.bam"

echo "Done."
echo "Tumor subset:  ${OUTDIR}/tumor.${CHROM}.bam"
echo "Normal subset: ${OUTDIR}/normal.${CHROM}.bam"
echo "Suggested quick checks:"
echo "  samtools idxstats ${OUTDIR}/tumor.${CHROM}.bam | head"
echo "  samtools idxstats ${OUTDIR}/normal.${CHROM}.bam | head"
