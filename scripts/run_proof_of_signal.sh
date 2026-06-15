#!/usr/bin/env bash
set -euo pipefail

# End-to-end proof-of-signal pipeline (chr22-oriented starter flow):
# 1) extract split + discordant evidence
# 2) build candidate loci with junk-region flags
# 3) annotate MEI family/subfamily support from clipped reads
#
# This wrapper runs the same commands used in interactive development so the
# workflow can be reproduced with one command.

TUMOR_BAM=""
NORMAL_BAM=""
MEI_FASTA=""
REFERENCE_FASTA=""
G1K_MEI_BED=""
G1K_MEI_VCF=""
OUTDIR="results/mei_step1_hg38_chr22"
REGION="chr22"
WINDOW_SIZE="200"
G1K_SPLIT_PADDING_BP="200"
G1K_DPE_PADDING_MIN_BP="200"
G1K_DPE_PADDING_MAX_BP="200"
G1K_DPE_PADDING_TLEN_FACTOR="0"
PYTHON_BIN="${PYTHON_BIN:-python}"
RUN_IN_ENV="${RUN_IN_ENV:-0}" # set RUN_IN_ENV=1 to use `micromamba run -n rtm-miner ...`

SEG_DUP_BED="data/public/annotation/hg38/segdup/genomicSuperDups.bed"
MAPPABILITY_BEDGRAPH="data/public/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.bedGraph"
GAP_BED="data/public/annotation/hg38/masks/gap.txt.gz"
BLACKLIST_BED="data/public/annotation/hg38/blacklist/ENCFF356LFX.bed.gz"

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
    --mei-fasta)
      MEI_FASTA="$2"
      shift 2
      ;;
    --reference-fasta)
      REFERENCE_FASTA="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --g1k-mei-bed)
      G1K_MEI_BED="$2"
      shift 2
      ;;
    --g1k-mei-vcf)
      G1K_MEI_VCF="$2"
      shift 2
      ;;
    --region)
      REGION="$2"
      shift 2
      ;;
    --window-size)
      WINDOW_SIZE="$2"
      shift 2
      ;;
    --g1k-split-padding-bp)
      G1K_SPLIT_PADDING_BP="$2"
      shift 2
      ;;
    --g1k-dpe-padding-min-bp)
      G1K_DPE_PADDING_MIN_BP="$2"
      shift 2
      ;;
    --g1k-dpe-padding-max-bp)
      G1K_DPE_PADDING_MAX_BP="$2"
      shift 2
      ;;
    --g1k-dpe-padding-tlen-factor)
      G1K_DPE_PADDING_TLEN_FACTOR="$2"
      shift 2
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    --run-in-env)
      RUN_IN_ENV="1"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  if [[ "${PYTHON_BIN}" == "python" ]] && command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  else
    echo "ERROR: python executable not found: ${PYTHON_BIN}" >&2
    exit 1
  fi
fi

for required in "${TUMOR_BAM}" "${NORMAL_BAM}" "${MEI_FASTA}"; do
  if [[ -z "${required}" ]]; then
    echo "ERROR: missing required args: --tumor-bam, --normal-bam, --mei-fasta" >&2
    exit 1
  fi
done

for f in "${TUMOR_BAM}" "${NORMAL_BAM}" "${MEI_FASTA}" "${SEG_DUP_BED}" "${MAPPABILITY_BEDGRAPH}" "${GAP_BED}" "${BLACKLIST_BED}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required file not found: ${f}" >&2
    exit 1
  fi
done
if [[ -n "${REFERENCE_FASTA}" ]] && [[ ! -f "${REFERENCE_FASTA}" ]]; then
  echo "ERROR: reference FASTA not found: ${REFERENCE_FASTA}" >&2
  exit 1
fi
if [[ -n "${G1K_MEI_BED}" ]] && [[ ! -f "${G1K_MEI_BED}" ]]; then
  echo "ERROR: 1000G/MELT MEI BED not found: ${G1K_MEI_BED}" >&2
  exit 1
fi
if [[ -n "${G1K_MEI_VCF}" ]] && [[ ! -f "${G1K_MEI_VCF}" ]]; then
  echo "ERROR: 1000G/MELT MEI VCF not found: ${G1K_MEI_VCF}" >&2
  exit 1
fi
if [[ -n "${G1K_MEI_BED}" ]] && [[ -n "${G1K_MEI_VCF}" ]]; then
  echo "ERROR: provide only one of --g1k-mei-bed or --g1k-mei-vcf" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

run_cli() {
  local cmd=("$@")
  if [[ "${RUN_IN_ENV}" == "1" ]]; then
    micromamba run -n rtm-miner env PYTHONPATH=src "${PYTHON_BIN}" -m retro_miner.cli "${cmd[@]}"
  else
    PYTHONPATH=src "${PYTHON_BIN}" -m retro_miner.cli "${cmd[@]}"
  fi
}

echo "[proof-of-signal] stage=extract-split-evidence region=${REGION}"
run_cli extract-split-evidence \
  --tumor-bam "${TUMOR_BAM}" \
  --normal-bam "${NORMAL_BAM}" \
  --outdir "${OUTDIR}" \
  --region "${REGION}" \
  --min-mapq 20 \
  --min-mapq-discordant 0 \
  --min-clip-len 20

echo "[proof-of-signal] stage=build-candidate-loci window_size=${WINDOW_SIZE}"
run_cli build-candidate-loci \
  --evidence-dir "${OUTDIR}" \
  --window-size "${WINDOW_SIZE}" \
  --pseudocount 1.0 \
  --segdup-bed "${SEG_DUP_BED}" \
  --segdup-min-fraction 0.1 \
  --mappability-bedgraph "${MAPPABILITY_BEDGRAPH}" \
  --mappability-low-threshold 0.5 \
  --mappability-min-fraction 0.5 \
  --gap-bed "${GAP_BED}" \
  --gap-min-fraction 0.1 \
  --encode-blacklist-bed "${BLACKLIST_BED}" \
  --encode-blacklist-min-fraction 0.1

echo "[proof-of-signal] stage=annotate-mei-support"
annotate_cmd=(
  annotate-mei-support
  --evidence-dir "${OUTDIR}"
  --candidate-loci "${OUTDIR}/candidate_loci.tsv"
  --mei-fasta "${MEI_FASTA}"
  --tumor-bam-depth "${TUMOR_BAM}"
  --normal-bam-depth "${NORMAL_BAM}"
)
if [[ -n "${REFERENCE_FASTA}" ]]; then
  annotate_cmd+=(--reference-fasta "${REFERENCE_FASTA}")
fi
if [[ -n "${G1K_MEI_BED}" ]]; then
  annotate_cmd+=(
    --g1k-mei-bed "${G1K_MEI_BED}"
    --g1k-split-padding-bp "${G1K_SPLIT_PADDING_BP}"
    --g1k-dpe-padding-min-bp "${G1K_DPE_PADDING_MIN_BP}"
    --g1k-dpe-padding-max-bp "${G1K_DPE_PADDING_MAX_BP}"
    --g1k-dpe-padding-tlen-factor "${G1K_DPE_PADDING_TLEN_FACTOR}"
  )
fi
if [[ -n "${G1K_MEI_VCF}" ]]; then
  annotate_cmd+=(
    --g1k-mei-vcf "${G1K_MEI_VCF}"
    --g1k-split-padding-bp "${G1K_SPLIT_PADDING_BP}"
    --g1k-dpe-padding-min-bp "${G1K_DPE_PADDING_MIN_BP}"
    --g1k-dpe-padding-max-bp "${G1K_DPE_PADDING_MAX_BP}"
    --g1k-dpe-padding-tlen-factor "${G1K_DPE_PADDING_TLEN_FACTOR}"
  )
fi
annotate_cmd+=(--out-tsv "${OUTDIR}/candidate_loci.mei.tsv")
run_cli "${annotate_cmd[@]}"

echo "[proof-of-signal] done"
echo "  ${OUTDIR}/split_evidence.summary.tsv"
echo "  ${OUTDIR}/candidate_loci.tsv"
echo "  ${OUTDIR}/candidate_loci.mei.tsv"
