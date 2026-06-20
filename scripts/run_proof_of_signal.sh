#!/usr/bin/env bash
set -euo pipefail

# End-to-end proof-of-signal pipeline (chr22-oriented starter flow):
# 1) extract split + discordant evidence
# 2) build candidate loci with junk-region flags
# 3) annotate MEI family/subfamily support from clipped reads
#
# This wrapper runs the same commands used in interactive development so the
# workflow can be reproduced with one command.

DISEASE_BAM=""
CONTROL_BAM=""
MEI_FASTA=""
REFERENCE_FASTA=""
RMSK_TABLE="data/public/annotation/hg38/repeats/rmsk.txt.gz"
G1K_MEI_VCF=""
LR_MEI_VCF=""
OUTDIR="results/mei_step1_hg38_chr22"
REGION="chr22"
WINDOW_SIZE="200"
G1K_SPLIT_PADDING_BP="200"
G1K_DPE_PADDING_MIN_BP="200"
G1K_DPE_PADDING_MAX_BP="200"
G1K_DPE_PADDING_TLEN_FACTOR="0"
EMPIRICAL_RANDOM_WINDOWS="1000"
EMPIRICAL_RANDOM_SCOPE="chromosome"
EMPIRICAL_RANDOM_SEED="13"
EMPIRICAL_HIGHCONF_BED=""
LOCAL_ASSEMBLY="1"
PYTHON_BIN="${PYTHON_BIN:-python}"
RUN_IN_ENV="${RUN_IN_ENV:-0}" # set RUN_IN_ENV=1 to use `micromamba run -n rtm-miner ...`

SEG_DUP_BED="data/public/annotation/hg38/segdup/genomicSuperDups.bed"
MAPPABILITY_BEDGRAPH="data/public/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.bedGraph"
MAPPABILITY_LOW_BED="data/public/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
GAP_BED="data/public/annotation/hg38/masks/gap.txt.gz"
BLACKLIST_BED="data/public/annotation/hg38/blacklist/ENCFF356LFX.bed.gz"
JUNK_MERGED_BED="data/public/annotation/hg38/junk/junk_exclusion_merged.bed"

now_epoch() {
  date +%s
}

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
    --rmsk-table)
      RMSK_TABLE="$2"
      shift 2
      ;;
    --g1k-mei-vcf)
      G1K_MEI_VCF="$2"
      shift 2
      ;;
    --lr-mei-vcf)
      LR_MEI_VCF="$2"
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
    --empirical-random-windows)
      EMPIRICAL_RANDOM_WINDOWS="$2"
      shift 2
      ;;
    --empirical-random-scope)
      EMPIRICAL_RANDOM_SCOPE="$2"
      shift 2
      ;;
    --empirical-random-seed)
      EMPIRICAL_RANDOM_SEED="$2"
      shift 2
      ;;
    --empirical-highconf-bed)
      EMPIRICAL_HIGHCONF_BED="$2"
      shift 2
      ;;
    --local-assembly)
      LOCAL_ASSEMBLY="1"
      shift 1
      ;;
    --no-local-assembly)
      LOCAL_ASSEMBLY="0"
      shift 1
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

for required in "${DISEASE_BAM}" "${CONTROL_BAM}" "${MEI_FASTA}"; do
  if [[ -z "${required}" ]]; then
    echo "ERROR: missing required args: --disease-bam, --control-bam, --mei-fasta" >&2
    exit 1
  fi
done

for f in "${DISEASE_BAM}" "${CONTROL_BAM}" "${MEI_FASTA}" "${SEG_DUP_BED}" "${MAPPABILITY_BEDGRAPH}" "${GAP_BED}" "${BLACKLIST_BED}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required file not found: ${f}" >&2
    exit 1
  fi
done
if [[ ! -f "${MAPPABILITY_LOW_BED}" ]]; then
  echo "ERROR: required low-mappability BED not found: ${MAPPABILITY_LOW_BED}" >&2
  echo "Re-run: python3 scripts/download_public_data.py --outdir data/public" >&2
  exit 1
fi
if [[ ! -f "${JUNK_MERGED_BED}" ]]; then
  echo "ERROR: required merged junk BED not found: ${JUNK_MERGED_BED}" >&2
  echo "Re-run: python3 scripts/download_public_data.py --outdir data/public" >&2
  exit 1
fi
if [[ -n "${REFERENCE_FASTA}" ]] && [[ ! -f "${REFERENCE_FASTA}" ]]; then
  echo "ERROR: reference FASTA not found: ${REFERENCE_FASTA}" >&2
  exit 1
fi
if [[ -n "${RMSK_TABLE}" ]] && [[ ! -f "${RMSK_TABLE}" ]]; then
  echo "ERROR: RepeatMasker table not found: ${RMSK_TABLE}" >&2
  exit 1
fi
if [[ -n "${G1K_MEI_VCF}" ]] && [[ ! -f "${G1K_MEI_VCF}" ]]; then
  echo "ERROR: 1000G/MELT MEI VCF not found: ${G1K_MEI_VCF}" >&2
  exit 1
fi
if [[ -n "${LR_MEI_VCF}" ]] && [[ ! -f "${LR_MEI_VCF}" ]]; then
  echo "ERROR: long-read SVAN MEI VCF not found: ${LR_MEI_VCF}" >&2
  exit 1
fi
if [[ -n "${EMPIRICAL_HIGHCONF_BED}" ]] && [[ ! -f "${EMPIRICAL_HIGHCONF_BED}" ]]; then
  echo "ERROR: empirical high-confidence BED not found: ${EMPIRICAL_HIGHCONF_BED}" >&2
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

stage_t0=$(now_epoch)
echo "[proof-of-signal] stage=extract-split-evidence region=${REGION}"
run_cli extract-split-evidence \
  --disease-bam "${DISEASE_BAM}" \
  --control-bam "${CONTROL_BAM}" \
  --outdir "${OUTDIR}" \
  --region "${REGION}" \
  --min-mapq 20 \
  --min-mapq-discordant 20 \
  --min-clip-len 20
echo "[proof-of-signal] stage=extract-split-evidence done elapsed=$(( $(now_epoch) - stage_t0 ))s"

stage_t0=$(now_epoch)
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
echo "[proof-of-signal] stage=build-candidate-loci done elapsed=$(( $(now_epoch) - stage_t0 ))s"

stage_t0=$(now_epoch)
echo "[proof-of-signal] stage=annotate-mei-support"
annotate_cmd=(
  annotate-mei-support
  --evidence-dir "${OUTDIR}"
  --candidate-loci "${OUTDIR}/candidate_loci.tsv"
  --mei-fasta "${MEI_FASTA}"
  --disease-bam-depth "${DISEASE_BAM}"
  --control-bam-depth "${CONTROL_BAM}"
  --empirical-exclude-merged-bed "${JUNK_MERGED_BED}"
  --empirical-exclude-segdup-bed "${SEG_DUP_BED}"
  --empirical-exclude-mappability-bedgraph "${MAPPABILITY_LOW_BED}"
  --empirical-exclude-mappability-threshold 0.5
  --empirical-exclude-gap-bed "${GAP_BED}"
  --empirical-exclude-blacklist-bed "${BLACKLIST_BED}"
  --empirical-random-windows "${EMPIRICAL_RANDOM_WINDOWS}"
  --empirical-random-scope "${EMPIRICAL_RANDOM_SCOPE}"
  --empirical-random-seed "${EMPIRICAL_RANDOM_SEED}"
)
if [[ -n "${EMPIRICAL_HIGHCONF_BED}" ]]; then
  annotate_cmd+=(--empirical-highconf-bed "${EMPIRICAL_HIGHCONF_BED}")
fi
if [[ -n "${REFERENCE_FASTA}" ]]; then
  annotate_cmd+=(--reference-fasta "${REFERENCE_FASTA}")
fi
if [[ -n "${RMSK_TABLE}" ]]; then
  annotate_cmd+=(--rmsk-table "${RMSK_TABLE}")
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
if [[ -n "${LR_MEI_VCF}" ]]; then
  annotate_cmd+=(--lr-mei-vcf "${LR_MEI_VCF}")
fi
if [[ "${LOCAL_ASSEMBLY}" == "1" ]]; then
  annotate_cmd+=(--local-assembly)
fi
annotate_cmd+=(--out-tsv "${OUTDIR}/candidate_loci.mei.tsv")
run_cli "${annotate_cmd[@]}"
echo "[proof-of-signal] stage=annotate-mei-support done elapsed=$(( $(now_epoch) - stage_t0 ))s"

echo "[proof-of-signal] done"
echo "  ${OUTDIR}/split_evidence.summary.tsv"
echo "  ${OUTDIR}/candidate_loci.tsv"
echo "  ${OUTDIR}/candidate_loci.mei.tsv"
echo "  ${OUTDIR}/candidate_loci.mei.gold_review.tsv"
echo "  ${OUTDIR}/candidate_loci.mei.gold_review.igv/ (when reference + BAMs provided)"
