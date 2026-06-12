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
OUTDIR="results/mei_step1_hg38_chr22"
REGION="chr22"
WINDOW_SIZE="200"
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
    --outdir)
      OUTDIR="$2"
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
run_cli annotate-mei-support \
  --evidence-dir "${OUTDIR}" \
  --candidate-loci "${OUTDIR}/candidate_loci.tsv" \
  --mei-fasta "${MEI_FASTA}" \
  --out-tsv "${OUTDIR}/candidate_loci.mei.tsv"

echo "[proof-of-signal] done"
echo "  ${OUTDIR}/split_evidence.summary.tsv"
echo "  ${OUTDIR}/candidate_loci.tsv"
echo "  ${OUTDIR}/candidate_loci.mei.tsv"
