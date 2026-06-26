#!/usr/bin/env bash
set -euo pipefail

# End-to-end candidate-discovery-and-annotation pipeline (chr22-oriented starter flow):
# 1) extract split + discordant evidence
# 2) build candidate loci with junk-region flags
# 3) annotate MEI family/subfamily support from clipped reads
#
# This wrapper runs the same commands used in interactive development so the
# workflow can be reproduced with one command.

RTM_WORKDIR="${RTM_WORKDIR:-${HOME}/retrotransposon-workdir}"
RTM_PUBLIC_DATA_DIR="${RTM_PUBLIC_DATA_DIR:-${RTM_WORKDIR}/data/public}"
RTM_RESULTS_DIR="${RTM_RESULTS_DIR:-${RTM_WORKDIR}/results}"

DISEASE_BAM=""
CONTROL_BAM=""
MEI_FASTA=""
REFERENCE_FASTA=""
RMSK_TABLE=""
G1K_MEI_VCF=""
LR_MEI_VCF=""
REFERENCE_BUILD="hg38"
OUTDIR=""
REGION="chr22"
CHR_ARG=""
CHR_CONCURRENCY="6"
CHR_ALL_MODE="0"
# Default resume-safe behavior: skip chromosomes already complete in this outdir.
SKIP_COMPLETE_EXISTING="1"
WINDOW_SIZE="200"
G1K_SPLIT_PADDING_BP="200"
G1K_DPE_PADDING_MIN_BP="200"
G1K_DPE_PADDING_MAX_BP="200"
G1K_DPE_PADDING_TLEN_FACTOR="0"
EMPIRICAL_RANDOM_WINDOWS="1000"
EMPIRICAL_RANDOM_SCOPE="chromosome"
EMPIRICAL_RANDOM_SEED="13"
EMPIRICAL_HIGHCONF_BED=""
# Keep empirical gating enabled by default so noisy repeat-context loci are
# filtered from gold-stage outputs unless explicitly disabled.
EMPIRICAL_STAGE="1"
LOCAL_ASSEMBLY="1"
PYTHON_BIN="${PYTHON_BIN:-python}"
RUN_IN_ENV="${RUN_IN_ENV:-0}" # set RUN_IN_ENV=1 to use `micromamba run -n rtm-miner ...`

SEG_DUP_BED="${SEG_DUP_BED:-}"
MAPPABILITY_BEDGRAPH="${MAPPABILITY_BEDGRAPH:-}"
MAPPABILITY_LOW_BED="${MAPPABILITY_LOW_BED:-}"
GAP_BED="${GAP_BED:-}"
BLACKLIST_BED="${BLACKLIST_BED:-}"
JUNK_MERGED_BED="${JUNK_MERGED_BED:-}"

now_epoch() {
  date +%s
}

normalize_chr_token() {
  local raw="$1"
  local t="${raw#chr}"
  t="${t#CHR}"
  case "${t}" in
    X|x) echo "chrX" ;;
    Y|y) echo "chrY" ;;
    [0-9]|1[0-9]|2[0-2]) echo "chr${t}" ;;
    *)
      echo "ERROR: invalid chromosome token '${raw}' for --chr. Use 1-22,X,Y or chr1-chr22,chrX,chrY." >&2
      return 1
      ;;
  esac
}

resolve_chr_list() {
  local chr_arg="$1"
  if [[ -z "${chr_arg}" ]]; then
    # Backward-compatible single-region mode.
    echo "${REGION}"
    return 0
  fi
  if [[ "${chr_arg,,}" == "all" ]]; then
    local i
    for i in $(seq 22 -1 1); do
      echo "chr${i}"
    done
    echo "chrX"
    echo "chrY"
    return 0
  fi
  IFS=',' read -r -a raw_tokens <<< "${chr_arg}"
  local tok=""
  local norm=""
  for tok in "${raw_tokens[@]}"; do
    tok="${tok// /}"
    if [[ -z "${tok}" ]]; then
      continue
    fi
    norm="$(normalize_chr_token "${tok}")" || return 1
    echo "${norm}"
  done
}

set_reference_build_defaults() {
  case "${REFERENCE_BUILD}" in
    hg38)
      [[ -n "${REFERENCE_FASTA}" ]] || REFERENCE_FASTA="${RTM_PUBLIC_DATA_DIR}/reference/hg38/Homo_sapiens_assembly38.fasta"
      [[ -n "${RMSK_TABLE}" ]] || RMSK_TABLE="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/repeats/rmsk.txt.gz"
      [[ -n "${G1K_MEI_VCF}" ]] || G1K_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz"
      [[ -n "${LR_MEI_VCF}" ]] || LR_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hg38/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz"
      [[ -n "${SEG_DUP_BED}" ]] || SEG_DUP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/segdup/genomicSuperDups.bed"
      [[ -n "${MAPPABILITY_BEDGRAPH}" ]] || MAPPABILITY_BEDGRAPH="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.bedGraph"
      [[ -n "${MAPPABILITY_LOW_BED}" ]] || MAPPABILITY_LOW_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
      [[ -n "${GAP_BED}" ]] || GAP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/masks/gap.txt.gz"
      [[ -n "${BLACKLIST_BED}" ]] || BLACKLIST_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/blacklist/ENCFF356LFX.bed.gz"
      [[ -n "${JUNK_MERGED_BED}" ]] || JUNK_MERGED_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg38/junk/junk_exclusion_merged.bed"
      ;;
    hg19)
      [[ -n "${REFERENCE_FASTA}" ]] || REFERENCE_FASTA="${RTM_PUBLIC_DATA_DIR}/reference/hg19/Homo_sapiens_assembly19.fasta"
      [[ -n "${RMSK_TABLE}" ]] || RMSK_TABLE="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/repeats/rmsk.txt.gz"
      [[ -n "${G1K_MEI_VCF}" ]] || G1K_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hg19/melt/nstd144.GRCh37.variant_call.vcf.gz"
      [[ -n "${LR_MEI_VCF}" ]] || LR_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hg19/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz"
      [[ -n "${SEG_DUP_BED}" ]] || SEG_DUP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/segdup/genomicSuperDups.bed"
      [[ -n "${MAPPABILITY_BEDGRAPH}" ]] || MAPPABILITY_BEDGRAPH="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/mappability/k100.Umap.MultiTrackMappability.bedGraph"
      [[ -n "${MAPPABILITY_LOW_BED}" ]] || MAPPABILITY_LOW_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
      [[ -n "${GAP_BED}" ]] || GAP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/masks/gap.txt.gz"
      [[ -n "${BLACKLIST_BED}" ]] || BLACKLIST_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/blacklist/ENCFF356LFX.hg19.bed"
      [[ -n "${JUNK_MERGED_BED}" ]] || JUNK_MERGED_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hg19/junk/junk_exclusion_merged.bed"
      ;;
    hs1)
      [[ -n "${REFERENCE_FASTA}" ]] || REFERENCE_FASTA="${RTM_PUBLIC_DATA_DIR}/reference/hs1/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa"
      [[ -n "${RMSK_TABLE}" ]] || RMSK_TABLE=""
      [[ -n "${G1K_MEI_VCF}" ]] || G1K_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hs1/melt/nstd144.hs1.variant_call.vcf.gz"
      [[ -n "${LR_MEI_VCF}" ]] || LR_MEI_VCF="${RTM_PUBLIC_DATA_DIR}/polymorphism/hs1/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz"
      [[ -n "${SEG_DUP_BED}" ]] || SEG_DUP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/segdup/genomicSuperDups.hs1.bed"
      [[ -n "${MAPPABILITY_BEDGRAPH}" ]] || MAPPABILITY_BEDGRAPH="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/mappability/k100.Umap.MultiTrackMappability.bedGraph"
      [[ -n "${MAPPABILITY_LOW_BED}" ]] || MAPPABILITY_LOW_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
      [[ -n "${GAP_BED}" ]] || GAP_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/masks/gap.bed"
      [[ -n "${BLACKLIST_BED}" ]] || BLACKLIST_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/blacklist/ENCFF356LFX.hs1.bed"
      [[ -n "${JUNK_MERGED_BED}" ]] || JUNK_MERGED_BED="${RTM_PUBLIC_DATA_DIR}/annotation/hs1/junk/junk_exclusion_merged.bed"
      ;;
    *)
      echo "ERROR: unsupported --reference-build '${REFERENCE_BUILD}'. Use hg19, hg38, or hs1." >&2
      exit 1
      ;;
  esac
}

validate_reference_path_consistency() {
  local label="$1"
  local path="$2"
  if [[ -z "${path}" ]]; then
    return 0
  fi
  local p
  p="$(echo "${path}" | tr '[:upper:]' '[:lower:]')"
  case "${REFERENCE_BUILD}" in
    hg38)
      if [[ "${p}" == *"/hg19/"* ]] || [[ "${p}" == *"/hs1/"* ]] || [[ "${p}" == *"grch37"* ]] || [[ "${p}" == *"chm13"* ]]; then
        echo "ERROR: ${label} appears to target a non-hg38 build while --reference-build=hg38: ${path}" >&2
        exit 1
      fi
      ;;
    hg19)
      if [[ "${p}" == *"/hg38/"* ]] || [[ "${p}" == *"/hs1/"* ]] || [[ "${p}" == *"grch38"* ]] || [[ "${p}" == *"assembly38"* ]] || [[ "${p}" == *"chm13"* ]]; then
        echo "ERROR: ${label} appears to target a non-hg19 build while --reference-build=hg19: ${path}" >&2
        exit 1
      fi
      ;;
    hs1)
      if [[ "${p}" == *"/hg38/"* ]] || [[ "${p}" == *"/hg19/"* ]] || [[ "${p}" == *"grch38"* ]] || [[ "${p}" == *"grch37"* ]] || [[ "${p}" == *"assembly38"* ]] || [[ "${p}" == *"human_g1k_v37"* ]]; then
        echo "ERROR: ${label} appears to target a non-hs1 build while --reference-build=hs1: ${path}" >&2
        exit 1
      fi
      ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr)
      CHR_ARG="$2"
      shift 2
      ;;
    --chr_concurrency|--chr-concurrency)
      CHR_CONCURRENCY="$2"
      shift 2
      ;;
    --skip-complete|--skip-complete-in-outdir|--skip-complete-today)
      # --skip-complete-today is kept as a backward-compatible alias.
      SKIP_COMPLETE_EXISTING="1"
      shift 1
      ;;
    --no-skip-complete)
      SKIP_COMPLETE_EXISTING="0"
      shift 1
      ;;
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
    --reference-build|--reference)
      REFERENCE_BUILD="${2,,}"
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
    --empirical-stage)
      EMPIRICAL_STAGE="1"
      shift 1
      ;;
    --no-empirical-stage)
      EMPIRICAL_STAGE="0"
      shift 1
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

set_reference_build_defaults
if [[ -z "${OUTDIR}" ]]; then
  OUTDIR="${RTM_RESULTS_DIR}/mei_step1_${REFERENCE_BUILD}_chr22"
fi
validate_reference_path_consistency "reference FASTA" "${REFERENCE_FASTA}"
validate_reference_path_consistency "RepeatMasker table" "${RMSK_TABLE}"
validate_reference_path_consistency "1000G/MELT VCF" "${G1K_MEI_VCF}"
validate_reference_path_consistency "long-read SVAN VCF" "${LR_MEI_VCF}"
validate_reference_path_consistency "segdup BED" "${SEG_DUP_BED}"
validate_reference_path_consistency "mappability track" "${MAPPABILITY_BEDGRAPH}"
validate_reference_path_consistency "low-mappability mask" "${MAPPABILITY_LOW_BED}"
validate_reference_path_consistency "gap track" "${GAP_BED}"
validate_reference_path_consistency "blacklist track" "${BLACKLIST_BED}"
validate_reference_path_consistency "merged junk BED" "${JUNK_MERGED_BED}"

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  if [[ "${PYTHON_BIN}" == "python" ]] && command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  else
    echo "ERROR: python executable not found: ${PYTHON_BIN}" >&2
    exit 1
  fi
fi

PYTHON_BIN="$(command -v "${PYTHON_BIN}")"

if [[ "${RUN_IN_ENV}" == "0" ]]; then
  PYTHON_BIN_DIR="$(cd "$(dirname "${PYTHON_BIN}")" && pwd)"
  case ":${PATH}:" in
    *":${PYTHON_BIN_DIR}:"*) ;;
    *)
      export PATH="${PYTHON_BIN_DIR}:${PATH}"
      echo "[candidate-pipeline] prepended python bin dir to PATH for subprocess tools: ${PYTHON_BIN_DIR}"
      ;;
  esac
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
  echo "Re-run: python3 scripts/download_public_data.py --references ${REFERENCE_BUILD} --outdir ${RTM_PUBLIC_DATA_DIR}" >&2
  exit 1
fi
if [[ ! -f "${JUNK_MERGED_BED}" ]]; then
  echo "ERROR: required merged junk BED not found: ${JUNK_MERGED_BED}" >&2
  echo "Re-run: python3 scripts/download_public_data.py --references ${REFERENCE_BUILD} --outdir ${RTM_PUBLIC_DATA_DIR}" >&2
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

required_runtime_bins=(bedtools samtools)
if [[ "${LOCAL_ASSEMBLY}" == "1" ]]; then
  required_runtime_bins+=(minimap2 spades.py)
fi
for bin in "${required_runtime_bins[@]}"; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    echo "ERROR: required runtime tool not found on PATH: ${bin}" >&2
    echo "Current PATH=${PATH}" >&2
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

consolidate_all_chrom_outputs() {
  local base_outdir="$1"
  shift
  local chr_list=("$@")
  local basename="candidate_loci.mei.gold_review.tsv"
  local inputs=()
  local chr=""
  for chr in "${chr_list[@]}"; do
    local p="${base_outdir}/${chr}/${basename}"
    if [[ -f "${p}" ]]; then
      inputs+=("${p}")
    fi
  done
  if [[ "${#inputs[@]}" -eq 0 ]]; then
    echo "[candidate-pipeline] no per-chrom gold review tables found to consolidate"
    return 0
  fi
  local out_path="${base_outdir}/${basename}"
  local inputs_joined
  inputs_joined="$(printf "%s\n" "${inputs[@]}")"
  export RTM_MERGE_INPUTS="${inputs_joined}"
  export RTM_MERGE_OUTPUT="${out_path}"
  if [[ "${RUN_IN_ENV}" == "1" ]]; then
    micromamba run -n rtm-miner env PYTHONPATH=src "${PYTHON_BIN}" - <<'PY'
import os
from pathlib import Path
import pandas as pd
from retro_miner.mei_support import _prioritize_mei_candidates

inputs = [p for p in os.environ.get("RTM_MERGE_INPUTS", "").splitlines() if p.strip()]
out_path = Path(os.environ["RTM_MERGE_OUTPUT"])
if not inputs:
    raise SystemExit(0)

frames = [pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for p in inputs]
merged = pd.concat(frames, ignore_index=True)
for col in merged.columns:
    # Allow prioritizer to perform numeric coercions while preserving empty-string semantics.
    merged[col] = merged[col].where(merged[col] != "", other=pd.NA)

sorted_out = _prioritize_mei_candidates(merged, stage_first=True)
for col in sorted_out.columns:
    if sorted_out[col].dtype == bool:
        sorted_out[col] = sorted_out[col].astype(int)
sorted_out = sorted_out.fillna("")
out_path.parent.mkdir(parents=True, exist_ok=True)
sorted_out.to_csv(out_path, sep="\t", index=False)
PY
  else
    PYTHONPATH=src "${PYTHON_BIN}" - <<'PY'
import os
from pathlib import Path
import pandas as pd
from retro_miner.mei_support import _prioritize_mei_candidates

inputs = [p for p in os.environ.get("RTM_MERGE_INPUTS", "").splitlines() if p.strip()]
out_path = Path(os.environ["RTM_MERGE_OUTPUT"])
if not inputs:
    raise SystemExit(0)

frames = [pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) for p in inputs]
merged = pd.concat(frames, ignore_index=True)
for col in merged.columns:
    # Allow prioritizer to perform numeric coercions while preserving empty-string semantics.
    merged[col] = merged[col].where(merged[col] != "", other=pd.NA)

sorted_out = _prioritize_mei_candidates(merged, stage_first=True)
for col in sorted_out.columns:
    if sorted_out[col].dtype == bool:
        sorted_out[col] = sorted_out[col].astype(int)
sorted_out = sorted_out.fillna("")
out_path.parent.mkdir(parents=True, exist_ok=True)
sorted_out.to_csv(out_path, sep="\t", index=False)
PY
  fi
  unset RTM_MERGE_INPUTS
  unset RTM_MERGE_OUTPUT
  echo "[candidate-pipeline] consolidated ${basename} -> ${out_path}"
}

run_single_pipeline() {
  local run_region="$1"
  local run_outdir="$2"
  local stage_t0
  mkdir -p "${run_outdir}"

  stage_t0=$(now_epoch)
  echo "[candidate-pipeline] stage=extract-split-evidence region=${run_region} outdir=${run_outdir}"
  run_cli extract-split-evidence \
    --disease-bam "${DISEASE_BAM}" \
    --control-bam "${CONTROL_BAM}" \
    --outdir "${run_outdir}" \
    --region "${run_region}" \
    --min-mapq 20 \
    --min-mapq-discordant 20 \
    --min-clip-len 20
  echo "[candidate-pipeline] stage=extract-split-evidence done region=${run_region} elapsed=$(( $(now_epoch) - stage_t0 ))s"

  stage_t0=$(now_epoch)
  echo "[candidate-pipeline] stage=build-candidate-loci region=${run_region} window_size=${WINDOW_SIZE}"
  run_cli build-candidate-loci \
    --evidence-dir "${run_outdir}" \
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
  echo "[candidate-pipeline] stage=build-candidate-loci done region=${run_region} elapsed=$(( $(now_epoch) - stage_t0 ))s"

  stage_t0=$(now_epoch)
  echo "[candidate-pipeline] stage=annotate-mei-support region=${run_region}"
  annotate_cmd=(
    annotate-mei-support
    --evidence-dir "${run_outdir}"
    --candidate-loci "${run_outdir}/candidate_loci.tsv"
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
  if [[ "${EMPIRICAL_STAGE}" == "1" ]]; then
    annotate_cmd+=(--empirical-stage)
  else
    annotate_cmd+=(--no-empirical-stage)
  fi
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
  annotate_cmd+=(--out-tsv "${run_outdir}/candidate_loci.mei.tsv")
  run_cli "${annotate_cmd[@]}"
  echo "[candidate-pipeline] stage=annotate-mei-support done region=${run_region} elapsed=$(( $(now_epoch) - stage_t0 ))s"
}

if ! [[ "${CHR_CONCURRENCY}" =~ ^[0-9]+$ ]] || [[ "${CHR_CONCURRENCY}" -lt 1 ]]; then
  echo "ERROR: --chr_concurrency must be a positive integer." >&2
  exit 1
fi

mapfile -t CHR_LIST < <(resolve_chr_list "${CHR_ARG}")
if [[ "${#CHR_LIST[@]}" -eq 0 ]]; then
  echo "ERROR: resolved empty chromosome list from --chr '${CHR_ARG}'." >&2
  exit 1
fi

if [[ "${SKIP_COMPLETE_EXISTING}" == "1" ]] && [[ "${#CHR_LIST[@]}" -gt 1 ]]; then
  filtered_chr_list=()
  skipped_count=0
  for chr in "${CHR_LIST[@]}"; do
    chr_outdir="${OUTDIR}/${chr}"
    chr_log="${OUTDIR}/logs/${chr}.log"
    chr_gold="${chr_outdir}/candidate_loci.mei.gold_review.tsv"
    if [[ -f "${chr_log}" ]] && [[ -f "${chr_gold}" ]]; then
      if grep -q "stage=annotate-mei-support done region=${chr}" "${chr_log}"; then
        echo "[candidate-pipeline] skip-complete skipping ${chr} (already complete in outdir)"
        skipped_count=$((skipped_count + 1))
        continue
      fi
    fi
    filtered_chr_list+=("${chr}")
  done
  if [[ "${skipped_count}" -gt 0 ]]; then
    echo "[candidate-pipeline] skip-complete skipped=${skipped_count} remaining=${#filtered_chr_list[@]}"
  fi
  CHR_LIST=("${filtered_chr_list[@]}")
  if [[ "${#CHR_LIST[@]}" -eq 0 ]]; then
    echo "[candidate-pipeline] all requested chromosomes already complete in outdir; nothing to run"
    exit 0
  fi
fi

if [[ -n "${CHR_ARG}" ]] && [[ "${CHR_ARG,,}" == "all" ]]; then
  CHR_ALL_MODE="1"
fi

if [[ "${#CHR_LIST[@]}" -gt 1 ]] && [[ "${LOCAL_ASSEMBLY}" == "1" ]] && [[ "${CHR_CONCURRENCY}" -gt 1 ]]; then
  echo "ERROR: --chr_concurrency > 1 is not allowed when local assembly is enabled." >&2
  echo "Use --chr_concurrency 1 or add --no-local-assembly." >&2
  exit 1
fi

if [[ "${#CHR_LIST[@]}" -eq 1 ]]; then
  REGION="${CHR_LIST[0]}"
  run_single_pipeline "${REGION}" "${OUTDIR}"
  echo "[candidate-pipeline] done"
  echo "  ${OUTDIR}/split_evidence.summary.tsv"
  echo "  ${OUTDIR}/candidate_loci.tsv"
  echo "  ${OUTDIR}/candidate_loci.mei.tsv"
  echo "  ${OUTDIR}/candidate_loci.mei.gold_review.tsv"
  echo "  ${OUTDIR}/candidate_loci.mei.gold_review.igv/ (when reference + BAMs provided)"
  exit 0
fi

BASE_OUTDIR="${OUTDIR}"
LOG_DIR="${BASE_OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

echo "[candidate-pipeline] multi-chrom run chr_count=${#CHR_LIST[@]} chr_concurrency=${CHR_CONCURRENCY} base_outdir=${BASE_OUTDIR}"

# Rolling worker queue:
# - launch chromosomes until CHR_CONCURRENCY slots are full
# - when any child exits, immediately free the slot and launch next chromosome
# This avoids chunk-wide synchronization barriers and keeps workers busy.
pids=()
labels=()
next_idx=0
total="${#CHR_LIST[@]}"

launch_next_chr() {
  local chr="$1"
  local chr_outdir="${BASE_OUTDIR}/${chr}"
  local chr_log="${LOG_DIR}/${chr}.log"
  echo "[candidate-pipeline] launch ${chr} -> ${chr_outdir} (log=${chr_log})"
  (
    run_single_pipeline "${chr}" "${chr_outdir}"
  ) >"${chr_log}" 2>&1 &
  pids+=("$!")
  labels+=("${chr}")
}

reap_finished_children() {
  local remaining_pids=()
  local remaining_labels=()
  local idx
  local pid
  local chr
  for idx in "${!pids[@]}"; do
    pid="${pids[$idx]}"
    chr="${labels[$idx]}"
    if kill -0 "${pid}" 2>/dev/null; then
      remaining_pids+=("${pid}")
      remaining_labels+=("${chr}")
      continue
    fi
    if wait "${pid}"; then
      echo "[candidate-pipeline] completed ${chr}"
    else
      echo "ERROR: chromosome run failed for ${chr}. See ${LOG_DIR}/${chr}.log" >&2
      exit 1
    fi
  done
  pids=("${remaining_pids[@]}")
  labels=("${remaining_labels[@]}")
}

while [[ "${next_idx}" -lt "${total}" ]] || [[ "${#pids[@]}" -gt 0 ]]; do
  while [[ "${next_idx}" -lt "${total}" ]] && [[ "${#pids[@]}" -lt "${CHR_CONCURRENCY}" ]]; do
    launch_next_chr "${CHR_LIST[$next_idx]}"
    next_idx=$((next_idx + 1))
  done

  reap_finished_children
  if [[ "${#pids[@]}" -gt 0 ]]; then
    sleep 2
  fi
done

echo "[candidate-pipeline] done multi-chrom"
echo "  per-chrom outputs under: ${BASE_OUTDIR}/chr*/"
echo "  logs under: ${LOG_DIR}/"
if [[ "${CHR_ALL_MODE}" == "1" ]]; then
  echo "[candidate-pipeline] consolidating --chr all outputs into ${BASE_OUTDIR}"
  consolidate_all_chrom_outputs "${BASE_OUTDIR}" "${CHR_LIST[@]}"
  echo "  consolidated gold review written to: ${BASE_OUTDIR}/candidate_loci.mei.gold_review.tsv"
fi
