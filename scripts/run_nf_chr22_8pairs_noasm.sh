#!/usr/bin/env bash
set -euo pipefail

# One-shot batch runner:
# - Streams unindexed S3 BAMs
# - Keeps only chr22 alignments locally
# - Indexes chr22 BAMs
# - Runs proof-of-signal with no local assembly
# - Runs at most 4 pairs in parallel

S3_PREFIX="${S3_PREFIX:-s3://l1tx-neurofibroma/synapse/syn4984617}"
MANIFEST="${MANIFEST:-/home/ec2-user/retrotransposon-miner/scripts/nf_chr22_8pairs_selected.tsv}"
WORKDIR="${WORKDIR:-/home/ec2-user/retrotransposon-workdir}"
REPO_DIR="${REPO_DIR:-/home/ec2-user/retrotransposon-miner}"
PAIR_CONCURRENCY="${PAIR_CONCURRENCY:-4}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"

ENV_BIN="${ENV_BIN:-/home/ec2-user/.local/share/mamba/envs/rtm-miner/bin}"
export PATH="${ENV_BIN}:$PATH"

RTM_PUBLIC_DATA_DIR="${RTM_PUBLIC_DATA_DIR:-${WORKDIR}/data/public}"
RTM_RESULTS_DIR="${RTM_RESULTS_DIR:-${WORKDIR}/results}"
RUNSTAMP="${RUNSTAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"

DATA_ROOT="${WORKDIR}/data/public/test_data/neurofibroma_chr22/${RUNSTAMP}"
OUT_ROOT="${WORKDIR}/results/mei_nf_chr22_8pairs_noasm_${RUNSTAMP}"
LOG_ROOT="${OUT_ROOT}/logs"
mkdir -p "${DATA_ROOT}" "${OUT_ROOT}" "${LOG_ROOT}"

MEI_FASTA="${MEI_FASTA:-${RTM_PUBLIC_DATA_DIR}/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta}"
REFERENCE_FASTA="${REFERENCE_FASTA:-${RTM_PUBLIC_DATA_DIR}/reference/hg19/Homo_sapiens_assembly19.fasta}"
RMSK_TABLE="${RMSK_TABLE:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/repeats/rmsk.txt.gz}"
G1K_MEI_VCF="${G1K_MEI_VCF:-${RTM_PUBLIC_DATA_DIR}/polymorphism/hg19/melt/nstd144.GRCh37.variant_call.vcf.gz}"
LR_MEI_VCF="${LR_MEI_VCF:-${RTM_PUBLIC_DATA_DIR}/polymorphism/hg19/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz}"
SEG_DUP_BED="${SEG_DUP_BED:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/segdup/genomicSuperDups.bed}"
MAPPABILITY_BEDGRAPH="${MAPPABILITY_BEDGRAPH:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/mappability/k100.Umap.MultiTrackMappability.bedGraph}"
MAPPABILITY_LOW_BED="${MAPPABILITY_LOW_BED:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed}"
GAP_BED="${GAP_BED:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/masks/gap.txt.gz}"
BLACKLIST_BED="${BLACKLIST_BED:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/blacklist/ENCFF356LFX.hg19.bed}"
JUNK_MERGED_BED="${JUNK_MERGED_BED:-${RTM_PUBLIC_DATA_DIR}/annotation/hg19/junk/junk_exclusion_merged.bed}"

for required_bin in aws samtools bash; do
  command -v "${required_bin}" >/dev/null 2>&1 || {
    echo "ERROR: missing required tool on PATH: ${required_bin}" >&2
    exit 1
  }
done

for required_file in \
  "${MANIFEST}" \
  "${MEI_FASTA}" \
  "${REFERENCE_FASTA}" \
  "${RMSK_TABLE}" \
  "${G1K_MEI_VCF}" \
  "${LR_MEI_VCF}" \
  "${SEG_DUP_BED}" \
  "${MAPPABILITY_BEDGRAPH}" \
  "${MAPPABILITY_LOW_BED}" \
  "${GAP_BED}" \
  "${BLACKLIST_BED}" \
  "${JUNK_MERGED_BED}"
do
  [[ -f "${required_file}" ]] || {
    echo "ERROR: required file not found: ${required_file}" >&2
    exit 1
  }
done

extract_chr22_from_s3_bam() {
  local s3_bam="$1"
  local out_bam="$2"
  local tmp_bam="${out_bam}.tmp"
  rm -f "${tmp_bam}"
  aws s3 cp "${s3_bam}" - \
    | samtools view -@ "${SAMTOOLS_THREADS}" -h - \
    | awk 'BEGIN{OFS="\t"} /^@/ || $3=="chr22"' \
    | samtools view -@ "${SAMTOOLS_THREADS}" -b -o "${tmp_bam}" -
  mv "${tmp_bam}" "${out_bam}"
}

run_one_pair() {
  local pair_id="$1"
  local disease_bam_name="$2"
  local control_bam_name="$3"
  local pair_data_dir="${DATA_ROOT}/${pair_id}"
  local pair_outdir="${OUT_ROOT}/${pair_id}"
  local pair_log="${LOG_ROOT}/${pair_id}.log"
  mkdir -p "${pair_data_dir}" "${pair_outdir}"

  {
    echo "[batch] start pair=${pair_id} ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "[batch] stream chr22 disease=${disease_bam_name}"
    extract_chr22_from_s3_bam "${S3_PREFIX}/${disease_bam_name}" "${pair_data_dir}/${disease_bam_name%.bam}.chr22.bam"
    echo "[batch] index disease chr22 bam"
    samtools index -@ "${SAMTOOLS_THREADS}" "${pair_data_dir}/${disease_bam_name%.bam}.chr22.bam"

    echo "[batch] stream chr22 control=${control_bam_name}"
    extract_chr22_from_s3_bam "${S3_PREFIX}/${control_bam_name}" "${pair_data_dir}/${control_bam_name%.bam}.chr22.bam"
    echo "[batch] index control chr22 bam"
    samtools index -@ "${SAMTOOLS_THREADS}" "${pair_data_dir}/${control_bam_name%.bam}.chr22.bam"

    echo "[batch] run proof-of-signal pair=${pair_id}"
    RUN_IN_ENV=1 \
    SEG_DUP_BED="${SEG_DUP_BED}" \
    MAPPABILITY_BEDGRAPH="${MAPPABILITY_BEDGRAPH}" \
    MAPPABILITY_LOW_BED="${MAPPABILITY_LOW_BED}" \
    GAP_BED="${GAP_BED}" \
    BLACKLIST_BED="${BLACKLIST_BED}" \
    JUNK_MERGED_BED="${JUNK_MERGED_BED}" \
    bash "${REPO_DIR}/scripts/run_proof_of_signal.sh" \
      --reference-build hg19 \
      --disease-bam "${pair_data_dir}/${disease_bam_name%.bam}.chr22.bam" \
      --control-bam "${pair_data_dir}/${control_bam_name%.bam}.chr22.bam" \
      --mei-fasta "${MEI_FASTA}" \
      --reference-fasta "${REFERENCE_FASTA}" \
      --rmsk-table "${RMSK_TABLE}" \
      --g1k-mei-vcf "${G1K_MEI_VCF}" \
      --lr-mei-vcf "${LR_MEI_VCF}" \
      --outdir "${pair_outdir}" \
      --chr chr22 \
      --no-local-assembly

    echo "[batch] done pair=${pair_id} ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "[batch] outputs=${pair_outdir}"
  } > "${pair_log}" 2>&1
}

mapfile -t PAIRS < <(awk -F'\t' 'NR>1 && NF>=3 {print $1"\t"$2"\t"$3}' "${MANIFEST}")
if [[ "${#PAIRS[@]}" -eq 0 ]]; then
  echo "ERROR: no pairs loaded from manifest: ${MANIFEST}" >&2
  exit 1
fi

echo "[batch] runstamp=${RUNSTAMP}"
echo "[batch] manifest=${MANIFEST}"
echo "[batch] pair_count=${#PAIRS[@]} pair_concurrency=${PAIR_CONCURRENCY}"
echo "[batch] data_root=${DATA_ROOT}"
echo "[batch] out_root=${OUT_ROOT}"
echo "[batch] logs=${LOG_ROOT}"

pids=()
labels=()

launch_pair() {
  local row="$1"
  IFS=$'\t' read -r pair_id disease_bam_name control_bam_name <<< "${row}"
  run_one_pair "${pair_id}" "${disease_bam_name}" "${control_bam_name}" &
  pids+=("$!")
  labels+=("${pair_id}")
  echo "[batch] launched pair=${pair_id} pid=${pids[-1]}"
}

reap_pairs() {
  local remaining_pids=()
  local remaining_labels=()
  local i pid label
  for i in "${!pids[@]}"; do
    pid="${pids[$i]}"
    label="${labels[$i]}"
    if kill -0 "${pid}" 2>/dev/null; then
      remaining_pids+=("${pid}")
      remaining_labels+=("${label}")
      continue
    fi
    if wait "${pid}"; then
      echo "[batch] completed pair=${label}"
    else
      echo "ERROR: pair failed=${label}; see ${LOG_ROOT}/${label}.log" >&2
      exit 1
    fi
  done
  pids=("${remaining_pids[@]}")
  labels=("${remaining_labels[@]}")
}

next_idx=0
total="${#PAIRS[@]}"
while [[ "${next_idx}" -lt "${total}" ]] || [[ "${#pids[@]}" -gt 0 ]]; do
  while [[ "${next_idx}" -lt "${total}" ]] && [[ "${#pids[@]}" -lt "${PAIR_CONCURRENCY}" ]]; do
    launch_pair "${PAIRS[$next_idx]}"
    next_idx=$((next_idx + 1))
  done
  reap_pairs
  if [[ "${#pids[@]}" -gt 0 ]]; then
    sleep 3
  fi
done

echo "[batch] all pairs complete"
echo "[batch] final_outputs=${OUT_ROOT}"
echo "[batch] final_logs=${LOG_ROOT}"
