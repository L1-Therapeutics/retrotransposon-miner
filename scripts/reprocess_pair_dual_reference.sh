#!/usr/bin/env bash
set -euo pipefail

# Re-map tumor/normal BAM pairs to both hg38 and hs1/T2T references.
# This is the assembly-harmonization path for having BAMs on both references.
#
# NOTE:
# - liftOver is suitable for intervals/coordinates, not read alignments in BAM.
# - BAM-to-other-assembly conversion should be done by re-alignment, which this script performs.
#
# Required tools: samtools, pigz, and one of: bwa or bwa-mem2
#
# Example:
#   bash scripts/reprocess_pair_dual_reference.sh \
#     --tumor-bam data/public/test_data/seqc2/WGS_EA_T_1.bwa.dedup.bam \
#     --normal-bam data/public/test_data/seqc2/WGS_EA_N_1.bwa.dedup.bam \
#     --hg38-fasta data/public/reference/hg38/hg38.fa \
#     --hs1-fasta data/public/reference/hs1/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna \
#     --prefix seqc2_chr22 \
#     --outdir results/reprocessed_bams \
#     --threads 16

TUMOR_BAM=""
NORMAL_BAM=""
HG38_FASTA=""
HS1_FASTA=""
OUTDIR="results/reprocessed_bams"
THREADS="${THREADS:-8}"
KEEP_FASTQ="0"
PREFIX=""

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
    --hg38-fasta)
      HG38_FASTA="$2"
      shift 2
      ;;
    --hs1-fasta)
      HS1_FASTA="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --prefix)
      PREFIX="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --keep-fastq)
      KEEP_FASTQ="1"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

make_output_stem() {
  local sample="$1"
  local ref_name="$2"
  if [[ -n "${PREFIX}" ]]; then
    echo "${PREFIX}.${sample}.${ref_name}"
  else
    echo "${sample}.${ref_name}"
  fi
}

for required in "${TUMOR_BAM}" "${NORMAL_BAM}" "${HG38_FASTA}" "${HS1_FASTA}"; do
  if [[ -z "${required}" ]]; then
    echo "ERROR: missing required arguments. See script header for usage." >&2
    exit 1
  fi
done

for exe in samtools pigz; do
  if ! command -v "${exe}" >/dev/null 2>&1; then
    echo "ERROR: required executable not found in PATH: ${exe}" >&2
    exit 1
  fi
done

if command -v bwa >/dev/null 2>&1; then
  ALIGNER="bwa"
elif command -v bwa-mem2 >/dev/null 2>&1; then
  # Fallback only. bwa-mem2 requires its own index format.
  ALIGNER="bwa-mem2"
else
  echo "ERROR: required executable not found in PATH: bwa or bwa-mem2" >&2
  exit 1
fi

for f in "${TUMOR_BAM}" "${NORMAL_BAM}" "${HG38_FASTA}" "${HS1_FASTA}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: file not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${OUTDIR}/fastq" "${OUTDIR}/hg38" "${OUTDIR}/hs1" "${OUTDIR}/logs"

extract_fastq() {
  local in_bam="$1"
  local sample="$2"
  local r1="${OUTDIR}/fastq/${sample}.R1.fastq.gz"
  local r2="${OUTDIR}/fastq/${sample}.R2.fastq.gz"

  if [[ -f "${r1}" && -f "${r2}" ]]; then
    echo "[fastq] ${sample}: existing FASTQs found, skipping extraction"
    return
  fi

  echo "[fastq] ${sample}: extracting paired FASTQ"
  samtools collate -@ "${THREADS}" -Ou "${in_bam}" \
    | samtools fastq -@ "${THREADS}" -n \
      -1 "${r1}" -2 "${r2}" -0 /dev/null -s /dev/null
}

align_sample() {
  local sample="$1"
  local ref_name="$2"
  local ref_fa="$3"
  local r1="${OUTDIR}/fastq/${sample}.R1.fastq.gz"
  local r2="${OUTDIR}/fastq/${sample}.R2.fastq.gz"
  local stem
  stem="$(make_output_stem "${sample}" "${ref_name}")"
  local out_bam="${OUTDIR}/${ref_name}/${stem}.bam"
  local log="${OUTDIR}/logs/${stem}.align.log"

  if [[ -f "${out_bam}" && -f "${out_bam}.bai" ]]; then
    echo "[align] ${sample} -> ${ref_name}: existing BAM+BAI found, skipping"
    return
  fi

  echo "[align] ${sample} -> ${ref_name} with ${ALIGNER}"
  if [[ "${ALIGNER}" == "bwa-mem2" ]]; then
    # If bwa-mem2 is used, ensure mem2 index exists for this exact reference.
    if [[ ! -f "${ref_fa}.bwt.2bit.64" && ! -f "${ref_fa}.0123" ]]; then
      echo "ERROR: bwa-mem2 selected but mem2 index not found for ${ref_fa}." >&2
      echo "Use bwa index files with classic bwa, or build bwa-mem2 indexes first." >&2
      exit 1
    fi
  fi

  "${ALIGNER}" mem -t "${THREADS}" "${ref_fa}" "${r1}" "${r2}" 2> "${log}" \
    | samtools sort -@ "${THREADS}" -o "${out_bam}"
  samtools index -@ "${THREADS}" "${out_bam}"
}

extract_fastq "${TUMOR_BAM}" "tumor"
extract_fastq "${NORMAL_BAM}" "normal"

align_sample "tumor" "hg38" "${HG38_FASTA}"
align_sample "normal" "hg38" "${HG38_FASTA}"
align_sample "tumor" "hs1" "${HS1_FASTA}"
align_sample "normal" "hs1" "${HS1_FASTA}"

if [[ "${KEEP_FASTQ}" != "1" ]]; then
  rm -f "${OUTDIR}/fastq/tumor.R1.fastq.gz" "${OUTDIR}/fastq/tumor.R2.fastq.gz"
  rm -f "${OUTDIR}/fastq/normal.R1.fastq.gz" "${OUTDIR}/fastq/normal.R2.fastq.gz"
fi

echo "Done."
echo "Outputs:"
echo "  ${OUTDIR}/hg38/$(make_output_stem tumor hg38).bam"
echo "  ${OUTDIR}/hg38/$(make_output_stem normal hg38).bam"
echo "  ${OUTDIR}/hs1/$(make_output_stem tumor hs1).bam"
echo "  ${OUTDIR}/hs1/$(make_output_stem normal hs1).bam"
