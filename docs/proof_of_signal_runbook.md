# Proof-of-Signal Runbook

This runbook reproduces the current end-to-end prototype for tumor-vs-normal
MEI signal discovery on chr22-style test data.

## Inputs

- Tumor BAM on hg38 (chr22-sliced in current test flow)
- Normal BAM on hg38 (chr22-sliced)
- MEI reference FASTA (LINE1/Alu/SVA subset), e.g.:
  - `data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta`

## Prerequisites

From repo root:

```bash
bash scripts/bootstrap_env.sh
eval "$($HOME/.local/bin/micromamba shell hook -s bash)"
micromamba activate rtm-miner
bash scripts/validate_environment.sh
```

## One-command pipeline

Run from repo root:

```bash
bash scripts/run_proof_of_signal.sh \
  --tumor-bam data/public/test_data/seqc2/chr22/tumor.chr22.hg38.bam \
  --normal-bam data/public/test_data/seqc2/chr22/normal.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --outdir results/mei_step1_hg38_chr22 \
  --region chr22 \
  --window-size 200
```

If you are not activating the shell environment and prefer command-level env
execution, add:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --tumor-bam data/public/test_data/seqc2/chr22/tumor.chr22.hg38.bam \
  --normal-bam data/public/test_data/seqc2/chr22/normal.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta
```

## Outputs

- `results/mei_step1_hg38_chr22/split_evidence.summary.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.mei.tsv`

## Rough high-confidence count

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} ($c["junk_flag_count"]+0)==0 && ($c["mate_junk_flag_count"]+0)==0 && ($c["tumor_mei_supported_reads"]+0)>=2 && ($c["mei_score_enrichment_ratio"]+0)>=2 && ($c["enrichment_ratio"]+0)>1 && ($c["tumor_insertion_mei_span"]+0)>=30 {n++} END{print n+0}' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv
```

## Notes

- This prototype intentionally avoids hard filtering in core pipeline stages.
- Junk regions are annotated as flags so downstream review can filter or rank.
- MEI support is derived from split-clip alignments to MEI reference FASTA.
