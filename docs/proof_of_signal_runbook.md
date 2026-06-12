# Proof-of-Signal Runbook

This runbook reproduces the current end-to-end prototype for tumor-vs-normal
MEI signal discovery on chr22-style test data.

## Inputs

- Tumor BAM on hg38 (chr22-sliced in current test flow)
- Normal BAM on hg38 (chr22-sliced)
- MEI reference FASTA (LINE1/Alu/SVA subset), e.g.:
  - `data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta`
- Reference genome FASTA (optional but recommended for TSD sequence and breakpoint
  context annotation), e.g.:
  - `data/public/reference/hg38/Homo_sapiens_assembly38.fasta`

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
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --outdir results/mei_step1_hg38_chr22 \
  --region chr22 \
  --window-size 200
```

### chr15 single-chromosome run

```bash
bash scripts/run_proof_of_signal.sh \
  --tumor-bam data/public/test_data/seqc2/WGS_EA_T_1.bwa.dedup.bam \
  --normal-bam data/public/test_data/seqc2/WGS_EA_N_1.bwa.dedup.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --outdir results/mei_step1_hg38_chr15 \
  --region chr15 \
  --window-size 200
```

If running via micromamba command-level execution:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --tumor-bam data/public/test_data/seqc2/WGS_EA_T_1.bwa.dedup.bam \
  --normal-bam data/public/test_data/seqc2/WGS_EA_N_1.bwa.dedup.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --outdir results/mei_step1_hg38_chr15 \
  --region chr15 \
  --window-size 200
```

If you are not activating the shell environment and prefer command-level env
execution, add:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --tumor-bam data/public/test_data/seqc2/chr22/tumor.chr22.hg38.bam \
  --normal-bam data/public/test_data/seqc2/chr22/normal.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta
```

## Outputs

- `results/mei_step1_hg38_chr22/split_evidence.summary.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.mei.tsv`

## Runtime Note (EC2)

- On `r6.xlarge`, `build-candidate-loci` now typically completes in a couple of
  minutes for the chr22 test flow after the clustering fix.
- The step now emits `[candidate-loci]` progress logs (load, cluster, assign,
  overlap annotation, write), so long silent periods generally indicate an issue
  with an older code version or environment contention.

## Rough high-confidence count

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} ($c["junk_flag_count"]+0)==0 && ($c["mate_junk_flag_count"]+0)==0 && ($c["tumor_mei_supported_reads"]+0)>=2 && ($c["mei_score_enrichment_ratio"]+0)>=2 && ($c["enrichment_ratio"]+0)>1 && ($c["tumor_insertion_mei_span"]+0)>=30 {n++} END{print n+0}' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv
```

## Condensed review table (key fields)

Use this to print a compact table for shared likely-real calls, including inferred
insertion position/size, polyA support, TSD fields, and L1-like breakpoint context:

```bash
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++) c[$i]=i
  print "chrom","window_start","window_end","tier","score","coherence", \
        "tumor_mei_supported_reads","normal_mei_supported_reads","mei_ratio", \
        "tumor_insertion_breakpoint_pos","tumor_insertion_mei_span","tumor_poly_at_max_run","tumor_poly_at_reads", \
        "tumor_poly_at_fraction_weighted","tsd_detected","tsd_len_estimate","tsd_seq", \
        "tumor_breakpoint_context_11bp_oriented","tumor_breakpoint_l1_en_pattern_yy_rrrr","tumor_breakpoint_l1_en_orientation_source", \
        "tumor_breakpoint_l1_en_best_motif", \
        "tumor_breakpoint_l1_en_motif_type","tumor_breakpoint_l1_en_best_match_seq","tumor_breakpoint_l1_en_best_match_offset","tumor_breakpoint_yyrrrr_logodds","tumor_breakpoint_yyrrrr_logodds_shift1_max","tumor_breakpoint_yyrrrr_logodds_shift1_mt_adj","tumor_breakpoint_yyrrrr_best_offset","tumor_breakpoint_l1_en_motif_like", \
        "tumor_subfamily","normal_subfamily"
  next
}
{
  t = $c["tumor_mei_supported_reads"]+0
  n = $c["normal_mei_supported_reads"]+0
  s = $c["insertion_model_score"]+0
  coh = $c["coherence_score"]+0
  if (t>=2 && n>=2 && s>=0.50 && coh>=0.45) {
    print $c["chrom"],$c["window_start"],$c["window_end"],$c["insertion_call_tier"],s,coh, \
          t,n,$c["mei_score_enrichment_ratio"], \
          $c["tumor_insertion_breakpoint_pos"],$c["tumor_insertion_mei_span"],$c["tumor_poly_at_max_run"],$c["tumor_poly_at_reads"], \
          $c["tumor_poly_at_fraction_weighted"],$c["tsd_detected"],$c["tsd_len_estimate"],$c["tsd_seq"], \
          $c["tumor_breakpoint_context_11bp_oriented"],$c["tumor_breakpoint_l1_en_pattern_yy_rrrr"],$c["tumor_breakpoint_l1_en_orientation_source"],$c["tumor_breakpoint_l1_en_best_motif"], \
          $c["tumor_breakpoint_l1_en_motif_type"],$c["tumor_breakpoint_l1_en_best_match_seq"],$c["tumor_breakpoint_l1_en_best_match_offset"],$c["tumor_breakpoint_yyrrrr_logodds"],$c["tumor_breakpoint_yyrrrr_logodds_shift1_max"],$c["tumor_breakpoint_yyrrrr_logodds_shift1_mt_adj"],$c["tumor_breakpoint_yyrrrr_best_offset"],$c["tumor_breakpoint_l1_en_motif_like"], \
          $c["tumor_R_mei_subfamily"],$c["normal_R_mei_subfamily"]
  }
}
' OFS='\t' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv | sort -t$'\t' -k5,5gr
```

For chr15, replace the input path with:

`results/mei_step1_hg38_chr15/candidate_loci.mei.tsv`

## Notes

- This prototype intentionally avoids hard filtering in core pipeline stages.
- Junk regions are annotated as flags so downstream review can filter or rank.
- MEI support is derived from split-clip alignments to MEI reference FASTA.
