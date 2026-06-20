# Proof-of-Signal Runbook

This runbook reproduces the current end-to-end prototype for disease-vs-control
MEI signal discovery on chr22-style test data.

## Inputs

- disease BAM on hg38 (chr22-sliced in current test flow)
- control BAM on hg38 (chr22-sliced)
- MEI reference FASTA (LINE1/Alu/SVA subset), e.g.:
  - `data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta`
- Reference genome FASTA (optional but recommended for TSD sequence and breakpoint
  context annotation), e.g.:
  - `data/public/reference/hg38/Homo_sapiens_assembly38.fasta`
- 1000G/MELT polymorphism VCF (optional; overlap + population frequency), e.g.:
  - `data/public/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz`
- 1000G ONT Vienna long-read SVAN polymorphism VCF (optional), e.g.:
  - `data/public/polymorphism/hg38/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz`
- RepeatMasker table (optional; nested insertion annotation), e.g.:
  - `data/public/annotation/hg38/repeats/rmsk.txt.gz`

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
  --disease-bam data/public/test_data/seqc2/chr22/disease.chr22.hg38.bam \
  --control-bam data/public/test_data/seqc2/chr22/control.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --outdir results/mei_step1_hg38_chr22 \
  --region chr22 \
  --window-size 200 \
  --local-assembly
```

### Add 1000G/MELT overlap annotation

Pass VCF directly (recommended). This is the preferred EC2 form using
repo-relative paths:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --disease-bam data/public/test_data/seqc2/chr22/disease.chr22.hg38.bam \
  --control-bam data/public/test_data/seqc2/chr22/control.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --g1k-mei-vcf data/public/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz \
  --outdir results/mei_step1_hg38_chr22_g1k \
  --region chr22 \
  --window-size 200 \
  --local-assembly
```

### Add both MELT + long-read cohort overlap annotation

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --disease-bam data/public/test_data/seqc2/chr22/disease.chr22.hg38.bam \
  --control-bam data/public/test_data/seqc2/chr22/control.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --g1k-mei-vcf data/public/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz \
  --lr-mei-vcf data/public/polymorphism/hg38/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz \
  --outdir results/mei_step1_hg38_chr22_g1k_lr \
  --region chr22 \
  --window-size 200 \
  --local-assembly
```

### EC2 command with empirical context scoring (N=1000)

Use this runbook command to include random-window empirical depth/read-quality
scoring during `annotate-mei-support`:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --disease-bam data/public/test_data/seqc2/chr22/disease.chr22.hg38.bam \
  --control-bam data/public/test_data/seqc2/chr22/control.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --g1k-mei-vcf data/public/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz \
  --outdir results/mei_step1_hg38_chr22_empirical \
  --region chr22 \
  --window-size 200 \
  --local-assembly
```

Defaults in the wrapper are:
- `empirical_random_windows=1000`
- `empirical_random_seed=13`
- `empirical_random_scope=chromosome`
- `local_assembly=on` (disable with `--no-local-assembly`)
- with `scope=chromosome`, this means `1000` random windows per chromosome in
  the run (for `scope=genome`, it is `1000` total)

Random windows are sampled outside candidate loci and outside the same junk
tracks used in candidate construction (segdup, low mappability, gap, ENCODE
blacklist), so no extra BED is required for the default flow.

Quick file check before running:

```bash
ls -lh data/public/polymorphism/hg38/melt/nstd144.GRCh38.variant_call.vcf.gz
```

### chr15 single-chromosome run

```bash
bash scripts/run_proof_of_signal.sh \
  --disease-bam data/public/test_data/seqc2/WGS_EA_T_1.bwa.dedup.bam \
  --control-bam data/public/test_data/seqc2/WGS_EA_N_1.bwa.dedup.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --outdir results/mei_step1_hg38_chr15 \
  --region chr15 \
  --window-size 200
```

If running via micromamba command-level execution:

```bash
RUN_IN_ENV=1 bash scripts/run_proof_of_signal.sh \
  --disease-bam data/public/test_data/seqc2/WGS_EA_T_1.bwa.dedup.bam \
  --control-bam data/public/test_data/seqc2/WGS_EA_N_1.bwa.dedup.bam \
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
  --disease-bam data/public/test_data/seqc2/chr22/disease.chr22.hg38.bam \
  --control-bam data/public/test_data/seqc2/chr22/control.chr22.hg38.bam \
  --mei-fasta data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta \
  --reference-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta
```

## Outputs

- `results/mei_step1_hg38_chr22/split_evidence.summary.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.tsv`
- `results/mei_step1_hg38_chr22/candidate_loci.mei.tsv`

Nested insertion columns in `candidate_loci.mei.tsv`:
- `nested_repeat_overlap`
- `nested_repeat_name`, `nested_repeat_class`, `nested_repeat_family`, `nested_repeat_strand`
- `nested_mei_family`, `nested_insertion_orientation`
- `nested_same_class`, `nested_same_orientation`, `nested_same_class_orientation`

## Runtime Note (EC2)

- On `r6.xlarge`, `build-candidate-loci` now typically completes in a couple of
  minutes for the chr22 test flow after the clustering fix.
- The step now emits `[candidate-loci]` progress logs (load, cluster, assign,
  overlap annotation, write), so long silent periods generally indicate an issue
  with an older code version or environment contention.

## Rough high-confidence count

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} ($c["junk_flag_count"]+0)==0 && ($c["mate_junk_flag_count"]+0)==0 && ($c["disease_mei_supported_reads"]+0)>=2 && ($c["mei_score_enrichment_ratio"]+0)>=2 && ($c["enrichment_ratio"]+0)>1 && ($c["disease_insertion_mei_span"]+0)>=30 {n++} END{print n+0}' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv
```

## Condensed review table (key fields)

Use this to print a compact table for shared likely-real calls, including inferred
insertion position/size, polyA support, TSD fields, and L1-like breakpoint context:

```bash
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++) c[$i]=i
  print "chrom","window_start","window_end","tier","score","coherence", \
        "disease_mei_supported_reads","control_mei_supported_reads","mei_ratio", \
        "insertion_breakpoint_pos","disease_insertion_mei_span","disease_poly_at_max_run","disease_poly_at_reads", \
        "disease_poly_at_fraction_weighted","tsd_detected","tsd_len_estimate","tsd_seq", \
        "breakpoint_context_11bp_oriented","breakpoint_l1_en_pattern_yy_rrrr","breakpoint_l1_en_orientation_source", \
        "breakpoint_l1_en_best_motif", \
        "breakpoint_l1_en_motif_type","breakpoint_l1_en_best_match_seq","breakpoint_l1_en_best_match_offset","breakpoint_yyrrrr_logodds","breakpoint_yyrrrr_logodds_shift1_max","breakpoint_yyrrrr_logodds_shift1_mt_adj","breakpoint_yyrrrr_best_offset","breakpoint_l1_en_motif_like", \
        "disease_subfamily","control_subfamily"
  next
}
{
  t = $c["disease_mei_supported_reads"]+0
  n = $c["control_mei_supported_reads"]+0
  s = $c["insertion_model_score"]+0
  coh = $c["coherence_score"]+0
  if (t>=2 && n>=2 && s>=0.50 && coh>=0.45) {
    print $c["chrom"],$c["window_start"],$c["window_end"],$c["insertion_call_tier"],s,coh, \
          t,n,$c["mei_score_enrichment_ratio"], \
          $c["insertion_breakpoint_pos"],$c["disease_insertion_mei_span"],$c["disease_poly_at_max_run"],$c["disease_poly_at_reads"], \
          $c["disease_poly_at_fraction_weighted"],$c["tsd_detected"],$c["tsd_len_estimate"],$c["tsd_seq"], \
          $c["breakpoint_context_11bp_oriented"],$c["breakpoint_l1_en_pattern_yy_rrrr"],$c["breakpoint_l1_en_orientation_source"],$c["breakpoint_l1_en_best_motif"], \
          $c["breakpoint_l1_en_motif_type"],$c["breakpoint_l1_en_best_match_seq"],$c["breakpoint_l1_en_best_match_offset"],$c["breakpoint_yyrrrr_logodds"],$c["breakpoint_yyrrrr_logodds_shift1_max"],$c["breakpoint_yyrrrr_logodds_shift1_mt_adj"],$c["breakpoint_yyrrrr_best_offset"],$c["breakpoint_l1_en_motif_like"], \
          $c["disease_R_mei_subfamily"],$c["control_R_mei_subfamily"]
  }
}
' OFS='\t' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv | sort -t$'\t' -k5,5gr
```

For chr15, replace the input path with:

`results/mei_step1_hg38_chr15/candidate_loci.mei.tsv`

## Strict disease-only DPE shortlist

Use this to keep only loci with bilateral disease DPE support (1+ each side),
perfect disease DPE family/strand consistency, no junk flag, and **zero** control
MEI support (split or DPE):

```bash
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++) c[$i]=i
  print "chrom","window_start","window_end","sample_status","tier", \
        "disease_split_total","control_split_total","disease_dpe_total","control_dpe_total", \
        "disease_dpe_left","disease_dpe_right","control_dpe_left","control_dpe_right", \
        "disease_dpe_family","disease_dpe_subfamily","disease_dpe_strand", \
        "disease_dpe_family_purity","disease_dpe_strand_purity", \
        "disease_anchor_mapq_mean","disease_anchor_mapq_min","control_anchor_mapq_mean","control_anchor_mapq_min"
  next
}
{
  tSplit = ($c["disease_L_mei_supported_reads"]+0) + ($c["disease_R_mei_supported_reads"]+0)
  nSplit = ($c["control_L_mei_supported_reads"]+0) + ($c["control_R_mei_supported_reads"]+0)
  tDpe = $c["disease_discordant_mei_supported_reads"]+0
  nDpe = $c["control_discordant_mei_supported_reads"]+0
  if (
    ($c["junk_flag_count"]+0)==0 &&
    ($c["disease_discordant_mei_left_supported_reads"]+0)>=1 &&
    ($c["disease_discordant_mei_right_supported_reads"]+0)>=1 &&
    ($c["disease_discordant_mei_family_purity"]+0)==1.0 &&
    ($c["disease_discordant_mei_strand_purity"]+0)==1.0 &&
    nSplit==0 &&
    nDpe==0
  ) {
    print $c["chrom"],$c["window_start"],$c["window_end"],$c["sample_status_label"],$c["insertion_call_tier"], \
          tSplit,nSplit,tDpe,nDpe, \
          $c["disease_discordant_mei_left_supported_reads"],$c["disease_discordant_mei_right_supported_reads"], \
          $c["control_discordant_mei_left_supported_reads"],$c["control_discordant_mei_right_supported_reads"], \
          $c["disease_discordant_mei_family"],$c["disease_discordant_mei_subfamily"],$c["disease_discordant_mei_strand"], \
          $c["disease_discordant_mei_family_purity"],$c["disease_discordant_mei_strand_purity"], \
          $c["discordant_disease_mapq_mean"],$c["discordant_disease_mapq_min"],$c["discordant_control_mapq_mean"],$c["discordant_control_mapq_min"]
  }
}
' OFS='\t' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv | sort -t$'\t' -k1,1 -k2,2n
```

## Symmetric DPE shortlist (keep shared, drop contradictory labels)

Use this to keep high-consistency bilateral DPE loci in clean regions while
removing contradictory single-sample labels:
- drop `disease_only` if any control support exists
- drop `control_only` if any disease support exists
- keep `shared`

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "chrom","window_start","window_end","sample_status","tier","disease_dpe_total","disease_dpe_left","disease_dpe_right","control_dpe_total","control_dpe_left","control_dpe_right","disease_dpe_family","control_dpe_family","disease_dpe_subfamily","control_dpe_subfamily","disease_dpe_strand","control_dpe_strand","disease_dpe_family_purity","control_dpe_family_purity","disease_dpe_strand_purity","control_dpe_strand_purity"; next} {tSplit=($c["disease_L_mei_supported_reads"]+0)+($c["disease_R_mei_supported_reads"]+0); nSplit=($c["control_L_mei_supported_reads"]+0)+($c["control_R_mei_supported_reads"]+0); tDpe=$c["disease_discordant_mei_supported_reads"]+0; nDpe=$c["control_discordant_mei_supported_reads"]+0; tAny=(tSplit+tDpe); nAny=(nSplit+nDpe); status=$c["sample_status_label"]; bilateral=((($c["disease_discordant_mei_left_supported_reads"]+0)>=1 && ($c["disease_discordant_mei_right_supported_reads"]+0)>=1 && ($c["disease_discordant_mei_family_purity"]+0)==1.0 && ($c["disease_discordant_mei_strand_purity"]+0)==1.0) || (($c["control_discordant_mei_left_supported_reads"]+0)>=1 && ($c["control_discordant_mei_right_supported_reads"]+0)>=1 && ($c["control_discordant_mei_family_purity"]+0)==1.0 && ($c["control_discordant_mei_strand_purity"]+0)==1.0)); drop_false_disease=((status=="disease_only") && nAny>=1); drop_false_control=((status=="control_only") && tAny>=1); if((($c["junk_flag_count"]+0)==0) && bilateral && !drop_false_disease && !drop_false_control){print $c["chrom"],$c["window_start"],$c["window_end"],status,$c["insertion_call_tier"],tDpe,$c["disease_discordant_mei_left_supported_reads"],$c["disease_discordant_mei_right_supported_reads"],nDpe,$c["control_discordant_mei_left_supported_reads"],$c["control_discordant_mei_right_supported_reads"],$c["disease_discordant_mei_family"],$c["control_discordant_mei_family"],$c["disease_discordant_mei_subfamily"],$c["control_discordant_mei_subfamily"],$c["disease_discordant_mei_strand"],$c["control_discordant_mei_strand"],$c["disease_discordant_mei_family_purity"],$c["control_discordant_mei_family_purity"],$c["disease_discordant_mei_strand_purity"],$c["control_discordant_mei_strand_purity"]}}' OFS='\t' results/mei_step1_hg38_chr22/candidate_loci.mei.tsv | sort -t$'\t' -k4,4 -k1,1 -k2,2n
```

## Notes

- This prototype intentionally avoids hard filtering in core pipeline stages.
- Junk regions are annotated as flags so downstream review can filter or rank.
- MEI support is derived from split-clip alignments to MEI reference FASTA.
