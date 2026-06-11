# retrotransposon-miner

Initial scaffold for a high-sensitivity retrotransposon insertion caller.

This first iteration focuses on reproducible environment setup for AWS Linux,
with tooling that is suitable for large BAM processing and permissive for
commercial use.

## Scope of this phase

- Reproducible conda environment with core bioinformatics + Python tools
- Command-line entry point scaffold for the future caller
- Environment validation script to verify required binaries and Python modules
- No calling logic yet (implemented in the next phase)

## Prerequisites

- Linux x86_64 VM (recommended for production runs)
- `curl` (bootstrap script installs `micromamba` automatically if needed)

## Create the environment (auto-installs package manager if missing)

```bash
bash scripts/bootstrap_env.sh
```

Then activate:

```bash
conda activate rtm-miner
```

If your shell only has micromamba available:

```bash
micromamba activate rtm-miner
```

## Validate installation

```bash
bash scripts/validate_environment.sh
```

This checks:
- Required command-line tools (`samtools`, `bedtools`, `bcftools`, `bwa-mem2`)
- Python runtime and key modules (`pysam`, `pandas`, `pyarrow`, `click`, `Bio`)
- Optional UCSC binary (`liftOver`)

## Install UCSC liftOver (optional but recommended)

```bash
bash scripts/install_ucsc_tools.sh
```

The installer tries to use your package manager first (preferred), then falls
back to standalone UCSC binaries in `~/.local/ucsc`.

If fallback binaries were used:

```bash
export PATH="$HOME/.local/ucsc:$PATH"
```

To persist PATH across sessions, add that export line to your shell profile.

If you see:

`liftOver: /lib64/libcurl.so.4: no version information available`

this usually comes from standalone binary/library mismatch on some Linux images.
It is often non-fatal, but package-manager installation is preferred to avoid it.

## Copy code to AWS VM and connect

From your local machine:

```bash
scp -i /Users/williambrandler/Dropbox/Mac/Documents/git/l1tx-data-pipelines/scripts/synapse-download-key.pem -r /Users/williambrandler/Dropbox/Mac/Documents/git/retrotransposon-miner ec2-user@44.204.61.87:~/
```

Connect to the VM:

```bash
ssh -i /Users/williambrandler/Dropbox/Mac/Documents/git/l1tx-data-pipelines/scripts/synapse-download-key.pem ec2-user@44.204.61.87
```

On the VM:

```bash
cd ~/retrotransposon-miner
bash scripts/bootstrap_env.sh
bash scripts/install_ucsc_tools.sh
export PATH="$HOME/.local/ucsc:$PATH"
conda activate rtm-miner || micromamba activate rtm-miner
bash scripts/validate_environment.sh
```

If `micromamba activate` complains about shell initialization, run this first in
that shell:

```bash
eval "$($HOME/.local/bin/micromamba shell hook -s bash)"
micromamba activate rtm-miner
```

Then add it to `~/.bashrc` (recommended on EC2) so future shells work directly:

```bash
echo 'eval "$($HOME/.local/bin/micromamba shell hook -s bash)"' >> ~/.bashrc
source ~/.bashrc
```

## Public tumor/normal BAM candidates for smoke tests

True matched tumor/normal WGS BAMs that are fully public in AWS S3 are uncommon
because many patient datasets are controlled-access. For immediate smoke testing,
use SEQC2 public benchmark BAMs (widely used, open access) from NCBI FTP:

- Tumor BAM:
  `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bam`
- Normal BAM:
  `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam`

These are large files; for quick iteration, subset one chromosome first.

## Smoke test (chr22 subset)

Recommended: download BAM + BAI locally first, then subset `chr22`.

```bash
bash scripts/smoke_test_chr.sh \
  --tumor-bam "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bam" \
  --normal-bam "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam" \
  --download-local \
  --chrom chr22 \
  --outdir smoke-test
```

Outputs:
- Downloaded inputs: `smoke-test/downloads/tumor.bam`, `smoke-test/downloads/normal.bam`
- Subsets: `smoke-test/tumor.chr22.bam`, `smoke-test/normal.chr22.bam`

## Download public reference/annotation resources

```bash
python3 scripts/download_public_data.py \
  --config resources/public_datasets.yaml \
  --outdir data/public
```

This initial dataset pack includes:
- hg38 reference FASTA (Broad AWS)
- T2T/hs1 reference FASTA (CHM13 analysis set on AWS)
- FASTA index artifacts generated locally (`.fai`, `.dict`)
- GENCODE annotation (v46, GRCh38)
- RepeatMasker + segmental duplication tables (hg38)
- UCSC 100bp mappability (`k100.Umap.MultiTrackMappability.bw`)
- gnomAD SV polymorphism track (includes MEI classes)
- hg38<->hs1 liftOver chains
- Dfam FamDB archive (required; extract human families with `famdb.py`)
- Public SEQC2 tumor/normal test BAMs (downloaded as chr22 slices)

Post-download processing is automatic and includes:
- remote slicing of SEQC2 tumor/normal BAMs to `chr22` (no full BAM local copy)
- converting selected hg38 resources to BED/BEDGRAPH
- extracting MEI-focused BED rows from gnomAD SV track
- lifting hg38 BED resources to hs1/T2T via `liftOver`
- downloading prebuilt BWA index files for hg38 and hs1 from AWS sources when available (fallback: local build)

Storage planning:
- Full whole-genome tumor+normal BAM remapping workflows can require ~300GB free disk.
- Default test-data download mode is chr22-sliced to keep storage much lower.

Download only selected categories:

```bash
python3 scripts/download_public_data.py \
  --categories reference liftover \
  --outdir data/public
```

Skip post-processing if you only want raw downloads:

```bash
python3 scripts/download_public_data.py \
  --skip-postprocess \
  --outdir data/public
```

Force local BWA index build (disable prebuilt-hg38 index download):

```bash
python3 scripts/download_public_data.py \
  --no-prebuilt-bwa-index \
  --outdir data/public
```

The dataset catalog now includes HPRC, 1KGP/MELT, GIAB index docs, ENCODE
blacklist, simple repeats, MANE, and hs1 add-ons as required downloads.

Notes on special cases:
- `dbRIP` and `euL1db` are currently portal-style/manual exports, not stable
  direct-file endpoints.
- `Repbase` requires a commercial license for proprietary use; keep it out of
  automated public-download workflows unless license terms are satisfied.

## Reprocess tumor/normal BAMs into both hg38 and hs1/T2T

To keep read alignments available on both assemblies, re-align the downloaded
test pair to each reference:

```bash
bash scripts/reprocess_pair_dual_reference.sh \
  --tumor-bam data/public/test_data/seqc2/chr22/tumor.chr22.hg38.bam \
  --normal-bam data/public/test_data/seqc2/chr22/normal.chr22.hg38.bam \
  --hg38-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --hs1-fasta data/public/reference/hs1/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa \
  --prefix seqc2_chr22 \
  --outdir results/reprocessed_bams \
  --threads 16
```

This uses `bwa-mem2` (or `bwa` fallback) for WGS read re-alignment.
It is a re-alignment step (not coordinate liftover), which is the appropriate
way to convert BAM evidence between assemblies.

### Reproducible chr22 test outputs

To regenerate the same chr22 test BAM set reproducibly from scratch:

```bash
python3 scripts/download_public_data.py --outdir data/public --threads 4

bash scripts/reprocess_pair_dual_reference.sh \
  --tumor-bam data/public/test_data/seqc2/chr22/tumor.chr22.hg38.bam \
  --normal-bam data/public/test_data/seqc2/chr22/normal.chr22.hg38.bam \
  --hg38-fasta data/public/reference/hg38/Homo_sapiens_assembly38.fasta \
  --hs1-fasta data/public/reference/hs1/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa \
  --prefix seqc2_chr22 \
  --outdir results/reprocessed_bams \
  --threads 16
```

Expected test outputs:
- `results/reprocessed_bams/hg38/seqc2_chr22.tumor.hg38.bam`
- `results/reprocessed_bams/hg38/seqc2_chr22.normal.hg38.bam`
- `results/reprocessed_bams/hs1/seqc2_chr22.tumor.hs1.bam`
- `results/reprocessed_bams/hs1/seqc2_chr22.normal.hs1.bam`

## Quick smoke test on a small chromosome

Use one chromosome (for example `chr22`) to iterate quickly before full-genome
runs:

```bash
# Example workflow shape only:
# 1) subset input BAM to chr22
# 2) index subset BAM
# 3) run downstream extraction on subset
samtools view -b INPUT.bam chr22 > subset.chr22.bam
samtools index subset.chr22.bam
```

## Project layout

```text
environment.yml
scripts/
  download_public_data.py
  smoke_test_chr.sh
  validate_environment.sh
resources/
  public_datasets.yaml
src/
  retro_miner/
    __init__.py
    cli.py
```

## Next step

After environment setup is confirmed, the next phase is to add:

- public reference/resource download + version manifest scripts
- discordant/split read extraction module
- table schema for candidate insertions and processed pseudogene events
