# retrotransposon-miner

Retrotransposon analysis toolkit with an EC2-first onboarding path for large whole-genome workflows.

This README shows a new user how to:
- create and manage an AWS EC2 instance,
- sync this repository to the VM,
- install required tools on the VM,
- run JupyterLab over an SSH tunnel,
- connect to the VM with plain SSH or Cursor Remote SSH.

For whole-genome runs, use at least `r6i.4xlarge`.

## Quick Start

From your local machine:

```bash
git clone https://github.com/<org>/retrotransposon-miner.git
cd retrotransposon-miner
chmod +x scripts/ec2_jlab.sh
./scripts/ec2_jlab.sh bootstrap
```

After bootstrap finishes:
- SSH to the instance: `ssh retro-ec2`
- Open JupyterLab locally: `http://127.0.0.1:8890/lab?token=<printed-token>`

## Prerequisites

### Local machine tools

- `aws` CLI v2 (`aws configure` completed)
- `ssh`
- `curl`
- `git`

### AWS IAM permissions

The IAM user/role you run this with needs:

- `ec2:Describe*`
- `ec2:RunInstances`
- `ec2:StartInstances`
- `ec2:StopInstances`
- `ec2:RebootInstances`
- `ec2:CreateTags`
- `ec2:CreateKeyPair`
- `ec2:DeleteKeyPair` (optional)
- `ec2:CreateSecurityGroup`
- `ec2:AuthorizeSecurityGroupIngress`
- `ec2:AllocateAddress`
- `ec2:AssociateAddress`
- `ec2:DescribeAddresses`
- `ec2:DescribeVpcs`
- `ec2:DescribeSubnets`
- `ssm:GetParameter` (for latest Amazon Linux AMI lookup)
- `iam:PassRole` (only if attaching an instance profile)

## What `scripts/ec2_jlab.sh` Does

The script:
- discovers or creates an instance by `Name` tag (`retrotransposon-miner-wgs`),
- creates/reuses an SSH key pair in `~/.ssh/`,
- creates/reuses a security group allowing SSH only from your current public IP,
- starts the instance and waits for health checks,
- allocates/associates an Elastic IP,
- writes SSH aliases to your local `~/.ssh/config`:
  - `retro-ec2` for normal SSH
  - `jlab` for Jupyter tunnel forwarding (`8890 -> 8888`)
- starts JupyterLab on the VM and the local tunnel.

## 1) Provision EC2 and Start Jupyter Tunnel

```bash
cd retrotransposon-miner
./scripts/ec2_jlab.sh bootstrap
```

Useful operations:

```bash
./scripts/ec2_jlab.sh status
./scripts/ec2_jlab.sh stop-instance
./scripts/ec2_jlab.sh start-instance
./scripts/ec2_jlab.sh reboot-instance
./scripts/ec2_jlab.sh start-jlab
./scripts/ec2_jlab.sh stop-jlab
./scripts/ec2_jlab.sh start-tunnel
```

## 2) SSH Into the VM

```bash
ssh retro-ec2
```

The default user is `ec2-user` (configured automatically in your SSH alias).

## 3) Sync This Repository on the VM

On the VM:

```bash
cd ~
git clone https://github.com/<org>/retrotransposon-miner.git
cd retrotransposon-miner
```

To update later:

```bash
cd ~/retrotransposon-miner
git pull
```

## 4) Install Project Tools on the VM

On the VM, from the repo root:

```bash
bash scripts/bootstrap_env.sh
bash scripts/install_ucsc_tools.sh
conda activate rtm-miner || micromamba activate rtm-miner
bash scripts/validate_environment.sh
```

If `micromamba activate` fails in a fresh shell:

```bash
eval "$($HOME/.local/bin/micromamba shell hook -s bash)"
micromamba activate rtm-miner
```

To persist that:

```bash
echo 'eval "$($HOME/.local/bin/micromamba shell hook -s bash)"' >> ~/.bashrc
source ~/.bashrc
```

## 5) Download Public Reference Data Needed by the Algorithm

Run this on the VM after environment setup.

For hg38 workflow inputs:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

For hs1 workflow inputs:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hs1 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

For dual-build work (recommended if you evaluate both):

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 hs1 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

Notes:
- `install_ucsc_tools.sh` belongs in environment setup because downloader post-processing needs `liftOver`, `bigWigToBedGraph`, and `bigBedToBed`.
- `--references` controls which reference-centered resources are materialized.
- SEQC2 test BAMs are included automatically as chr22-sliced datasets (remote slicing, not full BAM download).

## 6) Prepare Test BAM Inputs

`download_public_data.py` now includes SEQC2 chr22 test BAMs automatically in the output tree.

If you want a separate/manual smoke-test subset from URLs, you can still run:

```bash
bash scripts/smoke_test_chr.sh \
  --disease-bam "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_T_1.bwa.dedup.bam" \
  --control-bam "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam" \
  --download-local \
  --chrom chr22 \
  --outdir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/smoke-test"
```

If you need test BAMs aligned on both hg38 and hs1, re-align the chr22 pair:

```bash
bash scripts/reprocess_pair_dual_reference.sh \
  --disease-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/disease.chr22.hg38.bam" \
  --control-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/control.chr22.hg38.bam" \
  --hg38-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/reference/hg38/Homo_sapiens_assembly38.fasta" \
  --hs1-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/reference/hs1/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa" \
  --prefix seqc2_chr22 \
  --outdir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/reprocessed_bams" \
  --threads 16
```

## 7) Connect JupyterLab From Your Laptop

`bootstrap` already starts JupyterLab and the tunnel.

If needed, restart manually on your laptop:

```bash
./scripts/ec2_jlab.sh start-jlab
./scripts/ec2_jlab.sh start-tunnel
```

Open:

`http://127.0.0.1:8890/lab?token=<token>`

You can get the token with:

```bash
ssh retro-ec2 'cat ~/.jlab/token'
```

## 8) Connect With Cursor

Because `bootstrap` writes SSH config, Cursor can reuse it directly.

In Cursor:
1. Open Command Palette.
2. Run `Remote-SSH: Connect to Host...`
3. Choose `retro-ec2`.
4. Open `~/retrotransposon-miner` on the remote host.

You can still use normal terminal SSH at any time with `ssh retro-ec2`.

## Optional: Keep Large Data Outside Git Checkout

On the VM:

```bash
source scripts/use_external_workdir.sh
```

This keeps large intermediate/results files outside the repository directory.

## Sanity Check Workflow

After environment setup on VM:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
bash scripts/validate_environment.sh
```

Then run public-data provisioning and smoke-test prep:

```bash
python3 scripts/download_public_data.py \
  --references hg38 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

## Notes

- `scripts/ec2_jlab.sh` assumes default VPC/subnet behavior and Amazon Linux 2023 AMI lookup through SSM.
- SSH ingress is restricted to your current public IP at bootstrap time; if your IP changes, rerun `bootstrap`.
- For production, review security groups, key handling, tagging, and cost controls according to your org standards.
