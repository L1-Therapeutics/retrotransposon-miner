# Databricks notebook source
# DBTITLE 1,Bootstrap environment
import os
import subprocess
from pathlib import Path

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
result = subprocess.run(
    ["bash", "scripts/bootstrap_env.sh"],
    cwd=repo_dir,
    text=True,
    capture_output=True,
    check=False,
)
print(result.stdout)
if result.stderr:
    print(result.stderr)
if result.returncode != 0:
    raise RuntimeError(f"bootstrap_env.sh failed with exit code {result.returncode}")

# COMMAND ----------

# DBTITLE 1,Validate environment
import os
import subprocess
from pathlib import Path

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
micromamba = os.path.expanduser("~/.local/bin/micromamba")
cmd = [
    micromamba,
    "run",
    "-n",
    "rtm-miner",
    "bash",
    "scripts/validate_environment.sh",
]
result = subprocess.run(
    cmd,
    cwd=repo_dir,
    text=True,
    capture_output=True,
    check=False,
)
print(result.stdout)
if result.stderr:
    print(result.stderr)
if result.returncode != 0 and "All required tools detected." not in result.stdout:
    raise RuntimeError(f"validate_environment.sh failed with exit code {result.returncode}")

# COMMAND ----------

# DBTITLE 1,Configure external workdir
import os
import shutil
from pathlib import Path

from azure.core.exceptions import ResourceExistsError, ResourceNotFoundError
from azure.identity import ClientSecretCredential
from azure.storage.filedatalake import DataLakeServiceClient


storage_account = dbutils.widgets.get("rtm_storage_account").strip()
container = dbutils.widgets.get("rtm_container").strip() or "l1tx-rtm-miner-output"
tenant_id = dbutils.widgets.get("rtm_tenant_id_key").strip()
client_id = dbutils.widgets.get("rtm_client_id_key").strip()
client_secret = dbutils.widgets.get("rtm_client_secret_key").strip()
keep_intermediates = dbutils.widgets.get("rtm_keep_intermediates").strip() == "1"

missing = [
    name
    for name, value in [
        ("rtm_storage_account", storage_account),
        ("rtm_container", container),
        ("rtm_tenant_id_key", tenant_id),
        ("rtm_client_id_key", client_id),
        ("rtm_client_secret_key", client_secret),
    ]
    if not value
]
if missing:
    raise ValueError("Missing required parameters: " + ", ".join(missing))

account_fqdn = f"{storage_account}.dfs.core.windows.net"
account_url = f"https://{account_fqdn}"
workdir_prefix = "retrotransposon-workdir"
source_uri = f"abfss://{container}@{account_fqdn}/"
rtm_persistent_workdir = f"{source_uri}{workdir_prefix}"
rtm_public_data_dir = f"{rtm_persistent_workdir}/data/public"
rtm_results_dir = f"{rtm_persistent_workdir}/results"
local_workdir = Path("/local_disk0/tmp/retrotransposon-workdir")
local_public_data_dir = local_workdir / "data" / "public"
local_results_dir = local_workdir / "results"
local_logs_dir = local_workdir / "logs"


def ensure_local_dirs(root: Path) -> None:
    (root / "data" / "public").mkdir(parents=True, exist_ok=True)
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "logs").mkdir(parents=True, exist_ok=True)


def normalize_remote_path(path: str) -> str:
    return path.strip("/")


credential = ClientSecretCredential(tenant_id=tenant_id, client_id=client_id, client_secret=client_secret)
service_client = DataLakeServiceClient(account_url=account_url, credential=credential)
file_system_client = service_client.get_file_system_client(container)

spark.conf.set(f"fs.azure.account.auth.type.{account_fqdn}", "OAuth")
spark.conf.set(
    f"fs.azure.account.oauth.provider.type.{account_fqdn}",
    "org.apache.hadoop.fs.azurebfs.oauth2.ClientCredsTokenProvider",
)
spark.conf.set(f"fs.azure.account.oauth2.client.id.{account_fqdn}", client_id)
spark.conf.set(f"fs.azure.account.oauth2.client.secret.{account_fqdn}", client_secret)
spark.conf.set(
    f"fs.azure.account.oauth2.client.endpoint.{account_fqdn}",
    f"https://login.microsoftonline.com/{tenant_id}/oauth2/token",
)


for remote_dir in [workdir_prefix, f"{workdir_prefix}/data", f"{workdir_prefix}/data/public", f"{workdir_prefix}/results", f"{workdir_prefix}/logs"]:
    try:
        file_system_client.create_directory(normalize_remote_path(remote_dir))
    except ResourceExistsError:
        pass

ensure_local_dirs(local_workdir)


def _download_dir(remote_dir: str, local_dir: Path) -> bool:
    remote_dir = normalize_remote_path(remote_dir)
    try:
        paths = list(file_system_client.get_paths(path=remote_dir, recursive=True))
    except ResourceNotFoundError:
        return False
    if not paths:
        return False
    local_dir.mkdir(parents=True, exist_ok=True)
    for path_info in paths:
        rel_path = Path(path_info.name).relative_to(remote_dir)
        local_path = local_dir / rel_path
        if path_info.is_directory:
            local_path.mkdir(parents=True, exist_ok=True)
            continue
        local_path.parent.mkdir(parents=True, exist_ok=True)
        file_bytes = file_system_client.get_file_client(path_info.name).download_file().readall()
        local_path.write_bytes(file_bytes)
    return True


def _upload_file(local_file: Path, remote_file: str) -> None:
    remote_file = normalize_remote_path(remote_file)
    file_client = file_system_client.get_file_client(remote_file)
    with local_file.open("rb") as handle:
        file_client.upload_data(handle, overwrite=True)


def _upload_dir(local_dir: Path, remote_dir: str) -> None:
    remote_dir = normalize_remote_path(remote_dir)
    try:
        file_system_client.create_directory(remote_dir)
    except ResourceExistsError:
        pass
    for child in local_dir.rglob("*"):
        rel_path = child.relative_to(local_dir)
        remote_path = f"{remote_dir}/{rel_path.as_posix()}"
        if child.is_dir():
            try:
                file_system_client.create_directory(remote_path)
            except ResourceExistsError:
                pass
        else:
            _upload_file(child, remote_path)


def _local_relpath(path_like) -> Path | None:
    path = Path(path_like)
    try:
        return path.relative_to(local_workdir)
    except ValueError:
        return None


def persistent_path_for(path_like) -> str | None:
    rel = _local_relpath(path_like)
    if rel is None:
        return None
    return f"{rtm_persistent_workdir}/{rel.as_posix()}"


def sync_local_to_abfss(path_like) -> str | None:
    src = Path(path_like)
    rel = _local_relpath(src)
    if rel is None or not src.exists():
        return None
    remote_path = f"{workdir_prefix}/{rel.as_posix()}"
    if src.is_dir():
        _upload_dir(src, remote_path)
    else:
        _upload_file(src, remote_path)
    return f"{source_uri}{remote_path}"


def restore_abfss_dir(path_like) -> Path | None:
    dst = Path(path_like)
    rel = _local_relpath(dst)
    if rel is None:
        return None
    restored = _download_dir(f"{workdir_prefix}/{rel.as_posix()}", dst)
    return dst if restored else None


def sync_path_to_persistent(path_like) -> str | None:
    return sync_local_to_abfss(path_like)


def restore_path_from_persistent(path_like) -> Path | None:
    return restore_abfss_dir(path_like)


for local_dir, label in [(local_results_dir, "results"), (local_public_data_dir, "public data")]:
    if any(local_dir.iterdir()):
        continue
    restored = restore_abfss_dir(local_dir)
    if restored is not None:
        print(f"Restored {label} from ABFSS: {persistent_path_for(local_dir)}")

os.environ["RTM_STORAGE_ACCOUNT"] = storage_account
os.environ["RTM_CONTAINER"] = container
os.environ["RTM_PERSISTENT_WORKDIR"] = rtm_persistent_workdir
os.environ["RTM_WORKDIR"] = rtm_persistent_workdir
os.environ["RTM_PUBLIC_DATA_DIR"] = rtm_public_data_dir
os.environ["RTM_RESULTS_DIR"] = rtm_results_dir
os.environ["RTM_LOCAL_WORKDIR"] = str(local_workdir)
os.environ["RTM_LOCAL_PUBLIC_DATA_DIR"] = str(local_public_data_dir)
os.environ["RTM_LOCAL_RESULTS_DIR"] = str(local_results_dir)
os.environ["RTM_KEEP_INTERMEDIATES"] = "1" if keep_intermediates else "0"

print(f"RTM_PERSISTENT_WORKDIR={rtm_persistent_workdir}")
print(f"Blob source={source_uri}")
print(f"RTM_WORKDIR={os.environ['RTM_WORKDIR']}")
print(f"RTM_PUBLIC_DATA_DIR={os.environ['RTM_PUBLIC_DATA_DIR']}")
print(f"RTM_RESULTS_DIR={os.environ['RTM_RESULTS_DIR']}")
print(f"RTM_LOCAL_WORKDIR={os.environ['RTM_LOCAL_WORKDIR']}")
print(f"RTM_LOCAL_PUBLIC_DATA_DIR={os.environ['RTM_LOCAL_PUBLIC_DATA_DIR']}")
print(f"RTM_LOCAL_RESULTS_DIR={os.environ['RTM_LOCAL_RESULTS_DIR']}")
print(f"RTM_KEEP_INTERMEDIATES={keep_intermediates}")


# COMMAND ----------

# DBTITLE 1,Download test and reference data
import json
import os
import queue
import shutil
import subprocess
import tempfile
import threading
from datetime import datetime
from pathlib import Path

import yaml

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or os.environ.get("RTM_REPO_DIR", "") or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
micromamba = os.path.expanduser("~/.local/bin/micromamba")
outdir = Path(os.environ.get("RTM_LOCAL_PUBLIC_DATA_DIR", os.environ["RTM_PUBLIC_DATA_DIR"]))
persistent_outdir = os.environ["RTM_PUBLIC_DATA_DIR"]
config_path = repo_dir / "resources" / "public_datasets.yaml"
manifest_path = outdir / "manifest.json"


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _format_monitored_paths(paths):
    parts = []
    for path in paths:
        if path.exists():
            parts.append(f"{path.name}={path.stat().st_size / (1024 * 1024):.1f}MB")
        else:
            parts.append(f"{path.name}=missing")
    return ", ".join(parts)


def run_and_stream(cmd, stage_name, cwd=repo_dir, extra_env=None, monitored_paths=None, heartbeat_sec=30):
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if extra_env:
        env.update(extra_env)

    print(f"[{ts()}] START {stage_name}")
    print(" ".join(str(x) for x in cmd))
    proc = subprocess.Popen(
        [str(x) for x in cmd],
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None

    output_queue = queue.Queue()

    def _reader():
        assert proc.stdout is not None
        for line in proc.stdout:
            output_queue.put(line)
        proc.stdout.close()

    reader_thread = threading.Thread(target=_reader, daemon=True)
    reader_thread.start()

    monitored_paths = monitored_paths or []
    while True:
        try:
            line = output_queue.get(timeout=heartbeat_sec)
            print(line, end="")
            continue
        except queue.Empty:
            if proc.poll() is not None and output_queue.empty():
                break
            heartbeat = f"[{ts()}] heartbeat: {stage_name} still running"
            if monitored_paths:
                heartbeat += " | " + _format_monitored_paths(monitored_paths)
            print(heartbeat)

    reader_thread.join(timeout=5)
    rc = proc.wait()
    print(f"[{ts()}] END {stage_name} rc={rc}")
    if rc != 0:
        raise RuntimeError(f"{stage_name} failed with exit code {rc}")


def sync_public_data_stage(stage_name, target_path=outdir):
    if "sync_path_to_persistent" in globals():
        synced = sync_path_to_persistent(target_path)
        if synced is not None:
            print(f"[{ts()}] Synced public data after {stage_name}: {synced}")


def run_seqc2_remote_slice_with_mates(dataset_name, url, out_bam, region="chr22", threads=4):
    out_bam = Path(out_bam)
    out_bai = Path(f"{out_bam}.bai")
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if out_bam.exists() and out_bai.exists():
        print(f"[{ts()}] SKIP {dataset_name}: existing outputs present")
        print(f"[{ts()}] {_format_monitored_paths([out_bam, out_bai])}")
        return

    tmpdir = outdir / "_tmp_seqc2" / dataset_name
    if tmpdir.exists():
        shutil.rmtree(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    region_bam = tmpdir / "region.bam"
    region_bai = Path(f"{region_bam}.bai")
    names_tsv = tmpdir / "mate_qnames.txt"
    mates_bam = tmpdir / "mates.bam"
    merged_bam = tmpdir / "merged.bam"

    run_and_stream(
        [
            micromamba, "run", "-n", "rtm-miner",
            "samtools", "view", "-@", str(threads), "-b", url, region, "-o", str(region_bam),
        ],
        f"{dataset_name}:slice_region",
        monitored_paths=[region_bam],
    )
    run_and_stream(
        [micromamba, "run", "-n", "rtm-miner", "samtools", "index", "-@", str(threads), str(region_bam)],
        f"{dataset_name}:index_region",
        monitored_paths=[region_bam, region_bai],
    )

    collect_script = f"""
from pathlib import Path
import pysam

region = {region!r}
region_bam = {str(region_bam)!r}
names_tsv = {str(names_tsv)!r}
region_chrom = region.split(':', 1)[0]
qnames = set()
with pysam.AlignmentFile(region_bam, 'rb') as bam:
    for read in bam.fetch(region=region):
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_qcfail or read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if read.next_reference_id < 0:
            continue
        mate_chrom = bam.get_reference_name(read.next_reference_id)
        if mate_chrom != region_chrom:
            qnames.add(read.query_name)
Path(names_tsv).write_text(('\\n'.join(sorted(qnames)) + ('\\n' if qnames else '')), encoding='utf-8')
print(f'mate_qnames={{len(qnames)}}')
"""
    run_and_stream(
        [micromamba, "run", "-n", "rtm-miner", "python", "-u", "-c", collect_script],
        f"{dataset_name}:collect_mate_qnames",
        monitored_paths=[names_tsv],
    )

    mate_qname_count = 0
    if names_tsv.exists():
        with names_tsv.open("r", encoding="utf-8") as handle:
            mate_qname_count = sum(1 for line in handle if line.strip())
    print(f"[{ts()}] {dataset_name}: mate_qnames={mate_qname_count}")

    if mate_qname_count:
        run_and_stream(
            [
                micromamba, "run", "-n", "rtm-miner",
                "samtools", "view", "-@", str(threads), "-b", "-N", str(names_tsv), url, "-o", str(mates_bam),
            ],
            f"{dataset_name}:slice_discordant_mates",
            monitored_paths=[mates_bam],
        )
        run_and_stream(
            [
                micromamba, "run", "-n", "rtm-miner",
                "samtools", "merge", "-@", str(threads), "-f", str(merged_bam), str(region_bam), str(mates_bam),
            ],
            f"{dataset_name}:merge_region_and_mates",
            monitored_paths=[region_bam, mates_bam, merged_bam],
        )
        shutil.copy2(merged_bam, out_bam)
    else:
        shutil.copy2(region_bam, out_bam)

    run_and_stream(
        [micromamba, "run", "-n", "rtm-miner", "samtools", "index", "-@", str(threads), str(out_bam)],
        f"{dataset_name}:index_final_bam",
        monitored_paths=[out_bam, out_bai],
    )
    print(f"[{ts()}] COMPLETED {dataset_name} | {_format_monitored_paths([out_bam, out_bai])}")


with config_path.open("r", encoding="utf-8") as handle:
    cfg = yaml.safe_load(handle)

dataset_by_id = {ds["id"]: ds for ds in cfg["datasets"]}

# Work around HG00100 CRAM slicing on Azure/Databricks by handling that dataset
# separately with an explicit local reference (`samtools view -T`). Also run the
# two SEQC2 BAM slices as explicit stages so progress is visible.
filtered_cfg = {
    **cfg,
    "datasets": [
        ds for ds in cfg["datasets"]
        if ds.get("id") not in {
            "hg00100_shortread_highcov_cram",
            "seqc2_disease_bam",
            "seqc2_control_bam",
        }
    ],
}

with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tmp:
    yaml.safe_dump(filtered_cfg, tmp, sort_keys=False)
    filtered_config_path = tmp.name

seqc2_disease_bam = outdir / "test_data" / "seqc2" / "chr22" / "disease.chr22.hg38.bam"
seqc2_control_bam = outdir / "test_data" / "seqc2" / "chr22" / "control.chr22.hg38.bam"

main_download_cmd = [
    micromamba,
    "run",
    "-n",
    "rtm-miner",
    "python",
    "-u",
    "scripts/download_public_data.py",
    "--config",
    filtered_config_path,
    "--outdir",
    str(outdir),
    "--references",
    "hg38",
    "--threads",
    "4",
    "--download-workers",
    "4",
]
run_and_stream(main_download_cmd, "download_public_data_non_bam")
sync_public_data_stage("download_public_data_non_bam")

run_seqc2_remote_slice_with_mates(
    "seqc2_disease_bam",
    dataset_by_id["seqc2_disease_bam"]["url"],
    seqc2_disease_bam,
)
sync_public_data_stage("seqc2_disease_bam", seqc2_disease_bam.parent)
run_seqc2_remote_slice_with_mates(
    "seqc2_control_bam",
    dataset_by_id["seqc2_control_bam"]["url"],
    seqc2_control_bam,
)
sync_public_data_stage("seqc2_control_bam", seqc2_control_bam.parent)

if manifest_path.exists():
    manifest = json.loads(manifest_path.read_text())
    print(f"[{ts()}] manifest summary: {json.dumps(manifest.get('summary', {}), indent=2)}")

reference_fasta = outdir / "reference" / "hg38" / "Homo_sapiens_assembly38.fasta"
if not reference_fasta.exists():
    raise FileNotFoundError(f"Expected local hg38 reference FASTA not found: {reference_fasta}")

hg00100_url = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram"
hg00100_bam = outdir / "test_data" / "1kg_hg00100" / "chr22" / "hg00100.shortread.chr22.hg38.bam"
hg00100_bam.parent.mkdir(parents=True, exist_ok=True)

slice_cmd = [
    micromamba,
    "run",
    "-n",
    "rtm-miner",
    "samtools",
    "view",
    "-@",
    "4",
    "-T",
    str(reference_fasta),
    "-b",
    hg00100_url,
    "chr22",
    "-o",
    str(hg00100_bam),
]
run_and_stream(slice_cmd, "slice_hg00100_chr22", monitored_paths=[hg00100_bam])

index_cmd = [
    micromamba,
    "run",
    "-n",
    "rtm-miner",
    "samtools",
    "index",
    "-@",
    "4",
    str(hg00100_bam),
]
run_and_stream(index_cmd, "index_hg00100_chr22_bam", monitored_paths=[hg00100_bam, hg00100_bam.with_suffix(hg00100_bam.suffix + ".bai")])
sync_public_data_stage("hg00100_chr22", hg00100_bam.parent)

print(f"[{ts()}] Downloaded manual HG00100 test BAM: {hg00100_bam}")
sync_public_data_stage("all_public_data", outdir)

# COMMAND ----------

# DBTITLE 1,Run chr22 no assembly pipeline
import os
import queue
import subprocess
import threading
from datetime import datetime
from pathlib import Path

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or os.environ.get("RTM_REPO_DIR", "") or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
micromamba = os.path.expanduser("~/.local/bin/micromamba")
public_dir = Path(os.environ.get("RTM_LOCAL_PUBLIC_DATA_DIR", os.environ["RTM_PUBLIC_DATA_DIR"]))
results_dir = Path(os.environ.get("RTM_LOCAL_RESULTS_DIR", os.environ["RTM_RESULTS_DIR"]))
persistent_results_dir = os.environ["RTM_RESULTS_DIR"]
run_outdir = results_dir / "chr22_noasm_seqc2"
keep_intermediates = os.environ.get("RTM_KEEP_INTERMEDIATES", "1") == "1"


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def format_paths(paths):
    parts = []
    for path in paths:
        if path.exists():
            parts.append(f"{path.name}={path.stat().st_size / (1024 * 1024):.1f}MB")
        else:
            parts.append(f"{path.name}=missing")
    return ", ".join(parts)


def run_and_stream(
    cmd,
    stage_name,
    cwd=repo_dir,
    extra_env=None,
    monitored_paths=None,
    heartbeat_sec=60,
    fail_fast_markers=("ERROR:", "Traceback (most recent call last)"),
):
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if extra_env:
        env.update(extra_env)

    print(f"[{ts()}] START {stage_name}")
    print(" ".join(str(x) for x in cmd))
    proc = subprocess.Popen(
        [str(x) for x in cmd],
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None

    output_queue = queue.Queue()

    def _reader():
        assert proc.stdout is not None
        for line in proc.stdout:
            output_queue.put(line)
        proc.stdout.close()

    reader_thread = threading.Thread(target=_reader, daemon=True)
    reader_thread.start()

    monitored_paths = monitored_paths or []
    fail_fast_line = None
    while True:
        try:
            line = output_queue.get(timeout=heartbeat_sec)
            print(line, end="")
            if any(marker in line for marker in fail_fast_markers):
                fail_fast_line = line.rstrip()
                proc.terminate()
                continue
            continue
        except queue.Empty:
            if proc.poll() is not None and output_queue.empty():
                break
            heartbeat = f"[{ts()}] heartbeat: {stage_name} still running"
            if monitored_paths:
                heartbeat += " | " + format_paths(monitored_paths)
            print(heartbeat)

    reader_thread.join(timeout=5)
    rc = proc.wait()
    print(f"[{ts()}] END {stage_name} rc={rc}")
    if fail_fast_line is not None:
        raise RuntimeError(f"{stage_name} aborted after error output: {fail_fast_line}")
    if rc != 0:
        raise RuntimeError(f"{stage_name} failed with exit code {rc}")


disease_bam = public_dir / "test_data" / "seqc2" / "chr22" / "disease.chr22.hg38.bam"
control_bam = public_dir / "test_data" / "seqc2" / "chr22" / "control.chr22.hg38.bam"
mei_fasta = public_dir / "retrotransposon_db" / "dfam" / "dfam_human_mei_l1_alu_sva.fasta"
review_tsv = run_outdir / "candidate_loci.mei.gold_review.tsv"
run_log = run_outdir / "pipeline.log"

for required_path in [disease_bam, control_bam, mei_fasta]:
    if not required_path.exists():
        raise FileNotFoundError(f"Required input not found: {required_path}")

run_outdir.mkdir(parents=True, exist_ok=True)
print(f"Keep intermediate files: {keep_intermediates}")

cmd = [
    micromamba,
    "run",
    "-n",
    "rtm-miner",
    "bash",
    "scripts/run_candidate_discovery_and_annotation.sh",
    "--reference-build",
    "hg38",
    "--disease-bam",
    str(disease_bam),
    "--control-bam",
    str(control_bam),
    "--disease-mate-bam",
    str(disease_bam),
    "--control-mate-bam",
    str(control_bam),
    "--mei-fasta",
    str(mei_fasta),
    "--outdir",
    str(run_outdir),
    "--chr",
    "chr22",
    "--no-local-assembly",
]

extra_env = {
    "RTM_WORKDIR": os.environ.get("RTM_LOCAL_WORKDIR", os.environ["RTM_WORKDIR"]),
    "RTM_PUBLIC_DATA_DIR": str(public_dir),
    "RTM_RESULTS_DIR": str(results_dir),
    "RTM_PERSISTENT_WORKDIR": os.environ["RTM_PERSISTENT_WORKDIR"],
}

run_and_stream(
    cmd,
    "run_chr22_no_assembly",
    extra_env=extra_env,
    monitored_paths=[review_tsv, run_log],
)

if "sync_path_to_persistent" in globals():
    sync_target = results_dir if keep_intermediates else run_outdir
    synced = sync_path_to_persistent(sync_target)
    if synced is not None:
        print(f"Synced run output to persistent storage: {synced}")
        print(f"Persistent results root: {persistent_results_dir}")

print(f"Output directory: {run_outdir}")
print(f"Review table: {review_tsv}")

# COMMAND ----------

# DBTITLE 1,Compare gold output to README
import os
from pathlib import Path
import pandas as pd

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or os.environ.get("RTM_REPO_DIR", "") or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
readme_path = repo_dir / "README.md"
gold_path = Path(os.environ.get("RTM_LOCAL_RESULTS_DIR", os.environ["RTM_RESULTS_DIR"])) / "chr22_noasm_seqc2" / "candidate_loci.mei.gold_review.tsv"
if not gold_path.exists() and "restore_path_from_persistent" in globals():
    restored = restore_path_from_persistent(gold_path)
    if restored is not None:
        print(f"Restored review table from persistent storage: {restored}")
if not gold_path.exists():
    raise FileNotFoundError(
        f"Gold review table not found at {gold_path}. "
        "Local machine storage does not survive cluster restart; set RTM_PERSISTENT_WORKDIR "
        "to an abfss:// path and re-run cells 3-5 to make results durable."
    )

readme_lines = readme_path.read_text(encoding="utf-8").splitlines()
start = None
for i, line in enumerate(readme_lines):
    if line.strip().startswith("| chrom | consensus_insertion_breakpoint_pos |"):
        start = i
        break
if start is None:
    raise RuntimeError("Could not find README example table header")

rows = []
for line in readme_lines[start + 2:]:
    if not line.strip().startswith("|"):
        break
    rows.append([p.strip() for p in line.strip().strip("|").split("|")])

header = [p.strip() for p in readme_lines[start].strip().strip("|").split("|")]
readme_all_df = pd.DataFrame(rows, columns=header).copy()
readme_all_df["consensus_insertion_breakpoint_pos"] = readme_all_df["consensus_insertion_breakpoint_pos"].astype(int)

previous_top_n = min(100, len(readme_all_df))
current_top_n = 25
readme_df = readme_all_df.head(previous_top_n).copy()

gold_df = pd.read_csv(gold_path, sep="\t", dtype=str, keep_default_na=False).head(current_top_n).copy()
gold_df["consensus_insertion_breakpoint_pos"] = gold_df["consensus_insertion_breakpoint_pos"].astype(int)

locus_cols = [
    "chrom",
    "consensus_insertion_breakpoint_pos",
]
annotation_cols = [
    "sample_status_label",
    "consensus_mei_family",
    "consensus_mei_subfamily",
    "known_mei_polymorphism_id",
    "known_mei_polymorphism_source",
]

comparison = gold_df[locus_cols + annotation_cols].merge(
    readme_df[locus_cols + annotation_cols],
    on=locus_cols,
    how="left",
    suffixes=("_current", "_readme"),
    indicator=True,
)

matched_current = comparison[comparison["_merge"] == "both"].copy()
current_not_in_previous = comparison[comparison["_merge"] == "left_only"].copy()

annotation_changed = matched_current[
    (matched_current["sample_status_label_current"] != matched_current["sample_status_label_readme"])
    | (matched_current["consensus_mei_family_current"] != matched_current["consensus_mei_family_readme"])
    | (matched_current["consensus_mei_subfamily_current"] != matched_current["consensus_mei_subfamily_readme"])
    | (matched_current["known_mei_polymorphism_id_current"] != matched_current["known_mei_polymorphism_id_readme"])
    | (matched_current["known_mei_polymorphism_source_current"] != matched_current["known_mei_polymorphism_source_readme"])
].copy()

summary_df = pd.DataFrame([
    {
        "previous_top_n_checked": len(readme_df),
        "current_top_n_checked": len(gold_df),
        "current_loci_found_in_previous_top_n": len(matched_current),
        "current_loci_not_found_in_previous_top_n": len(current_not_in_previous),
        "fraction_current_loci_found": round(len(matched_current) / len(gold_df), 3),
        "matched_loci_with_annotation_changes": len(annotation_changed),
    }
])

display(summary_df)

matched_view = matched_current[locus_cols + [
    "sample_status_label_current",
    "sample_status_label_readme",
    "consensus_mei_subfamily_current",
    "consensus_mei_subfamily_readme",
]].copy()
matched_view.insert(0, "current_rank", matched_current.index + 1)
display(matched_view)

if len(annotation_changed):
    annotation_view = annotation_changed[locus_cols + [
        "sample_status_label_current",
        "sample_status_label_readme",
        "consensus_mei_family_current",
        "consensus_mei_family_readme",
        "consensus_mei_subfamily_current",
        "consensus_mei_subfamily_readme",
        "known_mei_polymorphism_id_current",
        "known_mei_polymorphism_id_readme",
        "known_mei_polymorphism_source_current",
        "known_mei_polymorphism_source_readme",
    ]].copy()
    annotation_view.insert(0, "current_rank", annotation_changed.index + 1)
    display(annotation_view)

if len(current_not_in_previous):
    missing_view = current_not_in_previous[locus_cols + [
        "sample_status_label_current",
        "consensus_mei_family_current",
        "consensus_mei_subfamily_current",
        "known_mei_polymorphism_id_current",
        "known_mei_polymorphism_source_current",
    ]].copy()
    missing_view.insert(0, "current_rank", current_not_in_previous.index + 1)
    display(missing_view)

    # --- Investigate WHY each absent locus is missing ---
    # Check if the same event exists in README at a slightly shifted coordinate (±500 bp)
    # and whether it exists anywhere in the full gold output (not just top 25).
    gold_full = pd.read_csv(gold_path, sep="\t", dtype=str, keep_default_na=False)
    gold_full["consensus_insertion_breakpoint_pos"] = gold_full["consensus_insertion_breakpoint_pos"].astype(int)
    gold_full["gold_rank"] = range(1, len(gold_full) + 1)

    proximity_bp = 500
    proximity_rows = []
    for _, row in current_not_in_previous.iterrows():
        chrom = row["chrom"]
        pos = row["consensus_insertion_breakpoint_pos"]
        current_rank = row.name + 1
        # Nearest README locus on same chrom
        readme_same_chrom = readme_all_df[readme_all_df["chrom"] == chrom].copy()
        readme_same_chrom["dist"] = (readme_same_chrom["consensus_insertion_breakpoint_pos"] - pos).abs()
        nearest = readme_same_chrom.nsmallest(1, "dist")
        if not nearest.empty:
            n = nearest.iloc[0]
            dist = int(n["dist"])
            near_readme_pos = int(n["consensus_insertion_breakpoint_pos"])
            near_readme_family = n["consensus_mei_family"]
            readme_rank = int(nearest.index[0]) + 1
        else:
            dist = None
            near_readme_pos = None
            near_readme_family = None
            readme_rank = None
        # Check if present in full gold (beyond top 25)
        gold_match = gold_full[
            (gold_full["chrom"] == chrom)
            & ((gold_full["consensus_insertion_breakpoint_pos"] - pos).abs() <= proximity_bp)
        ]
        gold_rank_str = ",".join(str(r) for r in gold_match["gold_rank"].tolist()) if not gold_match.empty else "absent"
        proximity_rows.append({
            "current_rank": current_rank,
            "chrom": chrom,
            "current_pos": pos,
            "nearest_readme_pos": near_readme_pos,
            "nearest_readme_dist_bp": dist,
            "coord_shifted": dist is not None and 0 < dist <= proximity_bp,
            "nearest_readme_family": near_readme_family,
            "readme_rank": readme_rank,
            "gold_rank_within_500bp": gold_rank_str,
        })

    proximity_df = pd.DataFrame(proximity_rows)
    display(proximity_df)

# COMMAND ----------

# DBTITLE 1,Show IGV plots for non-overlapping loci
import math
import os
import subprocess
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

notebook_path = dbutils.notebook.entry_point.getDbutils().notebook().getContext().notebookPath().get()
repo_dir_default = Path("/Workspace") / Path(notebook_path.lstrip("/")).parent
repo_dir = Path(dbutils.widgets.get("rtm_repo_dir").strip() or os.environ.get("RTM_REPO_DIR", "") or str(repo_dir_default))
os.environ["RTM_REPO_DIR"] = str(repo_dir)
micromamba = os.path.expanduser("~/.local/bin/micromamba")
public_dir = Path(os.environ.get("RTM_LOCAL_PUBLIC_DATA_DIR", os.environ["RTM_PUBLIC_DATA_DIR"]))
results_dir = Path(os.environ.get("RTM_LOCAL_RESULTS_DIR", os.environ["RTM_RESULTS_DIR"]))
run_outdir = results_dir / "chr22_noasm_seqc2"
snapshot_dir = run_outdir / "igv_non_overlap"
snapshot_dir.mkdir(parents=True, exist_ok=True)
index_path = snapshot_dir / "igv_snapshot_index.tsv"
variant_tsv = snapshot_dir / "non_overlap_variants.tsv"

if "current_not_in_previous" not in globals() or "gold_df" not in globals():
    raise RuntimeError("Run the comparison cell first so current_not_in_previous and gold_df are available.")

non_overlap_loci = current_not_in_previous[["chrom", "consensus_insertion_breakpoint_pos"]].drop_duplicates().copy()
if non_overlap_loci.empty:
    print("No non-overlapping loci found in the current top-ranked calls.")
else:
    igv_variants = gold_df.merge(
        non_overlap_loci,
        on=["chrom", "consensus_insertion_breakpoint_pos"],
        how="inner",
    ).copy()
    igv_variants.insert(0, "current_rank", range(1, len(igv_variants) + 1))
    igv_variants.to_csv(variant_tsv, sep="\t", index=False)

    disease_bam = public_dir / "test_data" / "seqc2" / "chr22" / "disease.chr22.hg38.bam"
    control_bam = public_dir / "test_data" / "seqc2" / "chr22" / "control.chr22.hg38.bam"
    reference_fasta = public_dir / "reference" / "hg38" / "Homo_sapiens_assembly38.fasta"
    required_paths = [
        reference_fasta,
        disease_bam,
        Path(f"{disease_bam}.bai"),
        control_bam,
        Path(f"{control_bam}.bai"),
    ]
    for path in required_paths:
        if not path.exists() and "restore_path_from_persistent" in globals():
            restore_path_from_persistent(path)
    if not index_path.exists() and "restore_path_from_persistent" in globals():
        restore_path_from_persistent(snapshot_dir)

    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required local files for IGV plots: " + ", ".join(missing))

    regenerate = True
    if index_path.exists():
        existing_index = pd.read_csv(index_path, sep="\t", dtype=str, keep_default_na=False)
        existing_pngs = [Path(path) for path in existing_index.get("snapshot_png", pd.Series(dtype=str)).tolist()]
        if len(existing_index) == len(igv_variants) and existing_pngs and all(path.exists() for path in existing_pngs):
            regenerate = False
            print(f"Using existing IGV snapshots from {snapshot_dir}")

    if regenerate:
        assembly_cache_dir = run_outdir / "assembly_cache"
        render_script = f"""
from pathlib import Path
import pandas as pd
from retro_miner.igv_plots import generate_gold_review_igv_plots

variants = pd.read_csv({str(variant_tsv)!r}, sep='\t', dtype=str, keep_default_na=False)
for column in [
    'current_rank',
    'consensus_insertion_breakpoint_pos',
    'window_start',
    'window_end',
    'discovery_window_start',
    'discovery_window_end',
    'insertion_breakpoint_pos',
]:
    if column in variants.columns:
        variants[column] = pd.to_numeric(variants[column], errors='coerce').fillna(0).astype(int)

index_path = generate_gold_review_igv_plots(
    variants,
    reference_fasta=Path({str(reference_fasta)!r}),
    disease_bam=Path({str(disease_bam)!r}),
    control_bam=Path({str(control_bam)!r}),
    snapshot_dir=Path({str(snapshot_dir)!r}),
    top_n=0,
    gold_only=False,
    assembly_cache_dir=Path({str(assembly_cache_dir)!r}),
)
print(index_path if index_path is not None else '')
"""
        result = subprocess.run(
            [micromamba, "run", "-n", "rtm-miner", "python", "-u", "-c", render_script],
            cwd=repo_dir,
            env={**os.environ, "PYTHONPATH": "src"},
            text=True,
            capture_output=True,
            check=False,
        )
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        if result.returncode != 0:
            raise RuntimeError(f"IGV plot generation failed with exit code {result.returncode}")
        if "sync_path_to_persistent" in globals():
            synced = sync_path_to_persistent(snapshot_dir)
            if synced is not None:
                print(f"Synced IGV snapshots to persistent storage: {synced}")

    index_df = pd.read_csv(index_path, sep="\t", dtype=str, keep_default_na=False)
    display(index_df[["plot_rank", "chrom", "discovery_window_start", "discovery_window_end", "snapshot_png"]])

    png_paths = [Path(path) for path in index_df["snapshot_png"].tolist() if Path(path).exists()]
    if not png_paths:
        raise FileNotFoundError(f"No IGV snapshot PNGs found in {snapshot_dir}")

    for png_path in png_paths:
        row = index_df.loc[index_df["snapshot_png"] == str(png_path)].iloc[0]
        fig, ax = plt.subplots(figsize=(18, 12))
        ax.imshow(mpimg.imread(png_path))
        ax.set_title(
            f"rank {row['plot_rank']} | {row['chrom']}:{row['discovery_window_start']}-{row['discovery_window_end']}",
            fontsize=16,
        )
        ax.axis("off")
        plt.tight_layout()
        display(fig)
        plt.close(fig)

