#!/usr/bin/env python3
"""Copy one mouse WGS run from the public SRA bucket and align it to mm39.

The FASTQs on ENA are the same reads as ``s3://sra-pub-run-odp``. The dump
uses the ENA defline, so the BAM QNAME, sequence, and quality match an
alignment of those FASTQs. Alignment is ``bwa-mem3`` against the mm39-mem3
index. Read order in the sorted BAM still follows whichever order that
emitted.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from retro_miner.mouse_public_data import (  # noqa: E402
    BWA_MEM3_BIN,
    BWA_MEM3_INDEX_SUFFIXES,
    NCBI_SRA_BUCKET,
    bwa_mem3_env,
    bwa_mem_command,
    fasterq_dump_command,
    ncbi_sra_s3_uri,
    s3_cp_command,
    s3_single_file_config_text,
    scratch_bytes_for_sra,
    validate_run_accession,
)

# A finished 30x mouse BAM is well above this. Short objects are failed uploads.
MIN_BAM_BYTES = 5 * 1024**3


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    print(f"{stamp} {message}", flush=True)


def load_mouse_libraries(config_path: Path) -> list[dict[str, str]]:
    with config_path.open(encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    rows = cfg.get("mouse_wgs_libraries") or []
    if not isinstance(rows, list):
        raise SystemExit(f"mouse_wgs_libraries is not a list in {config_path}")
    libraries = []
    for row in rows:
        libraries.append(
            {
                "id": str(row["id"]),
                "accession": validate_run_accession(str(row["accession"])),
                "paper": str(row["paper"]).strip(),
            }
        )
    return libraries


def write_s3_config(path: Path) -> dict[str, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s3_single_file_config_text(), encoding="utf-8")
    env = os.environ.copy()
    env["AWS_CONFIG_FILE"] = str(path)
    env["AWS_DEFAULT_REGION"] = "us-east-1"
    env["AWS_PAGER"] = ""
    return env


def run_cmd(cmd: list[str], env: dict[str, str] | None = None) -> None:
    log("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True, env=env)


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def s3_uri_parts(uri: str) -> tuple[str, str]:
    if not uri.startswith("s3://"):
        raise SystemExit(f"Not an S3 URI: {uri}")
    bucket, _, key = uri[5:].partition("/")
    if not bucket or not key:
        raise SystemExit(f"Not an S3 URI: {uri}")
    return bucket, key


def s3_size(uri: str, env: dict[str, str], *, anonymous: bool) -> int | None:
    bucket, key = s3_uri_parts(uri)
    cmd = [
        "aws",
        "s3api",
        "head-object",
        "--bucket",
        bucket,
        "--key",
        key,
        "--query",
        "ContentLength",
        "--output",
        "text",
    ]
    if anonymous:
        cmd.append("--no-sign-request")
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        return None
    text = proc.stdout.strip()
    if not text.isdigit():
        return None
    return int(text)


def wait_for_copy_pid(pidfile: Path, dest: Path) -> None:
    if not pidfile.exists():
        return
    try:
        pid = int(pidfile.read_text(encoding="utf-8").strip())
    except ValueError:
        return
    while pid_alive(pid):
        have = dest.stat().st_size if dest.exists() else 0
        log(f"s3 cp already running pid={pid} bytes={have}")
        time.sleep(15)


def ensure_local_sra(
    accession: str,
    dest: Path,
    env: dict[str, str],
    pidfile: Path,
) -> int:
    src = ncbi_sra_s3_uri(accession)
    expected = s3_size(src, env, anonymous=True)
    if not expected:
        raise SystemExit(f"Public SRA object is missing: {src}")
    wait_for_copy_pid(pidfile, dest)
    if dest.exists() and dest.stat().st_size == expected:
        log(f"SRA already complete {dest} bytes={expected}")
        return expected
    if dest.exists():
        log(f"removing short SRA {dest} bytes={dest.stat().st_size}")
        dest.unlink()
    dest.parent.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(dest.parent).free
    need = scratch_bytes_for_sra(expected)
    if free < need:
        raise SystemExit(
            f"Need {need} free bytes for {accession} and have {free} on {dest.parent}"
        )
    cmd = s3_cp_command(src, str(dest), anonymous=True)
    started = time.time()
    for attempt in range(1, 4):
        try:
            run_cmd(cmd, env)
        except subprocess.CalledProcessError:
            log(f"s3 cp attempt {attempt} failed for {accession}")
            if attempt == 3:
                raise
            time.sleep(5)
            continue
        break
    got = dest.stat().st_size if dest.exists() else 0
    if got != expected:
        raise SystemExit(f"SRA copy is {got} bytes, expected {expected}")
    elapsed = max(time.time() - started, 0.001)
    log(f"copied {accession} in {elapsed:.1f}s ({got / elapsed / 1e6:.0f} MB/s)")
    return expected


def dump_fastqs(accession: str, sra_path: Path, outdir: Path, threads: int) -> tuple[Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    tempdir = outdir.parent / "tmp"
    if tempdir.exists():
        shutil.rmtree(tempdir)
    tempdir.mkdir(parents=True)
    read1 = outdir / f"{accession}_1.fastq"
    read2 = outdir / f"{accession}_2.fastq"
    cmd = fasterq_dump_command(str(sra_path), str(outdir), str(tempdir), threads)
    run_cmd(cmd)
    if not read1.exists() or not read2.exists():
        raise SystemExit(f"fasterq-dump did not write {read1.name} and {read2.name}")
    with read1.open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n")
    log(f"first read header: {header}")
    if not header.startswith(f"@{accession}."):
        raise SystemExit(f"Unexpected FASTQ header for {accession}: {header}")
    shutil.rmtree(tempdir, ignore_errors=True)
    return read1, read2


def align_pair(
    reference: Path,
    read1: Path,
    read2: Path,
    bam_path: Path,
    threads: int,
    sample: str,
) -> None:
    if shutil.which(BWA_MEM3_BIN) is None:
        raise SystemExit(f"{BWA_MEM3_BIN} is not on PATH")
    for suffix in BWA_MEM3_INDEX_SUFFIXES:
        index = Path(str(reference) + suffix)
        if not index.is_file() or index.stat().st_size == 0:
            raise SystemExit(f"bwa-mem3 index file is missing: {index}")
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    bwa_cmd = bwa_mem_command(str(reference), str(read1), str(read2), threads, sample)
    sort_threads = max(2, min(4, threads))
    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(sort_threads),
        "-m",
        "2G",
        "-o",
        str(bam_path),
        "-",
    ]
    log("+ " + " ".join(bwa_cmd) + " | " + " ".join(sort_cmd))
    bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, env=bwa_mem3_env())
    assert bwa.stdout is not None
    sort = subprocess.Popen(sort_cmd, stdin=bwa.stdout)
    bwa.stdout.close()
    sort_rc = sort.wait()
    bwa_rc = bwa.wait()
    if bwa_rc != 0 or sort_rc != 0:
        bam_path.unlink(missing_ok=True)
        raise SystemExit(f"alignment failed bwa={bwa_rc} samtools={sort_rc}")
    run_cmd(["samtools", "index", "-@", str(sort_threads), str(bam_path)])
    run_cmd(["samtools", "quickcheck", "-v", str(bam_path)])


def upload_bam(bam_path: Path, dest_uri: str, env: dict[str, str]) -> None:
    bai_path = Path(str(bam_path) + ".bai")
    run_cmd(s3_cp_command(str(bam_path), dest_uri, anonymous=False), env)
    run_cmd(s3_cp_command(str(bai_path), dest_uri + ".bai", anonymous=False), env)
    if (s3_size(dest_uri, env, anonymous=False) or 0) < MIN_BAM_BYTES:
        raise SystemExit(f"Uploaded BAM is too small: {dest_uri}")


def bam_done(dest_uri: str, env: dict[str, str]) -> bool:
    bam_size = s3_size(dest_uri, env, anonymous=False) or 0
    bai_size = s3_size(dest_uri + ".bai", env, anonymous=False) or 0
    return bam_size >= MIN_BAM_BYTES and bai_size > 0


def align_library(
    library: dict[str, str],
    reference: Path,
    workdir: Path,
    s3_prefix: str,
    threads: int,
    env: dict[str, str],
    pidfile: Path,
) -> None:
    accession = library["accession"]
    dest_uri = f"{s3_prefix.rstrip('/')}/{library['paper']}/{accession}.bam"
    if bam_done(dest_uri, env):
        log(f"skip {accession}; {dest_uri} is already uploaded")
        return
    library_dir = workdir / accession
    library_dir.mkdir(parents=True, exist_ok=True)
    sra_path = library_dir / accession
    fastq_dir = library_dir / "fastq"
    read1 = fastq_dir / f"{accession}_1.fastq"
    read2 = fastq_dir / f"{accession}_2.fastq"
    bam_path = library_dir / f"{accession}.mm39.bam"
    stamp = fastq_dir / "dump.complete"
    if stamp.exists() and read1.exists() and read2.exists():
        log(f"FASTQs already present for {accession}")
    else:
        ensure_local_sra(accession, sra_path, env, pidfile)
        read1, read2 = dump_fastqs(accession, sra_path, fastq_dir, threads)
        stamp.write_text(accession + "\n", encoding="utf-8")
        if sra_path.exists():
            sra_path.unlink()
            log(f"removed {sra_path}")
    if not bam_path.exists():
        align_pair(reference, read1, read2, bam_path, threads, library["id"])
    upload_bam(bam_path, dest_uri, env)
    for path in (read1, read2, bam_path, Path(str(bam_path) + ".bai")):
        path.unlink(missing_ok=True)
    log(f"finished {accession} -> {dest_uri}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--s3-prefix", required=True)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--accession", action="append", default=[])
    parser.add_argument(
        "--copy-pidfile",
        type=Path,
        default=None,
        help="Pid of an s3 cp started ahead of this process. Wait for it before copying.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    index_marker = Path(str(args.reference) + ".bwt.2bit.64")
    if not args.reference.exists() and not index_marker.is_file():
        raise SystemExit(f"bwa-mem3 index prefix is missing: {args.reference}")
    libraries = load_mouse_libraries(args.config)
    wanted = {validate_run_accession(acc) for acc in args.accession}
    if wanted:
        libraries = [row for row in libraries if row["accession"] in wanted]
        missing = wanted - {row["accession"] for row in libraries}
        if missing:
            raise SystemExit(f"Accessions not in the mouse catalog: {', '.join(sorted(missing))}")
    if not libraries:
        raise SystemExit("No mouse libraries selected")
    args.workdir.mkdir(parents=True, exist_ok=True)
    env = write_s3_config(args.workdir / "s3-crt.config")
    env["TMPDIR"] = str(args.workdir)
    pidfile = args.copy_pidfile or (args.workdir / "s3-cp.pid")
    log(
        f"libraries={len(libraries)} threads={args.threads} "
        f"bucket={NCBI_SRA_BUCKET} prefix={args.s3_prefix}"
    )
    for library in libraries:
        log(f"start {library['id']} {library['accession']}")
        align_library(
            library,
            args.reference,
            args.workdir,
            args.s3_prefix,
            args.threads,
            env,
            pidfile,
        )
    log("all selected libraries aligned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
