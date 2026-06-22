#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Dataset:
    dataset_id: str
    category: str
    description: str
    source: str
    url: str
    target_path: str
    region: str | None = None
    required: bool = True


BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
DFAM_CURATED_CONSENSUS_0_URL = "https://www.dfam.org/releases/current/families/FamDB/dfam40.curated.consensus.0.h5.gz"


def _load_config(config_path: Path) -> list[Dataset]:
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    if not isinstance(cfg, dict) or "datasets" not in cfg:
        raise ValueError(f"Invalid config format: {config_path}")

    required_fields = ("id", "category", "url", "target_path")
    items = []
    for idx, row in enumerate(cfg["datasets"]):
        if not isinstance(row, dict):
            raise ValueError(f"Invalid dataset entry at index {idx}: expected mapping, got {type(row).__name__}")
        missing = [k for k in required_fields if k not in row]
        if missing:
            ds_id = row.get("id", f"index_{idx}")
            raise ValueError(
                f"Dataset '{ds_id}' missing required fields: {', '.join(missing)} in {config_path}"
            )

        items.append(
            Dataset(
                dataset_id=row["id"],
                category=row["category"],
                description=row.get("description", ""),
                source=row.get("source", "unknown"),
                url=row["url"],
                target_path=row["target_path"],
                region=row.get("region"),
                required=bool(row.get("required", True)),
            )
        )
    return items


def _download_file(url: str, out_path: Path, timeout_sec: int, force: bool) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        return {
            "status": "skipped_exists",
            "bytes": out_path.stat().st_size,
            "path": str(out_path),
            "url": url,
        }

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "retrotransposon-miner/0.1 (+https://github.com/L1-Therapeutics/retrotransposon-miner)"
        },
    )

    started = time.time()
    tmp_path = out_path.with_name(f"{out_path.name}.part")
    if tmp_path.exists():
        try:
            tmp_path.unlink()
        except OSError:
            pass
    try:
        with urllib.request.urlopen(request, timeout=timeout_sec) as resp, tmp_path.open("wb") as out_handle:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                out_handle.write(chunk)
        tmp_path.replace(out_path)
    except Exception:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass
        raise
    elapsed = time.time() - started
    size = out_path.stat().st_size
    return {
        "status": "downloaded",
        "bytes": size,
        "seconds": round(elapsed, 3),
        "path": str(out_path),
        "url": url,
    }


def _slice_remote_bam(url: str, region: str, out_bam: Path, threads: int, force: bool) -> dict[str, Any]:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if out_bam.exists() and Path(f"{out_bam}.bai").exists() and not force:
        return {"status": "skipped_exists", "path": str(out_bam), "bytes": out_bam.stat().st_size, "region": region}

    # Direct remote slicing avoids storing full-size BAM locally.
    _run_cmd(["samtools", "view", "-@", str(threads), "-b", url, region, "-o", str(out_bam)], required=True)
    _run_cmd(["samtools", "index", "-@", str(threads), str(out_bam)], required=True)
    return {"status": "sliced_remote_bam", "path": str(out_bam), "bytes": out_bam.stat().st_size, "region": region}


def _run_cmd(cmd: list[str], required: bool = True) -> tuple[bool, str]:
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True, proc.stdout.strip()
    except FileNotFoundError as err:
        msg = f"missing executable: {cmd[0]} ({err})"
        return (False, msg) if not required else (_raise_runtime(msg))
    except subprocess.CalledProcessError as err:
        msg = f"command failed: {' '.join(cmd)}\n{err.stderr}"
        return (False, msg) if not required else (_raise_runtime(msg))


def _ensure_postprocess_binaries() -> None:
    required_bins = ("bigWigToBedGraph", "bigBedToBed", "liftOver")
    missing = [name for name in required_bins if shutil.which(name) is None]
    if not missing:
        return
    raise RuntimeError(
        "Missing required postprocess executables: "
        f"{', '.join(missing)}. "
        "Install environment dependencies (including UCSC tools) and re-run, e.g. "
        "'bash scripts/bootstrap_env.sh' followed by 'bash scripts/validate_environment.sh'."
    )


def _raise_runtime(msg: str) -> tuple[bool, str]:
    raise RuntimeError(msg)


def _decompress_gzip(src: Path, dst: Path, force: bool) -> dict[str, Any]:
    if dst.exists() and not force:
        return {"status": "skipped_exists", "path": str(dst), "bytes": dst.stat().st_size}

    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as in_handle, dst.open("wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)
    return {"status": "decompressed", "path": str(dst), "bytes": dst.stat().st_size}


def _has_bwa_index(fasta_path: Path) -> bool:
    std = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if all(Path(f"{fasta_path}{s}").exists() for s in std):
        return True
    alt64 = [".64.amb", ".64.ann", ".64.bwt", ".64.pac", ".64.sa"]
    if all(Path(f"{fasta_path}{s}").exists() for s in alt64):
        return True
    return False


def _prebuilt_bwa_urls_from_fasta_url(fasta_url: str) -> dict[str, str]:
    parsed = urllib.parse.urlsplit(fasta_url)
    basename = os.path.basename(parsed.path)
    if basename.endswith(".gz"):
        return {}
    return {suffix: f"{fasta_url}{suffix}" for suffix in BWA_INDEX_SUFFIXES}


def _try_download_prebuilt_bwa_indexes(
    fasta_path: Path,
    fasta_url: str,
    ref_label: str,
    timeout_sec: int,
    force: bool,
) -> list[dict[str, Any]]:
    actions: list[dict[str, Any]] = []

    if _has_bwa_index(fasta_path) and not force:
        return [{"step": f"prebuilt_bwa_{ref_label}", "status": "skipped_exists", "path": str(fasta_path)}]

    source = _prebuilt_bwa_urls_from_fasta_url(fasta_url)
    if not source:
        actions.append(
            {
                "step": f"prebuilt_bwa_{ref_label}",
                "status": "skipped_unsupported_fasta_url",
                "path": str(fasta_path),
                "url": fasta_url,
            }
        )
        return actions

    print(f"[postprocess] attempting prebuilt {ref_label} BWA index download", file=sys.stderr)
    ok_all = True
    for suffix, url in source.items():
        out_path = Path(f"{fasta_path}{suffix}")
        try:
            res = _download_file(url, out_path, timeout_sec=timeout_sec, force=force)
            actions.append({"step": f"prebuilt_bwa_{ref_label}_file", "file": out_path.name, **res})
        except Exception as err:  # noqa: BLE001
            ok_all = False
            actions.append(
                {
                    "step": f"prebuilt_bwa_{ref_label}_file",
                    "file": out_path.name,
                    "status": "failed",
                    "error": str(err),
                    "url": url,
                }
            )
            break

    if ok_all and _has_bwa_index(fasta_path):
        actions.append({"step": f"prebuilt_bwa_{ref_label}", "status": "downloaded", "path": str(fasta_path)})
    else:
        actions.append({"step": f"prebuilt_bwa_{ref_label}", "status": "failed", "path": str(fasta_path)})
    return actions


def _build_fasta_indexes(
    fasta_path: Path, fasta_url: str, ref_label: str, force: bool, timeout_sec: int, prefer_prebuilt_bwa_index: bool
) -> list[dict[str, Any]]:
    actions: list[dict[str, Any]] = []

    fai_path = Path(f"{fasta_path}.fai")
    dict_path = fasta_path.with_suffix(".dict")

    if not fai_path.exists() or force:
        ok, msg = _run_cmd(["samtools", "faidx", str(fasta_path)], required=True)
        actions.append({"step": "samtools_faidx", "ok": ok, "message": msg, "path": str(fai_path)})
    else:
        actions.append({"step": "samtools_faidx", "ok": True, "status": "skipped_exists", "path": str(fai_path)})

    if not dict_path.exists() or force:
        ok, msg = _run_cmd(["samtools", "dict", "-o", str(dict_path), str(fasta_path)], required=True)
        actions.append({"step": "samtools_dict", "ok": ok, "message": msg, "path": str(dict_path)})
    else:
        actions.append({"step": "samtools_dict", "ok": True, "status": "skipped_exists", "path": str(dict_path)})

    # BWA index set for WGS remapping.
    # 1) Prefer prebuilt BWA index download when available from the same FASTA URL location.
    if prefer_prebuilt_bwa_index:
        for pre in _try_download_prebuilt_bwa_indexes(
            fasta_path, fasta_url=fasta_url, ref_label=ref_label, timeout_sec=timeout_sec, force=force
        ):
            actions.append(pre)
        if _has_bwa_index(fasta_path):
            actions.append({"step": "bwa_index", "ok": True, "status": "skipped_prebuilt", "path": str(fasta_path)})
            return actions

    # 2) Otherwise build locally (classic bwa first).
    if not _has_bwa_index(fasta_path) or force:
        bwa_msg = ""
        mem2_msg = ""
        if shutil.which("bwa"):
            ok, msg = _run_cmd(["bwa", "index", "-a", "bwtsw", str(fasta_path)], required=False)
            bwa_msg = msg
            actions.append(
                {
                    "step": "bwa_index",
                    "ok": ok,
                    "status": "indexed" if ok else "failed",
                    "message": msg,
                    "path": str(fasta_path),
                }
            )
            if ok:
                return actions

        if shutil.which("bwa-mem2"):
            ok, msg = _run_cmd(["bwa-mem2", "index", str(fasta_path)], required=False)
            mem2_msg = msg
            actions.append(
                {
                    "step": "bwa_mem2_index",
                    "ok": ok,
                    "status": "indexed" if ok else "failed",
                    "message": msg,
                    "path": str(fasta_path),
                }
            )
            if ok:
                return actions

        if not shutil.which("bwa") and not shutil.which("bwa-mem2"):
            raise RuntimeError(
                "Unable to build BWA indexes: neither 'bwa' nor 'bwa-mem2' was found in PATH. "
                "Install one of them or enable prebuilt index download."
            )
        raise RuntimeError("Unable to build BWA indexes with either bwa or bwa-mem2.\n" f"bwa: {bwa_msg}\n" f"bwa-mem2: {mem2_msg}")
    else:
        actions.append({"step": "bwa_index", "ok": True, "status": "skipped_exists", "path": str(fasta_path)})

    return actions


def _convert_ucsc_txt_to_bed(src_gz: Path, dst_bed: Path, force: bool) -> dict[str, Any]:
    if dst_bed.exists() and not force:
        return {"status": "skipped_exists", "path": str(dst_bed), "bytes": dst_bed.stat().st_size}

    dst_bed.parent.mkdir(parents=True, exist_ok=True)
    lines = 0
    with gzip.open(src_gz, "rt", encoding="utf-8", newline="") as in_handle, dst_bed.open(
        "w", encoding="utf-8"
    ) as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            # UCSC table dumps vary by schema:
            # - plain BED-like rows: chrom/start/end at columns 0/1/2
            # - many db-table dumps: `bin` + score fields before chrom/start/end
            # Locate the first `chr*` token and read the next two numeric columns.
            chrom = ""
            start = ""
            end = ""
            extra = "."
            for i, token in enumerate(row):
                if not token.startswith("chr"):
                    continue
                if i + 2 >= len(row):
                    continue
                try:
                    start_i = int(row[i + 1])
                    end_i = int(row[i + 2])
                except ValueError:
                    continue
                if end_i <= start_i:
                    continue
                chrom = token
                start = str(start_i)
                end = str(end_i)
                if i + 3 < len(row):
                    extra = row[i + 3] if row[i + 3] else "."
                break
            if not chrom:
                continue
            out_handle.write(f"{chrom}\t{start}\t{end}\t{extra}\n")
            lines += 1

    return {"status": "converted", "path": str(dst_bed), "rows": lines}


def _filter_mei_bed(src_bed: Path, dst_bed: Path, force: bool) -> dict[str, Any]:
    if dst_bed.exists() and not force:
        return {"status": "skipped_exists", "path": str(dst_bed), "bytes": dst_bed.stat().st_size}

    dst_bed.parent.mkdir(parents=True, exist_ok=True)
    patterns = ("mobile element insertion", "line1", "alu", "sva", "ins:me")
    kept = 0
    with src_bed.open("r", encoding="utf-8") as in_handle, dst_bed.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            low = line.lower()
            if any(p in low for p in patterns):
                parts = line.rstrip("\n").split("\t")
                chrom = ""
                start_i = -1
                end_i = -1
                extra = "."
                for i, token in enumerate(parts):
                    if not token.startswith("chr"):
                        continue
                    if i + 2 >= len(parts):
                        continue
                    try:
                        s = int(parts[i + 1])
                        e = int(parts[i + 2])
                    except ValueError:
                        continue
                    if e <= s:
                        continue
                    chrom = token
                    start_i = s
                    end_i = e
                    if i + 3 < len(parts) and parts[i + 3]:
                        extra = parts[i + 3]
                    break
                if not chrom:
                    continue
                out_handle.write(f"{chrom}\t{start_i}\t{end_i}\t{extra}\n")
                kept += 1
    return {"status": "filtered", "path": str(dst_bed), "rows": kept}


def _derive_low_mappability_bed(
    src_bedgraph: Path,
    dst_bed: Path,
    threshold: float,
    force: bool,
) -> dict[str, Any]:
    if dst_bed.exists() and not force:
        return {"status": "skipped_exists", "path": str(dst_bed), "bytes": dst_bed.stat().st_size}

    dst_bed.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with src_bedgraph.open("r", encoding="utf-8") as in_handle, dst_bed.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                score = float(parts[3])
            except ValueError:
                continue
            if score >= float(threshold):
                continue
            chrom = parts[0]
            start = parts[1]
            end = parts[2]
            out_handle.write(f"{chrom}\t{start}\t{end}\n")
            kept += 1
    return {"status": "derived", "path": str(dst_bed), "rows": kept, "threshold": float(threshold)}


def _open_textmaybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _collect_intervals_for_merge(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with _open_textmaybe_gz(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = ""
            start0 = -1
            end0 = -1
            try:
                # BED-like
                if len(parts) >= 3 and parts[0].startswith("chr"):
                    chrom = parts[0]
                    start0 = int(parts[1])
                    end0 = int(parts[2])
                # UCSC table with leading bin
                elif len(parts) >= 4 and parts[1].startswith("chr"):
                    chrom = parts[1]
                    start0 = int(parts[2])
                    end0 = int(parts[3])
            except ValueError:
                continue
            if not chrom or end0 <= start0:
                continue
            intervals.setdefault(chrom, []).append((start0, end0))
    return intervals


def _merge_junk_exclusion_bed(
    segdup_bed: Path,
    low_map_bed: Path,
    gap_bed_like: Path,
    blacklist_bed_like: Path,
    out_merged_bed: Path,
    out_gap_bed: Path,
    out_blacklist_bed: Path,
    force: bool,
) -> dict[str, Any]:
    if out_merged_bed.exists() and out_gap_bed.exists() and out_blacklist_bed.exists() and not force:
        return {"status": "skipped_exists", "path": str(out_merged_bed), "bytes": out_merged_bed.stat().st_size}

    # Normalize gap and blacklist to plain BED for easier downstream reuse.
    out_gap_bed.parent.mkdir(parents=True, exist_ok=True)
    out_blacklist_bed.parent.mkdir(parents=True, exist_ok=True)
    with _open_textmaybe_gz(gap_bed_like) as in_handle, out_gap_bed.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0].startswith("chr"):
                chrom, start, end = parts[0], parts[1], parts[2]
            elif len(parts) >= 4 and parts[1].startswith("chr"):
                chrom, start, end = parts[1], parts[2], parts[3]
            else:
                continue
            out_handle.write(f"{chrom}\t{start}\t{end}\n")
    with _open_textmaybe_gz(blacklist_bed_like) as in_handle, out_blacklist_bed.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0].startswith("chr"):
                out_handle.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")

    merged: dict[str, list[tuple[int, int]]] = {}
    for src in (segdup_bed, low_map_bed, out_gap_bed, out_blacklist_bed):
        for chrom, rows in _collect_intervals_for_merge(src).items():
            merged.setdefault(chrom, []).extend(rows)
    out_merged_bed.parent.mkdir(parents=True, exist_ok=True)
    total_intervals = 0
    with out_merged_bed.open("w", encoding="utf-8") as handle:
        for chrom in sorted(merged):
            rows = sorted(merged[chrom], key=lambda x: (x[0], x[1]))
            cur_s = -1
            cur_e = -1
            for s, e in rows:
                if cur_s < 0:
                    cur_s, cur_e = s, e
                    continue
                if s <= cur_e:
                    if e > cur_e:
                        cur_e = e
                else:
                    handle.write(f"{chrom}\t{cur_s}\t{cur_e}\n")
                    total_intervals += 1
                    cur_s, cur_e = s, e
            if cur_s >= 0:
                handle.write(f"{chrom}\t{cur_s}\t{cur_e}\n")
                total_intervals += 1
    return {
        "status": "merged",
        "path": str(out_merged_bed),
        "rows": total_intervals,
        "bytes": out_merged_bed.stat().st_size if out_merged_bed.exists() else 0,
        "gap_bed_path": str(out_gap_bed),
        "blacklist_bed_path": str(out_blacklist_bed),
    }


def _liftover_bed(in_bed: Path, chain_gz: Path, out_bed: Path, out_unmapped: Path, force: bool) -> dict[str, Any]:
    if out_bed.exists() and out_unmapped.exists() and not force:
        return {"status": "skipped_exists", "path": str(out_bed), "unmapped_path": str(out_unmapped)}

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    out_unmapped.parent.mkdir(parents=True, exist_ok=True)
    ok, msg = _run_cmd(
        ["liftOver", str(in_bed), str(chain_gz), str(out_bed), str(out_unmapped)],
        required=True,
    )
    return {"status": "lifted", "ok": ok, "message": msg, "path": str(out_bed), "unmapped_path": str(out_unmapped)}


def _ensure_famdb_repo(famdb_tools_dir: Path, timeout_sec: int) -> tuple[str | None, str]:
    famdb_py = famdb_tools_dir / "famdb.py"
    if famdb_py.exists():
        return str(famdb_py), "existing_repo"

    famdb_tools_dir.parent.mkdir(parents=True, exist_ok=True)
    ok, msg = _run_cmd(
        ["git", "clone", "https://github.com/Dfam-consortium/FamDB.git", str(famdb_tools_dir)],
        required=False,
    )
    if ok and famdb_py.exists():
        return str(famdb_py), "cloned_repo"

    system_famdb = shutil.which("famdb.py")
    if system_famdb:
        return system_famdb, "system_path"
    return None, f"missing_famdb_tool ({msg})"


def _export_human_mei_subset(curated_fasta: Path, mei_subset_fasta: Path) -> dict[str, Any]:
    header_re = re.compile(r"#SINE/Alu|#LINE/L1|(^|[^A-Za-z])SVA([^A-Za-z]|$)", re.IGNORECASE)
    kept = 0
    in_families = 0
    keep_seq = False
    with curated_fasta.open("r", encoding="utf-8") as in_handle, mei_subset_fasta.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if line.startswith(">"):
                in_families += 1
                keep_seq = bool(header_re.search(line))
                if keep_seq:
                    kept += 1
                    out_handle.write(line)
            elif keep_seq:
                out_handle.write(line)
    return {
        "input_families": in_families,
        "mei_families": kept,
        "output_path": str(mei_subset_fasta),
        "bytes": mei_subset_fasta.stat().st_size if mei_subset_fasta.exists() else 0,
    }


def _prepare_dfam_mei_library(outdir: Path, ds_map: dict[str, Dataset], timeout_sec: int, force: bool) -> dict[str, Any]:
    ds = ds_map.get("dfam_human_families_embl")
    if not ds:
        return {"status": "skipped_missing_dataset", "reason": "dfam_human_families_embl not in config"}

    archive_path = outdir / ds.target_path
    if not archive_path.exists():
        return {"status": "skipped_missing_archive", "path": str(archive_path)}

    decompressed_h5 = archive_path.with_suffix("")
    if decompressed_h5.suffix != ".h5":
        decompressed_h5 = archive_path.with_name("dfam40.0.h5")
    if archive_path.suffix == ".gz":
        _decompress_gzip(archive_path, decompressed_h5, force=force)

    dfam_dir = archive_path.parent
    db_dir = dfam_dir / "db"
    db_dir.mkdir(parents=True, exist_ok=True)

    db_base_h5 = db_dir / "dfam40.0.h5"
    if force or not db_base_h5.exists():
        shutil.copy2(decompressed_h5, db_base_h5)

    curated_consensus_h5 = db_dir / "dfam40.curated.consensus.0.h5"
    _download_file(DFAM_CURATED_CONSENSUS_0_URL, curated_consensus_h5.with_suffix(".h5.gz"), timeout_sec=timeout_sec, force=force)
    _decompress_gzip(curated_consensus_h5.with_suffix(".h5.gz"), curated_consensus_h5, force=force)

    famdb_tools_dir = dfam_dir / "tools" / "FamDB"
    famdb_py, tool_status = _ensure_famdb_repo(famdb_tools_dir, timeout_sec=timeout_sec)
    if famdb_py is None:
        return {"status": "failed", "error": tool_status, "db_dir": str(db_dir)}

    curated_fasta = dfam_dir / "dfam_human_curated.fasta"
    curated_missing_or_empty = (not curated_fasta.exists()) or (curated_fasta.stat().st_size <= 0)
    if force or curated_missing_or_empty:
        cmd = [
            sys.executable,
            famdb_py,
            "-i",
            str(db_dir),
            "families",
            "-f",
            "fasta_name",
            "--include-class-in-name",
            "-ad",
            "--curated",
            "9606",
        ]
        with curated_fasta.open("w", encoding="utf-8") as out_handle:
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"famdb export failed: {' '.join(cmd)}\n{proc.stderr}")
            out_handle.write(proc.stdout)
    curated_bytes = curated_fasta.stat().st_size if curated_fasta.exists() else 0
    if curated_bytes <= 0:
        raise RuntimeError(
            f"dfam curated FASTA is empty after export: {curated_fasta}. "
            "Re-run with --force and verify FamDB export output."
        )

    mei_subset_fasta = dfam_dir / "dfam_human_mei_l1_alu_sva.fasta"
    subset_stats = _export_human_mei_subset(curated_fasta, mei_subset_fasta)
    if int(subset_stats.get("bytes", 0)) <= 0 or int(subset_stats.get("mei_families", 0)) <= 0:
        raise RuntimeError(
            "dfam MEI subset FASTA is empty after export "
            f"(path={mei_subset_fasta}, families={subset_stats.get('mei_families', 0)})."
        )
    return {
        "status": "prepared",
        "famdb_tool": tool_status,
        "db_dir": str(db_dir),
        "curated_fasta": str(curated_fasta),
        "curated_fasta_bytes": curated_bytes,
        **subset_stats,
    }


def _postprocess(
    outdir: Path,
    datasets: list[Dataset],
    force: bool,
    timeout_sec: int,
    prefer_prebuilt_bwa_index: bool,
) -> list[dict[str, Any]]:
    ds_map = {d.dataset_id: d for d in datasets}
    steps: list[dict[str, Any]] = []

    def run_step(step_name: str, fn):
        print(f"[postprocess:start] {step_name}", file=sys.stderr)
        try:
            result = fn()
            status = result.get("status", "ok")
            print(f"[postprocess:done] {step_name} status={status}", file=sys.stderr)
            payload = {"step": step_name, **result}
            steps.append(payload)
            return payload
        except Exception as err:  # noqa: BLE001
            print(f"[postprocess:failed] {step_name} error={err}", file=sys.stderr)
            payload = {"step": step_name, "status": "failed", "error": str(err)}
            steps.append(payload)
            raise

    # 1) Reference FASTA prep: decompress + index files (.fai and .dict).
    for ref_id in ("hg38_reference_fasta", "hs1_t2t_reference_fasta"):
        ds = ds_map.get(ref_id)
        if not ds:
            continue
        source_path = outdir / ds.target_path
        if source_path.suffix == ".gz":
            fa_path = source_path.with_suffix("")
            run_step(f"{ref_id}_decompress", lambda: _decompress_gzip(source_path, fa_path, force=force))
        else:
            fa_path = source_path
            steps.append({"step": f"{ref_id}_decompress", "status": "skipped_uncompressed", "path": str(fa_path)})
            print(f"[postprocess:done] {ref_id}_decompress status=skipped_uncompressed", file=sys.stderr)
        run_step(
            f"{ref_id}_indexing",
            lambda: {
                "status": "ok",
                "results": _build_fasta_indexes(
                    fa_path,
                    fasta_url=ds.url,
                    ref_label=ref_id.replace("_reference_fasta", "").replace("_t2t", ""),
                    force=force,
                    timeout_sec=timeout_sec,
                    prefer_prebuilt_bwa_index=prefer_prebuilt_bwa_index,
                ),
            },
        )

    # 2) Convert hg38 annotations to BED.
    rmsk_gz = outdir / ds_map["hg38_repeatmasker_rmsk"].target_path
    rmsk_bed = outdir / "annotation/hg38/repeats/rmsk.bed"
    run_step("rmsk_to_bed", lambda: _convert_ucsc_txt_to_bed(rmsk_gz, rmsk_bed, force=force))

    segdup_gz = outdir / ds_map["hg38_segmental_duplications"].target_path
    segdup_bed = outdir / "annotation/hg38/segdup/genomicSuperDups.bed"
    run_step("segdup_to_bed", lambda: _convert_ucsc_txt_to_bed(segdup_gz, segdup_bed, force=force))

    # 3) Convert UCSC binary tracks (requires UCSC tools) + derive MEI-only BED.
    mappability_bw = outdir / ds_map["hg38_mappability_100bp_bw"].target_path
    mappability_bedgraph = outdir / "annotation/hg38/mappability/k100.Umap.MultiTrackMappability.bedGraph"
    if not mappability_bedgraph.exists() or force:
        run_step(
            "mappability_bw_to_bedgraph",
            lambda: {
                "status": "converted" if _run_cmd(["bigWigToBedGraph", str(mappability_bw), str(mappability_bedgraph)], required=True)[0] else "failed",
                "path": str(mappability_bedgraph),
            },
        )
    else:
        steps.append({"step": "mappability_bw_to_bedgraph", "status": "skipped_exists", "path": str(mappability_bedgraph)})
        print("[postprocess:done] mappability_bw_to_bedgraph status=skipped_exists", file=sys.stderr)
    low_map_bed = outdir / "annotation/hg38/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
    run_step(
        "mappability_low_to_bed",
        lambda: _derive_low_mappability_bed(
            src_bedgraph=mappability_bedgraph,
            dst_bed=low_map_bed,
            threshold=0.5,
            force=force,
        ),
    )
    gap_src = outdir / "annotation/hg38/masks/gap.txt.gz"
    blacklist_src = outdir / "annotation/hg38/blacklist/ENCFF356LFX.bed.gz"
    gap_bed = outdir / "annotation/hg38/masks/gap.bed"
    blacklist_bed = outdir / "annotation/hg38/blacklist/ENCFF356LFX.bed"
    merged_junk_bed = outdir / "annotation/hg38/junk/junk_exclusion_merged.bed"
    run_step(
        "merge_junk_exclusion_bed",
        lambda: _merge_junk_exclusion_bed(
            segdup_bed=segdup_bed,
            low_map_bed=low_map_bed,
            gap_bed_like=gap_src,
            blacklist_bed_like=blacklist_src,
            out_merged_bed=merged_junk_bed,
            out_gap_bed=gap_bed,
            out_blacklist_bed=blacklist_bed,
            force=force,
        ),
    )

    gnomad_bb = outdir / ds_map["gnomad_v41_sv_non_neuro_bb"].target_path
    gnomad_bed = outdir / "polymorphism/hg38/gnomad/gnomad.v4.1.sv.non_neuro_controls.sites.bed"
    if not gnomad_bed.exists() or force:
        run_step(
            "gnomad_bb_to_bed",
            lambda: {
                "status": "converted" if _run_cmd(["bigBedToBed", str(gnomad_bb), str(gnomad_bed)], required=True)[0] else "failed",
                "path": str(gnomad_bed),
            },
        )
    else:
        steps.append({"step": "gnomad_bb_to_bed", "status": "skipped_exists", "path": str(gnomad_bed)})
        print("[postprocess:done] gnomad_bb_to_bed status=skipped_exists", file=sys.stderr)

    mei_bed = outdir / "polymorphism/hg38/gnomad/gnomad.v4.1.sv.non_neuro_controls.mei.bed"
    run_step("gnomad_extract_mei", lambda: _filter_mei_bed(gnomad_bed, mei_bed, force=force))

    # 4) LiftOver hg38 resources to hs1 (T2T) where applicable.
    chain_hg38_to_hs1 = outdir / ds_map["hg38_to_hs1_chain"].target_path
    lift_inputs = [rmsk_bed, segdup_bed, mei_bed]
    for bed_in in lift_inputs:
        hs1_dir = Path(str(bed_in).replace("/hg38/", "/hs1/")).parent
        out_bed = hs1_dir / bed_in.name.replace(".bed", ".hs1.bed")
        out_unmapped = hs1_dir / bed_in.name.replace(".bed", ".hs1.unmapped.bed")
        run_step(
            f"liftover_{bed_in.name}_to_hs1",
            lambda bed_in=bed_in, out_bed=out_bed, out_unmapped=out_unmapped: _liftover_bed(
                bed_in, chain_hg38_to_hs1, out_bed, out_unmapped, force=force
            ),
        )

    # 5) Prepare Dfam-derived MEI FASTA library (human curated + LINE1/Alu/SVA subset).
    run_step(
        "dfam_prepare_human_mei_fasta",
        lambda: _prepare_dfam_mei_library(outdir=outdir, ds_map=ds_map, timeout_sec=timeout_sec, force=force),
    )

    return steps


def main() -> int:
    default_outdir = os.environ.get("RTM_PUBLIC_DATA_DIR", str(Path.home() / "retrotransposon-workdir" / "data" / "public"))
    parser = argparse.ArgumentParser(description="Download public references/annotations for retrotransposon-miner")
    parser.add_argument(
        "--config",
        default="resources/public_datasets.yaml",
        help="YAML dataset config path",
    )
    parser.add_argument(
        "--outdir",
        default=default_outdir,
        help="Base output directory",
    )
    parser.add_argument(
        "--categories",
        nargs="*",
        default=None,
        help="Optional category filter, e.g. reference annotation liftover",
    )
    parser.add_argument(
        "--include-optional",
        action="store_true",
        help="Include datasets marked required: false",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if target file exists",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=120,
        help="Per-request timeout",
    )
    parser.add_argument(
        "--manifest-path",
        default=None,
        help="Optional explicit manifest path (JSON). Defaults to <outdir>/manifest.json",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Thread count for slicing/indexing steps.",
    )
    parser.add_argument(
        "--skip-postprocess",
        action="store_true",
        help="Skip post-download preparation (fasta indexing, conversion, hg38->hs1 liftOver).",
    )
    parser.add_argument(
        "--no-prebuilt-bwa-index",
        action="store_true",
        help="Disable prebuilt hg38 BWA index download and always build indexes locally.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    outdir = Path(args.outdir).resolve()
    manifest_path = Path(args.manifest_path).resolve() if args.manifest_path else outdir / "manifest.json"

    datasets = _load_config(config_path)
    categories = set(args.categories) if args.categories else None

    outdir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, Any] = {
        "generated_at_unix": int(time.time()),
        "config": str(config_path),
        "outdir": str(outdir),
        "results": [],
        "postprocess": [],
        "summary": {},
    }

    print(
        "Storage note: full-BAM workflows (whole disease+control remap) may require ~300GB free disk. "
        "This downloader uses chr22 remote slicing for test BAMs to reduce footprint.",
        file=sys.stderr,
    )

    failures = 0
    for ds in datasets:
        if categories is not None and ds.category not in categories:
            continue
        if not ds.required and not args.include_optional:
            continue

        target = outdir / ds.target_path
        print(f"[download] {ds.dataset_id} -> {target}", file=sys.stderr)
        try:
            if ds.region and ds.url.endswith(".bam"):
                result = _slice_remote_bam(ds.url, ds.region, target, threads=args.threads, force=args.force)
            else:
                result = _download_file(ds.url, target, timeout_sec=args.timeout_sec, force=args.force)
            result.update(
                {
                    "id": ds.dataset_id,
                    "category": ds.category,
                    "description": ds.description,
                    "source": ds.source,
                    "region": ds.region,
                    "required": ds.required,
                }
            )
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as err:
            failures += 1
            result = {
                "id": ds.dataset_id,
                "category": ds.category,
                "description": ds.description,
                "source": ds.source,
                "required": ds.required,
                "status": "failed",
                "error": str(err),
                "url": ds.url,
                "path": str(target),
            }
            print(f"[error] {ds.dataset_id}: {err}", file=sys.stderr)

        manifest["results"].append(result)
        print(
            f"[download:{result.get('status', 'unknown')}] {ds.dataset_id} -> {target}",
            file=sys.stderr,
        )

    if not args.skip_postprocess:
        try:
            _ensure_postprocess_binaries()
        except Exception as err:  # noqa: BLE001
            manifest["postprocess"].append({"step": "postprocess_dependency_check", "status": "failed", "error": str(err)})
            manifest["summary"] = {
                "downloaded": sum(1 for r in manifest["results"] if r["status"] == "downloaded"),
                "skipped_exists": sum(1 for r in manifest["results"] if r["status"] == "skipped_exists"),
                "failed_required": sum(1 for r in manifest["results"] if r["status"] == "failed" and r["required"]),
                "failed_optional": sum(1 for r in manifest["results"] if r["status"] == "failed" and not r["required"]),
                "postprocess_failed": 1,
            }
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            with manifest_path.open("w", encoding="utf-8") as handle:
                json.dump(manifest, handle, indent=2)
                handle.write("\n")
            print(f"Wrote manifest: {manifest_path}")
            print(json.dumps(manifest["summary"], indent=2))
            return 1
        try:
            manifest["postprocess"] = _postprocess(
                outdir,
                datasets,
                force=args.force,
                timeout_sec=args.timeout_sec,
                prefer_prebuilt_bwa_index=not args.no_prebuilt_bwa_index,
            )
        except Exception as err:  # noqa: BLE001
            manifest["postprocess"].append({"step": "pipeline_postprocess", "status": "failed", "error": str(err)})

    downloaded = sum(1 for r in manifest["results"] if r["status"] == "downloaded")
    skipped = sum(1 for r in manifest["results"] if r["status"] == "skipped_exists")
    failed_required = sum(1 for r in manifest["results"] if r["status"] == "failed" and r["required"])
    failed_optional = sum(1 for r in manifest["results"] if r["status"] == "failed" and not r["required"])
    post_failed = sum(1 for r in manifest["postprocess"] if r.get("status") == "failed")

    manifest["summary"] = {
        "downloaded": downloaded,
        "skipped_exists": skipped,
        "failed_required": failed_required,
        "failed_optional": failed_optional,
        "postprocess_failed": post_failed,
    }

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    print(f"Wrote manifest: {manifest_path}")
    print(json.dumps(manifest["summary"], indent=2))

    return 1 if failed_required > 0 or post_failed > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
