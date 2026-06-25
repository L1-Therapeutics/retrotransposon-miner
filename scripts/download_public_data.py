#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from socket import timeout as SocketTimeoutError
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
UCSC_REPEATBROWSER_HG38REPS_URL = "https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38reps/hg38reps.fa"
SUPPORTED_REFERENCES = ("hg19", "hg38", "hs1")

# Datasets grouped by their native coordinate build.
REFERENCE_DATASET_IDS_BY_BUILD: dict[str, set[str]] = {
    "hg38": {
        "hg38_reference_fasta",
        "hg38_repeatmasker_rmsk",
        "hg38_segmental_duplications",
        "hg38_mappability_100bp_bw",
        "encode_blacklist_grch38",
        "hg38_gap",
        "gnomad_v41_sv_non_neuro_bb",
        "melt_1kg_mei_nstd144_vcf",
        "melt_1kg_mei_nstd144_vcf_tbi",
    },
    "hg19": {
        "hg19_reference_fasta",
        "hg19_repeatmasker_rmsk",
        "hg19_segmental_duplications",
        "hg19_mappability_100bp_bw",
        "hg19_gap",
    },
    "hs1": {
        "hs1_t2t_reference_fasta",
        "hs1_repeatmasker_out",
        "hs1_mappability_100bp_bw",
        "lr_1kg_ont_vienna_svan_mei_vcf",
        "lr_1kg_ont_vienna_svan_mei_vcf_tbi",
    },
}

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


def _resolve_selected_references(raw_refs: list[str] | None) -> tuple[str, ...]:
    if not raw_refs:
        return SUPPORTED_REFERENCES
    norm = []
    for ref in raw_refs:
        r = ref.strip().lower()
        if r not in SUPPORTED_REFERENCES:
            raise ValueError(f"Unsupported reference '{ref}'. Expected one of: {', '.join(SUPPORTED_REFERENCES)}")
        if r not in norm:
            norm.append(r)
    return tuple(norm)


def _select_dataset_ids_for_references(selected_refs: tuple[str, ...]) -> set[str]:
    ids: set[str] = set()

    # Build-native datasets for selected targets.
    for ref in selected_refs:
        ids.update(REFERENCE_DATASET_IDS_BY_BUILD.get(ref, set()))

    # MEI library + benchmark/helper resources are build-agnostic in this workflow.
    ids.update(
        {
            "dfam_human_families_embl",
            "ucsc_repeatbrowser_hg38reps_fasta",
        }
    )

    # Equivalent polymorphism outputs for every selected target require source datasets + chains.
    if "hg38" in selected_refs:
        ids.add("lr_1kg_ont_vienna_svan_mei_vcf")
        ids.add("lr_1kg_ont_vienna_svan_mei_vcf_tbi")
        ids.add("hs1_to_hg38_chain")

    if "hg19" in selected_refs:
        ids.update(
            {
                "lr_1kg_ont_vienna_svan_mei_vcf",
                "lr_1kg_ont_vienna_svan_mei_vcf_tbi",
                "hs1_to_hg19_chain",
                "hg38_to_hg19_chain",
                # MELT and gnomAD hg19 projections use hg38 sources.
                "melt_1kg_mei_nstd144_vcf",
                "melt_1kg_mei_nstd144_vcf_tbi",
                "gnomad_v41_sv_non_neuro_bb",
                "hg38_repeatmasker_rmsk",
                "hg38_segmental_duplications",
                "hg38_mappability_100bp_bw",
                "encode_blacklist_grch38",
                "hg38_gap",
            }
        )

    if "hs1" in selected_refs:
        ids.update(
            {
                # MELT and gnomAD hs1 projections use hg38 sources.
                "melt_1kg_mei_nstd144_vcf",
                "melt_1kg_mei_nstd144_vcf_tbi",
                "gnomad_v41_sv_non_neuro_bb",
                "hg38_repeatmasker_rmsk",
                "hg38_segmental_duplications",
                "hg38_mappability_100bp_bw",
                "encode_blacklist_grch38",
                "hg38_gap",
                "hg38_to_hs1_chain",
            }
        )

    return ids


def _download_file(url: str, out_path: Path, timeout_sec: int, force: bool) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        # Guard against previously interrupted/corrupted cached gzip files.
        if out_path.suffix == ".gz":
            try:
                with gzip.open(out_path, "rb") as gh:
                    _ = gh.read(1)
            except Exception:  # noqa: BLE001
                try:
                    out_path.unlink()
                except OSError:
                    pass
                # Continue to fresh download below.
            else:
                return {
                    "status": "skipped_exists",
                    "bytes": out_path.stat().st_size,
                    "path": str(out_path),
                    "url": url,
                }
        if out_path.exists():
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


def _download_dataset(
    ds: Dataset,
    target: Path,
    timeout_sec: int,
    threads: int,
    force: bool,
    retries: int,
    retry_backoff_sec: float,
) -> dict[str, Any]:
    attempts = max(1, int(retries))
    last_err: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            if ds.region and ds.url.endswith(".bam"):
                result = _slice_remote_bam(ds.url, ds.region, target, threads=threads, force=force)
            else:
                result = _download_file(ds.url, target, timeout_sec=timeout_sec, force=force)
            break
        except (
            urllib.error.URLError,
            urllib.error.HTTPError,
            TimeoutError,
            SocketTimeoutError,
            ConnectionError,
        ) as err:
            last_err = err
            if attempt >= attempts:
                raise
            sleep_s = float(retry_backoff_sec) * (2 ** (attempt - 1))
            time.sleep(max(0.0, sleep_s))
    else:
        if last_err is not None:
            raise last_err
        raise RuntimeError(f"Download failed unexpectedly for dataset {ds.dataset_id}")

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
    return result


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


def _extract_vcf_end_pos(pos1: int, ref: str, info_field: str) -> int:
    default_end = int(pos1) + max(1, len(ref or "")) - 1
    if not info_field or info_field == ".":
        return max(default_end, pos1)
    for token in info_field.split(";"):
        if token.startswith("END="):
            try:
                return max(int(token[4:]), pos1)
            except ValueError:
                return max(default_end, pos1)
    return max(default_end, pos1)


def _replace_vcf_end_pos(info_field: str, end1: int) -> str:
    target = f"END={int(end1)}"
    if not info_field or info_field == ".":
        return target
    parts = info_field.split(";")
    replaced = False
    for i, token in enumerate(parts):
        if token.startswith("END="):
            parts[i] = target
            replaced = True
            break
    if not replaced:
        parts.append(target)
    return ";".join(parts)


def _vcf_chrom_sort_key(chrom: str) -> tuple[int, str]:
    c = (chrom or "").strip()
    if c.startswith("chr"):
        core = c[3:]
    else:
        core = c
    if core.isdigit():
        n = int(core)
        if 1 <= n <= 22:
            return (n, c)
    if core == "X":
        return (23, c)
    if core == "Y":
        return (24, c)
    if core in {"M", "MT"}:
        return (25, c)
    return (1000, c)


def _liftover_vcf_to_target(
    source_vcf_gz: Path,
    chain_gz: Path,
    out_vcf_gz: Path,
    out_unmapped_bed: Path,
    force: bool,
    source_build: str,
    target_build: str,
    target_reference_meta: str,
) -> dict[str, Any]:
    if not source_vcf_gz.exists():
        return {"status": "skipped_missing_source", "source_path": str(source_vcf_gz)}
    if out_vcf_gz.exists() and out_unmapped_bed.exists() and not force:
        return {
            "status": "skipped_exists",
            "path": str(out_vcf_gz),
            "unmapped_path": str(out_unmapped_bed),
            "source_path": str(source_vcf_gz),
        }

    out_vcf_gz.parent.mkdir(parents=True, exist_ok=True)
    out_unmapped_bed.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"rtm_vcf_{source_build}_to_{target_build}_") as tmpdir:
        tmp = Path(tmpdir)
        query_bed = tmp / f"variants_{source_build}_intervals.bed"
        mapped_bed = tmp / f"variants_{source_build}_to_{target_build}.mapped.bed"
        unmapped_bed = tmp / f"variants_{source_build}_to_{target_build}.unmapped.bed"
        lifted_vcf = tmp / f"variants.{target_build}.vcf"

        metadata_headers: list[str] = []
        chrom_header = ""
        variants: list[list[str]] = []
        query_ids: list[str] = []

        with gzip.open(source_vcf_gz, "rt", encoding="utf-8") as ih, query_bed.open("w", encoding="utf-8") as bh:
            for line in ih:
                if line.startswith("##"):
                    metadata_headers.append(line)
                    continue
                if line.startswith("#CHROM"):
                    chrom_header = line
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                chrom = cols[0]
                try:
                    pos1 = int(cols[1])
                except ValueError:
                    continue
                ref = cols[3]
                end1 = _extract_vcf_end_pos(pos1=pos1, ref=ref, info_field=cols[7])
                start0 = max(0, pos1 - 1)
                end0 = max(start0 + 1, end1)
                rid = f"v{len(variants)}"
                bh.write(f"{chrom}\t{start0}\t{end0}\t{rid}\n")
                variants.append(cols)
                query_ids.append(rid)

        if not chrom_header:
            return {"status": "failed", "error": f"invalid_vcf_header: {source_vcf_gz}"}

        _run_cmd(
            ["liftOver", str(query_bed), str(chain_gz), str(mapped_bed), str(unmapped_bed)],
            required=True,
        )

        mapped_coords: dict[str, tuple[str, int, int]] = {}
        multi_mapped = 0
        with mapped_bed.open("r", encoding="utf-8") as mh:
            for line in mh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                rid = parts[3]
                try:
                    start0 = int(parts[1])
                    end0 = int(parts[2])
                except ValueError:
                    continue
                if rid in mapped_coords:
                    multi_mapped += 1
                    continue
                mapped_coords[rid] = (parts[0], start0, end0)

        with lifted_vcf.open("w", encoding="utf-8") as oh:
            saw_reference = False
            saw_liftover_meta = False
            for line in metadata_headers:
                if line.startswith("##reference="):
                    oh.write(f"##reference={target_reference_meta}\n")
                    saw_reference = True
                    continue
                if line.startswith("##rtm_liftover_source="):
                    saw_liftover_meta = True
                oh.write(line)
            if not saw_reference:
                oh.write(f"##reference={target_reference_meta}\n")
            if not saw_liftover_meta:
                oh.write(f"##rtm_liftover_source={source_build}_to_{target_build}_with_ucsc_liftOver\n")
            oh.write(chrom_header)

            lifted_rows: list[list[str]] = []
            for rid, cols in zip(query_ids, variants, strict=False):
                mapped = mapped_coords.get(rid)
                if mapped is None:
                    continue
                new_chrom, new_start0, new_end0 = mapped
                cols[0] = new_chrom
                cols[1] = str(max(1, new_start0 + 1))
                cols[7] = _replace_vcf_end_pos(cols[7], max(new_start0 + 1, new_end0))
                lifted_rows.append(cols)

            lifted_rows.sort(key=lambda row: (_vcf_chrom_sort_key(row[0]), int(row[1])))
            for cols in lifted_rows:
                oh.write("\t".join(cols) + "\n")
            written = len(lifted_rows)

        with out_vcf_gz.open("wb") as out_handle:
            proc = subprocess.run(
                ["bgzip", "-f", "-c", str(lifted_vcf)],
                check=True,
                stdout=out_handle,
                stderr=subprocess.PIPE,
            )
            _ = proc

        tabix_ok, tabix_msg = _run_cmd(["tabix", "-f", "-p", "vcf", str(out_vcf_gz)], required=False)
        shutil.copy2(unmapped_bed, out_unmapped_bed)

        return {
            "status": "lifted",
            "source_path": str(source_vcf_gz),
            "path": str(out_vcf_gz),
            "unmapped_path": str(out_unmapped_bed),
            "records_total": len(variants),
            "records_lifted": written,
            "records_unmapped": max(0, len(variants) - written),
            "multi_mapped_skipped": int(multi_mapped),
            "tabix_indexed": bool(tabix_ok),
            "tabix_message": tabix_msg,
        }


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


def _iter_fasta(path: Path):
    header = None
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts).upper()


def _mei_family_from_header(header: str) -> str:
    u = (header or "").upper()
    if "SVA" in u:
        return "SVA"
    if "ALU" in u:
        return "ALU"
    if "LINE/L1" in u or "L1" in u:
        return "LINE1"
    return "OTHER"


def _normalize_dfam_subfamily_for_full(header: str) -> str:
    base = (header or "").split("#", 1)[0]
    # Dfam subset contains fragment-oriented labels (3end/5end/orf2/short).
    # Convert to the corresponding full-length RepeatBrowser subfamily token.
    base = re.sub(r"_(3end|5end|orf2)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"_short_?$", "", base, flags=re.IGNORECASE)
    return base


def _build_mei_full_consensus_panel(outdir: Path, ds_map: dict[str, Dataset], timeout_sec: int, force: bool) -> dict[str, Any]:
    dfam_subset = outdir / "retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta"
    if not dfam_subset.exists():
        return {"status": "skipped_missing_dfam_subset", "path": str(dfam_subset)}

    ucsc_ds = ds_map.get("ucsc_repeatbrowser_hg38reps_fasta")
    ucsc_rel = "retrotransposon_db/ucsc_repeatbrowser/hg38reps.fa"
    ucsc_url = UCSC_REPEATBROWSER_HG38REPS_URL
    if ucsc_ds is not None:
        ucsc_rel = ucsc_ds.target_path
        ucsc_url = ucsc_ds.url
    ucsc_fa = outdir / ucsc_rel
    _download_file(ucsc_url, ucsc_fa, timeout_sec=timeout_sec, force=force)
    if not ucsc_fa.exists() or ucsc_fa.stat().st_size <= 0:
        return {"status": "failed", "error": "missing_ucsc_repeatbrowser_fasta", "path": str(ucsc_fa)}

    ucsc_sequences: dict[str, str] = {}
    for hdr, seq in _iter_fasta(ucsc_fa):
        if seq:
            ucsc_sequences[hdr] = seq

    out_dir = outdir / "retrotransposon_db/full_consensus"
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_fa = out_dir / "mei_full_canonical.panel.fa"
    panel_meta = out_dir / "mei_full_canonical.panel.tsv"

    min_line1_full_bp = 5000
    selected_records: list[tuple[str, str, str, str, int]] = []
    missing_records: list[tuple[str, str, str]] = []
    excluded_short_line1: list[tuple[str, str, int]] = []
    seen_norm: set[tuple[str, str]] = set()
    for hdr, _seq in _iter_fasta(dfam_subset):
        fam = _mei_family_from_header(hdr)
        if fam not in {"ALU", "SVA", "LINE1"}:
            continue
        cls = hdr.split("#", 1)[1] if "#" in hdr else ""
        norm = _normalize_dfam_subfamily_for_full(hdr)
        if not norm:
            continue
        key = (fam, norm)
        if key in seen_norm:
            continue
        seen_norm.add(key)
        full_seq = ucsc_sequences.get(norm, "")
        if not full_seq:
            missing_records.append((hdr, norm, fam))
            continue
        if fam == "LINE1" and len(full_seq) < min_line1_full_bp:
            excluded_short_line1.append((hdr, norm, len(full_seq)))
            continue
        out_hdr = f"{norm}_full#{cls}" if cls else f"{norm}_full"
        selected_records.append((out_hdr, norm, hdr, fam, len(full_seq)))

    if not selected_records:
        return {
            "status": "failed",
            "error": "no_full_consensus_records_selected",
            "dfam_subset": str(dfam_subset),
            "ucsc_fasta": str(ucsc_fa),
        }

    with panel_fa.open("w", encoding="utf-8") as out_handle:
        for out_hdr, norm, _src, _fam, _len in selected_records:
            seq = ucsc_sequences[norm]
            out_handle.write(f">{out_hdr}|source={norm}\n")
            for i in range(0, len(seq), 60):
                out_handle.write(seq[i : i + 60] + "\n")

    with panel_meta.open("w", encoding="utf-8", newline="") as meta_handle:
        writer = csv.writer(meta_handle, delimiter="\t")
        writer.writerow(["output_header", "source_repeatbrowser", "source_dfam_header", "family", "length_bp"])
        for rec in selected_records:
            writer.writerow(rec)
        writer.writerow([])
        writer.writerow([f"#excluded_short_line1_lt_{min_line1_full_bp}bp"])
        writer.writerow(["source_dfam_header", "normalized_repeatbrowser_name", "length_bp"])
        for src_hdr, norm, bp in excluded_short_line1:
            writer.writerow([src_hdr, norm, bp])
        writer.writerow([])
        writer.writerow(["#missing_dfam_subfamilies_not_found_in_repeatbrowser"])
        writer.writerow(["source_dfam_header", "normalized_repeatbrowser_name", "family"])
        for src_hdr, norm, fam in missing_records:
            writer.writerow([src_hdr, norm, fam])

    _run_cmd(["samtools", "faidx", str(panel_fa)], required=True)

    line1_lengths = [r[4] for r in selected_records if r[3] == "LINE1"]
    alu_lengths = [r[4] for r in selected_records if r[3] == "ALU"]
    sva_lengths = [r[4] for r in selected_records if r[3] == "SVA"]
    return {
        "status": "prepared",
        "panel_fasta": str(panel_fa),
        "panel_fai": str(panel_fa.with_suffix(panel_fa.suffix + ".fai")),
        "panel_metadata_tsv": str(panel_meta),
        "selected_total": len(selected_records),
        "selected_line1": len(line1_lengths),
        "selected_alu": len(alu_lengths),
        "selected_sva": len(sva_lengths),
        "missing_total": len(missing_records),
        "excluded_short_line1": len(excluded_short_line1),
        "line1_min_full_bp_threshold": min_line1_full_bp,
        "line1_length_min": min(line1_lengths) if line1_lengths else 0,
        "line1_length_max": max(line1_lengths) if line1_lengths else 0,
    }


def _postprocess(
    outdir: Path,
    datasets: list[Dataset],
    force: bool,
    timeout_sec: int,
    prefer_prebuilt_bwa_index: bool,
    selected_references: tuple[str, ...],
) -> list[dict[str, Any]]:
    ds_map = {d.dataset_id: d for d in datasets}
    steps: list[dict[str, Any]] = []
    want_hg38 = "hg38" in selected_references
    want_hg19 = "hg19" in selected_references
    want_hs1 = "hs1" in selected_references

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
    for ref_id in ("hg38_reference_fasta", "hg19_reference_fasta", "hs1_t2t_reference_fasta"):
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

    # 4) Prepare hg19 annotation/junk resources.
    # Prefer direct hg19 downloads when available in config; otherwise liftOver from hg38.
    hg38_to_hg19_ds = ds_map.get("hg38_to_hg19_chain")
    hg38_to_hg19_chain = outdir / hg38_to_hg19_ds.target_path if hg38_to_hg19_ds else None
    hg19_segdup_bed = outdir / "annotation/hg19/segdup/genomicSuperDups.bed"
    hg19_low_map_bed = outdir / "annotation/hg19/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
    if want_hg19:
        hg19_rmsk_ds = ds_map.get("hg19_repeatmasker_rmsk")
        hg19_rmsk_bed = outdir / "annotation/hg19/repeats/rmsk.bed"
        if hg19_rmsk_ds:
            hg19_rmsk_gz = outdir / hg19_rmsk_ds.target_path
            run_step("hg19_rmsk_to_bed", lambda: _convert_ucsc_txt_to_bed(hg19_rmsk_gz, hg19_rmsk_bed, force=force))
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_rmsk_to_hg19",
                lambda: _liftover_bed(
                    rmsk_bed,
                    hg38_to_hg19_chain,
                    hg19_rmsk_bed,
                    outdir / "annotation/hg19/repeats/rmsk.unmapped.bed",
                    force=force,
                ),
            )

        hg19_segdup_ds = ds_map.get("hg19_segmental_duplications")
        if hg19_segdup_ds:
            hg19_segdup_gz = outdir / hg19_segdup_ds.target_path
            run_step("hg19_segdup_to_bed", lambda: _convert_ucsc_txt_to_bed(hg19_segdup_gz, hg19_segdup_bed, force=force))
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_segdup_to_hg19",
                lambda: _liftover_bed(
                    segdup_bed,
                    hg38_to_hg19_chain,
                    hg19_segdup_bed,
                    outdir / "annotation/hg19/segdup/genomicSuperDups.unmapped.bed",
                    force=force,
                ),
            )

        hg19_mapp_bw_ds = ds_map.get("hg19_mappability_100bp_bw")
        hg19_mapp_bedgraph = outdir / "annotation/hg19/mappability/k100.Umap.MultiTrackMappability.bedGraph"
        if hg19_mapp_bw_ds:
            hg19_mapp_bw = outdir / hg19_mapp_bw_ds.target_path
            if hg19_mapp_bw.exists():
                if not hg19_mapp_bedgraph.exists() or force:
                    run_step(
                        "hg19_mappability_bw_to_bedgraph",
                        lambda: {
                            "status": "converted"
                            if _run_cmd(["bigWigToBedGraph", str(hg19_mapp_bw), str(hg19_mapp_bedgraph)], required=True)[0]
                            else "failed",
                            "path": str(hg19_mapp_bedgraph),
                        },
                    )
                else:
                    steps.append(
                        {
                            "step": "hg19_mappability_bw_to_bedgraph",
                            "status": "skipped_exists",
                            "path": str(hg19_mapp_bedgraph),
                        }
                    )
                    print("[postprocess:done] hg19_mappability_bw_to_bedgraph status=skipped_exists", file=sys.stderr)
            elif hg38_to_hg19_chain is not None:
                run_step(
                    "liftover_hg38_mappability_bedgraph_to_hg19",
                    lambda: _liftover_bed(
                        mappability_bedgraph,
                        hg38_to_hg19_chain,
                        hg19_mapp_bedgraph,
                        outdir / "annotation/hg19/mappability/k100.Umap.MultiTrackMappability.unmapped.bed",
                        force=force,
                    ),
                )
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_mappability_bedgraph_to_hg19",
                lambda: _liftover_bed(
                    mappability_bedgraph,
                    hg38_to_hg19_chain,
                    hg19_mapp_bedgraph,
                    outdir / "annotation/hg19/mappability/k100.Umap.MultiTrackMappability.unmapped.bed",
                    force=force,
                ),
            )

        if hg19_mapp_bedgraph.exists():
            run_step(
                "hg19_mappability_low_to_bed",
                lambda: _derive_low_mappability_bed(
                    src_bedgraph=hg19_mapp_bedgraph,
                    dst_bed=hg19_low_map_bed,
                    threshold=0.5,
                    force=force,
                ),
            )
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_low_mappability_to_hg19",
                lambda: _liftover_bed(
                    low_map_bed,
                    hg38_to_hg19_chain,
                    hg19_low_map_bed,
                    outdir / "annotation/hg19/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.unmapped.bed",
                    force=force,
                ),
            )

        hg19_gap_ds = ds_map.get("hg19_gap")
        hg19_gap_bed = outdir / "annotation/hg19/masks/gap.bed"
        hg19_gap_bed_like: Path = hg19_gap_bed
        if hg19_gap_ds:
            hg19_gap_bed_like = outdir / hg19_gap_ds.target_path
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_gap_to_hg19",
                lambda: _liftover_bed(
                    gap_bed,
                    hg38_to_hg19_chain,
                    hg19_gap_bed,
                    outdir / "annotation/hg19/masks/gap.unmapped.bed",
                    force=force,
                ),
            )

        hg19_blacklist_ds = ds_map.get("encode_blacklist_hg19")
        hg19_blacklist_bed = outdir / "annotation/hg19/blacklist/ENCFF356LFX.hg19.bed"
        hg19_blacklist_bed_like: Path = hg19_blacklist_bed
        if hg19_blacklist_ds:
            hg19_blacklist_bed_like = outdir / hg19_blacklist_ds.target_path
        elif hg38_to_hg19_chain is not None:
            run_step(
                "liftover_hg38_blacklist_to_hg19",
                lambda: _liftover_bed(
                    blacklist_bed,
                    hg38_to_hg19_chain,
                    hg19_blacklist_bed,
                    outdir / "annotation/hg19/blacklist/ENCFF356LFX.hg19.unmapped.bed",
                    force=force,
                ),
            )

        if hg19_segdup_bed.exists() and hg19_low_map_bed.exists():
            run_step(
                "merge_hg19_junk_exclusion_bed",
                lambda: _merge_junk_exclusion_bed(
                    segdup_bed=hg19_segdup_bed,
                    low_map_bed=hg19_low_map_bed,
                    gap_bed_like=hg19_gap_bed_like,
                    blacklist_bed_like=hg19_blacklist_bed_like,
                    out_merged_bed=outdir / "annotation/hg19/junk/junk_exclusion_merged.bed",
                    out_gap_bed=outdir / "annotation/hg19/masks/gap.bed",
                    out_blacklist_bed=outdir / "annotation/hg19/blacklist/ENCFF356LFX.hg19.bed",
                    force=force,
                ),
            )

    # 5) LiftOver hg38 resources to hs1 (T2T) where applicable.
    chain_hg38_to_hs1_ds = ds_map.get("hg38_to_hs1_chain")
    chain_hg38_to_hs1 = outdir / chain_hg38_to_hs1_ds.target_path if chain_hg38_to_hs1_ds else None
    hs1_segdup_bed = outdir / "annotation/hs1/segdup/genomicSuperDups.hs1.bed"
    hs1_low_map_bed = outdir / "annotation/hs1/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"
    if want_hs1 and chain_hg38_to_hs1 is not None:
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

        hs1_mapp_bw_ds = ds_map.get("hs1_mappability_100bp_bw")
        hs1_mapp_bedgraph = outdir / "annotation/hs1/mappability/k100.Umap.MultiTrackMappability.bedGraph"
        if hs1_mapp_bw_ds:
            hs1_mapp_bw = outdir / hs1_mapp_bw_ds.target_path
            if not hs1_mapp_bedgraph.exists() or force:
                run_step(
                    "hs1_mappability_bw_to_bedgraph",
                    lambda: {
                        "status": "converted"
                        if _run_cmd(["bigWigToBedGraph", str(hs1_mapp_bw), str(hs1_mapp_bedgraph)], required=True)[0]
                        else "failed",
                        "path": str(hs1_mapp_bedgraph),
                    },
                )
            else:
                steps.append(
                    {"step": "hs1_mappability_bw_to_bedgraph", "status": "skipped_exists", "path": str(hs1_mapp_bedgraph)}
                )
                print("[postprocess:done] hs1_mappability_bw_to_bedgraph status=skipped_exists", file=sys.stderr)
        else:
            run_step(
                "liftover_hg38_mappability_bedgraph_to_hs1",
                lambda: _liftover_bed(
                    mappability_bedgraph,
                    chain_hg38_to_hs1,
                    hs1_mapp_bedgraph,
                    outdir / "annotation/hs1/mappability/k100.Umap.MultiTrackMappability.unmapped.bed",
                    force=force,
                ),
            )
        if hs1_mapp_bedgraph.exists():
            run_step(
                "hs1_mappability_low_to_bed",
                lambda: _derive_low_mappability_bed(
                    src_bedgraph=hs1_mapp_bedgraph,
                    dst_bed=hs1_low_map_bed,
                    threshold=0.5,
                    force=force,
                ),
            )
        else:
            run_step(
                "liftover_hg38_low_mappability_to_hs1",
                lambda: _liftover_bed(
                    low_map_bed,
                    chain_hg38_to_hs1,
                    hs1_low_map_bed,
                    outdir / "annotation/hs1/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.unmapped.bed",
                    force=force,
                ),
            )

        hs1_gap_bed = outdir / "annotation/hs1/masks/gap.bed"
        hs1_blacklist_bed = outdir / "annotation/hs1/blacklist/ENCFF356LFX.hs1.bed"
        run_step(
            "liftover_hg38_gap_to_hs1",
            lambda: _liftover_bed(
                gap_bed,
                chain_hg38_to_hs1,
                hs1_gap_bed,
                outdir / "annotation/hs1/masks/gap.unmapped.bed",
                force=force,
            ),
        )
        run_step(
            "liftover_hg38_blacklist_to_hs1",
            lambda: _liftover_bed(
                blacklist_bed,
                chain_hg38_to_hs1,
                hs1_blacklist_bed,
                outdir / "annotation/hs1/blacklist/ENCFF356LFX.hs1.unmapped.bed",
                force=force,
            ),
        )
        if hs1_segdup_bed.exists() and hs1_low_map_bed.exists():
            run_step(
                "merge_hs1_junk_exclusion_bed",
                lambda: _merge_junk_exclusion_bed(
                    segdup_bed=hs1_segdup_bed,
                    low_map_bed=hs1_low_map_bed,
                    gap_bed_like=hs1_gap_bed,
                    blacklist_bed_like=hs1_blacklist_bed,
                    out_merged_bed=outdir / "annotation/hs1/junk/junk_exclusion_merged.bed",
                    out_gap_bed=hs1_gap_bed,
                    out_blacklist_bed=hs1_blacklist_bed,
                    force=force,
                ),
            )

    # 6) 1000G ONT Vienna SVAN is downloaded in hs1 coordinates; lift to hg38
    # and keep the historical hg38 filename for annotation compatibility.
    lr_ds = ds_map.get("lr_1kg_ont_vienna_svan_mei_vcf")
    hs1_to_hg38_chain = ds_map.get("hs1_to_hg38_chain")
    if want_hg38 and lr_ds and hs1_to_hg38_chain:
        lr_source_hs1_vcf = outdir / lr_ds.target_path
        lr_lifted_hg38_vcf = outdir / "polymorphism/hg38/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz"
        lr_lifted_unmapped = outdir / "polymorphism/hg38/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.unmapped.bed"
        run_step(
            "liftover_lr_1kg_ont_vienna_hs1_to_hg38_vcf",
            lambda: _liftover_vcf_to_target(
                source_vcf_gz=lr_source_hs1_vcf,
                chain_gz=outdir / hs1_to_hg38_chain.target_path,
                out_vcf_gz=lr_lifted_hg38_vcf,
                out_unmapped_bed=lr_lifted_unmapped,
                force=force,
                source_build="hs1",
                target_build="hg38",
                target_reference_meta="GCF_000001405.39",
            ),
        )

    # 6b) Also project the same Vienna SVAN callset to hg19 when chain is available.
    hs1_to_hg19_chain = ds_map.get("hs1_to_hg19_chain")
    if want_hg19 and lr_ds and hs1_to_hg19_chain:
        lr_source_hs1_vcf = outdir / lr_ds.target_path
        lr_lifted_hg19_vcf = outdir / "polymorphism/hg19/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.vcf.gz"
        lr_lifted_hg19_unmapped = (
            outdir / "polymorphism/hg19/long_read_1kg_ont_vienna/final-vcf.unphased.SVAN_1.3.unmapped.bed"
        )
        run_step(
            "liftover_lr_1kg_ont_vienna_hs1_to_hg19_vcf",
            lambda: _liftover_vcf_to_target(
                source_vcf_gz=lr_source_hs1_vcf,
                chain_gz=outdir / hs1_to_hg19_chain.target_path,
                out_vcf_gz=lr_lifted_hg19_vcf,
                out_unmapped_bed=lr_lifted_hg19_unmapped,
                force=force,
                source_build="hs1",
                target_build="hg19",
                target_reference_meta="GCF_000001405.25",
            ),
        )

    # 6c) Project MELT 1KG MEI callset from hg38 to hg19 when available.
    melt_ds = ds_map.get("melt_1kg_mei_nstd144_vcf")
    if want_hg19 and melt_ds and hg38_to_hg19_chain is not None:
        melt_hg38_vcf = outdir / melt_ds.target_path
        melt_hg19_vcf = outdir / "polymorphism/hg19/melt/nstd144.GRCh37.variant_call.vcf.gz"
        melt_hg19_unmapped = outdir / "polymorphism/hg19/melt/nstd144.GRCh37.variant_call.unmapped.bed"
        run_step(
            "liftover_melt_1kg_hg38_to_hg19_vcf",
            lambda: _liftover_vcf_to_target(
                source_vcf_gz=melt_hg38_vcf,
                chain_gz=hg38_to_hg19_chain,
                out_vcf_gz=melt_hg19_vcf,
                out_unmapped_bed=melt_hg19_unmapped,
                force=force,
                source_build="hg38",
                target_build="hg19",
                target_reference_meta="GCF_000001405.25",
            ),
        )

    if want_hs1 and melt_ds and chain_hg38_to_hs1 is not None:
        melt_hg38_vcf = outdir / melt_ds.target_path
        melt_hs1_vcf = outdir / "polymorphism/hs1/melt/nstd144.hs1.variant_call.vcf.gz"
        melt_hs1_unmapped = outdir / "polymorphism/hs1/melt/nstd144.hs1.variant_call.unmapped.bed"
        run_step(
            "liftover_melt_1kg_hg38_to_hs1_vcf",
            lambda: _liftover_vcf_to_target(
                source_vcf_gz=melt_hg38_vcf,
                chain_gz=chain_hg38_to_hs1,
                out_vcf_gz=melt_hs1_vcf,
                out_unmapped_bed=melt_hs1_unmapped,
                force=force,
                source_build="hg38",
                target_build="hs1",
                target_reference_meta="GCF_009914755.1",
            ),
        )

    if want_hg19 and hg38_to_hg19_chain is not None:
        gnomad_hg19_mei_bed = outdir / "polymorphism/hg19/gnomad/gnomad.v4.1.sv.non_neuro_controls.mei.bed"
        run_step(
            "liftover_gnomad_mei_hg38_to_hg19",
            lambda: _liftover_bed(
                mei_bed,
                hg38_to_hg19_chain,
                gnomad_hg19_mei_bed,
                outdir / "polymorphism/hg19/gnomad/gnomad.v4.1.sv.non_neuro_controls.mei.unmapped.bed",
                force=force,
            ),
        )

    # 7) Prepare Dfam-derived MEI FASTA library (human curated + LINE1/Alu/SVA subset).
    run_step(
        "dfam_prepare_human_mei_fasta",
        lambda: _prepare_dfam_mei_library(outdir=outdir, ds_map=ds_map, timeout_sec=timeout_sec, force=force),
    )
    run_step(
        "prepare_full_consensus_mei_panel",
        lambda: _build_mei_full_consensus_panel(outdir=outdir, ds_map=ds_map, timeout_sec=timeout_sec, force=force),
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
        "--references",
        nargs="+",
        default=None,
        help="Reference builds to materialize (any of: hg19 hg38 hs1). Default: all.",
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
        default=300,
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
        "--download-workers",
        type=int,
        default=6,
        help="Number of dataset downloads to run concurrently.",
    )
    parser.add_argument(
        "--download-retries",
        type=int,
        default=3,
        help="Retry attempts per dataset download on transient network failures.",
    )
    parser.add_argument(
        "--download-retry-backoff-sec",
        type=float,
        default=2.0,
        help="Initial retry backoff in seconds (exponential).",
    )
    parser.add_argument(
        "--skip-postprocess",
        action="store_true",
        help="Skip post-download preparation (indexing, track conversion, and cross-reference liftOver).",
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
    selected_references = _resolve_selected_references(args.references)
    selected_dataset_ids = _select_dataset_ids_for_references(selected_references)
    categories = set(args.categories) if args.categories else None
    active_datasets = [
        d
        for d in datasets
        if d.dataset_id in selected_dataset_ids
        and (categories is None or d.category in categories)
        and (d.required or args.include_optional)
    ]

    outdir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, Any] = {
        "generated_at_unix": int(time.time()),
        "config": str(config_path),
        "outdir": str(outdir),
        "selected_references": list(selected_references),
        "results": [],
        "postprocess": [],
        "summary": {},
    }

    print(
        "Storage note: full-BAM workflows (whole disease+control remap) may require ~300GB free disk. "
        "This downloader uses chr22 remote slicing for test BAMs to reduce footprint.",
        file=sys.stderr,
    )
    print(f"Selected references: {', '.join(selected_references)}", file=sys.stderr)

    failures = 0
    download_workers = max(1, int(args.download_workers))
    results_by_idx: list[dict[str, Any] | None] = [None] * len(active_datasets)
    future_to_meta: dict[concurrent.futures.Future[dict[str, Any]], tuple[int, Dataset, Path]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=download_workers) as executor:
        for idx, ds in enumerate(active_datasets):
            target = outdir / ds.target_path
            print(f"[download] {ds.dataset_id} -> {target}", file=sys.stderr)
            fut = executor.submit(
                _download_dataset,
                ds=ds,
                target=target,
                timeout_sec=args.timeout_sec,
                threads=args.threads,
                force=args.force,
                retries=args.download_retries,
                retry_backoff_sec=args.download_retry_backoff_sec,
            )
            future_to_meta[fut] = (idx, ds, target)

        for fut in concurrent.futures.as_completed(future_to_meta):
            idx, ds, target = future_to_meta[fut]
            try:
                result = fut.result()
            except Exception as err:  # noqa: BLE001
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

            results_by_idx[idx] = result
            print(
                f"[download:{result.get('status', 'unknown')}] {ds.dataset_id} -> {target}",
                file=sys.stderr,
            )

    for result in results_by_idx:
        if result is None:
            continue
        manifest["results"].append(result)
    if len(manifest["results"]) != len(active_datasets):
        raise RuntimeError(
            "Internal error: manifest result count mismatch after concurrent download "
            f"({len(manifest['results'])} vs {len(active_datasets)})."
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
                active_datasets,
                force=args.force,
                timeout_sec=args.timeout_sec,
                prefer_prebuilt_bwa_index=not args.no_prebuilt_bwa_index,
                selected_references=selected_references,
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
