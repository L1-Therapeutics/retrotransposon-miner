from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import pysam

_MINIMAP2_INDEX_CACHE: dict[str, Path] = {}
_MINIMAP2_INDEX_LOCK = threading.Lock()
_MEI_FASTA_LENGTH_CACHE: dict[str, dict[str, int]] = {}
_MEI_FASTA_LENGTH_CACHE_LOCK = threading.Lock()
_MIN_SIDE_ANCHOR_ALN_LEN = 30
_MIN_POLYA_RUN_FOR_FULL_3P_IMPUTE = 12
_ASSEMBLY_FEATURE_SCHEMA_VERSION = 4


def _safe_locus_id(chrom: str, start: int, end: int) -> str:
    chrom_safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(chrom))
    return f"{chrom_safe}_{int(start)}_{int(end)}"


def _window_locus_id_from_row(row: pd.Series) -> str:
    chrom = str(row.get("chrom", ""))
    window_start = int(row.get("window_start", 1))
    window_end = int(row.get("window_end", window_start))
    if window_end < window_start:
        window_start, window_end = window_end, window_start
    return _safe_locus_id(chrom, max(1, window_start), max(1, window_end))


def _interval_from_row(row: pd.Series, pad_bp: int) -> tuple[str, int, int]:
    chrom = str(row.get("chrom", ""))
    window_start = int(row.get("window_start", 1))
    window_end = int(row.get("window_end", window_start))
    center = int(row.get("insertion_breakpoint_pos", 0))
    if center <= 0:
        center = (window_start + window_end) // 2
    pad = max(1, int(pad_bp))
    start = min(window_start, center - pad)
    end = max(window_end, center + pad)
    return chrom, max(1, start), max(1, end)


def _canonical_read_name(name: str) -> str:
    n = str(name or "").strip()
    if not n:
        return ""
    if n.startswith("@"):
        n = n[1:]
    n = n.split()[0]
    if n.endswith("/1") or n.endswith("/2"):
        n = n[:-2]
    if re.search(r":[012]$", n):
        n = n[:-2]
    return n


def _preferred_read_lookup(preferred_read_names: set[str] | None) -> set[str]:
    if not preferred_read_names:
        return set()
    out: set[str] = set()
    for name in preferred_read_names:
        raw = str(name or "").strip()
        if not raw:
            continue
        out.add(raw)
        canon = _canonical_read_name(raw)
        if canon:
            out.add(canon)
    return out


def _is_primary_interval_read(read: pysam.AlignedSegment) -> bool:
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.query_sequence is None:
        return False
    return True


def _read_matches_preferred(read: pysam.AlignedSegment, preferred_lookup: set[str]) -> bool:
    if not preferred_lookup:
        return False
    qname = str(read.query_name or "")
    pair_suffix = "1" if read.is_read1 else "2" if read.is_read2 else "0"
    rid = f"{qname}:{pair_suffix}"
    candidates = {
        qname,
        rid,
        f"{qname}/1",
        f"{qname}/2",
        _canonical_read_name(qname),
        _canonical_read_name(rid),
    }
    return any(x and x in preferred_lookup for x in candidates)


def _write_interval_fastq(
    bam_path: Path,
    chrom: str,
    start_1based: int,
    end_1based: int,
    out_fastq_gz: Path,
    *,
    max_reads: int,
    non_perfect_only: bool = False,
    preferred_read_names: set[str] | None = None,
) -> tuple[int, set[str]]:
    start0 = max(0, int(start_1based) - 1)
    end0 = max(start0 + 1, int(end_1based))
    written = 0
    seen_ids: set[str] = set()
    emitted_qnames: set[str] = set()
    preferred_lookup = _preferred_read_lookup(preferred_read_names)
    out_fastq_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_fastq_gz, "wt", encoding="utf-8") as oh:
        def _emit_read(read: pysam.AlignedSegment) -> bool:
            nonlocal written
            rid = f"{read.query_name}:{'1' if read.is_read1 else '2' if read.is_read2 else '0'}"
            if rid in seen_ids:
                return False
            seq = str(read.query_sequence)
            qual = read.qual if read.qual is not None else ("I" * len(seq))
            oh.write(f"@{rid}\n{seq}\n+\n{qual}\n")
            seen_ids.add(rid)
            qname = str(read.query_name or "")
            if qname:
                emitted_qnames.add(qname)
            written += 1
            return True

        # Pass 0: force include locus-linked evidence read names first.
        if preferred_lookup:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam.fetch(chrom, start0, end0):
                    if not _is_primary_interval_read(read):
                        continue
                    if not _read_matches_preferred(read, preferred_lookup):
                        continue
                    _emit_read(read)
                    if written >= int(max_reads):
                        break

        # Pass 1: prioritize non-perfect/evidence-like reads.
        if written < int(max_reads):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam.fetch(chrom, start0, end0):
                    if not _is_primary_interval_read(read):
                        continue
                    if not _is_non_perfect_primary_alignment(read):
                        continue
                    _emit_read(read)
                    if written >= int(max_reads):
                        break

        # Pass 2: if full mode and under cap, backfill with clean/proper reads.
        if (not non_perfect_only) and written < int(max_reads):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam.fetch(chrom, start0, end0):
                    if not _is_primary_interval_read(read):
                        continue
                    if _is_non_perfect_primary_alignment(read):
                        continue
                    _emit_read(read)
                    if written >= int(max_reads):
                        break
    return written, emitted_qnames


def _read_recruited_qnames_from_fastq(path: Path) -> set[str]:
    if not path.exists():
        return set()
    out: set[str] = set()
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for idx, line in enumerate(handle):
                if idx % 4 != 0:
                    continue
                header = line.strip()
                if not header.startswith("@"):
                    continue
                rid = header[1:].split()[0]
                qname = rid.rsplit(":", 1)[0]
                qname = str(qname).strip()
                if not qname:
                    continue
                out.add(qname)
    except Exception:
        return set()
    return out


def _recruited_qnames_for_cache_hit(locus_dir: Path, sample: str, pad_bp: int, status: str) -> set[str]:
    sample_txt = str(sample)
    status_txt = str(status or "")
    if "highcap" in status_txt:
        path = locus_dir / f"{sample_txt}.pad{int(pad_bp)}.cap1000.fastq.gz"
        return _read_recruited_qnames_from_fastq(path)
    if status_txt.endswith("seed") or "_seed" in status_txt:
        path = locus_dir / f"{sample_txt}.pad{int(pad_bp)}.seedcap200.fastq.gz"
        return _read_recruited_qnames_from_fastq(path)
    path = locus_dir / f"{sample_txt}.pad{int(pad_bp)}.fastq.gz"
    return _read_recruited_qnames_from_fastq(path)


def _extract_microhomology_sequence(
    *,
    mei_hit: dict[str, object],
    left_hit: dict[str, object] | None,
    right_hit: dict[str, object] | None,
    seq_by_name: dict[str, str],
    max_len: int = 25,
    allow_homopolymer: bool = False,
) -> str:
    qname = str(mei_hit.get("qname", ""))
    seq = str(seq_by_name.get(qname, "") or "")
    if not qname or not seq:
        return ""
    m_qstart = int(mei_hit.get("qstart", 0))
    m_qend = int(mei_hit.get("qend", 0))
    if m_qstart < 0 or m_qend <= m_qstart:
        return ""

    candidates: list[str] = []
    if left_hit is not None and str(left_hit.get("qname", "")) == qname:
        left_qend = int(left_hit.get("qend", 0))
        ov = left_qend - m_qstart
        if ov > 0:
            start = max(0, m_qstart)
            end = min(len(seq), left_qend)
            if end > start:
                candidates.append(str(seq[start:end]).upper())
    if right_hit is not None and str(right_hit.get("qname", "")) == qname:
        right_qstart = int(right_hit.get("qstart", 0))
        ov = m_qend - right_qstart
        if ov > 0:
            start = max(0, right_qstart)
            end = min(len(seq), m_qend)
            if end > start:
                candidates.append(str(seq[start:end]).upper())
    if not candidates:
        return ""

    best = max(candidates, key=len)
    if len(best) > int(max_len):
        best = best[: int(max_len)]
    # Strict mode excludes pure homopolymer overlap.
    if (not bool(allow_homopolymer)) and len(set(best)) <= 1:
        return ""
    return best


def _is_non_perfect_primary_alignment(read: pysam.AlignedSegment) -> bool:
    # Non-perfect mapping heuristic for seed-first assembly:
    # soft clips/indels/skips, non-zero NM, low MAPQ, or pair-structure discordance.
    if read.mapping_quality < 60:
        return True
    if read.cigartuples:
        # 1=I, 2=D, 3=N, 4=S, 5=H
        if any(op in {1, 2, 3, 4, 5} and int(length) > 0 for op, length in read.cigartuples):
            return True
    try:
        if read.has_tag("NM") and int(read.get_tag("NM")) > 0:
            return True
    except Exception:
        pass
    if read.mate_is_unmapped or (not read.is_proper_pair):
        return True
    return False


def _run_spades(
    spades_exe: str,
    fastq_gz: Path,
    outdir: Path,
    *,
    threads: int,
    memory_gb: int,
) -> tuple[int, str]:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        spades_exe,
        "--only-assembler",
        "-s",
        str(fastq_gz),
        "-o",
        str(outdir),
        "-t",
        str(max(1, int(threads))),
        "-m",
        str(max(1, int(memory_gb))),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    stderr_tail = (proc.stderr or "")[-2000:]
    stdout_tail = (proc.stdout or "")[-2000:]
    combined = "\n".join([x for x in [stdout_tail, stderr_tail] if x]).strip()
    return int(proc.returncode), combined


def _resolve_spades_executable() -> str:
    for name in ("spades.py", "spades"):
        exe = shutil.which(name)
        if exe:
            return exe
    if "CONDA_PREFIX" in os.environ:
        p = Path(os.environ["CONDA_PREFIX"])
        for name in ("spades.py", "spades"):
            candidate = p / "bin" / name
            if candidate.exists():
                return str(candidate)
    raise FileNotFoundError(
        "SPAdes executable not found (spades.py/spades). Install with "
        "'micromamba install -n rtm-miner -c bioconda -c conda-forge spades' "
        "and ensure the environment is activated."
    )


def _summarize_contigs(contigs_fasta: Path) -> tuple[int, int]:
    if not contigs_fasta.exists():
        return 0, 0
    count = 0
    max_len = 0
    cur = 0
    with contigs_fasta.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if cur > 0:
                    max_len = max(max_len, cur)
                cur = 0
                count += 1
            else:
                cur += len(line.strip())
        if cur > 0:
            max_len = max(max_len, cur)
    return count, max_len


def _iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    out: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if name:
        out.append((name, "".join(seq_parts)))
    return out


def _run_minimap2_paf(
    query_fa: Path,
    target_fa: Path,
    *,
    preset: str = "asm5",
    threads: int = 1,
) -> list[dict[str, object]]:
    if not query_fa.exists() or not target_fa.exists():
        return []
    target_arg = str(target_fa)
    target_key = str(target_fa.resolve())
    with _MINIMAP2_INDEX_LOCK:
        cached_idx = _MINIMAP2_INDEX_CACHE.get(target_key)
        if cached_idx is not None and cached_idx.exists():
            target_arg = str(cached_idx)
        else:
            # Build once per process and reuse to avoid per-locus reference reindexing.
            preferred_idx = target_fa.with_suffix(target_fa.suffix + ".mmi")
            idx_path = preferred_idx
            if not preferred_idx.exists():
                idx_dir = Path(tempfile.gettempdir()) / "rtm_minimap2_indexes"
                idx_dir.mkdir(parents=True, exist_ok=True)
                digest = hashlib.sha1(target_key.encode("utf-8")).hexdigest()[:16]
                idx_path = idx_dir / f"{target_fa.stem}.{digest}.mmi"
            if not idx_path.exists():
                build = subprocess.run(
                    ["minimap2", "-d", str(idx_path), str(target_fa)],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if build.returncode == 0 and idx_path.exists():
                    _MINIMAP2_INDEX_CACHE[target_key] = idx_path
                    target_arg = str(idx_path)
            else:
                _MINIMAP2_INDEX_CACHE[target_key] = idx_path
                target_arg = str(idx_path)
    cmd = ["minimap2", "-x", preset, "--secondary=no", "-c", "-t", str(max(1, int(threads))), target_arg, str(query_fa)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        return []
    rows: list[dict[str, object]] = []
    for line in (proc.stdout or "").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        try:
            qname = parts[0]
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tstart = int(parts[7])
            tend = int(parts[8])
            nmatch = int(parts[9])
            alnlen = int(parts[10])
            mapq = int(parts[11])
            rows.append(
                {
                    "qname": qname,
                    "qstart": qstart,
                    "qend": qend,
                    "strand": strand,
                    "tname": tname,
                    "tstart": tstart,
                    "tend": tend,
                    "nmatch": nmatch,
                    "alnlen": alnlen,
                    "mapq": mapq,
                }
            )
        except Exception:
            continue
    return rows


def _family_from_target(target: str) -> str:
    t = (target or "").upper()
    if "ALU" in t:
        return "ALU"
    if "SVA" in t:
        return "SVA"
    if "LINE1" in t or "L1" in t:
        return "LINE1"
    if "HERV" in t or "ERV" in t:
        return "ERV"
    return ""


def _poly_at_max_run(seq: str) -> int:
    s = (seq or "").upper()
    best = 0
    cur = 0
    prev = ""
    for ch in s:
        if ch not in {"A", "T"}:
            cur = 0
            prev = ""
            continue
        if ch == prev:
            cur += 1
        else:
            cur = 1
            prev = ch
        if cur > best:
            best = cur
    return best


def _load_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    key = str(fasta_path.resolve())
    with _MEI_FASTA_LENGTH_CACHE_LOCK:
        cached = _MEI_FASTA_LENGTH_CACHE.get(key)
        if cached is not None:
            return cached
    lengths: dict[str, int] = {}
    try:
        with pysam.FastaFile(str(fasta_path)) as fa:
            for name, length in zip(fa.references, fa.lengths):
                lengths[str(name)] = int(length)
    except Exception:
        lengths = {}
    with _MEI_FASTA_LENGTH_CACHE_LOCK:
        _MEI_FASTA_LENGTH_CACHE[key] = lengths
    return lengths


def _sort_hit_key(hit: dict[str, object]) -> tuple[int, int, int]:
    return (int(hit.get("alnlen", 0)), int(hit.get("mapq", 0)), int(hit.get("nmatch", 0)))


def _polyA_impute_threshold_for_family(mei_family: str) -> int:
    fam = str(mei_family or "").upper()
    if fam == "ALU":
        return 15
    if fam in {"LINE1", "SVA"}:
        return 12
    return int(_MIN_POLYA_RUN_FOR_FULL_3P_IMPUTE)


def _infer_breakpoint_from_alignments(
    mei_hit: dict[str, object],
    ref_hits: list[dict[str, object]],
    *,
    default_pos: int,
) -> tuple[int, str, str, int]:
    if not mei_hit:
        return int(default_pos), "", "", 0
    contig = str(mei_hit.get("qname", ""))
    m_qstart = int(mei_hit.get("qstart", 0))
    m_qend = int(mei_hit.get("qend", 0))
    rh = [h for h in ref_hits if str(h.get("qname", "")) == contig]
    if not rh:
        return int(default_pos), "", "", 0

    left = [h for h in rh if int(h.get("qend", 0)) <= m_qstart + 20]
    right = [h for h in rh if int(h.get("qstart", 0)) >= m_qend - 20]
    left_hit = max(left, key=lambda h: int(h.get("qend", 0)), default=None)
    right_hit = min(right, key=lambda h: int(h.get("qstart", 0)), default=None)

    chrom = ""
    left_bp = 0
    right_bp = 0
    if left_hit is not None:
        chrom = str(left_hit.get("tname", ""))
        left_bp = int(left_hit.get("tend", 0))  # 1-based endpoint
    if right_hit is not None:
        chrom = chrom or str(right_hit.get("tname", ""))
        right_bp = int(right_hit.get("tstart", 0)) + 1  # 1-based start

    if left_bp > 0 and right_bp > 0 and chrom and str(left_hit.get("tname", "")) == str(right_hit.get("tname", "")):
        bp = int((left_bp + right_bp) // 2)
        tsd_len = int(right_bp - left_bp + 1) if right_bp >= left_bp else 0
        return bp, chrom, chrom, tsd_len
    if left_bp > 0:
        return left_bp, chrom, chrom, 0
    if right_bp > 0:
        return right_bp, chrom, chrom, 0
    return int(default_pos), "", "", 0


def _flank_ref_hits_for_mei(
    mei_hit: dict[str, object], ref_hits: list[dict[str, object]]
) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    if not mei_hit:
        return None, None
    contig = str(mei_hit.get("qname", ""))
    m_qstart = int(mei_hit.get("qstart", 0))
    m_qend = int(mei_hit.get("qend", 0))
    rh = [h for h in ref_hits if str(h.get("qname", "")) == contig]
    if not rh:
        return None, None
    left = [h for h in rh if int(h.get("qend", 0)) <= m_qstart + 20]
    right = [h for h in rh if int(h.get("qstart", 0)) >= m_qend - 20]
    left_hit = max(left, key=lambda h: int(h.get("qend", 0)), default=None)
    right_hit = min(right, key=lambda h: int(h.get("qstart", 0)), default=None)
    return left_hit, right_hit


def _breakpoint_side_status(mei_hit: dict[str, object], ref_hits: list[dict[str, object]]) -> str:
    left_hit, right_hit = _flank_ref_hits_for_mei(mei_hit, ref_hits)
    has_left = left_hit is not None
    has_right = right_hit is not None
    if has_left and has_right:
        return "both"
    if has_left:
        return "left_only"
    if has_right:
        return "right_only"
    return "none"


def _infer_non_mei_partner(
    *,
    left_hit: dict[str, object] | None,
    right_hit: dict[str, object] | None,
    side_status: str,
) -> tuple[str, int, str]:
    if side_status == "left_only" and left_hit is not None:
        return str(left_hit.get("tname", "")), int(left_hit.get("tend", 0)), "single_flank_ref_only"
    if side_status == "right_only" and right_hit is not None:
        return str(right_hit.get("tname", "")), int(right_hit.get("tstart", 0)) + 1, "single_flank_ref_only"
    if side_status != "both" or left_hit is None or right_hit is None:
        return "", -1, ""

    lchrom = str(left_hit.get("tname", ""))
    rchrom = str(right_hit.get("tname", ""))
    lstrand = str(left_hit.get("strand", ""))
    rstrand = str(right_hit.get("strand", ""))
    lpos = int(left_hit.get("tend", 0))
    rpos = int(right_hit.get("tstart", 0)) + 1
    if lchrom and rchrom and lchrom != rchrom:
        return rchrom, rpos, "interchrom_breakend"
    if lstrand in {"+", "-"} and rstrand in {"+", "-"} and lstrand != rstrand:
        return rchrom or lchrom, rpos, "inversion_breakend"
    if lchrom and rchrom and lchrom == rchrom and abs(rpos - lpos) > 1000:
        return rchrom, rpos, "large_gap_breakend"
    return "", -1, ""


def _summarize_complex_topology(
    records: list[tuple[str, str]],
    mei_hits: list[dict[str, object]],
    ref_hits: list[dict[str, object]],
    *,
    top_k: int = 3,
) -> dict[str, object]:
    seq_by_name = {name: seq for name, seq in records}
    mei_by_contig: dict[str, list[dict[str, object]]] = {}
    ref_by_contig: dict[str, list[dict[str, object]]] = {}
    for hit in mei_hits:
        mei_by_contig.setdefault(str(hit.get("qname", "")), []).append(hit)
    for hit in ref_hits:
        ref_by_contig.setdefault(str(hit.get("qname", "")), []).append(hit)

    scored: list[tuple[int, int, str]] = []
    for contig in set(list(seq_by_name.keys()) + list(mei_by_contig.keys()) + list(ref_by_contig.keys())):
        mei_best = max((int(h.get("alnlen", 0)) for h in mei_by_contig.get(contig, [])), default=0)
        ref_best = max((int(h.get("alnlen", 0)) for h in ref_by_contig.get(contig, [])), default=0)
        contig_len = len(seq_by_name.get(contig, ""))
        score = mei_best + ref_best
        scored.append((score, contig_len, contig))
    scored.sort(reverse=True)
    keep = [contig for _score, _len, contig in scored[: max(1, int(top_k))]]

    class_rank = {"unknown": 0, "simple_mei": 1, "mei_plus_sv": 2, "multi_junction": 3}
    side_rank = {"none": 0, "left_only": 1, "right_only": 1, "both": 2}
    best_class = "unknown"
    best_side_status = "none"
    partner_chrom = ""
    partner_pos = -1
    partner_type = ""
    best_sort_key = (-1, -1, -1, -1)

    for contig in keep:
        c_mei = sorted(
            mei_by_contig.get(contig, []),
            key=lambda h: (int(h.get("alnlen", 0)), int(h.get("mapq", 0)), int(h.get("nmatch", 0))),
            reverse=True,
        )
        c_ref = sorted(
            ref_by_contig.get(contig, []),
            key=lambda h: (int(h.get("alnlen", 0)), int(h.get("mapq", 0)), int(h.get("nmatch", 0))),
            reverse=True,
        )
        has_mei = len(c_mei) > 0
        strong_ref = [h for h in c_ref if int(h.get("alnlen", 0)) >= 80 and int(h.get("mapq", 0)) >= 20]
        ref_chroms = {str(h.get("tname", "")) for h in strong_ref if str(h.get("tname", ""))}
        ref_strands = {str(h.get("strand", "")) for h in strong_ref if str(h.get("strand", "")) in {"+", "-"}}
        ref_blocks = len(strong_ref)
        side_status = "none"
        this_partner_chrom = ""
        this_partner_pos = -1
        this_partner_type = ""
        this_class = "unknown"
        if has_mei:
            side_status = _breakpoint_side_status(c_mei[0], c_ref)
            left_hit, right_hit = _flank_ref_hits_for_mei(c_mei[0], c_ref)
            this_partner_chrom, this_partner_pos, this_partner_type = _infer_non_mei_partner(
                left_hit=left_hit,
                right_hit=right_hit,
                side_status=side_status,
            )
        if has_mei and ref_blocks >= 3:
            this_class = "multi_junction"
        elif has_mei and (this_partner_type or side_status in {"left_only", "right_only"}):
            this_class = "mei_plus_sv"
        elif has_mei and (len(ref_chroms) > 1 or len(ref_strands) > 1):
            this_class = "mei_plus_sv"
        elif has_mei:
            this_class = "simple_mei"
        elif ref_blocks >= 3 or len(ref_chroms) > 1 or len(ref_strands) > 1:
            this_class = "multi_junction"
        contig_mei_aln = int(c_mei[0].get("alnlen", 0)) if has_mei else 0
        sort_key = (
            class_rank.get(this_class, 0),
            side_rank.get(side_status, 0),
            1 if this_partner_type else 0,
            contig_mei_aln,
        )
        if sort_key > best_sort_key:
            best_sort_key = sort_key
            best_class = this_class
            best_side_status = side_status
            partner_chrom = this_partner_chrom
            partner_pos = int(this_partner_pos)
            partner_type = this_partner_type

    return {
        "complex_class": best_class,
        "breakpoint_side_status": best_side_status,
        "non_mei_partner_chrom": partner_chrom,
        "non_mei_partner_pos": int(partner_pos),
        "non_mei_partner_type": partner_type,
        "top_contigs": ",".join(keep),
    }


def _sample_complexity_key(feat: dict[str, object]) -> tuple[int, int, int, int]:
    class_rank = {"unknown": 0, "simple_mei": 1, "mei_plus_sv": 2, "multi_junction": 3}
    side_rank = {"none": 0, "left_only": 1, "right_only": 1, "both": 2}
    return (
        class_rank.get(str(feat.get("complex_class", "unknown")), 0),
        side_rank.get(str(feat.get("breakpoint_side_status", "none")), 0),
        1 if str(feat.get("non_mei_partner_type", "")) else 0,
        int(feat.get("mei_aln_len", 0)),
    )


def _choose_consensus_features(
    d_feat: dict[str, object], n_feat: dict[str, object]
) -> tuple[dict[str, object], str, dict[str, object], str]:
    def _coord_ok(feat: dict[str, object]) -> bool:
        s = int(feat.get("insertion_mei_start", -1))
        e = int(feat.get("insertion_mei_end", -1))
        return s > 0 and e > 0 and s >= e

    def _pick_score(feat: dict[str, object]) -> tuple[int, int, int, int]:
        left_aln = int(feat.get("left_support_mei_aln_len", 0))
        right_aln = int(feat.get("right_support_mei_aln_len", 0))
        mei_len = int(feat.get("mei_aln_len", 0))
        side_anchors = int(left_aln >= _MIN_SIDE_ANCHOR_ALN_LEN) + int(right_aln >= _MIN_SIDE_ANCHOR_ALN_LEN)
        coord_model = str(feat.get("coord_model", ""))
        is_informative_model = 0 if coord_model in {"", "no_contigs", "no_mei_hits"} else 1
        return (
            1 if _coord_ok(feat) else 0,
            side_anchors,
            max(left_aln, right_aln, mei_len),
            is_informative_model,
        )

    d_score = _pick_score(d_feat)
    n_score = _pick_score(n_feat)
    if d_score > (0, 0, 0, 0) and (n_score <= (0, 0, 0, 0) or d_score >= n_score):
        pick = d_feat
        pick_source = "disease"
    elif n_score > (0, 0, 0, 0):
        pick = n_feat
        pick_source = "control"
    else:
        pick = {}
        pick_source = ""
    if _sample_complexity_key(d_feat) >= _sample_complexity_key(n_feat):
        complex_pick = d_feat
        complex_source = "disease"
    else:
        complex_pick = n_feat
        complex_source = "control"
    return pick, pick_source, complex_pick, complex_source


def _should_escalate_read_cap(
    *,
    status: str,
    d_reads: int,
    n_reads: int,
    max_reads_per_sample: int,
    pick_source: str,
    d_cc: int,
    n_cc: int,
    escalated_max_reads: int,
) -> bool:
    if not str(status).startswith("assembled"):
        return False
    if int(max_reads_per_sample) >= int(escalated_max_reads):
        return False
    if int(d_reads) < int(max_reads_per_sample) or int(n_reads) < int(max_reads_per_sample):
        return False
    no_mei_source = str(pick_source or "") == ""
    low_complexity = int(d_cc) <= 1 and int(n_cc) <= 1
    return bool(no_mei_source or low_complexity)


def _extract_sample_assembly_features(
    *,
    contigs_fasta: Path,
    mei_fasta: Path,
    reference_fasta: Path | None,
    default_breakpoint_pos: int,
    minimap2_threads: int = 1,
) -> dict[str, object]:
    records = _iter_fasta_records(contigs_fasta)
    seq_by_name = {name: seq for name, seq in records}
    if not records:
        return {
            "primary_contig_id": "",
            "mei_subfamily": "",
            "mei_family": "",
            "mei_target_length": 0,
            "insertion_mei_start": -1,
            "insertion_mei_end": -1,
            "insertion_orientation": "",
            "insertion_length": 0,
            "insertion_length_observed": 0,
            "insertion_length_imputed": 0,
            "insertion_length_confidence_tier": "undetermined",
            "polyA_max_run": 0,
            "breakpoint_pos": int(default_breakpoint_pos),
            "breakpoint_chrom": "",
            "tsd_len": 0,
            "tsd_seq": "",
            "mei_aln_len": 0,
            "complex_class": "unknown",
            "non_mei_partner_chrom": "",
            "non_mei_partner_pos": -1,
            "non_mei_partner_type": "",
            "breakpoint_side_status": "none",
            "top_contigs": "",
            "mei_alignment_preset": "",
            "left_support_contig_id": "",
            "right_support_contig_id": "",
            "left_support_mei_start": -1,
            "left_support_mei_end": -1,
            "right_support_mei_start": -1,
            "right_support_mei_end": -1,
            "left_support_mei_aln_len": 0,
            "right_support_mei_aln_len": 0,
            "microhomology_sequence": "",
            "junction_overlap_sequence": "",
            "coord_model": "no_contigs",
            "coord_logic_version": int(_ASSEMBLY_FEATURE_SCHEMA_VERSION),
        }

    with tempfile.TemporaryDirectory(prefix="rtm_asm_parse_") as td:
        query_fa = Path(td) / "contigs.fa"
        with query_fa.open("w", encoding="utf-8") as oh:
            for name, seq in records:
                oh.write(f">{name}\n{seq}\n")
        mei_preset = "asm5"
        mei_hits = _run_minimap2_paf(query_fa, mei_fasta, preset=mei_preset, threads=minimap2_threads)
        if not mei_hits:
            # Contigs can be short/chimeric around breakpoints; sr preset is more
            # sensitive for recovering partial MEI-supporting alignments.
            mei_preset = "sr"
            mei_hits = _run_minimap2_paf(query_fa, mei_fasta, preset=mei_preset, threads=minimap2_threads)
        if reference_fasta is not None:
            ref_hits = _run_minimap2_paf(query_fa, reference_fasta, preset="asm5", threads=minimap2_threads)
        else:
            ref_hits = []

    if not mei_hits:
        poly = max((_poly_at_max_run(seq) for _n, seq in records), default=0)
        return {
            "primary_contig_id": "",
            "mei_subfamily": "",
            "mei_family": "",
            "mei_target_length": 0,
            "insertion_mei_start": -1,
            "insertion_mei_end": -1,
            "insertion_orientation": "",
            "insertion_length": 0,
            "insertion_length_observed": 0,
            "insertion_length_imputed": 0,
            "insertion_length_confidence_tier": "undetermined",
            "polyA_max_run": int(poly),
            "breakpoint_pos": int(default_breakpoint_pos),
            "breakpoint_chrom": "",
            "tsd_len": 0,
            "tsd_seq": "",
            "mei_aln_len": 0,
            "complex_class": "unknown",
            "non_mei_partner_chrom": "",
            "non_mei_partner_pos": -1,
            "non_mei_partner_type": "",
            "breakpoint_side_status": "none",
            "top_contigs": "",
            "mei_alignment_preset": mei_preset,
            "left_support_contig_id": "",
            "right_support_contig_id": "",
            "left_support_mei_start": -1,
            "left_support_mei_end": -1,
            "right_support_mei_start": -1,
            "right_support_mei_end": -1,
            "left_support_mei_aln_len": 0,
            "right_support_mei_aln_len": 0,
            "microhomology_sequence": "",
            "junction_overlap_sequence": "",
            "coord_model": "no_mei_hits",
            "coord_logic_version": int(_ASSEMBLY_FEATURE_SCHEMA_VERSION),
        }

    best = max(mei_hits, key=lambda h: (int(h["alnlen"]), int(h["mapq"]), int(h["nmatch"])))
    best_lo = int(best["tstart"]) + 1
    best_hi = int(best["tend"])
    mei_strand = str(best["strand"]) if str(best["strand"]) in {"+", "-"} else ""
    left_cands: list[dict[str, object]] = []
    right_cands: list[dict[str, object]] = []
    for hit in mei_hits:
        if int(hit.get("alnlen", 0)) < int(_MIN_SIDE_ANCHOR_ALN_LEN):
            continue
        side_status = _breakpoint_side_status(hit, ref_hits)
        if side_status in {"left_only", "both"}:
            left_cands.append(hit)
        if side_status in {"right_only", "both"}:
            right_cands.append(hit)
    left_hit = max(left_cands, key=_sort_hit_key, default=None)
    right_hit = max(right_cands, key=_sort_hit_key, default=None)

    left_lo = int(left_hit.get("tstart", -1)) + 1 if left_hit is not None else -1
    left_hi = int(left_hit.get("tend", -1)) if left_hit is not None else -1
    right_lo = int(right_hit.get("tstart", -1)) + 1 if right_hit is not None else -1
    right_hi = int(right_hit.get("tend", -1)) if right_hit is not None else -1

    poly_max = max((_poly_at_max_run(seq) for _n, seq in records), default=0)
    mei_lengths = _load_fasta_lengths(mei_fasta)
    target_name = str(best.get("tname", ""))
    target_len = int(mei_lengths.get(target_name, 0))
    target_family = _family_from_target(target_name)
    polyA_impute_threshold = _polyA_impute_threshold_for_family(target_family)
    allow_left_anchor_impute = mei_strand == "+"
    allow_right_anchor_impute = mei_strand == "-"
    coord_model = "single_best_contig"
    cx = _summarize_complex_topology(records, mei_hits, ref_hits, top_k=3)

    if left_hit is not None and right_hit is not None:
        fused_lo = min(left_lo, right_lo)
        fused_hi = max(left_hi, right_hi)
        if str(left_hit.get("qname", "")) == str(right_hit.get("qname", "")):
            coord_model = "single_contig_both_sides"
        else:
            coord_model = "multi_contig_fused_sides"
    elif left_hit is not None:
        fused_lo = left_lo
        fused_hi = left_hi
        coord_model = "left_anchor_only"
    elif right_hit is not None:
        fused_lo = right_lo
        fused_hi = right_hi
        coord_model = "right_anchor_only"
    else:
        fused_lo = best_lo
        fused_hi = best_hi
        coord_model = "single_best_no_side_anchor"

    observed_start = int(max(fused_lo, fused_hi))
    observed_end = int(min(fused_lo, fused_hi))
    observed_len = int(max(0, observed_start - observed_end + 1))

    imputation_applied = False
    imputation_reason = ""
    complex_class = str(cx.get("complex_class", "unknown"))
    complex_strict_poly_gate = poly_max >= int(max(polyA_impute_threshold, 18))
    complex_anchor_gate = (allow_left_anchor_impute and left_hit is not None) or (
        allow_right_anchor_impute and right_hit is not None
    )

    # Existing one-sided anchor + tail full-length imputation remains enabled.
    if coord_model == "left_anchor_only" and allow_left_anchor_impute and poly_max >= int(polyA_impute_threshold) and target_len > 0:
        fused_hi = max(fused_hi, target_len)
        coord_model = "left_anchor_plus_polyA_full_3p_impute"
        imputation_applied = True
        imputation_reason = "one_sided_polyA"
    elif (
        coord_model == "right_anchor_only"
        and allow_right_anchor_impute
        and poly_max >= int(polyA_impute_threshold)
        and target_len > 0
    ):
        fused_hi = max(fused_hi, target_len)
        coord_model = "right_anchor_plus_polyA_full_3p_impute"
        imputation_applied = True
        imputation_reason = "one_sided_polyA"
    # Strict complex-event imputation (explicitly without end-proximity gating).
    elif (
        complex_class == "mei_plus_sv"
        and target_len > 0
        and complex_strict_poly_gate
        and complex_anchor_gate
    ):
        fused_hi = max(fused_hi, target_len)
        coord_model = "complex_plus_polyA_full_3p_impute"
        imputation_applied = True
        imputation_reason = "complex_mei_plus_sv_polyA"

    # MEI target coordinates are on the consensus reference axis. Keep project
    # convention as 3'->5': start=high, end=low.
    ins_start = int(max(fused_lo, fused_hi))
    ins_end = int(min(fused_lo, fused_hi))
    ins_len = int(max(0, ins_start - ins_end + 1))
    imputed_len = int(max(ins_len, observed_len))
    breakpoint_resolved_model = coord_model in {
        "single_contig_both_sides",
        "multi_contig_fused_sides",
    }
    if target_len <= 0 or observed_len <= 0:
        length_confidence_tier = "undetermined"
    elif observed_len >= max(1, int(target_len * 0.9)):
        # A full-length MEI span on a single best contig without resolved two-sided
        # breakpoint geometry is informative but not junction-resolved confidence.
        # Reserve "high_conf" for two-sided breakpoint-resolved models.
        if breakpoint_resolved_model:
            length_confidence_tier = "full_length_high_conf"
        else:
            length_confidence_tier = "full_length_possible"
    elif imputation_applied and imputed_len >= target_len:
        if (
            imputation_reason == "complex_mei_plus_sv_polyA"
            and poly_max >= int(max(polyA_impute_threshold, 24))
            and breakpoint_resolved_model
        ):
            length_confidence_tier = "full_length_high_conf"
        else:
            length_confidence_tier = "full_length_possible"
    else:
        is_3p_truncated = bool(target_len > 0 and int(ins_start) < int(target_len))
        is_5p_truncated = bool(int(ins_end) > 1)
        if is_3p_truncated and is_5p_truncated:
            length_confidence_tier = "truncated_5p_3p"
        elif is_5p_truncated:
            length_confidence_tier = "truncated_5p"
        elif is_3p_truncated:
            length_confidence_tier = "truncated_3p"
        else:
            length_confidence_tier = "truncated"
    bp, bp_left_chrom, _bp_right_chrom, tsd_len = _infer_breakpoint_from_alignments(
        best, ref_hits, default_pos=default_breakpoint_pos
    )
    # Rescue missed TSDs with strong side-anchor flanks chosen for coordinate model
    # when the primary per-contig breakpoint inference does not yield a TSD.
    if tsd_len < 4 and left_hit is not None and right_hit is not None:
        try:
            lchrom = str(left_hit.get("tname", ""))
            rchrom = str(right_hit.get("tname", ""))
            if lchrom and rchrom and lchrom == rchrom:
                left_bp = int(left_hit.get("tend", 0))  # 1-based endpoint
                right_bp = int(right_hit.get("tstart", 0)) + 1  # 1-based start
                alt_tsd_len = int(right_bp - left_bp + 1) if right_bp >= left_bp else 0
                if alt_tsd_len >= 4:
                    bp = int((left_bp + right_bp) // 2)
                    bp_left_chrom = lchrom
                    tsd_len = alt_tsd_len
        except Exception:
            pass
    tsd_seq = ""
    if reference_fasta is not None and bp_left_chrom and tsd_len > 0 and tsd_len <= 50:
        try:
            with pysam.FastaFile(str(reference_fasta)) as ref:
                start0 = max(0, int(bp) - (tsd_len // 2) - 1)
                end0 = start0 + int(tsd_len)
                tsd_seq = ref.fetch(bp_left_chrom, start0, end0).upper()
        except Exception:
            tsd_seq = ""
    microhomology_sequence = _extract_microhomology_sequence(
        mei_hit=best,
        left_hit=left_hit,
        right_hit=right_hit,
        seq_by_name=seq_by_name,
        allow_homopolymer=False,
    )
    junction_overlap_sequence = _extract_microhomology_sequence(
        mei_hit=best,
        left_hit=left_hit,
        right_hit=right_hit,
        seq_by_name=seq_by_name,
        allow_homopolymer=True,
    )
    if (
        not microhomology_sequence
        and reference_fasta is not None
        and left_hit is not None
        and right_hit is not None
    ):
        try:
            lchrom = str(left_hit.get("tname", ""))
            rchrom = str(right_hit.get("tname", ""))
            if lchrom and rchrom and lchrom == rchrom:
                left_bp = int(left_hit.get("tend", 0))  # 1-based endpoint
                right_bp = int(right_hit.get("tstart", 0)) + 1  # 1-based start
                if right_bp <= left_bp:
                    mh_len = int(left_bp - right_bp + 1)
                    if 1 <= mh_len <= 25:
                        with pysam.FastaFile(str(reference_fasta)) as ref:
                            start0 = max(0, int(right_bp) - 1)
                            end0 = max(start0 + 1, int(left_bp))
                            mh_seq = ref.fetch(lchrom, start0, end0).upper()
                        if mh_seq and len(set(mh_seq)) > 1:
                            microhomology_sequence = str(mh_seq)
                        if mh_seq and not junction_overlap_sequence:
                            junction_overlap_sequence = str(mh_seq)
        except Exception:
            pass
    if not junction_overlap_sequence:
        qname = str(best.get("qname", ""))
        seq = str(seq_by_name.get(qname, "") or "")
        if seq:
            m_qstart = int(best.get("qstart", 0))
            m_qend = int(best.get("qend", 0))
            s0 = max(0, m_qstart - 8)
            e0 = min(len(seq), m_qstart + 17)
            seq_left = str(seq[s0:e0]).upper() if e0 > s0 else ""
            s1 = max(0, m_qend - 17)
            e1 = min(len(seq), m_qend + 8)
            seq_right = str(seq[s1:e1]).upper() if e1 > s1 else ""
            junction_overlap_sequence = seq_left if len(seq_left) >= len(seq_right) else seq_right
    return {
        "primary_contig_id": str(best["qname"]),
        "mei_subfamily": str(best["tname"]),
        "mei_family": _family_from_target(str(best["tname"])),
        "mei_target_length": int(target_len),
        "insertion_mei_start": int(ins_start),
        "insertion_mei_end": int(ins_end),
        "insertion_orientation": mei_strand,
        "insertion_length": int(imputed_len),
        "insertion_length_observed": int(observed_len),
        "insertion_length_imputed": int(imputed_len),
        "insertion_length_confidence_tier": str(length_confidence_tier),
        "polyA_max_run": int(poly_max),
        "breakpoint_pos": int(bp),
        "breakpoint_chrom": str(bp_left_chrom),
        "tsd_len": int(max(0, tsd_len)),
        "tsd_seq": str(tsd_seq),
        "mei_aln_len": int(best["alnlen"]),
        "complex_class": str(cx.get("complex_class", "simple_mei")),
        "non_mei_partner_chrom": str(cx.get("non_mei_partner_chrom", "")),
        "non_mei_partner_pos": int(cx.get("non_mei_partner_pos", -1)),
        "non_mei_partner_type": str(cx.get("non_mei_partner_type", "")),
        "breakpoint_side_status": str(cx.get("breakpoint_side_status", _breakpoint_side_status(best, ref_hits))),
        "top_contigs": str(cx.get("top_contigs", "")),
        "mei_alignment_preset": str(mei_preset),
        "left_support_contig_id": str(left_hit.get("qname", "") if left_hit is not None else ""),
        "right_support_contig_id": str(right_hit.get("qname", "") if right_hit is not None else ""),
        "left_support_mei_start": int(left_hi if left_hit is not None else -1),
        "left_support_mei_end": int(left_lo if left_hit is not None else -1),
        "right_support_mei_start": int(right_hi if right_hit is not None else -1),
        "right_support_mei_end": int(right_lo if right_hit is not None else -1),
        "left_support_mei_aln_len": int(left_hit.get("alnlen", 0) if left_hit is not None else 0),
        "right_support_mei_aln_len": int(right_hit.get("alnlen", 0) if right_hit is not None else 0),
        "microhomology_sequence": str(microhomology_sequence),
        "junction_overlap_sequence": str(junction_overlap_sequence),
        "coord_model": str(coord_model),
        "coord_logic_version": int(_ASSEMBLY_FEATURE_SCHEMA_VERSION),
    }


def _parse_existing_manifest(manifest_path: Path) -> dict[str, object] | None:
    if not manifest_path.exists():
        return None
    try:
        raw = json.loads(manifest_path.read_text(encoding="utf-8"))
        return raw if isinstance(raw, dict) else None
    except Exception:
        return None


def _has_assembled_cache(locus_dir: Path, manifest: dict[str, object] | None, *, interval_pad_bp: int) -> bool:
    if manifest is None:
        return False
    if not str(manifest.get("status", "")).startswith("assembled"):
        return False
    cached_pad = int(manifest.get("interval", {}).get("pad_bp", interval_pad_bp))
    d_out = locus_dir / f"disease.spades.pad{int(cached_pad)}"
    n_out = locus_dir / f"control.spades.pad{int(cached_pad)}"
    return bool((d_out / "contigs.fasta").exists() or (n_out / "contigs.fasta").exists())


def _process_single_locus(
    *,
    idx: int,
    total_loci: int,
    row_data: dict[str, object],
    disease_bam_path: Path,
    control_bam_path: Path,
    assembly_cache_dir: Path,
    mei_fasta: Path,
    reference_fasta: Path | None,
    interval_pad_bp: int,
    retry_pad_bp: int,
    max_reads_per_sample: int,
    spades_threads: int,
    spades_memory_gb: int,
    minimap2_threads: int,
    reuse_existing: bool,
    reuse_cache_only: bool,
    spades_exe: str | None,
    disease_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
    control_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
) -> tuple[dict[str, object], str]:
    def _has_sideaware_feature_schema(feat: object) -> bool:
        if not isinstance(feat, dict):
            return False
        required = {
            "coord_model",
            "left_support_contig_id",
            "right_support_contig_id",
            "left_support_mei_start",
            "left_support_mei_end",
            "right_support_mei_start",
            "right_support_mei_end",
            "left_support_mei_aln_len",
            "right_support_mei_aln_len",
            "mei_target_length",
            "insertion_length_observed",
            "insertion_length_imputed",
            "insertion_length_confidence_tier",
            "microhomology_sequence",
            "junction_overlap_sequence",
        }
        if not required.issubset(set(feat.keys())):
            return False
        version = int(feat.get("coord_logic_version", 0) or 0)
        return version >= int(_ASSEMBLY_FEATURE_SCHEMA_VERSION)

    as_row = pd.Series(row_data)
    locus_key = (
        str(as_row.get("chrom", "")),
        int(as_row.get("window_start", 0)),
        int(as_row.get("window_end", 0)),
    )
    disease_preferred_read_names = (
        disease_preferred_read_names_by_locus.get(locus_key, set()) if disease_preferred_read_names_by_locus else set()
    )
    control_preferred_read_names = (
        control_preferred_read_names_by_locus.get(locus_key, set()) if control_preferred_read_names_by_locus else set()
    )
    fallback_bp = int(as_row.get("insertion_breakpoint_pos", 0) or 0)
    chrom, i_start, i_end = _interval_from_row(as_row, interval_pad_bp)
    stable_locus_id = _window_locus_id_from_row(as_row)
    stable_locus_dir = assembly_cache_dir / stable_locus_id
    legacy_locus_dir = assembly_cache_dir / _safe_locus_id(chrom, i_start, i_end)
    locus_dir = stable_locus_dir
    manifest = None
    if reuse_existing:
        stable_manifest = _parse_existing_manifest(stable_locus_dir / "assembly_manifest.json")
        legacy_manifest = (
            _parse_existing_manifest(legacy_locus_dir / "assembly_manifest.json") if legacy_locus_dir != stable_locus_dir else None
        )
        if _has_assembled_cache(stable_locus_dir, stable_manifest, interval_pad_bp=interval_pad_bp):
            locus_dir = stable_locus_dir
            manifest = stable_manifest
        elif _has_assembled_cache(legacy_locus_dir, legacy_manifest, interval_pad_bp=interval_pad_bp):
            locus_dir = legacy_locus_dir
            manifest = legacy_manifest
    locus_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = locus_dir / "assembly_manifest.json"

    status = "assembled"
    retry_used = False
    err_msg = ""
    current_start = i_start
    current_end = i_end
    d_reads = 0
    n_reads = 0
    d_recruited_evidence_read_names: set[str] = set()
    n_recruited_evidence_read_names: set[str] = set()
    d_cc = 0
    n_cc = 0
    d_max = 0
    n_max = 0
    chosen_pad = int(interval_pad_bp)
    escalated_max_reads = 1000
    cache_hit = False

    if reuse_existing and _has_assembled_cache(locus_dir, manifest, interval_pad_bp=interval_pad_bp):
        cached_pad = int(manifest.get("interval", {}).get("pad_bp", interval_pad_bp))
        chosen_pad = cached_pad
        chrom, current_start, current_end = _interval_from_row(as_row, cached_pad)
        status = str(manifest.get("status", "assembled"))
        retry_used = bool(manifest.get("retry_used", False))
        d_reads = int(manifest.get("disease_reads_extracted", 0))
        n_reads = int(manifest.get("control_reads_extracted", 0))
        d_cc = int(manifest.get("disease_contigs", 0))
        n_cc = int(manifest.get("control_contigs", 0))
        d_max = int(manifest.get("disease_max_contig_len", 0))
        n_max = int(manifest.get("control_max_contig_len", 0))
        d_names_cached = manifest.get("disease_recruited_evidence_read_names", [])
        n_names_cached = manifest.get("control_recruited_evidence_read_names", [])
        if isinstance(d_names_cached, list):
            d_recruited_evidence_read_names = {str(x).strip() for x in d_names_cached if str(x).strip()}
        if isinstance(n_names_cached, list):
            n_recruited_evidence_read_names = {str(x).strip() for x in n_names_cached if str(x).strip()}
        if not d_recruited_evidence_read_names:
            d_recruited_evidence_read_names = _recruited_qnames_for_cache_hit(
                locus_dir=locus_dir,
                sample="disease",
                pad_bp=chosen_pad,
                status=status,
            )
        if not n_recruited_evidence_read_names:
            n_recruited_evidence_read_names = _recruited_qnames_for_cache_hit(
                locus_dir=locus_dir,
                sample="control",
                pad_bp=chosen_pad,
                status=status,
            )
        cache_hit = True
    elif manifest is None:
        manifest = _parse_existing_manifest(manifest_path)

    if manifest is None and reuse_cache_only:
        status = "assembly_cache_missing"
        err_msg = "reuse-cache-only enabled and no assembled cache found"
    elif manifest is None:
        seed_reads_cap = max(1, min(int(max_reads_per_sample), 200))
        assembled_from_attempt = False
        for attempt, pad in enumerate([interval_pad_bp, retry_pad_bp], start=1):
            chosen_pad = int(pad)
            chrom, current_start, current_end = _interval_from_row(as_row, pad)
            phase_specs: list[tuple[str, int, bool]] = [("seed", seed_reads_cap, True)]
            if int(seed_reads_cap) < int(max_reads_per_sample):
                phase_specs.append(("full", int(max_reads_per_sample), False))
            else:
                phase_specs = [("full", int(max_reads_per_sample), False)]
            for phase_name, phase_cap, non_perfect_only in phase_specs:
                suffix = "" if phase_name == "full" else f".seedcap{int(phase_cap)}"
                d_fq = locus_dir / f"disease.pad{int(pad)}{suffix}.fastq.gz"
                n_fq = locus_dir / f"control.pad{int(pad)}{suffix}.fastq.gz"
                d_reads, d_recruited_names_phase = _write_interval_fastq(
                    disease_bam_path,
                    chrom,
                    current_start,
                    current_end,
                    d_fq,
                    max_reads=int(phase_cap),
                    non_perfect_only=non_perfect_only,
                    preferred_read_names=disease_preferred_read_names,
                )
                n_reads, n_recruited_names_phase = _write_interval_fastq(
                    control_bam_path,
                    chrom,
                    current_start,
                    current_end,
                    n_fq,
                    max_reads=int(phase_cap),
                    non_perfect_only=non_perfect_only,
                    preferred_read_names=control_preferred_read_names,
                )
                if d_reads == 0 and n_reads == 0:
                    status = "assembly_no_reads"
                    err_msg = "no reads extracted in interval"
                    continue

                d_out_canonical = locus_dir / f"disease.spades.pad{int(pad)}"
                n_out_canonical = locus_dir / f"control.spades.pad{int(pad)}"
                d_out_try = d_out_canonical if phase_name == "full" else locus_dir / f"disease.spades.pad{int(pad)}{suffix}"
                n_out_try = n_out_canonical if phase_name == "full" else locus_dir / f"control.spades.pad{int(pad)}{suffix}"
                run_spades_exe = spades_exe or _resolve_spades_executable()
                d_rc, d_log = _run_spades(
                    run_spades_exe, d_fq, d_out_try, threads=spades_threads, memory_gb=spades_memory_gb
                )
                n_rc, n_log = _run_spades(
                    run_spades_exe, n_fq, n_out_try, threads=spades_threads, memory_gb=spades_memory_gb
                )
                if d_rc != 0 or n_rc != 0:
                    status = "assembly_failed"
                    err_msg = "; ".join([x for x in [d_log, n_log] if x]).strip()
                    continue

                d_cc, d_max = _summarize_contigs(d_out_try / "contigs.fasta")
                n_cc, n_max = _summarize_contigs(n_out_try / "contigs.fasta")
                if d_cc == 0 and n_cc == 0:
                    status = "assembly_no_contigs"
                    err_msg = "spades produced no contigs"
                    continue

                if phase_name == "seed":
                    d_feat_seed = _extract_sample_assembly_features(
                        contigs_fasta=d_out_try / "contigs.fasta",
                        mei_fasta=mei_fasta,
                        reference_fasta=reference_fasta,
                        default_breakpoint_pos=fallback_bp,
                        minimap2_threads=minimap2_threads,
                    )
                    n_feat_seed = _extract_sample_assembly_features(
                        contigs_fasta=n_out_try / "contigs.fasta",
                        mei_fasta=mei_fasta,
                        reference_fasta=reference_fasta,
                        default_breakpoint_pos=fallback_bp,
                        minimap2_threads=minimap2_threads,
                    )
                    _pick_seed, pick_source_seed, _cx_pick_seed, _cx_source_seed = _choose_consensus_features(
                        d_feat_seed, n_feat_seed
                    )
                    if str(pick_source_seed or "") == "":
                        # Seed contigs were not MEI-informative; escalate to full-read pass.
                        continue
                    if d_out_canonical.exists():
                        shutil.rmtree(d_out_canonical, ignore_errors=True)
                    if n_out_canonical.exists():
                        shutil.rmtree(n_out_canonical, ignore_errors=True)
                    if d_out_try.exists():
                        shutil.move(str(d_out_try), str(d_out_canonical))
                    if n_out_try.exists():
                        shutil.move(str(n_out_try), str(n_out_canonical))
                    retry_used = attempt > 1
                    status = "assembled_retry_seed" if retry_used else "assembled_seed"
                    err_msg = ""
                    d_recruited_evidence_read_names = set(d_recruited_names_phase)
                    n_recruited_evidence_read_names = set(n_recruited_names_phase)
                    assembled_from_attempt = True
                    break

                retry_used = attempt > 1
                status = "assembled_retry" if retry_used else "assembled"
                err_msg = ""
                d_recruited_evidence_read_names = set(d_recruited_names_phase)
                n_recruited_evidence_read_names = set(n_recruited_names_phase)
                assembled_from_attempt = True
                break
            if assembled_from_attempt:
                break

    d_out = locus_dir / f"disease.spades.pad{int(chosen_pad)}"
    n_out = locus_dir / f"control.spades.pad{int(chosen_pad)}"
    d_feat: dict[str, object]
    n_feat: dict[str, object]
    manifest_d_feat = manifest.get("disease_features", {}) if isinstance(manifest, dict) else {}
    manifest_n_feat = manifest.get("control_features", {}) if isinstance(manifest, dict) else {}
    if (
        cache_hit
        and isinstance(manifest_d_feat, dict)
        and isinstance(manifest_n_feat, dict)
        and _has_sideaware_feature_schema(manifest_d_feat)
        and _has_sideaware_feature_schema(manifest_n_feat)
    ):
        # Reuse precomputed feature payloads from cache manifests to avoid
        # per-locus minimap2 remapping when iterating on downstream annotation.
        d_feat = manifest_d_feat
        n_feat = manifest_n_feat
    else:
        d_feat = _extract_sample_assembly_features(
            contigs_fasta=d_out / "contigs.fasta",
            mei_fasta=mei_fasta,
            reference_fasta=reference_fasta,
            default_breakpoint_pos=fallback_bp,
            minimap2_threads=minimap2_threads,
        )
        n_feat = _extract_sample_assembly_features(
            contigs_fasta=n_out / "contigs.fasta",
            mei_fasta=mei_fasta,
            reference_fasta=reference_fasta,
            default_breakpoint_pos=fallback_bp,
            minimap2_threads=minimap2_threads,
        )
    pick, pick_source, complex_pick, complex_source = _choose_consensus_features(d_feat, n_feat)

    adaptive_rerun_used = False
    if (not cache_hit) and _should_escalate_read_cap(
        status=status,
        d_reads=d_reads,
        n_reads=n_reads,
        max_reads_per_sample=max_reads_per_sample,
        pick_source=pick_source,
        d_cc=d_cc,
        n_cc=n_cc,
        escalated_max_reads=escalated_max_reads,
    ):
        d_fq_hi = locus_dir / f"disease.pad{int(chosen_pad)}.cap{int(escalated_max_reads)}.fastq.gz"
        n_fq_hi = locus_dir / f"control.pad{int(chosen_pad)}.cap{int(escalated_max_reads)}.fastq.gz"
        d_reads_hi, d_recruited_names_hi = _write_interval_fastq(
            disease_bam_path,
            chrom,
            current_start,
            current_end,
            d_fq_hi,
            max_reads=int(escalated_max_reads),
            preferred_read_names=disease_preferred_read_names,
        )
        n_reads_hi, n_recruited_names_hi = _write_interval_fastq(
            control_bam_path,
            chrom,
            current_start,
            current_end,
            n_fq_hi,
            max_reads=int(escalated_max_reads),
            preferred_read_names=control_preferred_read_names,
        )
        if d_reads_hi > 0 or n_reads_hi > 0:
            d_out_hi = locus_dir / f"disease.spades.pad{int(chosen_pad)}.cap{int(escalated_max_reads)}"
            n_out_hi = locus_dir / f"control.spades.pad{int(chosen_pad)}.cap{int(escalated_max_reads)}"
            run_spades_exe = spades_exe or _resolve_spades_executable()
            d_rc_hi, d_log_hi = _run_spades(
                run_spades_exe, d_fq_hi, d_out_hi, threads=spades_threads, memory_gb=spades_memory_gb
            )
            n_rc_hi, n_log_hi = _run_spades(
                run_spades_exe, n_fq_hi, n_out_hi, threads=spades_threads, memory_gb=spades_memory_gb
            )
            if d_rc_hi == 0 and n_rc_hi == 0:
                d_cc_hi, d_max_hi = _summarize_contigs(d_out_hi / "contigs.fasta")
                n_cc_hi, n_max_hi = _summarize_contigs(n_out_hi / "contigs.fasta")
                if d_cc_hi > 0 or n_cc_hi > 0:
                    d_feat_hi = _extract_sample_assembly_features(
                        contigs_fasta=d_out_hi / "contigs.fasta",
                        mei_fasta=mei_fasta,
                        reference_fasta=reference_fasta,
                        default_breakpoint_pos=fallback_bp,
                        minimap2_threads=minimap2_threads,
                    )
                    n_feat_hi = _extract_sample_assembly_features(
                        contigs_fasta=n_out_hi / "contigs.fasta",
                        mei_fasta=mei_fasta,
                        reference_fasta=reference_fasta,
                        default_breakpoint_pos=fallback_bp,
                        minimap2_threads=minimap2_threads,
                    )
                    pick_hi, pick_source_hi, complex_pick_hi, complex_source_hi = _choose_consensus_features(d_feat_hi, n_feat_hi)
                    improved = (
                        str(pick_source_hi or "") != str(pick_source or "")
                        or int(pick_hi.get("mei_aln_len", 0) if pick_hi else 0) > int(pick.get("mei_aln_len", 0) if pick else 0)
                        or _sample_complexity_key(complex_pick_hi) > _sample_complexity_key(complex_pick)
                    )
                    if improved:
                        adaptive_rerun_used = True
                        d_reads = int(d_reads_hi)
                        n_reads = int(n_reads_hi)
                        d_cc = int(d_cc_hi)
                        n_cc = int(n_cc_hi)
                        d_max = int(d_max_hi)
                        n_max = int(n_max_hi)
                        d_feat = d_feat_hi
                        n_feat = n_feat_hi
                        d_recruited_evidence_read_names = set(d_recruited_names_hi)
                        n_recruited_evidence_read_names = set(n_recruited_names_hi)
                        pick = pick_hi
                        pick_source = pick_source_hi
                        complex_pick = complex_pick_hi
                        complex_source = complex_source_hi
                        if d_out.exists():
                            shutil.rmtree(d_out, ignore_errors=True)
                        if n_out.exists():
                            shutil.rmtree(n_out, ignore_errors=True)
                        if d_out_hi.exists():
                            shutil.move(str(d_out_hi), str(d_out))
                        if n_out_hi.exists():
                            shutil.move(str(n_out_hi), str(n_out))
                        status = "assembled_retry_highcap" if retry_used else "assembled_highcap"
                        err_msg = ""
            else:
                err_msg = "; ".join([x for x in [err_msg, d_log_hi, n_log_hi] if x]).strip()

    manifest_payload = {
        "status": status,
        "interval": {"chrom": chrom, "start": current_start, "end": current_end, "pad_bp": int(chosen_pad)},
        "retry_used": retry_used,
        "disease_reads_extracted": d_reads,
        "control_reads_extracted": n_reads,
        "disease_contigs": d_cc,
        "control_contigs": n_cc,
        "disease_max_contig_len": d_max,
        "control_max_contig_len": n_max,
        "error_message": err_msg,
        "disease_features": d_feat,
        "control_features": n_feat,
        "consensus_source": pick_source,
        "complexity_source": complex_source,
        "adaptive_read_cap_rerun": adaptive_rerun_used,
        "adaptive_read_cap_target": int(escalated_max_reads),
        "disease_recruited_evidence_read_names": sorted(d_recruited_evidence_read_names),
        "control_recruited_evidence_read_names": sorted(n_recruited_evidence_read_names),
    }
    manifest_path.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True), encoding="utf-8")

    result = {
        "chrom": str(as_row.get("chrom", "")),
        "window_start": int(as_row.get("window_start", 0)),
        "window_end": int(as_row.get("window_end", 0)),
        "asm_status": status,
        "asm_interval_start": int(current_start),
        "asm_interval_end": int(current_end),
        "asm_retry_used": bool(retry_used),
        "asm_disease_reads_extracted": int(d_reads),
        "asm_control_reads_extracted": int(n_reads),
        "asm_disease_contig_count": int(d_cc),
        "asm_control_contig_count": int(n_cc),
        "asm_disease_max_contig_len": int(d_max),
        "asm_control_max_contig_len": int(n_max),
        "asm_error_message": err_msg,
        "asm_disease_primary_contig_id": str(d_feat.get("primary_contig_id", "")),
        "asm_control_primary_contig_id": str(n_feat.get("primary_contig_id", "")),
        "asm_consensus_primary_contig_id": str(pick.get("primary_contig_id", "") if pick else ""),
        "asm_disease_breakpoint_pos": int(d_feat.get("breakpoint_pos", fallback_bp)),
        "asm_control_breakpoint_pos": int(n_feat.get("breakpoint_pos", fallback_bp)),
        "asm_consensus_breakpoint_pos": int(pick.get("breakpoint_pos", fallback_bp) if pick else fallback_bp),
        "asm_breakpoint_source": pick_source,
        "asm_tsd_seq": str(pick.get("tsd_seq", "") if pick else ""),
        "asm_tsd_len": int(pick.get("tsd_len", 0) if pick else 0),
        "asm_polyA_max_run": int(pick.get("polyA_max_run", 0) if pick else 0),
        "asm_insertion_length": int(pick.get("insertion_length", 0) if pick else 0),
        "asm_insertion_length_observed": int(pick.get("insertion_length_observed", 0) if pick else 0),
        "asm_insertion_length_imputed": int(pick.get("insertion_length_imputed", 0) if pick else 0),
        "asm_insertion_length_confidence_tier": str(pick.get("insertion_length_confidence_tier", "undetermined") if pick else "undetermined"),
        "asm_insertion_mei_start": int(pick.get("insertion_mei_start", -1) if pick else -1),
        "asm_insertion_mei_end": int(pick.get("insertion_mei_end", -1) if pick else -1),
        "asm_mei_family": str(pick.get("mei_family", "") if pick else ""),
        "asm_mei_subfamily": str(pick.get("mei_subfamily", "") if pick else ""),
        "asm_mei_target_length": int(pick.get("mei_target_length", 0) if pick else 0),
        "asm_insertion_orientation": str(pick.get("insertion_orientation", "") if pick else ""),
        "asm_complex_class": str(complex_pick.get("complex_class", "") if complex_pick else ""),
        "asm_non_mei_partner_chrom": str(complex_pick.get("non_mei_partner_chrom", "") if complex_pick else ""),
        "asm_non_mei_partner_pos": int(complex_pick.get("non_mei_partner_pos", -1) if complex_pick else -1),
        "asm_non_mei_partner_type": str(complex_pick.get("non_mei_partner_type", "") if complex_pick else ""),
        "asm_breakpoint_side_status": str(complex_pick.get("breakpoint_side_status", "none") if complex_pick else "none"),
        "asm_complexity_source": complex_source,
        "asm_top_contigs": str(complex_pick.get("top_contigs", "") if complex_pick else ""),
        "asm_mei_alignment_preset": str(pick.get("mei_alignment_preset", "") if pick else ""),
        "asm_left_support_contig_id": str(pick.get("left_support_contig_id", "") if pick else ""),
        "asm_right_support_contig_id": str(pick.get("right_support_contig_id", "") if pick else ""),
        "asm_left_support_mei_start": int(pick.get("left_support_mei_start", -1) if pick else -1),
        "asm_left_support_mei_end": int(pick.get("left_support_mei_end", -1) if pick else -1),
        "asm_right_support_mei_start": int(pick.get("right_support_mei_start", -1) if pick else -1),
        "asm_right_support_mei_end": int(pick.get("right_support_mei_end", -1) if pick else -1),
        "asm_left_support_mei_aln_len": int(pick.get("left_support_mei_aln_len", 0) if pick else 0),
        "asm_right_support_mei_aln_len": int(pick.get("right_support_mei_aln_len", 0) if pick else 0),
        "asm_microhomology_sequence": str(pick.get("microhomology_sequence", "") if pick else ""),
        "asm_junction_overlap_sequence": str(pick.get("junction_overlap_sequence", "") if pick else ""),
        "asm_coord_model": str(pick.get("coord_model", "") if pick else ""),
        "asm_coord_logic_version": int(pick.get("coord_logic_version", 0) if pick else 0),
        "asm_adaptive_read_cap_rerun": bool(adaptive_rerun_used),
        "asm_adaptive_read_cap_target": int(escalated_max_reads),
        "asm_disease_recruited_evidence_read_names": ",".join(sorted(d_recruited_evidence_read_names)),
        "asm_control_recruited_evidence_read_names": ",".join(sorted(n_recruited_evidence_read_names)),
    }
    log_line = (
        f"[local-assembly] locus {idx}/{total_loci} complete status={status} "
        f"reads(d={d_reads},n={n_reads}) contigs(d={d_cc},n={n_cc}) source={pick_source or 'none'}"
    )
    return result, log_line


def _detect_total_memory_gb() -> int:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return 0
    try:
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                parts = line.split()
                if len(parts) >= 2:
                    # /proc/meminfo reports kB.
                    return max(0, int(int(parts[1]) / (1024 * 1024)))
    except Exception:
        return 0
    return 0


def _resolve_locus_workers(
    *,
    requested_workers: int,
    total_loci: int,
    spades_threads: int,
    spades_memory_gb: int,
) -> int:
    if total_loci <= 0:
        return 1
    if requested_workers > 0:
        return max(1, min(int(requested_workers), int(total_loci)))

    cpu_count = max(1, int(os.cpu_count() or 1))
    cpu_headroom = 2
    cpu_limited_workers = max(1, cpu_count - cpu_headroom)
    total_mem_gb = _detect_total_memory_gb()
    reserve_mem_gb = 12
    mem_per_worker_gb = max(1, int(spades_memory_gb))
    if total_mem_gb > 0:
        mem_limited_workers = max(1, int((total_mem_gb - reserve_mem_gb) / mem_per_worker_gb))
    else:
        mem_limited_workers = cpu_count

    # Protect against oversubscription when callers increase SPAdes threads.
    spades_thread_limited_workers = max(1, int(cpu_count / max(1, int(spades_threads))))
    auto_workers = min(cpu_limited_workers, mem_limited_workers, spades_thread_limited_workers, int(total_loci))
    return max(1, int(auto_workers))


def annotate_silver_with_local_assembly(
    candidates: pd.DataFrame,
    *,
    disease_bam_path: Path,
    control_bam_path: Path,
    assembly_cache_dir: Path,
    mei_fasta: Path,
    reference_fasta: Path | None = None,
    interval_pad_bp: int = 250,
    retry_pad_bp: int = 600,
    max_reads_per_sample: int = 600,
    spades_threads: int = 1,
    spades_memory_gb: int = 8,
    minimap2_threads: int = 1,
    locus_workers: int = 0,
    reuse_existing: bool = True,
    reuse_cache_only: bool = False,
    disease_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
    control_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "window_start",
                "window_end",
                "asm_status",
                "asm_interval_start",
                "asm_interval_end",
                "asm_retry_used",
                "asm_disease_reads_extracted",
                "asm_control_reads_extracted",
                "asm_disease_contig_count",
                "asm_control_contig_count",
                "asm_disease_max_contig_len",
                "asm_control_max_contig_len",
                "asm_error_message",
                "asm_disease_primary_contig_id",
                "asm_control_primary_contig_id",
                "asm_disease_breakpoint_pos",
                "asm_control_breakpoint_pos",
                "asm_consensus_breakpoint_pos",
                "asm_breakpoint_source",
                "asm_tsd_seq",
                "asm_tsd_len",
                "asm_polyA_max_run",
                "asm_insertion_length",
                "asm_insertion_length_observed",
                "asm_insertion_length_imputed",
                "asm_insertion_length_confidence_tier",
                "asm_insertion_mei_start",
                "asm_insertion_mei_end",
                "asm_mei_family",
                "asm_mei_subfamily",
                "asm_mei_target_length",
                "asm_insertion_orientation",
                "asm_consensus_primary_contig_id",
                "asm_complex_class",
                "asm_non_mei_partner_chrom",
                "asm_non_mei_partner_pos",
                "asm_non_mei_partner_type",
                "asm_breakpoint_side_status",
                "asm_complexity_source",
                "asm_top_contigs",
                "asm_mei_alignment_preset",
                "asm_left_support_contig_id",
                "asm_right_support_contig_id",
                "asm_left_support_mei_start",
                "asm_left_support_mei_end",
                "asm_right_support_mei_start",
                "asm_right_support_mei_end",
                "asm_left_support_mei_aln_len",
                "asm_right_support_mei_aln_len",
                "asm_microhomology_sequence",
                "asm_junction_overlap_sequence",
                "asm_coord_model",
                "asm_coord_logic_version",
                "asm_adaptive_read_cap_rerun",
                "asm_adaptive_read_cap_target",
                "asm_disease_recruited_evidence_read_names",
                "asm_control_recruited_evidence_read_names",
            ]
        )

    assembly_cache_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    silver_stage = (
        candidates["silver_stage_pass"]
        if "silver_stage_pass" in candidates.columns
        else pd.Series(False, index=candidates.index)
    )
    silver = candidates.loc[silver_stage.fillna(False).astype(bool)].copy()
    total_loci = int(len(silver))
    if total_loci == 0:
        return pd.DataFrame(rows)
    spades_exe: str | None = None
    if not reuse_cache_only:
        spades_exe = _resolve_spades_executable()
    workers = _resolve_locus_workers(
        requested_workers=int(locus_workers),
        total_loci=total_loci,
        spades_threads=int(spades_threads),
        spades_memory_gb=int(spades_memory_gb),
    )
    print(
        f"[local-assembly] starting loci={total_loci} "
        f"spades={spades_exe or 'not-needed-unless-cache-miss'} "
        f"requested_locus_workers={int(locus_workers)} "
        f"minimap2_threads={max(1, int(minimap2_threads))} "
        f"locus_workers={workers}",
        flush=True,
    )
    row_records = silver.to_dict(orient="records")
    if workers == 1:
        for idx, row_data in enumerate(row_records, start=1):
            result, log_line = _process_single_locus(
                idx=idx,
                total_loci=total_loci,
                row_data=row_data,
                disease_bam_path=disease_bam_path,
                control_bam_path=control_bam_path,
                assembly_cache_dir=assembly_cache_dir,
                mei_fasta=mei_fasta,
                reference_fasta=reference_fasta,
                interval_pad_bp=interval_pad_bp,
                retry_pad_bp=retry_pad_bp,
                max_reads_per_sample=max_reads_per_sample,
                spades_threads=spades_threads,
                spades_memory_gb=spades_memory_gb,
                minimap2_threads=minimap2_threads,
                reuse_existing=reuse_existing,
                reuse_cache_only=reuse_cache_only,
                spades_exe=spades_exe,
                disease_preferred_read_names_by_locus=disease_preferred_read_names_by_locus,
                control_preferred_read_names_by_locus=control_preferred_read_names_by_locus,
            )
            rows.append(result)
            print(log_line, flush=True)
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    _process_single_locus,
                    idx=idx,
                    total_loci=total_loci,
                    row_data=row_data,
                    disease_bam_path=disease_bam_path,
                    control_bam_path=control_bam_path,
                    assembly_cache_dir=assembly_cache_dir,
                    mei_fasta=mei_fasta,
                    reference_fasta=reference_fasta,
                    interval_pad_bp=interval_pad_bp,
                    retry_pad_bp=retry_pad_bp,
                    max_reads_per_sample=max_reads_per_sample,
                    spades_threads=spades_threads,
                    spades_memory_gb=spades_memory_gb,
                    minimap2_threads=minimap2_threads,
                    reuse_existing=reuse_existing,
                    reuse_cache_only=reuse_cache_only,
                    spades_exe=spades_exe,
                    disease_preferred_read_names_by_locus=disease_preferred_read_names_by_locus,
                    control_preferred_read_names_by_locus=control_preferred_read_names_by_locus,
                )
                for idx, row_data in enumerate(row_records, start=1)
            ]
            for fut in as_completed(futures):
                result, log_line = fut.result()
                rows.append(result)
                print(log_line, flush=True)
    return pd.DataFrame(rows)
