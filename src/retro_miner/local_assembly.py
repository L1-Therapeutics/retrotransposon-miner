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


def _write_interval_fastq(
    bam_path: Path,
    chrom: str,
    start_1based: int,
    end_1based: int,
    out_fastq_gz: Path,
    *,
    max_reads: int,
    non_perfect_only: bool = False,
) -> int:
    start0 = max(0, int(start_1based) - 1)
    end0 = max(start0 + 1, int(end_1based))
    written = 0
    out_fastq_gz.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, gzip.open(out_fastq_gz, "wt", encoding="utf-8") as oh:
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_sequence is None:
                continue
            if non_perfect_only and not _is_non_perfect_primary_alignment(read):
                continue
            seq = str(read.query_sequence)
            qual = read.qual if read.qual is not None else ("I" * len(seq))
            rid = f"{read.query_name}:{'1' if read.is_read1 else '2' if read.is_read2 else '0'}"
            oh.write(f"@{rid}\n{seq}\n+\n{qual}\n")
            written += 1
            if written >= int(max_reads):
                break
    return written


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
    if int(d_feat.get("mei_aln_len", 0)) >= int(n_feat.get("mei_aln_len", 0)) and int(d_feat.get("mei_aln_len", 0)) > 0:
        pick = d_feat
        pick_source = "disease"
    elif int(n_feat.get("mei_aln_len", 0)) > 0:
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
    if not records:
        return {
            "primary_contig_id": "",
            "mei_subfamily": "",
            "mei_family": "",
            "insertion_mei_start": -1,
            "insertion_mei_end": -1,
            "insertion_orientation": "",
            "insertion_length": 0,
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
        }

    with tempfile.TemporaryDirectory(prefix="rtm_asm_parse_") as td:
        query_fa = Path(td) / "contigs.fa"
        with query_fa.open("w", encoding="utf-8") as oh:
            for name, seq in records:
                oh.write(f">{name}\n{seq}\n")
        mei_hits = _run_minimap2_paf(query_fa, mei_fasta, preset="asm5", threads=minimap2_threads)
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
            "insertion_mei_start": -1,
            "insertion_mei_end": -1,
            "insertion_orientation": "",
            "insertion_length": 0,
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
        }

    best = max(mei_hits, key=lambda h: (int(h["alnlen"]), int(h["mapq"]), int(h["nmatch"])))
    contig_seq = dict(records).get(str(best["qname"]), "")
    ins_len = max(0, int(best["qend"]) - int(best["qstart"]))
    bp, bp_left_chrom, _bp_right_chrom, tsd_len = _infer_breakpoint_from_alignments(
        best, ref_hits, default_pos=default_breakpoint_pos
    )
    tsd_seq = ""
    if reference_fasta is not None and bp_left_chrom and tsd_len > 0 and tsd_len <= 50:
        try:
            with pysam.FastaFile(str(reference_fasta)) as ref:
                start0 = max(0, int(bp) - (tsd_len // 2) - 1)
                end0 = start0 + int(tsd_len)
                tsd_seq = ref.fetch(bp_left_chrom, start0, end0).upper()
        except Exception:
            tsd_seq = ""
    cx = _summarize_complex_topology(records, mei_hits, ref_hits, top_k=3)
    return {
        "primary_contig_id": str(best["qname"]),
        "mei_subfamily": str(best["tname"]),
        "mei_family": _family_from_target(str(best["tname"])),
        "insertion_mei_start": int(best["tstart"]) + 1,
        "insertion_mei_end": int(best["tend"]),
        "insertion_orientation": str(best["strand"]) if str(best["strand"]) in {"+", "-"} else "",
        "insertion_length": int(ins_len),
        "polyA_max_run": int(_poly_at_max_run(contig_seq)),
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
) -> tuple[dict[str, object], str]:
    as_row = pd.Series(row_data)
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
    d_cc = 0
    n_cc = 0
    d_max = 0
    n_max = 0
    chosen_pad = int(interval_pad_bp)
    escalated_max_reads = 1000

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
                d_reads = _write_interval_fastq(
                    disease_bam_path,
                    chrom,
                    current_start,
                    current_end,
                    d_fq,
                    max_reads=int(phase_cap),
                    non_perfect_only=non_perfect_only,
                )
                n_reads = _write_interval_fastq(
                    control_bam_path,
                    chrom,
                    current_start,
                    current_end,
                    n_fq,
                    max_reads=int(phase_cap),
                    non_perfect_only=non_perfect_only,
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
                    assembled_from_attempt = True
                    break

                retry_used = attempt > 1
                status = "assembled_retry" if retry_used else "assembled"
                err_msg = ""
                assembled_from_attempt = True
                break
            if assembled_from_attempt:
                break

    d_out = locus_dir / f"disease.spades.pad{int(chosen_pad)}"
    n_out = locus_dir / f"control.spades.pad{int(chosen_pad)}"
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
    if _should_escalate_read_cap(
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
        d_reads_hi = _write_interval_fastq(
            disease_bam_path,
            chrom,
            current_start,
            current_end,
            d_fq_hi,
            max_reads=int(escalated_max_reads),
        )
        n_reads_hi = _write_interval_fastq(
            control_bam_path,
            chrom,
            current_start,
            current_end,
            n_fq_hi,
            max_reads=int(escalated_max_reads),
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
        "asm_insertion_mei_start": int(pick.get("insertion_mei_start", -1) if pick else -1),
        "asm_insertion_mei_end": int(pick.get("insertion_mei_end", -1) if pick else -1),
        "asm_mei_family": str(pick.get("mei_family", "") if pick else ""),
        "asm_mei_subfamily": str(pick.get("mei_subfamily", "") if pick else ""),
        "asm_insertion_orientation": str(pick.get("insertion_orientation", "") if pick else ""),
        "asm_complex_class": str(complex_pick.get("complex_class", "") if complex_pick else ""),
        "asm_non_mei_partner_chrom": str(complex_pick.get("non_mei_partner_chrom", "") if complex_pick else ""),
        "asm_non_mei_partner_pos": int(complex_pick.get("non_mei_partner_pos", -1) if complex_pick else -1),
        "asm_non_mei_partner_type": str(complex_pick.get("non_mei_partner_type", "") if complex_pick else ""),
        "asm_breakpoint_side_status": str(complex_pick.get("breakpoint_side_status", "none") if complex_pick else "none"),
        "asm_complexity_source": complex_source,
        "asm_top_contigs": str(complex_pick.get("top_contigs", "") if complex_pick else ""),
        "asm_adaptive_read_cap_rerun": bool(adaptive_rerun_used),
        "asm_adaptive_read_cap_target": int(escalated_max_reads),
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
    total_mem_gb = _detect_total_memory_gb()
    reserve_mem_gb = 12
    mem_per_worker_gb = max(1, int(spades_memory_gb))
    if total_mem_gb > 0:
        mem_limited_workers = max(1, int((total_mem_gb - reserve_mem_gb) / mem_per_worker_gb))
    else:
        mem_limited_workers = cpu_count

    # Protect against oversubscription when callers increase SPAdes threads.
    spades_thread_limited_workers = max(1, int(cpu_count / max(1, int(spades_threads))))
    auto_workers = min(cpu_count, mem_limited_workers, spades_thread_limited_workers, int(total_loci))
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
                "asm_insertion_mei_start",
                "asm_insertion_mei_end",
                "asm_mei_family",
                "asm_mei_subfamily",
                "asm_insertion_orientation",
                "asm_consensus_primary_contig_id",
                "asm_complex_class",
                "asm_non_mei_partner_chrom",
                "asm_non_mei_partner_pos",
                "asm_non_mei_partner_type",
                "asm_breakpoint_side_status",
                "asm_complexity_source",
                "asm_top_contigs",
                "asm_adaptive_read_cap_rerun",
                "asm_adaptive_read_cap_target",
            ]
        )

    assembly_cache_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    silver = candidates.loc[candidates.get("silver_stage_pass", False).fillna(False).astype(bool)].copy()
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
                )
                for idx, row_data in enumerate(row_records, start=1)
            ]
            for fut in as_completed(futures):
                result, log_line = fut.result()
                rows.append(result)
                print(log_line, flush=True)
    return pd.DataFrame(rows)
