#!/usr/bin/env python3
"""Benchmark MEI-consensus remappers on sampled consensus tips.

Samples random substrings (and reverse complements) from Alu / LINE-1 / SVA
Dfam consensus entries, then measures hit-rate and wall time for:

  - minimap2 -x sr          (current pipeline default)
  - minimap2 -x sr -k10 -w5
  - bwa mem                 (default)
  - bwa mem -k10 -T10
  - bwa aln + samse
  - bowtie2 --very-sensitive-local

Usage (from repo root, inside rtm-miner env):

  PYTHONPATH=src python scripts/benchmark_mei_aligners.py \\
    --mei-fasta $RTM_PUBLIC_DATA_DIR/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta
"""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


FAMILIES = ("ALU", "LINE1", "SVA")
DEFAULT_LENGTHS = (20, 22, 25, 30, 40, 60, 100)
DEFAULT_N_PER_CELL = 40  # per family × length × strand


def _family_of(name: str) -> str | None:
    u = name.upper()
    if "SINE/ALU" in u or name.startswith("Alu") or "#SINE/Alu" in name:
        return "ALU"
    if "LINE/L1" in u or name.startswith("L1") or "#LINE/L1" in name:
        return "LINE1"
    if "SVA" in u or "#Retroposon/SVA" in name:
        return "SVA"
    return None


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _load_fasta(path: Path) -> dict[str, list[tuple[str, str]]]:
    by_fam: dict[str, list[tuple[str, str]]] = {f: [] for f in FAMILIES}
    name = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    fam = _family_of(name)
                    seq = "".join(chunks).upper().replace("U", "T")
                    if fam and len(seq) >= 20:
                        by_fam[fam].append((name.split()[0], seq))
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        fam = _family_of(name)
        seq = "".join(chunks).upper().replace("U", "T")
        if fam and len(seq) >= 20:
            by_fam[fam].append((name.split()[0], seq))
    return by_fam


def _sample_queries(
    by_fam: dict[str, list[tuple[str, str]]],
    *,
    lengths: tuple[int, ...],
    n_per_cell: int,
    seed: int,
) -> list[tuple[str, str, str, int, str]]:
    """Return list of (qid, seq, family, length, strand)."""
    rng = random.Random(seed)
    out: list[tuple[str, str, str, int, str]] = []
    i = 0
    for fam in FAMILIES:
        entries = by_fam.get(fam) or []
        if not entries:
            continue
        for length in lengths:
            usable = [e for e in entries if len(e[1]) >= length]
            if not usable:
                continue
            for _ in range(n_per_cell):
                src_name, src = rng.choice(usable)
                start = rng.randint(0, len(src) - length)
                seq = src[start : start + length]
                strand = "+"
                if rng.random() < 0.5:
                    seq = _revcomp(seq)
                    strand = "-"
                qid = f"{fam}_{length}_{strand}_{i}"
                out.append((qid, seq, fam, length, strand))
                i += 1
    return out


def _write_fa(path: Path, queries: list[tuple[str, str, str, int, str]]) -> None:
    with path.open("w") as fh:
        for qid, seq, *_ in queries:
            fh.write(f">{qid}\n{seq}\n")


def _paf_mapped(stdout: str) -> set[str]:
    hits: set[str] = set()
    for line in stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        cols = line.split("\t")
        if len(cols) >= 6 and cols[0]:
            # PAF: qname ... tname
            hits.add(cols[0])
    return hits


def _sam_mapped(stdout: str) -> set[str]:
    hits: set[str] = set()
    for line in stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        if cols[2] != "*":
            hits.add(cols[0])
    return hits


def _time_cmd(cmd: list[str], *, parse_sam: bool = False) -> tuple[float, set[str], str]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    out = proc.stdout or ""
    hits = _sam_mapped(out) if parse_sam else _paf_mapped(out)
    err = (proc.stderr or "").strip().splitlines()
    note = err[-1] if err else ""
    return elapsed, hits, note


def _ensure_bowtie_index(mei_fasta: Path, work: Path) -> str:
    idx = work / "mei_bt2"
    # Build once into work dir (small MEI DB ~290 KB).
    subprocess.run(
        ["bowtie2-build", "--quiet", str(mei_fasta), str(idx)],
        check=True,
        capture_output=True,
    )
    return str(idx)


def run_aligners(
    mei_fasta: Path,
    queries: list[tuple[str, str, str, int, str]],
    work: Path,
) -> dict[str, tuple[float, set[str]]]:
    qfa = work / "queries.fa"
    _write_fa(qfa, queries)
    bt_idx = _ensure_bowtie_index(mei_fasta, work)
    results: dict[str, tuple[float, set[str]]] = {}

    specs: list[tuple[str, list[str], bool]] = [
        ("minimap2 -x sr", ["minimap2", "-x", "sr", "--secondary=yes", "-c", str(mei_fasta), str(qfa)], False),
        (
            "minimap2 -x sr -k10 -w5",
            ["minimap2", "-x", "sr", "-k", "10", "-w", "5", "--secondary=yes", "-c", str(mei_fasta), str(qfa)],
            False,
        ),
        ("bwa mem", ["bwa", "mem", "-t", "1", "-a", str(mei_fasta), str(qfa)], True),
        (
            "bwa mem -k10 -T10",
            ["bwa", "mem", "-t", "1", "-k", "10", "-T", "10", "-a", str(mei_fasta), str(qfa)],
            True,
        ),
        (
            "bowtie2 --very-sensitive-local",
            ["bowtie2", "-p", "1", "-f", "-x", bt_idx, "-U", str(qfa), "--very-sensitive-local"],
            True,
        ),
    ]
    for label, cmd, sam in specs:
        elapsed, hits, _note = _time_cmd(cmd, parse_sam=sam)
        results[label] = (elapsed, hits)

    # bwa aln is two-step; time both.
    sai = work / "q.sai"
    t0 = time.perf_counter()
    with sai.open("wb") as fh:
        subprocess.run(
            ["bwa", "aln", "-t", "1", "-n", "2", "-l", "10", str(mei_fasta), str(qfa)],
            stdout=fh,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    proc = subprocess.run(
        ["bwa", "samse", str(mei_fasta), str(sai), str(qfa)],
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    results["bwa aln -l10"] = (elapsed, _sam_mapped(proc.stdout or ""))
    return results


def _summarize(
    queries: list[tuple[str, str, str, int, str]],
    results: dict[str, tuple[float, set[str]]],
) -> None:
    meta = {qid: (fam, length) for qid, _seq, fam, length, _strand in queries}
    methods = list(results.keys())

    print(f"\n=== Overall ({len(queries)} queries) ===")
    print(f"{'method':<32} {'hit%':>7} {'hits':>6} {'time_s':>8} {'qps':>8}")
    for m in methods:
        elapsed, hits = results[m]
        n = len(queries)
        rate = 100.0 * len(hits) / max(n, 1)
        qps = n / elapsed if elapsed > 0 else float("inf")
        print(f"{m:<32} {rate:6.1f}% {len(hits):6d} {elapsed:8.3f} {qps:8.0f}")

    print("\n=== Hit-rate by family × length (%) ===")
    # header
    lengths = sorted({length for _q, _s, _f, length, _st in queries})
    for fam in FAMILIES:
        print(f"\n--- {fam} ---")
        header = f"{'method':<32}" + "".join(f"{L:>6}" for L in lengths)
        print(header)
        for m in methods:
            hits = results[m][1]
            cells = []
            for L in lengths:
                ids = [qid for qid, (f, length) in meta.items() if f == fam and length == L]
                if not ids:
                    cells.append(f"{'—':>6}")
                    continue
                n_hit = sum(1 for qid in ids if qid in hits)
                cells.append(f"{100.0 * n_hit / len(ids):5.0f}%")
            print(f"{m:<32}" + "".join(f"{c:>6}" for c in cells))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--mei-fasta",
        type=Path,
        default=Path(
            "/home/ec2-user/retrotransposon-workdir/data/public/retrotransposon_db/dfam/"
            "dfam_human_mei_l1_alu_sva.fasta"
        ),
    )
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-per-cell", type=int, default=DEFAULT_N_PER_CELL)
    ap.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=list(DEFAULT_LENGTHS),
    )
    ap.add_argument(
        "--runtime-n",
        type=int,
        default=5000,
        help="Extra same-length (40 bp) batch size for wall-time scaling",
    )
    args = ap.parse_args()

    by_fam = _load_fasta(args.mei_fasta)
    print("Consensus entries loaded:")
    for fam in FAMILIES:
        n = len(by_fam[fam])
        lens = sorted(len(s) for _, s in by_fam[fam])
        if lens:
            print(f"  {fam}: {n} seqs  len[{lens[0]}..{lens[-1]}]")
        else:
            print(f"  {fam}: 0 seqs")

    lengths = tuple(int(x) for x in args.lengths)
    queries = _sample_queries(
        by_fam, lengths=lengths, n_per_cell=int(args.n_per_cell), seed=int(args.seed)
    )
    print(
        f"\nSensitivity set: {len(queries)} queries "
        f"({args.n_per_cell}/family/length × {len(lengths)} lengths × ~50% RC)"
    )

    work = Path(tempfile.mkdtemp(prefix="rtm_mei_aln_bench_"))
    try:
        results = run_aligners(args.mei_fasta, queries, work)
        _summarize(queries, results)

        # Runtime scaling batch: fixed 40 bp mixed families (closer to SR clip length).
        if int(args.runtime_n) > 0:
            print(f"\n=== Runtime batch ({args.runtime_n} × 40 bp mixed) ===")
            rt_queries = _sample_queries(
                by_fam, lengths=(40,), n_per_cell=max(1, int(args.runtime_n) // 3), seed=int(args.seed) + 1
            )[: int(args.runtime_n)]
            # pad if needed
            while len(rt_queries) < int(args.runtime_n) and queries:
                rt_queries.append(queries[len(rt_queries) % len(queries)])
            rt_queries = rt_queries[: int(args.runtime_n)]
            rt_work = work / "runtime"
            rt_work.mkdir(exist_ok=True)
            rt = run_aligners(args.mei_fasta, rt_queries, rt_work)
            print(f"{'method':<32} {'time_s':>8} {'qps':>8}")
            for m, (elapsed, _hits) in rt.items():
                qps = len(rt_queries) / elapsed if elapsed > 0 else float("inf")
                print(f"{m:<32} {elapsed:8.3f} {qps:8.0f}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
