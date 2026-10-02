#!/usr/bin/env python3
"""Reference-genome opportunity, host annotation, and consensus projection.

Self-contained helpers for the ten-genome multi-sample analysis.  Nothing here
imports the five-genome ``dedup_samples`` module, so the ten-genome numbers do
not depend on that file's in-flight state.

Opportunity tiers implemented here follow the task specification:

* Tier 1 -- reference-genome opportunity: hg38 RepeatMasker host spans minus
  N-gaps in the reference assembly.
* Tier 2 -- labelled sensitivity: the same, additionally excluding reference
  hosts that overlap a known polymorphic MEI insertion.
* Tier 3 -- carriage-evidenced counts are reported descriptively only and are
  never used as a denominator, because conditioning on the calls under test is
  circular.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np

# Required wording for every Tier 1 / Tier 2 output.
REFERENCE_CAVEAT = (
    "conditional on reference-genome opportunity; reference host spans are "
    "assumed present in the cohort"
)
TIER2_LABEL = "excluding known-MEI-overlapping hosts"
TIER3_RATIONALE = (
    "carriage-evidenced counts are reported descriptively only and are never "
    "used as an opportunity denominator, because conditioning on the calls "
    "being tested is circular"
)

FAMILIES = ("ALU", "LINE1", "SVA")
PRIMARY_CHROMS: tuple[str, ...] = tuple(
    [f"chr{i}" for i in range(1, 23)] + ["chrX"]
)


# --------------------------------------------------------------------------
# Reference FASTA access
# --------------------------------------------------------------------------
def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    """Parse a samtools .fai into {name: (length, offset, line_bases, line_width)}."""
    fai = Path(str(path) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(f"FASTA index missing: {fai}")
    meta: dict[str, tuple[int, int, int, int]] = {}
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            meta[fields[0]] = (
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
                int(fields[4]),
            )
    return meta


def fasta_chunk(
    handle, meta: tuple[int, int, int, int], start: int, end: int
) -> bytes:
    """Return reference bases [start, end) using .fai byte arithmetic."""
    length, offset, line_bases, line_width = meta
    if not 0 <= start < end <= length:
        raise ValueError(f"bad FASTA slice {start}:{end} of {length}")
    line_index, col = divmod(start, line_bases)
    byte_offset = offset + line_index * line_width + col
    n_bases = end - start
    newline_bytes = ((col + n_bases - 1) // line_bases) * (line_width - line_bases)
    handle.seek(byte_offset)
    raw = handle.read(n_bases + newline_bytes)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")
    if len(seq) < n_bases:
        raise ValueError(f"short FASTA read at bases {start}:{end}")
    return seq[:n_bases].upper()


def scan_gaps(
    fasta: Path,
    meta: dict[str, tuple[int, int, int, int]],
    chroms: Sequence[str],
    *,
    min_run: int = 1,
    chunk_bp: int = 2_000_000,
) -> dict[str, np.ndarray]:
    """Find runs of N (or any non-ACGT) in the reference, per chromosome.

    Returns {chrom: int64 array of shape (n, 2)} with half-open [start, end)
    zero-based coordinates.  ``min_run`` keeps only runs of at least that many
    bases; the ten-genome run uses min_run=1 so that *any* uncallable base is
    removed from the opportunity denominator.
    """
    gaps: dict[str, np.ndarray] = {}
    with Path(fasta).open("rb") as handle:
        for chrom in chroms:
            if chrom not in meta:
                gaps[chrom] = np.zeros((0, 2), dtype=np.int64)
                continue
            length = meta[chrom][0]
            runs: list[tuple[int, int]] = []
            start0 = 0
            while start0 < length:
                end0 = min(start0 + chunk_bp, length)
                arr = np.frombuffer(
                    fasta_chunk(handle, meta[chrom], start0, end0), dtype=np.uint8
                )
                bad = ~np.isin(
                    arr,
                    np.asarray(
                        [ord("A"), ord("C"), ord("G"), ord("T"), ord("N")],
                        dtype=np.uint8,
                    ),
                )
                # treat N as a gap; any other non-ACGT byte is also uncallable
                bad |= arr == ord("N")
                idx = np.flatnonzero(np.diff(np.concatenate(([0], bad.view(np.int8), [0]))))
                for lo, hi in zip(idx[0::2], idx[1::2]):
                    gs, ge = start0 + int(lo), start0 + int(hi)
                    if ge - gs >= min_run:
                        runs.append((gs, ge))
                start0 = end0
            gaps[chrom] = (
                np.asarray(runs, dtype=np.int64).reshape(-1, 2)
                if runs
                else np.zeros((0, 2), dtype=np.int64)
            )
    return gaps


# --------------------------------------------------------------------------
# Interval algebra
# --------------------------------------------------------------------------
def merge(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for start, end in sorted(spans):
        if end <= start:
            continue
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(s, e) for s, e in out]


def subtract(
    spans: Sequence[tuple[int, int]], masks: Sequence[tuple[int, int]]
) -> list[tuple[int, int]]:
    """spans minus masks, both half-open and sorted-disjoint."""
    out: list[tuple[int, int]] = []
    for start, end in spans:
        cursor = start
        for m_start, m_end in masks:
            if m_end <= cursor:
                continue
            if m_start >= end:
                break
            if m_start > cursor:
                out.append((cursor, min(m_start, end)))
            cursor = max(cursor, m_end)
            if cursor >= end:
                break
        if cursor < end:
            out.append((cursor, end))
    return merge(out)


def overlap(spans: Sequence[tuple[int, int]], lo: int, hi: int) -> int:
    """Total bases shared between a sorted span list and [lo, hi)."""
    total = 0
    for start, end in spans:
        if start >= hi:
            break
        total += max(0, min(end, hi) - max(start, lo))
    return total


def bin_opportunity(
    spans: Sequence[tuple[int, int]],
    edges: Sequence[int],
    *,
    relative: bool = False,
) -> np.ndarray:
    """Callable bases per bin defined by consecutive ``edges``.

    With ``relative=True`` each span is rescaled to fill the bin range, so
    every span contributes exactly ``len(span)/n_bins`` to each bin regardless
    of its absolute length.  That is the relative-position sensitivity: it
    deliberately includes truncated hosts, and that truncation bias is stated
    in the report.
    """
    edges_arr = np.asarray(edges, dtype=np.int64)
    n_bins = max(0, edges_arr.size - 1)
    out = np.zeros(n_bins, dtype=np.int64)
    if n_bins == 0:
        return out
    if relative:
        per_bin: list[float] = [0.0] * n_bins
        for start, end in spans:
            span_len = end - start
            if span_len <= 0:
                continue
            share = span_len / float(n_bins)
            for b in range(n_bins):
                per_bin[b] += share
        return np.asarray([int(round(v)) for v in per_bin], dtype=np.int64)
    for start, end in spans:
        if end <= start:
            continue
        lo = max(0, int(np.searchsorted(edges_arr, start, side="right") - 1))
        hi = min(n_bins, int(np.searchsorted(edges_arr, end, side="left")))
        for b in range(lo, hi):
            out[b] += min(end, int(edges_arr[b + 1])) - max(start, int(edges_arr[b]))
    return out


# --------------------------------------------------------------------------
# RepeatMasker hosts
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Host:
    chrom: str
    start0: int
    end0: int
    strand: str
    name: str
    family: str

    @property
    def length(self) -> int:
        return self.end0 - self.start0

    @property
    def key(self) -> tuple[str, int, int, str, str]:
        return (self.chrom, self.start0, self.end0, self.name, self.strand)

    @property
    def host_id(self) -> str:
        return f"{self.chrom}:{self.start0}-{self.end0}:{self.name}:{self.strand}"


def normalize_family(token: str) -> str:
    value = (token or "").upper()
    if "ALU" in value:
        return "ALU"
    if "SVA" in value:
        return "SVA"
    if "LINE1" in value or "L1" in value:
        return "LINE1"
    return ""


def read_rmsk(
    path: Path, *, chroms: Sequence[str] = PRIMARY_CHROMS
) -> dict[tuple[str, str], list[Host]]:
    """Parse hg38 RepeatMasker into {(chrom, family): [Host, ...]}.

    Column layout is the UCSC rmsk.txt layout used by this project:
    5=chrom 6=start 7=end 9=strand 10=name 11=class 12=family.
    """
    keep = set(chroms)
    out: dict[tuple[str, str], list[Host]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13 or fields[5] not in keep:
                continue
            family = normalize_family(f"{fields[10]} {fields[11]} {fields[12]}")
            if family not in FAMILIES:
                continue
            try:
                start0, end0 = int(fields[6]), int(fields[7])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            out.setdefault((fields[5], family), []).append(
                Host(fields[5], start0, end0, fields[9], fields[10], family)
            )
    for hosts in out.values():
        hosts.sort(key=lambda h: (h.start0, h.end0, h.name, h.strand))
    return out


def host_segments(
    host: Host,
    gaps: np.ndarray,
    masks: Sequence[tuple[int, int]] | None = None,
) -> list[tuple[int, int]]:
    """Callable segments of one host after removing gaps (and optional masks)."""
    masked = [(int(s), int(e)) for s, e in gaps.tolist()] if gaps.size else []
    segs = subtract([(host.start0, host.end0)], merge(masked))
    if masks:
        segs = subtract(segs, merge(masks))
    return segs


def bin_host_opportunity(
    hosts: Sequence[Host],
    gaps_by_chrom: dict[str, np.ndarray],
    edges_by_chrom: dict[str, np.ndarray],
    *,
    excluded: frozenset[tuple[str, int, int, str, str]] = frozenset(),
    relative: bool = False,
) -> dict[str, np.ndarray]:
    """Opportunity per genomic bin, per chromosome, keyed by chromosome."""
    out: dict[str, np.ndarray] = {}
    by_chrom: dict[str, list[Host]] = {}
    for host in hosts:
        if host.key in excluded:
            continue
        by_chrom.setdefault(host.chrom, []).append(host)
    for chrom, chrom_hosts in by_chrom.items():
        if chrom not in edges_by_chrom:
            continue
        edges = np.asarray(edges_by_chrom[chrom], dtype=np.int64)
        segs: list[tuple[int, int]] = []
        for host in chrom_hosts:
            segs.extend(host_segments(host, gaps_by_chrom.get(chrom, _EMPTY)))
        out[chrom] = bin_opportunity(merge(segs), edges.tolist(), relative=relative)
    return out


_EMPTY = np.zeros((0, 2), dtype=np.int64)


# --------------------------------------------------------------------------
# L1 subfamily consensus (hg38reps.fa) and projection
# --------------------------------------------------------------------------
def load_consensus(path: Path, *, families: Sequence[str] = ("L1",)) -> dict[str, str]:
    """Load subfamily consensus sequences keyed by name.

    ``hg38reps.fa`` is a plain FASTA whose record names are RepeatMasker
    consensus identifiers such as ``L1PA2`` or ``AluJb``.
    """
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    if not families:
        return seqs
    return {k: v for k, v in seqs.items() if any(k.startswith(f) for f in families)}


def kmer_set(seq: str, k: int = 16) -> frozenset[str]:
    return frozenset(seq[i : i + k] for i in range(len(seq) - k + 1))


def assign_subfamily(
    seq: str, consensus_kmers: dict[str, frozenset[str]]
) -> tuple[str | None, float]:
    """Pick the consensus sharing the most 16-mers with ``seq``.

    A fast seed-containment score; the exact breakpoint projection is done
    afterwards by alignment for the hosts that actually carry a child.
    """
    best_name: str | None = None
    best_score = 0.0
    kmers = kmer_set(seq)
    if not kmers:
        return None, 0.0
    for name, ref in consensus_kmers.items():
        shared = len(kmers & ref)
        if shared:
            score = shared / float(len(kmers))
            if score > best_score:
                best_score, best_name = score, name
    return best_name, best_score


def project_to_consensus(
    host_seq: str,
    consensus: str,
    *,
    host_offset: int,
) -> tuple[int, str] | None:
    """Map a host-relative offset onto a consensus coordinate.

    ``host_offset`` is 0-based within ``host_seq`` (already strand-oriented so
    that index 0 is the 5' end of the host on its own strand).  Returns
    ``(consensus_offset, orientation)`` or None when the alignment is too poor
    to trust.  The consensus is reversed when the host aligns to it backwards,
    which is exactly the "oriented 5'-to-3' of the host strand" requirement.
    """
    from Bio import Align  # imported lazily: only Q4 needs it

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    def _map(query: str, target: str, orientation: str) -> tuple[int, str] | None:
        try:
            alignment = aligner.align(query, target)[0]
        except (ValueError, IndexError):
            return None
        query_to_target: dict[int, int] = {}
        for (q_start, q_end), (t_start, t_end) in zip(
            alignment.aligned[0], alignment.aligned[1]
        ):
            for offset in range(q_end - q_start):
                query_to_target[q_start + offset] = t_start + offset
        if host_offset not in query_to_target:
            return None
        return query_to_target[host_offset], orientation

    direct = _map(host_seq, consensus, "forward")
    if direct is not None:
        return direct
    rc = consensus[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    reverse = _map(host_seq, rc, "reverse")
    if reverse is not None:
        return reverse
    return None


def consensus_bin_edges(consensus_length: int, bin_bp: int = 500) -> list[int]:
    edges = list(range(0, consensus_length, bin_bp))
    if edges[-1] != consensus_length:
        edges.append(consensus_length)
    return edges


def iter_host_sequences(
    fasta: Path,
    meta: dict[str, tuple[int, int, int, int]],
    hosts: Iterable[Host],
) -> Iterator[tuple[Host, str]]:
    """Yield (host, strand-oriented sequence) for the given hosts."""
    by_chrom: dict[str, list[Host]] = {}
    for host in hosts:
        by_chrom.setdefault(host.chrom, []).append(host)
    with Path(fasta).open("rb") as handle:
        for chrom, chrom_hosts in by_chrom.items():
            if chrom not in meta:
                continue
            for host in sorted(chrom_hosts, key=lambda h: h.start0):
                seq = fasta_chunk(handle, meta[chrom], host.start0, host.end0)
                text = seq.decode("ascii", errors="replace")
                if host.strand == "-":
                    text = text[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
                yield host, text
