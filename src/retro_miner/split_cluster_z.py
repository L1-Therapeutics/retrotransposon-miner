"""Binomial z-score for split reads that share a junction cluster.

A read counts as clustered when it shares a clip side with other reads, its
junction is within 15 bp of the cluster median, and its junction-proximal
16-mer matches a cluster member with at most 2 mismatches. Only clusters of
3 or more reads count.

The null rate is the read-weighted cluster rate among silver calls on the
same chromosome that have at least 8 split reads in the discovery window.
The score is the binomial z of the clustered count:

    z = (k - n * p) / sqrt(n * p * (1 - p))

This uses the split-read count, so 0 clustered reads out of 60 is a stronger
departure than 0 out of 8. Gold drops calls below z = -2. Split reads in a
window are not independent draws, and silver windows are overdispersed
relative to a binomial, so the numerical z is not a calibrated tail probability.
"""

from __future__ import annotations

import math

import pandas as pd

POS_BP = 15
KMER = 16
MAX_MISMATCH = 2
MIN_CLUSTER = 3
MIN_BACKGROUND_SPLITS = 8

WINDOW_READS_COL = "split_cluster_window_reads"
CLUSTERED_READS_COL = "split_cluster_reads"
BINOMIAL_P_COL = "split_cluster_binomial_p"
BINOMIAL_Z_COL = "split_cluster_binomial_z"


def binomial_z(clustered: int, n_reads: int, rate: float) -> float:
    """Binomial z of ``clustered`` successes in ``n_reads`` at ``rate``."""
    if n_reads <= 0 or not math.isfinite(rate) or rate <= 0.0 or rate >= 1.0:
        return float("nan")
    mean = n_reads * rate
    sd = math.sqrt(n_reads * rate * (1.0 - rate))
    if sd <= 0.0:
        return float("nan")
    return (clustered - mean) / sd


def _junction_kmer(side: str, seq: str) -> str:
    if len(seq) < 8:
        return ""
    k = min(KMER, len(seq))
    return seq[-k:] if side == "L" else seq[:k]


def _hamming_ok(side: str, left: str, right: str) -> bool:
    n = min(len(left), len(right))
    if n < 8:
        return False
    pair = zip(left[-n:], right[-n:]) if side == "L" else zip(left[:n], right[:n])
    return sum(a != b for a, b in pair) <= MAX_MISMATCH


class _Cluster:
    __slots__ = ("side", "positions", "kmers", "n")

    def __init__(self, side: str, pos: int, kmer: str) -> None:
        self.side = side
        self.positions = [pos]
        self.kmers = [kmer]
        self.n = 1

    def accepts(self, pos: int, kmer: str) -> bool:
        med = sorted(self.positions)[len(self.positions) // 2]
        if abs(pos - med) > POS_BP:
            return False
        return any(_hamming_ok(self.side, kmer, prev) for prev in self.kmers[-6:])

    def add(self, pos: int, kmer: str) -> None:
        self.positions.append(pos)
        self.kmers.append(kmer)
        self.n += 1


def cluster_assigned_reads(pos: list[int], side: list[str], seq: list[str]) -> int:
    """Return how many reads fall in a junction cluster of size at least 3."""
    clusters: list[_Cluster] = []
    order = sorted(range(len(pos)), key=lambda i: (side[i], pos[i]))
    for i in order:
        sd = side[i]
        if sd not in {"L", "R"}:
            continue
        kmer = _junction_kmer(sd, seq[i])
        if not kmer:
            continue
        placed = False
        p = pos[i]
        for cluster in clusters:
            if cluster.side == sd and cluster.accepts(p, kmer):
                cluster.add(p, kmer)
                placed = True
                break
        if not placed:
            clusters.append(_Cluster(sd, p, kmer))
    return sum(cluster.n for cluster in clusters if cluster.n >= MIN_CLUSTER)


class _ChromSplits:
    def __init__(self, frame: pd.DataFrame) -> None:
        ordered = frame.sort_values("pos", kind="mergesort")
        self.pos = [int(v) for v in ordered["pos"].tolist()]
        self.side = ordered["clip_side"].tolist()
        self.seq = ordered["clip_seq"].tolist()

    def score(self, start: int, end: int) -> tuple[int, int]:
        import bisect

        i0 = bisect.bisect_left(self.pos, start)
        i1 = bisect.bisect_right(self.pos, end)
        n = i1 - i0
        if n <= 0:
            return 0, 0
        k = cluster_assigned_reads(self.pos[i0:i1], self.side[i0:i1], self.seq[i0:i1])
        return n, k


def _truthy(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).isin(["1", "True", "true", "TRUE"])


def _window_bounds(candidates: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    if "discovery_window_start" in candidates.columns and "discovery_window_end" in candidates.columns:
        start = pd.to_numeric(candidates["discovery_window_start"], errors="coerce")
        end = pd.to_numeric(candidates["discovery_window_end"], errors="coerce")
    else:
        start = pd.to_numeric(candidates.get("window_start"), errors="coerce")
        end = pd.to_numeric(candidates.get("window_end"), errors="coerce")
    return start.fillna(-1).astype(int), end.fillna(-1).astype(int)


def _prepare_splits(split_reads: pd.DataFrame) -> dict[str, _ChromSplits]:
    if split_reads is None or split_reads.empty:
        return {}
    required = {"chrom", "pos", "clip_side", "clip_seq"}
    missing = required.difference(split_reads.columns)
    if missing:
        raise ValueError(f"split reads missing columns: {sorted(missing)}")
    frame = split_reads.loc[:, ["chrom", "pos", "clip_side", "clip_seq"]].copy()
    frame["chrom"] = frame["chrom"].astype(str)
    frame["pos"] = pd.to_numeric(frame["pos"], errors="coerce")
    frame = frame.loc[frame["pos"].notna()].copy()
    frame["pos"] = frame["pos"].astype(int)
    frame["clip_side"] = frame["clip_side"].fillna("").astype(str).str.upper().str[:1]
    frame["clip_seq"] = frame["clip_seq"].fillna("").astype(str).str.upper()
    return {chrom: _ChromSplits(part) for chrom, part in frame.groupby("chrom", sort=False)}


def annotate_split_cluster_binomial_z(
    candidates: pd.DataFrame,
    split_reads: pd.DataFrame,
) -> pd.DataFrame:
    """Add discovery-window cluster counts and a silver-relative binomial z.

    ``split_reads`` is the disease split-evidence table. Germline runs clone
    control evidence from disease, so the control table is not included.
    """
    out = candidates.copy()
    n_rows = len(out)
    n_reads = [0] * n_rows
    n_clustered = [0] * n_rows
    rates = [float("nan")] * n_rows
    scores = [float("nan")] * n_rows
    if n_rows == 0:
        out[WINDOW_READS_COL] = pd.Series(dtype="int64")
        out[CLUSTERED_READS_COL] = pd.Series(dtype="int64")
        out[BINOMIAL_P_COL] = pd.Series(dtype="float64")
        out[BINOMIAL_Z_COL] = pd.Series(dtype="float64")
        return out

    indexes = _prepare_splits(split_reads)
    chrom = out["chrom"].astype(str) if "chrom" in out.columns else pd.Series([""] * n_rows, index=out.index)
    start, end = _window_bounds(out)
    silver = (
        _truthy(out["silver_stage_pass"])
        if "silver_stage_pass" in out.columns
        else pd.Series(False, index=out.index)
    )
    # Score each row once. Background rate is per chromosome.
    scored: list[tuple[str, int, int]] = []
    for i, (chrom_i, a, b) in enumerate(zip(chrom.tolist(), start.tolist(), end.tolist())):
        index = indexes.get(str(chrom_i))
        if index is None or a < 0 or b < a:
            scored.append((str(chrom_i), 0, 0))
            continue
        n, k = index.score(int(a), int(b))
        scored.append((str(chrom_i), n, k))

    background: dict[str, tuple[int, int]] = {}
    for i, (chrom_i, n, k) in enumerate(scored):
        if not bool(silver.iloc[i]) or n < MIN_BACKGROUND_SPLITS:
            continue
        prev_n, prev_k = background.get(chrom_i, (0, 0))
        background[chrom_i] = (prev_n + n, prev_k + k)

    for i, (chrom_i, n, k) in enumerate(scored):
        n_reads[i] = n
        n_clustered[i] = k
        total_n, total_k = background.get(chrom_i, (0, 0))
        if total_n <= 0:
            continue
        rate = total_k / total_n
        rates[i] = rate
        scores[i] = binomial_z(k, n, rate)

    out[WINDOW_READS_COL] = n_reads
    out[CLUSTERED_READS_COL] = n_clustered
    out[BINOMIAL_P_COL] = rates
    out[BINOMIAL_Z_COL] = scores
    return out
