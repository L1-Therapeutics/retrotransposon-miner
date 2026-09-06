"""Parity tests: vectorized helpers vs scalar reference implementations.

The production helpers in ``retro_miner`` favor vectorized/compiled primitives
over hand-written Python loops: NumPy ``diff``/``split`` for spatial clustering
(``candidate_loci._cluster_sorted_positions``), a compiled regex tokenizer for
CIGAR parsing (``mei_support``), and regex/two-pointer scans for polyA/T signal
(``_utils``).  Each helper is re-implemented here as a deliberately-scalar
reference and the two code paths are asserted to produce identical output over
hand-picked edge cases and seeded randomized inputs.
"""

from __future__ import annotations

import random
from types import SimpleNamespace

import pytest

from retro_miner._utils import _longest_poly_at_span, _poly_at_stats
from retro_miner.candidate_loci import _cluster_sorted_positions
from retro_miner.evidence_extract import _collect_soft_clips, _longest_soft_clip_from_read
from retro_miner.mei_support import _cigar_alignment_spans, _cigar_query_len


# ─────────────────────────────────────────────────────────────────────────────
# Scalar reference implementations
# ─────────────────────────────────────────────────────────────────────────────


def _ref_cluster_sorted_positions(positions: list[int], max_gap_bp: int) -> list[list[int]]:
    """Scalar reference for :func:`_cluster_sorted_positions` — one loop, no NumPy."""
    if not positions:
        return []
    gap = max(0, int(max_gap_bp))
    clusters: list[list[int]] = [[int(positions[0])]]
    for pos in positions[1:]:
        if int(pos) - int(clusters[-1][-1]) > gap:
            clusters.append([int(pos)])
        else:
            clusters[-1].append(int(pos))
    return clusters


def _scan_cigar(cigar: str) -> list[tuple[str, int]]:
    """Manual character-scan tokenizer — the scalar alternative to the production regex."""
    ops: list[tuple[str, int]] = []
    buf = 0
    for ch in cigar or "":
        if ch.isdigit():
            buf = buf * 10 + int(ch)
        else:
            ops.append((ch, buf))
            buf = 0
    return ops


def _ref_cigar_alignment_spans(cigar: str) -> tuple[int, int, int]:
    q_aln = 0
    r_span = 0
    alnlen = 0
    for op, n in _scan_cigar(cigar):
        if op in {"M", "=", "X"}:
            q_aln += n
            r_span += n
            alnlen += n
        elif op == "I":
            q_aln += n
            alnlen += n
        elif op in {"D", "N"}:
            r_span += n
            alnlen += n
    return q_aln, r_span, alnlen


def _ref_cigar_query_len(cigar: str) -> int:
    total = 0
    for op, n in _scan_cigar(cigar):
        if op in {"M", "I", "=", "X", "S", "H"}:
            total += n
    return total


def _ref_soft_clip_bounds(cigar: str) -> tuple[int, int]:
    """Scalar reference: ``(left_soft_clip_bp, right_soft_clip_bp)`` from a CIGAR.

    Mirrors ``_collect_soft_clips`` / ``_longest_soft_clip_from_read``, which only
    inspect the first and last CIGAR tuple: a soft clip is recorded only when the
    very first/last op is ``S`` — a hard clip (``H``) occupying that slot hides it.
    """
    ops = _scan_cigar(cigar)
    if not ops:
        return (0, 0)
    left = ops[0][1] if ops[0][0] == "S" else 0
    right = ops[-1][1] if ops[-1][0] == "S" else 0
    return left, right


# pysam CIGAR op codes (M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8)
_CIGAR_OP_CODES = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}


def _fake_read(
    cigar: str,
    query_seq: str = "",
    ref_start: int = 1000,
    ref_end: int = 2000,
) -> SimpleNamespace:
    """Minimal stand-in for a ``pysam.AlignedSegment`` with a parsed CIGAR.

    Only the attributes consumed by the evidence-extract soft-clip helpers are
    populated, so the tests exercise production code without needing a BAM.
    """
    ops = _scan_cigar(cigar)
    return SimpleNamespace(
        cigartuples=[(_CIGAR_OP_CODES[op], n) for op, n in ops] or None,
        query_sequence=query_seq,
        reference_start=ref_start,
        reference_end=ref_end,
    )


def _ref_poly_at_stats(seq: str) -> tuple[int, float, str]:
    """Scalar reference for :func:`_poly_at_stats` — counting loop, no regex."""
    s = (seq or "").upper()
    if not s:
        return (0, 0.0, "")
    n_a = s.count("A")
    n_t = s.count("T")
    if n_a <= 0 and n_t <= 0:
        return (0, 0.0, "")
    if n_a >= n_t:
        base, n_dom = "A", n_a
    else:
        base, n_dom = "T", n_t
    best = 0
    run = 0
    for ch in s:
        if ch == base:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return (best, float(n_dom) / float(len(s)), base)


def _ref_longest_poly_at_span(
    seq: str,
    *,
    min_frac: float = 0.90,
    min_len: int = 25,
) -> tuple[int, float, str, str]:
    """Scalar reference for :func:`_longest_poly_at_span` — no regex, explicit loops.

    Reproduces the exact production behavior: the greedy two-pointer that keeps a
    monotone ``left`` boundary and evaluates one window per right end, the same
    tie-breaking (length first, then fraction, base A before T), and the full-read
    override (best span within 2 bp of read length collapses to the whole read).
    Only the sequence cleaning differs: a character filter loop instead of the
    production ``re.sub``.
    """
    s = "".join(ch for ch in (seq or "").upper() if ch in "ACGT")
    n = len(s)
    if n < int(min_len):
        return (0, 0.0, "", "")
    thr = float(min_frac)
    best_len = 0
    best_frac = 0.0
    best_base = ""
    best_ij = (0, 0)
    for base in ("A", "T"):
        left = 0
        n_base = 0
        for right in range(n):
            if s[right] == base:
                n_base += 1
            while left <= right and (float(n_base) / float(right - left + 1)) < thr:
                if s[left] == base:
                    n_base -= 1
                left += 1
            cur_len = right - left + 1
            if cur_len >= int(min_len) and n_base > 0:
                frac = float(n_base) / float(cur_len)
                if cur_len > best_len or (cur_len == best_len and frac > best_frac):
                    best_len = cur_len
                    best_frac = frac
                    best_base = base
                    best_ij = (left, right + 1)
    if best_len <= 0:
        return (0, 0.0, "", "")
    span = s[best_ij[0] : best_ij[1]]
    if best_len >= 140 or best_len >= n - 2:
        best_len = n
        best_frac = float(span.count(best_base)) / float(len(span)) if span else best_frac
        span = s
    return (int(best_len), float(best_frac), best_base, span)


def _random_cigar_string(rng: random.Random) -> str:
    """Build a valid-ish random CIGAR with optional H/S groups at each end."""
    left_group = rng.choice(["", "5S", "3H4S", "4S", "2H"])
    right_group = rng.choice(["", "9S", "4S5H", "7S", "3H"])
    interior_ops = ["M", "I", "D", "N", "=", "X"]
    interior: list[str] = []
    last = None
    for _ in range(rng.randint(1, 6)):
        op = rng.choice(interior_ops)
        if op == last:
            continue
        last = op
        interior.append(f"{rng.randint(1, 60)}{op}")
    return left_group + "".join(interior) + right_group


# ─────────────────────────────────────────────────────────────────────────────
# _cluster_sorted_positions (NumPy diff/split) vs scalar loop
# ─────────────────────────────────────────────────────────────────────────────

_CLUSTER_EDGE_CASES = [
    ([], 100),  # empty list
    ([500], 100),  # single position
    ([100, 150, 200], 100),  # all within max_gap
    ([100, 200], 100),  # gap exactly at boundary stays together
    ([100, 202], 100),  # one over the gap splits
    ([100, 100, 100], 0),  # duplicate coordinates
    ([100, 100, 200], -5),  # negative max_gap clamped to 0
    ([100, 101], 0),  # zero gap splits adjacent positions
    ([42], -1),  # single position, clamped negative gap
    ([0, 1, 2, 1000, 1001, 2000, 100000], 500),
    (list(range(0, 1000, 10)), 25),  # dense run split at 25 bp
    (list(range(0, 1000, 10)), 9),
    (list(range(0, 1000, 10)), 10),
    ([1, 1, 1, 5, 6, 6, 100, 100, 100], 4),
    ([-(10**6), 0, 10**6], 1_000_000),  # large gaps
    ([0, 0, 10**9, 10**9 + 1], 10**9),
]


class TestClusterSortedPositionsParity:
    @pytest.mark.parametrize("positions,max_gap_bp", _CLUSTER_EDGE_CASES)
    def test_edge_case_parity(self, positions: list[int], max_gap_bp: int) -> None:
        assert _cluster_sorted_positions(positions, max_gap_bp=max_gap_bp) == _ref_cluster_sorted_positions(
            positions, max_gap_bp
        )

    def test_random_parity(self) -> None:
        rng = random.Random(20260829)
        for _ in range(300):
            n = rng.randint(0, 40)
            positions = sorted(rng.randint(rng.randint(-500, 0), rng.randint(0, 500)) for _ in range(n))
            if positions and rng.random() < 0.5:
                positions = sorted(positions + [rng.choice(positions) for _ in range(rng.randint(1, 5))])
            max_gap_bp = rng.choice([-10, 0, 1, 5, 17, 100, 500])
            assert _cluster_sorted_positions(positions, max_gap_bp=max_gap_bp) == _ref_cluster_sorted_positions(
                positions, max_gap_bp
            )


# ─────────────────────────────────────────────────────────────────────────────
# CIGAR spans / query length / soft-clip boundaries (regex) vs scalar scan
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEX_CIGARS = [
    "100M",  # no clips
    "10S90M",  # left soft clip only
    "90M10S",  # right soft clip only
    "10S50M10S",  # both-sided soft clip
    "6H5S15M2I10M3D20M100S",  # hard+soft, insertion, deletion, long right clip
    "10S60M10N40M10S",  # intron skip between clips
    "5S90M5H",  # left soft clip, right hard clip
    "5H95M5S",  # left hard clip, right soft clip
    "3S50M20X10=10S",  # X and = operators
    "10S10I30M4S",  # soft clips wrapping an insertion
    "5H100M10H",  # hard clips only
    "20S80M",  # one-sided left
    "80M20S",  # one-sided right
    "10S",  # fully soft-clipped
    "10S3H",  # soft then hard on the right edge
]


class TestCigarParity:
    @pytest.mark.parametrize("cigar", _COMPLEX_CIGARS + [None])
    def test_alignment_spans_parity(self, cigar: str) -> None:
        assert _cigar_alignment_spans(cigar) == _ref_cigar_alignment_spans(cigar)  # type: ignore[arg-type]

    @pytest.mark.parametrize("cigar", _COMPLEX_CIGARS + [None])
    def test_query_len_parity(self, cigar: str) -> None:
        assert _cigar_query_len(cigar) == _ref_cigar_query_len(cigar)  # type: ignore[arg-type]

    @pytest.mark.parametrize("cigar", _COMPLEX_CIGARS + [None])
    def test_soft_clip_boundaries_parity(self, cigar: str) -> None:
        read = _fake_read(str(cigar or ""), query_seq="Q" * _ref_cigar_query_len(cigar or ""))
        left, right = _ref_soft_clip_bounds(str(cigar or ""))
        expected = [("L", left)] * bool(left) + [("R", right)] * bool(right)
        assert sorted(_collect_soft_clips(read, min_clip_len=1)) == sorted(expected)

    @pytest.mark.parametrize("cigar", _COMPLEX_CIGARS + [None])
    def test_longest_soft_clip_parity(self, cigar: str) -> None:
        cigar_s = str(cigar or "")
        qlen = _ref_cigar_query_len(cigar_s)
        query = ("ATCG" * (1 + qlen // 4))[:qlen]
        read = _fake_read(cigar_s, query_seq=query, ref_start=1000, ref_end=2000)
        left, right = _ref_soft_clip_bounds(cigar_s)
        if left <= 0 and right <= 0:
            assert _longest_soft_clip_from_read(read) == ("", 0, 0, "")
        elif left >= right:
            assert _longest_soft_clip_from_read(read) == ("L", left, 1001, query[:left])
        else:
            assert _longest_soft_clip_from_read(read) == ("R", right, 2000, query[-right:])

    def test_multialignment_segments_parity(self) -> None:
        segments = ["10S90M5S", "5S80M100S", "6H12S48M4S12H"]
        for seg in segments:
            assert _cigar_alignment_spans(seg) == _ref_cigar_alignment_spans(seg)
            assert _cigar_query_len(seg) == _ref_cigar_query_len(seg)
        # Summing spans across chimeric segments is the whole-read alignment budget.
        spans = [_cigar_alignment_spans(seg) for seg in segments]
        assert sum(sp[0] for sp in spans) == 90 + 80 + 48
        assert sum(sp[1] for sp in spans) == 90 + 80 + 48
        assert sum(sp[2] for sp in spans) == 90 + 80 + 48
        # Per-segment soft-clip boundaries are preserved for supplementary records.
        # "6H12S48M4S12H" records no clip: H occupies the first/last CIGAR slot,
        # mirroring _collect_soft_clips which only inspects cigartuples[0]/[-1].
        for seg, (left, right) in zip(segments, [(10, 5), (5, 100), (0, 0)]):
            assert _ref_soft_clip_bounds(seg) == (left, right)
            read = _fake_read(seg, query_seq="Q" * _ref_cigar_query_len(seg))
            expected = [("L", left)] * bool(left) + [("R", right)] * bool(right)
            assert sorted(_collect_soft_clips(read, min_clip_len=1)) == sorted(expected)

    def test_random_cigar_parity(self) -> None:
        rng = random.Random(20260829)
        for _ in range(200):
            cigar = _random_cigar_string(rng)
            assert _cigar_alignment_spans(cigar) == _ref_cigar_alignment_spans(cigar)
            assert _cigar_query_len(cigar) == _ref_cigar_query_len(cigar)
            read = _fake_read(cigar, query_seq="Q" * _ref_cigar_query_len(cigar))
            left, right = _ref_soft_clip_bounds(cigar)
            expected = [("L", left)] * bool(left) + [("R", right)] * bool(right)
            assert sorted(_collect_soft_clips(read, min_clip_len=1)) == sorted(expected)


# ─────────────────────────────────────────────────────────────────────────────
# _poly_at_stats (regex run scan) vs scalar counting loop
# ─────────────────────────────────────────────────────────────────────────────


class TestPolyAtStatsParity:
    @pytest.mark.parametrize(
        "seq",
        [
            "",  # empty
            "AAAA",  # pure polyA
            "TTTTTT",  # pure polyT
            "A" * 40,
            "T" * 40,
            "AAAGAAA",  # mixed run, A dominant
            "TTTAA",  # T dominant
            "ATAT",  # tie -> A wins
            "GCGCGCGC",  # no A/T
            "aaaatttt",  # lowercase uppercased
            "ANANATTTT",  # non-ACGT chars counted in length only
        ],
    )
    def test_explicit_parity(self, seq: str) -> None:
        assert _poly_at_stats(seq) == _ref_poly_at_stats(seq)

    def test_random_parity(self) -> None:
        rng = random.Random(20260829)
        for _ in range(500):
            length = rng.randint(0, 60)
            alphabet = "ACGTN" if rng.random() < 0.2 else "ACGT"
            seq = "".join(rng.choice(alphabet) for _ in range(length))
            if rng.random() < 0.3:
                seq = seq.lower()
            assert _poly_at_stats(seq) == _ref_poly_at_stats(seq)


# ─────────────────────────────────────────────────────────────────────────────
# _longest_poly_at_span (regex + two-pointer) vs scalar two-pointer transcription
# ─────────────────────────────────────────────────────────────────────────────


class TestLongestPolyAtSpanParity:
    @pytest.mark.parametrize(
        "seq",
        [
            "A" * 30,  # pure polyA
            "T" * 30,  # pure polyT
            "A" * 80,  # pure polyA over the full-read override
            "T" * 80,
            "A" * 25,  # exactly min_len
            "A" * 20,  # below min_len
            "",  # empty
            "AG" * 40,  # 50% A -> below purity threshold
            "A" * 10 + "C" + "A" * 10 + "G" + "A" * 10 + "C",  # high-purity noisy
        ],
    )
    def test_explicit_parity(self, seq: str) -> None:
        assert _longest_poly_at_span(seq) == _ref_longest_poly_at_span(seq)

    @pytest.mark.parametrize(
        "seq,min_frac,min_len",
        [
            ("A" * 20 + "CC" + "A" * 20, 0.9, 8),  # noisy block between polyA runs
            ("A" * 15 + "G" + "A" * 15, 0.95, 8),  # ~96.8% A window
            ("T" * 12 + "AA" + "T" * 12, 0.9, 8),  # noisy A block inside polyT
            ("ACGT" * 10 + "A" * 30, 0.9, 8),  # noisy head + clean tail
            ("A" * 30 + "CGT" * 15, 0.9, 8),  # clean head + noisy tail
            ("N" * 20, 0.9, 8),  # all-N strips to empty
        ],
    )
    def test_noisy_sequences_parity(self, seq: str, min_frac: float, min_len: int) -> None:
        assert _longest_poly_at_span(seq, min_frac=min_frac, min_len=min_len) == _ref_longest_poly_at_span(
            seq, min_frac=min_frac, min_len=min_len
        )
        # Semantic invariants independent of the transcription under test.
        length, frac, base, span = _longest_poly_at_span(seq, min_frac=min_frac, min_len=min_len)
        if length > 0:
            assert base in {"A", "T"}
            assert length == len(span)
            assert length >= min_len
            assert frac >= min_frac

    def test_random_parity(self) -> None:
        rng = random.Random(20260829)
        for _ in range(200):
            length = rng.randint(0, 90)
            seq = "".join(rng.choice("ACGTNacgtn") for _ in range(length))
            min_len = rng.choice([1, 5, 8, 25])
            min_frac = rng.choice([0.5, 0.8, 0.9, 0.95])
            assert _longest_poly_at_span(seq, min_frac=min_frac, min_len=min_len) == _ref_longest_poly_at_span(
                seq, min_frac=min_frac, min_len=min_len
            )
