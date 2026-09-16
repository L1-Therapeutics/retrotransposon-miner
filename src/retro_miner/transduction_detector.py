"""L1 3' Transduction Detector.

Two complementary entry points operate on complementary inputs:

1. Read-level donor profiling (:func:`detect_3prime_transduction`), the TPRT /
   read-through detector.  Active L1 elements frequently bypass their weak
   canonical 3' polyadenylation signal, causing RNA polymerase II to transcribe
   unique genomic sequence downstream of the element; when the chimeric
   transcript integrates at a new locus it carries a unique "transduction
   fingerprint" that identifies the exact source (donor) element (Tubio et al.,
   2014, Science 345(6196):1251343).  Given read-level input:

   * canonical L1HS/Alu consensus sequence matches are stripped from the reads,
   * the ``3'`` poly(A)/poly(T) tract is trimmed,
   * the downstream sequence is kept only when it is highly complex
     (Shannon entropy ``H > 1.8``),
   * the retained unique sequence is seeded against a k-mer index of the
     reference genome to localise candidate donor loci.

2. Contig-level classification (:func:`detect_3prime_transduction_from_contig`).
   Given a De Bruijn breakpoint unitig structured ``[L1 body][poly(A) tail]
   [transduced flank]`` plus the coordinate where the L1 body alignment ends,
   classifies the ``3'`` sequence as ``3_PRIME_PARTNERED`` (identifiable source
   L1 via an L1HS 3' UTR polyadenylation signal), ``3_PRIME_ORPHAN`` (the
   source L1 lost its 3' UTR), or ``NONE``.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from dataclasses import dataclass
from functools import lru_cache

from pyfaidx import Fasta

from .tsd_refiner import calculate_sequence_entropy

POLY_TAIL_MIN_RUN = 5
POLY_TAIL_ENTROPY_THRESHOLD = 1.2
TRANSDUCTION_ENTROPY_THRESHOLD = 1.5
MIN_TRANSDUCTION_LEN = 20

_L1HS_3UTR_LOOKBACK_BP = 80
_L1HS_POLYA_SIGNALS = ("AATAAA", "TTTATT")

_CONF_PARTNERED_WITH_TAIL = 0.95
_CONF_ORPHAN_WITH_TAIL = 0.85
_CONF_PARTNERED_NO_TAIL = 0.80
_CONF_ORPHAN_NO_TAIL = 0.70

_TRANS_TYPE_PARTNERED = "3_PRIME_PARTNERED"
_TRANS_TYPE_ORPHAN = "3_PRIME_ORPHAN"
_TRANS_TYPE_NONE = "NONE"

# --- read-level donor profiler -------------------------------------------------
TAIL_MIN_LENGTH = 20
TAIL_MIN_ENTROPY = 1.8
SEED_K = 20
MIN_SEED_K = 8
DONOR_MAP_MIN_IDENTITY = 0.90
MAX_DONOR_CANDIDATES = 200
MEI_CONSENSUS_MIN_IDENTITY = 0.85

# Simplified consensus landmarks used to strip MEI sequence from read-level
# input before tail extraction.  The L1HS fragment is the polyadenylation
# signal-flanking 3' UTR core exercised by the transduction test fixtures; the
# Alu fragment is the conserved GC-rich SINE 5' head.  Both should be replaced
# with a curated consensus (e.g. RepeatMasker ``.cons`` output) before
# production use.
_L1HS_3UTR_CORE_CONSENSUS = "TGAGGGGCCTGTGCTTCTCCTTCAAGGAAATAAAATTAACCTTGAGCT"
_ALU_5PRIME_CONSENSUS = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
MEI_CONSENSUS_FRAGMENTS = (_L1HS_3UTR_CORE_CONSENSUS, _ALU_5PRIME_CONSENSUS)

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass(frozen=True)
class TransductionResult:
    """Outcome of 3' transduction detection on one assembled breakpoint contig.

    Attributes:
        has_transduction: Whether a transduction sequence passed all gates.
        transduction_seq: Transduced (source-flank) sequence beyond the ``3'``
            poly(A) tail, uppercased; ``""`` when no transduction.
        transduction_length: Length in bp of *transduction_seq*.
        transduction_type: One of ``3_PRIME_PARTNERED``, ``3_PRIME_ORPHAN`` or
            ``NONE``.
        confidence_score: 0..1 confidence in the reported classification.
    """

    has_transduction: bool
    transduction_seq: str
    transduction_length: int
    transduction_type: str
    confidence_score: float


@dataclass(frozen=True)
class DonorLocus:
    """A candidate genomic source (donor) interval for a transduced segment.

    Attributes:
        chrom: Reference contig name.
        donor_start: 0-based half-open start of the donor interval.
        donor_end: 0-based half-open end of the donor interval.
        strand: ``"+"`` if the tail mapped in the reference orientation,
            ``"-"`` if it mapped as the reverse complement.
        identity: Fraction of the tail matching the reference window.
        mismatches: Number of mismatching bases within the donor interval.
    """

    chrom: str
    donor_start: int
    donor_end: int
    strand: str
    identity: float
    mismatches: int


@dataclass(frozen=True)
class TransductionCall:
    """Outcome of read-level 3' transduction donor profiling.

    Attributes:
        has_transduction: True when at least one unique, high-complexity,
            non-MEI tail was recovered beyond a poly(A)/poly(T) tract.
        tail_sequence: Best recovered unique tail (uppercased).
        tail_length: Length in bp of *tail_sequence*.
        entropy: Shannon entropy of *tail_sequence*.
        poly_a_trimmed_bp: Total poly(A)/poly(T) bases trimmed in the profiled
            reads.
        source_read_count: Number of input reads that produced a qualifying
            tail.
        donor_loci: Candidate donor loci, sorted by decreasing identity.
        best_donor_locus: Highest-scoring donor locus, or ``None``.
    """

    has_transduction: bool
    tail_sequence: str
    tail_length: int
    entropy: float
    poly_a_trimmed_bp: int
    source_read_count: int
    donor_loci: tuple[DonorLocus, ...]
    best_donor_locus: DonorLocus | None


ReferenceIndex = tuple[tuple[tuple[str, str], ...], dict[str, list[tuple[int, int]]]]


def _leading_homopolymer_run(seq: str) -> tuple[str, int]:
    """Return ``(base, run_length)`` of the leading homopolymer in *seq*.

    Only a pure ``A`` or pure ``T`` run qualifies as a candidate poly(A) tail;
    any other leading base returns an empty tail.
    """
    if not seq or seq[0] not in {"A", "T"}:
        return "", 0
    base = seq[0]
    run = 0
    for char in seq:
        if char == base:
            run += 1
        else:
            break
    return base, run


def _has_l1hs_3utr_signal(body: str) -> bool:
    """True if the tail of the L1 body carries an L1HS 3' UTR polyadenylation signal.

    Looks for ``AATAAA`` (plus strand) or ``TTTATT`` (minus strand read on the
    contig) within the last ``_L1HS_3UTR_LOOKBACK_BP`` bases of the MEI body,
    i.e. immediately upstream of the L1 3' boundary.
    """
    window = body[-_L1HS_3UTR_LOOKBACK_BP:]
    return any(signal in window for signal in _L1HS_POLYA_SIGNALS)


def _no_transduction() -> TransductionResult:
    return TransductionResult(
        has_transduction=False,
        transduction_seq="",
        transduction_length=0,
        transduction_type=_TRANS_TYPE_NONE,
        confidence_score=0.0,
    )


def detect_3prime_transduction_from_contig(
    assembled_contig: str,
    mei_alignment_end: int | None = None,
    min_transduction_len: int = MIN_TRANSDUCTION_LEN,
) -> TransductionResult:
    """Detect a 3' transduction in an assembled breakpoint unitig.

    Args:
        assembled_contig: De Bruijn-assembled breakpoint unitig covering the
            MEI 3' junction: ``[L1 body][poly(A) tail][transduced flank]``.
        mei_alignment_end: 0-based index in *assembled_contig* where the MEI
            (L1) body alignment ends; bases at ``>= mei_alignment_end`` are
            candidate downstream (post-L1) sequence.  When omitted, the
            junction is auto-detected as the last pure ``A``/``T`` homopolymer
            run of at least ``POLY_TAIL_MIN_RUN`` bases; if none is found the
            call is conservatively ``NONE``.
        min_transduction_len: Minimum length in bp of the downstream sequence
            before it is considered a transduction.

    Returns:
        A :class:`TransductionResult` classifying the detected transduction.
    """
    if not assembled_contig:
        return _no_transduction()

    contig = assembled_contig.upper()
    if mei_alignment_end is None:
        auto_end = _auto_mei_alignment_end(contig)
        if auto_end is None:
            return _no_transduction()
        mei_alignment_end = auto_end
    end = max(0, min(int(mei_alignment_end), len(contig)))
    body = contig[:end]
    downstream = contig[end:]

    if len(downstream) < min_transduction_len:
        return _no_transduction()

    _, tail_len = _leading_homopolymer_run(downstream)
    if 0 < tail_len < POLY_TAIL_MIN_RUN:
        tail_len = 0
    core = downstream[tail_len:] if tail_len else downstream

    if len(core) < min_transduction_len:
        return _no_transduction()
    if calculate_sequence_entropy(core) < TRANSDUCTION_ENTROPY_THRESHOLD:
        return _no_transduction()

    partnered = _has_l1hs_3utr_signal(body)
    has_tail = int(bool(tail_len))

    if partnered:
        trans_type = _TRANS_TYPE_PARTNERED
        confidence = _CONF_PARTNERED_WITH_TAIL if has_tail else _CONF_PARTNERED_NO_TAIL
    else:
        trans_type = _TRANS_TYPE_ORPHAN
        confidence = _CONF_ORPHAN_WITH_TAIL if has_tail else _CONF_ORPHAN_NO_TAIL

    return TransductionResult(
        has_transduction=True,
        transduction_seq=core,
        transduction_length=len(core),
        transduction_type=trans_type,
        confidence_score=confidence,
    )


def _homopolymer_runs(
    seq: str, min_run: int = POLY_TAIL_MIN_RUN
) -> Iterator[tuple[int, int]]:
    """Yield ``(start, length)`` for pure ``A``/``T`` runs of length >= *min_run*."""
    seq_upper = seq.upper()
    i = 0
    while i < len(seq_upper):
        j = i + 1
        while j < len(seq_upper) and seq_upper[j] == seq_upper[i]:
            j += 1
        if j - i >= min_run and seq_upper[i] in "AT":
            yield i, j - i
        i = j


def _auto_mei_alignment_end(contig: str) -> int | None:
    """Best-effort MEI 3' junction: the last qualifying poly(A)/poly(T) run."""
    last_start: int | None = None
    for start, _run_len in _homopolymer_runs(contig):
        last_start = start
    return last_start


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _strip_mei_consensus_spans(
    seq: str, min_identity: float = MEI_CONSENSUS_MIN_IDENTITY
) -> str:
    """Remove spans matching the canonical L1HS/Alu consensus landmarks.

    Masking (rather than read rejection) preserves the transduced tail that a
    chimeric read carries downstream of its consensus-matching ``5'`` body.
    """
    if not seq:
        return ""
    seq_upper = seq.upper()
    mask = [False] * len(seq_upper)
    for fragment in MEI_CONSENSUS_FRAGMENTS:
        wlen = len(fragment)
        if wlen > len(seq_upper):
            continue
        probes = (fragment, _reverse_complement(fragment))
        for probe in probes:
            for start in range(len(seq_upper) - wlen + 1):
                window = seq_upper[start : start + wlen]
                matched = sum(a == b for a, b in zip(window, probe, strict=False))
                if matched / wlen >= min_identity:
                    mask[start : start + wlen] = [True] * wlen
    return "".join(ch for i, ch in enumerate(seq_upper) if not mask[i])


def _extract_transduction_tails(
    reads: list[str],
) -> tuple[list[str], int, int]:
    """Extract unique high-complexity tails extending beyond a poly(A) tract.

    For each read the poly(A)/poly(T) ``5'`` terminating run whose suffix is at
    least ``TAIL_MIN_LENGTH`` bp long with entropy ``> TAIL_MIN_ENTROPY`` is
    located in both orientations; the highest-entropy suffix wins.

    Returns ``(tails, covered_reads, trimmed_bp)`` where *tails* are deduplicated
    and order-preserving.
    """
    tails: list[str] = []
    covered_reads = 0
    trimmed_bp = 0
    seen: set[str] = set()
    for read in reads:
        best_tail: str | None = None
        best_trimmed = 0
        best_score = -1.0
        for probe in (read, _reverse_complement(read)):
            for start, run_len in _homopolymer_runs(probe):
                tail = probe[start + run_len :]
                if len(tail) < TAIL_MIN_LENGTH:
                    continue
                entropy = calculate_sequence_entropy(tail)
                if entropy <= TAIL_MIN_ENTROPY:
                    continue
                if entropy > best_score or (entropy == best_score and len(tail) > len(best_tail or "")):
                    best_tail = tail
                    best_trimmed = run_len
                    best_score = entropy
        if best_tail is not None:
            if best_tail not in seen:
                seen.add(best_tail)
                tails.append(best_tail)
            trimmed_bp += best_trimmed
            covered_reads += 1
    return tails, covered_reads, trimmed_bp


def _dedupe_sort_loci(results: list[DonorLocus]) -> list[DonorLocus]:
    best: dict[tuple[str, int, int, str], DonorLocus] = {}
    for locus in results:
        key = (locus.chrom, locus.donor_start, locus.donor_end, locus.strand)
        previous = best.get(key)
        if previous is None or locus.identity > previous.identity:
            best[key] = locus
    return sorted(
        best.values(),
        key=lambda locus: (-locus.identity, locus.mismatches, locus.chrom, locus.donor_start),
    )


def _load_reference_index(ref_fasta: str) -> ReferenceIndex:
    """Load (and cache) a 20-mer seed index over every contig of *ref_fasta*."""
    stat = os.stat(ref_fasta)
    return _build_reference_index(ref_fasta, stat.st_mtime_ns, stat.st_size)


@lru_cache(maxsize=8)
def _build_reference_index(path: str, mtime_ns: int, size: int) -> ReferenceIndex:
    contigs: list[tuple[str, str]] = []
    index: dict[str, list[tuple[int, int]]] = {}
    fasta = Fasta(path, as_raw=True)
    try:
        for chrom in fasta.keys():
            seq = str(fasta[chrom][:]).upper()
            contig_index = len(contigs)
            contigs.append((chrom, seq))
            if len(seq) < SEED_K:
                continue
            for pos in range(len(seq) - SEED_K + 1):
                index.setdefault(seq[pos : pos + SEED_K], []).append(
                    (contig_index, pos)
                )
    finally:
        fasta.close()
    return tuple(contigs), index


def _map_tail(tail: str, index: ReferenceIndex) -> list[DonorLocus]:
    """Seed-and-extend search of *tail* in both orientations against *index*."""
    contigs, kmers = index
    sequence = tail.upper()
    results: list[DonorLocus] = []
    for probe, strand in ((sequence, "+"), (_reverse_complement(sequence), "-")):
        k = min(SEED_K, len(probe))
        if k < MIN_SEED_K:
            continue
        windows: set[tuple[int, int]] = set()
        for offset in range(len(probe) - k + 1):
            seed = probe[offset : offset + k]
            for contig_index, pos in kmers.get(seed, ()):
                win_start = pos - offset
                win_end = win_start + len(probe)
                if 0 <= win_start <= win_end and win_end <= len(contigs[contig_index][1]):
                    windows.add((contig_index, win_start))
        verified = 0
        for contig_index, win_start in windows:
            if verified >= MAX_DONOR_CANDIDATES:
                break
            verified += 1
            ref_window = contigs[contig_index][1][win_start : win_start + len(probe)]
            mismatches = sum(a != b for a, b in zip(ref_window, probe, strict=False))
            identity = (len(probe) - mismatches) / len(probe)
            if identity >= DONOR_MAP_MIN_IDENTITY:
                end = win_start + len(probe)
                results.append(
                    DonorLocus(
                        chrom=contigs[contig_index][0],
                        donor_start=win_start,
                        donor_end=end,
                        strand=strand,
                        identity=round(identity, 4),
                        mismatches=mismatches,
                    )
                )
    return results


def detect_3prime_transduction(
    unmapped_mates: list[str],
    non_mei_clips: list[str],
    ref_fasta: str,
) -> TransductionCall:
    """Detect and localise 3' read-through transductions from read sequences.

    Chimeric reads produced by an L1 insertion that bypassed its canonical
    polyadenylation signal carry ``[L1 3' UTR][poly(A) tract][unique genomic
    tail]``.  The tail is the unique "transduction fingerprint" that maps back
    to the donor (source) element (Tubio et al., 2014, Science 345(6196):1251343).

    Args:
        unmapped_mates: Read sequences (unmapped mates / paired evidence) that
            failed to map against the reference.
        non_mei_clips: Soft-clip sequences already classified as non-MEI at the
            insertion junction.
        ref_fasta: Reference genome FASTA used to localise candidate donor
            loci.  Must exist; when it is missing the call is still made but
            without donor mapping.

    Returns:
        A :class:`TransductionCall` with the recovered unique tail and any
        candidate donor loci.
    """
    reads = _normalize_read_inputs(unmapped_mates, non_mei_clips)
    if not reads:
        return _no_transduction_call()
    stripped = [seq for seq in (_strip_mei_consensus_spans(seq) for seq in reads) if seq]
    if not stripped:
        return _no_transduction_call()
    tails, covered_reads, trimmed_bp = _extract_transduction_tails(stripped)
    if not tails:
        return _no_transduction_call(poly_a_trimmed_bp=trimmed_bp)

    donor_loci: list[DonorLocus] = []
    ref_path = str(ref_fasta)
    if ref_path and os.path.isfile(ref_path):
        index = _load_reference_index(ref_path)
        for tail in tails:
            donor_loci.extend(_map_tail(tail, index))
    donor_loci = _dedupe_sort_loci(donor_loci)

    best_tail = max(tails, key=lambda tail: (calculate_sequence_entropy(tail), len(tail)))
    best_entropy = calculate_sequence_entropy(best_tail)
    return TransductionCall(
        has_transduction=True,
        tail_sequence=best_tail,
        tail_length=len(best_tail),
        entropy=best_entropy,
        poly_a_trimmed_bp=trimmed_bp,
        source_read_count=covered_reads,
        donor_loci=tuple(donor_loci),
        best_donor_locus=donor_loci[0] if donor_loci else None,
    )


def _normalize_read_inputs(
    unmapped_mates: list[str], non_mei_clips: list[str]
) -> list[str]:
    """Uppercase, deduplicate and drop empties from the read-level inputs."""
    seen: set[str] = set()
    reads: list[str] = []
    for seq in (list(unmapped_mates) + list(non_mei_clips)):
        if not isinstance(seq, str) or not seq:
            continue
        upper = seq.upper()
        if upper not in seen:
            seen.add(upper)
            reads.append(upper)
    return reads


def _no_transduction_call(poly_a_trimmed_bp: int = 0) -> TransductionCall:
    return TransductionCall(
        has_transduction=False,
        tail_sequence="",
        tail_length=0,
        entropy=0.0,
        poly_a_trimmed_bp=poly_a_trimmed_bp,
        source_read_count=0,
        donor_loci=(),
        best_donor_locus=None,
    )