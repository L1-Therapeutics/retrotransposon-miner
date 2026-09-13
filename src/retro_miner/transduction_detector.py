"""L1 3' Transduction Detector.

Detects and classifies the 3' unique-coverage transductions of L1 retrotransposons
that are copied from the ``3'`` flank of a source (donor) L1 locus and dragged
along with the L1 body during a new insertion.  Each transduced segment is
flanked ``5'`` by the L1 poly(A) tail, so the engine trims the MEI body and the
poly(A) tail from an assembled breakpoint contig, keeps the downstream sequence,
and gates it on sequence complexity (Shannon entropy) before classifying:

* ``3_PRIME_PARTNERED`` - transduction attached to an intact L1HS ``3'`` UTR
  (identifiable polyadenylation signal upstream of the L1 end), so the segment
  can be linked back to a specific source locus.
* ``3_PRIME_ORPHAN`` - transduction dragged by a ``5'``-truncated L1 that lost
  its ``3'`` UTR; the segment has no identifiable source locus.
* ``NONE`` - no transduction (poly(A) extension only, low-complexity sequence,
  or truncated assembly).

Input is the De Bruijn micro-assembly breakpoint unitig produced by
:mod:`retro_miner.local_assembly` (``AssembledContig.sequence``) together with
the coordinate of the MEI ``3'`` end within the contig
(``mei_alignment_end``).

Literature Anchor: Tubio et al. (2015) LINE-1 3' transduction and its role in
corpome rearrangements, Science 347:1260408; Ewing & Kazazian (2010).
"""

from __future__ import annotations

from dataclasses import dataclass

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


def detect_3prime_transduction(
    assembled_contig: str,
    mei_alignment_end: int,
    min_transduction_len: int = MIN_TRANSDUCTION_LEN,
) -> TransductionResult:
    """Detect a 3' transduction in an assembled breakpoint unitig.

    Args:
        assembled_contig: De Bruijn-assembled breakpoint unitig covering the
            MEI 3' junction: ``[L1 body][poly(A) tail][transduced flank]``.
        mei_alignment_end: 0-based index in *assembled_contig* where the MEI
            (L1) body alignment ends; bases at ``>= mei_alignment_end`` are
            candidate downstream (post-L1) sequence.
        min_transduction_len: Minimum length in bp of the downstream sequence
            before it is considered a transduction.

    Returns:
        A :class:`TransductionResult` classifying the detected transduction.
    """
    if not assembled_contig:
        return _no_transduction()

    contig = assembled_contig.upper()
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