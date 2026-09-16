"""L1 5' Inversion (Twin-Priming) Breakpoint Analyzer.

Detects the characteristic ~20% of de novo LINE-1 insertions whose 5' terminus
carries an inverted, truncated copy of the L1 5' end attached to an upright
3' body.  Such structures arise when target-site second-strand cleavage primes
reverse transcription internally on the L1 RNA before plus-strand cDNA
synthesis completes, so that the cDNA is copied in the reverse orientation
upstream of the priming point (twin-priming; Ostertag & Kazazian, 2001, Genome
Research 11(4):459-465).

Structure modeled here (insertion plus-strand, 5' -> 3'):

    [5' flank] -- RC(C[a:b]) -- [junction micro-homology] -- C[p:] -- [poly(A)] -- [3' flank]

where ``C`` is the L1HS reference consensus, ``RC(C[a:b])`` is the reverse
complement of a 5'-proximal consensus interval (the inverted 5' copy), and
``C[p:]`` is the upright L1 body starting at internal consensus position
``p``.  The strand switch happened through a short micro-homology overlap of
``m`` bp where ``C[p:p+m] == RC(C[a:a+m])``.

Standard split-read callers miss these events because the 5' soft-clip maps in
the reverse orientation relative to the 3' clip.  This module therefore scans
5' soft-clips for *internal* strand switches: a prefix that aligns to the
consensus in reverse-complement orientation joined to a suffix that aligns
forward, with the micro-homology recovered from the consensus context at the
switch point.

Mirroring the reverse-complement, the copy's exact analytical picture is: the
linear junction reads ``RC(C[a+m:b]) ++ C[p:]`` -- the first ``m`` bases
``C[a:a+m]`` of the inverted interval (whose reverse complement is the
micro-homology) are withheld from the read and "absorbed" into the shared
block ``C[p:p+m]``.

Per-clip search:

* For every split position the clip is decomposed into ``prefix | suffix``.
* Both arms are aligned to ``C`` in the forward and reverse-complement
  directions.  Because the insertion can sit on either genomic strand, the
  inverted copy may map in either orientation; the four valid
  orientation assignments are each tested.
* For a tested assignment the inverted element is known up to its withheld
  head ``[inv_off - m, inv_off)`` and the upright element contributes head
  ``[upt_off, upt_off + m)``.  The micro-homology ``m`` is the longest
  ``m <= MAX_OVERLAP_MH`` with ``RC(C[inv_off-m:inv_off]) == C[upt_off:upt_off+m]``.
* ``inverted_length`` recovers the full inverted-copy consensus span
  ``b - a`` (the returned-overlap ``m`` plus the aligned arm length);
  ``inversion_breakpoint`` is the consensus position where the upright body
  (and therefore the strand switch) begins.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ._utils import _iter_fasta_records

MIN_CLIP_LEN = 20
MIN_SEG_LEN = 12
MIN_IDENTITY = 0.85
MIN_OVERLAP_MH = 2
MAX_OVERLAP_MH = 8

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass(frozen=True)
class TwinPrimingCall:
    """Outcome of twin-priming 5' inversion breakpoint analysis.

    Attributes:
        is_twin_primed: True when a 5' soft-clip exhibited a reverse-complement
            strand switch with a qualifying micro-homology junction.
        inversion_breakpoint: Consensus position ``p`` (0-based) where the
            upright L1 body begins after the inverted copy, i.e. the point of
            strand switching; ``-1`` when no inversion was detected.
        inverted_length: Length in bp of the inverted 5' copy
            ``RC(C[a:b])`` recovered from the clip; ``0`` when none detected.
        complementary_overlap_bp: Length in bp of the shared micro-homology at
            the strand-switch junction; ``0`` when none detected.
    """

    is_twin_primed: bool
    inversion_breakpoint: int = -1
    inverted_length: int = 0
    complementary_overlap_bp: int = 0


@dataclass(frozen=True)
class _Alignment:
    offset: int
    identity: float


@dataclass(frozen=True)
class _SwitchCandidate:
    score: float
    inverted_off: int
    breakpoint: int
    overlap_bp: int
    inverted_length: int


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _load_consensus(fasta_path: str) -> str:
    """Return the first non-empty L1 consensus sequence from *fasta_path*."""
    path = Path(fasta_path)
    if not path.is_file():
        return ""
    for _name, seq in _iter_fasta_records(path):
        if seq:
            return seq
    return ""


def _best_alignment(probe: str, consensus: str) -> _Alignment | None:
    """Best forward alignment of *probe* within *consensus*.

    Returns the 0-based consensus ``offset`` and identity of the highest
    scoring window whose identity reaches ``MIN_IDENTITY``, or ``None``.
    Windows are scanned with early abort once mismatches exceed the identity
    budget, so exact-match probes are validated in O(n) without full window
    evaluation.
    """
    m = len(probe)
    n = len(consensus)
    if m < MIN_SEG_LEN or n < m:
        return None
    allowed_mism = int(m * (1.0 - MIN_IDENTITY))
    best: _Alignment | None = None
    for offset in range(n - m + 1):
        window = consensus[offset : offset + m]
        mismatches = 0
        for i, ch in enumerate(probe):
            if ch != window[i]:
                mismatches += 1
                if mismatches > allowed_mism:
                    break
        else:
            identity = (m - mismatches) / m
            if identity >= MIN_IDENTITY:
                if best is None or identity > best.identity:
                    best = _Alignment(offset=offset, identity=identity)
    return best


def _longest_overlap(inv_off: int, upt_off: int, consensus: str) -> int:
    """Longest micro-homology ``m`` with ``RC(C[inv_off-m:inv_off]) == C[upt_off:upt_off+m]``.

    ``inv_off`` is the consensus offset of the *aligned* (known-head) arm of
    the inverted element; its withheld head is ``C[inv_off-m:inv_off]``.
    ``upt_off`` is the consensus head offset of the upright element.  Returns
    the longest exact qualifying ``m`` in
    ``[MIN_OVERLAP_MH, MAX_OVERLAP_MH]`` or ``0``.
    """
    n = len(consensus)
    for m in range(MAX_OVERLAP_MH, MIN_OVERLAP_MH - 1, -1):
        if inv_off - m < 0 or upt_off + m > n:
            continue
        withheld = consensus[inv_off - m : inv_off]
        upright_head = consensus[upt_off : upt_off + m]
        if _reverse_complement(withheld) == upright_head:
            return m
    return 0


def _evaluate_switch(
    inverted_off: int,
    inverted_arm_len: int,
    upright_off: int,
    score: float,
    consensus: str,
) -> _SwitchCandidate | None:
    """Score one inverted/upright arm pair sharing an ``m``-bp strand switch."""
    overlap_bp = _longest_overlap(inverted_off, upright_off, consensus)
    if overlap_bp < MIN_OVERLAP_MH:
        return None
    return _SwitchCandidate(
        score=score + overlap_bp,
        inverted_off=inverted_off,
        breakpoint=upright_off,
        overlap_bp=overlap_bp,
        inverted_length=inverted_arm_len + overlap_bp,
    )


def _detect_switch_in_clip(clip: str, consensus: str) -> _SwitchCandidate | None:
    """Search one soft-clip for the twin-priming reverse-complement switch.

    Both arms are aligned in each orientation and every valid inverted/upright
    assignment is scored; the strand of the insertion determines which
    assignment fires.  Score ties are broken toward the shorter inverted copy
    (the 5' inverted head is characteristically truncated), and exact
    parameter ties retain the first (lowest-split) candidate.
    """
    length = len(clip)
    best: _SwitchCandidate | None = None
    for split in range(MIN_SEG_LEN, length - MIN_SEG_LEN + 1):
        prefix = clip[:split]
        suffix = clip[split:]

        prefix_fwd = _best_alignment(prefix, consensus)
        prefix_rev = _best_alignment(_reverse_complement(prefix), consensus)
        suffix_fwd = _best_alignment(suffix, consensus)
        suffix_rev = _best_alignment(_reverse_complement(suffix), consensus)

        arm_score = len(prefix) + len(suffix)
        candidates: list[_SwitchCandidate | None] = []

        if prefix_rev is not None and suffix_fwd is not None:
            # Inverted copy aligned on the reverse-complement strand.
            candidates.append(
                _evaluate_switch(
                    prefix_rev.offset,
                    len(prefix),
                    suffix_fwd.offset,
                    arm_score * (prefix_rev.identity + suffix_fwd.identity),
                    consensus,
                )
            )
            # Upright body on the reverse-complement strand, inverted copy forward.
            candidates.append(
                _evaluate_switch(
                    suffix_fwd.offset,
                    len(suffix),
                    prefix_rev.offset,
                    arm_score * (prefix_rev.identity + suffix_fwd.identity),
                    consensus,
                )
            )
        if prefix_fwd is not None and suffix_rev is not None:
            candidates.append(
                _evaluate_switch(
                    prefix_fwd.offset,
                    len(prefix),
                    suffix_rev.offset,
                    arm_score * (prefix_fwd.identity + suffix_rev.identity),
                    consensus,
                )
            )
            candidates.append(
                _evaluate_switch(
                    suffix_rev.offset,
                    len(suffix),
                    prefix_fwd.offset,
                    arm_score * (prefix_fwd.identity + suffix_rev.identity),
                    consensus,
                )
            )

        for candidate in candidates:
            if candidate is None:
                continue
            better = best is None or candidate.score > best.score or (
                candidate.score == best.score and candidate.inverted_length < best.inverted_length
            )
            if better:
                best = candidate
    return best


def detect_twin_primed_inversion(soft_clips: list[str], l1_consensus_fasta: str) -> TwinPrimingCall:
    """Detect twin-primed 5' inversion breakpoints in *soft_clips*.

    Args:
        soft_clips: Soft-clip sequences at the MEI 5' junction (untrusted
            strand).  Sequences are uppercased and deduplicated before
            analysis.
        l1_consensus_fasta: Path to a FASTA file whose first record is the
            L1HS reference consensus.  The consensus is read 5' -> 3' and is
            used as the strand-switch reference.

    Returns:
        A :class:`TwinPrimingCall`.  ``is_twin_primed`` is True when any clip
        contained an inverted (reverse-complement) prefix joined through a
        micro-homology junction to an upright (forward) suffix.  When no
        consensus is available, no clip is usable, or no strand switch is
        found the call carries ``is_twin_primed=False`` with
        ``inversion_breakpoint=-1``, ``inverted_length=0`` and
        ``complementary_overlap_bp=0``.
    """
    consensus = _load_consensus(l1_consensus_fasta)
    if not consensus:
        return TwinPrimingCall(is_twin_primed=False)

    seen: set[str] = set()
    clips: list[str] = []
    for clip in soft_clips:
        if not isinstance(clip, str) or not clip:
            continue
        upper = clip.upper()
        if upper not in seen:
            seen.add(upper)
            clips.append(upper)
    if not clips:
        return TwinPrimingCall(is_twin_primed=False)

    best: _SwitchCandidate | None = None
    for clip in clips:
        if len(clip) < MIN_CLIP_LEN:
            continue
        candidate = _detect_switch_in_clip(clip, consensus)
        if candidate is not None and (best is None or candidate.score > best.score):
            best = candidate

    if best is None:
        return TwinPrimingCall(is_twin_primed=False)

    return TwinPrimingCall(
        is_twin_primed=True,
        inversion_breakpoint=best.breakpoint,
        inverted_length=best.inverted_length,
        complementary_overlap_bp=best.overlap_bp,
    )