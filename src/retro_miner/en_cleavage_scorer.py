"""L1 ORF2p Endonuclease (EN) Cleavage Site Preference Scorer.

LINE-1 ORF2p endonuclease makes the first-strand nick at a sequence-specific
consensus motif ``5'-TTTT/AA-3'`` (degenerates ``TTAA``, ``TTCA``, ``TTAG``).
Authentic TPRT retrotransposition insertions occur almost exclusively at these
motifs, so scoring the local reference context around a candidate breakpoint
provides a biological metric that penalizes chimera/assembly-artifact calls.

This module is deliberately *normalized* (PWM match score in ``[0.0, 1.0]``)
and window-centered (positions ``-4..+2`` around the nick), in contrast to the
log-odds scorer in :mod:`retro_miner.cleavage_motif`.

Literature Anchors: Flasch et al. (2019) NAR 47(21):11186-11201 / Monot et al. (2015).
"""

from __future__ import annotations

from dataclasses import dataclass

CANONICAL_EN_SITE_THRESHOLD: float = 0.70
DEGENERATE_EN_SITE_THRESHOLD: float = 0.40

MOTIF_LEN: int = 6

EN_CONSENSUS: str = "TTTTAA"

_EN_PWM: tuple[dict[str, float], ...] = (
    {"A": 0.06, "C": 0.04, "G": 0.04, "T": 0.86},
    {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    {"A": 0.78, "C": 0.10, "G": 0.08, "T": 0.04},
    {"A": 0.72, "C": 0.12, "G": 0.08, "T": 0.08},
)

CLASS_CANONICAL_EN = "CANONICAL_EN_SITE"
CLASS_DEGENERATE_EN = "DEGENERATE_EN_SITE"
CLASS_EN_INDEPENDENT = "EN_INDEPENDENT"


@dataclass(frozen=True)
class ElemCleavageScore:
    """PWM match of a breakpoint context to the L1 EN consensus motif.

    Attributes:
        is_canonical_en_site: True when ``raw_score >= 0.70``.
        raw_score: Normalized PWM match ``S_EN`` in ``[0.0, 1.0]``.
        motif: Six-base window spanning positions ``-4..+2`` around the nick.
        cleavage_class: Canonical, degenerate, or EN-independent annotation.
    """

    is_canonical_en_site: bool
    raw_score: float
    motif: str
    cleavage_class: str


def _extract_motif_window(ref_sequence_context: str) -> str:
    """Slice the informative ``5'-N_{-4}..N_{+2}-3'`` block from a nick-centered window.

    Assumes *ref_sequence_context* is a genomic window centered on the
    first-strand nick (e.g. the canonical 10-bp ``N_-5..N_+5`` extraction); the
    six informative columns sit immediately around the nick. Windows shorter
    than six bases are right-padded with ``N``.
    """
    seq = (ref_sequence_context or "").upper()
    half = len(seq) // 2
    start = min(max(half - 4, 0), max(len(seq) - MOTIF_LEN, 0))
    motif = seq[start : start + MOTIF_LEN]
    return motif.ljust(MOTIF_LEN, "N")


def _run_pwm_normalized(motif: str) -> float:
    """Score *motif* against the consensus PWM, normalized to ``[0.0, 1.0]``.

    Each informative position contributes the ratio of the observed base
    probability to the consensus base probability in that PWM column, averaged
    over the six columns. Unknown bases (``N``) contribute zero.
    """
    total = 0.0
    for index, base in enumerate(motif):
        column = _EN_PWM[index]
        if base not in column:
            continue
        consensus = column[EN_CONSENSUS[index]]
        total += column[base] / consensus
    return round(total / MOTIF_LEN, 4)


def score_en_cleavage_site(ref_sequence_context: str) -> ElemCleavageScore:
    """Score a reference breakpoint context against the L1 EN consensus motif.

    Args:
        ref_sequence_context: 10-bp genomic window centered on the breakpoint
            (``5'-N_-5..N_+5-3'``), i.e. the nick falls at the window midpoint.

    Returns:
        :class:`ElemCleavageScore` with normalized PWM score, extracted motif,
        and the binary canonical-EN flag (``S_EN >= 0.70``).
    """
    motif = _extract_motif_window(ref_sequence_context)
    raw_score = _run_pwm_normalized(motif)

    if raw_score >= CANONICAL_EN_SITE_THRESHOLD:
        cleavage_class = CLASS_CANONICAL_EN
        is_canonical = True
    elif raw_score >= DEGENERATE_EN_SITE_THRESHOLD:
        cleavage_class = CLASS_DEGENERATE_EN
        is_canonical = False
    else:
        cleavage_class = CLASS_EN_INDEPENDENT
        is_canonical = False

    return ElemCleavageScore(
        is_canonical_en_site=is_canonical,
        raw_score=raw_score,
        motif=motif,
        cleavage_class=cleavage_class,
    )