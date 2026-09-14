"""TPRT Endonuclease Cleavage Motif PWM Scorer for Mobile Element Insertions.

Quantifies the strength of the L1 EN cleavage signature (5'-TTTT/AA-3')
at insertion junctions using a position weight matrix (PWM) to support
Target-Primed Reverse Transcription (TPRT) mechanistic classification.

Literature Anchor: Jurka (1997) / Cost et al. (2002) Cell / Morrish et al. (2002) Nature.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class CleavageMotifResult:
    is_canonical_tprt: bool
    cleavage_class: str
    pwm_score: float
    breakpoint_motif: str


CLASS_CANONICAL_TPRT = "CANONICAL_TPRT"
CLASS_DEGENERATE_TPRT = "DEGENERATE_TPRT"
CLASS_EN_INDEPENDENT = "EN_INDEPENDENT"

CANONICAL_TPRT_THRESHOLD: float = 4.0
DEGENERATE_TPRT_THRESHOLD: float = 2.0

_PWM: dict[str, dict[str, float]] = {
    "1": {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    "2": {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    "3": {"A": 0.10, "C": 0.10, "G": 0.10, "T": 0.70},
    "4": {"A": 0.03, "C": 0.03, "G": 0.03, "T": 0.84},
    "5": {"A": 0.78, "C": 0.10, "G": 0.10, "T": 0.10},
    "6": {"A": 0.78, "C": 0.10, "G": 0.10, "T": 0.10},
    "7": {"A": 0.34, "C": 0.10, "G": 0.10, "T": 0.10},
    "8": {"A": 0.32, "C": 0.08, "G": 0.03, "T": 0.03},
}
_BG: float = 0.25


def _pwm_score(motif: str) -> float:
    score = 0.0
    for i, ch in enumerate(motif.upper()):
        pos = str(i + 1)
        p = _PWM.get(pos, {}).get(ch, 0.25)
        score += math.log2(max(p, 1e-6) / _BG)
    return score


def score_cleavage_motif(
    seq: str,
    *,
    breakpoint_offset: int = 4,
    window_size: int = 8,
) -> CleavageMotifResult:
    """Score TPRT endonuclease cleavage motif strength in a junction sequence.

    Extracts an 8-bp window centered on the breakpoint and scores it against
    a position weight matrix favoring the L1 EN consensus ``5'-TTTT/AA-3'``.

    Args:
        seq: Junction or soft-clip query sequence.
        breakpoint_offset: 0-based index of the nick site within *seq*.
        window_size: Width of the scored motif window.

    Returns:
        CleavageMotifResult with PWM score, class, and extracted motif.
    """
    s = (seq or "").upper().strip()
    if not s:
        return CleavageMotifResult(
            is_canonical_tprt=False,
            cleavage_class=CLASS_EN_INDEPENDENT,
            pwm_score=0.0,
            breakpoint_motif="N" * window_size,
        )

    start = max(0, breakpoint_offset - 4)
    end = min(len(s), start + window_size)
    motif = s[start:end]
    if len(motif) < window_size:
        motif = motif.ljust(window_size, "N")

    score = _pwm_score(motif)

    if score >= CANONICAL_TPRT_THRESHOLD:
        cls = CLASS_CANONICAL_TPRT
        is_canonical = True
    elif score >= DEGENERATE_TPRT_THRESHOLD:
        cls = CLASS_DEGENERATE_TPRT
        is_canonical = False
    else:
        cls = CLASS_EN_INDEPENDENT
        is_canonical = False

    return CleavageMotifResult(
        is_canonical_tprt=is_canonical,
        cleavage_class=cls,
        pwm_score=round(score, 4),
        breakpoint_motif=motif,
    )
