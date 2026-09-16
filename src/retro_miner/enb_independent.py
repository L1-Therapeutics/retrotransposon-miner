"""EN-Independent (DSB Repair) MEI Integration Site Classifier.

In DNA-repair-deficient or radiomimetic backgrounds, L1 retrotransposition can
proceed without ORF2p endonuclease (EN) cleavage.  These EN-independent
insertions are mechanistically distinct from canonical target-primed reverse
transcription (TPRT):

* no Target Site Duplication (``tsd_length == 0``),
* no canonical ``5'-TTTT/AA-3'`` cleavage motif (``en_score`` below threshold),
* a pre-existing double-strand break leaves a target-site genomic deletion
  (typically 10-1000 bp) at the insertion locus.

This classifier flags such events so the VCF can annotate them
(``EN_INDEPENDENT=1`` / ``TARGET_DEL_BP=N``).

Literature Anchor: Morrish et al. (2002), Nature Genetics 31(2):159-165.
"""

from __future__ import annotations

from dataclasses import dataclass

EN_MOTIF_SCORE_THRESHOLD: float = 0.30
MIN_TARGET_DELETION_BP: int = 1
DELETION_SCALE_BP: int = 1000


@dataclass(frozen=True)
class ENIndependentCall:
    """Classification of an insertion as an EN-independent (DSB-repair) event.

    Attributes:
        is_en_independent: True when no TSD, no canonical EN motif, and a
            target-site genomic deletion are all observed.
        deletion_size: Target-site genomic deletion size in bp (>= 0).
        cleavage_motif_score: PWM match score of the breakpoint context to the
            L1 EN consensus ``5'-TTTT/AA-3'`` in ``[0.0, 1.0]``.
        confidence: Calibrated confidence in the EN-independent call in
            ``[0.0, 1.0]`` (0.0 for non-EN-independent events).
    """

    is_en_independent: bool
    deletion_size: int
    cleavage_motif_score: float
    confidence: float


def classify_en_independent_event(
    tsd_length: int,
    en_score: float,
    span_deletion_bp: int,
) -> ENIndependentCall:
    """Classify an insertion as EN-independent integration at a pre-existing DSB.

    An event is flagged EN-independent when all three criteria hold:

    * ``tsd_length == 0`` (no target site duplication),
    * ``en_score < 0.30`` (no canonical ``5'-TTTT/AA-3'`` EN cleavage motif),
    * ``span_deletion_bp > 0`` (target-site genomic deletion present).

    Args:
        tsd_length: Resolved TSD length in bp (0 = no duplication).
        en_score: Normalized EN cleavage motif PWM score in ``[0.0, 1.0]``.
        span_deletion_bp: Size of the target-site genomic deletion in bp.

    Returns:
        :class:`ENIndependentCall` with the binary flag, deletion size, motif
        score, and calibrated confidence.
    """
    tsd_len = int(tsd_length)
    del_bp = max(int(span_deletion_bp), 0)
    motif_score = max(0.0, min(float(en_score), 1.0))

    is_en_independent = (
        tsd_len == 0
        and motif_score < EN_MOTIF_SCORE_THRESHOLD
        and del_bp >= MIN_TARGET_DELETION_BP
    )

    if is_en_independent:
        deletion_band = min(del_bp / DELETION_SCALE_BP, 1.0)
        motif_deficit = 1.0 - min(motif_score / EN_MOTIF_SCORE_THRESHOLD, 1.0)
        confidence = 0.5 + 0.25 * deletion_band + 0.25 * motif_deficit
    else:
        confidence = 0.0

    return ENIndependentCall(
        is_en_independent=is_en_independent,
        deletion_size=del_bp,
        cleavage_motif_score=round(motif_score, 4),
        confidence=round(min(max(confidence, 0.0), 1.0), 4),
    )