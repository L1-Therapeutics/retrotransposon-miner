"""L1 Endonuclease (ORF2p) Cleavage Site Preference Scorer.

Evaluates genomic target sequence context against the canonical 5'-TTTT/AA-3'
endonuclease cleavage motif.

Literature Anchor: Flasch et al. (2019) NAR / Monot et al. (2015).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ElemCleavageScore:
    motif: str
    match_score: float
    is_canonical_en_site: bool


# Position Weight Matrix for canonical 5'-TTTT/AA-3' integration preference (-4 to +2)
CANONICAL_MOTIFS: list[str] = [
    "TTTTAA",
    "TTAAAA",
    "TTTCAA",
    "TTTGAA",
    "TTCATT",
    "TTAAGA",
]


def score_en_cleavage_site(ref_sequence_context: str) -> ElemCleavageScore:
    """Score genomic sequence surrounding breakpoint for L1 EN cleavage compatibility.

    Args:
        ref_sequence_context: 10-14 bp genomic sequence surrounding breakpoint.

    Returns:
        ElemCleavageScore with raw motif score and binary canonical flag.
    """
    if not ref_sequence_context:
        return ElemCleavageScore(motif="", match_score=0.0, is_canonical_en_site=False)

    seq = ref_sequence_context.upper()
    best_score = 0.0
    best_motif = seq[:6] if len(seq) >= 6 else seq

    for canon in CANONICAL_MOTIFS:
        for i in range(len(seq) - len(canon) + 1):
            sub = seq[i : i + len(canon)]
            matches = sum(1 for a, b in zip(sub, canon, strict=False) if a == b)
            score = round(matches / len(canon), 4)
            if score > best_score:
                best_score = score
                best_motif = sub

    is_canonical = best_score >= 0.6667
    return ElemCleavageScore(
        motif=best_motif,
        match_score=best_score,
        is_canonical_en_site=is_canonical,
    )
