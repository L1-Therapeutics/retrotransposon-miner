"""Diagnostic Marker Subfamily Voting Engine for Mobile Element Insertions.

Implements k-mer log-likelihood classification to discriminate active,
polymorphic retrotransposons (L1HS, AluYa5, AluYb8, SVA-F) from fixed ancestral noise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

DIAGNOSTIC_KMERS: dict[str, list[str]] = {
    "L1HS": ["ACAGAG", "GACAGAG", "TGACAGAG", "AAGACAGAG"],
    "L1PA2": ["ACACAG", "GACACAG", "TGACACAG"],
    "AluYa5": ["CCATCT", "TCCATCT", "GCATCT"],
    "AluYb8": ["TGAGAC", "CTGAGAC", "ACTGAGAC"],
    "AluSx": ["CCACAC", "TCCACAC"],
    "SVA_F": ["CCTCTCT", "GCCTCTCT"],
}


@dataclass(frozen=True)
class SubfamilyCall:
    family: str
    top_subfamily: str
    log_likelihood_ratio: float
    posterior_prob: float
    matching_kmers: int


def classify_mei_subfamily(
    read_sequences: list[str],
    family_hint: str = "L1",
) -> SubfamilyCall:
    """Classify MEI subfamily based on diagnostic k-mer log-likelihood ratios.

    Args:
        read_sequences: List of soft-clipped or unmapped mate query sequences.
        family_hint: Target family hint ("L1", "Alu", or "SVA").

    Returns:
        SubfamilyCall containing top subfamily, LLR, and posterior probability.
    """
    if not read_sequences:
        return SubfamilyCall(family_hint, "UNKNOWN", 0.0, 0.3333, 0)

    combined_text = "$".join(s.upper() for s in read_sequences if s)
    scores: dict[str, int] = {}

    for subfam, kmers in DIAGNOSTIC_KMERS.items():
        if family_hint.upper() in subfam.upper() or subfam.startswith(family_hint.upper()):
            count = sum(combined_text.count(kmer) for kmer in kmers)
            scores[subfam] = count

    if not scores or max(scores.values()) == 0:
        default_sub = "L1HS" if family_hint.upper() == "L1" else "AluYa5" if family_hint.upper() == "ALU" else "SVA_F"
        return SubfamilyCall(family_hint, default_sub, 0.0, 0.50, 0)

    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    top_sub, top_count = sorted_scores[0]
    runner_up_count = sorted_scores[1][1] if len(sorted_scores) > 1 else 0

    # Log-likelihood ratio
    llr = round((top_count - runner_up_count) * math.log(2.0) + 0.1, 4)
    total_matches = sum(scores.values())
    prob = round(top_count / total_matches, 4) if total_matches > 0 else 0.50

    return SubfamilyCall(
        family=family_hint,
        top_subfamily=top_sub,
        log_likelihood_ratio=llr,
        posterior_prob=prob,
        matching_kmers=top_count,
    )
