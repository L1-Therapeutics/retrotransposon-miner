"""Scan the L1.3 consensus for strand-asymmetric RNA-processing signals.

Question this answers
---------------------
Belancio 2006 (NAR 34:1512) reports that L1 carries splice donor/acceptor sites
on BOTH strands, but that an ESE predictor found "four times as many ESEs in the
sense strand of L1.3 as in the antisense". Two mechanisms are therefore on the
table for why *sense*-oriented intronic L1 insertions are depleted (Medstrand
2002): premature polyadenylation, and aberrant splicing. They make different
strand predictions. This module measures the asymmetry directly.

Coordinates are 1-based on L1.3 (GenBank L19088, 6059 bp), the element Belancio
numbered their experimentally confirmed splice sites against.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


# --- poly(A) signals -------------------------------------------------------
# Canonical hexamer plus the single-nucleotide variants catalogued by
# Beaudoing et al. 2000 (Genome Res 10:1001) as functional in human 3' ends.
# AATAAA and ATTAAA together account for the large majority of real sites.
POLYA_HEXAMERS = [
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
    "AATATA", "AATACA", "AATAGA", "ACTAAA", "AAGAAA", "AATGAA",
]
POLYA_STRONG = {"AATAAA", "ATTAAA"}


@dataclass(frozen=True)
class Hit:
    strand: str        # '+' = L1 sense, '-' = L1 antisense
    kind: str          # 'polyA' | 'donor' | 'acceptor'
    l1_pos: int        # 1-based position on the L1.3 SENSE coordinate frame
    motif: str
    score: float


def scan_polya(seq: str) -> list[Hit]:
    """Find poly(A) hexamers on both strands, reported in L1 sense coordinates.

    A '-' strand hit means the hexamer reads 5'->3' on the antisense strand, i.e.
    it would terminate a transcript running antisense through the element.

    ``l1_pos`` is always the LEFTMOST (lowest) sense-frame coordinate of the
    motif footprint, on either strand, so window-containment tests are uniform.
    """
    n = len(seq)
    rc = revcomp(seq)
    hits: list[Hit] = []
    for hexamer in POLYA_HEXAMERS:
        weight = 1.0 if hexamer in POLYA_STRONG else 0.5
        for m in re.finditer(f"(?={hexamer})", seq):
            hits.append(Hit("+", "polyA", m.start() + 1, hexamer, weight))
        for m in re.finditer(f"(?={hexamer})", rc):
            # rc index i covers sense 1-based [n-i-L+1, n-i]; take the left edge.
            sense_pos = n - m.start() - len(hexamer) + 1
            hits.append(Hit("-", "polyA", sense_pos, hexamer, weight))
    return hits


# --- null models -----------------------------------------------------------
# A strand asymmetry in an A-rich motif is EXPECTED from base composition alone:
# L1's sense strand is A-rich, so its antisense strand is T-rich and cannot carry
# as many AATAAA-like hexamers. To claim biological enrichment we must exceed
# what composition predicts. Two nulls, increasing in strictness:
#   mono   - preserves exact base composition, destroys all order
#   markov - preserves base AND dinucleotide structure (first-order chain)


def shuffle_mono(seq: str, rng) -> str:
    a = list(seq)
    rng.shuffle(a)
    return "".join(a)


def fit_markov1(seq: str):
    """Return (states, start_probs, transition_matrix) for a first-order chain."""
    import numpy as np
    states = "ACGT"
    idx = {c: i for i, c in enumerate(states)}
    T = np.zeros((4, 4))
    for a, b in zip(seq, seq[1:]):
        T[idx[a], idx[b]] += 1
    T += 1e-9                       # avoid zero rows
    T /= T.sum(axis=1, keepdims=True)
    counts = np.array([seq.count(c) for c in states], dtype=float)
    return states, counts / counts.sum(), T


def simulate_markov1(n: int, states, p0, T, rng) -> str:
    import numpy as np
    k = rng.choice(4, p=p0)
    out = [states[k]]
    for _ in range(n - 1):
        k = rng.choice(4, p=T[k])
        out.append(states[k])
    return "".join(out)


def strand_counts(seq: str, canonical_only: bool = False) -> tuple[int, int]:
    """(sense, antisense) poly(A) hexamer counts for one sequence."""
    hits = scan_polya(seq)
    if canonical_only:
        hits = [h for h in hits if h.motif in POLYA_STRONG]
    s = sum(1 for h in hits if h.strand == "+")
    return s, len(hits) - s


def _self_test() -> None:
    """Verify the antisense coordinate mapping against a hand-checked case."""
    # AATAAA on the antisense strand == TTTATT on the sense strand.
    # Place it so the sense footprint is exactly positions 11..16 (1-based).
    seq = "G" * 10 + "TTTATT" + "C" * 10
    assert revcomp("TTTATT") == "AATAAA"
    anti = [h for h in scan_polya(seq) if h.strand == "-" and h.motif == "AATAAA"]
    assert len(anti) == 1, f"expected 1 antisense AATAAA, got {len(anti)}"
    assert anti[0].l1_pos == 11, f"expected left edge 11, got {anti[0].l1_pos}"
    # and the sense-strand case, same logic
    seq2 = "G" * 10 + "AATAAA" + "C" * 10
    sense = [h for h in scan_polya(seq2) if h.strand == "+" and h.motif == "AATAAA"]
    assert len(sense) == 1 and sense[0].l1_pos == 11, sense
    # Overlapping motifs must ALL be found (lookahead, non-consuming).
    # "AATAAATAAA" contains AATAAA at 1-based 1 and again at 5.
    ov = [h.l1_pos for h in scan_polya("AATAAATAAA")
          if h.strand == "+" and h.motif == "AATAAA"]
    assert sorted(ov) == [1, 5], f"overlap handling wrong: {ov}"


if __name__ == "__main__":
    _self_test()
    print("self-test passed")


def in_window(hits: list[Hit], start: int, end: int) -> list[Hit]:
    """Hits whose position falls inside a retained sub-span of the element."""
    return [h for h in hits if start <= h.l1_pos <= end]
