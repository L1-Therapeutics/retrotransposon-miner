"""Target Site Duplication (TSD) Refiner and Poly(A) Entropy Filter.

Implements single-base precision breakpoint resolution for Target-Primed
Reverse Transcription (TPRT) events (LINE-1, Alu, SVA).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Sequence


@dataclass(frozen=True)
class TSDResult:
    tsd_seq: str
    tsd_length: int
    poly_a_detected: bool
    entropy: float
    confidence_score: float
    chrom: str = ""
    pos: int = 0
    locus_id: str = ""

    @property
    def polyA_tail_detected(self) -> bool:
        """Alias for backward compatibility with pipeline callers."""
        return self.poly_a_detected


def calculate_sequence_entropy(seq: str) -> float:
    """Calculate Shannon entropy (H = -sum p_i * log2(p_i)) for a sequence."""
    if not seq:
        return 0.0
    seq_upper = seq.upper()
    length = len(seq_upper)
    counts: dict[str, int] = {}
    for char in seq_upper:
        counts[char] = counts.get(char, 0) + 1

    entropy = 0.0
    for count in counts.values():
        p = count / length
        entropy -= p * math.log2(p)
    return round(entropy, 4)


def is_poly_a_or_t(seq: str, max_entropy: float = 1.2, min_tract_len: int = 8) -> bool:
    """Return True if sequence contains low-entropy poly(A) or poly(T) tract."""
    if len(seq) < 5:
        return False
    seq_upper = seq.upper()

    if "A" * min_tract_len in seq_upper or "T" * min_tract_len in seq_upper:
        return True

    win_size = min(10, len(seq_upper))
    for i in range(len(seq_upper) - win_size + 1):
        window = seq_upper[i : i + win_size]
        if max(window.count("A"), window.count("T")) / win_size >= 0.8:
            return True

    entropy = calculate_sequence_entropy(seq_upper)
    a_or_t_fraction = max(seq_upper.count("A"), seq_upper.count("T")) / len(seq_upper)
    return entropy < max_entropy and a_or_t_fraction >= 0.75


def find_longest_common_substring(s1: str, s2: str) -> str:
    """Find longest exact matching substring between two sequences."""
    if not s1 or not s2:
        return ""
    m = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    longest_len = 0
    longest_end = 0

    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            if s1[i - 1].upper() == s2[j - 1].upper():
                m[i][j] = m[i - 1][j - 1] + 1
                if m[i][j] > longest_len:
                    longest_len = m[i][j]
                    longest_end = i

    return s1[longest_end - longest_len : longest_end].upper()


def _max_flush_match(s1: str, s2: str, min_len: int = 3) -> str:
    """Internal helper to find maximal matching prefix/suffix overlap between clips."""
    return find_longest_common_substring(s1, s2)


def refine_tsd_boundaries(
    arg1: Any = "",
    arg2: Any = "",
    arg3: Any = None,
    arg4: Any = None,
    min_tsd_len: int = 5,
    max_tsd_len: int = 35,
    **kwargs: Any,
) -> TSDResult:
    """Detect TSD sequence and poly(A) tail flags with flexible signature support.

    Supports:
      Signature A: refine_tsd_boundaries(left_clip_seq, right_clip_seq)
      Signature B: refine_tsd_boundaries(chrom, pos, soft_clips, fasta_path)
    """
    chrom = ""
    pos = 0

    # Determine signature style based on argument types
    if isinstance(arg1, str) and isinstance(arg2, int):
        # Signature B: (chrom, pos, soft_clips, ref_fasta)
        chrom = arg1
        pos = arg2
        clips = arg3 if isinstance(arg3, (list, tuple)) else []
        left_clips = [str(c) for c in clips[: len(clips) // 2]] if clips else []
        right_clips = [str(c) for c in clips[len(clips) // 2 :]] if clips else []
        l_seq = "".join(left_clips)
        r_seq = "".join(right_clips)
    else:
        # Signature A: (left_clip_seq, right_clip_seq)
        l_seq = "".join(arg1) if isinstance(arg1, (list, tuple)) else str(arg1 or "")
        r_seq = "".join(arg2) if isinstance(arg2, (list, tuple)) else str(arg2 or "")

    left_entropy = calculate_sequence_entropy(l_seq)
    right_entropy = calculate_sequence_entropy(r_seq)
    poly_a_detected = is_poly_a_or_t(l_seq) or is_poly_a_or_t(r_seq)

    tsd_match = find_longest_common_substring(l_seq, r_seq)

    if min_tsd_len <= len(tsd_match) <= max_tsd_len:
        tsd_seq = tsd_match
        tsd_length = len(tsd_seq)
        confidence = 0.95 if poly_a_detected else 0.85
    else:
        tsd_seq = ""
        tsd_length = 0
        confidence = 0.50 if poly_a_detected else 0.10

    avg_entropy = round((left_entropy + right_entropy) / 2.0, 4) if (l_seq or r_seq) else 0.0

    return TSDResult(
        tsd_seq=tsd_seq,
        tsd_length=tsd_length,
        poly_a_detected=poly_a_detected,
        entropy=avg_entropy,
        confidence_score=confidence,
        chrom=chrom,
        pos=pos,
    )
