"""Tests for the twin-priming 5' inversion breakpoint analyzer.

Builds a synthetic L1HS reference consensus and constructs inverted-truncated
junction sequences carrying a 4-bp micro-homology overlap (5'-ACTG-3'),
mirroring the twin-priming structures of Ostertag & Kazazian (2001).
"""

import random
from pathlib import Path

from retro_miner.twin_priming import TwinPrimingCall, detect_twin_primed_inversion

_CONSENSUS_LEN = 300
_BASE = "ACGT"
# Micro-homology at the strand-switch junction (shared 4-bp block).
_MICRO_HOMOLOGY = "ACTG"
# Withheld inverted-copy head: RC(C[14:18]) == "ACTG".
_HEAD_TETRAMER = "CAGT"
_HEAD_START = 14
_INVERTED_OFF = 18
_INVERTED_END = 74
_BREAKPOINT = 70
_UPRIGHT_LEN = 90


def build_consensus() -> str:
    rng = random.Random(7)
    bases = [rng.choice(_BASE) for _ in range(_CONSENSUS_LEN)]
    # Controlled flanks so the only exact micro-homology at the junction is the
    # 4-bp ACTG block (i.e. no RC-partner base extends it to m=5..8, and no
    # mirror split can win the score tie-break).
    bases[0:4] = list(_HEAD_TETRAMER)
    bases[4:6] = ["G", "G"]
    bases[10:14] = ["C", "T", "C", "A"]
    bases[_HEAD_START : _HEAD_START + 4] = list(_HEAD_TETRAMER)
    bases[18] = "A"
    bases[19] = "C"
    bases[68:70] = ["T", "T"]
    bases[69] = "G"
    bases[_BREAKPOINT : _BREAKPOINT + 4] = list(_MICRO_HOMOLOGY)
    bases[74:78] = ["C", "C", "C", "A"]
    return "".join(bases)


def _reversed_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _inverted_junction(consensus: str, *, noisy_bases: tuple[tuple[int, str], ...] = ()) -> str:
    """RC(C[18:74]) joined to C[70:160] across the ACTG micro-homology."""
    inverted_tail = _reversed_complement(consensus[_INVERTED_OFF:_INVERTED_END])
    upright = consensus[_BREAKPOINT : _BREAKPOINT + _UPRIGHT_LEN]
    junction = inverted_tail + upright
    seq = list(junction)
    for pos, base in noisy_bases:
        if 0 <= pos < len(seq) and base != seq[pos]:
            seq[pos] = base
    return "".join(seq)


def _write_consensus(tmp_path: Path, consensus: str) -> str:
    fasta = tmp_path / "l1hs_consensus.fa"
    fasta.write_text(f">L1HS\n{consensus}\n")
    return str(fasta)


class TestDetectInvertedJunction:
    def test_inverted_truncated_junction_is_twin_primed(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        junction = _inverted_junction(consensus)

        call = detect_twin_primed_inversion([junction], fasta)

        assert isinstance(call, TwinPrimingCall)
        assert call.is_twin_primed is True
        assert call.inversion_breakpoint == _BREAKPOINT
        assert call.inverted_length == _INVERTED_END - _HEAD_START
        assert call.complementary_overlap_bp == 4

    def test_micro_homology_is_actg_at_switch_point(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        junction = _inverted_junction(consensus)

        call = detect_twin_primed_inversion([junction], fasta)

        assert consensus[_BREAKPOINT : _BREAKPOINT + call.complementary_overlap_bp] == _MICRO_HOMOLOGY

    def test_accepts_mild_sequencing_noise(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        junction = _inverted_junction(
            consensus,
            noisy_bases=[(10, "T"), (120, "G")],
        )

        call = detect_twin_primed_inversion([junction], fasta)

        assert call.is_twin_primed is True
        assert call.complementary_overlap_bp == 4

    def test_output_is_strand_orientation_agnostic(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)

        plus_junction = _inverted_junction(consensus)
        minus_junction = _reversed_complement(plus_junction)

        plus_call = detect_twin_primed_inversion([plus_junction], fasta)
        minus_call = detect_twin_primed_inversion([minus_junction], fasta)

        assert plus_call.is_twin_primed is True
        assert minus_call.is_twin_primed is True
        assert plus_call.inversion_breakpoint == minus_call.inversion_breakpoint == _BREAKPOINT
        assert plus_call.inverted_length == minus_call.inverted_length == _INVERTED_END - _HEAD_START
        assert plus_call.complementary_overlap_bp == minus_call.complementary_overlap_bp == 4


class TestDetectNonInverted:
    def test_upright_l1_insertion_is_not_twin_primed(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        upright_clip = consensus[20:170]

        call = detect_twin_primed_inversion([upright_clip], fasta)

        assert call.is_twin_primed is False
        assert call.inversion_breakpoint == -1
        assert call.inverted_length == 0
        assert call.complementary_overlap_bp == 0

    def test_pure_inverted_copy_without_upright_partner_is_not_twin_primed(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        inverted_only = _reversed_complement(consensus[0:100])

        call = detect_twin_primed_inversion([inverted_only], fasta)

        assert call.is_twin_primed is False

    def test_non_mei_genomic_clip_is_not_twin_primed(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        rng = random.Random(99)
        genomic = "".join(rng.choice(_BASE) for _ in range(120))

        call = detect_twin_primed_inversion([genomic], fasta)

        assert call.is_twin_primed is False

    def test_short_clip_below_minimum_length_is_not_twin_primed(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        junction = _inverted_junction(consensus)[:15]

        call = detect_twin_primed_inversion([junction], fasta)

        assert call.is_twin_primed is False


class TestEdgeCases:
    def test_empty_clip_list_returns_negative_defaults(self, tmp_path):
        fasta = _write_consensus(tmp_path, build_consensus())

        call = detect_twin_primed_inversion([], fasta)

        assert call.is_twin_primed is False
        assert call.inversion_breakpoint == -1
        assert call.inverted_length == 0
        assert call.complementary_overlap_bp == 0

    def test_missing_consensus_fasta_returns_negative_defaults(self, tmp_path):
        call = detect_twin_primed_inversion(
            ["CTGCGGCCACTGACCGCTCAAACCGGCAGCGGTCGATTAATACGACTCACTATAGGGCCGC"],
            str(tmp_path / "missing.fa"),
        )

        assert call.is_twin_primed is False
        assert call.inversion_breakpoint == -1

    def test_clips_are_uppercased_and_deduplicated(self, tmp_path):
        consensus = build_consensus()
        fasta = _write_consensus(tmp_path, consensus)
        junction = _inverted_junction(consensus)
        lowercase = junction.lower()

        call = detect_twin_primed_inversion([junction, lowercase, junction], fasta)

        assert call.is_twin_primed is True