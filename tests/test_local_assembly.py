import pytest
from retro_miner.local_assembly import assemble_locus_clips, AssembledContig

def test_de_bruijn_assembly_overlapping_clips():
    # Two overlapping soft clips sharing 10-bp motif
    clips = [
        "ATTGCAGCTAGCTAG",
        "AGCTAGCTAGAAAAA",
        "ATTGCAGCTAGCTAG",
        "AGCTAGCTAGAAAAA",
    ]

    contigs = assemble_locus_clips(clips, k=8, min_coverage=2)

    assert len(contigs) > 0
    top = contigs[0]
    assert isinstance(top, AssembledContig)
    assert "ATTGCAGCTAGCTAG" in top.sequence
    assert top.length >= 15
    assert top.support_score >= 2.0

def test_noise_kmers_pruned_below_coverage_threshold():
    # True clips repeated twice, noise clip present once
    clips = [
        "CGATCGATCGATCG",
        "CGATCGATCGATCG",
        "AAAAAAAAAAAAAA", # Noise clip (coverage=1)
    ]

    contigs = assemble_locus_clips(clips, k=6, min_coverage=2)

    assert len(contigs) == 1
    assert "CGATCG" in contigs[0].sequence
    assert "AAAAAA" not in contigs[0].sequence

def test_empty_sequence_list():
    assert assemble_locus_clips([], k=15) == []
