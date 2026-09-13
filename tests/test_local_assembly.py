import pytest

from retro_miner.local_assembly import AssembledContig, assemble_locus_clips, build_kmer_graph


class TestBuildKmerGraph:
    def test_non_empty_sequences_produce_adjacency(self):
        adj, in_degrees, kmer_counts = build_kmer_graph(["ATTGCAGCTAGCTAG", "AGCTAGCTAGAAAAA"], k=4)
        assert isinstance(adj, dict)
        assert len(adj) > 0

    def test_empty_sequences_return_empty_graph(self):
        adj, in_degrees, kmer_counts = build_kmer_graph([], k=21)
        assert adj == {}
        assert in_degrees == {}
        assert kmer_counts == {}

    def test_singleton_noise_pruned(self):
        adj, in_degrees, kmer_counts = build_kmer_graph(["ACGTACGTACGT", "TTTTTTTT"], k=4)
        for neighbors in adj.values():
            for nxt in neighbors:
                assert nxt in adj


class TestAssembleLocusClips:
    def test_overlapping_soft_clips_assemble_contig(self):
        sequences = ["ATTGCAGCTAGCTAG", "AGCTAGCTAGAAAAA"]
        contigs = assemble_locus_clips(sequences, k=4, min_coverage=2)
        assert len(contigs) >= 1
        longest = contigs[0]
        assert isinstance(longest, AssembledContig)
        assert longest.length >= 1
        assert "AGCTAG" in longest.sequence

    def test_low_coverage_noise_filtered(self):
        clips = [
            "CGATCGATCGATCG",
            "CGATCGATCGATCG",
            "AAAAAAAAAAAAAA",
        ]
        contigs = assemble_locus_clips(clips, k=6, min_coverage=2)
        assert len(contigs) == 1
        assert "CGATCG" in contigs[0].sequence
        assert "AAAAAA" not in contigs[0].sequence

    def test_empty_input_returns_empty(self):
        contigs = assemble_locus_clips([], k=15)
        assert contigs == []

    def test_short_sequences_below_k_skipped(self):
        contigs = assemble_locus_clips(["ACG", "TTT"], k=5, min_coverage=2)
        assert contigs == []

    def test_contigs_ordered_by_length_desc(self):
        sequences = ["AAAAAA", "CCCC", "GGGGGGGG"]
        contigs = assemble_locus_clips(sequences, k=3, min_coverage=2)
        for i in range(len(contigs) - 1):
            assert contigs[i].length >= contigs[i + 1].length
