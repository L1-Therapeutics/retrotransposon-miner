"""De Bruijn Micro-Assembly Engine for Breakpoint Junction Reconstruction.

Assembles soft-clipped read fragments at candidate loci into linear unitigs
to resolve full-length MEI insertion junctions and target site duplications.

Literature Anchor: Cameron et al. (2017) GRIDSS / DeBruijn graph micro-assembly.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass


@dataclass(frozen=True)
class AssembledContig:
    sequence: str
    length: int
    kmer_count: int
    support_score: float


def build_kmer_graph(
    sequences: list[str],
    k: int = 15,
    min_coverage: int = 2,
) -> tuple[dict[str, list[str]], dict[str, int], dict[str, int]]:
    """Build a directed de Bruijn graph from soft-clipped sequences.

    Args:
        sequences: Soft-clipped query sequences.
        k: Length of k-mers (default 15).
        min_coverage: Minimum observation count to keep a k-mer (default 2).

    Returns:
        Tuple of (adjacency_list, in_degrees, kmer_counts).
    """
    kmer_counts: dict[str, int] = defaultdict(int)

    for seq in sequences:
        s = seq.upper()
        if len(s) < k:
            continue
        for i in range(len(s) - k + 1):
            kmer = s[i : i + k]
            kmer_counts[kmer] += 1

    adj: dict[str, list[str]] = defaultdict(list)
    in_degrees: dict[str, int] = defaultdict(int)

    for kmer, count in kmer_counts.items():
        if count < min_coverage:
            continue
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adj[prefix].append(suffix)
        in_degrees[suffix] += 1
        if prefix not in in_degrees:
            in_degrees[prefix] = 0

    return dict(adj), dict(in_degrees), dict(kmer_counts)


def assemble_locus_clips(
    sequences: list[str],
    k: int = 15,
    min_coverage: int = 2,
) -> list[AssembledContig]:
    """Reconstruct contiguous unitigs from soft-clipped sequences.

    Args:
        sequences: List of soft-clipped read sequences at insertion locus.
        k: k-mer size for graph construction.
        min_coverage: Minimum k-mer frequency threshold.

    Returns:
        List of AssembledContig objects sorted by sequence length and support.
    """
    if not sequences:
        return []

    adj, in_degrees, kmer_counts = build_kmer_graph(sequences, k=k, min_coverage=min_coverage)
    if not adj:
        return []

    # Identify source nodes (in-degree == 0 or branching)
    sources = [node for node in adj if in_degrees.get(node, 0) == 0]
    if not sources:
        sources = list(adj.keys())[:1]

    visited_edges: set[tuple[str, str]] = set()
    contigs: list[AssembledContig] = []

    for start_node in sources:
        curr = start_node
        path = [curr]

        while curr in adj and len(adj[curr]) > 0:
            next_node = adj[curr][0]
            edge = (curr, next_node)
            if edge in visited_edges:
                break
            visited_edges.add(edge)
            path.append(next_node[-1])
            curr = next_node

        contig_seq = "".join(path)
        if len(contig_seq) >= k:
            total_kmers = len(contig_seq) - k + 1
            support = sum(
                kmer_counts.get(contig_seq[i : i + k], 1)
                for i in range(total_kmers)
            ) / max(1, total_kmers)

            contigs.append(
                AssembledContig(
                    sequence=contig_seq,
                    length=len(contig_seq),
                    kmer_count=total_kmers,
                    support_score=round(support, 2),
                )
            )

    contigs.sort(key=lambda c: (c.length, c.support_score), reverse=True)
    return contigs
