"""De Bruijn Micro-Assembly Engine for Breakpoint Junction Reconstruction.

Assembles soft-clipped read fragments at candidate loci into linear unitigs
to resolve full-length MEI insertion junctions and target site duplications.

Also provides ``annotate_silver_with_local_assembly`` - the mei-annotate
integration entrypoint.  It adapts the De Bruijn engine (plus clip-anchor
breakpoint estimation, reference-assisted TSD resolution through
:mod:`retro_miner.tsd_refiner`, and poly(A) tail detection) onto the legacy
SPAdes-oriented signature and output schema expected by
:mod:`retro_miner.mei_support`, so silver-stage candidate refinements keep
working without the external assembler.

Literature Anchor: Cameron et al. (2017) GRIDSS / DeBruijn graph micro-assembly.
"""

from __future__ import annotations

import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

_ASM_CLIP_MIN_BP = 5


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


_ASM_OUTPUT_COLUMNS: list[str] = [
    "chrom",
    "window_start",
    "window_end",
    "asm_status",
    "asm_interval_start",
    "asm_interval_end",
    "asm_retry_used",
    "asm_disease_reads_extracted",
    "asm_control_reads_extracted",
    "asm_disease_contig_count",
    "asm_control_contig_count",
    "asm_disease_max_contig_len",
    "asm_control_max_contig_len",
    "asm_error_message",
    "asm_disease_primary_contig_id",
    "asm_control_primary_contig_id",
    "asm_disease_breakpoint_pos",
    "asm_control_breakpoint_pos",
    "asm_consensus_breakpoint_pos",
    "asm_breakpoint_source",
    "asm_tsd_seq",
    "asm_tsd_len",
    "asm_polyA_max_run",
    "asm_insertion_length",
    "asm_insertion_length_observed",
    "asm_insertion_length_imputed",
    "asm_insertion_length_confidence_tier",
    "asm_insertion_mei_start",
    "asm_insertion_mei_end",
    "asm_mei_family",
    "asm_mei_subfamily",
    "asm_mei_target_length",
    "asm_insertion_orientation",
    "asm_consensus_primary_contig_id",
    "asm_complex_class",
    "asm_non_mei_partner_chrom",
    "asm_non_mei_partner_pos",
    "asm_non_mei_partner_type",
    "asm_breakpoint_side_status",
    "asm_complexity_source",
    "asm_top_contigs",
    "asm_mei_alignment_preset",
    "asm_left_support_contig_id",
    "asm_right_support_contig_id",
    "asm_left_support_mei_start",
    "asm_left_support_mei_end",
    "asm_right_support_mei_start",
    "asm_right_support_mei_end",
    "asm_left_support_mei_aln_len",
    "asm_right_support_mei_aln_len",
    "asm_microhomology_sequence",
    "asm_junction_overlap_sequence",
    "asm_coord_model",
    "asm_coord_logic_version",
    "asm_adaptive_read_cap_rerun",
    "asm_adaptive_read_cap_target",
    "asm_disease_recruited_evidence_read_names",
    "asm_control_recruited_evidence_read_names",
]


def _clip_lens(read) -> tuple[int, int]:
    """Return ``(leading_clip_bp, trailing_clip_bp)`` for an aligned read."""
    if read.cigartuples is None:
        return 0, 0
    leading = read.cigartuples[0]
    trailing = read.cigartuples[-1]
    lead = leading[1] if leading[0] == 4 else 0
    trail = trailing[1] if trailing[0] == 4 else 0
    return lead, trail


def _median_or_zero(values: list[int]) -> int:
    if not values:
        return 0
    return int(round(statistics.median(values)))


def _longest_at_run(sequences: list[str]) -> int:
    """Longest run of consecutive ``A``/``T`` bases across *sequences*."""
    best = 0
    for seq in sequences:
        run = 0
        for base in seq.upper():
            if base in {"A", "T"}:
                run += 1
                best = max(best, run)
            else:
                run = 0
    return best


def _sample_clip_evidence(
    bam,
    chrom: str,
    center: int,
    pad_bp: int,
    max_reads: int,
    preferred_names: set[str] | None,
) -> tuple[list[int], list[int], list[str]]:
    """Collect clip-anchored junction evidence near *center*.

    Returns ``(left_anchor_positions, right_anchor_positions, clip_sequences,
    clip_reads)`` where left anchors (trailing soft clips) mark the 5' insertion
    junction and right anchors (leading soft clips) mark the 3' junction.
    """
    left_anchors: list[int] = []
    right_anchors: list[int] = []
    clip_sequences: list[str] = []
    clip_reads: list = []

    start = max(0, center - pad_bp - 1)
    stop = center + pad_bp + 1
    for read in bam.fetch(chrom, start, stop):
        if read is None or read.is_unmapped or read.query_sequence is None:
            continue
        qname = read.query_name
        if preferred_names is not None and qname is not None and qname not in preferred_names:
            continue
        if len(clip_sequences) >= max_reads:
            break
        lead, trail = _clip_lens(read)
        qs = read.query_sequence
        if trail >= _ASM_CLIP_MIN_BP and read.reference_end is not None:
            left_anchors.append(int(read.reference_end))
            clip_sequences.append(qs[len(qs) - trail :].upper())
            clip_reads.append(read)
        elif lead >= _ASM_CLIP_MIN_BP and read.reference_start is not None:
            right_anchors.append(int(read.reference_start))
            clip_sequences.append(qs[:lead].upper())
            clip_reads.append(read)

    return left_anchors, right_anchors, clip_sequences, clip_reads


def annotate_silver_with_local_assembly(
    candidates,
    *,
    disease_bam_path: Path,
    control_bam_path: Path,
    assembly_cache_dir: Path,
    mei_fasta: Path,
    reference_fasta: Path | None = None,
    interval_pad_bp: int = 250,
    retry_pad_bp: int = 600,
    max_reads_per_sample: int = 600,
    spades_threads: int = 1,
    spades_memory_gb: int = 8,
    minimap2_threads: int = 1,
    locus_workers: int = 0,
    reuse_existing: bool = True,
    reuse_cache_only: bool = False,
    disease_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
    control_preferred_read_names_by_locus: dict[tuple[str, int, int], set[str]] | None = None,
):
    """Annotate silver-stage candidate loci with De Bruijn micro-assembly.

    Args:
        candidates: Candidate locus DataFrame (chrom/window_start/window_end).
        disease_bam_path: Disease-sample BAM with soft-clip evidence.
        control_bam_path: Control-sample BAM with soft-clip evidence.
        assembly_cache_dir: Unused cache directory retained for API parity.
        mei_fasta: Unused MEI reference FASTA retained for API parity.
        reference_fasta: Optional reference genome FASTA (.fai indexed) used to
            resolve TSD boundaries precisely.
        interval_pad_bp: Flanking window radius used to recruit reads.
        retry_pad_bp: Unused retry radius retained for API parity.
        max_reads_per_sample: Cap on recruited soft-clip reads per sample.
        spades_threads: Unused SPAdes option retained for API parity.
        spades_memory_gb: Unused SPAdes option retained for API parity.
        minimap2_threads: Unused minimap2 option retained for API parity.
        locus_workers: Unused parallelization option retained for API parity.
        reuse_existing: Unused cache option retained for API parity.
        reuse_cache_only: Unused cache option retained for API parity.
        disease_preferred_read_names_by_locus: Preferred read names per locus.
        control_preferred_read_names_by_locus: Preferred read names per locus.
    """
    import pandas as pd
    import pysam

    if candidates is None or candidates.empty:
        return pd.DataFrame(columns=_ASM_OUTPUT_COLUMNS)

    silver_stage = (
        candidates["silver_stage_pass"]
        if "silver_stage_pass" in candidates.columns
        else pd.Series(False, index=candidates.index)
    )
    silver = candidates.loc[silver_stage.fillna(False).astype(bool)].copy()
    if silver.empty:
        return pd.DataFrame(columns=_ASM_OUTPUT_COLUMNS)

    assembly_cache_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    with pysam.AlignmentFile(str(disease_bam_path)) as dis_bam, pysam.AlignmentFile(
        str(control_bam_path)
    ) as ctl_bam:
        for record in silver.to_dict(orient="records"):
            chrom = str(record.get("chrom", "chr1"))
            ws = int(record.get("window_start", 0))
            we = int(record.get("window_end", ws))
            center = (ws + we) // 2
            key = (chrom, ws, we)

            left_d, right_d, clips_d, reads_d = _sample_clip_evidence(
                dis_bam,
                chrom,
                center,
                interval_pad_bp,
                max_reads_per_sample,
                disease_preferred_read_names_by_locus.get(key)
                if disease_preferred_read_names_by_locus
                else None,
            )
            left_c, right_c, clips_c, reads_c = _sample_clip_evidence(
                ctl_bam,
                chrom,
                center,
                interval_pad_bp,
                max_reads_per_sample,
                control_preferred_read_names_by_locus.get(key)
                if control_preferred_read_names_by_locus
                else None,
            )

            disease_bp = _median_or_zero(left_d) or _median_or_zero(right_d)
            control_bp = _median_or_zero(left_c) or _median_or_zero(right_c)

            disease_contigs = assemble_locus_clips(clips_d) if clips_d else []
            control_contigs = assemble_locus_clips(clips_c) if clips_c else []

            support_by_sample = {"disease": len(clips_d), "control": len(clips_c)}
            best_source = max(support_by_sample, key=lambda s: support_by_sample[s])
            consensus_bp = max(disease_bp, control_bp) if (disease_bp or control_bp) else 0
            source = best_source if support_by_sample[best_source] > 0 and consensus_bp > 0 else ""

            tsd_seq = ""
            tsd_len = 0
            clips_for_tsd = reads_d if consensus_bp == disease_bp else reads_c
            if consensus_bp > 0 and clips_for_tsd and reference_fasta is not None:
                from retro_miner.tsd_refiner import refine_tsd_boundaries

                tsd_result = refine_tsd_boundaries(
                    chrom,
                    consensus_bp,
                    clips_for_tsd,
                    reference_fasta,
                    min_tsd_len=4,
                    max_tsd_len=40,
                )
                tsd_seq = tsd_result.tsd_seq
                tsd_len = tsd_result.tsd_length

            rows.append(
                {
                    "chrom": chrom,
                    "window_start": ws,
                    "window_end": we,
                    "asm_status": "assembled" if consensus_bp > 0 else "no_clips",
                    "asm_interval_start": center - interval_pad_bp,
                    "asm_interval_end": center + interval_pad_bp,
                    "asm_retry_used": False,
                    "asm_disease_reads_extracted": len(clips_d),
                    "asm_control_reads_extracted": len(clips_c),
                    "asm_disease_contig_count": len(disease_contigs),
                    "asm_control_contig_count": len(control_contigs),
                    "asm_disease_max_contig_len": (
                        max(c.length for c in disease_contigs) if disease_contigs else 0
                    ),
                    "asm_control_max_contig_len": (
                        max(c.length for c in control_contigs) if control_contigs else 0
                    ),
                    "asm_error_message": "",
                    "asm_disease_primary_contig_id": (
                        f"contig_{0}" if disease_contigs else ""
                    ),
                    "asm_control_primary_contig_id": (
                        f"contig_{0}" if control_contigs else ""
                    ),
                    "asm_disease_breakpoint_pos": disease_bp,
                    "asm_control_breakpoint_pos": control_bp,
                    "asm_consensus_breakpoint_pos": consensus_bp,
                    "asm_breakpoint_source": source,
                    "asm_tsd_seq": tsd_seq,
                    "asm_tsd_len": int(tsd_len),
                    "asm_polyA_max_run": _longest_at_run(clips_d + clips_c),
                    "asm_insertion_length": 0,
                    "asm_insertion_length_observed": False,
                    "asm_insertion_length_imputed": False,
                    "asm_insertion_length_confidence_tier": "",
                    "asm_insertion_mei_start": 0,
                    "asm_insertion_mei_end": 0,
                    "asm_mei_family": "",
                    "asm_mei_subfamily": "",
                    "asm_mei_target_length": 0,
                    "asm_insertion_orientation": "",
                    "asm_consensus_primary_contig_id": "",
                    "asm_complex_class": "",
                    "asm_non_mei_partner_chrom": "",
                    "asm_non_mei_partner_pos": 0,
                    "asm_non_mei_partner_type": "",
                    "asm_breakpoint_side_status": "",
                    "asm_complexity_source": "",
                    "asm_top_contigs": "",
                    "asm_mei_alignment_preset": "",
                    "asm_left_support_contig_id": "",
                    "asm_right_support_contig_id": "",
                    "asm_left_support_mei_start": 0,
                    "asm_left_support_mei_end": 0,
                    "asm_right_support_mei_start": 0,
                    "asm_right_support_mei_end": 0,
                    "asm_left_support_mei_aln_len": 0,
                    "asm_right_support_mei_aln_len": 0,
                    "asm_microhomology_sequence": "",
                    "asm_junction_overlap_sequence": "",
                    "asm_coord_model": "",
                    "asm_coord_logic_version": 0,
                    "asm_adaptive_read_cap_rerun": False,
                    "asm_adaptive_read_cap_target": 0,
                    "asm_disease_recruited_evidence_read_names": "",
                    "asm_control_recruited_evidence_read_names": "",
                }
            )

    if not rows:
        return pd.DataFrame(columns=_ASM_OUTPUT_COLUMNS)
    return pd.DataFrame(rows)