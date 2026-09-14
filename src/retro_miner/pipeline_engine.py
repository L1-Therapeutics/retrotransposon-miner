"""Unified Scientific MEI Discovery Pipeline Engine.

Integrates streaming BAM candidate extraction, De Bruijn micro-assembly,
Shannon-entropy TSD refinement, k-mer subfamily voting, Bayesian diploid
genotyping, and VCF v4.3 output generation into an end-to-end workflow.

Literature Anchors: Gardner et al. (2017) / Layer et al. (2014) / Cameron et al. (2017) / Li et al. (2011).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from retro_miner.genotyper import calculate_mei_genotype
from retro_miner.local_assembly import assemble_locus_clips
from retro_miner.subfamily_voter import classify_mei_subfamily
from retro_miner.tsd_refiner import refine_tsd_boundaries
from retro_miner.vcf_export import write_mei_vcf


def process_candidate_locus(
    chrom: str,
    pos: int,
    left_clips: list[str],
    right_clips: list[str],
    k_alt: int,
    k_ref: int,
    family_hint: str = "L1",
) -> dict[str, Any]:
    """Perform full scientific refinement and genotyping for a single candidate locus.

    Args:
        chrom: Reference chromosome.
        pos: Breakpoint coordinate.
        left_clips: 5' soft-clipped sequences.
        right_clips: 3' soft-clipped sequences.
        k_alt: Supporting non-reference evidence read count.
        k_ref: Spanning concordant reference read count.
        family_hint: Targeted MEI family ("L1", "Alu", or "SVA").

    Returns:
        Enriched locus dictionary ready for VCF v4.3 export.
    """
    all_clips = left_clips + right_clips
    l_clip_concat = "".join(left_clips)
    r_clip_concat = "".join(right_clips)

    # 1. Target Site Duplication (TSD) & Poly(A) Refinement
    tsd_res = refine_tsd_boundaries(l_clip_concat, r_clip_concat)

    # 2. Subfamily Voting (k-mer LLR)
    subfam_call = classify_mei_subfamily(all_clips, family_hint=family_hint)

    # 3. De Bruijn Micro-Assembly
    assembled_contigs = assemble_locus_clips(all_clips, k=15, min_coverage=2)
    top_contig = assembled_contigs[0].sequence if assembled_contigs else ""

    # 4. Bayesian Diploid Genotyping
    gt_call = calculate_mei_genotype(k_alt=k_alt, k_ref=k_ref)

    return {
        "chrom": chrom,
        "pos": pos,
        "ref_base": "N",
        "family": family_hint,
        "subfamily": subfam_call.top_subfamily,
        "genotype": gt_call.genotype,
        "genotype_quality": gt_call.genotype_quality,
        "vaf": gt_call.vaf,
        "k_ref": k_ref,
        "k_alt": k_alt,
        "tsd_seq": tsd_res.tsd_seq,
        "tsd_len": tsd_res.tsd_length,
        "poly_a_detected": tsd_res.polyA_tail_detected,
        "mei_llr": subfam_call.log_likelihood_ratio,
        "assembled_contig": top_contig,
    }


def run_scientific_mei_pipeline(
    raw_candidate_loci: list[dict[str, Any]],
    output_vcf_path: str | Path,
    sample_name: str = "SAMPLE",
) -> Path:
    """Execute end-to-end candidate locus refinement and export to VCF v4.3.

    Args:
        raw_candidate_loci: List of unrefined candidate dictionary loci.
        output_vcf_path: Destination path for generated VCF v4.3 file.
        sample_name: Individual sample header identifier.

    Returns:
        Path object pointing to finalized VCF file.
    """
    refined_records = []
    for loc in raw_candidate_loci:
        refined = process_candidate_locus(
            chrom=loc.get("chrom", "chr1"),
            pos=int(loc.get("pos", 10000)),
            left_clips=loc.get("left_clips", []),
            right_clips=loc.get("right_clips", []),
            k_alt=int(loc.get("k_alt", 5)),
            k_ref=int(loc.get("k_ref", 5)),
            family_hint=loc.get("family", "L1"),
        )
        refined_records.append(refined)

    return write_mei_vcf(refined_records, output_vcf_path, sample_name=sample_name)
