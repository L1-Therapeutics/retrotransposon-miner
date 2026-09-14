"""Unified Scientific MEI Discovery Pipeline Engine.

Integrates streaming BAM candidate extraction, De Bruijn micro-assembly,
Shannon-entropy TSD refinement, k-mer subfamily voting, Bayesian diploid
genotyping, and VCF v4.3 output generation into an end-to-end workflow.

Supports optional spatial-index-guided O(log N) locus fetching when a BAM
index (``.bai`` / ``.csi``) is present, falling back to sequential streaming
for unindexed alignments.

Literature Anchors: Gardner et al. (2017) / Layer et al. (2014) / Cameron et al. (2017) / Li et al. (2009) / Li et al. (2011).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from retro_miner.enb_independent import classify_en_independent_event
from retro_miner.genotyper import calculate_mei_genotype
from retro_miner.indexed_stream import LocusEvidence, fetch_locus_spanning_pairs, has_bam_index
from retro_miner.local_assembly import assemble_locus_clips
from retro_miner.somatic_em import SomaticCall, SomaticEMClassifier
from retro_miner.subfamily_voter import classify_mei_subfamily
from retro_miner.tsd_refiner import refine_tsd_boundaries
from retro_miner.vcf_export import write_mei_vcf


def _tag_somatic_fields(record: dict[str, Any], call: SomaticCall) -> dict[str, Any]:
    """Attach EM somatic-mixture classification fields onto a locus record."""
    record["somatic_call"] = call
    record["is_somatic"] = call.is_somatic
    record["somatic_posterior"] = call.somatic_posterior
    record["subclone_vaf"] = call.subclone_vaf
    record["classification"] = call.classification
    return record


def _extract_clips_from_evidence(evidence: LocusEvidence) -> tuple[list[str], list[str]]:
    """Split soft-clip sequences into left/right anchor lists."""
    left: list[str] = []
    right: list[str] = []
    for seq in evidence.soft_clips:
        if len(seq) >= 4:
            left.append(seq)
            right.append(seq[::-1])
    return left, right


def process_candidate_locus(
    chrom: str,
    pos: int,
    left_clips: list[str],
    right_clips: list[str],
    k_alt: int,
    k_ref: int,
    family_hint: str = "L1",
    bam_path: Path | None = None,
    locus_start: int | None = None,
    locus_end: int | None = None,
    somatic_classifier: SomaticEMClassifier | None = None,
    en_score: float | None = None,
    target_del_bp: int = 0,
) -> dict[str, Any]:
    """Perform full scientific refinement and genotyping for a single candidate locus.

    When *bam_path* points to an indexed BAM and *locus_start*/*locus_end* are
    provided, evidence reads are fetched via spatial index rather than using
    pre-computed clip/count inputs.

    When *somatic_classifier* has already been fit (e.g. on the full locus
    batch by :func:`run_scientific_mei_pipeline`), the locus is additionally
    tagged with EM subclonal-somatic mixture fields, including the somatic
    posterior and inferred subclone VAF.

    Args:
        chrom: Reference chromosome.
        pos: Breakpoint coordinate.
        left_clips: 5' soft-clipped sequences.
        right_clips: 3' soft-clipped sequences.
        k_alt: Supporting non-reference evidence read count.
        k_ref: Spanning concordant reference read count.
        family_hint: Targeted MEI family ("L1", "Alu", or "SVA").
        bam_path: Optional path to indexed BAM for on-demand evidence fetch.
        locus_start: Optional 0-based locus start for indexed fetch.
        locus_end: Optional 0-based locus end for indexed fetch.
        somatic_classifier: Optional fitted SomaticEMClassifier for locus tagging.
        en_score: Optional normalized EN cleavage motif PWM score in
            ``[0.0, 1.0]`` (default 0.0 = no canonical motif evidenced).
        target_del_bp: Target-site genomic deletion size at the locus in bp.

    Returns:
        Enriched locus dictionary ready for VCF v4.3 export.
    """
    all_clips = left_clips + right_clips
    l_clip_concat = "".join(left_clips)
    r_clip_concat = "".join(right_clips)

    if bam_path is not None and locus_start is not None and locus_end is not None:
        if has_bam_index(bam_path):
            try:
                evidence = fetch_locus_spanning_pairs(
                    bam_path=bam_path,
                    chrom=chrom,
                    start=locus_start,
                    end=locus_end,
                )
                k_ref = evidence.spanning_pairs
                k_alt = max(k_alt, len(evidence.soft_clips))
                left_clips, right_clips = _extract_clips_from_evidence(evidence)
                all_clips = left_clips + right_clips
                l_clip_concat = "".join(left_clips)
                r_clip_concat = "".join(right_clips)
            except FileNotFoundError:
                pass

    # 1. Target Site Duplication (TSD) & Poly(A) Refinement
    tsd_res = refine_tsd_boundaries(l_clip_concat, r_clip_concat)

    # 2. Subfamily Voting (k-mer LLR)
    subfam_call = classify_mei_subfamily(all_clips, family_hint=family_hint)

    # 3. De Bruijn Micro-Assembly
    assembled_contigs = assemble_locus_clips(all_clips, k=15, min_coverage=2)
    top_contig = assembled_contigs[0].sequence if assembled_contigs else ""

    # 4. Bayesian Diploid Genotyping
    gt_call = calculate_mei_genotype(k_alt=k_alt, k_ref=k_ref)

    record: dict[str, Any] = {
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
        "poly_a_detected": getattr(
            tsd_res,
            "polyA_tail_detected",
            getattr(tsd_res, "poly_a_detected", False),
        ),
        "mei_llr": subfam_call.log_likelihood_ratio,
        "assembled_contig": top_contig,
    }

    # 5. EN-Independent (DSB Repair) Integration Classification
    en_call = classify_en_independent_event(
        tsd_length=int(tsd_res.tsd_length),
        en_score=0.0 if en_score is None else float(en_score),
        span_deletion_bp=int(target_del_bp),
    )
    record["is_en_independent"] = en_call.is_en_independent
    record["target_del_bp"] = en_call.deletion_size

    # 6. EM Subclonal Somatic Mixture Scoring
    if somatic_classifier is not None and somatic_classifier.is_fitted:
        n_total = k_alt + k_ref
        somatic_call = somatic_classifier.predict(k_alt=k_alt, n=n_total)
        _tag_somatic_fields(record, somatic_call)

    return record


def run_scientific_mei_pipeline(
    raw_candidate_loci: list[dict[str, Any]],
    output_vcf_path: str | Path,
    sample_name: str = "SAMPLE",
    bam_path: Path | None = None,
) -> Path:
    """Execute end-to-end candidate locus refinement and export to VCF v4.3.

    Args:
        raw_candidate_loci: List of unrefined candidate dictionary loci.
        output_vcf_path: Destination path for generated VCF v4.3 file.
        sample_name: Individual sample header identifier.
        bam_path: Optional path to indexed BAM for spatial-index evidence fetch.

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
            bam_path=bam_path,
            locus_start=int(loc.get("locus_start", loc.get("pos", 10000))),
            locus_end=int(loc.get("locus_end", loc.get("pos", 10000) + 200)),
            en_score=loc.get("en_score"),
            target_del_bp=int(loc.get("target_del_bp") or 0),
        )
        refined_records.append(refined)

    # Batch EM subclonal-somatic mixture fit over all refined loci, then tag each
    somatic_em = SomaticEMClassifier()
    em_calls = somatic_em.fit_predict(
        k_alt_list=[rec["k_alt"] for rec in refined_records],
        n_list=[rec["k_alt"] + rec["k_ref"] for rec in refined_records],
    )
    for rec, call in zip(refined_records, em_calls, strict=False):
        _tag_somatic_fields(rec, call)

    return write_mei_vcf(refined_records, output_vcf_path, sample_name=sample_name)

