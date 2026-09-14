"""Read-Spanning Haplotype Linkage & Phase Estimator for Mobile Element Insertions.

Physically links supporting evidence reads (MEI-insertion reads vs concordant
reference reads) to nearby heterozygous germline SNPs within a user-specified
window around the insertion locus.  Reads overlapping both the MEI breakpoint
and a SNP establish phase blocks (``H1`` vs ``H2``), enabling haplotype-aware
VAF calculation and parent/chromosomal-lineage-of-origin determination for
somatic and germline MEIs.

Phase confidence is reported Phred-scaled:

    PQ = -10 * log10(P_error)

where ``P_error`` is the two-sided Fisher exact p-value of the ``2 x 2``
contingency of alternate-vs-reference reads per hapotype group.  Without any
informative heterozygous SNP, the locus is left ``UNPHASED`` (``PQ = 0.0``).

Literature Anchors: Delaneau et al. (2019) Nat Genet 51(5):747-751 /
Martin et al. (2016) Genome Biol 17(1):230.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pysam

PHASE_BLOCK_H1 = "H1"
PHASE_BLOCK_H2 = "H2"
PHASE_UNPHASED = "UNPHASED"

MIN_INFORMATIVE_READS = 2
PHASE_QUALITY_PVALUE_THRESHOLD = 0.05


@dataclass(frozen=True)
class PhaseLinkageResult:
    """Haplotype linkage of an MEI locus to nearby heterozygous SNPs.

    Attributes:
        haplotype_assigned: ``"H1"``, ``"H2"``, or ``"UNPHASED"``.
        linked_snps_count: Number of heterozygous SNPs with phasing evidence.
        phase_confidence: Phred-scaled phase quality ``PQ`` (0.0 when unphased).
        phase_block_id: Genomic coordinates of the phased window
            (e.g. ``"chr1:100450-102450"``).
    """

    haplotype_assigned: str
    linked_snps_count: int
    phase_confidence: float
    phase_block_id: str


def _unphased_result(window_start: int, window_end: int) -> PhaseLinkageResult:
    return PhaseLinkageResult(
        haplotype_assigned=PHASE_UNPHASED,
        linked_snps_count=0,
        phase_confidence=0.0,
        phase_block_id=f"{window_start}-{window_end}",
    )


def _log_hypergeom(x: int, row1: int, row2: int, col1: int, total: int) -> float:
    """Log-probability of a ``2 x 2`` table under fixed margins (Fisher exact)."""
    if x < max(0, col1 - row2) or x > min(row1, col1) or x < 0:
        return -math.inf
    col2 = total - col1
    if col1 - x < 0 or col1 - x > row2 or row1 - x < 0:
        return -math.inf
    log_num = (
        math.lgamma(row1 + 1)
        - math.lgamma(x + 1)
        - math.lgamma(row1 - x + 1)
        + math.lgamma(row2 + 1)
        - math.lgamma(col1 - x + 1)
        - math.lgamma(row2 - (col1 - x) + 1)
    )
    log_den = math.lgamma(total + 1) - math.lgamma(col1 + 1) - math.lgamma(col2 + 1)
    return log_num - log_den


def _fisher_two_sided_pvalue(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for table ``[[a, b], [c, d]]``."""
    total = a + b + c + d
    if total == 0:
        return 1.0
    row1 = a + b
    col1 = a + c
    log_p_obs = _log_hypergeom(a, row1, c + d, col1, total)
    p_value = 0.0
    for x in range(max(0, col1 - (c + d)), min(row1, col1) + 1):
        lp = _log_hypergeom(x, row1, c + d, col1, total)
        if lp <= log_p_obs:
            p_value += math.exp(lp)
    return min(p_value, 1.0)


def _read_spans_breakpoint(read, pos: int) -> bool:
    """True when the read's aligned reference span covers the breakpoint."""
    start = read.reference_start
    end = read.reference_end
    if start is None or end is None:
        return False
    return start <= pos < end


def _read_carries_mei_evidence(read) -> bool:
    """True when the read carries an insertion/soft-clip consistent with an MEI."""
    if read.cigartuples is None:
        return False
    return any(op in (1, 4, 5) for op, _ in read.cigartuples)


def _allele_on_read(read, snp_pos: int) -> str | None:
    """Return the read's base (uppercased) at *snp_pos*, or None if deleted/absent."""
    query_seq = read.query_sequence
    if not query_seq:
        return None
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if ref_pos == snp_pos:
            if query_pos is None:
                return None
            return query_seq[query_pos].upper()
    return None


def _collect_phase_counts(
    reads,
    pos: int,
    het_snps: list[tuple[int, str, str]],
) -> tuple[list[int], list[tuple[int, int, int, int]]]:
    """Return ``(informative_snp_count, oriented_2x2_tables)`` per SNP.

    Each table is ``(mei_alt, mei_ref, ref_alt, ref_ref)`` oriented so that
    ``alt`` is the alternate allele of that SNP.
    """
    informative: list[int] = []
    tables: list[tuple[int, int, int, int]] = []
    for snp_pos, snp_ref, snp_alt in het_snps:
        mei_alt = mei_ref = ref_alt = ref_ref = 0
        for read in reads:
            if not _read_spans_breakpoint(read, pos):
                continue
            allele = _allele_on_read(read, snp_pos)
            if allele not in (snp_ref.upper(), snp_alt.upper()):
                continue
            if _read_carries_mei_evidence(read):
                if allele == snp_alt.upper():
                    mei_alt += 1
                else:
                    mei_ref += 1
            else:
                if allele == snp_alt.upper():
                    ref_alt += 1
                else:
                    ref_ref += 1
        if mei_alt + mei_ref >= MIN_INFORMATIVE_READS:
            informative.append(snp_pos)
            tables.append((mei_alt, mei_ref, ref_alt, ref_ref))
    return informative, tables


def phase_mei_locus(
    bam_path: str,
    chrom: str,
    pos: int,
    het_snps: list[tuple[int, str, str]],
    window_bp: int = 1000,
) -> PhaseLinkageResult:
    """Phase an MEI breakpoint against nearby heterozygous SNPs.

    Reads spanning ``[pos - window_bp, pos + window_bp]`` are collected with
    :class:`pysam.AlignmentFile`.  Reads that both overlap the MEI breakpoint
    and carry one of the SNP alleles contribute to a per-SNP ``2 x 2``
    contingency of MEI-insertion vs reference reads per allele.  Tables are
    pooled across consistently-oriented SNPs and scored with Fisher's exact
    test to yield the Phred-scaled phase quality ``PQ``.

    Args:
        bam_path: Path to an indexed BAM file of the sample.
        chrom: Reference chromosome of the insertion.
        pos: MEI insertion breakpoint coordinate (0-based).
        het_snps: Heterozygous SNP loci as ``(position, ref_base, alt_base)``.
        window_bp: Half-width of the phasing window in bp around *pos*.

    Returns:
        :class:`PhaseLinkageResult` with haplotype assignment, linked SNP
        count, Phred-scaled phase quality, and phase block coordinates.
    """
    window_start = max(0, pos - window_bp)
    window_end = pos + window_bp
    blocked = _unphased_result(window_start, window_end)

    bam_file = Path(bam_path) if bam_path else None
    if bam_file is None or not bam_file.exists():
        return blocked
    if not het_snps:
        return _unphased_result(window_start, window_end)

    try:
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            reads = list(bam.fetch(chrom, max(0, window_start), window_end))
    except (OSError, ValueError):
        return blocked

    informative_snps, tables = _collect_phase_counts(reads, pos, het_snps)
    if not informative_snps:
        return _unphased_result(window_start, window_end)

    mei_alt = sum(t[0] for t in tables)
    mei_ref = sum(t[1] for t in tables)
    ref_alt = sum(t[2] for t in tables)
    ref_ref = sum(t[3] for t in tables)

    if mei_alt + mei_ref == 0 or mei_alt == mei_ref:
        return _unphased_result(window_start, window_end)

    p_value = _fisher_two_sided_pvalue(mei_alt, mei_ref, ref_alt, ref_ref)
    pq = -10.0 * math.log10(max(p_value, 1e-100))
    phase_block_id = f"{chrom}:{window_start}-{window_end}"

    if p_value > PHASE_QUALITY_PVALUE_THRESHOLD:
        return PhaseLinkageResult(
            haplotype_assigned=PHASE_UNPHASED,
            linked_snps_count=len(informative_snps),
            phase_confidence=round(pq, 2),
            phase_block_id=phase_block_id,
        )

    haplotype = PHASE_BLOCK_H1 if mei_alt > mei_ref else PHASE_BLOCK_H2
    return PhaseLinkageResult(
        haplotype_assigned=haplotype,
        linked_snps_count=len(informative_snps),
        phase_confidence=round(pq, 2),
        phase_block_id=phase_block_id,
    )