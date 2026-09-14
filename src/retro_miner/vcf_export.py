"""VCF v4.3 Exporter for Mobile Element Insertions (MEIs).

Serializes candidate MEI loci enriched with Bayesian genotype calls, TSD
refinement metrics, and diagnostic subfamily classifications into standard
VCF v4.3 structural variant callsets.

Candidate records are Python dicts that may carry either the refined typed
objects (``GenotypeCall``, ``TSDResult``, ``SubfamilyCall``) or the legacy
flat keys (``genotype``/``vaf``/``genotype_quality``, ``tsd_seq``/``tsd_len``/
``poly_a_detected``, ``subfamily``/``mei_llr``).  Typed objects take
precedence; flat keys are used as a fallback so both interfaces produce
identical output.

Literature Anchors: Li et al. (2011) / VCF v4.3 Spec / Gardner et al. (2017) MELT format.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from retro_miner.genotyper import GenotypeCall
from retro_miner.subfamily_voter import SubfamilyCall
from retro_miner.tsd_refiner import TSDResult

VCF_META_HEADER = "##fileformat=VCFv4.3"
VCF_FILE_DATE = "##fileDate=20260913"
VCF_SOURCE = "##source=retrotransposon-miner-v1.0"

ALT_DEFINITIONS: dict[str, str] = {
    "L1HS": '##ALT=<ID=INS:MEI:L1HS,Description="Line-1 Mobile Element Insertion (L1HS)">',
    "ALU": '##ALT=<ID=INS:MEI:ALU,Description="Alu Mobile Element Insertion">',
    "SVA": '##ALT=<ID=INS:MEI:SVA,Description="SVA Mobile Element Insertion">',
}

INFO_HEADER_LINES: list[str] = [
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=MEI_TYPE,Number=1,Type=String,Description="Mobile element insertion subfamily">',
    '##INFO=<ID=TSD,Number=1,Type=String,Description="Target Site Duplication sequence">',
    '##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description="Target Site Duplication length in bp">',
    '##INFO=<ID=POLYA,Number=1,Type=Integer,Description="Flag indicating 3\' poly(A) tail detection (1=True, 0=False)">',
    '##INFO=<ID=MEI_LLR,Number=1,Type=Float,Description="Log-likelihood ratio for subfamily call">',
    '##INFO=<ID=SUBFAM,Number=1,Type=String,Description="Classified mobile element subfamily">',
    '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of supporting non-reference reads">',
    '##INFO=<ID=TRANSDUCTION_TYPE,Number=1,Type=String,Description="3\' transduction classification (3_PRIME_PARTNERED, 3_PRIME_ORPHAN, NONE)">',
    '##INFO=<ID=TRANSDUCTION_SEQ,Number=1,Type=String,Description="Transduced genomic sequence barcode">',
    '##INFO=<ID=TRANSDUCTION_LENGTH,Number=1,Type=Integer,Description="Length of transduced sequence in bp">',
    '##INFO=<ID=TPRT_MOTIF_SCORE,Number=1,Type=Float,Description="TPRT endonuclease cleavage motif confidence score">',
]

FILTER_HEADER_LINES: list[str] = [
    '##FILTER=<ID=PASS,Description="High-confidence structural variant call">',
    '##FILTER=<ID=LowQual,Description="Genotype quality below threshold (GQ < 20.0) or insufficient alt support (< 3 reads)">',
]

FORMAT_HEADER_LINES: list[str] = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">',
    '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for reference and alternate alleles">',
]

MIN_PASS_GQ = 20.0
MIN_PASS_SUPPORT = 3


def _alt_label(family: str) -> str:
    """Map a family hint ("L1", "Alu", "SVA") to the symbolic ALT label."""
    fam = (family or "").upper()
    if "ALU" in fam:
        return "ALU"
    if "SVA" in fam:
        return "SVA"
    return "L1HS"


def _extract_genotype_fields(rec: dict[str, Any]) -> tuple[str, float, float]:
    """Return ``(genotype, vaf, genotype_quality)`` from a typed or flat record.

    In the flat fallback, a present-but-``None`` key (as produced by pandas
    ``DataFrame.to_dict("records")`` on empty cells) is treated as absent and
    falls back to the default.  Explicit zeros (``vaf=0.0``, ``gq=0.0``) are
    preserved via an ``is None`` check rather than falsiness.
    """
    call = rec.get("genotype_call")
    if isinstance(call, GenotypeCall):
        return call.genotype, float(call.vaf), float(call.genotype_quality)
    genotype = rec.get("genotype", "0/1")
    vaf = rec.get("vaf", 0.50)
    gq = rec.get("genotype_quality", 30.0)
    if genotype is None:
        genotype = "0/1"
    if vaf is None:
        vaf = 0.50
    if gq is None:
        gq = 30.0
    return str(genotype), float(vaf), float(gq)


def _extract_subfamily_fields(rec: dict[str, Any]) -> tuple[str, float]:
    """Return ``(subfamily, log_likelihood_ratio)`` from a typed or flat record."""
    call = rec.get("subfamily_call")
    if isinstance(call, SubfamilyCall):
        return call.top_subfamily, float(call.log_likelihood_ratio)
    subfamily = rec.get("subfamily", "L1HS")
    llr = rec.get("mei_llr", 3.0)
    if subfamily is None:
        subfamily = "L1HS"
    if llr is None:
        llr = 3.0
    return str(subfamily), float(llr)


def _extract_tsd_fields(rec: dict[str, Any]) -> tuple[str, int, int]:
    """Return ``(tsd_seq, tsd_length, poly_a_flag)`` from a typed or flat record."""
    result = rec.get("tsd_result")
    if isinstance(result, TSDResult):
        tsd_seq = result.tsd_seq
        tsd_len = int(result.tsd_length)
        poly_a = 1 if result.polyA_tail_detected else 0
    else:
        tsd_seq = rec.get("tsd_seq", "")
        tsd_len = rec.get("tsd_len", len(tsd_seq) if tsd_seq else 0)
        poly_a_raw = rec.get("poly_a_detected")
        if tsd_seq is None:
            tsd_seq = ""
        if tsd_len is None:
            tsd_len = len(tsd_seq)
        if poly_a_raw is None:
            poly_a_raw = False
        tsd_seq = str(tsd_seq)
        tsd_len = int(tsd_len)
        poly_a = 1 if bool(poly_a_raw) else 0
    return tsd_seq, tsd_len, poly_a


def _extract_transduction_fields(rec: dict[str, Any]) -> tuple[str, str, int]:
    """Return ``(transduction_type, transduction_seq, transduction_length)``."""
    return (
        str(rec.get("transduction_type", "NONE")),
        str(rec.get("transduction_seq", "")),
        int(rec.get("transduction_length", 0)),
    )


def _extract_tprt_motif_score(rec: dict[str, Any]) -> float:
    """Return TPRT cleavage motif score."""
    return float(rec.get("tprt_motif_score", 0.0))


def _extract_support(rec: dict[str, Any]) -> int:
    """Alt-supporting read count, from ``support`` or the legacy ``k_alt`` key."""
    value = rec.get("support")
    if value is None:
        value = rec.get("k_alt", 0)
    if value is None:
        value = 0
    return int(value)


def _extract_depth(rec: dict[str, Any]) -> tuple[int, int]:
    """Return ``(k_ref, k_alt)`` allele depths used to populate AD."""
    k_ref = rec.get("k_ref", 10)
    k_alt = rec.get("k_alt", 10)
    if k_ref is None:
        k_ref = 10
    if k_alt is None:
        k_alt = 10
    return int(k_ref), int(k_alt)


def _build_header_lines(records: list[dict[str, Any]], sample_name: str) -> list[str]:
    """Assemble the VCF header, including per-record contigs and used ALT symbols."""
    chrom_seen: set[str] = set()
    alt_seen: set[str] = set()
    contig_lines: list[str] = []
    alt_lines: list[str] = []
    for rec in records:
        chrom = str(rec.get("chrom", "chr1"))
        if chrom not in chrom_seen:
            chrom_seen.add(chrom)
            contig_lines.append(f"##contig=<ID={chrom}>")
        label = _alt_label(str(rec.get("family", "L1")))
        if label not in alt_seen:
            alt_seen.add(label)
            alt_lines.append(ALT_DEFINITIONS[label])

    lines = [VCF_META_HEADER, VCF_FILE_DATE, VCF_SOURCE]
    lines.extend(contig_lines)
    lines.extend(alt_lines)
    lines.extend(INFO_HEADER_LINES)
    lines.extend(FILTER_HEADER_LINES)
    lines.extend(FORMAT_HEADER_LINES)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name)
    return lines


def write_mei_vcf(
    candidate_records: list[dict[str, Any]],
    output_path: str | Path,
    sample_name: str = "SAMPLE",
) -> Path:
    """Export candidate MEI records to a VCF v4.3 file.

    Each record may carry typed ``GenotypeCall``/``TSDResult``/``SubfamilyCall``
    objects (preferred) or the legacy flat keys.  A record passes ``FILTER=PASS``
    only when genotype quality is >= 20.0 *and* at least 3 supporting
    non-reference reads are called; otherwise it is ``LowQual``.

    Args:
        candidate_records: List of dicts with locus, GT, TSD, and subfamily calls.
        output_path: Path to output VCF file.
        sample_name: Individual sample identifier for VCF column header.

    Returns:
        Path object of written VCF file.
    """
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    lines = _build_header_lines(candidate_records, sample_name)

    for rec in candidate_records:
        chrom = str(rec.get("chrom", "chr1"))
        pos = int(rec.get("pos", 10000))
        mei_id = f"MEI_{chrom}_{pos}"
        ref_base = rec.get("ref_base", "N")
        if ref_base is None:
            ref_base = "N"
        ref_base = str(ref_base)[0].upper()
        family = str(rec.get("family", "L1"))

        subfamily, llr = _extract_subfamily_fields(rec)
        gt, vaf, gq = _extract_genotype_fields(rec)
        tsd_seq, tsd_len, poly_a = _extract_tsd_fields(rec)
        support = _extract_support(rec)
        k_ref, k_alt = _extract_depth(rec)
        tr_type, tr_seq, tr_len = _extract_transduction_fields(rec)
        tprt_score = _extract_tprt_motif_score(rec)

        alt_symbol = f"<INS:MEI:{_alt_label(family)}>"
        qual_str = f"{gq:.1f}" if gq > 0 else "."
        filter_str = "PASS" if gq >= MIN_PASS_GQ and support >= MIN_PASS_SUPPORT else "LowQual"

        info_fields = [
            "SVTYPE=INS",
            f"MEI_TYPE={subfamily}",
            f"TSD={tsd_seq if tsd_seq else 'NONE'}",
            f"TSDLEN={tsd_len}",
            f"POLYA={poly_a}",
            f"MEI_LLR={llr:.2f}",
            f"SUBFAM={subfamily}",
            f"SUPPORT={support}",
            f"TRANSDUCTION_TYPE={tr_type}",
            f"TRANSDUCTION_SEQ={tr_seq if tr_seq else 'NONE'}",
            f"TRANSDUCTION_LENGTH={tr_len}",
            f"TPRT_MOTIF_SCORE={tprt_score:.2f}",
        ]
        info_str = ";".join(info_fields)

        sample_vals = f"{gt}:{gq:.1f}:{vaf:.2f}:{k_ref},{k_alt}"
        line = (
            f"{chrom}\t{pos}\t{mei_id}\t{ref_base}\t{alt_symbol}\t{qual_str}\t"
            f"{filter_str}\t{info_str}\tGT:GQ:VAF:AD\t{sample_vals}"
        )
        lines.append(line)

    out_file.write_text("\n".join(lines) + "\n")
    return out_file