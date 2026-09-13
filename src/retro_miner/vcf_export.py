"""VCF v4.3 Exporter for Mobile Element Insertions (MEIs).

Serializes candidate MEI loci enriched with Bayesian genotype calls, TSD
refinement metrics, and diagnostic subfamily classifications into standard
VCF v4.3 structural variant callsets.

Literature Anchors: Li et al. (2011) / VCF v4.3 Spec / Gardner et al. (2017) MELT format.
"""

from __future__ import annotations

import re
from pathlib import Path

VCF_HEADER_TEMPLATE = """##fileformat=VCFv4.3
##fileDate=20260913
##source=retrotransposon-miner-v1.0
##ALT=<ID=INS:MEI:L1HS,Description="Line-1 Mobile Element Insertion (L1HS)">
##ALT=<ID=INS:MEI:ALU,Description="Alu Mobile Element Insertion">
##ALT=<ID=INS:MEI:SVA,Description="SVA Mobile Element Insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MEI_TYPE,Number=1,Type=String,Description="Mobile element insertion subfamily">
##INFO=<ID=TSD,Number=1,Type=String,Description="Target Site Duplication sequence">
##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description="Target Site Duplication length in bp">
##INFO=<ID=POLYA,Number=1,Type=Integer,Description="Flag indicating 3' poly(A) tail detection (1=True, 0=False)">
##INFO=<ID=MEI_LLR,Number=1,Type=Float,Description="Log-likelihood ratio for subfamily call">
##FILTER=<ID=PASS,Description="High-confidence structural variant call">
##FILTER=<ID=LowQual,Description="Genotype quality below threshold (GQ < 20.0)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for reference and alternate alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}
"""


def write_mei_vcf(
    candidate_records: list[dict],
    output_path: str | Path,
    sample_name: str = "SAMPLE",
) -> Path:
    """Export candidate MEI records to a VCF v4.3 file.

    Args:
        candidate_records: List of dictionaries containing locus metrics, GT, TSD, and Subfamily calls.
        output_path: Path to output VCF file.
        sample_name: Individual sample identifier for VCF column header.

    Returns:
        Path object of written VCF file.
    """
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    header_text = VCF_HEADER_TEMPLATE.format(sample_name=sample_name)

    lines = [header_text.strip()]

    for rec in candidate_records:
        chrom = rec.get("chrom", "chr1")
        pos = rec.get("pos", 10000)
        mei_id = f"MEI_{chrom}_{pos}"
        ref_base = rec.get("ref_base", "N")

        family = rec.get("family", "L1").upper()
        subfamily = rec.get("subfamily", "L1HS")

        alt_symbol = f"<INS:MEI:{'L1HS' if 'L1' in family else 'ALU' if 'ALU' in family else 'SVA'}>"

        gq = float(rec.get("genotype_quality", 30.0))
        qual_str = f"{gq:.1f}" if gq > 0 else "."
        filter_str = "PASS" if gq >= 20.0 else "LowQual"

        tsd_seq = rec.get("tsd_seq", "")
        tsd_len = int(rec.get("tsd_len", len(tsd_seq)))
        poly_a = 1 if rec.get("poly_a_detected", False) else 0
        llr = float(rec.get("mei_llr", 3.0))

        info_fields = [
            "SVTYPE=INS",
            f"MEI_TYPE={subfamily}",
            f"TSD={tsd_seq if tsd_seq else 'NONE'}",
            f"TSDLEN={tsd_len}",
            f"POLYA={poly_a}",
            f"MEI_LLR={llr:.2f}",
        ]
        info_str = ";".join(info_fields)

        gt = rec.get("genotype", "0/1")
        vaf = float(rec.get("vaf", 0.50))
        k_ref = int(rec.get("k_ref", 10))
        k_alt = int(rec.get("k_alt", 10))

        format_keys = "GT:GQ:VAF:AD"
        sample_vals = f"{gt}:{gq:.1f}:{vaf:.2f}:{k_ref},{k_alt}"

        line = f"{chrom}\t{pos}\t{mei_id}\t{ref_base}\t{alt_symbol}\t{qual_str}\t{filter_str}\t{info_str}\t{format_keys}\t{sample_vals}"
        lines.append(line)

    out_file.write_text("\n".join(lines) + "\n")
    return out_file
