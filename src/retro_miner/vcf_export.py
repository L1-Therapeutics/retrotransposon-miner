"""Export annotated MEI candidate loci to VCF v4.3.

Input is the annotated TSV produced by ``annotate-mei-support``
(``--out-tsv``); its columns are a superset of
:data:`retro_miner.readme_example_table.TABLE_COLS`, which is the
canonical, tested column set for a finished locus record.

Genotype fields are intentionally left blank (``GT=./.``, ``GQ=.``).
This table reports pooled disease-vs-control read support across a
tumor/normal (or disease/control) cohort comparison -- there is no
per-individual diploid sample here, so there is nothing to genotype in
the classical sense. A read-count-ratio genotyper was attempted on an
earlier branch (PR#39) and rejected in review: it compared two
structurally different read populations (split/discordant "alt"
evidence in a 200bp window vs. proper-pair "ref" evidence in a +/-50bp
window) as if they were directly comparable at a 50/50 het ratio, which
they are not, and the model additionally ignored ploidy variation (chrX),
somatic context, and mosaicism. Rather than carry that logic forward,
GT/GQ are left as VCF-standard missing values so the file is still
valid and consumable by downstream tools (bcftools, IGV, vcf-validator),
and genotyping can be added later with a model that's actually correct
for MEI evidence (e.g. a pangenome/graph-based genotyper).
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

VCF_FILEFORMAT = "VCFv4.3"

#: Reserved VCF INFO/text characters that must not appear raw in a value.
#: We substitute with printable, unambiguous replacements rather than
#: percent-encoding, since these values are also meant to be human-readable
#: when a VCF is opened directly.
_VCF_UNSAFE_CHARS = {
    ";": "_",
    "=": ":",
    ",": "|",
    " ": "_",
    "\t": "_",
    "\n": "_",
}

#: Known MEI family -> symbolic ALT allele ID, following the 1000 Genomes
#: MEI VCF convention (e.g. MELT, PMC7017486). Unrecognized families fall
#: back to the generic ``<INS:ME>`` symbolic allele.
_MEI_FAMILY_ALT_ID = {
    "ALU": "INS:ME:ALU",
    "LINE1": "INS:ME:LINE1",
    "SVA": "INS:ME:SVA",
}

#: ID-only contig declarations (no lengths -- this tool does not carry
#: genome-build metadata, and fabricating lengths would be worse than
#: omitting them). Covers the standard GRCh38-style "chr"-prefixed
#: autosomes/sex/mito contigs that candidate loci are called against.
_STANDARD_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

VCF_HEADER_LINES = [
    f"##fileformat={VCF_FILEFORMAT}",
    "##source=retrotransposon-miner:vcf_export",
    *[f"##contig=<ID={c}>" for c in _STANDARD_CONTIGS],
    '##ALT=<ID=INS:ME:ALU,Description="Insertion of an Alu mobile element">',
    '##ALT=<ID=INS:ME:LINE1,Description="Insertion of a LINE-1 mobile element">',
    '##ALT=<ID=INS:ME:SVA,Description="Insertion of an SVA mobile element">',
    '##ALT=<ID=INS:ME,Description="Insertion of a mobile element of unresolved family">',
    '##INFO=<ID=WINDOWSTART,Number=1,Type=Integer,Description="Candidate locus window start">',
    '##INFO=<ID=WINDOWEND,Number=1,Type=Integer,Description="Candidate locus window end">',
    '##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="Consensus MEI family (e.g. ALU, LINE1, SVA)">',
    '##INFO=<ID=MEISUBFAMILY,Number=1,Type=String,Description="Consensus MEI subfamily call">',
    '##INFO=<ID=TSD,Number=1,Type=String,Description="Consensus target site duplication sequence">',
    '##INFO=<ID=POLYA_MIN_BP,Number=1,Type=Integer,Description="Minimum observed poly-A/T tail length in bp">',
    '##INFO=<ID=ORIENT,Number=1,Type=String,Description="Consensus insertion orientation relative to reference (+/-)">',
    '##INFO=<ID=NESTED,Number=1,Type=String,Description="Whether the insertion is nested within another copy of the same MEI">',
    '##INFO=<ID=MEI_SPAN,Number=1,Type=Integer,Description="Full-length span (bp) of the consensus MEI alignment">',
    '##INFO=<ID=MEI_5P,Number=1,Type=Integer,Description="5-prime coordinate of the consensus MEI alignment on the full-length reference">',
    '##INFO=<ID=MEI_3P,Number=1,Type=Integer,Description="3-prime coordinate of the consensus MEI alignment on the full-length reference">',
    '##INFO=<ID=KNOWN_ID,Number=1,Type=String,Description="Matching known MEI polymorphism identifier, if any">',
    '##INFO=<ID=KNOWN_SRC,Number=1,Type=String,Description="Source database for KNOWN_ID (e.g. melt_1kg, long_read_1kg_ont_vienna)">',
    '##INFO=<ID=SAMPLE_STATUS,Number=1,Type=String,Description="shared / disease_only / control_only support classification">',
    '##INFO=<ID=CTRL_SUPPORT,Number=1,Type=String,Description="Raw control-sample supporting-read evidence string (pipe-delimited key=value pairs)">',
    '##INFO=<ID=DISEASE_SUPPORT,Number=1,Type=String,Description="Raw disease-sample supporting-read evidence string (pipe-delimited key=value pairs)">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (unset: this table reports pooled disease-vs-control read support, not per-individual diploid genotypes)">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality (unset: no validated genotyping model for MEI evidence yet -- see module docstring)">',
]

VCF_COLUMN_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"


def _vcf_safe(value: Any) -> str:
    """Sanitize a value for use inside a single VCF INFO/text field."""
    text = "" if value is None else str(value)
    text = text.strip()
    if text == "" or text.lower() == "nan":
        return ""
    for bad, good in _VCF_UNSAFE_CHARS.items():
        text = text.replace(bad, good)
    return text


def _alt_allele_for_family(family: str) -> str:
    alt_id = _MEI_FAMILY_ALT_ID.get(family.strip().upper(), "INS:ME")
    return f"<{alt_id}>"


def _info_field(key: str, value: str) -> str | None:
    """Return ``KEY=value`` or None if value is empty (omit rather than fabricate)."""
    if value == "":
        return None
    return f"{key}={value}"


def build_vcf_record(row: dict[str, Any]) -> str:
    """Build one VCF data line from a single annotated candidate-loci row.

    ``row`` is expected to have (at least) the columns in
    :data:`retro_miner.readme_example_table.TABLE_COLS`. Missing/blank
    optional columns are simply omitted from INFO rather than filled with
    a placeholder value.
    """
    chrom = _vcf_safe(row.get("chrom")) or "."
    pos_raw = row.get("consensus_insertion_breakpoint_pos")
    try:
        pos = str(int(float(pos_raw)))
    except (TypeError, ValueError):
        pos = "."

    known_id = _vcf_safe(row.get("known_mei_polymorphism_id"))
    vcf_id = known_id if known_id else "."

    ref = "N"  # No reference-genome lookup performed; POS marks the breakpoint, not a called base.
    family = _vcf_safe(row.get("consensus_mei_family"))
    alt = _alt_allele_for_family(family) if family else "<INS:ME>"

    qual = "."
    filt = "."  # Not evaluated -- no validated pass/fail model without genotyping.

    info_parts = [
        _info_field("WINDOWSTART", _vcf_safe(row.get("window_start"))),
        _info_field("WINDOWEND", _vcf_safe(row.get("window_end"))),
        _info_field("MEIFAMILY", family),
        _info_field("MEISUBFAMILY", _vcf_safe(row.get("consensus_mei_subfamily"))),
        _info_field("TSD", _vcf_safe(row.get("consensus_tsd_seq"))),
        _info_field("POLYA_MIN_BP", _vcf_safe(row.get("consensus_poly_at_min_bp"))),
        _info_field("ORIENT", _vcf_safe(row.get("consensus_insertion_orientation"))),
        _info_field("NESTED", _vcf_safe(row.get("nested_in_same_MEI"))),
        _info_field("MEI_SPAN", _vcf_safe(row.get("consensus_insertion_mei_span_full"))),
        _info_field("MEI_5P", _vcf_safe(row.get("consensus_insertion_mei_5p_coord_full"))),
        _info_field("MEI_3P", _vcf_safe(row.get("consensus_insertion_mei_3p_coord_full"))),
        _info_field("KNOWN_ID", known_id),
        _info_field("KNOWN_SRC", _vcf_safe(row.get("known_mei_polymorphism_source"))),
        _info_field("SAMPLE_STATUS", _vcf_safe(row.get("sample_status_label"))),
        _info_field("CTRL_SUPPORT", _vcf_safe(row.get("control_supporting_reads"))),
        _info_field("DISEASE_SUPPORT", _vcf_safe(row.get("disease_supporting_reads"))),
    ]
    info = ";".join(p for p in info_parts if p is not None)
    if info == "":
        info = "."

    fmt = "GT:GQ"
    sample_value = "./.:."  # Intentionally blank -- see module docstring.

    return "\t".join([chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, sample_value])


#: Contig -> ordinal, matching the order contigs are declared in the header.
#: VCF requires records be sorted by contig (in header-declaration order)
#: then by position; otherwise ``bcftools index`` refuses the file and most
#: downstream annotation tools cannot consume it.
_CONTIG_ORDER = {name: i for i, name in enumerate(_STANDARD_CONTIGS)}


def _sort_key(row: dict[str, Any]) -> tuple[int, str, int]:
    """Coordinate sort key: (contig ordinal, contig name, position).

    Contigs not in the standard set sort after all standard ones, then
    alphabetically by name, so an unexpected contig (e.g. a decoy or alt
    scaffold) is still emitted deterministically rather than dropped.
    """
    chrom = str(row.get("chrom", "")).strip()
    ordinal = _CONTIG_ORDER.get(chrom, len(_CONTIG_ORDER))
    try:
        pos = int(float(row.get("consensus_insertion_breakpoint_pos")))
    except (TypeError, ValueError):
        pos = 0
    return (ordinal, chrom, pos)


def export_vcf(
    rows: list[dict[str, Any]],
    out_path: Path,
    *,
    sample_name: str = "SAMPLE",
    sort: bool = True,
) -> int:
    """Write ``rows`` (annotated candidate loci) to ``out_path`` as VCF v4.3.

    Records are coordinate-sorted by default. The upstream candidate-loci
    table is sorted by ``enrichment_ratio`` (most interesting first), which
    is useful for human review but is not valid VCF ordering -- ``bcftools
    index`` rejects it with "Unsorted positions", and an unindexable VCF
    cannot be fed to most downstream annotation tools. Pass ``sort=False``
    only if the caller has already guaranteed coordinate order.

    Returns the number of records written.
    """
    header = list(VCF_HEADER_LINES)
    header.append(VCF_COLUMN_HEADER.format(sample=sample_name))

    ordered = sorted(rows, key=_sort_key) if sort else list(rows)

    with open(out_path, "w", newline="") as fh:
        for line in header:
            fh.write(line + "\n")
        n = 0
        for row in ordered:
            fh.write(build_vcf_record(row) + "\n")
            n += 1
    return n


def load_tsv_rows(tsv_path: Path) -> list[dict[str, Any]]:
    """Load an annotated candidate-loci TSV into a list of row dicts."""
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def export_vcf_from_tsv(
    tsv_path: Path,
    out_path: Path,
    *,
    sample_name: str = "SAMPLE",
    sort: bool = True,
) -> int:
    """Convenience wrapper: read an annotated TSV, write it out as VCF."""
    rows = load_tsv_rows(tsv_path)
    return export_vcf(rows, out_path, sample_name=sample_name, sort=sort)
