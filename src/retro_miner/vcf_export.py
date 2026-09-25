"""Export MEI candidate loci to VCF v4.3.

Two inputs are accepted:

* The genome-wide gold review table (``candidate_loci.mei.gold_review.tsv``),
  or the annotated TSV from ``annotate-mei-support``. Both carry
  ``consensus_insertion_breakpoint_pos`` and ``consensus_mei_family``.
* A classifier-ranked table (``chrom``, ``window_start``, ``window_end``,
  ``mei_family``, ``gold_score``, ``classifier_rank``). That table has no
  breakpoint, so pass the gold review table as ``breakpoint_tsv``. Rows are
  joined on ``chrom`` + ``window_start`` + ``window_end``, and the classifier
  score is written to INFO.

``SVLEN`` is intentionally omitted. Ensembl VEP treats ``POS + SVLEN`` as a
reference span for symbolic insertions and ignores ``END``. The insertion
length stays in ``MEI_SPAN``. ``END`` is set equal to ``POS``.

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
import re
from decimal import ROUND_HALF_UP, Context, Decimal, InvalidOperation
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

#: Labels used in gold tables (ALU/LINE1/SVA) and in the classifier table
#: (Alu/L1/SVA). Keys are compared after uppercasing and stripping ``-``/``_``.
_FAMILY_ALIASES = {
    "ALU": "ALU",
    "L1": "LINE1",
    "LINE": "LINE1",
    "LINE1": "LINE1",
    "SVA": "SVA",
}

#: Extra INFO fields copied when the source column is present and non-blank.
#: (INFO id, source column, VCF Number/Type, description)
_EXTRA_INFO_FIELDS: tuple[tuple[str, str, str, str], ...] = (
    (
        "L1TXGOLDSCORE",
        "gold_score",
        "Sig4",
        "Classifier probability that the gold call is a true insertion, 4 significant figures",
    ),
    (
        "L1TXRANKPCT",
        "l1tx_rank_pct",
        "Integer",
        "Percentile of classifier_rank in the input table; 100 is the best call, rounded to the nearest integer",
    ),
    (
        "L1TXGOLDRANKPCT",
        "l1tx_gold_rank_pct",
        "Integer",
        "Percentile of gold_rank when the table has no classifier rank; 100 is the best call, rounded to the nearest integer",
    ),
    (
        "INSERTIONSCORE",
        "insertion_model_score",
        "Float",
        "Pipeline insertion-model heuristic score",
    ),
    (
        "PEAKDEPTHZ",
        "local_bam_peak_depth_z",
        "Float",
        "Local BAM peak-depth z score",
    ),
    (
        "SPLITCLUSTERZ",
        "split_cluster_binomial_z",
        "Float",
        "Split-cluster binomial z score",
    ),
    (
        "CALLTIER",
        "insertion_call_tier",
        "String",
        "Insertion call tier from the gold review table",
    ),
    (
        "KNOWNMEI",
        "known_mei_polymorphism",
        "String",
        "Whether the call overlaps a known MEI polymorphism",
    ),
)

#: ID-only contig declarations (no lengths -- this tool does not carry
#: genome-build metadata, and fabricating lengths would be worse than
#: omitting them). Covers the standard GRCh38-style "chr"-prefixed
#: autosomes/sex/mito contigs that candidate loci are called against.
_STANDARD_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

#: Pipeline ``--reference-build`` value -> VCF contig assembly tag.
_REFERENCE_ASSEMBLY = {
    "hg38": "GRCh38",
    "hg19": "GRCh37",
    "hs1": "T2T-CHM13v2.0",
}

#: Written by the candidate pipeline next to the run directory.
PIPELINE_PARAMS_NAME = "pipeline_params.env"

#: Internal kinds that are not VCF Types. Sig4 is a Float written at 4 significant figures.
_VCF_TYPE = {"Sig4": "Float"}

_KNOWN_ID_SPLIT = re.compile(r"[|;]")
_NSSV_TOKEN = re.compile(r"nssv\d+")

VCF_HEADER_LINES = [
    f"##fileformat={VCF_FILEFORMAT}",
    "##source=retrotransposon-miner:vcf_export",
    *[f"##contig=<ID={c}>" for c in _STANDARD_CONTIGS],
    '##ALT=<ID=INS:ME:ALU,Description="Insertion of an Alu mobile element">',
    '##ALT=<ID=INS:ME:LINE1,Description="Insertion of a LINE-1 mobile element">',
    '##ALT=<ID=INS:ME:SVA,Description="Insertion of an SVA mobile element">',
    '##ALT=<ID=INS:ME,Description="Insertion of a mobile element of unresolved family">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position; equal to POS for a symbolic insertion, which consumes no reference bases">',
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
    '##INFO=<ID=L1TXNSSV,Number=1,Type=String,Description="dbVar nssv accession of an overlapping known MEI">',
    '##INFO=<ID=L1TXG1K,Number=1,Type=String,Description="1000 Genomes MELT identifier of an overlapping known MEI">',
    '##INFO=<ID=L1TXLR,Number=1,Type=String,Description="Long-read SVAN identifier of an overlapping known MEI">',
    '##INFO=<ID=KNOWN_SRC,Number=1,Type=String,Description="Source database for the known-MEI overlap (e.g. melt_1kg, long_read_1kg_ont_vienna)">',
    '##INFO=<ID=SAMPLE_STATUS,Number=1,Type=String,Description="shared / disease_only / control_only support classification">',
    '##INFO=<ID=CTRL_SUPPORT,Number=1,Type=String,Description="Raw control-sample supporting-read evidence string (pipe-delimited key=value pairs)">',
    '##INFO=<ID=DISEASE_SUPPORT,Number=1,Type=String,Description="Raw disease-sample supporting-read evidence string (pipe-delimited key=value pairs)">',
    *[
        f'##INFO=<ID={info_id},Number=1,Type={_VCF_TYPE.get(kind, kind)},Description="{desc}">'
        for info_id, _column, kind, desc in _EXTRA_INFO_FIELDS
    ],
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


def _row_value(row: dict[str, Any], *keys: str) -> str:
    """First non-blank value among ``keys``."""
    for key in keys:
        if key not in row:
            continue
        text = _vcf_safe(row.get(key))
        if text:
            return text
    return ""


def _canonical_family(raw: str) -> str:
    """Map gold (ALU/LINE1/SVA) and classifier (Alu/L1/SVA) labels onto one family."""
    key = raw.strip().upper().replace("-", "").replace("_", "")
    return _FAMILY_ALIASES.get(key, raw.strip())


def _alt_allele_for_family(family: str) -> str:
    alt_id = _MEI_FAMILY_ALT_ID.get(_canonical_family(family).upper(), "INS:ME")
    return f"<{alt_id}>"


def _info_field(key: str, value: str) -> str | None:
    """Return ``KEY=value`` or None if value is empty (omit rather than fabricate)."""
    if value == "":
        return None
    return f"{key}={value}"


def _info_int(key: str, raw: Any) -> str | None:
    text = _vcf_safe(raw)
    if text == "":
        return None
    try:
        return f"{key}={int(float(text))}"
    except ValueError:
        return None


def _info_float(key: str, raw: Any) -> str | None:
    text = _vcf_safe(raw)
    if text == "":
        return None
    try:
        float(text)
    except ValueError:
        return None
    return f"{key}={text}"


def format_sigfigs(raw: Any, n: int = 4) -> str:
    """Round ``raw`` to ``n`` significant figures and keep a fixed-point rendering."""
    text = _vcf_safe(raw)
    if text == "":
        raise ValueError("empty")
    ctx = Context(prec=n, rounding=ROUND_HALF_UP)
    try:
        dec = ctx.create_decimal(text)
    except InvalidOperation as exc:
        raise ValueError(text) from exc
    return format(dec, "f")


def _info_sigfigs(key: str, raw: Any, n: int = 4) -> str | None:
    try:
        return f"{key}={format_sigfigs(raw, n)}"
    except ValueError:
        return None


def split_known_mei_ids(raw: Any) -> dict[str, str]:
    """Split a combined polymorphism id into nssv, g1k, and long-read ids.

    Combined rows look like ``g1k:nssv14058545|lr:chr5-43795330-INS->...``.
    A bare ``nssv...`` value is the 1000 Genomes MELT id. A bare non-nssv
    value is the long-read SVAN id.
    """
    text = "" if raw is None else str(raw).strip()
    if text == "" or text.lower() == "nan":
        return {}
    found: dict[str, str] = {}
    tagged = False
    for part in _KNOWN_ID_SPLIT.split(text):
        piece = part.strip()
        if piece.startswith("g1k:"):
            tagged = True
            found["g1k"] = piece[4:]
        elif piece.startswith("lr:"):
            tagged = True
            found["lr"] = piece[3:]
    if not tagged:
        if text.startswith("nssv"):
            found["g1k"] = text
        else:
            found["lr"] = text
    nssv = _NSSV_TOKEN.search(found.get("g1k", "")) or _NSSV_TOKEN.search(text)
    if nssv is not None:
        found["nssv"] = nssv.group(0)
    return found


def _has_value(rows: list[dict[str, Any]], column: str) -> bool:
    return any(_vcf_safe(row.get(column)) != "" for row in rows)


def _fill_percentile(rows: list[dict[str, Any]], rank_column: str, out_column: str) -> None:
    """Rank 1 is best. Percentile is ``100 * (N - rank + 1) / N``, nearest integer."""
    ranks: list[int] = []
    for row in rows:
        text = _vcf_safe(row.get(rank_column))
        if text == "":
            continue
        try:
            ranks.append(int(float(text)))
        except ValueError:
            continue
    if not ranks:
        return
    n = max(ranks)
    for row in rows:
        text = _vcf_safe(row.get(rank_column))
        if text == "":
            continue
        try:
            rank = int(float(text))
        except ValueError:
            continue
        pct = int(round(100.0 * (n - rank + 1) / n))
        row[out_column] = str(max(0, min(100, pct)))


def attach_rank_percentiles(rows: list[dict[str, Any]]) -> None:
    """Percentile from the classifier rank, or from ``gold_rank`` when that is absent.

    A classifier table also carries ``gold_rank``. That absolute rank is not
    written once ``L1TXRANKPCT`` is available.
    """
    if _has_value(rows, "classifier_rank"):
        if not _has_value(rows, "l1tx_rank_pct"):
            _fill_percentile(rows, "classifier_rank", "l1tx_rank_pct")
        return
    if not _has_value(rows, "l1tx_gold_rank_pct"):
        _fill_percentile(rows, "gold_rank", "l1tx_gold_rank_pct")


def parse_pipeline_params(path: Path) -> dict[str, str]:
    """Read ``key=value`` lines written by the candidate pipeline."""
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if text == "" or text.startswith("#") or "=" not in text:
            continue
        key, value = text.split("=", 1)
        values[key.strip()] = value.strip().strip('"').strip("'")
    return values


def lookup_reference_build(start: Path) -> tuple[str, str]:
    """Return ``(reference_build, reference_fasta)`` from ``pipeline_params.env``.

    The file is the ``--reference-build`` choice that selected the FASTA
    (``hg38`` -> ``Homo_sapiens_assembly38.fasta``). Search starts at ``start``
    and walks up to the run directory.
    """
    current = start.resolve()
    if current.is_file():
        current = current.parent
    for _ in range(6):
        params = current / PIPELINE_PARAMS_NAME
        if params.is_file():
            values = parse_pipeline_params(params)
            return values.get("reference_build", ""), values.get("reference_fasta", "")
        if current.parent == current:
            break
        current = current.parent
    return "", ""


def vcf_header_lines(reference_build: str | None = None) -> list[str]:
    """Header lines, with ``##reference`` and contig assembly when the build is known."""
    build = (reference_build or "").strip()
    if build == "":
        return list(VCF_HEADER_LINES)
    assembly = _REFERENCE_ASSEMBLY.get(build, "")
    lines: list[str] = []
    for line in VCF_HEADER_LINES:
        lines.append(line)
        if line.startswith("##source="):
            lines.append(f"##reference={build}")
        elif assembly and line.startswith("##contig=<ID=") and line.endswith(">"):
            lines[-1] = line[:-1] + f",assembly={assembly}>"
    return lines


def locus_id(chrom: str, pos: str, family: str) -> str:
    """``L1TX-<chrom>-<pos>-<family>`` using the canonical family token."""
    fam = family if family in _MEI_FAMILY_ALT_ID else "ME"
    return f"L1TX-{chrom}-{pos}-{fam}"


def breakpoint_pos(row: dict[str, Any]) -> str:
    """1-based insertion breakpoint, or ``.`` when the row has none."""
    raw = _row_value(row, "consensus_insertion_breakpoint_pos")
    try:
        return str(int(float(raw)))
    except (TypeError, ValueError):
        return "."


def _window_key(row: dict[str, Any]) -> tuple[str, str, str]:
    chrom = str(row.get("chrom", "")).strip()

    def norm(value: Any) -> str:
        text = "" if value is None else str(value).strip()
        try:
            return str(int(float(text)))
        except (TypeError, ValueError):
            return text

    return (chrom, norm(row.get("window_start")), norm(row.get("window_end")))


def fill_breakpoint_columns(
    rows: list[dict[str, Any]],
    locus_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Copy blank columns from the gold table onto classifier rows.

    Join key is ``chrom``, ``window_start``, ``window_end``. Existing
    classifier columns (score, rank, family) are left as they are.
    """
    index: dict[tuple[str, str, str], dict[str, Any]] = {}
    for locus in locus_rows:
        index[_window_key(locus)] = locus

    merged: list[dict[str, Any]] = []
    missed = 0
    for row in rows:
        locus = index.get(_window_key(row))
        if locus is None:
            missed += 1
            merged.append(dict(row))
            continue
        out = dict(row)
        for column, value in locus.items():
            if _vcf_safe(out.get(column)) == "":
                out[column] = value
        merged.append(out)
    if missed:
        raise ValueError(
            f"{missed} of {len(rows)} rows did not match the breakpoint table "
            "on chrom, window_start, window_end"
        )
    return merged


def filter_min_score(
    rows: list[dict[str, Any]],
    min_score: float,
    column: str = "gold_score",
) -> list[dict[str, Any]]:
    """Keep rows whose ``column`` is a number >= ``min_score``."""
    if rows and column not in rows[0]:
        raise ValueError(f"score column {column!r} is not in the input table")
    kept: list[dict[str, Any]] = []
    for row in rows:
        try:
            score = float(row.get(column, ""))
        except (TypeError, ValueError):
            continue
        if score >= min_score:
            kept.append(row)
    return kept


def resolve_sample_name(rows: list[dict[str, Any]], sample_name: str) -> str:
    """Use a uniform ``sample`` column when the caller left the default label."""
    if sample_name != "SAMPLE":
        return sample_name
    values = {_vcf_safe(row.get("sample")) for row in rows}
    values.discard("")
    if len(values) == 1:
        return values.pop()
    return sample_name


def build_vcf_record(row: dict[str, Any], *, mark_pass: bool = False) -> str:
    """Build one VCF data line from a gold-review or classifier row.

    Missing optional columns are omitted from INFO rather than filled with
    a placeholder. ``GT`` and ``GQ`` stay missing. ``SVLEN`` is not written.
    """
    chrom = _vcf_safe(row.get("chrom")) or "."
    pos = breakpoint_pos(row)

    ref = "N"  # No reference-genome lookup performed; POS marks the breakpoint, not a called base.
    family_raw = _row_value(row, "consensus_mei_family", "mei_family")
    family = _canonical_family(family_raw) if family_raw else ""
    alt = _alt_allele_for_family(family) if family else "<INS:ME>"
    vcf_id = locus_id(chrom, pos, family) if chrom != "." and pos != "." else "."
    known = split_known_mei_ids(row.get("known_mei_polymorphism_id"))

    qual = "."
    filt = "PASS" if mark_pass else "."

    info_parts = [
        _info_field("SVTYPE", "INS"),
        _info_int("END", pos) if pos != "." else None,
        _info_int("WINDOWSTART", row.get("window_start")),
        _info_int("WINDOWEND", row.get("window_end")),
        _info_field("MEIFAMILY", _vcf_safe(family)),
        _info_field("MEISUBFAMILY", _row_value(row, "consensus_mei_subfamily", "mei_subfamily")),
        _info_field("TSD", _vcf_safe(row.get("consensus_tsd_seq"))),
        _info_int("POLYA_MIN_BP", row.get("consensus_poly_at_min_bp")),
        _info_field("ORIENT", _vcf_safe(row.get("consensus_insertion_orientation"))),
        _info_field("NESTED", _vcf_safe(row.get("nested_in_same_MEI"))),
        _info_int("MEI_SPAN", row.get("consensus_insertion_mei_span_full")),
        _info_int("MEI_5P", row.get("consensus_insertion_mei_5p_coord_full")),
        _info_int("MEI_3P", row.get("consensus_insertion_mei_3p_coord_full")),
        _info_field("L1TXNSSV", _vcf_safe(known.get("nssv"))),
        _info_field("L1TXG1K", _vcf_safe(known.get("g1k"))),
        _info_field("L1TXLR", _vcf_safe(known.get("lr"))),
        _info_field("KNOWN_SRC", _vcf_safe(row.get("known_mei_polymorphism_source"))),
        _info_field("SAMPLE_STATUS", _vcf_safe(row.get("sample_status_label"))),
        _info_field("CTRL_SUPPORT", _vcf_safe(row.get("control_supporting_reads"))),
        _info_field("DISEASE_SUPPORT", _vcf_safe(row.get("disease_supporting_reads"))),
    ]
    for info_id, column, kind, _desc in _EXTRA_INFO_FIELDS:
        raw = row.get(column)
        if kind == "Integer":
            info_parts.append(_info_int(info_id, raw))
        elif kind == "Float":
            info_parts.append(_info_float(info_id, raw))
        elif kind == "Sig4":
            info_parts.append(_info_sigfigs(info_id, raw, 4))
        else:
            info_parts.append(_info_field(info_id, _vcf_safe(raw)))
    info = ";".join(part for part in info_parts if part is not None)
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
    pos_text = breakpoint_pos(row)
    pos = int(pos_text) if pos_text != "." else 0
    return (ordinal, chrom, pos)


def export_vcf(
    rows: list[dict[str, Any]],
    out_path: Path,
    *,
    sample_name: str = "SAMPLE",
    sort: bool = True,
    mark_pass: bool = False,
    reference_build: str | None = None,
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
    attach_rank_percentiles(rows)
    header = vcf_header_lines(reference_build)
    header.append(VCF_COLUMN_HEADER.format(sample=sample_name))

    ordered = sorted(rows, key=_sort_key) if sort else list(rows)

    with open(out_path, "w", newline="") as fh:
        for line in header:
            fh.write(line + "\n")
        n = 0
        for row in ordered:
            fh.write(build_vcf_record(row, mark_pass=mark_pass) + "\n")
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
    min_score: float | None = None,
    score_column: str = "gold_score",
    breakpoint_tsv: Path | None = None,
    reference_build: str | None = None,
) -> int:
    """Read a gold-review or classifier TSV and write VCF.

    ``breakpoint_tsv`` is the genome-wide gold table used to fill the
    insertion breakpoint when ``tsv_path`` is a classifier ranking. When
    ``min_score`` is set, only rows with ``score_column >= min_score`` are
    written and FILTER is PASS.
    """
    rows = load_tsv_rows(tsv_path)
    if breakpoint_tsv is not None:
        rows = fill_breakpoint_columns(rows, load_tsv_rows(breakpoint_tsv))
    elif rows and not any(breakpoint_pos(row) != "." for row in rows):
        raise ValueError(
            "input has no consensus_insertion_breakpoint_pos; pass breakpoint_tsv "
            "pointing at the genome-wide gold review table"
        )
    # Percentile uses the full input, including rows the score filter will drop.
    attach_rank_percentiles(rows)
    if min_score is not None:
        rows = filter_min_score(rows, min_score, score_column)
    sample = resolve_sample_name(rows, sample_name)
    build = (reference_build or "").strip()
    if build == "":
        build, _fasta = lookup_reference_build(tsv_path)
    return export_vcf(
        rows,
        out_path,
        sample_name=sample,
        sort=sort,
        mark_pass=min_score is not None,
        reference_build=build or None,
    )
