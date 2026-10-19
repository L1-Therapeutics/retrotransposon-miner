"""Gene / functional-consequence annotation of MEI calls via local snpEff.

Consumes a VCF of mobile-element insertions (symbolic ``<INS:ME:*>`` ALT
alleles, as written by ``vcf_export``) and appends gene-level and
transcript-aware splicing consequence fields to each record's INFO column.

``CONSEQUENCE_SEVERITY`` reproduces Ensembl's "Calculated variant
consequences" table
(https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html),
most severe first. snpEff does not order ``ANN`` entries that way, so
:func:`parse_snpeff_ann` applies the table itself, prioritizing MANE Select /
canonical and protein-coding transcripts over retained-intron or non-coding
isoforms. Terms are Sequence Ontology identifiers.

References
----------
Cingolani P. et al. (2012) A program for annotating and predicting the
    effects of single nucleotide polymorphisms, SnpEff. Fly 6(2):80-92.
    doi:10.4161/fly.19695
Eilbeck K. et al. (2005) The Sequence Ontology: a tool for the unification
    of genome annotations. Genome Biology 6:R44. doi:10.1186/gb-2005-6-5-r44
"""

from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

# Ensembl "Calculated variant consequences" table, most severe first.
CONSEQUENCE_SEVERITY: tuple[str, ...] = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
)
_SEVERITY_RANK = {term: i for i, term in enumerate(CONSEQUENCE_SEVERITY)}

# INFO keys this module writes, with their VCF header definitions.
ANNOTATION_INFO_HEADERS: tuple[str, ...] = (
    '##INFO=<ID=GENE,Number=.,Type=String,Description="Gene symbols overlapped by the insertion breakpoint (snpEff)">',
    '##INFO=<ID=GENEID,Number=.,Type=String,Description="Ensembl gene IDs overlapped by the insertion breakpoint">',
    '##INFO=<ID=GENE_STRAND,Number=.,Type=String,Description="Transcribed strand (+ or -) of each GENEID, in the same order. An insertion is sense when ORIENT matches this strand">',
    '##INFO=<ID=MEI_GENE_ORIENT,Number=.,Type=String,Description="Orientation of the insertion (ORIENT) relative to each GENEID transcribed strand, in the same order: sense, antisense, or . when ORIENT or the strand is unknown">',
    '##INFO=<ID=CSQ,Number=1,Type=String,Description="Most severe Sequence Ontology consequence term across overlapping transcripts (snpEff)">',
    '##INFO=<ID=CSQ_TERMS,Number=.,Type=String,Description="All distinct consequence terms across overlapping transcripts">',
    '##INFO=<ID=CSQ_NTX,Number=1,Type=Integer,Description="Number of transcripts overlapping the insertion breakpoint">',
    '##INFO=<ID=SPLICE_DIST,Number=1,Type=Integer,Description="Distance in bp to the nearest annotated splice donor or acceptor site">',
    '##INFO=<ID=DIST_TO_SPLICE_DONOR,Number=1,Type=Integer,Description="Distance in bp to the upstream splice donor junction in the primary overlapping intron">',
    '##INFO=<ID=DIST_TO_SPLICE_ACCEPTOR,Number=1,Type=Integer,Description="Distance in bp to the downstream splice acceptor junction in the primary overlapping intron">',
    '##INFO=<ID=SPLICE_SIDE,Number=1,Type=String,Description="Nearest splice junction type: donor, acceptor, or exonic">',
    '##INFO=<ID=SPLICE_REGION,Number=0,Type=Flag,Description="Insertion breakpoint falls within canonical splice site or region (<=8bp)">',
    '##INFO=<ID=MEI_EXONIZATION_RISK,Number=1,Type=String,Description="Risk of alternative pseudo-exon inclusion or transcript truncation (HIGH, MODERATE, NONE)">',
    '##INFO=<ID=ALU_EXONIZATION_RISK,Number=1,Type=String,Description="Antisense Alu cryptic 3-prime splice acceptor exonization risk (HIGH, MODERATE, NONE)">',
    '##INFO=<ID=L1_TRUNCATION_RISK,Number=1,Type=String,Description="Sense LINE-1 premature polyadenylation or internal donor splicing truncation risk (HIGH, MODERATE, NONE)">',
    '##INFO=<ID=SVA_EXON_TRAPPING_RISK,Number=1,Type=String,Description="SVA cryptic splice site activation and exon trapping risk in sense or antisense orientation (HIGH, MODERATE, NONE)">',
    '##INFO=<ID=CRYPTIC_DONOR_COUNT,Number=1,Type=Integer,Description="Count of strong cryptic 5-prime splice donor (9-mer GT) motifs on the transcribed pre-mRNA strand of the assembled insertion sequence">',
    '##INFO=<ID=CRYPTIC_ACCEPTOR_COUNT,Number=1,Type=Integer,Description="Count of strong cryptic 3-prime splice acceptor (polypyrimidine + AG 23-mer) motifs on the transcribed pre-mRNA strand of the assembled insertion sequence">',
    '##INFO=<ID=CRYPTIC_PSEUDOEXON_PAIRS,Number=1,Type=Integer,Description="Count of candidate internal pseudo-exon pairs (3-prime acceptor followed 40-350 bp downstream by a 5-prime donor) in the transcribed insertion sequence">',
    '##INFO=<ID=CRYPTIC_SPLICE_MAX_SCORE,Number=1,Type=Float,Description="Maximum normalized consensus splice motif score (0.00-1.00) across scanned cryptic 5-prime donor and 3-prime acceptor sites">',
    '##INFO=<ID=TWIN_PRIMING_INV,Number=0,Type=Flag,Description="Insertion exhibits a twin-priming 5-prime inversion junction (5-prime inverted segment joined to 3-prime non-inverted segment)">',
    '##INFO=<ID=INV_JUNCTION_SPLICE_COUNT,Number=1,Type=Integer,Description="Count of chimeric cryptic 5-prime donor and 3-prime acceptor splice motifs spanning or adjacent to the 5-prime inversion junction">',
    '##INFO=<ID=INV_JUNCTION_SPLICE_MAX_SCORE,Number=1,Type=Float,Description="Maximum normalized consensus splice motif score (0.00-1.00) across the twin-priming 5-prime inversion junction window">',
    '##INFO=<ID=TX_ID,Number=.,Type=String,Description="Transcript IDs impacted by the insertion breakpoint">',
    '##INFO=<ID=TX_BIOTYPE,Number=.,Type=String,Description="Biotypes of affected transcripts">',
    '##INFO=<ID=TX_RANK,Number=.,Type=String,Description="Exon or intron rank (e.g. 3/10) in the primary affected transcript">',
    '##INFO=<ID=EXON_RANK,Number=1,Type=String,Description="Exon rank (e.g. 3/10) when the primary affected transcript consequence is exonic">',
    '##INFO=<ID=INTRON_RANK,Number=1,Type=String,Description="Intron rank (e.g. 3/10) when the primary affected transcript consequence is intronic">',
)
ANNOTATION_INFO_KEYS: tuple[str, ...] = (
    "GENE",
    "GENEID",
    "GENE_STRAND",
    "MEI_GENE_ORIENT",
    "CSQ",
    "CSQ_TERMS",
    "CSQ_NTX",
    "SPLICE_DIST",
    "DIST_TO_SPLICE_DONOR",
    "DIST_TO_SPLICE_ACCEPTOR",
    "SPLICE_SIDE",
    "SPLICE_REGION",
    "MEI_EXONIZATION_RISK",
    "ALU_EXONIZATION_RISK",
    "L1_TRUNCATION_RISK",
    "SVA_EXON_TRAPPING_RISK",
    "CRYPTIC_DONOR_COUNT",
    "CRYPTIC_ACCEPTOR_COUNT",
    "CRYPTIC_PSEUDOEXON_PAIRS",
    "CRYPTIC_SPLICE_MAX_SCORE",
    "TWIN_PRIMING_INV",
    "INV_JUNCTION_SPLICE_COUNT",
    "INV_JUNCTION_SPLICE_MAX_SCORE",
    "TX_ID",
    "TX_BIOTYPE",
    "TX_RANK",
    "EXON_RANK",
    "INTRON_RANK",
)

# Standard Sequence Ontology splice-site interval parameters for snpEff:
#   - Core splice donor/acceptor dinucleotide: 2 bp (-spliceSiteSize 2)
#   - Exonic splice region edge: 3 bp into exon (-spliceRegionExonSize 3)
#   - Intronic splice region edge: 3-8 bp into intron (-spliceRegionIntronMin 3, -spliceRegionIntronMax 8)
SNPEFF_SPLICE_ARGS: tuple[str, ...] = (
    "-spliceSiteSize", "2",
    "-spliceRegionExonSize", "3",
    "-spliceRegionIntronMin", "3",
    "-spliceRegionIntronMax", "8",
)

SPLICE_REGION_MAX_BP = 8
DEEP_INTRONIC_EXONIZATION_BP = 500

_INTRONIC_OR_SPLICE_TERMS: frozenset[str] = frozenset({
    "intron_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
})

_EXONIC_TERMS: frozenset[str] = frozenset({
    "coding_sequence_variant",
    "inframe_insertion",
    "inframe_deletion",
    "frameshift_variant",
    "missense_variant",
    "protein_altering_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "incomplete_terminal_codon_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
})

_UP_DOWN_INTERGENIC_TERMS: frozenset[str] = frozenset({
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_variant",
})


@dataclass
class VcfRecord:
    """One data line of a VCF, split into fields; INFO parsed to an ordered dict."""

    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: dict[str, str | None]
    rest: list[str] = field(default_factory=list)  # FORMAT + sample columns, verbatim

    def to_line(self) -> str:
        parts: list[str] = []
        for k, v in self.info.items():
            if v is None:
                parts.append(k)
            elif v != "" and v != ".":
                parts.append(f"{k}={v}")
        info = ";".join(parts) or "."
        cols = [self.chrom, str(self.pos), self.id, self.ref, self.alt, self.qual, self.filter, info, *self.rest]
        return "\t".join(cols)


@dataclass(frozen=True)
class TranscriptAnn:
    """Parsed per-transcript entry from a snpEff ``ANN`` subfield."""

    allele: str
    terms: tuple[str, ...]
    most_severe: str
    impact: str
    gene_name: str
    gene_id: str
    feature_type: str
    transcript_id: str
    biotype: str
    rank: str


@dataclass(frozen=True)
class ExonRecord:
    """Single annotated GTF exon feature (1-based inclusive coordinates)."""

    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str
    exon_number: int = 0
    biotype: str = ""
    is_mane_or_canonical: bool = False


@dataclass(frozen=True)
class SpliceMetrics:
    """Distance and junction context from the native GTF splice-site engine."""

    splice_dist: int | None
    splice_side: str | None  # "donor" | "acceptor" | "exonic"
    splice_region: bool
    inferred_csq: str | None
    transcript_id: str = ""
    gene_id: str = ""
    rank: str = ""
    dist_to_donor: int | None = None
    dist_to_acceptor: int | None = None


@dataclass(frozen=True)
class GeneAnnotation:
    """Gene- and transcript-level summary collapsed from a snpEff ANN field."""

    gene_symbols: tuple[str, ...]
    gene_ids: tuple[str, ...]
    most_severe: str
    terms: tuple[str, ...]
    n_transcripts: int
    transcript_ids: tuple[str, ...] = ()
    transcript_biotypes: tuple[str, ...] = ()
    primary_rank: str = ""
    exon_rank: str = ""
    intron_rank: str = ""
    transcripts: tuple[TranscriptAnn, ...] = ()

    def info_fields(
        self,
        gene_strands: Mapping[str, str] | None = None,
        orient: str | None = None,
    ) -> dict[str, str]:
        out: dict[str, str] = {}
        if self.gene_symbols:
            out["GENE"] = ",".join(self.gene_symbols)
        if self.gene_ids:
            out["GENEID"] = ",".join(self.gene_ids)
            strands = [
                gene_strands.get(_gene_id_key(gid)) if gene_strands is not None else None
                for gid in self.gene_ids
            ]
            if gene_strands is not None:
                out["GENE_STRAND"] = ",".join(s or "." for s in strands)
            # One entry per gene ID, so this stays index-aligned with GENEID even
            # when a strand or the insertion orientation is unknown.
            out["MEI_GENE_ORIENT"] = ",".join(
                _gene_orientation(orient, s) for s in strands
            )
        out["CSQ"] = self.most_severe
        if self.terms:
            out["CSQ_TERMS"] = ",".join(self.terms)
        out["CSQ_NTX"] = str(self.n_transcripts)
        if self.transcript_ids:
            out["TX_ID"] = ",".join(self.transcript_ids)
        if self.transcript_biotypes:
            out["TX_BIOTYPE"] = ",".join(self.transcript_biotypes)
        if self.primary_rank:
            out["TX_RANK"] = self.primary_rank
        if self.exon_rank:
            out["EXON_RANK"] = self.exon_rank
        if self.intron_rank:
            out["INTRON_RANK"] = self.intron_rank
        return out


# --------------------------------------------------------------------------- VCF I/O


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "rt")


def parse_info(info: str) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    if info in ("", "."):
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = None  # flag
    return out


def read_vcf(path: str | Path) -> tuple[list[str], list[VcfRecord]]:
    """Return (header_lines, records). Header lines keep their trailing '#CHROM' line last."""
    headers: list[str] = []
    records: list[VcfRecord] = []
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            f = line.split("\t")
            if len(f) < 8:
                raise ValueError(f"VCF data line has {len(f)} columns (<8): {line[:80]!r}")
            records.append(
                VcfRecord(
                    chrom=f[0], pos=int(f[1]), id=f[2], ref=f[3], alt=f[4],
                    qual=f[5], filter=f[6], info=parse_info(f[7]), rest=f[8:],
                )
            )
    if not headers or not headers[-1].startswith("#CHROM"):
        raise ValueError(f"{path}: missing #CHROM header line")
    return headers, records


def write_vcf(path: str | Path, headers: Sequence[str], records: Iterable[VcfRecord]) -> None:
    with open(path, "w") as fh:
        for h in headers:
            fh.write(h + "\n")
        for r in records:
            fh.write(r.to_line() + "\n")


def add_info_headers(headers: Sequence[str]) -> list[str]:
    """Insert this module's ##INFO lines before #CHROM, skipping any already present."""
    present = {h.split("ID=", 1)[1].split(",", 1)[0] for h in headers if h.startswith("##INFO=<ID=")}
    new = [h for h in ANNOTATION_INFO_HEADERS if h.split("ID=", 1)[1].split(",", 1)[0] not in present]
    return [*headers[:-1], *new, headers[-1]]


# --------------------------------------------------------------------------- consequence ranking


def most_severe_term(terms: Iterable[str]) -> str:
    """Pick the most severe SO term per Ensembl's ranking; unknown terms rank last."""
    terms = list(terms)
    if not terms:
        return "intergenic_variant"
    return min(terms, key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))


def _record_key(rec: VcfRecord) -> tuple[str, int, str]:
    return (rec.chrom, rec.pos, rec.alt)


@dataclass
class AnnotationStats:
    n_records: int
    n_annotated: int
    n_with_gene: int
    n_unmatched: int
    seconds: float


# --------------------------------------------------------------------------- snpEff backend
#
# snpEff (Cingolani et al. 2012) runs locally against a downloaded Ensembl or
# RefSeq database and writes an ``ANN`` INFO field (one comma-separated entry
# per transcript, 16 pipe-separated columns).
# Default JVM heap OOMs on GRCh38 -> pass -Xmx8g.
# ``-noHgvs`` is always set: HGVS protein notation is not part of this module's
# output, and computing it makes snpEff load the reference sequence into heap.
# Standard Sequence Ontology splice intervals (-spliceSiteSize 2,
# -spliceRegionExonSize 3, -spliceRegionIntronMin 3, -spliceRegionIntronMax 8)
# are enabled so insertions at splice donors, acceptors, and splice regions are
# accurately classified.
#
# ANN entries are NOT ordered by Ensembl severity (snpEff lists e.g.
# upstream_gene_variant ahead of intron_variant for the same gene), so the
# parser below applies CONSEQUENCE_SEVERITY and transcript biotype/MANE
# prioritization itself. ``intergenic_region`` is snpEff's name for SO
# ``intergenic_variant``.

SNPEFF_TERM_ALIASES: dict[str, str] = {
    "intergenic_region": "intergenic_variant",
    "gene_variant": "coding_transcript_variant",
    "transcript_variant": "coding_transcript_variant",
}

ANN_ALLELE = 0
ANN_ANNOTATION = 1
ANN_IMPACT = 2
ANN_GENE_NAME = 3
ANN_GENE_ID = 4
ANN_FEATURE_TYPE = 5
ANN_FEATURE_ID = 6
ANN_BIOTYPE = 7
ANN_RANK = 8

# snpEff ANN has no strand column. This table is Ensembl GRCh38.115 gene
# features (gene_id, transcribed strand), bundled so annotate-genes writes
# GENE_STRAND without a second database load.
_BUNDLED_GENE_STRANDS = Path(__file__).resolve().parent / "data" / "grch38.115.gene_strand.tsv.gz"


def _gene_id_key(gene_id: str) -> str:
    """Drop an Ensembl version suffix so ``ENSG00000093072.20`` matches the table."""
    return gene_id.split(".", 1)[0]


def _tx_id_key(tx_id: str) -> str:
    """Drop an Ensembl/RefSeq version suffix so ``ENST00000399837.8`` matches GTF keys."""
    return tx_id.split(".", 1)[0]


def _gene_orientation(orient: str | None, strand: str | None) -> str:
    """Return ``sense``/``antisense`` for an insertion on ``orient`` overlapping a gene on ``strand``.

    ``"."`` is returned when either side is unknown (missing, or not ``+``/``-``) so
    ``MEI_GENE_ORIENT`` stays index-aligned with ``GENEID`` rather than dropping a position.
    For twin-primed 5'-inverted calls whose ``ORIENT`` encodes both segments (e.g. ``+-`` or
    ``-+``), the 3' poly(A)-bearing segment (last ``+``/``-`` character) defines the primary
    insertion orientation relative to the reference.
    """
    if not orient or strand not in ("+", "-"):
        return "."
    eff_orient = orient.strip()
    if eff_orient not in ("+", "-"):
        signs = [ch for ch in eff_orient if ch in ("+", "-")]
        eff_orient = signs[-1] if signs else ""
    if eff_orient not in ("+", "-"):
        return "."
    return "sense" if eff_orient == strand else "antisense"


_GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')


def _parse_gtf_attributes(attr_text: str) -> tuple[dict[str, str], set[str]]:
    """Parse GTF attribute column into a dict plus the set of ``tag`` values."""
    attrs: dict[str, str] = {}
    tags: set[str] = set()
    for match in _GTF_ATTR_RE.finditer(attr_text):
        key, val = match.group(1), match.group(2)
        if key == "tag":
            tags.add(val)
        elif key not in attrs:
            attrs[key] = val
    return attrs, tags


def load_gene_strands(path: str | Path) -> dict[str, str]:
    """Read gene_id -> ``+``/``-`` from a two-column TSV or an Ensembl/GENCODE GTF.

    GTF input uses ``gene`` features (falling back to ``exon``/``transcript``
    features if a minimal GTF omits ``gene`` lines). Version suffixes on gene IDs
    are stripped.
    """
    path = Path(path)
    opener = gzip.open if path.name.endswith(".gz") else open
    strands: dict[str, str] = {}
    fallback_strands: dict[str, str] = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) == 2 and cols[1] in {"+", "-"}:
                strands[_gene_id_key(cols[0])] = cols[1]
                continue
            if len(cols) < 9 or cols[6] not in {"+", "-"}:
                continue
            gid = ""
            for attr in cols[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id "):
                    gid = attr.split(" ", 1)[1].strip().strip('"')
                    break
            if not gid:
                continue
            if cols[2] == "gene":
                strands[_gene_id_key(gid)] = cols[6]
            elif _gene_id_key(gid) not in fallback_strands:
                fallback_strands[_gene_id_key(gid)] = cols[6]
    if not strands and fallback_strands:
        strands = fallback_strands
    if not strands:
        raise ValueError(f"no gene strands in {path}")
    return strands


@lru_cache(maxsize=1)
def bundled_gene_strands() -> dict[str, str]:
    """Ensembl GRCh38.115 transcribed strand for each gene ID."""
    if not _BUNDLED_GENE_STRANDS.is_file():
        raise FileNotFoundError(f"bundled gene-strand table missing: {_BUNDLED_GENE_STRANDS}")
    return load_gene_strands(_BUNDLED_GENE_STRANDS)


def _resolve_gene_strands(
    gene_strands: Mapping[str, str] | None,
    gene_gtf: str | Path | None,
) -> Mapping[str, str]:
    if gene_gtf is not None:
        return load_gene_strands(gene_gtf)
    if gene_strands is not None:
        return gene_strands
    return bundled_gene_strands()


# --------------------------------------------------------------------------- GTF exon & splice-site engine


@dataclass(frozen=True)
class TranscriptModel:
    """Ordered exons and metadata for a single transcript from a GTF."""

    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str
    is_mane_or_canonical: bool
    exons: tuple[ExonRecord, ...]  # sorted in ascending genomic coordinate order
    tx_start: int
    tx_end: int


@dataclass(frozen=True)
class GtfSpliceIndex:
    """Cached GTF transcript/exon structures and MANE/canonical sets."""

    by_chrom: Mapping[str, tuple[TranscriptModel, ...]]
    mane_transcripts: frozenset[str]
    transcript_biotypes: Mapping[str, str]
    gene_strands: Mapping[str, str]


_MANE_TAGS: frozenset[str] = frozenset({
    "MANE_Select",
    "MANE_Plus_Clinical",
    "Ensembl_canonical",
    "canonical",
})


@lru_cache(maxsize=8)
def load_gtf_exons(path: str | Path) -> tuple[ExonRecord, ...]:
    """Parse and cache all ``exon`` features from a GTF file."""
    p = Path(path)
    opener = gzip.open if p.name.endswith(".gz") else open
    exons: list[ExonRecord] = []
    with opener(p, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs, tags = _parse_gtf_attributes(cols[8])
            tx_id = attrs.get("transcript_id", "")
            gene_id = attrs.get("gene_id", "")
            if not tx_id:
                continue
            exon_num_raw = attrs.get("exon_number", "0")
            try:
                exon_num = int(exon_num_raw)
            except ValueError:
                exon_num = 0
            biotype = (
                attrs.get("transcript_type")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or ""
            )
            is_mane = bool(tags & _MANE_TAGS) or attrs.get("mane_select", "").lower() == "true"
            exons.append(
                ExonRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    transcript_id=tx_id,
                    gene_id=gene_id,
                    exon_number=exon_num,
                    biotype=biotype,
                    is_mane_or_canonical=is_mane,
                )
            )
    return tuple(exons)


@lru_cache(maxsize=8)
def build_gtf_splice_index(path: str | Path) -> GtfSpliceIndex:
    """Build per-chromosome transcript models and MANE/canonical lookup sets from a GTF."""
    p = Path(path)
    opener = gzip.open if p.name.endswith(".gz") else open
    tx_exons: dict[str, list[ExonRecord]] = {}
    tx_meta: dict[str, tuple[str, str, str, str, str, bool]] = {}
    gene_strands: dict[str, str] = {}

    with opener(p, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feat = cols[2]
            if feat not in ("gene", "transcript", "exon"):
                continue
            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs, tags = _parse_gtf_attributes(cols[8])
            gid = attrs.get("gene_id", "")
            gname = attrs.get("gene_name", "") or _gene_id_key(gid)
            if gid and strand in ("+", "-"):
                if feat == "gene" or _gene_id_key(gid) not in gene_strands:
                    gene_strands[_gene_id_key(gid)] = strand
            if feat == "gene":
                continue
            tx_id = attrs.get("transcript_id", "")
            if not tx_id:
                continue
            biotype = (
                attrs.get("transcript_type")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or ""
            )
            is_mane = bool(tags & _MANE_TAGS) or attrs.get("mane_select", "").lower() == "true"
            prev = tx_meta.get(tx_id)
            if prev is None:
                tx_meta[tx_id] = (gid, gname, chrom, strand, biotype, is_mane)
            else:
                tx_meta[tx_id] = (
                    prev[0] or gid,
                    prev[1] or gname,
                    prev[2] or chrom,
                    prev[3] or strand,
                    prev[4] or biotype,
                    prev[5] or is_mane,
                )
            if feat == "exon":
                exon_num_raw = attrs.get("exon_number", "0")
                try:
                    exon_num = int(exon_num_raw)
                except ValueError:
                    exon_num = 0
                tx_exons.setdefault(tx_id, []).append(
                    ExonRecord(
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        transcript_id=tx_id,
                        gene_id=gid,
                        exon_number=exon_num,
                        biotype=biotype,
                        is_mane_or_canonical=is_mane,
                    )
                )

    by_chrom_list: dict[str, list[TranscriptModel]] = {}
    mane_set: set[str] = set()
    biotype_map: dict[str, str] = {}

    for tx_id, raw_exons in tx_exons.items():
        if not raw_exons:
            continue
        gid, gname, chrom, strand, biotype, is_mane = tx_meta[tx_id]
        sorted_exons = tuple(sorted(raw_exons, key=lambda e: (e.start, e.end)))
        tx_model = TranscriptModel(
            transcript_id=tx_id,
            gene_id=gid,
            gene_name=gname,
            chrom=chrom,
            strand=strand,
            biotype=biotype,
            is_mane_or_canonical=is_mane,
            exons=sorted_exons,
            tx_start=sorted_exons[0].start,
            tx_end=sorted_exons[-1].end,
        )
        by_chrom_list.setdefault(chrom, []).append(tx_model)
        if not chrom.startswith("chr"):
            by_chrom_list.setdefault(f"chr{chrom}", []).append(tx_model)
        else:
            by_chrom_list.setdefault(chrom[3:], []).append(tx_model)
        if is_mane:
            mane_set.add(tx_id)
            mane_set.add(_tx_id_key(tx_id))
        if biotype:
            biotype_map[tx_id] = biotype
            biotype_map[_tx_id_key(tx_id)] = biotype

    by_chrom = {c: tuple(models) for c, models in by_chrom_list.items()}
    return GtfSpliceIndex(
        by_chrom=by_chrom,
        mane_transcripts=frozenset(mane_set),
        transcript_biotypes=biotype_map,
        gene_strands=gene_strands,
    )


def _biotype_priority(biotype: str) -> int:
    """Lower rank means higher priority when selecting primary transcript consequence."""
    bt = (biotype or "").lower()
    if bt == "protein_coding":
        return 0
    if bt in ("nonsense_mediated_decay", "non_stop_decay", "polymorphic_pseudogene"):
        return 1
    if bt in ("lncrna", "lincrna", "antisense", "processed_transcript"):
        return 2
    if bt == "retained_intron":
        return 3
    if "pseudogene" in bt:
        return 4
    if bt:
        return 2
    return 5


def compute_transcript_splice_metrics(tx: TranscriptModel, pos: int) -> SpliceMetrics | None:
    """Compute splice distance and side for breakpoint ``pos`` within ``tx``.

    Returns None when ``pos`` lies outside ``[tx.tx_start, tx.tx_end]``.
    For ``+`` strand transcripts:
      * Distance to splice donor (upstream exon 3' end): ``pos - exon_end``
      * Distance to splice acceptor (downstream exon 5' start): ``next_exon_start - pos``
    For ``-`` strand transcripts:
      * The donor is the 3' end of the transcriptionally upstream exon (which sits
        at higher genomic coordinate ``next_exon_start``), so ``next_exon_start - pos``.
      * The acceptor is the 5' start of the transcriptionally downstream exon (which
        sits at lower genomic coordinate ``exon_end``), so ``pos - exon_end``.
    """
    if pos < tx.tx_start or pos > tx.tx_end:
        return None
    exons = tx.exons
    n_exons = len(exons)
    for idx, ex in enumerate(exons):
        if ex.start <= pos <= ex.end:
            exon_rank_num = (idx + 1) if tx.strand != "-" else (n_exons - idx)
            rank_str = f"{exon_rank_num}/{n_exons}"
            # Distance to nearest internal splice junction from inside the exon
            dists: list[tuple[int, str]] = []
            if idx > 0:
                d_left = pos - ex.start
                side_left = "acceptor" if tx.strand != "-" else "donor"
                dists.append((d_left, side_left))
            if idx + 1 < n_exons:
                d_right = ex.end - pos
                side_right = "donor" if tx.strand != "-" else "acceptor"
                dists.append((d_right, side_right))
            min_edge = min((d for d, _ in dists), default=None)
            in_region = min_edge is not None and min_edge <= 3
            inferred = "splice_region_variant" if in_region else (
                "coding_sequence_variant" if tx.biotype == "protein_coding" else "non_coding_transcript_exon_variant"
            )
            return SpliceMetrics(
                splice_dist=0,
                splice_side="exonic",
                splice_region=in_region,
                inferred_csq=inferred,
                transcript_id=tx.transcript_id,
                gene_id=tx.gene_id,
                rank=rank_str,
                dist_to_donor=0 if any(s == "donor" and d == 0 for d, s in dists) else None,
                dist_to_acceptor=0 if any(s == "acceptor" and d == 0 for d, s in dists) else None,
            )

    # Intronic between exons[idx] and exons[idx + 1]
    n_introns = max(1, n_exons - 1)
    for idx in range(n_exons - 1):
        left_ex = exons[idx]
        right_ex = exons[idx + 1]
        if left_ex.end < pos < right_ex.start:
            if tx.strand != "-":
                dist_donor = pos - left_ex.end
                dist_acceptor = right_ex.start - pos
                intron_num = idx + 1
            else:
                dist_donor = right_ex.start - pos
                dist_acceptor = pos - left_ex.end
                intron_num = n_introns - idx
            rank_str = f"{intron_num}/{n_introns}"
            if dist_donor <= dist_acceptor:
                splice_dist = dist_donor
                splice_side = "donor"
            else:
                splice_dist = dist_acceptor
                splice_side = "acceptor"
            in_region = splice_dist <= SPLICE_REGION_MAX_BP
            if splice_dist <= 2:
                inferred = "splice_donor_variant" if splice_side == "donor" else "splice_acceptor_variant"
            elif in_region:
                inferred = "splice_region_variant"
            else:
                inferred = "intron_variant"
            return SpliceMetrics(
                splice_dist=splice_dist,
                splice_side=splice_side,
                splice_region=in_region,
                inferred_csq=inferred,
                transcript_id=tx.transcript_id,
                gene_id=tx.gene_id,
                rank=rank_str,
                dist_to_donor=dist_donor,
                dist_to_acceptor=dist_acceptor,
            )
    return None


def compute_locus_splice_metrics(
    gtf_index: GtfSpliceIndex | None,
    chrom: str,
    pos: int,
    *,
    candidate_tx_ids: Sequence[str] = (),
) -> SpliceMetrics | None:
    """Compute locus-level splice distance and side across overlapping GTF transcripts.

    Prioritizes MANE Select / canonical transcripts first, then ``protein_coding``
    biotype, then minimum ``splice_dist``.
    """
    if gtf_index is None:
        return None
    models = gtf_index.by_chrom.get(chrom, ())
    if not models:
        return None
    tx_filter = {_tx_id_key(t) for t in candidate_tx_ids if t}
    hits: list[tuple[tuple[int, int, int, int], SpliceMetrics]] = []
    for tx in models:
        m = compute_transcript_splice_metrics(tx, pos)
        if m is None or m.splice_dist is None:
            continue
        in_cand = 0 if (not tx_filter or _tx_id_key(tx.transcript_id) in tx_filter) else 1
        mane_rank = 0 if tx.is_mane_or_canonical else 1
        bio_rank = _biotype_priority(tx.biotype)
        key = (mane_rank, bio_rank, in_cand, m.splice_dist)
        hits.append((key, m))
    if not hits:
        return None
    hits.sort(key=lambda item: item[0])
    return hits[0][1]


# --------------------------------------------------------------------------- ANN parsing & isoform prioritization


def parse_snpeff_ann(
    ann: str,
    *,
    mane_transcripts: Iterable[str] = (),
) -> GeneAnnotation:
    """Parse a snpEff ``ANN`` value while retaining transcript isoforms, ranks, and MANE priority.

    Genes and transcripts are taken from ``Feature_Type == transcript`` entries
    only, so intergenic entries -- whose Gene_Name is ``GENE1-GENE2`` -- do not
    leak flanking genes into ``GENE`` or ``TX_ID``.

    When selecting the primary consequence (``CSQ``) and primary rank
    (``TX_RANK`` / ``EXON_RANK`` / ``INTRON_RANK``), overlapping transcripts
    (intronic/exonic/UTR/splice) are prioritized over flanking
    ``upstream_gene_variant`` / ``downstream_gene_variant`` entries, and
    MANE Select / canonical or ``protein_coding`` transcripts take precedence
    over ``retained_intron`` or non-coding isoforms.
    """
    mane_set = {t for t in mane_transcripts} | {_tx_id_key(t) for t in mane_transcripts}
    gene_ids: list[str] = []
    gid_to_symbol: dict[str, str] = {}
    terms: list[str] = []
    tx_ids: list[str] = []
    tx_biotypes: list[str] = []
    tx_entries: list[TranscriptAnn] = []
    n_tx = 0

    for entry in ann.split(","):
        cols = entry.split("|")
        if len(cols) <= ANN_FEATURE_TYPE:
            continue
        raw_terms = [
            SNPEFF_TERM_ALIASES.get(raw, raw)
            for raw in cols[ANN_ANNOTATION].split("&")
            if raw
        ]
        for term in raw_terms:
            if term not in terms:
                terms.append(term)
        if cols[ANN_FEATURE_TYPE] == "transcript":
            n_tx += 1
            sym = cols[ANN_GENE_NAME] if len(cols) > ANN_GENE_NAME else ""
            gid = cols[ANN_GENE_ID] if len(cols) > ANN_GENE_ID else ""
            tx_id = cols[ANN_FEATURE_ID] if len(cols) > ANN_FEATURE_ID else ""
            biotype = cols[ANN_BIOTYPE] if len(cols) > ANN_BIOTYPE else ""
            rank = cols[ANN_RANK] if len(cols) > ANN_RANK else ""
            impact = cols[ANN_IMPACT] if len(cols) > ANN_IMPACT else ""
            if gid and gid not in gene_ids:
                gene_ids.append(gid)
            if gid and sym and gid not in gid_to_symbol:
                gid_to_symbol[gid] = sym
            if tx_id and tx_id not in tx_ids:
                tx_ids.append(tx_id)
            if biotype and biotype not in tx_biotypes:
                tx_biotypes.append(biotype)
            entry_severe = most_severe_term(raw_terms)
            tx_entries.append(
                TranscriptAnn(
                    allele=cols[ANN_ALLELE],
                    terms=tuple(raw_terms),
                    most_severe=entry_severe,
                    impact=impact,
                    gene_name=sym,
                    gene_id=gid,
                    feature_type=cols[ANN_FEATURE_TYPE],
                    transcript_id=tx_id,
                    biotype=biotype,
                    rank=rank,
                )
            )

    symbols = [gid_to_symbol.get(gid, gid) for gid in gene_ids]
    ranked_terms = sorted(terms, key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))

    primary_csq = most_severe_term(terms)
    primary_rank = ""
    exon_rank = ""
    intron_rank = ""

    if tx_entries:
        def _tx_sort_key(t: TranscriptAnn) -> tuple[int, int, int, int]:
            is_flank = 1 if t.most_severe in _UP_DOWN_INTERGENIC_TERMS else 0
            is_mane = 0 if (t.transcript_id in mane_set or _tx_id_key(t.transcript_id) in mane_set) else 1
            bio_pri = _biotype_priority(t.biotype)
            sev_pri = _SEVERITY_RANK.get(t.most_severe, len(CONSEQUENCE_SEVERITY))
            # If no MANE set is active and biotypes are equal, severity wins;
            # when MANE or protein_coding vs retained_intron differs among body-overlapping
            # transcripts, prioritize MANE / protein_coding first, then severity.
            return (is_flank, is_mane, bio_pri, sev_pri)

        best_tx = min(tx_entries, key=_tx_sort_key)
        # Use best_tx's consequence if either (a) it is as severe as the global
        # term, or (b) the only more-severe transcript was a non-coding /
        # retained-intron / non-MANE isoform while a MANE or protein-coding
        # body-overlapping transcript exists.
        if best_tx.most_severe not in _UP_DOWN_INTERGENIC_TERMS:
            primary_csq = best_tx.most_severe
        primary_rank = best_tx.rank
        if primary_rank:
            if primary_csq in _EXONIC_TERMS:
                exon_rank = primary_rank
            elif primary_csq in _INTRONIC_OR_SPLICE_TERMS:
                intron_rank = primary_rank

    return GeneAnnotation(
        gene_symbols=tuple(symbols),
        gene_ids=tuple(gene_ids),
        most_severe=primary_csq,
        terms=tuple(ranked_terms),
        n_transcripts=n_tx,
        transcript_ids=tuple(tx_ids),
        transcript_biotypes=tuple(tx_biotypes),
        primary_rank=primary_rank,
        exon_rank=exon_rank,
        intron_rank=intron_rank,
        transcripts=tuple(tx_entries),
    )


# --------------------------------------------------------------------------- Sequence-level cryptic splice motif scanner


_COMPLEMENT_TABLE = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")

# Shapiro & Senapathy (1987) / MaxEnt-style position weights for the 9-mer 5' splice donor
# [-3, -2, -1 | +1(G), +2(T), +3, +4, +5, +6] with invariant GT at indices 3..4.
# Consensus: (C/A) A G | G T (A/G) A G T
_DONOR_POS_WEIGHTS: tuple[dict[str, float], ...] = (
    {"C": 0.38, "A": 0.34, "G": 0.15, "T": 0.13},  # -3
    {"A": 0.62, "C": 0.14, "G": 0.12, "T": 0.12},  # -2
    {"G": 0.80, "A": 0.09, "T": 0.06, "C": 0.05},  # -1
    {"G": 1.00},                                     # +1 (invariant G)
    {"T": 1.00},                                     # +2 (invariant T)
    {"A": 0.54, "G": 0.38, "C": 0.04, "T": 0.04},  # +3 (purine R)
    {"A": 0.70, "G": 0.12, "T": 0.10, "C": 0.08},  # +4
    {"G": 0.82, "A": 0.08, "T": 0.06, "C": 0.04},  # +5
    {"T": 0.48, "C": 0.20, "A": 0.17, "G": 0.15},  # +6
)
_DONOR_MIN_RAW = sum(min(d.values()) for d in _DONOR_POS_WEIGHTS)
_DONOR_MAX_RAW = sum(max(d.values()) for d in _DONOR_POS_WEIGHTS)

SPLICE_MOTIF_SCORE_THRESHOLD = 0.70
PSEUDOEXON_MIN_LEN_BP = 40
PSEUDOEXON_MAX_LEN_BP = 350


@dataclass(frozen=True)
class CrypticSpliceSite:
    """Single scored 5' donor or 3' acceptor motif hit on a transcribed RNA strand (0-based)."""

    kind: str          # "donor" | "acceptor"
    junction_pos: int  # 0-based exon/intron boundary coordinate in transcribed RNA
    score: float       # normalized [0.0, 1.0]
    motif: str         # 9-mer (donor) or 23-mer (acceptor)


@dataclass(frozen=True)
class CrypticSpliceScanResult:
    """Summary of cryptic 5' donor, 3' acceptor, pseudo-exon pairs, and poly(A) signals."""

    sequence_len: int
    donor_sites: tuple[CrypticSpliceSite, ...]
    acceptor_sites: tuple[CrypticSpliceSite, ...]
    donor_count: int
    acceptor_count: int
    pseudoexon_pairs: int
    max_donor_score: float
    max_acceptor_score: float
    max_score: float
    polya_signal_count: int


def reverse_complement_dna(seq: str) -> str:
    """Return the DNA reverse-complement of ``seq`` (U converted to T)."""
    return seq.translate(_COMPLEMENT_TABLE)[::-1].upper()


def score_donor_9mer(mer9: str) -> float:
    """Score a 9-bp 5' splice donor window (3 exonic + GT + 4 intronic) in ``[0.0, 1.0]``.

    Returns ``0.0`` if the dinucleotide at positions +1..+2 (indices 3..5) is not ``GT``
    or if non-ACGT bases are present.
    """
    s = mer9.upper().replace("U", "T")
    if len(s) != 9 or s[3:5] != "GT":
        return 0.0
    raw = 0.0
    for idx, ch in enumerate(s):
        w = _DONOR_POS_WEIGHTS[idx].get(ch)
        if w is None:
            return 0.0
        raw += w
    norm = (raw - _DONOR_MIN_RAW) / (_DONOR_MAX_RAW - _DONOR_MIN_RAW)
    return round(max(0.0, min(1.0, norm)), 4)


def score_acceptor_23mer(mer23: str) -> float:
    """Score a 23-bp 3' splice acceptor window (18-bp PPT/branch + C/T + AG + 2 exonic) in ``[0.0, 1.0]``.

    Layout of ``mer23``:
      * ``mer23[0:17]`` (positions -20..-4): polypyrimidine tract (C/T density, optimal >= 75%,
        penalized for premature AG dinucleotides inside the last 12 bp of the PPT).
      * ``mer23[17]`` (position -3): pyrimidine preference (C > T >> A > G).
      * ``mer23[18:20]`` (positions -2..-1): invariant ``AG`` acceptor dinucleotide.
      * ``mer23[20]`` (position +1, first exonic base): purine preference (G > A > T/C).
      * ``mer23[21:23]`` (positions +2..+3): minor exonic context.
    """
    s = mer23.upper().replace("U", "T")
    if len(s) != 23 or s[18:20] != "AG":
        return 0.0
    for ch in s:
        if ch not in ("A", "C", "G", "T"):
            return 0.0
    ppt = s[0:17]
    pyr_count = sum(1 for ch in ppt if ch in ("C", "T"))
    pyr_frac = pyr_count / 17.0
    # Normalize pyrimidine fraction so 50% -> 0.25 and >=82% -> 1.0
    ppt_score = min(1.0, max(0.0, (pyr_frac - 0.30) / 0.52))
    # AG exclusion zone in positions -15..-4 (indices 5..17)
    ag_penalty = 0.18 if "AG" in s[5:17] else 0.0

    pos_m3_weights = {"C": 1.0, "T": 0.90, "A": 0.20, "G": 0.05}
    pos_p1_weights = {"G": 1.0, "A": 0.70, "T": 0.35, "C": 0.30}
    m3_score = pos_m3_weights[s[17]]
    p1_score = pos_p1_weights[s[20]]

    combined = 0.60 * ppt_score + 0.25 * m3_score + 0.15 * p1_score - ag_penalty
    return round(max(0.0, min(1.0, combined)), 4)


def scan_cryptic_splice_motifs(
    seq: str,
    *,
    orient_to_gene: str = "sense",
    min_score: float = SPLICE_MOTIF_SCORE_THRESHOLD,
    pseudoexon_min_bp: int = PSEUDOEXON_MIN_LEN_BP,
    pseudoexon_max_bp: int = PSEUDOEXON_MAX_LEN_BP,
) -> CrypticSpliceScanResult:
    """Scan an assembled insertion sequence for cryptic 5' donor and 3' acceptor motifs.

    Parameters
    ----------
    seq : str
        Assembled mobile-element insertion sequence (in the element's orientation).
    orient_to_gene : {"sense", "antisense"}
        When ``"antisense"``, ``seq`` is reverse-complemented before scanning so
        motifs are evaluated on the host gene's transcribed pre-mRNA strand.
    min_score : float
        Minimum normalized score (default ``0.70``) for counting a cryptic site.
    """
    cleaned = re.sub(r"[^ACGTUacgtuNn]", "", seq or "").upper().replace("U", "T")
    if orient_to_gene == "antisense":
        rna = reverse_complement_dna(cleaned)
    else:
        rna = cleaned
    n = len(rna)
    donors: list[CrypticSpliceSite] = []
    acceptors: list[CrypticSpliceSite] = []
    max_d = 0.0
    max_a = 0.0

    if n >= 9:
        for i in range(n - 8):
            if rna[i + 3 : i + 5] != "GT":
                continue
            sc = score_donor_9mer(rna[i : i + 9])
            if sc > max_d:
                max_d = sc
            if sc >= min_score:
                # Junction is after the 3rd exonic base (index i + 3)
                donors.append(CrypticSpliceSite("donor", i + 3, sc, rna[i : i + 9]))

    if n >= 23:
        for i in range(n - 22):
            if rna[i + 18 : i + 20] != "AG":
                continue
            sc = score_acceptor_23mer(rna[i : i + 23])
            if sc > max_a:
                max_a = sc
            if sc >= min_score:
                # Exon starts at index i + 20 (immediately after AG)
                acceptors.append(CrypticSpliceSite("acceptor", i + 20, sc, rna[i : i + 23]))

    # Count internal pseudo-exon pairs: a 3' acceptor at A_pos followed by a 5' donor at D_pos
    # such that pseudoexon_min_bp <= (D_pos - A_pos) <= pseudoexon_max_bp.
    pairs = 0
    for acc in acceptors:
        for don in donors:
            exon_len = don.junction_pos - acc.junction_pos
            if pseudoexon_min_bp <= exon_len <= pseudoexon_max_bp:
                pairs += 1

    polya_hits = len(re.findall(r"A[AT]TAAA", rna))
    max_overall = round(max(max_d, max_a), 4)
    return CrypticSpliceScanResult(
        sequence_len=n,
        donor_sites=tuple(donors),
        acceptor_sites=tuple(acceptors),
        donor_count=len(donors),
        acceptor_count=len(acceptors),
        pseudoexon_pairs=pairs,
        max_donor_score=round(max_d, 4),
        max_acceptor_score=round(max_a, 4),
        max_score=max_overall,
        polya_signal_count=polya_hits,
    )


def _extract_assembled_sequence(info: Mapping[str, str | None]) -> str:
    """Return the assembled insertion nucleotide sequence from ``info`` if present."""
    for key in ("CONSENSUS_SEQ", "SEQ", "INS_SEQ", "ASSEMBLED_SEQ"):
        raw = (info.get(key) or "").strip()
        if raw and raw != ".":
            cleaned = re.sub(r"[^ACGTUacgtuNn]", "", raw)
            if len(cleaned) >= 20:
                return cleaned
    return ""


@lru_cache(maxsize=8)
def load_panel_sequences(path: str | Path) -> dict[str, tuple[str, str]]:
    """Load FASTA sequences from a Dfam/Repbase or full-consensus MEI panel.

    Returns a dictionary mapping both the raw header label (e.g. ``AluY#SINE/Alu``)
    and its stripped subfamily name (e.g. ``AluY``, ``AluY_short``) to
    ``(mei_class, sequence_uppercase)``.
    """
    p = Path(path)
    opener = gzip.open if p.name.endswith(".gz") else open
    seqs: dict[str, tuple[str, str]] = {}
    raw_head: str | None = None
    name: str | None = None
    mei_class = ""
    chunks: list[str] = []

    def _store() -> None:
        if not name:
            return
        seq_str = "".join(chunks).upper()
        seqs[name] = (mei_class, seq_str)
        stripped = name.rstrip("_")
        if stripped and stripped not in seqs:
            seqs[stripped] = (mei_class, seq_str)
        if raw_head:
            seqs[raw_head] = (mei_class, seq_str)

    with opener(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _store()
                raw_head = line[1:].strip().split()[0].split("|", 1)[0]
                name, _, mei_class = raw_head.partition("#")
                chunks = []
            elif name is not None:
                chunks.append(re.sub(r"[^ACGTUacgtuNn]", "", line.strip()))
    _store()
    return seqs


def slice_panel_consensus_span(
    info: Mapping[str, str | None],
    *,
    family: str = "",
    panel_paths: Sequence[str | Path] | None = None,
) -> str:
    """Slice the 1-based inclusive ``[MEI_5P..MEI_3P]`` span from the Dfam/Repbase consensus panel.

    Resolution order:
      1. Full-length canonical axis (``full_axis_name(subfamily)``, e.g. ``L1HS_full``,
         ``AluY_full``, ``SVA_F_full``) because ``MEI_5P`` and ``MEI_3P`` are reported
         on the full consensus axis.
      2. Exact subfamily match in the full or fragment panel (e.g. ``AluY``, ``SVA_D``,
         ``L1HS_5end``).
      3. Family canonical full record (``FAMILY_CANONICAL_FULL[family]``).
    """
    from retro_miner.mei_panel_index import (
        FAMILY_CANONICAL_FULL,
        find_panel_paths,
        full_axis_name,
    )

    raw_sub = (info.get("MEISUBFAMILY") or "").strip()
    if not raw_sub and info.get("MEINFO"):
        raw_sub = (info.get("MEINFO") or "").split(",", 1)[0].strip()
    sub_name = raw_sub.split("|", 1)[0].partition("#")[0].strip()
    sub_clean = sub_name.rstrip("_")

    if not sub_clean and not family:
        return ""

    if panel_paths is not None:
        resolved_paths: list[Path] = [Path(p) for p in panel_paths if p and Path(p).is_file()]
    else:
        resolved_paths = [p for _kind, p in find_panel_paths()]

    if not resolved_paths:
        return ""

    merged: dict[str, tuple[str, str]] = {}
    # Load fragment panels first, then full panels so full-axis records take precedence on key collisions
    for p in resolved_paths:
        merged.update(load_panel_sequences(p))

    candidates: list[str] = []
    if sub_clean:
        candidates.append(full_axis_name(sub_clean))
    if sub_name:
        candidates.append(full_axis_name(sub_name))
        candidates.append(sub_name)
    if sub_clean:
        candidates.append(sub_clean)
    if raw_sub:
        candidates.append(raw_sub)
    if family and family in FAMILY_CANONICAL_FULL:
        candidates.append(FAMILY_CANONICAL_FULL[family])
        base_fam = FAMILY_CANONICAL_FULL[family].removesuffix("_full")
        candidates.append(base_fam)

    target_seq = ""
    for cand in candidates:
        if cand in merged:
            target_seq = merged[cand][1]
            break

    if not target_seq:
        return ""

    seq_len = len(target_seq)
    try:
        start_1b = int(float(info.get("MEI_5P") or "1"))
    except ValueError:
        start_1b = 1
    try:
        end_1b = int(float(info.get("MEI_3P") or str(seq_len)))
    except ValueError:
        end_1b = seq_len

    if start_1b > end_1b:
        start_1b, end_1b = end_1b, start_1b
    start_idx = max(0, min(seq_len, start_1b - 1))
    end_idx = max(start_idx, min(seq_len, end_1b))
    sliced = target_seq[start_idx:end_idx]
    return sliced if len(sliced) >= 20 else ""


def _resolve_sequence_for_motif_scan(
    info: Mapping[str, str | None],
    *,
    family: str = "",
    panel_paths: Sequence[str | Path] | None = None,
) -> str:
    """Return assembled ``CONSENSUS_SEQ`` if present, otherwise slice ``(MEISUBFAMILY, MEI_5P, MEI_3P)`` from panel."""
    assembled = _extract_assembled_sequence(info)
    if assembled:
        return assembled
    return slice_panel_consensus_span(info, family=family, panel_paths=panel_paths)


# --------------------------------------------------------------------------- Twin-priming 5' inversion junction modeling


INV_JUNCTION_WINDOW_RADIUS_BP = 25


@dataclass(frozen=True)
class TwinPrimingInversionSpec:
    """Descriptor of a twin-priming 5' inversion on a LINE-1 or SVA insertion."""

    is_inverted: bool
    inv_5p: int | None = None
    inv_3p: int | None = None
    noninv_5p: int | None = None
    noninv_3p: int | None = None
    junction_offset: int | None = None


@dataclass(frozen=True)
class InversionJunctionScanResult:
    """Summary of chimeric cryptic splice sites spanning or adjacent to a 5' inversion junction."""

    is_inverted: bool
    chimeric_sequence: str
    junction_offset: int
    junction_donor_count: int
    junction_acceptor_count: int
    junction_total_count: int
    junction_max_score: float
    full_scan: CrypticSpliceScanResult | None = None


def detect_twin_priming_5p_inversion(info: Mapping[str, str | None]) -> TwinPrimingInversionSpec:
    """Detect whether an insertion carries a twin-priming 5' inversion and extract segment coordinates.

    Recognizes:
      1. Explicit boolean/flag fields: ``TWIN_PRIMING_INV``, ``MEI_5P_INV``, ``INV_5P``
         (as a flag or ``1``/``true``/``yes``/``5p_inv``), or ``ORIENT`` tokens containing
         both strands (``+-``, ``-+``, ``+|-``, ``-|+``, ``5p_inv``).
      2. Explicit inverted + non-inverted segment coordinate fields:
         ``(INV_5P_START / INV_5P_END, NONINV_5P_START / NONINV_3P_END)`` or
         ``(MEI_INV_5P, MEI_INV_3P)``.
      3. Inverted consensus coordinates ``MEI_5P > MEI_3P`` (when ``MEI_3P > 0``),
         or an explicit ``INV_JUNCTION_POS`` / ``INV_BREAKPOINT_OFFSET`` within the
         assembled sequence.
    """
    def _opt_int(*keys: str) -> int | None:
        for k in keys:
            raw = (info.get(k) or "").strip()
            if raw and raw != ".":
                try:
                    return int(float(raw))
                except ValueError:
                    continue
        return None

    junc_offset = _opt_int("INV_JUNCTION_POS", "INV_BREAKPOINT_OFFSET", "INV_JUNCTION_OFFSET")
    inv_5p = _opt_int("MEI_INV_5P", "INV_5P_START", "INV_START")
    inv_3p = _opt_int("MEI_INV_3P", "INV_5P_END", "INV_END")
    noninv_5p = _opt_int("NONINV_5P", "NONINV_5P_START", "DIRECT_5P")
    noninv_3p = _opt_int("NONINV_3P", "NONINV_3P_END", "DIRECT_3P")
    mei_5p = _opt_int("MEI_5P")
    mei_3p = _opt_int("MEI_3P")

    flag_present = False
    for flag_key in ("TWIN_PRIMING_INV", "MEI_5P_INV"):
        if flag_key in info:
            val = info.get(flag_key)
            if val is None or val.strip().lower() in ("1", "true", "yes", "y", "5p_inv", "twin_priming"):
                flag_present = True
    if "INV_5P" in info and inv_5p is None:
        val = info.get("INV_5P")
        if val is None or val.strip().lower() in ("1", "true", "yes", "y", "5p_inv"):
            flag_present = True

    orient_raw = (info.get("ORIENT") or "").strip().lower()
    if any(tok in orient_raw for tok in ("+-", "-+", "+|-", "-|+", "+/-", "-/+", "5p_inv", "inv", "twin")):
        flag_present = True

    svtype_raw = (info.get("SVTYPE") or "").strip().upper()
    meinfo_raw = (info.get("MEINFO") or "").strip().lower()
    if "INV" in svtype_raw or "5p_inv" in meinfo_raw or "twin_prim" in meinfo_raw:
        flag_present = True

    if inv_5p is not None and inv_3p is not None:
        flag_present = True
        if noninv_5p is None and mei_5p is not None and mei_3p is not None and mei_5p <= mei_3p:
            noninv_5p, noninv_3p = mei_5p, mei_3p

    if junc_offset is not None and junc_offset > 0:
        flag_present = True

    # If MEI_5P > MEI_3P (and both positive), treat as a 5'-inverted coordinate span
    if not flag_present and mei_5p is not None and mei_3p is not None and mei_5p > mei_3p > 0:
        flag_present = True
        inv_5p, inv_3p = mei_3p, mei_5p

    if not flag_present:
        return TwinPrimingInversionSpec(is_inverted=False)

    # Populate default twin-priming segment bounds when coordinates are given via MEI_5P..MEI_3P
    if inv_5p is None and inv_3p is None and mei_5p is not None and mei_3p is not None:
        lo, hi = min(mei_5p, mei_3p), max(mei_5p, mei_3p)
        mid = (lo + hi) // 2
        if hi - lo >= 40:
            inv_5p, inv_3p = lo, mid
            if noninv_5p is None:
                noninv_5p, noninv_3p = mid + 1, hi

    return TwinPrimingInversionSpec(
        is_inverted=True,
        inv_5p=inv_5p,
        inv_3p=inv_3p,
        noninv_5p=noninv_5p,
        noninv_3p=noninv_3p,
        junction_offset=junc_offset,
    )


def build_twin_priming_chimeric_sequence(
    info: Mapping[str, str | None],
    spec: TwinPrimingInversionSpec,
    *,
    family: str = "",
    panel_paths: Sequence[str | Path] | None = None,
) -> tuple[str, int]:
    """Construct the chimeric ``[RC(5'_inverted_segment) + 3'_non_inverted_segment]`` sequence
    and return ``(chimeric_seq_in_element_3p_orientation, junction_offset_0based)``.

    In twin priming (Ostertag & Kazazian 2001), the 3' segment integrates in the
    canonical orientation (defining the insertion's 3' poly(A) strand ``ORIENT``),
    while the 5' segment is reverse-complemented relative to the 3' segment,
    forming an internal head-to-head / tail-to-tail inversion junction.
    """
    if not spec.is_inverted:
        return "", 0

    # 1. If an assembled sequence is provided along with an explicit junction offset,
    #    the assembled sequence is already chimeric around ``junction_offset``.
    assembled = _extract_assembled_sequence(info)
    if assembled and spec.junction_offset is not None and 0 < spec.junction_offset < len(assembled):
        return assembled, spec.junction_offset

    # 2. If explicit inverted and non-inverted segments are provided as sequences:
    inv_seq_raw = re.sub(r"[^ACGTUacgtuNn]", "", (info.get("INV_5P_SEQ") or "").upper())
    noninv_seq_raw = re.sub(r"[^ACGTUacgtuNn]", "", (info.get("NONINV_3P_SEQ") or "").upper())
    if len(inv_seq_raw) >= 10 and len(noninv_seq_raw) >= 10:
        rc_inv = reverse_complement_dna(inv_seq_raw)
        return rc_inv + noninv_seq_raw, len(rc_inv)

    # 3. Slice the 5' inverted and 3' non-inverted segments from the consensus panel
    if spec.inv_5p is not None and spec.inv_3p is not None:
        inv_info = dict(info, MEI_5P=str(min(spec.inv_5p, spec.inv_3p)), MEI_3P=str(max(spec.inv_5p, spec.inv_3p)))
        inv_slice = slice_panel_consensus_span(inv_info, family=family, panel_paths=panel_paths)
        noninv_slice = ""
        if spec.noninv_5p is not None and spec.noninv_3p is not None:
            noninv_info = dict(
                info,
                MEI_5P=str(min(spec.noninv_5p, spec.noninv_3p)),
                MEI_3P=str(max(spec.noninv_5p, spec.noninv_3p)),
            )
            noninv_slice = slice_panel_consensus_span(noninv_info, family=family, panel_paths=panel_paths)
        if inv_slice and noninv_slice:
            rc_inv = reverse_complement_dna(inv_slice)
            return rc_inv + noninv_slice, len(rc_inv)
        if inv_slice and assembled:
            rc_inv = reverse_complement_dna(inv_slice)
            return rc_inv + assembled, len(rc_inv)

    # 4. If assembled sequence is present without an explicit junction offset, split at midpoint
    #    (or construct the twin-primed chimera by inverting its 5' half if INV_APPLY_RC=1)
    if assembled:
        mid = len(assembled) // 2
        if (info.get("INV_RAW_UNINVERTED") or "").strip().lower() in ("1", "true", "yes"):
            rc_inv = reverse_complement_dna(assembled[:mid])
            return rc_inv + assembled[mid:], len(rc_inv)
        return assembled, mid

    return "", 0


def scan_twin_priming_inversion_junction(
    info: Mapping[str, str | None],
    *,
    family: str = "",
    orient_to_gene: str = "sense",
    panel_paths: Sequence[str | Path] | None = None,
    min_score: float = SPLICE_MOTIF_SCORE_THRESHOLD,
    junction_radius_bp: int = INV_JUNCTION_WINDOW_RADIUS_BP,
) -> InversionJunctionScanResult:
    """Detect a twin-priming 5' inversion and scan its chimeric junction for cryptic splice sites.

    In a twin-primed insertion, the 5' segment is inverted relative to the 3'
    segment:
      * When the insertion (3' poly-A segment) is ``sense`` to the host gene, the
        5' inverted segment is transcribed in the **antisense** orientation of
        the consensus element (creating antisense-L1 / antisense-SVA cryptic
        splice sites) and joins the sense 3' segment at the inversion breakpoint.
      * When the insertion (3' poly-A segment) is ``antisense`` to the host gene,
        the 5' inverted segment is transcribed in the **sense** orientation of
        the consensus element (exposing sense L1 5'UTR SD 54 / SD 97 or sense SVA
        Alu-like/VNTR acceptors even in an nominally antisense call!).
      * Furthermore, the sequence juxtaposition across the inversion breakpoint
        itself frequently creates de novo chimeric 5' donor (``GT``) or 3'
        acceptor (``PPT + AG``) splice motifs within ``±junction_radius_bp`` of
        the breakpoint.
    """
    spec = detect_twin_priming_5p_inversion(info)
    if not spec.is_inverted:
        return InversionJunctionScanResult(
            is_inverted=False,
            chimeric_sequence="",
            junction_offset=0,
            junction_donor_count=0,
            junction_acceptor_count=0,
            junction_total_count=0,
            junction_max_score=0.0,
            full_scan=None,
        )

    chimeric_seq, junc_offset = build_twin_priming_chimeric_sequence(
        info, spec, family=family, panel_paths=panel_paths
    )
    if not chimeric_seq:
        return InversionJunctionScanResult(
            is_inverted=True,
            chimeric_sequence="",
            junction_offset=0,
            junction_donor_count=0,
            junction_acceptor_count=0,
            junction_total_count=0,
            junction_max_score=0.0,
            full_scan=None,
        )

    # Scan the chimeric sequence on the host gene's transcribed RNA strand.
    # Because chimeric_seq is in the 3' segment's orientation:
    #   - if orient_to_gene == "sense", transcribed RNA is chimeric_seq (so 5' segment is RC(5') and 3' is forward).
    #   - if orient_to_gene == "antisense", transcribed RNA is RC(chimeric_seq), which places the
    #     junction at (len(chimeric_seq) - junc_offset) and restores the 5' segment to forward sense!
    full_scan = scan_cryptic_splice_motifs(chimeric_seq, orient_to_gene=orient_to_gene, min_score=min_score)
    seq_len = full_scan.sequence_len
    rna_junc = (seq_len - junc_offset) if orient_to_gene == "antisense" else junc_offset

    lo = max(0, rna_junc - junction_radius_bp)
    hi = min(seq_len, rna_junc + junction_radius_bp)

    j_donors = [d for d in full_scan.donor_sites if lo <= d.junction_pos <= hi]
    j_acceptors = [a for a in full_scan.acceptor_sites if lo <= a.junction_pos <= hi]

    # Also compute the maximum raw donor/acceptor motif score across any window overlapping rna_junc
    rna_seq = reverse_complement_dna(chimeric_seq) if orient_to_gene == "antisense" else chimeric_seq.upper()
    max_j_score = max(
        [d.score for d in j_donors] + [a.score for a in j_acceptors],
        default=0.0,
    )
    win_start = max(0, rna_junc - junction_radius_bp - 22)
    win_end = min(len(rna_seq), rna_junc + junction_radius_bp + 22)
    sub_rna = rna_seq[win_start:win_end]
    if len(sub_rna) >= 9:
        for i in range(len(sub_rna) - 8):
            sc = score_donor_9mer(sub_rna[i : i + 9])
            if sc > max_j_score:
                max_j_score = sc
    if len(sub_rna) >= 23:
        for i in range(len(sub_rna) - 22):
            sc = score_acceptor_23mer(sub_rna[i : i + 23])
            if sc > max_j_score:
                max_j_score = sc

    return InversionJunctionScanResult(
        is_inverted=True,
        chimeric_sequence=chimeric_seq,
        junction_offset=rna_junc,
        junction_donor_count=len(j_donors),
        junction_acceptor_count=len(j_acceptors),
        junction_total_count=len(j_donors) + len(j_acceptors),
        junction_max_score=round(max_j_score, 4),
        full_scan=full_scan,
    )


# --------------------------------------------------------------------------- Exonization & truncation risk modeling


def _normalize_family(info: Mapping[str, str | None], alt: str) -> str:
    raw = (info.get("MEIFAMILY") or "").strip().upper()
    if not raw and alt.startswith("<INS:ME:") and alt.endswith(">"):
        raw = alt[len("<INS:ME:"):-1].strip().upper()
    if raw in ("ALU",):
        return "ALU"
    if raw in ("LINE1", "L1", "LINE-1", "LINE"):
        return "LINE1"
    if raw in ("SVA",):
        return "SVA"
    return raw


def _l1_preserves_polya_or_donor(info: Mapping[str, str | None]) -> bool:
    """True when a LINE-1 insertion preserves its 5'UTR splice donor or 3' poly(A)/tail signal."""
    seq_fields = (info.get("SEQ") or "", info.get("CONSENSUS_SEQ") or "", info.get("MEINFO") or "")
    for s in seq_fields:
        upper = s.upper()
        if "AATAAA" in upper or "ATTAAA" in upper:
            return True
    try:
        polya_bp = int(float(info.get("POLYA_MIN_BP") or "0"))
        if polya_bp >= 10:
            return True
    except ValueError:
        pass
    try:
        mei_5p = int(float(info.get("MEI_5P") or "0"))
        if 0 < mei_5p <= 100:  # Preserves SD 54 / SD 97 in L1 5'UTR
            return True
    except ValueError:
        pass
    try:
        mei_3p = int(float(info.get("MEI_3P") or "0"))
        if mei_3p >= 5900:     # Reaches 3'UTR poly(A) signal
            return True
    except ValueError:
        pass
    # Default for intronic sense LINE-1 when truncation coordinates are not restrictive:
    # if neither MEI_5P nor MEI_3P nor POLYA_MIN_BP is provided, assume intact poly(A)/donor risk.
    if not any(k in info for k in ("MEI_5P", "MEI_3P", "POLYA_MIN_BP", "SEQ", "CONSENSUS_SEQ")):
        return True
    return False


def _sva_preserves_sense_trapping_domains(info: Mapping[str, str | None]) -> bool:
    """True when a sense SVA insertion preserves its 5' Alu-like/VNTR splice acceptors or SINE-R poly(A) signal.

    Full-length SVA (~1.3-1.5 kb) comprises:
      1. 5' (CCCTCT)n hexamer + antisense Alu-like domain (~1..400 bp) with strong
         cryptic 3' splice acceptors and branchpoint/polypyrimidine tracts.
      2. GC-rich VNTR (~400..900 bp) with dense cryptic 3'SA / 5'SD pairs.
      3. 3' SINE-R (HERV-K10 LTR-derived, ~900..1366 bp) + canonical AATAAA poly(A)
         signal + poly(A) tail that drives exon trapping / premature cleavage.
    """
    seq_fields = (info.get("SEQ") or "", info.get("CONSENSUS_SEQ") or "", info.get("MEINFO") or "")
    for s in seq_fields:
        upper = s.upper()
        if "AATAAA" in upper or "ATTAAA" in upper or "CCCTCT" in upper:
            return True
    try:
        polya_bp = int(float(info.get("POLYA_MIN_BP") or "0"))
        if polya_bp >= 10:
            return True
    except ValueError:
        pass
    try:
        mei_5p = int(float(info.get("MEI_5P") or "0"))
        if 0 < mei_5p <= 450:  # Preserves 5' Alu-like / VNTR cryptic splice acceptors
            return True
    except ValueError:
        pass
    try:
        mei_3p = int(float(info.get("MEI_3P") or "0"))
        if mei_3p >= 1100:     # Reaches SINE-R 3' LTR poly(A) signal
            return True
    except ValueError:
        pass
    try:
        mei_span = int(float(info.get("MEI_SPAN") or info.get("SVLEN") or "0"))
        if mei_span >= 500:
            return True
    except ValueError:
        pass
    if not any(k in info for k in ("MEI_5P", "MEI_3P", "MEI_SPAN", "SVLEN", "POLYA_MIN_BP", "SEQ", "CONSENSUS_SEQ")):
        return True
    return False


def _sva_preserves_antisense_cryptic_splice_domains(info: Mapping[str, str | None]) -> bool:
    """True when an antisense SVA insertion preserves the 5' Alu-like (transcribed as sense Alu)
    or VNTR cryptic splice site cluster, or the 5' poly(U) tract from the inverted poly(A) tail.
    """
    seq_fields = (info.get("SEQ") or "", info.get("CONSENSUS_SEQ") or "", info.get("MEINFO") or "")
    for s in seq_fields:
        upper = s.upper()
        if "AGAGGG" in upper or "TTTTTT" in upper:
            return True
    try:
        mei_5p = int(float(info.get("MEI_5P") or "0"))
        if 0 < mei_5p <= 700:  # Includes Alu-like or VNTR domain
            return True
    except ValueError:
        pass
    try:
        mei_span = int(float(info.get("MEI_SPAN") or info.get("SVLEN") or "0"))
        if mei_span >= 400:
            return True
    except ValueError:
        pass
    if not any(k in info for k in ("MEI_5P", "MEI_3P", "MEI_SPAN", "SVLEN", "POLYA_MIN_BP", "SEQ", "CONSENSUS_SEQ")):
        return True
    return False


def evaluate_splicing_and_exonization_risk(
    *,
    family: str,
    mei_gene_orient: str | None,
    csq: str,
    terms: Sequence[str],
    splice_dist: int | None,
    info: Mapping[str, str | None],
    panel_paths: Sequence[str | Path] | None = None,
) -> dict[str, str]:
    """Compute repeat-family-specific alternative splicing and truncation risk fields.

    * Antisense intronic Alu: cryptic 3' splice acceptors (positions 279 and 435)
      and the poly-U tract from the inverted poly(A) tail trigger pseudo-exon
      exonization -> ``HIGH`` within 500 bp of an exon (or when distance is
      unmeasured), ``MODERATE`` when deep intronic (> 500 bp).
    * Sense intronic LINE-1: premature polyadenylation (canonical AATAAA/ATTAAA)
      or 5'UTR splice donor (SD 54 / SD 97) capture -> ``HIGH`` when poly(A) or
      donor coordinates are preserved, ``MODERATE`` otherwise.
    * Intronic SVA (sense & antisense exon trapping / cryptic splice activation):
      - Sense SVA: 5' antisense Alu-like + VNTR cryptic 3' splice acceptors paired
        with the 3' SINE-R canonical AATAAA polyadenylation signal cause potent
        exon trapping and premature transcript truncation -> ``HIGH`` when 5'
        acceptor/VNTR or 3' SINE-R poly(A) domains are preserved (or within
        500 bp of a splice junction), ``MODERATE`` when severely truncated and
        deep intronic.
      - Antisense SVA: inverted 3' poly(A) tail forms a 5' polypyrimidine tract
        and theVNTR / inverted Alu-like domains provide cryptic 3'SA and 5'SD
        pairs that induce pseudo-exon inclusion -> ``HIGH`` within 500 bp of an
        annotated splice junction when cryptic splice domains are preserved (or
        when distance is unmeasured and domains are preserved), ``MODERATE``
        when deep intronic (> 500 bp) or truncated.
    * Sequence-level motif refinement: when an assembled insertion sequence
      (``CONSENSUS_SEQ``, ``SEQ``, ``INS_SEQ``, or ``ASSEMBLED_SEQ``) is present,
      or when a ``(MEISUBFAMILY, MEI_5P, MEI_3P)`` span can be sliced from the
      Dfam/Repbase consensus panel via :func:`slice_panel_consensus_span`,
      :func:`scan_cryptic_splice_motifs` scans the transcribed pre-mRNA strand
      (reverse-complemented when ``MEI_GENE_ORIENT == "antisense"``) for cryptic
      5' donor 9-mers, 3' acceptor + polypyrimidine 23-mers, internal pseudo-exon
      3'SS->5'SD pairs (40-350 bp), and canonical poly(A) hexamers, emitting
      ``CRYPTIC_DONOR_COUNT``, ``CRYPTIC_ACCEPTOR_COUNT``,
      ``CRYPTIC_PSEUDOEXON_PAIRS``, and ``CRYPTIC_SPLICE_MAX_SCORE`` and refining
      the risk tier accordingly.
    """
    orients = [o.strip() for o in (mei_gene_orient or "").split(",") if o.strip() in ("sense", "antisense")]
    all_terms = {csq, *terms}
    is_intronic = bool(all_terms & _INTRONIC_OR_SPLICE_TERMS)
    near_exon = splice_dist is None or splice_dist <= DEEP_INTRONIC_EXONIZATION_BP

    assembled_seq = _resolve_sequence_for_motif_scan(info, family=family, panel_paths=panel_paths)
    inv_spec = detect_twin_priming_5p_inversion(info) if family in ("LINE1", "SVA") else TwinPrimingInversionSpec(is_inverted=False)
    inv_by_orient: dict[str, InversionJunctionScanResult] = {}
    best_inv: InversionJunctionScanResult | None = None

    scan_by_orient: dict[str, CrypticSpliceScanResult] = {}
    best_scan: CrypticSpliceScanResult | None = None
    if orients:
        for o in dict.fromkeys(orients):
            if inv_spec.is_inverted:
                inv_res = scan_twin_priming_inversion_junction(
                    info, family=family, orient_to_gene=o, panel_paths=panel_paths
                )
                inv_by_orient[o] = inv_res
                if best_inv is None or (
                    inv_res.junction_total_count,
                    inv_res.junction_max_score,
                ) > (
                    best_inv.junction_total_count,
                    best_inv.junction_max_score,
                ):
                    best_inv = inv_res
                if inv_res.full_scan is not None:
                    scan_by_orient[o] = inv_res.full_scan
            if o not in scan_by_orient and assembled_seq:
                scan_by_orient[o] = scan_cryptic_splice_motifs(assembled_seq, orient_to_gene=o)
            res = scan_by_orient.get(o)
            if res is not None and (
                best_scan is None
                or (
                    res.pseudoexon_pairs,
                    res.acceptor_count + res.donor_count,
                    res.max_score,
                )
                > (
                    best_scan.pseudoexon_pairs,
                    best_scan.acceptor_count + best_scan.donor_count,
                    best_scan.max_score,
                )
            ):
                best_scan = res

    alu_risk = "NONE"
    l1_risk = "NONE"
    sva_risk = "NONE"
    mei_risk = "NONE"

    if is_intronic and family == "ALU":
        if "antisense" in orients:
            as_scan = scan_by_orient.get("antisense")
            if as_scan is not None:
                if as_scan.acceptor_count == 0 and as_scan.pseudoexon_pairs == 0:
                    # Sequence lacks functional cryptic 3' acceptors -> downgrade
                    alu_risk = "MODERATE" if near_exon else "NONE"
                elif as_scan.pseudoexon_pairs >= 1 and as_scan.max_score >= 0.80:
                    # Strong complete pseudo-exon (3'SA + 5'SD pair) -> HIGH even in deeper introns
                    alu_risk = "HIGH"
                else:
                    alu_risk = "MODERATE" if (splice_dist is not None and splice_dist > DEEP_INTRONIC_EXONIZATION_BP) else "HIGH"
            else:
                if splice_dist is not None and splice_dist > DEEP_INTRONIC_EXONIZATION_BP:
                    alu_risk = "MODERATE"
                else:
                    alu_risk = "HIGH"
            mei_risk = alu_risk
        elif "sense" in orients:
            s_scan = scan_by_orient.get("sense")
            if s_scan is not None and s_scan.pseudoexon_pairs >= 1 and near_exon:
                alu_risk = "MODERATE"
            else:
                alu_risk = "NONE"
            mei_risk = alu_risk

    elif is_intronic and family == "LINE1":
        l1_candidates: list[str] = []
        if "sense" in orients:
            s_scan = scan_by_orient.get("sense")
            s_inv = inv_by_orient.get("sense")
            if s_inv is not None and s_inv.is_inverted and (
                s_inv.junction_total_count >= 1 or (s_scan is not None and s_scan.pseudoexon_pairs >= 1)
            ):
                l1_candidates.append("HIGH")
            elif s_scan is not None:
                if (
                    s_scan.polya_signal_count >= 1
                    or s_scan.donor_count >= 1
                    or s_scan.pseudoexon_pairs >= 1
                    or _l1_preserves_polya_or_donor(info)
                ):
                    l1_candidates.append("HIGH")
                elif s_scan.acceptor_count >= 1 or (s_inv is not None and s_inv.is_inverted):
                    l1_candidates.append("MODERATE")
                else:
                    l1_candidates.append("NONE" if (splice_dist is not None and splice_dist > DEEP_INTRONIC_EXONIZATION_BP) else "MODERATE")
            else:
                if _l1_preserves_polya_or_donor(info) or (s_inv is not None and s_inv.is_inverted and near_exon):
                    l1_candidates.append("HIGH")
                else:
                    l1_candidates.append("MODERATE")
        if "antisense" in orients:
            as_scan = scan_by_orient.get("antisense")
            as_inv = inv_by_orient.get("antisense")
            if as_inv is not None and as_inv.is_inverted:
                # Antisense twin-primed L1 flips its 5' segment into sense orientation on the
                # pre-mRNA strand and creates chimeric splice junctions at the inversion breakpoint.
                if (
                    as_inv.junction_total_count >= 1
                    or (as_scan is not None and (as_scan.pseudoexon_pairs >= 1 or as_scan.donor_count >= 1))
                ):
                    l1_candidates.append("HIGH" if (near_exon or as_inv.junction_total_count >= 1 or (as_scan is not None and as_scan.pseudoexon_pairs >= 1)) else "MODERATE")
                else:
                    l1_candidates.append("MODERATE")
            elif as_scan is not None:
                if as_scan.pseudoexon_pairs >= 1 and near_exon and as_scan.max_score >= 0.80:
                    l1_candidates.append("HIGH")
                elif as_scan.pseudoexon_pairs >= 1 or (as_scan.acceptor_count >= 1 and as_scan.donor_count >= 1):
                    l1_candidates.append("MODERATE")
                else:
                    l1_candidates.append("NONE")
            else:
                l1_candidates.append("NONE")
        if "HIGH" in l1_candidates:
            l1_risk = "HIGH"
        elif "MODERATE" in l1_candidates:
            l1_risk = "MODERATE"
        mei_risk = l1_risk

    elif is_intronic and family == "SVA":
        candidate_risks: list[str] = []
        if "sense" in orients:
            s_scan = scan_by_orient.get("sense")
            s_inv = inv_by_orient.get("sense")
            if s_inv is not None and s_inv.is_inverted and (
                s_inv.junction_total_count >= 1 or (s_scan is not None and s_scan.pseudoexon_pairs >= 1)
            ):
                candidate_risks.append("HIGH")
            elif s_scan is not None:
                if (
                    s_scan.pseudoexon_pairs >= 1
                    or (s_scan.acceptor_count >= 1 and s_scan.polya_signal_count >= 1)
                    or (_sva_preserves_sense_trapping_domains(info) and (s_scan.acceptor_count + s_scan.donor_count + s_scan.polya_signal_count) > 0)
                ):
                    candidate_risks.append("HIGH")
                elif s_scan.acceptor_count == 0 and s_scan.donor_count == 0 and s_scan.polya_signal_count == 0 and not near_exon and not (s_inv is not None and s_inv.is_inverted):
                    candidate_risks.append("NONE")
                else:
                    candidate_risks.append("MODERATE")
            else:
                if _sva_preserves_sense_trapping_domains(info) or (splice_dist is not None and splice_dist <= DEEP_INTRONIC_EXONIZATION_BP):
                    candidate_risks.append("HIGH")
                else:
                    candidate_risks.append("MODERATE")
        if "antisense" in orients:
            as_scan = scan_by_orient.get("antisense")
            as_inv = inv_by_orient.get("antisense")
            if as_inv is not None and as_inv.is_inverted and (
                as_inv.junction_total_count >= 1 or (as_scan is not None and as_scan.pseudoexon_pairs >= 1)
            ):
                candidate_risks.append("HIGH")
            elif as_scan is not None:
                if as_scan.pseudoexon_pairs >= 1 and (near_exon or as_scan.max_score >= 0.80):
                    candidate_risks.append("HIGH")
                elif as_scan.acceptor_count == 0 and as_scan.donor_count == 0 and not near_exon and not (as_inv is not None and as_inv.is_inverted):
                    candidate_risks.append("NONE")
                elif near_exon and _sva_preserves_antisense_cryptic_splice_domains(info) and (as_scan.acceptor_count + as_scan.donor_count) > 0:
                    candidate_risks.append("HIGH")
                else:
                    candidate_risks.append("MODERATE")
            else:
                if near_exon and _sva_preserves_antisense_cryptic_splice_domains(info):
                    candidate_risks.append("HIGH")
                else:
                    candidate_risks.append("MODERATE")
        if "HIGH" in candidate_risks:
            sva_risk = "HIGH"
        elif "MODERATE" in candidate_risks:
            sva_risk = "MODERATE"
        mei_risk = sva_risk

    elif is_intronic and best_scan is not None:
        if best_scan.pseudoexon_pairs >= 1 and near_exon:
            mei_risk = "HIGH"
        elif best_scan.pseudoexon_pairs >= 1 or (best_scan.acceptor_count >= 1 and best_scan.donor_count >= 1):
            mei_risk = "MODERATE"

    out: dict[str, str | None] = {}
    if is_intronic and orients:
        out["MEI_EXONIZATION_RISK"] = mei_risk
        if family == "ALU":
            out["ALU_EXONIZATION_RISK"] = alu_risk
        elif family == "LINE1":
            out["L1_TRUNCATION_RISK"] = l1_risk
        elif family == "SVA":
            out["SVA_EXON_TRAPPING_RISK"] = sva_risk
        if best_scan is not None:
            out["CRYPTIC_DONOR_COUNT"] = str(best_scan.donor_count)
            out["CRYPTIC_ACCEPTOR_COUNT"] = str(best_scan.acceptor_count)
            out["CRYPTIC_PSEUDOEXON_PAIRS"] = str(best_scan.pseudoexon_pairs)
            out["CRYPTIC_SPLICE_MAX_SCORE"] = f"{best_scan.max_score:.4f}"
        if inv_spec.is_inverted:
            out["TWIN_PRIMING_INV"] = None
            if best_inv is not None and best_inv.chimeric_sequence:
                out["INV_JUNCTION_SPLICE_COUNT"] = str(best_inv.junction_total_count)
                out["INV_JUNCTION_SPLICE_MAX_SCORE"] = f"{best_inv.junction_max_score:.4f}"
    return out  # type: ignore[return-value]


# --------------------------------------------------------------------------- snpEff execution


def snpeff_command(
    in_path: str | Path, genome: str, *, snpeff_bin: str = "snpEff",
    config: str | Path | None = None, xmx: str = "8g",
) -> list[str]:
    cmd = [snpeff_bin, f"-Xmx{xmx}", "-noStats", "-noHgvs", *SNPEFF_SPLICE_ARGS]
    if config is not None:
        cmd += ["-c", str(config)]
    cmd += [genome, str(in_path)]
    return cmd


def _is_gtf_file(path: str | Path | None) -> bool:
    if path is None:
        return False
    p = Path(path)
    if not p.is_file():
        return False
    name = p.name.lower()
    if name.endswith((".gtf", ".gtf.gz", ".gff", ".gff3", ".gff.gz")):
        return True
    opener = gzip.open if name.endswith(".gz") else open
    try:
        with opener(p, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                return len(line.rstrip("\n").split("\t")) >= 9
    except OSError:
        return False
    return False


def annotate_records(
    records: Sequence[VcfRecord],
    headers: Sequence[str],
    genome: str,
    *,
    snpeff_bin: str = "snpEff",
    config: str | Path | None = None,
    xmx: str = "8g",
    gene_strands: Mapping[str, str] | None = None,
    gene_gtf: str | Path | None = None,
    panel_fasta: str | Path | Sequence[str | Path] | None = None,
    run: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> tuple[list[VcfRecord], AnnotationStats]:
    """Run snpEff on ``records`` and translate its ``ANN`` field + GTF splice metrics into INFO fields.

    ``run`` is injectable (defaults to ``subprocess.run``) so tests can substitute
    a canned snpEff output. The raw ``ANN`` field is dropped; run snpEff directly
    for the full per-transcript table. ``GENE_STRAND`` is filled from the bundled
    Ensembl GRCh38.115 table unless ``gene_strands`` or ``gene_gtf`` is passed.
    When ``gene_gtf`` points to a GTF with ``exon`` features, native splice-site
    distances (``SPLICE_DIST``, ``SPLICE_SIDE``, ``SPLICE_REGION``) and MANE
    Select prioritization are computed directly from the GTF exon intervals.
    """
    if run is subprocess.run and shutil.which(snpeff_bin) is None:
        raise FileNotFoundError(f"snpEff launcher {snpeff_bin!r} not found on PATH")
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="rtm_snpeff_") as td:
        tmp_in = Path(td) / "in.vcf"
        write_vcf(tmp_in, headers, records)
        proc = run(
            snpeff_command(tmp_in, genome, snpeff_bin=snpeff_bin, config=config, xmx=xmx),
            capture_output=True, text=True, cwd=td,
        )
        if proc.returncode != 0:
            stderr_tail = (proc.stderr or "")[-2000:]
            stdout_tail = (proc.stdout or "")[-1000:]
            combined = f"{stderr_tail}\n{stdout_tail}"
            oom_indicators = ("OutOfMemoryError", "Java heap space", "GC overhead limit exceeded", "Cannot allocate memory")
            if proc.returncode == 1 or any(tok in combined for tok in oom_indicators):
                raise RuntimeError(
                    f"snpEff exited {proc.returncode} (with -Xmx{xmx}): {stderr_tail or stdout_tail}. "
                    "If snpEff ran out of JVM memory building the splice interval forest, "
                    "increase --snpeff-xmx (e.g. --snpeff-xmx 12g or 16g) or batch the input VCF."
                )
            raise RuntimeError(f"snpEff exited {proc.returncode}: {stderr_tail}")
        _, out_records = _read_vcf_text(proc.stdout)

    gtf_index: GtfSpliceIndex | None = None
    if _is_gtf_file(gene_gtf) and gene_gtf is not None:
        gtf_index = build_gtf_splice_index(gene_gtf)

    if gene_strands is not None:
        strands: Mapping[str, str] = gene_strands
    elif gtf_index is not None and gtf_index.gene_strands:
        strands = gtf_index.gene_strands
    else:
        strands = _resolve_gene_strands(gene_strands, gene_gtf)

    mane_tx = gtf_index.mane_transcripts if gtf_index is not None else ()
    by_key = {_record_key(r): r for r in out_records}
    out: list[VcfRecord] = []
    n_annot = n_gene = n_unmatched = 0
    for rec in records:
        hit = by_key.get(_record_key(rec))
        new_info = dict(rec.info)
        # VCF 4.3/4.4 compliance: for symbolic insertions, END must equal POS.
        if rec.alt.startswith("<INS") and "END" in new_info:
            new_info["END"] = str(rec.pos)
        ann = hit.info.get("ANN") if hit is not None else None
        if not ann:
            n_unmatched += 1
        else:
            summary = parse_snpeff_ann(ann, mane_transcripts=mane_tx)
            orient = rec.info.get("ORIENT")
            info_updates = summary.info_fields(strands, orient=orient)

            splice_metrics = compute_locus_splice_metrics(
                gtf_index,
                rec.chrom,
                rec.pos,
                candidate_tx_ids=summary.transcript_ids,
            )
            effective_csq = summary.most_severe
            effective_terms = list(summary.terms)
            if splice_metrics is not None and splice_metrics.splice_dist is not None:
                info_updates["SPLICE_DIST"] = str(splice_metrics.splice_dist)
                if splice_metrics.dist_to_donor is not None:
                    info_updates["DIST_TO_SPLICE_DONOR"] = str(splice_metrics.dist_to_donor)
                if splice_metrics.dist_to_acceptor is not None:
                    info_updates["DIST_TO_SPLICE_ACCEPTOR"] = str(splice_metrics.dist_to_acceptor)
                if splice_metrics.splice_side:
                    info_updates["SPLICE_SIDE"] = splice_metrics.splice_side
                if splice_metrics.splice_region:
                    new_info["SPLICE_REGION"] = None
                else:
                    new_info.pop("SPLICE_REGION", None)
                if splice_metrics.inferred_csq:
                    if splice_metrics.inferred_csq not in effective_terms:
                        effective_terms.append(splice_metrics.inferred_csq)
                        effective_terms.sort(key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))
                    if _SEVERITY_RANK.get(splice_metrics.inferred_csq, 999) < _SEVERITY_RANK.get(effective_csq, 999):
                        effective_csq = splice_metrics.inferred_csq
                    info_updates["CSQ"] = effective_csq
                    info_updates["CSQ_TERMS"] = ",".join(effective_terms)
                if not info_updates.get("TX_RANK") and splice_metrics.rank:
                    info_updates["TX_RANK"] = splice_metrics.rank
                    if splice_metrics.splice_side == "exonic":
                        info_updates["EXON_RANK"] = splice_metrics.rank
                    else:
                        info_updates["INTRON_RANK"] = splice_metrics.rank
            else:
                # Infer SPLICE_SIDE / SPLICE_REGION from snpEff SO terms when no GTF is provided
                if effective_csq == "splice_donor_variant" or "splice_donor_variant" in effective_terms:
                    info_updates["SPLICE_SIDE"] = "donor"
                    new_info["SPLICE_REGION"] = None
                elif effective_csq == "splice_acceptor_variant" or "splice_acceptor_variant" in effective_terms:
                    info_updates["SPLICE_SIDE"] = "acceptor"
                    new_info["SPLICE_REGION"] = None
                elif effective_csq in {
                    "splice_region_variant",
                    "splice_donor_5th_base_variant",
                    "splice_donor_region_variant",
                    "splice_polypyrimidine_tract_variant",
                }:
                    new_info["SPLICE_REGION"] = None
                elif effective_csq in _EXONIC_TERMS:
                    info_updates["SPLICE_SIDE"] = "exonic"

            new_info.update(info_updates)
            family = _normalize_family(new_info, rec.alt)
            if panel_fasta is None:
                resolved_panels: Sequence[str | Path] | None = None
            elif isinstance(panel_fasta, (str, Path)):
                resolved_panels = (panel_fasta,)
            else:
                resolved_panels = tuple(panel_fasta)
            risk_fields = evaluate_splicing_and_exonization_risk(
                family=family,
                mei_gene_orient=new_info.get("MEI_GENE_ORIENT"),
                csq=effective_csq,
                terms=effective_terms,
                splice_dist=splice_metrics.splice_dist if splice_metrics is not None else None,
                info=new_info,
                panel_paths=resolved_panels,
            )
            new_info.update(risk_fields)
            n_annot += 1
            n_gene += bool(summary.gene_symbols)
        out.append(
            VcfRecord(
                rec.chrom, rec.pos, rec.id, rec.ref, rec.alt,
                rec.qual, rec.filter, new_info, list(rec.rest),
            )
        )
    return out, AnnotationStats(len(records), n_annot, n_gene, n_unmatched, time.time() - t0)


def _read_vcf_text(text: str) -> tuple[list[str], list[VcfRecord]]:
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
        fh.write(text)
        name = fh.name
    try:
        return read_vcf(name)
    finally:
        Path(name).unlink(missing_ok=True)


def annotate_vcf(
    in_path: str | Path,
    out_path: str | Path,
    genome: str,
    *,
    tsv_path: str | Path | None = None,
    snpeff_bin: str = "snpEff",
    config: str | Path | None = None,
    xmx: str = "8g",
    gene_strands: Mapping[str, str] | None = None,
    gene_gtf: str | Path | None = None,
    panel_fasta: str | Path | Sequence[str | Path] | None = None,
    run: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> AnnotationStats:
    """Read ``in_path``, annotate with snpEff, write ``out_path`` (and optional flat TSV)."""
    headers, records = read_vcf(in_path)
    annotated, stats = annotate_records(
        records, headers, genome,
        snpeff_bin=snpeff_bin, config=config, xmx=xmx,
        gene_strands=gene_strands, gene_gtf=gene_gtf,
        panel_fasta=panel_fasta, run=run,
    )
    write_vcf(out_path, add_info_headers(headers), annotated)
    if tsv_path is not None:
        write_annotation_tsv(tsv_path, annotated)
    return stats


TSV_COLUMNS: tuple[str, ...] = ("chrom", "pos", "id", "alt", "svlen", "meinfo", *ANNOTATION_INFO_KEYS)


def _tsv_cell_value(info: Mapping[str, str | None], key: str) -> str:
    if key in ("SPLICE_REGION", "TWIN_PRIMING_INV"):
        return "True" if key in info else ""
    return info.get(key) or ""


def write_annotation_tsv(path: str | Path, records: Iterable[VcfRecord]) -> None:
    """One row per record: coordinates, MEI descriptors, and the gene annotation fields."""
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in records:
            row = [
                r.chrom, str(r.pos), r.id, r.alt,
                r.info.get("SVLEN") or "", r.info.get("MEINFO") or "",
                *[_tsv_cell_value(r.info, k) for k in ANNOTATION_INFO_KEYS],
            ]
            fh.write("\t".join(row) + "\n")
