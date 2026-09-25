"""Gene / functional-consequence annotation of MEI calls via local snpEff.

Consumes a VCF of mobile-element insertions (symbolic ``<INS:ME:*>`` ALT
alleles, as written by ``vcf_export``) and appends gene-level consequence
fields to each record's INFO column.

``CONSEQUENCE_SEVERITY`` reproduces Ensembl's "Calculated variant
consequences" table
(https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html),
most severe first. snpEff does not order ``ANN`` entries that way, so
:func:`parse_snpeff_ann` applies the table itself. Terms are Sequence
Ontology identifiers.

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
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

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
    '##INFO=<ID=CSQ,Number=1,Type=String,Description="Most severe Sequence Ontology consequence term across overlapping transcripts (snpEff)">',
    '##INFO=<ID=CSQ_TERMS,Number=.,Type=String,Description="All distinct consequence terms across overlapping transcripts">',
    '##INFO=<ID=CSQ_NTX,Number=1,Type=Integer,Description="Number of transcripts overlapping the insertion breakpoint">',
)
ANNOTATION_INFO_KEYS: tuple[str, ...] = ("GENE", "GENEID", "CSQ", "CSQ_TERMS", "CSQ_NTX")


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
        info = ";".join(k if v is None else f"{k}={v}" for k, v in self.info.items()) or "."
        cols = [self.chrom, str(self.pos), self.id, self.ref, self.alt, self.qual, self.filter, info, *self.rest]
        return "\t".join(cols)


@dataclass(frozen=True)
class GeneAnnotation:
    """Gene-level summary collapsed from a snpEff ANN field."""

    gene_symbols: tuple[str, ...]
    gene_ids: tuple[str, ...]
    most_severe: str
    terms: tuple[str, ...]
    n_transcripts: int

    def info_fields(self) -> dict[str, str]:
        out: dict[str, str] = {}
        if self.gene_symbols:
            out["GENE"] = ",".join(self.gene_symbols)
        if self.gene_ids:
            out["GENEID"] = ",".join(self.gene_ids)
        out["CSQ"] = self.most_severe
        if self.terms:
            out["CSQ_TERMS"] = ",".join(self.terms)
        out["CSQ_NTX"] = str(self.n_transcripts)
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
# per transcript, 16 pipe-separated columns). Measured 2026-09-18 on the same
# 30 chr22 calls with ``snpEff 5.4c -Xmx8g GRCh38.99``: 57 s wall, of which
# ~55 s is fixed database load and ~1 s annotation; 30/30 annotated; the
# 6 kb L1 at chr22:19223382 evaluated as a point (intron_variant, CLTCL1).
# ``chr22`` input was mapped to database contig ``22``
# and echoed back as ``chr22``. Default JVM heap OOMs on GRCh38 -> pass -Xmx.
# ``-noHgvs`` is always set: an insertion has no codon change, and HGVS
# protein notation is not part of this module's output. Computing it makes
# snpEff load the reference sequence into the heap.
# Splice-site sizes are set to 0. ``Exon.createSpliceSiteRegionEnd`` returns
# before allocating when the size is <= 0; the default sizes allocate a
# splice interval on every exon and threw ``OutOfMemoryError`` in
# ``buildForest`` at ``-Xmx3g``. An insertion is not a splice-site variant.
#
# ANN entries are NOT ordered by Ensembl severity (snpEff lists e.g.
# upstream_gene_variant ahead of intron_variant for the same gene), so the
# parser below applies CONSEQUENCE_SEVERITY itself. ``intergenic_region`` is
# snpEff's name for SO ``intergenic_variant``.

SNPEFF_TERM_ALIASES: dict[str, str] = {
    "intergenic_region": "intergenic_variant",
    "gene_variant": "coding_transcript_variant",
    "transcript_variant": "coding_transcript_variant",
}

ANN_ANNOTATION, ANN_GENE_NAME, ANN_GENE_ID, ANN_FEATURE_TYPE = 1, 3, 4, 5

# Size 0 makes snpEff skip splice-site interval allocation. See module comment above.
SNPEFF_NO_SPLICE: tuple[str, ...] = (
    "-spliceSiteSize", "0",
    "-spliceRegionExonSize", "0",
    "-spliceRegionIntronMin", "0",
    "-spliceRegionIntronMax", "0",
)


def parse_snpeff_ann(ann: str) -> GeneAnnotation:
    """Collapse a snpEff ``ANN`` value to gene symbols, gene IDs, and a consequence summary.

    Genes are taken from transcript entries only (``Feature_Type == transcript``),
    so intergenic entries -- whose Gene_Name is ``GENE1-GENE2`` -- do not leak
    flanking genes into ``GENE``.
    """
    symbols: list[str] = []
    gene_ids: list[str] = []
    terms: list[str] = []
    n_tx = 0
    for entry in ann.split(","):
        cols = entry.split("|")
        if len(cols) <= ANN_FEATURE_TYPE:
            continue
        for raw in cols[ANN_ANNOTATION].split("&"):
            term = SNPEFF_TERM_ALIASES.get(raw, raw)
            if term and term not in terms:
                terms.append(term)
        if cols[ANN_FEATURE_TYPE] == "transcript":
            n_tx += 1
            sym, gid = cols[ANN_GENE_NAME], cols[ANN_GENE_ID]
            if sym and sym not in symbols:
                symbols.append(sym)
            if gid and gid not in gene_ids:
                gene_ids.append(gid)
    ranked = sorted(terms, key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))
    return GeneAnnotation(
        gene_symbols=tuple(symbols), gene_ids=tuple(gene_ids),
        most_severe=most_severe_term(terms), terms=tuple(ranked), n_transcripts=n_tx,
    )


def snpeff_command(
    in_path: str | Path, genome: str, *, snpeff_bin: str = "snpEff",
    config: str | Path | None = None, xmx: str = "8g",
) -> list[str]:
    cmd = [snpeff_bin, f"-Xmx{xmx}", "-noStats", "-noHgvs", *SNPEFF_NO_SPLICE]
    if config is not None:
        cmd += ["-c", str(config)]
    cmd += [genome, str(in_path)]
    return cmd


def annotate_records(
    records: Sequence[VcfRecord],
    headers: Sequence[str],
    genome: str,
    *,
    snpeff_bin: str = "snpEff",
    config: str | Path | None = None,
    xmx: str = "8g",
    run: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> tuple[list[VcfRecord], AnnotationStats]:
    """Run snpEff on ``records`` and translate its ``ANN`` field into the module's INFO fields.

    ``run`` is injectable (defaults to ``subprocess.run``) so tests can substitute
    a canned snpEff output. The raw ``ANN`` field is dropped; run snpEff directly
    for the full per-transcript table.
    """
    if run is subprocess.run and shutil.which(snpeff_bin) is None:
        raise FileNotFoundError(f"snpEff launcher {snpeff_bin!r} not found on PATH")
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="rtm_snpeff_") as td:
        tmp_in = Path(td) / "in.vcf"
        write_vcf(tmp_in, headers, records)
        proc = run(snpeff_command(tmp_in, genome, snpeff_bin=snpeff_bin, config=config, xmx=xmx),
                   capture_output=True, text=True, cwd=td)
        if proc.returncode != 0:
            raise RuntimeError(f"snpEff exited {proc.returncode}: {proc.stderr[-2000:]}")
        _, out_records = _read_vcf_text(proc.stdout)

    by_key = {_record_key(r): r for r in out_records}
    out: list[VcfRecord] = []
    n_annot = n_gene = n_unmatched = 0
    for rec in records:
        hit = by_key.get(_record_key(rec))
        new_info = dict(rec.info)
        ann = hit.info.get("ANN") if hit is not None else None
        if not ann:
            n_unmatched += 1
        else:
            summary = parse_snpeff_ann(ann)
            new_info.update(summary.info_fields())
            n_annot += 1
            n_gene += bool(summary.gene_symbols)
        out.append(VcfRecord(rec.chrom, rec.pos, rec.id, rec.ref, rec.alt, rec.qual, rec.filter, new_info, list(rec.rest)))
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
    run: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> AnnotationStats:
    """Read ``in_path``, annotate with snpEff, write ``out_path`` (and optional flat TSV)."""
    headers, records = read_vcf(in_path)
    annotated, stats = annotate_records(records, headers, genome, snpeff_bin=snpeff_bin, config=config, xmx=xmx, run=run)
    write_vcf(out_path, add_info_headers(headers), annotated)
    if tsv_path is not None:
        write_annotation_tsv(tsv_path, annotated)
    return stats


TSV_COLUMNS: tuple[str, ...] = ("chrom", "pos", "id", "alt", "svlen", "meinfo", *ANNOTATION_INFO_KEYS)


def write_annotation_tsv(path: str | Path, records: Iterable[VcfRecord]) -> None:
    """One row per record: coordinates, MEI descriptors, and the gene annotation fields."""
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for r in records:
            row = [
                r.chrom, str(r.pos), r.id, r.alt,
                r.info.get("SVLEN") or "", r.info.get("MEINFO") or "",
                *[(r.info.get(k) or "") for k in ANNOTATION_INFO_KEYS],
            ]
            fh.write("\t".join(row) + "\n")
