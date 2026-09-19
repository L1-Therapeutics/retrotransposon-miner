"""Gene / functional-consequence annotation of MEI calls via Ensembl VEP.

Consumes a VCF of mobile-element insertions (symbolic ``<INS:ME:*>`` ALT
alleles, as written by ``vcf_export``) and appends gene-level consequence
fields to each record's INFO column using the Ensembl VEP REST endpoint
``POST /vep/:species/region``.

Why REST and not a local install
--------------------------------
Measured on the 30 chr22 GRCh38 calls in ``tests/data/chr22_mei.vcf``
(2026-09-18, rest.ensembl.org): 9.4 s for 30 variants as submitted
(~312 ms/variant), 5.4 s with SVLEN stripped (~180 ms/variant). No local
database, no JVM, no heap tuning. A ~1000-record chromosome-scale callset
therefore annotates in roughly 3-5 minutes.

The SVLEN pitfall (important, verified)
---------------------------------------
For a symbolic insertion VEP recomputes the affected reference span as
``POS + SVLEN`` and *ignores the END field supplied in INFO*. An L1
insertion with ``SVLEN=6018`` at chr22:19223382 was therefore evaluated as
a 6,018 bp *span* (VEP reported start=19223383, end=19229400) and returned
``coding_sequence_variant`` on 38 CLTCL1 transcripts. The insertion is a
point event in reference coordinates (VCF 4.4 §5.6: END == POS for a
symbolic INS); with SVLEN removed the same locus is ``intron_variant`` on
every transcript. Across the 30-record chr22 set, 2/30 loci changed
most-severe consequence because of this. :func:`vep_region_string`
therefore drops SVLEN by default; SVLEN is preserved untouched in the
output VCF because it is correct VCF, just misread by VEP.

Consequence severity ranking
----------------------------
``CONSEQUENCE_SEVERITY`` reproduces the order in Ensembl's
"Calculated variant consequences" table
(https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html),
most severe first. Terms are Sequence Ontology (SO) identifiers, so a
downstream reader can look each one up unambiguously.

References
----------
McLaren W. et al. (2016) The Ensembl Variant Effect Predictor.
    Genome Biology 17:122. doi:10.1186/s13059-016-0974-4
Eilbeck K. et al. (2005) The Sequence Ontology: a tool for the unification
    of genome annotations. Genome Biology 6:R44. doi:10.1186/gb-2005-6-5-r44
VCF specification v4.4, §1.4.7 (INFO reserved keys) and §5.6 (symbolic
    structural variant alleles). https://samtools.github.io/hts-specs/
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

VEP_REST_URL = "https://rest.ensembl.org/vep/human/region"

# rest.ensembl.org rejects POST bodies with more than 200 variants.
VEP_MAX_BATCH = 200

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
    '##INFO=<ID=GENE,Number=.,Type=String,Description="Gene symbols overlapped by the insertion breakpoint (Ensembl VEP)">',
    '##INFO=<ID=GENEID,Number=.,Type=String,Description="Ensembl gene IDs overlapped by the insertion breakpoint">',
    '##INFO=<ID=CSQ,Number=1,Type=String,Description="Most severe Sequence Ontology consequence term across overlapping transcripts (Ensembl VEP)">',
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
    """Gene-level summary of one VEP response record."""

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


# --------------------------------------------------------------------------- VEP request


def vep_region_string(rec: VcfRecord, *, keep_svlen: bool = False) -> str:
    """Render one record in VEP's whitespace-separated VCF-line input format.

    SVLEN is dropped by default: VEP widens a symbolic insertion to a
    POS+SVLEN reference span and ignores END, producing spurious exonic
    consequences (see module docstring). END is forced to POS for symbolic
    insertions so the breakpoint is evaluated as a point.
    """
    info_keys = ("SVTYPE", "END") + (("SVLEN",) if keep_svlen else ())
    parts: list[str] = []
    for k in info_keys:
        if k == "END" and rec.alt.startswith("<INS"):
            parts.append(f"END={rec.pos}")
        elif k in rec.info and rec.info[k] is not None:
            parts.append(f"{k}={rec.info[k]}")
    info = ";".join(parts) or "."
    return f"{rec.chrom} {rec.pos} . {rec.ref} {rec.alt} . . {info}"


def _default_post(url: str, body: bytes, timeout: float) -> tuple[int, bytes, dict[str, str]]:
    """POST JSON; return (status, body, lowercase headers). HTTP errors are returned, not raised."""
    req = urllib.request.Request(
        url, data=body, method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read(), {k.lower(): v for k, v in resp.headers.items()}
    except urllib.error.HTTPError as e:
        return e.code, e.read(), {k.lower(): v for k, v in e.headers.items()}


PostFn = Callable[[str, bytes, float], tuple[int, bytes, dict[str, str]]]


def query_vep(
    variants: Sequence[str],
    *,
    url: str = VEP_REST_URL,
    batch_size: int = VEP_MAX_BATCH,
    timeout: float = 180.0,
    max_retries: int = 5,
    post: PostFn = _default_post,
    sleep: Callable[[float], None] = time.sleep,
) -> list[dict]:
    """POST variants to VEP in batches of <= ``batch_size``; return the concatenated JSON records.

    Honours HTTP 429 ``Retry-After``. ``post``/``sleep`` are injectable so
    unit tests run offline.
    """
    if batch_size < 1 or batch_size > VEP_MAX_BATCH:
        raise ValueError(f"batch_size must be in 1..{VEP_MAX_BATCH}, got {batch_size}")
    out: list[dict] = []
    for start in range(0, len(variants), batch_size):
        chunk = list(variants[start : start + batch_size])
        body = json.dumps({"variants": chunk}).encode()
        for attempt in range(max_retries + 1):
            status, payload, headers = post(url, body, timeout)
            if status == 200:
                data = json.loads(payload)
                if not isinstance(data, list):
                    raise RuntimeError(f"VEP returned non-list JSON for batch at {start}")
                out.extend(data)
                break
            if status == 429 and attempt < max_retries:
                sleep(float(headers.get("retry-after", "1")))
                continue
            raise RuntimeError(f"VEP HTTP {status} for batch at {start}: {payload[:200]!r}")
    return out


# --------------------------------------------------------------------------- response parsing


def most_severe_term(terms: Iterable[str]) -> str:
    """Pick the most severe SO term per Ensembl's ranking; unknown terms rank last."""
    terms = list(terms)
    if not terms:
        return "intergenic_variant"
    return min(terms, key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))


def summarize_vep_record(rec: dict) -> GeneAnnotation:
    """Collapse one VEP JSON record to gene symbols / IDs / consequence summary.

    Only ``transcript_consequences`` contribute genes; regulatory and
    intergenic blocks are folded into the term list only.
    """
    symbols: list[str] = []
    gene_ids: list[str] = []
    terms: list[str] = []
    txs = rec.get("transcript_consequences") or []
    for t in txs:
        sym = t.get("gene_symbol")
        gid = t.get("gene_id")
        if sym and sym not in symbols:
            symbols.append(sym)
        if gid and gid not in gene_ids:
            gene_ids.append(gid)
        for c in t.get("consequence_terms", []):
            if c not in terms:
                terms.append(c)
    for block in ("regulatory_feature_consequences", "motif_feature_consequences", "intergenic_consequences"):
        for t in rec.get(block) or []:
            for c in t.get("consequence_terms", []):
                if c not in terms:
                    terms.append(c)
    severe = rec.get("most_severe_consequence") or most_severe_term(terms)
    return GeneAnnotation(
        gene_symbols=tuple(symbols),
        gene_ids=tuple(gene_ids),
        most_severe=severe,
        terms=tuple(sorted(terms, key=lambda t: _SEVERITY_RANK.get(t, len(CONSEQUENCE_SEVERITY)))),
        n_transcripts=len(txs),
    )


def _record_key(rec: VcfRecord) -> tuple[str, int, str]:
    return (rec.chrom, rec.pos, rec.alt)


def _vep_input_key(vep_rec: dict) -> tuple[str, int, str] | None:
    """Recover (chrom, pos, alt) from the echoed ``input`` string."""
    parts = str(vep_rec.get("input", "")).split()
    if len(parts) < 5:
        return None
    try:
        return (parts[0], int(parts[1]), parts[4])
    except ValueError:
        return None


# --------------------------------------------------------------------------- driver


@dataclass
class AnnotationStats:
    n_records: int
    n_annotated: int
    n_with_gene: int
    n_unmatched: int
    seconds: float


def annotate_records(
    records: Sequence[VcfRecord],
    *,
    keep_svlen: bool = False,
    batch_size: int = VEP_MAX_BATCH,
    post: PostFn = _default_post,
    sleep: Callable[[float], None] = time.sleep,
) -> tuple[list[VcfRecord], AnnotationStats]:
    """Return records with GENE/GENEID/CSQ/CSQ_TERMS/CSQ_NTX added to INFO.

    Records VEP did not echo back are left unannotated and counted in
    ``n_unmatched`` rather than silently dropped.
    """
    t0 = time.time()
    variants = [vep_region_string(r, keep_svlen=keep_svlen) for r in records]
    responses = query_vep(variants, batch_size=batch_size, post=post, sleep=sleep)
    by_key: dict[tuple[str, int, str], GeneAnnotation] = {}
    for resp in responses:
        key = _vep_input_key(resp)
        if key is not None:
            by_key[key] = summarize_vep_record(resp)

    out: list[VcfRecord] = []
    n_annot = n_gene = n_unmatched = 0
    for rec in records:
        ann = by_key.get(_record_key(rec))
        new_info = dict(rec.info)
        if ann is None:
            n_unmatched += 1
        else:
            new_info.update(ann.info_fields())
            n_annot += 1
            n_gene += bool(ann.gene_symbols)
        out.append(VcfRecord(rec.chrom, rec.pos, rec.id, rec.ref, rec.alt, rec.qual, rec.filter, new_info, list(rec.rest)))
    return out, AnnotationStats(len(records), n_annot, n_gene, n_unmatched, time.time() - t0)


def annotate_vcf(
    in_path: str | Path,
    out_path: str | Path,
    *,
    tsv_path: str | Path | None = None,
    keep_svlen: bool = False,
    batch_size: int = VEP_MAX_BATCH,
    post: PostFn = _default_post,
    sleep: Callable[[float], None] = time.sleep,
) -> AnnotationStats:
    """Read ``in_path``, annotate via VEP, write ``out_path`` (and optional flat TSV)."""
    headers, records = read_vcf(in_path)
    annotated, stats = annotate_records(records, keep_svlen=keep_svlen, batch_size=batch_size, post=post, sleep=sleep)
    write_vcf(out_path, add_info_headers(headers), annotated)
    if tsv_path is not None:
        write_annotation_tsv(tsv_path, annotated)
    return stats


# --------------------------------------------------------------------------- snpEff backend
#
# snpEff (Cingolani et al. 2012) runs locally against a downloaded Ensembl or
# RefSeq database and writes an ``ANN`` INFO field (one comma-separated entry
# per transcript, 16 pipe-separated columns). Measured 2026-09-18 on the same
# 30 chr22 calls with ``snpEff 5.4c -Xmx8g GRCh38.99``: 57 s wall, of which
# ~55 s is fixed database load and ~1 s annotation; 30/30 annotated; the
# 6 kb L1 at chr22:19223382 evaluated as a point (intron_variant, CLTCL1) with
# no SVLEN workaround. ``chr22`` input was mapped to database contig ``22``
# and echoed back as ``chr22``. Default JVM heap OOMs on GRCh38 -> pass -Xmx.
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


def parse_snpeff_ann(ann: str) -> GeneAnnotation:
    """Collapse a snpEff ``ANN`` value to the same GeneAnnotation shape as :func:`summarize_vep_record`.

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
    cmd = [snpeff_bin, f"-Xmx{xmx}", "-noStats"]
    if config is not None:
        cmd += ["-c", str(config)]
    cmd += [genome, str(in_path)]
    return cmd


def annotate_records_snpeff(
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
    a canned snpEff output. The raw ``ANN`` field is dropped from the output so
    both backends produce the same INFO contract; use snpEff directly if you
    need the full per-transcript table.
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


def annotate_vcf_snpeff(
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
    """snpEff counterpart of :func:`annotate_vcf`; identical output contract."""
    headers, records = read_vcf(in_path)
    annotated, stats = annotate_records_snpeff(records, headers, genome, snpeff_bin=snpeff_bin, config=config, xmx=xmx, run=run)
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
