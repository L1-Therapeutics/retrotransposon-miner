"""Annotate an MEI VCF with each insertion's position relative to the nearest exon.

The pipeline's gene annotation records WHICH gene an insertion lands in, never
WHERE in that gene relative to the splicing machinery. Exonization -- the route
by which an intronic Alu becomes a new exon and causes Alport syndrome, Sly
syndrome or OAT deficiency -- depends on an antisense element sitting just
upstream of an exon (Gal-Mark 2008; Lev-Maor 2008). These two fields record that
position. They make no prediction about splicing; they record where the element
is so a later analysis can ask.

snpEff cannot supply this. It annotates the breakpoint in reference coordinates
and never reads the inserted sequence, and its splice-region window is ~3 bp
exonic / 3-8 bp intronic -- which captures 1 of 2,985 intronic Alu sites in our
ten-genome callset.

Usage:
  python exon_distance_annotate.py in.vcf gencode.gtf.gz out.vcf
"""
from __future__ import annotations

import argparse
import gzip
from collections import defaultdict

import numpy as np

EXON_DIST_HEADER = (
    '##INFO=<ID=EXON_DIST,Number=1,Type=Integer,Description="Distance in bp from the insertion '
    'breakpoint to the nearest annotated exon boundary of the overlapping gene. 0 when the '
    'breakpoint is inside an exon; -1 when it is not inside exactly one single-strand gene">'
)
EXON_SIDE_HEADER = (
    '##INFO=<ID=EXON_SIDE,Number=1,Type=String,Description="Whether the breakpoint lies upstream '
    'or downstream of that nearest exon in the direction of transcription. \'.\' when EXON_DIST '
    'is 0 or -1">'
)
INFO_HEADERS = (EXON_DIST_HEADER, EXON_SIDE_HEADER)
INFO_KEYS = ("EXON_DIST", "EXON_SIDE")

NO_GENE = -1
IN_EXON = 0


def load_gencode(path: str) -> tuple[dict, dict]:
    """Per-chromosome gene records and merged exon intervals, 0-based half-open."""
    genes: dict[str, list] = defaultdict(list)
    exon_rows: dict[str, list] = defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if f[2] == "gene":
                if strand in ("+", "-"):
                    genes[chrom].append((start, end, strand))
            else:
                exon_rows[chrom].append((start, end))
    exons = {c: merge_intervals(np.array(sorted(r), dtype=np.int64)) for c, r in exon_rows.items()}
    return dict(genes), exons


def merge_intervals(a: np.ndarray) -> np.ndarray:
    """Merge overlapping or touching [start, end) rows. Input sorted by start."""
    if len(a) == 0:
        return a
    out = [a[0].copy()]
    for s, e in a[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append(np.array([s, e]))
    return np.array(out, dtype=np.int64)


def gene_strand_at(gene_list, pos: int):
    """Strand of the gene(s) covering pos; None when none or mixed strands."""
    hits = {st for s, e, st in gene_list if s <= pos < e}
    return hits.pop() if len(hits) == 1 else None


def exon_distance(exons: np.ndarray, pos: int, strand: str) -> tuple[int, str]:
    """(distance, side) for pos against merged exon intervals of its gene."""
    if len(exons) == 0:
        return NO_GENE, "."
    i = int(np.searchsorted(exons[:, 0], pos, side="right")) - 1
    if 0 <= i < len(exons) and exons[i][0] <= pos < exons[i][1]:
        return IN_EXON, "."
    best_d, best_side = None, "."
    for j in (i, i + 1):
        if not (0 <= j < len(exons)):
            continue
        s, e = exons[j]
        if pos < s:
            d, genomic_left = int(s - pos), True
        else:
            d, genomic_left = int(pos - (e - 1)), False
        if best_d is None or d < best_d:
            best_d = d
            if strand == "+":
                best_side = "upstream" if genomic_left else "downstream"
            else:
                best_side = "downstream" if genomic_left else "upstream"
    if best_d is None:
        return NO_GENE, "."
    return best_d, best_side


def annotate_info(info: str, dist: int, side: str) -> str:
    """Append the two keys, replacing any value already present."""
    parts = [p for p in info.split(";") if p and p.split("=", 1)[0] not in INFO_KEYS]
    parts.append(f"EXON_DIST={dist}")
    parts.append(f"EXON_SIDE={side}")
    return ";".join(parts)


def insert_headers(headers: list[str]) -> list[str]:
    """Add this module's ##INFO lines, skipping ones already present."""
    present = {
        h.split("ID=", 1)[1].split(",", 1)[0]
        for h in headers
        if h.startswith("##INFO=<ID=")
    }
    return [h for h in INFO_HEADERS if h.split("ID=", 1)[1].split(",", 1)[0] not in present]


def annotate_vcf(vcf_in: str, gtf: str, vcf_out: str) -> dict:
    genes, exons = load_gencode(gtf)
    stats = {"records": 0, "in_gene": 0, "in_exon": 0, "no_gene": 0}
    headers: list[str] = []
    out_lines: list[str] = []
    with open(vcf_in) as fh:
        for line in fh:
            if line.startswith("##"):
                headers.append(line.rstrip("\n"))
                out_lines.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                out_lines.extend(insert_headers(headers))
                out_lines.append(line.rstrip("\n"))
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                out_lines.append(line.rstrip("\n"))
                continue
            stats["records"] += 1
            chrom, pos = f[0], int(f[1])
            strand = gene_strand_at(genes.get(chrom, []), pos)
            if strand is None:
                dist, side = NO_GENE, "."
                stats["no_gene"] += 1
            else:
                dist, side = exon_distance(
                    exons.get(chrom, np.empty((0, 2), np.int64)), pos, strand
                )
                stats["in_gene"] += 1
                if dist == IN_EXON:
                    stats["in_exon"] += 1
            f[7] = annotate_info(f[7], dist, side)
            out_lines.append("\t".join(f))
    with open(vcf_out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    return stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf_in")
    ap.add_argument("gtf")
    ap.add_argument("vcf_out")
    a = ap.parse_args(argv)
    s = annotate_vcf(a.vcf_in, a.gtf, a.vcf_out)
    print(
        f"[exon-distance] records={s['records']} in_gene={s['in_gene']} "
        f"in_exon={s['in_exon']} no_gene={s['no_gene']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
