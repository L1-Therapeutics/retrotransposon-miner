"""Per-call, overlap-only gene annotation from a GENCODE GTF, compared with the snpEff gene fields.

For every call in a directory of snpEff-annotated MEI VCFs this writes one row with:
  * the genes whose body contains the insertion point (not the +/-5 kb snpEff window),
    with symbol, biotype and strand, and the exon/intron/intergenic class;
  * the snpEff GENEID list, and whether it contains genes the breakpoint does not overlap
    (``snpeff_neighbour_genes``);
  * whether GENE and GENEID differ in length (``register_mismatch``) and a ``proposed_GENE``
    rebuilt from the GTF symbols in GENEID order, so it is index-aligned with GENEID and GENE_STRAND.

An insertion after POS lies inside a feature when start <= POS < end (1-based, inclusive GTF coordinates).

Usage:
    python gene_overlap_table.py VCF_DIR GENCODE.annotation.gtf.gz --out-dir OUT
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_gtf_check import inside, union_intervals  # noqa: E402
from gene_orientation_summary import load_calls  # noqa: E402


def load_gene_model(path: str) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    cols = ["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"]
    g = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols, dtype={"chrom": str, "attr": str},
                    usecols=[0, 2, 3, 4, 6, 8], compression="infer")
    g = g[g.feat.isin(("gene", "exon"))].copy()
    g["gene_id"] = g.attr.str.extract(r'gene_id "([^"]+)"')[0].str.split(".").str[0]
    g["gene_name"] = g.attr.str.extract(r'gene_name "([^"]+)"')[0]
    g["gene_type"] = g.attr.str.extract(r'gene_type "([^"]+)"')[0]
    genes = g[g.feat == "gene"].drop(columns=["attr", "feat"]).reset_index(drop=True)
    exons = {c: union_intervals(x) for c, x in g[g.feat == "exon"].groupby("chrom")}
    return genes, exons


def overlap_annotation(positions: pd.DataFrame, genes: pd.DataFrame, exons: dict[str, np.ndarray]) -> pd.DataFrame:
    """One row per (chrom, pos): overlapping gene ids/names/types/strands and the class."""
    by_chrom = {c: x for c, x in genes.groupby("chrom")}
    rows = []
    for chrom, x in positions.groupby("chrom"):
        pos = x.pos.values
        body = by_chrom.get(chrom)
        if body is None:
            rows += [(chrom, int(p), [], [], [], [], False) for p in pos]
            continue
        hit = (body.start.values[None, :] <= pos[:, None]) & (pos[:, None] < body.end.values[None, :])
        in_exon = inside(exons[chrom], pos) if chrom in exons else np.zeros(len(pos), bool)
        for k, p in enumerate(pos):
            idx = np.where(hit[k])[0]
            rows.append((chrom, int(p), body.gene_id.values[idx].tolist(),
                         [n if isinstance(n, str) else i for n, i in zip(body.gene_name.values[idx], body.gene_id.values[idx])],
                         body.gene_type.values[idx].tolist(), body.strand.values[idx].tolist(), bool(in_exon[k])))
    out = pd.DataFrame(rows, columns=["chrom", "pos", "ov_ids", "ov_names", "ov_types", "ov_strands", "in_exon"])
    out["gtf_class"] = np.where(out.ov_ids.map(len) == 0, "intergenic", np.where(out.in_exon, "exonic", "intronic"))
    return out


def build_table(calls: pd.DataFrame, genes: pd.DataFrame, exons: dict[str, np.ndarray]) -> pd.DataFrame:
    ov = overlap_annotation(calls[["chrom", "pos"]].drop_duplicates(), genes, exons)
    t = calls.merge(ov, on=["chrom", "pos"], how="left")
    name_of = dict(zip(genes.gene_id, genes.gene_name.where(genes.gene_name.notna(), genes.gene_id)))
    t["register_mismatch"] = t.gene_l.map(len) != t.geneid_l.map(len)
    t["snpeff_neighbour_genes"] = [[g for g in ids if g not in set(ov_ids)] for ids, ov_ids in zip(t.geneid_l, t.ov_ids)]
    t["proposed_GENE"] = [[name_of.get(g, g) for g in ids] for ids in t.geneid_l]
    t["overlap_strand_class"] = t.ov_strands.map(lambda s: "none" if not s else ("single" if len(set(s)) == 1 else "mixed"))
    cols = {"sample": t["sample"], "chrom": t.chrom, "pos": t.pos, "family": t.MEIFAMILY, "orient": t.ORIENT,
            "snpeff_csq": t.CSQ, "gtf_class": t.gtf_class,
            "overlap_geneids": t.ov_ids.map(",".join), "overlap_genes": t.ov_names.map(",".join),
            "overlap_biotypes": t.ov_types.map(",".join), "overlap_strands": t.ov_strands.map(",".join),
            "overlap_strand_class": t.overlap_strand_class,
            "snpeff_geneids": t.geneid_l.map(",".join), "snpeff_GENE": t.gene_l.map(",".join),
            "snpeff_neighbour_genes": t.snpeff_neighbour_genes.map(",".join),
            "register_mismatch": t.register_mismatch, "proposed_GENE": t.proposed_GENE.map(",".join)}
    return pd.DataFrame(cols)


def summarise(t: pd.DataFrame) -> dict:
    has_snp = t.snpeff_geneids != ""
    return {
        "calls": int(len(t)),
        "gtf_class": t.gtf_class.value_counts().to_dict(),
        "snpeff_gene_list_nonempty": int(has_snp.sum()),
        "calls_with_neighbour_genes_in_snpeff_list": int((t.snpeff_neighbour_genes != "").sum()),
        "register_mismatch_calls": int(t.register_mismatch.sum()),
        "overlap_strand_class": t.overlap_strand_class.value_counts().to_dict(),
        "intronic_calls_single_strand_by_overlap": int(((t.gtf_class == "intronic") & (t.overlap_strand_class == "single")).sum()),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("gtf")
    ap.add_argument("--out-dir", default="gene_overlap_out")
    args = ap.parse_args(argv)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes, exons = load_gene_model(args.gtf)
    table = build_table(load_calls(args.vcf_dir), genes, exons)
    table.to_csv(out / "gene_overlap_per_call.tsv", sep="\t", index=False)
    table[table.register_mismatch].to_csv(out / "register_mismatch_calls.tsv", sep="\t", index=False)
    summary = summarise(table)
    (out / "gene_overlap_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
