"""Tests for scripts/gene_overlap_table.py.

Covers the overlap rule (an insertion after POS lies inside a feature when start <= POS < end), the
exonic/intronic/intergenic class, the neighbour-gene flag, the register-mismatch flag with a proposed
GENE rebuilt from GENEID order, and one end-to-end run on a hand-built callset and GTF.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "gene_overlap_table.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("gene_overlap_table", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gene_overlap_table"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


got = _load_module()


def _genes():
    rows = [("chr1", 100, 500, "+", "A", "GENEA", "protein_coding"),
            ("chr1", 400, 900, "-", "B", "GENEB", "lncRNA"),
            ("chr1", 2000, 2500, "+", "C", "GENEA", "protein_coding")]       # same symbol as A, different gene ID
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])


def _exons():
    return {"chr1": np.array([[100, 200]])}


def test_overlap_annotation_classes_and_strands():
    pos = pd.DataFrame({"chrom": ["chr1"] * 5 + ["chr2"], "pos": [150, 300, 450, 700, 1500, 100]})
    out = got.overlap_annotation(pos, _genes(), _exons()).set_index(["chrom", "pos"])
    assert out.loc[("chr1", 150), "gtf_class"] == "exonic"
    assert out.loc[("chr1", 300), "gtf_class"] == "intronic" and out.loc[("chr1", 300), "ov_strands"] == ["+"]
    assert sorted(out.loc[("chr1", 450), "ov_strands"]) == ["+", "-"]       # overlapping genes on opposite strands
    assert out.loc[("chr1", 700), "ov_ids"] == ["B"]
    assert out.loc[("chr1", 1500), "gtf_class"] == "intergenic"
    assert out.loc[("chr2", 100), "gtf_class"] == "intergenic"               # chromosome absent from the GTF


def test_the_last_base_of_a_gene_is_outside_it():
    pos = pd.DataFrame({"chrom": ["chr1", "chr1"], "pos": [499, 500]})
    out = got.overlap_annotation(pos, _genes(), _exons()).set_index("pos")
    assert out.loc[499, "ov_ids"] == ["A", "B"] and out.loc[500, "ov_ids"] == ["B"]


def _calls(rows):
    base = {"sample": "S1", "chrom": "chr1", "MEIFAMILY": "ALU", "ORIENT": "+", "CSQ": "intron_variant",
            "gene_l": [], "geneid_l": []}
    return pd.DataFrame([{**base, **r} for r in rows])


def test_build_table_flags_neighbour_genes_and_leaves_overlap_genes_alone():
    calls = _calls([{"pos": 300, "gene_l": ["GENEA", "GENEB"], "geneid_l": ["A", "B"]}])   # B starts at 400
    t = got.build_table(calls, _genes(), _exons()).iloc[0]
    assert t.overlap_geneids == "A" and t.snpeff_neighbour_genes == "B"
    assert t.overlap_strand_class == "single"


def test_build_table_rebuilds_an_index_aligned_gene_list_when_symbols_repeat():
    # snpEff collapsed the repeated symbol: 3 gene IDs but only 2 symbols (the real chr6:25124585 pattern)
    calls = _calls([{"pos": 2200, "gene_l": ["GENEA", "ENSX"], "geneid_l": ["A", "C", "ENSX"]}])
    t = got.build_table(calls, _genes(), _exons()).iloc[0]
    assert bool(t.register_mismatch)
    assert t.proposed_GENE == "GENEA,GENEA,ENSX"                    # unknown ID falls back to itself
    assert len(t.proposed_GENE.split(",")) == len(t.snpeff_geneids.split(","))


def test_build_table_marks_mixed_strand_overlap_and_intergenic_calls():
    calls = _calls([{"pos": 450, "gene_l": ["GENEA", "GENEB"], "geneid_l": ["A", "B"]},
                    {"pos": 1500, "CSQ": "intergenic_variant"}])
    t = got.build_table(calls, _genes(), _exons()).set_index("pos")
    assert t.loc[450, "overlap_strand_class"] == "mixed" and t.loc[450, "snpeff_neighbour_genes"] == ""
    assert t.loc[1500, "gtf_class"] == "intergenic" and t.loc[1500, "overlap_strand_class"] == "none"


def _write_gtf(path: Path) -> None:
    def line(feat, s, e, strand, gid, name):
        return f'chr1\tHAVANA\t{feat}\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{gid}.1"; gene_type "protein_coding"; gene_name "{name}";'
    with gzip.open(path, "wt") as fh:
        fh.write("##test\n" + "\n".join([line("gene", 100, 500, "+", "A", "GENEA"), line("exon", 100, 200, "+", "A", "GENEA"),
                                          line("gene", 2000, 2500, "+", "C", "GENEA")]) + "\n")


def test_load_gene_model_strips_versions_and_keeps_symbols(tmp_path):
    gtf = tmp_path / "t.gtf.gz"
    _write_gtf(gtf)
    genes, exons = got.load_gene_model(str(gtf))
    assert genes.gene_id.tolist() == ["A", "C"] and genes.gene_name.tolist() == ["GENEA", "GENEA"]
    assert exons["chr1"].tolist() == [[100, 200]]


def test_main_writes_the_per_call_table_the_mismatch_table_and_a_summary(tmp_path, monkeypatch):
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir()
    (vcf_dir / "S1_x.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t300\t.\tN\t<INS:ME:ALU>\t.\tPASS\tMEIFAMILY=ALU;ORIENT=+;CSQ=intron_variant;GENE=GENEA;GENEID=A;GENE_STRAND=+\n"
        "chr1\t2200\t.\tN\t<INS:ME:ALU>\t.\tPASS\tMEIFAMILY=ALU;ORIENT=-;CSQ=intron_variant;GENE=GENEA;GENEID=A,C;GENE_STRAND=+,+\n")
    gtf = tmp_path / "t.gtf.gz"
    _write_gtf(gtf)
    monkeypatch.syspath_prepend(str(SCRIPT_PATH.parent))
    out = tmp_path / "out"
    assert got.main([str(vcf_dir), str(gtf), "--out-dir", str(out)]) == 0
    table = pd.read_csv(out / "gene_overlap_per_call.tsv", sep="\t", keep_default_na=False)
    assert len(table) == 2 and table.register_mismatch.tolist() == [False, True]
    assert (out / "register_mismatch_calls.tsv").read_text().count("\n") == 2     # header + one row
    summary = json.loads((out / "gene_overlap_summary.json").read_text())
    assert summary["calls"] == 2 and summary["register_mismatch_calls"] == 1
