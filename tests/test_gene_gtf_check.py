"""Tests for scripts/gene_gtf_check.py.

Covers the interval logic the GENCODE cross-check depends on (union of exons, the "insertion lies inside a
feature" rule, gene-body strand sets, intronic/exonic/intergenic classes), the GTF reader, the sense
summary, and one end-to-end run of ``main`` on a hand-built two-genome VCF directory and a tiny GTF.
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

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "gene_gtf_check.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("gene_gtf_check", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gene_gtf_check"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


ggc = _load_module()


def test_union_intervals_merges_overlapping_and_touching_intervals():
    x = pd.DataFrame({"start": [10, 1, 3, 5], "end": [12, 5, 8, 9]})
    assert ggc.union_intervals(x).tolist() == [[1, 9], [10, 12]]


def test_inside_requires_the_feature_to_cover_both_flanking_bases():
    iv = np.array([[10, 20]])
    pos = np.array([9, 10, 15, 19, 20, 25])
    assert ggc.inside(iv, pos).tolist() == [False, True, True, True, False, False]


def test_inside_is_false_for_positions_before_the_first_interval():
    iv = np.array([[100, 200], [300, 400]])
    assert ggc.inside(iv, np.array([1, 250, 350])).tolist() == [False, False, True]


def _write_gtf(path: Path, lines: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("##description: test\n")
        fh.write("\n".join(lines) + "\n")


def _gtf_line(chrom, feat, start, end, strand, gene_id):
    return f'{chrom}\tHAVANA\t{feat}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_id}"; gene_type "x";'


def test_load_gtf_keeps_genes_and_exons_and_strips_the_version(tmp_path):
    path = tmp_path / "t.gtf.gz"
    _write_gtf(path, [_gtf_line("chr1", "gene", 100, 500, "+", "ENSG1.7"),
                      _gtf_line("chr1", "exon", 100, 200, "+", "ENSG1.7"),
                      _gtf_line("chr1", "transcript", 100, 500, "+", "ENSG1.7"),
                      _gtf_line("chr1", "CDS", 120, 180, "+", "ENSG1.7")])
    out = ggc.load_gtf(str(path))
    assert sorted(out.feat) == ["exon", "gene"]
    assert set(out.gene_id) == {"ENSG1"}


def _gtf_frame():
    rows = [("chr1", "gene", 100, 500, "+", "A"), ("chr1", "exon", 100, 200, "+", "A"),
            ("chr1", "gene", 400, 900, "-", "B")]
    return pd.DataFrame(rows, columns=["chrom", "feat", "start", "end", "strand", "gene_id"])


def test_overlay_gtf_classifies_exonic_intronic_mixed_and_intergenic_sites():
    sites = pd.DataFrame({"site": [1, 2, 3, 4, 5, 6], "chrom": ["chr1"] * 5 + ["chr2"],
                          "pos": [150, 300, 450, 700, 5000, 100]})
    out = ggc.overlay_gtf(sites, _gtf_frame()).set_index("site")
    assert out.gtf_class.to_dict() == {1: "exonic", 2: "intronic", 3: "intronic", 4: "intronic",
                                        5: "intergenic", 6: "intergenic"}
    assert out.loc[2, "g_strands"] == frozenset({"+"})
    assert out.loc[3, "g_strands"] == frozenset({"+", "-"})      # overlapping genes on opposite strands
    assert out.loc[4, "g_strands"] == frozenset({"-"})
    assert out.loc[6, "n_gene"] == 0                              # chromosome absent from the GTF


def test_overlay_gtf_treats_the_last_base_of_a_gene_as_outside():
    sites = pd.DataFrame({"site": [1, 2], "chrom": ["chr1", "chr1"], "pos": [499, 500]})
    gtf = _gtf_frame()
    gtf = gtf[gtf.gene_id == "A"]
    out = ggc.overlay_gtf(sites, gtf).set_index("site")
    assert out.loc[1, "n_gene"] == 1 and out.loc[2, "n_gene"] == 0


def test_sense_summary_reports_counts_fractions_and_the_strand_split():
    z = pd.DataFrame({
        "family": ["LINE1"] * 8,
        "gene_strand": ["+"] * 4 + ["-"] * 4,
        "orient": ["+", "-", "-", "-", "-", "+", "+", "+"],
    })
    out = ggc.sense_summary(z, "orient")["LINE1"]
    assert (out["n"], out["n_sense"]) == (8, 2)                   # one sense call in + genes, one in - genes
    assert out["plus_genes"] == [1, 4] and out["minus_genes"] == [1, 4]
    assert out["sense_frac"] == pytest.approx(2 / 8)


def _vcf(path: Path, rows: list[str]) -> None:
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(rows) + "\n")


def _row(pos, orient, csq, genes, strands, side):
    info = (f"MEIFAMILY=ALU;ORIENT={orient};CSQ={csq};CSQ_TERMS={csq};GENEID={genes};GENE_STRAND={strands};"
            f"SUPPORT=SR_L:1|SR_R:1|polyA_side:{side}")
    return f"chr1\t{pos}\t.\tN\t<INS:ME:ALU>\t.\tPASS\t{info}"


def test_main_end_to_end_on_a_tiny_callset(tmp_path, monkeypatch):
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir()
    # site 1 (pos 300): intronic gene A (+), ORIENT + and poly(A) side L (=> +): sense, both measures agree
    # site 2 (pos 700): intronic gene B (-), ORIENT + but side L (=> +): antisense, measures agree
    # site 3 (pos 5000): intergenic
    # site 4 (pos 150): exonic in gene A
    _vcf(vcf_dir / "S1_x.vcf", [_row(300, "+", "intron_variant", "A", "+", "L"),
                                _row(700, "+", "intron_variant", "B", "-", "L"),
                                _row(5000, "-", "intergenic_variant", "", "", "R"),
                                _row(150, "+", "non_coding_transcript_exon_variant", "A", "+", "L")])
    _vcf(vcf_dir / "S2_x.vcf", [_row(310, "+", "intron_variant", "A", "+", "R")])    # poly(A) side disagrees
    gtf = tmp_path / "t.gtf.gz"
    _write_gtf(gtf, [_gtf_line("chr1", "gene", 100, 500, "+", "A.1"), _gtf_line("chr1", "exon", 100, 200, "+", "A.1"),
                     _gtf_line("chr1", "gene", 600, 900, "-", "B.1")])
    monkeypatch.syspath_prepend(str(SCRIPT_PATH.parent))
    out_dir = tmp_path / "out"
    assert ggc.main([str(vcf_dir), str(gtf), "--out-dir", str(out_dir)]) == 0
    report = json.loads((out_dir / "gene_gtf_report.json").read_text())
    assert report["gene_ids_in_callset"] == 2 and report["gene_ids_found_in_gtf"] == 2
    assert report["strand_resolution"]["intronic_sites"] == 2
    assert report["strand_resolution"]["strand_agrees"] == 2
    assert report["analysed_sites"]["n"] == 2
    # site at pos 300/310 has members with the same ORIENT; poly(A) sides L and R disagree -> no poly(A) orientation
    assert report["sense_by_ORIENT"]["ALU"]["n"] == 2
    assert report["sense_by_polyA"]["ALU"]["n"] == 1
    assert report["orientation_agreement_calls"]["ALU"]["n"] == 5          # 4 agree, 1 (pos 310: ORIENT + vs side R) does not
    assert report["orientation_agreement_calls"]["ALU"]["agreement"] == pytest.approx(4 / 5)
